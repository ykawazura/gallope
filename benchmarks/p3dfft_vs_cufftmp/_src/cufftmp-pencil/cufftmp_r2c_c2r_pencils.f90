!
! This sample illustrates the use of cufftXtSetDistribution
! to support arbitrary user data distributions in the case of an R2C - C2R transform
! 
! It performs
! - Forward FFT
! - Element wise kernel (normalization)
! - Inverse FFT
!
! Modified to match cufftmp_r2c.f90 style:
! - Direct nx, ny, nz, nrepeat parameter setting
! - Warmup loop
! - Timing measurement
! - Result validation

module cufft_required
    integer :: planr2c, planc2r
    integer :: local_rshape_in(3), local_rshape_out(3), local_cshape_out(3)

end module cufft_required


program cufftmp_r2c_c2r_pencils
    use iso_c_binding
    use cudafor
    use cufftXt
    use cufft
    use openacc
    use mpi
    use cufft_required
    implicit none

    integer :: size, rank, ndevices, ierr
    integer :: nrepeat, nx, ny, nz ! nx slowest
    integer :: params(4), iounit
    logical :: have_params
    integer :: nz_complex, nz_real_padded
    integer :: i, j, k, m, nranks1d
    integer :: my_nx, my_ny, my_nz, ranks_cutoff, whichgpu(1)
    integer :: local_nx, local_ny
    real(8), dimension(:, :, :), allocatable :: u, ref
    complex(8), dimension(:,:,:), allocatable :: u_permuted
    real(8) :: max_norm, max_diff, max_norm_global, max_diff_global
    real(8) :: Nglob, factor, prec
    real(8) :: rtime1, rtime2, rtime0

    integer :: hostname_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
    character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: all_hostnames(:)
    integer :: unique_nodes
    logical :: is_unique

    ! cufft-related
    integer(c_size_t) :: worksize(1)
    type(cudaLibXtDesc), pointer :: u_desc
    type(cudaXtDesc), pointer    :: u_descptr
    complex(8), pointer, device     :: u_dptr(:,:,:)
    integer(kind=cuda_stream_kind) :: stream

    type :: Box3D
        integer(c_long_long), dimension(0:2) :: lower
        integer(c_long_long), dimension(0:2) :: upper
        integer(c_long_long), dimension(0:2) :: strides
    end type Box3D
    type(Box3D), dimension(:), allocatable :: input_boxes, output_boxes

    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD,size,ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)

    call checkCuda(cudaGetDeviceCount(ndevices))
    nranks1d = nint(sqrt(real(size)))
    if (nranks1d * nranks1d .ne. size) then
        if (rank == 0) print*,"The number of MPI ranks should be a perfect square"
        call mpi_finalize(ierr)
        stop
    end if  

    call checkCuda(cudaSetDevice(mod(rank, ndevices)))
    whichgpu(1) = mod(rank, ndevices)

!vvvvv   Setting parameters   vvvvv
    nx      = 2048
    ny      = 2048
    nz      = 2048
    nrepeat = 10
    ! A file "params.in" in the run directory, holding "nx ny nz nrepeat" on a
    ! single line, overrides the defaults above.  The weak-scaling series changes
    ! the grid at every point, and without this the benchmark would need one
    ! binary per shape.  Absent the file the run is bit-for-bit the original one.
    if (rank == 0) then
        inquire(file='params.in', exist=have_params)
        if (have_params) then
            open(newunit=iounit, file='params.in', status='old', action='read')
            read(iounit, *) nx, ny, nz, nrepeat
            close(iounit)
        end if
    end if
    params = [nx, ny, nz, nrepeat]
    call mpi_bcast(params, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    nx = params(1); ny = params(2); nz = params(3); nrepeat = params(4)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! Validate that nx and ny are divisible by nranks1d
    if (mod(nx, nranks1d) /= 0 .or. mod(ny, nranks1d) /= 0) then
        if (rank == 0) then
            print*, "Error: nx and ny must be divisible by nranks1d (sqrt of nproc)"
            print*, "nx =", nx, ", ny =", ny, ", nranks1d =", nranks1d
        endif
        call mpi_finalize(ierr)
        stop
    endif

    ! Calculate local sizes
    local_nx = nx / nranks1d
    local_ny = ny / nranks1d
    nz_complex = nz/2+1
    nz_real_padded = 2 * nz_complex

    ! Get number of nodes and number of total GPUs
    if (rank.eq.0) then
        print *, 'cuFFTMp pencil test, random input'
        write(*, '(" nx = ", I5, ", ny = ", I5, ", nz = ", I5, ", nproc = ", I6, ", ndevices = ", I5, ", nrepeat = ", I5)') &
              nx, ny, nz, size, ndevices, nrepeat
    endif

    if (rank == 0) then
        allocate(all_hostnames(size))
    endif
   
    call MPI_GET_PROCESSOR_NAME(hostname, hostname_len, ierr)
    call MPI_Gather(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                    all_hostnames, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                    0, MPI_COMM_WORLD, ierr)
   
    if (rank == 0) then
        ! Count unique node number
        unique_nodes = 1
        do i = 2, size
            is_unique = .true.
            do j = 1, i-1
                if (trim(all_hostnames(i)) == trim(all_hostnames(j))) then
                    is_unique = .false.
                    exit
                endif
            enddo
            if (is_unique) unique_nodes = unique_nodes + 1
        enddo
      
        write(*, '(" ==========================================")') 
        write(*, '(" Total MPI processes: ", I6)') size
        write(*, '(" Unique nodes:        ", I6)') unique_nodes
        write(*, '(" GPUs per node:       ", I6)') ndevices
        write(*, '(" Total GPUs used:     ", I6)') unique_nodes * ndevices
        write(*, '(" nranks1d:            ", I6)') nranks1d
      
        if (unique_nodes < size) then
            print*, "WARNING: Multiple processes per node detected!"
            print*, "This may cause GPU contention."
        endif
      
        ! Show list of hostnames
        !do i = 1, size
        !   write(*, '(" Rank = ", I6, " on ", A)') i - 1, trim(all_hostnames(i))
        !enddo
        write(*, '(" ==========================================")') 
      
        deallocate(all_hostnames)
    endif

    ! Input data is real, distributed in X and Y, contiguous at Z
    ! Output data is complex, distributed in X and Y (picked arbitrarily)
    local_rshape_in = [nz_real_padded, local_ny, local_nx]
    local_cshape_out = [nz_complex, local_ny, local_nx]  
    local_rshape_out = [nz_real_padded, local_ny, local_nx]  

    if (rank == 0) then
        write(*,*) "local_rshape_in  :", local_rshape_in(1), local_rshape_in(2), local_rshape_in(3)
        write(*,*) "local_cshape_out :", local_cshape_out(1), local_cshape_out(2), local_cshape_out(3)
    end if

    ! cufftBox3d boxes are defined as { {x_start, y_start, z_start}, {x_end, y_end, z_end}, {x%strides, y%strides, z%strides} }
    ! where
    ! - {x/y/z}_{start/end} are the lower and upper corner of the boxes relatived to the global 3D box (of size nx * ny * nz)
    ! - {x/y/z} strides are the local strides. Note the use of the padded stride at z
    allocate(input_boxes(0:size-1), output_boxes(0:size-1) )
    do i = 0 , nranks1d - 1
        do j = 0, nranks1d - 1
            ! input data are pencils in X and Y, contiguous at Z
            input_boxes( i*nranks1d + j )%lower = [i*local_nx, j*local_ny, 0]
            input_boxes( i*nranks1d + j )%upper = [(i+1)*local_nx, (j+1)*local_ny, nz]
            input_boxes( i*nranks1d + j )%strides = [local_ny * nz_real_padded, nz_real_padded, 1]
            ! output data are pencils in X and Y, picked arbitrarily
            output_boxes( i*nranks1d + j )%lower = [i*local_nx, j*local_ny, 0]
            output_boxes( i*nranks1d + j )%upper = [(i+1)*local_nx, (j+1)*local_ny, nz_complex]
            output_boxes( i*nranks1d + j )%strides = [local_ny * nz_complex, nz_complex, 1]
        end do  
    end do

    ! Generate local, distributed data
    allocate(u(local_rshape_in(1), local_rshape_in(2), local_rshape_in(3)))
    allocate(u_permuted(local_cshape_out(1), local_cshape_out(2), local_cshape_out(3)))
    allocate(ref(local_rshape_in(1), local_rshape_in(2), local_rshape_in(3)))

    ! Initialize with random data
    call generate_random(nz, local_rshape_in(1), local_rshape_in(2), local_rshape_in(3), u)
    ref = u
    u_permuted = (0.0d0, 0.0d0)

    ! Calculate normalization factor
    Nglob = real(nx, 8) * real(ny, 8) * real(nz, 8)
    factor = 1.0d0 / Nglob

    ! Create cuFFT plans
    call checkCufft(cufftCreate(planr2c))
    call checkCufft(cufftCreate(planc2r))
    call checkCuda(cudaStreamCreate(stream))
    call checkCufft(cufftMpAttachComm(planr2c, CUFFT_COMM_MPI, MPI_COMM_WORLD), 'cufftMpAttachComm error')
    call checkCufft(cufftMpAttachComm(planc2r, CUFFT_COMM_MPI, MPI_COMM_WORLD), 'cufftMpAttachComm error')

    call checkCufft(cufftSetStream(planr2c, stream), 'cufftSetStream error')
    call checkCufft(cufftSetStream(planc2r, stream), 'cufftSetStream error')
    call checkCufft(cufftXtSetDistribution(planr2c, 3, input_boxes(rank)%lower, input_boxes(rank)%upper, &
                    output_boxes(rank)%lower, output_boxes(rank)%upper, &
                    input_boxes(rank)%strides, output_boxes(rank)%strides), 'cufftXtSetDistribution error')
    call checkCufft(cufftXtSetDistribution(planc2r, 3, input_boxes(rank)%lower, input_boxes(rank)%upper, &
                    output_boxes(rank)%lower, output_boxes(rank)%upper, &
                    input_boxes(rank)%strides, output_boxes(rank)%strides), 'cufftXtSetDistribution error')

    call checkCufft(cufftMakePlan3d(planr2c, nx, ny, nz, CUFFT_D2Z, worksize), 'cufftMakePlan3d r2c error')
    call checkCufft(cufftMakePlan3d(planc2r, nx, ny, nz, CUFFT_Z2D, worksize), 'cufftMakePlan3d c2r error')

    call checkCufft(cufftXtMalloc(planr2c, u_desc, CUFFT_XT_FORMAT_DISTRIBUTED_INPUT), 'cufftXtMalloc error')
    call cufft_memcpyH2D(u_desc, u, CUFFT_XT_FORMAT_INPLACE, .true.)

    ! Warming up
    do m = 1, 3
        ! Forward transform
        call checkCufft(cufftXtExecDescriptor(planr2c, u_desc, u_desc, CUFFT_FORWARD),'forward fft failed')
        call checkCuda(cudaStreamSynchronize(stream))

        ! Normalize - scale the output
        call c_f_pointer(u_desc%descriptor, u_descptr)
        call c_f_pointer(u_descptr%data(1), u_dptr, [local_cshape_out(1), local_cshape_out(2), local_cshape_out(3)])
        !$cuf kernel do (3)
        do k = 1, local_cshape_out(3)
            do j = 1, local_cshape_out(2)
                do i = 1, local_cshape_out(1)
                    u_dptr(i,j,k) = u_dptr(i,j,k) * factor
                end do
            end do
        end do
        call checkCuda(cudaDeviceSynchronize())

        ! Inverse transform
        call checkCufft(cufftXtExecDescriptor(planc2r, u_desc, u_desc, CUFFT_INVERSE), 'inverse fft failed')
        call checkCuda(cudaStreamSynchronize(stream))
    end do

    ! Main loop: repeat FFT transforms nrepeat times
    rtime1 = 0.0d0
    rtime0 = 0.0d0
    do m = 1, nrepeat
        ! Barrier for correct timing
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        ! Forward transform
        rtime1 = rtime1 - MPI_wtime()
        call checkCufft(cufftXtExecDescriptor(planr2c, u_desc, u_desc, CUFFT_FORWARD),'forward fft failed')
        call checkCuda(cudaStreamSynchronize(stream))
        rtime1 = rtime1 + MPI_wtime()

        ! Normalize - scale the output
        call c_f_pointer(u_desc%descriptor, u_descptr)
        call c_f_pointer(u_descptr%data(1), u_dptr, [local_cshape_out(1), local_cshape_out(2), local_cshape_out(3)])
        !$cuf kernel do (3)
        do k = 1, local_cshape_out(3)
            do j = 1, local_cshape_out(2)
                do i = 1, local_cshape_out(1)
                    u_dptr(i,j,k) = u_dptr(i,j,k) * factor
                end do
            end do
        end do
        call checkCuda(cudaDeviceSynchronize())

        ! Barrier for correct timing
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        ! Inverse transform
        rtime1 = rtime1 - MPI_wtime()
        call checkCufft(cufftXtExecDescriptor(planc2r, u_desc, u_desc, CUFFT_INVERSE), 'inverse fft failed')
        call checkCuda(cudaStreamSynchronize(stream))
        rtime1 = rtime1 + MPI_wtime()

        if(rank .eq. 0) then
            print *,'Iteration = ',m, ', time = ', rtime1 - rtime0
        endif
        rtime0 = rtime1

    end do

    ! Copy result back to CPU
    call cufft_memcpyD2H(u, u_desc, CUFFT_XT_FORMAT_INPLACE, .true.)

    ! Free resources
    call checkCufft(cufftXtFree(u_desc))
    call checkCufft(cufftDestroy(planr2c))
    call checkCufft(cufftDestroy(planc2r))
    call checkCuda(cudaStreamDestroy(stream))

    ! Check results (local)
    call checkNormDiff(nz, local_rshape_in(1), local_rshape_in(2), local_rshape_in(3), u, ref, max_norm, max_diff)

    ! Reduce to get global max
    call MPI_Reduce(max_norm, max_norm_global, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(max_diff, max_diff_global, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    if(rank .eq. 0) then
        prec = 1e-5
        if(max_diff_global .gt. prec * Nglob * 0.25) then
            print *,'Results are incorrect'
        else
            print *,'Results are correct'
        endif
        write (6,*) 'max diff =', max_diff_global
        write (6,*) 'Relative Linf = ', max_diff_global/max_norm_global
    endif

    ! Process timing statistics
    call MPI_Reduce(rtime1, rtime2, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    if (rank.eq.0) then
        write(6,*)'proc_id = ', rank, ', cpu time per loop = ', rtime2/dble(nrepeat)
    endif

    deallocate(u)
    deallocate(ref)
    deallocate(u_permuted)
    deallocate(input_boxes)
    deallocate(output_boxes)
    
    call mpi_finalize(ierr)

    
contains 
    subroutine checkCuda(istat, message)
        implicit none
        integer, intent(in)                   :: istat
        character(len=*),intent(in), optional :: message
        if (istat /= cudaSuccess) then
            write(*,"('Error code: ',I0, ': ')") istat
            write(*,*) cudaGetErrorString(istat)
            if(present(message)) write(*,*) message
            call mpi_finalize(ierr)
        endif
    end subroutine checkCuda

    subroutine checkCufft(istat, message)
        implicit none
        integer, intent(in)                   :: istat
        character(len=*),intent(in), optional :: message
        if (istat /= CUFFT_SUCCESS) then
            write(*,"('Error code: ',I0, ': ')") istat
            write(*,*) cudaGetErrorString(istat)
            if(present(message)) write(*,*) message
            call mpi_finalize(ierr)
        endif
    end subroutine checkCufft

    subroutine generate_random(nz1, nz, ny, nx, data)
        implicit none
        integer, intent(in) :: nx, ny, nz, nz1
        real(8), dimension(nz, ny, nx), intent(out) :: data
        real(8) :: rand(1)
        integer :: i,j,k
        do k = 1, nx
            do j = 1, ny
                do i = 1, nz1
                    call random_number(rand)
                    data(i,j,k) = rand(1)
                end do
            end do
        end do
        
    end subroutine generate_random

    subroutine checkNorm(nz1, nz, ny, nx, data, max_norm)
        implicit none
        integer, intent(in)  :: nx, ny, nz, nz1
        real(8), dimension(nz, ny, nx), intent(in) :: data
        real(8) :: max_norm
        integer :: i, j, k
        max_norm = 0
        do k = 1, nx
            do j = 1, ny
                do i = 1, nz1
                    max_norm = max(max_norm, abs(data(i,j,k)))
                end do
            end do
        end do
    end subroutine checkNorm

    subroutine checkNormComplex(nz, ny, nx, data, max_norm)
        implicit none
        integer, intent(in)  :: nx, ny, nz
        complex(8), dimension(nz, ny, nx), intent(in) :: data
        real(8) :: max_norm, max_diff
        integer :: i,j,k
        max_norm = 0
        do k = 1, nx
            do j = 1, ny
                do i = 1, nz
                    max_norm = max(max_norm, abs(data(i,j,k)%re))
                    max_norm = max(max_norm, abs(data(i,j,k)%im))
                end do
            end do
        end do
    end subroutine checkNormComplex

    subroutine checkNormDiff(nz1, nz, ny, nx, data, ref, max_norm, max_diff)
        implicit none
        integer, intent(in)  :: nx, ny, nz, nz1
        real(8), dimension(nz, ny, nx), intent(in) :: data, ref
        real(8) :: max_norm, max_diff
        integer :: i, j, k
        max_norm = 0
        max_diff = 0
        do k = 1, nx
            do j = 1, ny
                do i = 1, nz1
                    max_norm = max(max_norm, abs(data(i,j,k)))
                    max_diff = max(max_diff, abs(ref(i,j,k)-data(i,j,k)))
                end do
            end do
        end do
    end subroutine checkNormDiff


subroutine cufft_memcpyH2D(ulibxt, u_h, data_format, ismemcpy)
    implicit none
    type(cudaLibXtDesc), pointer, intent(out) :: ulibxt
    real(8), dimension(*), intent(in)          :: u_h
    integer, intent(in)                         :: data_format
    logical, intent(in)                         :: ismemcpy
    type(cudaXtDesc), pointer  :: uxt
    real(8), dimension(:,:,:), device, pointer :: u_d

    if(data_format == CUFFT_XT_FORMAT_INPLACE_SHUFFLED) then
      if (ismemcpy == .false.) then
        call checkCufft(cufftXtMemcpy(planc2r, ulibxt, u_h, CUFFT_COPY_HOST_TO_DEVICE), "cufft_memcpyHToD pinv Error")
      else
        call c_f_pointer(ulibxt%descriptor, uxt)
        call c_f_pointer(uxt%data(1), u_d, local_rshape_out)
        call checkCuda(cudaMemcpy(u_d, u_h, product(int(local_rshape_out,kind=8))), "cudamemcpy H2D Error")
        nullify(u_d, uxt)
      endif
    endif 

    if (data_format == CUFFT_XT_FORMAT_INPLACE) then
      if (ismemcpy == .false.) then
        call checkCufft(cufftXtMemcpy(planr2c, ulibxt, u_h, CUFFT_COPY_HOST_TO_DEVICE), "cufft_memcpyHToD pfor Error")
      else 
        call c_f_pointer(ulibxt%descriptor, uxt) 
        call c_f_pointer(uxt%data(1), u_d, local_rshape_in)
        call checkCuda(cudaMemcpy(u_d, u_h, product(int(local_rshape_in,kind=8))), "cudamemcpy H2D Error")
        nullify(u_d, uxt)
      endif
    endif 
end subroutine cufft_memcpyH2D


subroutine cufft_memcpyD2H(u_h, ulibxt, data_format,ismemcpy)
    implicit none
    type(cudaLibXtDesc), pointer, intent(in) :: ulibxt
    real(8), dimension(*), intent(out)      :: u_h
    integer, intent(in)                      :: data_format
    logical, intent(in)                      :: ismemcpy
    type(cudaXtDesc), pointer  :: uxt
    real(8), dimension(:,:,:), device, pointer :: u_d

    if(data_format == CUFFT_XT_FORMAT_INPLACE_SHUFFLED) then
      if (ismemcpy == .false.) then
        call checkCufft(cufftXtMemcpy(planr2c, u_h, ulibxt, CUFFT_COPY_DEVICE_TO_HOST), "cufft_memcpyDToH pfor Error")
      else
        call c_f_pointer(ulibxt%descriptor, uxt)
        call c_f_pointer(uxt%data(1), u_d, local_rshape_out)
        call checkCuda(cudaMemcpy(u_h, u_d, product(int(local_rshape_out,kind=8))), "cudamemcpy D2H Error")
        nullify(u_d, uxt)
      endif 
    endif
    
    if (data_format == CUFFT_XT_FORMAT_INPLACE) then
      if (ismemcpy == .false.) then
        call checkCufft(cufftXtMemcpy(planc2r, u_h, ulibxt, CUFFT_COPY_DEVICE_TO_HOST), "cufft_memcpyDToH pinv Error")
      else
        call c_f_pointer(ulibxt%descriptor, uxt)
        call c_f_pointer(uxt%data(1), u_d, local_rshape_in)
        call checkCufft(cudamemcpy(u_h, u_d, product(int(local_rshape_in,kind=8))), "cufft_memcpyD2H error")
        nullify(u_d, uxt)
      endif
    endif 
end subroutine cufft_memcpyD2H


end program cufftmp_r2c_c2r_pencils
