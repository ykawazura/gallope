!-----------------------------------------------!
!> @author  YK
!! @brief   CUDA and cuFFTmp setting
!-----------------------------------------------!
module cufftmp
  use iso_c_binding
  use cudafor
  use cufftXt
  use cufft
  use openacc
  use mpi
  implicit none

  public  init_cuda, init_cuFFTmp, finish_cuFFTmp
  public  ftran_r2c, btran_c2r
  public  ndevice, device_id

  integer :: ndevice, device_id
  integer :: planr2c, planc2r

  type(cudaLibXtDesc), pointer :: r2c_desc => null()
  type(cudaLibXtDesc), pointer :: c2r_desc => null()
  type(cudaLibXtDesc), pointer :: r2c_desc_tmp => null()
  type(cudaLibXtDesc), pointer :: c2r_desc_tmp => null()
contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of CUDA
!-----------------------------------------------!
  subroutine init_cuda
    use mp, only: proc_id, proc0, nproc
    implicit none

    integer :: hostname_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
    character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: all_hostnames(:)
    integer :: nnode
    logical :: is_unique
    integer, allocatable :: all_device_id(:)
    integer :: i, j, ierr
    integer :: node_comm, node_rank

    call check_cuda(cudaGetDeviceCount(ndevice))
    ! Select the GPU by intra-node rank so that several ranks sharing a node map
    ! to distinct devices. This stays machine-independent: it does not assume one
    ! GPU per node (Miyabi), and generalises to multi-GPU nodes via mpiprocs.
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                             MPI_INFO_NULL, node_comm, ierr)
    call MPI_Comm_rank(node_comm, node_rank, ierr)
    call MPI_Comm_free(node_comm, ierr)
    device_id = mod(node_rank, ndevice)
    call check_cuda(cudaSetDevice(device_id))

   ! Get number of nodes and number of total GPUs
   if (proc0) then
      allocate(all_hostnames(nproc))
      allocate(all_device_id(nproc))
   endif
   
   call MPI_GET_PROCESSOR_NAME(hostname, hostname_len, ierr)
   call MPI_Gather(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                   all_hostnames, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                   0, MPI_COMM_WORLD, ierr)
   call MPI_Gather(device_id, 1, MPI_INT, &
                   all_device_id, 1, MPI_INT, &
                   0, MPI_COMM_WORLD, ierr)
   
   if (proc0) then
      ! Count unique node number
      nnode = 1
      do i = 2, nproc
         is_unique = .true.
         do j = 1, i-1
            if (trim(all_hostnames(i)) == trim(all_hostnames(j))) then
               is_unique = .false.
               exit
            endif
         enddo
         if (is_unique) nnode = nnode + 1
      enddo
      
      write(*, '(" ==========================================")') 
      write(*, '(" MPI processes: ", I6)') nproc
      write(*, '(" Nodes:         ", I6)') nnode
      write(*, '(" GPUs per node: ", I6)') ndevice
      
      if (nproc > nnode * ndevice) then
         print*
         print*, "WARNING: # of processes is larger than # of GPUs!"
         print*
         stop
      endif
      
#ifdef DEBG
      ! Show list of hostnames
      do i = 1, nproc
        write(*, '(" proc_id = ", I6, ", gpu_id = ", I6, " on ", A)') i - 1, all_device_id(i), trim(all_hostnames(i))
      enddo
#endif
      write(*, '(" ==========================================")') 
      
      deallocate(all_hostnames)
      deallocate(all_device_id)
   endif
  end subroutine init_cuda


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of cuFFTmp
!-----------------------------------------------!
  subroutine init_cuFFTmp
    use grid, only: nlx, nly, nlz
    use mp, only: comm_fft
    implicit none
    integer(c_size_t) :: worksize(1)

    ! Attach each plan to comm_fft (not the world): the distributed 3D FFT runs
    ! independently within every FFT group. With P_m=P_s=1, comm_fft spans all
    ! ranks, so this is identical to the pre-refactor MPI_COMM_WORLD attach.
    call check_cufft(cufftCreate(planr2c))
    call check_cufft(cufftCreate(planc2r))
    call check_cufft(cufftMpAttachComm(planr2c, CUFFT_COMM_MPI, comm_fft), 'cufftMpAttachComm error')
    call check_cufft(cufftMpAttachComm(planc2r, CUFFT_COMM_MPI, comm_fft), 'cufftMpAttachComm error')

    call check_cufft(cufftMakePlan3d(planr2c, nlx, nly, nlz, CUFFT_D2Z, worksize), 'cufftMakePlan3d r2c error')
    call check_cufft(cufftMakePlan3d(planc2r, nlx, nly, nlz, CUFFT_Z2D, worksize), 'cufftMakePlan3d c2r error')

    call check_cufft(cufftXtMalloc(planr2c, r2c_desc, CUFFT_XT_FORMAT_INPLACE), 'r2c descriptor allocation failed')
    call check_cufft(cufftXtMalloc(planc2r, c2r_desc, CUFFT_XT_FORMAT_INPLACE_SHUFFLED), 'c2r descriptor allocation failed')
    call check_cufft(cufftXtMalloc(planr2c, r2c_desc_tmp, CUFFT_XT_FORMAT_INPLACE), 'r2c descriptor allocation failed')
    call check_cufft(cufftXtMalloc(planc2r, c2r_desc_tmp, CUFFT_XT_FORMAT_INPLACE_SHUFFLED), 'c2r descriptor allocation failed')

  end subroutine init_cuFFTmp


!-----------------------------------------------!
!> @author  YK
!! @brief   Forward FFT via cuFFTmp
!-----------------------------------------------!
  subroutine ftran_r2c(ur, uk)
    use grid, only: nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use time_stamp, only: put_time_stamp, timer_fft, timer_fft_prepost
    use mp, only: proc0
    implicit none
    real   (8), dimension(:,:,:), intent(in ) :: ur
    complex(8), dimension(:,:,:), intent(out) :: uk

    type(cudaXtDesc), pointer  :: uxt
    real   (8), dimension(:,:,:), device, pointer :: tmpr
    complex(8), dimension(:,:,:), device, pointer :: tmpk

    integer :: i, j, k, ierr

    if (proc0) call put_time_stamp(timer_fft_prepost)
    call c_f_pointer(r2c_desc%descriptor, uxt)
    call c_f_pointer(uxt%data(1), tmpr, [nlz_padded, nly, nlx_local])

    !$acc parallel loop collapse(3) present(ur) deviceptr(tmpr)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          tmpr(k,j,i) = ur(k,j,i)
        enddo
      enddo
    enddo
    if (proc0) call put_time_stamp(timer_fft_prepost)

    ! Forward transform
    if (proc0) call put_time_stamp(timer_fft)
    call check_cufft(cufftXtExecDescriptor(planr2c, r2c_desc, r2c_desc, CUFFT_FORWARD), 'forward fft failed')
    if (proc0) call put_time_stamp(timer_fft)
    ierr = cudaDeviceSynchronize()

    if (proc0) call put_time_stamp(timer_fft_prepost)
    call c_f_pointer(r2c_desc%descriptor, uxt)
    call c_f_pointer(uxt%data(1), tmpk, [nkz, nky_local, nkx])

    !$acc parallel loop collapse(3) present(uk) deviceptr(tmpk)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          uk(k,j,i) = tmpk(k,j,i)
        enddo
      enddo
    enddo
  
    r2c_desc = r2c_desc_tmp 

    nullify(uxt, tmpr, tmpk)
    if (proc0) call put_time_stamp(timer_fft_prepost)

  end subroutine ftran_r2c


!-----------------------------------------------!
!> @author  YK
!! @brief   Backward FFT via cuFFTmp
!-----------------------------------------------!
  subroutine btran_c2r(uk, ur)
    use grid, only: nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use time_stamp, only: put_time_stamp, timer_fft, timer_fft_prepost
    use mp, only: proc0
    implicit none
    complex(8), dimension(:,:,:), intent(in ) :: uk
    real   (8), dimension(:,:,:), intent(out) :: ur

    type(cudaXtDesc), pointer  :: uxt
    real   (8), dimension(:,:,:), device, pointer :: tmpr
    complex(8), dimension(:,:,:), device, pointer :: tmpk

    integer :: i, j, k, ierr

    if (proc0) call put_time_stamp(timer_fft_prepost)
    call c_f_pointer(c2r_desc%descriptor, uxt)
    call c_f_pointer(uxt%data(1), tmpk, [nkz, nky_local, nkx])

    !$acc parallel loop collapse(3) present(uk) deviceptr(tmpk)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          tmpk(k,j,i) = uk(k,j,i)
        enddo
      enddo
    enddo
    if (proc0) call put_time_stamp(timer_fft_prepost)

    ! Backward transform
    if (proc0) call put_time_stamp(timer_fft)
    call check_cufft(cufftXtExecDescriptor(planc2r, c2r_desc, c2r_desc, CUFFT_INVERSE), 'backward fft failed')
    if (proc0) call put_time_stamp(timer_fft)
    ierr = cudaDeviceSynchronize()

    if (proc0) call put_time_stamp(timer_fft_prepost)
    call c_f_pointer(c2r_desc%descriptor, uxt)
    call c_f_pointer(uxt%data(1), tmpr, [nlz_padded, nly, nlx_local])

    !$acc parallel loop collapse(3) present(ur) deviceptr(tmpr)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          ur(k,j,i) = tmpr(k,j,i)
        enddo
      enddo
    enddo
  
    c2r_desc = c2r_desc_tmp 
    if (proc0) call put_time_stamp(timer_fft_prepost)

    nullify(uxt, tmpr, tmpk)

  end subroutine btran_c2r


!-----------------------------------------------!
!> @author  YK
!! @brief   Finalization of cuFFTmp
!-----------------------------------------------!
  subroutine finish_cuFFTmp
    implicit none

    if (associated(r2c_desc)) then
      ! call check_cufft(cufftXtFree(r2c_desc))
      nullify(r2c_desc)
    endif
    if (associated(c2r_desc)) then
      call check_cufft(cufftXtFree(c2r_desc))
      nullify(c2r_desc)
    endif
    ! if (associated(r2c_desc_tmp)) then
    !   call check_cufft(cufftXtFree(r2c_desc_tmp))
    !   nullify(r2c_desc_tmp)
    ! endif
    ! if (associated(c2r_desc_tmp)) then
    !   call check_cufft(cufftXtFree(c2r_desc_tmp))
    !   nullify(c2r_desc_tmp)
    ! endif

    call check_cufft(cufftDestroy(planr2c))
    call check_cufft(cufftDestroy(planc2r))
  end subroutine finish_cuFFTmp


  subroutine check_cuda(istat, message)
    implicit none
    integer, intent(in)                   :: istat
    character(len=*),intent(in), optional :: message
    integer :: ierr

    if (istat /= cudaSuccess) then
      write(*,"('Error code: ',I0, ': ')") istat
      write(*,*) cudaGetErrorString(istat)
      if(present(message)) write(*,*) message
      call mpi_finalize(ierr)
    endif
  end subroutine check_cuda


  subroutine check_cufft(istat, message)
    implicit none
    integer, intent(in)                   :: istat
    character(len=*),intent(in), optional :: message
    integer :: ierr

    if (istat /= CUFFT_SUCCESS) then
      write(*,"('Error code: ',I0, ': ')") istat
      write(*,*) cudaGetErrorString(istat)
      if(present(message)) write(*,*) message
      call mpi_finalize(ierr)
    endif
  end subroutine check_cufft

end module cufftmp

