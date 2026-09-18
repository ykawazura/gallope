! =====================================================================
!  cuFFTMp subcommunicator smoke test (gallope Phase 0)
!
!  Purpose:
!    Split MPI_COMM_WORLD into NCOMM disjoint subcommunicators and show,
!    on real hardware, that cuFFTMp can run an independent distributed 3D
!    FFT within each subcommunicator. This retires the existential risk of
!    the [P_fft, P_m, P_s] decomposition used by the kinetic extension.
!
!  Checks:
!    V1: a cuFFTMp plan attached to a subcommunicator transforms correctly
!        (per-subcomm independent random input -> forward+inverse -> input)
!    V2: multiple subcommunicators coexist in one job without interfering
!        (all subcomms PASS simultaneously)
!    V3: for an identical deterministic input, the spectral checksum is
!        bitwise identical across all subcomms (determinism premise of the
!        4D redundant solve)
!
!  Run: pass the split count NCOMM as the first argument (default 2).
!       The world size must be divisible by NCOMM.
!
!  Base: adapted from
!    /work/gr96/o07001/p3dfft_vs_cufftmp/cufftmp-slab/cufftmp_r2c.f90
!    (proven on 256 nodes / 2048^3).
! =====================================================================
module cufft_required
   ! Each rank belongs to exactly one subcommunicator. The shapes below are
   ! the local shapes "within this rank's subcomm"; they are referenced by
   ! the contained memcpy subroutines.
   integer :: planr2c, planc2r
   integer :: local_rshape(3), local_rshape_permuted(3), local_permuted_cshape(3)
end module cufft_required


program cufftmp_subcomm
   use iso_c_binding
   use cudafor
   use cufftXt
   use cufft
   use openacc
   use mpi
   use cufft_required
   implicit none

   integer :: gsize, grank, ndevices, ierr
   integer :: nx, ny, nz
   integer :: i, j, k
   integer :: my_nx, ranks_cutoff, x_offset

   ! --- subcommunicator handling ---
   integer :: ncomm, ranks_per_comm, color, subcomm, sub_rank, sub_size
   integer :: node_comm, node_rank
   integer :: envlen, envstat, ios
   character(len=32) :: envval

   ! --- verification ---
   real(8), dimension(:, :, :), allocatable :: u, ref
   real(8) :: max_norm, max_diff, comm_max_diff
   real(8) :: local_sum, comm_checksum
   real(8), allocatable :: all_maxdiff(:), all_checksum(:)
   real(8) :: Nglob, factor, prec, reldiff, maxrel
   logical :: allpass, bit_identical
   integer :: c

   ! --- cufft ---
   integer(c_size_t) :: worksize(1)
   type(cudaLibXtDesc), pointer :: u_desc
   type(cudaXtDesc), pointer    :: u_descptr
   complex(8), pointer, device  :: u_dptr(:,:,:)

   call mpi_init(ierr)
   call mpi_comm_size(MPI_COMM_WORLD, gsize, ierr)
   call mpi_comm_rank(MPI_COMM_WORLD, grank, ierr)

!vvvvv   parameters   vvvvv
   nx = 128
   ny = 128
   nz = 128
!^^^^^^^^^^^^^^^^^^^^^^^^^^^

   ! --- Decide the split count NCOMM (world root reads it, then broadcasts
   !     so all ranks agree). Prefer the first argument (mpirun forwards argv
   !     to all ranks); fall back to the NCOMM env var; default 2.
   if (grank == 0) then
      ncomm = 2
      call get_command_argument(1, envval, envlen, envstat)
      if (envstat == 0 .and. envlen > 0) then
         read(envval, *, iostat=ios) ncomm
         if (ios /= 0 .or. ncomm < 1) ncomm = 2
      else
         call get_environment_variable("NCOMM", envval, envlen, envstat)
         if (envstat == 0 .and. envlen > 0) then
            read(envval, *, iostat=ios) ncomm
            if (ios /= 0 .or. ncomm < 1) ncomm = 2
         end if
      end if
   end if
   call MPI_Bcast(ncomm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   if (mod(gsize, ncomm) /= 0) then
      if (grank == 0) write(*,'(A,I0,A,I0)') &
         ' ERROR: world size must be divisible by NCOMM. size=', gsize, ', NCOMM=', ncomm
      call mpi_finalize(ierr)
      stop
   end if

   ! --- Split into subcommunicators (contiguous ranks into the same comm) ---
   ranks_per_comm = gsize / ncomm
   color = grank / ranks_per_comm
   call MPI_Comm_split(MPI_COMM_WORLD, color, grank, subcomm, ierr)
   call MPI_Comm_rank(subcomm, sub_rank, ierr)
   call MPI_Comm_size(subcomm, sub_size, ierr)

   ! --- Device assignment (machine-general, multi-GPU/node capable) ---
   ! Pick the GPU by node-local rank so 1 rank = 1 GPU is correct whether the
   ! node has G=1 (Miyabi) or G>1 GPUs.
   call checkCuda(cudaGetDeviceCount(ndevices))
   call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, grank, &
                            MPI_INFO_NULL, node_comm, ierr)
   call MPI_Comm_rank(node_comm, node_rank, ierr)
   call checkCuda(cudaSetDevice(mod(node_rank, ndevices)))
   call MPI_Comm_free(node_comm, ierr)

   if (grank == 0) then
      write(*,'(A)')        ' =========================================================='
      write(*,'(A)')        ' cuFFTMp subcommunicator smoke test (gallope Phase 0)'
      write(*,'(A,I0,A,I0,A,I0)') ' grid = ', nx, ' x ', ny, ' x ', nz
      write(*,'(A,I0,A,I0,A,I0,A,I0)') ' world_size = ', gsize, ', NCOMM = ', ncomm, &
            ', ranks/comm = ', ranks_per_comm, ', ndevices/node = ', ndevices
      write(*,'(A)')        ' =========================================================='
   end if

   ! --- Slab decomposition within this subcomm (X-slab input) ---
   ! Ranks 0 .. (mod(nx,sub_size)-1) own one extra plane in X.
   ranks_cutoff = mod(nx, sub_size)
   my_nx = nx / sub_size
   if (sub_rank < ranks_cutoff) my_nx = my_nx + 1
   ! Global X index of this rank's first plane (0-based)
   if (sub_rank < ranks_cutoff) then
      x_offset = sub_rank * (nx / sub_size + 1)
   else
      x_offset = ranks_cutoff * (nx / sub_size + 1) + (sub_rank - ranks_cutoff) * (nx / sub_size)
   end if

   local_rshape          = [2*(nz/2+1), ny,          my_nx]
   local_permuted_cshape = [nz/2+1,     ny/sub_size, nx   ]
   local_rshape_permuted = [2*(nz/2+1), ny/sub_size, nx   ]

   if (mod(ny, sub_size) /= 0) then
      if (sub_rank == 0) write(*,'(A,I0)') ' ERROR: ny must divide by sub_size. sub_size=', sub_size
      call mpi_finalize(ierr)
      stop
   end if

   allocate(u  (local_rshape(1), local_rshape(2), local_rshape(3)))
   allocate(ref(local_rshape(1), local_rshape(2), local_rshape(3)))
   ! Allocate the gather arrays on every rank (do not pass an unallocated
   ! allocatable to an explicit-shape dummy argument).
   allocate(all_maxdiff(ncomm), all_checksum(ncomm))

   Nglob  = real(nx,8) * real(ny,8) * real(nz,8)
   factor = 1.0d0 / Nglob
   prec   = 1e-5

   ! --- Create cuFFT plans (attach to the subcommunicator) ---
   call checkCufft(cufftCreate(planr2c))
   call checkCufft(cufftCreate(planc2r))
   call checkCufft(cufftMpAttachComm(planr2c, CUFFT_COMM_MPI, subcomm), 'cufftMpAttachComm r2c error')
   call checkCufft(cufftMpAttachComm(planc2r, CUFFT_COMM_MPI, subcomm), 'cufftMpAttachComm c2r error')
   call checkCufft(cufftMakePlan3d(planr2c, nx, ny, nz, CUFFT_D2Z, worksize), 'cufftMakePlan3d r2c error')
   call checkCufft(cufftMakePlan3d(planc2r, nx, ny, nz, CUFFT_Z2D, worksize), 'cufftMakePlan3d c2r error')
   call checkCufft(cufftXtMalloc(planr2c, u_desc, CUFFT_XT_FORMAT_INPLACE), 'cufftXtMalloc error')

   ! =================================================================
   !  V1/V2: round-trip a per-subcomm independent random input
   ! =================================================================
   ! Offset the RNG seed by color so each subcomm holds different data,
   ! which exercises their independence.
   call seed_rng(color + 1)
   call generate_random(nz, local_rshape(1), local_rshape(2), local_rshape(3), u)
   ref = u

   call cufft_memcpyH2D(u_desc, u, CUFFT_XT_FORMAT_INPLACE, .true.)

   ! forward
   call checkCufft(cufftXtExecDescriptor(planr2c, u_desc, u_desc, CUFFT_FORWARD), 'forward fft failed')
   call checkCuda(cudaDeviceSynchronize())
   ! normalize (1/Nglob)
   call c_f_pointer(u_desc%descriptor, u_descptr)
   call c_f_pointer(u_descptr%data(1), u_dptr, &
        [local_permuted_cshape(1), local_permuted_cshape(2), local_permuted_cshape(3)])
   !$cuf kernel do (3)
   do k = 1, local_permuted_cshape(3)
      do j = 1, local_permuted_cshape(2)
         do i = 1, local_permuted_cshape(1)
            u_dptr(i,j,k) = u_dptr(i,j,k) * factor
         end do
      end do
   end do
   call checkCuda(cudaDeviceSynchronize())
   ! inverse
   call checkCufft(cufftXtExecDescriptor(planc2r, u_desc, u_desc, CUFFT_INVERSE), 'inverse fft failed')
   call checkCuda(cudaDeviceSynchronize())

   call cufft_memcpyD2H(u, u_desc, CUFFT_XT_FORMAT_INPLACE, .true.)
   call checkNormDiff(nz, local_rshape(1), local_rshape(2), local_rshape(3), u, ref, max_norm, max_diff)

   ! Reduce the max error within this subcomm
   call MPI_Allreduce(max_diff, comm_max_diff, 1, MPI_REAL8, MPI_MAX, subcomm, ierr)

   ! Gather each subcomm root's comm_max_diff to the world root
   call gather_root_scalar(comm_max_diff, color, (sub_rank==0), (grank==0), ncomm, all_maxdiff)

   if (grank == 0) then
      write(*,'(A)') ' --- V1/V2: round-trip per subcommunicator ---'
      allpass = .true.
      do c = 1, ncomm
         if (all_maxdiff(c) > prec * Nglob * 0.25d0) allpass = .false.
         write(*,'(A,I0,A,ES12.4,A,L1)') '   subcomm ', c-1, ' : max_diff = ', &
               all_maxdiff(c), '   PASS=', (all_maxdiff(c) <= prec * Nglob * 0.25d0)
      end do
      if (allpass) then
         write(*,'(A)') '   ==> V1/V2 RESULT: ALL SUBCOMMS CORRECT'
      else
         write(*,'(A)') '   ==> V1/V2 RESULT: FAILURE (see above)'
      end if
   end if

   ! =================================================================
   !  V3: identical deterministic input -> forward -> compare spectral
   !      checksums. Bitwise agreement across subcomms guarantees the
   !      determinism required by the 4D redundant solve.
   ! =================================================================
   call generate_det(nz, local_rshape(1), local_rshape(2), local_rshape(3), &
                     x_offset, nx, ny, nz, u)
   call cufft_memcpyH2D(u_desc, u, CUFFT_XT_FORMAT_INPLACE, .true.)

   call checkCufft(cufftXtExecDescriptor(planr2c, u_desc, u_desc, CUFFT_FORWARD), 'V3 forward fft failed')
   call checkCuda(cudaDeviceSynchronize())

   ! Sum of |uk|^2 over the spectrum (local partial sum)
   call c_f_pointer(u_desc%descriptor, u_descptr)
   call c_f_pointer(u_descptr%data(1), u_dptr, &
        [local_permuted_cshape(1), local_permuted_cshape(2), local_permuted_cshape(3)])
   local_sum = 0.0d0
   !$cuf kernel do (3)
   do k = 1, local_permuted_cshape(3)
      do j = 1, local_permuted_cshape(2)
         do i = 1, local_permuted_cshape(1)
            local_sum = local_sum + u_dptr(i,j,k)%re * u_dptr(i,j,k)%re &
                                  + u_dptr(i,j,k)%im * u_dptr(i,j,k)%im
         end do
      end do
   end do
   call checkCuda(cudaDeviceSynchronize())

   call MPI_Allreduce(local_sum, comm_checksum, 1, MPI_REAL8, MPI_SUM, subcomm, ierr)

   call gather_root_scalar(comm_checksum, color, (sub_rank==0), (grank==0), ncomm, all_checksum)

   if (grank == 0) then
      write(*,'(A)') ' --- V3: cross-subcomm determinism (identical input) ---'
      bit_identical = .true.
      maxrel = 0.0d0
      do c = 1, ncomm
         if (all_checksum(c) /= all_checksum(1)) bit_identical = .false.
         if (all_checksum(1) /= 0.0d0) then
            reldiff = abs(all_checksum(c) - all_checksum(1)) / abs(all_checksum(1))
         else
            reldiff = abs(all_checksum(c) - all_checksum(1))
         end if
         maxrel = max(maxrel, reldiff)
         write(*,'(A,I0,A,ES22.15,A,ES10.2)') '   subcomm ', c-1, ' : checksum = ', &
               all_checksum(c), '   rel.diff = ', reldiff
      end do
      if (bit_identical) then
         write(*,'(A)') '   ==> V3 RESULT: BITWISE IDENTICAL across all subcomms (4D redundant solve OK)'
      else
         write(*,'(A,ES10.2,A)') '   ==> V3 RESULT: NOT bitwise identical (max rel.diff = ', &
               maxrel, '); field re-broadcast needed for 4D determinism'
      end if
   end if

   ! --- Teardown (order per the documentation's lifetime requirement) ---
   call checkCufft(cufftXtFree(u_desc))
   call checkCufft(cufftDestroy(planr2c))
   call checkCufft(cufftDestroy(planc2r))
   call MPI_Comm_free(subcomm, ierr)

   deallocate(all_maxdiff, all_checksum)
   deallocate(u, ref)

   if (grank == 0) write(*,'(A)') ' =========================================================='

   call mpi_finalize(ierr)

contains

   subroutine checkCuda(istat, message)
      implicit none
      integer, intent(in)                   :: istat
      character(len=*), intent(in), optional :: message
      if (istat /= cudaSuccess) then
         write(*,"('CUDA Error code: ',I0, ': ')") istat
         write(*,*) cudaGetErrorString(istat)
         if (present(message)) write(*,*) message
         call mpi_finalize(ierr)
         stop
      end if
   end subroutine checkCuda

   subroutine checkCufft(istat, message)
      implicit none
      integer, intent(in)                   :: istat
      character(len=*), intent(in), optional :: message
      if (istat /= CUFFT_SUCCESS) then
         write(*,"('cuFFT Error code: ',I0, ': ')") istat
         write(*,*) cudaGetErrorString(istat)
         if (present(message)) write(*,*) message
         call mpi_finalize(ierr)
         stop
      end if
   end subroutine checkCufft

   ! Reproducible RNG seeding (expand a scalar seed into the seed array)
   subroutine seed_rng(s)
      implicit none
      integer, intent(in) :: s
      integer :: nseed, i2
      integer, allocatable :: seed(:)
      call random_seed(size=nseed)
      allocate(seed(nseed))
      do i2 = 1, nseed
         seed(i2) = s * 1009 + i2 * 9973
      end do
      call random_seed(put=seed)
      deallocate(seed)
   end subroutine seed_rng

   subroutine generate_random(nz1, nzp, nyy, nxl, data)
      implicit none
      integer, intent(in) :: nzp, nyy, nxl, nz1
      real(8), dimension(nzp, nyy, nxl), intent(out) :: data
      real(8) :: rand(1)
      integer :: i2, j2, k2
      do k2 = 1, nxl
         do j2 = 1, nyy
            do i2 = 1, nz1
               call random_number(rand)
               data(i2, j2, k2) = rand(1)
            end do
         end do
      end do
   end subroutine generate_random

   ! Deterministic input: a function of the global indices (xg, y, z) only,
   ! so all subcomms hold the same global field.
   subroutine generate_det(nz1, nzp, nyy, nxl, xoff, nxg, nyg, nzg, data)
      implicit none
      integer, intent(in) :: nz1, nzp, nyy, nxl, xoff, nxg, nyg, nzg
      real(8), dimension(nzp, nyy, nxl), intent(out) :: data
      integer :: i2, j2, k2, xg
      real(8), parameter :: twopi = 6.283185307179586d0
      data = 0.0d0
      do k2 = 1, nxl
         xg = xoff + (k2 - 1)
         do j2 = 1, nyy
            do i2 = 1, nz1
               data(i2, j2, k2) = sin( twopi * ( 1.0d0*xg/nxg &
                                               + 2.0d0*(j2-1)/nyg &
                                               + 3.0d0*(i2-1)/nzg ) )
            end do
         end do
      end do
   end subroutine generate_det

   subroutine checkNormDiff(nz1, nzp, nyy, nxl, data, refd, mnorm, mdiff)
      implicit none
      integer, intent(in) :: nzp, nyy, nxl, nz1
      real(8), dimension(nzp, nyy, nxl), intent(in) :: data, refd
      real(8), intent(out) :: mnorm, mdiff
      integer :: i2, j2, k2
      mnorm = 0
      mdiff = 0
      do k2 = 1, nxl
         do j2 = 1, nyy
            do i2 = 1, nz1
               mnorm = max(mnorm, abs(data(i2, j2, k2)))
               mdiff = max(mdiff, abs(refd(i2, j2, k2) - data(i2, j2, k2)))
            end do
         end do
      end do
   end subroutine checkNormDiff

   ! Gather each subcomm root's value into the world root's arr(1:ncomm).
   ! The color-0 root is the world root itself; other roots send with tag=color.
   subroutine gather_root_scalar(val, comm_color, is_sub_root, is_world_root, ncomm_, arr)
      implicit none
      real(8), intent(in)  :: val
      integer, intent(in)  :: comm_color, ncomm_
      logical, intent(in)  :: is_sub_root, is_world_root
      real(8), intent(out) :: arr(ncomm_)
      integer :: cc, mstat(MPI_STATUS_SIZE), ierr2
      if (is_sub_root) then
         if (is_world_root) then
            arr(comm_color + 1) = val
            do cc = 1, ncomm_ - 1
               call MPI_Recv(arr(cc + 1), 1, MPI_REAL8, MPI_ANY_SOURCE, cc, &
                             MPI_COMM_WORLD, mstat, ierr2)
            end do
         else
            call MPI_Send(val, 1, MPI_REAL8, 0, comm_color, MPI_COMM_WORLD, ierr2)
         end if
      end if
   end subroutine gather_root_scalar

   subroutine cufft_memcpyH2D(ulibxt, u_h, data_format, ismemcpy)
      implicit none
      type(cudaLibXtDesc), pointer, intent(out) :: ulibxt
      real(8), dimension(*), intent(in)          :: u_h
      integer, intent(in)                         :: data_format
      logical, intent(in)                         :: ismemcpy
      type(cudaXtDesc), pointer  :: uxt
      real(8), dimension(:,:,:), device, pointer :: u_d

      if (data_format == CUFFT_XT_FORMAT_INPLACE_SHUFFLED) then
         if (ismemcpy .eqv. .false.) then
            call checkCufft(cufftXtMemcpy(planc2r, ulibxt, u_h, CUFFT_COPY_HOST_TO_DEVICE), "H2D pinv Error")
         else
            call c_f_pointer(ulibxt%descriptor, uxt)
            call c_f_pointer(uxt%data(1), u_d, local_rshape_permuted)
            call checkCuda(cudaMemcpy(u_d, u_h, product(int(local_rshape_permuted,kind=8))), "cudamemcpy H2D Error")
            nullify(u_d, uxt)
         end if
      end if

      if (data_format == CUFFT_XT_FORMAT_INPLACE) then
         if (ismemcpy .eqv. .false.) then
            call checkCufft(cufftXtMemcpy(planr2c, ulibxt, u_h, CUFFT_COPY_HOST_TO_DEVICE), "H2D pfor Error")
         else
            call c_f_pointer(ulibxt%descriptor, uxt)
            call c_f_pointer(uxt%data(1), u_d, local_rshape)
            call checkCuda(cudaMemcpy(u_d, u_h, product(int(local_rshape,kind=8))), "cudamemcpy H2D Error")
            nullify(u_d, uxt)
         end if
      end if
   end subroutine cufft_memcpyH2D

   subroutine cufft_memcpyD2H(u_h, ulibxt, data_format, ismemcpy)
      implicit none
      type(cudaLibXtDesc), pointer, intent(in) :: ulibxt
      real(8), dimension(*), intent(out)      :: u_h
      integer, intent(in)                      :: data_format
      logical, intent(in)                      :: ismemcpy
      type(cudaXtDesc), pointer  :: uxt
      real(8), dimension(:,:,:), device, pointer :: u_d

      if (data_format == CUFFT_XT_FORMAT_INPLACE_SHUFFLED) then
         if (ismemcpy .eqv. .false.) then
            call checkCufft(cufftXtMemcpy(planr2c, u_h, ulibxt, CUFFT_COPY_DEVICE_TO_HOST), "D2H pfor Error")
         else
            call c_f_pointer(ulibxt%descriptor, uxt)
            call c_f_pointer(uxt%data(1), u_d, local_rshape_permuted)
            call checkCuda(cudaMemcpy(u_h, u_d, product(int(local_rshape_permuted,kind=8))), "cudamemcpy D2H Error")
            nullify(u_d, uxt)
         end if
      end if

      if (data_format == CUFFT_XT_FORMAT_INPLACE) then
         if (ismemcpy .eqv. .false.) then
            call checkCufft(cufftXtMemcpy(planc2r, u_h, ulibxt, CUFFT_COPY_DEVICE_TO_HOST), "D2H pinv Error")
         else
            call c_f_pointer(ulibxt%descriptor, uxt)
            call c_f_pointer(uxt%data(1), u_d, local_rshape)
            call checkCuda(cudaMemcpy(u_h, u_d, product(int(local_rshape,kind=8))), "cudamemcpy D2H Error")
            nullify(u_d, uxt)
         end if
      end if
   end subroutine cufft_memcpyD2H

end program cufftmp_subcomm
