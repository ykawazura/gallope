!-----------------------------------------------!
!> @author  YK
!! @brief   Field setting for RMHD
!-----------------------------------------------!
module fields
  implicit none

  public :: init_fields, finish_fields
  public :: phi, omg, psi
  public :: phi_old, omg_old, psi_old
  public :: g
  public :: nfields
  public :: iomg, ipsi

  private

  complex(8), allocatable, dimension(:,:,:) :: phi, omg, psi
  complex(8), allocatable, dimension(:,:,:) :: phi_old, omg_old, psi_old
  ! Hermite-moment field g(nkz,nky_local,nkx,0:nm_local+1). The v_parallel axis
  ! (index 4) is decomposed over comm_m; ghosts at 0 and nm_local+1 hold the
  ! m+/-1 neighbours filled by halo_exchange_m. Interior 1:nm_local is real data.
  complex(8), allocatable, dimension(:,:,:,:) :: g
  character(100) :: init_type
  integer :: init_seed

  ! Field index
  integer, parameter :: nfields = 2
  integer, parameter :: iomg = 1, ipsi = 2

  !$acc declare create(phi, omg, psi)
  !$acc declare create(phi_old, omg_old, psi_old)
  !$acc declare create(g)
  !$acc declare copyin(nfields)
  !$acc declare copyin(iomg, ipsi)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of fields
!-----------------------------------------------!
  subroutine init_fields
    use grid, only: nkx, nky_local, nkz
    use grid, only: nm_local
    use params, only: inputfile
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    allocate(src(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(phi    , source=src)
    allocate(omg    , source=src)
    allocate(psi    , source=src)
    allocate(phi_old, source=src)
    allocate(omg_old, source=src)
    allocate(psi_old, source=src)
    !$acc update device(phi, omg, psi)
    !$acc update device(phi_old, omg_old, psi_old)
    deallocate(src)

    allocate(g(nkz, nky_local, nkx, 0:nm_local+1), source=(0.d0, 0.d0))
    !$acc update device(g)

    call read_parameters(inputfile)

    if(init_type == 'zero') then
    endif
    if(init_type == 'single_mode') then
      call init_single_mode
    endif
    if(init_type == 'OT2') then
      call init_OT2
    endif
    ! if(init_type == 'OT3') then
    !   call init_OT3
    ! endif
    ! if(init_type == 'KH') then
    !   call init_KH
    ! endif
    if(init_type == 'random') then
      call init_random
    endif
    if(init_type == 'restart') then
      call restart
    endif
    ! if(init_type == 'convert_restartfiles_to_real') then
    !   call convert_restartfiles_to_real
    ! endif

    ! Hermite moments g. restart reads g inside restart(); the other init types
    ! set g here, per-moment distinct so a mis-mapped m_offset is caught by the
    ! P_m-invariance / restart XOR-fold checks.
    if(init_type /= 'restart') call init_g

  end subroutine init_fields


!-----------------------------------------------!
!> @author  YK
!! @brief   Read inputfile for initial condition
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /initial_condition/ init_type, init_seed

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

    init_type = 'zero'
    ! init_seed < 0 (default): non-reproducible random IC (system clock).
    ! init_seed >= 0        : reproducible IC, required for redundant-solve
    !                         consistency across comm_m/comm_s and for regression.
    init_seed = -1
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=initial_condition,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @brief   Single mode initialization
!-----------------------------------------------!
  subroutine init_single_mode
    use grid, only: nkx, nky_local, nkz
    use grid, only: kprp2, kx, ky
    use mp, only: proc0, proc_id
    implicit none
    integer :: i, j, k

    if(proc0) then
      print '("Single mode initialization")'
    endif

    ! Seed the single global mode (kx=ky=0, kz=2pi/lz) by WAVENUMBER, not by
    ! local index: ky is split across comm_fft, so global ky=0 lives only on
    ! iproc_fft=0 at j=1. A local "j==1" test would fire on every fft rank and
    ! seed one spurious mode per rank.
    !$acc data present(phi, psi, kprp2, kx, ky)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if (kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. k == 2) then
            phi(k, j, i) = 1.d0
            psi(k, j, i) = 1.d0
          endif
          omg(k, j, i) = phi(k, j, i)*(-kprp2(k, j, i))
        enddo
      enddo
    enddo
    !$acc end data
    
    !$acc data present(phi, omg, psi, phi_old, omg_old, psi_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi_old(k, j, i) = phi(k, j, i)
          omg_old(k, j, i) = omg(k, j, i)
          psi_old(k, j, i) = psi(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine init_single_mode
!
!
!-----------------------------------------------!
!> @author  YK
!! @brief   Random initialization
!-----------------------------------------------!
  subroutine init_random
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky, nky_local, nkz
    use grid, only: kx, ky, kz, kprp2
    use grid, only: ntot
    use mp, only: proc0, iproc_fft, broadcast
    use time, only: microsleep
    use cuFFTmp, only: ftran_r2c, btran_c2r
    implicit none
    real(8), allocatable, dimension(:,:,:) :: phi_r, psi_r
    integer :: seedsize
    integer, allocatable :: seed(:)
    integer :: i, j, k

    if(proc0) then
      print *, 'Random initialization'
    endif

    allocate(phi_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(psi_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(phi_r)
    !$acc enter data create(psi_r)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             create random number            v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    call random_seed(size=seedsize)
    allocate(seed(seedsize))

    ! Build the base seed once on proc0 and broadcast it over MPI_COMM_WORLD so
    ! every rank shares it. Offsetting by (iproc_fft+1) makes ranks holding the
    ! same slab (same iproc_fft) across redundant comm_m/comm_s groups draw an
    ! identical random field, while distinct slabs still differ.
    if (proc0) then
      if (init_seed < 0) then
        ! Legacy behaviour: non-reproducible seed from the system clock.
        do i = 1, seedsize
          call system_clock(count=seed(i))
          call microsleep(1000)
          call system_clock(count=seed(i))
        end do
      else
        ! Reproducible seed for redundant-solve / regression testing.
        do i = 1, seedsize
          seed(i) = init_seed + i
        end do
      endif
    endif
    call broadcast(seed)
    call random_seed(put=(iproc_fft+1)*seed(:))

    call random_number(phi_r)
    call random_number(psi_r)

    !$acc update device(phi_r)
    !$acc update device(psi_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             compute r2c transform           v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !$acc data present(psi_r, phi_r, phi, psi)
    call ftran_r2c(phi_r, phi)
    call ftran_r2c(psi_r, psi)
    !$acc end data
    
    !$acc data present(phi, psi)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi(k,j,i)/ntot
          psi(k,j,i) = psi(k,j,i)/ntot
        enddo
      enddo
    enddo
    !$acc end data

    !$acc data present(phi, psi, kx, ky, kz)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if (kx(i)**2 + ky(j)**2 + kz(k)**2 /= 0.d0) then
            phi(k, j, i) = phi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            psi(k, j, i) = psi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
          endif
          if (kx(i)**2 + ky(j)**2 + kz(k)**2 >= (0.8*maxval(kx))**2) then
            phi(k, j, i) = 0.d0
            psi(k, j, i) = 0.d0
          endif
        enddo
      enddo
    enddo
    !$acc end data

    !$acc data present(phi, omg, kprp2)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          omg(k, j, i) = phi(k, j, i)*(-kprp2(k, j, i))
        enddo
      enddo
    enddo
    !$acc end data
    
    !$acc data present(phi, omg, psi, phi_old, omg_old, psi_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi_old(k, j, i) = phi(k, j, i)
          omg_old(k, j, i) = omg(k, j, i)
          psi_old(k, j, i) = psi(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data
    
    !$acc exit data delete(phi_r)
    !$acc exit data delete(psi_r)
    deallocate(phi_r)
    deallocate(psi_r)

  end subroutine init_random


!-----------------------------------------------!
!> @author  YK
!! @brief   Hermite-moment (g) initialization for non-restart runs.
!!          Each moment is seeded distinctly by GLOBAL m so its field is
!!          independent of how the m axis is split over comm_m (P_m), which is
!!          what the 1-A P_m-invariance / restart XOR-fold checks rely on.
!-----------------------------------------------!
  subroutine init_g
    use grid, only: nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot, nm_local, m_offset
    use grid, only: kx, ky
    use mp, only: iproc_fft
    use cuFFTmp, only: ftran_r2c
    implicit none
    real(8), allocatable, dimension(:,:,:) :: g_r
    integer, allocatable :: seed(:)
    integer :: seedsize, i, j, k, mm, m_glob

    if(init_type == 'single_mode') then
      ! Single (k,m) mode: only g_0 is excited (unit amplitude) at the single
      ! parallel mode (kx=0, ky=0, kz(2)); every higher moment starts at zero.
      ! Under nonlinear=F this is the parallel free-streaming IC for the linear
      ! Landau-damping / Hermite-recurrence test. Only the rank holding the
      ! global m=0 moment writes, so the field stays P_m-invariant.
      ! Seed by WAVENUMBER (kx=ky=0 lives only on iproc_fft=0 at the first local
      ! ky), not by local index: a local "j==1" test fires on every fft rank and
      ! would seed one spurious mode per rank, breaking the single-mode IC.
      !$acc data present(g, kx, ky)
      !$acc parallel loop collapse(4)
      do mm = 1, nm_local
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              if (m_offset+mm-1 == 0 .and. kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. k == 2) then
                g(k,j,i,mm) = (1.d0, 0.d0)
              else
                g(k,j,i,mm) = (0.d0, 0.d0)
              endif
            enddo
          enddo
        enddo
      enddo
      !$acc end data
      return
    endif

    if(init_type == 'random') then
      allocate(g_r(nlz_padded, nly, nlx_local), source=0.d0)
      !$acc enter data create(g_r)
      call random_seed(size=seedsize)
      allocate(seed(seedsize))

      do mm = 1, nm_local
        m_glob = m_offset + mm - 1
        ! Deterministic, P_m-independent seed keyed on the global moment index
        ! and the slab (iproc_fft), mirroring the RMHD random-field seeding.
        do i = 1, seedsize
          seed(i) = abs(init_seed) + 1 + i + 97*(m_glob + 1)
        end do
        call random_seed(put=(iproc_fft+1)*seed(:))
        call random_number(g_r)
        !$acc update device(g_r)

        !$acc data present(g_r, g)
        call ftran_r2c(g_r, g(:,:,:,mm))
        !$acc end data
        !$acc data present(g)
        !$acc parallel loop collapse(3)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              g(k,j,i,mm) = g(k,j,i,mm)/ntot
            enddo
          enddo
        enddo
        !$acc end data
      enddo

      !$acc exit data delete(g_r)
      deallocate(g_r)
      deallocate(seed)
      return
    endif

  end subroutine init_g


!-----------------------------------------------!
!> @author  YK
!! @brief   2D Orszag Tang problem initialization
!-----------------------------------------------!
  subroutine init_OT2
    use grid, only: xx, yy
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0
    use params, only: zi, pi
    use cuFFTmp, only: ftran_r2c
    implicit none
    real(8), allocatable, dimension(:,:,:) :: phi_r, psi_r
    real(8), allocatable, dimension(:,:,:) :: src
    integer :: i, j, k

    if(proc0) then
      print *, 'OT2 initialization'
    endif

    allocate(src(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(phi_r, source=src)
    allocate(psi_r, source=src)
    deallocate(src)
    !$acc enter data create(phi_r)
    !$acc enter data create(psi_r)

    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          phi_r(k, j, i) = -cos(xx(i)) - cos(yy(j))
          psi_r(k, j, i) = -0.5d0*cos(2.d0*(xx(i))) - cos(yy(j))
        end do
      end do
    end do
  
    ! compute r2c transform
    call ftran_r2c(phi_r, phi)
    call ftran_r2c(psi_r, psi)
  
    !$acc data present(phi, psi)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi(k,j,i)/ntot
          psi(k,j,i) = psi(k,j,i)/ntot
        enddo
      enddo
    enddo
    !$acc end data

    !$acc data present(phi, omg, psi, phi_old, omg_old, psi_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi_old(k, j, i) = phi(k, j, i)
          omg_old(k, j, i) = omg(k, j, i)
          psi_old(k, j, i) = psi(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data

    !$acc exit data delete(phi_r)
    !$acc exit data delete(psi_r)
    deallocate(phi_r)
    deallocate(psi_r)

  end subroutine init_OT2


!-----------------------------------------------!
!> @author  YK
!! @brief   Restart
!-----------------------------------------------!
  subroutine restart
    use mp, only: proc0, iproc_fft, comm_fft, comm_fm
    use time, only: tt, dt
    use grid, only: nkx, nky, nky_local, nkz
    use grid, only: nm, nm_local, m_offset
    use params, only: restart_dir
    use file, only: open_input_file, close_file
    use mpiio, only: mpiio_read_one
    use shearing_box, only: tsc
    use MPI
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts
    integer, dimension(4) :: sizes4, subsizes4, starts4

    if(proc0) then
      print *, 'Restart'
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                   Read time                 v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    call open_input_file (time_unit, trim(restart_dir)//'time.dat')
    read (unit=time_unit, fmt=*)
    read (unit=time_unit, fmt="(100es30.21)") tt, tsc
    call close_file (time_unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    !$acc update device(tt, tsc)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v               Read Binary file              v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    sizes(1) = nkz
    sizes(2) = nky
    sizes(3) = nkx
    subsizes(1) = nkz
    subsizes(2) = nky_local
    subsizes(3) = nkx
    starts(1) = 0
    ! Offset by iproc_fft (field is decomposed along comm_fft only). Reads are
    ! NOT gated: every comm_fft group reads the same file so each redundant
    ! group loads an identical initial condition.
    starts(2) = nky_local*iproc_fft
    starts(3) = 0

    call mpiio_read_one(phi, sizes, subsizes, starts, trim(restart_dir)//'phi.dat', comm_fft)
    call mpiio_read_one(omg, sizes, subsizes, starts, trim(restart_dir)//'omg.dat', comm_fft)
    call mpiio_read_one(psi, sizes, subsizes, starts, trim(restart_dir)//'psi.dat', comm_fft)

    phi_old = phi
    omg_old = omg
    psi_old = psi

    !$acc update device(phi, omg, psi)
    !$acc update device(phi_old, omg_old, psi_old)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v         Read Hermite moments g (rank-4)      v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    ! g is distributed over the (comm_fft x comm_m) plane = comm_fm. The read is
    ! NOT gated: every comm_s group reads the same file (redundant identical load).
    sizes4    = (/ nkz, nky      , nkx, nm       /)
    subsizes4 = (/ nkz, nky_local, nkx, nm_local /)
    starts4   = (/ 0  , nky_local*iproc_fft, 0, m_offset /)

    call mpiio_read_one(g(:,:,:,1:nm_local), sizes4, subsizes4, starts4, &
                        trim(restart_dir)//'g.dat', comm_fm)
    !$acc update device(g)
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @brief   Finalization of Fields
!-----------------------------------------------!
  subroutine finish_fields
    implicit none

    deallocate(phi)
    deallocate(omg)
    deallocate(psi)
    deallocate(phi_old)
    deallocate(omg_old)
    deallocate(psi_old)
    deallocate(g)

  end subroutine finish_fields

end module fields

