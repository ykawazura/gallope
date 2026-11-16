!-----------------------------------------------!
!> @author  YK
!! @date    25 Feb 2021
!! @brief   Field setting for MHD_COMP_ISOTH
!!          Gallope port: real space f(nlz_padded, nly, nlx_local),
!!          spectral f(nkz, nky_local, nkx) (kz >= 0 half).
!!          Initial conditions are built on the host and
!!          transformed on the device; the conversions
!!          (m_to_u, u_to_m, rho_to_lnrho, lnrho_to_rho)
!!          work on device data.
!-----------------------------------------------!
module fields
  implicit none

  public :: init_fields, finish_fields, m_to_u, u_to_m
  public :: lnrho_to_rho, rho_to_lnrho
  public :: lnrho
  public :: mx, my, mz
  public :: ux, uy, uz
  public :: bx, by, bz
  public :: lnrho_old
  public :: mx_old, my_old, mz_old
  public :: bx_old, by_old, bz_old
  public :: ux_old, uy_old, uz_old
  public :: nfields
  public :: ilnrho
  public :: imx, imy, imz
  public :: ibx, iby, ibz

  private

  complex(8), allocatable, dimension(:,:,:) :: lnrho ! state variable s = ln(rho); rho = exp(s) > 0 by construction
  complex(8), allocatable, dimension(:,:,:) :: mx, my, mz
  complex(8), allocatable, dimension(:,:,:) :: ux, uy, uz
  complex(8), allocatable, dimension(:,:,:) :: bx, by, bz
  complex(8), allocatable, dimension(:,:,:) :: lnrho_old
  complex(8), allocatable, dimension(:,:,:) :: mx_old, my_old, mz_old
  complex(8), allocatable, dimension(:,:,:) :: ux_old, uy_old, uz_old
  complex(8), allocatable, dimension(:,:,:) :: bx_old, by_old, bz_old
  character(100) :: init_type
  integer :: init_seed
  real   (8) :: b0(3)

  ! Field index
  integer, parameter :: nfields = 7
  integer, parameter :: ilnrho = 1
  integer, parameter :: imx  = 2, imy = 3, imz = 4
  integer, parameter :: ibx  = 5, iby = 6, ibz = 7

  !$acc declare create(lnrho, mx, my, mz, ux, uy, uz, bx, by, bz)
  !$acc declare create(lnrho_old, mx_old, my_old, mz_old, ux_old, uy_old, uz_old, bx_old, by_old, bz_old)
  !$acc declare copyin(nfields)
  !$acc declare copyin(ilnrho, imx, imy, imz, ibx, iby, ibz)

contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialization of fields
!-----------------------------------------------!
  subroutine init_fields
    use grid, only: nkx, nky_local, nkz
    use params, only: inputfile
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    allocate(src(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(lnrho  , source=src)
    allocate(mx     , source=src)
    allocate(my     , source=src)
    allocate(mz     , source=src)
    allocate(ux     , source=src)
    allocate(uy     , source=src)
    allocate(uz     , source=src)
    allocate(bx     , source=src)
    allocate(by     , source=src)
    allocate(bz     , source=src)
    allocate(lnrho_old, source=src)
    allocate(mx_old , source=src)
    allocate(my_old , source=src)
    allocate(mz_old , source=src)
    allocate(ux_old , source=src)
    allocate(uy_old , source=src)
    allocate(uz_old , source=src)
    allocate(bx_old , source=src)
    allocate(by_old , source=src)
    allocate(bz_old , source=src)
    deallocate(src)
    !$acc update device(lnrho, mx, my, mz, ux, uy, uz, bx, by, bz)
    !$acc update device(lnrho_old, mx_old, my_old, mz_old, ux_old, uy_old, uz_old, bx_old, by_old, bz_old)

    call read_parameters(inputfile)

    if(init_type == 'zero') then
      call init_zero
    endif
    if(init_type == 'single_mode') then
      call init_single_mode
    endif
    if(init_type == 'random') then
      call init_random
    endif
    if(init_type == 'sound_wave' .or. init_type == 'simple_wave') then
      call init_wave1d
    endif
    if(init_type == 'riemann') then
      call init_riemann
    endif
    if(init_type == 'profile1d') then
      call init_profile1d
    endif
    if(init_type == 'restart') then
      call restart
    endif

    if(       init_type /= 'zero'       .and. init_type /= 'single_mode' .and. init_type /= 'random' &
        .and. init_type /= 'sound_wave' .and. init_type /= 'simple_wave' .and. init_type /= 'riemann' &
        .and. init_type /= 'profile1d'  .and. init_type /= 'restart') then
      print *, 'Unknown init_type: ', trim(init_type)
      stop
    endif

  end subroutine init_fields


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
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
    ! init_seed >= 0        : reproducible IC, required for regression.
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
!! @date    16 Feb 2021
!! @brief   Zero initialization
!-----------------------------------------------!
  subroutine init_zero
    use mp, only: proc0
    implicit none

    if(proc0) then
      print *, 'Zero initialization'
    endif

    call zero_c(lnrho) ! rho = 1

    call zero_c(mx)
    call zero_c(my)
    call zero_c(mz)
    call zero_c(bx)
    call zero_c(by)
    call zero_c(bz)

    call copy_to_old

  end subroutine init_zero


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Single mode initialization
!!          Gallope: the fields are built on the host.
!!          Calliope's (i, k, j) = (kx, kz, ky) mode indices
!!          are kept and placed by put_mode, which also
!!          reproduces the real field that Calliope's c2r
!!          makes of a mode on its ky = 0 plane.
!-----------------------------------------------!
  subroutine init_single_mode
    use mp, only: nproc, proc0
    use grid, only: nlx_local, nly, nlz, nlz_padded, kx, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: inputfile, cs2va2
    use file, only: get_unused_unit
    implicit none
    character(20) :: mode_type
    integer :: mode(3)
    real(8) :: u1(3), b1(3), rho0, rho1
    real(8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r, rho_r
    complex(8), allocatable, dimension(:,:,:) :: rho
    real(8) :: kprp, kpar, omega, alpha
    complex(8) :: val
    integer :: i, j, k

    integer  :: unit, ierr
    namelist /initial_condition_params/ mode_type, mode, b0, b1, u1, rho0, rho1

    if(proc0) then
      print '("Single mode initialization")'
      print *
    endif

    if(nproc /= 1 .and. proc0) then
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      print *, '!              Error!              !'
      print *, '!  nproc must one for single_mode  !'
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      stop
    endif


    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                read inputfile               v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    mode_type = 'arbitrary'
    mode = (/1, 1, 1/)
    b0  = (/0.d0, 0.d0, 0.d0/)
    b1   = 0.d0
    u1   = 0.d0
    rho0 = 1.d0
    rho1 = 0.d0

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')

    read(unit,nml=initial_condition_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if(trim(mode_type) /= 'Alfven'    .and. &
       trim(mode_type) /= 'slow'      .and. &
       trim(mode_type) /= 'fast'      .and. &
       trim(mode_type) /= 'arbitrary' .and. proc0) then
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      print *, '!              Error!              !'
      print *, '!  mode_type must be either        !'
      print *, '!  "Alfven", "slow", "fast".       !'
      print *, '!  or "arbitrary".                 !'
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      stop
    endif

    ! Host arrays; the padding rows k > nlz stay zero
    allocate(rho_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(bx_r , source=rho_r)
    allocate(by_r , source=rho_r)
    allocate(bz_r , source=rho_r)
    allocate(rho(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    !$acc enter data copyin(rho_r, bx_r, by_r, bz_r, rho)

    ! Calliope's mode indices (i, j, k) = (kx, ky, kz) = (2, 1, 2)
    i = 2; j = 1; k = 2

    ! Alfven eigen mode
    if(trim(mode_type) == 'Alfven') then
      rho_r(1:nlz, :, :) = 1.d0
      bz_r (1:nlz, :, :) = 1.d0
      !$acc update device(rho_r, bz_r)
      call ftran_norm(rho_r, rho)
      call ftran_norm(bz_r , bz )
      !$acc update self(rho, bz)

      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      bx = 0.d0
      by = 0.d0

      val = 1.d-1
      call put_mode(uy, i, j, k,  val) ! uy(i, k, j) = 1.d-1
      call put_mode(by, i, j, k, -val) ! by(i, k, j) = -uy(i, k, j)

      ! uy(i, nkz, j) = dconjg(uy(i, k, j)) ! reality condition
      ! by(i, nkz, j) = dconjg(by(i, k, j)) ! reality condition
    endif

    ! slow eigen mode
    if(trim(mode_type) == 'slow') then
      kprp  = kx(2)
      kpar  = kz(2)
      omega = dsqrt( &
                      0.5d0*(kprp**2 + kpar**2)* &
                      (1.d0 + cs2va2 - dsqrt((1.d0 + cs2va2)**2 - 4.d0*cs2va2*kpar**2/(kprp**2 + kpar**2))) &
                   )
      alpha = 1.d0 - (kprp**2 + kpar**2)/omega**2

      rho_r(1:nlz, :, :) = 1.d0
      bz_r (1:nlz, :, :) = 1.d0
      !$acc update device(rho_r, bz_r)
      call ftran_norm(rho_r, rho)
      call ftran_norm(bz_r , bz )
      !$acc update self(rho, bz)

      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      bx = 0.d0
      by = 0.d0

      ! val is Calliope's ux(i, k, j)
      val = 1.d-3
      call put_mode(ux , i, j, k, val)
      call put_mode(uz , i, j, k, alpha*kpar/kprp*val)
      call put_mode(bx , i, j, k, -kpar/omega*val)
      call put_mode(bz , i, j, k,  kprp/omega*val)
      call put_mode(rho, i, j, k, alpha*omega/(kprp*cs2va2)*val)

    endif

    ! fast eigen mode
    if(trim(mode_type) == 'fast') then
      kprp  = kx(2)
      kpar  = kz(2)
      omega = dsqrt( &
                      0.5d0*(kprp**2 + kpar**2)* &
                      (1.d0 + cs2va2 + dsqrt((1.d0 + cs2va2)**2 - 4.d0*cs2va2*kpar**2/(kprp**2 + kpar**2))) &
                   )
      alpha = 1.d0 - (kprp**2 + kpar**2)/omega**2

      rho_r(1:nlz, :, :) = 1.d0
      bz_r (1:nlz, :, :) = 1.d0
      !$acc update device(rho_r, bz_r)
      call ftran_norm(rho_r, rho)
      call ftran_norm(bz_r , bz )
      !$acc update self(rho, bz)

      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      bx = 0.d0
      by = 0.d0

      ! val is Calliope's ux(i, k, j)
      val = 1.d-3
      call put_mode(ux , i, j, k, val)
      call put_mode(uz , i, j, k, alpha*kpar/kprp*val)
      call put_mode(bx , i, j, k, -kpar/omega*val)
      call put_mode(bz , i, j, k,  kprp/omega*val)
      call put_mode(rho, i, j, k, alpha*omega/(kprp*dsqrt(cs2va2))*val)

    endif

    ! Arbitrary mode
    if(trim(mode_type) == 'arbitrary') then
      rho_r(1:nlz, :, :) = rho0
      bx_r (1:nlz, :, :) = b0(1)
      by_r (1:nlz, :, :) = b0(2)
      bz_r (1:nlz, :, :) = b0(3)
      !$acc update device(rho_r, bx_r, by_r, bz_r)
      call ftran_norm(rho_r, rho)
      call ftran_norm(bx_r , bx )
      call ftran_norm(by_r , by )
      call ftran_norm(bz_r , bz )
      !$acc update self(rho, bx, by, bz)

      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      ! Calliope sets the mode (i, j, k) = mode(1:3) if it is on the grid
      call put_mode(rho, mode(1), mode(2), mode(3), dcmplx(rho1 , 0.d0))
      call put_mode(ux , mode(1), mode(2), mode(3), dcmplx(u1(1), 0.d0))
      call put_mode(uy , mode(1), mode(2), mode(3), dcmplx(u1(2), 0.d0))
      call put_mode(uz , mode(1), mode(2), mode(3), dcmplx(u1(3), 0.d0))
      call put_mode(bx , mode(1), mode(2), mode(3), dcmplx(b1(1), 0.d0))
      call put_mode(by , mode(1), mode(2), mode(3), dcmplx(b1(2), 0.d0))
      call put_mode(bz , mode(1), mode(2), mode(3), dcmplx(b1(3), 0.d0))
    endif

    !$acc update device(rho, ux, uy, uz, bx, by, bz)

    call u_to_m(rho, ux, uy, uz, mx, my, mz)
    call rho_to_lnrho(rho, lnrho)

    call copy_to_old

    !$acc exit data delete(rho_r, bx_r, by_r, bz_r, rho)
    deallocate(rho_r)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)
    deallocate(rho)

  end subroutine init_single_mode


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Random initialization
!!          Gallope: the random numbers are drawn in
!!          Calliope's (nproc = 1) global order and
!!          placed by global index (draw_global), so the
!!          field equals Calliope's for the same seed.
!!          Calliope's local ux, uy, uz shadow the module
!!          arrays; here they are ux_, uy_, uz_, and the
!!          module ux, uy, uz stay zero as in Calliope.
!-----------------------------------------------!
  subroutine init_random
    use grid, only: kx, ky, kz, iky, nlx, nly, nlz
    use grid, only: nlx_local, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: kx_max, ky_max, kz_max
    use mp, only: proc0, broadcast, sum_allreduce, comm_fft
    use params, only: zi, cs2va2
    use params, only: inputfile
    use file, only: get_unused_unit
    use time, only: microsleep
    use cuFFTmp, only: btran_c2r
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: phi, psi, ux_, uy_, uz_, rho
    real   (8), allocatable, dimension(:,:,:) :: phi_r, psi_r, rho_r
    real   (8), allocatable, dimension(:,:,:) :: m_r
    real   (8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r
    real   (8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r
    integer :: seedsize
    integer, allocatable :: seed(:)
    real(8) :: b_rms, u_rms, rho_rms, kmin(3), kmax(3)
    real(8) :: u1(3), b1(3), rho0, rho1
    real(8) :: smach_rms, amach_rms, beta_rms
    real(8) :: kxs, kys
    integer :: i, j, k
    character(*), parameter :: fmt_mach = "('rms of (1) sound mach num = ', es10.3, ',  (2) Alfven mach num = ', es10.3, " &
                                       // "',  (3) beta = ', es10.3)"

    integer  :: unit, ierr
    namelist /initial_condition_params/ kmin, kmax, b0, b1, u1, rho0, rho1

    if(proc0) then
      print *, 'Random initialization'
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                read inputfile               v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    kmin = (/0.d0, 0.d0, 0.d0/)
    ! Calliope: (/maxval(kx), maxval(ky), maxval(kz)/), the Nyquist wavenumbers
    kmax = (/kx_max, ky_max, kz_max/)
    b0  = (/0.d0, 0.d0, 0.d0/)
    b1   = 0.d0
    u1   = 0.d0
    rho0 = 1.d0
    rho1 = 0.d0

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')

    read(unit,nml=initial_condition_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    ! Host arrays; the padding rows k > nlz stay zero
    allocate(phi_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(psi_r, source=phi_r)
    allocate(rho_r, source=phi_r)

    allocate(phi(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(psi, source=phi)
    allocate(rho, source=phi)
    !$acc enter data create(phi, psi, rho)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             create random number             v
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    call random_seed(size=seedsize)
    allocate(seed(seedsize))

    ! Calliope seeds every rank with (proc_id+1)*seed and fills its own
    ! slab. Here proc0 builds the seed, every rank puts it unchanged
    ! (Calliope's proc 0 factor 1) and draws the whole grid in Calliope's
    ! nproc = 1 order, keeping its own x slab.
    if (proc0) then
      if (init_seed < 0) then
        ! Non-reproducible seed from the system clock (as Calliope).
        do i = 1, seedsize
          call system_clock(count=seed(i))
          call microsleep(1000)
          call system_clock(count=seed(i))
        end do
      else
        do i = 1, seedsize
          seed(i) = init_seed + i
        end do
      endif
    endif
    call broadcast(seed)
    call random_seed(put=seed(:))

    call draw_global(phi_r)
    call draw_global(psi_r)
    call draw_global(rho_r)
    deallocate(seed)
    !$acc enter data copyin(phi_r, psi_r, rho_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    ! compute r2c transform
    call ftran_norm(phi_r, phi)
    call ftran_norm(psi_r, psi)
    call ftran_norm(rho_r, rho)

    ! phi_r and psi_r are not used any more
    !$acc exit data delete(phi_r, psi_r)
    deallocate(phi_r)
    deallocate(psi_r)

    if(nlz <= 2) then
      ! 2D case
      !$acc parallel loop collapse(3) present(phi, psi, rho, kx, ky, kz) copyin(kmin, kmax)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            if (kx(i)**2 + ky(j)**2 + kz(k)**2 /= 0.d0) then
              phi(k, j, i) = phi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
              psi(k, j, i) = psi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
              rho(k, j, i) = rho(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            endif

            if (kx(i)**2 < kmin(1)**2 .or. ky(j)**2 < kmin(2)**2 .or. &
                kx(i)**2 > kmax(1)**2 .or. ky(j)**2 > kmax(2)**2 .or. &
                kz(k) /= 0.d0) then
              phi(k, j, i) = 0.d0
              psi(k, j, i) = 0.d0
              rho(k, j, i) = 0.d0
            endif
          end do
        end do
      end do
    else
      ! 3D case
      !$acc parallel loop collapse(3) present(phi, psi, rho, kx, ky, kz) copyin(kmin, kmax)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            if (kx(i)**2 + ky(j)**2 + kz(k)**2 /= 0.d0) then
              phi(k, j, i) = phi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
              psi(k, j, i) = psi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
              rho(k, j, i) = rho(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            endif

            if (kx(i)**2 < kmin(1)**2 .or. ky(j)**2 < kmin(2)**2 .or. kz(k)**2 < kmin(3)**2 .or. &
                kx(i)**2 > kmax(1)**2 .or. ky(j)**2 > kmax(2)**2 .or. kz(k)**2 > kmax(3)**2) then
              phi(k, j, i) = 0.d0
              psi(k, j, i) = 0.d0
              rho(k, j, i) = 0.d0
            endif
          end do
        end do
      end do
    endif

    allocate(ux_(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(uy_, source=ux_)
    allocate(uz_, source=ux_)
    !$acc enter data create(ux_, uy_, uz_)

    ! u = z x grad(phi)
    ! Calliope applies ky and kx on its ky >= 0 half, and its c2r keeps
    ! only the Hermitian part of the ky = 0 and ky Nyquist planes. kys and
    ! kxs give the same real field on Gallope's full ky range:
    !   ky Nyquist plane           : kys = 0
    !   kx Nyquist (i = nlx/2 + 1) : kxs = 0 on ky = 0 and the ky Nyquist plane,
    !                                kxs = kx(i) for ky > 0, -kx(i) for ky < 0
    !$acc parallel loop collapse(3) present(phi, psi, ux_, uy_, uz_, bx, by, bz, kx, ky, iky) private(kxs, kys)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          kys = ky(j)
          if (mod(nly, 2) == 0 .and. iky(j) == nly/2) kys = 0.d0
          kxs = kx(i)
          if (mod(nlx, 2) == 0 .and. i == nlx/2 + 1) then
            if (iky(j) == 0 .or. (mod(nly, 2) == 0 .and. iky(j) == nly/2)) then
              kxs = 0.d0
            elseif (iky(j) < 0) then
              kxs = -kx(i)
            endif
          endif
          ux_(k, j, i) = -zi*kys*phi(k, j, i)
          uy_(k, j, i) =  zi*kxs*phi(k, j, i)
          uz_(k, j, i) =  0.d0
          bx (k, j, i) = -zi*kys*psi(k, j, i)
          by (k, j, i) =  zi*kxs*psi(k, j, i)
          bz (k, j, i) =  0.d0
        enddo
      enddo
    enddo

    ! phi and psi are not used any more
    !$acc exit data delete(phi, psi)
    deallocate(phi)
    deallocate(psi)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                   set rho                   v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    ! normalize rho by rms of |rho|
    call btran_c2r(rho, rho_r)
    rho_rms = 0.d0
    !$acc parallel loop collapse(3) present(rho_r) reduction(+:rho_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_rms = rho_rms + rho_r(k, j, i)**2
        end do
      end do
    end do
    call sum_allreduce(rho_rms, comm=comm_fft)
    rho_rms = sqrt(rho_rms/nlx/nly/nlz)

    !$acc parallel loop collapse(3) present(rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_r(k, j, i) = rho1*rho_r(k, j, i) + rho0
        end do
      end do
    end do

    ! recompute rho_rms after adding mean field
    rho_rms = 0.d0
    !$acc parallel loop collapse(3) present(rho_r) reduction(+:rho_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_rms = rho_rms + rho_r(k, j, i)**2
        end do
      end do
    end do
    call sum_allreduce(rho_rms, comm=comm_fft)
    rho_rms = sqrt(rho_rms/nlx/nly/nlz)

    call ftran_norm(rho_r, rho)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!


    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                    set u                    v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    allocate(ux_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(uy_r, source=ux_r)
    allocate(uz_r, source=ux_r)
    allocate(m_r , source=ux_r)
    !$acc enter data create(ux_r, uy_r, uz_r, m_r)

    ! normalize ux, uy, uz by rms of |u|
    call btran_c2r(ux_, ux_r)
    call btran_c2r(uy_, uy_r)
    call btran_c2r(uz_, uz_r)
    u_rms = 0.d0
    !$acc parallel loop collapse(3) present(ux_r, uy_r, uz_r) reduction(+:u_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          u_rms = u_rms + (ux_r(k, j, i)**2 + uy_r(k, j, i)**2 + uz_r(k, j, i)**2)
        end do
      end do
    end do
    call sum_allreduce(u_rms, comm=comm_fft)
    u_rms = sqrt(u_rms/nlx/nly/nlz)

    !$acc parallel loop collapse(3) present(ux_r, uy_r, uz_r) copyin(u1)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          ux_r(k, j, i) = u1(1)*ux_r(k, j, i)/u_rms
          uy_r(k, j, i) = u1(2)*uy_r(k, j, i)/u_rms
          uz_r(k, j, i) = u1(3)*uz_r(k, j, i)/u_rms
        end do
      end do
    end do

    ! recompute u_rms
    u_rms = 0.d0
    !$acc parallel loop collapse(3) present(ux_r, uy_r, uz_r) reduction(+:u_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          u_rms = u_rms + (ux_r(k, j, i)**2 + uy_r(k, j, i)**2 + uz_r(k, j, i)**2)
        end do
      end do
    end do
    call sum_allreduce(u_rms, comm=comm_fft)
    u_rms = sqrt(u_rms/nlx/nly/nlz)

    ! m = rho u. Calliope also transforms ux_r, uy_r, uz_r back into its
    ! local ux, uy, uz, which are never used; that is skipped here.
    !$acc parallel loop collapse(3) present(m_r, rho_r, ux_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          m_r(k, j, i) = rho_r(k, j, i)*ux_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(m_r, mx)

    !$acc parallel loop collapse(3) present(m_r, rho_r, uy_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          m_r(k, j, i) = rho_r(k, j, i)*uy_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(m_r, my)

    !$acc parallel loop collapse(3) present(m_r, rho_r, uz_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          m_r(k, j, i) = rho_r(k, j, i)*uz_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(m_r, mz)

    !$acc exit data delete(ux_r, uy_r, uz_r, m_r, ux_, uy_, uz_)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
    deallocate(m_r)
    deallocate(ux_)
    deallocate(uy_)
    deallocate(uz_)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!


    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                    set b                    v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    allocate(bx_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(by_r, source=bx_r)
    allocate(bz_r, source=bx_r)
    !$acc enter data create(bx_r, by_r, bz_r)

    ! normalize bx, by, bz by rms of |b|
    call btran_c2r(bx, bx_r)
    call btran_c2r(by, by_r)
    call btran_c2r(bz, bz_r)
    b_rms = 0.d0
    !$acc parallel loop collapse(3) present(bx_r, by_r, bz_r) reduction(+:b_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          b_rms = b_rms + (bx_r(k, j, i)**2 + by_r(k, j, i)**2 + bz_r(k, j, i)**2)
        end do
      end do
    end do
    call sum_allreduce(b_rms, comm=comm_fft)
    b_rms = sqrt(b_rms/nlx/nly/nlz)

    !$acc parallel loop collapse(3) present(bx_r, by_r, bz_r) copyin(b1, b0)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          bx_r(k, j, i) = b1(1)*bx_r(k, j, i)/b_rms + b0(1)
          by_r(k, j, i) = b1(2)*by_r(k, j, i)/b_rms + b0(2)
          bz_r(k, j, i) = b1(3)*bz_r(k, j, i)/b_rms + b0(3)
        end do
      end do
    end do

    ! recompute b_rms after adding mean field
    b_rms = 0.d0
    !$acc parallel loop collapse(3) present(bx_r, by_r, bz_r) reduction(+:b_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          b_rms = b_rms + (bx_r(k, j, i)**2 + by_r(k, j, i)**2 + bz_r(k, j, i)**2)
        end do
      end do
    end do
    call sum_allreduce(b_rms, comm=comm_fft)
    b_rms = sqrt(b_rms/nlx/nly/nlz)

    call ftran_norm(bx_r, bx)
    call ftran_norm(by_r, by)
    call ftran_norm(bz_r, bz)

    !$acc exit data delete(bx_r, by_r, bz_r)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    ! check div free of u and b
    call is_div_free('b', bx, by, bz)

    smach_rms  = u_rms/sqrt(cs2va2)
    amach_rms  = sqrt(rho_rms)*u_rms/b_rms
    beta_rms   = 2.d0*cs2va2*rho_rms/b_rms**2
    if(proc0) write(*, fmt_mach) smach_rms, amach_rms, beta_rms

    call rho_to_lnrho(rho, lnrho)

    call copy_to_old

    !$acc exit data delete(rho_r, rho)
    deallocate(rho_r)
    deallocate(rho )
  end subroutine init_random


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   Fill the host array f_r(nlz_padded, nly, nlx_local)
!!          with random numbers in Calliope's nproc = 1 order:
!!          random_number(f(nly, nlz, nlx)) runs over y, then z,
!!          then x, so one (nly, nlz) plane is drawn per global
!!          x index and only the planes of the local slab are
!!          kept. Every rank draws the whole grid. The padding
!!          rows k > nlz are not touched.
!-----------------------------------------------!
  subroutine draw_global(f_r)
    use grid, only: nlx, nly, nlz, nlx_local
    use mp, only: iproc_fft
    implicit none
    real(8), dimension(:,:,:), intent(inout) :: f_r
    real(8), allocatable, dimension(:,:) :: plane
    integer :: ig, il, j, k

    allocate(plane(nly, nlz))
    do ig = 1, nlx
      call random_number(plane)
      il = ig - nlx_local*iproc_fft
      if (il >= 1 .and. il <= nlx_local) then
        do j = 1, nly
          do k = 1, nlz
            f_r(k, j, il) = plane(j, k)
          end do
        end do
      endif
    end do
    deallocate(plane)

  end subroutine draw_global


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   Set one Fourier mode on the host copy of
!!          f(nkz, nky_local, nkx) as Calliope's
!!          f(ic, kc, jc) = val does, where ic, jc, kc
!!          are Calliope's kx, ky (ky >= 0 half) and kz
!!          indices (jc is also Gallope's global ky index).
!!          Calliope's c2r makes the real field of the
!!          stored mode and its conjugate at (-kx, -ky, -kz);
!!          on its ky = 0 and ky Nyquist planes it keeps only
!!          the Hermitian part, so the mode and its mirror get
!!          val/2 and conjg(val)/2 there (the mirror slot of
!!          Calliope holds zero in every caller).
!!          Gallope stores the kz >= 0 half: a slot with
!!          kz < 0 is skipped (its mirror carries it), and a
!!          slot on another rank's ky range is skipped.
!!          Out-of-range indices set nothing, as Calliope.
!-----------------------------------------------!
  subroutine put_mode(f, ic, jc, kc, val)
    use grid, only: nlx, nly, nlz, nkz, nky_local
    use mp, only: iproc_fft
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: f
    integer, intent(in) :: ic, jc, kc
    complex(8), intent(in) :: val
    integer :: im, jm, km
    logical :: herm

    if (ic < 1 .or. ic > nlx .or. jc < 1 .or. jc > nly/2 + 1 .or. kc < 1 .or. kc > nlz) return

    im = mirror(ic, nlx)
    jm = mirror(jc, nly)
    km = mirror(kc, nlz)
    herm = jc == 1 .or. (mod(nly, 2) == 0 .and. jc == nly/2 + 1)

    if (herm) then
      call set_slot(kc, jc, ic, (0.d0, 0.d0), .false.)
      call set_slot(km, jm, im, (0.d0, 0.d0), .false.)
      call set_slot(kc, jc, ic, val/2.d0, .true.)
      call set_slot(km, jm, im, conjg(val)/2.d0, .true.)
    else
      call set_slot(kc, jc, ic, val, .false.)
      call set_slot(km, jm, im, conjg(val), .false.)
    endif

  contains

    !> global index of -k for the FFT index i of an n-point grid
    integer function mirror(i, n)
      integer, intent(in) :: i, n
      if (i == 1) then
        mirror = 1
      else
        mirror = n + 2 - i
      endif
    end function mirror

    !> f(kk, jg, ii) = v (add = .false.) or += v (add = .true.),
    !! if the slot is stored on this rank (jg is the global ky index)
    subroutine set_slot(kk, jg, ii, v, add)
      integer, intent(in) :: kk, jg, ii
      complex(8), intent(in) :: v
      logical, intent(in) :: add
      integer :: j_off

      if (kk > nkz) return
      j_off = nky_local*iproc_fft
      if (jg <= j_off .or. jg > j_off + nky_local) return
      if (add) then
        f(kk, jg - j_off, ii) = f(kk, jg - j_off, ii) + v
      else
        f(kk, jg - j_off, ii) = v
      endif
    end subroutine set_slot

  end subroutine put_mode


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Restart
!!          Gallope: EMHD layout (kz, ky, kx), ky split over comm_fft
!-----------------------------------------------!
  subroutine restart
    use mp, only: proc0, iproc_fft, comm_fft
    use time, only: tt
    use grid, only: nkx, nky, nky_local, nkz
    use params, only: restart_dir
    use file, only: open_input_file, close_file
    use mpiio, only: mpiio_read_one
    use shearing_box, only: tsc
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts

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
    starts(2) = nky_local*iproc_fft
    starts(3) = 0

    call mpiio_read_one(lnrho, sizes, subsizes, starts, trim(restart_dir)//'lnrho.dat', comm_fft)
    call mpiio_read_one(mx , sizes, subsizes, starts, trim(restart_dir)//'mx.dat', comm_fft)
    call mpiio_read_one(my , sizes, subsizes, starts, trim(restart_dir)//'my.dat', comm_fft)
    call mpiio_read_one(mz , sizes, subsizes, starts, trim(restart_dir)//'mz.dat', comm_fft)
    call mpiio_read_one(bx , sizes, subsizes, starts, trim(restart_dir)//'bx.dat', comm_fft)
    call mpiio_read_one(by , sizes, subsizes, starts, trim(restart_dir)//'by.dat', comm_fft)
    call mpiio_read_one(bz , sizes, subsizes, starts, trim(restart_dir)//'bz.dat', comm_fft)

    lnrho_old = lnrho
    mx_old  = mx
    my_old  = my
    mz_old  = mz
    bx_old  = bx
    by_old  = by
    bz_old  = bz

    !$acc update device(lnrho, mx, my, mz, bx, by, bz)
    !$acc update device(lnrho_old, mx_old, my_old, mz_old, bx_old, by_old, bz_old)
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @date    15 Dec 2020
!! @brief   m to u (u = m/exp(lnrho))
!!          Gallope: device data in and out. The FFTs keep
!!          their input, so Calliope's copies are not needed.
!-----------------------------------------------!
  subroutine m_to_u(lnrho, mx, my, mz, ux, uy, uz)
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use cuFFTmp, only: btran_c2r
    implicit none
    complex(8), dimension(:,:,:), intent(in ) :: lnrho, mx, my, mz
    complex(8), dimension(:,:,:), intent(out) :: ux, uy, uz
    real   (8), allocatable, dimension(:,:,:) :: rho_r, w_r
    integer :: i, j, k

    allocate(rho_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(w_r  , source=rho_r)
    !$acc enter data create(rho_r, w_r)

    call btran_c2r(lnrho, rho_r)

    !$acc parallel loop collapse(3) present(rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_r(k, j, i) = exp(rho_r(k, j, i))
        end do
      end do
    end do

    ! ux_r = mx_r/rho_r
    call btran_c2r(mx, w_r)
    !$acc parallel loop collapse(3) present(w_r, rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          w_r(k, j, i) = w_r(k, j, i)/rho_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(w_r, ux)

    ! uy_r = my_r/rho_r
    call btran_c2r(my, w_r)
    !$acc parallel loop collapse(3) present(w_r, rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          w_r(k, j, i) = w_r(k, j, i)/rho_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(w_r, uy)

    ! uz_r = mz_r/rho_r
    call btran_c2r(mz, w_r)
    !$acc parallel loop collapse(3) present(w_r, rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          w_r(k, j, i) = w_r(k, j, i)/rho_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(w_r, uz)

    !$acc exit data delete(rho_r, w_r)
    deallocate(rho_r)
    deallocate(w_r)
  end subroutine m_to_u


!-----------------------------------------------!
!> @author  YK
!! @date    15 Dec 2020
!! @brief   u to m (rho is the density, not ln(rho); used only for initial conditions)
!!          Gallope: device data in and out.
!-----------------------------------------------!
  subroutine u_to_m(rho, ux, uy, uz, mx, my, mz)
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use cuFFTmp, only: btran_c2r
    implicit none
    complex(8), dimension(:,:,:), intent(in ) :: rho, ux, uy, uz
    complex(8), dimension(:,:,:), intent(out) :: mx, my, mz
    real   (8), allocatable, dimension(:,:,:) :: rho_r, w_r
    integer :: i, j, k

    allocate(rho_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(w_r  , source=rho_r)
    !$acc enter data create(rho_r, w_r)

    call btran_c2r(rho, rho_r)

    ! mx_r = ux_r*rho_r
    call btran_c2r(ux, w_r)
    !$acc parallel loop collapse(3) present(w_r, rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          w_r(k, j, i) = w_r(k, j, i)*rho_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(w_r, mx)

    ! my_r = uy_r*rho_r
    call btran_c2r(uy, w_r)
    !$acc parallel loop collapse(3) present(w_r, rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          w_r(k, j, i) = w_r(k, j, i)*rho_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(w_r, my)

    ! mz_r = uz_r*rho_r
    call btran_c2r(uz, w_r)
    !$acc parallel loop collapse(3) present(w_r, rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          w_r(k, j, i) = w_r(k, j, i)*rho_r(k, j, i)
        end do
      end do
    end do
    call ftran_norm(w_r, mz)

    !$acc exit data delete(rho_r, w_r)
    deallocate(rho_r)
    deallocate(w_r)
  end subroutine u_to_m


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   rho to ln(rho) (Fourier to Fourier).
!!          Stops if rho <= 0 at any grid point; no floor is applied.
!!          Gallope: device data in and out.
!-----------------------------------------------!
  subroutine rho_to_lnrho(rho, lnrho)
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use mp, only: proc0, min_allreduce, comm_fft
    use cuFFTmp, only: btran_c2r
    implicit none
    complex(8), dimension(:,:,:), intent(in ) :: rho
    complex(8), dimension(:,:,:), intent(out) :: lnrho
    real   (8), allocatable, dimension(:,:,:) :: rho_r
    real   (8) :: rho_min_r
    integer :: i, j, k

    allocate(rho_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(rho_r)

    call btran_c2r(rho, rho_r)

    rho_min_r = huge(1.d0)
    !$acc parallel loop collapse(3) present(rho_r) reduction(min:rho_min_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_min_r = min(rho_min_r, rho_r(k, j, i))
        end do
      end do
    end do
    call min_allreduce(rho_min_r, comm=comm_fft)
    if(rho_min_r <= 0.d0) then
      if(proc0) then
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        print *, '!              Error!              !'
        print *, '!  initial rho <= 0 somewhere      !'
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        print '("  min(rho) = ", es12.4e3)', rho_min_r
      endif
      stop
    endif

    !$acc parallel loop collapse(3) present(rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_r(k, j, i) = log(rho_r(k, j, i))
        end do
      end do
    end do
    call ftran_norm(rho_r, lnrho)

    !$acc exit data delete(rho_r)
    deallocate(rho_r)
  end subroutine rho_to_lnrho


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   ln(rho) to rho (Fourier to Fourier)
!!          Gallope: device data in and out.
!-----------------------------------------------!
  subroutine lnrho_to_rho(lnrho, rho)
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use cuFFTmp, only: btran_c2r
    implicit none
    complex(8), dimension(:,:,:), intent(in ) :: lnrho
    complex(8), dimension(:,:,:), intent(out) :: rho
    real   (8), allocatable, dimension(:,:,:) :: rho_r
    integer :: i, j, k

    allocate(rho_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(rho_r)

    call btran_c2r(lnrho, rho_r)
    !$acc parallel loop collapse(3) present(rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_r(k, j, i) = exp(rho_r(k, j, i))
        end do
      end do
    end do
    call ftran_norm(rho_r, rho)

    !$acc exit data delete(rho_r)
    deallocate(rho_r)
  end subroutine lnrho_to_rho


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   1D wave along dir ('x', 'y', or 'z') on a uniform background (rho0, b0)
!!          'sound_wave' : standing wave,  ln(rho) = ln(rho0) + amp*cos(k xi), u = 0
!!          'simple_wave': right-going simple wave, ln(rho) = ln(rho0) + amp*sin(k xi),
!!                         u_dir = c_s ln(rho/rho0) (Riemann invariant u - c_s ln(rho) is uniform).
!!                         Exact solution of isothermal Euler (B parallel to dir) until breaking.
!!          Gallope: built on the host (xx is the local x slab).
!-----------------------------------------------!
  subroutine init_wave1d
    use mp, only: proc0
    use grid, only: nlx_local, nly, nlz, nlz_padded, lx, ly, lz, xx, yy, zz
    use params, only: inputfile, cs2va2, pi
    use file, only: get_unused_unit
    implicit none
    character(1) :: dir
    integer  :: mode
    real(8) :: amp, rho0
    real(8) :: kk, xi, cs
    real(8), allocatable, dimension(:,:,:) :: s_r, u_r, m_r, b_r
    integer :: i, j, k
    integer :: unit, ierr
    namelist /wave1d_params/ dir, mode, amp, rho0, b0

    dir  = 'x'
    mode = 1
    amp  = 1.d-3
    rho0 = 1.d0
    b0   = (/0.d0, 0.d0, 0.d0/)

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')
    read(unit,nml=wave1d_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading wave1d_params failed"
    close(unit)

    if(proc0) then
      print '("1D wave initialization: ", A, ", dir = ", A, ", mode = ", I4, ", amp = ", es12.4e3)', &
            trim(init_type), dir, mode, amp
    endif
    if(rho0 <= 0.d0 .or. (dir /= 'x' .and. dir /= 'y' .and. dir /= 'z')) then
      if(proc0) print *, 'init_wave1d: rho0 must be positive and dir must be x, y, or z'
      stop
    endif

    ! Host arrays; the padding rows k > nlz stay zero
    allocate(s_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(u_r, source=s_r)
    allocate(m_r, source=s_r)
    allocate(b_r, source=s_r)

    cs = dsqrt(cs2va2)
    if(dir == 'x') kk = 2.d0*pi*mode/lx
    if(dir == 'y') kk = 2.d0*pi*mode/ly
    if(dir == 'z') kk = 2.d0*pi*mode/lz

    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          if(dir == 'x') xi = xx(i)
          if(dir == 'y') xi = yy(j)
          if(dir == 'z') xi = zz(k)
          if(init_type == 'sound_wave') then
            s_r(k, j, i) = log(rho0) + amp*cos(kk*xi)
            u_r(k, j, i) = 0.d0
          else
            s_r(k, j, i) = log(rho0) + amp*sin(kk*xi)
            u_r(k, j, i) = cs*amp*sin(kk*xi)
          endif
        end do
      end do
    end do
    m_r(1:nlz, :, :) = exp(s_r(1:nlz, :, :))*u_r(1:nlz, :, :)

    !$acc enter data copyin(s_r, m_r, b_r)
    call ftran_norm(s_r, lnrho)
    call zero_c(mx); call zero_c(my); call zero_c(mz)
    if(dir == 'x') then
      call ftran_norm(m_r, mx)
    endif
    if(dir == 'y') then
      call ftran_norm(m_r, my)
    endif
    if(dir == 'z') then
      call ftran_norm(m_r, mz)
    endif
    b_r(1:nlz, :, :) = b0(1)
    !$acc update device(b_r)
    call ftran_norm(b_r, bx)
    b_r(1:nlz, :, :) = b0(2)
    !$acc update device(b_r)
    call ftran_norm(b_r, by)
    b_r(1:nlz, :, :) = b0(3)
    !$acc update device(b_r)
    call ftran_norm(b_r, bz)

    call copy_to_old

    !$acc exit data delete(s_r, m_r, b_r)
    deallocate(s_r)
    deallocate(u_r)
    deallocate(m_r)
    deallocate(b_r)
  end subroutine init_wave1d


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Periodic Riemann problem along dir ('x', 'y', or 'z').
!!          The R state occupies L/4 < xi < 3L/4 and the L state the rest:
!!            H(xi) = [tanh((xi - L/4)/w) - tanh((xi - 3L/4)/w)]/2,  w = width*dl
!!            (width <= 0: step, H = 1/2 at the interfaces),
!!            (rho, u, B) = (1 - H)*(L state) + H*(R state).
!!          The interface at L/4 is the Riemann problem (L | R); the one at 3L/4 is its mirror (R | L).
!!          B_dir must be the same on both sides (div B = 0).
!!          Gallope: built on the host (xx is the local x slab).
!-----------------------------------------------!
  subroutine init_riemann
    use mp, only: proc0
    use grid, only: nlx_local, nly, nlz, nlz_padded, lx, ly, lz, xx, yy, zz, dlx, dly, dlz
    use params, only: inputfile
    use file, only: get_unused_unit
    implicit none
    character(1) :: dir
    real(8) :: rho_L, rho_R, u_L(3), u_R(3), b_L(3), b_R(3), width
    real(8) :: ll, ww, xi, hh
    real(8), allocatable, dimension(:,:,:) :: s_r, den_r, hh_r, f_r ! density is den_r: rho_r is the same name as rho_R
    integer :: i, j, k, idir
    integer :: unit, ierr
    namelist /riemann_params/ dir, rho_L, rho_R, u_L, u_R, b_L, b_R, width

    dir   = 'x'
    rho_L = 1.d0
    rho_R = 1.d0
    u_L   = 0.d0
    u_R   = 0.d0
    b_L   = 0.d0
    b_R   = 0.d0
    width = 1.d0

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')
    read(unit,nml=riemann_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading riemann_params failed"
    close(unit)

    if(dir == 'x') then
      idir = 1; ll = lx; ww = width*dlx
    elseif(dir == 'y') then
      idir = 2; ll = ly; ww = width*dly
    elseif(dir == 'z') then
      idir = 3; ll = lz; ww = width*dlz
    else
      if(proc0) print *, 'init_riemann: dir must be x, y, or z'
      stop
    endif
    if(rho_L <= 0.d0 .or. rho_R <= 0.d0) then
      if(proc0) print *, 'init_riemann: rho_L and rho_R must be positive'
      stop
    endif
    if(b_L(idir) /= b_R(idir)) then
      if(proc0) print *, 'init_riemann: b_L(dir) must be equal to b_R(dir) (div B = 0)'
      stop
    endif

    if(proc0) then
      print '("Riemann initialization: dir = ", A, ", width = ", f6.2, " cells")', dir, width
      print '("  L: rho = ", es12.4e3, ", u = ", 3es12.4e3, ", b = ", 3es12.4e3)', rho_L, u_L, b_L
      print '("  R: rho = ", es12.4e3, ", u = ", 3es12.4e3, ", b = ", 3es12.4e3)', rho_R, u_R, b_R
    endif

    ! Host arrays; only the rows k <= nlz are set (no log of the padding)
    allocate(s_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(den_r, source=s_r)
    allocate(hh_r , source=s_r)
    allocate(f_r  , source=s_r)

    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          if(dir == 'x') xi = xx(i)
          if(dir == 'y') xi = yy(j)
          if(dir == 'z') xi = zz(k)
          if(ww > 0.d0) then
            hh = 0.5d0*(tanh((xi - 0.25d0*ll)/ww) - tanh((xi - 0.75d0*ll)/ww))
          else
            hh = 0.d0
            if(xi > 0.25d0*ll .and. xi < 0.75d0*ll) hh = 1.d0
            if(xi == 0.25d0*ll .or. xi == 0.75d0*ll) hh = 0.5d0
          endif
          hh_r(k, j, i) = hh
        end do
      end do
    end do

    !$acc enter data copyin(s_r, f_r)

    den_r(1:nlz, :, :) = (1.d0 - hh_r(1:nlz, :, :))*rho_L + hh_r(1:nlz, :, :)*rho_R
    s_r  (1:nlz, :, :) = log(den_r(1:nlz, :, :))
    !$acc update device(s_r)
    call ftran_norm(s_r, lnrho)

    f_r(1:nlz, :, :) = den_r(1:nlz, :, :)*((1.d0 - hh_r(1:nlz, :, :))*u_L(1) + hh_r(1:nlz, :, :)*u_R(1))
    !$acc update device(f_r)
    call ftran_norm(f_r, mx)
    f_r(1:nlz, :, :) = den_r(1:nlz, :, :)*((1.d0 - hh_r(1:nlz, :, :))*u_L(2) + hh_r(1:nlz, :, :)*u_R(2))
    !$acc update device(f_r)
    call ftran_norm(f_r, my)
    f_r(1:nlz, :, :) = den_r(1:nlz, :, :)*((1.d0 - hh_r(1:nlz, :, :))*u_L(3) + hh_r(1:nlz, :, :)*u_R(3))
    !$acc update device(f_r)
    call ftran_norm(f_r, mz)
    f_r(1:nlz, :, :) =                     (1.d0 - hh_r(1:nlz, :, :))*b_L(1) + hh_r(1:nlz, :, :)*b_R(1)
    !$acc update device(f_r)
    call ftran_norm(f_r, bx)
    f_r(1:nlz, :, :) =                     (1.d0 - hh_r(1:nlz, :, :))*b_L(2) + hh_r(1:nlz, :, :)*b_R(2)
    !$acc update device(f_r)
    call ftran_norm(f_r, by)
    f_r(1:nlz, :, :) =                     (1.d0 - hh_r(1:nlz, :, :))*b_L(3) + hh_r(1:nlz, :, :)*b_R(3)
    !$acc update device(f_r)
    call ftran_norm(f_r, bz)

    call copy_to_old

    !$acc exit data delete(s_r, f_r)
    deallocate(s_r)
    deallocate(den_r)
    deallocate(hh_r)
    deallocate(f_r)
  end subroutine init_riemann


!-----------------------------------------------!
!> @author  YK
!! @date    2 Nov 2026
!! @brief   1D profile along dir ('x', 'y', or 'z') read from the ASCII table 'table', one row per grid point:
!!            xi, rho, u_x, u_y, u_z, b_x, b_y, b_z
!!          Blank lines and lines starting with '#' are skipped. Every rank reads the whole table and takes its
!!          own range. The rows must be the grid points xi = (l - 1)*L/n in order, rho > 0, and b_dir must be
!!          uniform (div B = 0).
!!          Gallope: built on the host; along x the row of the local
!!          point i is the global index i + nlx_local*iproc_fft.
!-----------------------------------------------!
  subroutine init_profile1d
    use mp, only: proc0, iproc_fft
    use grid, only: nlx, nly, nlz, nlx_local, nlz_padded, lx, ly, lz, xx_global, yy, zz
    use params, only: inputfile
    use file, only: get_unused_unit
    implicit none
    character(1)    :: dir
    character(256)  :: table
    character(1024) :: line
    real(8) :: ll
    real(8), allocatable :: tab(:,:), xi(:)
    real(8), allocatable, dimension(:,:,:) :: s_r, den_r, f_r ! density is den_r, as in init_riemann
    integer :: i, j, k, l, n, nrow, idir
    integer :: unit, ierr
    namelist /profile1d_params/ dir, table

    dir   = 'x'
    table = 'profile1d.dat'

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')
    read(unit,nml=profile1d_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading profile1d_params failed"
    close(unit)

    if(dir == 'x') then
      idir = 1; ll = lx; n = nlx
    elseif(dir == 'y') then
      idir = 2; ll = ly; n = nly
    elseif(dir == 'z') then
      idir = 3; ll = lz; n = nlz
    else
      if(proc0) print *, 'init_profile1d: dir must be x, y, or z'
      stop
    endif
    allocate(xi(n))
    if(idir == 1) xi = xx_global
    if(idir == 2) xi = yy
    if(idir == 3) xi = zz

    allocate(tab(8, n), source=0.d0)
    call get_unused_unit (unit)
    open(unit=unit,file=trim(table),status='old',action='read',iostat=ierr)
    if(ierr /= 0) then
      if(proc0) print *, 'init_profile1d: cannot open ', trim(table)
      stop
    endif
    nrow = 0
    do
      read(unit,'(A)',iostat=ierr) line
      if(ierr /= 0) exit
      line = adjustl(line)
      if(len_trim(line) == 0 .or. line(1:1) == '#') cycle
      nrow = nrow + 1
      if(nrow > n) exit
      read(line,*,iostat=ierr) tab(:, nrow)
      if(ierr /= 0) then
        if(proc0) print '("init_profile1d: cannot read 8 numbers from data row ", I8, " of ", A)', nrow, trim(table)
        stop
      endif
    end do
    close(unit)

    if(nrow /= n) then
      if(proc0) print '("init_profile1d: ", A, " must have ", I8, " data rows (one per grid point along ", A, ")")', &
                      trim(table), n, dir
      stop
    endif
    if(maxval(abs(tab(1, :) - xi)) > 1.d-10*ll) then
      if(proc0) print '("init_profile1d: the first column of ", A, " is not the grid (max error ", es10.2e3, ")")', &
                      trim(table), maxval(abs(tab(1, :) - xi))
      stop
    endif
    if(minval(tab(2, :)) <= 0.d0) then
      if(proc0) print *, 'init_profile1d: rho must be positive'
      stop
    endif
    if(maxval(tab(5 + idir, :)) /= minval(tab(5 + idir, :))) then
      if(proc0) print *, 'init_profile1d: b_dir must be uniform (div B = 0)'
      stop
    endif

    if(proc0) then
      print '("1D profile initialization: dir = ", A, ", table = ", A)', dir, trim(table)
      print '("  rho in [", es12.4e3, ", ", es12.4e3, "], b_dir = ", es12.4e3)', &
            minval(tab(2, :)), maxval(tab(2, :)), tab(5 + idir, 1)
    endif

    ! Host arrays; only the rows k <= nlz are set (no log of the padding)
    allocate(s_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(den_r, source=s_r)
    allocate(f_r  , source=s_r)
    !$acc enter data copyin(s_r, f_r)

    call fill(2); den_r = f_r
    s_r(1:nlz, :, :) = log(den_r(1:nlz, :, :))
    !$acc update device(s_r)
    call ftran_norm(s_r, lnrho)

    call fill(3); f_r(1:nlz, :, :) = den_r(1:nlz, :, :)*f_r(1:nlz, :, :)
    !$acc update device(f_r)
    call ftran_norm(f_r, mx)
    call fill(4); f_r(1:nlz, :, :) = den_r(1:nlz, :, :)*f_r(1:nlz, :, :)
    !$acc update device(f_r)
    call ftran_norm(f_r, my)
    call fill(5); f_r(1:nlz, :, :) = den_r(1:nlz, :, :)*f_r(1:nlz, :, :)
    !$acc update device(f_r)
    call ftran_norm(f_r, mz)
    call fill(6)
    !$acc update device(f_r)
    call ftran_norm(f_r, bx)
    call fill(7)
    !$acc update device(f_r)
    call ftran_norm(f_r, by)
    call fill(8)
    !$acc update device(f_r)
    call ftran_norm(f_r, bz)

    call copy_to_old

    !$acc exit data delete(s_r, f_r)
    deallocate(s_r)
    deallocate(den_r)
    deallocate(f_r)
    deallocate(tab)
    deallocate(xi)

  contains

    !> f_r = column col of the table at the local grid points (host)
    subroutine fill(col)
      integer, intent(in) :: col
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            if(idir == 1) l = i + nlx_local*iproc_fft
            if(idir == 2) l = j
            if(idir == 3) l = k
            f_r(k, j, i) = tab(col, l)
          end do
        end do
      end do
    end subroutine fill

  end subroutine init_profile1d


!-----------------------------------------------!
!> @author  YK
!! @date    3 Oct 2020
!! @brief   Check if divergence free is satisfied
!!          Gallope: summed on the device over the stored
!!          kz >= 0 half (Calliope sums over its ky >= 0 half),
!!          so the printed ratio is not identical to Calliope's.
!-----------------------------------------------!
  subroutine is_div_free(name, fx, fy, fz)
    use grid, only: kx, ky, kz, k2
    use grid, only: nkx, nky_local, nkz
    use mp, only: proc0
    use mp, only: sum_reduce, comm_fft
    implicit none
    character(*) :: name
    integer :: i, j, k
    complex(8), dimension(:,:,:), intent(in) :: fx, fy, fz
    real(8) :: abs_div_f_sum, abs_f_sum, abs_k_sum

    abs_div_f_sum = 0.d0
    abs_f_sum     = 0.d0
    abs_k_sum     = 0.d0
    ! abs_f_sum is the real part of the sum of sqrt(fx**2 + fy**2 + fz**2), as Calliope
    !$acc parallel loop collapse(3) present(fx, fy, fz, kx, ky, kz, k2) &
    !$acc reduction(+:abs_div_f_sum, abs_f_sum, abs_k_sum)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          abs_div_f_sum = abs_div_f_sum + abs(kx(i)*fx(k, j, i) + ky(j)*fy(k, j, i) + kz(k)*fz(k, j, i))
          abs_f_sum     = abs_f_sum     + dble(sqrt(fx(k, j, i)**2 + fy(k, j, i)**2 + fz(k, j, i)**2))
          abs_k_sum     = abs_k_sum     + sqrt(k2(k, j, i))
        end do
      end do
    end do
    call sum_reduce(abs_div_f_sum, 0, comm=comm_fft)
    call sum_reduce(abs_f_sum, 0, comm=comm_fft)
    call sum_reduce(abs_k_sum, 0, comm=comm_fft)

    if(proc0) write(*, "('<|k.',A2,'|>/(<|k|><|',A2,'|>) = ', es10.3)") &
                    trim(name), trim(name), abs_div_f_sum/(abs_f_sum*abs_k_sum)

  end subroutine is_div_free


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   fk = r2c(fr)/nlx/nly/nlz on the device
!!          (Calliope: call p3dfft_ftran_r2c(fr, fk, 'fft'); fk = fk/nlx/nly/nlz)
!-----------------------------------------------!
  subroutine ftran_norm(fr, fk)
    use grid, only: nlx, nly, nlz
    use grid, only: nkx, nky_local, nkz
    use cuFFTmp, only: ftran_r2c
    implicit none
    real   (8), dimension(:,:,:), intent(in ) :: fr
    complex(8), dimension(:,:,:), intent(out) :: fk
    integer :: i, j, k

    call ftran_r2c(fr, fk)

    !$acc parallel loop collapse(3) present(fk)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fk(k, j, i) = fk(k, j, i)/nlx/nly/nlz
        end do
      end do
    end do

  end subroutine ftran_norm


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   f = 0 on the device
!-----------------------------------------------!
  subroutine zero_c(f)
    use grid, only: nkx, nky_local, nkz
    implicit none
    complex(8), dimension(:,:,:), intent(out) :: f
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(f)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          f(k, j, i) = 0.d0
        end do
      end do
    end do

  end subroutine zero_c


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   *_old = * on the device
!!          (lnrho, m and b, as Calliope; u_old is not set)
!-----------------------------------------------!
  subroutine copy_to_old
    use grid, only: nkx, nky_local, nkz
    implicit none
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(lnrho, mx, my, mz, bx, by, bz) &
    !$acc present(lnrho_old, mx_old, my_old, mz_old, bx_old, by_old, bz_old)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          lnrho_old(k, j, i) = lnrho(k, j, i)
          mx_old   (k, j, i) = mx   (k, j, i)
          my_old   (k, j, i) = my   (k, j, i)
          mz_old   (k, j, i) = mz   (k, j, i)
          bx_old   (k, j, i) = bx   (k, j, i)
          by_old   (k, j, i) = by   (k, j, i)
          bz_old   (k, j, i) = bz   (k, j, i)
        enddo
      enddo
    enddo

  end subroutine copy_to_old


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Finalization of Fields
!-----------------------------------------------!
  subroutine finish_fields
    use shock_viscosity, only: report_timers, finish_shock_viscosity
    implicit none

    call report_timers
    call finish_shock_viscosity

    deallocate(lnrho)
    deallocate(mx)
    deallocate(my)
    deallocate(mz)
    deallocate(ux)
    deallocate(uy)
    deallocate(uz)
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
    deallocate(lnrho_old)
    deallocate(mx_old)
    deallocate(my_old)
    deallocate(mz_old)
    deallocate(ux_old)
    deallocate(uy_old)
    deallocate(uz_old)
    deallocate(bx_old)
    deallocate(by_old)
    deallocate(bz_old)

  end subroutine finish_fields

end module fields
