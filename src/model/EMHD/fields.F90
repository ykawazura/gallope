!-----------------------------------------------!
!> @author  YK
!! @brief   Field setting for EMHD
!!          Evolved field: b = (bx, by, bz), divergence free.
!!          The guide field b0 is the k = 0 mode of b.
!-----------------------------------------------!
module fields
  implicit none

  public :: init_fields, finish_fields
  public :: bx, by, bz
  public :: bx_old, by_old, bz_old
  public :: nfields
  public :: ibx, iby, ibz
  public :: b0
  public :: div_free, div_b_rel

  private

  complex(8), allocatable, dimension(:,:,:) :: bx, by, bz
  complex(8), allocatable, dimension(:,:,:) :: bx_old, by_old, bz_old
  character(100) :: init_type
  integer :: init_seed
  ! Guide field, set as the k = 0 mode by every initial condition but restart
  real(8) :: b0(3)
  ! random: rms of the fluctuation |b| and the band kmin(d) <= |k_d| <= kmax(d)
  !         (wavenumbers, i.e. 2 pi x mode number/box length)
  real(8) :: b1
  real(8) :: kmin(3), kmax(3)
  ! single_mode: global mode numbers (kx, ky, kz in units of 2 pi/l),
  !              amplitude of the Fourier coefficient and helicity sign
  integer :: mode_kx, mode_ky, mode_kz
  real(8) :: mode_amp
  integer :: mode_hel

  ! Field index
  integer, parameter :: nfields = 3
  integer, parameter :: ibx = 1, iby = 2, ibz = 3

  !$acc declare create(bx, by, bz)
  !$acc declare create(bx_old, by_old, bz_old)
  !$acc declare copyin(nfields)
  !$acc declare copyin(ibx, iby, ibz)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of fields
!-----------------------------------------------!
  subroutine init_fields
    use grid, only: nkx, nky_local, nkz
    use params, only: inputfile
    use mp, only: proc0
    use MPI
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src
    real(8) :: divb
    integer :: ierr

    allocate(src(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(bx    , source=src)
    allocate(by    , source=src)
    allocate(bz    , source=src)
    allocate(bx_old, source=src)
    allocate(by_old, source=src)
    allocate(bz_old, source=src)
    !$acc update device(bx, by, bz)
    !$acc update device(bx_old, by_old, bz_old)
    deallocate(src)

    call read_parameters(inputfile)

    select case (trim(init_type))
    case ('zero')
      call set_guide_field
      call copy_to_old
    case ('single_mode')
      call init_single_mode
    case ('OT3')
      call init_OT3
    case ('random')
      call init_random
    case ('restart')
      call restart
    case default
      if (proc0) print *, 'Unknown init_type for EMHD: ', trim(init_type)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end select

    ! div_b_rel is collective: evaluate it on every rank
    divb = div_b_rel(bx, by, bz)
    if (proc0) print '("  max_k |k.b|/(|k||b|) of the initial field = ", es10.3)', divb

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
    namelist /initial_condition_params/ b0, b1, kmin, kmax, &
                                        mode_kx, mode_ky, mode_kz, mode_amp, mode_hel

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    init_type = 'zero'
    ! init_seed < 0 (default): non-reproducible random IC (system clock).
    ! init_seed >= 0        : reproducible IC, required for regression.
    init_seed = -1
    b0   = 0.d0
    b1   = 0.d0
    kmin = 0.d0
    kmax = huge(1.d0)
    mode_kx  = 0
    mode_ky  = 0
    mode_kz  = 1
    mode_amp = 1.d0
    mode_hel = 1
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')
    read(unit,nml=initial_condition,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')
    read(unit,nml=initial_condition_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition_params failed"
    close(unit)

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @brief   Set the k = 0 mode of b to the guide
!!          field b0 (on the device)
!-----------------------------------------------!
  subroutine set_guide_field
    use mp, only: iproc_fft
    implicit none

    ! ky(1) = 0 lives on iproc_fft = 0 only
    if (iproc_fft /= 0) return

    !$acc serial present(bx, by, bz) copyin(b0)
    bx(1, 1, 1) = dcmplx(b0(1), 0.d0)
    by(1, 1, 1) = dcmplx(b0(2), 0.d0)
    bz(1, 1, 1) = dcmplx(b0(3), 0.d0)
    !$acc end serial

  end subroutine set_guide_field


!-----------------------------------------------!
!> @author  YK
!! @brief   Single helical mode
!!
!!          b = b0 + mode_amp*e_s exp(i k.x) + c.c.,
!!          k = 2 pi (mode_kx/lx, mode_ky/ly, mode_kz/lz),
!!          k x e_s = i s |k| e_s, s = mode_hel = +-1,
!!          e_s = (e1 - i s e2)/sqrt(2), e2 = k/|k| x e1.
!!          This is divergence free and an exact solution of
!!          the full (nonlinear) equation: u = -i k x b is
!!          parallel to b, so b* x u has no k /= 0 part
!!          except b0 x u, and
!!            b_k(t) = b_k(0) exp(i s (k.b0)|k| t/(1 + de^2 k^2))
!!          times the resistive decay.
!!          The mode is given by its global wavenumber, so the
!!          result does not depend on the decomposition.
!-----------------------------------------------!
  subroutine init_single_mode
    use grid, only: nkx, nky, nky_local, nkz
    use grid, only: lx, ly, lz
    use params, only: pi, zi
    use mp, only: proc0, iproc_fft
    use MPI
    implicit none
    integer :: im, jm, km, j_off
    integer :: sgn, ierr
    real(8) :: kv(3), kh(3), a(3), e1(3), e2(3)
    complex(8) :: es(3)

    if(proc0) then
      print '("Single helical mode: (kx, ky, kz) mode number = (", i0, ", ", i0, ", ", i0, "), s = ", i0)', &
            mode_kx, mode_ky, mode_kz, mode_hel
    endif

    if ((mode_kx == 0 .and. mode_ky == 0 .and. mode_kz == 0) .or. abs(mode_hel) /= 1) then
      if (proc0) print *, 'single_mode: k = 0 or mode_hel /= +-1'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif

    ! The r2c axis is kz, so only kz >= 0 is stored. A mode given with kz < 0
    ! is represented by its -k partner; e_s is built for the stored k, and
    ! the helicity sign is the same for k and -k.
    sgn = 1
    if (mode_kz < 0) sgn = -1
    im = modulo(sgn*mode_kx, nkx) + 1
    jm = modulo(sgn*mode_ky, nky) + 1
    km = sgn*mode_kz + 1
    if (km > nkz .or. abs(mode_kx) > nkx/2 .or. abs(mode_ky) > nky/2) then
      if (proc0) print *, 'single_mode: mode number exceeds the grid'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif

    kv = 2.d0*pi*sgn*[dble(mode_kx)/lx, dble(mode_ky)/ly, dble(mode_kz)/lz]
    kh = kv/norm2(kv)
    if (abs(kh(3)) < 0.9d0) then
      a = [0.d0, 0.d0, 1.d0]
    else
      a = [1.d0, 0.d0, 0.d0]
    endif
    e1 = a - dot_product(a, kh)*kh
    e1 = e1/norm2(e1)
    e2 = [kh(2)*e1(3) - kh(3)*e1(2), kh(3)*e1(1) - kh(1)*e1(3), kh(1)*e1(2) - kh(2)*e1(1)]
    es = mode_amp*(e1 - zi*mode_hel*e2)/sqrt(2.d0)

    ! Fields are zero on the host and the device at this point, so the mode is
    ! set on the host and copied over.
    j_off = nky_local*iproc_fft
    call set_mode(im, jm, es)
    ! On the kz = 0 plane the -k partner is stored too and carries conjg.
    if (km == 1) call set_mode(modulo(-sgn*mode_kx, nkx) + 1, modulo(-sgn*mode_ky, nky) + 1, conjg(es))
    !$acc update device(bx, by, bz)

    call set_guide_field
    call copy_to_old

  contains

    subroutine set_mode(ii, jj, v)
      implicit none
      integer, intent(in) :: ii, jj
      complex(8), intent(in) :: v(3)

      if (jj > j_off .and. jj <= j_off + nky_local) then
        bx(km, jj - j_off, ii) = v(1)
        by(km, jj - j_off, ii) = v(2)
        bz(km, jj - j_off, ii) = v(3)
      endif
    end subroutine set_mode

  end subroutine init_single_mode


!-----------------------------------------------!
!> @author  YK
!! @brief   Random initialization
!!
!!          Uniform random numbers are transformed, the
!!          modes outside the band kmin(d) <= |k_d| <= kmax(d)
!!          and k = 0 are removed, b is projected onto
!!          div b = 0 and scaled so that the rms of the
!!          fluctuation |b| is b1; then b0 is added.
!!          (Calliope scales each component by b1(d), which
!!          breaks div b = 0, and its kmin test removes every
!!          mode with a zero component; neither is copied.)
!-----------------------------------------------!
  subroutine init_random
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: kx, ky, kz
    use grid, only: ntot
    use mp, only: proc0, iproc_fft, broadcast, sum_allreduce, comm_fft
    use time, only: microsleep
    use cuFFTmp, only: ftran_r2c
    use MPI
    implicit none
    real(8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r
    integer :: seedsize
    integer, allocatable :: seed(:)
    real(8) :: rms, fac
    integer :: i, j, k, ierr

    if(proc0) then
      print *, 'Random initialization'
    endif

    allocate(bx_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(by_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(bz_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(bx_r, by_r, bz_r)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             create random number            v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    call random_seed(size=seedsize)
    allocate(seed(seedsize))

    ! Build the base seed once on proc0 and broadcast it, then offset by
    ! (iproc_fft+1): distinct slabs differ, and the field does not depend on
    ! anything but init_seed and the decomposition.
    if (proc0) then
      if (init_seed < 0) then
        ! Non-reproducible seed from the system clock.
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
    call random_seed(put=(iproc_fft+1)*seed(:))

    call random_number(bx_r)
    call random_number(by_r)
    call random_number(bz_r)
    ! in-place r2c padding
    bx_r(nlz+1:, :, :) = 0.d0
    by_r(nlz+1:, :, :) = 0.d0
    bz_r(nlz+1:, :, :) = 0.d0

    !$acc update device(bx_r, by_r, bz_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !$acc data present(bx_r, by_r, bz_r, bx, by, bz)
    call ftran_r2c(bx_r, bx)
    call ftran_r2c(by_r, by)
    call ftran_r2c(bz_r, bz)
    !$acc end data

    !$acc parallel loop collapse(3) present(bx, by, bz, kx, ky, kz) copyin(kmin, kmax)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if ( abs(kx(i)) < kmin(1) .or. abs(kx(i)) > kmax(1) .or. &
               abs(ky(j)) < kmin(2) .or. abs(ky(j)) > kmax(2) .or. &
               abs(kz(k)) < kmin(3) .or. abs(kz(k)) > kmax(3) .or. &
               (i == 1 .and. ky(j) == 0.d0 .and. k == 1) ) then
            bx(k, j, i) = 0.d0
            by(k, j, i) = 0.d0
            bz(k, j, i) = 0.d0
          else
            bx(k, j, i) = bx(k, j, i)/ntot
            by(k, j, i) = by(k, j, i)/ntot
            bz(k, j, i) = bz(k, j, i)/ntot
          endif
        enddo
      enddo
    enddo

    call div_free(bx, by, bz)

    ! rms of |b| by Parseval: <|b|^2> = sum over all k of |b_k|^2
    rms = 0.d0
    !$acc parallel loop collapse(3) present(bx, by, bz) private(fac) reduction(+:rms)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fac = 1.d0
          if (k /= 1) fac = 2.d0
          rms = rms + fac*( dble(bx(k, j, i))**2 + aimag(bx(k, j, i))**2 &
                          + dble(by(k, j, i))**2 + aimag(by(k, j, i))**2 &
                          + dble(bz(k, j, i))**2 + aimag(bz(k, j, i))**2 )
        enddo
      enddo
    enddo
    call sum_allreduce(rms, comm=comm_fft)
    rms = sqrt(rms)
    if (rms == 0.d0) then
      if (proc0) print *, 'random: no mode left in the band kmin..kmax'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif

    fac = b1/rms
    !$acc parallel loop collapse(3) present(bx, by, bz) copyin(fac)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          bx(k, j, i) = fac*bx(k, j, i)
          by(k, j, i) = fac*by(k, j, i)
          bz(k, j, i) = fac*bz(k, j, i)
        enddo
      enddo
    enddo

    call set_guide_field
    call copy_to_old

    !$acc exit data delete(bx_r, by_r, bz_r)
    deallocate(bx_r, by_r, bz_r)

  end subroutine init_random


!-----------------------------------------------!
!> @author  YK
!! @brief   3D Orszag Tang type initialization
!!            bx = -sin(y/y0 + z/z0)
!!            by =  sin(2x/x0 + z/z0)
!!            bz =  0
!!          plus b0 (x0 = lx/(2 pi) etc.). Divergence free.
!!          Calliope's EMHD has no OT3; the same field is
!!          added there by tests/EMHD/calliope_OT3.patch.
!-----------------------------------------------!
  subroutine init_OT3
    use grid, only: lx, ly, lz, xx, yy, zz
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0
    use params, only: pi
    use cuFFTmp, only: ftran_r2c
    implicit none
    real(8), allocatable, dimension(:,:,:) :: bx_r, by_r
    real(8) :: x0, y0, z0
    integer :: i, j, k

    if(proc0) then
      print *, 'OT3 initialization'
    endif

    x0 = lx/(2.d0*pi)
    y0 = ly/(2.d0*pi)
    z0 = lz/(2.d0*pi)

    allocate(bx_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(by_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(bx_r, by_r)

    ! zz is dimensioned nlz; the nlz_padded - nlz padding slots are zeroed.
    !$acc parallel loop collapse(3) present(bx_r, by_r, xx, yy, zz) copyin(x0, y0, z0)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          if (k <= nlz) then
            bx_r(k, j, i) = -sin(yy(j)/y0 + zz(k)/z0)
            by_r(k, j, i) =  sin(2.d0*xx(i)/x0 + zz(k)/z0)
          else
            bx_r(k, j, i) = 0.d0
            by_r(k, j, i) = 0.d0
          endif
        end do
      end do
    end do

    call ftran_r2c(bx_r, bx)
    call ftran_r2c(by_r, by)

    !$acc parallel loop collapse(3) present(bx, by, bz)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          bx(k, j, i) = bx(k, j, i)/ntot
          by(k, j, i) = by(k, j, i)/ntot
          bz(k, j, i) = 0.d0
        enddo
      enddo
    enddo

    call set_guide_field
    call copy_to_old

    !$acc exit data delete(bx_r, by_r)
    deallocate(bx_r, by_r)

  end subroutine init_OT3


!-----------------------------------------------!
!> @author  YK
!! @brief   *_old = * on the device
!-----------------------------------------------!
  subroutine copy_to_old
    use grid, only: nkx, nky_local, nkz
    implicit none
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(bx, by, bz, bx_old, by_old, bz_old)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          bx_old(k, j, i) = bx(k, j, i)
          by_old(k, j, i) = by(k, j, i)
          bz_old(k, j, i) = bz(k, j, i)
        enddo
      enddo
    enddo

  end subroutine copy_to_old


!-----------------------------------------------!
!> @author  YK
!! @brief   Restart (b0 is part of the k = 0 mode
!!          in the files; the input b0 is not used)
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

    call mpiio_read_one(bx, sizes, subsizes, starts, trim(restart_dir)//'bx.dat', comm_fft)
    call mpiio_read_one(by, sizes, subsizes, starts, trim(restart_dir)//'by.dat', comm_fft)
    call mpiio_read_one(bz, sizes, subsizes, starts, trim(restart_dir)//'bz.dat', comm_fft)

    bx_old = bx
    by_old = by
    bz_old = bz

    !$acc update device(bx, by, bz)
    !$acc update device(bx_old, by_old, bz_old)
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @brief   Projection onto divergence-free fields,
!!          w <- w - k (k.w)/k^2 (on the device)
!-----------------------------------------------!
  subroutine div_free(wx, wy, wz)
    use grid, only: nkx, nky_local, nkz
    use grid, only: kx, ky, kz, k2inv
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: wx, wy, wz
    complex(8) :: kw
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(wx, wy, wz, kx, ky, kz, k2inv) private(kw)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          kw = (kx(i)*wx(k, j, i) + ky(j)*wy(k, j, i) + kz(k)*wz(k, j, i))*k2inv(k, j, i)
          wx(k, j, i) = wx(k, j, i) - kx(i)*kw
          wy(k, j, i) = wy(k, j, i) - ky(j)*kw
          wz(k, j, i) = wz(k, j, i) - kz(k)*kw
        enddo
      enddo
    enddo

  end subroutine div_free


!-----------------------------------------------!
!> @author  YK
!! @brief   max over modes of |k.w|/(|k||w|),
!!          a round-off level check of div w = 0
!!          (reduced over comm_fft, on every rank)
!-----------------------------------------------!
  function div_b_rel(wx, wy, wz) result(r)
    use grid, only: nkx, nky_local, nkz
    use grid, only: kx, ky, kz, k2
    use mp, only: max_allreduce, comm_fft
    implicit none
    complex(8), dimension(:,:,:), intent(in) :: wx, wy, wz
    real(8) :: r, w2
    integer :: i, j, k

    r = 0.d0
    !$acc parallel loop collapse(3) present(wx, wy, wz, kx, ky, kz, k2) private(w2) reduction(max:r)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          w2 = abs(wx(k, j, i))**2 + abs(wy(k, j, i))**2 + abs(wz(k, j, i))**2
          if (k2(k, j, i)*w2 > 0.d0) then
            r = max(r, abs(kx(i)*wx(k, j, i) + ky(j)*wy(k, j, i) + kz(k)*wz(k, j, i)) &
                       /sqrt(k2(k, j, i)*w2))
          endif
        enddo
      enddo
    enddo
    call max_allreduce(r, comm=comm_fft)

  end function div_b_rel


!-----------------------------------------------!
!> @author  YK
!! @brief   Finalization of Fields
!-----------------------------------------------!
  subroutine finish_fields
    implicit none

    deallocate(bx, by, bz)
    deallocate(bx_old, by_old, bz_old)

  end subroutine finish_fields

end module fields
