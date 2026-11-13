!-----------------------------------------------!
!> @author  YK
!! @brief   Field setting for HRMHD
!!          Evolved fields: omg (vorticity), psi (flux function),
!!          upa (parallel flow), bpa (parallel field); phi is derived from omg.
!-----------------------------------------------!
module fields
  implicit none

  public :: init_fields, finish_fields
  public :: phi, omg, psi, upa, bpa
  public :: phi_old, omg_old, psi_old, upa_old, bpa_old
  public :: nfields
  public :: iomg, ipsi, iupa, ibpa

  private

  complex(8), allocatable, dimension(:,:,:) :: phi, omg, psi, upa, bpa
  complex(8), allocatable, dimension(:,:,:) :: phi_old, omg_old, psi_old, upa_old, bpa_old
  character(100) :: init_type
  integer :: init_seed
  ! single_mode: mode numbers (kx, ky, kz in units of 2 pi/l) and amplitudes
  integer :: mode_kx, mode_ky, mode_kz
  real(8) :: amp_phi, amp_psi, amp_upa, amp_bpa

  ! Field index
  integer, parameter :: nfields = 4
  integer, parameter :: iomg = 1, ipsi = 2, iupa = 3, ibpa = 4

  !$acc declare create(phi, omg, psi, upa, bpa)
  !$acc declare create(phi_old, omg_old, psi_old, upa_old, bpa_old)
  !$acc declare copyin(nfields)
  !$acc declare copyin(iomg, ipsi, iupa, ibpa)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of fields
!-----------------------------------------------!
  subroutine init_fields
    use grid, only: nkx, nky_local, nkz
    use params, only: inputfile
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    allocate(src(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(phi    , source=src)
    allocate(omg    , source=src)
    allocate(psi    , source=src)
    allocate(upa    , source=src)
    allocate(bpa    , source=src)
    allocate(phi_old, source=src)
    allocate(omg_old, source=src)
    allocate(psi_old, source=src)
    allocate(upa_old, source=src)
    allocate(bpa_old, source=src)
    !$acc update device(phi, omg, psi, upa, bpa)
    !$acc update device(phi_old, omg_old, psi_old, upa_old, bpa_old)
    deallocate(src)

    call read_parameters(inputfile)

    if(init_type == 'zero') then
    endif
    if(init_type == 'single_mode') then
      call init_single_mode
    endif
    if(init_type == 'OT2') then
      call init_OT2
    endif
    if(init_type == 'OT3') then
      call init_OT3
    endif
    if(init_type == 'random') then
      call init_random
    endif
    if(init_type == 'restart') then
      call restart
    endif

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

    namelist /initial_condition/ init_type, init_seed, &
                                 mode_kx, mode_ky, mode_kz, &
                                 amp_phi, amp_psi, amp_upa, amp_bpa

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

    init_type = 'zero'
    ! init_seed < 0 (default): non-reproducible random IC (system clock).
    ! init_seed >= 0        : reproducible IC, required for redundant-solve
    !                         consistency across comm_m/comm_s and for regression.
    init_seed = -1
    mode_kx = 1
    mode_ky = 0
    mode_kz = 1
    amp_phi = 1.d0
    amp_psi = 1.d0
    amp_upa = 1.d0
    amp_bpa = 1.d0
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
!!          Sets the Fourier mode (mode_kx, mode_ky, mode_kz) of
!!          phi, psi, upa, bpa to amp_phi, ..., amp_bpa. The real-space
!!          field is then 2*amp*cos(k.x). The mode is given by its global
!!          wavenumber, so the result does not depend on the decomposition.
!-----------------------------------------------!
  subroutine init_single_mode
    use grid, only: nkx, nky, nky_local, nkz
    use grid, only: kprp2
    use mp, only: proc0, iproc_fft
    use MPI
    implicit none
    integer :: i, j, k, im, jm, km, j_off
    integer :: sgn, ierr

    if(proc0) then
      print '("Single mode initialization: (kx, ky, kz) mode number = (", i0, ", ", i0, ", ", i0, ")")', &
            mode_kx, mode_ky, mode_kz
    endif

    ! The r2c axis is kz, so only kz >= 0 is stored. A mode given with kz < 0
    ! is represented by its -k partner (the amplitudes are real).
    sgn = 1
    if (mode_kz < 0) sgn = -1
    im = modulo(sgn*mode_kx, nkx) + 1
    jm = modulo(sgn*mode_ky, nky) + 1
    km = sgn*mode_kz + 1
    if (km > nkz) then
      if (proc0) print *, 'single_mode: |mode_kz| exceeds the grid'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif

    ! Fields are zero on the host and the device at this point, so the mode is
    ! set on the host and copied over.
    j_off = nky_local*iproc_fft
    call set_mode(im, jm)
    ! On the kz = 0 plane the -k partner is stored too; it carries the complex
    ! conjugate, which equals the amplitude since the amplitudes are real.
    if (km == 1) call set_mode(modulo(-sgn*mode_kx, nkx) + 1, modulo(-sgn*mode_ky, nky) + 1)
    !$acc update device(phi, psi, upa, bpa)

    !$acc data present(phi, omg, psi, upa, bpa, phi_old, omg_old, psi_old, upa_old, bpa_old, kprp2)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          omg(k, j, i) = phi(k, j, i)*(-kprp2(k, j, i))
          phi_old(k, j, i) = phi(k, j, i)
          omg_old(k, j, i) = omg(k, j, i)
          psi_old(k, j, i) = psi(k, j, i)
          upa_old(k, j, i) = upa(k, j, i)
          bpa_old(k, j, i) = bpa(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data

  contains

    subroutine set_mode(ii, jj)
      implicit none
      integer, intent(in) :: ii, jj

      if (jj > j_off .and. jj <= j_off + nky_local) then
        phi(km, jj - j_off, ii) = amp_phi
        psi(km, jj - j_off, ii) = amp_psi
        upa(km, jj - j_off, ii) = amp_upa
        bpa(km, jj - j_off, ii) = amp_bpa
      endif
    end subroutine set_mode

  end subroutine init_single_mode


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
    real(8), allocatable, dimension(:,:,:) :: phi_r, psi_r, upa_r, bpa_r
    integer :: seedsize
    integer, allocatable :: seed(:)
    integer :: i, j, k

    if(proc0) then
      print *, 'Random initialization'
    endif

    allocate(phi_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(psi_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(upa_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(bpa_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(phi_r, psi_r, upa_r, bpa_r)

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
    call random_number(upa_r)
    call random_number(bpa_r)

    !$acc update device(phi_r, psi_r, upa_r, bpa_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             compute r2c transform           v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !$acc data present(phi_r, psi_r, upa_r, bpa_r, phi, psi, upa, bpa)
    call ftran_r2c(phi_r, phi)
    call ftran_r2c(psi_r, psi)
    call ftran_r2c(upa_r, upa)
    call ftran_r2c(bpa_r, bpa)
    !$acc end data

    !$acc data present(phi, omg, psi, upa, bpa, kx, ky, kz, kprp2)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi(k,j,i)/ntot
          psi(k,j,i) = psi(k,j,i)/ntot
          upa(k,j,i) = upa(k,j,i)/ntot
          bpa(k,j,i) = bpa(k,j,i)/ntot
          if (kx(i)**2 + ky(j)**2 + kz(k)**2 /= 0.d0) then
            phi(k, j, i) = phi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            psi(k, j, i) = psi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            upa(k, j, i) = upa(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            bpa(k, j, i) = bpa(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
          endif
          if (kx(i)**2 + ky(j)**2 + kz(k)**2 >= (0.8*maxval(kx))**2) then
            phi(k, j, i) = 0.d0
            psi(k, j, i) = 0.d0
            upa(k, j, i) = 0.d0
            bpa(k, j, i) = 0.d0
          endif
          omg(k, j, i) = phi(k, j, i)*(-kprp2(k, j, i))
        enddo
      enddo
    enddo
    !$acc end data

    call copy_to_old

    !$acc exit data delete(phi_r, psi_r, upa_r, bpa_r)
    deallocate(phi_r, psi_r, upa_r, bpa_r)

  end subroutine init_random


!-----------------------------------------------!
!> @author  YK
!! @brief   2D Orszag Tang problem initialization
!!          (upa = bpa = 0; the Hall terms generate them)
!-----------------------------------------------!
  subroutine init_OT2
    use grid, only: xx, yy, kprp2
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0
    use cuFFTmp, only: ftran_r2c
    implicit none
    real(8), allocatable, dimension(:,:,:) :: phi_r, psi_r
    integer :: i, j, k

    if(proc0) then
      print *, 'OT2 initialization'
    endif

    allocate(phi_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(psi_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(phi_r, psi_r)

    !$acc parallel loop collapse(3) present(phi_r, psi_r, xx, yy)
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

    !$acc data present(phi, omg, psi, upa, bpa, kprp2)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi(k,j,i)/ntot
          psi(k,j,i) = psi(k,j,i)/ntot
          upa(k,j,i) = 0.d0
          bpa(k,j,i) = 0.d0
          omg(k,j,i) = phi(k,j,i)*(-kprp2(k,j,i))
        enddo
      enddo
    enddo
    !$acc end data

    call copy_to_old

    !$acc exit data delete(phi_r, psi_r)
    deallocate(phi_r, psi_r)

  end subroutine init_OT2


!-----------------------------------------------!
!> @author  YK
!! @brief   3D Orszag Tang problem initialization
!!          (upa = bpa = 0; the Hall terms generate them)
!-----------------------------------------------!
  subroutine init_OT3
    use grid, only: lx, ly, lz, xx, yy, zz, kprp2
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0
    use params, only: pi
    use cuFFTmp, only: ftran_r2c
    implicit none
    real(8), allocatable, dimension(:,:,:) :: phi_r, psi_r
    real(8) :: x0, y0, z0
    integer :: i, j, k

    if(proc0) then
      print *, 'OT3 initialization'
    endif

    x0 = lx/(2.d0*pi)
    y0 = ly/(2.d0*pi)
    z0 = lz/(2.d0*pi)

    allocate(phi_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(psi_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(phi_r, psi_r)

    ! Same fields as Calliope's HRMHD init_OT3
    ! (calliope_dev/src/model/HRMHD/fields.F90), including the yy(j)/x0 in the
    ! second psi term, where y0 was evidently intended; it is kept so that the
    ! two codes start from identical fields (the comparison runs with lx = ly).
    ! Calliope's version then overwrites psi with the transform of upa_r = 0;
    ! that is not reproduced here.
    !
    ! zz is dimensioned nlz, so the fill must stop at nlz; the remaining
    ! nlz_padded - nlz slots are the in-place r2c padding and are zeroed.
    !$acc parallel loop collapse(3) present(phi_r, psi_r, xx, yy, zz)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          if (k <= nlz) then
            phi_r(k, j, i) = -(cos(xx(i)/x0 + 1.4d0 + zz(k)/z0) &
                             + cos(yy(j)/y0 + 0.5d0 + zz(k)/z0))
            psi_r(k, j, i) = -(0.5d0*cos(2.d0*(xx(i)/x0) + 2.3d0 + zz(k)/z0) &
                             + cos(yy(j)/x0 + 4.1d0 + zz(k)/z0))
          else
            phi_r(k, j, i) = 0.d0
            psi_r(k, j, i) = 0.d0
          endif
        end do
      end do
    end do

    ! compute r2c transform
    call ftran_r2c(phi_r, phi)
    call ftran_r2c(psi_r, psi)

    !$acc data present(phi, omg, psi, upa, bpa, kprp2)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi(k,j,i)/ntot
          psi(k,j,i) = psi(k,j,i)/ntot
          upa(k,j,i) = 0.d0
          bpa(k,j,i) = 0.d0
          omg(k,j,i) = phi(k,j,i)*(-kprp2(k,j,i))
        enddo
      enddo
    enddo
    !$acc end data

    call copy_to_old

    !$acc exit data delete(phi_r, psi_r)
    deallocate(phi_r, psi_r)

  end subroutine init_OT3


!-----------------------------------------------!
!> @author  YK
!! @brief   *_old = * on the device
!-----------------------------------------------!
  subroutine copy_to_old
    use grid, only: nkx, nky_local, nkz
    implicit none
    integer :: i, j, k

    !$acc data present(phi, omg, psi, upa, bpa, phi_old, omg_old, psi_old, upa_old, bpa_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi_old(k, j, i) = phi(k, j, i)
          omg_old(k, j, i) = omg(k, j, i)
          psi_old(k, j, i) = psi(k, j, i)
          upa_old(k, j, i) = upa(k, j, i)
          bpa_old(k, j, i) = bpa(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine copy_to_old


!-----------------------------------------------!
!> @author  YK
!! @brief   Restart
!-----------------------------------------------!
  subroutine restart
    use mp, only: proc0, iproc_fft, comm_fft
    use time, only: tt, dt
    use grid, only: nkx, nky, nky_local, nkz
    use params, only: restart_dir
    use file, only: open_input_file, close_file
    use mpiio, only: mpiio_read_one
    use shearing_box, only: tsc
    use MPI
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
    ! Offset by iproc_fft (field is decomposed along comm_fft only). Reads are
    ! NOT gated: every comm_fft group reads the same file so each redundant
    ! group loads an identical initial condition.
    starts(2) = nky_local*iproc_fft
    starts(3) = 0

    call mpiio_read_one(phi, sizes, subsizes, starts, trim(restart_dir)//'phi.dat', comm_fft)
    call mpiio_read_one(omg, sizes, subsizes, starts, trim(restart_dir)//'omg.dat', comm_fft)
    call mpiio_read_one(psi, sizes, subsizes, starts, trim(restart_dir)//'psi.dat', comm_fft)
    call mpiio_read_one(upa, sizes, subsizes, starts, trim(restart_dir)//'upa.dat', comm_fft)
    call mpiio_read_one(bpa, sizes, subsizes, starts, trim(restart_dir)//'bpa.dat', comm_fft)

    phi_old = phi
    omg_old = omg
    psi_old = psi
    upa_old = upa
    bpa_old = bpa

    !$acc update device(phi, omg, psi, upa, bpa)
    !$acc update device(phi_old, omg_old, psi_old, upa_old, bpa_old)
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @brief   Finalization of Fields
!-----------------------------------------------!
  subroutine finish_fields
    implicit none

    deallocate(phi, omg, psi, upa, bpa)
    deallocate(phi_old, omg_old, psi_old, upa_old, bpa_old)

  end subroutine finish_fields

end module fields
