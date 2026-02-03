!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics for RMHD
!-----------------------------------------------!
module diagnostics_common
  implicit none

  public :: write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, write_intvl_SF2, SF2_nsample
  public :: write_intvl_nltrans
  public :: series_output, n_series_modes, series_modes
  public :: init_polar_spectrum_2d, init_polar_spectrum_3d, init_series_modes, finish_diagnostics
  public :: get_polar_spectrum_2d, get_polar_spectrum_3d, write_polar_spectrum_2d_in_3d, cut_2d_r, cut_2d_k
  public :: nkpolar, kpbin
  public :: nkpolar_log, kpolar_log_sep, kpbin_log
  public :: series_modes_unit
  public :: read_parameters

  private

  real(8)              :: write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, write_intvl_SF2, SF2_nsample
  real(8)              :: write_intvl_nltrans
  integer              :: n_series_modes
  integer, allocatable :: series_modes(:, :)
  logical              :: series_output = .false.

  integer :: nkpolar !    = min(nkx, nky) for 2D bin
                     ! or = min(nkx, nky, nkz) for 3D bin
  real(8), dimension(:), allocatable :: kpbin

  ! kpbin_log = 0, k0, A*k0, A^2*k0, ..., A^nkpolar_log*k0 ~= k_max
  ! This A is set in inputfile
  integer :: nkpolar_log
  real(8) :: kpolar_log_sep
  real(8), dimension(:), allocatable :: kpbin_log

  integer :: series_modes_unit

  !$acc declare create(nkpolar, kpbin, nkpolar_log, kpolar_log_sep, kpbin_log)
  !$acc declare create(series_modes)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Read inputfile for various parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none
    integer :: unit
    
    character(len=100), intent(in) :: filename
    integer :: output_modes(100)
    integer  :: i, ierr

    namelist /diagnostics_parameters/ write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, &
                                      write_intvl_SF2, SF2_nsample, write_intvl_nltrans, output_modes, kpolar_log_sep

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    write_intvl         = 0.d0
    write_intvl_2D      = 0.d0
    write_intvl_3D      = 0.d0
    write_intvl_kpar    = 0.d0
    write_intvl_SF2     = 0.d0
    SF2_nsample         = 100000
    write_intvl_nltrans = 0.d0
    output_modes(:)  = 0
    kpolar_log_sep   = 2.d0**(1.d0/3.d0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=diagnostics_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading diagnostics_parameters failed"
    close(unit)

    if (write_intvl         == 0.d0) write_intvl         = -dble(huge(1))
    if (write_intvl_2D      == 0.d0) write_intvl_2D      = -dble(huge(1))
    if (write_intvl_3D      == 0.d0) write_intvl_3D      = -dble(huge(1))
    if (write_intvl_kpar    == 0.d0) write_intvl_kpar    = -dble(huge(1))
    if (write_intvl_SF2     == 0.d0) write_intvl_SF2     = -dble(huge(1))
    if (write_intvl_nltrans == 0.d0) write_intvl_nltrans = -dble(huge(1))

    ! Time series output of modes
    if(count(output_modes /= 0) > 0) then
      if(mod(count(output_modes /= 0), 3) /= 0) then
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        print *, '!              Error!              !'
        print *, '!  output modes must be multiple   !'
        print *, '!  of 3; set of (ikx, iky, ikz)    !'
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        stop
      endif

      n_series_modes = count(output_modes /= 0) / 3
      allocate(series_modes(n_series_modes, 3))

      do i = 1, n_series_modes
        series_modes(i, 1) = output_modes((i-1)*3 + 1)
        series_modes(i, 2) = output_modes((i-1)*3 + 2)
        series_modes(i, 3) = output_modes((i-1)*3 + 3)
      enddo

      series_output = .true.
    endif

    !$acc update device(series_modes)
  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of polar spectrum
!!          variables, based on AGK subroutine
!-----------------------------------------------!
  ! bin over (x, y)
  subroutine init_polar_spectrum_2d
    use grid, only: kx, ky_global
    use params, only: pi
    implicit none
    real(8) :: dkp
    integer :: i

    dkp = max(kx(2), ky_global(2))
    nkpolar = int(max(maxval(kx), maxval(ky_global))/dkp)

    allocate (kpbin(1:nkpolar))

    do i = 1, nkpolar
      kpbin(i) = (i - 1)*dkp
    enddo

    !$acc update device(nkpolar, kpbin)
  end subroutine init_polar_spectrum_2d

  ! bin over (x, y, z)
  subroutine init_polar_spectrum_3d
    use grid, only: kx, ky_global, kz
    use params, only: pi
    implicit none
    real(8) :: dkp, k0, kmax
    integer :: i

    ! set kpbin
    dkp = max(kx(2), ky_global(2), kz(2))
    nkpolar = int(max(maxval(kx), maxval(ky_global), maxval(kz))/dkp)

    allocate (kpbin(1:nkpolar))

    do i = 1, nkpolar
      kpbin(i) = (i - 1)*dkp
    enddo

    ! set kpbin_log
    k0   = min(kx(2), ky_global(2), kz(2))
    kmax = max(maxval(kx), maxval(ky_global), maxval(kz))
    nkpolar_log = floor( dlog(kmax/k0)/dlog(kpolar_log_sep) ) + 2

    allocate (kpbin_log(1:nkpolar_log))

    kpbin_log(1) = 0.d0
    do i = 1, nkpolar_log - 1
      kpbin_log(i+1) = k0*kpolar_log_sep**(i-1)
    enddo

    !$acc update device(nkpolar, kpbin, nkpolar_log, kpolar_log_sep, kpbin_log)
  end subroutine init_polar_spectrum_3d


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization series output of modes
!-----------------------------------------------!
  subroutine init_series_modes
    use file, only: open_output_file
    use params, only: runname
    use mp, only: proc0
    implicit none

    if(series_output .and. proc0) call open_output_file (series_modes_unit, trim(runname)//'.modes.out')

  end subroutine init_series_modes


!-----------------------------------------------!
!> @author  YK
!! @brief   Return kpolar binned and log averaged
!           spectrum, based on AGK subroutine
!-----------------------------------------------!
  ! bin over (x, y)
  subroutine get_polar_spectrum_2d(ee, ebin)
    use params, only: pi
    use mp, only: sum_reduce
    use grid, only: nkx, nky_local, nkz
    use grid, only: kx, ky
    implicit none
    real(8), dimension (:,:,:), intent (in)  :: ee   !variable ee by (kx,ky) mode
    real(8), dimension (1:nkpolar, nkz), intent (out) :: ebin 
    real(8) :: k2
    integer :: idx, i, j, k

    ! Initialize ebin on device
    !$acc parallel loop collapse(2) default(present)
    do k = 1, nkz
      do idx = 1, nkpolar
        ebin(idx, k) = 0.d0
      end do
    end do
    !$acc end parallel loop

    ! Loop through all modes and bin by polar wavenumber
    ! Use atomic update for thread-safe accumulation into ebin
    !$acc parallel loop collapse(3) gang vector default(present) private(k2, idx)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          k2 = kx(i)**2 + ky(j)**2
          ! Find the appropriate bin and accumulate
          do idx = 1, nkpolar - 1
            if (k2 >= kpbin(idx)**2 .and. k2 < kpbin(idx+1)**2) then
              !$acc atomic update
              ebin(idx, k) = ebin(idx, k) + ee(k, j, i)
              exit  ! Found the bin, no need to continue
            endif
          enddo
        enddo
      enddo
    enddo
    !$acc end parallel loop

    ! Transfer ebin to host for MPI reduction
    !$acc update host(ebin)

    ! MPI reduction across all processes
    call sum_reduce(ebin, 0)

    ! Transfer reduced result back to device
    !$acc update device(ebin)
  end subroutine get_polar_spectrum_2d

  ! bin over (x, y, z)
  subroutine get_polar_spectrum_3d(ee, ebin)
    use params, only: pi
    use mp, only: sum_reduce
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: shear, q
    use shearing_box, only: shear_flg, tsc
    implicit none
    real(8), dimension (:,:,:), intent (in)  :: ee   ! Energy density by (kz, ky, kx) mode
    real(8), dimension (1:nkpolar)  , intent (out) :: ebin 
    real(8) :: kxt, k2
    integer :: idx, i, j, k

    ! Initialize ebin on device
    !$acc parallel loop default(present)
    do idx = 1, nkpolar
      ebin(idx) = 0.d0
    end do
    !$acc end parallel loop

    ! Loop through all modes and bin by polar wavenumber
    ! Use atomic update for thread-safe accumulation into ebin
    !$acc parallel loop collapse(3) gang vector default(present) private(kxt, k2, idx)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if(shear) then
            kxt = kx(i) + q*shear_flg*tsc*ky(j)
          else
            kxt = kx(i)
          endif
          k2 = kxt**2 + ky(j)**2 + kz(k)**2

          ! Find the appropriate bin and accumulate
          do idx = 1, nkpolar - 1
            if (k2 >= kpbin(idx)**2 .and. k2 < kpbin(idx+1)**2) then
              !$acc atomic update
              ebin(idx) = ebin(idx) + ee(k, j, i)
              exit  ! Found the bin, no need to continue
            endif
          enddo
        enddo
      enddo
    enddo
    !$acc end parallel loop

    ! Transfer ebin to host for MPI reduction
    !$acc update host(ebin)

    ! MPI reduction across all processes
    call sum_reduce(ebin, 0)

    ! Transfer reduced result back to device
    !$acc update device(ebin)

  end subroutine get_polar_spectrum_3d

  ! bin over (x, y) leaving z
  !       or (x, z) leaving y
  ! no x leaving because kx will be kxt when shear is on.
  subroutine write_polar_spectrum_2d_in_3d(ee, direction, unit)
    use mp, only: proc0, sum_reduce
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky, nky_local, nkz
    use params, only: shear, q
    use shearing_box, only: shear_flg, tsc
    implicit none
    real(8), dimension (:,:,:), intent (in)  :: ee   !variable ee by (kx,ky) mode
    character(len=1), intent(in) :: direction
    integer, intent(in) :: unit

    real(8), dimension (1:nkpolar) :: ebin 
    real(8) :: kxt
    integer :: idx, i, j, k

    if(trim(direction) == 'z') then
      do k = 1, nkz
        ebin = 0.d0

        do idx = 1, nkpolar - 1
          do i = 1, nkx
            do j = 1, nky_local
              if(shear) then
                kxt = kx(i) + q*shear_flg*tsc*ky(j)
              else
                kxt = kx(i)
              endif

              if (      kxt**2 + ky(j)**2 >= kpbin(idx)**2 &
                  .and. kxt**2 + ky(j)**2 <  kpbin(idx+1)**2) then
                ebin(idx) = ebin(idx) + ee(k, j, i)
              endif

            enddo
          enddo
        enddo
        call sum_reduce(ebin, 0)

        if(proc0) write (unit=unit) ebin

      enddo
    elseif(direction == 'y') then
      do j = 1, nky
        ebin = 0.d0

        do idx = 1, nkpolar - 1
          if(j >= 1 .and. j <= nky_local) then
            
            do k = 1, nkz
              do i = 1, nkx
                if(shear) then
                  kxt = kx(i) + q*shear_flg*tsc*ky(j)
                else
                  kxt = kx(i)
                endif

                if (      kxt**2 + kz(k)**2 >= kpbin(idx)**2 &
                    .and. kxt**2 + kz(k)**2 <  kpbin(idx+1)**2) then
                  ebin(idx) = ebin(idx) + ee(k, j, i)
                endif

              enddo
            enddo

          endif
        enddo
        call sum_reduce(ebin, 0)

        if(proc0) write (unit=unit) ebin

      enddo
    endif

    if(proc0) call flush(unit) 

  end subroutine write_polar_spectrum_2d_in_3d


!-----------------------------------------------!
!> @author  YK
!! @brief   2D cut of real field 
!!                           along z=0, x=0, y=0
!-----------------------------------------------!
  subroutine cut_2d_r(f, fr_z0, fr_x0, fr_y0)
    use grid, only: xx, yy, zz
    use grid, only: nlx_local, nly, nlz
    implicit none
    real(8), dimension (:,:,:), intent(in) :: f
    real(8), dimension (1:nlx_local, 1:nly), intent(out) :: fr_z0
    real(8), dimension (1:nly      , 1:nlz), intent(out) :: fr_x0
    real(8), dimension (1:nlx_local, 1:nlz), intent(out) :: fr_y0

    integer :: i, j, k

    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          if(zz(k) == 0.d0) fr_z0(i, j) = f(k, j, i)
          if(xx(i) == 0.d0) fr_x0(j, k) = f(k, j, i)
          if(yy(j) == 0.d0) fr_y0(i, k) = f(k, j, i)
        end do
      end do
    end do

  end subroutine cut_2d_r


!-----------------------------------------------!
!> @author  YK
!! @brief   2D cut of spectral field 
!!                           along kz=0, kx=0, ky=0
!-----------------------------------------------!
  subroutine cut_2d_k(fk, fk_kxy, fk_kyz, fk_kxz)
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    implicit none
    real(8), dimension (:,:,:), intent(in) :: fk
    real(8), dimension (1:nkx      , 1:nky_local), intent(out) :: fk_kxy
    real(8), dimension (1:nky_local, 1:nkz      ), intent(out) :: fk_kyz
    real(8), dimension (1:nkx      , 1:nkz      ), intent(out) :: fk_kxz
    real(8), allocatable :: message(:)

    integer :: i, j, k
    integer :: proc_j1, proc_k1

    fk_kxy(:, :) = 0.d0
    fk_kyz(:, :) = 0.d0
    fk_kxz(:, :) = 0.d0

    do i =1, nkx
      do j = 1, nky_local
        do k = 1, nkz
         if(kz(k) == 0.d0) fk_kxy(i, j) = fk_kxy(i, j) + fk(i, k, j)
         if(kx(i) == 0.d0) fk_kyz(j, k) = fk_kyz(j, k) + fk(i, k, j)
         if(ky(j) == 0.d0) fk_kxz(i, k) = fk_kxz(i, k) + fk(i, k, j)
        end do
      end do
    end do

  end subroutine cut_2d_k


!-----------------------------------------------!
!> @author  YK
!! @brief   Finalization of diagnostics
!-----------------------------------------------!
  subroutine finish_diagnostics
    use io, only: finish_io
    implicit none

    deallocate(kpbin)

    call finish_io
  end subroutine finish_diagnostics

end module diagnostics_common

