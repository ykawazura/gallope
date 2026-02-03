!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Diagnostics for MHD_INCOMP
!-----------------------------------------------!
module diagnostics
  use diagnostics_common
  implicit none

  public :: init_diagnostics, finish_diagnostics
  public :: loop_diagnostics, loop_diagnostics_2D, loop_diagnostics_kpar, loop_diagnostics_SF2
  public :: loop_diagnostics_nltrans

  private

  integer  :: nl, nsamp_r0, nsamp_ang
  real(8), allocatable :: ll(:), lpar(:), lper(:)
  integer  :: unit_u2kxy_vs_kz, unit_u2kxz_vs_ky
  integer  :: unit_b2kxy_vs_kz, unit_b2kxz_vs_ky
contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of diagnostics
!-----------------------------------------------!
  subroutine init_diagnostics
    use params, only: inputfile
    use diagnostics_common, only: read_parameters
    use diagnostics_common, only: init_polar_spectrum_2d, init_polar_spectrum_3d
    use diagnostics_common, only: init_series_modes
    use io, only: init_io 
    use grid, only: nlz
    implicit none

    call read_parameters(inputfile)

    if(nlz == 2) then
      call init_polar_spectrum_2d
    else
      call init_polar_spectrum_3d
    endif
    call init_SF2
    call init_io(nkpolar, kpbin, nkpolar_log, kpbin_log, nl, lpar, lper)
    call init_series_modes
    call set_unit_for_polar_spectrum_3d_in_2d
  end subroutine init_diagnostics

  subroutine set_unit_for_polar_spectrum_3d_in_2d
    use file, only: open_output_file_binary
    use mp, only: proc0
    implicit none

    if(proc0) then
      call open_output_file_binary (unit_u2kxy_vs_kz, 'out2d/u2_kxy_vs_kz.dat')
      call open_output_file_binary (unit_u2kxz_vs_ky, 'out2d/u2_kxz_vs_ky.dat')
      call open_output_file_binary (unit_b2kxy_vs_kz, 'out2d/b2_kxy_vs_kz.dat')
      call open_output_file_binary (unit_b2kxz_vs_ky, 'out2d/b2_kxz_vs_ky.dat')
    endif
  end subroutine set_unit_for_polar_spectrum_3d_in_2d


!-----------------------------------------------!
!> @author  YK
!! @date    15 Apr 2020
!! @brief   Diagnostics in loop
!-----------------------------------------------!
  subroutine loop_diagnostics
    use io, only: loop_io
    use mp, only: proc0
    use grid, only: kx, ky, kz, k2_max
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx, nly, nlz_padded
    use grid, only: dlx, dly, dlz
    use grid, only:  lx,  ly,  lz
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use fields, only: ux_old, uy_old, uz_old
    use fields, only: bx_old, by_old, bz_old
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, nu, nu_h, nu_h_exp, eta, eta_h, eta_h_exp, shear, q
    use force, only: fux, fuy, fuz, fux_old, fuy_old, fuz_old
    use force, only: fbx, fby, fbz, fbx_old, fby_old, fbz_old
    use shearing_box, only: shear_flg, tsc, nremap, k2t
    implicit none
    integer :: i, j, k

    real(8), allocatable, dimension(:,:,:) :: u2, ux2, uy2, uz2
    real(8), allocatable, dimension(:,:,:) :: b2, bx2, by2, bz2
    real(8), allocatable, dimension(:,:,:) :: u2old, b2old
    real(8), allocatable, dimension(:,:,:) :: u2dissip, b2dissip
    real(8), allocatable, dimension(:,:,:) :: p_ext_ene, p_ext_xhl, p_re, p_ma
    real(8), allocatable, dimension(:,:,:) :: zp2, zm2
    real(8), allocatable, dimension(:,:,:) :: src

    real(8) :: u2_sum, b2_sum
    real(8) :: u2dot_sum, b2dot_sum
    real(8) :: u2dissip_sum, b2dissip_sum
    real(8) :: p_ext_ene_sum, p_ext_xhl_sum, p_re_sum, p_ma_sum
    real(8) :: zp2_sum, zm2_sum
    real(8) :: bx0, by0, bz0 ! mean magnetic field

    real(8), dimension(:), allocatable :: u2_bin, ux2_bin, uy2_bin, uz2_bin
    real(8), dimension(:), allocatable :: b2_bin, bx2_bin, by2_bin, bz2_bin
    real(8), dimension(:), allocatable :: zp2_bin, zm2_bin
    real(8), dimension(:), allocatable :: u2dissip_bin, b2dissip_bin
    real(8), dimension(:), allocatable :: p_re_bin, p_ma_bin
    complex(8) ::  ux_mid,  uy_mid,  uz_mid
    complex(8) ::  bx_mid,  by_mid,  bz_mid
    complex(8) :: fux_mid, fuy_mid, fuz_mid
    complex(8) :: fbx_mid, fby_mid, fbz_mid

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(src(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(u2      , source=src)
    allocate(ux2     , source=src)
    allocate(uy2     , source=src)
    allocate(uz2     , source=src)
    allocate(b2      , source=src)
    allocate(bx2     , source=src)
    allocate(by2     , source=src)
    allocate(bz2     , source=src)
    allocate(u2old   , source=src)
    allocate(b2old   , source=src)
    allocate(u2dissip, source=src)
    allocate(b2dissip, source=src)
    allocate(p_ext_ene, source=src)
    allocate(p_ext_xhl, source=src)
    allocate(p_re    , source=src)
    allocate(p_ma    , source=src)
    allocate(zp2     , source=src)
    allocate(zm2     , source=src)
    deallocate(src)

    allocate(  u2_bin(1:nkpolar), source=0.d0)
    allocate( ux2_bin(1:nkpolar), source=0.d0)
    allocate( uy2_bin(1:nkpolar), source=0.d0)
    allocate( uz2_bin(1:nkpolar), source=0.d0)
    allocate(  b2_bin(1:nkpolar), source=0.d0)
    allocate( bx2_bin(1:nkpolar), source=0.d0)
    allocate( by2_bin(1:nkpolar), source=0.d0)
    allocate( bz2_bin(1:nkpolar), source=0.d0)
    allocate( zp2_bin(1:nkpolar), source=0.d0)
    allocate( zm2_bin(1:nkpolar), source=0.d0)
    allocate(u2dissip_bin(1:nkpolar), source=0.d0)
    allocate(b2dissip_bin(1:nkpolar), source=0.d0)
    allocate(p_re_bin(1:nkpolar), source=0.d0)
    allocate(p_ma_bin(1:nkpolar), source=0.d0)

    ! Create work arrays on device
    !$acc enter data create(u2, ux2, uy2, uz2, b2, bx2, by2, bz2)
    !$acc enter data create(u2old, b2old, u2dissip, b2dissip)
    !$acc enter data create(p_ext_ene, p_ext_xhl, p_re, p_ma, zp2, zm2)
    !$acc enter data create(u2_bin, ux2_bin, uy2_bin, uz2_bin)
    !$acc enter data create(b2_bin, bx2_bin, by2_bin, bz2_bin)
    !$acc enter data create(zp2_bin, zm2_bin, u2dissip_bin, b2dissip_bin)
    !$acc enter data create(p_re_bin, p_ma_bin)

    ! Initialize reduction variables
    bx0           = 0.d0
    by0           = 0.d0
    bz0           = 0.d0
    u2_sum        = 0.d0
    b2_sum        = 0.d0
    u2dot_sum     = 0.d0
    b2dot_sum     = 0.d0
    u2dissip_sum  = 0.d0
    b2dissip_sum  = 0.d0
    p_ext_ene_sum = 0.d0
    p_ext_xhl_sum = 0.d0
    p_re_sum      = 0.d0
    p_ma_sum      = 0.d0
    zp2_sum       = 0.d0
    zm2_sum       = 0.d0

    !$acc parallel loop collapse(3) gang vector default(present) &
    !$acc& private(ux_mid, uy_mid, uz_mid, bx_mid, by_mid, bz_mid) &
    !$acc& private(fux_mid, fuy_mid, fuz_mid, fbx_mid, fby_mid, fbz_mid) &
    !$acc& reduction(+:u2_sum, b2_sum, u2dot_sum, b2dot_sum) &
    !$acc& reduction(+:u2dissip_sum, b2dissip_sum) &
    !$acc& reduction(+:p_ext_ene_sum, p_ext_xhl_sum, p_re_sum, p_ma_sum) &
    !$acc& reduction(+:zp2_sum, zm2_sum)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
            bx0 = abs(bx(k,j,i))
            by0 = abs(by(k,j,i))
            bz0 = abs(bz(k,j,i))
          endif
          u2   (k, j, i) = 0.5d0*(abs(ux    (k, j, i))**2 + abs(uy    (k, j, i))**2 + abs(uz    (k, j, i))**2)
          ux2  (k, j, i) = 0.5d0*(abs(ux    (k, j, i))**2)
          uy2  (k, j, i) = 0.5d0*(abs(uy    (k, j, i))**2)
          uz2  (k, j, i) = 0.5d0*(abs(uz    (k, j, i))**2)
          b2   (k, j, i) = 0.5d0*(abs(bx    (k, j, i))**2 + abs(by    (k, j, i))**2 + abs(bz    (k, j, i))**2)
          bx2  (k, j, i) = 0.5d0*(abs(bx    (k, j, i))**2)
          by2  (k, j, i) = 0.5d0*(abs(by    (k, j, i))**2)
          bz2  (k, j, i) = 0.5d0*(abs(bz    (k, j, i))**2)
          u2old(k, j, i) = 0.5d0*(abs(ux_old(k, j, i))**2 + abs(uy_old(k, j, i))**2 + abs(uz_old(k, j, i))**2)
          b2old(k, j, i) = 0.5d0*(abs(bx_old(k, j, i))**2 + abs(by_old(k, j, i))**2 + abs(bz_old(k, j, i))**2)
          zp2  (k, j, i) = abs(ux(k,j,i) + bx(k,j,i))**2 + abs(uy(k,j,i) + by(k,j,i))**2 + abs(uz(k,j,i) + bz(k,j,i))**2
          zm2  (k, j, i) = abs(ux(k,j,i) - bx(k,j,i))**2 + abs(uy(k,j,i) - by(k,j,i))**2 + abs(uz(k,j,i) - bz(k,j,i))**2

           ux_mid  = 0.5d0*( ux(k, j, i) +  ux_old(k, j, i))
           uy_mid  = 0.5d0*( uy(k, j, i) +  uy_old(k, j, i))
           uz_mid  = 0.5d0*( uz(k, j, i) +  uz_old(k, j, i))
           bx_mid  = 0.5d0*( bx(k, j, i) +  bx_old(k, j, i))
           by_mid  = 0.5d0*( by(k, j, i) +  by_old(k, j, i))
           bz_mid  = 0.5d0*( bz(k, j, i) +  bz_old(k, j, i))
          fux_mid  = 0.5d0*(fux(k, j, i) + fux_old(k, j, i))
          fuy_mid  = 0.5d0*(fuy(k, j, i) + fuy_old(k, j, i))
          fuz_mid  = 0.5d0*(fuz(k, j, i) + fuz_old(k, j, i))
          fbx_mid  = 0.5d0*(fbx(k, j, i) + fbx_old(k, j, i))
          fby_mid  = 0.5d0*(fby(k, j, i) + fby_old(k, j, i))
          fbz_mid  = 0.5d0*(fbz(k, j, i) + fbz_old(k, j, i))

          u2dissip(k, j, i) = (nu *(k2t(k, j, i)/k2_max) + nu_h *(k2t(k, j, i)/k2_max)**nu_h_exp ) &
                                *(abs(ux_mid)**2 + abs(uy_mid)**2 + abs(uz_mid)**2)
          b2dissip(k, j, i) = (eta*(k2t(k, j, i)/k2_max) + eta_h*(k2t(k, j, i)/k2_max)**eta_h_exp) &
                                *(abs(bx_mid)**2 + abs(by_mid)**2 + abs(bz_mid)**2)
          p_ext_ene(k, j, i) = 0.5d0*( &
                                  (fux_mid*conjg(ux_mid) + conjg(fux_mid)*ux_mid) &
                                + (fuy_mid*conjg(uy_mid) + conjg(fuy_mid)*uy_mid) &
                                + (fuz_mid*conjg(uz_mid) + conjg(fuz_mid)*uz_mid) &
                                + (fbx_mid*conjg(bx_mid) + conjg(fbx_mid)*bx_mid) &
                                + (fby_mid*conjg(by_mid) + conjg(fby_mid)*by_mid) &
                                + (fbz_mid*conjg(bz_mid) + conjg(fbz_mid)*bz_mid) &
                              )
          p_ext_xhl(k, j, i) = 0.5d0*( &
                                  (fux_mid*conjg(bx_mid) + conjg(fux_mid)*bx_mid) &
                                + (fuy_mid*conjg(by_mid) + conjg(fuy_mid)*by_mid) &
                                + (fuz_mid*conjg(bz_mid) + conjg(fuz_mid)*bz_mid) &
                                + (fbx_mid*conjg(ux_mid) + conjg(fbx_mid)*ux_mid) &
                                + (fby_mid*conjg(uy_mid) + conjg(fby_mid)*uy_mid) &
                                + (fbz_mid*conjg(uz_mid) + conjg(fbz_mid)*uz_mid) &
                              )
          p_re    (k, j, i) = + 0.5d0*q*shear_flg*(ux_mid*conjg(uy_mid) + conjg(ux_mid)*uy_mid)
          p_ma    (k, j, i) = - 0.5d0*q*shear_flg*(bx_mid*conjg(by_mid) + conjg(bx_mid)*by_mid)

          ! The reason for the following treatment for kz == 0 mode is the following. Compile it with LaTeX.
          !-----------------------------------------------------------------------------------------------------------------------------------
          ! The volume integral of a quadratic function is
          ! \int \mathrm{d}^3\mathbf{r}\, f(x,y,z)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2}
          ! \sum_{k_z = -n_{k_z}/2}^{n_{k_z}/2}|f_{k_x, k_y, k_z}|^2 = \left( \sum_{k_z = -n_{k_z}/2}^{-1}\sum_{k_x, k_y} + \sum_{k_z = 0}
          ! \sum_{k_x, k_y} + \sum_{k_z = 1}^{n_{k_z}/2}\sum_{k_x, k_y} \right) |f_{k_x, k_y, k_z}|^2          !
          ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
          !-----------------------------------------------------------------------------------------------------------------------------------
          if (k /= 1) then
            u2       (k, j, i) = 2.0d0*u2       (k, j, i)
            ux2      (k, j, i) = 2.0d0*ux2      (k, j, i)
            uy2      (k, j, i) = 2.0d0*uy2      (k, j, i)
            uz2      (k, j, i) = 2.0d0*uz2      (k, j, i)
            b2       (k, j, i) = 2.0d0*b2       (k, j, i)
            bx2      (k, j, i) = 2.0d0*bx2      (k, j, i)
            by2      (k, j, i) = 2.0d0*by2      (k, j, i)
            bz2      (k, j, i) = 2.0d0*bz2      (k, j, i)
            u2old    (k, j, i) = 2.0d0*u2old    (k, j, i)
            b2old    (k, j, i) = 2.0d0*b2old    (k, j, i)
            u2dissip (k, j, i) = 2.0d0*u2dissip (k, j, i)
            b2dissip (k, j, i) = 2.0d0*b2dissip (k, j, i)
            p_ext_ene(k, j, i) = 2.0d0*p_ext_ene(k, j, i)
            p_ext_xhl(k, j, i) = 2.0d0*p_ext_xhl(k, j, i)
            p_re     (k, j, i) = 2.0d0*p_re     (k, j, i)
            p_ma     (k, j, i) = 2.0d0*p_ma     (k, j, i)
            zp2      (k, j, i) = 2.0d0*zp2      (k, j, i)
            zm2      (k, j, i) = 2.0d0*zm2      (k, j, i)
          endif

          u2_sum        = u2_sum        + u2       (k, j, i)
          b2_sum        = b2_sum        + b2       (k, j, i)
          u2dot_sum     = u2dot_sum     + (u2(k, j, i) - u2old(k, j, i))/dt
          b2dot_sum     = b2dot_sum     + (b2(k, j, i) - b2old(k, j, i))/dt
          u2dissip_sum  = u2dissip_sum  + u2dissip (k, j, i)
          b2dissip_sum  = b2dissip_sum  + b2dissip (k, j, i)
          p_ext_ene_sum = p_ext_ene_sum + p_ext_ene(k, j, i)
          p_ext_xhl_sum = p_ext_xhl_sum + p_ext_xhl(k, j, i)
          p_re_sum      = p_re_sum      + p_re     (k, j, i)
          p_ma_sum      = p_ma_sum      + p_ma     (k, j, i)
          zp2_sum       = zp2_sum       + zp2      (k, j, i)
          zm2_sum       = zm2_sum       + zm2      (k, j, i)

        end do
      end do
    end do
    !$acc end parallel loop
    
    !$acc parallel loop collapse(3) default(present) &
    !$acc& reduction(+:bx0, by0, bz0)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
            bx0 = abs(bx(k,j,i))
            by0 = abs(by(k,j,i))
            bz0 = abs(bz(k,j,i))
          endif
        end do
      end do
    end do
    !$acc end parallel loop

    ! MPI reduction for summed quantities (operates on host scalars)
    call sum_reduce(u2_sum       , 0)
    call sum_reduce(b2_sum       , 0)
    call sum_reduce(u2dot_sum    , 0)
    call sum_reduce(b2dot_sum    , 0)
    call sum_reduce(u2dissip_sum , 0)
    call sum_reduce(b2dissip_sum , 0)
    call sum_reduce(p_ext_ene_sum, 0)
    call sum_reduce(p_ext_xhl_sum, 0)
    call sum_reduce(p_re_sum     , 0)
    call sum_reduce(p_ma_sum     , 0)
    call sum_reduce(zp2_sum      , 0)
    call sum_reduce(zm2_sum      , 0)

    ! Bin spectra over (kx, ky, kz) on device
    call get_polar_spectrum_3d(  u2,   u2_bin)
    call get_polar_spectrum_3d( ux2,  ux2_bin)
    call get_polar_spectrum_3d( uy2,  uy2_bin)
    call get_polar_spectrum_3d( uz2,  uz2_bin)
    call get_polar_spectrum_3d(  b2,   b2_bin)
    call get_polar_spectrum_3d( bx2,  bx2_bin)
    call get_polar_spectrum_3d( by2,  by2_bin)
    call get_polar_spectrum_3d( bz2,  bz2_bin)
    call get_polar_spectrum_3d( zp2,  zp2_bin)
    call get_polar_spectrum_3d( zm2,  zm2_bin)
    call get_polar_spectrum_3d(u2dissip, u2dissip_bin)
    call get_polar_spectrum_3d(b2dissip, b2dissip_bin)
    call get_polar_spectrum_3d(p_re, p_re_bin)
    call get_polar_spectrum_3d(p_ma, p_ma_bin)

    ! Transfer binned spectra from device to host for I/O
    !$acc update host(u2_bin, ux2_bin, uy2_bin, uz2_bin)
    !$acc update host(b2_bin, bx2_bin, by2_bin, bz2_bin)
    !$acc update host(zp2_bin, zm2_bin, u2dissip_bin, b2dissip_bin)
    !$acc update host(p_re_bin, p_ma_bin)
  
    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  u2_sum, b2_sum, &
                  u2dot_sum, b2dot_sum, &
                  u2dissip_sum, b2dissip_sum, &
                  p_ext_ene_sum, p_ext_xhl_sum, &
                  p_re_sum, p_ma_sum, &
                  zp2_sum, zm2_sum, &
                  bx0, by0, bz0, &
                  !
                  nkpolar, &
                  u2_bin, ux2_bin, uy2_bin, uz2_bin, &
                  b2_bin, bx2_bin, by2_bin, bz2_bin, &
                  zp2_bin, zm2_bin, &
                  u2dissip_bin, b2dissip_bin, &
                  p_re_bin, p_ma_bin &
                )

    !$acc exit data delete(u2, ux2, uy2, uz2, b2, bx2, by2, bz2)
    !$acc exit data delete(u2old, b2old, u2dissip, b2dissip)
    !$acc exit data delete(p_ext_ene, p_ext_xhl, p_re, p_ma, zp2, zm2)
    !$acc exit data delete(u2_bin, ux2_bin, uy2_bin, uz2_bin)
    !$acc exit data delete(b2_bin, bx2_bin, by2_bin, bz2_bin)
    !$acc exit data delete(zp2_bin, zm2_bin, u2dissip_bin, b2dissip_bin)
    !$acc exit data delete(p_re_bin, p_ma_bin)
    deallocate(u2)
    deallocate(ux2)
    deallocate(uy2)
    deallocate(uz2)
    deallocate(b2)
    deallocate(bx2)
    deallocate(by2)
    deallocate(bz2)
    deallocate(u2old)
    deallocate(b2old)
    deallocate(u2dissip)
    deallocate(b2dissip)
    deallocate(p_ext_ene)
    deallocate(p_ext_xhl)
    deallocate(p_ma)
    deallocate(zp2)
    deallocate(zm2)

    deallocate(  u2_bin    )
    deallocate( ux2_bin    )
    deallocate( uy2_bin    )
    deallocate( uz2_bin    )
    deallocate(  b2_bin    )
    deallocate( bx2_bin    )
    deallocate( by2_bin    )
    deallocate( bz2_bin    )
    deallocate( zp2_bin    )
    deallocate( zm2_bin    )
    deallocate(u2dissip_bin)
    deallocate(b2dissip_bin)
    deallocate(p_re_bin    )
    deallocate(p_ma_bin    )
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Diagnostics for cross section of fileds
!-----------------------------------------------!
  subroutine loop_diagnostics_2D
    use io, only: loop_io_2D
    use utils, only: curl
    use mp, only: proc0
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use cuFFTmp, only: btran_c2r
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: shear
    use shearing_box, only: to_non_shearing_coordinate, tsc, nremap
    implicit none

    complex(8), allocatable, dimension(:,:,:) :: f
    real(8)   , allocatable, dimension(:,:,:) :: fr
    real(8)   , allocatable, dimension(:,:)   :: ux_r_z0, ux_r_x0, ux_r_y0
    real(8)   , allocatable, dimension(:,:)   :: uy_r_z0, uy_r_x0, uy_r_y0
    real(8)   , allocatable, dimension(:,:)   :: uz_r_z0, uz_r_x0, uz_r_y0
    real(8)   , allocatable, dimension(:,:)   :: wx_r_z0, wx_r_x0, wx_r_y0
    real(8)   , allocatable, dimension(:,:)   :: wy_r_z0, wy_r_x0, wy_r_y0
    real(8)   , allocatable, dimension(:,:)   :: wz_r_z0, wz_r_x0, wz_r_y0
    real(8)   , allocatable, dimension(:,:)   :: bx_r_z0, bx_r_x0, bx_r_y0
    real(8)   , allocatable, dimension(:,:)   :: by_r_z0, by_r_x0, by_r_y0
    real(8)   , allocatable, dimension(:,:)   :: bz_r_z0, bz_r_x0, bz_r_y0
    real(8)   , allocatable, dimension(:,:)   :: jx_r_z0, jx_r_x0, jx_r_y0
    real(8)   , allocatable, dimension(:,:)   :: jy_r_z0, jy_r_x0, jy_r_y0
    real(8)   , allocatable, dimension(:,:)   :: jz_r_z0, jz_r_x0, jz_r_y0
    complex(8), allocatable, dimension(:,:,:) :: wx, wy, wz, jx, jy, jz

    real(8)   , allocatable, dimension(:,:,:) :: u2, b2
    real(8)   , allocatable, dimension(:,:)   :: u2_kxy, u2_kyz, u2_kxz
    real(8)   , allocatable, dimension(:,:)   :: b2_kxy, b2_kyz, b2_kxz

    real(8)   , allocatable, dimension(:,:)   :: src1, src2, src3
    complex(8), allocatable, dimension(:,:,:) :: src4

    integer :: i, j, k

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(f (nkz       , nky_local, nkx      )); f   = 0.d0
    allocate(fr(nlz_padded, nly      , nlx_local)); fr  = 0.d0
    !$acc enter data create(f )
    !$acc enter data create(fr)

    allocate(src1(nlx_local, nly), source=0.d0) 
    allocate(src2(nly      , nlz), source=0.d0) 
    allocate(src3(nlx_local, nlz), source=0.d0)
    allocate(ux_r_z0, source=src1)
    allocate(ux_r_x0, source=src2)
    allocate(ux_r_y0, source=src3)

    allocate(uy_r_z0, source=src1)
    allocate(uy_r_x0, source=src2)
    allocate(uy_r_y0, source=src3)

    allocate(uz_r_z0, source=src1)
    allocate(uz_r_x0, source=src2)
    allocate(uz_r_y0, source=src3)

    allocate(wx_r_z0, source=src1)
    allocate(wx_r_x0, source=src2)
    allocate(wx_r_y0, source=src3)

    allocate(wy_r_z0, source=src1)
    allocate(wy_r_x0, source=src2)
    allocate(wy_r_y0, source=src3)

    allocate(wz_r_z0, source=src1)
    allocate(wz_r_x0, source=src2)
    allocate(wz_r_y0, source=src3)

    allocate(bx_r_z0, source=src1)
    allocate(bx_r_x0, source=src2)
    allocate(bx_r_y0, source=src3)

    allocate(by_r_z0, source=src1)
    allocate(by_r_x0, source=src2)
    allocate(by_r_y0, source=src3)

    allocate(bz_r_z0, source=src1)
    allocate(bz_r_x0, source=src2)
    allocate(bz_r_y0, source=src3)

    allocate(jx_r_z0, source=src1)
    allocate(jx_r_x0, source=src2)
    allocate(jx_r_y0, source=src3)

    allocate(jy_r_z0, source=src1)
    allocate(jy_r_x0, source=src2)
    allocate(jy_r_y0, source=src3)

    allocate(jz_r_z0, source=src1)
    allocate(jz_r_x0, source=src2)
    allocate(jz_r_y0, source=src3)
    deallocate(src1, src2, src3)

    allocate(src4(nkz, nky_local, nkx), source=(0.d0,0.d0))
    allocate(wx, source=src4)
    allocate(wy, source=src4)
    allocate(wz, source=src4)
    allocate(jx, source=src4)
    allocate(jy, source=src4)
    allocate(jz, source=src4)
    !$acc enter data create(wx)
    !$acc enter data create(wy)
    !$acc enter data create(wz)
    !$acc enter data create(jx)
    !$acc enter data create(jy)
    !$acc enter data create(jz)
    deallocate(src4)

    allocate(u2 (nkz, nky_local, nkx), source=0.d0)
    allocate(b2 (nkz, nky_local, nkx), source=0.d0)

    allocate(src1(nkx      , nky_local), source=0.d0); 
    allocate(src2(nky_local, nkz      ), source=0.d0); 
    allocate(src3(nkx      , nkz      ), source=0.d0)
    allocate(u2_kxy, source=src1)
    allocate(u2_kyz, source=src2)
    allocate(u2_kxz, source=src3)
    allocate(b2_kxy, source=src1)
    allocate(b2_kyz, source=src2)
    allocate(b2_kxz, source=src3)
    deallocate(src1, src2, src3)

    !vvvvvvvvvvvvvvvvvv         2D cut of fields          vvvvvvvvvvvvvvvvvv!
    call curl(ux, uy, uz, wx, wy, wz)
    call curl(bx, by, bz, jx, jy, jz)

    !$acc kernels
    f = ux 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, ux_r_z0, ux_r_x0, ux_r_y0)

    !$acc kernels
    f = uy 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, uy_r_z0, uy_r_x0, uy_r_y0)

    !$acc kernels
    f = uz 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, uz_r_z0, uz_r_x0, uz_r_y0)

    !$acc kernels
    f = bx 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, bx_r_z0, bx_r_x0, bx_r_y0)

    !$acc kernels
    f = by 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, by_r_z0, by_r_x0, by_r_y0)

    !$acc kernels
    f = bz 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, bz_r_z0, bz_r_x0, bz_r_y0)

    !$acc kernels
    f = wx 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, wx_r_z0, wx_r_x0, wx_r_y0)

    !$acc kernels
    f = wy 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, wy_r_z0, wy_r_x0, wy_r_y0)

    !$acc kernels
    f = wz 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, wz_r_z0, wz_r_x0, wz_r_y0)

    !$acc kernels
    f = jx 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, jx_r_z0, jx_r_x0, jx_r_y0)

    !$acc kernels
    f = jy 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, jy_r_z0, jy_r_x0, jy_r_y0)

    !$acc kernels
    f = jz 
    !$acc end kernels
    call btran_c2r(f, fr) 
    if(shear) call to_non_shearing_coordinate(fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, jz_r_z0, jz_r_x0, jz_r_y0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv      2D spectra (avg over 1D)     vvvvvvvvvvvvvvvvvv!
    !$acc update host(ux, uy, uz, bx, by, bz)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          u2(k, j, i) = 0.5d0*(abs(ux(k, j, i))**2 + abs(uy(k, j, i))**2 + abs(uz(k, j, i))**2)
          b2(k, j, i) = 0.5d0*(abs(bx(k, j, i))**2 + abs(by(k, j, i))**2 + abs(bz(k, j, i))**2)
          ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
          !-----------------------------------------------------------------------------------------------------------------------------------
          ! The volume integral of a quadratic function is
          ! \int \mathrm{d}^3\mathbf{r}\, f(x,y,z)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2}
          ! \sum_{k_z = -n_{k_z}/2}^{n_{k_z}/2}|f_{k_x, k_y, k_z}|^2 = \left( \sum_{k_z = -n_{k_z}/2}^{-1}\sum_{k_x, k_y} + \sum_{k_z = 0}
          ! \sum_{k_x, k_y} + \sum_{k_z = 1}^{n_{k_z}/2}\sum_{k_x, k_y} \right) |f_{k_x, k_y, k_z}|^2          !
          ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
          !-----------------------------------------------------------------------------------------------------------------------------------
          if (k /= 1) then
            u2(k, j, i) = 2.0d0*u2(k, j, i)
            b2(k, j, i) = 2.0d0*b2(k, j, i)
          endif
        end do
      end do
    end do
    call cut_2d_k(u2, u2_kxy, u2_kyz, u2_kxz)
    call cut_2d_k(b2, b2_kxy, b2_kyz, b2_kxz)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv       bin over (kx, ky) vs kz     vvvvvvvvvvvvvvvvvv!
    call write_polar_spectrum_2d_in_3d(u2, 'z', unit_u2kxy_vs_kz)
    call write_polar_spectrum_2d_in_3d(b2, 'z', unit_b2kxy_vs_kz)
    !vvvvvvvvvvvvvvvvvv       bin over (kx, kz) vs ky     vvvvvvvvvvvvvvvvvv!
    call write_polar_spectrum_2d_in_3d(u2, 'y', unit_u2kxz_vs_ky)
    call write_polar_spectrum_2d_in_3d(b2, 'y', unit_b2kxz_vs_ky)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io_2D( &
                  ux_r_z0, ux_r_x0, ux_r_y0, &
                  uy_r_z0, uy_r_x0, uy_r_y0, &
                  uz_r_z0, uz_r_x0, uz_r_y0, &
                  !
                  wx_r_z0, wx_r_x0, wx_r_y0, &
                  wy_r_z0, wy_r_x0, wy_r_y0, &
                  wz_r_z0, wz_r_x0, wz_r_y0, &
                  !
                  bx_r_z0, bx_r_x0, bx_r_y0, &
                  by_r_z0, by_r_x0, by_r_y0, &
                  bz_r_z0, bz_r_x0, bz_r_y0, &
                  !
                  jx_r_z0, jx_r_x0, jx_r_y0, &
                  jy_r_z0, jy_r_x0, jy_r_y0, &
                  jz_r_z0, jz_r_x0, jz_r_y0, &
                  !
                  u2_kxy, u2_kyz, u2_kxz, &
                  b2_kxy, b2_kyz, b2_kxz  &
                )

    !$acc exit data delete(f )
    !$acc exit data delete(fr)
    deallocate(f)
    deallocate(fr)

    deallocate(ux_r_z0)
    deallocate(ux_r_x0)
    deallocate(ux_r_y0)

    deallocate(uy_r_z0)
    deallocate(uy_r_x0)
    deallocate(uy_r_y0)

    deallocate(uz_r_z0)
    deallocate(uz_r_x0)
    deallocate(uz_r_y0)

    deallocate(wx_r_z0)
    deallocate(wx_r_x0)
    deallocate(wx_r_y0)

    deallocate(wy_r_z0)
    deallocate(wy_r_x0)
    deallocate(wy_r_y0)

    deallocate(wz_r_z0)
    deallocate(wz_r_x0)
    deallocate(wz_r_y0)

    deallocate(bx_r_z0)
    deallocate(bx_r_x0)
    deallocate(bx_r_y0)

    deallocate(by_r_z0)
    deallocate(by_r_x0)
    deallocate(by_r_y0)

    deallocate(bz_r_z0)
    deallocate(bz_r_x0)
    deallocate(bz_r_y0)

    deallocate(jx_r_z0)
    deallocate(jx_r_x0)
    deallocate(jx_r_y0)

    deallocate(jy_r_z0)
    deallocate(jy_r_x0)
    deallocate(jy_r_y0)

    deallocate(jz_r_z0)
    deallocate(jz_r_x0)
    deallocate(jz_r_y0)

    !$acc exit data delete(wx)
    !$acc exit data delete(wy)
    !$acc exit data delete(wz)
    !$acc exit data delete(jx)
    !$acc exit data delete(jy)
    !$acc exit data delete(jz)
    deallocate(wx)
    deallocate(wy)
    deallocate(wz)
    deallocate(jx)
    deallocate(jy)
    deallocate(jz)

    deallocate(u2)
    deallocate(b2)

    deallocate(u2_kxy)
    deallocate(u2_kyz)
    deallocate(u2_kxz)
    deallocate(b2_kxy)
    deallocate(b2_kyz)
    deallocate(b2_kxz)
  end subroutine loop_diagnostics_2D


!-----------------------------------------------!
!> @author  YK
!! @date    4 Jul 2021
!! @brief   Initialize order structure function
!-----------------------------------------------!
  subroutine init_SF2
    use grid, only: lx, ly, lz, dlx, dly, dlz
    use mp, only: nproc
    implicit none
    real(8) :: l, d
    integer :: il

    l = 0.5d0*min(lx, ly, lz)
    d = max(dlx, dly, dlz)

    nl = int(l/d)

    allocate(ll(nl))

    do il = 1, nl
      ! ll(il) = 10.d0**((dlog10(l) - dlog10(d))/(nl - 1)*(il - 1) + dlog10(d))
      ll(il) = (l - d)/(nl - 1)*(il - 1) + d
    enddo

    allocate(lpar(nl), source=ll)
    allocate(lper(nl), source=ll)

    nsamp_r0  = int(sqrt(real(SF2_nsample)/nproc))
    nsamp_ang = int(sqrt(real(SF2_nsample)/nproc))

  end subroutine init_SF2


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Second order structure function
!-----------------------------------------------!
  subroutine loop_diagnostics_SF2
    use io, only: loop_io_SF2
    use mp, only: proc0, proc_id, nproc
    use mp, only: sum_allreduce
    use grid, only: xx, yy, zz
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use fields, only: bx, by, bz, ux, uy, uz
    use cuFFTmp, only: btran_c2r
    use params, only: pi
    use utils, only: ranf
    use time, only: microsleep
    use time_stamp, only: put_time_stamp, timer_diagnostics_total, timer_diagnostics_SF2
    implicit none
    include 'mpif.h'
    real(8) :: ll_x, ll_y, ll_z, theta, phi
    real(8) :: x0, y0, z0, x1, y1, z1
    integer :: i0, j0, k0, i1, j1, k1
    integer :: il, isamp_r0, isamp_ang , iproc
    integer :: demand(nproc, 3)
    real(8) :: supply_b(nproc, 3), supply_u(nproc, 3)
    real(8) :: b0x, b0y, b0z, b1x, b1y, b1z
    real(8) :: u0x, u0y, u0z, u1x, u1y, u1z
    real(8) :: bloc_x, bloc_y, bloc_z, lpar_, lper_
    real(8) :: sf2b(nl, nl), sf2u(nl, nl), count(nl, nl)
    integer :: ilpar, ilper
    complex(8), allocatable, dimension(:,:,:) :: f
    real(8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r, ux_r, uy_r, uz_r
    real(8), allocatable, dimension(:,:,:) :: src

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_SF2)

    allocate(f  (nkz, nky_local, nkx)); f    = 0.d0
    !$acc enter data create(f)

    allocate(src(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(bx_r, source=src)
    allocate(by_r, source=src)
    allocate(bz_r, source=src)
    allocate(ux_r, source=src)
    allocate(uy_r, source=src)
    allocate(uz_r, source=src)

    !$acc enter data create(bx_r)
    !$acc enter data create(by_r)
    !$acc enter data create(bz_r)
    !$acc enter data create(ux_r)
    !$acc enter data create(uy_r)
    !$acc enter data create(uz_r)
    deallocate(src)

    !$acc kernels
    f = bx
    !$acc end kernels
    call btran_c2r(f, bx_r)

    !$acc kernels
    f = by
    !$acc end kernels
    call btran_c2r(f, by_r)
    
    !$acc kernels
    f = bz
    !$acc end kernels
    call btran_c2r(f, bz_r)

    !$acc kernels
    f = ux
    !$acc end kernels
    call btran_c2r(f, ux_r)

    !$acc kernels
    f = uy
    !$acc end kernels
    call btran_c2r(f, uy_r)

    !$acc kernels
    f = uz
    !$acc end kernels
    call btran_c2r(f, uz_r)

    !$acc update host(bx_r, by_r, bz_r, ux_r, uy_r, uz_r)

    sf2b (:, :) = 0.d0
    sf2u (:, :) = 0.d0
    count(:, :) = 0
    do il = 1, nl
      if(proc0) print *, il, '/', nl
      do isamp_r0 = 1, nsamp_r0
        call microsleep(1000)

        ! Pick a random initial point
        i0 = 1 + floor((nlx_local)*ranf()); x0 = xx(i0)
        j0 = 1 + floor((nly      )*ranf()); y0 = yy(j0)
        k0 = 1 + floor((nlz      )*ranf()); z0 = zz(k0)

        b0x = bx_r(k0, j0, i0)
        b0y = by_r(k0, j0, i0)
        b0z = bz_r(k0, j0, i0)
        u0x = ux_r(k0, j0, i0)
        u0y = uy_r(k0, j0, i0)
        u0z = uz_r(k0, j0, i0)

        do isamp_ang = 1, nsamp_ang
          theta = pi*(ranf() - 0.5d0)
          phi   = 2.d0*pi*ranf() 

          ll_x = ll(il)*dsin(theta)*dcos(phi)
          ll_y = ll(il)*dsin(theta)*dsin(phi)
          ll_z = ll(il)*dcos(theta)

          x1 = x0 + ll_x
          y1 = y0 + ll_y
          z1 = z0 + ll_z

          ! When the separation vector crosses the boundary, flip the direction of the separation
          ! so that the separation stays inside the domain.
          ! if(x1 > maxval(xx) .or. x1 < minval(xx)) phi   = pi - phi
          ! if(y1 > maxval(yy) .or. y1 < minval(yy)) phi   =    - phi
          ! if(z1 > maxval(zz) .or. z1 < minval(zz)) theta = pi - theta
          ! x1 = x0 + ll(il)*dsin(theta)*dcos(phi)
          ! y1 = y0 + ll(il)*dsin(theta)*dsin(phi)
          ! z1 = z0 + ll(il)*dcos(theta)
          ! Ignore if the separation goes beyond the boundary
          ! if(x1 > maxval(xx) .or. x1 < minval(xx) .or. &
             ! y1 > maxval(yy) .or. y1 < minval(yy) .or. &
             ! z1 > maxval(zz) .or. z1 < minval(zz)      &
            ! ) then

            ! x1 = x0
            ! y1 = y0
            ! z1 = z0
          ! endif

          i1 = minloc(abs(xx - x1), 1); x1 = xx(i1)
          j1 = minloc(abs(yy - y1), 1); y1 = yy(j1)
          k1 = minloc(abs(zz - z1), 1); z1 = zz(k1)

          ll_x = (x1 - x0)
          ll_y = (y1 - y0)
          ll_z = (z1 - z0)

          ! Periodic boundary
          if(x1 > maxval(xx)) x1 = minval(xx) + x1 - maxval(xx)
          if(y1 > maxval(yy)) y1 = minval(yy) + y1 - maxval(yy) 
          if(z1 > maxval(zz)) z1 = minval(zz) + z1 - maxval(zz)
          if(x1 < minval(xx)) x1 = maxval(xx) + x1 - minval(xx)
          if(y1 < minval(yy)) y1 = maxval(yy) + y1 - minval(yy) 
          if(z1 < minval(zz)) z1 = maxval(zz) + z1 - minval(zz)
          i1 = minloc(abs(xx - x1), 1)
          j1 = minloc(abs(yy - y1), 1)
          k1 = minloc(abs(zz - z1), 1)

          ! Set demand array; "proc_id" process needs a value at (i1, j1, k1) = demand(proc_id+1, :) 
          demand(:, :) = 0
          demand(proc_id+1, :) = [i1, j1, k1]

          call sum_allreduce(demand)

          ! Get B and u vector at the demanded location
          supply_b(:, :) = 0
          supply_u(:, :) = 0
          do iproc = 1, nproc
            i1 = demand(iproc, 1)
            j1 = demand(iproc, 2)
            k1 = demand(iproc, 3)

            if(i1 >= 1 .and. i1 <= nkx .and. &
               j1 >= 1 .and. j1 <= nky_local .and. &
               k1 >= 1 .and. k1 <= nkz) then

               supply_b(iproc, 1) = bx_r(k1, j1, i1)
               supply_b(iproc, 2) = by_r(k1, j1, i1)
               supply_b(iproc, 3) = bz_r(k1, j1, i1)
               supply_u(iproc, 1) = ux_r(k1, j1, i1)
               supply_u(iproc, 2) = uy_r(k1, j1, i1)
               supply_u(iproc, 3) = uz_r(k1, j1, i1)
            endif
          enddo

          call sum_allreduce(supply_b)
          call sum_allreduce(supply_u)

          ! Get B vector at the demanding location
          b1x = supply_b(proc_id+1, 1)
          b1y = supply_b(proc_id+1, 2)
          b1z = supply_b(proc_id+1, 3)
          u1x = supply_u(proc_id+1, 1)
          u1y = supply_u(proc_id+1, 2)
          u1z = supply_u(proc_id+1, 3)

          ! Local mean magnetic field
          bloc_x = (b0x + b1x)/2.d0
          bloc_y = (b0y + b1y)/2.d0
          bloc_z = (b0z + b1z)/2.d0

          ! Parallel and perpendicular components of the separation vector
          if(ll_x**2 + ll_y**2 + ll_z**2 > 0.d0) then
            lpar_ = dabs((ll_x*bloc_x + ll_y*bloc_y + ll_z*bloc_z))/dsqrt((bloc_x**2 + bloc_y**2 + bloc_z**2))
            lper_ = dsqrt((   (bloc_y*ll_z - bloc_z*ll_y)**2 &
                            + (bloc_z*ll_x - bloc_x*ll_z)**2 &
                            + (bloc_x*ll_y - bloc_y*ll_x)**2)/(bloc_x**2 + bloc_y**2 + bloc_z**2))
            ilpar = minloc(abs(lpar_ - lpar), 1)
            ilper = minloc(abs(lper_ - lper), 1)

            sf2b(ilpar, ilper) = sf2b(ilpar, ilper) &
                                + (b1x - b0x)**2 + (b1y - b0y)**2 + (b1z - b0z)**2
            sf2u(ilpar, ilper) = sf2u(ilpar, ilper) &
                                + (u1x - u0x)**2 + (u1y - u0y)**2 + (u1z - u0z)**2
            count(ilpar, ilper) = count(ilpar, ilper) + 1
          endif
        enddo
      enddo
    enddo

    call sum_allreduce(sf2b)
    call sum_allreduce(sf2u)
    call sum_allreduce(count)

    ! Average the structure function
    do ilper = 1, nl
      do ilpar = 1, nl
        if(count(ilpar, ilper) == 0) then
          sf2b(ilpar, ilper) = 0.d0
          sf2u(ilpar, ilper) = 0.d0
        else
          sf2b(ilpar, ilper) = sf2b(ilpar, ilper)/count(ilpar, ilper)
          sf2u(ilpar, ilper) = sf2u(ilpar, ilper)/count(ilpar, ilper)
        endif
      enddo
    enddo

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_SF2)

    call loop_io_SF2(nl, sf2b, sf2u)

    !$acc exit data delete(f)
    !$acc exit data delete(bx_r)
    !$acc exit data delete(by_r)
    !$acc exit data delete(bz_r)
    !$acc exit data delete(ux_r)
    !$acc exit data delete(uy_r)
    !$acc exit data delete(uz_r)
    deallocate(f)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
  end subroutine loop_diagnostics_SF2


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Calculate kpar(k) & delta b/b0
!           k_\|(k) = \left(\frac{\langle|\mathbf{b}_{0,k} \cdot\nabla \delta\mathbf{b}_k|^2\rangle}
!           {\langle b_{0,k}^2\rangle\langle \delta b_k^2\rangle}\right)^{1/2}  \\ 
!           \mathbf{b}_{0,k}(\mathbf{x}) = \calF^{-1}\sum_{|\bm{k}|' \le k/2} \mathbf{b}_{\mathbf{k}'} \\  
!           \delta\mathbf{b}_{k}(\mathbf{x}) = \calF^{-1}\sum_{k/2 \le |\bm{k}|' \le 2k} \mathbf{b}_{\mathbf{k}'}
!-----------------------------------------------!
  subroutine loop_diagnostics_kpar
    use fields, only: bx, by, bz, ux, uy, uz
    use mp, only: proc0, sum_allreduce
    use grid, only: nlx, nlx_local, nly, nlz
    use grid, only: ky, kz, nkx, nky_local, nkz
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx_local, nly, nlz_padded
    use cuFFTmp, only: btran_c2r
    use params, only: zi
    use shearing_box, only: k2t, kxt
    use io, only: loop_io_kpar
    use time_stamp, only: put_time_stamp, timer_diagnostics_total, timer_diagnostics_kpar
    implicit none
    integer :: ii, i, j, k

    real   (8), dimension(:), allocatable :: kpar_b, kpar_u, b1_ovr_b0

    complex(8), allocatable, dimension(:,:,:) :: bx0, by0, bz0                ! local mean field
    complex(8), allocatable, dimension(:,:,:) :: bx1, by1, bz1                ! local fluctuating field
    complex(8), allocatable, dimension(:,:,:) :: ux1, uy1, uz1                ! local fluctuating field
    complex(8), allocatable, dimension(:,:,:) :: dbx1_dx, dby1_dx, dbz1_dx    
    complex(8), allocatable, dimension(:,:,:) :: dbx1_dy, dby1_dy, dbz1_dy    
    complex(8), allocatable, dimension(:,:,:) :: dbx1_dz, dby1_dz, dbz1_dz    
    complex(8), allocatable, dimension(:,:,:) :: dux1_dx, duy1_dx, duz1_dx    
    complex(8), allocatable, dimension(:,:,:) :: dux1_dy, duy1_dy, duz1_dy    
    complex(8), allocatable, dimension(:,:,:) :: dux1_dz, duy1_dz, duz1_dz    
    real   (8), allocatable, dimension(:,:,:) :: bx0r, by0r, bz0r             
    real   (8), allocatable, dimension(:,:,:) :: bx1r, by1r, bz1r             
    real   (8), allocatable, dimension(:,:,:) :: ux1r, uy1r, uz1r             
    real   (8), allocatable, dimension(:,:,:) :: dbx1_dxr, dby1_dxr, dbz1_dxr 
    real   (8), allocatable, dimension(:,:,:) :: dbx1_dyr, dby1_dyr, dbz1_dyr 
    real   (8), allocatable, dimension(:,:,:) :: dbx1_dzr, dby1_dzr, dbz1_dzr 
    real   (8), allocatable, dimension(:,:,:) :: dux1_dxr, duy1_dxr, duz1_dxr 
    real   (8), allocatable, dimension(:,:,:) :: dux1_dyr, duy1_dyr, duz1_dyr 
    real   (8), allocatable, dimension(:,:,:) :: dux1_dzr, duy1_dzr, duz1_dzr 

    real   (8), allocatable, dimension(:,:,:) :: b0_gradb1_sq, b0_gradu1_sq, b0sq, b1sq, u1sq
    real   (8) :: b0_gradb1_sq_avg, b0_gradu1_sq_avg, b0sq_avg, b1sq_avg, u1sq_avg

    real   (8), allocatable, dimension(:,:,:) :: bx0hat, by0hat, bz0hat ! local mean field unit vector
    real   (8), allocatable, dimension(:,:,:) :: b1par, b1prpx, b1prpy  ! b1par : projection of b1 to b0 => Pseudo AW
                                                                        ! b1per : b1 - b1par*b0hat       => Shear AW
    real   (8), allocatable, dimension(:,:,:) :: u1par, u1prpx, u1prpy  ! u1par : projection of u1 to b0 => Pseudo AW
                                                                        ! u1per : u1 - u1par*b0hat       => Shear AW
    real   (8), allocatable, dimension(:) :: b1par2, b1prp2
    real   (8), allocatable, dimension(:) :: u1par2, u1prp2

    complex(8), allocatable, dimension(:,:,:) :: src_c
    real   (8), allocatable, dimension(:,:,:) :: src_r

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_kpar)

    allocate(kpar_b   (nkpolar_log), source=0.d0)
    allocate(kpar_u   (nkpolar_log), source=0.d0)
    allocate(b1_ovr_b0(nkpolar_log), source=0.d0)
    allocate(b1par2   (nkpolar_log), source=0.d0)
    allocate(b1prp2   (nkpolar_log), source=0.d0)
    allocate(u1par2   (nkpolar_log), source=0.d0)
    allocate(u1prp2   (nkpolar_log), source=0.d0)

    allocate(src_c(nkz, nky_local, nkx), source=(0.d0,0.d0))
    allocate(bx0    , source=src_c)
    allocate(by0    , source=src_c)
    allocate(bz0    , source=src_c)
    allocate(bx1    , source=src_c)
    allocate(by1    , source=src_c)
    allocate(bz1    , source=src_c)
    allocate(ux1    , source=src_c)
    allocate(uy1    , source=src_c)
    allocate(uz1    , source=src_c)
    allocate(dbx1_dx, source=src_c)
    allocate(dby1_dx, source=src_c)
    allocate(dbz1_dx, source=src_c)
    allocate(dbx1_dy, source=src_c)
    allocate(dby1_dy, source=src_c)
    allocate(dbz1_dy, source=src_c)
    allocate(dbx1_dz, source=src_c)
    allocate(dby1_dz, source=src_c)
    allocate(dbz1_dz, source=src_c)
    allocate(dux1_dx, source=src_c)
    allocate(duy1_dx, source=src_c)
    allocate(duz1_dx, source=src_c)
    allocate(dux1_dy, source=src_c)
    allocate(duy1_dy, source=src_c)
    allocate(duz1_dy, source=src_c)
    allocate(dux1_dz, source=src_c)
    allocate(duy1_dz, source=src_c)
    allocate(duz1_dz, source=src_c)

    !$acc enter data create(bx0    )
    !$acc enter data create(by0    )
    !$acc enter data create(bz0    )
    !$acc enter data create(bx1    )
    !$acc enter data create(by1    )
    !$acc enter data create(bz1    )
    !$acc enter data create(ux1    )
    !$acc enter data create(uy1    )
    !$acc enter data create(uz1    )
    !$acc enter data create(dbx1_dx)
    !$acc enter data create(dby1_dx)
    !$acc enter data create(dbz1_dx)
    !$acc enter data create(dbx1_dy)
    !$acc enter data create(dby1_dy)
    !$acc enter data create(dbz1_dy)
    !$acc enter data create(dbx1_dz)
    !$acc enter data create(dby1_dz)
    !$acc enter data create(dbz1_dz)
    !$acc enter data create(dux1_dx)
    !$acc enter data create(duy1_dx)
    !$acc enter data create(duz1_dx)
    !$acc enter data create(dux1_dy)
    !$acc enter data create(duy1_dy)
    !$acc enter data create(duz1_dy)
    !$acc enter data create(dux1_dz)
    !$acc enter data create(duy1_dz)
    !$acc enter data create(duz1_dz)
    deallocate(src_c)

    allocate(src_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(bx0r    , source=src_r)
    allocate(by0r    , source=src_r)
    allocate(bz0r    , source=src_r)
    allocate(bx1r    , source=src_r)
    allocate(by1r    , source=src_r)
    allocate(bz1r    , source=src_r)
    allocate(ux1r    , source=src_r)
    allocate(uy1r    , source=src_r)
    allocate(uz1r    , source=src_r)
    allocate(dbx1_dxr, source=src_r)
    allocate(dby1_dxr, source=src_r)
    allocate(dbz1_dxr, source=src_r)
    allocate(dbx1_dyr, source=src_r)
    allocate(dby1_dyr, source=src_r)
    allocate(dbz1_dyr, source=src_r)
    allocate(dbx1_dzr, source=src_r)
    allocate(dby1_dzr, source=src_r)
    allocate(dbz1_dzr, source=src_r)
    allocate(dux1_dxr, source=src_r)
    allocate(duy1_dxr, source=src_r)
    allocate(duz1_dxr, source=src_r)
    allocate(dux1_dyr, source=src_r)
    allocate(duy1_dyr, source=src_r)
    allocate(duz1_dyr, source=src_r)
    allocate(dux1_dzr, source=src_r)
    allocate(duy1_dzr, source=src_r)
    allocate(duz1_dzr, source=src_r)

    !$acc enter data create(bx0r    )
    !$acc enter data create(by0r    )
    !$acc enter data create(bz0r    )
    !$acc enter data create(bx1r    )
    !$acc enter data create(by1r    )
    !$acc enter data create(bz1r    )
    !$acc enter data create(ux1r    )
    !$acc enter data create(uy1r    )
    !$acc enter data create(uz1r    )
    !$acc enter data create(dbx1_dxr)
    !$acc enter data create(dby1_dxr)
    !$acc enter data create(dbz1_dxr)
    !$acc enter data create(dbx1_dyr)
    !$acc enter data create(dby1_dyr)
    !$acc enter data create(dbz1_dyr)
    !$acc enter data create(dbx1_dzr)
    !$acc enter data create(dby1_dzr)
    !$acc enter data create(dbz1_dzr)
    !$acc enter data create(dux1_dxr)
    !$acc enter data create(duy1_dxr)
    !$acc enter data create(duz1_dxr)
    !$acc enter data create(dux1_dyr)
    !$acc enter data create(duy1_dyr)
    !$acc enter data create(duz1_dyr)
    !$acc enter data create(dux1_dzr)
    !$acc enter data create(duy1_dzr)
    !$acc enter data create(duz1_dzr)

    allocate(b0_gradb1_sq, source=src_r)
    allocate(b0_gradu1_sq, source=src_r)
    allocate(b0sq        , source=src_r)
    allocate(b1sq        , source=src_r)
    allocate(u1sq        , source=src_r)

    !$acc enter data create(b0_gradb1_sq)
    !$acc enter data create(b0_gradu1_sq)
    !$acc enter data create(b0sq        )
    !$acc enter data create(b1sq        )
    !$acc enter data create(u1sq        )

    allocate(bx0hat      , source=src_r)
    allocate(by0hat      , source=src_r)
    allocate(bz0hat      , source=src_r)
    allocate(b1par       , source=src_r)
    allocate(b1prpx      , source=src_r)
    allocate(b1prpy      , source=src_r)
    allocate(u1par       , source=src_r)
    allocate(u1prpx      , source=src_r)
    allocate(u1prpy      , source=src_r)
    deallocate(src_r)

    ! get kpar for each kprp_log(ii)
    do ii = 1, nkpolar_log

      ! filter out to get local mean field and local fluctuating field
      !$acc kernels
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            ! smaller than kprp/2
            if(k2t(k, j, i) < (0.5d0*kpbin_log(ii))**2) then
              bx0(k, j, i) = bx(k, j, i)
              by0(k, j, i) = by(k, j, i)
              bz0(k, j, i) = bz(k, j, i)
            else
              bx0(k, j, i) = 0.d0
              by0(k, j, i) = 0.d0
              bz0(k, j, i) = 0.d0
            endif

            ! larger than kprp/2 and smaller than 2*kprp
            if(k2t(k, j, i) >= (0.5d0*kpbin_log(ii))**2 .and. k2t(k, j, i) < (2.0d0*kpbin_log(ii))**2) then
              bx1(k, j, i) = bx(k, j, i)
              by1(k, j, i) = by(k, j, i)
              bz1(k, j, i) = bz(k, j, i)
              ux1(k, j, i) = ux(k, j, i)
              uy1(k, j, i) = uy(k, j, i)
              uz1(k, j, i) = uz(k, j, i)
            else
              bx1(k, j, i) = 0.d0
              by1(k, j, i) = 0.d0
              bz1(k, j, i) = 0.d0
              ux1(k, j, i) = 0.d0
              uy1(k, j, i) = 0.d0
              uz1(k, j, i) = 0.d0
            endif

            if(abs(bx0(k,j,i)) < epsilon(1.d0) .or. abs(bx0(k,j,i)) > 1.d0/epsilon(1.d0)) bx0(k,j,i) = 0.d0
            if(abs(by0(k,j,i)) < epsilon(1.d0) .or. abs(by0(k,j,i)) > 1.d0/epsilon(1.d0)) by0(k,j,i) = 0.d0
            if(abs(bz0(k,j,i)) < epsilon(1.d0) .or. abs(bz0(k,j,i)) > 1.d0/epsilon(1.d0)) bz0(k,j,i) = 0.d0
            if(abs(bx1(k,j,i)) < epsilon(1.d0) .or. abs(bx1(k,j,i)) > 1.d0/epsilon(1.d0)) bx1(k,j,i) = 0.d0
            if(abs(by1(k,j,i)) < epsilon(1.d0) .or. abs(by1(k,j,i)) > 1.d0/epsilon(1.d0)) by1(k,j,i) = 0.d0
            if(abs(bz1(k,j,i)) < epsilon(1.d0) .or. abs(bz1(k,j,i)) > 1.d0/epsilon(1.d0)) bz1(k,j,i) = 0.d0
            if(abs(ux1(k,j,i)) < epsilon(1.d0) .or. abs(ux1(k,j,i)) > 1.d0/epsilon(1.d0)) ux1(k,j,i) = 0.d0
            if(abs(uy1(k,j,i)) < epsilon(1.d0) .or. abs(uy1(k,j,i)) > 1.d0/epsilon(1.d0)) uy1(k,j,i) = 0.d0
            if(abs(uz1(k,j,i)) < epsilon(1.d0) .or. abs(uz1(k,j,i)) > 1.d0/epsilon(1.d0)) uz1(k,j,i) = 0.d0

            dbx1_dx(k, j, i) = zi*kxt(i,j)*bx1(k, j, i)
            dby1_dx(k, j, i) = zi*kxt(i,j)*by1(k, j, i)
            dbz1_dx(k, j, i) = zi*kxt(i,j)*bz1(k, j, i)

            dbx1_dy(k, j, i) = zi*ky(j)   *bx1(k, j, i)
            dby1_dy(k, j, i) = zi*ky(j)   *by1(k, j, i)
            dbz1_dy(k, j, i) = zi*ky(j)   *bz1(k, j, i)
                                          
            dbx1_dz(k, j, i) = zi*kz(k)   *bx1(k, j, i)
            dby1_dz(k, j, i) = zi*kz(k)   *by1(k, j, i)
            dbz1_dz(k, j, i) = zi*kz(k)   *bz1(k, j, i)
                                          
            dux1_dx(k, j, i) = zi*kxt(i,j)*ux1(k, j, i)
            duy1_dx(k, j, i) = zi*kxt(i,j)*uy1(k, j, i)
            duz1_dx(k, j, i) = zi*kxt(i,j)*uz1(k, j, i)
                                          
            dux1_dy(k, j, i) = zi*ky(j)   *ux1(k, j, i)
            duy1_dy(k, j, i) = zi*ky(j)   *uy1(k, j, i)
            duz1_dy(k, j, i) = zi*ky(j)   *uz1(k, j, i)
                                          
            dux1_dz(k, j, i) = zi*kz(k)   *ux1(k, j, i)
            duy1_dz(k, j, i) = zi*kz(k)   *uy1(k, j, i)
            duz1_dz(k, j, i) = zi*kz(k)   *uz1(k, j, i)
          end do
        end do
      end do
      !$acc end kernels

      call btran_c2r(bx0, bx0r)
      call btran_c2r(by0, by0r)
      call btran_c2r(bz0, bz0r)
      call btran_c2r(bx1, bx1r)
      call btran_c2r(by1, by1r)
      call btran_c2r(bz1, bz1r)
      call btran_c2r(ux1, ux1r)
      call btran_c2r(uy1, uy1r)
      call btran_c2r(uz1, uz1r)

      call btran_c2r(dbx1_dx, dbx1_dxr)
      call btran_c2r(dby1_dx, dby1_dxr)
      call btran_c2r(dbz1_dx, dbz1_dxr)
      call btran_c2r(dbx1_dy, dbx1_dyr)
      call btran_c2r(dby1_dy, dby1_dyr)
      call btran_c2r(dbz1_dy, dbz1_dyr)
      call btran_c2r(dbx1_dz, dbx1_dzr)
      call btran_c2r(dby1_dz, dby1_dzr)
      call btran_c2r(dbz1_dz, dbz1_dzr)

      call btran_c2r(dux1_dx, dux1_dxr)
      call btran_c2r(duy1_dx, duy1_dxr)
      call btran_c2r(duz1_dx, duz1_dxr)
      call btran_c2r(dux1_dy, dux1_dyr)
      call btran_c2r(duy1_dy, duy1_dyr)
      call btran_c2r(duz1_dy, duz1_dyr)
      call btran_c2r(dux1_dz, dux1_dzr)
      call btran_c2r(duy1_dz, duy1_dzr)
      call btran_c2r(duz1_dz, duz1_dzr)

      !$acc update host(bx0r)
      !$acc update host(by0r)
      !$acc update host(bz0r)
      !$acc update host(bx1r)
      !$acc update host(by1r)
      !$acc update host(bz1r)
      !$acc update host(ux1r)
      !$acc update host(uy1r)
      !$acc update host(uz1r)

      !$acc update host(dbx1_dxr)
      !$acc update host(dby1_dxr)
      !$acc update host(dbz1_dxr)
      !$acc update host(dbx1_dyr)
      !$acc update host(dby1_dyr)
      !$acc update host(dbz1_dyr)
      !$acc update host(dbx1_dzr)
      !$acc update host(dby1_dzr)
      !$acc update host(dbz1_dzr)
                                
      !$acc update host(dux1_dxr)
      !$acc update host(duy1_dxr)
      !$acc update host(duz1_dxr)
      !$acc update host(dux1_dyr)
      !$acc update host(duy1_dyr)
      !$acc update host(duz1_dyr)
      !$acc update host(dux1_dzr)
      !$acc update host(duy1_dzr)
      !$acc update host(duz1_dzr)

      ! Get kpar & delta b/b0
      b0_gradb1_sq =   (bx0r*dbx1_dxr + by0r*dbx1_dyr + bz0r*dbx1_dzr)**2 &
                     + (bx0r*dby1_dxr + by0r*dby1_dyr + bz0r*dby1_dzr)**2 &
                     + (bx0r*dbz1_dxr + by0r*dbz1_dyr + bz0r*dbz1_dzr)**2 
      b0_gradu1_sq =   (bx0r*dux1_dxr + by0r*dux1_dyr + bz0r*dux1_dzr)**2 &
                     + (bx0r*duy1_dxr + by0r*duy1_dyr + bz0r*duy1_dzr)**2 &
                     + (bx0r*duz1_dxr + by0r*duz1_dyr + bz0r*duz1_dzr)**2 
      b0sq = bx0r**2 + by0r**2 + bz0r**2
      b1sq = bx1r**2 + by1r**2 + bz1r**2
      u1sq = ux1r**2 + uy1r**2 + uz1r**2

      b0_gradb1_sq_avg = sum(b0_gradb1_sq(:, :, :nlz)); call sum_allreduce(b0_gradb1_sq_avg); b0_gradb1_sq_avg = b0_gradb1_sq_avg/nlx/nly/nlz
      b0_gradu1_sq_avg = sum(b0_gradu1_sq(:, :, :nlz)); call sum_allreduce(b0_gradu1_sq_avg); b0_gradu1_sq_avg = b0_gradu1_sq_avg/nlx/nly/nlz
      b0sq_avg         = sum(b0sq        (:, :, :nlz)); call sum_allreduce(b0sq_avg        ); b0sq_avg         = b0sq_avg        /nlx/nly/nlz
      b1sq_avg         = sum(b1sq        (:, :, :nlz)); call sum_allreduce(b1sq_avg        ); b1sq_avg         = b1sq_avg        /nlx/nly/nlz
      u1sq_avg         = sum(u1sq        (:, :, :nlz)); call sum_allreduce(u1sq_avg        ); u1sq_avg         = u1sq_avg        /nlx/nly/nlz

      if (b0sq_avg /= 0.d0 .and. b1sq_avg /= 0.d0 .and. u1sq_avg /= 0.d0) then
        kpar_b   (ii) = dsqrt( b0_gradb1_sq_avg/(b1sq_avg*b0sq_avg) )
        kpar_u   (ii) = dsqrt( b0_gradu1_sq_avg/(u1sq_avg*b0sq_avg) )
        b1_ovr_b0(ii) = dsqrt( b1sq_avg/b0sq_avg )
      else
        kpar_b   (ii) = 0.d0
        kpar_u   (ii) = 0.d0
        b1_ovr_b0(ii) = 0.d0
      endif

      ! Get Shear AWs and pseudo AWs
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            if(bx0r(k,j,i)**2 + by0r(k,j,i)**2 + bz0r(k,j,i)**2 /= 0.d0) then 
              bx0hat(k,j,i) = bx0r(k,j,i)/dsqrt(bx0r(k,j,i)**2 + by0r(k,j,i)**2 + bz0r(k,j,i)**2)
              by0hat(k,j,i) = by0r(k,j,i)/dsqrt(bx0r(k,j,i)**2 + by0r(k,j,i)**2 + bz0r(k,j,i)**2)
              bz0hat(k,j,i) = bz0r(k,j,i)/dsqrt(bx0r(k,j,i)**2 + by0r(k,j,i)**2 + bz0r(k,j,i)**2)
            else 
              bx0hat(k,j,i) = 0.d0 
              by0hat(k,j,i) = 0.d0 
              bz0hat(k,j,i) = 0.d0 
            endif

            b1par (k,j,i) = bx1r(k,j,i)*bx0hat(k,j,i) + by1r(k,j,i)*by0hat(k,j,i) + bz1r(k,j,i)*bz0hat(k,j,i)
            b1prpx(k,j,i) = bx1r(k,j,i) - b1par(k,j,i)*bx0hat(k,j,i)
            b1prpy(k,j,i) = by1r(k,j,i) - b1par(k,j,i)*by0hat(k,j,i)

            u1par (k,j,i) = ux1r(k,j,i)*bx0hat(k,j,i) + uy1r(k,j,i)*by0hat(k,j,i) + uz1r(k,j,i)*bz0hat(k,j,i)
            u1prpx(k,j,i) = ux1r(k,j,i) - u1par(k,j,i)*bx0hat(k,j,i)
            u1prpy(k,j,i) = uy1r(k,j,i) - u1par(k,j,i)*by0hat(k,j,i)
          enddo
        enddo
      enddo

      b1par2(ii) = sum(b1par**2)             ; call sum_allreduce(b1par2(ii)); b1par2(ii) = b1par2(ii)/nlx/nly/nlz
      b1prp2(ii) = sum(b1prpx**2 + b1prpy**2); call sum_allreduce(b1prp2(ii)); b1prp2(ii) = b1prp2(ii)/nlx/nly/nlz
      u1par2(ii) = sum(u1par**2)             ; call sum_allreduce(u1par2(ii)); u1par2(ii) = u1par2(ii)/nlx/nly/nlz
      u1prp2(ii) = sum(u1prpx**2 + u1prpy**2); call sum_allreduce(u1prp2(ii)); u1prp2(ii) = u1prp2(ii)/nlx/nly/nlz
    enddo

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_kpar)

    call loop_io_kpar(nkpolar_log, kpar_b, kpar_u, b1_ovr_b0, &
                      b1par2, b1prp2, u1par2, u1prp2)

    !$acc exit data delete(bx0    )
    !$acc exit data delete(by0    )
    !$acc exit data delete(bz0    )
    !$acc exit data delete(bx1    )
    !$acc exit data delete(by1    )
    !$acc exit data delete(bz1    )
    !$acc exit data delete(ux1    )
    !$acc exit data delete(uy1    )
    !$acc exit data delete(uz1    )
    !$acc exit data delete(dbx1_dx)
    !$acc exit data delete(dby1_dx)
    !$acc exit data delete(dbz1_dx)
    !$acc exit data delete(dbx1_dy)
    !$acc exit data delete(dby1_dy)
    !$acc exit data delete(dbz1_dy)
    !$acc exit data delete(dbx1_dz)
    !$acc exit data delete(dby1_dz)
    !$acc exit data delete(dbz1_dz)
    !$acc exit data delete(dux1_dx)
    !$acc exit data delete(duy1_dx)
    !$acc exit data delete(duz1_dx)
    !$acc exit data delete(dux1_dy)
    !$acc exit data delete(duy1_dy)
    !$acc exit data delete(duz1_dy)
    !$acc exit data delete(dux1_dz)
    !$acc exit data delete(duy1_dz)
    !$acc exit data delete(duz1_dz)

    !$acc exit data delete(bx0r    )
    !$acc exit data delete(by0r    )
    !$acc exit data delete(bz0r    )
    !$acc exit data delete(bx1r    )
    !$acc exit data delete(by1r    )
    !$acc exit data delete(bz1r    )
    !$acc exit data delete(ux1r    )
    !$acc exit data delete(uy1r    )
    !$acc exit data delete(uz1r    )
    !$acc exit data delete(dbx1_dxr)
    !$acc exit data delete(dby1_dxr)
    !$acc exit data delete(dbz1_dxr)
    !$acc exit data delete(dbx1_dyr)
    !$acc exit data delete(dby1_dyr)
    !$acc exit data delete(dbz1_dyr)
    !$acc exit data delete(dbx1_dzr)
    !$acc exit data delete(dby1_dzr)
    !$acc exit data delete(dbz1_dzr)
    !$acc exit data delete(dux1_dxr)
    !$acc exit data delete(duy1_dxr)
    !$acc exit data delete(duz1_dxr)
    !$acc exit data delete(dux1_dyr)
    !$acc exit data delete(duy1_dyr)
    !$acc exit data delete(duz1_dyr)
    !$acc exit data delete(dux1_dzr)
    !$acc exit data delete(duy1_dzr)
    !$acc exit data delete(duz1_dzr)

    !$acc exit data delete(b0_gradb1_sq)
    !$acc exit data delete(b0_gradu1_sq)
    !$acc exit data delete(b0sq        )
    !$acc exit data delete(b1sq        )
    !$acc exit data delete(u1sq        )

    deallocate(kpar_b)
    deallocate(kpar_u)
    deallocate(b1_ovr_b0)
    deallocate(bx0, by0, bz0)
    deallocate(bx1, by1, bz1)
    deallocate(ux1, uy1, uz1)
    deallocate(dbx1_dx, dby1_dx, dbz1_dx)
    deallocate(dbx1_dy, dby1_dy, dbz1_dy)
    deallocate(dbx1_dz, dby1_dz, dbz1_dz)
    deallocate(dux1_dx, duy1_dx, duz1_dx)
    deallocate(dux1_dy, duy1_dy, duz1_dy)
    deallocate(dux1_dz, duy1_dz, duz1_dz)
    deallocate(bx0r, by0r, bz0r)
    deallocate(bx1r, by1r, bz1r)
    deallocate(ux1r, uy1r, uz1r)
    deallocate(dbx1_dxr, dby1_dxr, dbz1_dxr)
    deallocate(dbx1_dyr, dby1_dyr, dbz1_dyr)
    deallocate(dbx1_dzr, dby1_dzr, dbz1_dzr)
    deallocate(dux1_dxr, duy1_dxr, duz1_dxr)
    deallocate(dux1_dyr, duy1_dyr, duz1_dyr)
    deallocate(dux1_dzr, duy1_dzr, duz1_dzr)
    deallocate(b0_gradb1_sq, b0_gradu1_sq, b0sq, b1sq, u1sq)
    deallocate(bx0hat)
    deallocate(by0hat)
    deallocate(bz0hat)
    deallocate(b1par )
    deallocate(b1prpx)
    deallocate(b1prpy)
    deallocate(u1par )
    deallocate(u1prpx)
    deallocate(u1prpy)
    deallocate(b1par2)
    deallocate(b1prp2)
    deallocate(u1par2)
    deallocate(u1prp2)

  end subroutine loop_diagnostics_kpar


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Calculate shell-to-shell transfer
!!          trans_uu(K, Q) : u(Q) -> u(K)
!!          trans_bb(K, Q) : b(Q) -> b(K)
!!          trans_ub(K, Q) : u(Q) -> b(K)
!!          trans_bu(K, Q) : b(Q) -> u(K)
!-----------------------------------------------!
  subroutine loop_diagnostics_nltrans
    use fields, only: nfields
    use fields, only: iux, iuy, iuz
    use fields, only: ibx, iby, ibz
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use params, only: zi
    use grid, only: ky, kz
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use cuFFTmp, only: ftran_r2c, btran_c2r
    use mp, only: proc0, sum_reduce
    use shearing_box, only: k2t, kxt
    use io, only: loop_io_nltrans
    use time_stamp, only: put_time_stamp, timer_diagnostics_total, timer_diagnostics_nltrans
    implicit none


    ! Forward FFT variables
    integer :: nftran = 36
    integer :: iflx_uu_xx = 1 , iflx_uu_xy = 2 , iflx_uu_xz = 3  !  
    integer :: iflx_uu_yx = 4 , iflx_uu_yy = 5 , iflx_uu_yz = 6  ! uu^Q (u^Q is filtered on |k| = Q) 
    integer :: iflx_uu_zx = 7 , iflx_uu_zy = 8 , iflx_uu_zz = 9  !  

    integer :: iflx_bb_xx = 10, iflx_bb_xy = 11, iflx_bb_xz = 12 !  
    integer :: iflx_bb_yx = 13, iflx_bb_yy = 14, iflx_bb_yz = 15 ! bb^Q (b^Q is filtered on |k| = Q) 
    integer :: iflx_bb_zx = 16, iflx_bb_zy = 17, iflx_bb_zz = 18 !  

    integer :: iflx_ub_xx = 19, iflx_ub_xy = 20, iflx_ub_xz = 21 !  
    integer :: iflx_ub_yx = 22, iflx_ub_yy = 23, iflx_ub_yz = 24 ! ub^Q (b^Q is filtered on |k| = Q) 
    integer :: iflx_ub_zx = 25, iflx_ub_zy = 26, iflx_ub_zz = 27 !  

    integer :: iflx_bu_xx = 28, iflx_bu_xy = 29, iflx_bu_xz = 30 !  
    integer :: iflx_bu_yx = 31, iflx_bu_yy = 32, iflx_bu_yz = 33 ! bu^Q (u^Q is filtered on |k| = Q) 
    integer :: iflx_bu_zx = 34, iflx_bu_zy = 35, iflx_bu_zz = 36 !  

    integer :: nnonlin = 12
    integer :: inonlin_uu_x = 1 , inonlin_uu_y = 2 , inonlin_uu_z = 3 
    integer :: inonlin_bb_x = 4 , inonlin_bb_y = 5 , inonlin_bb_z = 6 
    integer :: inonlin_ub_x = 7 , inonlin_ub_y = 8 , inonlin_ub_z = 9 
    integer :: inonlin_bu_x = 10, inonlin_bu_y = 11, inonlin_bu_z = 12


    integer :: i, j, k, ii, jj, idx
    real(8) :: filter

    real   (8), dimension(:,:), allocatable :: trans_uu, trans_bb, trans_ub, trans_bu
    complex(8), allocatable, dimension(:,:,:,:) :: w, wtmp 
    complex(8), allocatable, dimension(:,:,:,:) :: w_filtd
    complex(8), allocatable, dimension(:,:,:,:) :: flx
    complex(8), allocatable, dimension(:,:,:,:) :: nonlin

    real   (8), allocatable, dimension(:,:,:,:) :: w_r
    real   (8), allocatable, dimension(:,:,:,:) :: w_filtd_r
    real   (8), allocatable, dimension(:,:,:,:) :: flx_r

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_nltrans)

    allocate(trans_uu (nkpolar_log, nkpolar_log), source=0.d0)
    allocate(trans_bb (nkpolar_log, nkpolar_log), source=0.d0)
    allocate(trans_ub (nkpolar_log, nkpolar_log), source=0.d0)
    allocate(trans_bu (nkpolar_log, nkpolar_log), source=0.d0)

    if(.not. allocated(w   )) then
      allocate(w   (nkz, nky_local, nkx, nfields), source=(0.d0,0.d0))
      !$acc enter data create(w)
    endif
    if(.not. allocated(wtmp)) then
      allocate(wtmp(nkz, nky_local, nkx, nfields), source=(0.d0,0.d0))
      !$acc enter data create(wtmp)
    endif
    if(.not. allocated(w_r )) then
      allocate(w_r (nlz_padded, nly, nlx_local, nfields), source=0.d0)
      !$acc enter data create(w_r)
    endif

    !$acc parallel loop collapse(3)
    do i =1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          w(k,j,i,iux) = ux(k,j,i)
          w(k,j,i,iuy) = uy(k,j,i)
          w(k,j,i,iuz) = uz(k,j,i)
          w(k,j,i,ibx) = bx(k,j,i)
          w(k,j,i,iby) = by(k,j,i)
          w(k,j,i,ibz) = bz(k,j,i)
        enddo
      enddo
    enddo

    ! get unfiltered fields in real space
    !$acc kernels
    wtmp = w
    !$acc end kernels
    
    do i = 1, nfields
      call btran_c2r(w(:,:,:,i), w_r(:,:,:,i))
    enddo

    !$acc kernels
    w = wtmp
    !$acc end kernels
    
    if(allocated(wtmp)) then
      !$acc exit data delete(wtmp)
      deallocate(wtmp)
    endif


    ! get nonlinear transfer for each kprp_log(ii)
    do ii = 1, nkpolar_log
      if(proc0) print *, ii, '/', nkpolar_log

      if(.not. allocated(w_filtd)) then
        allocate(w_filtd(nkz, nky_local, nkx, nfields), source=(0.d0,0.d0))
        !$acc enter data create(w_filtd)
      endif

      ! filter out where kprp != kprp_log(ii)
      !$acc kernels
      do i =1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            if(ii == nkpolar_log) then
              if(k2t(k, j, i) >= (kpbin_log(ii))**2) then
                filter = 1.d0
              else
                filter = 0.d0
              endif
            else
              if(k2t(k, j, i) >= (kpbin_log(ii))**2 .and. k2t(k, j, i) < (kpbin_log(ii + 1))**2) then
                filter = 1.d0
              else
                filter = 0.d0
              endif
            endif

            do idx = 1, nfields
              w_filtd(k, j, i, idx) = filter*w(k, j, i, idx)
            enddo

          end do
        end do
      end do
      !$acc end kernels

      if(.not. allocated(w_filtd_r)) then
        allocate(w_filtd_r(nlz_padded, nly, nlx_local, nfields), source=0.d0)
        !$acc enter data create(w_filtd_r)
      endif

      ! get filtered fields in real space
      do i = 1, nfields
        call btran_c2r(w_filtd(:,:,:,i), w_filtd_r(:,:,:,i))
      enddo

      if(allocated(w_filtd)) then
        !$acc exit data delete(w_filtd)
        deallocate(w_filtd)
      endif

      if(.not. allocated(flx_r)) then
        allocate(flx_r(nlz_padded, nly, nlx_local, nftran ), source=0.d0)
        !$acc enter data create(flx_r)
      endif

      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz_padded
            ! uu^Q (u^Q is filtered on |k| = Q) 
            flx_r(k,j,i,iflx_uu_xx) = w_r(k,j,i,iux)*w_filtd_r(k,j,i,iux)
            flx_r(k,j,i,iflx_uu_xy) = w_r(k,j,i,iux)*w_filtd_r(k,j,i,iuy)
            flx_r(k,j,i,iflx_uu_xz) = w_r(k,j,i,iux)*w_filtd_r(k,j,i,iuz)

            flx_r(k,j,i,iflx_uu_yx) = w_r(k,j,i,iuy)*w_filtd_r(k,j,i,iux)
            flx_r(k,j,i,iflx_uu_yy) = w_r(k,j,i,iuy)*w_filtd_r(k,j,i,iuy)
            flx_r(k,j,i,iflx_uu_yz) = w_r(k,j,i,iuy)*w_filtd_r(k,j,i,iuz)

            flx_r(k,j,i,iflx_uu_zx) = w_r(k,j,i,iuz)*w_filtd_r(k,j,i,iux)
            flx_r(k,j,i,iflx_uu_zy) = w_r(k,j,i,iuz)*w_filtd_r(k,j,i,iuy)
            flx_r(k,j,i,iflx_uu_zz) = w_r(k,j,i,iuz)*w_filtd_r(k,j,i,iuz)

            ! bb^Q (b^Q is filtered on |k| = Q) 
            flx_r(k,j,i,iflx_bb_xx) = w_r(k,j,i,ibx)*w_filtd_r(k,j,i,ibx)
            flx_r(k,j,i,iflx_bb_xy) = w_r(k,j,i,ibx)*w_filtd_r(k,j,i,iby)
            flx_r(k,j,i,iflx_bb_xz) = w_r(k,j,i,ibx)*w_filtd_r(k,j,i,ibz)

            flx_r(k,j,i,iflx_bb_yx) = w_r(k,j,i,iby)*w_filtd_r(k,j,i,ibx)
            flx_r(k,j,i,iflx_bb_yy) = w_r(k,j,i,iby)*w_filtd_r(k,j,i,iby)
            flx_r(k,j,i,iflx_bb_yz) = w_r(k,j,i,iby)*w_filtd_r(k,j,i,ibz)

            flx_r(k,j,i,iflx_bb_zx) = w_r(k,j,i,ibz)*w_filtd_r(k,j,i,ibx)
            flx_r(k,j,i,iflx_bb_zy) = w_r(k,j,i,ibz)*w_filtd_r(k,j,i,iby)
            flx_r(k,j,i,iflx_bb_zz) = w_r(k,j,i,ibz)*w_filtd_r(k,j,i,ibz)

            ! ub^Q (b^Q is filtered on |k| = Q) 
            flx_r(k,j,i,iflx_ub_xx) = w_r(k,j,i,iux)*w_filtd_r(k,j,i,ibx)
            flx_r(k,j,i,iflx_ub_xy) = w_r(k,j,i,iux)*w_filtd_r(k,j,i,iby)
            flx_r(k,j,i,iflx_ub_xz) = w_r(k,j,i,iux)*w_filtd_r(k,j,i,ibz)

            flx_r(k,j,i,iflx_ub_yx) = w_r(k,j,i,iuy)*w_filtd_r(k,j,i,ibx)
            flx_r(k,j,i,iflx_ub_yy) = w_r(k,j,i,iuy)*w_filtd_r(k,j,i,iby)
            flx_r(k,j,i,iflx_ub_yz) = w_r(k,j,i,iuy)*w_filtd_r(k,j,i,ibz)

            flx_r(k,j,i,iflx_ub_zx) = w_r(k,j,i,iuz)*w_filtd_r(k,j,i,ibx)
            flx_r(k,j,i,iflx_ub_zy) = w_r(k,j,i,iuz)*w_filtd_r(k,j,i,iby)
            flx_r(k,j,i,iflx_ub_zz) = w_r(k,j,i,iuz)*w_filtd_r(k,j,i,ibz)

            ! bu^Q (u^Q is filtered on |k| = Q) 
            flx_r(k,j,i,iflx_bu_xx) = w_r(k,j,i,ibx)*w_filtd_r(k,j,i,iux)
            flx_r(k,j,i,iflx_bu_xy) = w_r(k,j,i,ibx)*w_filtd_r(k,j,i,iuy)
            flx_r(k,j,i,iflx_bu_xz) = w_r(k,j,i,ibx)*w_filtd_r(k,j,i,iuz)

            flx_r(k,j,i,iflx_bu_yx) = w_r(k,j,i,iby)*w_filtd_r(k,j,i,iux)
            flx_r(k,j,i,iflx_bu_yy) = w_r(k,j,i,iby)*w_filtd_r(k,j,i,iuy)
            flx_r(k,j,i,iflx_bu_yz) = w_r(k,j,i,iby)*w_filtd_r(k,j,i,iuz)

            flx_r(k,j,i,iflx_bu_zx) = w_r(k,j,i,ibz)*w_filtd_r(k,j,i,iux)
            flx_r(k,j,i,iflx_bu_zy) = w_r(k,j,i,ibz)*w_filtd_r(k,j,i,iuy)
            flx_r(k,j,i,iflx_bu_zz) = w_r(k,j,i,ibz)*w_filtd_r(k,j,i,iuz)

          enddo
        enddo
      enddo

      if(allocated(w_filtd_r)) then
        !$acc exit data delete(w_filtd_r)
        deallocate(w_filtd_r)
      endif

      if(.not. allocated(flx)) then
        allocate(flx(nkz, nky_local, nkx, nftran ), source=(0.d0,0.d0))
        !$acc enter data create(flx)
      endif

      do i = 1, nftran
        call ftran_r2c(flx_r(:,:,:,i), flx(:,:,:,i))
      enddo

      if(allocated(flx_r)) then
        !$acc exit data delete(flx_r)
        deallocate(flx_r)
      endif

      !$acc kernels
      flx = flx/ntot
      !$acc end kernels

      if(.not. allocated(nonlin)) then
        allocate(nonlin(nkz, nky_local, nkx, nnonlin), source=(0.d0,0.d0))
        !$acc enter data create(nonlin)
      endif

      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            nonlin(k,j,i,inonlin_uu_x) = zi*(  kxt(i,j)*flx(k,j,i,iflx_uu_xx) &
                                             + ky (j)  *flx(k,j,i,iflx_uu_yx) &
                                             + kz (k)  *flx(k,j,i,iflx_uu_zx) ) 
            nonlin(k,j,i,inonlin_uu_y) = zi*(  kxt(i,j)*flx(k,j,i,iflx_uu_xy) &
                                             + ky (j)  *flx(k,j,i,iflx_uu_yy) &
                                             + kz (k)  *flx(k,j,i,iflx_uu_zy) ) 
            nonlin(k,j,i,inonlin_uu_z) = zi*(  kxt(i,j)*flx(k,j,i,iflx_uu_xz) &
                                             + ky (j)  *flx(k,j,i,iflx_uu_yz) &
                                             + kz (k)  *flx(k,j,i,iflx_uu_zz) ) 

            nonlin(k,j,i,inonlin_bb_x) = zi*(  kxt(i,j)*flx(k,j,i,iflx_bb_xx) &
                                             + ky (j)  *flx(k,j,i,iflx_bb_yx) &
                                             + kz (k)  *flx(k,j,i,iflx_bb_zx) ) 
            nonlin(k,j,i,inonlin_bb_y) = zi*(  kxt(i,j)*flx(k,j,i,iflx_bb_xy) &
                                             + ky (j)  *flx(k,j,i,iflx_bb_yy) &
                                             + kz (k)  *flx(k,j,i,iflx_bb_zy) ) 
            nonlin(k,j,i,inonlin_bb_z) = zi*(  kxt(i,j)*flx(k,j,i,iflx_bb_xz) &
                                             + ky (j)  *flx(k,j,i,iflx_bb_yz) &
                                             + kz (k)  *flx(k,j,i,iflx_bb_zz) ) 

            nonlin(k,j,i,inonlin_ub_x) = zi*(  kxt(i,j)*flx(k,j,i,iflx_ub_xx) &
                                             + ky (j)  *flx(k,j,i,iflx_ub_yx) &
                                             + kz (k)  *flx(k,j,i,iflx_ub_zx) ) 
            nonlin(k,j,i,inonlin_ub_y) = zi*(  kxt(i,j)*flx(k,j,i,iflx_ub_xy) &
                                             + ky (j)  *flx(k,j,i,iflx_ub_yy) &
                                             + kz (k)  *flx(k,j,i,iflx_ub_zy) ) 
            nonlin(k,j,i,inonlin_ub_z) = zi*(  kxt(i,j)*flx(k,j,i,iflx_ub_xz) &
                                             + ky (j)  *flx(k,j,i,iflx_ub_yz) &
                                             + kz (k)  *flx(k,j,i,iflx_ub_zz) ) 

            nonlin(k,j,i,inonlin_bu_x) = zi*(  kxt(i,j)*flx(k,j,i,iflx_bu_xx) &
                                             + ky (j)  *flx(k,j,i,iflx_bu_yx) &
                                             + kz (k)  *flx(k,j,i,iflx_bu_zx) ) 
            nonlin(k,j,i,inonlin_bu_y) = zi*(  kxt(i,j)*flx(k,j,i,iflx_bu_xy) &
                                             + ky (j)  *flx(k,j,i,iflx_bu_yy) &
                                             + kz (k)  *flx(k,j,i,iflx_bu_zy) ) 
            nonlin(k,j,i,inonlin_bu_z) = zi*(  kxt(i,j)*flx(k,j,i,iflx_bu_xz) &
                                             + ky (j)  *flx(k,j,i,iflx_bu_yz) &
                                             + kz (k)  *flx(k,j,i,iflx_bu_zz) ) 

          enddo
        enddo
      enddo

      if(allocated(flx)) then
        !$acc exit data delete(flx)
        deallocate(flx)
      endif

      ! get nonlinear transfer for each kprp_log(jj)
      do jj = 1, nkpolar_log
        !$acc parallel loop collapse(3)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz

              if(jj == nkpolar_log) then
                if(k2t(k, j, i) >= (kpbin_log(jj))**2) then
                  filter = 1.d0
                else
                  filter = 0.d0
                endif
              else
                if(k2t(k, j, i) >= (kpbin_log(jj))**2 .and. k2t(k, j, i) < (kpbin_log(jj + 1))**2) then
                  filter = 1.d0
                else
                  filter = 0.d0
                endif
              endif

              ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
              !-----------------------------------------------------------------------------------------------------------------------------------
              ! The volume integral of a quadratic function is
              ! \int \mathrm{d}^3\mathbf{r}\, f(x,y,z)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2}
              ! \sum_{k_z = -n_{k_z}/2}^{n_{k_z}/2}|f_{k_x, k_y, k_z}|^2 = \left( \sum_{k_z = -n_{k_z}/2}^{-1}\sum_{k_x, k_y} + \sum_{k_z = 0}
              ! \sum_{k_x, k_y} + \sum_{k_z = 1}^{n_{k_z}/2}\sum_{k_x, k_y} \right) |f_{k_x, k_y, k_z}|^2          !
              ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
              !-----------------------------------------------------------------------------------------------------------------------------------
              if (k /= 1) filter = filter * 2

              trans_uu(jj, ii) = trans_uu(jj, ii) - filter*dble( &
                                    w(k,j,i,iux)*conjg(nonlin(k,j,i,inonlin_uu_x)) &
                                  + w(k,j,i,iuy)*conjg(nonlin(k,j,i,inonlin_uu_y)) &
                                  + w(k,j,i,iuz)*conjg(nonlin(k,j,i,inonlin_uu_z)) &
                                 )                                                                                                 
                                                                                                                                   
              trans_bb(jj, ii) = trans_bb(jj, ii) - filter*dble( &                                                                 
                                    w(k,j,i,ibx)*conjg(nonlin(k,j,i,inonlin_ub_x)) &
                                  + w(k,j,i,iby)*conjg(nonlin(k,j,i,inonlin_ub_y)) &
                                  + w(k,j,i,ibz)*conjg(nonlin(k,j,i,inonlin_ub_z)) &
                                 )                                                                                                 
                                                                                                                                   
              trans_ub(jj, ii) = trans_ub(jj, ii) + filter*dble( &                                                                 
                                    w(k,j,i,ibx)*conjg(nonlin(k,j,i,inonlin_bu_x)) &
                                  + w(k,j,i,iby)*conjg(nonlin(k,j,i,inonlin_bu_y)) &
                                  + w(k,j,i,ibz)*conjg(nonlin(k,j,i,inonlin_bu_z)) &
                                 )                                                                                                 
                                                                                                                                   
              trans_bu(jj, ii) = trans_bu(jj, ii) + filter*dble( &                                                                 
                                    w(k,j,i,iux)*conjg(nonlin(k,j,i,inonlin_bb_x)) &
                                  + w(k,j,i,iuy)*conjg(nonlin(k,j,i,inonlin_bb_y)) &
                                  + w(k,j,i,iuz)*conjg(nonlin(k,j,i,inonlin_bb_z)) &
                                 )

            enddo
          enddo
        enddo
      enddo

      if(allocated(nonlin)) then
        !$acc exit data delete(nonlin)
        deallocate(nonlin)
      endif


    enddo

    call sum_reduce(trans_uu, 0)
    call sum_reduce(trans_bb, 0)
    call sum_reduce(trans_ub, 0)
    call sum_reduce(trans_bu, 0)

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_nltrans)

    call loop_io_nltrans(nkpolar_log, trans_uu, trans_bb, trans_ub, trans_bu)

    deallocate(trans_uu )
    deallocate(trans_bb )
    deallocate(trans_ub )
    deallocate(trans_bu )
    if(allocated(w)) then
      !$acc exit data delete(w)
      deallocate(w  )
    endif
    if(allocated(w_r)) then
      !$acc exit data delete(w_r)
      deallocate(w_r)
    endif

  end subroutine loop_diagnostics_nltrans

end module diagnostics




