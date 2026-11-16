!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Diagnostics for MHD_COMP_ISOTH
!!
!!          Gallope port: spectral data are (kz, ky, kx)
!!          with kz the r2c half axis and ky split over
!!          comm_fft, so the Fourier sums over the kz >= 0
!!          half are doubled for kz /= 0 (index k /= 1);
!!          real data are (z, y, x) with x split over
!!          comm_fft. The model runs with P_m = P_s = 1
!!          (model_specific), so comm_fft spans all ranks
!!          and proc0 is iproc_fft = 0, which holds the
!!          k = 0 mode, the ky = 0 plane and the x = 0
!!          plane.
!-----------------------------------------------!
module diagnostics
  use diagnostics_common
  implicit none

  public :: init_diagnostics, finish_diagnostics
  public :: loop_diagnostics, loop_diagnostics_2D, loop_diagnostics_kpar, loop_diagnostics_SF2
  public :: loop_diagnostics_nltrans

  private

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
    call init_io(nkpolar, kpbin)
    call init_series_modes
  end subroutine init_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @date    15 Apr 2020
!! @brief   Diagnostics in loop
!-----------------------------------------------!
  subroutine loop_diagnostics
    use io, only: loop_io
    use mp, only: proc0, iproc_fft, comm_fft
    use grid, only: kx, ky, kz, k2_max
    use grid, only: nlx, nly, nlz
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx_local, nlz_padded
    use fields, only: lnrho, lnrho_to_rho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use fields, only: lnrho_old
    use fields, only: mx_old, my_old, mz_old
    use fields, only: bx_old, by_old, bz_old
    use fields, only: m_to_u
    use mp, only: sum_reduce
    use cuFFTmp, only: btran_c2r
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, nu, nu_h, nu_h_exp, eta, eta_h, eta_h_exp, lmd, lmd_h, lmd_h_exp, cs2va2, shear, q
    use force, only: fmx, fmy, fmz, fmx_old, fmy_old, fmz_old
    use shearing_box, only: shear_flg, tsc, nremap, k2t
    implicit none
    integer :: i, j, k

    ! Quantities that are only summed over k are per-mode scalars of the
    ! kernel below (Calliope stores them as arrays and sums them afterwards);
    ! the quantities binned over |k| stay arrays.
    real(8) :: wkin, wmag
    real(8) :: wkin_old, wmag_old
    real(8) :: wkin_dissip, wmag_dissip, wrho_dissip
    real(8) :: wrho
    real(8) :: wrho_old
    real(8)   , allocatable, dimension(:,:,:) :: u2, ux2, uy2, uz2
    real(8)   , allocatable, dimension(:,:,:) :: b2, bx2, by2, bz2
    real(8) :: p_ext, p_re, p_ma
    real(8)   , allocatable, dimension(:,:,:) :: zp2, zm2
    real(8)   , allocatable, dimension(:,:,:) :: rho_abs ! abs(rho), binned unweighted as in Calliope
    real(8)   , allocatable, dimension(:,:,:) :: src_r
    complex(8), allocatable, dimension(:,:,:) :: ux, uy, uz
    complex(8), allocatable, dimension(:,:,:) :: ux_old, uy_old, uz_old
    complex(8), allocatable, dimension(:,:,:) :: rho, rho_old
    complex(8), allocatable, dimension(:,:,:) :: u2half, u2half_old
    complex(8), allocatable, dimension(:,:,:) :: src_c
    real(8)   , allocatable, dimension(:,:,:) :: rho_r, ux_r, uy_r, uz_r, bx_r, by_r, bz_r

    real(8) :: wkin_sum, wmag_sum
    real(8) :: wkin_dot_sum, wmag_dot_sum
    real(8) :: wkin_dissip_sum, wmag_dissip_sum, wrho_dissip_sum
    real(8) :: wrho_sum
    real(8) :: wrho_dot_sum
    real(8) :: p_ext_sum, p_re_sum, p_ma_sum
    real(8) :: zp2_sum, zm2_sum
    real(8) :: rho_rms, u_rms, b_rms
    real(8) :: smach_rms, amach_rms, beta_rms
    real(8) :: bx0, by0, bz0 ! mean magnetic field

    real(8), dimension(:), allocatable :: u2_bin, ux2_bin, uy2_bin, uz2_bin
    real(8), dimension(:), allocatable :: b2_bin, bx2_bin, by2_bin, bz2_bin
    real(8), dimension(:), allocatable :: rho_bin
    real(8), dimension(:), allocatable :: zp2_bin, zm2_bin

    complex(8) ::  mx_mid,    my_mid,  mz_mid
    complex(8) ::  ux_mid,    uy_mid,  uz_mid
    complex(8) ::  bx_mid,    by_mid,  bz_mid
    complex(8) :: rho_mid, lnrho_mid,  u2half_mid
    complex(8) :: fmx_mid, fmy_mid, fmz_mid

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(src_c(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(ux        , source=src_c)
    allocate(uy        , source=src_c)
    allocate(uz        , source=src_c)
    allocate(ux_old    , source=src_c)
    allocate(uy_old    , source=src_c)
    allocate(uz_old    , source=src_c)
    allocate(rho       , source=src_c)
    allocate(rho_old   , source=src_c)
    allocate(u2half    , source=src_c)
    allocate(u2half_old, source=src_c)
    deallocate(src_c)
    !$acc enter data copyin(ux, uy, uz, ux_old, uy_old, uz_old)
    !$acc enter data copyin(rho, rho_old, u2half, u2half_old)

    allocate(src_r(nkz, nky_local, nkx), source=0.d0)
    allocate(u2     , source=src_r)
    allocate(ux2    , source=src_r)
    allocate(uy2    , source=src_r)
    allocate(uz2    , source=src_r)
    allocate(b2     , source=src_r)
    allocate(bx2    , source=src_r)
    allocate(by2    , source=src_r)
    allocate(bz2    , source=src_r)
    allocate(zp2    , source=src_r)
    allocate(zm2    , source=src_r)
    allocate(rho_abs, source=src_r)
    deallocate(src_r)
    !$acc enter data create(u2, ux2, uy2, uz2, b2, bx2, by2, bz2)
    !$acc enter data create(zp2, zm2, rho_abs)

    allocate(src_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(rho_r, source=src_r)
    allocate( ux_r, source=src_r)
    allocate( uy_r, source=src_r)
    allocate( uz_r, source=src_r)
    allocate( bx_r, source=src_r)
    allocate( by_r, source=src_r)
    allocate( bz_r, source=src_r)
    deallocate(src_r)
    !$acc enter data create(rho_r, ux_r, uy_r, uz_r, bx_r, by_r, bz_r)

    allocate( u2_bin(1:nkpolar), source=0.d0)
    allocate(ux2_bin(1:nkpolar), source=0.d0)
    allocate(uy2_bin(1:nkpolar), source=0.d0)
    allocate(uz2_bin(1:nkpolar), source=0.d0)
    allocate( b2_bin(1:nkpolar), source=0.d0)
    allocate(bx2_bin(1:nkpolar), source=0.d0)
    allocate(by2_bin(1:nkpolar), source=0.d0)
    allocate(bz2_bin(1:nkpolar), source=0.d0)
    allocate(rho_bin(1:nkpolar), source=0.d0)
    allocate(zp2_bin(1:nkpolar), source=0.d0)
    allocate(zm2_bin(1:nkpolar), source=0.d0)
    !$acc enter data create(u2_bin, ux2_bin, uy2_bin, uz2_bin)
    !$acc enter data create(b2_bin, bx2_bin, by2_bin, bz2_bin)
    !$acc enter data create(rho_bin, zp2_bin, zm2_bin)

    call lnrho_to_rho(lnrho    , rho    )
    call lnrho_to_rho(lnrho_old, rho_old)

    call m_to_u(lnrho    , mx    , my    , mz    , ux    , uy    , uz    )
    call m_to_u(lnrho_old, mx_old, my_old, mz_old, ux_old, uy_old, uz_old)

    call get_u2half(ux    , uy    , uz    , u2half    )
    call get_u2half(ux_old, uy_old, uz_old, u2half_old)

    wkin_sum        = 0.d0
    wmag_sum        = 0.d0
    wrho_sum        = 0.d0
    wkin_dot_sum    = 0.d0
    wmag_dot_sum    = 0.d0
    wrho_dot_sum    = 0.d0
    wkin_dissip_sum = 0.d0
    wmag_dissip_sum = 0.d0
    wrho_dissip_sum = 0.d0
    p_ext_sum       = 0.d0
    p_re_sum        = 0.d0
    p_ma_sum        = 0.d0
    zp2_sum         = 0.d0
    zm2_sum         = 0.d0

    !$acc parallel loop collapse(3) gang vector default(present) &
    !$acc& private(wkin, wmag, wkin_old, wmag_old, wrho, wrho_old) &
    !$acc& private(wkin_dissip, wmag_dissip, wrho_dissip, p_ext, p_re, p_ma) &
    !$acc& private(mx_mid, my_mid, mz_mid, ux_mid, uy_mid, uz_mid, bx_mid, by_mid, bz_mid) &
    !$acc& private(fmx_mid, fmy_mid, fmz_mid, rho_mid, lnrho_mid, u2half_mid) &
    !$acc& reduction(+:wkin_sum, wmag_sum, wrho_sum, wkin_dot_sum, wmag_dot_sum, wrho_dot_sum) &
    !$acc& reduction(+:wkin_dissip_sum, wmag_dissip_sum, wrho_dissip_sum) &
    !$acc& reduction(+:p_ext_sum, p_re_sum, p_ma_sum, zp2_sum, zm2_sum)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ! The mean field (kx = ky = kz = 0) is taken after this kernel.

          wkin     = 0.25d0*( &
                               mx    (k, j, i)*conjg(ux    (k, j, i)) + ux    (k, j, i)*conjg(mx    (k, j, i)) &
                             + my    (k, j, i)*conjg(uy    (k, j, i)) + uy    (k, j, i)*conjg(my    (k, j, i)) &
                             + mz    (k, j, i)*conjg(uz    (k, j, i)) + uz    (k, j, i)*conjg(mz    (k, j, i)) &
                            )
          wkin_old = 0.25d0*( &
                               mx_old(k, j, i)*conjg(ux_old(k, j, i)) + ux_old(k, j, i)*conjg(mx_old(k, j, i)) &
                             + my_old(k, j, i)*conjg(uy_old(k, j, i)) + uy_old(k, j, i)*conjg(my_old(k, j, i)) &
                             + mz_old(k, j, i)*conjg(uz_old(k, j, i)) + uz_old(k, j, i)*conjg(mz_old(k, j, i)) &
                            )

          wmag     = 0.5d0*( &
                               bx    (k, j, i)*conjg(bx    (k, j, i)) &
                             + by    (k, j, i)*conjg(by    (k, j, i)) &
                             + bz    (k, j, i)*conjg(bz    (k, j, i)) &
                            )
          wmag_old = 0.5d0*( &
                               bx_old(k, j, i)*conjg(bx_old(k, j, i)) &
                             + by_old(k, j, i)*conjg(by_old(k, j, i)) &
                             + bz_old(k, j, i)*conjg(bz_old(k, j, i)) &
                            )

          u2      (k, j, i) = 0.5d0*(abs(ux(k, j, i))**2 + abs(uy(k, j, i))**2 + abs(uz(k, j, i))**2)
          ux2     (k, j, i) = 0.5d0*(abs(ux(k, j, i))**2)
          uy2     (k, j, i) = 0.5d0*(abs(uy(k, j, i))**2)
          uz2     (k, j, i) = 0.5d0*(abs(uz(k, j, i))**2)

          b2      (k, j, i) = 0.5d0*(abs(bx(k, j, i))**2 + abs(by(k, j, i))**2 + abs(bz(k, j, i))**2)
          bx2     (k, j, i) = 0.5d0*(abs(bx(k, j, i))**2)
          by2     (k, j, i) = 0.5d0*(abs(by(k, j, i))**2)
          bz2     (k, j, i) = 0.5d0*(abs(bz(k, j, i))**2)

          wrho     = 0.5d0*cs2va2*( &
                          rho    (k, j, i)*conjg(lnrho    (k, j, i)) + lnrho    (k, j, i)*conjg(rho    (k, j, i)) &
                        )
          wrho_old = 0.5d0*cs2va2*( &
                          rho_old(k, j, i)*conjg(lnrho_old(k, j, i)) + lnrho_old(k, j, i)*conjg(rho_old(k, j, i)) &
                        )

          zp2     (k, j, i) = abs(ux(k,j,i) + bx(k,j,i))**2 + abs(uy(k,j,i) + by(k,j,i))**2 + abs(uz(k,j,i) + bz(k,j,i))**2
          zm2     (k, j, i) = abs(ux(k,j,i) - bx(k,j,i))**2 + abs(uy(k,j,i) - by(k,j,i))**2 + abs(uz(k,j,i) - bz(k,j,i))**2

           mx_mid    = 0.5d0*( mx   (k, j, i) +  mx_old   (k, j, i))
           my_mid    = 0.5d0*( my   (k, j, i) +  my_old   (k, j, i))
           mz_mid    = 0.5d0*( mz   (k, j, i) +  mz_old   (k, j, i))
           ux_mid    = 0.5d0*( ux   (k, j, i) +  ux_old   (k, j, i))
           uy_mid    = 0.5d0*( uy   (k, j, i) +  uy_old   (k, j, i))
           uz_mid    = 0.5d0*( uz   (k, j, i) +  uz_old   (k, j, i))
           bx_mid    = 0.5d0*( bx   (k, j, i) +  bx_old   (k, j, i))
           by_mid    = 0.5d0*( by   (k, j, i) +  by_old   (k, j, i))
           bz_mid    = 0.5d0*( bz   (k, j, i) +  bz_old   (k, j, i))
          fmx_mid    = 0.5d0*(fmx   (k, j, i) + fmx_old   (k, j, i))
          fmy_mid    = 0.5d0*(fmy   (k, j, i) + fmy_old   (k, j, i))
          fmz_mid    = 0.5d0*(fmz   (k, j, i) + fmz_old   (k, j, i))
          rho_mid    = 0.5d0*(rho   (k, j, i) + rho_old   (k, j, i))
          lnrho_mid  = 0.5d0*(lnrho (k, j, i) + lnrho_old (k, j, i))
          u2half_mid = 0.5d0*(u2half(k, j, i) + u2half_old(k, j, i))

          wkin_dissip = (nu *(k2t(k, j, i)/k2_max) + nu_h *(k2t(k, j, i)/k2_max)**nu_h_exp )*0.5d0*( &
                                                             mx_mid*conjg(ux_mid) + ux_mid*conjg(mx_mid) &
                                                           + my_mid*conjg(uy_mid) + uy_mid*conjg(my_mid) &
                                                           + mz_mid*conjg(uz_mid) + uz_mid*conjg(mz_mid) &
                                                          )
          wmag_dissip = (eta*(k2t(k, j, i)/k2_max) + eta_h*(k2t(k, j, i)/k2_max)**eta_h_exp) &
                          *(abs(bx_mid)**2 + abs(by_mid)**2 + abs(bz_mid)**2)
          wrho_dissip = (lmd*(k2t(k, j, i)/k2_max) + lmd_h *(k2t(k, j, i)/k2_max)**lmd_h_exp)*0.5d0*( &
                                                               (-u2half_mid + cs2va2*lnrho_mid)*conjg(rho_mid) &
                                                             + rho_mid*conjg(-u2half_mid + cs2va2*lnrho_mid) &
                                                           )

          p_ext = 0.5d0*( &
                           (fmx_mid*conjg(ux_mid) + conjg(fmx_mid)*ux_mid) &
                         + (fmy_mid*conjg(uy_mid) + conjg(fmy_mid)*uy_mid) &
                         + (fmz_mid*conjg(uz_mid) + conjg(fmz_mid)*uz_mid) &
                       )
          p_re  = + 0.5d0*q*shear_flg*(ux_mid*conjg(my_mid) + conjg(ux_mid)*my_mid)
          p_ma  = - 0.5d0*q*shear_flg*(bx_mid*conjg(by_mid) + conjg(bx_mid)*by_mid)

          ! The reason for the following treatment for kz == 0 mode is the following. Compile it with LaTeX.
          !-----------------------------------------------------------------------------------------------------------------------------------
          ! The volume integral of a quadratic function is
          ! \int \mathrm{d}^3\mathbf{r}\, f(x,y,z)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2}
          ! \sum_{k_z = -n_{k_z}/2}^{n_{k_z}/2}|f_{k_x, k_y, k_z}|^2 = \left( \sum_{k_z = -n_{k_z}/2}^{-1}\sum_{k_x, k_y} + \sum_{k_z = 0}
          ! \sum_{k_x, k_y} + \sum_{k_z = 1}^{n_{k_z}/2}\sum_{k_x, k_y} \right) |f_{k_x, k_y, k_z}|^2
          ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
          !-----------------------------------------------------------------------------------------------------------------------------------
          if (k /= 1) then
            wkin              = 2.0d0*wkin
            wkin_old          = 2.0d0*wkin_old
            wkin_dissip       = 2.0d0*wkin_dissip
            wmag              = 2.0d0*wmag
            wmag_old          = 2.0d0*wmag_old
            wmag_dissip       = 2.0d0*wmag_dissip
            wrho              = 2.0d0*wrho
            wrho_old          = 2.0d0*wrho_old
            wrho_dissip       = 2.0d0*wrho_dissip
            u2       (k, j, i) = 2.0d0*u2       (k, j, i)
            ux2      (k, j, i) = 2.0d0*ux2      (k, j, i)
            uy2      (k, j, i) = 2.0d0*uy2      (k, j, i)
            uz2      (k, j, i) = 2.0d0*uz2      (k, j, i)
            b2       (k, j, i) = 2.0d0*b2       (k, j, i)
            bx2      (k, j, i) = 2.0d0*bx2      (k, j, i)
            by2      (k, j, i) = 2.0d0*by2      (k, j, i)
            bz2      (k, j, i) = 2.0d0*bz2      (k, j, i)
            p_ext             = 2.0d0*p_ext
            p_re              = 2.0d0*p_re
            p_ma              = 2.0d0*p_ma
            zp2      (k, j, i) = 2.0d0*zp2      (k, j, i)
            zm2      (k, j, i) = 2.0d0*zm2      (k, j, i)
          endif

          ! abs(rho) is binned without the kz weight, as in Calliope
          rho_abs(k, j, i) = abs(rho(k, j, i))

          ! local part of the sums over (kx, ky, kz) below
          wkin_sum        = wkin_sum        + wkin
          wmag_sum        = wmag_sum        + wmag
          wrho_sum        = wrho_sum        + wrho
          wkin_dot_sum    = wkin_dot_sum    + (wkin - wkin_old)/dt
          wmag_dot_sum    = wmag_dot_sum    + (wmag - wmag_old)/dt
          wrho_dot_sum    = wrho_dot_sum    + (wrho - wrho_old)/dt
          wkin_dissip_sum = wkin_dissip_sum + wkin_dissip
          wmag_dissip_sum = wmag_dissip_sum + wmag_dissip
          wrho_dissip_sum = wrho_dissip_sum + wrho_dissip
          p_ext_sum       = p_ext_sum       + p_ext
          p_re_sum        = p_re_sum        + p_re
          p_ma_sum        = p_ma_sum        + p_ma
          zp2_sum         = zp2_sum         + zp2(k, j, i)
          zm2_sum         = zm2_sum         + zm2(k, j, i)
        end do
      end do
    end do
    !$acc end parallel loop

    ! Mean field: the kx = ky = kz = 0 mode is (1, 1, 1) on iproc_fft = 0,
    ! which is proc0 (the only rank loop_io writes from)
    bx0 = 0.d0
    by0 = 0.d0
    bz0 = 0.d0
    if (iproc_fft == 0) then
      !$acc serial default(present) copyout(bx0, by0, bz0)
      bx0 = abs(bx(1, 1, 1))
      by0 = abs(by(1, 1, 1))
      bz0 = abs(bz(1, 1, 1))
      !$acc end serial
    endif

    !vvvvvvvvvvvvvvvvvv     integrate over kx, ky, kz     vvvvvvvvvvvvvvvvvv!
    ! the local sums are accumulated in the kernel above
    call sum_reduce(wkin_sum, 0, comm=comm_fft)
    call sum_reduce(wmag_sum, 0, comm=comm_fft)
    call sum_reduce(wrho_sum, 0, comm=comm_fft)

    call sum_reduce(wkin_dot_sum, 0, comm=comm_fft)
    call sum_reduce(wmag_dot_sum, 0, comm=comm_fft)
    call sum_reduce(wrho_dot_sum, 0, comm=comm_fft)

    call sum_reduce(wkin_dissip_sum, 0, comm=comm_fft)
    call sum_reduce(wmag_dissip_sum, 0, comm=comm_fft)
    call sum_reduce(wrho_dissip_sum, 0, comm=comm_fft)

    call sum_reduce(p_ext_sum, 0, comm=comm_fft)
    call sum_reduce(p_re_sum , 0, comm=comm_fft)
    call sum_reduce(p_ma_sum , 0, comm=comm_fft)

    call sum_reduce(zp2_sum, 0, comm=comm_fft)
    call sum_reduce(zm2_sum, 0, comm=comm_fft)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv       bin over (kx, ky, kz)       vvvvvvvvvvvvvvvvvv!
    ! each call reduces over comm_fft; the host copy of the bin holds the
    ! reduced spectrum on proc0
    call get_polar_spectrum_3d( u2,  u2_bin)
    call get_polar_spectrum_3d(ux2, ux2_bin)
    call get_polar_spectrum_3d(uy2, uy2_bin)
    call get_polar_spectrum_3d(uz2, uz2_bin)
    call get_polar_spectrum_3d( b2,  b2_bin)
    call get_polar_spectrum_3d(bx2, bx2_bin)
    call get_polar_spectrum_3d(by2, by2_bin)
    call get_polar_spectrum_3d(bz2, bz2_bin)
    call get_polar_spectrum_3d(rho_abs, rho_bin)
    call get_polar_spectrum_3d(zp2, zp2_bin)
    call get_polar_spectrum_3d(zm2, zm2_bin)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv         maximum mach values       vvvvvvvvvvvvvvvvvv!
    ! btran_c2r preserves its input, so Calliope's copy into f is not needed
    call btran_c2r(rho, rho_r)
    call btran_c2r(ux , ux_r )
    call btran_c2r(uy , uy_r )
    call btran_c2r(uz , uz_r )
    call btran_c2r(bx , bx_r )
    call btran_c2r(by , by_r )
    call btran_c2r(bz , bz_r )

    ! local sums over the physical points (k = 1, nlz; the padding is excluded)
    rho_rms = 0.d0
    u_rms   = 0.d0
    b_rms   = 0.d0
    !$acc parallel loop collapse(3) gang vector default(present) &
    !$acc& reduction(+:rho_rms, u_rms, b_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_rms = rho_rms + rho_r(k, j, i)**2
          u_rms   = u_rms   + (ux_r(k, j, i)**2 + uy_r(k, j, i)**2 + uz_r(k, j, i)**2)
          b_rms   = b_rms   + (bx_r(k, j, i)**2 + by_r(k, j, i)**2 + bz_r(k, j, i)**2)
        end do
      end do
    end do
    !$acc end parallel loop

    call sum_reduce(rho_rms, 0, comm=comm_fft)
    rho_rms = sqrt(rho_rms/nlx/nly/nlz)

    call sum_reduce(u_rms, 0, comm=comm_fft)
    u_rms = sqrt(u_rms/nlx/nly/nlz)

    call sum_reduce(b_rms, 0, comm=comm_fft)
    b_rms = sqrt(b_rms/nlx/nly/nlz)

    smach_rms = u_rms/sqrt(cs2va2)
    amach_rms = sqrt(rho_rms)*u_rms/b_rms
    beta_rms  = 2.d0*cs2va2*rho_rms/b_rms**2
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  wkin_sum, wmag_sum, wrho_sum, &
                  wkin_dot_sum, wmag_dot_sum, wrho_dot_sum, &
                  wkin_dissip_sum, wmag_dissip_sum, wrho_dissip_sum, &
                  p_ext_sum, p_re_sum, p_ma_sum, &
                  smach_rms, amach_rms, beta_rms, &
                  zp2_sum, zm2_sum, &
                  bx0, by0, bz0, &
                  !
                  nkpolar, &
                  rho_bin, &
                  u2_bin, ux2_bin, uy2_bin, uz2_bin, &
                  b2_bin, bx2_bin, by2_bin, bz2_bin, &
                  zp2_bin, zm2_bin &
                  !
                )

    !$acc exit data delete(ux, uy, uz, ux_old, uy_old, uz_old)
    !$acc exit data delete(rho, rho_old, u2half, u2half_old)
    !$acc exit data delete(u2, ux2, uy2, uz2, b2, bx2, by2, bz2)
    !$acc exit data delete(zp2, zm2, rho_abs)
    !$acc exit data delete(rho_r, ux_r, uy_r, uz_r, bx_r, by_r, bz_r)
    !$acc exit data delete(u2_bin, ux2_bin, uy2_bin, uz2_bin)
    !$acc exit data delete(b2_bin, bx2_bin, by2_bin, bz2_bin)
    !$acc exit data delete(rho_bin, zp2_bin, zm2_bin)

    deallocate(ux        )
    deallocate(uy        )
    deallocate(uz        )
    deallocate(ux_old    )
    deallocate(uy_old    )
    deallocate(uz_old    )
    deallocate(rho       )
    deallocate(rho_old   )
    deallocate(u2half    )
    deallocate(u2half_old)

    deallocate(u2         )
    deallocate(ux2        )
    deallocate(uy2        )
    deallocate(uz2        )
    deallocate(b2         )
    deallocate(bx2        )
    deallocate(by2        )
    deallocate(bz2        )
    deallocate(zp2        )
    deallocate(zm2        )
    deallocate(rho_abs    )

    deallocate( u2_bin)
    deallocate(ux2_bin)
    deallocate(uy2_bin)
    deallocate(uz2_bin)
    deallocate( b2_bin)
    deallocate(bx2_bin)
    deallocate(by2_bin)
    deallocate(bz2_bin)
    deallocate(rho_bin)
    deallocate(zp2_bin)
    deallocate(zm2_bin)

    deallocate(rho_r)
    deallocate( ux_r)
    deallocate( uy_r)
    deallocate( uz_r)
    deallocate( bx_r)
    deallocate( by_r)
    deallocate( bz_r)
  end subroutine loop_diagnostics


!---------------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Diagnostics for cross section of fileds
!---------------------------------------------------!
  subroutine loop_diagnostics_2D
    use io, only: loop_io_2D
    use utils, only: curl
    use mp, only: proc0, broadcast
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use fields, only: lnrho, lnrho_to_rho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use fields, only: m_to_u
    use cuFFTmp, only: btran_c2r
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: shear
    use shearing_box, only: to_non_shearing_coordinate, tsc, nremap
    implicit none

    complex(8), allocatable, dimension(:,:,:) :: rho
    real(8)   , allocatable, dimension(:,:,:) :: fr
    complex(8), allocatable, dimension(:,:,:) :: ux, uy, uz
    complex(8), allocatable, dimension(:,:,:) :: wx, wy, wz, jx, jy, jz
    real(8)   , allocatable, dimension(:,:)   :: rho_r_z0, rho_r_x0, rho_r_y0
    real(8)   , allocatable, dimension(:,:)   :: mx_r_z0, mx_r_x0, mx_r_y0
    real(8)   , allocatable, dimension(:,:)   :: my_r_z0, my_r_x0, my_r_y0
    real(8)   , allocatable, dimension(:,:)   :: mz_r_z0, mz_r_x0, mz_r_y0
    real(8)   , allocatable, dimension(:,:)   :: wx_r_z0, wx_r_x0, wx_r_y0
    real(8)   , allocatable, dimension(:,:)   :: wy_r_z0, wy_r_x0, wy_r_y0
    real(8)   , allocatable, dimension(:,:)   :: wz_r_z0, wz_r_x0, wz_r_y0
    real(8)   , allocatable, dimension(:,:)   :: bx_r_z0, bx_r_x0, bx_r_y0
    real(8)   , allocatable, dimension(:,:)   :: by_r_z0, by_r_x0, by_r_y0
    real(8)   , allocatable, dimension(:,:)   :: bz_r_z0, bz_r_x0, bz_r_y0
    real(8)   , allocatable, dimension(:,:)   :: jx_r_z0, jx_r_x0, jx_r_y0
    real(8)   , allocatable, dimension(:,:)   :: jy_r_z0, jy_r_x0, jy_r_y0
    real(8)   , allocatable, dimension(:,:)   :: jz_r_z0, jz_r_x0, jz_r_y0

    real(8)   , allocatable, dimension(:,:,:) :: u2, b2
    real(8)   , allocatable, dimension(:,:,:) :: rho_abs
    real(8)   , allocatable, dimension(:,:)   :: rho_kxy, rho_kyz, rho_kxz
    real(8)   , allocatable, dimension(:,:)   ::  u2_kxy,  u2_kyz,  u2_kxz
    real(8)   , allocatable, dimension(:,:)   ::  b2_kxy,  b2_kyz,  b2_kxz

    real(8)   , allocatable, dimension(:,:)   :: src1, src2, src3
    complex(8), allocatable, dimension(:,:,:) :: src4
    real(8)   , allocatable, dimension(:,:,:) :: src5

    integer :: i, j, k

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(fr(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(fr)

    allocate(src1(nlx_local, nly), source=0.d0)
    allocate(src2(nly      , nlz), source=0.d0)
    allocate(src3(nlx_local, nlz), source=0.d0)
    allocate(rho_r_z0, source=src1)
    allocate(rho_r_x0, source=src2)
    allocate(rho_r_y0, source=src3)

    allocate(mx_r_z0, source=src1)
    allocate(mx_r_x0, source=src2)
    allocate(mx_r_y0, source=src3)

    allocate(my_r_z0, source=src1)
    allocate(my_r_x0, source=src2)
    allocate(my_r_y0, source=src3)

    allocate(mz_r_z0, source=src1)
    allocate(mz_r_x0, source=src2)
    allocate(mz_r_y0, source=src3)

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

    allocate(src4(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(rho, source=src4)
    allocate(ux , source=src4)
    allocate(uy , source=src4)
    allocate(uz , source=src4)
    allocate(wx , source=src4)
    allocate(wy , source=src4)
    allocate(wz , source=src4)
    allocate(jx , source=src4)
    allocate(jy , source=src4)
    allocate(jz , source=src4)
    deallocate(src4)
    !$acc enter data create(rho, ux, uy, uz)
    !$acc enter data create(wx, wy, wz, jx, jy, jz)

    allocate(src5(nkz, nky_local, nkx), source=0.d0)
    allocate(u2     , source=src5)
    allocate(b2     , source=src5)
    allocate(rho_abs, source=src5)
    deallocate(src5)
    !$acc enter data create(u2, b2, rho_abs)

    allocate(src1(nkx      , nky_local), source=0.d0)
    allocate(src2(nky_local, nkz      ), source=0.d0)
    allocate(src3(nkx      , nkz      ), source=0.d0)
    allocate(rho_kxy, source=src1)
    allocate(rho_kyz, source=src2)
    allocate(rho_kxz, source=src3)
    allocate( u2_kxy, source=src1)
    allocate( u2_kyz, source=src2)
    allocate( u2_kxz, source=src3)
    allocate( b2_kxy, source=src1)
    allocate( b2_kyz, source=src2)
    allocate( b2_kxz, source=src3)
    deallocate(src1, src2, src3)

    !vvvvvvvvvvvvvvvvvv         2D cut of fields          vvvvvvvvvvvvvvvvvv!
    call lnrho_to_rho(lnrho, rho)
    call m_to_u(lnrho, mx, my, mz, ux, uy, uz)
    call curl(ux, uy, uz, wx, wy, wz)
    call curl(bx, by, bz, jx, jy, jz)

    ! btran_c2r preserves its input, so Calliope's copy into f is not needed;
    ! cut_2d_r runs on the host
    call btran_c2r(rho, fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, rho_r_z0, rho_r_x0, rho_r_y0)
    call btran_c2r(mx , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, mx_r_z0, mx_r_x0, mx_r_y0)
    call btran_c2r(my , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, my_r_z0, my_r_x0, my_r_y0)
    call btran_c2r(mz , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, mz_r_z0, mz_r_x0, mz_r_y0)
    call btran_c2r(bx , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, bx_r_z0, bx_r_x0, bx_r_y0)
    call btran_c2r(by , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, by_r_z0, by_r_x0, by_r_y0)
    call btran_c2r(bz , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, bz_r_z0, bz_r_x0, bz_r_y0)
    call btran_c2r(wx , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, wx_r_z0, wx_r_x0, wx_r_y0)
    call btran_c2r(wy , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, wy_r_z0, wy_r_x0, wy_r_y0)
    call btran_c2r(wz , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, wz_r_z0, wz_r_x0, wz_r_y0)
    call btran_c2r(jx , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, jx_r_z0, jx_r_x0, jx_r_y0)
    call btran_c2r(jy , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, jy_r_z0, jy_r_x0, jy_r_y0)
    call btran_c2r(jz , fr); if(shear) call to_non_shearing_coordinate(fr)
    !$acc update host(fr)
    call cut_2d_r(fr, jz_r_z0, jz_r_x0, jz_r_y0)

    ! x is split in real space: only iproc_fft = 0 (= proc0, the root of
    ! broadcast) has x = 0. The other ranks write one dummy element of the
    ! x = 0 plane in loop_io_2D; broadcasting the plane makes that element
    ! equal to the true value, so the overlapping MPI-IO write is harmless.
    call broadcast(rho_r_x0)
    call broadcast( mx_r_x0)
    call broadcast( my_r_x0)
    call broadcast( mz_r_x0)
    call broadcast( bx_r_x0)
    call broadcast( by_r_x0)
    call broadcast( bz_r_x0)
    call broadcast( wx_r_x0)
    call broadcast( wy_r_x0)
    call broadcast( wz_r_x0)
    call broadcast( jx_r_x0)
    call broadcast( jy_r_x0)
    call broadcast( jz_r_x0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv      2D spectra (avg over 1D)     vvvvvvvvvvvvvvvvvv!
    !$acc parallel loop collapse(3) gang vector default(present)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          u2(k, j, i) = 0.5d0*(abs(ux(k, j, i))**2 + abs(uy(k, j, i))**2 + abs(uz(k, j, i))**2)
          b2(k, j, i) = 0.5d0*(abs(bx(k, j, i))**2 + abs(by(k, j, i))**2 + abs(bz(k, j, i))**2)
          ! The reason for the following treatment for kz == 0 mode is the following. Compile it with LaTeX.
          !-----------------------------------------------------------------------------------------------------------------------------------
          ! The volume integral of a quadratic function is
          ! \int \mathrm{d}^3\mathbf{r}\, f(x,y,z)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2}
          ! \sum_{k_z = -n_{k_z}/2}^{n_{k_z}/2}|f_{k_x, k_y, k_z}|^2 = \left( \sum_{k_z = -n_{k_z}/2}^{-1}\sum_{k_x, k_y} + \sum_{k_z = 0}
          ! \sum_{k_x, k_y} + \sum_{k_z = 1}^{n_{k_z}/2}\sum_{k_x, k_y} \right) |f_{k_x, k_y, k_z}|^2
          ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
          !-----------------------------------------------------------------------------------------------------------------------------------
          if (k /= 1) then
            u2(k, j, i) = 2.0d0*u2(k, j, i)
            b2(k, j, i) = 2.0d0*b2(k, j, i)
          endif
          ! abs(rho) is summed without the kz weight, as in Calliope
          rho_abs(k, j, i) = abs(rho(k, j, i))
        end do
      end do
    end do
    !$acc end parallel loop

    ! sum_2d_k runs on the host
    !$acc update host(u2, b2, rho_abs)
    call sum_2d_k(rho_abs, rho_kxy, rho_kyz, rho_kxz)
    call sum_2d_k(     u2,  u2_kxy,  u2_kyz,  u2_kxz)
    call sum_2d_k(     b2,  b2_kxy,  b2_kyz,  b2_kxz)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io_2D( &
                  rho_r_z0, rho_r_x0, rho_r_y0, &
                  !
                  mx_r_z0, mx_r_x0, mx_r_y0, &
                  my_r_z0, my_r_x0, my_r_y0, &
                  mz_r_z0, mz_r_x0, mz_r_y0, &
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
                  rho_kxy, rho_kyz, rho_kxz, &
                   u2_kxy,  u2_kyz,  u2_kxz, &
                   b2_kxy,  b2_kyz,  b2_kxz  &
                )

    !$acc exit data delete(fr)
    !$acc exit data delete(rho, ux, uy, uz)
    !$acc exit data delete(wx, wy, wz, jx, jy, jz)
    !$acc exit data delete(u2, b2, rho_abs)

    deallocate(fr)
    deallocate(rho)
    deallocate(ux)
    deallocate(uy)
    deallocate(uz)
    deallocate(wx)
    deallocate(wy)
    deallocate(wz)
    deallocate(jx)
    deallocate(jy)
    deallocate(jz)

    deallocate(rho_r_z0)
    deallocate(rho_r_x0)
    deallocate(rho_r_y0)

    deallocate(mx_r_z0)
    deallocate(mx_r_x0)
    deallocate(mx_r_y0)

    deallocate(my_r_z0)
    deallocate(my_r_x0)
    deallocate(my_r_y0)

    deallocate(mz_r_z0)
    deallocate(mz_r_x0)
    deallocate(mz_r_y0)

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

    deallocate(u2)
    deallocate(b2)
    deallocate(rho_abs)

    deallocate(rho_kxy)
    deallocate(rho_kyz)
    deallocate(rho_kxz)
    deallocate( u2_kxy)
    deallocate( u2_kyz)
    deallocate( u2_kxz)
    deallocate( b2_kxy)
    deallocate( b2_kyz)
    deallocate( b2_kxz)
  end subroutine loop_diagnostics_2D


!-----------------------------------------------!
!> @author  YK
!! @date    8 Sep 2021
!! @brief   sum of spectrum over kz, kx, and ky
!!
!!          Gallope port: local copy of Calliope's
!!          diagnostics_common::sum_2d_k (Gallope's
!!          diagnostics_common has only cut_2d_k). Runs on
!!          the host on the (kz, ky, kx) layout. kz and kx
!!          are not split, so the sums over them are
!!          complete on every rank; ky is split over
!!          comm_fft, so the sum over ky is reduced to
!!          iproc_fft = 0 (which holds ky = 0 and writes the
!!          plane) and broadcast, so that the dummy element
!!          the other ranks write in loop_io_2D holds the
!!          same value.
!-----------------------------------------------!
  subroutine sum_2d_k(fk, fk_kxy, fk_kyz, fk_kxz)
    use grid, only: nkx, nky_local, nkz
    use mp, only: sum_reduce, broadcast, comm_fft
    implicit none
    real(8), dimension (:,:,:), intent(in) :: fk
    real(8), dimension (1:nkx      , 1:nky_local), intent(out) :: fk_kxy
    real(8), dimension (1:nky_local, 1:nkz      ), intent(out) :: fk_kyz
    real(8), dimension (1:nkx      , 1:nkz      ), intent(out) :: fk_kxz

    integer :: i, j, k

    fk_kxy(:, :) = 0.d0
    fk_kyz(:, :) = 0.d0
    fk_kxz(:, :) = 0.d0

    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fk_kxy(i, j) = fk_kxy(i, j) + fk(k, j, i)
          fk_kyz(j, k) = fk_kyz(j, k) + fk(k, j, i)
          fk_kxz(i, k) = fk_kxz(i, k) + fk(k, j, i)
        end do
      end do
    end do

    call sum_reduce(fk_kxz, 0, comm=comm_fft)
    call broadcast(fk_kxz)

  end subroutine sum_2d_k


!-----------------------------------------------!
!> @author  YK
!! @date    24 Mar 2021
!! @brief   get u2/2
!-----------------------------------------------!
  subroutine get_u2half(ux, uy, uz, u2half)
    use grid, only: nlx, nly, nlz
    use grid, only: nlx_local, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use cuFFTmp, only: btran_c2r, ftran_r2c
    implicit none
    complex(8), dimension (:,:,:), intent(in ) :: ux, uy, uz
    complex(8), dimension (:,:,:), intent(out) :: u2half
    real   (8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r, u2half_r
    real   (8), allocatable, dimension(:,:,:) :: src_r
    integer :: i, j, k

    allocate(src_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(ux_r    , source=src_r)
    allocate(uy_r    , source=src_r)
    allocate(uz_r    , source=src_r)
    allocate(u2half_r, source=src_r)
    deallocate(src_r)
    ! u2half_r is copied in so that its padding is zero
    !$acc enter data create(ux_r, uy_r, uz_r)
    !$acc enter data copyin(u2half_r)

    ! btran_c2r preserves its input, so Calliope's copies ux_ = ux, ... are not needed
    call btran_c2r(ux, ux_r)
    call btran_c2r(uy, uy_r)
    call btran_c2r(uz, uz_r)

    !$acc parallel loop collapse(3) gang vector default(present)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          u2half_r(k, j, i) = 0.5d0*(ux_r(k, j, i)**2 + uy_r(k, j, i)**2 + uz_r(k, j, i)**2)
        end do
      end do
    end do
    !$acc end parallel loop

    call ftran_r2c(u2half_r, u2half)

    !$acc parallel loop collapse(3) gang vector default(present)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          u2half(k, j, i) = u2half(k, j, i)/nlx/nly/nlz
        end do
      end do
    end do
    !$acc end parallel loop

    !$acc exit data delete(ux_r, uy_r, uz_r, u2half_r)
    deallocate(ux_r, uy_r, uz_r, u2half_r)

  end subroutine get_u2half


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Second order structure function
!-----------------------------------------------!
  subroutine loop_diagnostics_SF2
  ! under development...
  end subroutine loop_diagnostics_SF2


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Return kpar(k) & delta b/b0
!           k_\|(k) = \left(\frac{\langle|\mathbf{b}_{0,k} \cdot\nabla \delta\mathbf{b}_k|^2\rangle}
!           {\langle b_{0,k}^2\rangle\langle \delta b_k^2\rangle}\right)^{1/2}  \\
!           \mathbf{b}_{0,k}(\mathbf{x}) = \calF^{-1}\sum_{|\bm{k}|' \le k/2} \mathbf{b}_{\mathbf{k}'} \\
!           \delta\mathbf{b}_{k}(\mathbf{x}) = \calF^{-1}\sum_{k/2 \le |\bm{k}|' \le 2k} \mathbf{b}_{\mathbf{k}'}
!-----------------------------------------------!
  subroutine loop_diagnostics_kpar
  ! under development...
  end subroutine loop_diagnostics_kpar


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Calculate shell-to-shell transfer
!-----------------------------------------------!
  subroutine loop_diagnostics_nltrans
  ! under development...
  end subroutine loop_diagnostics_nltrans

end module diagnostics
