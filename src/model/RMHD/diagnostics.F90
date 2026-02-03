!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics for RMHD
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
!! @brief   Initialization of diagnostics
!-----------------------------------------------!
  subroutine init_diagnostics
    use params, only: inputfile
    use diagnostics_common, only: read_parameters
    use diagnostics_common, only: init_polar_spectrum_2d
    use io, only: init_io 
    implicit none

    call read_parameters(inputfile)

    call init_polar_spectrum_2d
    call init_io(nkpolar, kpbin)
  end subroutine init_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics in loop
!-----------------------------------------------!
  subroutine loop_diagnostics
    use io, only: loop_io
    use mp, only: proc0
    use grid, only: kprp2, kz2, kprp2_max, kz2_max
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx, nly, nlz_padded
    use fields, only: phi, psi
    use fields, only: phi_old, psi_old
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp
    use force, only: fphi, fpsi, fphi_old, fpsi_old
    use utils, only: cabs2
    implicit none
    integer :: i, j, k

    real(8), allocatable, dimension(:,:,:) :: upe2, bpe2
    real(8), allocatable, dimension(:,:,:) :: upe2old, bpe2old
    real(8), allocatable, dimension(:,:,:) :: upe2dissip_x, upe2dissip_z
    real(8), allocatable, dimension(:,:,:) :: bpe2dissip_x, bpe2dissip_z
    real(8), allocatable, dimension(:,:,:) :: p_phi, p_psi, p_xhl
    real(8), allocatable, dimension(:,:,:) :: zppe2, zmpe2
    real(8), allocatable, dimension(:,:,:) :: src

    real(8) :: upe2_sum, bpe2_sum
    real(8) :: upe2dot_sum, bpe2dot_sum
    real(8) :: upe2dissip_sum, bpe2dissip_sum
    real(8) :: zppe2_sum, zmpe2_sum
    real(8) :: p_phi_sum, p_psi_sum, p_xhl_sum

    real(8), dimension(:, :), allocatable :: upe2_bin , bpe2_bin      ! [kprp, kz]
    real(8), dimension(:, :), allocatable :: zppe2_bin, zmpe2_bin
    complex(8) :: phi_mid, psi_mid, jpa_mid, fphi_mid, fpsi_mid
    complex(8) :: zppe_mid, zmpe_mid, fzppe_mid, fzmpe_mid

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(src(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(upe2        , source=src)
    allocate(bpe2        , source=src)
    allocate(upe2old     , source=src)
    allocate(bpe2old     , source=src)
    allocate(upe2dissip_x, source=src)
    allocate(upe2dissip_z, source=src)
    allocate(bpe2dissip_x, source=src)
    allocate(bpe2dissip_z, source=src)
    allocate(p_phi       , source=src)
    allocate(p_psi       , source=src)
    allocate(p_xhl       , source=src)
    allocate(zppe2       , source=src)
    allocate(zmpe2       , source=src)
    deallocate(src)

    allocate (upe2_bin (1:nkpolar, nkz)); upe2_bin  = 0.d0
    allocate (bpe2_bin (1:nkpolar, nkz)); bpe2_bin  = 0.d0
    allocate (zppe2_bin(1:nkpolar, nkz)); zppe2_bin = 0.d0
    allocate (zmpe2_bin(1:nkpolar, nkz)); zmpe2_bin = 0.d0

    ! Create work arrays on device
    !$acc enter data create(upe2, bpe2, upe2old, bpe2old)
    !$acc enter data create(upe2dissip_x, upe2dissip_z, bpe2dissip_x, bpe2dissip_z)
    !$acc enter data create(p_phi, p_psi, p_xhl, zppe2, zmpe2)
    !$acc enter data create(upe2_bin, bpe2_bin)
    !$acc enter data create(zppe2_bin, zmpe2_bin)

    ! Initialize reduction variables
    upe2_sum       = 0.d0
    bpe2_sum       = 0.d0
    upe2dot_sum    = 0.d0
    bpe2dot_sum    = 0.d0
    upe2dissip_sum = 0.d0
    bpe2dissip_sum = 0.d0
    zppe2_sum      = 0.d0
    zmpe2_sum      = 0.d0
    p_phi_sum      = 0.d0
    p_psi_sum      = 0.d0
    p_xhl_sum      = 0.d0

    !$acc parallel loop collapse(3) gang vector default(present) &
    !$acc& private(phi_mid, psi_mid, jpa_mid, fphi_mid, fpsi_mid) &
    !$acc& private(zppe_mid, zmpe_mid, fzppe_mid, fzmpe_mid) &
    !$acc& reduction(+:upe2_sum, bpe2_sum, upe2dot_sum, bpe2dot_sum) &
    !$acc& reduction(+:upe2dissip_sum, bpe2dissip_sum) &
    !$acc& reduction(+:zppe2_sum, zmpe2_sum) &
    !$acc& reduction(+:p_phi_sum, p_psi_sum, p_xhl_sum)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          upe2   (k, j, i) = 0.5d0*cabs2(phi(k, j, i))*kprp2(k, j, i)
          bpe2   (k, j, i) = 0.5d0*cabs2(psi(k, j, i))*kprp2(k, j, i)

          upe2old(k, j, i) = 0.5d0*cabs2(phi_old(k, j, i))*kprp2(k, j, i)
          bpe2old(k, j, i) = 0.5d0*cabs2(psi_old(k, j, i))*kprp2(k, j, i)

          phi_mid  = 0.5d0*(phi (k, j, i) + phi_old (k, j, i))
          psi_mid  = 0.5d0*(psi (k, j, i) + psi_old (k, j, i))
          jpa_mid  = -kprp2(k, j, i)*psi_mid
          fphi_mid = 0.5d0*(fphi(k, j, i) + fphi_old(k, j, i))
          fpsi_mid = 0.5d0*(fpsi(k, j, i) + fpsi_old(k, j, i))

          upe2dissip_x(k, j, i) =  nupe_x*(kprp2(k, j, i)/kprp2_max)** nupe_x_exp*cabs2(phi_mid)*kprp2(k, j, i)
          upe2dissip_z(k, j, i) =  nupe_z*(kz2  (k)      /kz2_max  )** nupe_z_exp*cabs2(phi_mid)*kprp2(k, j, i)
          bpe2dissip_x(k, j, i) = etape_x*(kprp2(k, j, i)/kprp2_max)**etape_x_exp*cabs2(psi_mid)*kprp2(k, j, i)
          bpe2dissip_z(k, j, i) = etape_z*(kz2  (k)      /kz2_max  )**etape_z_exp*cabs2(psi_mid)*kprp2(k, j, i)
          p_phi       (k, j, i) = - 0.5d0*(-kprp2(k, j, i)*fphi_mid*conjg(phi_mid) + conjg(-kprp2(k, j, i)*fphi_mid)*phi_mid) 
          p_psi       (k, j, i) = - 0.5d0*(fpsi_mid*conjg(jpa_mid) + conjg(fpsi_mid)*jpa_mid) 
          zppe_mid  =  phi_mid +  psi_mid
          zmpe_mid  =  phi_mid -  psi_mid
          fzppe_mid = fphi_mid + fpsi_mid
          fzmpe_mid = fphi_mid - fpsi_mid
          zppe2       (k, j, i) = 0.5d0*cabs2(phi_mid + psi_mid)*kprp2(k, j, i)
          zmpe2       (k, j, i) = 0.5d0*cabs2(phi_mid - psi_mid)*kprp2(k, j, i)

          p_xhl       (k, j, i) = 0.25d0*( &
                                            zppe_mid*conjg(kprp2(k, j, i)*fzppe_mid) + conjg(zppe_mid)*kprp2(k, j, i)*fzppe_mid &
                                          - zmpe_mid*conjg(kprp2(k, j, i)*fzmpe_mid) - conjg(zmpe_mid)*kprp2(k, j, i)*fzmpe_mid &
                                        ) 

          ! The reason for the following treatment for kz == 0 mode is the following. Compile it with LaTeX.
          !-----------------------------------------------------------------------------------------------------------------------------------
          ! The volume integral of a quadratic function is
          ! \int \mathrm{d}^3\mathbf{r}\, f(x,y,z)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2}
          ! \sum_{k_z = -n_{k_z}/2}^{n_{k_z}/2}|f_{k_x, k_y, k_z}|^2 = \left( \sum_{k_z = -n_{k_z}/2}^{-1}\sum_{k_x, k_y} + \sum_{k_z = 0}
          ! \sum_{k_x, k_y} + \sum_{k_z = 1}^{n_{k_z}/2}\sum_{k_x, k_y} \right) |f_{k_x, k_y, k_z}|^2          !
          ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
          !-----------------------------------------------------------------------------------------------------------------------------------
          if (k /= 1) then
            upe2   (k, j, i) = 2.0d0*upe2   (k, j, i)
            bpe2   (k, j, i) = 2.0d0*bpe2   (k, j, i)

            upe2old(k, j, i) = 2.0d0*upe2old(k, j, i)
            bpe2old(k, j, i) = 2.0d0*bpe2old(k, j, i)

            upe2dissip_x(k, j, i) = 2.0d0*upe2dissip_x(k, j, i)
            upe2dissip_z(k, j, i) = 2.0d0*upe2dissip_z(k, j, i)
            bpe2dissip_x(k, j, i) = 2.0d0*bpe2dissip_x(k, j, i)
            bpe2dissip_z(k, j, i) = 2.0d0*bpe2dissip_z(k, j, i)
            zppe2       (k, j, i) = 2.0d0*zppe2       (k, j, i)
            zmpe2       (k, j, i) = 2.0d0*zmpe2       (k, j, i)
            p_phi       (k, j, i) = 2.0d0*p_phi       (k, j, i)
            p_psi       (k, j, i) = 2.0d0*p_psi       (k, j, i)
            p_xhl       (k, j, i) = 2.0d0*p_xhl       (k, j, i)
          endif

          upe2_sum       = upe2_sum       + upe2       (k, j, i)
          bpe2_sum       = bpe2_sum       + bpe2       (k, j, i)
          upe2dot_sum    = upe2dot_sum    + (upe2(k, j, i) - upe2old(k, j, i))/dt
          bpe2dot_sum    = bpe2dot_sum    + (bpe2(k, j, i) - bpe2old(k, j, i))/dt
          upe2dissip_sum = upe2dissip_sum + upe2dissip_x (k, j, i) + upe2dissip_z (k, j, i)
          bpe2dissip_sum = bpe2dissip_sum + bpe2dissip_x (k, j, i) + bpe2dissip_z (k, j, i)
          zppe2_sum      = zppe2_sum      + zppe2      (k, j, i)
          zmpe2_sum      = zmpe2_sum      + zmpe2      (k, j, i)
          p_phi_sum      = p_phi_sum      + p_phi(k, j, i)
          p_psi_sum      = p_psi_sum      + p_psi(k, j, i)
          p_xhl_sum      = p_xhl_sum      + p_xhl(k, j, i)

        end do
      end do
    end do
    !$acc end parallel loop

    ! MPI reduction for summed quantities (operates on host scalars)
    call sum_reduce(upe2_sum       , 0)
    call sum_reduce(bpe2_sum       , 0)
    call sum_reduce(upe2dot_sum    , 0)
    call sum_reduce(bpe2dot_sum    , 0)
    call sum_reduce(upe2dissip_sum , 0)
    call sum_reduce(bpe2dissip_sum , 0)
    call sum_reduce(zppe2_sum      , 0)
    call sum_reduce(zmpe2_sum      , 0)
    call sum_reduce(p_phi_sum      , 0)
    call sum_reduce(p_psi_sum      , 0)
    call sum_reduce(p_xhl_sum      , 0)

    ! Bin spectra over kprp on device
    call get_polar_spectrum_2d(upe2 , upe2_bin )
    call get_polar_spectrum_2d(bpe2 , bpe2_bin )
    call get_polar_spectrum_2d(zppe2, zppe2_bin)
    call get_polar_spectrum_2d(zmpe2, zmpe2_bin)

    ! Transfer binned spectra from device to host for I/O
    !$acc update host(upe2_bin, bpe2_bin)
    !$acc update host(zppe2_bin, zmpe2_bin)
  
    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  upe2_sum, bpe2_sum, &
                  upe2dot_sum, bpe2dot_sum, &
                  upe2dissip_sum, bpe2dissip_sum, &
                  p_phi_sum, p_psi_sum, p_xhl_sum, &
                  zppe2_sum, zmpe2_sum, &
                  !
                  nkpolar, &
                  upe2_bin , bpe2_bin , &
                  zppe2_bin, zmpe2_bin  &
                )

    !$acc exit data delete(upe2, bpe2)
    !$acc exit data delete(upe2old, bpe2old)
    !$acc exit data delete(upe2dissip_x, upe2dissip_z, bpe2dissip_x, bpe2dissip_z)
    !$acc exit data delete(p_phi, p_psi, p_xhl)
    !$acc exit data delete(zppe2, zmpe2)
    !$acc exit data delete(upe2_bin, bpe2_bin)
    !$acc exit data delete(zppe2_bin, zmpe2_bin)
    deallocate(upe2)
    deallocate(bpe2)
    deallocate(upe2old)
    deallocate(bpe2old)
    deallocate(upe2dissip_x)
    deallocate(upe2dissip_z)
    deallocate(bpe2dissip_x)
    deallocate(bpe2dissip_z)
    deallocate(p_phi)
    deallocate(p_psi)
    deallocate(p_xhl)
    deallocate(zppe2)
    deallocate(zmpe2)

    deallocate (upe2_bin)
    deallocate (bpe2_bin)
    deallocate (zppe2_bin)
    deallocate (zmpe2_bin)
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics for cross section of fileds
!-----------------------------------------------!
  subroutine loop_diagnostics_2D
    use io, only: loop_io_2D
    use mp, only: proc0
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: kx, ky, kprp2
    use grid, only: nkx, nky_local, nkz
    use fields, only: phi, psi
    use cuFFTmp, only: btran_c2r
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi
    implicit none
    integer :: i, j, k

    complex(8), allocatable, dimension(:,:,:) :: f
    real(8)   , allocatable, dimension(:,:,:) :: fr
    real(8)   , allocatable, dimension(:,:)   :: phi_r_z0, phi_r_x0, phi_r_y0
    real(8)   , allocatable, dimension(:,:)   :: psi_r_z0, psi_r_x0, psi_r_y0
    real(8)   , allocatable, dimension(:,:)   :: omg_r_z0, omg_r_x0, omg_r_y0
    real(8)   , allocatable, dimension(:,:)   :: jpa_r_z0, jpa_r_x0, jpa_r_y0
    real(8)   , allocatable, dimension(:,:)   ::  ux_r_z0,  ux_r_x0,  ux_r_y0
    real(8)   , allocatable, dimension(:,:)   ::  uy_r_z0,  uy_r_x0,  uy_r_y0
    real(8)   , allocatable, dimension(:,:)   ::  bx_r_z0,  bx_r_x0,  bx_r_y0
    real(8)   , allocatable, dimension(:,:)   ::  by_r_z0,  by_r_x0,  by_r_y0
    complex(8), allocatable, dimension(:,:,:) :: omg, jpa, ux, uy, bx, by

    real(8)   , allocatable, dimension(:,:)   :: src1, src2, src3
    complex(8), allocatable, dimension(:,:,:) :: src4

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(f (nkz       , nky_local, nkx      )); f   = 0.d0
    allocate(fr(nlz_padded, nly      , nlx_local)); fr  = 0.d0
    !$acc enter data create(f )
    !$acc enter data create(fr)

    allocate(src1(nlx_local, nly), source=0.d0) 
    allocate(src2(nly      , nlz), source=0.d0) 
    allocate(src3(nlx_local, nlz), source=0.d0)
    allocate(phi_r_z0, source=src1)
    allocate(phi_r_x0, source=src2)
    allocate(phi_r_y0, source=src3)

    allocate(psi_r_z0, source=src1)
    allocate(psi_r_x0, source=src2)
    allocate(psi_r_y0, source=src3)

    allocate(omg_r_z0, source=src1)
    allocate(omg_r_x0, source=src2)
    allocate(omg_r_y0, source=src3)

    allocate(jpa_r_z0, source=src1)
    allocate(jpa_r_x0, source=src2)
    allocate(jpa_r_y0, source=src3)

    allocate(ux_r_z0, source=src1)
    allocate(ux_r_x0, source=src2)
    allocate(ux_r_y0, source=src3)

    allocate(uy_r_z0, source=src1)
    allocate(uy_r_x0, source=src2)
    allocate(uy_r_y0, source=src3)

    allocate(bx_r_z0, source=src1)
    allocate(bx_r_x0, source=src2)
    allocate(bx_r_y0, source=src3)

    allocate(by_r_z0, source=src1)
    allocate(by_r_x0, source=src2)
    allocate(by_r_y0, source=src3)
    deallocate(src1, src2, src3)

    allocate(src4(nkz, nky_local, nkx), source=(0.d0,0.d0))
    allocate(omg, source=src4)
    allocate(jpa, source=src4)
    allocate(ux , source=src4)
    allocate(uy , source=src4)
    allocate(bx , source=src4)
    allocate(by , source=src4)
    !$acc enter data create(omg)
    !$acc enter data create(jpa)
    !$acc enter data create(ux )
    !$acc enter data create(uy )
    !$acc enter data create(bx )
    !$acc enter data create(by )
    deallocate(src4)

    !vvvvvvvvvvvvvvvvvv         2D cut of fields          vvvvvvvvvvvvvvvvvv!
    !$acc data present(phi, psi, omg, jpa, ux, uy, bx, kx, ky, kprp2)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          omg(k, j, i) = -kprp2(k, j, i)*phi(k, j, i)
          jpa(k, j, i) = -kprp2(k, j, i)*psi(k, j, i)
          ux (k, j, i) = -zi*ky(j)*phi(k, j, i)
          uy (k, j, i) =  zi*kx(i)*phi(k, j, i)
          bx (k, j, i) = -zi*ky(j)*psi(k, j, i)
          by (k, j, i) =  zi*kx(i)*psi(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data

    !$acc kernels
    f = phi 
    !$acc end kernels
    call btran_c2r(f, fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, phi_r_z0, phi_r_x0, phi_r_y0)

    !$acc kernels
    f = psi 
    !$acc end kernels
    call btran_c2r(f, fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, psi_r_z0, psi_r_x0, psi_r_y0)

    !$acc kernels
    f = omg 
    !$acc end kernels
    call btran_c2r(f, fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, omg_r_z0, omg_r_x0, omg_r_y0)

    !$acc kernels
    f = jpa 
    !$acc end kernels
    call btran_c2r(f, fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, jpa_r_z0, jpa_r_x0, jpa_r_y0)

    !$acc kernels
    f = ux 
    !$acc end kernels
    call btran_c2r(f, fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, ux_r_z0, ux_r_x0, ux_r_y0)

    !$acc kernels
    f = uy 
    !$acc end kernels
    call btran_c2r(f, fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, uy_r_z0, uy_r_x0, uy_r_y0)

    !$acc kernels
    f = bx 
    !$acc end kernels
    call btran_c2r(f, fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, bx_r_z0, bx_r_x0, bx_r_y0)

    !$acc kernels
    f = by 
    !$acc end kernels
    call btran_c2r(f, fr) 
    !$acc update host(fr)
    call cut_2d_r(fr, by_r_z0, by_r_x0, by_r_y0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io_2D( &
                  phi_r_z0, phi_r_x0, phi_r_y0, &
                  psi_r_z0, psi_r_x0, psi_r_y0, &
                  omg_r_z0, omg_r_x0, omg_r_y0, &
                  jpa_r_z0, jpa_r_x0, jpa_r_y0, &
                   ux_r_z0,  ux_r_x0,  ux_r_y0, &
                   uy_r_z0,  uy_r_x0,  uy_r_y0, &
                   bx_r_z0,  bx_r_x0,  bx_r_y0, &
                   by_r_z0,  by_r_x0,  by_r_y0  &
                )

    !$acc exit data delete(f )
    !$acc exit data delete(fr)
    deallocate(f)
    deallocate(fr)

    deallocate (phi_r_z0)
    deallocate (phi_r_x0)
    deallocate (phi_r_y0)

    deallocate (psi_r_z0)
    deallocate (psi_r_x0)
    deallocate (psi_r_y0)

    deallocate (omg_r_z0)
    deallocate (omg_r_x0)
    deallocate (omg_r_y0)

    deallocate (jpa_r_z0)
    deallocate (jpa_r_x0)
    deallocate (jpa_r_y0)

    deallocate ( ux_r_z0)
    deallocate ( ux_r_x0)
    deallocate ( ux_r_y0)

    deallocate ( uy_r_z0)
    deallocate ( uy_r_x0)
    deallocate ( uy_r_y0)

    deallocate ( bx_r_z0)
    deallocate ( bx_r_x0)
    deallocate ( bx_r_y0)

    deallocate ( by_r_z0)
    deallocate ( by_r_x0)
    deallocate ( by_r_y0)

    !$acc exit data delete(omg)
    !$acc exit data delete(jpa)
    !$acc exit data delete(ux )
    !$acc exit data delete(uy )
    !$acc exit data delete(bx )
    !$acc exit data delete(by )
    deallocate (omg)
    deallocate (jpa)
    deallocate (ux)
    deallocate (uy)
    deallocate (bx)
    deallocate (by)
  end subroutine loop_diagnostics_2D


!-----------------------------------------------!
!> @author  YK
!! @brief   Second order structure function
!-----------------------------------------------!
  subroutine loop_diagnostics_SF2
  ! under development...
  end subroutine loop_diagnostics_SF2


!-----------------------------------------------!
!> @author  YK
!! @brief   Calculate kpar(k) & delta b/b0
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
!! @brief   Calculate shell-to-shell transfer
!-----------------------------------------------!
  subroutine loop_diagnostics_nltrans
  ! under development...
  end subroutine loop_diagnostics_nltrans

end module diagnostics




