!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics for HRMHD
!!
!!          Invariants (ideal, unforced):
!!            W = 1/2 <|grad phi|^2 + |grad psi|^2 + upa^2 + bpa^2>
!!            X = <grad phi.grad psi + upa*bpa/sgm + (rho/sgm)*upa*lap phi>
!!
!!          Fourier sums over the kz >= 0 half are doubled
!!          for kz /= 0 (index k /= 1) to recover the full sum.
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
!!
!!          Energies, the Elsasser energies and X are
!!          taken at the current time; dissipation and
!!          injection at the midpoint of the last step,
!!          so that d/dt(energy) = injection - dissipation
!!          holds to O(dt^2).
!!
!!          KAW/ICW: W-orthogonal split of (phi, bpa)
!!          and (psi, upa) onto the large-k_perp*rho
!!          eigenvectors (b = rho*kprp2*phi, u = sgm*psi/rho);
!!          ICW = total - KAW. The split is only defined
!!          for rho > 0 and sgm > 0; otherwise its outputs
!!          are zero.
!-----------------------------------------------!
  subroutine loop_diagnostics
    use io, only: loop_io, nseries, nbins
    use io, only: is_upe2, is_bpe2, is_upa2, is_bpa2
    use io, only: is_upe2dot, is_bpe2dot, is_upa2dot, is_bpa2dot
    use io, only: is_upe2dissip, is_bpe2dissip, is_upa2dissip, is_bpa2dissip
    use io, only: is_p_aw, is_p_compr, is_p_xhl
    use io, only: is_zppe2, is_zmpe2, is_zppa2, is_zmpa2
    use io, only: is_upe2_KAW_dissip, is_bpe2_KAW_dissip, is_upa2_KAW_dissip, is_bpa2_KAW_dissip
    use io, only: is_upe2_ICW_dissip, is_bpe2_ICW_dissip, is_upa2_ICW_dissip, is_bpa2_ICW_dissip
    use io, only: is_KAW_dissip_prp, is_KAW_dissip_par, is_ICW_dissip_prp, is_ICW_dissip_par
    use io, only: is_xhl
    use io, only: ib_upe2, ib_bpe2, ib_upa2, ib_bpa2
    use io, only: ib_ux2, ib_uy2, ib_bx2, ib_by2
    use io, only: ib_p_aw, ib_p_compr
    use io, only: ib_dissip_aw, ib_dissip_compr, ib_dissip_KAW, ib_dissip_ICW
    use io, only: ib_zppe2, ib_zmpe2, ib_zppa2, ib_zmpa2
    use io, only: ib_upe2_KAW, ib_bpe2_KAW, ib_upa2_KAW, ib_bpa2_KAW
    use io, only: ib_upe2_ICW, ib_bpe2_ICW, ib_upa2_ICW, ib_bpa2_ICW
    use io, only: ib_flux_W, ib_flux_X
    use mp, only: proc0
    use grid, only: kprp2, kz2, kprp2_max, kz2_max
    use grid, only: kx, ky
    use grid, only: nkx, nky_local, nkz
    use fields, only: phi, psi, upa, bpa
    use fields, only: phi_old, psi_old, upa_old, bpa_old
    use mp, only: sum_reduce, comm_fft, nproc_m, nproc_s
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: rho, sgm, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp, &
                      nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                      etapa_x, etapa_x_exp, etapa_z, etapa_z_exp
    use force, only: fphi, fpsi, fupa, fbpa
    use force, only: fphi_old, fpsi_old, fupa_old, fbpa_old
    use utils, only: cabs2
    implicit none
    ! Binned quantities computed per mode here; the two fluxes (the last
    ! two bins) come from get_flux_WX.
    integer, parameter :: nq = ib_flux_W - 1
    integer :: i, j, k, l
    logical :: split_ok

    real(8) :: series(nseries)
    real(8), allocatable, dimension(:,:,:)   :: bins   ! (kpbin, kz, bin) on host
    real(8), allocatable, dimension(:,:,:,:) :: qq     ! (kz, ky, kx, bin) per-mode work
    real(8), allocatable, dimension(:,:)     :: ebin   ! (kpbin, kz)

    ! per-mode scalars
    real(8)    :: fac, kp2
    real(8)    :: rx_upe, rz_upe, rx_bpe, rz_bpe, rx_upa, rz_upa, rx_bpa, rz_bpa
    real(8)    :: e_upe, e_bpe, e_upa, e_bpa, eo_upe, eo_bpe, eo_upa, eo_bpa
    real(8)    :: dx_upe, dz_upe, dx_bpe, dz_bpe, dx_upa, dz_upa, dx_bpa, dz_bpa
    real(8)    :: p_aw, p_compr
    complex(8) :: phi_m, psi_m, upa_m, bpa_m
    complex(8) :: fphi_m, fpsi_m, fupa_m, fbpa_m
    complex(8) :: a, c
    complex(8) :: phi_K, psi_K, upa_K, bpa_K
    complex(8) :: phi_I, psi_I, upa_I, bpa_I
    real(8)    :: dKx_upe, dKz_upe, dKx_bpe, dKz_bpe, dKx_upa, dKz_upa, dKx_bpa, dKz_bpa
    real(8)    :: dIx_upe, dIz_upe, dIx_bpe, dIz_bpe, dIx_upa, dIz_upa, dIx_bpa, dIz_bpa

    ! reduction variables
    real(8) :: upe2_sum, bpe2_sum, upa2_sum, bpa2_sum
    real(8) :: upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum
    real(8) :: upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum
    real(8) :: p_aw_sum, p_compr_sum, p_xhl_sum
    real(8) :: zppe2_sum, zmpe2_sum, zppa2_sum, zmpa2_sum
    real(8) :: xhl_sum
    real(8) :: upe2_KAW_dissip_sum, bpe2_KAW_dissip_sum, upa2_KAW_dissip_sum, bpa2_KAW_dissip_sum
    real(8) :: upe2_ICW_dissip_sum, bpe2_ICW_dissip_sum, upa2_ICW_dissip_sum, bpa2_ICW_dissip_sum
    real(8) :: KAW_dissip_prp_sum, KAW_dissip_par_sum, ICW_dissip_prp_sum, ICW_dissip_par_sum

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(bins(1:nkpolar, nkz, nbins), source=0.d0)

    ! Fluxes of W and X first: they need the nonlinear-term buffers, which
    ! are released before the per-mode work array below is allocated.
    call get_flux_WX(bins(:,:,ib_flux_W), bins(:,:,ib_flux_X))

    allocate(qq(nkz, nky_local, nkx, nq))
    allocate(ebin(1:nkpolar, nkz))
    !$acc enter data create(qq, ebin)

    ! Slots that a kernel below does not write (KAW/ICW when the split is
    ! undefined) stay zero.
    !$acc kernels present(qq)
    qq = 0.d0
    !$acc end kernels

    upe2_sum       = 0.d0; bpe2_sum       = 0.d0; upa2_sum       = 0.d0; bpa2_sum       = 0.d0
    upe2dot_sum    = 0.d0; bpe2dot_sum    = 0.d0; upa2dot_sum    = 0.d0; bpa2dot_sum    = 0.d0
    upe2dissip_sum = 0.d0; bpe2dissip_sum = 0.d0; upa2dissip_sum = 0.d0; bpa2dissip_sum = 0.d0
    p_aw_sum       = 0.d0; p_compr_sum    = 0.d0; p_xhl_sum      = 0.d0
    zppe2_sum      = 0.d0; zmpe2_sum      = 0.d0; zppa2_sum      = 0.d0; zmpa2_sum      = 0.d0
    xhl_sum        = 0.d0

    !vvvvvvvvvvvvvvvv  energies, rates, injection, Elsasser, X  vvvvvvvvvvvvvvvv!
    !$acc parallel loop collapse(3) gang vector default(present) &
    !$acc& private(fac, kp2, rx_upe, rz_upe, rx_bpe, rz_bpe, rx_upa, rz_upa, rx_bpa, rz_bpa) &
    !$acc& private(e_upe, e_bpe, e_upa, e_bpa, eo_upe, eo_bpe, eo_upa, eo_bpa) &
    !$acc& private(dx_upe, dz_upe, dx_bpe, dz_bpe, dx_upa, dz_upa, dx_bpa, dz_bpa) &
    !$acc& private(p_aw, p_compr) &
    !$acc& private(phi_m, psi_m, upa_m, bpa_m, fphi_m, fpsi_m, fupa_m, fbpa_m) &
    !$acc& reduction(+:upe2_sum, bpe2_sum, upa2_sum, bpa2_sum) &
    !$acc& reduction(+:upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum) &
    !$acc& reduction(+:upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum) &
    !$acc& reduction(+:p_aw_sum, p_compr_sum, p_xhl_sum) &
    !$acc& reduction(+:zppe2_sum, zmpe2_sum, zppa2_sum, zmpa2_sum, xhl_sum)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ! kz /= 0 modes stand for their kz < 0 conjugates as well
          fac = 1.d0
          if (k /= 1) fac = 2.d0
          kp2 = kprp2(k, j, i)

          ! energy-dissipation coefficients, same as advance::get_imp_rates
          rx_upe = nupe_x *(kp2/kprp2_max)**nupe_x_exp ; rz_upe = nupe_z *(kz2(k)/kz2_max)**nupe_z_exp
          rx_bpe = etape_x*(kp2/kprp2_max)**etape_x_exp; rz_bpe = etape_z*(kz2(k)/kz2_max)**etape_z_exp
          rx_upa = nupa_x *(kp2/kprp2_max)**nupa_x_exp ; rz_upa = nupa_z *(kz2(k)/kz2_max)**nupa_z_exp
          rx_bpa = etapa_x*(kp2/kprp2_max)**etapa_x_exp; rz_bpa = etapa_z*(kz2(k)/kz2_max)**etapa_z_exp

          e_upe  = 0.5d0*cabs2(phi    (k, j, i))*kp2
          e_bpe  = 0.5d0*cabs2(psi    (k, j, i))*kp2
          e_upa  = 0.5d0*cabs2(upa    (k, j, i))
          e_bpa  = 0.5d0*cabs2(bpa    (k, j, i))
          eo_upe = 0.5d0*cabs2(phi_old(k, j, i))*kp2
          eo_bpe = 0.5d0*cabs2(psi_old(k, j, i))*kp2
          eo_upa = 0.5d0*cabs2(upa_old(k, j, i))
          eo_bpa = 0.5d0*cabs2(bpa_old(k, j, i))

          phi_m  = 0.5d0*(phi (k, j, i) + phi_old (k, j, i))
          psi_m  = 0.5d0*(psi (k, j, i) + psi_old (k, j, i))
          upa_m  = 0.5d0*(upa (k, j, i) + upa_old (k, j, i))
          bpa_m  = 0.5d0*(bpa (k, j, i) + bpa_old (k, j, i))
          fphi_m = 0.5d0*(fphi(k, j, i) + fphi_old(k, j, i))
          fpsi_m = 0.5d0*(fpsi(k, j, i) + fpsi_old(k, j, i))
          fupa_m = 0.5d0*(fupa(k, j, i) + fupa_old(k, j, i))
          fbpa_m = 0.5d0*(fbpa(k, j, i) + fbpa_old(k, j, i))

          ! no factor 2: d/dt (1/2 |f|^2) = -r |f|^2
          dx_upe = rx_upe*cabs2(phi_m)*kp2; dz_upe = rz_upe*cabs2(phi_m)*kp2
          dx_bpe = rx_bpe*cabs2(psi_m)*kp2; dz_bpe = rz_bpe*cabs2(psi_m)*kp2
          dx_upa = rx_upa*cabs2(upa_m)    ; dz_upa = rz_upa*cabs2(upa_m)
          dx_bpa = rx_bpa*cabs2(bpa_m)    ; dz_bpa = rz_bpa*cabs2(bpa_m)

          p_aw    = kp2*(dble(fphi_m*conjg(phi_m)) + dble(fpsi_m*conjg(psi_m)))
          p_compr =      dble(fupa_m*conjg(upa_m)) + dble(fbpa_m*conjg(bpa_m))

          upe2_sum       = upe2_sum       + fac*e_upe
          bpe2_sum       = bpe2_sum       + fac*e_bpe
          upa2_sum       = upa2_sum       + fac*e_upa
          bpa2_sum       = bpa2_sum       + fac*e_bpa
          upe2dot_sum    = upe2dot_sum    + fac*(e_upe - eo_upe)/dt
          bpe2dot_sum    = bpe2dot_sum    + fac*(e_bpe - eo_bpe)/dt
          upa2dot_sum    = upa2dot_sum    + fac*(e_upa - eo_upa)/dt
          bpa2dot_sum    = bpa2dot_sum    + fac*(e_bpa - eo_bpa)/dt
          upe2dissip_sum = upe2dissip_sum + fac*(dx_upe + dz_upe)
          bpe2dissip_sum = bpe2dissip_sum + fac*(dx_bpe + dz_bpe)
          upa2dissip_sum = upa2dissip_sum + fac*(dx_upa + dz_upa)
          bpa2dissip_sum = bpa2dissip_sum + fac*(dx_bpa + dz_bpa)
          p_aw_sum       = p_aw_sum       + fac*p_aw
          p_compr_sum    = p_compr_sum    + fac*p_compr

          ! Elsasser energies (no 1/2, as in Calliope HRMHD) and the
          ! perpendicular part of the X injection,
          ! 1/4 (z+ conj(kprp2 fz+) + c.c. - z- conj(kprp2 fz-) - c.c.)
          zppe2_sum = zppe2_sum + fac*cabs2(phi_m + psi_m)*kp2
          zmpe2_sum = zmpe2_sum + fac*cabs2(phi_m - psi_m)*kp2
          zppa2_sum = zppa2_sum + fac*cabs2(upa_m + bpa_m)
          zmpa2_sum = zmpa2_sum + fac*cabs2(upa_m - bpa_m)
          p_xhl_sum = p_xhl_sum + fac*0.5d0*kp2*( dble((phi_m + psi_m)*conjg(fphi_m + fpsi_m)) &
                                                - dble((phi_m - psi_m)*conjg(fphi_m - fpsi_m)) )

          ! X = kprp2 Re(phi psi*) + Re(upa bpa*)/sgm - (rho/sgm) kprp2 Re(upa phi*)
          if (sgm > 0.d0) then
            xhl_sum = xhl_sum + fac*( kp2*dble(phi(k, j, i)*conjg(psi(k, j, i))) &
                                    + ( dble(upa(k, j, i)*conjg(bpa(k, j, i))) &
                                      - rho*kp2*dble(upa(k, j, i)*conjg(phi(k, j, i))) )/sgm )
          endif

          qq(k, j, i, ib_upe2)         = fac*e_upe
          qq(k, j, i, ib_bpe2)         = fac*e_bpe
          qq(k, j, i, ib_upa2)         = fac*e_upa
          qq(k, j, i, ib_bpa2)         = fac*e_bpa
          qq(k, j, i, ib_ux2)          = fac*0.5d0*ky(j)**2*cabs2(phi(k, j, i))
          qq(k, j, i, ib_uy2)          = fac*0.5d0*kx(i)**2*cabs2(phi(k, j, i))
          qq(k, j, i, ib_bx2)          = fac*0.5d0*ky(j)**2*cabs2(psi(k, j, i))
          qq(k, j, i, ib_by2)          = fac*0.5d0*kx(i)**2*cabs2(psi(k, j, i))
          qq(k, j, i, ib_p_aw)         = fac*p_aw
          qq(k, j, i, ib_p_compr)      = fac*p_compr
          qq(k, j, i, ib_dissip_aw)    = fac*(dx_upe + dz_upe + dx_bpe + dz_bpe)
          qq(k, j, i, ib_dissip_compr) = fac*(dx_upa + dz_upa + dx_bpa + dz_bpa)
          qq(k, j, i, ib_zppe2)        = fac*cabs2(phi_m + psi_m)*kp2
          qq(k, j, i, ib_zmpe2)        = fac*cabs2(phi_m - psi_m)*kp2
          qq(k, j, i, ib_zppa2)        = fac*cabs2(upa_m + bpa_m)
          qq(k, j, i, ib_zmpa2)        = fac*cabs2(upa_m - bpa_m)
        end do
      end do
    end do
    !$acc end parallel loop
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    upe2_KAW_dissip_sum = 0.d0; bpe2_KAW_dissip_sum = 0.d0
    upa2_KAW_dissip_sum = 0.d0; bpa2_KAW_dissip_sum = 0.d0
    upe2_ICW_dissip_sum = 0.d0; bpe2_ICW_dissip_sum = 0.d0
    upa2_ICW_dissip_sum = 0.d0; bpa2_ICW_dissip_sum = 0.d0
    KAW_dissip_prp_sum  = 0.d0; KAW_dissip_par_sum  = 0.d0
    ICW_dissip_prp_sum  = 0.d0; ICW_dissip_par_sum  = 0.d0

    split_ok = rho > 0.d0 .and. sgm > 0.d0

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvv  KAW/ICW split  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    ! Energies are split at the current time, dissipation at the midpoint
    ! (the same times as the totals above). The W-orthogonal split makes
    ! the KAW and ICW energies add up to the total; the dissipation does
    ! not add up in general, since the dissipation coefficients of the four
    ! fields differ.
    if (split_ok) then
      !$acc parallel loop collapse(3) gang vector default(present) &
      !$acc& private(fac, kp2, rx_upe, rz_upe, rx_bpe, rz_bpe, rx_upa, rz_upa, rx_bpa, rz_bpa) &
      !$acc& private(a, c, phi_K, psi_K, upa_K, bpa_K, phi_I, psi_I, upa_I, bpa_I) &
      !$acc& private(phi_m, psi_m, upa_m, bpa_m) &
      !$acc& private(dKx_upe, dKz_upe, dKx_bpe, dKz_bpe, dKx_upa, dKz_upa, dKx_bpa, dKz_bpa) &
      !$acc& private(dIx_upe, dIz_upe, dIx_bpe, dIz_bpe, dIx_upa, dIz_upa, dIx_bpa, dIz_bpa) &
      !$acc& reduction(+:upe2_KAW_dissip_sum, bpe2_KAW_dissip_sum, upa2_KAW_dissip_sum, bpa2_KAW_dissip_sum) &
      !$acc& reduction(+:upe2_ICW_dissip_sum, bpe2_ICW_dissip_sum, upa2_ICW_dissip_sum, bpa2_ICW_dissip_sum) &
      !$acc& reduction(+:KAW_dissip_prp_sum, KAW_dissip_par_sum, ICW_dissip_prp_sum, ICW_dissip_par_sum)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            fac = 1.d0
            if (k /= 1) fac = 2.d0
            kp2 = kprp2(k, j, i)

            rx_upe = nupe_x *(kp2/kprp2_max)**nupe_x_exp ; rz_upe = nupe_z *(kz2(k)/kz2_max)**nupe_z_exp
            rx_bpe = etape_x*(kp2/kprp2_max)**etape_x_exp; rz_bpe = etape_z*(kz2(k)/kz2_max)**etape_z_exp
            rx_upa = nupa_x *(kp2/kprp2_max)**nupa_x_exp ; rz_upa = nupa_z *(kz2(k)/kz2_max)**nupa_z_exp
            rx_bpa = etapa_x*(kp2/kprp2_max)**etapa_x_exp; rz_bpa = etapa_z*(kz2(k)/kz2_max)**etapa_z_exp

            ! current time: energies
            a = (phi(k, j, i) + rho*bpa(k, j, i))/(1.d0 + rho**2*kp2)
            c = (rho*kp2*psi(k, j, i) + sgm*upa(k, j, i))/(sgm**2 + rho**2*kp2)
            phi_K = a    ; bpa_K = rho*kp2*a
            psi_K = rho*c; upa_K = sgm*c
            phi_I = phi(k, j, i) - phi_K; psi_I = psi(k, j, i) - psi_K
            upa_I = upa(k, j, i) - upa_K; bpa_I = bpa(k, j, i) - bpa_K

            qq(k, j, i, ib_upe2_KAW) = fac*0.5d0*cabs2(phi_K)*kp2
            qq(k, j, i, ib_bpe2_KAW) = fac*0.5d0*cabs2(psi_K)*kp2
            qq(k, j, i, ib_upa2_KAW) = fac*0.5d0*cabs2(upa_K)
            qq(k, j, i, ib_bpa2_KAW) = fac*0.5d0*cabs2(bpa_K)
            qq(k, j, i, ib_upe2_ICW) = fac*0.5d0*cabs2(phi_I)*kp2
            qq(k, j, i, ib_bpe2_ICW) = fac*0.5d0*cabs2(psi_I)*kp2
            qq(k, j, i, ib_upa2_ICW) = fac*0.5d0*cabs2(upa_I)
            qq(k, j, i, ib_bpa2_ICW) = fac*0.5d0*cabs2(bpa_I)

            ! midpoint: dissipation
            phi_m = 0.5d0*(phi(k, j, i) + phi_old(k, j, i))
            psi_m = 0.5d0*(psi(k, j, i) + psi_old(k, j, i))
            upa_m = 0.5d0*(upa(k, j, i) + upa_old(k, j, i))
            bpa_m = 0.5d0*(bpa(k, j, i) + bpa_old(k, j, i))
            a = (phi_m + rho*bpa_m)/(1.d0 + rho**2*kp2)
            c = (rho*kp2*psi_m + sgm*upa_m)/(sgm**2 + rho**2*kp2)
            phi_K = a    ; bpa_K = rho*kp2*a
            psi_K = rho*c; upa_K = sgm*c
            phi_I = phi_m - phi_K; psi_I = psi_m - psi_K
            upa_I = upa_m - upa_K; bpa_I = bpa_m - bpa_K

            dKx_upe = rx_upe*cabs2(phi_K)*kp2; dKz_upe = rz_upe*cabs2(phi_K)*kp2
            dKx_bpe = rx_bpe*cabs2(psi_K)*kp2; dKz_bpe = rz_bpe*cabs2(psi_K)*kp2
            dKx_upa = rx_upa*cabs2(upa_K)    ; dKz_upa = rz_upa*cabs2(upa_K)
            dKx_bpa = rx_bpa*cabs2(bpa_K)    ; dKz_bpa = rz_bpa*cabs2(bpa_K)
            dIx_upe = rx_upe*cabs2(phi_I)*kp2; dIz_upe = rz_upe*cabs2(phi_I)*kp2
            dIx_bpe = rx_bpe*cabs2(psi_I)*kp2; dIz_bpe = rz_bpe*cabs2(psi_I)*kp2
            dIx_upa = rx_upa*cabs2(upa_I)    ; dIz_upa = rz_upa*cabs2(upa_I)
            dIx_bpa = rx_bpa*cabs2(bpa_I)    ; dIz_bpa = rz_bpa*cabs2(bpa_I)

            upe2_KAW_dissip_sum = upe2_KAW_dissip_sum + fac*(dKx_upe + dKz_upe)
            bpe2_KAW_dissip_sum = bpe2_KAW_dissip_sum + fac*(dKx_bpe + dKz_bpe)
            upa2_KAW_dissip_sum = upa2_KAW_dissip_sum + fac*(dKx_upa + dKz_upa)
            bpa2_KAW_dissip_sum = bpa2_KAW_dissip_sum + fac*(dKx_bpa + dKz_bpa)
            upe2_ICW_dissip_sum = upe2_ICW_dissip_sum + fac*(dIx_upe + dIz_upe)
            bpe2_ICW_dissip_sum = bpe2_ICW_dissip_sum + fac*(dIx_bpe + dIz_bpe)
            upa2_ICW_dissip_sum = upa2_ICW_dissip_sum + fac*(dIx_upa + dIz_upa)
            bpa2_ICW_dissip_sum = bpa2_ICW_dissip_sum + fac*(dIx_bpa + dIz_bpa)
            KAW_dissip_prp_sum  = KAW_dissip_prp_sum  + fac*(dKx_upe + dKx_bpe + dKx_upa + dKx_bpa)
            KAW_dissip_par_sum  = KAW_dissip_par_sum  + fac*(dKz_upe + dKz_bpe + dKz_upa + dKz_bpa)
            ICW_dissip_prp_sum  = ICW_dissip_prp_sum  + fac*(dIx_upe + dIx_bpe + dIx_upa + dIx_bpa)
            ICW_dissip_par_sum  = ICW_dissip_par_sum  + fac*(dIz_upe + dIz_bpe + dIz_upa + dIz_bpa)

            qq(k, j, i, ib_dissip_KAW) = fac*( dKx_upe + dKz_upe + dKx_bpe + dKz_bpe &
                                             + dKx_upa + dKz_upa + dKx_bpa + dKz_bpa )
            qq(k, j, i, ib_dissip_ICW) = fac*( dIx_upe + dIz_upe + dIx_bpe + dIz_bpe &
                                             + dIx_upa + dIz_upa + dIx_bpa + dIz_bpa )
          end do
        end do
      end do
      !$acc end parallel loop
    endif
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    series(is_upe2           ) = upe2_sum
    series(is_bpe2           ) = bpe2_sum
    series(is_upa2           ) = upa2_sum
    series(is_bpa2           ) = bpa2_sum
    series(is_upe2dot        ) = upe2dot_sum
    series(is_bpe2dot        ) = bpe2dot_sum
    series(is_upa2dot        ) = upa2dot_sum
    series(is_bpa2dot        ) = bpa2dot_sum
    series(is_upe2dissip     ) = upe2dissip_sum
    series(is_bpe2dissip     ) = bpe2dissip_sum
    series(is_upa2dissip     ) = upa2dissip_sum
    series(is_bpa2dissip     ) = bpa2dissip_sum
    series(is_p_aw           ) = p_aw_sum
    series(is_p_compr        ) = p_compr_sum
    series(is_p_xhl          ) = p_xhl_sum
    series(is_zppe2          ) = zppe2_sum
    series(is_zmpe2          ) = zmpe2_sum
    series(is_zppa2          ) = zppa2_sum
    series(is_zmpa2          ) = zmpa2_sum
    series(is_upe2_KAW_dissip) = upe2_KAW_dissip_sum
    series(is_bpe2_KAW_dissip) = bpe2_KAW_dissip_sum
    series(is_upa2_KAW_dissip) = upa2_KAW_dissip_sum
    series(is_bpa2_KAW_dissip) = bpa2_KAW_dissip_sum
    series(is_upe2_ICW_dissip) = upe2_ICW_dissip_sum
    series(is_bpe2_ICW_dissip) = bpe2_ICW_dissip_sum
    series(is_upa2_ICW_dissip) = upa2_ICW_dissip_sum
    series(is_bpa2_ICW_dissip) = bpa2_ICW_dissip_sum
    series(is_KAW_dissip_prp ) = KAW_dissip_prp_sum
    series(is_KAW_dissip_par ) = KAW_dissip_par_sum
    series(is_ICW_dissip_prp ) = ICW_dissip_prp_sum
    series(is_ICW_dissip_par ) = ICW_dissip_par_sum
    series(is_xhl            ) = xhl_sum

    ! Spectral partial sums: reduce on comm_fft; the global proc0 is the
    ! root of comm_fft group 0, so the result lands where loop_io writes.
    call sum_reduce(series, 0, comm=comm_fft)

    ! Bin over kprp (each call reduces on comm_fft; the host copy of ebin
    ! holds the reduced bins on proc0)
    do l = 1, nq
      call get_polar_spectrum_2d(qq(:,:,:,l), ebin)
      bins(:,:,l) = ebin
    enddo

    !$acc exit data delete(qq, ebin)
    deallocate(qq, ebin)

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io(series, nkpolar, bins)

    deallocate(bins)

    ! With P_m>1 (and/or P_s>1) every comm_fft group redundantly solves the
    ! same problem; verify the fields stay bitwise identical across groups.
    if (nproc_m*nproc_s > 1) call check_redundant_consistency
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @brief   Fluxes of W and X through kprp,
!!          from the nonlinear terms of the
!!          midpoint fields
!!
!!          tw, tx: rates of change of the W and X
!!          densities under the nonlinear terms;
!!          flux(kpbin(idx)) = - sum of the transfer
!!          over the bins below kpbin(idx).
!!          Zero if not nonlinear; flux_X is zero for
!!          sgm = 0 (X is undefined).
!!          flux_W/flux_X are set on proc0 only.
!-----------------------------------------------!
  subroutine get_flux_WX(flux_W, flux_X)
    use params, only: nonlinear, rho, sgm
    use mp, only: proc0
    use grid, only: kprp2
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx_local, nly, nlz_padded
    use fields, only: phi, psi, upa, bpa
    use fields, only: phi_old, psi_old, upa_old, bpa_old
    use fields, only: nfields, iomg, ipsi, iupa, ibpa
    use model_specific, only: nbtran, gradients_to_real, nonlinear_to_spectral
    implicit none
    real(8), dimension(1:nkpolar, nkz), intent(out) :: flux_W, flux_X
    integer :: i, j, k, idx
    real(8) :: fac, kp2
    logical :: x_ok

    complex(8), allocatable, dimension(:,:,:)   :: phi_m, psi_m, upa_m, bpa_m
    complex(8), allocatable, dimension(:,:,:,:) :: w, nonlin
    real(8)   , allocatable, dimension(:,:,:,:) :: w_r, nonlin_r
    real(8)   , allocatable, dimension(:,:,:)   :: tw, tx
    real(8)   , allocatable, dimension(:,:)     :: t_bin

    flux_W = 0.d0
    flux_X = 0.d0
    if (.not. nonlinear) return

    x_ok = sgm > 0.d0

    allocate(phi_m   (nkz       , nky_local, nkx      ))
    allocate(psi_m   (nkz       , nky_local, nkx      ))
    allocate(upa_m   (nkz       , nky_local, nkx      ))
    allocate(bpa_m   (nkz       , nky_local, nkx      ))
    allocate(w       (nkz       , nky_local, nkx      , nbtran ))
    allocate(nonlin  (nkz       , nky_local, nkx      , nfields))
    allocate(w_r     (nlz_padded, nly      , nlx_local, nbtran ))
    allocate(nonlin_r(nlz_padded, nly      , nlx_local, nfields))
    !$acc enter data create(phi_m, psi_m, upa_m, bpa_m, w, nonlin, w_r, nonlin_r)

    !$acc parallel loop collapse(3) default(present)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi_m(k, j, i) = 0.5d0*(phi(k, j, i) + phi_old(k, j, i))
          psi_m(k, j, i) = 0.5d0*(psi(k, j, i) + psi_old(k, j, i))
          upa_m(k, j, i) = 0.5d0*(upa(k, j, i) + upa_old(k, j, i))
          bpa_m(k, j, i) = 0.5d0*(bpa(k, j, i) + bpa_old(k, j, i))
        end do
      end do
    end do
    !$acc end parallel loop

    ! the same nonlinear terms as the time integration (model_specific)
    call gradients_to_real(phi_m, psi_m, upa_m, bpa_m, w, w_r)
    call nonlinear_to_spectral(w_r, nonlin_r, nonlin)

    !$acc exit data delete(w, w_r, nonlin_r)
    deallocate(w, w_r, nonlin_r)

    allocate(tw(nkz, nky_local, nkx))
    allocate(tx(nkz, nky_local, nkx))
    allocate(t_bin(1:nkpolar, nkz))
    !$acc enter data create(tw, tx, t_bin)

    !$acc parallel loop collapse(3) gang vector default(present) private(fac, kp2)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fac = 1.d0
          if (k /= 1) fac = 2.d0
          kp2 = kprp2(k, j, i)

          tw(k, j, i) = fac*( - dble(conjg(phi_m(k, j, i))*nonlin(k, j, i, iomg)) &
                              + kp2*dble(conjg(psi_m(k, j, i))*nonlin(k, j, i, ipsi)) &
                              + dble(conjg(upa_m(k, j, i))*nonlin(k, j, i, iupa)) &
                              + dble(conjg(bpa_m(k, j, i))*nonlin(k, j, i, ibpa)) )

          if (x_ok) then
            tx(k, j, i) = fac*( - dble(conjg(nonlin(k, j, i, iomg))*psi_m(k, j, i)) &
                                + kp2*dble(conjg(phi_m(k, j, i))*nonlin(k, j, i, ipsi)) &
                                + ( dble(conjg(nonlin(k, j, i, iupa))*bpa_m(k, j, i)) &
                                  + dble(conjg(upa_m(k, j, i))*nonlin(k, j, i, ibpa)) )/sgm &
                                + rho/sgm*( dble(conjg(upa_m(k, j, i))*nonlin(k, j, i, iomg)) &
                                          - kp2*dble(conjg(nonlin(k, j, i, iupa))*phi_m(k, j, i)) ) )
          else
            tx(k, j, i) = 0.d0
          endif
        end do
      end do
    end do
    !$acc end parallel loop

    call get_polar_spectrum_2d(tw, t_bin)
    if (proc0) then
      flux_W(1, :) = 0.d0
      do idx = 2, nkpolar
        flux_W(idx, :) = flux_W(idx-1, :) - t_bin(idx-1, :)
      enddo
    endif

    if (x_ok) then
      call get_polar_spectrum_2d(tx, t_bin)
      if (proc0) then
        flux_X(1, :) = 0.d0
        do idx = 2, nkpolar
          flux_X(idx, :) = flux_X(idx-1, :) - t_bin(idx-1, :)
        enddo
      endif
    endif

    !$acc exit data delete(phi_m, psi_m, upa_m, bpa_m, nonlin, tw, tx, t_bin)
    deallocate(phi_m, psi_m, upa_m, bpa_m, nonlin, tw, tx, t_bin)
  end subroutine get_flux_WX


!-----------------------------------------------!
!> @author  YK
!! @brief   Runtime assertion that the redundantly-solved fields
!!          (phi, psi, upa, bpa) are bitwise identical across all
!!          comm_m/comm_s groups. A fixed loop order makes the scalar
!!          checksum reproducible, so max==min over a redundant axis
!!          holds iff the field data agree. Active only when
!!          nproc_m*nproc_s > 1.
!-----------------------------------------------!
  subroutine check_redundant_consistency
    use MPI
    use mp, only: proc0
    use mp, only: comm_m, comm_s, nproc_m, nproc_s
    use mp, only: max_allreduce, min_allreduce, sum_allreduce
    use grid, only: nkx, nky_local, nkz
    use fields, only: phi, psi, upa, bpa
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    integer :: i, j, k, nbad, ierr
    real(8) :: chk, chk_max, chk_min
    integer(int64) :: h, h_or, h_and

    ! (1) EXACT bitwise-identity test. XOR-fold the raw IEEE bit patterns of every
    ! field component in a fixed host-side order. XOR is associative AND
    ! commutative, so -- unlike a floating-point sum -- the fold is immune to the
    ! GPU reduction ORDER, which is not guaranteed identical across the distinct
    ! GPUs of a redundant comm_m group. The fold changes iff the field data differ
    ! bit-for-bit, so it tests true redundant-solve identity rather than
    ! reduction-order identity (a plain FP sum yields false ~1 ULP mismatches once
    ! the field grows enough for the summation order to matter).
    !$acc update host(phi, psi, upa, bpa)
    h = 0_int64
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          call fold(h, phi(k,j,i))
          call fold(h, psi(k,j,i))
          call fold(h, upa(k,j,i))
          call fold(h, bpa(k,j,i))
        end do
      end do
    end do

    ! (2) floating-point checksum, kept ONLY to report the magnitude of any
    ! disagreement. This sum is order-sensitive and must never drive PASS/FAIL.
    chk = 0.d0
    !$acc data present(phi, psi, upa, bpa)
    !$acc parallel loop collapse(3) reduction(+:chk)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          chk = chk + dble(phi(k,j,i)) + aimag(phi(k,j,i)) &
                    + dble(psi(k,j,i)) + aimag(psi(k,j,i)) &
                    + dble(upa(k,j,i)) + aimag(upa(k,j,i)) &
                    + dble(bpa(k,j,i)) + aimag(bpa(k,j,i))
        end do
      end do
    end do
    !$acc end data

    ! A group is bitwise identical iff every rank's hash agrees, i.e. the bitwise
    ! OR across the group equals the bitwise AND across the group.
    nbad = 0
    chk_max = chk; chk_min = chk
    if (nproc_m > 1) then
      call mpi_allreduce(h, h_or,  1, MPI_INTEGER8, MPI_BOR,  comm_m, ierr)
      call mpi_allreduce(h, h_and, 1, MPI_INTEGER8, MPI_BAND, comm_m, ierr)
      if (h_or /= h_and) nbad = nbad + 1
      call max_allreduce(chk_max, comm=comm_m)
      call min_allreduce(chk_min, comm=comm_m)
    endif
    if (nproc_s > 1) then
      call mpi_allreduce(h, h_or,  1, MPI_INTEGER8, MPI_BOR,  comm_s, ierr)
      call mpi_allreduce(h, h_and, 1, MPI_INTEGER8, MPI_BAND, comm_s, ierr)
      if (h_or /= h_and) nbad = nbad + 1
      call max_allreduce(chk_max, comm=comm_s)
      call min_allreduce(chk_min, comm=comm_s)
    endif

    call sum_allreduce(nbad)
    if (proc0) then
      if (nbad == 0) then
        write (*, '("[check_redundant_consistency] PASS: comm_m/comm_s groups bitwise identical")')
      else
        write (*, '("[check_redundant_consistency] FAIL: ", i0, &
          &" group mismatch(es); FP chk span=", es12.4e3, " rel=", es10.2e3)') &
          nbad, chk_max - chk_min, (chk_max - chk_min)/max(abs(chk_max), 1.d-300)
      endif
    endif

  contains

    subroutine fold(h, z)
      integer(int64), intent(inout) :: h
      complex(8), intent(in) :: z
      h = ieor(h, transfer(dble (z), 0_int64))
      h = ieor(h, transfer(aimag(z), 0_int64))
    end subroutine fold

  end subroutine check_redundant_consistency


!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics for cross section of fileds
!!          (fields in the order of io::i2d_*)
!-----------------------------------------------!
  subroutine loop_diagnostics_2D
    use io, only: loop_io_2D, n2d
    use io, only: i2d_phi, i2d_psi, i2d_upa, i2d_bpa, i2d_omg, i2d_jpa
    use io, only: i2d_ux, i2d_uy, i2d_bx, i2d_by
    use mp, only: proc0
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: kx, ky, kprp2
    use grid, only: nkx, nky_local, nkz
    use fields, only: phi, psi, upa, bpa
    use cuFFTmp, only: btran_c2r
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi
    implicit none
    integer :: i, j, k, l

    complex(8), allocatable, dimension(:,:,:) :: f
    real(8)   , allocatable, dimension(:,:,:) :: fr
    real(8)   , allocatable, dimension(:,:,:) :: cut_z0, cut_x0, cut_y0

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(f (nkz       , nky_local, nkx      ), source=(0.d0, 0.d0))
    allocate(fr(nlz_padded, nly      , nlx_local), source=0.d0)
    !$acc enter data create(f, fr)

    allocate(cut_z0(nlx_local, nly, n2d), source=0.d0)
    allocate(cut_x0(nly      , nlz, n2d), source=0.d0)
    allocate(cut_y0(nlx_local, nlz, n2d), source=0.d0)

    do l = 1, n2d
      !$acc parallel loop collapse(3) present(phi, psi, upa, bpa, f, kx, ky, kprp2)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            if      (l == i2d_phi) then
              f(k, j, i) = phi(k, j, i)
            else if (l == i2d_psi) then
              f(k, j, i) = psi(k, j, i)
            else if (l == i2d_upa) then
              f(k, j, i) = upa(k, j, i)
            else if (l == i2d_bpa) then
              f(k, j, i) = bpa(k, j, i)
            else if (l == i2d_omg) then
              f(k, j, i) = -kprp2(k, j, i)*phi(k, j, i)
            else if (l == i2d_jpa) then
              f(k, j, i) = -kprp2(k, j, i)*psi(k, j, i)
            else if (l == i2d_ux ) then
              f(k, j, i) = -zi*ky(j)*phi(k, j, i)
            else if (l == i2d_uy ) then
              f(k, j, i) =  zi*kx(i)*phi(k, j, i)
            else if (l == i2d_bx ) then
              f(k, j, i) = -zi*ky(j)*psi(k, j, i)
            else if (l == i2d_by ) then
              f(k, j, i) =  zi*kx(i)*psi(k, j, i)
            endif
          enddo
        enddo
      enddo
      !$acc end parallel loop

      call btran_c2r(f, fr)
      !$acc update host(fr)
      call cut_2d_r(fr, cut_z0(:,:,l), cut_x0(:,:,l), cut_y0(:,:,l))
    enddo

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io_2D(cut_z0, cut_x0, cut_y0)

    !$acc exit data delete(f, fr)
    deallocate(f, fr)
    deallocate(cut_z0, cut_x0, cut_y0)
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
