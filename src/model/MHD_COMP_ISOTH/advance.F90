!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Time stepping for MHD_COMP_ISOTH
!!
!!          State (Fourier): s = ln(rho), m = rho*u, B
!!            ds/dt = -(div_N m)/rho
!!            dm/dt = -div_N [m u + (cs^2 rho + B^2/2 - q) I - B B]
!!            dB/dt = -curl_N (B x u)
!!          q = rho*zeta*D, D = div_N u (see shock_viscosity).
!!          div_N = F^-1 i k_d M F with k_d = 0 at Nyquist and M the
!!          filter of the differentiated field.
!!
!!          rho = exp(s) is positive by construction; no floor is
!!          applied. A non-finite value or |s| > s_abort stops the run.
!!          Only eSSPIFRK3 is implemented. Its effective weights are
!!          b = (1/4, 3/16, 9/16) at c = (0, 2/3, 2/3), which are used
!!          for the energy budget.
!!
!!          Gallope (OpenACC + cuFFTMp): the layout is f(k,j,i) with
!!          (nlz_padded, nly, nlx_local) in real space and
!!          (nkz, nky_local, nkx) in Fourier space. Every work array
!!          is allocated once (allocate_advance) and stays on the
!!          device; only scalars reach the host inside a step. Real-
!!          space loops run over k = 1, nlz (the FFT padding rows
!!          k > nlz are not physical points).
!-----------------------------------------------!
module advance
  use fields, only: nfields
  use fields, only: ilnrho
  use fields, only: imx, imy, imz
  use fields, only: ibx, iby, ibz
  implicit none

  public solve, is_allocated, allocate_advance, deallocate_advance

  logical :: is_allocated = .false.

  integer :: counter = 0
  integer :: istep_adv = 0 ! number of steps completed in this run
  complex(8), allocatable, dimension(:,:,:)   :: lnrho_new
  complex(8), allocatable, dimension(:,:,:)   ::    mx_new
  complex(8), allocatable, dimension(:,:,:)   ::    my_new
  complex(8), allocatable, dimension(:,:,:)   ::    mz_new
  complex(8), allocatable, dimension(:,:,:)   ::    bx_new
  complex(8), allocatable, dimension(:,:,:)   ::    by_new
  complex(8), allocatable, dimension(:,:,:)   ::    bz_new
  complex(8), allocatable, dimension(:,:,:,:) :: w
  real   (8), allocatable, dimension(:,:,:,:) :: w_r, flx_r
  real   (8), allocatable, dimension(:,:,:)   :: rho_r, ux_r, uy_r, uz_r
  real   (8), allocatable, dimension(:,:,:)   :: div_r, zeta_r, q_r ! D, zeta, q = rho*zeta*D
  real   (8), allocatable, dimension(:,:,:)   :: eta_r ! shock resistivity eta (C_eta > 0)
  real   (8), allocatable, dimension(:,:,:,:) :: ej_r  ! eta J, J = curl_N B (C_eta > 0)
  real   (8), allocatable, dimension(:,:,:)   :: svv_r ! SVV term of the s equation, (d rho/dt)_svv/rho
  complex(8), allocatable, dimension(:,:,:,:) :: svv_mh ! SVV term of m ('entropy'), filtered and normalized as the state
  complex(8), allocatable, dimension(:,:,:,:) :: svv_uh ! F u/(nlx nly nlz) of the current stage ('entropy')
  complex(8), allocatable, dimension(:,:,:,:) :: flx, exp_terms
  real   (8), allocatable, dimension(:,:)     :: kxt ! sheared kx
  real   (8), allocatable, dimension(:,:)     :: kxd ! kx for derivatives: kx_d (shear = F) or kxt (shear = T)

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(8), allocatable, dimension(:,:,:)   :: lnrho_tmp
  complex(8), allocatable, dimension(:,:,:)   ::    mx_tmp
  complex(8), allocatable, dimension(:,:,:)   ::    my_tmp
  complex(8), allocatable, dimension(:,:,:)   ::    mz_tmp
  complex(8), allocatable, dimension(:,:,:)   ::    bx_tmp
  complex(8), allocatable, dimension(:,:,:)   ::    by_tmp
  complex(8), allocatable, dimension(:,:,:)   ::    bz_tmp
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms0
  real(8), parameter :: b_rk(3) = (/0.25d0, 0.1875d0, 0.5625d0/) ! effective weights of eSSPIFRK3

  ! per-field settings: filter after each stage, linear dissipation (lmd for s, nu for m, eta for B)
  logical :: filt_on  (nfields)
  real(8) :: dcoef    (nfields), dcoef_h(nfields)
  integer :: dcoef_exp(nfields)

  integer :: cfl_unit, cons_unit, budget_unit
  integer :: nline1d = 0
  integer :: nplanes = 0 ! number of records in out2d/plane_*.dat

  ! Forward FFT variables
  integer, parameter :: nftran = 10
  integer, parameter :: iflx_s   = 1                               ! -(div_N m)/rho
  integer, parameter :: iflx_mxx = 2, iflx_mxy = 3, iflx_mxz = 4   ! m u - B B + (cs^2 rho + B^2/2 - q) I (tensor)
  integer, parameter ::               iflx_myy = 5, iflx_myz = 6   ! m u - B B + (cs^2 rho + B^2/2 - q) I (tensor)
  integer, parameter ::                             iflx_mzz = 7   ! m u - B B + (cs^2 rho + B^2/2 - q) I (tensor)
  integer, parameter :: iflx_bx  = 8, iflx_by  = 9, iflx_bz  = 10  ! B x u

  ! Volume integrals of the state (conservation.dat)
  integer, parameter :: nmon = 7
  integer, parameter :: imon_M  = 1                         ! mass
  integer, parameter :: imon_px = 2, imon_py = 3, imon_pz = 4 ! momentum
  integer, parameter :: imon_ek = 5, imon_em = 6, imon_ei = 7 ! kinetic, magnetic, internal energy
  real(8) :: mon_n(nmon)          ! at the beginning of the current step (after the mass correction)
  real(8) :: mass0 = 0.d0         ! reference mass = mass at the first step of this run
  real(8) :: dmass_corr     = 0.d0 ! ln(M0/M) applied at the current step
  real(8) :: dmass_corr_tot = 0.d0

  ! Of the latest RHS evaluation (volume integrals)
  real(8) :: dMdt_raw  ! -sum div_N m (zero in exact arithmetic)
  real(8) :: dMdt_abs  !  sum |div_N m| (scale of dMdt_raw)
  real(8) :: visc_diss !  sum rho zeta D^2 (>= 0)
  real(8) :: eta_diss = 0.d0 !  sum eta |J|^2 (>= 0)

  ! SVV coefficients eps_d of the current step (set at stage 1 and fixed over the step)
  real(8) :: eps_svv(3) = 0.d0
  ! 'entropy': D_svv = sum [a t^2 + b |g_m|^2 + kap rho g_h^2] dV (>= 0) of the latest RHS evaluation
  real(8) :: svv_diss = 0.d0
  ! dt before reset_dt of the current step (for the check of dt*lambda_svv)
  real(8) :: dt_prev = 0.d0

  ! Energy and mass budget accumulated over the RK stages with weights b_rk (budget.dat)
  integer, parameter :: nbgt = 19
  integer, parameter :: ib_dMdt_raw = 1, ib_dMdt_proj = 2, ib_dMdt_abs = 3
  integer, parameter :: ib_eps_N    = 4, ib_eps_shock = 5, ib_visc     = 6, ib_R_zeta = 7
  integer, parameter :: ib_eps_force = 8, ib_eps_res  = 9, ib_eps_nu   = 10, ib_eps_lmd = 11
  integer, parameter :: ib_eps_svv  = 12, ib_eps_svv_b = 13 ! SVV of (rho, m) and of B
  integer, parameter :: ib_D_svv    = 14, ib_R_svv     = 15 ! D_svv and R_svv = eps_svv + D_svv ('entropy')
  integer, parameter :: ib_A_svv    = 16 ! sum |dE/dU.U_svv| dV over (rho, m): rounding scale of eps_svv
  integer, parameter :: ib_eps_eta  = 17, ib_eta_diss = 18, ib_R_eta = 19 ! shock resistivity: R_eta = eps_eta + eta_diss
  real(8) :: bgt(nbgt) = 0.d0
  real(8) :: R_zeta_ratio_max = 0.d0
  real(8) :: R_svv_ratio_max  = 0.d0
  real(8) :: R_svv_scale_max  = 0.d0
  real(8) :: R_eta_ratio_max  = 0.d0
  real(8) :: t_bgt, dt_bgt
  logical :: bgt_pending = .false.

  ! Work arrays of stage_budget (budget = T), allocated once in allocate_advance
  ! (Calliope allocates them at every call)
  complex(8), allocatable, dimension(:,:,:)   :: wk, qh
  complex(8), allocatable, dimension(:,:,:,:) :: ejh
  real   (8), allocatable, dimension(:,:,:)   :: wk_r, dEds_r

  ! Device buffers of write_line (lbuf(nline, ncol)) and write_planes (x-z and x-y planes,
  ! pbuf_y(nlx_local, nlz, nq_plane) and pbuf_z(nlx_local, nly, nq_plane))
  integer, parameter :: nq_plane = 10 ! arrays per plane record (nq of write_planes)
  real   (8), allocatable, dimension(:,:)     :: lbuf
  real   (8), allocatable, dimension(:,:,:)   :: pbuf_y, pbuf_z

  ! Device buffer of output_series_modes, smodes(n_series_modes, nfields)
  complex(8), allocatable, dimension(:,:)     :: smodes

  !$acc declare create(kxt, kxd)
  !$acc declare create(filt_on, dcoef, dcoef_h, dcoef_exp)
  !$acc declare create(eps_svv)

contains


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Solve the time evolution (eSSPIFRK3)
!-----------------------------------------------!
  subroutine solve
    use fields, only: lnrho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use fields, only: lnrho_old
    use fields, only: mx_old, my_old, mz_old
    use fields, only: bx_old, by_old, bz_old
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use time, only: dt, tt, nstep, t_end_reached
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: shear, budget
    use shearing_box, only: k2t, k2t_inv, tsc, tremap
    use force, only: fmx, fmz, fmy, fmx_old, fmy_old, fmz_old, driven, update_force, get_force
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use diagnostics_common, only: series_output
    use shock_viscosity, only: tmr_on, tmr_off, itmr_step, itmr_budget
    implicit none
    real(8) :: imp_terms_tintg0(nfields), imp_terms_tintg1(nfields)
    real(8) :: imp_terms_tintg2(nfields), imp_terms_tintg3(nfields)
    real(8) :: t0, t2, t3 ! tsc, tsc + 2/3 dt and tsc + dt (evaluated on the host)
    integer :: i, j, k, n

    if (proc0) call put_time_stamp(timer_advance)
    call tmr_on(itmr_step)

    ! initialize tmp fields
    if(counter == 0) then
      call init_work_fields
      counter = 1
    endif

    ! Calcualte force terms
    if (driven) then
      ! at n
      call get_force('mx', fmx)
      call get_force('my', fmy)
      call get_force('mz', fmz)
    endif

    !---------------  RK 1st step  ---------------
    ! Calcualte nonlinear terms (also CFL, monitor, and mass correction of the real-space s)
    call get_nonlinear_terms(lnrho, mx, my, mz, bx, by, bz, 1)
    if(dmass_corr /= 0.d0) call shift_mean(lnrho, dmass_corr)

    call get_explicit_terms(mx, my, mz, bx, by, bz)
    if(budget) then
      call tmr_on(itmr_budget)
      call stage_budget(1, lnrho, mx, my, mz, bx, by, bz)
      call tmr_off(itmr_budget)
    endif

    ! dt is final after the stage-1 RHS (reset_dt, dt_adjust)
    t0 = tsc
    t2 = tsc + 2.d0/3.d0*dt

    !$acc data present(lnrho, mx, my, mz, bx, by, bz, &
    !$acc              lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, &
    !$acc              exp_terms, kx, ky, kz, dcoef, dcoef_h, dcoef_exp)
    !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg1)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ! Calculate time integral of implicit terms
          !$acc loop seq
          do n = 1, nfields
            call get_imp_terms_tintg(imp_terms_tintg0(n), t0, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
            call get_imp_terms_tintg(imp_terms_tintg1(n), t2, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
          enddo

          ! update s
          call eSSPIFRK1(lnrho_tmp(k,j,i), lnrho(k,j,i), &
             exp_terms(k,j,i,ilnrho), &
             imp_terms_tintg1(ilnrho), imp_terms_tintg0(ilnrho) &
          )

          ! update m
          call eSSPIFRK1(mx_tmp(k,j,i), mx(k,j,i), &
             exp_terms(k,j,i,imx), &
             imp_terms_tintg1(imx), imp_terms_tintg0(imx) &
          )
          call eSSPIFRK1(my_tmp(k,j,i), my(k,j,i), &
             exp_terms(k,j,i,imy), &
             imp_terms_tintg1(imy), imp_terms_tintg0(imy) &
          )
          call eSSPIFRK1(mz_tmp(k,j,i), mz(k,j,i), &
             exp_terms(k,j,i,imz), &
             imp_terms_tintg1(imz), imp_terms_tintg0(imz) &
          )

          ! update b
          call eSSPIFRK1(bx_tmp(k,j,i), bx(k,j,i), &
             exp_terms(k,j,i,ibx), &
             imp_terms_tintg1(ibx), imp_terms_tintg0(ibx) &
          )
          call eSSPIFRK1(by_tmp(k,j,i), by(k,j,i), &
             exp_terms(k,j,i,iby), &
             imp_terms_tintg1(iby), imp_terms_tintg0(iby) &
          )
          call eSSPIFRK1(bz_tmp(k,j,i), bz(k,j,i), &
             exp_terms(k,j,i,ibz), &
             imp_terms_tintg1(ibz), imp_terms_tintg0(ibz) &
          )
        enddo
      enddo
    enddo
    !$acc end data

    ! save explicit terms at the previous step
    !$acc parallel loop collapse(4) present(exp_terms, exp_terms0)
    do n = 1, nfields
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            exp_terms0(k,j,i,n) = exp_terms(k,j,i,n)
          enddo
        enddo
      enddo
    enddo

    ! Dealiasing
    call apply_filter(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)

    !---------------  RK 2nd step  ---------------
    ! Calcualte force terms
    if (driven) then
      ! at n + 2/3
      call update_force(2.d0/3.d0*dt)
      call get_force('mx', fmx)
      call get_force('my', fmy)
      call get_force('mz', fmz)

      ! go to n + 1
      call update_force(1.d0/3.d0*dt)
    endif

    ! Calcualte kxt at n + 2/3
    if(shear) call set_kx(tsc + 2.d0/3.d0*dt)

    ! Calcualte nonlinear terms
    call get_nonlinear_terms(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, 2)

    call get_explicit_terms(mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)
    if(budget) then
      call tmr_on(itmr_budget)
      call stage_budget(2, lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)
      call tmr_off(itmr_budget)
    endif

    !$acc data present(lnrho, mx, my, mz, bx, by, bz, &
    !$acc              lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, &
    !$acc              exp_terms, kx, ky, kz, dcoef, dcoef_h, dcoef_exp)
    !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg2)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ! Calculate time integral of implicit terms
          !$acc loop seq
          do n = 1, nfields
            call get_imp_terms_tintg(imp_terms_tintg0(n), t0, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
            call get_imp_terms_tintg(imp_terms_tintg2(n), t2, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
          enddo

          ! update s
          call eSSPIFRK2(lnrho_tmp(k,j,i), lnrho_tmp(k,j,i), lnrho(k,j,i), &
             exp_terms(k,j,i,ilnrho), &
             imp_terms_tintg2(ilnrho), imp_terms_tintg0(ilnrho) &
          )

          ! update m
          call eSSPIFRK2(mx_tmp(k,j,i), mx_tmp(k,j,i), mx(k,j,i), &
             exp_terms(k,j,i,imx), &
             imp_terms_tintg2(imx), imp_terms_tintg0(imx) &
          )
          call eSSPIFRK2(my_tmp(k,j,i), my_tmp(k,j,i), my(k,j,i), &
             exp_terms(k,j,i,imy), &
             imp_terms_tintg2(imy), imp_terms_tintg0(imy) &
          )
          call eSSPIFRK2(mz_tmp(k,j,i), mz_tmp(k,j,i), mz(k,j,i), &
             exp_terms(k,j,i,imz), &
             imp_terms_tintg2(imz), imp_terms_tintg0(imz) &
          )

          ! update b
          call eSSPIFRK2(bx_tmp(k,j,i), bx_tmp(k,j,i), bx(k,j,i), &
             exp_terms(k,j,i,ibx), &
             imp_terms_tintg2(ibx), imp_terms_tintg0(ibx) &
          )
          call eSSPIFRK2(by_tmp(k,j,i), by_tmp(k,j,i), by(k,j,i), &
             exp_terms(k,j,i,iby), &
             imp_terms_tintg2(iby), imp_terms_tintg0(iby) &
          )
          call eSSPIFRK2(bz_tmp(k,j,i), bz_tmp(k,j,i), bz(k,j,i), &
             exp_terms(k,j,i,ibz), &
             imp_terms_tintg2(ibz), imp_terms_tintg0(ibz) &
          )
        enddo
      enddo
    enddo
    !$acc end data

    ! Dealiasing
    call apply_filter(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)

    !---------------  RK 3rd step  ---------------
    ! Calcualte nonlinear terms
    call get_nonlinear_terms(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, 3)

    call get_explicit_terms(mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)
    if(budget) then
      call tmr_on(itmr_budget)
      call stage_budget(3, lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)
      call tmr_off(itmr_budget)
    endif

    t3 = tsc + dt

    !$acc data present(lnrho, mx, my, mz, bx, by, bz, &
    !$acc              lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, &
    !$acc              lnrho_new, mx_new, my_new, mz_new, bx_new, by_new, bz_new, &
    !$acc              exp_terms, exp_terms0, kx, ky, kz, dcoef, dcoef_h, dcoef_exp)
    !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg2, imp_terms_tintg3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ! Calculate time integral of implicit terms
          !$acc loop seq
          do n = 1, nfields
            call get_imp_terms_tintg(imp_terms_tintg0(n), t0, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
            call get_imp_terms_tintg(imp_terms_tintg2(n), t2, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
            call get_imp_terms_tintg(imp_terms_tintg3(n), t3, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
          enddo

          ! update s
          call eSSPIFRK3(lnrho_new(k,j,i), lnrho_tmp(k,j,i), lnrho(k,j,i), &
             exp_terms(k,j,i,ilnrho), exp_terms0(k,j,i,ilnrho), &
             imp_terms_tintg3(ilnrho), imp_terms_tintg2(ilnrho), imp_terms_tintg0(ilnrho) &
          )

          ! update m
          call eSSPIFRK3(mx_new(k,j,i), mx_tmp(k,j,i), mx(k,j,i), &
             exp_terms(k,j,i,imx), exp_terms0(k,j,i,imx), &
             imp_terms_tintg3(imx), imp_terms_tintg2(imx), imp_terms_tintg0(imx) &
          )
          call eSSPIFRK3(my_new(k,j,i), my_tmp(k,j,i), my(k,j,i), &
             exp_terms(k,j,i,imy), exp_terms0(k,j,i,imy), &
             imp_terms_tintg3(imy), imp_terms_tintg2(imy), imp_terms_tintg0(imy) &
          )
          call eSSPIFRK3(mz_new(k,j,i), mz_tmp(k,j,i), mz(k,j,i), &
             exp_terms(k,j,i,imz), exp_terms0(k,j,i,imz), &
             imp_terms_tintg3(imz), imp_terms_tintg2(imz), imp_terms_tintg0(imz) &
          )

          ! update b
          call eSSPIFRK3(bx_new(k,j,i), bx_tmp(k,j,i), bx(k,j,i), &
             exp_terms(k,j,i,ibx), exp_terms0(k,j,i,ibx), &
             imp_terms_tintg3(ibx), imp_terms_tintg2(ibx), imp_terms_tintg0(ibx) &
          )
          call eSSPIFRK3(by_new(k,j,i), by_tmp(k,j,i), by(k,j,i), &
             exp_terms(k,j,i,iby), exp_terms0(k,j,i,iby), &
             imp_terms_tintg3(iby), imp_terms_tintg2(iby), imp_terms_tintg0(iby) &
          )
          call eSSPIFRK3(bz_new(k,j,i), bz_tmp(k,j,i), bz(k,j,i), &
             exp_terms(k,j,i,ibz), exp_terms0(k,j,i,ibz), &
             imp_terms_tintg3(ibz), imp_terms_tintg2(ibz), imp_terms_tintg0(ibz) &
          )
        enddo
      enddo
    enddo
    !$acc end data

    ! save fields at the previous step
    !$acc parallel loop collapse(3) present(lnrho, mx, my, mz, bx, by, bz, &
    !$acc   lnrho_old, mx_old, my_old, mz_old, bx_old, by_old, bz_old)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          lnrho_old(k,j,i) = lnrho(k,j,i)
             mx_old(k,j,i) =    mx(k,j,i)
             my_old(k,j,i) =    my(k,j,i)
             mz_old(k,j,i) =    mz(k,j,i)
             bx_old(k,j,i) =    bx(k,j,i)
             by_old(k,j,i) =    by(k,j,i)
             bz_old(k,j,i) =    bz(k,j,i)
        enddo
      enddo
    enddo

    if (driven) then
      !$acc parallel loop collapse(3) present(fmx, fmy, fmz, fmx_old, fmy_old, fmz_old)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            fmx_old(k,j,i) = fmx(k,j,i)
            fmy_old(k,j,i) = fmy(k,j,i)
            fmz_old(k,j,i) = fmz(k,j,i)
          enddo
        enddo
      enddo
    endif

    ! Dealiasing
    call apply_filter(lnrho_new, mx_new, my_new, mz_new, bx_new, by_new, bz_new)

    !$acc parallel loop collapse(3) present(lnrho, mx, my, mz, bx, by, bz, &
    !$acc   lnrho_new, mx_new, my_new, mz_new, bx_new, by_new, bz_new)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          lnrho(k,j,i) = lnrho_new(k,j,i)
             mx(k,j,i) =    mx_new(k,j,i)
             my(k,j,i) =    my_new(k,j,i)
             mz(k,j,i) =    mz_new(k,j,i)
             bx(k,j,i) =    bx_new(k,j,i)
             by(k,j,i) =    by_new(k,j,i)
             bz(k,j,i) =    bz_new(k,j,i)
        enddo
      enddo
    enddo

    tt  = tt  + dt
    tsc = tsc + dt
    !$acc update device(tt, tsc)

    if(shear .and. tsc > tremap) call remap

    if(shear) then
      call set_kx(tsc)
      !$acc parallel loop collapse(3) present(k2t, k2t_inv, kxt, ky, kz)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            k2t(k,j,i) = kxt(i,j)**2 + ky(j)**2 + kz(k)**2

            if(k2t(k,j,i) == 0.d0) then
              k2t_inv(k,j,i) = 0.d0
            else
              k2t_inv(k,j,i) = 1.d0/k2t(k,j,i)
            endif
          enddo
        enddo
      enddo
    endif

    ! Div b cleaing
    call div_free(bx, by, bz)

    if(series_output) call output_series_modes

    istep_adv = istep_adv + 1

    ! Monitor and line output of the final state (the RHS is not evaluated): after nstep steps or at t_end
    if(istep_adv == nstep .or. t_end_reached()) call get_nonlinear_terms(lnrho, mx, my, mz, bx, by, bz, 0)

    call tmr_off(itmr_step)
    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Check the run parameters, allocate
!!          fields used only here, and open the
!!          output files of this module
!-----------------------------------------------!
  subroutine init_work_fields
    use fields, only: lnrho, mx, my, mz, bx, by, bz
    use grid, only: ky, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: nonlinear, shear, time_step_scheme
    use params, only: nu, nu_h, nu_h_exp, eta, eta_h, eta_h_exp, lmd, lmd_h, lmd_h_exp
    use params, only: formulation, visc_form, filter_lnrho, filter_mom, filter_mag
    use params, only: C_shock, s_abort, budget, write_intvl_cons, svv_eps0, svv_form, cfl_svv, C_eta, svv_mass
    use shearing_box, only: tsc, k2t, k2t_inv
    use shock_viscosity, only: init_shock_viscosity, abort_run
    use file, only: open_output_file
    use mp, only: proc0
    implicit none
    integer :: i, j, k

    if(time_step_scheme /= 'eSSPIFRK3') &
      call abort_run('MHD_COMP_ISOTH: time_step_scheme = '//trim(time_step_scheme)//' is not implemented (use eSSPIFRK3)')
    if(trim(formulation) /= 'C') &
      call abort_run('MHD_COMP_ISOTH: formulation = '//trim(formulation)//' is not implemented (use C)')
    if(trim(visc_form) /= 'flux') &
      call abort_run('MHD_COMP_ISOTH: visc_form = '//trim(visc_form)//' is not implemented yet (use flux)')
    if(.not. nonlinear) &
      call abort_run('MHD_COMP_ISOTH: nonlinear = F is not supported (mass flux and pressure are nonlinear terms)')
    if(shear .and. C_shock > 0.d0) &
      call abort_run('MHD_COMP_ISOTH: shear = T with C_shock > 0 has not been verified')
    if(shear .and. C_eta > 0.d0) &
      call abort_run('MHD_COMP_ISOTH: shear = T with C_eta > 0 has not been verified')
    if(shear .and. svv_eps0 > 0.d0) &
      call abort_run('MHD_COMP_ISOTH: shear = T with svv_eps0 > 0 is not implemented (the SVV kernel uses the unsheared kx)')
    if(trim(svv_form) /= 'entropy' .and. trim(svv_form) /= 'linear') &
      call abort_run('MHD_COMP_ISOTH: svv_form = '//trim(svv_form)//' is not implemented (use entropy or linear)')
    if(trim(svv_form) == 'entropy' .and. cfl_svv <= 0.d0) &
      call abort_run('MHD_COMP_ISOTH: cfl_svv must be > 0')
    if(.not. svv_mass .and. trim(svv_form) /= 'entropy') &
      call abort_run('MHD_COMP_ISOTH: svv_mass = F is implemented for svv_form = entropy only')
    ! Calliope's check of pruned = T is dropped: Gallope has no pruned FFT
    if(s_abort > 7.d2) &
      call abort_run('MHD_COMP_ISOTH: s_abort must be <= 700 so that exp(s) does not overflow')
    if(proc0 .and. (lmd /= 0.d0 .or. lmd_h /= 0.d0)) then
      print *
      print *, 'WARNING: lmd /= 0 diffuses ln(rho); this term does not conserve mass'
      print *
    endif
    if(proc0 .and. svv_eps0 > 0.d0 .and. trim(svv_form) == 'entropy' .and. filter_lnrho) then
      print *
      print *, 'WARNING: filter_lnrho = T filters the SVV term of ln(rho); dE/dt <= 0 of the entropy-stable SVV'
      print *, '         is then not guaranteed (see R_svv in budget.dat)'
      print *
    endif

    filt_on(ilnrho) = filter_lnrho
    filt_on(imx:imz) = filter_mom
    filt_on(ibx:ibz) = filter_mag
    dcoef(ilnrho) = lmd; dcoef_h(ilnrho) = lmd_h; dcoef_exp(ilnrho) = lmd_h_exp
    dcoef(imx:imz) = nu ; dcoef_h(imx:imz) = nu_h ; dcoef_exp(imx:imz) = nu_h_exp
    dcoef(ibx:ibz) = eta; dcoef_h(ibx:ibz) = eta_h; dcoef_exp(ibx:ibz) = eta_h_exp
    !$acc update device(filt_on, dcoef, dcoef_h, dcoef_exp)
    !$acc update device(eps_svv)

    call init_shock_viscosity

    allocate(kxt(nkx, nky_local))
    allocate(kxd(nkx, nky_local))
    call set_kx(tsc)

    !$acc parallel loop collapse(3) present(k2t, k2t_inv, kxt, ky, kz)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          k2t(k, j, i) = kxt(i, j)**2 + ky(j)**2 + kz(k)**2
          if(k2t(k, j, i) == 0.d0) then
            k2t_inv(k, j, i) = 0.d0
          else
            k2t_inv(k, j, i) = 1.0d0/k2t(k, j, i)
          endif
        enddo
      enddo
    enddo

    call allocate_advance

    ! The RK stages keep the filtered fields; start from the filtered initial condition
    if(proc0) then
      print *
      print '(" MHD_COMP_ISOTH: the initial condition is filtered (filter_lnrho, filter_mom, filter_mag = ", 3L2, ")")', &
            filter_lnrho, filter_mom, filter_mag
      print *
    endif
    call apply_filter(lnrho, mx, my, mz, bx, by, bz)

    if(proc0) then
      call open_output_file (cfl_unit, 'cfl.dat')
      write(cfl_unit, '(a)') '# 1:t 2:dt 3:dt_adv 4:dt_visc 5:max(|ux|+cf_x) 6:max(|uy|+cf_y) 7:max(|uz|+cf_z) ' &
                           //'8:zeta_max 9:max(rho*zeta)/(rho_min*zeta_max) 10:eps_svv_x 11:eps_svv_y 12:eps_svv_z ' &
                           //'13:lam_svv 14:F_svv(=13/sum_d eps_d max k_d^2 Q) 15:dt*lam_svv 16:n_iter 17:eta_max'
      write(cfl_unit, '(a)') '# 13-16: power method of the entropy-stable SVV (0 for svv_form = linear); ' &
                           //'17: shock resistivity (0 for C_eta = 0)'

      if(write_intvl_cons >= 0.d0) then
        call open_output_file (cons_unit, 'conservation.dat')
        write(cons_unit, '(a)') '# 1:t 2:dt 3:M 4:(M-M0)/M0 5:Px 6:Py 7:Pz 8:E_kin 9:E_mag 10:E_int 11:E ' &
                              //'12:rho_min 13:rho_max 14:zeta_max 15:D_min 16:sum(rho*zeta*D^2)dV ' &
                              //'17:dmass_corr 18:dmass_corr_tot 19:rms(div_N B)/(k_N rms|B|) 20:max|div_N B|/(k_N max|B|)'
        write(cons_unit, '(a)') '# M and E are before the mass correction of this step; M0 = M at the first step of this run'
        write(cons_unit, '(a)') '# 19, 20: k_N = sqrt(sum_d max k_d^2) (19, 20 <= 1)'
      endif

      if(budget) then
        call open_output_file (budget_unit, 'budget.dat')
        write(budget_unit, '(a)') '# one line per step [t, t+dt]; rates are sum_i b_i X_i over the RK stages, b = (1/4,3/16,9/16)'
        write(budget_unit, '(a)') '# 1:t 2:dt 3:dMdt_actual 4:dMdt_raw 5:dMdt_proj 6:dMdt_abs 7:dEdt_actual 8:dEdt_sum ' &
                                //'9:err_t(=7-8) 10:eps_ideal(=eps_N-11-15-21-27) 11:eps_shock 12:-visc_diss 13:R_zeta ' &
                                //'14:max|R_zeta|/visc_diss 15:eps_force 16:eps_res 17:eps_nu 18:eps_lmd 19:E_n 20:E_n+1 ' &
                                //'21:eps_svv(all) 22:eps_svv_B 23:-D_svv 24:R_svv(=21-22-23) 25:max|R_svv|/D_svv ' &
                                //'26:max|R_svv|/sum|dE/dU.U_svv|dV 27:eps_eta 28:-eta_diss 29:R_eta(=27-28) ' &
                                //'30:max|R_eta|/eta_diss'
        write(budget_unit, '(a)') '# 23-26: entropy-stable SVV of (rho, m) (0 for svv_form = linear); 26 is relative to the ' &
                                //'rounding scale of 21-22 (25 is not when D_svv is itself at the rounding level)'
        write(budget_unit, '(a)') '# 27-30: shock resistivity, eta_diss = sum eta |J|^2 dV (0 for C_eta = 0)'
      endif
    endif

  end subroutine init_work_fields


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   kxt = kx + q*shear_flg*t*ky and the
!!          x-wavenumber for derivatives
!!
!!          Computed on the host (nkx*nky_local
!!          values) and copied to the device
!-----------------------------------------------!
  subroutine set_kx(t)
    use grid, only: kx, ky
    use grid, only: nkx, nky_local
    use params, only: q, shear
    use shearing_box, only: shear_flg
    use shock_viscosity, only: kx_d
    implicit none
    real(8), intent(in) :: t
    integer :: i, j

    do j = 1, nky_local
      do i = 1, nkx
        kxt(i,j) = kx(i) + q*shear_flg*t*ky(j)
        if(shear) then
          kxd(i,j) = kxt(i,j)
        else
          kxd(i,j) = kx_d(i)
        endif
      enddo
    enddo
    !$acc update device(kxt, kxd)

  end subroutine set_kx


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Calculate nonlinear terms via
!!          1. Inverse FFT of (s, m, B)
!!          2. Validity check, rho = exp(s),
!!             monitor and mass correction
!!             (istage <= 1), u = m/rho
!!          3. D, zeta, q; SVV coefficients and
!!             CFL (istage = 1)
!!          4. Fluxes in real space, SVV of ln rho
!!          5. Forward FFT
!!          istage = 1, 2, 3: RK stage
!!          istage = 0      : steps 1-3 only, for
!!                            the final state
!!
!!          Gallope: the FFTs of cuFFTMp keep their
!!          input and take timer_fft themselves.
!!          Grid sums and maxima are device
!!          reductions over k = 1, nlz.
!-----------------------------------------------!
  subroutine get_nonlinear_terms(lnrho, mx, my, mz, bx, by, bz, istage)
    use grid, only: dlx, dly, dlz
    use grid, only: nlx, nly, nlz
    use grid, only: nlx_local, nkx, nky_local, nkz
    use params, only: zi, cs2va2, C_shock, cfl_visc, filter_mom, mass_correction, debug_nan_step
    use params, only: write_intvl_cons, write_intvl_1D, cfl_svv, svv_pi_check, C_eta, filter_mag
    use params, only: write_intvl_planes, svv_mass, svv_pi_init
    use dealias, only: filter
    use mp, only: proc0, max_allreduce, sum_allreduce, comm_fft
    use time, only: cfl, dt, tt, reset_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms
    use cufftmp, only: ftran_r2c, btran_c2r
    use advance_common, only: dt_adjust_while_running
    use shock_viscosity, only: get_div_u, get_zeta, ky_d, kz_d, kd2_max
    use shock_viscosity, only: active_x, active_y, active_z, is_nonfinite, abort_run
    use shock_viscosity, only: svv_on, kq2x, kq2y, kq2z, svv_lam_m, svv_lam_s, svv_lam_b, get_svv_eps
    use shock_viscosity, only: svv_es, get_svv_es, svv_power, svv_lam, svv_seed_nomass
    use shock_viscosity, only: tmr_on, tmr_off, itmr_zeta, itmr_svv_es, itmr_svv_pi, itmr_planes
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    complex(8), dimension (:,:,:), intent(in) :: lnrho, mx, my, mz, bx, by, bz
    integer, intent(in) :: istage

    real(8) :: mon_pre(nmon), sums(3), maxs(9), divb_err(2)
    real(8) :: vad2_r(3), va2_r, cf2_r, b2h_r, ptot_r
    real(8) :: dV, mfilt, bfilt, dt_adv, dt_visc, dt_cfl, rate_visc, F_svv
    integer :: i, j, k, n, d, pi_nit, pi_nit0
    logical :: first, pi_check, pi_conv, pi_seeded
    real(8), parameter :: rk3_lim = 2.51d0 ! real-axis stability limit of SSPRK3 (2.5127)
    character(len=256) :: msg
    ! device reduction variables and host scalars of the kernels
    real(8) :: s1, s2, s3, m1, m2, m3, m4, m5, m6, m7, m8, m9, e1, nan_val
    logical :: zeta_on

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    dV    = dlx*dly*dlz
    first = istage == 1 .and. istep_adv == 0

    ! 1. Inverse FFT
    ! for some reason c2r_many doesn't work at Fugaku
    ! (btran_c2r keeps its input, so the state is transformed without the copy into w)
    call btran_c2r(lnrho, w_r(:,:,:,ilnrho))
    call btran_c2r(mx   , w_r(:,:,:,imx   ))
    call btran_c2r(my   , w_r(:,:,:,imy   ))
    call btran_c2r(mz   , w_r(:,:,:,imz   ))
    call btran_c2r(bx   , w_r(:,:,:,ibx   ))
    call btran_c2r(by   , w_r(:,:,:,iby   ))
    call btran_c2r(bz   , w_r(:,:,:,ibz   ))

    ! 2. Validity check (no floor), rho, monitor, mass correction
    if(istage == 1 .and. istep_adv == debug_nan_step .and. proc0) then
      print *, 'debug_nan_step: NaN is injected into mx at the first grid point of proc 0'
      nan_val = ieee_value(1.d0, ieee_quiet_nan)
      !$acc serial present(w_r)
      w_r(1, 1, 1, imx) = nan_val
      !$acc end serial
    endif
    call check_state(istage)

    !$acc parallel loop collapse(3) present(w_r, rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          rho_r(k,j,i) = exp(w_r(k,j,i,ilnrho))
        enddo
      enddo
    enddo

    if(istage <= 1) then
      call get_monitor_sums(mon_pre)
      if(first) mass0 = mon_pre(imon_M)

      ! close the budget of the previous step with the state before the mass correction
      if(bgt_pending) call write_budget(mon_pre)

      dmass_corr = 0.d0
      if(istage == 1 .and. mass_correction .and. .not. first) then
        ! uniform shift s <- s + ln(M0/M); this does not restore the local conservation
        dmass_corr     = log(mass0/mon_pre(imon_M))
        dmass_corr_tot = dmass_corr_tot + dmass_corr
        !$acc parallel loop collapse(3) present(w_r, rho_r)
        do i = 1, nlx_local
          do j = 1, nly
            do k = 1, nlz
              w_r(k,j,i,ilnrho) = w_r(k,j,i,ilnrho) + dmass_corr
              rho_r(k,j,i) = exp(w_r(k,j,i,ilnrho))
            enddo
          enddo
        enddo
        call get_monitor_sums(mon_n)
      else
        mon_n = mon_pre
      endif
    endif

    !$acc parallel loop collapse(3) present(w_r, rho_r, ux_r, uy_r, uz_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          ux_r(k,j,i) = w_r(k,j,i,imx)/rho_r(k,j,i)
          uy_r(k,j,i) = w_r(k,j,i,imy)/rho_r(k,j,i)
          uz_r(k,j,i) = w_r(k,j,i,imz)/rho_r(k,j,i)
        enddo
      enddo
    enddo

    ! 3. Shock viscosity: D = div_N u, zeta, q = rho*zeta*D; shock resistivity eta
    call tmr_on(itmr_zeta)
    if(svv_es) then
      call get_div_u(ux_r, uy_r, uz_r, kxd, div_r, svv_uh)
    elseif(C_shock > 0.d0 .or. C_eta > 0.d0 .or. istage <= 1) then
      call get_div_u(ux_r, uy_r, uz_r, kxd, div_r)
    endif
    if(C_eta > 0.d0) then
      call get_zeta(div_r, zeta_r, eta_r)
    elseif(C_shock > 0.d0) then
      call get_zeta(div_r, zeta_r)
    endif
    ! q = rho*zeta*D if C_eta > 0 or C_shock > 0; zeta = q = 0 otherwise
    zeta_on = C_eta > 0.d0 .or. C_shock > 0.d0
    !$acc parallel loop collapse(3) present(rho_r, zeta_r, div_r, q_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          if(zeta_on) then
            q_r(k,j,i) = rho_r(k,j,i)*zeta_r(k,j,i)*div_r(k,j,i)
          else
            zeta_r(k,j,i) = 0.d0
            q_r   (k,j,i) = 0.d0
          endif
        enddo
      enddo
    enddo
    call tmr_off(itmr_zeta)

    ! div_N m (w is used as a work array; its content is destroyed by btran anyway)
    !$acc parallel loop collapse(3) present(w, mx, my, mz, filter, kxd, ky_d, kz_d) private(mfilt)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          mfilt = merge(filter(k,j,i), 1.d0, filter_mom)
          w(k,j,i,1) = zi*(kxd(i,j)*mx(k,j,i) + ky_d(j)*my(k,j,i) + kz_d(k)*mz(k,j,i))*mfilt
        enddo
      enddo
    enddo
    call btran_c2r(w(:,:,:,1), flx_r(:,:,:,iflx_s))

    s1 = 0.d0
    s2 = 0.d0
    s3 = 0.d0
    !$acc parallel loop collapse(3) present(flx_r, q_r, div_r) reduction(+:s1, s2, s3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          s1 = s1 + flx_r(k,j,i,iflx_s)
          s2 = s2 + abs(flx_r(k,j,i,iflx_s))
          s3 = s3 + q_r(k,j,i)*div_r(k,j,i)
        enddo
      enddo
    enddo
    sums(1) = -s1
    sums(2) =  s2
    sums(3) =  s3
    call sum_allreduce(sums, comm=comm_fft)
    dMdt_raw  = sums(1)*dV
    dMdt_abs  = sums(2)*dV
    visc_diss = sums(3)*dV

    if(istage <= 1) then
      ! maxs: max rho, -min rho, max zeta, -min D, max rho*zeta, max(|u_d| + cf_d) for d = x, y, z, max eta
      m1 = -huge(1.d0)
      m2 = -huge(1.d0)
      m3 = -huge(1.d0)
      m4 = -huge(1.d0)
      m5 = -huge(1.d0)
      m6 = 0.d0
      m7 = 0.d0
      m8 = 0.d0
      m9 = 0.d0
      !$acc parallel loop collapse(3) present(w_r, rho_r, zeta_r, div_r, ux_r, uy_r, uz_r) &
      !$acc   private(vad2_r, va2_r, cf2_r) reduction(max:m1, m2, m3, m4, m5, m6, m7, m8)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            m1 = max(m1,  rho_r(k,j,i))
            m2 = max(m2, -rho_r(k,j,i))
            m3 = max(m3,  zeta_r(k,j,i))
            m4 = max(m4, -div_r(k,j,i))
            m5 = max(m5,  rho_r(k,j,i)*zeta_r(k,j,i))

            vad2_r(1) = w_r(k,j,i,ibx)**2/rho_r(k,j,i)
            vad2_r(2) = w_r(k,j,i,iby)**2/rho_r(k,j,i)
            vad2_r(3) = w_r(k,j,i,ibz)**2/rho_r(k,j,i)
            va2_r     = vad2_r(1) + vad2_r(2) + vad2_r(3)
            ! d = 1, 2, 3 of Calliope's loop, unrolled
            cf2_r = 0.5d0*(cs2va2 + va2_r + sqrt(max(0.d0, (cs2va2 + va2_r)**2 - 4.d0*cs2va2*vad2_r(1))))
            m6 = max(m6, abs(ux_r(k,j,i)) + sqrt(cf2_r))
            cf2_r = 0.5d0*(cs2va2 + va2_r + sqrt(max(0.d0, (cs2va2 + va2_r)**2 - 4.d0*cs2va2*vad2_r(2))))
            m7 = max(m7, abs(uy_r(k,j,i)) + sqrt(cf2_r))
            cf2_r = 0.5d0*(cs2va2 + va2_r + sqrt(max(0.d0, (cs2va2 + va2_r)**2 - 4.d0*cs2va2*vad2_r(3))))
            m8 = max(m8, abs(uz_r(k,j,i)) + sqrt(cf2_r))
          enddo
        enddo
      enddo
      if(C_eta > 0.d0) then
        m9 = -huge(1.d0)
        !$acc parallel loop collapse(3) present(eta_r) reduction(max:m9)
        do i = 1, nlx_local
          do j = 1, nly
            do k = 1, nlz
              m9 = max(m9, eta_r(k,j,i))
            enddo
          enddo
        enddo
      endif
      maxs = (/m1, m2, m3, m4, m5, m6, m7, m8, m9/)
      call max_allreduce(maxs, comm=comm_fft)
    endif

    ! SVV: eps_d (fixed over the step) and, for svv_form = 'entropy', the terms of rho and m
    ! and the power method for lambda_svv (one iteration per step; until converged at the check steps)
    if(istage == 1) then
      call get_svv_eps(maxs(6:8), eps_svv)
      !$acc update device(eps_svv)
    endif
    pi_nit   = 0
    pi_check = .false.
    pi_conv  = .false.
    if(svv_es .and. istage >= 1) then
      call tmr_on(itmr_svv_es)
      call get_svv_es(w_r(:,:,:,ilnrho), rho_r, ux_r, uy_r, uz_r, eps_svv, svv_r, svv_mh, svv_diss, svv_uh)
      call tmr_off(itmr_svv_es)
      if(istage == 1) then
        call tmr_on(itmr_svv_pi)
        ! svv_mass = F: svv_seed_nomass sets the initial vector (also in m) before svv_power would (s only), and
        ! svv_pi_init - 1 iterations here and the call below make the svv_pi_init of svv_power's first call
        if(.not. svv_mass .and. .not. all(eps_svv == 0.d0)) then
          pi_nit0 = 0
          call svv_seed_nomass(pi_seeded)
          if(pi_seeded) call svv_power(w_r(:,:,:,ilnrho), rho_r, w_r(:,:,:,imx), w_r(:,:,:,imy), w_r(:,:,:,imz), &
                                       eps_svv, svv_r, svv_mh, svv_pi_init - 1, 0.d0, pi_nit0, pi_conv)
        endif
        pi_check = svv_pi_check > 0 .and. istep_adv > 0 .and. mod(istep_adv, svv_pi_check) == 0
        if(pi_check) then
          call svv_power(w_r(:,:,:,ilnrho), rho_r, w_r(:,:,:,imx), w_r(:,:,:,imy), w_r(:,:,:,imz), eps_svv, &
                         svv_r, svv_mh, 200, 1.d-3, pi_nit, pi_conv)
          if(.not. pi_conv .and. proc0) &
            print '(" WARNING: the power method for lambda_svv did not converge in ", i4, " iterations")', pi_nit
        else
          call svv_power(w_r(:,:,:,ilnrho), rho_r, w_r(:,:,:,imx), w_r(:,:,:,imy), w_r(:,:,:,imz), eps_svv, &
                         svv_r, svv_mh, 1, 0.d0, pi_nit, pi_conv)
        endif
        call tmr_off(itmr_svv_pi)
        ! svv_mass = F: J has no ln rho row, so the initial vector is also in m (J v = 0 is then not expected);
        ! lambda_svv = 0 would drop the SVV limit on dt
        if(.not. svv_mass .and. .not. all(eps_svv == 0.d0)) then
          pi_nit = pi_nit + pi_nit0
          if(svv_lam == 0.d0) &
            call abort_run('get_nonlinear_terms: svv_mass = F and lambda_svv = 0; the SVV would not limit dt')
        endif
      endif
    endif

    ! CFL: advection by fast waves over the active directions, and the explicit viscous terms
    if(istage == 1) then
      dt_prev = dt

      ! When the file 'dt_adjust' including a float number is created,
      ! dt will be manually adjusted to that value while running.
      call dt_adjust_while_running()
      ! dt_adjust_while_running sets the host dt only
      !$acc update device(dt)

      dt_adv = 0.d0
      if(active_x) dt_adv = dt_adv + maxs(6)/dlx
      if(active_y) dt_adv = dt_adv + maxs(7)/dly
      if(active_z) dt_adv = dt_adv + maxs(8)/dlz
      dt_adv = cfl/dt_adv

      if(svv_es) then
        ! 1/dt_visc: shock viscosity/cfl_visc + lambda_svv/cfl_svv in (rho, m), (SVV + eta) in B/cfl_visc
        rate_visc = max(maxs(3)*kd2_max/cfl_visc + svv_lam/cfl_svv, (sum(eps_svv*svv_lam_b) + maxs(9)*kd2_max)/cfl_visc)
        if(rate_visc > 0.d0) then
          dt_visc = 1.d0/rate_visc
        else
          dt_visc = huge(1.d0)
        endif
      else
        ! largest damping rate: shock viscosity + SVV in m, SVV in ln rho and SVV + eta in B
        rate_visc = max(maxs(3)*kd2_max + sum(eps_svv*svv_lam_m), sum(eps_svv*svv_lam_s), &
                        sum(eps_svv*svv_lam_b) + maxs(9)*kd2_max)
        if(rate_visc > 0.d0) then
          dt_visc = cfl_visc/rate_visc
        else
          dt_visc = huge(1.d0)
        endif
      endif

      dt_cfl = min(dt_adv, dt_visc)
      if(is_nonfinite(dt_cfl) .or. dt_cfl <= 0.d0) call abort_run('get_nonlinear_terms: invalid dt_cfl')

      ! reset_method = 'track' raises dt only with the converged lambda_svv
      call reset_dt(dt_cfl, counter, (pi_check .and. pi_conv) .or. .not. svv_es)
      n = 0
      do while(dt_cfl < dt .and. n < 100)
        call reset_dt(dt_cfl, counter)
        n = n + 1
      enddo
      if(dt_cfl < dt) call abort_run('get_nonlinear_terms: dt could not be reduced below dt_cfl (check reset_method)')

      ! check of the explicit SVV step with the converged lambda_svv: dt of the previous and of this step
      if(pi_check .and. max(dt_prev, dt)*svv_lam > rk3_lim) then
        write(msg, '("get_nonlinear_terms: dt*lambda_svv = ", es12.4, " > ", f5.2, " (dt_prev = ", es12.4, &
                    &", dt = ", es12.4, ", lambda_svv = ", es12.4, ", step = ", i10, ")")') &
              max(dt_prev, dt)*svv_lam, rk3_lim, dt_prev, dt, svv_lam, istep_adv
        call abort_run(trim(msg))
      endif

      F_svv = 0.d0
      if(svv_es .and. sum(eps_svv*svv_lam_s) > 0.d0) F_svv = svv_lam/sum(eps_svv*svv_lam_s)
      if(proc0) then
        write (unit=cfl_unit, fmt="(100es26.17e3)") tt, dt, dt_adv, dt_visc, maxs(6:8), maxs(3), &
                                                  maxs(5)/max(-maxs(2)*maxs(3), tiny(1.d0)), eps_svv, &
                                                  merge(svv_lam, 0.d0, svv_es), F_svv, &
                                                  merge(dt*svv_lam, 0.d0, svv_es), dble(pi_nit), maxs(9)
        call flush(cfl_unit)
      endif
    endif

    ! Output of the state at tt
    if(istage <= 1) then
      if(write_intvl_cons == 0.d0 .or. (write_intvl_cons > 0.d0 .and. &
         (first .or. istage == 0 .or. write_due(write_intvl_cons)))) then
        call get_divb_err(bx, by, bz, divb_err)
        if(proc0) then
          write (unit=cons_unit, fmt="(100es26.17e3)") tt, dt, mon_pre(imon_M), (mon_pre(imon_M) - mass0)/mass0, &
                mon_pre(imon_px:imon_pz), mon_pre(imon_ek:imon_ei), sum(mon_pre(imon_ek:imon_ei)), &
                -maxs(2), maxs(1), maxs(3), -maxs(4), visc_diss, dmass_corr, dmass_corr_tot, divb_err
          call flush(cons_unit)
        endif
      endif

      if(write_intvl_1D > 0.d0 .and. (first .or. istage == 0 .or. write_due(write_intvl_1D))) call write_line
      if(write_intvl_planes > 0.d0 .and. (first .or. istage == 0 .or. write_due(write_intvl_planes))) then
        call tmr_on(itmr_planes)
        call write_planes
        call tmr_off(itmr_planes)
      endif
    endif

    if(istage == 0) then
      if (proc0) call put_time_stamp(timer_nonlinear_terms)
      return
    endif

    ! 4. Calculate nonlinear terms in real space
    !$acc parallel loop collapse(3) present(flx_r, rho_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          flx_r(k,j,i,iflx_s) = -flx_r(k,j,i,iflx_s)/rho_r(k,j,i)
        enddo
      enddo
    enddo

    ! SVV of the mass equation: ds/dt += (d rho/dt)_svv/rho, so that sum rho ds/dt = 0 is kept
    if(svv_es) then
      ! 'entropy': (d rho/dt)_svv from get_svv_es
      !$acc parallel loop collapse(3) present(flx_r, svv_r, rho_r)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            svv_r(k,j,i) = svv_r(k,j,i)/rho_r(k,j,i)
            flx_r(k,j,i,iflx_s) = flx_r(k,j,i,iflx_s) + svv_r(k,j,i)
          enddo
        enddo
      enddo
    elseif(svv_on) then
      ! 'linear': d rho/dt = F^-1(-lambda F rho) (lambda(0) = 0)
      ! (ftran_r2c keeps its input, so rho_r is transformed without the copy into svv_r)
      call ftran_r2c(rho_r, w(:,:,:,1))
      !$acc parallel loop collapse(3) present(w, kq2x, kq2y, kq2z, eps_svv)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            w(k,j,i,1) = -(eps_svv(1)*kq2x(i) + eps_svv(2)*kq2y(j) + eps_svv(3)*kq2z(k))*w(k,j,i,1)/nlx/nly/nlz
          enddo
        enddo
      enddo
      call btran_c2r(w(:,:,:,1), svv_r)
      !$acc parallel loop collapse(3) present(flx_r, svv_r, rho_r)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            svv_r(k,j,i) = svv_r(k,j,i)/rho_r(k,j,i)
            flx_r(k,j,i,iflx_s) = flx_r(k,j,i,iflx_s) + svv_r(k,j,i)
          enddo
        enddo
      enddo
    endif

    ! Shock resistivity: J = curl_N B = F^-1 i k_d x M B^ (M: filter if filter_mag, as in D),
    ! eta_diss = sum eta |J|^2 dV, and eta J is added to the flux B x u of the induction equation
    if(C_eta > 0.d0) then
      !$acc parallel loop collapse(3) present(w, bx, by, bz, filter, kxd, ky_d, kz_d) private(bfilt)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            bfilt = merge(filter(k,j,i), 1.d0, filter_mag)
            w(k,j,i,1) = zi*(ky_d(j)*bz(k,j,i) - kz_d(k)*by(k,j,i))*bfilt
            w(k,j,i,2) = zi*(kz_d(k)*bx(k,j,i) - kxd(i,j)*bz(k,j,i))*bfilt
            w(k,j,i,3) = zi*(kxd(i,j)*by(k,j,i) - ky_d(j)*bx(k,j,i))*bfilt
          enddo
        enddo
      enddo
      do d = 1, 3
        call btran_c2r(w(:,:,:,d), ej_r(:,:,:,d))
      enddo
      e1 = 0.d0
      !$acc parallel loop collapse(3) present(eta_r, ej_r) reduction(+:e1)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            e1 = e1 + eta_r(k,j,i)*(ej_r(k,j,i,1)**2 + ej_r(k,j,i,2)**2 + ej_r(k,j,i,3)**2)
          enddo
        enddo
      enddo
      eta_diss = e1
      call sum_allreduce(eta_diss, comm=comm_fft)
      eta_diss = eta_diss*dV
      !$acc parallel loop collapse(3) present(eta_r, ej_r)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            ej_r(k,j,i,1) = eta_r(k,j,i)*ej_r(k,j,i,1)
            ej_r(k,j,i,2) = eta_r(k,j,i)*ej_r(k,j,i,2)
            ej_r(k,j,i,3) = eta_r(k,j,i)*ej_r(k,j,i,3)
          enddo
        enddo
      enddo
    endif

    !$acc parallel loop collapse(3) present(w_r, flx_r, rho_r, ux_r, uy_r, uz_r, q_r) private(b2h_r, ptot_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          b2h_r  = 0.5d0*(w_r(k,j,i,ibx)**2 + w_r(k,j,i,iby)**2 + w_r(k,j,i,ibz)**2)
          ptot_r = cs2va2*rho_r(k,j,i) + b2h_r - q_r(k,j,i)

          flx_r(k,j,i,iflx_mxx) = rho_r(k,j,i)*ux_r(k,j,i)*ux_r(k,j,i) - w_r(k,j,i,ibx)*w_r(k,j,i,ibx) + ptot_r
          flx_r(k,j,i,iflx_mxy) = rho_r(k,j,i)*ux_r(k,j,i)*uy_r(k,j,i) - w_r(k,j,i,ibx)*w_r(k,j,i,iby)
          flx_r(k,j,i,iflx_mxz) = rho_r(k,j,i)*ux_r(k,j,i)*uz_r(k,j,i) - w_r(k,j,i,ibx)*w_r(k,j,i,ibz)
          flx_r(k,j,i,iflx_myy) = rho_r(k,j,i)*uy_r(k,j,i)*uy_r(k,j,i) - w_r(k,j,i,iby)*w_r(k,j,i,iby) + ptot_r
          flx_r(k,j,i,iflx_myz) = rho_r(k,j,i)*uy_r(k,j,i)*uz_r(k,j,i) - w_r(k,j,i,iby)*w_r(k,j,i,ibz)
          flx_r(k,j,i,iflx_mzz) = rho_r(k,j,i)*uz_r(k,j,i)*uz_r(k,j,i) - w_r(k,j,i,ibz)*w_r(k,j,i,ibz) + ptot_r

          flx_r(k,j,i,iflx_bx ) = w_r(k,j,i,iby)*uz_r(k,j,i) - w_r(k,j,i,ibz)*uy_r(k,j,i)
          flx_r(k,j,i,iflx_by ) = w_r(k,j,i,ibz)*ux_r(k,j,i) - w_r(k,j,i,ibx)*uz_r(k,j,i)
          flx_r(k,j,i,iflx_bz ) = w_r(k,j,i,ibx)*uy_r(k,j,i) - w_r(k,j,i,iby)*ux_r(k,j,i)
        enddo
      enddo
    enddo
    if(C_eta > 0.d0) then
      !$acc parallel loop collapse(3) present(flx_r, ej_r)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            flx_r(k,j,i,iflx_bx) = flx_r(k,j,i,iflx_bx) + ej_r(k,j,i,1)
            flx_r(k,j,i,iflx_by) = flx_r(k,j,i,iflx_by) + ej_r(k,j,i,2)
            flx_r(k,j,i,iflx_bz) = flx_r(k,j,i,iflx_bz) + ej_r(k,j,i,3)
          enddo
        enddo
      enddo
    endif

    ! 5. Forward FFT
    ! for some reason r2c_many doesn't work at Fugaku
    do n = 1, nftran
      call ftran_r2c(flx_r(:,:,:,n), flx(:,:,:,n))
    enddo
    !$acc parallel loop collapse(4) present(flx)
    do n = 1, nftran
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            flx(k,j,i,n) = flx(k,j,i,n)/nlx/nly/nlz
          enddo
        enddo
      enddo
    enddo

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Stop if a field is non-finite or
!!          |ln rho| > s_abort (no correction)
!!
!!          Gallope: the device counts the bad
!!          points; w_r is copied to the host only
!!          when the state is invalid
!-----------------------------------------------!
  subroutine check_state(istage)
    use grid, only: nlx_local, nly, nlz
    use params, only: s_abort
    use mp, only: proc_id, iproc_fft
    use time, only: tt
    use shock_viscosity, only: is_nonfinite, abort_run
    implicit none
    integer, intent(in) :: istage
    character(len=5), parameter :: fname(nfields) = (/'lnrho', 'mx   ', 'my   ', 'mz   ', 'bx   ', 'by   ', 'bz   '/)
    integer :: i, j, k, n
    integer :: nbad
    real(8) :: smax

    nbad = 0
    !$acc parallel loop collapse(4) present(w_r) reduction(+:nbad)
    do n = 1, nfields
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            if(is_nonfinite(w_r(k,j,i,n))) nbad = nbad + 1
          enddo
        enddo
      enddo
    enddo
    smax = 0.d0
    !$acc parallel loop collapse(3) present(w_r) reduction(max:smax)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          smax = max(smax, abs(w_r(k,j,i,ilnrho)))
        enddo
      enddo
    enddo

    if(.not. (nbad > 0 .or. smax > s_abort)) return

    !$acc update self(w_r)
    do n = 1, nfields
      do i = 1, nlx_local
        do k = 1, nlz
          do j = 1, nly
            if(is_nonfinite(w_r(k,j,i,n)) .or. (n == ilnrho .and. abs(w_r(k,j,i,n)) > s_abort)) then
              write(*, '("[proc ", i6, "] invalid state at t = ", es24.16, ", step ", i10, ", stage ", i2, &
                       &": ", a, "(iy, iz, ix) = (", i8, ",", i8, ",", i8, ") = ", es24.16)') &
                    proc_id, tt, istep_adv + 1, istage, trim(fname(n)), j, k, i + nlx_local*iproc_fft, w_r(k,j,i,n)
              call abort_run('non-finite value or |ln rho| > s_abort')
            endif
          enddo
        enddo
      enddo
    enddo

  end subroutine check_state


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Volume integrals of the state in
!!          w_r and rho_r (grid sums times dV)
!-----------------------------------------------!
  subroutine get_monitor_sums(mon)
    use grid, only: dlx, dly, dlz
    use grid, only: nlx_local, nly, nlz
    use params, only: cs2va2
    use mp, only: sum_allreduce, comm_fft
    implicit none
    real(8), intent(out) :: mon(nmon)
    integer :: i, j, k
    real(8) :: sM, spx, spy, spz, sek, sem, sei

    sM  = 0.d0
    spx = 0.d0
    spy = 0.d0
    spz = 0.d0
    sek = 0.d0
    sem = 0.d0
    sei = 0.d0
    !$acc parallel loop collapse(3) present(w_r, rho_r) reduction(+:sM, spx, spy, spz, sek, sem, sei)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          sM  = sM  + rho_r(k,j,i)
          spx = spx + w_r(k,j,i,imx)
          spy = spy + w_r(k,j,i,imy)
          spz = spz + w_r(k,j,i,imz)
          sek = sek + 0.5d0*(w_r(k,j,i,imx)**2 + w_r(k,j,i,imy)**2 + w_r(k,j,i,imz)**2)/rho_r(k,j,i)
          sem = sem + 0.5d0*(w_r(k,j,i,ibx)**2 + w_r(k,j,i,iby)**2 + w_r(k,j,i,ibz)**2)
          sei = sei + cs2va2*rho_r(k,j,i)*(w_r(k,j,i,ilnrho) - 1.d0)
        enddo
      enddo
    enddo
    mon(imon_M ) = sM
    mon(imon_px) = spx
    mon(imon_py) = spy
    mon(imon_pz) = spz
    mon(imon_ek) = sek
    mon(imon_em) = sem
    mon(imon_ei) = sei
    call sum_allreduce(mon, comm=comm_fft)
    mon = mon*dlx*dly*dlz

  end subroutine get_monitor_sums


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Same criterion as check_write_now
!!          in calliope.F90
!-----------------------------------------------!
  function write_due(intvl)
    use time, only: tt, dt
    implicit none
    real(8), intent(in) :: intvl
    logical :: write_due
    integer, parameter :: i16 = selected_int_kind(16)

    write_due = floor(tt/intvl, i16) /= floor((tt+dt)/intvl, i16)

  end function write_due


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Add a to the k = 0 mode of f
!!          (uniform shift in real space)
!-----------------------------------------------!
  subroutine shift_mean(f, a)
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    implicit none
    complex(8), dimension (:,:,:), intent(inout) :: f
    real(8), intent(in) :: a
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(f, kx, ky, kz)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) f(k,j,i) = f(k,j,i) + a
        enddo
      enddo
    enddo

  end subroutine shift_mean


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Multiply each field by the dealiasing
!!          filter if filter_* of the field is T
!-----------------------------------------------!
  subroutine apply_filter(lnrho, mx, my, mz, bx, by, bz)
    use grid, only: nkx, nky_local, nkz
    use dealias, only: filter
    implicit none
    complex(8), dimension (:,:,:), intent(inout) :: lnrho, mx, my, mz, bx, by, bz
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(lnrho, mx, my, mz, bx, by, bz, filter, filt_on)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          lnrho(k,j,i) = lnrho(k,j,i)*merge(filter(k,j,i), 1.d0, filt_on(ilnrho))
          mx   (k,j,i) = mx   (k,j,i)*merge(filter(k,j,i), 1.d0, filt_on(imx   ))
          my   (k,j,i) = my   (k,j,i)*merge(filter(k,j,i), 1.d0, filt_on(imy   ))
          mz   (k,j,i) = mz   (k,j,i)*merge(filter(k,j,i), 1.d0, filt_on(imz   ))
          bx   (k,j,i) = bx   (k,j,i)*merge(filter(k,j,i), 1.d0, filt_on(ibx   ))
          by   (k,j,i) = by   (k,j,i)*merge(filter(k,j,i), 1.d0, filt_on(iby   ))
          bz   (k,j,i) = bz   (k,j,i)*merge(filter(k,j,i), 1.d0, filt_on(ibz   ))
        enddo
      enddo
    enddo

  end subroutine apply_filter


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Explicit terms of all modes from flx
!!          (filled by get_nonlinear_terms) and
!!          the SVV of m and B: -lambda(k) B^, and
!!          svv_mh ('entropy') or -lambda(k) m^
!!          ('linear')
!!
!!          Gallope: the SVV is added by a second
!!          kernel (svv_mh exists only for
!!          'entropy'); the operations per mode
!!          are those of Calliope
!-----------------------------------------------!
  subroutine get_explicit_terms(mx, my, mz, bx, by, bz)
    use grid, only: nkx, nky_local, nkz
    use force, only: fmx, fmy, fmz
    use shock_viscosity, only: ky_d, kz_d
    use shock_viscosity, only: svv_on, svv_es, kq2x, kq2y, kq2z
    implicit none
    complex(8), dimension (:,:,:), intent(in) :: mx, my, mz, bx, by, bz
    real(8) :: lam
    integer :: i, j, k, n
    complex(8) :: fl(nftran), ex(nfields)

    !$acc parallel loop collapse(3) present(exp_terms, flx, mx, my, mz, bx, by, bz, fmx, fmy, fmz, &
    !$acc   kxd, ky_d, kz_d) private(fl, ex)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          !$acc loop seq
          do n = 1, nftran
            fl(n) = flx(k,j,i,n)
          enddo
          call get_ext_terms(ex, &
                             mx(k,j,i), my(k,j,i), mz(k,j,i), bx(k,j,i), by(k,j,i), bz(k,j,i), &
                             fl, &
                             fmx(k,j,i), fmy(k,j,i), fmz(k,j,i), &
                             kxd(i,j), ky_d(j), kz_d(k))
          !$acc loop seq
          do n = 1, nfields
            exp_terms(k,j,i,n) = ex(n)
          enddo
        enddo
      enddo
    enddo

    if(svv_on) then
      if(svv_es) then
        !$acc parallel loop collapse(3) present(exp_terms, svv_mh, bx, by, bz, kq2x, kq2y, kq2z, eps_svv) &
        !$acc   private(lam)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              lam = eps_svv(1)*kq2x(i) + eps_svv(2)*kq2y(j) + eps_svv(3)*kq2z(k)
              exp_terms(k,j,i,imx) = exp_terms(k,j,i,imx) + svv_mh(k,j,i,1)
              exp_terms(k,j,i,imy) = exp_terms(k,j,i,imy) + svv_mh(k,j,i,2)
              exp_terms(k,j,i,imz) = exp_terms(k,j,i,imz) + svv_mh(k,j,i,3)
              exp_terms(k,j,i,ibx) = exp_terms(k,j,i,ibx) - lam*bx(k,j,i)
              exp_terms(k,j,i,iby) = exp_terms(k,j,i,iby) - lam*by(k,j,i)
              exp_terms(k,j,i,ibz) = exp_terms(k,j,i,ibz) - lam*bz(k,j,i)
            enddo
          enddo
        enddo
      else
        !$acc parallel loop collapse(3) present(exp_terms, mx, my, mz, bx, by, bz, kq2x, kq2y, kq2z, eps_svv) &
        !$acc   private(lam)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              lam = eps_svv(1)*kq2x(i) + eps_svv(2)*kq2y(j) + eps_svv(3)*kq2z(k)
              exp_terms(k,j,i,imx) = exp_terms(k,j,i,imx) - lam*mx(k,j,i)
              exp_terms(k,j,i,imy) = exp_terms(k,j,i,imy) - lam*my(k,j,i)
              exp_terms(k,j,i,imz) = exp_terms(k,j,i,imz) - lam*mz(k,j,i)
              exp_terms(k,j,i,ibx) = exp_terms(k,j,i,ibx) - lam*bx(k,j,i)
              exp_terms(k,j,i,iby) = exp_terms(k,j,i,iby) - lam*by(k,j,i)
              exp_terms(k,j,i,ibz) = exp_terms(k,j,i,ibz) - lam*bz(k,j,i)
            enddo
          enddo
        enddo
      endif
    endif

  end subroutine get_explicit_terms


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms, &
                           mx, my, mz, bx, by, bz, &
                           flx, &
                           fmx, fmy, fmz, &
                           kxd, kyd, kzd)
    !$acc routine seq
    use params, only: zi, q
    use shearing_box, only: shear_flg
    implicit none
    complex(8), intent(out) :: exp_terms(nfields)
    complex(8), intent(in ) :: mx, my, mz, bx, by, bz
    complex(8), intent(in ) :: flx(nftran)
    complex(8), intent(in ) :: fmx, fmy, fmz
    real(8)   , intent(in ) :: kxd, kyd, kzd
    complex(8)              :: nl(nfields)

    ! -(div_N m)/rho (computed in real space)
    nl(ilnrho) = flx(iflx_s)

    ! div {m u - B B + [cs^2 rho + B^2/2 - q]*I}
    nl(imx ) = -zi*( kxd*flx(iflx_mxx) + kyd*flx(iflx_mxy) + kzd*flx(iflx_mxz) )
    nl(imy ) = -zi*( kxd*flx(iflx_mxy) + kyd*flx(iflx_myy) + kzd*flx(iflx_myz) )
    nl(imz ) = -zi*( kxd*flx(iflx_mxz) + kyd*flx(iflx_myz) + kzd*flx(iflx_mzz) )

    ! curl (b x u)
    nl(ibx ) = -zi*( kyd*flx(iflx_bz) - kzd*flx(iflx_by) )
    nl(iby ) = -zi*( kzd*flx(iflx_bx) - kxd*flx(iflx_bz) )
    nl(ibz ) = -zi*( kxd*flx(iflx_by) - kyd*flx(iflx_bx) )

    exp_terms(ilnrho) = nl(ilnrho)

    exp_terms(imx ) = nl(imx) + fmx + 2.d0*shear_flg*my
    exp_terms(imy ) = nl(imy) + fmy - (2.d0 - q)*shear_flg*mx
    exp_terms(imz ) = nl(imz) + fmz

    exp_terms(ibx ) = nl(ibx)
    exp_terms(iby ) = nl(iby) - q*shear_flg*bx
    exp_terms(ibz ) = nl(ibz)

  end subroutine get_ext_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Time integral of hyperdissipation
!-----------------------------------------------!
  subroutine get_imp_terms_tintg(imp_terms_tintg, t, kx, ky, kz, coeff, coeff_h, nexp)
    !$acc routine seq
    use grid, only: k2_max
    use params, only: shear
    use shearing_box, only: get_imp_terms_tintg_with_shear
    implicit none
    real(8), intent(out) :: imp_terms_tintg
    real(8), intent(in) :: t, kx, ky, kz, coeff, coeff_h
    integer, intent(in) :: nexp

    if(shear) then
      call get_imp_terms_tintg_with_shear(imp_terms_tintg, t, kx, ky, kz, coeff, coeff_h, nexp )
    else
      imp_terms_tintg = -(coeff*((kx**2 + ky**2 + kz**2)/k2_max) + coeff_h*((kx**2 + ky**2 + kz**2)/k2_max)**nexp)*t
    endif
  end subroutine get_imp_terms_tintg


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Energy and mass rates of one RK stage
!!          dE/dt = sum [dE/ds Pds + u.Pdm + B.PdB] dV
!!          with E = m^2/(2 rho) + B^2/2
!!                 + cs^2 rho (ln rho - 1),
!!          dE/ds = rho (cs^2 s - u^2/2), and P the
!!          filter applied after the stage.
!!          eps_shock = sum u.F^-1(P i k_d q^) dV;
!!          R_zeta = eps_shock + sum rho zeta D^2 dV
!!          vanishes in exact arithmetic.
!!          eps_svv, eps_svv_b: energy rates of the SVV
!!          of (rho, m) and of B (included in eps_N).
!!          eps_svv_b <= 0. 'linear': eps_svv is not
!!          sign-definite. 'entropy': R_svv = eps_svv
!!          + D_svv vanishes in exact arithmetic and
!!          D_svv >= 0 (a sum of squares). A_svv =
!!          sum |dE/dU.U_svv| dV bounds the rounding
!!          error of eps_svv.
!!          The linear terms (integrating factor) are
!!          evaluated with the stage state.
!!
!!          Gallope: wk, wk_r, dEds_r, qh and ejh are
!!          module work arrays (allocate_advance).
!!          Each grid sum is a device reduction over
!!          k = 1, nlz (sum_prod_r, sum_abs_prod_r);
!!          only the local sums loc reach the host.
!-----------------------------------------------!
  subroutine stage_budget(istage, lnrho, mx, my, mz, bx, by, bz)
    use grid, only: dlx, dly, dlz, nlx, nly, nlz, k2_max
    use grid, only: ky, kz
    use grid, only: nlx_local, nkx, nky_local, nkz
    use params, only: zi, cs2va2, C_shock, filter_mom, C_eta, filter_mag
    use dealias, only: filter
    use force, only: driven, fmx, fmy, fmz
    use mp, only: sum_allreduce, comm_fft
    use time, only: tt, dt
    use cufftmp, only: ftran_r2c, btran_c2r
    use shock_viscosity, only: ky_d, kz_d
    use shock_viscosity, only: svv_on, svv_es, kq2x, kq2y, kq2z
    implicit none
    integer, intent(in) :: istage
    complex(8), dimension (:,:,:), intent(in) :: lnrho, mx, my, mz, bx, by, bz
    real(8) :: loc(nbgt), kd, mfilt, bfilt, k2r, lrate
    integer  :: i, j, k, n, d

    !$acc parallel loop collapse(3) present(dEds_r, rho_r, w_r, ux_r, uy_r, uz_r)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          dEds_r(k,j,i) = rho_r(k,j,i)*(cs2va2*w_r(k,j,i,ilnrho) &
                        - 0.5d0*(ux_r(k,j,i)**2 + uy_r(k,j,i)**2 + uz_r(k,j,i)**2))
        enddo
      enddo
    enddo

    loc = 0.d0

    ! Rates of the explicit terms after the filter of each field
    do n = 1, nfields
      !$acc parallel loop collapse(3) present(wk, exp_terms, filter, filt_on)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            wk(k,j,i) = exp_terms(k,j,i,n)*merge(filter(k,j,i), 1.d0, filt_on(n))
          enddo
        enddo
      enddo
      call btran_wk
      select case(n)
        case(ilnrho)
          loc(ib_eps_N    ) = loc(ib_eps_N) + sum_prod_r(dEds_r, wk_r)
          loc(ib_dMdt_proj) = sum_prod_r(rho_r, wk_r)
        case(imx); loc(ib_eps_N) = loc(ib_eps_N) + sum_prod_r(ux_r, wk_r)
        case(imy); loc(ib_eps_N) = loc(ib_eps_N) + sum_prod_r(uy_r, wk_r)
        case(imz); loc(ib_eps_N) = loc(ib_eps_N) + sum_prod_r(uz_r, wk_r)
        case(ibx); loc(ib_eps_N) = loc(ib_eps_N) + sum_prod_r(w_r(:,:,:,ibx), wk_r)
        case(iby); loc(ib_eps_N) = loc(ib_eps_N) + sum_prod_r(w_r(:,:,:,iby), wk_r)
        case(ibz); loc(ib_eps_N) = loc(ib_eps_N) + sum_prod_r(w_r(:,:,:,ibz), wk_r)
      end select
    enddo

    ! Shock viscosity: u.F^-1(P i k_d q^)
    if(C_shock > 0.d0) then
      ! (ftran_r2c keeps its input, so q_r is transformed without the copy into wk_r)
      call ftran_r2c(q_r, qh)
      !$acc parallel loop collapse(3) present(qh)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            qh(k,j,i) = qh(k,j,i)/nlx/nly/nlz
          enddo
        enddo
      enddo
      do d = 1, 3
        !$acc parallel loop collapse(3) present(wk, qh, filter, kxd, ky_d, kz_d) private(kd, mfilt)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              if(d == 1) then
                kd = kxd(i,j)
              elseif(d == 2) then
                kd = ky_d(j)
              else
                kd = kz_d(k)
              endif
              mfilt = merge(filter(k,j,i), 1.d0, filter_mom)
              wk(k,j,i) = zi*kd*qh(k,j,i)*mfilt
            enddo
          enddo
        enddo
        call btran_wk
        select case(d)
          case(1); loc(ib_eps_shock) = loc(ib_eps_shock) + sum_prod_r(ux_r, wk_r)
          case(2); loc(ib_eps_shock) = loc(ib_eps_shock) + sum_prod_r(uy_r, wk_r)
          case(3); loc(ib_eps_shock) = loc(ib_eps_shock) + sum_prod_r(uz_r, wk_r)
        end select
      enddo
    endif

    ! Shock resistivity: B.F^-1(P(-i k_d x F(eta J)))
    if(C_eta > 0.d0) then
      do d = 1, 3
        call ftran_r2c(ej_r(:,:,:,d), ejh(:,:,:,d))
      enddo
      !$acc parallel loop collapse(4) present(ejh)
      do d = 1, 3
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              ejh(k,j,i,d) = ejh(k,j,i,d)/nlx/nly/nlz
            enddo
          enddo
        enddo
      enddo
      do d = 1, 3
        !$acc parallel loop collapse(3) present(wk, ejh, filter, kxd, ky_d, kz_d) private(bfilt)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              bfilt = merge(filter(k,j,i), 1.d0, filter_mag)
              if(d == 1) then
                wk(k,j,i) = -zi*(ky_d(j)*ejh(k,j,i,3) - kz_d(k)*ejh(k,j,i,2))*bfilt
              elseif(d == 2) then
                wk(k,j,i) = -zi*(kz_d(k)*ejh(k,j,i,1) - kxd(i,j)*ejh(k,j,i,3))*bfilt
              else
                wk(k,j,i) = -zi*(kxd(i,j)*ejh(k,j,i,2) - ky_d(j)*ejh(k,j,i,1))*bfilt
              endif
            enddo
          enddo
        enddo
        call btran_wk
        loc(ib_eps_eta) = loc(ib_eps_eta) + sum_prod_r(w_r(:,:,:,ibx+d-1), wk_r)
      enddo
    endif

    ! SVV: dE/ds.F^-1(P F svv_r) + u.F^-1(P svv_m^) (eps_svv) + B.F^-1(P(-lambda B^)) (eps_svv_b),
    ! svv_m^ = svv_mh ('entropy') or -lambda m^ ('linear')
    if(svv_on) then
      call ftran_r2c(svv_r, wk)
      !$acc parallel loop collapse(3) present(wk, filter, filt_on)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            wk(k,j,i) = wk(k,j,i)/nlx/nly/nlz*merge(filter(k,j,i), 1.d0, filt_on(ilnrho))
          enddo
        enddo
      enddo
      call btran_wk
      loc(ib_eps_svv) = loc(ib_eps_svv) + sum_prod_r(dEds_r, wk_r)
      loc(ib_A_svv  ) = loc(ib_A_svv  ) + sum_abs_prod_r(dEds_r, wk_r)
      do n = imx, ibz
        ! Gallope: svv_mh exists only for 'entropy', so its branch is a kernel of its own
        if(n <= imz .and. svv_es) then
          !$acc parallel loop collapse(3) present(wk, svv_mh, filter, filt_on)
          do i = 1, nkx
            do j = 1, nky_local
              do k = 1, nkz
                wk(k,j,i) = svv_mh(k,j,i,n-imx+1)*merge(filter(k,j,i), 1.d0, filt_on(n))
              enddo
            enddo
          enddo
        else
          !$acc parallel loop collapse(3) present(wk, mx, my, mz, bx, by, bz, filter, filt_on, &
          !$acc   kq2x, kq2y, kq2z, eps_svv) private(lrate)
          do i = 1, nkx
            do j = 1, nky_local
              do k = 1, nkz
                lrate = (eps_svv(1)*kq2x(i) + eps_svv(2)*kq2y(j) + eps_svv(3)*kq2z(k)) &
                       *merge(filter(k,j,i), 1.d0, filt_on(n))
                if(n == imx) then
                  wk(k,j,i) = -lrate*mx(k,j,i)
                elseif(n == imy) then
                  wk(k,j,i) = -lrate*my(k,j,i)
                elseif(n == imz) then
                  wk(k,j,i) = -lrate*mz(k,j,i)
                elseif(n == ibx) then
                  wk(k,j,i) = -lrate*bx(k,j,i)
                elseif(n == iby) then
                  wk(k,j,i) = -lrate*by(k,j,i)
                else
                  wk(k,j,i) = -lrate*bz(k,j,i)
                endif
              enddo
            enddo
          enddo
        endif
        call btran_wk
        select case(n)
          case(imx)
            loc(ib_eps_svv) = loc(ib_eps_svv) + sum_prod_r(ux_r, wk_r)
            loc(ib_A_svv  ) = loc(ib_A_svv  ) + sum_abs_prod_r(ux_r, wk_r)
          case(imy)
            loc(ib_eps_svv) = loc(ib_eps_svv) + sum_prod_r(uy_r, wk_r)
            loc(ib_A_svv  ) = loc(ib_A_svv  ) + sum_abs_prod_r(uy_r, wk_r)
          case(imz)
            loc(ib_eps_svv) = loc(ib_eps_svv) + sum_prod_r(uz_r, wk_r)
            loc(ib_A_svv  ) = loc(ib_A_svv  ) + sum_abs_prod_r(uz_r, wk_r)
          case(ibx); loc(ib_eps_svv_b) = loc(ib_eps_svv_b) + sum_prod_r(w_r(:,:,:,ibx), wk_r)
          case(iby); loc(ib_eps_svv_b) = loc(ib_eps_svv_b) + sum_prod_r(w_r(:,:,:,iby), wk_r)
          case(ibz); loc(ib_eps_svv_b) = loc(ib_eps_svv_b) + sum_prod_r(w_r(:,:,:,ibz), wk_r)
        end select
      enddo
      ! svv_diss is a local sum (get_svv_es); it is reduced with loc
      if(svv_es) loc(ib_D_svv) = svv_diss
    endif

    ! Forcing: u.F^-1(P f^)
    if(driven) then
      do d = 1, 3
        !$acc parallel loop collapse(3) present(wk, fmx, fmy, fmz, filter)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              if(d == 1) then
                wk(k,j,i) = fmx(k,j,i)*merge(filter(k,j,i), 1.d0, filter_mom)
              elseif(d == 2) then
                wk(k,j,i) = fmy(k,j,i)*merge(filter(k,j,i), 1.d0, filter_mom)
              else
                wk(k,j,i) = fmz(k,j,i)*merge(filter(k,j,i), 1.d0, filter_mom)
              endif
            enddo
          enddo
        enddo
        call btran_wk
        select case(d)
          case(1); loc(ib_eps_force) = loc(ib_eps_force) + sum_prod_r(ux_r, wk_r)
          case(2); loc(ib_eps_force) = loc(ib_eps_force) + sum_prod_r(uy_r, wk_r)
          case(3); loc(ib_eps_force) = loc(ib_eps_force) + sum_prod_r(uz_r, wk_r)
        end select
      enddo
    endif

    ! Linear dissipation: -L(k) w^ with L = coeff*k^2/k2_max + coeff_h*(k^2/k2_max)^nexp
    do n = 1, nfields
      if(dcoef(n) == 0.d0 .and. dcoef_h(n) == 0.d0) cycle
      !$acc parallel loop collapse(3) present(wk, lnrho, mx, my, mz, bx, by, bz, kxt, ky, kz, &
      !$acc   dcoef, dcoef_h, dcoef_exp) private(k2r, lrate)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            k2r   = (kxt(i,j)**2 + ky(j)**2 + kz(k)**2)/k2_max
            lrate = dcoef(n)*k2r + dcoef_h(n)*k2r**dcoef_exp(n)
            if(n == ilnrho) then
              wk(k,j,i) = -lrate*lnrho(k,j,i)
            elseif(n == imx) then
              wk(k,j,i) = -lrate*mx   (k,j,i)
            elseif(n == imy) then
              wk(k,j,i) = -lrate*my   (k,j,i)
            elseif(n == imz) then
              wk(k,j,i) = -lrate*mz   (k,j,i)
            elseif(n == ibx) then
              wk(k,j,i) = -lrate*bx   (k,j,i)
            elseif(n == iby) then
              wk(k,j,i) = -lrate*by   (k,j,i)
            else
              wk(k,j,i) = -lrate*bz   (k,j,i)
            endif
          enddo
        enddo
      enddo
      call btran_wk
      select case(n)
        case(ilnrho); loc(ib_eps_lmd) = loc(ib_eps_lmd) + sum_prod_r(dEds_r, wk_r)
        case(imx   ); loc(ib_eps_nu ) = loc(ib_eps_nu ) + sum_prod_r(ux_r, wk_r)
        case(imy   ); loc(ib_eps_nu ) = loc(ib_eps_nu ) + sum_prod_r(uy_r, wk_r)
        case(imz   ); loc(ib_eps_nu ) = loc(ib_eps_nu ) + sum_prod_r(uz_r, wk_r)
        case(ibx   ); loc(ib_eps_res) = loc(ib_eps_res) + sum_prod_r(w_r(:,:,:,ibx), wk_r)
        case(iby   ); loc(ib_eps_res) = loc(ib_eps_res) + sum_prod_r(w_r(:,:,:,iby), wk_r)
        case(ibz   ); loc(ib_eps_res) = loc(ib_eps_res) + sum_prod_r(w_r(:,:,:,ibz), wk_r)
      end select
    enddo

    call sum_allreduce(loc, comm=comm_fft)
    loc = loc*dlx*dly*dlz

    ! already reduced in get_nonlinear_terms
    loc(ib_dMdt_raw) = dMdt_raw
    loc(ib_dMdt_abs) = dMdt_abs
    loc(ib_visc    ) = visc_diss
    loc(ib_R_zeta  ) = loc(ib_eps_shock) + visc_diss
    if(svv_es) loc(ib_R_svv) = loc(ib_eps_svv) + loc(ib_D_svv)
    if(C_eta > 0.d0) then
      loc(ib_eta_diss) = eta_diss
      loc(ib_R_eta   ) = loc(ib_eps_eta) + eta_diss
    endif

    if(istage == 1) then
      t_bgt  = tt
      dt_bgt = dt
    endif
    bgt = bgt + b_rk(istage)*loc
    if(visc_diss > 0.d0) R_zeta_ratio_max = max(R_zeta_ratio_max, abs(loc(ib_R_zeta))/visc_diss)
    if(loc(ib_D_svv) > 0.d0) R_svv_ratio_max = max(R_svv_ratio_max, abs(loc(ib_R_svv))/loc(ib_D_svv))
    if(svv_es .and. loc(ib_A_svv) > 0.d0) R_svv_scale_max = max(R_svv_scale_max, abs(loc(ib_R_svv))/loc(ib_A_svv))
    if(loc(ib_eta_diss) > 0.d0) R_eta_ratio_max = max(R_eta_ratio_max, abs(loc(ib_R_eta))/loc(ib_eta_diss))
    if(istage == 3) bgt_pending = .true.

  contains

    ! cuFFTMp takes timer_fft itself
    subroutine btran_wk
      call btran_c2r(wk, wk_r)
    end subroutine btran_wk

  end subroutine stage_budget


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   sum a*b over the physical grid points
!!          (k = 1, nlz) of this rank, as a device
!!          reduction (Calliope: sum(a*b))
!-----------------------------------------------!
  function sum_prod_r(a, b)
    use grid, only: nlx_local, nly, nlz
    implicit none
    real(8), dimension(:,:,:), intent(in) :: a, b
    real(8) :: sum_prod_r
    real(8) :: s
    integer :: i, j, k

    s = 0.d0
    !$acc parallel loop collapse(3) present(a, b) reduction(+:s)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          s = s + a(k,j,i)*b(k,j,i)
        enddo
      enddo
    enddo
    sum_prod_r = s

  end function sum_prod_r


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   sum |a*b| over the physical grid
!!          points (k = 1, nlz) of this rank, as a
!!          device reduction (Calliope:
!!          sum(abs(a*b)))
!-----------------------------------------------!
  function sum_abs_prod_r(a, b)
    use grid, only: nlx_local, nly, nlz
    implicit none
    real(8), dimension(:,:,:), intent(in) :: a, b
    real(8) :: sum_abs_prod_r
    real(8) :: s
    integer :: i, j, k

    s = 0.d0
    !$acc parallel loop collapse(3) present(a, b) reduction(+:s)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          s = s + abs(a(k,j,i)*b(k,j,i))
        enddo
      enddo
    enddo
    sum_abs_prod_r = s

  end function sum_abs_prod_r


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Write the budget of the step that has
!!          just finished. mon_np1: state at the
!!          end of the step (before the mass
!!          correction). E_n = mon_n.
!-----------------------------------------------!
  subroutine write_budget(mon_np1)
    use mp, only: proc0
    implicit none
    real(8), intent(in) :: mon_np1(nmon)
    real(8) :: E_n, E_np1, dMdt_act, dEdt_act, dEdt_sum, eps_ideal

    E_n   = sum(mon_n  (imon_ek:imon_ei))
    E_np1 = sum(mon_np1(imon_ek:imon_ei))
    dMdt_act  = (mon_np1(imon_M) - mon_n(imon_M))/dt_bgt
    dEdt_act  = (E_np1 - E_n)/dt_bgt
    dEdt_sum  = bgt(ib_eps_N) + bgt(ib_eps_res) + bgt(ib_eps_nu) + bgt(ib_eps_lmd)
    eps_ideal = bgt(ib_eps_N) - bgt(ib_eps_shock) - bgt(ib_eps_force) - bgt(ib_eps_svv) - bgt(ib_eps_svv_b) &
              - bgt(ib_eps_eta)

    if(proc0) then
      write (unit=budget_unit, fmt="(100es26.17e3)") t_bgt, dt_bgt, &
            dMdt_act, bgt(ib_dMdt_raw), bgt(ib_dMdt_proj), bgt(ib_dMdt_abs), &
            dEdt_act, dEdt_sum, dEdt_act - dEdt_sum, &
            eps_ideal, bgt(ib_eps_shock), -bgt(ib_visc), bgt(ib_R_zeta), R_zeta_ratio_max, &
            bgt(ib_eps_force), bgt(ib_eps_res), bgt(ib_eps_nu), bgt(ib_eps_lmd), E_n, E_np1, &
            bgt(ib_eps_svv) + bgt(ib_eps_svv_b), bgt(ib_eps_svv_b), -bgt(ib_D_svv), bgt(ib_R_svv), R_svv_ratio_max, &
            R_svv_scale_max, bgt(ib_eps_eta), -bgt(ib_eta_diss), bgt(ib_R_eta), R_eta_ratio_max
      call flush(budget_unit)
    endif

    bgt = 0.d0
    R_zeta_ratio_max = 0.d0
    R_svv_ratio_max  = 0.d0
    R_svv_scale_max  = 0.d0
    R_eta_ratio_max  = 0.d0
    bgt_pending = .false.

  end subroutine write_budget


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   Direction and number of points of the
!!          line of write_line: the first active
!!          direction (Gallope: shared with
!!          allocate_advance, which allocates lbuf)
!-----------------------------------------------!
  subroutine line_dir(dir, nl)
    use grid, only: nlx, nly, nlz
    use shock_viscosity, only: active_x, active_y
    implicit none
    character(len=1), intent(out) :: dir
    integer, intent(out) :: nl

    if(active_x) then
      dir = 'x'; nl = nlx
    elseif(active_y) then
      dir = 'y'; nl = nly
    else
      dir = 'z'; nl = nlz
    endif

  end subroutine line_dir


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   ASCII output of rho, u, B, zeta, D
!!          (and eta for C_eta > 0) on the line
!!          along the first active direction
!!          through grid index 1 of the other
!!          directions (out1d/line_NNNNN.dat)
!!
!!          Gallope: x slabs in real space. An x
!!          line (y = z = 1) has points on every
!!          rank; a y or z line (x = 1) lies on the
!!          rank with iproc_fft = 0. The line is
!!          built on the device in lbuf (allocated
!!          in allocate_advance) and copied to the
!!          host.
!-----------------------------------------------!
  subroutine write_line
    use grid, only: nlx_local, nly, nlz, xx, yy, zz
    use params, only: C_eta
    use mp, only: proc0, sum_reduce, iproc_fft
    use time, only: tt
    use file, only: get_unused_unit
    use shock_viscosity, only: abort_run
    implicit none
    character(len=1)   :: dir
    character(len=100) :: fn
    integer :: i, j, k, n, nl, nc, unit, ierr, ioff

    call line_dir(dir, nl)
    nc = merge(11, 10, C_eta > 0.d0)

    !$acc parallel loop collapse(2) present(lbuf)
    do n = 1, nc
      do i = 1, nl
        lbuf(i,n) = 0.d0
      enddo
    enddo

    select case(dir)
      case('x')
        ioff = nlx_local*iproc_fft
        !$acc parallel loop present(lbuf, xx, rho_r, ux_r, uy_r, uz_r, w_r, zeta_r, div_r)
        do i = 1, nlx_local
          lbuf(i+ioff, 1) = xx(i)
          lbuf(i+ioff, 2) = rho_r(1,1,i)
          lbuf(i+ioff, 3) = ux_r(1,1,i)
          lbuf(i+ioff, 4) = uy_r(1,1,i)
          lbuf(i+ioff, 5) = uz_r(1,1,i)
          lbuf(i+ioff, 6) = w_r(1,1,i,ibx)
          lbuf(i+ioff, 7) = w_r(1,1,i,iby)
          lbuf(i+ioff, 8) = w_r(1,1,i,ibz)
          lbuf(i+ioff, 9) = zeta_r(1,1,i)
          lbuf(i+ioff,10) = div_r(1,1,i)
        enddo
        if(nc == 11) then
          !$acc parallel loop present(lbuf, eta_r)
          do i = 1, nlx_local
            lbuf(i+ioff,11) = eta_r(1,1,i)
          enddo
        endif
      case('y')
        if(iproc_fft == 0) then
          !$acc parallel loop present(lbuf, yy, rho_r, ux_r, uy_r, uz_r, w_r, zeta_r, div_r)
          do j = 1, nly
            lbuf(j, 1) = yy(j)
            lbuf(j, 2) = rho_r(1,j,1)
            lbuf(j, 3) = ux_r(1,j,1)
            lbuf(j, 4) = uy_r(1,j,1)
            lbuf(j, 5) = uz_r(1,j,1)
            lbuf(j, 6) = w_r(1,j,1,ibx)
            lbuf(j, 7) = w_r(1,j,1,iby)
            lbuf(j, 8) = w_r(1,j,1,ibz)
            lbuf(j, 9) = zeta_r(1,j,1)
            lbuf(j,10) = div_r(1,j,1)
          enddo
          if(nc == 11) then
            !$acc parallel loop present(lbuf, eta_r)
            do j = 1, nly
              lbuf(j,11) = eta_r(1,j,1)
            enddo
          endif
        endif
      case('z')
        if(iproc_fft == 0) then
          !$acc parallel loop present(lbuf, zz, rho_r, ux_r, uy_r, uz_r, w_r, zeta_r, div_r)
          do k = 1, nlz
            lbuf(k, 1) = zz(k)
            lbuf(k, 2) = rho_r(k,1,1)
            lbuf(k, 3) = ux_r(k,1,1)
            lbuf(k, 4) = uy_r(k,1,1)
            lbuf(k, 5) = uz_r(k,1,1)
            lbuf(k, 6) = w_r(k,1,1,ibx)
            lbuf(k, 7) = w_r(k,1,1,iby)
            lbuf(k, 8) = w_r(k,1,1,ibz)
            lbuf(k, 9) = zeta_r(k,1,1)
            lbuf(k,10) = div_r(k,1,1)
          enddo
          if(nc == 11) then
            !$acc parallel loop present(lbuf, eta_r)
            do k = 1, nlz
              lbuf(k,11) = eta_r(k,1,1)
            enddo
          endif
        endif
    end select
    !$acc update self(lbuf)

    ! each point is owned by one rank; the other ranks add 0, which is exact
    call sum_reduce(lbuf, 0)

    if(proc0) then
      write(fn, '("out1d/line_", i5.5, ".dat")') nline1d
      call get_unused_unit(unit)
      open(unit=unit, file=trim(fn), status='replace', iostat=ierr)
      if(ierr /= 0) call abort_run('write_line: cannot open '//trim(fn)//' (create the directory out1d)')
      write(unit, '("# t = ", es26.17e3, ", dir = ", a1, ", step = ", i10)') tt, dir, istep_adv
      if(nc == 11) then
        write(unit, '(a)') '# 1:'//dir//' 2:rho 3:ux 4:uy 5:uz 6:bx 7:by 8:bz 9:zeta 10:D 11:eta'
      else
        write(unit, '(a)') '# 1:'//dir//' 2:rho 3:ux 4:uy 5:uz 6:bx 7:by 8:bz 9:zeta 10:D'
      endif
      do i = 1, nl
        write(unit, '(11es26.17e3)') lbuf(i,:)
      enddo
      close(unit)
    endif
    nline1d = nline1d + 1

  end subroutine write_line


!-----------------------------------------------!
!> @author  YK
!! @date    3 Nov 2026
!! @brief   Error of div B of the state at tt
!!            err(1) = rms(div_N B)/(k_N rms|B|)
!!            err(2) = max|div_N B|/(k_N max|B|)
!!          with div_N = F^-1 i k_d M F as in
!!          curl_N of the induction equation and
!!          k_N = sqrt(kd2_max), so err <= 1.
!!          B in real space must be in w_r.
!!          w(:,:,:,1) and flx_r(:,:,:,iflx_bx)
!!          are used as work arrays.
!-----------------------------------------------!
  subroutine get_divb_err(bx, by, bz, err)
    use grid, only: nlx_local, nly, nlz, nkx, nky_local, nkz
    use params, only: zi, filter_mag
    use dealias, only: filter
    use mp, only: sum_allreduce, max_allreduce, comm_fft
    use cufftmp, only: btran_c2r
    use shock_viscosity, only: ky_d, kz_d, kd2_max
    implicit none
    complex(8), dimension (:,:,:), intent(in) :: bx, by, bz
    real(8), intent(out) :: err(2)
    real(8) :: sums(2), maxs(2), bfilt, kn
    real(8) :: s1, s2, m1, m2 ! device reduction variables
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(w, bx, by, bz, filter, kxd, ky_d, kz_d) private(bfilt)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          bfilt = merge(filter(k,j,i), 1.d0, filter_mag)
          w(k,j,i,1) = zi*(kxd(i,j)*bx(k,j,i) + ky_d(j)*by(k,j,i) + kz_d(k)*bz(k,j,i))*bfilt
        enddo
      enddo
    enddo
    call btran_c2r(w(:,:,:,1), flx_r(:,:,:,iflx_bx))

    ! the maxima are of non-negative values, so 0 is a valid initial value
    s1 = 0.d0
    s2 = 0.d0
    m1 = 0.d0
    m2 = 0.d0
    !$acc parallel loop collapse(3) present(flx_r, w_r) reduction(+:s1, s2) reduction(max:m1, m2)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          s1 = s1 + flx_r(k,j,i,iflx_bx)**2
          s2 = s2 + (w_r(k,j,i,ibx)**2 + w_r(k,j,i,iby)**2 + w_r(k,j,i,ibz)**2)
          m1 = max(m1, abs(flx_r(k,j,i,iflx_bx)))
          m2 = max(m2, w_r(k,j,i,ibx)**2 + w_r(k,j,i,iby)**2 + w_r(k,j,i,ibz)**2)
        enddo
      enddo
    enddo
    sums(1) = s1
    sums(2) = s2
    maxs(1) = m1
    maxs(2) = sqrt(m2)
    call sum_allreduce(sums, comm=comm_fft)
    call max_allreduce(maxs, comm=comm_fft)

    kn = sqrt(max(kd2_max, tiny(1.d0)))
    err(1) = sqrt(sums(1)/max(sums(2), tiny(1.d0)))/kn
    err(2) = maxs(1)/max(maxs(2), tiny(1.d0))/kn

  end subroutine get_divb_err


!-----------------------------------------------!
!> @author  YK
!! @date    3 Nov 2026
!! @brief   Binary output of rho, ux, uy, uz, bx,
!!          by, bz, zeta, D, eta (0 for C_eta = 0)
!!          of the state at tt on the planes
!!            y = yy(1): out2d/plane_y0.dat,
!!                       arrays (nlx, nlz)
!!            z = zz(1): out2d/plane_z0.dat,
!!                       arrays (nlx, nly), nly > 2 only
!!          Each call appends one record of the 10
!!          arrays (real*8, x fastest, MPI-IO);
!!          record number, t and step are written
!!          to out2d/plane_time.dat.
!!
!!          Gallope: x slabs in real space, so every
!!          rank owns its x slab of both planes. The
!!          planes are built on the device in pbuf_y
!!          and pbuf_z (allocated in
!!          allocate_advance) and copied to the host.
!-----------------------------------------------!
  subroutine write_planes
    use grid, only: nlx, nly, nlz, nlx_local
    use params, only: C_eta
    use mp, only: proc0, iproc_fft
    use time, only: tt
    use file, only: get_unused_unit
    use shock_viscosity, only: active_y, abort_run
    implicit none
    integer :: i, j, k, unit, ierr

    ! y = yy(1)
    !$acc parallel loop collapse(2) present(pbuf_y, rho_r, ux_r, uy_r, uz_r, w_r, zeta_r, div_r)
    do k = 1, nlz
      do i = 1, nlx_local
        pbuf_y(i,k, 1) = rho_r(k,1,i)
        pbuf_y(i,k, 2) = ux_r(k,1,i)
        pbuf_y(i,k, 3) = uy_r(k,1,i)
        pbuf_y(i,k, 4) = uz_r(k,1,i)
        pbuf_y(i,k, 5) = w_r(k,1,i,ibx)
        pbuf_y(i,k, 6) = w_r(k,1,i,iby)
        pbuf_y(i,k, 7) = w_r(k,1,i,ibz)
        pbuf_y(i,k, 8) = zeta_r(k,1,i)
        pbuf_y(i,k, 9) = div_r(k,1,i)
        pbuf_y(i,k,10) = 0.d0
      enddo
    enddo
    if(C_eta > 0.d0) then
      !$acc parallel loop collapse(2) present(pbuf_y, eta_r)
      do k = 1, nlz
        do i = 1, nlx_local
          pbuf_y(i,k,10) = eta_r(k,1,i)
        enddo
      enddo
    endif
    !$acc update self(pbuf_y)
    call put_plane('out2d/plane_y0.dat', (/nlx, nlz/), (/nlx_local, nlz/), &
                   (/nlx_local*iproc_fft, 0/), .true., pbuf_y)

    ! z = zz(1)
    if(active_y) then
      !$acc parallel loop collapse(2) present(pbuf_z, rho_r, ux_r, uy_r, uz_r, w_r, zeta_r, div_r)
      do j = 1, nly
        do i = 1, nlx_local
          pbuf_z(i,j, 1) = rho_r(1,j,i)
          pbuf_z(i,j, 2) = ux_r(1,j,i)
          pbuf_z(i,j, 3) = uy_r(1,j,i)
          pbuf_z(i,j, 4) = uz_r(1,j,i)
          pbuf_z(i,j, 5) = w_r(1,j,i,ibx)
          pbuf_z(i,j, 6) = w_r(1,j,i,iby)
          pbuf_z(i,j, 7) = w_r(1,j,i,ibz)
          pbuf_z(i,j, 8) = zeta_r(1,j,i)
          pbuf_z(i,j, 9) = div_r(1,j,i)
          pbuf_z(i,j,10) = 0.d0
        enddo
      enddo
      if(C_eta > 0.d0) then
        !$acc parallel loop collapse(2) present(pbuf_z, eta_r)
        do j = 1, nly
          do i = 1, nlx_local
            pbuf_z(i,j,10) = eta_r(1,j,i)
          enddo
        enddo
      endif
      !$acc update self(pbuf_z)
      call put_plane('out2d/plane_z0.dat', (/nlx, nly/), (/nlx_local, nly/), &
                     (/nlx_local*iproc_fft, 0/), .true., pbuf_z)
    endif

    if(proc0) then
      call get_unused_unit(unit)
      if(nplanes == 0) then
        open(unit=unit, file='out2d/plane_time.dat', status='replace', iostat=ierr)
        if(ierr == 0) write(unit, '(a)') '# 1:record (from 0) 2:t 3:step; 10 arrays per record: ' &
                                       //'rho ux uy uz bx by bz zeta D eta'
      else
        open(unit=unit, file='out2d/plane_time.dat', status='old', position='append', iostat=ierr)
      endif
      if(ierr /= 0) call abort_run('write_planes: cannot open out2d/plane_time.dat')
      write(unit, '(i8, es26.17e3, i10)') nplanes, tt, istep_adv
      close(unit)
    endif
    nplanes = nplanes + 1

  end subroutine write_planes


!-----------------------------------------------!
!> @author  YK
!! @date    3 Nov 2026
!! @brief   Write the record nplanes of the plane
!!          file fn (collective). Ranks with
!!          own = F write nothing (count 0).
!!          The file is truncated at the first
!!          record of the run.
!-----------------------------------------------!
  subroutine put_plane(fn, sizes, subsizes, starts, own, buf)
    use MPI
    use shock_viscosity, only: abort_run
    implicit none
    character(len=*), intent(in) :: fn
    integer, intent(in) :: sizes(2), subsizes(2), starts(2)
    logical, intent(in) :: own
    real(8), intent(in) :: buf(:,:,:)
    integer :: fh, newtype, cnt, q, ierr
    integer(kind=MPI_OFFSET_KIND) :: disp

    call MPI_FILE_OPEN(MPI_COMM_WORLD, fn, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    if(ierr /= MPI_SUCCESS) call abort_run('write_planes: cannot open '//fn//' (create the directory out2d)')
    if(nplanes == 0) call MPI_FILE_SET_SIZE(fh, 0_MPI_OFFSET_KIND, ierr)

    ! the view of a rank without data is a valid 1 x 1 subarray; it writes 0 elements
    if(own) then
      call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_REAL8, newtype, ierr)
      cnt = subsizes(1)*subsizes(2)
    else
      call MPI_TYPE_CREATE_SUBARRAY(2, sizes, (/1, 1/), (/0, 0/), MPI_ORDER_FORTRAN, MPI_REAL8, newtype, ierr)
      cnt = 0
    endif
    call MPI_TYPE_COMMIT(newtype, ierr)
    do q = 1, size(buf, 3)
      disp = (int(nplanes, MPI_OFFSET_KIND)*size(buf, 3) + (q - 1))*sizes(1)*sizes(2)*8
      call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8, newtype, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(fh, buf(:,:,q), cnt, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
    enddo
    call MPI_TYPE_FREE(newtype, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

  end subroutine put_plane


!-----------------------------------------------!
!> @author  YK
!! @date    15 Jul 2021
!! @brief   Output series modes
!!
!!          Gallope: series_modes(n, 2) is the
!!          global ky index; the rank with that ky
!!          in its slab copies the mode into smodes
!!          on the device (allocated in
!!          allocate_advance). ky is labelled with
!!          ky_global.
!-----------------------------------------------!
  subroutine output_series_modes
    use fields, only: lnrho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use grid, only: kx, ky_global, kz
    use grid, only: nkx, nky_local, nkz
    use mp, only: proc0, sum_reduce, iproc_fft
    use time, only: tt
    use diagnostics_common, only: n_series_modes, series_modes
    use diagnostics_common, only: series_modes_unit
    implicit none
    integer :: n, i, j, k, jl, joff

    joff = nky_local*iproc_fft
    !$acc parallel loop present(smodes, series_modes, lnrho, mx, my, mz, bx, by, bz) private(i, j, k, jl)
    do n = 1, n_series_modes
      i  = series_modes(n, 1)
      j  = series_modes(n, 2)
      k  = series_modes(n, 3)
      jl = j - joff

      smodes(n, ilnrho) = (0.d0, 0.d0)
      smodes(n, imx   ) = (0.d0, 0.d0)
      smodes(n, imy   ) = (0.d0, 0.d0)
      smodes(n, imz   ) = (0.d0, 0.d0)
      smodes(n, ibx   ) = (0.d0, 0.d0)
      smodes(n, iby   ) = (0.d0, 0.d0)
      smodes(n, ibz   ) = (0.d0, 0.d0)

      if(       (i  >= 1 .and. i  <= nkx      ) &
          .and. (jl >= 1 .and. jl <= nky_local) &
          .and. (k  >= 1 .and. k  <= nkz      ) &
        ) then

        smodes(n, ilnrho) = lnrho(k, jl, i)
        smodes(n, imx   ) =    mx(k, jl, i)
        smodes(n, imy   ) =    my(k, jl, i)
        smodes(n, imz   ) =    mz(k, jl, i)

        smodes(n, ibx   ) =    bx(k, jl, i)
        smodes(n, iby   ) =    by(k, jl, i)
        smodes(n, ibz   ) =    bz(k, jl, i)

      endif
    enddo
    !$acc update self(smodes)

    call sum_reduce(smodes(:, ilnrho), 0)
    call sum_reduce(smodes(:, imx   ), 0)
    call sum_reduce(smodes(:, imy   ), 0)
    call sum_reduce(smodes(:, imz   ), 0)
    call sum_reduce(smodes(:, ibx   ), 0)
    call sum_reduce(smodes(:, iby   ), 0)
    call sum_reduce(smodes(:, ibz   ), 0)

    do n = 1, n_series_modes
      if(proc0) then
        i = series_modes(n, 1)
        j = series_modes(n, 2)
        k = series_modes(n, 3)
999 format(es30.21, A6, 5es30.21e3)
        write (unit=series_modes_unit, fmt=999) tt, 'lnrho', kx(i), ky_global(j), kz(k), smodes(n, ilnrho)
        write (unit=series_modes_unit, fmt=999) tt,    'mx', kx(i), ky_global(j), kz(k), smodes(n, imx   )
        write (unit=series_modes_unit, fmt=999) tt,    'my', kx(i), ky_global(j), kz(k), smodes(n, imy   )
        write (unit=series_modes_unit, fmt=999) tt,    'mz', kx(i), ky_global(j), kz(k), smodes(n, imz   )
        write (unit=series_modes_unit, fmt=999) tt,    'bx', kx(i), ky_global(j), kz(k), smodes(n, ibx   )
        write (unit=series_modes_unit, fmt=999) tt,    'by', kx(i), ky_global(j), kz(k), smodes(n, iby   )
        write (unit=series_modes_unit, fmt=999) tt,    'bz', kx(i), ky_global(j), kz(k), smodes(n, ibz   )
        call flush(series_modes_unit)

      endif
    enddo

  end subroutine output_series_modes


!-----------------------------------------------!
!> @author  YK
!! @date    3 Mar 2021
!! @brief   Remap
!!
!!          Gallope: shear = T stops the run in
!!          init_model_specific, so the remap of
!!          Calliope is not ported; this stub keeps
!!          the call in solve.
!-----------------------------------------------!
  subroutine remap
    use shock_viscosity, only: abort_run
    implicit none

    call abort_run('remap: shear = T is not supported in Gallope')

  end subroutine remap


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Enforce div free for (wx, wy, wz)
!!
!!          Gallope: nabla^-2 (div w) is a scalar
!!          of each mode (Calliope: a work array)
!-----------------------------------------------!
  subroutine div_free(wx, wy, wz)
    use grid, only: ky, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: zi
    use shearing_box, only: k2t_inv
    implicit none
    complex(8), dimension (:,:,:), intent(inout) :: wx, wy, wz
    complex(8) :: nbl2inv_div_w ! nabla^-2 (div w)
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(wx, wy, wz, kxt, ky, kz, k2t_inv) private(nbl2inv_div_w)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          nbl2inv_div_w = -zi*(   kxt(i,j)*wx(k,j,i) &
                                + ky(j)   *wy(k,j,i) &
                                + kz(k)   *wz(k,j,i) )*k2t_inv(k,j,i)

          wx(k,j,i) = wx(k,j,i) - zi*kxt(i,j)*nbl2inv_div_w
          wy(k,j,i) = wy(k,j,i) - zi*ky(j)   *nbl2inv_div_w
          wz(k,j,i) = wz(k,j,i) - zi*kz(k)   *nbl2inv_div_w
        enddo
      enddo
    enddo

  end subroutine div_free


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Allocates arrays for this module
!!
!!          Gallope: every array is also created on
!!          the device (enter data). The work arrays
!!          that Calliope allocates at every call of
!!          stage_budget, write_line, write_planes
!!          and output_series_modes are allocated
!!          here, when the output is on.
!-----------------------------------------------!
  subroutine allocate_advance
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use params, only: C_eta, C_shock, budget, write_intvl_1D, write_intvl_planes
    use shock_viscosity, only: svv_es, active_y
    use diagnostics_common, only: series_output, n_series_modes
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src
    real   (8), allocatable, dimension(:,:,:) :: src_r
    character(len=1) :: dir
    integer :: nl

    if (is_allocated) then
      return
    else
      is_allocated = .true.
    endif

    allocate(src  (nkz, nky_local, nkx), source=(0.d0,0.d0))
    allocate(src_r(nlz_padded, nly, nlx_local), source=0.d0)

    allocate(lnrho_new, source=src)
    allocate(   mx_new, source=src)
    allocate(   my_new, source=src)
    allocate(   mz_new, source=src)
    allocate(   bx_new, source=src)
    allocate(   by_new, source=src)
    allocate(   bz_new, source=src)

    allocate(w        (nkz       , nky_local, nkx      , nfields), source=(0.d0, 0.d0))
    allocate(flx      (nkz       , nky_local, nkx      , nftran ), source=(0.d0, 0.d0))
    allocate(w_r      (nlz_padded, nly      , nlx_local, nfields), source=0.d0)
    allocate(flx_r    (nlz_padded, nly      , nlx_local, nftran ), source=0.d0)
    allocate(exp_terms(nkz       , nky_local, nkx      , nfields), source=(0.d0, 0.d0))

    allocate(rho_r , source=src_r)
    allocate(ux_r  , source=src_r)
    allocate(uy_r  , source=src_r)
    allocate(uz_r  , source=src_r)
    allocate(div_r , source=src_r)
    allocate(zeta_r, source=src_r)
    allocate(q_r   , source=src_r)
    allocate(svv_r , source=src_r)
    if(C_eta > 0.d0) then
      allocate(eta_r, source=src_r)
      allocate(ej_r(nlz_padded, nly, nlx_local, 3), source=0.d0)
    endif
    if(svv_es) then
      allocate(svv_mh(nkz, nky_local, nkx, 3), source=(0.d0, 0.d0))
      allocate(svv_uh(nkz, nky_local, nkx, 3), source=(0.d0, 0.d0))
    endif

    !$acc enter data copyin(lnrho_new, mx_new, my_new, mz_new, bx_new, by_new, bz_new)
    !$acc enter data copyin(w, flx, w_r, flx_r, exp_terms)
    !$acc enter data copyin(rho_r, ux_r, uy_r, uz_r, div_r, zeta_r, q_r, svv_r)
    if(C_eta > 0.d0) then
      !$acc enter data copyin(eta_r, ej_r)
    endif
    if(svv_es) then
      !$acc enter data copyin(svv_mh, svv_uh)
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    allocate(lnrho_tmp, source=src)
    allocate(   mx_tmp, source=src)
    allocate(   my_tmp, source=src)
    allocate(   mz_tmp, source=src)
    allocate(   bx_tmp, source=src)
    allocate(   by_tmp, source=src)
    allocate(   bz_tmp, source=src)

    allocate(exp_terms0(nkz, nky_local, nkx, nfields)); exp_terms0  = 0.d0

    !$acc enter data copyin(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)
    !$acc enter data copyin(exp_terms0)

    ! Work arrays of stage_budget
    if(budget) then
      allocate(wk    , source=src)
      allocate(wk_r  , source=src_r)
      allocate(dEds_r, source=src_r)
      !$acc enter data copyin(wk, wk_r, dEds_r)
      if(C_shock > 0.d0) then
        allocate(qh, source=src)
        !$acc enter data copyin(qh)
      endif
      if(C_eta > 0.d0) then
        allocate(ejh(nkz, nky_local, nkx, 3), source=(0.d0, 0.d0))
        !$acc enter data copyin(ejh)
      endif
    endif

    ! Line buffer of write_line
    if(write_intvl_1D > 0.d0) then
      call line_dir(dir, nl)
      allocate(lbuf(nl, merge(11, 10, C_eta > 0.d0)), source=0.d0)
      !$acc enter data copyin(lbuf)
    endif

    ! Plane buffers of write_planes
    if(write_intvl_planes > 0.d0) then
      allocate(pbuf_y(nlx_local, nlz, nq_plane), source=0.d0)
      !$acc enter data copyin(pbuf_y)
      if(active_y) then
        allocate(pbuf_z(nlx_local, nly, nq_plane), source=0.d0)
        !$acc enter data copyin(pbuf_z)
      endif
    endif

    ! Mode buffer of output_series_modes
    if(series_output) then
      allocate(smodes(n_series_modes, nfields), source=(0.d0, 0.d0))
      !$acc enter data copyin(smodes)
    endif

    deallocate(src, src_r)

  end subroutine allocate_advance


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Deallocates arrays for this module
!!
!!          Gallope: kxt and kxd (allocated in
!!          init_work_fields) are kept, as in
!!          Calliope.
!-----------------------------------------------!
  subroutine deallocate_advance
    use fields, only: lnrho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use time, only: nstep
    use terminate, only: terminated
    implicit none

    if (is_allocated) then
      is_allocated = .false.
    else
      return
    endif

    ! Monitor and line output of the final state when check_terminate stopped the run before nstep
    ! (solve does it at istep_adv = nstep); this also writes the budget row of the last step
    if(terminated .and. istep_adv < nstep) call get_nonlinear_terms(lnrho, mx, my, mz, bx, by, bz, 0)

    !$acc exit data delete(lnrho_new, mx_new, my_new, mz_new, bx_new, by_new, bz_new)
    !$acc exit data delete(w, flx, w_r, flx_r, exp_terms)
    !$acc exit data delete(rho_r, ux_r, uy_r, uz_r, div_r, zeta_r, q_r, svv_r)

    deallocate(lnrho_new)
    deallocate(   mx_new)
    deallocate(   my_new)
    deallocate(   mz_new)
    deallocate(   bx_new)
    deallocate(   by_new)
    deallocate(   bz_new)

    deallocate(w        )
    deallocate(flx      )
    deallocate(w_r      )
    deallocate(flx_r    )
    deallocate(exp_terms)

    deallocate(rho_r )
    deallocate(ux_r  )
    deallocate(uy_r  )
    deallocate(uz_r  )
    deallocate(div_r )
    deallocate(zeta_r)
    deallocate(q_r   )
    deallocate(svv_r )
    if(allocated(eta_r )) then
      !$acc exit data delete(eta_r)
      deallocate(eta_r )
    endif
    if(allocated(ej_r  )) then
      !$acc exit data delete(ej_r)
      deallocate(ej_r  )
    endif
    if(allocated(svv_mh)) then
      !$acc exit data delete(svv_mh)
      deallocate(svv_mh)
    endif
    if(allocated(svv_uh)) then
      !$acc exit data delete(svv_uh)
      deallocate(svv_uh)
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !$acc exit data delete(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)
    !$acc exit data delete(exp_terms0)

    deallocate(lnrho_tmp)
    deallocate(   mx_tmp)
    deallocate(   my_tmp)
    deallocate(   mz_tmp)
    deallocate(   bx_tmp)
    deallocate(   by_tmp)
    deallocate(   bz_tmp)

    deallocate(exp_terms0)

    ! Work arrays of stage_budget, write_line, write_planes and output_series_modes
    if(allocated(wk    )) then
      !$acc exit data delete(wk)
      deallocate(wk    )
    endif
    if(allocated(wk_r  )) then
      !$acc exit data delete(wk_r)
      deallocate(wk_r  )
    endif
    if(allocated(dEds_r)) then
      !$acc exit data delete(dEds_r)
      deallocate(dEds_r)
    endif
    if(allocated(qh    )) then
      !$acc exit data delete(qh)
      deallocate(qh    )
    endif
    if(allocated(ejh   )) then
      !$acc exit data delete(ejh)
      deallocate(ejh   )
    endif
    if(allocated(lbuf  )) then
      !$acc exit data delete(lbuf)
      deallocate(lbuf  )
    endif
    if(allocated(pbuf_y)) then
      !$acc exit data delete(pbuf_y)
      deallocate(pbuf_y)
    endif
    if(allocated(pbuf_z)) then
      !$acc exit data delete(pbuf_z)
      deallocate(pbuf_z)
    endif
    if(allocated(smodes)) then
      !$acc exit data delete(smodes)
      deallocate(smodes)
    endif

  end subroutine deallocate_advance

end module advance
