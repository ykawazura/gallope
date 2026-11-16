!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Local shock viscosity for MHD_COMP_ISOTH
!!
!!          zeta = C_shock*Delta^2*S[max(0, -D)],  D = div_N u
!!
!!          S = (3-point local max)^shock_nmax followed by
!!              ([1,2,1]/4 smoothing)^shock_nsmooth.
!!          Each pass is applied direction by direction in the fixed
!!          order y, z, x, only in the active directions (nl > 2).
!!          Both operations have non-negative weights, so zeta >= 0,
!!          and zeta = 0 farther than shock_nmax + shock_nsmooth cells
!!          from compression.
!!
!!          y and z are local in real space (periodic wrap inside the
!!          rank); the x halo (width 1) is exchanged with halo_exchange_x.
!!          At every grid point each pass uses the same operands in the
!!          same order, irrespective of the MPI decomposition, so zeta
!!          for a given D does not depend on the decomposition (bitwise).
!!
!!          Derivative: div_N = F^-1 i k_d M F, where k_d is set to 0 at
!!          the Nyquist mode (skew-adjointness) and M is the dealiasing
!!          filter applied to the momentum (identity if filter_mom = F).
!!
!!          Spectral vanishing viscosity (Tadmor 1989; Maday, Ould Kaber
!!          & Tadmor 1993) is also set up here. Its kernel and coefficients are
!!            Q = 0                                for |k_d| <= m_N,
!!                exp(-((|k_d| - k_c)/(|k_d| - m_N))^2) for m_N < |k_d| < k_c,
!!                1 (svv_above_kc = T) or 0 (F)    for |k_d| >= k_c,
!!            eps_d = svv_eps0*dl_d*max(|u_d| + cf_d)   (active d only),
!!          with k_c,d = 2/3 of the Nyquist wavenumber (the edge of the
!!          dealiasing filter) and m_N = svv_mN*k_c,d. Here k_d is the
!!          wavenumber itself (the Nyquist mode is damped).
!!          B: dB^/dt = -lambda(k) B^, lambda(k) = sum_d eps_d k_d^2 Q(|k_d|).
!!          rho and m, svv_form = 'linear': the same -lambda(k) on m^ and,
!!            in the mass-conserving form, on rho^ (the energy rate is not
!!            sign-definite).
!!          rho and m, svv_form = 'entropy': flux form in the entropy
!!            variables w = (cs^2 s - u^2/2, u) (get_svv_es),
!!              A_L,d = F^-1 sig_d M F,  A_H,d = F^-1 sig_d (1 - M) F,
!!              sig_d = i k_d sqrt(Q) (|k_d| sqrt(Q) at Nyquist, so that A
!!                      maps a real field to a real field),
!!              g_r = A_L,d w_r, g_m = A_L,d u, g_h = A_H,d w_r,
!!              f_r = a t, t = g_r + u.g_m, f_m = u f_r + b g_m, f_h = kap rho g_h,
!!              a = eps_d rho/cs^2, b = eps_d rho, kap = eps_d/(cs^2 + max|u|^2),
!!              d rho/dt = -sum_d (A_L,d^* f_r + A_H,d^* f_h),
!!              dm/dt    = -sum_d  A_L,d^* f_m,
!!            where M is the filter of m (identity if filter_mom = F; then
!!            A_H = 0; M is 0 or 1, so the filter of m after each stage
!!            leaves dm/dt unchanged). Then sum d rho/dt = sum dm/dt = 0 (sig_d = 0 at
!!            k_d = 0) and dE/dt = -D_svv with
!!              D_svv = sum_d sum [a t^2 + b |g_m|^2 + kap rho g_h^2] dV >= 0
!!            for any state, if the s tendency is not filtered
!!            (filter_lnrho = F). The explicit step is limited by the
!!            largest |eigenvalue| of the Jacobian of this term, estimated by
!!            the power method (svv_power).
!!          svv_mass = F (test only) drops the mass fluxes, f_r = f_h = 0, so
!!            that d rho/dt has no SVV term and f_m = b g_m. Then
!!              D_svv = sum_d sum b |g_m|^2 dV >= 0,
!!            mass and momentum are still conserved, but ln rho is not
!!            dissipated by the SVV.
!!
!!          Gallope (OpenACC + cuFFTMp): real-space arrays are
!!          (nlz_padded, nly, nlx_local) and spectral arrays (nkz, nky_local, nkx);
!!          real-space loops skip the FFT padding k = nlz+1..nlz_padded.
!!          Array arguments are device-present on entry and their results
!!          are valid on the device on exit. Every work array is allocated
!!          once in init_shock_viscosity (the SVV ones only for
!!          svv_form = 'entropy') and freed in finish_shock_viscosity.
!-----------------------------------------------!
module shock_viscosity
  implicit none
  private
  include 'mpif.h'

  public init_shock_viscosity, finish_shock_viscosity
  public get_div_u, get_zeta
  public kx_d, ky_d, kz_d, kd2_max
  public delta_shock, active_x, active_y, active_z
  public is_nonfinite, abort_run
  public svv_on, kq2x, kq2y, kq2z, svv_lam_m, svv_lam_s, svv_lam_b, get_svv_eps
  public svv_es, get_svv_es, svv_power, svv_lam, svv_seed_nomass
  public tmr_on, tmr_off, report_timers
  public itmr_step, itmr_zeta, itmr_svv_es, itmr_svv_pi, itmr_budget, itmr_planes

  logical  :: initialized = .false.
  real(8), allocatable :: kx_d(:), ky_d(:), kz_d(:) ! wavenumbers for derivatives (0 at Nyquist)
  real(8) :: kd2_max     ! sum_d (max |k_d| over the modes kept by M)^2; used for the viscous dt
  real(8) :: delta_shock ! Delta in zeta = C_shock*Delta^2*S[...]
  logical  :: active_x, active_y, active_z
  !$acc declare create(kx_d, ky_d, kz_d, kd2_max, delta_shock, active_x, active_y, active_z)

  ! spectral vanishing viscosity
  logical  :: svv_on = .false.
  real(8), allocatable :: kq2x(:), kq2y(:), kq2z(:) ! k_d^2 Q(|k_d|); 0 in inactive directions
  real(8) :: svv_lam_m(3), svv_lam_s(3), svv_lam_b(3) ! max_k k_d^2 Q over the modes kept in m, ln rho and B
  logical  :: svv_es = .false.                        ! svv_form = 'entropy'
  complex(8), allocatable :: sgx(:), sgy(:), sgz(:)   ! sig_d = i k_d sqrt(Q) (|k_d| sqrt(Q) at Nyquist); 0 if inactive
  real(8) :: svv_lam = 0.d0                           ! power-method estimate of the largest |eigenvalue| (svv_power)
  !$acc declare create(svv_on, kq2x, kq2y, kq2z, svv_lam_m, svv_lam_s, svv_lam_b, svv_es, sgx, sgy, sgz, svv_lam)
  real(8), allocatable :: pi_v(:,:,:,:)               ! power-method vector (s, mx, my, mz) in real space
  logical  :: pi_set = .false.                        ! pi_v holds a vector (Calliope: allocated(pi_v))

  ! work arrays (allocated once in alloc_work; enter data)
  complex(8), allocatable :: dv_u(:,:,:,:), dv_div(:,:,:)      ! get_div_u: F u (not normalized), F D
  real(8)   , allocatable :: zw(:,:,:)                         ! get_zeta: output of a pass (ping-pong with zeta_r)
  real(8)   , allocatable :: halo_lo(:,:), halo_hi(:,:)        ! pass_1d: x planes below/above the slab
  complex(8), allocatable :: es_uh(:,:,:,:)                    ! get_svv_es: F u/(nlx nly nlz) if uh is absent
  real(8)   , allocatable :: es_gr(:,:,:), es_gh(:,:,:), es_gm(:,:,:,:) ! svv_es_core: g_r, g_h, g_m (then fluxes)
  complex(8), allocatable :: es_wh(:,:,:), es_fh(:,:,:), es_rdh(:,:,:)  ! svv_es_core: F w_r, work, F d rho/dt
  real(8)   , allocatable :: G_rho(:,:,:), G_ux(:,:,:), G_uy(:,:,:), G_uz(:,:,:) ! svv_G: rho, u
  complex(8), allocatable :: G_mdh(:,:,:,:)                    ! svv_G: dm/dt (Fourier)
  real(8)   , allocatable :: pw_G0s(:,:,:), pw_Gs(:,:,:), pw_ps(:,:,:)       ! svv_power: G(U), G(U + h v), U + h v (s)
  real(8)   , allocatable :: pw_G0m(:,:,:,:), pw_Gm(:,:,:,:), pw_pm(:,:,:,:) ! svv_power: the same for m

  ! wall-clock time [s] of parts of the step, summed over the run on each rank (see report_timers)
  integer, parameter :: ntmr = 6
  integer, parameter :: itmr_step = 1, itmr_zeta = 2, itmr_svv_es = 3, itmr_svv_pi = 4, itmr_budget = 5, itmr_planes = 6
  character(len=28), parameter :: tmr_name(ntmr) = [character(len=28) :: &
    'step (solve)', 'div u, zeta, eta (halo)', 'get_svv_es', 'svv_power', 'stage_budget', 'write_planes']
  real(8) :: tmr(ntmr) = 0.d0, tmr_t0(ntmr) = 0.d0

contains


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Initialization: derivative wavenumbers,
!!          Delta, SVV kernel, work arrays. Must be
!!          called after init_filter.
!-----------------------------------------------!
  subroutine init_shock_viscosity
    use grid, only: nlx, nly, nlz, nkx, nky_local, nkz, dlx, dly, dlz
    use grid, only: kx, ky, kz, ikx, iky, ikz
    use grid, only: kx_max, ky_max, kz_max
    use dealias, only: filter
    use mp, only: proc0, max_allreduce, comm_fft
    use params, only: C_shock, shock_dl, shock_nmax, shock_nsmooth, filter_mom, shock_selftest, C_eta
    use params, only: svv_eps0, svv_mN, svv_above_kc, filter_lnrho, filter_mag, svv_form, cfl_svv, svv_mass
    implicit none
    real(8) :: kmx2, kmy2, kmz2, mfilt, sfilt, bfilt
    integer :: i, j, k

    if(initialized) return
    initialized = .true.

    ! wavenumbers for derivatives: Nyquist mode set to 0
    allocate(kx_d, source=kx)
    allocate(ky_d, source=ky)
    allocate(kz_d, source=kz)
    do i = 1, size(kx)
      if(abs(ikx(i)) == nlx/2) kx_d(i) = 0.d0
    enddo
    do j = 1, size(ky)
      if(abs(iky(j)) == nly/2) ky_d(j) = 0.d0
    enddo
    do k = 1, size(kz)
      if(abs(ikz(k)) == nlz/2) kz_d(k) = 0.d0
    enddo

    ! largest |k_d| among the modes kept in the momentum
    kmx2 = 0.d0; kmy2 = 0.d0; kmz2 = 0.d0
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          mfilt = merge(filter(k,j,i), 1.d0, filter_mom)
          if(mfilt > 0.d0) then
            kmx2 = max(kmx2, kx_d(i)**2)
            kmy2 = max(kmy2, ky_d(j)**2)
            kmz2 = max(kmz2, kz_d(k)**2)
          endif
        enddo
      enddo
    enddo
    call max_allreduce(kmx2, comm=comm_fft)
    call max_allreduce(kmy2, comm=comm_fft)
    call max_allreduce(kmz2, comm=comm_fft)
    kd2_max = kmx2 + kmy2 + kmz2

    ! active directions
    active_x = nlx > 2
    active_y = nly > 2
    active_z = nlz > 2

    ! Delta
    if(shock_dl > 0.d0) then
      delta_shock = shock_dl
    else
      delta_shock = huge(1.d0)
      if(active_x) delta_shock = min(delta_shock, dlx)
      if(active_y) delta_shock = min(delta_shock, dly)
      if(active_z) delta_shock = min(delta_shock, dlz)
      if(.not. (active_x .or. active_y .or. active_z)) delta_shock = 0.d0
    endif

    ! spectral vanishing viscosity: k_d^2 Q per direction and its maxima over the kept modes
    svv_on = svv_eps0 > 0.d0
    allocate(kq2x, mold=kx)
    allocate(kq2y, mold=ky)
    allocate(kq2z, mold=kz)
    call svv_kernel(kx, kx_max, active_x, kq2x)
    call svv_kernel(ky, ky_max, active_y, kq2y)
    call svv_kernel(kz, kz_max, active_z, kq2z)
    svv_es = svv_on .and. trim(svv_form) == 'entropy'
    allocate(sgx(size(kx)), sgy(size(ky)), sgz(size(kz)))
    call svv_symbol(kx, ikx, nlx, kq2x, sgx)
    call svv_symbol(ky, iky, nly, kq2y, sgy)
    call svv_symbol(kz, ikz, nlz, kq2z, sgz) ! kz is the r2c axis: its Nyquist mode is k = nkz
    svv_lam_m = 0.d0
    svv_lam_s = 0.d0
    svv_lam_b = 0.d0
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          mfilt = merge(filter(k,j,i), 1.d0, filter_mom)
          sfilt = merge(filter(k,j,i), 1.d0, filter_lnrho)
          bfilt = merge(filter(k,j,i), 1.d0, filter_mag)
          if(mfilt > 0.d0) svv_lam_m = max(svv_lam_m, (/kq2x(i), kq2y(j), kq2z(k)/))
          if(sfilt > 0.d0) svv_lam_s = max(svv_lam_s, (/kq2x(i), kq2y(j), kq2z(k)/))
          if(bfilt > 0.d0) svv_lam_b = max(svv_lam_b, (/kq2x(i), kq2y(j), kq2z(k)/))
        enddo
      enddo
    enddo
    call max_allreduce(svv_lam_m, comm=comm_fft)
    call max_allreduce(svv_lam_s, comm=comm_fft)
    call max_allreduce(svv_lam_b, comm=comm_fft)
    svv_lam = 0.d0

    !$acc update device(kx_d, ky_d, kz_d, kd2_max, delta_shock, active_x, active_y, active_z)
    !$acc update device(svv_on, kq2x, kq2y, kq2z, svv_lam_m, svv_lam_s, svv_lam_b, svv_es, sgx, sgy, sgz, svv_lam)

    if(proc0) then
      print *
      print '("  shock viscosity: C_shock = ", es12.4, ", Delta = ", es12.4, &
             &", n_max = ", i3, ", n_smooth = ", i3)', C_shock, delta_shock, shock_nmax, shock_nsmooth
      print '("                   active (x, y, z) = ", 3l2, ", sum_d k_d,max^2 = ", es12.4)', &
             active_x, active_y, active_z, kd2_max
      if(C_eta > 0.d0) print '("  shock resistivity: C_eta = ", es12.4, " (eta = C_eta*Delta^2*S[max(0, -D)])")', C_eta
      if(svv_on) then
        print '("  spectral vanishing viscosity: eps0 = ", es12.4, ", m_N/k_c = ", es12.4, &
               &", Q(|k| >= k_c) = ", i2)', svv_eps0, svv_mN, merge(1, 0, svv_above_kc)
        print '("                   max k_d^2 Q (x, y, z) in m: ", 3es12.4, ", in ln rho: ", 3es12.4)', &
               svv_lam_m, svv_lam_s
        if(svv_es) then
          print '("                   rho and m: entropy-stable flux form, cfl_svv = ", es12.4)', cfl_svv
          if(.not. svv_mass) &
            print '("                   svv_mass = F: no mass fluxes (no SVV in ln rho, f_m = eps rho A_L u)")'
        else
          print '("                   rho and m: linear form, -lambda(k)")'
        endif
      else
        print '("  spectral vanishing viscosity: off")'
      endif
      print *
    endif

    call alloc_work

    if(shock_selftest) call selftest

  end subroutine init_shock_viscosity


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   Allocate the work arrays once (host and
!!          device); the SVV ones only if svv_es
!-----------------------------------------------!
  subroutine alloc_work
    use grid, only: nlx_local, nly, nlz_padded, nkx, nky_local, nkz
    implicit none

    allocate(dv_u   (nkz, nky_local, nkx, 3), source=(0.d0, 0.d0))
    allocate(dv_div (nkz, nky_local, nkx)   , source=(0.d0, 0.d0))
    allocate(zw     (nlz_padded, nly, nlx_local), source=0.d0)
    allocate(halo_lo(nlz_padded, nly), source=0.d0)
    allocate(halo_hi(nlz_padded, nly), source=0.d0)
    !$acc enter data copyin(dv_u, dv_div, zw, halo_lo, halo_hi)

    if(.not. svv_es) return

    allocate(es_uh (nkz, nky_local, nkx, 3), source=(0.d0, 0.d0))
    allocate(es_gr (nlz_padded, nly, nlx_local), source=0.d0)
    allocate(es_gh (nlz_padded, nly, nlx_local), source=0.d0)
    allocate(es_gm (nlz_padded, nly, nlx_local, 3), source=0.d0)
    allocate(es_wh (nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(es_fh (nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(es_rdh(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    !$acc enter data copyin(es_uh, es_gr, es_gh, es_gm, es_wh, es_fh, es_rdh)

    allocate(G_rho(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(G_ux (nlz_padded, nly, nlx_local), source=0.d0)
    allocate(G_uy (nlz_padded, nly, nlx_local), source=0.d0)
    allocate(G_uz (nlz_padded, nly, nlx_local), source=0.d0)
    allocate(G_mdh(nkz, nky_local, nkx, 3), source=(0.d0, 0.d0))
    !$acc enter data copyin(G_rho, G_ux, G_uy, G_uz, G_mdh)

    allocate(pw_G0s(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(pw_Gs (nlz_padded, nly, nlx_local), source=0.d0)
    allocate(pw_ps (nlz_padded, nly, nlx_local), source=0.d0)
    allocate(pw_G0m(nlz_padded, nly, nlx_local, 3), source=0.d0)
    allocate(pw_Gm (nlz_padded, nly, nlx_local, 3), source=0.d0)
    allocate(pw_pm (nlz_padded, nly, nlx_local, 3), source=0.d0)
    allocate(pi_v  (nlz_padded, nly, nlx_local, 4), source=0.d0)
    !$acc enter data copyin(pw_G0s, pw_Gs, pw_ps, pw_G0m, pw_Gm, pw_pm, pi_v)
    pi_set = .false.

  end subroutine alloc_work


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   Finalization: free the work arrays and
!!          the module arrays. Call once at the end
!!          of the run (after report_timers).
!-----------------------------------------------!
  subroutine finish_shock_viscosity
    implicit none

    if(.not. initialized) return

    !$acc exit data delete(dv_u, dv_div, zw, halo_lo, halo_hi)
    deallocate(dv_u, dv_div, zw, halo_lo, halo_hi)

    if(svv_es) then
      !$acc exit data delete(es_uh, es_gr, es_gh, es_gm, es_wh, es_fh, es_rdh)
      deallocate(es_uh, es_gr, es_gh, es_gm, es_wh, es_fh, es_rdh)
      !$acc exit data delete(G_rho, G_ux, G_uy, G_uz, G_mdh)
      deallocate(G_rho, G_ux, G_uy, G_uz, G_mdh)
      !$acc exit data delete(pw_G0s, pw_Gs, pw_ps, pw_G0m, pw_Gm, pw_pm, pi_v)
      deallocate(pw_G0s, pw_Gs, pw_ps, pw_G0m, pw_Gm, pw_pm, pi_v)
      pi_set = .false.
    endif

    deallocate(kx_d, ky_d, kz_d, kq2x, kq2y, kq2z, sgx, sgy, sgz)
    initialized = .false.

  end subroutine finish_shock_viscosity


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   kq2 = k^2 Q(|k|) of the SVV kernel for
!!          one direction (0 if the direction is
!!          inactive); kmax_noprune = Nyquist k
!-----------------------------------------------!
  subroutine svv_kernel(k, kmax_noprune, active, kq2)
    use params, only: svv_mN, svv_above_kc
    implicit none
    real(8), intent(in ) :: k(:), kmax_noprune
    logical, intent(in ) :: active
    real(8), intent(out) :: kq2(:)
    real(8) :: kc, mN, a, r, Q
    integer :: i

    kc = kmax_noprune*2.d0/3.d0 ! same threshold as init_filter
    mN = svv_mN*kc
    do i = 1, size(k)
      a = abs(k(i))
      if(a >= kc) then
        Q = merge(1.d0, 0.d0, svv_above_kc)
      elseif(a > mN) then
        r = (kc - a)/(a - mN)
        Q = merge(exp(-r**2), 0.d0, r < 30.d0) ! no underflow next to m_N
      else
        Q = 0.d0
      endif
      kq2(i) = merge(k(i)**2*Q, 0.d0, active)
    enddo

  end subroutine svv_kernel


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   sig = i k sqrt(Q) of the entropy-stable
!!          SVV for one direction, from kq2 = k^2 Q
!!          (real |k| sqrt(Q) at Nyquist; 0 if the
!!          direction is inactive, where kq2 = 0)
!-----------------------------------------------!
  subroutine svv_symbol(k, ik, nl, kq2, sg)
    use params, only: zi
    implicit none
    real(8)   , intent(in ) :: k(:), kq2(:)
    integer   , intent(in ) :: ik(:), nl
    complex(8), intent(out) :: sg(:)
    integer :: i

    do i = 1, size(k)
      if(abs(ik(i)) == nl/2) then
        sg(i) = sqrt(kq2(i))
      else
        sg(i) = zi*sign(sqrt(kq2(i)), k(i))
      endif
    enddo

  end subroutine svv_symbol


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   SVV coefficients eps_d = svv_eps0*dl_d
!!          *vmax_d (active d), vmax_d = max(|u_d| + cf_d)
!-----------------------------------------------!
  subroutine get_svv_eps(vmax, eps)
    use grid, only: dlx, dly, dlz
    use params, only: svv_eps0
    implicit none
    real(8), intent(in ) :: vmax(3)
    real(8), intent(out) :: eps(3)

    eps = 0.d0
    if(.not. svv_on) return
    if(active_x) eps(1) = svv_eps0*dlx*vmax(1)
    if(active_y) eps(2) = svv_eps0*dly*vmax(2)
    if(active_z) eps(3) = svv_eps0*dlz*vmax(3)

  end subroutine get_svv_eps


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   D = div_N u = F^-1 i k_d M F u
!!          kxd(ikx, iky) is the x-wavenumber
!!          (sheared when shear = T).
!!          uh (optional): F u/(nlx nly nlz)
!-----------------------------------------------!
  subroutine get_div_u(ux_r, uy_r, uz_r, kxd, div_r, uh)
    use grid, only: nlx, nly, nlz, nkx, nky_local, nkz
    use dealias, only: filter
    use params, only: zi, filter_mom
    use cufftmp, only: ftran_r2c, btran_c2r
    implicit none
    real(8), dimension(:,:,:), intent(in ) :: ux_r, uy_r, uz_r
    real(8), dimension(:,:)  , intent(in ) :: kxd
    real(8), dimension(:,:,:), intent(out) :: div_r
    complex(8), dimension(:,:,:,:), intent(out), optional :: uh
    real(8) :: mfilt
    logical :: lmom
    integer :: i, j, k

    lmom = filter_mom

    call ftran_r2c(ux_r, dv_u(:,:,:,1))
    call ftran_r2c(uy_r, dv_u(:,:,:,2))
    call ftran_r2c(uz_r, dv_u(:,:,:,3))

    if(present(uh)) then
      !$acc data present(uh, dv_u)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            uh(k,j,i,1) = dv_u(k,j,i,1)/nlx/nly/nlz
            uh(k,j,i,2) = dv_u(k,j,i,2)/nlx/nly/nlz
            uh(k,j,i,3) = dv_u(k,j,i,3)/nlx/nly/nlz
          enddo
        enddo
      enddo
      !$acc end data
    endif

    !$acc data present(kxd, dv_u, dv_div, filter)
    !$acc parallel loop collapse(3) private(mfilt)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if(lmom) then
            mfilt = filter(k,j,i)
          else
            mfilt = 1.d0
          endif
          dv_div(k,j,i) = zi*(kxd(i,j)*dv_u(k,j,i,1) + ky_d(j)*dv_u(k,j,i,2) + kz_d(k)*dv_u(k,j,i,3)) &
                          *mfilt/nlx/nly/nlz
        enddo
      enddo
    enddo
    !$acc end data

    call btran_c2r(dv_div, div_r)

  end subroutine get_div_u


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Entropy-stable SVV of rho and m for the
!!          coefficients eps (module header):
!!          rdot_r = d rho/dt (real space),
!!          mdh    = dm/dt (Fourier, normalized as
!!                   the state; filtered by M),
!!          dsum   = local sum of a t^2 + b |g_m|^2
!!                   + kap rho g_h^2 (>= 0) over the
!!                   grid points and the directions
!!                   (this rank only; the caller reduces it).
!!          uh (optional): F u/(nlx nly nlz); it is
!!          computed here if absent.
!!          Requires svv_es (work arrays).
!-----------------------------------------------!
  subroutine get_svv_es(s_r, rho_r, ux_r, uy_r, uz_r, eps, rdot_r, mdh, dsum, uh)
    use grid, only: nlx, nly, nlz, nkx, nky_local, nkz
    use cufftmp, only: ftran_r2c
    implicit none
    real(8), dimension(:,:,:), intent(in ) :: s_r, rho_r, ux_r, uy_r, uz_r
    real(8), intent(in ) :: eps(3)
    real(8), dimension(:,:,:), intent(out) :: rdot_r
    complex(8), dimension(:,:,:,:), intent(out) :: mdh
    real(8), intent(out) :: dsum
    complex(8), dimension(:,:,:,:), intent(in), optional :: uh
    integer :: i, j, k

    if(.not. svv_es) call abort_run('get_svv_es: called without svv_form = entropy (no work arrays)')

    if(present(uh)) then
      call svv_es_core(s_r, rho_r, ux_r, uy_r, uz_r, uh, eps, rdot_r, mdh, dsum)
    else
      call ftran_r2c(ux_r, es_uh(:,:,:,1))
      call ftran_r2c(uy_r, es_uh(:,:,:,2))
      call ftran_r2c(uz_r, es_uh(:,:,:,3))
      !$acc data present(es_uh)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            es_uh(k,j,i,1) = es_uh(k,j,i,1)/nlx/nly/nlz
            es_uh(k,j,i,2) = es_uh(k,j,i,2)/nlx/nly/nlz
            es_uh(k,j,i,3) = es_uh(k,j,i,3)/nlx/nly/nlz
          enddo
        enddo
      enddo
      !$acc end data
      call svv_es_core(s_r, rho_r, ux_r, uy_r, uz_r, es_uh, eps, rdot_r, mdh, dsum)
    endif

  end subroutine get_svv_es


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Body of get_svv_es. Per active
!!          direction: 4 (5 with A_H) inverse and
!!          4 (5) forward FFTs; in addition 1
!!          forward (w_r) and 1 inverse (d rho/dt).
!!          svv_mass = F: 3 inverse and 3 forward
!!          per direction, nothing else.
!-----------------------------------------------!
  subroutine svv_es_core(s_r, rho_r, ux_r, uy_r, uz_r, uh, eps, rdot_r, mdh, dsum)
    use grid, only: nlx, nly, nlz, nlx_local, nkx, nky_local, nkz
    use params, only: cs2va2, filter_mom, svv_mass
    use mp, only: max_allreduce, comm_fft
    use cufftmp, only: ftran_r2c, btran_c2r
    implicit none
    real(8), dimension(:,:,:), intent(in ) :: s_r, rho_r, ux_r, uy_r, uz_r
    complex(8), dimension(:,:,:,:), intent(in ) :: uh
    real(8), intent(in ) :: eps(3)
    real(8), dimension(:,:,:), intent(out) :: rdot_r
    complex(8), dimension(:,:,:,:), intent(out) :: mdh
    real(8), intent(out) :: dsum
    real(8) :: umax2, a, b, kap, t, fr, g1, g2, g3, cs2, epsd, sd, sh
    logical :: lmass, lmom
    integer :: i, j, k, c, d

    cs2   = cs2va2
    lmass = svv_mass
    lmom  = filter_mom

    umax2 = 0.d0
    !$acc data present(ux_r, uy_r, uz_r)
    !$acc parallel loop collapse(3) reduction(max:umax2)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          umax2 = max(umax2, ux_r(k,j,i)**2 + uy_r(k,j,i)**2 + uz_r(k,j,i)**2)
        enddo
      enddo
    enddo
    !$acc end data
    call max_allreduce(umax2, comm=comm_fft)

    ! w_r = cs^2 s - u^2/2 (not needed without the mass fluxes)
    if(lmass) then
      !$acc data present(s_r, ux_r, uy_r, uz_r, es_gr)
      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            es_gr(k,j,i) = cs2*s_r(k,j,i) - 0.5d0*(ux_r(k,j,i)**2 + uy_r(k,j,i)**2 + uz_r(k,j,i)**2)
          enddo
        enddo
      enddo
      !$acc end data
      call ftran_r2c(es_gr, es_wh)
      !$acc data present(es_wh)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            es_wh(k,j,i) = es_wh(k,j,i)/nlx/nly/nlz
          enddo
        enddo
      enddo
      !$acc end data
    endif

    !$acc data present(es_rdh, mdh)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          es_rdh(k,j,i) = (0.d0, 0.d0)
          mdh(k,j,i,1)  = (0.d0, 0.d0)
          mdh(k,j,i,2)  = (0.d0, 0.d0)
          mdh(k,j,i,3)  = (0.d0, 0.d0)
        enddo
      enddo
    enddo
    !$acc end data
    dsum = 0.d0
    do d = 1, 3
      if(eps(d) == 0.d0) cycle ! inactive direction or svv off
      epsd = eps(d)

      ! gradients g_r = A_L w_r, g_m = A_L u, g_h = A_H w_r
      if(lmass) call svv_grad(es_wh, d, .true., es_gr)
      do c = 1, 3
        call svv_grad(uh(:,:,:,c), d, .true., es_gm(:,:,:,c))
      enddo
      if(lmass .and. lmom) call svv_grad(es_wh, d, .false., es_gh)

      ! fluxes (overwrite the gradients) and the dissipation density
      sd = 0.d0
      if(lmass) then
        kap = epsd/(cs2 + umax2)
        !$acc data present(rho_r, ux_r, uy_r, uz_r, es_gr, es_gm)
        !$acc parallel loop collapse(3) private(a, b, t, fr, g1, g2, g3) reduction(+:sd)
        do i = 1, nlx_local
          do j = 1, nly
            do k = 1, nlz
              a  = epsd*rho_r(k,j,i)/cs2
              b  = epsd*rho_r(k,j,i)
              g1 = es_gm(k,j,i,1)
              g2 = es_gm(k,j,i,2)
              g3 = es_gm(k,j,i,3)
              t  = es_gr(k,j,i) + ux_r(k,j,i)*g1 + uy_r(k,j,i)*g2 + uz_r(k,j,i)*g3
              sd = sd + a*t**2 + b*(g1**2 + g2**2 + g3**2)
              fr = a*t
              es_gr(k,j,i)   = fr
              es_gm(k,j,i,1) = ux_r(k,j,i)*fr + b*g1
              es_gm(k,j,i,2) = uy_r(k,j,i)*fr + b*g2
              es_gm(k,j,i,3) = uz_r(k,j,i)*fr + b*g3
            enddo
          enddo
        enddo
        !$acc end data
        dsum = dsum + sd
        if(lmom) then
          sh = 0.d0
          !$acc data present(rho_r, es_gh)
          !$acc parallel loop collapse(3) reduction(+:sh)
          do i = 1, nlx_local
            do j = 1, nly
              do k = 1, nlz
                sh = sh + kap*rho_r(k,j,i)*es_gh(k,j,i)**2
                es_gh(k,j,i) = kap*rho_r(k,j,i)*es_gh(k,j,i)
              enddo
            enddo
          enddo
          !$acc end data
          dsum = dsum + sh
        endif
      else
        ! f_r = f_h = 0, f_m = b g_m
        !$acc data present(rho_r, es_gm)
        !$acc parallel loop collapse(3) private(b, g1, g2, g3) reduction(+:sd)
        do i = 1, nlx_local
          do j = 1, nly
            do k = 1, nlz
              b  = epsd*rho_r(k,j,i)
              g1 = es_gm(k,j,i,1)
              g2 = es_gm(k,j,i,2)
              g3 = es_gm(k,j,i,3)
              sd = sd + b*(g1**2 + g2**2 + g3**2)
              es_gm(k,j,i,1) = b*g1
              es_gm(k,j,i,2) = b*g2
              es_gm(k,j,i,3) = b*g3
            enddo
          enddo
        enddo
        !$acc end data
        dsum = dsum + sd
      endif

      ! adjoints: d rho/dt -= A_L^* f_r + A_H^* f_h, dm/dt -= A_L^* f_m
      if(lmass) then
        call svv_grad_adj(es_gr, d, .true.)
        call sub_fh(es_rdh)
        if(lmom) then
          call svv_grad_adj(es_gh, d, .false.)
          call sub_fh(es_rdh)
        endif
      endif
      do c = 1, 3
        call svv_grad_adj(es_gm(:,:,:,c), d, .true.)
        call sub_fh(mdh(:,:,:,c))
      enddo
    enddo

    if(lmass) then
      call btran_c2r(es_rdh, rdot_r)
    else
      !$acc data present(rdot_r)
      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            rdot_r(k,j,i) = 0.d0
          enddo
        enddo
      enddo
      !$acc end data
    endif

  end subroutine svv_es_core


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   Helpers of svv_es_core (its contained
!!          grad, grad_adj and mult_sig in Calliope);
!!          es_fh is the work array, d the direction.
!!          svv_grad:     y  = F^-1 sig_d mf x,
!!                        mf = M (low = T) or 1 - M (low = F)
!!          svv_grad_adj: es_fh = conj(sig_d) mf F y/(nlx nly nlz)
!!          svv_mult_sig: es_fh <- sig_d mf es_fh
!!                        (conj(sig_d) for adj = T)
!!          sub_fh:       x <- x - es_fh
!-----------------------------------------------!
  subroutine svv_grad(x, d, low, y)
    use grid, only: nkx, nky_local, nkz
    use cufftmp, only: btran_c2r
    implicit none
    complex(8), dimension(:,:,:), intent(in ) :: x
    integer, intent(in) :: d
    logical, intent(in) :: low
    real(8), dimension(:,:,:), intent(out) :: y
    integer :: i, j, k

    !$acc data present(x, es_fh)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          es_fh(k,j,i) = x(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data
    call svv_mult_sig(d, low, .false.)
    call btran_c2r(es_fh, y)

  end subroutine svv_grad


  subroutine svv_grad_adj(y, d, low)
    use grid, only: nlx, nly, nlz, nkx, nky_local, nkz
    use cufftmp, only: ftran_r2c
    implicit none
    real(8), dimension(:,:,:), intent(in) :: y
    integer, intent(in) :: d
    logical, intent(in) :: low
    integer :: i, j, k

    call ftran_r2c(y, es_fh)
    !$acc data present(es_fh)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          es_fh(k,j,i) = es_fh(k,j,i)/nlx/nly/nlz
        enddo
      enddo
    enddo
    !$acc end data
    call svv_mult_sig(d, low, .true.)

  end subroutine svv_grad_adj


  subroutine svv_mult_sig(d, low, adj)
    use grid, only: nkx, nky_local, nkz
    use dealias, only: filter
    use params, only: filter_mom
    implicit none
    integer, intent(in) :: d
    logical, intent(in) :: low, adj
    complex(8) :: sg
    real(8) :: mf
    integer :: i, j, k, ld
    logical :: llow, ladj, lmom

    ld   = d
    llow = low
    ladj = adj
    lmom = filter_mom
    !$acc data present(es_fh, filter)
    !$acc parallel loop collapse(3) private(sg, mf)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if(ld == 1) then
            sg = sgx(i)
          elseif(ld == 2) then
            sg = sgy(j)
          else
            sg = sgz(k)
          endif
          if(ladj) sg = conjg(sg)
          if(lmom) then
            mf = filter(k,j,i)
          else
            mf = 1.d0
          endif
          if(.not. llow) mf = 1.d0 - mf
          es_fh(k,j,i) = sg*mf*es_fh(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine svv_mult_sig


  subroutine sub_fh(x)
    use grid, only: nkx, nky_local, nkz
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: x
    integer :: i, j, k

    !$acc data present(x, es_fh)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          x(k,j,i) = x(k,j,i) - es_fh(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine sub_fh


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   SVV tendency of the state (s, m) with
!!          eps fixed: Gs = (d rho/dt)/rho,
!!          Gm = dm/dt (real space)
!-----------------------------------------------!
  subroutine svv_G(s_r, mx_r, my_r, mz_r, eps, Gs, Gm)
    use grid, only: nlx_local, nly, nlz
    use cufftmp, only: btran_c2r
    implicit none
    real(8), dimension(:,:,:), intent(in ) :: s_r, mx_r, my_r, mz_r
    real(8), intent(in ) :: eps(3)
    real(8), dimension(:,:,:)  , intent(out) :: Gs
    real(8), dimension(:,:,:,:), intent(out) :: Gm
    real(8) :: dsum
    integer :: i, j, k, c

    !$acc data present(s_r, mx_r, my_r, mz_r, G_rho, G_ux, G_uy, G_uz)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          G_rho(k,j,i) = exp(s_r(k,j,i))
          G_ux(k,j,i)  = mx_r(k,j,i)/G_rho(k,j,i)
          G_uy(k,j,i)  = my_r(k,j,i)/G_rho(k,j,i)
          G_uz(k,j,i)  = mz_r(k,j,i)/G_rho(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data
    call get_svv_es(s_r, G_rho, G_ux, G_uy, G_uz, eps, Gs, G_mdh, dsum)
    !$acc data present(Gs, G_rho)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          Gs(k,j,i) = Gs(k,j,i)/G_rho(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data
    do c = 1, 3
      call btran_c2r(G_mdh(:,:,:,c), Gm(:,:,:,c))
    enddo

  end subroutine svv_G


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Power method for svv_lam, the largest
!!          |eigenvalue| of the Jacobian J of
!!          G(s, m) = ((d rho/dt)/rho, dm/dt) of the
!!          entropy-stable SVV with eps fixed:
!!            J v = (G(U + h v) - G(U))/h, h = 1e-7,
!!            v <- J v/|v|, svv_lam = |J v|/|v|,
!!          |v|^2 = sum (v_s^2 + |v_m|^2)/(nlx nly nlz).
!!          The vector is kept between calls (not
!!          across a restart). Its initial value is a
!!          fixed function of the global grid index,
!!          so that it does not depend on the
!!          decomposition. niter iterations at most
!!          (max(niter, svv_pi_init) when the vector is
!!          initialized); with tol > 0 it stops when
!!          two successive estimates of this call agree
!!          within tol (relative). nit: iterations
!!          done; conv: the tol test was met.
!!          G(U) is given by rdot0_r (d rho/dt) and
!!          mdh0 (dm/dt, Fourier) of get_svv_es.
!!          svv_lam is set on the host and the device.
!-----------------------------------------------!
  subroutine svv_power(s_r, rho_r, mx_r, my_r, mz_r, eps, rdot0_r, mdh0, niter, tol, nit, conv)
    use grid, only: nlx_local, nly, nlz
    use params, only: svv_pi_init
    use cufftmp, only: btran_c2r
    implicit none
    real(8), dimension(:,:,:), intent(in) :: s_r, rho_r, mx_r, my_r, mz_r, rdot0_r
    complex(8), dimension(:,:,:,:), intent(in) :: mdh0
    real(8), intent(in ) :: eps(3), tol
    integer, intent(in ) :: niter
    integer, intent(out) :: nit
    logical, intent(out) :: conv
    real(8), parameter :: h = 1.d-7
    real(8) :: nv, lam, lam_old, v1, v2, v3, v4
    integer :: i, j, k, c, it, nmax

    nit  = 0
    conv = .false.
    if(all(eps == 0.d0)) then
      svv_lam = 0.d0
      !$acc update device(svv_lam)
      return
    endif
    if(.not. svv_es) call abort_run('svv_power: called without svv_form = entropy (no work arrays)')

    nmax = niter
    if(.not. pi_set) then
      nmax = max(niter, svv_pi_init)
      ! (-1)^(sum_d (i_d - 1)) (1 + 0.3 sum_d cos(2 pi (i_d - 1)/nl_d + 1)) in s, 0 in m, sums over the
      ! active directions d: with nl_d = 2 the modes (-1)^(i_d - 1) of an inactive direction are not coupled
      ! by J to the fields uniform in it, and their stiffness is lower (M = 0 there, so only A_H acts)
      call seed_pi_v(.false.)
      pi_set = .true.
    endif

    ! G(U)
    !$acc data present(rdot0_r, rho_r, pw_G0s)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          pw_G0s(k,j,i) = rdot0_r(k,j,i)/rho_r(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data
    do c = 1, 3
      call btran_c2r(mdh0(:,:,:,c), pw_G0m(:,:,:,c))
    enddo

    lam_old = 0.d0
    lam     = 0.d0
    do it = 1, nmax
      nv = pi_norm()
      if(nv == 0.d0) then
        ! J v = 0 (degenerate); restart from the initial vector at the next call
        pi_set = .false.
        lam = 0.d0
        exit
      endif
      !$acc data present(s_r, mx_r, my_r, mz_r, pi_v, pw_ps, pw_pm)
      !$acc parallel loop collapse(3) private(v1, v2, v3, v4)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            v1 = pi_v(k,j,i,1)/nv
            v2 = pi_v(k,j,i,2)/nv
            v3 = pi_v(k,j,i,3)/nv
            v4 = pi_v(k,j,i,4)/nv
            pi_v(k,j,i,1) = v1
            pi_v(k,j,i,2) = v2
            pi_v(k,j,i,3) = v3
            pi_v(k,j,i,4) = v4
            pw_ps(k,j,i)   = s_r(k,j,i)  + h*v1
            pw_pm(k,j,i,1) = mx_r(k,j,i) + h*v2
            pw_pm(k,j,i,2) = my_r(k,j,i) + h*v3
            pw_pm(k,j,i,3) = mz_r(k,j,i) + h*v4
          enddo
        enddo
      enddo
      !$acc end data
      call svv_G(pw_ps, pw_pm(:,:,:,1), pw_pm(:,:,:,2), pw_pm(:,:,:,3), eps, pw_Gs, pw_Gm)
      !$acc data present(pi_v, pw_Gs, pw_G0s, pw_Gm, pw_G0m)
      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            pi_v(k,j,i,1) = (pw_Gs(k,j,i)   - pw_G0s(k,j,i)  )/h
            pi_v(k,j,i,2) = (pw_Gm(k,j,i,1) - pw_G0m(k,j,i,1))/h
            pi_v(k,j,i,3) = (pw_Gm(k,j,i,2) - pw_G0m(k,j,i,2))/h
            pi_v(k,j,i,4) = (pw_Gm(k,j,i,3) - pw_G0m(k,j,i,3))/h
          enddo
        enddo
      enddo
      !$acc end data
      lam  = pi_norm()
      nit  = it
      conv = tol > 0.d0 .and. it > 1 .and. abs(lam - lam_old) <= tol*lam
      lam_old = lam
      if(conv) exit
    enddo
    svv_lam = lam
    !$acc update device(svv_lam)

  end subroutine svv_power


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   |pi_v| = sqrt(sum (v_s^2 + |v_m|^2)
!!          /(nlx nly nlz)) over all ranks
!!          (norm() contained in svv_power in Calliope)
!-----------------------------------------------!
  real(8) function pi_norm()
    use grid, only: nlx, nly, nlz, nlx_local
    use mp, only: sum_allreduce, comm_fft
    implicit none
    real(8) :: ss
    integer :: i, j, k

    ss = 0.d0
    !$acc data present(pi_v)
    !$acc parallel loop collapse(3) reduction(+:ss)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          ss = ss + pi_v(k,j,i,1)**2 + pi_v(k,j,i,2)**2 + pi_v(k,j,i,3)**2 + pi_v(k,j,i,4)**2
        enddo
      enddo
    enddo
    !$acc end data
    call sum_allreduce(ss, comm=comm_fft)
    pi_norm = sqrt(ss/(dble(nlx)*dble(nly)*dble(nlz)))

  end function pi_norm


!-----------------------------------------------!
!> @author  YK
!! @date    14 Nov 2026
!! @brief   Initial vector of the power method, on
!!          the host from the global grid index
!!          (x: i + nlx_local*iproc_fft), then copied
!!          to the device. s part: the vector of
!!          svv_power; m part (with_m = T, from
!!          svv_seed_nomass) or 0. Calliope writes the
!!          s formula in both routines; here it is
!!          written once, with the same arithmetic.
!-----------------------------------------------!
  subroutine seed_pi_v(with_m)
    use grid, only: nlx, nly, nlz, nlx_local
    use params, only: pi
    use mp, only: iproc_fft
    implicit none
    logical, intent(in) :: with_m
    real(8) :: amp
    integer :: i, j, k, ig, nsg

    pi_v = 0.d0
    do i = 1, nlx_local
      ig = i + nlx_local*iproc_fft
      do j = 1, nly
        do k = 1, nlz
          amp = 1.d0
          nsg = 0
          if(active_x) then
            amp = amp + 0.3d0*cos(2.d0*pi*dble(ig - 1)/dble(nlx) + 1.d0)
            nsg = nsg + ig - 1
          endif
          if(active_y) then
            amp = amp + 0.3d0*cos(2.d0*pi*dble(j - 1)/dble(nly) + 1.d0)
            nsg = nsg + j - 1
          endif
          if(active_z) then
            amp = amp + 0.3d0*cos(2.d0*pi*dble(k - 1)/dble(nlz) + 1.d0)
            nsg = nsg + k - 1
          endif
          pi_v(k,j,i,1) = merge(1.d0, -1.d0, mod(nsg, 2) == 0)*amp

          if(with_m) then
            amp = 1.d0
            if(active_x) amp = amp*cos(2.d0*pi*dble(((nlx - 1)/3)*(ig - 1))/dble(nlx) + 1.d0)
            if(active_y) amp = amp*cos(2.d0*pi*dble(((nly - 1)/3)*(j - 1))/dble(nly) + 1.d0)
            if(active_z) amp = amp*cos(2.d0*pi*dble(((nlz - 1)/3)*(k - 1))/dble(nlz) + 1.d0)
            pi_v(k,j,i,2:4) = amp
          endif
        enddo
      enddo
    enddo
    !$acc update device(pi_v)

  end subroutine seed_pi_v


!-----------------------------------------------!
!> @author  YK
!! @date    9 Nov 2026
!! @brief   Initial vector of svv_power for
!!          svv_mass = F
!! @details If pi_v holds no vector, puts the vector
!!          below in it and returns seeded = .true.;
!!          otherwise leaves it and returns
!!          seeded = .false.
!!          Called for svv_mass = F only, before
!!          svv_power, which then finds pi_v set.
!!          In s the vector is the initial vector of
!!          svv_power, and in m it is
!!            prod_d cos(2 pi K_d (i_d - 1)/nl_d + 1)
!!          over the active d, K_d = (nl_d - 1)/3 the
!!          top index kept by the 2/3 filter on m.
!-----------------------------------------------!
  subroutine svv_seed_nomass(seeded)
    implicit none
    logical, intent(out) :: seeded

    if(.not. svv_es) call abort_run('svv_seed_nomass: called without svv_form = entropy (no work arrays)')
    seeded = .not. pi_set
    if(.not. seeded) return

    ! svv_mass = F: J = [[0, 0], [J_ms, J_mm]] has no ln rho row, and J v = 0 for the vector of svv_power (s only)
    ! if A_L(u v_s) = 0 (e.g. s = 0 and u only at |k_d| <= svv_mN k_c,d); the iteration would then start from
    ! round-off, which depends on the decomposition. The m part: for uniform rho it is an eigenvector of J_mm with
    ! the largest eigenvalue on the mask, and in general J_mm v_m = 0 only if A_L,d(v_m/rho) = 0 for all d
    ! (<v_m/rho, J_mm v_m> = -sum_d eps_d sum W |A_L,d(v_m/rho)|^2).
    ! The s part is the same formula as in svv_power (seed_pi_v).
    call seed_pi_v(.true.)
    pi_set = .true.

  end subroutine svv_seed_nomass


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   zeta = C_shock*Delta^2*S[max(0, -D)]
!!          and (optional) the shock resistivity
!!          eta = C_eta*Delta^2*S[max(0, -D)]
!!          The passes alternate between zeta_r and
!!          the work array zw (in_z: S is in zeta_r).
!-----------------------------------------------!
  subroutine get_zeta(div_r, zeta_r, eta_r)
    use grid, only: nlx_local, nly, nlz
    use params, only: C_shock, shock_nmax, shock_nsmooth, C_eta
    implicit none
    real(8), dimension(:,:,:), intent(in ) :: div_r
    real(8), dimension(:,:,:), intent(out) :: zeta_r
    real(8), dimension(:,:,:), intent(out), optional :: eta_r
    real(8) :: ce, cz
    logical :: in_z
    integer :: i, j, k, n

    !$acc data present(div_r, zeta_r)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          zeta_r(k,j,i) = max(0.d0, -div_r(k,j,i))
        enddo
      enddo
    enddo
    !$acc end data
    in_z = .true.

    do n = 1, shock_nmax
      if(active_y) call pass('y', 'm')
      if(active_z) call pass('z', 'm')
      if(active_x) call pass('x', 'm')
    enddo
    do n = 1, shock_nsmooth
      if(active_y) call pass('y', 's')
      if(active_z) call pass('z', 's')
      if(active_x) call pass('x', 's')
    enddo

    if(.not. in_z) then
      !$acc data present(zeta_r, zw)
      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            zeta_r(k,j,i) = zw(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data
    endif

    ! the products are evaluated as in Calliope, (C*Delta**2)*S
    ce = C_eta*delta_shock**2
    cz = C_shock*delta_shock**2
    if(present(eta_r)) then
      !$acc data present(zeta_r, eta_r)
      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            eta_r(k,j,i)  = ce*zeta_r(k,j,i)
            zeta_r(k,j,i) = cz*zeta_r(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data
    else
      !$acc data present(zeta_r)
      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            zeta_r(k,j,i) = cz*zeta_r(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data
    endif

  contains

    ! one pass from the array that holds S into the other one
    subroutine pass(dir, op)
      character(len=1), intent(in) :: dir, op
      if(in_z) then
        call pass_1d(zeta_r, zw, dir, op)
      else
        call pass_1d(zw, zeta_r, dir, op)
      endif
      in_z = .not. in_z
    end subroutine pass

  end subroutine get_zeta


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   One 3-point pass in direction dir, g = P f
!!          op = 'm': f_j <- max(f_{j-1}, f_j, f_{j+1})
!!          op = 's': f_j <- (f_{j-1} + 2 f_j + f_{j+1})/4,
!!                    evaluated as 0.5*(f_j + 0.5*(f_{j-1} + f_{j+1}))
!!                    so that a uniform field is kept exactly
!!          y, z: periodic wrap inside the rank;
!!          x: halo planes from halo_exchange_x
!!          (also for nproc_fft = 1, where it copies)
!-----------------------------------------------!
  subroutine pass_1d(f, g, dir, op)
    use grid, only: nlx_local, nly, nlz
    use mp, only: halo_exchange_x
    implicit none
    real(8), dimension(:,:,:), intent(in ) :: f
    real(8), dimension(:,:,:), intent(out) :: g
    character(len=1), intent(in) :: dir, op
    real(8) :: fm, fc, fp
    integer :: i, j, k, idir, nx, ny, nz
    logical :: lmax

    select case(dir)
    case('y')
      idir = 2
    case('z')
      idir = 3
    case default
      idir = 1
    end select
    lmax = op == 'm'
    nx = nlx_local
    ny = nly
    nz = nlz

    if(idir == 1) call halo_exchange_x(f, halo_lo, halo_hi)

    !$acc data present(f, g, halo_lo, halo_hi)
    !$acc parallel loop collapse(3) private(fm, fc, fp)
    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          fc = f(k,j,i)
          if(idir == 2) then
            if(j == 1) then
              fm = f(k,ny,i)
            else
              fm = f(k,j-1,i)
            endif
            if(j == ny) then
              fp = f(k,1,i)
            else
              fp = f(k,j+1,i)
            endif
          elseif(idir == 3) then
            if(k == 1) then
              fm = f(nz,j,i)
            else
              fm = f(k-1,j,i)
            endif
            if(k == nz) then
              fp = f(1,j,i)
            else
              fp = f(k+1,j,i)
            endif
          else
            if(i == 1) then
              fm = halo_lo(k,j)
            else
              fm = f(k,j,i-1)
            endif
            if(i == nx) then
              fp = halo_hi(k,j)
            else
              fp = f(k,j,i+1)
            endif
          endif
          if(lmax) then
            g(k,j,i) = max(fm, fc, fp)
          else
            g(k,j,i) = 0.5d0*(fc + 0.5d0*(fm + fp))
          endif
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine pass_1d


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Self test (shock_selftest = T):
!!          zeta of an analytic D on the grid is
!!          gathered to proc0 and written to
!!          shock_selftest.dat (es25.17). Files from
!!          different decompositions must be identical.
!-----------------------------------------------!
  subroutine selftest
    use grid, only: nlx, nly, nlz, nlx_local, nlz_padded, xx, yy, zz, lx, ly, lz
    use params, only: pi
    use mp, only: proc0, sum_reduce, comm_fft, iproc_fft
    use file, only: open_output_file
    implicit none
    real(8), allocatable, dimension(:,:,:) :: d_r, z_r, z_all
    integer :: i, j, k, unit

    allocate(d_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(z_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(z_all(nly, nlz, nlx), source=0.d0)

    ! xx is local: xx(i) is at the global index i + nlx_local*iproc_fft
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          d_r(k,j,i) = sin(2.d0*pi*xx(i)/lx + 0.3d0)*cos(4.d0*pi*yy(j)/ly + 1.1d0) &
                     + 0.7d0*sin(6.d0*pi*zz(k)/lz + 0.5d0) - 0.2d0
        enddo
      enddo
    enddo

    !$acc enter data copyin(d_r, z_r)
    call get_zeta(d_r, z_r)
    !$acc update self(z_r)
    !$acc exit data delete(d_r, z_r)

    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz
          z_all(j, k, i + nlx_local*iproc_fft) = z_r(k,j,i)
        enddo
      enddo
    enddo
    call sum_reduce(z_all, 0, comm=comm_fft) ! exact: other ranks contribute 0

    if(proc0) then
      call open_output_file(unit, 'shock_selftest.dat')
      do i = 1, nlx
        do k = 1, nlz
          do j = 1, nly
            write(unit, '(3i6, es25.17)') i, j, k, z_all(j,k,i)
          enddo
        enddo
      enddo
      close(unit)
      print '("  shock_selftest.dat written: min/max zeta = ", 2es12.4)', minval(z_all), maxval(z_all)
    endif

    deallocate(d_r, z_r, z_all)

  end subroutine selftest


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   True for Inf/NaN. Tests the exponent
!!          bits so that it is not removed by
!!          value-unsafe optimizations (-O3, fast fp).
!!          Callable in kernels (acc routine seq).
!!          -DISNONFINITE_IEEE: ieee_is_finite instead
!!          of the bit test, in case transfer is not
!!          accepted in device code.
!-----------------------------------------------!
  elemental logical function is_nonfinite(x)
    !$acc routine seq
#ifdef ISNONFINITE_IEEE
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8), intent(in) :: x

    is_nonfinite = .not. ieee_is_finite(x)
#else
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    real(8), intent(in) :: x
    integer(int64) :: bits

    bits = transfer(x, 0_int64)
    is_nonfinite = iand(ishft(bits, -52), 2047_int64) == 2047_int64
#endif

  end function is_nonfinite


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Print a message with the rank and
!!          abort all ranks (non-zero exit code).
!!          Does not require the other ranks to
!!          reach the same point.
!-----------------------------------------------!
  subroutine abort_run(msg)
    use mp, only: proc_id
    implicit none
    character(len=*), intent(in) :: msg
    integer :: ierr

    write(*, '("[proc ", i6, "] ABORT: ", a)') proc_id, trim(msg)
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  end subroutine abort_run


!-----------------------------------------------!
!> @author  YK
!! @date    3 Nov 2026
!! @brief   Start and stop the timer i (MPI_Wtime
!!          on each rank; no synchronization)
!-----------------------------------------------!
  subroutine tmr_on(i)
    implicit none
    integer, intent(in) :: i

    tmr_t0(i) = MPI_Wtime()

  end subroutine tmr_on

  subroutine tmr_off(i)
    implicit none
    integer, intent(in) :: i

    tmr(i) = tmr(i) + (MPI_Wtime() - tmr_t0(i))

  end subroutine tmr_off


!-----------------------------------------------!
!> @author  YK
!! @date    3 Nov 2026
!! @brief   Print the min and max over the ranks
!!          of the timers (collective)
!-----------------------------------------------!
  subroutine report_timers
    use mp, only: proc0, max_reduce, min_reduce
    implicit none
    real(8) :: tmin(ntmr), tmax(ntmr)
    integer :: i

    tmin = tmr
    tmax = tmr
    call min_reduce(tmin, 0)
    call max_reduce(tmax, 0)
    if(proc0) then
      print '(/, " MHD_COMP_ISOTH timers [s] (min, max over ranks; % of max step)")'
      do i = 1, ntmr
        print '(3x, a28, 2es12.4, f8.2, " %")', tmr_name(i), tmin(i), tmax(i), &
              100.d0*tmax(i)/max(tmax(itmr_step), tiny(1.d0))
      enddo
    endif

  end subroutine report_timers

end module shock_viscosity
