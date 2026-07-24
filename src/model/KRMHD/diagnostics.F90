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
    use grid, only: nkx, nky_local, nkz, nm
    use grid, only: nlx, nly, nlz_padded
    use fields, only: phi, psi
    use fields, only: phi_old, psi_old
    use mp, only: sum_reduce, comm_fft, nproc_m, nproc_s
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, alpha, write_hermite_flux, &
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
    real(8), dimension(:, :, :), allocatable :: g2bin                 ! [kprp, m, kz]
    real(8), dimension(:, :, :), allocatable :: gam_bin               ! [kprp, m, kz] Hermite flux (eq 9)
    real(8), dimension(:), allocatable :: Wm                          ! per-m free energy
    real(8), dimension(:), allocatable :: gam_kint                    ! k-integrated Hermite flux Gamma(m)
    real(8) :: W_free                                                 ! Meyrand 2019 eq (6)
    real(8) :: p_g_sum, Dg_sum                                        ! g power balance: dW_free/dt = P_g - D_g
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

    allocate (g2bin(1:nkpolar, nm, nkz)); g2bin = 0.d0
    allocate (gam_bin(1:nkpolar, nm, nkz)); gam_bin = 0.d0
    allocate (Wm(nm)); Wm = 0.d0
    allocate (gam_kint(nm)); gam_kint = 0.d0

    ! Create work arrays on device
    !$acc enter data create(upe2, bpe2, upe2old, bpe2old)
    !$acc enter data create(upe2dissip_x, upe2dissip_z, bpe2dissip_x, bpe2dissip_z)
    !$acc enter data create(p_phi, p_psi, p_xhl, zppe2, zmpe2)
    !$acc enter data create(upe2_bin, bpe2_bin)
    !$acc enter data create(zppe2_bin, zmpe2_bin)
    !$acc enter data create(g2bin)
    !$acc enter data create(gam_bin)

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

    ! MPI reduction for summed quantities (operates on host scalars).
    ! These are grid/spectral partial sums, so reduce on comm_fft. The global
    ! proc0 is the root of comm_fft group 0 (contiguous rank layout), so the
    ! reduced value still lands there and the existing "if(proc0) write" holds.
    call sum_reduce(upe2_sum       , 0, comm=comm_fft)
    call sum_reduce(bpe2_sum       , 0, comm=comm_fft)
    call sum_reduce(upe2dot_sum    , 0, comm=comm_fft)
    call sum_reduce(bpe2dot_sum    , 0, comm=comm_fft)
    call sum_reduce(upe2dissip_sum , 0, comm=comm_fft)
    call sum_reduce(bpe2dissip_sum , 0, comm=comm_fft)
    call sum_reduce(zppe2_sum      , 0, comm=comm_fft)
    call sum_reduce(zmpe2_sum      , 0, comm=comm_fft)
    call sum_reduce(p_phi_sum      , 0, comm=comm_fft)
    call sum_reduce(p_psi_sum      , 0, comm=comm_fft)
    call sum_reduce(p_xhl_sum      , 0, comm=comm_fft)

    ! Bin spectra over kprp on device
    call get_polar_spectrum_2d(upe2 , upe2_bin )
    call get_polar_spectrum_2d(bpe2 , bpe2_bin )
    call get_polar_spectrum_2d(zppe2, zppe2_bin)
    call get_polar_spectrum_2d(zmpe2, zmpe2_bin)

    ! Transfer binned spectra from device to host for I/O
    !$acc update host(upe2_bin, bpe2_bin)
    !$acc update host(zppe2_bin, zmpe2_bin)

    ! Per-m Hermite spectrum + exact per-m free energy (assembled on proc0 via
    ! comm_fm inside the helper: it already updates host and reduces).
    call get_g_spectrum_2d(g2bin, Wm)
    ! Meyrand 2019 eq (6) free energy: sum_m W_m + alpha*W_{m=0} (Wm(1)=m=0 rung).
    ! Only proc0 holds the fully comm_fm-reduced Wm; other ranks keep partials.
    W_free = 0.d0
    if (proc0) W_free = sum(Wm) + alpha*Wm(1)

    ! g free-energy power balance scalars P_g (injection) and D_g (dissipation),
    ! assembled on proc0 via comm_fm inside the helper. Post-processing checks
    ! dW_free/dt = P_g - D_g (dW_free/dt from finite-differencing W_free(t)).
    call get_g_power_balance(p_g_sum, Dg_sum)

    ! Hermite free-energy flux Gamma_m(k_perp) (Meyrand 2019 eq 9). Costly (extra
    ! FFTs to reconstruct grad_par S_m), so it is gated behind write_hermite_flux.
    ! When disabled, gam_bin stays zero and is not written to NetCDF.
    if (write_hermite_flux) call get_hermite_flux_2d(gam_bin, gam_kint)

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
                  zppe2_bin, zmpe2_bin, &
                  !
                  nm, g2bin, Wm, W_free, p_g_sum, Dg_sum, gam_bin, gam_kint &
                )

    !$acc exit data delete(upe2, bpe2)
    !$acc exit data delete(upe2old, bpe2old)
    !$acc exit data delete(upe2dissip_x, upe2dissip_z, bpe2dissip_x, bpe2dissip_z)
    !$acc exit data delete(p_phi, p_psi, p_xhl)
    !$acc exit data delete(zppe2, zmpe2)
    !$acc exit data delete(upe2_bin, bpe2_bin)
    !$acc exit data delete(zppe2_bin, zmpe2_bin)
    !$acc exit data delete(g2bin)
    !$acc exit data delete(gam_bin)
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
    deallocate (g2bin)
    deallocate (gam_bin)
    deallocate (Wm)
    deallocate (gam_kint)

    ! With P_m>1 (and/or P_s>1) every comm_fft group redundantly solves the
    ! same problem; verify the fields stay bitwise identical across groups.
    if (nproc_m*nproc_s > 1) call check_redundant_consistency
    ! One-time derived-alpha banner, plus a g redundancy assert across comm_s
    ! when P_s>1 (silent no-op otherwise).
    call check_g_consistency
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @brief   Per-m Hermite spectrum (1/2)|g_m|^2 binned over
!!          kperp, plus the exact per-m free energy Wm.
!!          Mirrors the RMHD polar binning (get_polar_spectrum_2d)
!!          but (a) writes into the global Hermite slot
!!          mg = m_offset + mm and (b) assembles kperp (comm_fft)
!!          and m (comm_m) together with a single comm_fm
!!          reduction, so no per-axis double reduce is needed.
!!          The kz/=1 factor-of-2 follows the RMHD energy
!!          convention (compensates the unstored negative-kz half).
!-----------------------------------------------!
  subroutine get_g_spectrum_2d(g2bin, Wm)
    use mp, only: sum_reduce, comm_fm
    use grid, only: nkx, nky_local, nkz
    use grid, only: nm, nm_local, m_offset
    use grid, only: kx, ky
    use fields, only: g
    use utils, only: cabs2
    implicit none
    real(8), dimension(1:nkpolar, nm, nkz), intent(out) :: g2bin  ! device (present)
    real(8), dimension(nm),                 intent(out) :: Wm     ! host
    real(8) :: k2, g2, Wsum
    integer :: i, j, k, mm, mg, idx

    ! Zero the binned accumulator on device and the host per-m energy.
    !$acc parallel loop collapse(3) default(present)
    do k = 1, nkz
      do mg = 1, nm
        do idx = 1, nkpolar
          g2bin(idx, mg, k) = 0.d0
        end do
      end do
    end do
    !$acc end parallel loop
    Wm = 0.d0

    ! Bin (1/2)|g_m|^2 into the owned global slot mg; the unowned slots stay 0
    ! so the comm_fm reduction assembles the full (kperp x m) spectrum exactly.
    do mm = 1, nm_local
      mg   = m_offset + mm
      Wsum = 0.d0
      !$acc parallel loop collapse(3) gang vector default(present) &
      !$acc& private(k2, g2, idx) reduction(+:Wsum)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            g2 = 0.5d0*cabs2(g(k, j, i, mm))
            ! compensate the unstored negative-kz half (RMHD convention)
            if (k /= 1) g2 = 2.0d0*g2
            Wsum = Wsum + g2
            k2 = kx(i)**2 + ky(j)**2
            do idx = 1, nkpolar - 1
              if (k2 >= kpbin(idx)**2 .and. k2 < kpbin(idx+1)**2) then
                !$acc atomic update
                g2bin(idx, mg, k) = g2bin(idx, mg, k) + g2
                exit  ! found the bin
              endif
            enddo
          end do
        end do
      end do
      !$acc end parallel loop
      Wm(mg) = Wsum
    end do

    ! Bring the binned spectrum to host and assemble both axes (kperp over
    ! comm_fft, m over comm_m) in one reduction on the comm_fm plane. g is
    ! redundant across comm_s, so comm_fm covers the whole field exactly once.
    !$acc update host(g2bin)
    call sum_reduce(g2bin, 0, comm=comm_fm)
    call sum_reduce(Wm,    0, comm=comm_fm)
  end subroutine get_g_spectrum_2d


!-----------------------------------------------!
!> @author  YK
!! @brief   Compressive free-energy power balance scalars (Meyrand 2019):
!!            dW_free/dt = P_g - D_g.
!!          P_g  = sum_k Re[conjg(g_1)*fg]                 (g_1 forcing injection)
!!          D_g  = sum_{k,m} (1 + alpha*delta_{m,0}) * D_op(m,k) * |g_m|^2
!!            D_op = mu_hyper_perp*(kprp2/kprp2_max)^nexp_perp
!!                 + nu_hyper_m*(m/nm)^nexp_m              (mirrors get_imp_terms_tintg_g)
!!          Both carry the kz/=1 factor-of-2 (unstored negative-kz half, RMHD
!!          convention) and NO kprp2 weight (W_free is 1/2|g_m|^2, not the
!!          1/2 kprp2|phi|^2 Alfven energy). The (1+alpha) weight on the m=0 rung
!!          matches the alpha-weighted free energy (eq 6); P_g sees no alpha since
!!          it injects into m=1. Reduced on comm_fm (kperp over comm_fft, m over
!!          comm_m; g is redundant across comm_s) so the total lands on proc0.
!!          dW_free/dt itself is NOT computed here: it is formed in post-processing
!!          as a finite difference of the stored W_free(t) time series.
!-----------------------------------------------!
  subroutine get_g_power_balance(p_g_sum, Dg_sum)
    use mp, only: sum_reduce, comm_fm
    use grid, only: nkx, nky_local, nkz
    use grid, only: nm, nm_local, m_offset
    use grid, only: kprp2, kprp2_max
    use fields, only: g
    use force, only: fg, driven, is_forced
    use params, only: alpha, mu_hyper_perp, nu_hyper_m, nexp_perp, nexp_m
    use utils, only: cabs2
    implicit none
    real(8), intent(out) :: p_g_sum, Dg_sum
    real(8) :: Dsum, Psum, Dop, wt, dloc, ploc
    integer :: i, j, k, mm, m_phys, mm1
    logical :: force_g, owns_m1

    ! --- D_g : free-energy dissipation, summed over all local (k, m) ---
    Dsum = 0.d0
    do mm = 1, nm_local
      m_phys = m_offset + mm - 1
      wt = 1.d0
      if (m_phys == 0) wt = 1.d0 + alpha   ! (1+alpha) weight on the m=0 rung (eq 6)
      !$acc parallel loop collapse(3) gang vector default(present) &
      !$acc& private(Dop, dloc) reduction(+:Dsum)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            Dop = mu_hyper_perp*(kprp2(k, j, i)/kprp2_max)**nexp_perp &
                + nu_hyper_m*(dble(m_phys)/dble(nm))**nexp_m
            dloc = wt*Dop*cabs2(g(k, j, i, mm))
            ! compensate the unstored negative-kz half (RMHD convention)
            if (k /= 1) dloc = 2.0d0*dloc
            Dsum = Dsum + dloc
          end do
        end do
      end do
      !$acc end parallel loop
    end do
    Dg_sum = Dsum

    ! --- P_g : g_1 forcing injection Re[conjg(g_1)*fg] ---
    ! Only the comm_m rank owning global m=1 contributes (fg is allocated there
    ! alone, and only when g1 is forced); the comm_fm reduction below carries the
    ! total to proc0. Guard mirrors advance.F90 solve (force_g .and. owns_m1).
    Psum    = 0.d0
    force_g = driven .and. is_forced('g1')
    mm1     = 2 - m_offset
    owns_m1 = (mm1 >= 1 .and. mm1 <= nm_local)
    if (force_g .and. owns_m1) then
      !$acc parallel loop collapse(3) gang vector default(present) &
      !$acc& private(ploc) reduction(+:Psum)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            ! Re[conjg(g_1)*fg], symmetrized as in force.F90:normalize_force_g
            ploc = dble(0.5d0*( g(k, j, i, mm1)*conjg(fg(k, j, i)) &
                              + conjg(g(k, j, i, mm1))*fg(k, j, i) ))
            if (k /= 1) ploc = 2.0d0*ploc
            Psum = Psum + ploc
          end do
        end do
      end do
      !$acc end parallel loop
    endif
    p_g_sum = Psum

    ! g is m-distributed (comm_m) x kperp-distributed (comm_fft), redundant across
    ! comm_s, so a single comm_fm reduction assembles the full sums on proc0.
    call sum_reduce(Dg_sum , 0, comm=comm_fm)
    call sum_reduce(p_g_sum, 0, comm=comm_fm)
  end subroutine get_g_power_balance


!-----------------------------------------------!
!> @author  YK
!! @brief   Hermite free-energy flux Gamma_m(k_perp) (Meyrand 2019 eq 9):
!!            Gamma_m = sum_{m'=0}^{m} (1 + alpha*delta_{m',0})
!!                                     * v_th * Re[ conjg(g_{m'}) . grad_par S_{m'} ]
!!          binned over k_perp, with the parallel-streaming operator
!!            grad_par S_m = zi*kz*S_m + {Psi, S_m},  S_m = cp*g_{m+1} + cm*g_{m-1}.
!!          Two outputs are produced from the same per-rung source src_m:
!!            gam_bin (k_perp,m) : polar-binned flux spectrum (echo diagnostic);
!!            gam_kint(m)        : k-integrated flux Gamma(m) = sum_k src (all modes).
!!          Telescoping: the LINEAR streaming zi*kz*S_m is k-diagonal, so sum_m src_m
!!          = 0 at EVERY fixed k (weighted antisymmetry w_m*cp_m = w_{m+1}*cm_{m+1}
!!          + ghost-0 boundaries). The {Psi,S_m} bracket, however, redistributes free
!!          energy across k_perp (a perpendicular transfer), so sum_m src_m(k) /= 0
!!          per k -- it cancels only AFTER the k-sum (int a{Psi,b} = -int b{Psi,a}).
!!          Hence gam_kint(nm-1) = sum_{m,k} src_m ~ 0 is the correct telescoping
!!          invariant (= -dW/dt at mu=nu=0); the binned gam_bin top slot is NON-zero
!!          per k for the bracket and must NOT be used as a telescoping check. The
!!          polar binning also drops the k_perp>k_max,1D Cartesian corners, which the
!!          unbinned gam_kint retains -- another reason gam_kint carries the invariant.
!!
!!          The streaming coefficients (cp, cm) replicate advance.F90:s_coeffs
!!          (canonical) inline: the Makefile compiles diagnostics.F90 BEFORE
!!          advance.F90 (no dependency tracking), so `use advance` is unavailable
!!          here. Any drift from the canonical coefficients is caught by the
!!          telescoping (Stage-0) check |gam_kint(nm-1)|/W << 1.
!!
!!          The {Psi, S_m} bracket is gated by `nonlinear` to match the operator
!!          the run actually integrates (get_nonlinear_terms_g is nonlinear-guarded;
!!          the linear -v_th*d_z S_m is unconditional). It is reconstructed with the
!!          code's own grad_par (btran gradients -> real bracket -> ftran /ntot, same
!!          sign/normalization as advance.F90) so diagnostic and dynamics never diverge.
!!
!!          Assembled onto proc0 with the same single comm_fm reduction idiom as
!!          get_g_spectrum_2d (global slot mg = m_offset + mm); the cumulative sum
!!          over m is then done on proc0 (the only rank holding the full spectrum).
!-----------------------------------------------!
  subroutine get_hermite_flux_2d(gam_bin, gam_kint)
    use mp, only: sum_reduce, comm_fm, proc0, halo_exchange_m
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx_local, nly, nlz_padded, ntot
    use grid, only: nm, nm_local, m_offset
    use grid, only: kx, ky, kz
    use fields, only: g, psi
    use params, only: zi, v_th, alpha, nonlinear
    use cuFFTmp, only: btran_c2r, ftran_r2c
    implicit none
    real(8), dimension(1:nkpolar, nm, nkz), intent(out) :: gam_bin  ! device (present)
    real(8), dimension(nm),                 intent(out) :: gam_kint ! host (k-integrated)
    complex(8), allocatable, dimension(:,:,:) :: sk, gx, gy, brk         ! spectral scratch
    real(8),    allocatable, dimension(:,:,:) :: pdx, pdy, sdx, sdy, brr ! real scratch
    real(8) :: k2, src, cp, cm, wt, Ssum
    integer :: i, j, k, mm, mg, m_glob, idx

    ! g_{m+-1} live in the m-halo; refresh it (nproc_m==1 -> early return, ghosts
    ! stay 0 = global Hermite boundary). Safe: solve() re-halos before it reads
    ! ghosts, and get_g_spectrum_2d only touches interior moments.
    call halo_exchange_m(g)

    allocate(sk (nkz, nky_local, nkx))
    allocate(gx (nkz, nky_local, nkx))
    allocate(gy (nkz, nky_local, nkx))
    allocate(brk(nkz, nky_local, nkx))
    allocate(pdx(nlz_padded, nly, nlx_local))
    allocate(pdy(nlz_padded, nly, nlx_local))
    allocate(sdx(nlz_padded, nly, nlx_local))
    allocate(sdy(nlz_padded, nly, nlx_local))
    allocate(brr(nlz_padded, nly, nlx_local))
    !$acc enter data create(sk, gx, gy, brk, pdx, pdy, sdx, sdy, brr)

    ! Zero the binned accumulator (owned + unowned slots) so the comm_fm reduction
    ! assembles the full (k_perp x m) flux spectrum exactly.
    !$acc parallel loop collapse(3) default(present)
    do k = 1, nkz
      do mg = 1, nm
        do idx = 1, nkpolar
          gam_bin(idx, mg, k) = 0.d0
        end do
      end do
    end do
    !$acc end parallel loop
    gam_kint = 0.d0   ! host k-integrated flux accumulator (unowned m slots stay 0)

    ! Zero the bracket once; when nonlinear=F it stays 0 for every moment (grad_par
    ! reduces to zi*kz*S_m), when nonlinear=T it is overwritten per moment below.
    !$acc parallel loop collapse(3) default(present)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          brk(k,j,i) = (0.d0, 0.d0)
        end do
      end do
    end do
    !$acc end parallel loop

    ! Perpendicular gradients of Psi (once), needed only for the {Psi, S_m} bracket.
    if (nonlinear) then
      !$acc parallel loop collapse(3) default(present)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            gx(k,j,i) = zi*kx(i)*psi(k,j,i)
            gy(k,j,i) = zi*ky(j)*psi(k,j,i)
          end do
        end do
      end do
      !$acc end parallel loop
      call btran_c2r(gx, pdx)   ! d_x Psi
      call btran_c2r(gy, pdy)   ! d_y Psi
    endif

    do mm = 1, nm_local
      m_glob = m_offset + mm - 1
      mg     = m_offset + mm          ! global 1-based slot (physical m = mg-1)
      ! streaming stencil coefficients (canonical: advance.F90:s_coeffs)
      cp = sqrt(dble(m_glob + 1)/2.d0)
      cm = sqrt(dble(m_glob    )/2.d0)         ! = 0 at m=0 (boundary g_{-1}=0)
      if (m_glob == 1) cm = cm*(1.d0 + alpha)
      ! per-rung free-energy weight (1 + alpha*delta_{m,0})
      wt = 1.d0
      if (m_glob == 0) wt = 1.d0 + alpha

      ! S_m in Fourier space and its perpendicular gradients (gx=d_x, gy=d_y).
      !$acc parallel loop collapse(3) default(present)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            sk(k,j,i) = cp*g(k,j,i,mm+1) + cm*g(k,j,i,mm-1)
            gx(k,j,i) = zi*kx(i)*sk(k,j,i)
            gy(k,j,i) = zi*ky(j)*sk(k,j,i)
          end do
        end do
      end do
      !$acc end parallel loop

      if (nonlinear) then
        ! Real-space bracket {Psi, S_m} = d_x Psi d_y S_m - d_y Psi d_x S_m,
        ! reconstructed with the code's own grad_par (same sign as advance).
        call btran_c2r(gx, sdx)   ! d_x S_m
        call btran_c2r(gy, sdy)   ! d_y S_m
        !$acc parallel loop collapse(3) default(present)
        do i = 1, nlx_local
          do j = 1, nly
            do k = 1, nlz_padded
              brr(k,j,i) = pdx(k,j,i)*sdy(k,j,i) - pdy(k,j,i)*sdx(k,j,i)
            end do
          end do
        end do
        !$acc end parallel loop
        call ftran_r2c(brr, brk)  ! {Psi, S_m} in Fourier space (needs /ntot below)
      endif

      ! Per-rung source src = v_th * wt * Re[ conjg(g_m) . grad_par S_m ],
      !   grad_par S_m = zi*kz*S_m + {Psi, S_m}/ntot   (bracket = 0 if nonlinear=F),
      ! binned over k_perp into the owned global slot mg. The kz/=1 factor-of-2
      ! matches the W_m convention (compensates the unstored negative-kz half); it
      ! does not break telescoping since sum_m src_m = 0 at every k regardless.
      ! Ssum accumulates src over ALL modes (no polar-bin drop of the k_perp>k_max,1D
      ! Cartesian corners), giving the k-integrated source; unlike the binned
      ! gam_bin it preserves the exact global sum, so sum_m gam_kint(m) = 0 is the
      ! binning-immune telescoping invariant (= -dW/dt with mu=nu=0). The {Psi,S_m}
      ! bracket telescopes only after the k-sum (perp transfer is non-zero per k).
      Ssum = 0.d0
      !$acc parallel loop collapse(3) gang vector default(present) &
      !$acc& private(k2, src, idx) reduction(+:Ssum)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            src = v_th*wt*dble( conjg(g(k,j,i,mm)) &
                                * (zi*kz(k)*sk(k,j,i) + brk(k,j,i)/ntot) )
            if (k /= 1) src = 2.0d0*src
            Ssum = Ssum + src
            k2 = kx(i)**2 + ky(j)**2
            do idx = 1, nkpolar - 1
              if (k2 >= kpbin(idx)**2 .and. k2 < kpbin(idx+1)**2) then
                !$acc atomic update
                gam_bin(idx, mg, k) = gam_bin(idx, mg, k) + src
                exit  ! found the bin
              endif
            enddo
          end do
        end do
      end do
      !$acc end parallel loop
      gam_kint(mg) = Ssum
    end do

    ! Assemble k_perp (comm_fft) and m (comm_m) onto proc0 with one reduction; g is
    ! redundant across comm_s, so comm_fm covers the whole field exactly once. Here
    ! gam_bin holds the PER-RUNG source src_m at each global slot mg.
    !$acc update host(gam_bin)
    call sum_reduce(gam_bin, 0, comm=comm_fm)
    call sum_reduce(gam_kint, 0, comm=comm_fm)  ! assemble k_perp (comm_fft) + m (comm_m)

    ! Cumulative sum over m on proc0 (the only rank with the full spectrum) turns
    ! the per-rung source into the flux Gamma_m = sum_{m'<=m} src_{m'}. The same
    ! prefix sum on the k-integrated source gives Gamma(m); its top slot m=nm-1 is
    ! the global telescoping residual (~0) and equals -dW/dt when mu=nu=0.
    if (proc0) then
      do k = 1, nkz
        do idx = 1, nkpolar
          do mg = 2, nm
            gam_bin(idx, mg, k) = gam_bin(idx, mg, k) + gam_bin(idx, mg-1, k)
          end do
        end do
      end do
      do mg = 2, nm
        gam_kint(mg) = gam_kint(mg) + gam_kint(mg-1)
      end do
    endif

    !$acc exit data delete(sk, gx, gy, brk, pdx, pdy, sdx, sdy, brr)
    deallocate(sk, gx, gy, brk, pdx, pdy, sdx, sdy, brr)
  end subroutine get_hermite_flux_2d


!-----------------------------------------------!
!> @author  YK
!! @brief   Runtime assertion that the redundantly-solved fields (phi, psi)
!!          are bitwise identical across all comm_m/comm_s groups. A fixed
!!          loop order makes the scalar checksum reproducible, so max==min
!!          over a redundant axis holds iff the field data agree. Active only
!!          when nproc_m*nproc_s > 1.
!-----------------------------------------------!
  subroutine check_redundant_consistency
    use MPI
    use mp, only: proc0
    use mp, only: comm_m, comm_s, nproc_m, nproc_s
    use mp, only: max_allreduce, min_allreduce, sum_allreduce
    use grid, only: nkx, nky_local, nkz
    use fields, only: phi, psi
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    integer :: i, j, k, nbad, ierr
    real(8) :: chk, chk_max, chk_min, re, im
    integer(int64) :: h, h_or, h_and

    ! (1) EXACT bitwise-identity test. XOR-fold the raw IEEE bit patterns of every
    ! phi/psi component in a fixed host-side order. XOR is associative AND
    ! commutative, so -- unlike a floating-point sum -- the fold is immune to the
    ! GPU reduction ORDER, which is not guaranteed identical across the distinct
    ! GPUs of a redundant comm_m group. The fold changes iff the field data differ
    ! bit-for-bit, so it tests true redundant-solve identity rather than
    ! reduction-order identity (a plain FP sum yields false ~1 ULP mismatches once
    ! the field grows enough for the summation order to matter).
    !$acc update host(phi, psi)
    h = 0_int64
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          re = dble(phi(k,j,i)); im = aimag(phi(k,j,i))
          h = ieor(h, transfer(re, 0_int64)); h = ieor(h, transfer(im, 0_int64))
          re = dble(psi(k,j,i)); im = aimag(psi(k,j,i))
          h = ieor(h, transfer(re, 0_int64)); h = ieor(h, transfer(im, 0_int64))
        end do
      end do
    end do

    ! (2) floating-point checksum, kept ONLY to report the magnitude of any
    ! disagreement. This sum is order-sensitive and must never drive PASS/FAIL.
    chk = 0.d0
    !$acc data present(phi, psi)
    !$acc parallel loop collapse(3) reduction(+:chk)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          chk = chk + dble(phi(k,j,i)) + aimag(phi(k,j,i)) &
                    + dble(psi(k,j,i)) + aimag(psi(k,j,i))
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
  end subroutine check_redundant_consistency


!-----------------------------------------------!
!> @author  YK
!! @brief   Global XOR-fold checksum of the Hermite-moment field g over the
!!          (P_m x P_fft) plane comm_fm. XOR is associative and commutative, so
!!          folding every local interior g component and reducing with MPI_BXOR
!!          over comm_fm yields a checksum of the WHOLE g field that is
!!          independent of how the Hermite axis m is split across comm_m. Two
!!          runs with identical init/steps but different P_m must therefore print
!!          the same value iff the m-split, m-halo and rank-4 I/O are consistent.
!!          g is redundant across comm_s, so the plane checksum is additionally
!!          required to agree across comm_s (BOR == BAND), never folded over it
!!          (XOR would cancel the identical redundant copies).
!-----------------------------------------------!
  subroutine check_g_consistency
    use MPI
    use mp, only: proc0, comm_fm, comm_s, nproc_s
    use grid, only: nkx, nky_local, nkz, nm_local
    use fields, only: g
    use params, only: v_th, alpha, beta_i, tau, Zcharge, alpha_root
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    integer :: i, j, k, mm, ierr
    real(8) :: re, im
    integer(int64) :: h, h_fm, h_or, h_and
    ! one-time echo of the derived ion-sound coupling alpha for the unit check
    logical, save :: alpha_banner = .true.

    ! Provenance: echo the derived ion-sound coupling once for the unit check.
    if (proc0 .and. alpha_banner) then
      write (*, '("[alpha_check] beta_i tau Zcharge alpha_root = ", &
        &3es24.16e3, i3)') beta_i, tau, Zcharge, alpha_root
      write (*, '("[alpha_check] v_th alpha = ", 2es24.16e3)') v_th, alpha
      alpha_banner = .false.
    endif

    ! g is stored redundantly across comm_s (split over the m x fft plane).
    ! With P_s>1, assert the per-plane XOR checksum is identical on every
    ! s-group; with P_s==1 there is no redundant copy to compare, so skip the
    ! whole thing (no device copy, no reduction, no stdout).
    if (nproc_s > 1) then
      !$acc update host(g)
      h = 0_int64
      do mm = 1, nm_local
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              re = dble(g(k,j,i,mm)); im = aimag(g(k,j,i,mm))
              h = ieor(h, transfer(re, 0_int64)); h = ieor(h, transfer(im, 0_int64))
            end do
          end do
        end do
      end do
      ! P_m-invariant global g checksum over the (m x fft) plane, then compare
      ! it across the redundant s-groups.
      call mpi_allreduce(h, h_fm, 1, MPI_INTEGER8, MPI_BXOR, comm_fm, ierr)
      call mpi_allreduce(h_fm, h_or,  1, MPI_INTEGER8, MPI_BOR,  comm_s, ierr)
      call mpi_allreduce(h_fm, h_and, 1, MPI_INTEGER8, MPI_BAND, comm_s, ierr)
      if (proc0 .and. h_or /= h_and) write (*, &
        '("[check_g_consistency] FAIL: comm_s redundancy mismatch")')
    endif
  end subroutine check_g_consistency


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




