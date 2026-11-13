!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics for RRMHD
!!
!!          Energy (conserved for q = kappa_p = 0 without
!!          dissipation; kappa_b is energy-neutral):
!!            E = 1/2 <|grad phi|^2 + |grad psi|^2 + upa^2 + A bpa^2>,
!!            A = va2cs2_plus_1 = 1 + 2/(beta*gamma)
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
!!          Energies and Elsasser energies are taken at the
!!          current time; the dissipation and the linear
!!          sources at the
!!          midpoint of the last step, so that
!!            d/dt E = p_aw_rot + p_compr_rot
!!                   + p_aw_grd + p_compr_grd - dissipation
!!          holds to O(dt^2) (RRMHD has no forcing).
!!
!!          Linear sources (the rotation terms with q and
!!          the background-gradient terms with kappa_b,
!!          kappa_p), in explicit real form; the subscript m
!!          is the midpoint and c = cs2va2:
!!            p_aw_rot    =  2 ky Im(phi_m conj(upa_m))
!!            p_compr_rot = -ky [ q Im(psi_m conj(bpa_m))
!!                              + (2 - q) Im(phi_m conj(upa_m)) ]
!!            p_aw_grd    = -2 ky (kappa_b + c kappa_p) Im(phi_m conj(bpa_m))
!!            p_compr_grd =  ky [ (2 kappa_b + (c - 1) kappa_p) Im(phi_m conj(bpa_m))
!!                              + c kappa_p Im(psi_m conj(upa_m)) ]
!!
!!          Elsasser energies at the current time without a
!!          factor 1/2, as in Calliope RRMHD:
!!            zpep2 = kprp2 |phi + psi|^2
!!            zpem2 = kprp2 |phi - psi|^2
!!            zpap2 = |upa + sqrt(A) bpa|^2
!!            zpam2 = |upa - sqrt(A) bpa|^2
!!          so that, to round-off, in the series and per bin,
!!            zpep2 + zpem2 = 4 (upe2 + bpe2)
!!            zpap2 + zpam2 = 4 (upa2 + bpa2)
!!          (bpa2 includes the factor A).
!-----------------------------------------------!
  subroutine loop_diagnostics
    use io, only: loop_io, nseries, nbins
    use io, only: is_upe2, is_bpe2, is_upa2, is_bpa2
    use io, only: is_upe2dot, is_bpe2dot, is_upa2dot, is_bpa2dot
    use io, only: is_upe2dissip, is_bpe2dissip, is_upa2dissip, is_bpa2dissip
    use io, only: is_p_aw_rot, is_p_compr_rot, is_p_aw_grd, is_p_compr_grd
    use io, only: is_zpep2, is_zpem2, is_zpap2, is_zpam2
    use io, only: ib_upe2, ib_bpe2, ib_upa2, ib_bpa2
    use io, only: ib_ux2, ib_uy2, ib_bx2, ib_by2
    use io, only: ib_p_aw_rot, ib_p_compr_rot, ib_p_aw_grd, ib_p_compr_grd
    use io, only: ib_dissip_aw, ib_dissip_compr
    use io, only: ib_zpep2, ib_zpem2, ib_zpap2, ib_zpam2
    use mp, only: proc0
    use grid, only: kprp2, kz2, kprp2_max, kz2_max
    use grid, only: kx, ky
    use grid, only: nkx, nky_local, nkz
    use fields, only: phi, psi, upa, bpa
    use fields, only: phi_old, psi_old, upa_old, bpa_old
    use mp, only: sum_reduce, comm_fft, nproc_m, nproc_s
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: q, kappa_b, kappa_p, va2cs2_plus_1, cs2va2, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp, &
                      nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                      etapa_x, etapa_x_exp, etapa_z, etapa_z_exp
    use utils, only: cabs2
    implicit none
    ! every bin is computed per mode
    integer, parameter :: nq = nbins
    integer :: i, j, k, l

    real(8) :: series(nseries)
    real(8), allocatable, dimension(:,:,:)   :: bins   ! (kpbin, kz, bin) on host
    real(8), allocatable, dimension(:,:,:,:) :: qq     ! (kz, ky, kx, bin) per-mode work
    real(8), allocatable, dimension(:,:)     :: ebin   ! (kpbin, kz)

    ! coefficients, copied to the device as firstprivate scalars
    real(8)    :: aa, sqa, c_rot_u, c_aw_grd, c_compr_grd_b, c_compr_grd_u

    ! per-mode scalars
    real(8)    :: fac, kp2, kyj
    real(8)    :: rx_upe, rz_upe, rx_bpe, rz_bpe, rx_upa, rz_upa, rx_bpa, rz_bpa
    real(8)    :: e_upe, e_bpe, e_upa, e_bpa, eo_upe, eo_bpe, eo_upa, eo_bpa
    real(8)    :: dx_upe, dz_upe, dx_bpe, dz_bpe, dx_upa, dz_upa, dx_bpa, dz_bpa
    real(8)    :: im_phi_upa, im_psi_bpa, im_phi_bpa, im_psi_upa
    real(8)    :: p_aw_rot, p_compr_rot, p_aw_grd, p_compr_grd
    real(8)    :: zpep2, zpem2, zpap2, zpam2
    complex(8) :: phi_m, psi_m, upa_m, bpa_m

    ! reduction variables
    real(8) :: upe2_sum, bpe2_sum, upa2_sum, bpa2_sum
    real(8) :: upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum
    real(8) :: upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum
    real(8) :: p_aw_rot_sum, p_compr_rot_sum, p_aw_grd_sum, p_compr_grd_sum
    real(8) :: zpep2_sum, zpem2_sum, zpap2_sum, zpam2_sum

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    aa            = va2cs2_plus_1
    sqa           = sqrt(va2cs2_plus_1)
    c_rot_u       = 2.d0 - q
    c_aw_grd      = 2.d0*(kappa_b + cs2va2*kappa_p)
    c_compr_grd_b = 2.d0*kappa_b + (cs2va2 - 1.d0)*kappa_p
    c_compr_grd_u = cs2va2*kappa_p

    allocate(bins(1:nkpolar, nkz, nbins), source=0.d0)
    allocate(qq(nkz, nky_local, nkx, nq))
    allocate(ebin(1:nkpolar, nkz))
    !$acc enter data create(qq, ebin)

    upe2_sum       = 0.d0; bpe2_sum        = 0.d0; upa2_sum     = 0.d0; bpa2_sum        = 0.d0
    upe2dot_sum    = 0.d0; bpe2dot_sum     = 0.d0; upa2dot_sum  = 0.d0; bpa2dot_sum     = 0.d0
    upe2dissip_sum = 0.d0; bpe2dissip_sum  = 0.d0; upa2dissip_sum = 0.d0; bpa2dissip_sum = 0.d0
    p_aw_rot_sum   = 0.d0; p_compr_rot_sum = 0.d0; p_aw_grd_sum = 0.d0; p_compr_grd_sum = 0.d0
    zpep2_sum      = 0.d0; zpem2_sum       = 0.d0; zpap2_sum    = 0.d0; zpam2_sum       = 0.d0

    !$acc parallel loop collapse(3) gang vector default(present) &
    !$acc& private(fac, kp2, kyj, rx_upe, rz_upe, rx_bpe, rz_bpe, rx_upa, rz_upa, rx_bpa, rz_bpa) &
    !$acc& private(e_upe, e_bpe, e_upa, e_bpa, eo_upe, eo_bpe, eo_upa, eo_bpa) &
    !$acc& private(dx_upe, dz_upe, dx_bpe, dz_bpe, dx_upa, dz_upa, dx_bpa, dz_bpa) &
    !$acc& private(im_phi_upa, im_psi_bpa, im_phi_bpa, im_psi_upa) &
    !$acc& private(p_aw_rot, p_compr_rot, p_aw_grd, p_compr_grd) &
    !$acc& private(zpep2, zpem2, zpap2, zpam2) &
    !$acc& private(phi_m, psi_m, upa_m, bpa_m) &
    !$acc& reduction(+:upe2_sum, bpe2_sum, upa2_sum, bpa2_sum) &
    !$acc& reduction(+:upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum) &
    !$acc& reduction(+:upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum) &
    !$acc& reduction(+:p_aw_rot_sum, p_compr_rot_sum, p_aw_grd_sum, p_compr_grd_sum) &
    !$acc& reduction(+:zpep2_sum, zpem2_sum, zpap2_sum, zpam2_sum)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ! kz /= 0 modes stand for their kz < 0 conjugates as well
          fac = 1.d0
          if (k /= 1) fac = 2.d0
          kp2 = kprp2(k, j, i)
          kyj = ky(j)

          ! Energy-dissipation coefficients, the same as advance::get_imp_rates
          ! up to the 1/A of bpa: the bpa rate there is r/A and the bpa energy
          ! weight is A, so the bpa loss is r |bpa|^2 without A.
          rx_upe = nupe_x *(kp2/kprp2_max)**nupe_x_exp ; rz_upe = nupe_z *(kz2(k)/kz2_max)**nupe_z_exp
          rx_bpe = etape_x*(kp2/kprp2_max)**etape_x_exp; rz_bpe = etape_z*(kz2(k)/kz2_max)**etape_z_exp
          rx_upa = nupa_x *(kp2/kprp2_max)**nupa_x_exp ; rz_upa = nupa_z *(kz2(k)/kz2_max)**nupa_z_exp
          rx_bpa = etapa_x*(kp2/kprp2_max)**etapa_x_exp; rz_bpa = etapa_z*(kz2(k)/kz2_max)**etapa_z_exp

          e_upe  = 0.5d0*cabs2(phi    (k, j, i))*kp2
          e_bpe  = 0.5d0*cabs2(psi    (k, j, i))*kp2
          e_upa  = 0.5d0*cabs2(upa    (k, j, i))
          e_bpa  = 0.5d0*cabs2(bpa    (k, j, i))*aa
          eo_upe = 0.5d0*cabs2(phi_old(k, j, i))*kp2
          eo_bpe = 0.5d0*cabs2(psi_old(k, j, i))*kp2
          eo_upa = 0.5d0*cabs2(upa_old(k, j, i))
          eo_bpa = 0.5d0*cabs2(bpa_old(k, j, i))*aa

          phi_m  = 0.5d0*(phi(k, j, i) + phi_old(k, j, i))
          psi_m  = 0.5d0*(psi(k, j, i) + psi_old(k, j, i))
          upa_m  = 0.5d0*(upa(k, j, i) + upa_old(k, j, i))
          bpa_m  = 0.5d0*(bpa(k, j, i) + bpa_old(k, j, i))

          ! no factor 2: d/dt (1/2 |f|^2) = -r |f|^2
          dx_upe = rx_upe*cabs2(phi_m)*kp2; dz_upe = rz_upe*cabs2(phi_m)*kp2
          dx_bpe = rx_bpe*cabs2(psi_m)*kp2; dz_bpe = rz_bpe*cabs2(psi_m)*kp2
          dx_upa = rx_upa*cabs2(upa_m)    ; dz_upa = rz_upa*cabs2(upa_m)
          dx_bpa = rx_bpa*cabs2(bpa_m)    ; dz_bpa = rz_bpa*cabs2(bpa_m)

          im_phi_upa = aimag(phi_m*conjg(upa_m))
          im_psi_bpa = aimag(psi_m*conjg(bpa_m))
          im_phi_bpa = aimag(phi_m*conjg(bpa_m))
          im_psi_upa = aimag(psi_m*conjg(upa_m))

          p_aw_rot    =  2.d0*kyj*im_phi_upa
          p_compr_rot = -kyj*(q*im_psi_bpa + c_rot_u*im_phi_upa)
          p_aw_grd    = -kyj*c_aw_grd*im_phi_bpa
          p_compr_grd =  kyj*(c_compr_grd_b*im_phi_bpa + c_compr_grd_u*im_psi_upa)

          ! state quantities: current fields, no 1/2 (see the header)
          zpep2 = cabs2(phi(k, j, i) + psi(k, j, i))*kp2
          zpem2 = cabs2(phi(k, j, i) - psi(k, j, i))*kp2
          zpap2 = cabs2(upa(k, j, i) + sqa*bpa(k, j, i))
          zpam2 = cabs2(upa(k, j, i) - sqa*bpa(k, j, i))

          upe2_sum        = upe2_sum        + fac*e_upe
          bpe2_sum        = bpe2_sum        + fac*e_bpe
          upa2_sum        = upa2_sum        + fac*e_upa
          bpa2_sum        = bpa2_sum        + fac*e_bpa
          upe2dot_sum     = upe2dot_sum     + fac*(e_upe - eo_upe)/dt
          bpe2dot_sum     = bpe2dot_sum     + fac*(e_bpe - eo_bpe)/dt
          upa2dot_sum     = upa2dot_sum     + fac*(e_upa - eo_upa)/dt
          bpa2dot_sum     = bpa2dot_sum     + fac*(e_bpa - eo_bpa)/dt
          upe2dissip_sum  = upe2dissip_sum  + fac*(dx_upe + dz_upe)
          bpe2dissip_sum  = bpe2dissip_sum  + fac*(dx_bpe + dz_bpe)
          upa2dissip_sum  = upa2dissip_sum  + fac*(dx_upa + dz_upa)
          bpa2dissip_sum  = bpa2dissip_sum  + fac*(dx_bpa + dz_bpa)
          p_aw_rot_sum    = p_aw_rot_sum    + fac*p_aw_rot
          p_compr_rot_sum = p_compr_rot_sum + fac*p_compr_rot
          p_aw_grd_sum    = p_aw_grd_sum    + fac*p_aw_grd
          p_compr_grd_sum = p_compr_grd_sum + fac*p_compr_grd
          zpep2_sum       = zpep2_sum       + fac*zpep2
          zpem2_sum       = zpem2_sum       + fac*zpem2
          zpap2_sum       = zpap2_sum       + fac*zpap2
          zpam2_sum       = zpam2_sum       + fac*zpam2

          qq(k, j, i, ib_upe2)         = fac*e_upe
          qq(k, j, i, ib_bpe2)         = fac*e_bpe
          qq(k, j, i, ib_upa2)         = fac*e_upa
          qq(k, j, i, ib_bpa2)         = fac*e_bpa
          qq(k, j, i, ib_ux2)          = fac*0.5d0*kyj**2  *cabs2(phi(k, j, i))
          qq(k, j, i, ib_uy2)          = fac*0.5d0*kx(i)**2*cabs2(phi(k, j, i))
          qq(k, j, i, ib_bx2)          = fac*0.5d0*kyj**2  *cabs2(psi(k, j, i))
          qq(k, j, i, ib_by2)          = fac*0.5d0*kx(i)**2*cabs2(psi(k, j, i))
          qq(k, j, i, ib_p_aw_rot)     = fac*p_aw_rot
          qq(k, j, i, ib_p_compr_rot)  = fac*p_compr_rot
          qq(k, j, i, ib_p_aw_grd)     = fac*p_aw_grd
          qq(k, j, i, ib_p_compr_grd)  = fac*p_compr_grd
          qq(k, j, i, ib_dissip_aw)    = fac*(dx_upe + dz_upe + dx_bpe + dz_bpe)
          qq(k, j, i, ib_dissip_compr) = fac*(dx_upa + dz_upa + dx_bpa + dz_bpa)
          qq(k, j, i, ib_zpep2)        = fac*zpep2
          qq(k, j, i, ib_zpem2)        = fac*zpem2
          qq(k, j, i, ib_zpap2)        = fac*zpap2
          qq(k, j, i, ib_zpam2)        = fac*zpam2
        end do
      end do
    end do
    !$acc end parallel loop

    series(is_upe2       ) = upe2_sum
    series(is_bpe2       ) = bpe2_sum
    series(is_upa2       ) = upa2_sum
    series(is_bpa2       ) = bpa2_sum
    series(is_upe2dot    ) = upe2dot_sum
    series(is_bpe2dot    ) = bpe2dot_sum
    series(is_upa2dot    ) = upa2dot_sum
    series(is_bpa2dot    ) = bpa2dot_sum
    series(is_upe2dissip ) = upe2dissip_sum
    series(is_bpe2dissip ) = bpe2dissip_sum
    series(is_upa2dissip ) = upa2dissip_sum
    series(is_bpa2dissip ) = bpa2dissip_sum
    series(is_p_aw_rot   ) = p_aw_rot_sum
    series(is_p_compr_rot) = p_compr_rot_sum
    series(is_p_aw_grd   ) = p_aw_grd_sum
    series(is_p_compr_grd) = p_compr_grd_sum
    series(is_zpep2      ) = zpep2_sum
    series(is_zpem2      ) = zpem2_sum
    series(is_zpap2      ) = zpap2_sum
    series(is_zpam2      ) = zpam2_sum

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
