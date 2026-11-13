!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics for EMHD
!!
!!          Invariants (ideal, unforced):
!!            W  = 1/2 <(1 + de^2 k^2)|b|^2> = 1/2 <b.b*>
!!            H* = <A*.b*>, b* = (1 - de^2 lap) b, curl A* = b*
!!          W is conserved with the guide field; H* only
!!          for b0 = 0 (the uniform field has no vector
!!          potential in a periodic box).
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
    use diagnostics_common, only: init_polar_spectrum_3d
    use io, only: init_io
    implicit none

    call read_parameters(inputfile)

    call init_polar_spectrum_3d
    call init_io(nkpolar, kpbin)
  end subroutine init_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @brief   Diagnostics in loop
!!
!!          W and H* are taken at the current time;
!!          dissipation and injection at the midpoint
!!          of the last step, so that
!!            d/dt W = p_ext - dissipation
!!          holds to O(dt^2). The guide field is the
!!          real part of the k = 0 mode (Calliope
!!          writes its absolute value).
!-----------------------------------------------!
  subroutine loop_diagnostics
    use io, only: loop_io, nseries, nbins
    use io, only: is_wmag, is_wmag_dot, is_wmag_dissip, is_p_ext, is_hel, is_divb_rel
    use io, only: ib_b2, ib_bx2, ib_by2, ib_bz2, ib_wmag
    use mp, only: proc0, iproc_fft
    use mp, only: sum_reduce, comm_fft
    use grid, only: kx, ky, kz, k2, k2inv, k2_max
    use grid, only: nkx, nky_local, nkz
    use fields, only: bx, by, bz
    use fields, only: bx_old, by_old, bz_old
    use fields, only: div_b_rel
    use force, only: fbx, fby, fbz
    use force, only: fbx_old, fby_old, fbz_old
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, eta, eta_exp, de2
    use utils, only: cabs2
    implicit none
    integer :: i, j, k, l

    real(8) :: series(nseries)
    real(8) :: b0vec(3)
    real(8) :: divb
    real(8), allocatable, dimension(:,:)     :: bins   ! (kpbin, bin) on host
    real(8), allocatable, dimension(:,:,:,:) :: qq     ! (kz, ky, kx, bin) per-mode work
    real(8), allocatable, dimension(:)       :: ebin   ! (kpbin)

    ! per-mode scalars
    real(8)    :: fac, d, r, w, wo
    complex(8) :: bx_m, by_m, bz_m, fbx_m, fby_m, fbz_m
    complex(8) :: ax, ay, az

    ! reduction variables
    real(8) :: wmag_sum, wmag_dot_sum, wmag_dissip_sum, p_ext_sum, hel_sum

    ! div_b_rel is collective: every rank calls it
    divb = div_b_rel(bx, by, bz)

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(bins(1:nkpolar, nbins), source=0.d0)
    allocate(qq(nkz, nky_local, nkx, nbins))
    allocate(ebin(1:nkpolar))
    !$acc enter data create(qq, ebin)

    wmag_sum = 0.d0; wmag_dot_sum = 0.d0; wmag_dissip_sum = 0.d0
    p_ext_sum = 0.d0; hel_sum = 0.d0

    !$acc parallel loop collapse(3) gang vector default(present) &
    !$acc& private(fac, d, r, w, wo, bx_m, by_m, bz_m, fbx_m, fby_m, fbz_m, ax, ay, az) &
    !$acc& reduction(+:wmag_sum, wmag_dot_sum, wmag_dissip_sum, p_ext_sum, hel_sum)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ! kz /= 0 modes stand for their kz < 0 conjugates as well
          fac = 1.d0
          if (k /= 1) fac = 2.d0

          d = 1.d0 + de2*k2(k, j, i)
          ! energy-dissipation coefficient: advance::imp_rate times
          ! (1 + de^2 k^2), since dW/dt = -sum r |b|^2. Zero at k = 0
          ! so that eta_exp = 0 does not dissipate the guide field.
          r = 0.d0
          if (k2(k, j, i) > 0.d0) r = eta*(k2(k, j, i)/k2_max)**eta_exp

          w  = 0.5d0*d*(cabs2(bx    (k, j, i)) + cabs2(by    (k, j, i)) + cabs2(bz    (k, j, i)))
          wo = 0.5d0*d*(cabs2(bx_old(k, j, i)) + cabs2(by_old(k, j, i)) + cabs2(bz_old(k, j, i)))

          bx_m  = 0.5d0*(bx (k, j, i) + bx_old (k, j, i))
          by_m  = 0.5d0*(by (k, j, i) + by_old (k, j, i))
          bz_m  = 0.5d0*(bz (k, j, i) + bz_old (k, j, i))
          fbx_m = 0.5d0*(fbx(k, j, i) + fbx_old(k, j, i))
          fby_m = 0.5d0*(fby(k, j, i) + fby_old(k, j, i))
          fbz_m = 0.5d0*(fbz(k, j, i) + fbz_old(k, j, i))

          ! A = i k x b/k^2 (zero at k = 0)
          ax = zi*(ky(j)*bz(k, j, i) - kz(k)*by(k, j, i))*k2inv(k, j, i)
          ay = zi*(kz(k)*bx(k, j, i) - kx(i)*bz(k, j, i))*k2inv(k, j, i)
          az = zi*(kx(i)*by(k, j, i) - ky(j)*bx(k, j, i))*k2inv(k, j, i)

          wmag_sum        = wmag_sum        + fac*w
          wmag_dot_sum    = wmag_dot_sum    + fac*(w - wo)/dt
          wmag_dissip_sum = wmag_dissip_sum + fac*r*(cabs2(bx_m) + cabs2(by_m) + cabs2(bz_m))
          p_ext_sum       = p_ext_sum       + fac*( dble(fbx_m*conjg(bx_m)) &
                                                  + dble(fby_m*conjg(by_m)) &
                                                  + dble(fbz_m*conjg(bz_m)) )
          hel_sum         = hel_sum         + fac*d*d*dble( ax*conjg(bx(k, j, i)) &
                                                          + ay*conjg(by(k, j, i)) &
                                                          + az*conjg(bz(k, j, i)) )

          qq(k, j, i, ib_bx2 ) = fac*0.5d0*cabs2(bx(k, j, i))
          qq(k, j, i, ib_by2 ) = fac*0.5d0*cabs2(by(k, j, i))
          qq(k, j, i, ib_bz2 ) = fac*0.5d0*cabs2(bz(k, j, i))
          qq(k, j, i, ib_b2  ) = qq(k, j, i, ib_bx2) + qq(k, j, i, ib_by2) + qq(k, j, i, ib_bz2)
          qq(k, j, i, ib_wmag) = fac*w
        end do
      end do
    end do
    !$acc end parallel loop

    series(is_wmag       ) = wmag_sum
    series(is_wmag_dot   ) = wmag_dot_sum
    series(is_wmag_dissip) = wmag_dissip_sum
    series(is_p_ext      ) = p_ext_sum
    series(is_hel        ) = hel_sum
    series(is_divb_rel   ) = 0.d0

    ! Spectral partial sums: reduce on comm_fft; the global proc0 is the
    ! root of comm_fft group 0, so the result lands where loop_io writes.
    call sum_reduce(series, 0, comm=comm_fft)
    ! already reduced (max over all ranks)
    series(is_divb_rel) = divb

    ! The k = 0 mode lives on iproc_fft = 0, which holds the global proc0
    b0vec = 0.d0
    if (iproc_fft == 0) then
      !$acc serial present(bx, by, bz) copyout(b0vec)
      b0vec(1) = dble(bx(1, 1, 1))
      b0vec(2) = dble(by(1, 1, 1))
      b0vec(3) = dble(bz(1, 1, 1))
      !$acc end serial
    endif

    ! Bin over |k| (each call reduces on comm_fft; the host copy of ebin
    ! holds the reduced bins on proc0)
    do l = 1, nbins
      call get_polar_spectrum_3d(qq(:,:,:,l), ebin)
      bins(:,l) = ebin
    enddo

    !$acc exit data delete(qq, ebin)
    deallocate(qq, ebin)

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io(series, b0vec, nkpolar, bins)

    deallocate(bins)
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @brief   2D cuts of b and j = curl b
!-----------------------------------------------!
  subroutine loop_diagnostics_2D
    use io, only: loop_io_2D, n2d
    use io, only: i2d_bx, i2d_by, i2d_bz, i2d_jx, i2d_jy, i2d_jz
    use mp, only: proc0
    use grid, only: nlx_local, nly, nlz
    use grid, only: nlz_padded
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use fields, only: bx, by, bz
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
      !$acc parallel loop collapse(3) present(bx, by, bz, f, kx, ky, kz)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            if      (l == i2d_bx) then
              f(k, j, i) = bx(k, j, i)
            else if (l == i2d_by) then
              f(k, j, i) = by(k, j, i)
            else if (l == i2d_bz) then
              f(k, j, i) = bz(k, j, i)
            else if (l == i2d_jx) then
              f(k, j, i) = zi*(ky(j)*bz(k, j, i) - kz(k)*by(k, j, i))
            else if (l == i2d_jy) then
              f(k, j, i) = zi*(kz(k)*bx(k, j, i) - kx(i)*bz(k, j, i))
            else if (l == i2d_jz) then
              f(k, j, i) = zi*(kx(i)*by(k, j, i) - ky(j)*bx(k, j, i))
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
