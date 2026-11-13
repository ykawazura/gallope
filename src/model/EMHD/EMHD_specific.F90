!-----------------------------------------------!
!> @author  YK
!! @brief   Model specific subroutines for EMHD
!!
!!          The nonlinear terms live here because both
!!          advance (time stepping) and diagnostics
!!          need them; one copy keeps the two from
!!          diverging.
!!
!!   (1 + de^2 k^2) d_t b = curl(u x b*) + fb
!!   with b* = (1 - de^2 nabla^2) b and the electron
!!   velocity u = -curl b. The guide field is the
!!   k = 0 mode of b, so it enters through b*.
!-----------------------------------------------!
module model_specific
  use fields, only: nfields, ibx, iby, ibz
  implicit none

  public init_model_specific
  public fields_to_real, nonlinear_to_spectral

  private

  ! Fields transformed to real space
  integer, parameter, public :: nbtran = 6
  integer, parameter, public :: ibsx = 1, ibsy = 2, ibsz = 3  !< b*
  integer, parameter, public :: iux  = 4, iuy  = 5, iuz  = 6  !< u = -curl b
  ! Fluxes transformed back to Fourier space
  integer, parameter, public :: nftran = 3

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialize model specific subroutines
!!
!!          EMHD runs on the P_fft decomposition only.
!!          The forcing must act on b directly: the
!!          div-free projection couples bx, by, bz, so
!!          the Elsasser forcing and other field names
!!          have no meaning here.
!-----------------------------------------------!
  subroutine init_model_specific
    use params, only: P_m, P_s
    use force, only: init_force, alloc_force
    use force_common, only: driven, elsasser, field_names
    use mp, only: proc0
    use MPI
    implicit none
    integer :: l, ierr

    if (P_m /= 1 .or. P_s /= 1) then
      if (proc0) print *, 'EMHD supports P_m = P_s = 1 only'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif

    call init_force
    call alloc_force

    if (driven) then
      if (elsasser) then
        if (proc0) print *, 'EMHD: elsasser forcing is not defined'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      endif
      do l = 1, size(field_names)
        if (trim(field_names(l)) /= 'bx' .and. trim(field_names(l)) /= 'by' &
                                         .and. trim(field_names(l)) /= 'bz') then
          if (proc0) print *, 'EMHD: unknown forced field ', trim(field_names(l))
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        endif
      enddo
    endif

  end subroutine init_model_specific


!-----------------------------------------------!
!> @author  YK
!! @brief   b* = (1 + de^2 k^2) b and u = -i k x b
!!          in real space
!!          (w: spectral work, w_r: result)
!-----------------------------------------------!
  subroutine fields_to_real(bx, by, bz, w, w_r)
    use grid, only: kx, ky, kz, k2
    use grid, only: nkx, nky_local, nkz
    use params, only: zi, de2
    use cuFFTmp, only: btran_c2r
    implicit none
    complex(8), dimension(:,:,:)  , intent(in)    :: bx, by, bz
    complex(8), dimension(:,:,:,:), intent(inout) :: w
    real   (8), dimension(:,:,:,:), intent(inout) :: w_r
    real   (8) :: d
    integer :: i, j, k, l

    !$acc data present(bx, by, bz, w)
    !$acc parallel loop collapse(3) private(d)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          d = 1.d0 + de2*k2(k,j,i)
          w(k,j,i,ibsx) = d*bx(k,j,i)
          w(k,j,i,ibsy) = d*by(k,j,i)
          w(k,j,i,ibsz) = d*bz(k,j,i)

          w(k,j,i,iux) = -zi*(ky(j)*bz(k,j,i) - kz(k)*by(k,j,i))
          w(k,j,i,iuy) = -zi*(kz(k)*bx(k,j,i) - kx(i)*bz(k,j,i))
          w(k,j,i,iuz) = -zi*(kx(i)*by(k,j,i) - ky(j)*bx(k,j,i))
        enddo
      enddo
    enddo
    !$acc end data

    do l = 1, nbtran
      call btran_c2r(w(:,:,:,l), w_r(:,:,:,l))
    enddo

  end subroutine fields_to_real


!-----------------------------------------------!
!> @author  YK
!! @brief   Nonlinear terms curl(u x b*) from the
!!          real-space fields, returned in Fourier
!!          space (nonlin_r: real-space work)
!!
!!   flx = b* x u (real space),
!!   N   = -i k x flx = curl(u x b*)
!-----------------------------------------------!
  subroutine nonlinear_to_spectral(w_r, nonlin_r, nonlin)
    use grid, only: nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: kx, ky, kz
    use grid, only: ntot
    use params, only: zi
    use cuFFTmp, only: ftran_r2c
    implicit none
    real   (8), dimension(:,:,:,:), intent(in)    :: w_r
    real   (8), dimension(:,:,:,:), intent(inout) :: nonlin_r
    complex(8), dimension(:,:,:,:), intent(inout) :: nonlin
    complex(8) :: fx, fy, fz
    integer :: i, j, k, l

    !$acc data present(w_r, nonlin_r)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          nonlin_r(k,j,i,ibx) = w_r(k,j,i,ibsy)*w_r(k,j,i,iuz) - w_r(k,j,i,ibsz)*w_r(k,j,i,iuy)
          nonlin_r(k,j,i,iby) = w_r(k,j,i,ibsz)*w_r(k,j,i,iux) - w_r(k,j,i,ibsx)*w_r(k,j,i,iuz)
          nonlin_r(k,j,i,ibz) = w_r(k,j,i,ibsx)*w_r(k,j,i,iuy) - w_r(k,j,i,ibsy)*w_r(k,j,i,iux)
        enddo
      enddo
    enddo
    !$acc end data

    do l = 1, nftran
      call ftran_r2c(nonlin_r(:,:,:,l), nonlin(:,:,:,l))
    enddo

    !$acc data present(nonlin)
    !$acc parallel loop collapse(3) private(fx, fy, fz)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fx = nonlin(k,j,i,ibx)/ntot
          fy = nonlin(k,j,i,iby)/ntot
          fz = nonlin(k,j,i,ibz)/ntot
          nonlin(k,j,i,ibx) = -zi*(ky(j)*fz - kz(k)*fy)
          nonlin(k,j,i,iby) = -zi*(kz(k)*fx - kx(i)*fz)
          nonlin(k,j,i,ibz) = -zi*(kx(i)*fy - ky(j)*fx)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine nonlinear_to_spectral

end module model_specific
