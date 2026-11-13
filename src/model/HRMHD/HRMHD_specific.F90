!-----------------------------------------------!
!> @author  YK
!! @brief   Model specific subroutines
!!
!!          The nonlinear terms live here because both
!!          advance (time stepping) and diagnostics
!!          (transfer functions of W and X) need them;
!!          one copy keeps the two from diverging.
!-----------------------------------------------!
module model_specific
  use force_common
  use fields, only: nfields, iomg, ipsi, iupa, ibpa
  implicit none

  public init_model_specific
  public gradients_to_real, nonlinear_to_spectral

  ! Gradients transformed to real space
  integer, parameter, public :: nbtran = 12
  integer, parameter, public :: idphi_dx =  1, idphi_dy =  2
  integer, parameter, public :: idomg_dx =  3, idomg_dy =  4
  integer, parameter, public :: idpsi_dx =  5, idpsi_dy =  6
  integer, parameter, public :: idjpa_dx =  7, idjpa_dy =  8
  integer, parameter, public :: idupa_dx =  9, idupa_dy = 10
  integer, parameter, public :: idbpa_dx = 11, idbpa_dy = 12

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialize model specific subroutines
!-----------------------------------------------!
  subroutine init_model_specific
    use force, only :init_force, alloc_force
    implicit none

    call init_force
    call alloc_force

  end subroutine init_model_specific


!-----------------------------------------------!
!> @author  YK
!! @brief   x and y derivatives of phi, omg, psi,
!!          jpa, upa, bpa in real space
!!          (w: spectral work, w_r: result)
!-----------------------------------------------!
  subroutine gradients_to_real(phi, psi, upa, bpa, w, w_r)
    use grid, only: kprp2, kx, ky
    use grid, only: nkx, nky_local, nkz
    use params, only: zi
    use cuFFTmp, only: btran_c2r
    implicit none
    complex(8), dimension(:,:,:)  , intent(in)    :: phi, psi, upa, bpa
    complex(8), dimension(:,:,:,:), intent(inout) :: w
    real   (8), dimension(:,:,:,:), intent(inout) :: w_r
    integer :: i, j, k, l

    !$acc data present(phi, psi, upa, bpa, w)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          w(k,j,i,idphi_dx) = zi*kx(i)                *phi(k,j,i)
          w(k,j,i,idphi_dy) = zi*ky(j)                *phi(k,j,i)
          w(k,j,i,idomg_dx) = zi*kx(i)*(-kprp2(k,j,i))*phi(k,j,i)
          w(k,j,i,idomg_dy) = zi*ky(j)*(-kprp2(k,j,i))*phi(k,j,i)

          w(k,j,i,idpsi_dx) = zi*kx(i)                *psi(k,j,i)
          w(k,j,i,idpsi_dy) = zi*ky(j)                *psi(k,j,i)
          w(k,j,i,idjpa_dx) = zi*kx(i)*(-kprp2(k,j,i))*psi(k,j,i)
          w(k,j,i,idjpa_dy) = zi*ky(j)*(-kprp2(k,j,i))*psi(k,j,i)

          w(k,j,i,idupa_dx) = zi*kx(i)                *upa(k,j,i)
          w(k,j,i,idupa_dy) = zi*ky(j)                *upa(k,j,i)
          w(k,j,i,idbpa_dx) = zi*kx(i)                *bpa(k,j,i)
          w(k,j,i,idbpa_dy) = zi*ky(j)                *bpa(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

    do l = 1, nbtran
      call btran_c2r(w(:,:,:,l), w_r(:,:,:,l))
    enddo

  end subroutine gradients_to_real


!-----------------------------------------------!
!> @author  YK
!! @brief   Nonlinear terms from the real-space
!!          gradients, returned in Fourier space
!!          (nonlin_r: real-space work)
!!
!!   N_omg = -{phi,omg} + {psi,jpa}
!!   N_psi =  {psi, phi + rho*bpa}
!!   N_upa = -{phi,upa} + sgm*{psi,bpa}
!!   N_bpa = -{phi,bpa} + {psi, sgm*upa - rho*jpa}
!!
!!   with {f,g} = f_x g_y - f_y g_x
!-----------------------------------------------!
  subroutine nonlinear_to_spectral(w_r, nonlin_r, nonlin)
    use grid, only: nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use params, only: rho, sgm
    use cuFFTmp, only: ftran_r2c
    implicit none
    real   (8), dimension(:,:,:,:), intent(in)    :: w_r
    real   (8), dimension(:,:,:,:), intent(inout) :: nonlin_r
    complex(8), dimension(:,:,:,:), intent(inout) :: nonlin
    integer :: i, j, k, l

    !$acc data present(w_r, nonlin_r)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          nonlin_r(k,j,i,iomg) = - w_r(k,j,i,idphi_dx)*w_r(k,j,i,idomg_dy) &
                                 + w_r(k,j,i,idphi_dy)*w_r(k,j,i,idomg_dx) &
                                 + w_r(k,j,i,idpsi_dx)*w_r(k,j,i,idjpa_dy) &
                                 - w_r(k,j,i,idpsi_dy)*w_r(k,j,i,idjpa_dx)
          nonlin_r(k,j,i,ipsi) =   w_r(k,j,i,idpsi_dx)*(w_r(k,j,i,idphi_dy) + rho*w_r(k,j,i,idbpa_dy)) &
                                 - w_r(k,j,i,idpsi_dy)*(w_r(k,j,i,idphi_dx) + rho*w_r(k,j,i,idbpa_dx))
          nonlin_r(k,j,i,iupa) = - w_r(k,j,i,idphi_dx)*w_r(k,j,i,idupa_dy) &
                                 + w_r(k,j,i,idphi_dy)*w_r(k,j,i,idupa_dx) &
                                 + sgm*( w_r(k,j,i,idpsi_dx)*w_r(k,j,i,idbpa_dy) &
                                       - w_r(k,j,i,idpsi_dy)*w_r(k,j,i,idbpa_dx) )
          nonlin_r(k,j,i,ibpa) = - w_r(k,j,i,idphi_dx)*w_r(k,j,i,idbpa_dy) &
                                 + w_r(k,j,i,idphi_dy)*w_r(k,j,i,idbpa_dx) &
                                 + w_r(k,j,i,idpsi_dx)*(sgm*w_r(k,j,i,idupa_dy) - rho*w_r(k,j,i,idjpa_dy)) &
                                 - w_r(k,j,i,idpsi_dy)*(sgm*w_r(k,j,i,idupa_dx) - rho*w_r(k,j,i,idjpa_dx))
        enddo
      enddo
    enddo
    !$acc end data

    do l = 1, nfields
      call ftran_r2c(nonlin_r(:,:,:,l), nonlin(:,:,:,l))
    enddo

    !$acc data present(nonlin)
    !$acc parallel loop collapse(4)
    do l = 1, nfields
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            nonlin(k,j,i,l) = nonlin(k,j,i,l)/ntot
          enddo
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine nonlinear_to_spectral

end module model_specific
