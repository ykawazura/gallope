!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../force_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Forcing setting specific to RMHD
!!          'force_common' is inherited
!-----------------------------------------------!
module force
  use force_common
  implicit none

  public get_force

  complex(8), dimension(:,:,:), allocatable :: fphi, fphi_old
  complex(8), dimension(:,:,:), allocatable :: fpsi, fpsi_old
  complex(8), dimension(:,:,:), allocatable :: fzppe
  complex(8), dimension(:,:,:), allocatable :: fzmpe

  !$acc declare create(fphi, fphi_old)
  !$acc declare create(fpsi, fpsi_old)
  !$acc declare create(fzppe)
  !$acc declare create(fzmpe)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialize model specific subroutines
!-----------------------------------------------!
  subroutine alloc_force
    use grid, only: nkx, nky_local, nkz
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    allocate(src(nkz, nky_local, nkx), source=(0.d0,0.d0))
    allocate(fphi    , source=src)
    allocate(fphi_old, source=src)
    allocate(fpsi    , source=src)
    allocate(fpsi_old, source=src)
    if(elsasser) then
      allocate(fzppe, source=src)
      allocate(fzmpe, source=src)
    endif
    !$acc update device(fphi, fphi_old)
    !$acc update device(fpsi, fpsi_old)
    !$acc update device(fzppe)
    !$acc update device(fzmpe)
    deallocate(src)

  end subroutine alloc_force


!-----------------------------------------------!
!> @author  YK
!! @brief   Compute the power of forcing 
!!          and normalize with it 
!-----------------------------------------------!
  subroutine normalize_force(fphi, fpsi)
    use grid, only: nkx, nky_local, nkz, kprp2
    use fields, only: phi, psi
    use mp, only: sum_allreduce
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: fphi, fpsi
    real(8) :: phi_dot_nbl2_fphi_sum, psi_dot_nbl2_fpsi_sum
    real(8), dimension(:,:,:), allocatable :: phi_dot_nbl2_fphi, psi_dot_nbl2_fpsi

    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_force)

    allocate(phi_dot_nbl2_fphi(nkz, nky_local, nkx), source=0.d0)
    allocate(psi_dot_nbl2_fpsi(nkz, nky_local, nkx), source=0.d0)
    
    if (fix_power) then
      phi_dot_nbl2_fphi_sum = 0.d0
      psi_dot_nbl2_fpsi_sum = 0.d0

      !$acc data present(phi, psi, fphi, fpsi) create(phi_dot_nbl2_fphi, psi_dot_nbl2_fpsi)
      !$acc parallel loop collapse(3) reduction(+:phi_dot_nbl2_fphi_sum, psi_dot_nbl2_fpsi_sum)
      do i =1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            phi_dot_nbl2_fphi(k, j, i) = - 0.5d0*( &
                                    phi(k, j, i)*conjg(kprp2(k, j, i)*fphi(k, j, i)) &
                                  + conjg(phi(k, j, i))*kprp2(k, j, i)*fphi(k, j, i) &
                                ) 
            psi_dot_nbl2_fpsi(k, j, i) = - 0.5d0*( &
                                    psi(k, j, i)*conjg(kprp2(k, j, i)*fpsi(k, j, i)) &
                                  + conjg(psi(k, j, i))*kprp2(k, j, i)*fpsi(k, j, i) &
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
              phi_dot_nbl2_fphi(k, j, i) = 2.0d0*phi_dot_nbl2_fphi(k, j, i)
              psi_dot_nbl2_fpsi(k, j, i) = 2.0d0*psi_dot_nbl2_fpsi(k, j, i)
            endif

            phi_dot_nbl2_fphi_sum = phi_dot_nbl2_fphi_sum + phi_dot_nbl2_fphi(k, j, i)
            psi_dot_nbl2_fpsi_sum = psi_dot_nbl2_fpsi_sum + psi_dot_nbl2_fpsi(k, j, i)

          end do
        end do
      end do
      !$acc end data

      call sum_allreduce(phi_dot_nbl2_fphi_sum)
      call sum_allreduce(psi_dot_nbl2_fpsi_sum)

      if(abs(phi_dot_nbl2_fphi_sum) > 1.0d-6 .and. abs(psi_dot_nbl2_fpsi_sum) > 1.0d-6) then
        !$acc parallel loop present(fphi, fpsi) copyin(phi_dot_nbl2_fphi_sum, psi_dot_nbl2_fpsi_sum)
        do i = 1, nkx
           do j = 1, nky_local
              do k = 1, nkz
                 fphi(k, j, i) = 0.5d0*ene_inj*(1.d0 + res_inj)*fphi(k, j, i)/abs(phi_dot_nbl2_fphi_sum)*(-sign(1.d0, phi_dot_nbl2_fphi_sum))
                 fpsi(k, j, i) = 0.5d0*ene_inj*(1.d0 - res_inj)*fpsi(k, j, i)/abs(psi_dot_nbl2_fpsi_sum)*(-sign(1.d0, psi_dot_nbl2_fpsi_sum))
              end do
           end do
        end do
      endif
    endif

    deallocate(phi_dot_nbl2_fphi, psi_dot_nbl2_fpsi)

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force


!-----------------------------------------------!
!> @author  YK
!! @brief   Compute the power and helicity of forcing 
!!          and normalize with it 
!-----------------------------------------------!
  subroutine normalize_force_els(fzppe, fzmpe)
    use grid, only: nkx, nky_local, nkz, kprp2
    use fields, only: phi, psi
    use mp, only: sum_allreduce
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: fzppe, fzmpe
    complex(8) :: zppe, zmpe
    real(8) :: zppe_dot_nbl2_fzppe_sum, zmpe_dot_nbl2_fzmpe_sum
    real(8), dimension(:,:,:), allocatable :: zppe_dot_nbl2_fzppe, zmpe_dot_nbl2_fzmpe

    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_force)

    allocate(zppe_dot_nbl2_fzppe(nkz, nky_local, nkx), source=0.d0)
    allocate(zmpe_dot_nbl2_fzmpe(nkz, nky_local, nkx), source=0.d0)
    
    if (fix_power) then
      zppe_dot_nbl2_fzppe_sum = 0.d0
      zmpe_dot_nbl2_fzmpe_sum = 0.d0

      !$acc data present(phi, psi, fphi, fpsi) create(zppe_dot_nbl2_fzppe, zmpe_dot_nbl2_fzmpe)
      !$acc parallel loop collapse(3) reduction(+:zppe_dot_nbl2_fzppe_sum, zmpe_dot_nbl2_fzmpe_sum)
      do i =1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            zppe = phi(k, j, i) + psi(k, j, i)
            zmpe = phi(k, j, i) - psi(k, j, i)
            zppe_dot_nbl2_fzppe(k, j, i) = - 0.5d0*( &
                                    zppe*conjg(kprp2(k, j, i)*fzppe(k, j, i)) &
                                  + conjg(zppe)*kprp2(k, j, i)*fzppe(k, j, i) &
                                ) 
            zmpe_dot_nbl2_fzmpe(k, j, i) = - 0.5d0*( &
                                    zmpe*conjg(kprp2(k, j, i)*fzmpe(k, j, i)) &
                                  + conjg(zmpe)*kprp2(k, j, i)*fzmpe(k, j, i) &
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
              zppe_dot_nbl2_fzppe(k, j, i) = 2.0d0*zppe_dot_nbl2_fzppe(k, j, i)
              zmpe_dot_nbl2_fzmpe(k, j, i) = 2.0d0*zmpe_dot_nbl2_fzmpe(k, j, i)
            endif

            zppe_dot_nbl2_fzppe_sum = zppe_dot_nbl2_fzppe_sum + zppe_dot_nbl2_fzppe(k, j, i)
            zmpe_dot_nbl2_fzmpe_sum = zmpe_dot_nbl2_fzmpe_sum + zmpe_dot_nbl2_fzmpe(k, j, i)

          end do
        end do
      end do
      !$acc end data

      call sum_allreduce(zppe_dot_nbl2_fzppe_sum)
      call sum_allreduce(zmpe_dot_nbl2_fzmpe_sum)

      if(abs(zppe_dot_nbl2_fzppe_sum) > 1.0d-6 .and. abs(zmpe_dot_nbl2_fzmpe_sum) > 1.0d-6) then
        !$acc parallel loop present(fphi, fpsi) copyin(zppe_dot_nbl2_fzppe_sum, zmpe_dot_nbl2_fzmpe_sum)
        do i = 1, nkx
           do j = 1, nky_local
              do k = 1, nkz
                 fzppe(k, j, i) = 0.5d0*ene_inj*(1.d0 + xhl_inj)*fzppe(k, j, i)/abs(zppe_dot_nbl2_fzppe_sum)*(-sign(1.d0, zppe_dot_nbl2_fzppe_sum))
                 fzmpe(k, j, i) = 0.5d0*ene_inj*(1.d0 - xhl_inj)*fzmpe(k, j, i)/abs(zmpe_dot_nbl2_fzmpe_sum)*(-sign(1.d0, zmpe_dot_nbl2_fzmpe_sum))
              end do
           end do
        end do
      endif
    endif

    deallocate(zppe_dot_nbl2_fzppe, zmpe_dot_nbl2_fzmpe)

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force_els

end module force

