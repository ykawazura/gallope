!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../force_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2021
!! @brief   Forcing setting specific to MHD_INCOMP
!!          'force_common' is inherited
!-----------------------------------------------!
module force
  use force_common
  implicit none

  public get_force

  complex(8), dimension(:,:,:), allocatable :: fux, fux_old
  complex(8), dimension(:,:,:), allocatable :: fuy, fuy_old
  complex(8), dimension(:,:,:), allocatable :: fuz, fuz_old
  complex(8), dimension(:,:,:), allocatable :: fbx, fbx_old
  complex(8), dimension(:,:,:), allocatable :: fby, fby_old
  complex(8), dimension(:,:,:), allocatable :: fbz, fbz_old
  complex(8), dimension(:,:,:), allocatable :: fzpx
  complex(8), dimension(:,:,:), allocatable :: fzpy
  complex(8), dimension(:,:,:), allocatable :: fzpz
  complex(8), dimension(:,:,:), allocatable :: fzmx
  complex(8), dimension(:,:,:), allocatable :: fzmy
  complex(8), dimension(:,:,:), allocatable :: fzmz

  !$acc declare create(fux, fuy, fuz)
  !$acc declare create(fbx, fby, fbz)
  !$acc declare create(fux_old, fuy_old, fuz_old)
  !$acc declare create(fbx_old, fby_old, fbz_old)
  !$acc declare create(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)

contains


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2021
!! @brief   Initialize model specific subroutines
!-----------------------------------------------!
  subroutine alloc_force
    use grid, only: nkx, nky_local, nkz
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    allocate(src(nkz, nky_local, nkx), source=(0.d0,0.d0))
    allocate(fux    , source=src)
    allocate(fux_old, source=src)
    allocate(fuy    , source=src)
    allocate(fuy_old, source=src)
    allocate(fuz    , source=src)
    allocate(fuz_old, source=src)
    allocate(fbx    , source=src)
    allocate(fbx_old, source=src)
    allocate(fby    , source=src)
    allocate(fby_old, source=src)
    allocate(fbz    , source=src)
    allocate(fbz_old, source=src)
    if(elsasser) then
      allocate(fzpx, source=src)
      allocate(fzpy, source=src)
      allocate(fzpz, source=src)
      allocate(fzmx, source=src)
      allocate(fzmy, source=src)
      allocate(fzmz, source=src)
    endif
    !$acc update device(fux, fuy, fuz)
    !$acc update device(fbx, fby, fbz)
    !$acc update device(fux_old, fuy_old, fuz_old)
    !$acc update device(fbx_old, fby_old, fbz_old)
    !$acc update device(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
    deallocate(src)

  end subroutine alloc_force


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2021
!! @brief   Compute the power of forcing 
!!          and normalize with it 
!-----------------------------------------------!
  subroutine normalize_force(fux, fuy, fuz, fbx, fby, fbz)
    use grid, only: nkx, nky_local, nkz
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use mp, only: sum_allreduce
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: fux, fuy, fuz, fbx, fby, fbz
    real(8) :: p_ext_ene_sum
    real(8), dimension(:,:,:), allocatable :: p_ext_ene

    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_force)

    allocate(p_ext_ene(nkz, nky_local, nkx), source=0.d0)
    
    if (fix_power) then
      p_ext_ene_sum = 0.d0

      !$acc data present(ux, uy, uz, bx, by, bz, fux, fuy, fuz, fbx, fby, fbz) create(p_ext_ene)
      !$acc parallel loop collapse(3) reduction(+:p_ext_ene_sum)
      do i =1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            p_ext_ene(k, j, i) = 0.5d0*( &
                                    (fux(k, j, i)*conjg(ux(k, j, i)) + conjg(fux(k, j, i))*ux(k, j, i)) &
                                  + (fuy(k, j, i)*conjg(uy(k, j, i)) + conjg(fuy(k, j, i))*uy(k, j, i)) &
                                  + (fuz(k, j, i)*conjg(uz(k, j, i)) + conjg(fuz(k, j, i))*uz(k, j, i)) &
                                  + (fbx(k, j, i)*conjg(bx(k, j, i)) + conjg(fbx(k, j, i))*bx(k, j, i)) &
                                  + (fby(k, j, i)*conjg(by(k, j, i)) + conjg(fby(k, j, i))*by(k, j, i)) &
                                  + (fbz(k, j, i)*conjg(bz(k, j, i)) + conjg(fbz(k, j, i))*bz(k, j, i)) &
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
              p_ext_ene(k, j, i) = 2.0d0*p_ext_ene(k, j, i)
            endif

            p_ext_ene_sum = p_ext_ene_sum + p_ext_ene(k, j, i)

          end do
        end do
      end do
      !$acc end data

      call sum_allreduce(p_ext_ene_sum)

      if(p_ext_ene_sum > 1.0d-6) then
        !$acc parallel loop present(fux, fuy, fuz, fbx, fby, fbz) copyin(p_ext_ene_sum)
        do i = 1, nkx
           do j = 1, nky_local
              do k = 1, nkz
                 fux(k, j, i) = ene_inj*fux(k, j, i)/p_ext_ene_sum
                 fuy(k, j, i) = ene_inj*fuy(k, j, i)/p_ext_ene_sum
                 fuz(k, j, i) = ene_inj*fuz(k, j, i)/p_ext_ene_sum
                 fbx(k, j, i) = ene_inj*fbx(k, j, i)/p_ext_ene_sum
                 fby(k, j, i) = ene_inj*fby(k, j, i)/p_ext_ene_sum
                 fbz(k, j, i) = ene_inj*fbz(k, j, i)/p_ext_ene_sum
              end do
           end do
        end do
      endif
    endif

    deallocate(p_ext_ene)

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2021
!! @brief   Compute the power and helicity of forcing 
!!          and normalize with it 
!-----------------------------------------------!
  subroutine normalize_force_els(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
    use grid, only: nkx, nky_local, nkz
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use mp, only: sum_allreduce
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(8), dimension (:,:,:), intent(inout) :: fzpx, fzpy, fzpz, fzmx, fzmy, fzmz
    complex(8) :: zpx, zpy, zpz, zmx, zmy, zmz
    real(8)    :: zp_dot_fzp_sum, zm_dot_fzm_sum
    real(8), allocatable, dimension(:,:,:) :: zp_dot_fzp, zm_dot_fzm

    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_force)

    allocate(zp_dot_fzp(nkz, nky_local, nkx), source=0.d0)
    allocate(zm_dot_fzm(nkz, nky_local, nkx), source=0.d0)
    
    if (fix_power) then
      zp_dot_fzp_sum = 0.d0
      zm_dot_fzm_sum = 0.d0

      !$acc data present(ux, uy, uz, bx, by, bz, fzpx, fzpy, fzpz, fzmx, fzmy, fzmz) create(zp_dot_fzp, zm_dot_fzm)
      !$acc parallel loop collapse(3) reduction(+:zp_dot_fzp_sum, zm_dot_fzm_sum)
      do i =1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            zpx = ux(k, j, i) + bx(k, j, i)
            zpy = uy(k, j, i) + by(k, j, i)
            zpz = uz(k, j, i) + bz(k, j, i)
            zmx = ux(k, j, i) - bx(k, j, i)
            zmy = uy(k, j, i) - by(k, j, i)
            zmz = uz(k, j, i) - bz(k, j, i)

            zp_dot_fzp(k, j, i) = 0.5d0*( &
                                    (fzpx(k, j, i)*conjg(zpx) + conjg(fzpx(k, j, i))*zpx) &
                                  + (fzpy(k, j, i)*conjg(zpy) + conjg(fzpy(k, j, i))*zpy) &
                                  + (fzpz(k, j, i)*conjg(zpz) + conjg(fzpz(k, j, i))*zpz) &
                                )
            zm_dot_fzm(k, j, i) = 0.5d0*( &
                                    (fzmx(k, j, i)*conjg(zmx) + conjg(fzmx(k, j, i))*zmx) &
                                  + (fzmy(k, j, i)*conjg(zmy) + conjg(fzmy(k, j, i))*zmy) &
                                  + (fzmz(k, j, i)*conjg(zmz) + conjg(fzmz(k, j, i))*zmz) &
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
              zp_dot_fzp(k, j, i) = 2.0d0*zp_dot_fzp(k, j, i)
              zm_dot_fzm(k, j, i) = 2.0d0*zm_dot_fzm(k, j, i)
            endif

            zp_dot_fzp_sum = zp_dot_fzp_sum + zp_dot_fzp(k, j, i)
            zm_dot_fzm_sum = zm_dot_fzm_sum + zm_dot_fzm(k, j, i)

          end do
        end do
      end do
      !$acc end data

      call sum_allreduce(zp_dot_fzp_sum)
      call sum_allreduce(zm_dot_fzm_sum)

      if(abs(zp_dot_fzp_sum) > 1.0d-6 .and. abs(zm_dot_fzm_sum) > 1.0d-6) then
        !$acc parallel loop present(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz) copyin(zp_dot_fzp_sum, zm_dot_fzm_sum)
        do i = 1, nkx
           do j = 1, nky_local
              do k = 1, nkz
                fzpx(k, j, i) = 0.5d0*ene_inj*(1.d0 + xhl_inj)*fzpx(k, j, i)/abs(zp_dot_fzp_sum)
                fzpy(k, j, i) = 0.5d0*ene_inj*(1.d0 + xhl_inj)*fzpy(k, j, i)/abs(zp_dot_fzp_sum)
                fzpz(k, j, i) = 0.5d0*ene_inj*(1.d0 + xhl_inj)*fzpz(k, j, i)/abs(zp_dot_fzp_sum)
                fzmx(k, j, i) = 0.5d0*ene_inj*(1.d0 - xhl_inj)*fzmx(k, j, i)/abs(zm_dot_fzm_sum)
                fzmy(k, j, i) = 0.5d0*ene_inj*(1.d0 - xhl_inj)*fzmy(k, j, i)/abs(zm_dot_fzm_sum)
                fzmz(k, j, i) = 0.5d0*ene_inj*(1.d0 - xhl_inj)*fzmz(k, j, i)/abs(zm_dot_fzm_sum)
              end do
           end do
        end do
      endif
    endif

    deallocate(zp_dot_fzp)
    deallocate(zm_dot_fzm)

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force_els

end module force

