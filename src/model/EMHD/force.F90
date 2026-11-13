!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../force_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Forcing setting specific to EMHD
!!          'force_common' is inherited
!!          Forced fields: 'bx', 'by', 'bz'
!-----------------------------------------------!
module force
  use force_common
  implicit none

  public get_force
  public alloc_force, normalize_force
  public fbx, fby, fbz, fbx_old, fby_old, fbz_old

  complex(8), dimension(:,:,:), allocatable :: fbx, fbx_old
  complex(8), dimension(:,:,:), allocatable :: fby, fby_old
  complex(8), dimension(:,:,:), allocatable :: fbz, fbz_old

  !$acc declare create(fbx, fbx_old)
  !$acc declare create(fby, fby_old)
  !$acc declare create(fbz, fbz_old)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Allocate the forcing arrays
!-----------------------------------------------!
  subroutine alloc_force
    use grid, only: nkx, nky_local, nkz
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    allocate(src(nkz, nky_local, nkx), source=(0.d0,0.d0))
    allocate(fbx    , source=src)
    allocate(fbx_old, source=src)
    allocate(fby    , source=src)
    allocate(fby_old, source=src)
    allocate(fbz    , source=src)
    allocate(fbz_old, source=src)
    !$acc update device(fbx, fbx_old)
    !$acc update device(fby, fby_old)
    !$acc update device(fbz, fbz_old)
    deallocate(src)

  end subroutine alloc_force


!-----------------------------------------------!
!> @author  YK
!! @brief   Scale the forcing so that its power
!!          at the current field equals ene_inj
!!
!!          The forcing enters as (1 + de^2 k^2) db/dt = fb,
!!          so with W = 1/2 <(1 + de^2 k^2)|b|^2> it injects
!!            p_ext = sum_k Re(fb . conjg(b)).
!!          Calliope reads fix_power for EMHD but never
!!          applies it; here it takes effect.
!-----------------------------------------------!
  subroutine normalize_force(fbx, fby, fbz)
    use grid, only: nkx, nky_local, nkz
    use fields, only: bx, by, bz
    use mp, only: sum_allreduce
    use mp, only: proc0, comm_fft
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: fbx, fby, fbz
    real(8) :: p_ext, p_ext_sum, fac
    integer :: i, j, k

    if (.not. fix_power) return

    if (proc0) call put_time_stamp(timer_force)

    p_ext_sum = 0.d0
    !$acc parallel loop collapse(3) present(bx, by, bz, fbx, fby, fbz) &
    !$acc& private(p_ext) reduction(+:p_ext_sum)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          p_ext = dble(fbx(k, j, i)*conjg(bx(k, j, i))) &
                + dble(fby(k, j, i)*conjg(by(k, j, i))) &
                + dble(fbz(k, j, i)*conjg(bz(k, j, i)))
          ! kz > 0 modes stand for their kz < 0 partners too
          if (k /= 1) p_ext = 2.d0*p_ext
          p_ext_sum = p_ext_sum + p_ext
        end do
      end do
    end do

    call sum_allreduce(p_ext_sum, comm=comm_fft)

    if (abs(p_ext_sum) > 1.0d-6) then
      fac = ene_inj/p_ext_sum
      !$acc parallel loop collapse(3) present(fbx, fby, fbz) copyin(fac)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            fbx(k, j, i) = fac*fbx(k, j, i)
            fby(k, j, i) = fac*fby(k, j, i)
            fbz(k, j, i) = fac*fbz(k, j, i)
          end do
        end do
      end do
    endif

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force

end module force
