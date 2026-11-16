!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../force_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2021
!! @brief   Forcing setting specific to MHD_COMP_ISOTH
!!          'force_common' is inherited
!!
!!          driven = T aborts in init_model_specific,
!!          so fmx, fmy, fmz stay zero; they exist
!!          because advance and diagnostics keep the
!!          forcing terms of Calliope.
!-----------------------------------------------!
module force
  use force_common
  implicit none

  public get_force
  public alloc_force
  public fmx, fmy, fmz, fmx_old, fmy_old, fmz_old

  complex(8), allocatable, dimension(:,:,:) :: fmx, fmx_old
  complex(8), allocatable, dimension(:,:,:) :: fmy, fmy_old
  complex(8), allocatable, dimension(:,:,:) :: fmz, fmz_old

  !$acc declare create(fmx, fmx_old)
  !$acc declare create(fmy, fmy_old)
  !$acc declare create(fmz, fmz_old)

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
    allocate(fmx    , source=src)
    allocate(fmx_old, source=src)
    allocate(fmy    , source=src)
    allocate(fmy_old, source=src)
    allocate(fmz    , source=src)
    allocate(fmz_old, source=src)
    !$acc update device(fmx, fmx_old)
    !$acc update device(fmy, fmy_old)
    !$acc update device(fmz, fmz_old)
    deallocate(src)
  end subroutine alloc_force

end module force
