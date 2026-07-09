!-----------------------------------------------!
!> @author  YK
!! @brief   Model specific subroutines
!-----------------------------------------------!
module model_specific
  use force_common
  implicit none

  public init_model_specific

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
end module model_specific

