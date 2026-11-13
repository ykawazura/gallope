!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../force_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Forcing setting specific to RRMHD
!!          'force_common' is inherited
!!
!!          RRMHD has no forcing. Calliope reads &force,
!!          allocates unused forcing arrays and silently
!!          ignores driven = .true.; here &force is still
!!          read (so common inputs work) and driven = .true.
!!          stops the run (check_no_forcing).
!-----------------------------------------------!
module force
  use force_common
  implicit none

  public check_no_forcing

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Stop when forcing is requested
!-----------------------------------------------!
  subroutine check_no_forcing
    use MPI
    use mp, only: proc0
    implicit none
    integer :: ierr

    if (driven) then
      if (proc0) print *, 'RRMHD has no forcing: set driven = .false. in &force'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif

  end subroutine check_no_forcing

end module force
