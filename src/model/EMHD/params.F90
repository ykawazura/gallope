!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../params_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Parameter setting specific to EMHD
!!          (electron MHD with electron inertia)
!!          'params_common' is inherited
!-----------------------------------------------!
module params
  use params_common
  implicit none

  public  init_params
  public  nonlinear
  public  eta, eta_exp, de2
  public  shear, q
  private read_parameters

  logical :: nonlinear
  real(8) :: eta      !< coefficient of the (hyper-)resistivity eta*(k^2/k2_max)^eta_exp
  integer :: eta_exp  !< exponent of the (hyper-)resistivity
  real(8) :: de       !< electron inertial length in units of the box
  real(8) :: de2      !< de**2
  logical, parameter  :: shear = .false.
  real(8), parameter :: q = 0.d0
  !$acc declare create(eta, eta_exp, de2)
  !$acc declare create(shear, q)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of run parameters,
!!          followed by input file reading
!-----------------------------------------------!
  subroutine init_params
    implicit none

    call init_params_common
    call read_parameters(inputfile)

  end subroutine init_params


!-----------------------------------------------!
!> @author  YK
!! @brief   Read inputfile for various parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none

    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /operation_parameters/ nonlinear
    namelist /physical_parameters/ eta, eta_exp, de

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    nonlinear = .false.
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=operation_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading operation_parameters failed"
    close(unit)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    eta     = 0.d0
    eta_exp = 1
    de      = 0.d0
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=physical_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading physical_parameters failed"
    close(unit)

    de2 = de**2

    !$acc update device(eta, eta_exp, de2)

  end subroutine read_parameters

end module params
