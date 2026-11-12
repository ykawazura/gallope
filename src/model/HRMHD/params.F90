!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../params_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Parameter setting specific to HRMHD (Hall reduced MHD)
!!          'params_common' is inherited
!-----------------------------------------------!
module params
  use params_common
  implicit none

  public  init_params
  public  nonlinear
  public  nupe_x , nupe_x_exp , nupe_z , nupe_z_exp
  public  etape_x, etape_x_exp, etape_z, etape_z_exp
  public  nupa_x , nupa_x_exp , nupa_z , nupa_z_exp
  public  etapa_x, etapa_x_exp, etapa_z, etapa_z_exp
  public  rho, sgm
  public  shear, q
  private read_parameters

  logical :: nonlinear
  real(8) :: nupe_x , nupe_z
  real(8) :: etape_x, etape_z
  integer :: nupe_x_exp , nupe_z_exp
  integer :: etape_x_exp, etape_z_exp
  real(8) :: nupa_x , nupa_z
  real(8) :: etapa_x, etapa_z
  integer :: nupa_x_exp , nupa_z_exp
  integer :: etapa_x_exp, etapa_z_exp
  real(8) :: rho  !< Hall length (ion sound/inertial scale) in units of the box
  real(8) :: sgm  !< parallel phase-speed ratio of the compressive (slow) branch
  logical, parameter  :: shear = .false.
  real(8), parameter :: q = 0.d0
  !$acc declare create(nupe_x , nupe_x_exp , nupe_z , nupe_z_exp)
  !$acc declare create(etape_x, etape_x_exp, etape_z, etape_z_exp)
  !$acc declare create(nupa_x , nupa_x_exp , nupa_z , nupa_z_exp)
  !$acc declare create(etapa_x, etapa_x_exp, etapa_z, etapa_z_exp)
  !$acc declare create(rho, sgm)
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
    namelist /physical_parameters/ nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                                   etape_x, etape_x_exp, etape_z, etape_z_exp, &
                                   nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                                   etapa_x, etapa_x_exp, etapa_z, etapa_z_exp, &
                                   rho, sgm

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
    nupe_x  = 0.d0
    nupe_z  = 0.d0
    etape_x = 0.d0
    etape_z = 0.d0
    nupa_x  = 0.d0
    nupa_z  = 0.d0
    etapa_x = 0.d0
    etapa_z = 0.d0
    nupe_x_exp  = 1; nupe_z_exp  = 1
    etape_x_exp = 1; etape_z_exp = 1
    nupa_x_exp  = 1; nupa_z_exp  = 1
    etapa_x_exp = 1; etapa_z_exp = 1
    rho = 0.d0
    sgm = 1.d0
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=physical_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading physical_parameters failed"
    close(unit)

    !$acc update device(nupe_x , nupe_x_exp , nupe_z , nupe_z_exp)
    !$acc update device(etape_x, etape_x_exp, etape_z, etape_z_exp)
    !$acc update device(nupa_x , nupa_x_exp , nupa_z , nupa_z_exp)
    !$acc update device(etapa_x, etapa_x_exp, etapa_z, etapa_z_exp)
    !$acc update device(rho, sgm)

  end subroutine read_parameters

end module params
