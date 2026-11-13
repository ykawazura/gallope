!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../params_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Parameter setting specific to RRMHD (rotating reduced
!!          MHD with parallel flow and compressive field)
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
  public  beta, gamma, kappa_b, kappa_p
  public  va2cs2_plus_1, cs2va2
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
  real(8) :: beta     !< plasma beta
  real(8) :: gamma    !< adiabatic index
  real(8) :: q        !< shear parameter of the rotating frame (3/2: Keplerian)
  real(8) :: kappa_b  !< background magnetic-field gradient
  real(8) :: kappa_p  !< background pressure gradient
  real(8) :: va2cs2_plus_1  !< A = 1 + vA^2/cs^2 = 1 + 2/(beta*gamma)
  real(8) :: cs2va2         !< c = cs^2/vA^2 = beta*gamma/2
  ! The shearing-box remap is not used: rotation enters only through the
  ! linear q and Coriolis terms of the equations (as in Calliope RRMHD).
  logical, parameter  :: shear = .false.
  !$acc declare create(nupe_x , nupe_x_exp , nupe_z , nupe_z_exp)
  !$acc declare create(etape_x, etape_x_exp, etape_z, etape_z_exp)
  !$acc declare create(nupa_x , nupa_x_exp , nupa_z , nupa_z_exp)
  !$acc declare create(etapa_x, etapa_x_exp, etapa_z, etapa_z_exp)
  !$acc declare create(beta, gamma, q, kappa_b, kappa_p)
  !$acc declare create(va2cs2_plus_1, cs2va2)
  !$acc declare create(shear)

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
                                   beta, gamma, q, kappa_b, kappa_p

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
    ! Calliope leaves the eight exponents undefined; they get defaults here
    nupe_x_exp  = 1; nupe_z_exp  = 1
    etape_x_exp = 1; etape_z_exp = 1
    nupa_x_exp  = 1; nupa_z_exp  = 1
    etapa_x_exp = 1; etapa_z_exp = 1
    beta    = 1.d0
    gamma   = 5.d0/3.d0
    q       = 0.d0
    kappa_b = 0.d0
    kappa_p = 0.d0
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
    ! Derived quantities; they must reach the device after they are set
    ! (A = 0 on the device would divide every bpa term by zero).
    va2cs2_plus_1 = 2.d0/(beta*gamma) + 1.d0
    cs2va2        = (beta*gamma)/2.d0

    !$acc update device(beta, gamma, q, kappa_b, kappa_p)
    !$acc update device(va2cs2_plus_1, cs2va2)

  end subroutine read_parameters

end module params
