!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../params_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Parameter setting specific to RMHD
!!          'params_common' is inherited
!-----------------------------------------------!
module params
  use params_common
  implicit none

  public  init_params
  public  nonlinear
  public  write_hermite_flux
  public  nupe_x , nupe_x_exp , nupe_z , nupe_z_exp
  public  etape_x, etape_x_exp, etape_z, etape_z_exp
  public  shear, q
  ! KRMHD (Hermite g) parameters
  public  v_th, alpha
  public  mu_hyper_perp, nexp_perp, nu_hyper_m, nexp_m
  public  beta_i, tau, Zcharge, alpha_root
  private read_parameters

  logical :: nonlinear
  logical :: write_hermite_flux   ! gate the (costly) Hermite flux Gamma_m diagnostic (eq 9)
  real(8) :: nupe_x , nupe_z
  real(8) :: etape_x, etape_z
  integer :: nupe_x_exp , nupe_z_exp
  integer :: etape_x_exp, etape_z_exp
  logical, parameter  :: shear = .false.
  real(8), parameter :: q = 0.d0
  !$acc declare create(nupe_x , nupe_x_exp , nupe_z , nupe_z_exp)
  !$acc declare create(etape_x, etape_x_exp, etape_z, etape_z_exp)
  !$acc declare create(shear, q)

  ! KRMHD (Hermite v_parallel moments)
  real(8) :: v_th                          ! thermal speed (streaming coefficient, 1-B)
  real(8) :: alpha                         ! ion-sound coupling; computed from beta_i,tau,Z (1-B)
  real(8) :: mu_hyper_perp, nu_hyper_m     ! perp hyperviscosity / Hermite hypercollision
  integer :: nexp_perp, nexp_m             ! their exponents: mu*(kp2/kp2max)^nexp_perp, nu*(m/nm)^nexp_m
  real(8) :: beta_i, tau, Zcharge          ! ion beta, T_i/T_e, q_i/|e| (alpha inputs)
  integer :: alpha_root                    ! +1/-1: which +/-kappa root of alpha to select
  !$acc declare create(v_th, alpha)
  !$acc declare create(mu_hyper_perp, nu_hyper_m, nexp_perp, nexp_m)

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
    real(8)  :: kappa

    namelist /operation_parameters/ nonlinear, write_hermite_flux
    namelist /physical_parameters/ nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                                   etape_x, etape_x_exp, etape_z, etape_z_exp, &
                                   v_th, mu_hyper_perp, nexp_perp, nu_hyper_m, nexp_m, &
                                   beta_i, tau, Zcharge, alpha_root

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    nonlinear          = .false.
    write_hermite_flux = .false.
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

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
    ! KRMHD defaults: mu=nu=0 -> identity integrating factor (pure passive
    ! advection in 1-A); v_th unused until the streaming term is wired in 1-B.
    v_th          = 1.d0
    mu_hyper_perp = 0.d0
    nexp_perp     = 4        ! mu*k_perp^8 = mu*(k_perp^2)^4
    nu_hyper_m    = 0.d0
    nexp_m        = 6        ! nu*m^6
    beta_i        = 1.d0
    tau           = 1.d0
    Zcharge       = 1.d0
    alpha_root    = 1
    alpha         = 0.d0     ! computed from beta_i,tau,Z in milestone 1-B
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=physical_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading physical_parameters failed"
    close(unit)

    ! ion-sound coupling alpha from (beta_i, tau, Zcharge) [Meyrand et al. 2019].
    !   kappa = sqrt((1 + tau/Z)^2 + 1/beta_i^2)
    !   alpha = 1 / (tau/Z - 1/beta_i + alpha_root*kappa)
    ! The +/-kappa roots are the two compressive modes; alpha_root (+1/-1) picks
    ! one for the single g-hierarchy of phase 1. For physical inputs (beta_i,
    ! tau, Z > 0) the denominator is bounded away from zero for either root.
    if (beta_i > 0.d0) then
      kappa = sqrt((1.d0 + tau/Zcharge)**2 + 1.d0/beta_i**2)
      alpha = 1.d0/(tau/Zcharge - 1.d0/beta_i + dble(alpha_root)*kappa)
    else
      alpha = 0.d0
    endif

    !$acc update device(nupe_x , nupe_x_exp , nupe_z , nupe_z_exp)
    !$acc update device(etape_x, etape_x_exp, etape_z, etape_z_exp)
    !$acc update device(v_th, alpha)
    !$acc update device(mu_hyper_perp, nu_hyper_m, nexp_perp, nexp_m)

  end subroutine read_parameters

end module params
