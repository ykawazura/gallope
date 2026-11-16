!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Model specific subroutines
!-----------------------------------------------!
module model_specific
  use force_common
  implicit none

  public init_model_specific

  private abort_specific

contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialize model specific subroutines
!!
!!          Gallope runs MHD_COMP_ISOTH only with
!!          shear = F, driven = F, formulation = 'C',
!!          visc_form = 'flux', time_step_scheme =
!!          'eSSPIFRK3' and P_m = P_s = 1, and needs
!!          nlx and nly divisible by P_fft (x slabs in
!!          real space, ky slabs in Fourier space).
!!          Any other setting stops the run here.
!!
!!          Calliope's report of k0*vA/Omega,
!!          L/lmd_MRI and L/H for shear = T is not
!!          ported, because shear = T stops the run.
!-----------------------------------------------!
  subroutine init_model_specific
    use mp, only: nproc_fft
    use grid, only: nlx, nly
    use force, only :init_force, alloc_force
    use shearing_box, only: init_shearing_time
    use params, only: shear, formulation, visc_form
    use params, only: time_step_scheme, P_m, P_s
    implicit none
    character(len=200) :: msg

    if(shear) &
      call abort_specific('MHD_COMP_ISOTH: shear = T is not supported in Gallope')
    if(trim(formulation) /= 'C') &
      call abort_specific('MHD_COMP_ISOTH: formulation = '//trim(formulation)//' is not implemented (use C)')
    if(trim(visc_form) /= 'flux') &
      call abort_specific('MHD_COMP_ISOTH: visc_form = '//trim(visc_form)//' is not implemented yet (use flux)')
    if(trim(time_step_scheme) /= 'eSSPIFRK3') &
      call abort_specific('MHD_COMP_ISOTH: time_step_scheme = '//trim(time_step_scheme)// &
                          ' is not implemented (use eSSPIFRK3)')
    if(P_m /= 1 .or. P_s /= 1) &
      call abort_specific('MHD_COMP_ISOTH supports P_m = P_s = 1 only')
    if(mod(nlx, nproc_fft) /= 0 .or. mod(nly, nproc_fft) /= 0) then
      write(msg, '("MHD_COMP_ISOTH: nlx and nly must be divisible by P_fft (nlx = ", i0, &
                 & ", nly = ", i0, ", P_fft = ", i0, ")")') nlx, nly, nproc_fft
      call abort_specific(trim(msg))
    endif

    call init_shearing_time

    ! if(proc0 .and. maxval(kx**2) /= maxval(ky**2)) then
      ! print *, maxval(kx**2) , maxval(ky**2)
      ! print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      ! print *, '!              Error!              !'
      ! print *, '!  kx_max must be equal to ky_max  !'
      ! print *, '!  i.e., lx/nlx == ly/nly          !'
      ! print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      ! stop
    ! endif

    call init_force
    call alloc_force

    ! driven is read by init_force
    if(driven) &
      call abort_specific('MHD_COMP_ISOTH: driven = T is not supported in Gallope (fmx, fmy, fmz are zero)')

  end subroutine init_model_specific


!-----------------------------------------------!
!> @author  YK
!! @brief   Print msg on proc0 and abort all ranks.
!!          Every rank evaluates the same global
!!          condition, so the barrier lets proc0
!!          print before any rank aborts.
!-----------------------------------------------!
  subroutine abort_specific(msg)
    use mp, only: proc0
    use MPI
    implicit none
    character(len=*), intent(in) :: msg
    integer :: ierr

    if (proc0) then
      print *, trim(msg)
      flush(6)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  end subroutine abort_specific

end module model_specific
