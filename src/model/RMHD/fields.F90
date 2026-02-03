!-----------------------------------------------!
!> @author  YK
!! @brief   Field setting for RMHD
!-----------------------------------------------!
module fields
  implicit none

  public :: init_fields, finish_fields
  public :: phi, omg, psi
  public :: phi_old, omg_old, psi_old
  public :: nfields
  public :: iomg, ipsi

  private

  complex(8), allocatable, dimension(:,:,:) :: phi, omg, psi
  complex(8), allocatable, dimension(:,:,:) :: phi_old, omg_old, psi_old
  character(100) :: init_type

  ! Field index
  integer, parameter :: nfields = 2
  integer, parameter :: iomg = 1, ipsi = 2

  !$acc declare create(phi, omg, psi)
  !$acc declare create(phi_old, omg_old, psi_old)
  !$acc declare copyin(nfields)
  !$acc declare copyin(iomg, ipsi)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of fields
!-----------------------------------------------!
  subroutine init_fields
    use grid, only: nkx, nky_local, nkz
    use params, only: inputfile
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    allocate(src(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(phi    , source=src)
    allocate(omg    , source=src)
    allocate(psi    , source=src)
    allocate(phi_old, source=src)
    allocate(omg_old, source=src)
    allocate(psi_old, source=src)
    !$acc update device(phi, omg, psi)
    !$acc update device(phi_old, omg_old, psi_old)
    deallocate(src)

    call read_parameters(inputfile)

    if(init_type == 'zero') then
    endif
    if(init_type == 'single_mode') then
      call init_single_mode
    endif
    if(init_type == 'OT2') then
      call init_OT2
    endif
    ! if(init_type == 'OT3') then
    !   call init_OT3
    ! endif
    ! if(init_type == 'KH') then
    !   call init_KH
    ! endif
    if(init_type == 'random') then
      call init_random
    endif
    if(init_type == 'restart') then
      call restart
    endif
    ! if(init_type == 'convert_restartfiles_to_real') then
    !   call convert_restartfiles_to_real
    ! endif

  end subroutine init_fields


!-----------------------------------------------!
!> @author  YK
!! @brief   Read inputfile for initial condition
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /initial_condition/ init_type

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

    init_type = 'zero'
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=initial_condition,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @brief   Single mode initialization
!-----------------------------------------------!
  subroutine init_single_mode
    use grid, only: nkx, nky_local, nkz
    use grid, only: kprp2
    use mp, only: proc0, proc_id
    implicit none
    integer :: i, j, k

    if(proc0) then
      print '("Single mode initialization")'
    endif

    !$acc data present(phi, psi, kprp2)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if (i == 1 .and. j == 1 .and. k == 2) then
            phi(k, j, i) = 1.d0
            psi(k, j, i) = 1.d0
          endif
          omg(k, j, i) = phi(k, j, i)*(-kprp2(k, j, i))
        enddo
      enddo
    enddo
    !$acc end data
    
    !$acc data present(phi, omg, psi, phi_old, omg_old, psi_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi_old(k, j, i) = phi(k, j, i)
          omg_old(k, j, i) = omg(k, j, i)
          psi_old(k, j, i) = psi(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine init_single_mode
!
!
!-----------------------------------------------!
!> @author  YK
!! @brief   Random initialization
!-----------------------------------------------!
  subroutine init_random
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky, nky_local, nkz
    use grid, only: kx, ky, kz, kprp2
    use grid, only: ntot
    use mp, only: proc0, proc_id
    use time, only: microsleep
    use cuFFTmp, only: ftran_r2c, btran_c2r
    implicit none
    real(8), allocatable, dimension(:,:,:) :: phi_r, psi_r
    integer :: seedsize
    integer, allocatable :: seed(:)
    integer :: i, j, k

    if(proc0) then
      print *, 'Random initialization'
    endif

    allocate(phi_r(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(psi_r(nlz_padded, nly, nlx_local), source=0.d0)
    !$acc enter data create(phi_r)
    !$acc enter data create(psi_r)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             create random number            v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    call random_seed(size=seedsize)
    allocate(seed(seedsize)) 

    do i = 1, seedsize
      call system_clock(count=seed(i))
      call microsleep(1000)
      call system_clock(count=seed(i))
    end do
    call random_seed(put=(proc_id+1)*seed(:)) 

    call random_number(phi_r)
    call random_number(psi_r)
    !$acc update device(phi_r)
    !$acc update device(psi_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             compute r2c transform           v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !$acc data present(psi_r, phi_r, phi, psi)
    call ftran_r2c(phi_r, phi)
    call ftran_r2c(psi_r, psi)
    !$acc end data
    
    !$acc data present(phi, psi)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi(k,j,i)/ntot
          psi(k,j,i) = psi(k,j,i)/ntot
        enddo
      enddo
    enddo
    !$acc end data

    !$acc data present(phi, psi, kx, ky, kz)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if (kx(i)**2 + ky(j)**2 + kz(k)**2 /= 0.d0) then
            phi(k, j, i) = phi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            psi(k, j, i) = psi(k, j, i)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
          endif
          if (kx(i)**2 + ky(j)**2 + kz(k)**2 >= (0.8*maxval(kx))**2) then
            phi(k, j, i) = 0.d0
            psi(k, j, i) = 0.d0
          endif
        enddo
      enddo
    enddo
    !$acc end data

    !$acc data present(phi, omg, kprp2)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          omg(k, j, i) = phi(k, j, i)*(-kprp2(k, j, i))
        enddo
      enddo
    enddo
    !$acc end data
    
    !$acc data present(phi, omg, psi, phi_old, omg_old, psi_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi_old(k, j, i) = phi(k, j, i)
          omg_old(k, j, i) = omg(k, j, i)
          psi_old(k, j, i) = psi(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data
    
    !$acc exit data delete(phi_r)
    !$acc exit data delete(psi_r)
    deallocate(phi_r)
    deallocate(psi_r)

  end subroutine init_random


!-----------------------------------------------!
!> @author  YK
!! @brief   2D Orszag Tang problem initialization
!-----------------------------------------------!
  subroutine init_OT2
    use grid, only: xx, yy
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0
    use params, only: zi, pi
    use cuFFTmp, only: ftran_r2c
    implicit none
    real(8), allocatable, dimension(:,:,:) :: phi_r, psi_r
    real(8), allocatable, dimension(:,:,:) :: src
    integer :: i, j, k

    if(proc0) then
      print *, 'OT2 initialization'
    endif

    allocate(src(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(phi_r, source=src)
    allocate(psi_r, source=src)
    deallocate(src)
    !$acc enter data create(phi_r)
    !$acc enter data create(psi_r)

    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          phi_r(k, j, i) = -cos(xx(i)) - cos(yy(j))
          psi_r(k, j, i) = -0.5d0*cos(2.d0*(xx(i))) - cos(yy(j))
        end do
      end do
    end do
  
    ! compute r2c transform
    call ftran_r2c(phi_r, phi)
    call ftran_r2c(psi_r, psi)
  
    !$acc data present(phi, psi)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi(k,j,i)/ntot
          psi(k,j,i) = psi(k,j,i)/ntot
        enddo
      enddo
    enddo
    !$acc end data

    !$acc data present(phi, omg, psi, phi_old, omg_old, psi_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi_old(k, j, i) = phi(k, j, i)
          omg_old(k, j, i) = omg(k, j, i)
          psi_old(k, j, i) = psi(k, j, i)
        enddo
      enddo
    enddo
    !$acc end data

    !$acc exit data delete(phi_r)
    !$acc exit data delete(psi_r)
    deallocate(phi_r)
    deallocate(psi_r)

  end subroutine init_OT2


!-----------------------------------------------!
!> @author  YK
!! @brief   Restart
!-----------------------------------------------!
  subroutine restart
    use mp, only: proc0, proc_id
    use time, only: tt, dt
    use grid, only: nkx, nky, nky_local, nkz
    use params, only: restart_dir
    use file, only: open_input_file, close_file
    use mpiio, only: mpiio_read_one
    use shearing_box, only: tsc
    use MPI
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts

    if(proc0) then
      print *, 'Restart'
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                   Read time                 v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    call open_input_file (time_unit, trim(restart_dir)//'time.dat')
    read (unit=time_unit, fmt=*)
    read (unit=time_unit, fmt="(100es30.21)") tt, tsc
    call close_file (time_unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    !$acc update device(tt, tsc)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v               Read Binary file              v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    sizes(1) = nkz
    sizes(2) = nky
    sizes(3) = nkx
    subsizes(1) = nkz
    subsizes(2) = nky_local
    subsizes(3) = nkx
    starts(1) = 0
    starts(2) = nky_local*proc_id
    starts(3) = 0

    call mpiio_read_one(phi, sizes, subsizes, starts, trim(restart_dir)//'phi.dat')
    call mpiio_read_one(omg, sizes, subsizes, starts, trim(restart_dir)//'omg.dat')
    call mpiio_read_one(psi, sizes, subsizes, starts, trim(restart_dir)//'psi.dat')

    phi_old = phi
    omg_old = omg
    psi_old = psi

    !$acc update device(phi, omg, psi)
    !$acc update device(phi_old, omg_old, psi_old)
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @brief   Finalization of Fields
!-----------------------------------------------!
  subroutine finish_fields
    implicit none

    deallocate(phi)
    deallocate(omg)
    deallocate(psi)
    deallocate(phi_old)
    deallocate(omg_old)
    deallocate(psi_old)

  end subroutine finish_fields

end module fields

