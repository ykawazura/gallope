!-----------------------------------------------!
!> @author  YK
!! @date    25 Feb 2021
!! @brief   Field setting for MHD_INCOMP
!-----------------------------------------------!
module fields
  implicit none

  public :: init_fields, finish_fields
  public :: ux, uy, uz
  public :: bx, by, bz
  public :: ux_old, uy_old, uz_old
  public :: bx_old, by_old, bz_old
  public :: nfields
  public :: iux, iuy, iuz
  public :: ibx, iby, ibz

  private

  complex(8), dimension(:,:,:), allocatable :: ux, uy, uz
  complex(8), dimension(:,:,:), allocatable :: bx, by, bz
  complex(8), dimension(:,:,:), allocatable :: ux_old, uy_old, uz_old
  complex(8), dimension(:,:,:), allocatable :: bx_old, by_old, bz_old
  character(100) :: init_type
  real   (8) :: b0(3)

  ! Field index
  integer, parameter :: nfields = 6
  integer, parameter :: iux = 1, iuy = 2, iuz = 3
  integer, parameter :: ibx = 4, iby = 5, ibz = 6

  !$acc declare create(ux, uy, uz)
  !$acc declare create(bx, by, bz)
  !$acc declare create(ux_old, uy_old, uz_old)
  !$acc declare create(bx_old, by_old, bz_old)
  !$acc declare copyin(nfields)
  !$acc declare copyin(iux, iuy, iuz)
  !$acc declare copyin(ibx, iby, ibz)

contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialization of fields
!-----------------------------------------------!
  subroutine init_fields
    use grid, only: nkx, nky_local, nkz
    use params, only: inputfile
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    allocate(src(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    allocate(ux    , source=src)
    allocate(uy    , source=src)
    allocate(uz    , source=src)
    allocate(bx    , source=src)
    allocate(by    , source=src)
    allocate(bz    , source=src)
    allocate(ux_old, source=src)
    allocate(uy_old, source=src)
    allocate(uz_old, source=src)
    allocate(bx_old, source=src)
    allocate(by_old, source=src)
    allocate(bz_old, source=src)
    !$acc update device(ux, uy, uz)
    !$acc update device(bx, by, bz)
    !$acc update device(ux_old, uy_old, uz_old)
    !$acc update device(bx_old, by_old, bz_old)
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
    if(init_type == 'OT3') then
      call init_OT3
    endif
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
!! @date    29 Dec 2018
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
!! @date    16 Feb 2021
!! @brief   Single mode initialization
!-----------------------------------------------!
  subroutine init_single_mode
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky, nky_local, nkz
    use grid, only: xx, yy, zz
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0, proc_id
    use params, only: zi
    use mp, only: sum_allreduce, sum_reduce
    use params, only: inputfile
    use file, only: get_unused_unit
    use time, only: microsleep
    use cuFFTmp, only: ftran_r2c, btran_c2r
    implicit none
    real(8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r
    real(8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r
    real(8), allocatable, dimension(:,:,:) :: src
    integer :: seedsize
    integer, allocatable :: seed(:)
    real(8) :: rms, kmin(3), kmax(3)
    real(8) :: u1(3), b1(3)
    character(100) :: nf_or_znf
    integer :: i, j, k
    character(len=20)  :: filename

    integer  :: unit, ierr
    namelist /initial_condition_params/ kmin, kmax, b0, b1, u1, nf_or_znf

    if(proc0) then
      print '("Single mode initialization")'
    endif

    allocate(src(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(ux_r, source=src)
    allocate(uy_r, source=src)
    allocate(uz_r, source=src)
    allocate(bx_r, source=src)
    allocate(by_r, source=src)
    allocate(bz_r, source=src)
    !$acc enter data create(ux_r)
    !$acc enter data create(uy_r)
    !$acc enter data create(uz_r)
    !$acc enter data create(bx_r)
    !$acc enter data create(by_r)
    !$acc enter data create(bz_r)
    deallocate(src)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v       Set real first then R2C and C2R        v
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

    !$acc data present(xx, ux_r, uy_r, uz_r, bx_r, by_r, bz_r)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          ux_r(k, j, i) = 1*dcos(xx(i))
          uy_r(k, j, i) = 2*dcos(xx(i))
          uz_r(k, j, i) = 3*dcos(xx(i))
          bx_r(k, j, i) = 4*dcos(xx(i))
          by_r(k, j, i) = 5*dcos(xx(i))
          bz_r(k, j, i) = 6*dcos(xx(i))
        end do
      end do
    end do
    !$acc end data

    !$acc update host(ux_r, uy_r, uz_r, bx_r, by_r, bz_r)
    write(filename, '(I0, "_r1.dat")') proc_id
    open(unit=10, file=filename, status="replace", action="write", form="formatted")
    do i = 1, nlx_local
      do j = 1, nly
        write(10, '(100E22.6)') xx(i), yy(j), ux_r(1, j, i), uy_r(1, j, i), uz_r(1, j, i), bx_r(1, j, i), by_r(1, j, i), bz_r(1, j, i)
      end do
    end do

    !$acc data present(ux_r, uy_r, uz_r, bx_r, by_r, bz_r, ux, uy, uz, bx, by, bz)
    call ftran_r2c(ux_r, ux)
    call ftran_r2c(uy_r, uy)
    call ftran_r2c(uz_r, uz)
    call ftran_r2c(bx_r, bx)
    call ftran_r2c(by_r, by)
    call ftran_r2c(bz_r, bz)
    !$acc end data

    !$acc data present(ux, uy, uz, bx, by, bz)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux(k,j,i) = ux(k,j,i)/ntot
          uy(k,j,i) = uy(k,j,i)/ntot
          uz(k,j,i) = uz(k,j,i)/ntot
          bx(k,j,i) = bx(k,j,i)/ntot
          by(k,j,i) = by(k,j,i)/ntot
          bz(k,j,i) = bz(k,j,i)/ntot
        enddo
      enddo
    enddo
    !$acc end data

    !$acc update host(ux, uy, uz, bx, by, bz)
    write(filename, '(I0, "_k.dat")') proc_id
    open(unit=10, file=filename, status="replace", action="write", form="formatted")
    do i = 1, nkx
      do j = 1, nky_local
        write(10, '(100E22.6)') kx(i), ky(j), abs(ux(1, j, i)), abs(uy(1, j, i)), abs(uz(1, j, i)), abs(bx(1, j, i)), abs(by(1, j, i)), abs(bz(1, j, i))
      end do
    end do

    !$acc data present(ux_r, uy_r, uz_r, bx_r, by_r, bz_r, ux, uy, uz, bx, by, bz)
    call btran_c2r(ux, ux_r)
    call btran_c2r(uy, uy_r)
    call btran_c2r(uz, uz_r)
    call btran_c2r(bx, bx_r)
    call btran_c2r(by, by_r)
    call btran_c2r(bz, bz_r)
    !$acc end data

    !$acc update host(ux_r, uy_r, uz_r, bx_r, by_r, bz_r)
    write(filename, '(I0, "_r2.dat")') proc_id
    open(unit=10, file=filename, status="replace", action="write", form="formatted")
    do i = 1, nlx_local
      do j = 1, nly
        write(10, '(100E22.6)') xx(i), yy(j), ux_r(1, j, i), uy_r(1, j, i), uz_r(1, j, i), bx_r(1, j, i), by_r(1, j, i), bz_r(1, j, i)
      end do
    end do
    !
    ! !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    ! !v      Set complex first then C2R and R2C      v
    ! !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !
    ! !$acc data present(ux, uy, uz, bx, by, bz)
    ! !$acc parallel loop collapse(3)
    ! do i = 1, nkx
    !   do j = 1, nky_local
    !     do k = 1, nkz
    !       if(kx(i) == 0.d0 .and. ky(j) == 1.d0 .and. kz(k) == 0) then
    !         ux(k, j, i) = 1.d0
    !         uy(k, j, i) = 2.d0
    !         uz(k, j, i) = 3.d0
    !         bx(k, j, i) = 4.d0
    !         by(k, j, i) = 5.d0
    !         bz(k, j, i) = 6.d0
    !       endif
    !     end do
    !   end do
    ! end do
    !
    ! !$acc update host(ux, uy, uz, bx, by, bz)
    ! write(filename, '(I0, "_k1.dat")') proc_id
    ! open(unit=10, file=filename, status="replace", action="write", form="formatted")
    ! do i = 1, nkx
    !   do j = 1, nky_local
    !     write(10, '(100E22.6)') kx(i), ky(j), abs(ux(1, j, i)), abs(uy(1, j, i)), abs(uz(1, j, i)), abs(bx(1, j, i)), abs(by(1, j, i)), abs(bz(1, j, i))
    !   end do
    ! end do
    ! !$acc end data
    !
    ! !$acc data present(ux_r, uy_r, uz_r, bx_r, by_r, bz_r, ux, uy, uz, bx, by, bz)
    ! call btran_c2r(ux, ux_r)
    ! call btran_c2r(uy, uy_r)
    ! call btran_c2r(uz, uz_r)
    ! call btran_c2r(bx, bx_r)
    ! call btran_c2r(by, by_r)
    ! call btran_c2r(bz, bz_r)
    ! !$acc end data
    !
    ! !$acc update host(ux_r, uy_r, uz_r, bx_r, by_r, bz_r)
    ! write(filename, '(I0, "_r1.dat")') proc_id
    ! open(unit=10, file=filename, status="replace", action="write", form="formatted")
    ! do i = 1, nlx_local
    !   do j = 1, nly
    !     write(10, '(100E22.6)') xx(i), yy(j), ux_r(1, j, i), uy_r(1, j, i), uz_r(1, j, i), bx_r(1, j, i), by_r(1, j, i), bz_r(1, j, i)
    !   end do
    ! end do
    !
    ! !$acc data present(ux_r, uy_r, uz_r, bx_r, by_r, bz_r, ux, uy, uz, bx, by, bz)
    ! call ftran_r2c(ux_r, ux)
    ! call ftran_r2c(uy_r, uy)
    ! call ftran_r2c(uz_r, uz)
    ! call ftran_r2c(bx_r, bx)
    ! call ftran_r2c(by_r, by)
    ! call ftran_r2c(bz_r, bz)
    ! !$acc end data
    !
    ! !$acc data present(ux, uy, uz, bx, by, bz)
    ! !$acc kernels
    ! ux = ux/ntot
    ! uy = uy/ntot
    ! uz = uz/ntot
    ! bx = bx/ntot
    ! by = by/ntot
    ! bz = bz/ntot
    ! !$acc end kernels
    ! !$acc end data
    !
    ! !$acc update host(ux, uy, uz, bx, by, bz)
    ! write(filename, '(I0, "_k2.dat")') proc_id
    ! open(unit=10, file=filename, status="replace", action="write", form="formatted")
    ! do i = 1, nkx
    !   do j = 1, nky_local
    !     write(10, '(100E22.6)') kx(i), ky(j), abs(ux(1, j, i)), abs(uy(1, j, i)), abs(uz(1, j, i)), abs(bx(1, j, i)), abs(by(1, j, i)), abs(bz(1, j, i))
    !   end do
    ! end do
    !
    ! !$acc data present(ux_r, uy_r, uz_r, bx_r, by_r, bz_r, ux, uy, uz, bx, by, bz)
    ! call btran_c2r(ux, ux_r)
    ! call btran_c2r(uy, uy_r)
    ! call btran_c2r(uz, uz_r)
    ! call btran_c2r(bx, bx_r)
    ! call btran_c2r(by, by_r)
    ! call btran_c2r(bz, bz_r)
    ! !$acc end data
    !
    ! !$acc update host(ux_r, uy_r, uz_r, bx_r, by_r, bz_r)
    ! write(filename, '(I0, "_r2.dat")') proc_id
    ! open(unit=10, file=filename, status="replace", action="write", form="formatted")
    ! do i = 1, nlx_local
    !   do j = 1, nly
    !     write(10, '(100E22.6)') xx(i), yy(j), ux_r(1, j, i), uy_r(1, j, i), uz_r(1, j, i), bx_r(1, j, i), by_r(1, j, i), bz_r(1, j, i)
    !   end do
    ! end do


    !$acc exit data delete(ux_r)
    !$acc exit data delete(uy_r)
    !$acc exit data delete(uz_r)
    !$acc exit data delete(bx_r)
    !$acc exit data delete(by_r)
    !$acc exit data delete(bz_r)
    deallocate(ux_r )
    deallocate(uy_r )
    deallocate(uz_r )
    deallocate(bx_r )
    deallocate(by_r )
    deallocate(bz_r )
  end subroutine init_single_mode
!
!
!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Random initialization
!-----------------------------------------------!
  subroutine init_random
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky, nky_local, nkz
    use grid, only: xx
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0, proc_id
    use params, only: zi
    use mp, only: sum_allreduce, sum_reduce
    use params, only: inputfile
    use file, only: get_unused_unit
    use time, only: microsleep
    use cuFFTmp, only: ftran_r2c, btran_c2r
    implicit none
    real(8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r
    real(8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r
    real(8), allocatable, dimension(:,:,:) :: src
    integer :: seedsize
    integer, allocatable :: seed(:)
    real(8) :: rms, kmin(3), kmax(3)
    real(8) :: u1(3), b1(3)
    character(100) :: nf_or_znf
    integer :: i, j, k

    integer  :: unit, ierr
    namelist /initial_condition_params/ kmin, kmax, b0, b1, u1, nf_or_znf

    if(proc0) then
      print *, 'Random initialization'
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                read inputfile               v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    kmin = (/0.d0, 0.d0, 0.d0/)
    kmax = (/maxval(kx), maxval(ky), maxval(kz)/)
    b0  = (/0.d0, 0.d0, 0.d0/)
    b1   = 0.d0
    u1   = 0.d0
    nf_or_znf = 'nf'

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')

    read(unit,nml=initial_condition_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    !$acc enter data copyin(kmin, kmax, b0, b1, u1)

    allocate(src(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(ux_r, source=src)
    allocate(uy_r, source=src)
    allocate(uz_r, source=src)
    allocate(bx_r, source=src)
    allocate(by_r, source=src)
    allocate(bz_r, source=src)
    !$acc enter data create(ux_r)
    !$acc enter data create(uy_r)
    !$acc enter data create(uz_r)
    !$acc enter data create(bx_r)
    !$acc enter data create(by_r)
    !$acc enter data create(bz_r)
    deallocate(src)

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

    call random_number(ux_r)
    call random_number(uy_r)
    call random_number(uz_r)
    call random_number(bx_r)
    call random_number(by_r)
    call random_number(bz_r)
    !$acc update device(ux_r)
    !$acc update device(uy_r)
    !$acc update device(uz_r)
    !$acc update device(bx_r)
    !$acc update device(by_r)
    !$acc update device(bz_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             compute r2c transform           v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !$acc data present(ux_r, uy_r, uz_r, bx_r, by_r, bz_r, ux, uy, uz, bx, by, bz)
    call ftran_r2c(ux_r, ux)
    call ftran_r2c(uy_r, uy)
    call ftran_r2c(uz_r, uz)
    call ftran_r2c(bx_r, bx)
    call ftran_r2c(by_r, by)
    call ftran_r2c(bz_r, bz)
    !$acc end data

    !$acc data present(ux, uy, uz, bx, by, bz)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux(k,j,i) = ux(k,j,i)/ntot
          uy(k,j,i) = uy(k,j,i)/ntot
          uz(k,j,i) = uz(k,j,i)/ntot
          bx(k,j,i) = bx(k,j,i)/ntot
          by(k,j,i) = by(k,j,i)/ntot
          bz(k,j,i) = bz(k,j,i)/ntot
        enddo
      enddo
    enddo
    !$acc end data

    !$acc data present(nkx, nky_local, nkz, kx, ky, kz, kmin, kmax, ux, uy, uz, bx, by, bz)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ! remove mean fields
          if (kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
            ux(k, j, i) = 0.d0
            uy(k, j, i) = 0.d0
            uz(k, j, i) = 0.d0
            bx(k, j, i) = 0.d0
            by(k, j, i) = 0.d0
            bz(k, j, i) = 0.d0
          endif

          if(nkx > 1)then
            if (kx(i)**2 <= kmin(1)**2) then
              ux(k, j, i) = 0.d0
              uy(k, j, i) = 0.d0
              uz(k, j, i) = 0.d0
              bx(k, j, i) = 0.d0
              by(k, j, i) = 0.d0
              bz(k, j, i) = 0.d0
            endif
            if (kx(i)**2 >= kmax(2)**2) then
              ux(k, j, i) = 0.d0
              uy(k, j, i) = 0.d0
              uz(k, j, i) = 0.d0
              bx(k, j, i) = 0.d0
              by(k, j, i) = 0.d0
              bz(k, j, i) = 0.d0
            endif
          endif

          if(nky > 1)then
            if (ky(j)**2 <= kmin(1)**2) then
              ux(k, j, i) = 0.d0
              uy(k, j, i) = 0.d0
              uz(k, j, i) = 0.d0
              bx(k, j, i) = 0.d0
              by(k, j, i) = 0.d0
              bz(k, j, i) = 0.d0
            endif
            if (ky(j)**2 >= kmax(3)**2) then
              ux(k, j, i) = 0.d0
              uy(k, j, i) = 0.d0
              uz(k, j, i) = 0.d0
              bx(k, j, i) = 0.d0
              by(k, j, i) = 0.d0
              bz(k, j, i) = 0.d0
            endif
          endif

          if(nkz > 1)then
            if (kz(k)**2 <= kmin(2)**2) then
              ux(k, j, i) = 0.d0
              uy(k, j, i) = 0.d0
              uz(k, j, i) = 0.d0
              bx(k, j, i) = 0.d0
              by(k, j, i) = 0.d0
              bz(k, j, i) = 0.d0
            endif

            if (kz(k)**2 >= kmax(3)**2) then
              ux(k, j, i) = 0.d0
              uy(k, j, i) = 0.d0
              uz(k, j, i) = 0.d0
              bx(k, j, i) = 0.d0
              by(k, j, i) = 0.d0
              bz(k, j, i) = 0.d0
            endif
          endif

        end do
      end do
    end do
    !$acc end data

    !$acc update host(ux, uy, uz, bx, by, bz)
    call div_free(ux, uy, uz)
    call div_free(bx, by, bz)
    !$acc update device(ux, uy, uz, bx, by, bz)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v      normalize ux, uy, uz by rms of |u|     v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !$acc data present(ux_r, uy_r, uz_r, ux, uy, uz)
    call btran_c2r(ux, ux_r)
    call btran_c2r(uy, uy_r)
    call btran_c2r(uz, uz_r)
    !$acc end data
    
    !$acc update host(ux_r, uy_r, uz_r)
    rms = sum(ux_r**2 + uy_r**2 + uz_r**2)
    call sum_allreduce(rms)
    rms = sqrt(rms/ntot)
    !$acc enter data copyin(rms)

    !$acc data present(ux_r, uy_r, uz_r, u1, rms)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          ux_r(k,j,i) = u1(1)*ux_r(k,j,i)/rms
          uy_r(k,j,i) = u1(2)*uy_r(k,j,i)/rms
          uz_r(k,j,i) = u1(3)*uz_r(k,j,i)/rms
        end do
      end do
    end do
    !$acc end data
    
    !$acc data present(ux_r, uy_r, uz_r, ux, uy, uz)
    call ftran_r2c(ux_r, ux)
    call ftran_r2c(uy_r, uy)
    call ftran_r2c(uz_r, uz)
    !$acc end data

    !$acc data present(ux, uy, uz, ntot)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux(k,j,i) = ux(k,j,i)/ntot
          uy(k,j,i) = uy(k,j,i)/ntot
          uz(k,j,i) = uz(k,j,i)/ntot
        enddo
      enddo
    enddo
    !$acc end data

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v      normalize bx, by, bz by rms of |b|     v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !$acc data present(bx_r, by_r, bz_r, bx, by, bz)
    call btran_c2r(bx, bx_r)
    call btran_c2r(by, by_r)
    call btran_c2r(bz, bz_r)
    !$acc end data
    
    !$acc update host(bx_r, by_r, bz_r)
    rms = sum(bx_r**2 + by_r**2 + bz_r**2)
    call sum_allreduce(rms)
    rms = sqrt(rms/ntot)
    !$acc update device(rms)

    if (nf_or_znf == 'nf') then
      !$acc data present(bx_r, by_r, bz_r, b1, b0, rms)
      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz_padded
            bx_r(k,j,i) = b1(1)*bx_r(k,j,i)/rms + b0(1)
            by_r(k,j,i) = b1(2)*by_r(k,j,i)/rms + b0(2)
            bz_r(k,j,i) = b1(3)*bz_r(k,j,i)/rms + b0(3)
          end do
        end do
      end do
      !$acc end data
    elseif (nf_or_znf == 'znf') then
      !$acc data present(bx_r, by_r, bz_r, b1, b0, rms, nlx_local, nly, nlz_padded, xx)
      !$acc parallel loop collapse(3)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz_padded
            bx_r(k, j, i) = b1(1)*bx_r(k, j, i)/rms
            by_r(k, j, i) = b1(2)*by_r(k, j, i)/rms
            bz_r(k, j, i) = b1(3)*bz_r(k, j, i)/rms + b0(3)*sin(xx(i))
          end do
        end do
      end do
      !$acc end data
    endif

    !$acc data present(bx_r, by_r, bz_r, bx, by, bz)
    call ftran_r2c(bx_r, bx)
    call ftran_r2c(by_r, by)
    call ftran_r2c(bz_r, bz)
    !$acc end data
    
    !$acc data present(bx, by, bz, ntot)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          bx(k,j,i) = bx(k,j,i)/ntot
          by(k,j,i) = by(k,j,i)/ntot
          bz(k,j,i) = bz(k,j,i)/ntot
        enddo
      enddo
    enddo
    !$acc end data

    call is_div_free('u', ux, uy, uz)
    call is_div_free('b', bx, by, bz)

    !$acc exit data delete(ux_r)
    !$acc exit data delete(uy_r)
    !$acc exit data delete(uz_r)
    !$acc exit data delete(bx_r)
    !$acc exit data delete(by_r)
    !$acc exit data delete(bz_r)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)

    !$acc exit data delete(kmin, kmax, b0, b1, u1)
    !$acc exit data delete(rms)
  end subroutine init_random


!-----------------------------------------------!
!> @author  YK
!! @date    18 Feb 2021
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
    real(8), allocatable, dimension(:,:,:) :: ux_r, uy_r, bx_r, by_r
    real(8), allocatable, dimension(:,:,:) :: src
    integer :: i, j, k

    if(proc0) then
      print *, 'OT2 initialization'
    endif

    allocate(src(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(ux_r, source=src)
    allocate(uy_r, source=src)
    allocate(bx_r, source=src)
    allocate(by_r, source=src)
    deallocate(src)
    !$acc enter data create(ux_r)
    !$acc enter data create(uy_r)
    !$acc enter data create(bx_r)
    !$acc enter data create(by_r)

    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          ux_r(k, j, i) = -sin(yy(j))
          uy_r(k, j, i) =  sin(xx(i))
          bx_r(k, j, i) = -sin(yy(j))
          by_r(k, j, i) =  sin(2.d0*xx(i))
        end do
      end do
    end do
  !
  !   ! compute r2c transform
    call ftran_r2c(ux_r, ux)
    call ftran_r2c(uy_r, uy)
    call ftran_r2c(bx_r, bx)
    call ftran_r2c(by_r, by)
  !
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux(k, j, i) = ux(k, j, i) / ntot
          uy(k, j, i) = uy(k, j, i) / ntot
          bx(k, j, i) = bx(k, j, i) / ntot
          by(k, j, i) = by(k, j, i) / ntot
        end do
      end do
    end do

    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux_old(k,j,i) = ux(k,j,i)
          uy_old(k,j,i) = uy(k,j,i)
          uz_old(k,j,i) = uz(k,j,i)
          bx_old(k,j,i) = bx(k,j,i)
          by_old(k,j,i) = by(k,j,i)
          bz_old(k,j,i) = bz(k,j,i)
        enddo
      enddo
    enddo

    !$acc exit data delete(ux_r)
    !$acc exit data delete(uy_r)
    !$acc exit data delete(bx_r)
    !$acc exit data delete(by_r)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(bx_r)
    deallocate(by_r)

    ! check div free of u and b
    call is_div_free('u', ux, uy, uz)
    call is_div_free('b', bx, by, bz)

  end subroutine init_OT2


!-----------------------------------------------!
!> @author  YK
!! @date    18 Feb 2021
!! @brief   3D Orszag Tang problem initialization
!-----------------------------------------------!
  subroutine init_OT3
    use grid, only: xx, yy, zz
    use grid, only: nlx_local, nly, nlz, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0
    use cuFFTmp, only: ftran_r2c
    implicit none
    real(8), allocatable, dimension(:,:,:) :: ux_r, uy_r, bx_r, by_r
    real(8), allocatable, dimension(:,:,:) :: src
    integer :: i, j, k

    if(proc0) then
      print *, 'OT3 initialization'
    endif

    allocate(src(nlz_padded, nly, nlx_local), source=0.d0)
    allocate(ux_r, source=src)
    allocate(uy_r, source=src)
    allocate(bx_r, source=src)
    allocate(by_r, source=src)
    deallocate(src)
    !$acc enter data create(ux_r)
    !$acc enter data create(uy_r)
    !$acc enter data create(bx_r)
    !$acc enter data create(by_r)

    ! zz is dimensioned nlz, so the fill must stop at nlz; the remaining
    ! nlz_padded - nlz slots are the in-place r2c padding and are zeroed.
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          if (k <= nlz) then
            ux_r(k, j, i) = -sin(yy(j) + zz(k))
            uy_r(k, j, i) =  sin(xx(i) + zz(k))
            bx_r(k, j, i) = -sin(yy(j) + zz(k))
            by_r(k, j, i) =  sin(2.d0*xx(i) + zz(k))
          else
            ux_r(k, j, i) = 0.d0
            uy_r(k, j, i) = 0.d0
            bx_r(k, j, i) = 0.d0
            by_r(k, j, i) = 0.d0
          endif
        end do
      end do
    end do

    ! compute r2c transform
    call ftran_r2c(ux_r, ux)
    call ftran_r2c(uy_r, uy)
    call ftran_r2c(bx_r, bx)
    call ftran_r2c(by_r, by)

    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux(k, j, i) = ux(k, j, i) / ntot
          uy(k, j, i) = uy(k, j, i) / ntot
          bx(k, j, i) = bx(k, j, i) / ntot
          by(k, j, i) = by(k, j, i) / ntot
        end do
      end do
    end do

    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux_old(k,j,i) = ux(k,j,i)
          uy_old(k,j,i) = uy(k,j,i)
          uz_old(k,j,i) = uz(k,j,i)
          bx_old(k,j,i) = bx(k,j,i)
          by_old(k,j,i) = by(k,j,i)
          bz_old(k,j,i) = bz(k,j,i)
        enddo
      enddo
    enddo

    !$acc exit data delete(ux_r)
    !$acc exit data delete(uy_r)
    !$acc exit data delete(bx_r)
    !$acc exit data delete(by_r)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(bx_r)
    deallocate(by_r)

    ! check div free of u and b
    call is_div_free('u', ux, uy, uz)
    call is_div_free('b', bx, by, bz)

  end subroutine init_OT3


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Restart
!-----------------------------------------------!
  subroutine restart
    use mp, only: proc0, proc_id, comm_fft
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

    call mpiio_read_one(ux, sizes, subsizes, starts, trim(restart_dir)//'ux.dat', comm_fft)
    call mpiio_read_one(uy, sizes, subsizes, starts, trim(restart_dir)//'uy.dat', comm_fft)
    call mpiio_read_one(uz, sizes, subsizes, starts, trim(restart_dir)//'uz.dat', comm_fft)
    call mpiio_read_one(bx, sizes, subsizes, starts, trim(restart_dir)//'bx.dat', comm_fft)
    call mpiio_read_one(by, sizes, subsizes, starts, trim(restart_dir)//'by.dat', comm_fft)
    call mpiio_read_one(bz, sizes, subsizes, starts, trim(restart_dir)//'bz.dat', comm_fft)

    ux_old = ux
    uy_old = uy
    uz_old = uz
    bx_old = bx
    by_old = by
    bz_old = bz

    !$acc update device(ux, uy, uz, bx, by, bz)
    !$acc update device(ux_old, uy_old, uz_old, bx_old, by_old, bz_old)
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Enforce div free for (wx, wy, wz)
!-----------------------------------------------!
  subroutine div_free(wx, wy, wz)
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: zi
    implicit none
    complex(8), dimension (:,:,:), intent(inout) :: wx, wy, wz
    complex(8), allocatable, dimension(:,:,:)   :: nbl2inv_div_w ! nabla^-2 (div w)
    complex(8) :: k2
    integer :: i, j, k

    allocate(nbl2inv_div_w(nkz, nky_local, nkx), source=(0.d0,0.d0))

    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          k2 = kx(i)**2 + ky(j)**2 + kz(k)**2
          if(k2 == 0.d0) then
            nbl2inv_div_w(k,j,i) = 0.d0
          else
            nbl2inv_div_w(k,j,i) = -zi*(   kx(i)*wx(k,j,i) &
                                         + ky(j)*wy(k,j,i) &
                                         + kz(k)*wz(k,j,i) )/k2
          endif

          wx(k,j,i) = wx(k,j,i) - zi*kx(i)*nbl2inv_div_w(k,j,i)
          wy(k,j,i) = wy(k,j,i) - zi*ky(j)*nbl2inv_div_w(k,j,i)
          wz(k,j,i) = wz(k,j,i) - zi*kz(k)*nbl2inv_div_w(k,j,i)
        enddo
      enddo
    enddo

    deallocate(nbl2inv_div_w)

  end subroutine div_free

!
!-----------------------------------------------!
!> @author  YK
!! @date    3 Oct 2020
!! @brief   Check if divergence free is satisfied
!-----------------------------------------------!
  subroutine is_div_free(name, fx, fy, fz)
    use grid, only: kx, ky, kz, k2
    use grid, only: nkx, nky_local, nkz
    use mp, only: proc0
    use params, only: zi
    use mp, only: sum_reduce
    implicit none
    character(*) :: name
    integer :: i, j, k
    complex(8), dimension(:,:,:), intent(in) :: fx, fy, fz
    real(8) :: abs_div_f_sum, abs_f_sum, abs_k_sum

    ! The loop below runs on the host, but the fields live on the device
    ! between kernels.  Refresh here rather than at the call sites: two of
    ! the three callers used to omit it and printed NaN from uninitialized
    ! host memory.
    !$acc update host(fx, fy, fz)

    abs_div_f_sum  = 0.d0
    abs_f_sum      = 0.d0

    do i =1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          abs_div_f_sum = abs_div_f_sum + abs(kx(i)*fx(k, j, i) + ky(j)*fy(k, j, i) + kz(k)*fz(k, j, i))
          abs_f_sum     = abs_f_sum     + sqrt(abs(fx(k, j, i))**2 + abs(fy(k, j, i))**2 + abs(fz(k, j, i))**2)
        end do
      end do
    end do

    call sum_reduce(abs_div_f_sum, 0)
    call sum_reduce(abs_f_sum, 0)
    abs_k_sum = sum(sqrt(k2)); call sum_reduce(abs_k_sum, 0)

    if(proc0) write(*, "('<|k.',A2,'|>/(<|k|><|',A2,'|>) = ', es10.3)") trim(name), trim(name), abs_div_f_sum/(abs_f_sum*abs_k_sum)

  end subroutine is_div_free


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Finalization of Fields
!-----------------------------------------------!
  subroutine finish_fields
    implicit none

    deallocate(ux)
    deallocate(uy)
    deallocate(uz)
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
    deallocate(ux_old)
    deallocate(uy_old)
    deallocate(uz_old)
    deallocate(bx_old)
    deallocate(by_old)
    deallocate(bz_old)

  end subroutine finish_fields

end module fields

