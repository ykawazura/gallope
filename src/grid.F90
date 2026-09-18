!-----------------------------------------------!
!> @author  YK
!! @brief   Grid setting
!-----------------------------------------------!
module grid
  implicit none

  public  init_grid, finish_grid
  public  lx, ly, lz
  public  nlx, nly, nlz, nlx_local, nlz_padded
  public  nkx, nky, nkz, nky_local
  public  nm, nm_local, m_offset
  public  ntot
  public  xx, yy, zz, xx_global
  public  kx, ky, kz, ky_global
  public  kx_max, ky_max, kz_max
  public  k2, k2inv, k2_max
  public  kprp2, kprp2inv, kz2, kprp2_max, kz2_max
  public  ikx, iky, ikz, iky_global
  public  dlx, dly, dlz, dkx, dky, dkz
  private read_parameters

  real(8) :: lx, ly, lz
  integer :: nlx, nly, nlz, nlx_local, nlz_padded
  integer :: nkx, nky, nkz, nky_local
  ! nm      : total number of Hermite (v_parallel) moments (KRMHD and up; nm=1 for RMHD)
  ! nm_local: moments held by this rank after the comm_m split
  ! m_offset: global index of this rank's first local moment (global m runs 0..nm-1)
  integer :: nm, nm_local, m_offset
  integer :: ntot
  real(8), allocatable :: xx(:), yy(:), zz(:), xx_global(:)
  real(8), allocatable :: kx(:), ky(:), kz(:), ky_global(:)
  real(8), allocatable :: k2(:, :, :), k2inv(:, :, :), kprp2(:, :, :), kprp2inv(:, :, :), kz2(:)
  integer, allocatable :: ikx(:), iky(:), ikz(:), iky_global(:)
  real(8) :: dlx, dly, dlz, dkx, dky, dkz
  real(8) :: kx_max, ky_max, kz_max
  real(8) :: kprp2_max, kz2_max, k2_max
  !$acc declare create(lx, ly, lz)
  !$acc declare create(nlx, nly, nlz, nlx_local, nlz_padded)
  !$acc declare create(nkx, nky, nkz, nky_local)
  !$acc declare create(nm, nm_local, m_offset)
  !$acc declare create(ntot)
  !$acc declare create(xx, yy, zz)
  !$acc declare create(kx, ky, kz)
  !$acc declare create(kx_max, ky_max, kz_max)
  !$acc declare create(k2, k2inv, k2_max)
  !$acc declare create(kprp2, kprp2inv, kz2, kprp2_max, kz2_max)
  !$acc declare create(ikx, iky, ikz)
  !$acc declare create(dlx, dly, dlz, dkx, dky, dkz)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of grid
!-----------------------------------------------!
  subroutine init_grid
    use mp, only: proc0, sum_allreduce
    use mp, only: proc_id, nproc
    ! the real/spectral slabs are decomposed over comm_fft, not the world.
    ! With P_m=P_s=1, nproc_fft==nproc and iproc_fft==proc_id (bitwise regression).
    use mp, only: iproc_fft, nproc_fft
    ! Hermite (v_parallel) moments are decomposed over comm_m.
    use mp, only: iproc_m, nproc_m
    use params, only: pi, inputfile
    implicit none
    include 'mpif.h'
    integer :: ranks_cutoff
    integer :: i, j, k, ierr

    call read_parameters(inputfile)

    nlz_padded = 2*(nlz/2 + 1)

    nkx = nlx
    nky = nly
    nkz = nlz/2 + 1

    ! We start with X-Slabs
    ! Ranks 0 ... (nlx % nproc_fft - 1) have 1 more element in the X dimension
    ! and every rank own all elements in the Y and Z dimensions.
    ranks_cutoff = mod(nlx, nproc_fft)
    nlx_local = nlx/nproc_fft
    if (iproc_fft < ranks_cutoff) nlx_local = nlx_local + 1

    nky_local =  nly / nproc_fft;

    ! Hermite-moment decomposition over comm_m (mirror of the nky_local split).
    ! Require an even divide so m_offset = nm_local*iproc_m is exact, which keeps
    ! the comm_m halo and the collective rank-4 g I/O consistent across ranks.
    if (mod(nm, nproc_m) > 0) then
      if (proc0) print*," nm has to divide evenly by the comm_m rank count (P_m)"
      call mpi_finalize(ierr)
      stop
    end if
    nm_local = nm / nproc_m
    m_offset = nm_local * iproc_m

    allocate(xx(nlx_local))
    allocate(xx_global(nlx))
    allocate(yy(nly))
    allocate(zz(nlz))
    allocate(kx(nkx))
    allocate(ky(nky_local))
    allocate(ky_global(nky))
    allocate(kz(nkz))
    allocate(ikx(nkx))
    allocate(iky(nky_local))
    allocate(iky_global(nky))
    allocate(ikz(nkz))
    allocate(kz2(nkz))

    do i = 1, nlx_local
      xx(i) = lx*dble(1.d0/nlx*(i + nlx_local*iproc_fft - 1))
    enddo
    do i = 1, nlx
      xx_global(i) = lx*dble(1.d0/nlx*(i - 1))
    enddo
    do j = 1, nly
      yy(j) = ly*dble(1.d0/nly*(j-1))
    enddo
    do k = 1, nlz
      zz(k) = lz*dble(1.d0/nlz*(k-1))
    enddo

    do i = 1, nkx
      if (i <= nkx/2 + 1) then
        ikx(i) = i - 1
      else
        ikx(i) = i - nkx - 1
      endif

      kx(i) = 2.d0*pi*dble(ikx(i))/lx
    enddo
    do j = 1, nky_local
      if (j + nky_local*iproc_fft <= nky/2 + 1) then
        iky(j) = j + nky_local*iproc_fft - 1
      else
        iky(j) = j + nky_local*iproc_fft - nky - 1
      endif

      ky(j) = 2.d0*pi*dble(iky(j))/ly
    enddo
    do j = 1, nky
      if (j <= nky/2 + 1) then
        iky_global(j) = j - 1
      else
        iky_global(j) = j - nky - 1
      endif

      ky_global(j) = 2.d0*pi*dble(iky_global(j))/ly
    enddo
    do k = 1, nkz
      ikz(k) = k - 1
      kz(k) = 2.d0*pi*dble(ikz(k))/lz
      kz2(k) = kz(k)**2
    enddo

    if (mod(nly, nproc_fft) > 0) then
      print*," nly has to divide evenly by the comm_fft rank count"
      call mpi_finalize(ierr)
    end if

    dlx = abs(xx_global(2) - xx_global(1))
    dly = abs(yy(2) - yy(1))
    dlz = abs(zz(2) - zz(1))
    dkx = abs(kx(2) - kx(1))
    dky = abs(ky_global(2) - ky_global(1))
    dkz = abs(kz(2) - kz(1))

    kx_max = 2.d0*pi*dble(nlx/2)/lx
    ky_max = 2.d0*pi*dble(nly/2)/ly
    kz_max = 2.d0*pi*dble(nlz/2)/lz

    kprp2_max = maxval(abs(kx))**2 + maxval(abs(ky_global))**2
    kz2_max   = maxval(kz2)
    k2_max    = maxval(abs(kx))**2 + maxval(abs(ky_global))**2 + maxval(abs(kz))**2

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(proc0) then
      print '("  nlx  = ", i6, ",    nly  = ", i6, ",    nlz  = ", i6)', nlx , nly , nlz
      print *
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call sleep(1)

!
#ifdef DEBG
    do i = 0, nproc
      if(i == proc_id) then
        print '("  proc_id = ", i6, ", nlx_local = ", i6)', proc_id, nlx_local
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call sleep(1)
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(proc0) print *
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    do i = 0, nproc
      if(i == proc_id) then
        print '("  proc_id = ", i6, ", nky_local = ", i6)', proc_id, nky_local
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call sleep(1)
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(proc0) print *
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call sleep(1)
#endif

    allocate(k2      (nkz, nky_local, nkx))
    allocate(k2inv   (nkz, nky_local, nkx))
    allocate(kprp2   (nkz, nky_local, nkx))
    allocate(kprp2inv(nkz, nky_local, nkx))

    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          kprp2(k, j, i) = kx(i)**2 + ky(j)**2
          if(kprp2(k, j, i) == 0.d0) then
            kprp2inv(k, j, i) = 0.d0
          else
            kprp2inv(k, j, i) = 1.0d0/kprp2(k, j, i)
          endif

          k2(k, j, i) = kx(i)**2 + ky(j)**2 + kz(k)**2
          if(k2(k, j, i) == 0.d0) then
            k2inv(k, j, i) = 0.d0
          else
            k2inv(k, j, i) = 1.0d0/k2(k, j, i)
          endif
        enddo
      enddo
    enddo

    !$acc update device(lx, ly, lz)
    !$acc update device(nlx, nly, nlz, nlx_local, nlz_padded)
    !$acc update device(nkx, nky, nkz, nky_local)
    !$acc update device(nm, nm_local, m_offset)
    !$acc update device(ntot)
    !$acc update device(xx, yy, zz)
    !$acc update device(kx, ky, kz)
    !$acc update device(kx_max, ky_max, kz_max)
    !$acc update device(k2, k2inv, k2_max)
    !$acc update device(kprp2, kprp2inv, kz2, kprp2_max, kz2_max)
    !$acc update device(ikx, iky, ikz)
    !$acc update device(dlx, dly, dlz, dkx, dky, dkz)
  end subroutine init_grid


!-----------------------------------------------!
!> @author  YK
!! @brief   Read inputfile for box parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use mp, only: nproc, proc0, proc_id
    use file, only: get_unused_unit
    use params, only: pi
    implicit none
    
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /box_parameters/ lx, ly, lz, nlx, nly, nlz, nm

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    lx = 1.d0
    ly = 1.d0
    lz = 1.d0

    nlx = 32
    nly = 32
    nlz = 32

    nm  = 1     ! single moment (RMHD limit); KRMHD sets nm>1 in the input file
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=box_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading box_parameters failed"
    close(unit)

    lx = 2.d0*pi*lx
    ly = 2.d0*pi*ly
    lz = 2.d0*pi*lz

    ntot = nlx*nly*nlz

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @brief   Finalization of grid
!-----------------------------------------------!
  subroutine finish_grid
    implicit none

    deallocate(xx)
    deallocate(xx_global)
    deallocate(yy)
    deallocate(zz)
    deallocate(kx)
    deallocate(ky)
    deallocate(ky_global)
    deallocate(kz)
    deallocate(ikx)
    deallocate(iky)
    deallocate(iky_global)
    deallocate(ikz)
    deallocate(kz2)

    deallocate(k2)
    deallocate(k2inv)
    deallocate(kprp2)
    deallocate(kprp2inv)
  end subroutine finish_grid

end module grid

