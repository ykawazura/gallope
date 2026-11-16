!-----------------------------------------------!
!> @brief  Gate G1 of the MHD_COMP_ISOTH tests: cuFFTMp with a size-2 dimension
!!
!! Link with the objects of any Gallope build (everything but gallope.o) and
!! run with the same input file conventions as gallope:
!!   mpirun -np P ./fft_roundtrip RUNNAME     (reads RUNNAME.in)
!! Only &box_parameters matters; P must divide nlx and nly.
!!
!! Checks, with the global indices x = i + nlx_local*iproc_fft and
!! ky = j + nky_local*iproc_fft:
!!   1. round trip: btran_c2r(ftran_r2c(f))/(nlx*nly*nlz) = f on k <= nlz,
!!      with garbage in the padding rows k > nlz;
!!   2. axis mapping: for f = cos(2 pi m.x/n) + sin(2 pi m'.x/n)/2 every stored
!!      mode equals the exact DFT, N/2 [d(k-m) + d(k+m)] + N/(4i) [d(k-m') - d(k+m')]
!!      (indices mod n; the forward transform is exp(-i k.x)), for m, m' that
!!      include the Nyquist modes of every direction.
!! Prints PASS or FAIL per case; the exit status of the job script greps it.
!-----------------------------------------------!
program fft_roundtrip
  use mp, only: init_mp, init_comm_layout, finish_mp, proc0, iproc_fft, nproc_fft
  use mp, only: max_allreduce, comm_fft
  use cufftmp, only: init_cuda, init_cuFFTmp, finish_cuFFTmp, ftran_r2c, btran_c2r
  use params, only: init_params
  use grid, only: init_grid, nlx, nly, nlz, nlx_local, nlz_padded, nkx, nky_local, nkz
  implicit none
  real(8), parameter :: pi = 3.14159265358979323846d0, tol = 1.d-13
  real(8), allocatable :: ur(:,:,:), vr(:,:,:)
  complex(8), allocatable :: uk(:,:,:)
  integer, parameter :: ncase = 9
  integer :: m(3, ncase), m2(3, ncase)
  real(8) :: err, scale, ntot
  integer :: i, j, k, ic, ig, jg
  logical :: ok_all

  call init_params
  call init_mp
  call init_comm_layout(1, 1)
  call init_cuda
  call init_grid
  call init_cuFFTmp

  if (mod(nlx, nproc_fft) /= 0 .or. mod(nly, nproc_fft) /= 0) then
    if (proc0) print *, 'FAIL: the process count must divide nlx and nly'
    call finish_mp
    stop
  endif
  ntot = dble(nlx)*dble(nly)*dble(nlz)

  allocate(ur(nlz_padded, nly, nlx_local), vr(nlz_padded, nly, nlx_local))
  allocate(uk(nkz, nky_local, nkx))
  !$acc enter data create(ur, vr, uk)
  ok_all = .true.
  if (proc0) write(*, '("grid nlx, nly, nlz = ", 3i6, "   processes = ", i4)') nlx, nly, nlz, nproc_fft

  ! 1. round trip of a generic field
  do i = 1, nlx_local
    ig = i + nlx_local*iproc_fft
    do j = 1, nly
      do k = 1, nlz_padded
        if (k <= nlz) then
          ur(k, j, i) = sin(1.3d0*ig + 0.7d0*j*j + 0.37d0*k) + 0.1d0*ig/nlx
        else
          ur(k, j, i) = 1.d30   ! padding: must be ignored by the r2c transform
        endif
      enddo
    enddo
  enddo
  !$acc update device(ur)
  call ftran_r2c(ur, uk)
  call btran_c2r(uk, vr)
  !$acc update self(vr)
  err = 0.d0
  scale = 0.d0
  do i = 1, nlx_local
    do j = 1, nly
      do k = 1, nlz
        err = max(err, abs(vr(k, j, i)/ntot - ur(k, j, i)))
        scale = max(scale, abs(ur(k, j, i)))
      enddo
    enddo
  enddo
  call max_allreduce(err, comm_fft)
  call max_allreduce(scale, comm_fft)
  call report('round trip', err/scale)

  ! 2. axis mapping; m(:, ic) carries the cos mode, m2(:, ic) the sin mode
  m (:, 1) = [1, 0, 0];             m2(:, 1) = [0, 1, 1]
  m (:, 2) = [0, 1, 0];             m2(:, 2) = [1, 0, 1]
  m (:, 3) = [0, 0, 1];             m2(:, 3) = [1, 1, 0]
  m (:, 4) = [nlx/2, 0, 0];         m2(:, 4) = [1, 1, 1]
  m (:, 5) = [0, nly/2, 0];         m2(:, 5) = [nlx/2, 1, 1]
  m (:, 6) = [0, 0, nlz/2];         m2(:, 6) = [1, nly/2, 1]
  m (:, 7) = [nlx/2, nly/2, nlz/2]; m2(:, 7) = [1, 1, nlz/2]
  m (:, 8) = [1, nly/2, 1];         m2(:, 8) = [nlx/2, 1, nlz/2]
  m (:, 9) = [nlx/2, 1, nlz/2];     m2(:, 9) = [nlx - 1, nly - 1, 1]
  do ic = 1, ncase
    do i = 1, nlx_local
      ig = i + nlx_local*iproc_fft
      do j = 1, nly
        do k = 1, nlz_padded
          if (k <= nlz) then
            ur(k, j, i) = cos(phase(m(:, ic), ig - 1, j - 1, k - 1)) &
                        + 0.5d0*sin(phase(m2(:, ic), ig - 1, j - 1, k - 1))
          else
            ur(k, j, i) = -1.d30
          endif
        enddo
      enddo
    enddo
    !$acc update device(ur)
    call ftran_r2c(ur, uk)
    !$acc update self(uk)
    err = 0.d0
    do i = 1, nkx
      do j = 1, nky_local
        jg = j + nky_local*iproc_fft
        do k = 1, nkz
          err = max(err, abs(uk(k, j, i)/ntot - expected(i - 1, jg - 1, k - 1, ic)))
        enddo
      enddo
    enddo
    call max_allreduce(err, comm_fft)
    block
      character(len=80) :: label
      write(label, '("mapping m=(", 3i4, ") m2=(", 3i4, ")")') m(:, ic), m2(:, ic)
      call report(trim(label), err)
    end block
  enddo

  if (proc0) then
    if (ok_all) then
      print *, 'ALL PASS'
    else
      print *, 'SOME FAIL'
    endif
  endif

  !$acc exit data delete(ur, vr, uk)
  deallocate(ur, vr, uk)
  call finish_cuFFTmp
  call finish_mp

contains

  real(8) function phase(mm, x, y, z)
    integer, intent(in) :: mm(3), x, y, z
    phase = 2.d0*pi*(dble(mm(1))*x/nlx + dble(mm(2))*y/nly + dble(mm(3))*z/nlz)
  end function phase

  logical function hit(kk, mm)
    integer, intent(in) :: kk(3), mm(3)
    hit = mod(kk(1) - mm(1), nlx) == 0 .and. mod(kk(2) - mm(2), nly) == 0 .and. &
          mod(kk(3) - mm(3), nlz) == 0
  end function hit

  !> exact DFT / ntot of cos(m.x) + sin(m2.x)/2 at the stored mode (kx, ky, kz) >= 0
  complex(8) function expected(kx, ky, kz, ic)
    integer, intent(in) :: kx, ky, kz, ic
    integer :: kk(3)
    kk = [kx, ky, kz]
    expected = (0.d0, 0.d0)
    if (hit(kk,  m(:, ic))) expected = expected + 0.5d0
    if (hit(kk, -m(:, ic))) expected = expected + 0.5d0
    if (hit(kk,  m2(:, ic))) expected = expected + 0.5d0/(2.d0*(0.d0, 1.d0))
    if (hit(kk, -m2(:, ic))) expected = expected - 0.5d0/(2.d0*(0.d0, 1.d0))
  end function expected

  subroutine report(label, e)
    character(len=*), intent(in) :: label
    real(8), intent(in) :: e
    if (e > tol) ok_all = .false.
    if (proc0) write(*, '(a, t48, es10.2, 2x, a)') label, e, merge('PASS', 'FAIL', e <= tol)
  end subroutine report

end program fft_roundtrip
