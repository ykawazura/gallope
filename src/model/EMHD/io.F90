!-----------------------------------------------!
!> @author  YK
!! @brief   IO for EMHD
!!
!!          The time series, the binned spectra and
!!          the 2D cuts are handled as lists: each
!!          entry has an index parameter (is_*, ib_*,
!!          i2d_*) and a name. The diagnostics fill
!!          the arrays through the index parameters,
!!          so a quantity cannot be written under the
!!          wrong name. Names follow Calliope EMHD
!!          where Calliope has the quantity.
!-----------------------------------------------!
module io
  use netcdf
  use MPI
  implicit none

  public :: init_io, finish_io, loop_io, loop_io_2D, loop_io_3D, save_restart

  private

  !--------------------------------------------------!
  ! Time series <runname>.out.nc: '<name>' (tt)
  !--------------------------------------------------!
  integer, parameter, public :: is_wmag   = 1, is_wmag_dot = 2, is_wmag_dissip = 3, &
                                is_p_ext  = 4, is_hel      = 5, is_divb_rel    = 6
  integer, parameter, public :: nseries = 6

  character(len=20), parameter :: series_names(nseries) = [character(len=20) :: &
    'wmag_sum', 'wmag_dot_sum', 'wmag_dissip_sum', 'p_ext_sum', 'hel_sum', 'divb_rel' ]

  !--------------------------------------------------!
  ! Isotropic spectra <runname>.out.nc: '<name>' (kpbin, tt)
  !--------------------------------------------------!
  integer, parameter, public :: ib_b2 = 1, ib_bx2 = 2, ib_by2 = 3, ib_bz2 = 4, ib_wmag = 5
  integer, parameter, public :: nbins = 5

  character(len=20), parameter :: bin_names(nbins) = [character(len=20) :: &
    'b2_bin', 'bx2_bin', 'by2_bin', 'bz2_bin', 'wmag_bin' ]

  !--------------------------------------------------!
  ! 2D cuts out2d/<name>_r_{z0,x0,y0}.dat
  !--------------------------------------------------!
  integer, parameter, public :: i2d_bx = 1, i2d_by = 2, i2d_bz = 3, &
                                i2d_jx = 4, i2d_jy = 5, i2d_jz = 6
  integer, parameter, public :: n2d = 6

  character(len=2), parameter :: cut2d_names(n2d) = ['bx', 'by', 'bz', 'jx', 'jy', 'jz']

  ! cut planes, in the order of the last index of the fh_2d/disp_2d arrays
  integer, parameter :: iz0 = 1, ix0 = 2, iy0 = 3
  character(len=2), parameter :: plane_names(3) = ['z0', 'x0', 'y0']

  !--------------------------------------------------!
  ! 3D fields out3d/<name>.dat and restart files
  !--------------------------------------------------!
  integer, parameter :: n3d = 3
  character(len=2), parameter :: out3d_names(n3d) = ['bx', 'by', 'bz']

  ! MPIIO for 3D
  integer :: fh_3d(n3d)
  integer (kind=MPI_OFFSET_KIND) :: disp_3d(n3d)
  integer :: out3d_time_unit

  ! MPIIO for 2D
  integer :: fh_2d(n2d, 3)
  integer (kind=MPI_OFFSET_KIND) :: disp_2d(n2d, 3)
  integer :: out2d_time_unit

  ! NETCDF for regular output file
  character(len=100) :: filename
  integer :: status
  integer, parameter :: kind_nf = kind (NF90_NOERR)
  integer (kind_nf) :: ncid
  integer :: run_id
  integer (kind_nf) :: char10_dim
  ! parameter
  integer, parameter :: nparams = 3
  integer :: param_id(nparams)
  ! coordinate
  integer :: xx_id, yy_id, zz_id, kx_id, ky_id, kz_id, kpbin_id, tt_id
  ! time series, mean field and spectra
  integer :: series_id(nseries), b0_id, bin_id(nbins)

  integer (kind_nf) :: xx_dim, yy_dim, zz_dim, kx_dim, ky_dim, kz_dim, kpbin_dim, tt_dim, vec_dim
  integer, dimension (2) :: bin_dim, b0_dim

  integer :: nout

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of IO
!-----------------------------------------------!
  subroutine init_io(nkpolar, kpbin)
    implicit none
    integer, intent(in) :: nkpolar
    real(8), intent(in) :: kpbin(1:nkpolar)

    call init_io_decomp_2d
    call init_io_decomp_3d
    call init_io_netcdf(nkpolar, kpbin)
  end subroutine init_io


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of MPIIO for 3D
!-----------------------------------------------!
  subroutine init_io_decomp_3d
    use mp, only: proc0, comm_fft
    use file, only: open_output_file
    implicit none
    integer :: l

    ! EMHD runs with P_m = P_s = 1 (model_specific), so comm_fft spans all ranks.
    do l = 1, n3d
      call set_file_handle('out3d/'//trim(out3d_names(l))//'.dat', fh_3d(l), disp_3d(l), comm_fft)
    enddo

    if(proc0) then
      call open_output_file (out3d_time_unit, 'out3d/time.dat')
    endif
  end subroutine init_io_decomp_3d


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of MPIIO for 2D
!-----------------------------------------------!
  subroutine init_io_decomp_2d
    use mp, only: proc0, comm_fft
    use file, only: open_output_file
    implicit none
    integer :: l, p

    do l = 1, n2d
      do p = 1, 3
        call set_file_handle('out2d/'//trim(cut2d_names(l))//'_r_'//plane_names(p)//'.dat', &
                             fh_2d(l, p), disp_2d(l, p), comm_fft)
      enddo
    enddo

    if(proc0) then
      call open_output_file (out2d_time_unit, 'out2d/time.dat')
    endif
  end subroutine init_io_decomp_2d


!-----------------------------------------------!
!> @author  YK
!! @brief   Set file handle for MPIIO
!-----------------------------------------------!
  subroutine set_file_handle(fn, fh, disp, comm)
    implicit none
    character(*) :: fn
    integer :: fh
    integer (kind=MPI_OFFSET_KIND) :: disp
    integer, intent(in) :: comm
    integer :: ierr

    call MPI_FILE_OPEN(comm, trim(fn), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_FILE_SET_SIZE(fh, 0_MPI_OFFSET_KIND, ierr)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
  end subroutine set_file_handle


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of NETCDF
!-----------------------------------------------!
  subroutine init_io_netcdf(nkpolar, kpbin)
    use grid, only: xx_global, yy, zz, kx, ky_global, kz
    use mp, only: proc0
    use params, only: runname, eta, eta_exp, de
    implicit none
    integer, intent(in) :: nkpolar
    real(8), intent(in) :: kpbin(1:nkpolar)
    character(len=12), parameter :: param_names(nparams) = [character(len=12) :: &
      'eta', 'eta_exp', 'de' ]
    real(8) :: param_values(nparams)
    integer :: l

    if(proc0) then
      ! exponents are written as doubles, as in Calliope
      param_values = [ eta, dble(eta_exp), de ]

      !--------------------------------------------------!
      ! Output for parameters, time history, and spectra
      !--------------------------------------------------!
      filename = trim(runname)//'.out.nc' ! File name
      status = nf90_create (filename, NF90_CLOBBER, ncid)

      status = nf90_put_att (ncid, NF90_GLOBAL, 'title', 'calliope simulation data')
      status = nf90_def_dim (ncid, 'char10', 10, char10_dim)
      status = nf90_def_var (ncid, 'run_info', NF90_CHAR, char10_dim, run_id)
      status = nf90_put_att (ncid, run_id, 'model', _MODEL_)

      status = nf90_def_dim (ncid, 'yy', size(yy), yy_dim)
      status = nf90_def_dim (ncid, 'xx', size(xx_global), xx_dim)
      status = nf90_def_dim (ncid, 'zz', size(zz), zz_dim)
      status = nf90_def_dim (ncid, 'kx', size(kx), kx_dim)
      status = nf90_def_dim (ncid, 'ky', size(ky_global), ky_dim)
      status = nf90_def_dim (ncid, 'kz', size(kz), kz_dim)
      status = nf90_def_dim (ncid, 'kpbin', size(kpbin), kpbin_dim)
      status = nf90_def_dim (ncid, 'tt', NF90_UNLIMITED, tt_dim)
      status = nf90_def_dim (ncid, '3', 3, vec_dim)

      do l = 1, nparams
        status = nf90_def_var (ncid, trim(param_names(l)), NF90_DOUBLE, param_id(l))
      enddo

      status = nf90_def_var (ncid, 'xx', NF90_DOUBLE, xx_dim, xx_id)
      status = nf90_def_var (ncid, 'yy', NF90_DOUBLE, yy_dim, yy_id)
      status = nf90_def_var (ncid, 'zz', NF90_DOUBLE, zz_dim, zz_id)
      status = nf90_def_var (ncid, 'kx', NF90_DOUBLE, kx_dim, kx_id)
      status = nf90_def_var (ncid, 'ky', NF90_DOUBLE, ky_dim, ky_id)
      status = nf90_def_var (ncid, 'kz', NF90_DOUBLE, kz_dim, kz_id)
      status = nf90_def_var (ncid, 'kpbin', NF90_DOUBLE, kpbin_dim, kpbin_id)
      status = nf90_def_var (ncid, 'tt', NF90_DOUBLE, tt_dim, tt_id)

      ! time series
      do l = 1, nseries
        status = nf90_def_var (ncid, trim(series_names(l)), NF90_DOUBLE, tt_dim, series_id(l))
      enddo

      ! mean (guide) field = the k = 0 mode of b
      b0_dim (1) = vec_dim
      b0_dim (2) = tt_dim
      status = nf90_def_var (ncid, 'b0', NF90_DOUBLE, b0_dim, b0_id)

      ! isotropic spectrum
      bin_dim (1) = kpbin_dim
      bin_dim (2) = tt_dim
      do l = 1, nbins
        status = nf90_def_var (ncid, trim(bin_names(l)), NF90_DOUBLE, bin_dim, bin_id(l))
      enddo

      status = nf90_enddef (ncid)  ! out of definition mode

      do l = 1, nparams
        status = nf90_put_var (ncid, param_id(l), param_values(l))
      enddo

      status = nf90_put_var (ncid, xx_id, xx_global)
      status = nf90_put_var (ncid, yy_id, yy)
      status = nf90_put_var (ncid, zz_id, zz)
      status = nf90_put_var (ncid, kx_id, kx)
      status = nf90_put_var (ncid, ky_id, ky_global)
      status = nf90_put_var (ncid, kz_id, kz)
      status = nf90_put_var (ncid, kpbin_id, kpbin)

      nout = 1
    endif
  end subroutine init_io_netcdf


!-----------------------------------------------!
!> @author  YK
!! @brief   Append variables to NETCDF
!           for time history & spectra
!!          (series, b0 and bins are filled on the
!!           global proc0 only)
!-----------------------------------------------!
  subroutine loop_io(series, b0, nkpolar, bins)
    use time, only: tt
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_io_total
    implicit none
    real(8), intent(in) :: series(nseries)
    real(8), intent(in) :: b0(3)
    integer, intent(in) :: nkpolar
    real(8), intent(in) :: bins(1:nkpolar, nbins)
    integer :: l

    if (proc0) call put_time_stamp(timer_io_total)

    ! output via NETCDF
    if(proc0) then
      status = nf90_put_var (ncid, tt_id, tt, start=(/nout/))

      do l = 1, nseries
        status = nf90_put_var (ncid, series_id(l), series(l), start=(/nout/))
      enddo

      status = nf90_put_var (ncid, b0_id, b0, start=(/1, nout/), count=(/3, 1/))

      do l = 1, nbins
        status = nf90_put_var (ncid, bin_id(l), bins(:,l), &
                               start=(/1, nout/), count=(/nkpolar, 1/))
      enddo

      status = nf90_sync (ncid)

      nout = nout + 1
    endif

    if (proc0) call put_time_stamp(timer_io_total)
  end subroutine loop_io


!-----------------------------------------------!
!> @author  YK
!! @brief   Append cross section via MPIIO
!!          cut_z0(:,:,l), cut_x0(:,:,l), cut_y0(:,:,l)
!!          are the z = 0, x = 0, y = 0 planes of the
!!          field l = i2d_*
!-----------------------------------------------!
  subroutine loop_io_2D(cut_z0, cut_x0, cut_y0)
    use grid, only: xx, yy, zz
    use grid, only: nlx, nly, nlz
    use grid, only: nlx_local
    use time, only: tt
    use mp, only: proc0, iproc_fft
    use mpiio, only: mpiio_write_var_2d
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_2D
    implicit none
    real(8), intent(in) :: cut_z0(1:nlx_local, 1:nly, n2d), &
                           cut_x0(1:nly      , 1:nlz, n2d), &
                           cut_y0(1:nlx_local, 1:nlz, n2d)

    integer, dimension(2) :: sizes, subsizes, starts
    integer :: l

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_2D)

    !--------------------------------------------------!
    !                    z = 0 cut
    !--------------------------------------------------!
    sizes = [nlx, nly]
    if(any(zz == 0.0)) then ! only the processes that have z = 0 write
      subsizes = [nlx_local, nly]
      starts   = [nlx_local*iproc_fft, 0]
    else
      subsizes = [1, 1]
      starts   = [0, 0]
    endif
    do l = 1, n2d
      call mpiio_write_var_2d(fh_2d(l, iz0), disp_2d(l, iz0), sizes, subsizes, starts, cut_z0(:,:,l))
    enddo

    !--------------------------------------------------!
    !                    x = 0 cut
    !--------------------------------------------------!
    sizes = [nly, nlz]
    if(any(xx == 0.0)) then ! only the processes that have x = 0 write
      subsizes = [nly, nlz]
      starts   = [0, 0]
    else
      subsizes = [1, 1]
      starts   = [0, 0]
    endif
    do l = 1, n2d
      call mpiio_write_var_2d(fh_2d(l, ix0), disp_2d(l, ix0), sizes, subsizes, starts, cut_x0(:,:,l))
    enddo

    !--------------------------------------------------!
    !                    y = 0 cut
    !--------------------------------------------------!
    sizes = [nlx, nlz]
    if(any(yy == 0.0)) then ! only the processes that have y = 0 write
      subsizes = [nlx_local, nlz]
      starts   = [nlx_local*iproc_fft, 0]
    else
      subsizes = [1, 1]
      starts   = [0, 0]
    endif
    do l = 1, n2d
      call mpiio_write_var_2d(fh_2d(l, iy0), disp_2d(l, iy0), sizes, subsizes, starts, cut_y0(:,:,l))
    enddo

    if(proc0) then
      write (unit=out2d_time_unit, fmt="(100es30.21)") tt
      flush (out2d_time_unit)
    endif

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_2D)
  end subroutine loop_io_2D


!-----------------------------------------------!
!> @author  YK
!! @brief   Append 3D field via MPIIO
!-----------------------------------------------!
  subroutine loop_io_3D
    use fields, only: bx, by, bz
    use mp, only: proc0, iproc_fft
    use time, only: tt
    use grid, only: nkx, nky, nky_local, nkz
    use mpiio, only: mpiio_write_var
    use shearing_box, only: tsc
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_3D
    implicit none
    integer, dimension(3) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_3D)

    sizes    = [nkz, nky      , nkx]
    subsizes = [nkz, nky_local, nkx]
    starts   = [0  , nky_local*iproc_fft, 0]

    !$acc update host(bx, by, bz)

    call mpiio_write_var(fh_3d(1), disp_3d(1), sizes, subsizes, starts, bx)
    call mpiio_write_var(fh_3d(2), disp_3d(2), sizes, subsizes, starts, by)
    call mpiio_write_var(fh_3d(3), disp_3d(3), sizes, subsizes, starts, bz)

    if(proc0) then
      write (unit=out3d_time_unit, fmt="(100es30.21)") tt, tsc
      flush (out3d_time_unit)
    endif

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_3D)
  end subroutine loop_io_3D


!-----------------------------------------------!
!> @author  YK
!! @brief   Save restart file via MPIIO
!!          (read back by fields::restart in the
!!           same order: bx, by, bz)
!-----------------------------------------------!
  subroutine save_restart
    use fields, only: bx, by, bz
    use mp, only: proc0, iproc_fft, comm_fft
    use grid, only: nkx, nky, nky_local, nkz
    use time, only: tt
    use params, only: restart_dir
    use file, only: open_output_file, close_file
    use mpiio, only: mpiio_write_one
    use shearing_box, only: tsc
    use time_stamp, only: put_time_stamp, timer_save_restart
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_save_restart)

    sizes    = [nkz, nky      , nkx]
    subsizes = [nkz, nky_local, nkx]
    starts   = [0  , nky_local*iproc_fft, 0]

    !$acc update host(bx, by, bz)

    call mpiio_write_one(bx, sizes, subsizes, starts, trim(restart_dir)//'bx.dat', comm_fft)
    call mpiio_write_one(by, sizes, subsizes, starts, trim(restart_dir)//'by.dat', comm_fft)
    call mpiio_write_one(bz, sizes, subsizes, starts, trim(restart_dir)//'bz.dat', comm_fft)

    if(proc0) then
      call open_output_file (time_unit, trim(restart_dir)//'time.dat')
      write (unit=time_unit, fmt="(3X, 'tt', 28X, 'tst')")
      write (unit=time_unit, fmt="(100es30.21)") tt, tsc
      call close_file (time_unit)
    endif

    if (proc0) call put_time_stamp(timer_save_restart)
  end subroutine save_restart


!-----------------------------------------------!
!> @author  YK
!! @brief   Finalization of NETCDF
!-----------------------------------------------!
  subroutine finish_io
    use mp, only: proc0
    use file, only: close_file
    implicit none
    integer :: ierr, l, p

    do l = 1, n3d
      call MPI_FILE_CLOSE(fh_3d(l), ierr)
    enddo
    do l = 1, n2d
      do p = 1, 3
        call MPI_FILE_CLOSE(fh_2d(l, p), ierr)
      enddo
    enddo

    if(proc0) then
      call close_file (out2d_time_unit)
      call close_file (out3d_time_unit)
    endif

    if(proc0) then
      status = nf90_close (ncid)
    endif

  end subroutine finish_io

end module io
