!-----------------------------------------------!
!> @author  YK
!! @brief   IO for HRMHD
!!
!!          The time series, the binned spectra and
!!          the 2D cuts are handled as lists: each
!!          entry has an index parameter (is_*, ib_*,
!!          i2d_*) and a name. The diagnostics fill
!!          the arrays through the index parameters,
!!          so a quantity cannot be written under the
!!          wrong name.
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
  integer, parameter, public :: is_upe2            =  1, is_bpe2            =  2, &
                                is_upa2            =  3, is_bpa2            =  4
  integer, parameter, public :: is_upe2dot         =  5, is_bpe2dot         =  6, &
                                is_upa2dot         =  7, is_bpa2dot         =  8
  integer, parameter, public :: is_upe2dissip      =  9, is_bpe2dissip      = 10, &
                                is_upa2dissip      = 11, is_bpa2dissip      = 12
  integer, parameter, public :: is_p_aw            = 13, is_p_compr         = 14, &
                                is_p_xhl           = 15
  integer, parameter, public :: is_zppe2           = 16, is_zmpe2           = 17, &
                                is_zppa2           = 18, is_zmpa2           = 19
  integer, parameter, public :: is_upe2_KAW_dissip = 20, is_bpe2_KAW_dissip = 21, &
                                is_upa2_KAW_dissip = 22, is_bpa2_KAW_dissip = 23
  integer, parameter, public :: is_upe2_ICW_dissip = 24, is_bpe2_ICW_dissip = 25, &
                                is_upa2_ICW_dissip = 26, is_bpa2_ICW_dissip = 27
  integer, parameter, public :: is_KAW_dissip_prp  = 28, is_KAW_dissip_par  = 29, &
                                is_ICW_dissip_prp  = 30, is_ICW_dissip_par  = 31
  integer, parameter, public :: is_xhl             = 32
  integer, parameter, public :: nseries = 32

  character(len=20), parameter :: series_names(nseries) = [character(len=20) :: &
    'upe2_sum'           , 'bpe2_sum'           , 'upa2_sum'           , 'bpa2_sum'           , &
    'upe2dot_sum'        , 'bpe2dot_sum'        , 'upa2dot_sum'        , 'bpa2dot_sum'        , &
    'upe2dissip_sum'     , 'bpe2dissip_sum'     , 'upa2dissip_sum'     , 'bpa2dissip_sum'     , &
    'p_aw_sum'           , 'p_compr_sum'        , 'p_xhl_sum'          ,                        &
    'zppe2_sum'          , 'zmpe2_sum'          , 'zppa2_sum'          , 'zmpa2_sum'          , &
    'upe2_KAW_dissip_sum', 'bpe2_KAW_dissip_sum', 'upa2_KAW_dissip_sum', 'bpa2_KAW_dissip_sum', &
    'upe2_ICW_dissip_sum', 'bpe2_ICW_dissip_sum', 'upa2_ICW_dissip_sum', 'bpa2_ICW_dissip_sum', &
    'KAW_dissip_prp_sum' , 'KAW_dissip_par_sum' , 'ICW_dissip_prp_sum' , 'ICW_dissip_par_sum' , &
    'xhl_sum'                                                                                     &
  ]

  !--------------------------------------------------!
  ! Polar spectra <runname>.out.nc: '<name>' (kpbin, kz, tt)
  !--------------------------------------------------!
  integer, parameter, public :: ib_upe2         =  1, ib_bpe2          =  2, &
                                ib_upa2         =  3, ib_bpa2          =  4
  integer, parameter, public :: ib_ux2          =  5, ib_uy2           =  6, &
                                ib_bx2          =  7, ib_by2           =  8
  integer, parameter, public :: ib_p_aw         =  9, ib_p_compr       = 10
  integer, parameter, public :: ib_dissip_aw    = 11, ib_dissip_compr  = 12, &
                                ib_dissip_KAW   = 13, ib_dissip_ICW    = 14
  integer, parameter, public :: ib_zppe2        = 15, ib_zmpe2         = 16, &
                                ib_zppa2        = 17, ib_zmpa2         = 18
  integer, parameter, public :: ib_upe2_KAW     = 19, ib_bpe2_KAW      = 20, &
                                ib_upa2_KAW     = 21, ib_bpa2_KAW      = 22
  integer, parameter, public :: ib_upe2_ICW     = 23, ib_bpe2_ICW      = 24, &
                                ib_upa2_ICW     = 25, ib_bpa2_ICW      = 26
  integer, parameter, public :: ib_flux_W       = 27, ib_flux_X        = 28
  integer, parameter, public :: nbins = 28

  character(len=20), parameter :: bin_names(nbins) = [character(len=20) :: &
    'upe2_bin'        , 'bpe2_bin'        , 'upa2_bin'        , 'bpa2_bin'        , &
    'ux2_bin'         , 'uy2_bin'         , 'bx2_bin'         , 'by2_bin'         , &
    'p_aw_bin'        , 'p_compr_bin'     ,                                         &
    'dissip_aw_bin'   , 'dissip_compr_bin', 'dissip_KAW_bin'  , 'dissip_ICW_bin'  , &
    'zppe2_bin'       , 'zmpe2_bin'       , 'zppa2_bin'       , 'zmpa2_bin'       , &
    'upe2_KAW_bin'    , 'bpe2_KAW_bin'    , 'upa2_KAW_bin'    , 'bpa2_KAW_bin'    , &
    'upe2_ICW_bin'    , 'bpe2_ICW_bin'    , 'upa2_ICW_bin'    , 'bpa2_ICW_bin'    , &
    'flux_W_bin'      , 'flux_X_bin'                                                &
  ]

  !--------------------------------------------------!
  ! 2D cuts out2d/<name>_r_{z0,x0,y0}.dat
  !--------------------------------------------------!
  integer, parameter, public :: i2d_phi = 1, i2d_psi = 2, i2d_upa = 3, i2d_bpa = 4, &
                                i2d_omg = 5, i2d_jpa = 6, i2d_ux  = 7, i2d_uy  = 8, &
                                i2d_bx  = 9, i2d_by  = 10
  integer, parameter, public :: n2d = 10

  character(len=3), parameter :: cut2d_names(n2d) = [character(len=3) :: &
    'phi', 'psi', 'upa', 'bpa', 'omg', 'jpa', 'ux', 'uy', 'bx', 'by']

  ! cut planes, in the order of the last index of the fh_2d/disp_2d arrays
  integer, parameter :: iz0 = 1, ix0 = 2, iy0 = 3
  character(len=2), parameter :: plane_names(3) = ['z0', 'x0', 'y0']

  !--------------------------------------------------!
  ! 3D fields out3d/<name>.dat and restart files
  !--------------------------------------------------!
  integer, parameter :: n3d = 4
  character(len=3), parameter :: out3d_names(n3d) = ['phi', 'psi', 'upa', 'bpa']

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
  integer, parameter :: nparams = 18
  integer :: param_id(nparams)
  ! coordinate
  integer :: xx_id, yy_id, zz_id, kx_id, ky_id, kz_id, kpbin_id, tt_id
  ! time series and polar spectra
  integer :: series_id(nseries), bin_id(nbins)

  integer (kind_nf) :: xx_dim, yy_dim, zz_dim, kx_dim, ky_dim, kz_dim, kpbin_dim, tt_dim
  integer, dimension (3) :: bin_dim

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
    use mp, only: proc0, proc0_m, proc0_s, comm_fft
    use file, only: open_output_file
    implicit none
    integer :: l

    ! Only comm_fft group 0 performs MPI-IO for the redundantly-solved fields;
    ! every group holds an identical copy, so opening/writing on comm_fft from
    ! group 0 alone avoids cross-group file corruption.
    if (proc0_m .and. proc0_s) then
      do l = 1, n3d
        call set_file_handle('out3d/'//trim(out3d_names(l))//'.dat', fh_3d(l), disp_3d(l), comm_fft)
      enddo
    endif

    if(proc0) then
      call open_output_file (out3d_time_unit, 'out3d/time.dat')
    endif
  end subroutine init_io_decomp_3d


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of MPIIO for 2D
!-----------------------------------------------!
  subroutine init_io_decomp_2d
    use mp, only: proc0, proc0_m, proc0_s, comm_fft
    use file, only: open_output_file
    implicit none
    integer :: l, p

    ! Only comm_fft group 0 performs the 2D-cut MPI-IO (see init_io_decomp_3d).
    if (proc0_m .and. proc0_s) then
      do l = 1, n2d
        do p = 1, 3
          call set_file_handle('out2d/'//trim(cut2d_names(l))//'_r_'//plane_names(p)//'.dat', &
                               fh_2d(l, p), disp_2d(l, p), comm_fft)
        enddo
      enddo
    endif

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
    use params, only: runname, sgm, rho, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp, &
                      nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                      etapa_x, etapa_x_exp, etapa_z, etapa_z_exp
    implicit none
    integer, intent(in) :: nkpolar
    real(8), intent(in) :: kpbin(1:nkpolar)
    character(len=12), parameter :: param_names(nparams) = [character(len=12) :: &
      'sgm'    , 'rho'        ,                                         &
      'nupe_x' , 'nupe_x_exp' , 'nupe_z' , 'nupe_z_exp' ,               &
      'nupa_x' , 'nupa_x_exp' , 'nupa_z' , 'nupa_z_exp' ,               &
      'etape_x', 'etape_x_exp', 'etape_z', 'etape_z_exp',               &
      'etapa_x', 'etapa_x_exp', 'etapa_z', 'etapa_z_exp'                &
    ]
    real(8) :: param_values(nparams)
    integer :: l

    if(proc0) then
      ! exponents are written as doubles, as in Calliope
      param_values = [ sgm, rho, &
                       nupe_x , dble(nupe_x_exp) , nupe_z , dble(nupe_z_exp) , &
                       nupa_x , dble(nupa_x_exp) , nupa_z , dble(nupa_z_exp) , &
                       etape_x, dble(etape_x_exp), etape_z, dble(etape_z_exp), &
                       etapa_x, dble(etapa_x_exp), etapa_z, dble(etapa_z_exp) ]

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

      ! polar spectrum
      bin_dim (1) = kpbin_dim
      bin_dim (2) = kz_dim
      bin_dim (3) = tt_dim
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
!!          (series and bins are filled on the
!!           global proc0 only)
!-----------------------------------------------!
  subroutine loop_io(series, nkpolar, bins)
    use time, only: tt
    use grid, only: nkz
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_io_total
    implicit none
    real(8), intent(in) :: series(nseries)
    integer, intent(in) :: nkpolar
    real(8), intent(in) :: bins(1:nkpolar, nkz, nbins)
    integer :: l

    if (proc0) call put_time_stamp(timer_io_total)

    ! output via NETCDF
    if(proc0) then
      status = nf90_put_var (ncid, tt_id, tt, start=(/nout/))

      do l = 1, nseries
        status = nf90_put_var (ncid, series_id(l), series(l), start=(/nout/))
      enddo

      do l = 1, nbins
        status = nf90_put_var (ncid, bin_id(l), bins(:,:,l), &
                               start=(/1, 1, nout/), count=(/nkpolar, nkz, 1/))
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
    use mp, only: proc0, iproc_fft, proc0_m, proc0_s
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

    ! Only comm_fft group 0 opened the 2D-cut files, so only it writes the cuts
    ! (collective over comm_fft). Other groups hold identical redundant copies.
    if (proc0_m .and. proc0_s) then

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

    endif ! proc0_m .and. proc0_s

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
    use fields, only: phi, psi, upa, bpa
    use mp, only: proc0, iproc_fft, proc0_m, proc0_s
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

    !$acc update host(phi, psi, upa, bpa)

    ! Only comm_fft group 0 opened the 3D files (see init_io_decomp_3d).
    if (proc0_m .and. proc0_s) then
      call mpiio_write_var(fh_3d(1), disp_3d(1), sizes, subsizes, starts, phi)
      call mpiio_write_var(fh_3d(2), disp_3d(2), sizes, subsizes, starts, psi)
      call mpiio_write_var(fh_3d(3), disp_3d(3), sizes, subsizes, starts, upa)
      call mpiio_write_var(fh_3d(4), disp_3d(4), sizes, subsizes, starts, bpa)
    endif

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
!!          (read back by fields::init_restart in
!!           the same order: phi, omg, psi, upa, bpa)
!-----------------------------------------------!
  subroutine save_restart
    use fields, only: phi, omg, psi, upa, bpa
    use mp, only: proc0, iproc_fft, proc0_m, proc0_s, comm_fft
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

    !$acc update host(phi, omg, psi, upa, bpa)

    ! Only comm_fft group 0 writes the checkpoint; every group holds an identical
    ! copy, so writing from all groups would corrupt the same files.
    if (proc0_m .and. proc0_s) then
      call mpiio_write_one(phi, sizes, subsizes, starts, trim(restart_dir)//'phi.dat', comm_fft)
      call mpiio_write_one(omg, sizes, subsizes, starts, trim(restart_dir)//'omg.dat', comm_fft)
      call mpiio_write_one(psi, sizes, subsizes, starts, trim(restart_dir)//'psi.dat', comm_fft)
      call mpiio_write_one(upa, sizes, subsizes, starts, trim(restart_dir)//'upa.dat', comm_fft)
      call mpiio_write_one(bpa, sizes, subsizes, starts, trim(restart_dir)//'bpa.dat', comm_fft)
    endif

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
    use mp, only: proc0, proc0_m, proc0_s
    use file, only: close_file
    implicit none
    integer :: ierr, l, p

    ! Only comm_fft group 0 opened the MPI-IO handles (see init_io_decomp_*).
    if (proc0_m .and. proc0_s) then
      do l = 1, n3d
        call MPI_FILE_CLOSE(fh_3d(l), ierr)
      enddo
      do l = 1, n2d
        do p = 1, 3
          call MPI_FILE_CLOSE(fh_2d(l, p), ierr)
        enddo
      enddo
    endif

    if(proc0) then
      call close_file (out2d_time_unit)
      call close_file (out3d_time_unit)
    endif

    if(proc0) then
      status = nf90_close (ncid)
    endif

  end subroutine finish_io

end module io
