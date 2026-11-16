!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   IO for MHD_COMP_ISOTH
!!
!!          Gallope port: spectral data are (kz, ky, kx)
!!          with kz the r2c half axis and ky split over
!!          comm_fft; real data are (z, y, x) with x
!!          split over comm_fft. The model runs with
!!          P_m = P_s = 1 (model_specific), so comm_fft
!!          spans all ranks and proc0 is iproc_fft = 0.
!-----------------------------------------------!
module io
  use netcdf
  use MPI
  implicit none

  public :: init_io, finish_io, loop_io, loop_io_2D, loop_io_3D, save_restart

  private

  ! MPIIO for 3D
  integer :: fh_rho, fh_mx, fh_my, fh_mz, fh_bx, fh_by, fh_bz
  character(len=100) :: filename
  integer (kind=MPI_OFFSET_KIND) :: disp_rho
  integer (kind=MPI_OFFSET_KIND) :: disp_mx
  integer (kind=MPI_OFFSET_KIND) :: disp_my
  integer (kind=MPI_OFFSET_KIND) :: disp_mz
  integer (kind=MPI_OFFSET_KIND) :: disp_bx
  integer (kind=MPI_OFFSET_KIND) :: disp_by
  integer (kind=MPI_OFFSET_KIND) :: disp_bz
  integer :: out3d_time_unit

  ! MPIIO for 2D
  integer :: fh_rho_r_z0, fh_rho_r_x0, fh_rho_r_y0
  integer :: fh_mx_r_z0 , fh_mx_r_x0 , fh_mx_r_y0, &
             fh_my_r_z0 , fh_my_r_x0 , fh_my_r_y0, &
             fh_mz_r_z0 , fh_mz_r_x0 , fh_mz_r_y0
  integer :: fh_wx_r_z0 , fh_wx_r_x0 , fh_wx_r_y0, &
             fh_wy_r_z0 , fh_wy_r_x0 , fh_wy_r_y0, &
             fh_wz_r_z0 , fh_wz_r_x0 , fh_wz_r_y0
  integer :: fh_bx_r_z0 , fh_bx_r_x0 , fh_bx_r_y0, &
             fh_by_r_z0 , fh_by_r_x0 , fh_by_r_y0, &
             fh_bz_r_z0 , fh_bz_r_x0 , fh_bz_r_y0
  integer :: fh_jx_r_z0 , fh_jx_r_x0 , fh_jx_r_y0, &
             fh_jy_r_z0 , fh_jy_r_x0 , fh_jy_r_y0, &
             fh_jz_r_z0 , fh_jz_r_x0 , fh_jz_r_y0
  integer :: fh_rho_kxy, fh_rho_kyz, fh_rho_kxz, &
             fh_u2_kxy , fh_u2_kyz , fh_u2_kxz, &
             fh_b2_kxy , fh_b2_kyz , fh_b2_kxz

  integer (kind=MPI_OFFSET_KIND) :: disp_rho_r_z0, disp_rho_r_x0, disp_rho_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_mx_r_z0 , disp_mx_r_x0 , disp_mx_r_y0 , &
                                    disp_my_r_z0 , disp_my_r_x0 , disp_my_r_y0 , &
                                    disp_mz_r_z0 , disp_mz_r_x0 , disp_mz_r_y0 
  integer (kind=MPI_OFFSET_KIND) :: disp_wx_r_z0 , disp_wx_r_x0 , disp_wx_r_y0 , &
                                    disp_wy_r_z0 , disp_wy_r_x0 , disp_wy_r_y0 , &
                                    disp_wz_r_z0 , disp_wz_r_x0 , disp_wz_r_y0 
  integer (kind=MPI_OFFSET_KIND) :: disp_bx_r_z0 , disp_bx_r_x0 , disp_bx_r_y0 , &
                                    disp_by_r_z0 , disp_by_r_x0 , disp_by_r_y0 , &
                                    disp_bz_r_z0 , disp_bz_r_x0 , disp_bz_r_y0 
  integer (kind=MPI_OFFSET_KIND) :: disp_jx_r_z0 , disp_jx_r_x0 , disp_jx_r_y0 , &
                                    disp_jy_r_z0 , disp_jy_r_x0 , disp_jy_r_y0 , &
                                    disp_jz_r_z0 , disp_jz_r_x0 , disp_jz_r_y0 
  integer (kind=MPI_OFFSET_KIND) :: disp_rho_kxy, disp_rho_kyz, disp_rho_kxz, &
                                    disp_u2_kxy , disp_u2_kyz , disp_u2_kxz , &
                                    disp_b2_kxy , disp_b2_kyz , disp_b2_kxz 
  integer :: out2d_time_unit

  ! NETCDF for regular output file
  integer :: status
  integer, parameter :: kind_nf = kind (NF90_NOERR)
  integer (kind_nf) :: ncid
  integer :: run_id
  integer (kind_nf) :: char10_dim
  ! parameter
  integer :: beta0_id
  integer :: shear_flg_id, q_id
  integer :: nu_id , nu_h_id , nu_h_exp_id
  integer :: eta_id, eta_h_id, eta_h_exp_id
  integer :: lmd_id, lmd_h_id, lmd_h_exp_id
  ! coordinate
  integer :: xx_id, yy_id, zz_id, kx_id, ky_id, kz_id, kpbin_id, tt_id
  ! total energy and max history
  integer :: wkin_sum_id, wmag_sum_id, wrho_sum_id
  integer :: wkin_dot_sum_id, wmag_dot_sum_id, wrho_dot_sum_id
  integer :: wkin_dissip_sum_id, wmag_dissip_sum_id, wrho_dissip_sum_id
  integer :: p_ext_sum_id, p_re_sum_id, p_ma_sum_id
  integer :: smach_rms_id, amach_rms_id, beta_rms_id
  integer :: zp2_sum_id, zm2_sum_id
  ! mean magnetic field
  integer :: b0_id
  ! polar spectrum
  integer :: rho_bin_id
  integer :: u2_bin_id, ux2_bin_id, uy2_bin_id, uz2_bin_id
  integer :: b2_bin_id, bx2_bin_id, by2_bin_id, bz2_bin_id
  integer :: zp2_bin_id, zm2_bin_id

  integer (kind_nf) :: xx_dim, yy_dim, zz_dim, kx_dim, ky_dim, kz_dim, kpbin_dim, tt_dim, vec_dim
  integer, dimension (2) :: mean_fld_dim
  integer, dimension (2) :: bin_dim

  integer :: nout

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
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
!! @date    29 Dec 2018
!! @brief   Initialization of MPIIO for 3D
!-----------------------------------------------!
  subroutine init_io_decomp_3d
    use mp, only: proc0, comm_fft
    use file, only: open_output_file
    implicit none

    call set_file_handle('out3d/rho.dat', fh_rho, disp_rho, comm_fft)
    call set_file_handle('out3d/mx.dat' , fh_mx , disp_mx , comm_fft)
    call set_file_handle('out3d/my.dat' , fh_my , disp_my , comm_fft)
    call set_file_handle('out3d/mz.dat' , fh_mz , disp_mz , comm_fft)
    call set_file_handle('out3d/bx.dat' , fh_bx , disp_bx , comm_fft)
    call set_file_handle('out3d/by.dat' , fh_by , disp_by , comm_fft)
    call set_file_handle('out3d/bz.dat' , fh_bz , disp_bz , comm_fft)

    if(proc0) then
      call open_output_file (out3d_time_unit, 'out3d/time.dat')
    endif
  end subroutine init_io_decomp_3d


!-----------------------------------------------!
!> @author  YK
!! @date    7 May 2022
!! @brief   Initialization of MPIIO for 2D
!-----------------------------------------------!
  subroutine init_io_decomp_2d
    use mp, only: proc0, comm_fft
    use file, only: open_output_file
    implicit none

    !--------------------------------------------------!
    !                       rho
    !--------------------------------------------------!
    call set_file_handle('out2d/rho_r_z0.dat', fh_rho_r_z0, disp_rho_r_z0, comm_fft)
    call set_file_handle('out2d/rho_r_x0.dat', fh_rho_r_x0, disp_rho_r_x0, comm_fft)
    call set_file_handle('out2d/rho_r_y0.dat', fh_rho_r_y0, disp_rho_r_y0, comm_fft)

    !--------------------------------------------------!
    !                        m
    !--------------------------------------------------!
    call set_file_handle('out2d/mx_r_z0.dat', fh_mx_r_z0, disp_mx_r_z0, comm_fft)
    call set_file_handle('out2d/mx_r_x0.dat', fh_mx_r_x0, disp_mx_r_x0, comm_fft)
    call set_file_handle('out2d/mx_r_y0.dat', fh_mx_r_y0, disp_mx_r_y0, comm_fft)

    call set_file_handle('out2d/my_r_z0.dat', fh_my_r_z0, disp_my_r_z0, comm_fft)
    call set_file_handle('out2d/my_r_x0.dat', fh_my_r_x0, disp_my_r_x0, comm_fft)
    call set_file_handle('out2d/my_r_y0.dat', fh_my_r_y0, disp_my_r_y0, comm_fft)

    call set_file_handle('out2d/mz_r_z0.dat', fh_mz_r_z0, disp_mz_r_z0, comm_fft)
    call set_file_handle('out2d/mz_r_x0.dat', fh_mz_r_x0, disp_mz_r_x0, comm_fft)
    call set_file_handle('out2d/mz_r_y0.dat', fh_mz_r_y0, disp_mz_r_y0, comm_fft)

    !--------------------------------------------------!
    !                        w
    !--------------------------------------------------!
    call set_file_handle('out2d/wx_r_z0.dat', fh_wx_r_z0, disp_wx_r_z0, comm_fft)
    call set_file_handle('out2d/wx_r_x0.dat', fh_wx_r_x0, disp_wx_r_x0, comm_fft)
    call set_file_handle('out2d/wx_r_y0.dat', fh_wx_r_y0, disp_wx_r_y0, comm_fft)

    call set_file_handle('out2d/wy_r_z0.dat', fh_wy_r_z0, disp_wy_r_z0, comm_fft)
    call set_file_handle('out2d/wy_r_x0.dat', fh_wy_r_x0, disp_wy_r_x0, comm_fft)
    call set_file_handle('out2d/wy_r_y0.dat', fh_wy_r_y0, disp_wy_r_y0, comm_fft)

    call set_file_handle('out2d/wz_r_z0.dat', fh_wz_r_z0, disp_wz_r_z0, comm_fft)
    call set_file_handle('out2d/wz_r_x0.dat', fh_wz_r_x0, disp_wz_r_x0, comm_fft)
    call set_file_handle('out2d/wz_r_y0.dat', fh_wz_r_y0, disp_wz_r_y0, comm_fft)

    !--------------------------------------------------!
    !                        b
    !--------------------------------------------------!
    call set_file_handle('out2d/bx_r_z0.dat', fh_bx_r_z0, disp_bx_r_z0, comm_fft)
    call set_file_handle('out2d/bx_r_x0.dat', fh_bx_r_x0, disp_bx_r_x0, comm_fft)
    call set_file_handle('out2d/bx_r_y0.dat', fh_bx_r_y0, disp_bx_r_y0, comm_fft)

    call set_file_handle('out2d/by_r_z0.dat', fh_by_r_z0, disp_by_r_z0, comm_fft)
    call set_file_handle('out2d/by_r_x0.dat', fh_by_r_x0, disp_by_r_x0, comm_fft)
    call set_file_handle('out2d/by_r_y0.dat', fh_by_r_y0, disp_by_r_y0, comm_fft)

    call set_file_handle('out2d/bz_r_z0.dat', fh_bz_r_z0, disp_bz_r_z0, comm_fft)
    call set_file_handle('out2d/bz_r_x0.dat', fh_bz_r_x0, disp_bz_r_x0, comm_fft)
    call set_file_handle('out2d/bz_r_y0.dat', fh_bz_r_y0, disp_bz_r_y0, comm_fft)

    !--------------------------------------------------!
    !                        j
    !--------------------------------------------------!
    call set_file_handle('out2d/jx_r_z0.dat', fh_jx_r_z0, disp_jx_r_z0, comm_fft)
    call set_file_handle('out2d/jx_r_x0.dat', fh_jx_r_x0, disp_jx_r_x0, comm_fft)
    call set_file_handle('out2d/jx_r_y0.dat', fh_jx_r_y0, disp_jx_r_y0, comm_fft)

    call set_file_handle('out2d/jy_r_z0.dat', fh_jy_r_z0, disp_jy_r_z0, comm_fft)
    call set_file_handle('out2d/jy_r_x0.dat', fh_jy_r_x0, disp_jy_r_x0, comm_fft)
    call set_file_handle('out2d/jy_r_y0.dat', fh_jy_r_y0, disp_jy_r_y0, comm_fft)

    call set_file_handle('out2d/jz_r_z0.dat', fh_jz_r_z0, disp_jz_r_z0, comm_fft)
    call set_file_handle('out2d/jz_r_x0.dat', fh_jz_r_x0, disp_jz_r_x0, comm_fft)
    call set_file_handle('out2d/jz_r_y0.dat', fh_jz_r_y0, disp_jz_r_y0, comm_fft)

    !--------------------------------------------------!
    !                      rho_k
    !--------------------------------------------------!
    call set_file_handle('out2d/rho_kxy_sum_kz.dat', fh_rho_kxy, disp_rho_kxy, comm_fft)
    call set_file_handle('out2d/rho_kyz_sum_kx.dat', fh_rho_kyz, disp_rho_kyz, comm_fft)
    call set_file_handle('out2d/rho_kxz_sum_ky.dat', fh_rho_kxz, disp_rho_kxz, comm_fft)

    !--------------------------------------------------!
    !                      u^2_k
    !--------------------------------------------------!
    call set_file_handle('out2d/u2_kxy_sum_kz.dat', fh_u2_kxy, disp_u2_kxy, comm_fft)
    call set_file_handle('out2d/u2_kyz_sum_kx.dat', fh_u2_kyz, disp_u2_kyz, comm_fft)
    call set_file_handle('out2d/u2_kxz_sum_ky.dat', fh_u2_kxz, disp_u2_kxz, comm_fft)

    !--------------------------------------------------!
    !                      b^2_k
    !--------------------------------------------------!
    call set_file_handle('out2d/b2_kxy_sum_kz.dat', fh_b2_kxy, disp_b2_kxy, comm_fft)
    call set_file_handle('out2d/b2_kyz_sum_kx.dat', fh_b2_kyz, disp_b2_kyz, comm_fft)
    call set_file_handle('out2d/b2_kxz_sum_ky.dat', fh_b2_kxz, disp_b2_kxz, comm_fft)

    if(proc0) then
      call open_output_file (out2d_time_unit, 'out2d/time.dat')
    endif
  end subroutine init_io_decomp_2d


!-----------------------------------------------!
!> @author  YK
!! @date    20 Jun 2022
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
!! @date    29 Dec 2018
!! @brief   Initialization of NETCDF
!-----------------------------------------------!
  subroutine init_io_netcdf(nkpolar, kpbin)
    use grid, only: nlx, nly, nlz 
    use grid, only: xx_global, yy, zz, kx, ky_global, kz
    use mp, only: proc0
    use params, only: runname, nu, nu_h, nu_h_exp, eta, eta_h, eta_h_exp, lmd, lmd_h, lmd_h_exp, beta0, q
    use shearing_box, only: shear_flg
    implicit none
    integer , intent(in) :: nkpolar
    real(8), intent(in) :: kpbin(1:nkpolar)

    if(proc0) then
      !--------------------------------------------------!
      ! Output for parameters, time history, and spectra
      !--------------------------------------------------!
      filename = trim(runname)//'.out.nc' ! File name
      status = nf90_create (filename, NF90_CLOBBER, ncid)

      status = nf90_put_att (ncid, NF90_GLOBAL, 'title', 'calliope simulation data')
      status = nf90_def_dim (ncid, 'char10', 10, char10_dim)
      status = nf90_def_var (ncid, 'run_info', NF90_CHAR, char10_dim, run_id)
      status = nf90_put_att (ncid, run_id, 'model', _MODEL_)

      status = nf90_def_dim (ncid, 'xx', size(xx_global), xx_dim)
      status = nf90_def_dim (ncid, 'yy', size(yy), yy_dim)
      status = nf90_def_dim (ncid, 'zz', size(zz), zz_dim)
      status = nf90_def_dim (ncid, 'kx', size(kx), kx_dim)
      status = nf90_def_dim (ncid, 'ky', size(ky_global), ky_dim)
      status = nf90_def_dim (ncid, 'kz', size(kz), kz_dim)
      status = nf90_def_dim (ncid, 'kpbin', size(kpbin), kpbin_dim)
      status = nf90_def_dim (ncid, 'tt', NF90_UNLIMITED, tt_dim)
      status = nf90_def_dim (ncid, '3', 3, vec_dim)

      status = nf90_def_var (ncid, 'beta0', NF90_DOUBLE, beta0_id)
      status = nf90_def_var (ncid, 'shear_flg', NF90_INT, shear_flg_id)
      status = nf90_def_var (ncid, 'q', NF90_DOUBLE, q_id)
      status = nf90_def_var (ncid, 'nu'   , NF90_DOUBLE, nu_id   )
      status = nf90_def_var (ncid, 'nu_h' , NF90_DOUBLE, nu_h_id )
      status = nf90_def_var (ncid, 'eta'  , NF90_DOUBLE, eta_id  )
      status = nf90_def_var (ncid, 'eta_h', NF90_DOUBLE, eta_h_id)
      status = nf90_def_var (ncid, 'lmd'  , NF90_DOUBLE, lmd_id  )
      status = nf90_def_var (ncid, 'lmd_h', NF90_DOUBLE, lmd_h_id)
      status = nf90_def_var (ncid, 'nu_h_exp' , NF90_DOUBLE, nu_h_exp_id )
      status = nf90_def_var (ncid, 'eta_h_exp', NF90_DOUBLE, eta_h_exp_id)
      status = nf90_def_var (ncid, 'lmd_h_exp', NF90_DOUBLE, lmd_h_exp_id)

      status = nf90_def_var (ncid, 'xx', NF90_DOUBLE, xx_dim, xx_id)
      status = nf90_def_var (ncid, 'yy', NF90_DOUBLE, yy_dim, yy_id)
      status = nf90_def_var (ncid, 'zz', NF90_DOUBLE, zz_dim, zz_id)
      status = nf90_def_var (ncid, 'kx', NF90_DOUBLE, kx_dim, kx_id)
      status = nf90_def_var (ncid, 'ky', NF90_DOUBLE, ky_dim, ky_id)
      status = nf90_def_var (ncid, 'kz', NF90_DOUBLE, kz_dim, kz_id)
      status = nf90_def_var (ncid, 'kpbin', NF90_DOUBLE, kpbin_dim, kpbin_id)
      status = nf90_def_var (ncid, 'tt', NF90_DOUBLE, tt_dim, tt_id)
      ! energy and max values history
      status = nf90_def_var (ncid, 'wkin_sum', NF90_DOUBLE, tt_dim, wkin_sum_id)
      status = nf90_def_var (ncid, 'wmag_sum', NF90_DOUBLE, tt_dim, wmag_sum_id)
      status = nf90_def_var (ncid, 'wrho_sum', NF90_DOUBLE, tt_dim, wrho_sum_id)
      status = nf90_def_var (ncid, 'wkin_dot_sum', NF90_DOUBLE, tt_dim, wkin_dot_sum_id)
      status = nf90_def_var (ncid, 'wmag_dot_sum', NF90_DOUBLE, tt_dim, wmag_dot_sum_id)
      status = nf90_def_var (ncid, 'wrho_dot_sum', NF90_DOUBLE, tt_dim, wrho_dot_sum_id)
      status = nf90_def_var (ncid, 'wkin_dissip_sum', NF90_DOUBLE, tt_dim, wkin_dissip_sum_id)
      status = nf90_def_var (ncid, 'wmag_dissip_sum', NF90_DOUBLE, tt_dim, wmag_dissip_sum_id)
      status = nf90_def_var (ncid, 'wrho_dissip_sum', NF90_DOUBLE, tt_dim, wrho_dissip_sum_id)
      status = nf90_def_var (ncid, 'p_ext_sum'   , NF90_DOUBLE, tt_dim, p_ext_sum_id  )
      status = nf90_def_var (ncid, 'p_re_sum'    , NF90_DOUBLE, tt_dim, p_re_sum_id   )
      status = nf90_def_var (ncid, 'p_ma_sum'    , NF90_DOUBLE, tt_dim, p_ma_sum_id   )
      status = nf90_def_var (ncid, 'smach_rms'   , NF90_DOUBLE, tt_dim, smach_rms_id )
      status = nf90_def_var (ncid, 'amach_rms'   , NF90_DOUBLE, tt_dim, amach_rms_id )
      status = nf90_def_var (ncid, 'beta_rms'    , NF90_DOUBLE, tt_dim, beta_rms_id  )
      status = nf90_def_var (ncid, 'zp2_sum', NF90_DOUBLE, tt_dim, zp2_sum_id)
      status = nf90_def_var (ncid, 'zm2_sum', NF90_DOUBLE, tt_dim, zm2_sum_id)
      ! polar spectrum
      bin_dim (1) = kpbin_dim
      bin_dim (2) = tt_dim
      status = nf90_def_var (ncid, 'rho_bin', NF90_DOUBLE, bin_dim, rho_bin_id)
      status = nf90_def_var (ncid,  'u2_bin', NF90_DOUBLE, bin_dim,  u2_bin_id)
      status = nf90_def_var (ncid, 'ux2_bin', NF90_DOUBLE, bin_dim, ux2_bin_id)
      status = nf90_def_var (ncid, 'uy2_bin', NF90_DOUBLE, bin_dim, uy2_bin_id)
      status = nf90_def_var (ncid, 'uz2_bin', NF90_DOUBLE, bin_dim, uz2_bin_id)
      status = nf90_def_var (ncid,  'b2_bin', NF90_DOUBLE, bin_dim,  b2_bin_id)
      status = nf90_def_var (ncid, 'bx2_bin', NF90_DOUBLE, bin_dim, bx2_bin_id)
      status = nf90_def_var (ncid, 'by2_bin', NF90_DOUBLE, bin_dim, by2_bin_id)
      status = nf90_def_var (ncid, 'bz2_bin', NF90_DOUBLE, bin_dim, bz2_bin_id)
      status = nf90_def_var (ncid, 'zp2_bin', NF90_DOUBLE, bin_dim, zp2_bin_id)
      status = nf90_def_var (ncid, 'zm2_bin', NF90_DOUBLE, bin_dim, zm2_bin_id)
      ! mean magnetic field
      mean_fld_dim (1) = vec_dim
      mean_fld_dim (2) = tt_dim
      status = nf90_def_var (ncid, 'b0', NF90_DOUBLE, mean_fld_dim, b0_id)

      status = nf90_enddef (ncid)  ! out of definition mode

      status = nf90_put_var (ncid, beta0_id, beta0)
      status = nf90_put_var (ncid, shear_flg_id, shear_flg)
      status = nf90_put_var (ncid, q_id, q)
      status = nf90_put_var (ncid, nu_id   , nu   )
      status = nf90_put_var (ncid, nu_h_id , nu_h )
      status = nf90_put_var (ncid, eta_id  , eta  )
      status = nf90_put_var (ncid, eta_h_id, eta_h)
      status = nf90_put_var (ncid, lmd_id  , lmd  )
      status = nf90_put_var (ncid, lmd_h_id, lmd_h)
      status = nf90_put_var (ncid, nu_h_exp_id , dble(nu_h_exp ))
      status = nf90_put_var (ncid, eta_h_exp_id, dble(eta_h_exp))
      status = nf90_put_var (ncid, lmd_h_exp_id, dble(lmd_h_exp))

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
!! @date    29 Dec 2018
!! @brief   Append variables to NETCDF
!           for params, time history & spectra
!-----------------------------------------------!
  subroutine loop_io( &
                      wkin_sum, wmag_sum, wrho_sum, &
                      wkin_dot_sum, wmag_dot_sum, wrho_dot_sum, &
                      wkin_dissip_sum, wmag_dissip_sum, wrho_dissip_sum, &
                      p_ext_sum, p_re_sum, p_ma_sum, &
                      smach_rms, amach_rms, beta_rms, &
                      zp2_sum, zm2_sum, &
                      bx0, by0, bz0, &
                      !
                      nkpolar, &
                      rho_bin, &
                      u2_bin, ux2_bin, uy2_bin, uz2_bin, &
                      b2_bin, bx2_bin, by2_bin, bz2_bin, &
                      zp2_bin, zm2_bin &
                    )
    use time, only: tt
    use grid, only: nlx, nly, nlz
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_io_total
    implicit none
    real(8), intent(in) :: wkin_sum, wmag_sum, wrho_sum
    real(8), intent(in) :: wkin_dot_sum, wmag_dot_sum, wrho_dot_sum
    real(8), intent(in) :: wkin_dissip_sum, wmag_dissip_sum, wrho_dissip_sum
    real(8), intent(in) :: p_ext_sum, p_re_sum, p_ma_sum
    real(8), intent(in) :: zp2_sum, zm2_sum
    real(8), intent(in) :: smach_rms, amach_rms, beta_rms
    real(8), intent(in) :: bx0, by0, bz0

    integer , intent(in) :: nkpolar
    real(8), intent(in) :: rho_bin(1:nkpolar)
    real(8), intent(in) :: u2_bin (1:nkpolar), ux2_bin(1:nkpolar), uy2_bin(1:nkpolar), uz2_bin(1:nkpolar)
    real(8), intent(in) :: b2_bin (1:nkpolar), bx2_bin(1:nkpolar), by2_bin(1:nkpolar), bz2_bin(1:nkpolar)
    real(8), intent(in) :: zp2_bin(1:nkpolar), zm2_bin(1:nkpolar)

    integer, dimension (2) :: start2, count2

    if (proc0) call put_time_stamp(timer_io_total)

    ! output via NETCDF
    if(proc0) then
      ! energy and max values history
      status = nf90_put_var (ncid, tt_id, tt, start=(/nout/))
      status = nf90_put_var (ncid, wkin_sum_id, wkin_sum, start=(/nout/))
      status = nf90_put_var (ncid, wmag_sum_id, wmag_sum, start=(/nout/))
      status = nf90_put_var (ncid, wrho_sum_id, wrho_sum, start=(/nout/))
      status = nf90_put_var (ncid, wkin_dot_sum_id, wkin_dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, wmag_dot_sum_id, wmag_dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, wrho_dot_sum_id, wrho_dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, wkin_dissip_sum_id, wkin_dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, wmag_dissip_sum_id, wmag_dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, wrho_dissip_sum_id, wrho_dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_ext_sum_id, p_ext_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_re_sum_id , p_re_sum , start=(/nout/))
      status = nf90_put_var (ncid, p_ma_sum_id , p_ma_sum , start=(/nout/))
      status = nf90_put_var (ncid, zp2_sum_id, zp2_sum, start=(/nout/))
      status = nf90_put_var (ncid, zm2_sum_id, zm2_sum, start=(/nout/))
      status = nf90_put_var (ncid, smach_rms_id, smach_rms, start=(/nout/))
      status = nf90_put_var (ncid, amach_rms_id, amach_rms, start=(/nout/))
      status = nf90_put_var (ncid, beta_rms_id , beta_rms , start=(/nout/))
      ! mean magnetic field
      start2(1) = 1
      start2(2) = nout

      count2(1) = 3
      count2(2) = 1
      status = nf90_put_var (ncid, b0_id, (/bx0, by0, bz0/), start=start2, count=count2)
      ! polar spectrum
      start2(1) = 1
      start2(2) = nout

      count2(1) = nkpolar
      count2(2) = 1
      status = nf90_put_var (ncid, rho_bin_id, rho_bin, start=start2, count=count2)
      status = nf90_put_var (ncid,  u2_bin_id,  u2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, ux2_bin_id, ux2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, uy2_bin_id, uy2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, uz2_bin_id, uz2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid,  b2_bin_id,  b2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, bx2_bin_id, bx2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, by2_bin_id, by2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, bz2_bin_id, bz2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, zp2_bin_id, zp2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, zm2_bin_id, zm2_bin, start=start2, count=count2)

      status = nf90_sync (ncid)

      nout = nout + 1
    endif

    if (proc0) call put_time_stamp(timer_io_total)
  end subroutine loop_io


!-----------------------------------------------!
!> @author  YK
!! @date    28 Jun 2021
!! @brief   Append cross section via MPIIO
!-----------------------------------------------!
  subroutine loop_io_2D( &
                      rho_r_z0, rho_r_x0, rho_r_y0, &
                      !
                      mx_r_z0, mx_r_x0, mx_r_y0, &
                      my_r_z0, my_r_x0, my_r_y0, &
                      mz_r_z0, mz_r_x0, mz_r_y0, &
                      !
                      wx_r_z0, wx_r_x0, wx_r_y0, &
                      wy_r_z0, wy_r_x0, wy_r_y0, &
                      wz_r_z0, wz_r_x0, wz_r_y0, &
                      !
                      bx_r_z0, bx_r_x0, bx_r_y0, &
                      by_r_z0, by_r_x0, by_r_y0, &
                      bz_r_z0, bz_r_x0, bz_r_y0, &
                      !
                      jx_r_z0, jx_r_x0, jx_r_y0, &
                      jy_r_z0, jy_r_x0, jy_r_y0, &
                      jz_r_z0, jz_r_x0, jz_r_y0, &
                      !
                      rho_kxy, rho_kyz, rho_kxz, &
                       u2_kxy,  u2_kyz,  u2_kxz, &
                       b2_kxy,  b2_kyz,  b2_kxz  &
                    )
    use grid, only: nlx, nly, nlz, nkx, nky, nkz
    use grid, only: nlx_local, nky_local
    use grid, only: xx, yy, zz, kx, ky, kz
    use time, only: tt
    use shearing_box, only: tsc
    use mp, only: proc0, iproc_fft
    use mpiio, only: mpiio_write_var_2d
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_2D
    implicit none

    real(8), intent(in) :: rho_r_z0(1:nlx_local, 1:nly      ), &
                            rho_r_x0(1:nly      , 1:nlz      ), &
                            rho_r_y0(1:nlx_local, 1:nlz      )

    real(8), intent(in) ::  mx_r_z0(1:nlx_local, 1:nly      ), &
                             mx_r_x0(1:nly      , 1:nlz      ), &
                             mx_r_y0(1:nlx_local, 1:nlz      )
    real(8), intent(in) ::  my_r_z0(1:nlx_local, 1:nly      ), &
                             my_r_x0(1:nly      , 1:nlz      ), &
                             my_r_y0(1:nlx_local, 1:nlz      )
    real(8), intent(in) ::  mz_r_z0(1:nlx_local, 1:nly      ), &
                             mz_r_x0(1:nly      , 1:nlz      ), &
                             mz_r_y0(1:nlx_local, 1:nlz      )
                                                                                   
    real(8), intent(in) ::  wx_r_z0(1:nlx_local, 1:nly      ), &
                             wx_r_x0(1:nly      , 1:nlz      ), &
                             wx_r_y0(1:nlx_local, 1:nlz      )
    real(8), intent(in) ::  wy_r_z0(1:nlx_local, 1:nly      ), &
                             wy_r_x0(1:nly      , 1:nlz      ), &
                             wy_r_y0(1:nlx_local, 1:nlz      )
    real(8), intent(in) ::  wz_r_z0(1:nlx_local, 1:nly      ), &
                             wz_r_x0(1:nly      , 1:nlz      ), &
                             wz_r_y0(1:nlx_local, 1:nlz      )
                                                                            
    real(8), intent(in) ::  bx_r_z0(1:nlx_local, 1:nly      ), &
                             bx_r_x0(1:nly      , 1:nlz      ), &
                             bx_r_y0(1:nlx_local, 1:nlz      )
    real(8), intent(in) ::  by_r_z0(1:nlx_local, 1:nly      ), &
                             by_r_x0(1:nly      , 1:nlz      ), &
                             by_r_y0(1:nlx_local, 1:nlz      )
    real(8), intent(in) ::  bz_r_z0(1:nlx_local, 1:nly      ), &
                             bz_r_x0(1:nly      , 1:nlz      ), &
                             bz_r_y0(1:nlx_local, 1:nlz      )
                                                                            
    real(8), intent(in) ::  jx_r_z0(1:nlx_local, 1:nly      ), &
                             jx_r_x0(1:nly      , 1:nlz      ), &
                             jx_r_y0(1:nlx_local, 1:nlz      )
    real(8), intent(in) ::  jy_r_z0(1:nlx_local, 1:nly      ), &
                             jy_r_x0(1:nly      , 1:nlz      ), &
                             jy_r_y0(1:nlx_local, 1:nlz      )
    real(8), intent(in) ::  jz_r_z0(1:nlx_local, 1:nly      ), &
                             jz_r_x0(1:nly      , 1:nlz      ), &
                             jz_r_y0(1:nlx_local, 1:nlz      )

    real(8), intent(in) :: rho_kxy(1:nkx      , 1:nky_local), &
                            rho_kyz(1:nky_local, 1:nkz      ), &
                            rho_kxz(1:nkx      , 1:nkz      )
    real(8), intent(in) ::  u2_kxy(1:nkx      , 1:nky_local), &
                             u2_kyz(1:nky_local, 1:nkz      ), &
                             u2_kxz(1:nkx      , 1:nkz      )
    real(8), intent(in) ::  b2_kxy(1:nkx      , 1:nky_local), &
                             b2_kyz(1:nky_local, 1:nkz      ), &
                             b2_kxz(1:nkx      , 1:nkz      )

    integer, dimension(2) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_2D)

    !--------------------------------------------------!
    !                    z = 0 cut
    !--------------------------------------------------!
    ! z is not split in real space: every rank has z = 0 and writes its x slab
    if(any(zz == 0.d0)) then ! only the processe that has z = 0 write
      sizes(1) = nlx
      sizes(2) = nly
      subsizes(1) = nlx_local
      subsizes(2) = nly
      starts(1) = nlx_local*iproc_fft
      starts(2) = 0
    else
      sizes(1) = nlx
      sizes(2) = nly
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_rho_r_z0, disp_rho_r_z0, sizes, subsizes, starts, rho_r_z0)

    call mpiio_write_var_2d(fh_mx_r_z0 , disp_mx_r_z0 , sizes, subsizes, starts, mx_r_z0 )
    call mpiio_write_var_2d(fh_my_r_z0 , disp_my_r_z0 , sizes, subsizes, starts, my_r_z0 )
    call mpiio_write_var_2d(fh_mz_r_z0 , disp_mz_r_z0 , sizes, subsizes, starts, mz_r_z0 )
                                                                                         
    call mpiio_write_var_2d(fh_wx_r_z0 , disp_wx_r_z0 , sizes, subsizes, starts, wx_r_z0 )
    call mpiio_write_var_2d(fh_wy_r_z0 , disp_wy_r_z0 , sizes, subsizes, starts, wy_r_z0 )
    call mpiio_write_var_2d(fh_wz_r_z0 , disp_wz_r_z0 , sizes, subsizes, starts, wz_r_z0 )
                                                                                         
    call mpiio_write_var_2d(fh_bx_r_z0 , disp_bx_r_z0 , sizes, subsizes, starts, bx_r_z0 )
    call mpiio_write_var_2d(fh_by_r_z0 , disp_by_r_z0 , sizes, subsizes, starts, by_r_z0 )
    call mpiio_write_var_2d(fh_bz_r_z0 , disp_bz_r_z0 , sizes, subsizes, starts, bz_r_z0 )
                                                                                         
    call mpiio_write_var_2d(fh_jx_r_z0 , disp_jx_r_z0 , sizes, subsizes, starts, jx_r_z0 )
    call mpiio_write_var_2d(fh_jy_r_z0 , disp_jy_r_z0 , sizes, subsizes, starts, jy_r_z0 )
    call mpiio_write_var_2d(fh_jz_r_z0 , disp_jz_r_z0 , sizes, subsizes, starts, jz_r_z0 )

    !--------------------------------------------------!
    !                    x = 0 cut
    !--------------------------------------------------!
    ! x is split in real space: only iproc_fft = 0 has x = 0. The other ranks
    ! write one dummy element at (0, 0), as in Calliope and Gallope EMHD.
    if(any(xx == 0.d0)) then ! only the processe that has x = 0 write
      sizes(1) = nly
      sizes(2) = nlz
      subsizes(1) = nly
      subsizes(2) = nlz
      starts(1) = 0
      starts(2) = 0
    else
      sizes(1) = nly
      sizes(2) = nlz
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_rho_r_x0, disp_rho_r_x0, sizes, subsizes, starts, rho_r_x0)

    call mpiio_write_var_2d(fh_mx_r_x0 , disp_mx_r_x0 , sizes, subsizes, starts, mx_r_x0 )
    call mpiio_write_var_2d(fh_my_r_x0 , disp_my_r_x0 , sizes, subsizes, starts, my_r_x0 )
    call mpiio_write_var_2d(fh_mz_r_x0 , disp_mz_r_x0 , sizes, subsizes, starts, mz_r_x0 )
                                                                                         
    call mpiio_write_var_2d(fh_wx_r_x0 , disp_wx_r_x0 , sizes, subsizes, starts, wx_r_x0 )
    call mpiio_write_var_2d(fh_wy_r_x0 , disp_wy_r_x0 , sizes, subsizes, starts, wy_r_x0 )
    call mpiio_write_var_2d(fh_wz_r_x0 , disp_wz_r_x0 , sizes, subsizes, starts, wz_r_x0 )
                                                                                         
    call mpiio_write_var_2d(fh_bx_r_x0 , disp_bx_r_x0 , sizes, subsizes, starts, bx_r_x0 )
    call mpiio_write_var_2d(fh_by_r_x0 , disp_by_r_x0 , sizes, subsizes, starts, by_r_x0 )
    call mpiio_write_var_2d(fh_bz_r_x0 , disp_bz_r_x0 , sizes, subsizes, starts, bz_r_x0 )
                                                                                         
    call mpiio_write_var_2d(fh_jx_r_x0 , disp_jx_r_x0 , sizes, subsizes, starts, jx_r_x0 )
    call mpiio_write_var_2d(fh_jy_r_x0 , disp_jy_r_x0 , sizes, subsizes, starts, jy_r_x0 )
    call mpiio_write_var_2d(fh_jz_r_x0 , disp_jz_r_x0 , sizes, subsizes, starts, jz_r_x0 )

    !--------------------------------------------------!
    !                    y = 0 cut
    !--------------------------------------------------!
    ! y is not split in real space: every rank has y = 0 and writes its x slab
    if(any(yy == 0.d0)) then ! only the processe that has y = 0 write
      sizes(1) = nlx
      sizes(2) = nlz
      subsizes(1) = nlx_local
      subsizes(2) = nlz
      starts(1) = nlx_local*iproc_fft
      starts(2) = 0
    else
      sizes(1) = nlx
      sizes(2) = nlz
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_rho_r_y0, disp_rho_r_y0, sizes, subsizes, starts, rho_r_y0)

    call mpiio_write_var_2d(fh_mx_r_y0 , disp_mx_r_y0 , sizes, subsizes, starts, mx_r_y0 )
    call mpiio_write_var_2d(fh_my_r_y0 , disp_my_r_y0 , sizes, subsizes, starts, my_r_y0 )
    call mpiio_write_var_2d(fh_mz_r_y0 , disp_mz_r_y0 , sizes, subsizes, starts, mz_r_y0 )
                                                                                         
    call mpiio_write_var_2d(fh_wx_r_y0 , disp_wx_r_y0 , sizes, subsizes, starts, wx_r_y0 )
    call mpiio_write_var_2d(fh_wy_r_y0 , disp_wy_r_y0 , sizes, subsizes, starts, wy_r_y0 )
    call mpiio_write_var_2d(fh_wz_r_y0 , disp_wz_r_y0 , sizes, subsizes, starts, wz_r_y0 )
                                                                                         
    call mpiio_write_var_2d(fh_bx_r_y0 , disp_bx_r_y0 , sizes, subsizes, starts, bx_r_y0 )
    call mpiio_write_var_2d(fh_by_r_y0 , disp_by_r_y0 , sizes, subsizes, starts, by_r_y0 )
    call mpiio_write_var_2d(fh_bz_r_y0 , disp_bz_r_y0 , sizes, subsizes, starts, bz_r_y0 )
                                                                                         
    call mpiio_write_var_2d(fh_jx_r_y0 , disp_jx_r_y0 , sizes, subsizes, starts, jx_r_y0 )
    call mpiio_write_var_2d(fh_jy_r_y0 , disp_jy_r_y0 , sizes, subsizes, starts, jy_r_y0 )
    call mpiio_write_var_2d(fh_jz_r_y0 , disp_jz_r_y0 , sizes, subsizes, starts, jz_r_y0 )

    !--------------------------------------------------!
    !                      kz sum
    !--------------------------------------------------!
    ! kz is not split in Fourier space: the sum over kz is local and
    ! every rank writes its ky slab
    if(any(kz == 0.d0)) then ! only the processe that has kz = 0 write
      sizes(1) = nkx
      sizes(2) = nky
      subsizes(1) = nkx
      subsizes(2) = nky_local
      starts(1) = 0
      starts(2) = nky_local*iproc_fft
    else
      sizes(1) = nkx
      sizes(2) = nky
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_rho_kxy, disp_rho_kxy, sizes, subsizes, starts, rho_kxy)
    call mpiio_write_var_2d(fh_u2_kxy , disp_u2_kxy , sizes, subsizes, starts, u2_kxy )
    call mpiio_write_var_2d(fh_b2_kxy , disp_b2_kxy , sizes, subsizes, starts, b2_kxy )

    !--------------------------------------------------!
    !                      kx sum
    !--------------------------------------------------!
    ! kx is not split in Fourier space: the sum over kx is local and
    ! every rank writes its ky slab
    if(any(kx == 0.d0)) then ! only the processe that has kx = 0 write
      sizes(1) = nky
      sizes(2) = nkz
      subsizes(1) = nky_local
      subsizes(2) = nkz
      starts(1) = nky_local*iproc_fft
      starts(2) = 0
    else
      sizes(1) = nky
      sizes(2) = nkz
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_rho_kyz, disp_rho_kyz, sizes, subsizes, starts, rho_kyz)
    call mpiio_write_var_2d(fh_u2_kyz , disp_u2_kyz , sizes, subsizes, starts, u2_kyz )
    call mpiio_write_var_2d(fh_b2_kyz , disp_b2_kyz , sizes, subsizes, starts, b2_kyz )

    !--------------------------------------------------!
    !                      ky sum
    !--------------------------------------------------!
    ! ky is split in Fourier space: the sum over ky is reduced over comm_fft
    ! and broadcast by the diagnostics, so that iproc_fft = 0 (ky = 0) writes
    ! the whole plane and the dummy element (0, 0) of the other ranks holds
    ! the same value
    if(any(ky == 0.d0)) then ! only the processe that has ky = 0 write
      sizes(1) = nkx
      sizes(2) = nkz
      subsizes(1) = nkx
      subsizes(2) = nkz
      starts(1) = 0
      starts(2) = 0
    else
      sizes(1) = nkx
      sizes(2) = nkz
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_rho_kxz, disp_rho_kxz, sizes, subsizes, starts, rho_kxz)
    call mpiio_write_var_2d(fh_u2_kxz , disp_u2_kxz , sizes, subsizes, starts, u2_kxz )
    call mpiio_write_var_2d(fh_b2_kxz , disp_b2_kxz , sizes, subsizes, starts, b2_kxz )

    if(proc0) then
      write (unit=out2d_time_unit, fmt="(100es30.21)") tt, tsc
      flush (out2d_time_unit)
    endif

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_2D)
  end subroutine loop_io_2D


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Append 3D field via MPIIO
!-----------------------------------------------!
  subroutine loop_io_3D
    use fields, only: lnrho, lnrho_to_rho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use mp, only: proc0, iproc_fft
    use time, only: tt
    use grid, only: nkx, nky, nkz
    use grid, only: nky_local
    use mpiio, only: mpiio_write_var
    use shearing_box, only: tsc
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_3D
    implicit none
    integer, dimension(3) :: sizes, subsizes, starts
    complex(8), allocatable, dimension(:,:,:) :: rho

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_3D)

    ! out3d/rho.dat keeps the density (not ln(rho)) so that the output format is unchanged
    allocate(rho(nkz, nky_local, nkx), source=(0.d0, 0.d0))
    !$acc enter data copyin(rho)
    call lnrho_to_rho(lnrho, rho)

    ! Gallope EMHD layout: (kz, ky, kx), ky split over comm_fft
    sizes(1) = nkz
    sizes(2) = nky
    sizes(3) = nkx
    subsizes(1) = nkz
    subsizes(2) = nky_local
    subsizes(3) = nkx
    starts(1) = 0
    starts(2) = nky_local*iproc_fft
    starts(3) = 0

    ! MPI-IO writes the host copies
    !$acc update host(rho, mx, my, mz, bx, by, bz)

    call mpiio_write_var(fh_rho, disp_rho, sizes, subsizes, starts, rho)
    call mpiio_write_var(fh_mx , disp_mx , sizes, subsizes, starts, mx )
    call mpiio_write_var(fh_my , disp_my , sizes, subsizes, starts, my )
    call mpiio_write_var(fh_mz , disp_mz , sizes, subsizes, starts, mz )
    call mpiio_write_var(fh_bx , disp_bx , sizes, subsizes, starts, bx )
    call mpiio_write_var(fh_by , disp_by , sizes, subsizes, starts, by )
    call mpiio_write_var(fh_bz , disp_bz , sizes, subsizes, starts, bz )
    !$acc exit data delete(rho)
    deallocate(rho)

    if(proc0) then
      write (unit=out3d_time_unit, fmt="(100es30.21)") tt, tsc
      flush (out3d_time_unit)
    endif

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_3D)
  end subroutine loop_io_3D


!-----------------------------------------------!
!> @author  YK
!! @date    16 Jan 2019
!! @brief   Save restart file via MPIIO
!-----------------------------------------------!
  subroutine save_restart
    use fields, only: lnrho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use mp, only: proc0, iproc_fft, comm_fft
    use grid, only: nkx, nky, nkz
    use grid, only: nky_local
    use time, only: tt, dt
    use params, only: restart_dir
    use file, only: open_output_file, close_file
    use mpiio, only: mpiio_write_one
    use shearing_box, only: tsc
    use time_stamp, only: put_time_stamp, timer_save_restart
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_save_restart)

    ! Gallope EMHD layout: (kz, ky, kx), ky split over comm_fft
    ! (read back by fields::restart with the same sizes)
    sizes(1) = nkz
    sizes(2) = nky
    sizes(3) = nkx
    subsizes(1) = nkz
    subsizes(2) = nky_local
    subsizes(3) = nkx
    starts(1) = 0
    starts(2) = nky_local*iproc_fft
    starts(3) = 0

    ! MPI-IO writes the host copies
    !$acc update host(lnrho, mx, my, mz, bx, by, bz)

    call mpiio_write_one(lnrho, sizes, subsizes, starts, trim(restart_dir)//'lnrho.dat', comm_fft)
    call mpiio_write_one(mx , sizes, subsizes, starts, trim(restart_dir)//'mx.dat' , comm_fft)
    call mpiio_write_one(my , sizes, subsizes, starts, trim(restart_dir)//'my.dat' , comm_fft)
    call mpiio_write_one(mz , sizes, subsizes, starts, trim(restart_dir)//'mz.dat' , comm_fft)
    call mpiio_write_one(bx , sizes, subsizes, starts, trim(restart_dir)//'bx.dat' , comm_fft)
    call mpiio_write_one(by , sizes, subsizes, starts, trim(restart_dir)//'by.dat' , comm_fft)
    call mpiio_write_one(bz , sizes, subsizes, starts, trim(restart_dir)//'bz.dat' , comm_fft)

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
!! @date    29 Dec 2018
!! @brief   Finalization of NETCDF
!-----------------------------------------------!
  subroutine finish_io
    use mp, only: proc0
    use file, only: close_file
    implicit none
    integer :: ierr

    !3D
    call MPI_FILE_CLOSE(fh_rho,ierr)
    call MPI_FILE_CLOSE(fh_mx ,ierr)
    call MPI_FILE_CLOSE(fh_my ,ierr)
    call MPI_FILE_CLOSE(fh_mz ,ierr)
    call MPI_FILE_CLOSE(fh_bx ,ierr)
    call MPI_FILE_CLOSE(fh_by ,ierr)
    call MPI_FILE_CLOSE(fh_bz ,ierr)


    !2D
    call MPI_FILE_CLOSE(fh_rho_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_rho_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_rho_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_mx_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_mx_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_mx_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_my_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_my_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_my_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_mz_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_mz_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_mz_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_wx_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_wx_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_wx_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_wy_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_wy_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_wy_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_wz_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_wz_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_wz_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_bx_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_bx_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_bx_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_by_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_by_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_by_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_bz_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_bz_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_bz_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_jx_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_jx_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_jx_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_jy_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_jy_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_jy_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_jz_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_jz_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_jz_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_rho_kxy,ierr)
    call MPI_FILE_CLOSE(fh_rho_kyz,ierr)
    call MPI_FILE_CLOSE(fh_rho_kxz,ierr)

    call MPI_FILE_CLOSE(fh_u2_kxy,ierr)
    call MPI_FILE_CLOSE(fh_u2_kyz,ierr)
    call MPI_FILE_CLOSE(fh_u2_kxz,ierr)

    call MPI_FILE_CLOSE(fh_b2_kxy,ierr)
    call MPI_FILE_CLOSE(fh_b2_kyz,ierr)
    call MPI_FILE_CLOSE(fh_b2_kxz,ierr)

    if(proc0) then
      call close_file (out2d_time_unit)
      call close_file (out3d_time_unit)
    endif

    if(proc0) then
      status = nf90_close (ncid)
    endif

  end subroutine finish_io

end module io

