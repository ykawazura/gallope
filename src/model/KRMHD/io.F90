!-----------------------------------------------!
!> @author  YK
!! @brief   IO for RMHD
!-----------------------------------------------!
module io
  use netcdf
  use MPI
  implicit none

  public :: init_io, finish_io, loop_io, loop_io_2D, loop_io_3D, save_restart

  private

  ! MPIIO for 3D
  integer :: fh_phi, fh_psi
  character(len=100) :: filename
  integer (kind=MPI_OFFSET_KIND) :: disp_phi
  integer (kind=MPI_OFFSET_KIND) :: disp_psi
  integer :: out3d_time_unit

  ! MPIIO for 2D
  integer :: fh_phi_r_z0, fh_phi_r_x0, fh_phi_r_y0
  integer :: fh_psi_r_z0, fh_psi_r_x0, fh_psi_r_y0
  integer :: fh_omg_r_z0, fh_omg_r_x0, fh_omg_r_y0
  integer :: fh_jpa_r_z0, fh_jpa_r_x0, fh_jpa_r_y0
  integer :: fh_ux_r_z0 , fh_ux_r_x0 , fh_ux_r_y0, &
             fh_uy_r_z0 , fh_uy_r_x0 , fh_uy_r_y0
  integer :: fh_bx_r_z0 , fh_bx_r_x0 , fh_bx_r_y0, &
             fh_by_r_z0 , fh_by_r_x0 , fh_by_r_y0

  integer (kind=MPI_OFFSET_KIND) :: disp_phi_r_z0, disp_phi_r_x0, disp_phi_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_psi_r_z0, disp_psi_r_x0, disp_psi_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_omg_r_z0, disp_omg_r_x0, disp_omg_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_jpa_r_z0, disp_jpa_r_x0, disp_jpa_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_ux_r_z0, disp_ux_r_x0, disp_ux_r_y0, &
                                    disp_uy_r_z0, disp_uy_r_x0, disp_uy_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_bx_r_z0, disp_bx_r_x0, disp_bx_r_y0, &
                                    disp_by_r_z0, disp_by_r_x0, disp_by_r_y0
  integer :: out2d_time_unit

  ! NETCDF for regular output file
  integer :: status
  integer, parameter :: kind_nf = kind (NF90_NOERR)
  integer (kind_nf) :: ncid
  integer :: run_id
  integer (kind_nf) :: char10_dim
  ! parameter
  integer :: nupe_x_id, nupe_x_exp_id, nupe_z_id, nupe_z_exp_id
  integer :: etape_x_id, etape_x_exp_id, etape_z_id, etape_z_exp_id
  ! coordinate
  integer :: xx_id, yy_id, zz_id, kx_id, ky_id, kz_id, kpbin_id, tt_id
  ! total energy
  integer :: upe2_sum_id , bpe2_sum_id
  integer :: zppe2_sum_id, zmpe2_sum_id
  integer :: upe2dot_sum_id, bpe2dot_sum_id
  integer :: upe2dissip_sum_id, bpe2dissip_sum_id
  integer :: p_phi_sum_id, p_psi_sum_id, p_xhl_sum_id
  ! polar spectrum
  integer :: upe2_bin_id , bpe2_bin_id
  integer :: zppe2_bin_id, zmpe2_bin_id
  ! Hermite (g) spectrum + free energy (Meyrand 2019)
  integer :: mm_id, W_free_id, W_m_id, g2_bin_id
  ! g free-energy power balance: injection P_g and dissipation D_g (dW_free/dt = P_g - D_g)
  integer :: p_g_sum_id, Dg_sum_id
  ! Hermite free-energy flux Gamma_m (Meyrand 2019 eq 9; gated by write_hermite_flux)
  integer :: Gamma_m_id
  ! k-integrated Hermite flux Gamma(m) (telescoping invariant; gated likewise)
  integer :: Gamma_m_kint_id

  integer (kind_nf) :: xx_dim, yy_dim, zz_dim, kx_dim, ky_dim, kz_dim, kpbin_dim, tt_dim
  integer (kind_nf) :: mm_dim
  integer, dimension (3) :: bin_dim
  integer, dimension (4) :: bin4_dim

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

    ! Only comm_fft group 0 performs MPI-IO for the redundantly-solved fields;
    ! every group holds an identical copy, so opening/writing on comm_fft from
    ! group 0 alone avoids cross-group file corruption. At P_m=P_s=1 the gate is
    ! always true and comm_fft==MPI_COMM_WORLD, so behaviour is unchanged.
    if (proc0_m .and. proc0_s) then
      call set_file_handle('out3d/phi.dat', fh_phi, disp_phi, comm_fft)
      call set_file_handle('out3d/psi.dat', fh_psi, disp_psi, comm_fft)
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

    ! Only comm_fft group 0 performs the 2D-cut MPI-IO (see init_io_decomp_3d).
    if (proc0_m .and. proc0_s) then
    call set_file_handle('out2d/phi_r_z0.dat', fh_phi_r_z0, disp_phi_r_z0, comm_fft)
    call set_file_handle('out2d/phi_r_x0.dat', fh_phi_r_x0, disp_phi_r_x0, comm_fft)
    call set_file_handle('out2d/phi_r_y0.dat', fh_phi_r_y0, disp_phi_r_y0, comm_fft)

    call set_file_handle('out2d/psi_r_z0.dat', fh_psi_r_z0, disp_psi_r_z0, comm_fft)
    call set_file_handle('out2d/psi_r_x0.dat', fh_psi_r_x0, disp_psi_r_x0, comm_fft)
    call set_file_handle('out2d/psi_r_y0.dat', fh_psi_r_y0, disp_psi_r_y0, comm_fft)

    call set_file_handle('out2d/omg_r_z0.dat', fh_omg_r_z0, disp_omg_r_z0, comm_fft)
    call set_file_handle('out2d/omg_r_x0.dat', fh_omg_r_x0, disp_omg_r_x0, comm_fft)
    call set_file_handle('out2d/omg_r_y0.dat', fh_omg_r_y0, disp_omg_r_y0, comm_fft)

    call set_file_handle('out2d/jpa_r_z0.dat', fh_jpa_r_z0, disp_jpa_r_z0, comm_fft)
    call set_file_handle('out2d/jpa_r_x0.dat', fh_jpa_r_x0, disp_jpa_r_x0, comm_fft)
    call set_file_handle('out2d/jpa_r_y0.dat', fh_jpa_r_y0, disp_jpa_r_y0, comm_fft)

    call set_file_handle('out2d/ux_r_z0.dat' , fh_ux_r_z0 , disp_ux_r_z0 , comm_fft)
    call set_file_handle('out2d/ux_r_x0.dat' , fh_ux_r_x0 , disp_ux_r_x0 , comm_fft)
    call set_file_handle('out2d/ux_r_y0.dat' , fh_ux_r_y0 , disp_ux_r_y0 , comm_fft)
                                                                         
    call set_file_handle('out2d/uy_r_z0.dat' , fh_uy_r_z0 , disp_uy_r_z0 , comm_fft)
    call set_file_handle('out2d/uy_r_x0.dat' , fh_uy_r_x0 , disp_uy_r_x0 , comm_fft)
    call set_file_handle('out2d/uy_r_y0.dat' , fh_uy_r_y0 , disp_uy_r_y0 , comm_fft)
                                                                         
    call set_file_handle('out2d/bx_r_z0.dat' , fh_bx_r_z0 , disp_bx_r_z0 , comm_fft)
    call set_file_handle('out2d/bx_r_x0.dat' , fh_bx_r_x0 , disp_bx_r_x0 , comm_fft)
    call set_file_handle('out2d/bx_r_y0.dat' , fh_bx_r_y0 , disp_bx_r_y0 , comm_fft)
                                                                         
    call set_file_handle('out2d/by_r_z0.dat' , fh_by_r_z0 , disp_by_r_z0 , comm_fft)
    call set_file_handle('out2d/by_r_x0.dat' , fh_by_r_x0 , disp_by_r_x0 , comm_fft)
    call set_file_handle('out2d/by_r_y0.dat' , fh_by_r_y0 , disp_by_r_y0 , comm_fft)
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
    use grid, only: nlx, nly, nlz
    use grid, only: xx_global, yy, zz, kx, ky_global, kz
    use grid, only: nm
    use mp, only: proc0
    use params, only: runname, write_hermite_flux, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp, &
                      etape_x, etape_x_exp, etape_z, etape_z_exp
    implicit none
    integer, intent(in) :: nkpolar
    real(8), intent(in) :: kpbin(1:nkpolar)
    integer :: im

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

      status = nf90_def_dim (ncid, 'yy', size(yy), yy_dim)
      status = nf90_def_dim (ncid, 'xx', size(xx_global), xx_dim)
      status = nf90_def_dim (ncid, 'zz', size(zz), zz_dim)
      status = nf90_def_dim (ncid, 'kx', size(kx), kx_dim)
      status = nf90_def_dim (ncid, 'ky', size(ky_global), ky_dim)
      status = nf90_def_dim (ncid, 'kz', size(kz), kz_dim)
      status = nf90_def_dim (ncid, 'kpbin', size(kpbin), kpbin_dim)
      status = nf90_def_dim (ncid, 'mm', nm, mm_dim)
      status = nf90_def_dim (ncid, 'tt', NF90_UNLIMITED, tt_dim)

      status = nf90_def_var (ncid, 'nupe_x', NF90_DOUBLE, nupe_x_id)
      status = nf90_def_var (ncid, 'nupe_x_exp', NF90_DOUBLE, nupe_x_exp_id)
      status = nf90_def_var (ncid, 'nupe_z', NF90_DOUBLE, nupe_z_id)
      status = nf90_def_var (ncid, 'nupe_z_exp', NF90_DOUBLE, nupe_z_exp_id)
      status = nf90_def_var (ncid, 'etape_x', NF90_DOUBLE, etape_x_id)
      status = nf90_def_var (ncid, 'etape_x_exp', NF90_DOUBLE, etape_x_exp_id)
      status = nf90_def_var (ncid, 'etape_z', NF90_DOUBLE, etape_z_id)
      status = nf90_def_var (ncid, 'etape_z_exp', NF90_DOUBLE, etape_z_exp_id)

      status = nf90_def_var (ncid, 'xx', NF90_DOUBLE, xx_dim, xx_id)
      status = nf90_def_var (ncid, 'yy', NF90_DOUBLE, yy_dim, yy_id)
      status = nf90_def_var (ncid, 'zz', NF90_DOUBLE, zz_dim, zz_id)
      status = nf90_def_var (ncid, 'kx', NF90_DOUBLE, kx_dim, kx_id)
      status = nf90_def_var (ncid, 'ky', NF90_DOUBLE, ky_dim, ky_id)
      status = nf90_def_var (ncid, 'kz', NF90_DOUBLE, kz_dim, kz_id)
      status = nf90_def_var (ncid, 'kpbin', NF90_DOUBLE, kpbin_dim, kpbin_id)
      status = nf90_def_var (ncid, 'mm', NF90_DOUBLE, mm_dim, mm_id)
      status = nf90_def_var (ncid, 'tt', NF90_DOUBLE, tt_dim, tt_id)
      ! total energy
      status = nf90_def_var (ncid, 'upe2_sum' , NF90_DOUBLE, tt_dim, upe2_sum_id )
      status = nf90_def_var (ncid, 'bpe2_sum' , NF90_DOUBLE, tt_dim, bpe2_sum_id )
      status = nf90_def_var (ncid, 'zppe2_sum', NF90_DOUBLE, tt_dim, zppe2_sum_id)
      status = nf90_def_var (ncid, 'zmpe2_sum', NF90_DOUBLE, tt_dim, zmpe2_sum_id)
      status = nf90_def_var (ncid, 'upe2dot_sum', NF90_DOUBLE, tt_dim, upe2dot_sum_id)
      status = nf90_def_var (ncid, 'bpe2dot_sum', NF90_DOUBLE, tt_dim, bpe2dot_sum_id)
      status = nf90_def_var (ncid, 'upe2dissip_sum', NF90_DOUBLE, tt_dim, upe2dissip_sum_id)
      status = nf90_def_var (ncid, 'bpe2dissip_sum', NF90_DOUBLE, tt_dim, bpe2dissip_sum_id)
      status = nf90_def_var (ncid, 'p_phi_sum'     , NF90_DOUBLE, tt_dim, p_phi_sum_id   )
      status = nf90_def_var (ncid, 'p_psi_sum'     , NF90_DOUBLE, tt_dim, p_psi_sum_id   )
      status = nf90_def_var (ncid, 'p_xhl_sum'     , NF90_DOUBLE, tt_dim, p_xhl_sum_id   )
      ! polar spectrum
      bin_dim (1) = kpbin_dim
      bin_dim (2) = kz_dim
      bin_dim (3) = tt_dim
      status = nf90_def_var (ncid, 'upe2_bin' , NF90_DOUBLE, bin_dim, upe2_bin_id )
      status = nf90_def_var (ncid, 'bpe2_bin' , NF90_DOUBLE, bin_dim, bpe2_bin_id )
      status = nf90_def_var (ncid, 'zppe2_bin', NF90_DOUBLE, bin_dim, zppe2_bin_id)
      status = nf90_def_var (ncid, 'zmpe2_bin', NF90_DOUBLE, bin_dim, zmpe2_bin_id)
      ! Hermite (g) spectrum + free energy (Meyrand 2019 eq 6)
      bin4_dim (1) = kpbin_dim
      bin4_dim (2) = mm_dim
      bin4_dim (3) = kz_dim
      bin4_dim (4) = tt_dim
      status = nf90_def_var (ncid, 'W_free', NF90_DOUBLE, tt_dim, W_free_id)
      status = nf90_def_var (ncid, 'W_m'   , NF90_DOUBLE, (/mm_dim, tt_dim/), W_m_id)
      status = nf90_def_var (ncid, 'g2_bin', NF90_DOUBLE, bin4_dim, g2_bin_id)
      ! g free-energy power balance scalars (Meyrand 2019: dW_free/dt = P_g - D_g)
      status = nf90_def_var (ncid, 'p_g_sum', NF90_DOUBLE, tt_dim, p_g_sum_id)
      status = nf90_def_var (ncid, 'Dg_sum' , NF90_DOUBLE, tt_dim, Dg_sum_id )
      ! Hermite free-energy flux Gamma_m(kpbin,mm,kz,tt) (Meyrand 2019 eq 9)
      if (write_hermite_flux) then
        status = nf90_def_var (ncid, 'Gamma_m', NF90_DOUBLE, bin4_dim, Gamma_m_id)
        ! k-integrated flux Gamma(mm,tt): top slot = telescoping residual (~0)
        status = nf90_def_var (ncid, 'Gamma_m_kint', NF90_DOUBLE, &
                               (/mm_dim, tt_dim/), Gamma_m_kint_id)
      endif

      status = nf90_enddef (ncid)  ! out of definition mode

      status = nf90_put_var (ncid, nupe_x_id, nupe_x)
      status = nf90_put_var (ncid, nupe_x_exp_id, dble(nupe_x_exp))
      status = nf90_put_var (ncid, nupe_z_id, nupe_z)
      status = nf90_put_var (ncid, nupe_z_exp_id, dble(nupe_z_exp))
      status = nf90_put_var (ncid, etape_x_id, etape_x)
      status = nf90_put_var (ncid, etape_x_exp_id, dble(etape_x_exp))
      status = nf90_put_var (ncid, etape_z_id, etape_z)
      status = nf90_put_var (ncid, etape_z_exp_id, dble(etape_z_exp))

      status = nf90_put_var (ncid, xx_id, xx_global)
      status = nf90_put_var (ncid, yy_id, yy)
      status = nf90_put_var (ncid, zz_id, zz)
      status = nf90_put_var (ncid, kx_id, kx)
      status = nf90_put_var (ncid, ky_id, ky_global)
      status = nf90_put_var (ncid, kz_id, kz)
      status = nf90_put_var (ncid, kpbin_id, kpbin)
      ! Hermite index coordinate: physical m = 0, 1, ..., nm-1
      status = nf90_put_var (ncid, mm_id, (/(dble(im), im=0,nm-1)/))

      nout = 1
    endif
  end subroutine init_io_netcdf


!-----------------------------------------------!
!> @author  YK
!! @brief   Append variables to NETCDF
!           for params, time history & spectra
!-----------------------------------------------!
  subroutine loop_io( &
                      upe2_sum, bpe2_sum, &
                      upe2dot_sum, bpe2dot_sum, &
                      upe2dissip_sum, bpe2dissip_sum, &
                      p_phi_sum, p_psi_sum, p_xhl_sum, &
                      zppe2_sum, zmpe2_sum, &
                      !
                      nkpolar, &
                      upe2_bin , bpe2_bin , &
                      zppe2_bin, zmpe2_bin, &
                      !
                      nm, g2bin, Wm, W_free, p_g_sum, Dg_sum, Gamma_m, Gamma_m_kint &
                    )
    use time, only: tt
    use grid, only: nlx, nly, nlz, nkz
    use mp, only: proc0
    use params, only: write_hermite_flux
    use time_stamp, only: put_time_stamp, timer_io_total
    implicit none
    real(8), intent(in) :: upe2_sum , bpe2_sum
    real(8), intent(in) :: zppe2_sum, zmpe2_sum
    real(8), intent(in) :: upe2dot_sum, bpe2dot_sum
    real(8), intent(in) :: upe2dissip_sum, bpe2dissip_sum
    real(8), intent(in) :: p_phi_sum, p_psi_sum, p_xhl_sum

    integer, intent(in) :: nkpolar
    real(8), intent(in) :: upe2_bin (1:nkpolar, nkz), bpe2_bin (1:nkpolar, nkz)
    real(8), intent(in) :: zppe2_bin(1:nkpolar, nkz), zmpe2_bin(1:nkpolar, nkz)

    integer, intent(in) :: nm
    real(8), intent(in) :: g2bin(1:nkpolar, nm, nkz)
    real(8), intent(in) :: Wm(nm)
    real(8), intent(in) :: W_free
    real(8), intent(in) :: p_g_sum, Dg_sum
    real(8), intent(in) :: Gamma_m(1:nkpolar, nm, nkz)
    real(8), intent(in) :: Gamma_m_kint(nm)

    integer, dimension (3) :: start3, count3
    integer, dimension (4) :: start4, count4

    if (proc0) call put_time_stamp(timer_io_total)

    ! output via NETCDF
    if(proc0) then
      ! total energy
      status = nf90_put_var (ncid, tt_id, tt, start=(/nout/))
      status = nf90_put_var (ncid, upe2_sum_id , upe2_sum , start=(/nout/))
      status = nf90_put_var (ncid, bpe2_sum_id , bpe2_sum , start=(/nout/))
      status = nf90_put_var (ncid, zppe2_sum_id, zppe2_sum, start=(/nout/))
      status = nf90_put_var (ncid, zmpe2_sum_id, zmpe2_sum, start=(/nout/))
      status = nf90_put_var (ncid, upe2dot_sum_id, upe2dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, bpe2dot_sum_id, bpe2dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, upe2dissip_sum_id, upe2dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, bpe2dissip_sum_id, bpe2dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_phi_sum_id, p_phi_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_psi_sum_id, p_psi_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_xhl_sum_id, p_xhl_sum, start=(/nout/))
      ! mean magnetic field
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout

      count3(1) = nkpolar
      count3(2) = nkz
      count3(3) = 1
      status = nf90_put_var (ncid, upe2_bin_id , upe2_bin , start=start3, count=count3)
      status = nf90_put_var (ncid, bpe2_bin_id , bpe2_bin , start=start3, count=count3)
      status = nf90_put_var (ncid, zppe2_bin_id, zppe2_bin, start=start3, count=count3)
      status = nf90_put_var (ncid, zmpe2_bin_id, zmpe2_bin, start=start3, count=count3)
      ! Hermite (g) spectrum + free energy (Meyrand 2019 eq 6)
      status = nf90_put_var (ncid, W_free_id, W_free, start=(/nout/))
      status = nf90_put_var (ncid, p_g_sum_id, p_g_sum, start=(/nout/))
      status = nf90_put_var (ncid, Dg_sum_id , Dg_sum , start=(/nout/))
      status = nf90_put_var (ncid, W_m_id, Wm, start=(/1, nout/), count=(/nm, 1/))
      start4(1) = 1;       start4(2) = 1;  start4(3) = 1;   start4(4) = nout
      count4(1) = nkpolar; count4(2) = nm; count4(3) = nkz; count4(4) = 1
      status = nf90_put_var (ncid, g2_bin_id, g2bin, start=start4, count=count4)
      ! Hermite free-energy flux Gamma_m (eq 9), same rank-4 layout as g2_bin
      if (write_hermite_flux) then
        status = nf90_put_var (ncid, Gamma_m_id, Gamma_m, start=start4, count=count4)
        status = nf90_put_var (ncid, Gamma_m_kint_id, Gamma_m_kint, &
                               start=(/1, nout/), count=(/nm, 1/))
      endif

      status = nf90_sync (ncid)

      nout = nout + 1
    endif

    if (proc0) call put_time_stamp(timer_io_total)
  end subroutine loop_io


!-----------------------------------------------!
!> @author  YK
!! @brief   Append cross section via MPIIO
!-----------------------------------------------!
  subroutine loop_io_2D( &
                      phi_r_z0, phi_r_x0, phi_r_y0, &
                      psi_r_z0, psi_r_x0, psi_r_y0, &
                      omg_r_z0, omg_r_x0, omg_r_y0, &
                      jpa_r_z0, jpa_r_x0, jpa_r_y0, &
                       ux_r_z0,  ux_r_x0,  ux_r_y0, &
                       uy_r_z0,  uy_r_x0,  uy_r_y0, &
                       bx_r_z0,  bx_r_x0,  bx_r_y0, &
                       by_r_z0,  by_r_x0,  by_r_y0  &
                    )
    use grid, only: xx, yy, zz, kx, ky, kz
    use grid, only: nlx, nly, nlz, nkx, nky, nkz
    use grid, only: nlx_local, nky_local
    use time, only: tt
    use mp, only: proc0, iproc_fft, proc0_m, proc0_s
    use mpiio, only: mpiio_write_var_2d
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_2D
    implicit none

    real(8), intent(in) :: phi_r_z0(1:nlx_local, 1:nly), &
                           phi_r_x0(1:nly      , 1:nlz), &
                           phi_r_y0(1:nlx_local, 1:nlz)

    real(8), intent(in) :: psi_r_z0(1:nlx_local, 1:nly), &
                           psi_r_x0(1:nly      , 1:nlz), &
                           psi_r_y0(1:nlx_local, 1:nlz)

    real(8), intent(in) :: omg_r_z0(1:nlx_local, 1:nly), &
                           omg_r_x0(1:nly      , 1:nlz), &
                           omg_r_y0(1:nlx_local, 1:nlz)

    real(8), intent(in) :: jpa_r_z0(1:nlx_local, 1:nly), &
                           jpa_r_x0(1:nly      , 1:nlz), &
                           jpa_r_y0(1:nlx_local, 1:nlz)

    real(8), intent(in) :: ux_r_z0(1:nlx_local, 1:nly), &
                           ux_r_x0(1:nly      , 1:nlz), &
                           ux_r_y0(1:nlx_local, 1:nlz)
    real(8), intent(in) :: uy_r_z0(1:nlx_local, 1:nly), &
                           uy_r_x0(1:nly      , 1:nlz), &
                           uy_r_y0(1:nlx_local, 1:nlz)

    real(8), intent(in) :: bx_r_z0(1:nlx_local, 1:nly), &
                           bx_r_x0(1:nly      , 1:nlz), &
                           bx_r_y0(1:nlx_local, 1:nlz)
    real(8), intent(in) :: by_r_z0(1:nlx_local, 1:nly), &
                           by_r_x0(1:nly      , 1:nlz), &
                           by_r_y0(1:nlx_local, 1:nlz)

    integer, dimension(2) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_2D)

    ! Only comm_fft group 0 opened the 2D-cut files, so only it writes the cuts
    ! (collective over comm_fft). Other groups hold identical redundant copies.
    if (proc0_m .and. proc0_s) then

    !--------------------------------------------------!
    !                    z = 0 cut
    !--------------------------------------------------!
    if(any(zz == 0.0)) then ! only the processe that has z = 0 write
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

    call mpiio_write_var_2d(fh_phi_r_z0, disp_phi_r_z0, sizes, subsizes, starts, phi_r_z0)
    call mpiio_write_var_2d(fh_psi_r_z0, disp_psi_r_z0, sizes, subsizes, starts, psi_r_z0)
    call mpiio_write_var_2d(fh_omg_r_z0, disp_omg_r_z0, sizes, subsizes, starts, omg_r_z0)
    call mpiio_write_var_2d(fh_jpa_r_z0, disp_jpa_r_z0, sizes, subsizes, starts, jpa_r_z0)

    call mpiio_write_var_2d(fh_ux_r_z0 , disp_ux_r_z0 , sizes, subsizes, starts, ux_r_z0 )
    call mpiio_write_var_2d(fh_uy_r_z0 , disp_uy_r_z0 , sizes, subsizes, starts, uy_r_z0 )
    call mpiio_write_var_2d(fh_bx_r_z0 , disp_bx_r_z0 , sizes, subsizes, starts, bx_r_z0 )
    call mpiio_write_var_2d(fh_by_r_z0 , disp_by_r_z0 , sizes, subsizes, starts, by_r_z0 )

    !--------------------------------------------------!
    !                    x = 0 cut
    !--------------------------------------------------!
    if(any(xx == 0.0)) then ! only the processe that has z = 0 write
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

    call mpiio_write_var_2d(fh_phi_r_x0, disp_phi_r_x0, sizes, subsizes, starts, phi_r_x0)
    call mpiio_write_var_2d(fh_psi_r_x0, disp_psi_r_x0, sizes, subsizes, starts, psi_r_x0)
    call mpiio_write_var_2d(fh_omg_r_x0, disp_omg_r_x0, sizes, subsizes, starts, omg_r_x0)
    call mpiio_write_var_2d(fh_jpa_r_x0, disp_jpa_r_x0, sizes, subsizes, starts, jpa_r_x0)

    call mpiio_write_var_2d(fh_ux_r_x0 , disp_ux_r_x0 , sizes, subsizes, starts, ux_r_x0 )
    call mpiio_write_var_2d(fh_uy_r_x0 , disp_uy_r_x0 , sizes, subsizes, starts, uy_r_x0 )
    call mpiio_write_var_2d(fh_bx_r_x0 , disp_bx_r_x0 , sizes, subsizes, starts, bx_r_x0 )
    call mpiio_write_var_2d(fh_by_r_x0 , disp_by_r_x0 , sizes, subsizes, starts, by_r_x0 )

    !--------------------------------------------------!
    !                    y = 0 cut
    !--------------------------------------------------!
    if(any(yy == 0.0)) then ! only the processe that has z = 0 write
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

    call mpiio_write_var_2d(fh_phi_r_y0, disp_phi_r_y0, sizes, subsizes, starts, phi_r_y0)
    call mpiio_write_var_2d(fh_psi_r_y0, disp_psi_r_y0, sizes, subsizes, starts, psi_r_y0)
    call mpiio_write_var_2d(fh_omg_r_y0, disp_omg_r_y0, sizes, subsizes, starts, omg_r_y0)
    call mpiio_write_var_2d(fh_jpa_r_y0, disp_jpa_r_y0, sizes, subsizes, starts, jpa_r_y0)

    call mpiio_write_var_2d(fh_ux_r_y0 , disp_ux_r_y0 , sizes, subsizes, starts, ux_r_y0 )
    call mpiio_write_var_2d(fh_uy_r_y0 , disp_uy_r_y0 , sizes, subsizes, starts, uy_r_y0 )
    call mpiio_write_var_2d(fh_bx_r_y0 , disp_bx_r_y0 , sizes, subsizes, starts, bx_r_y0 )
    call mpiio_write_var_2d(fh_by_r_y0 , disp_by_r_y0 , sizes, subsizes, starts, by_r_y0 )

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
    use fields, only: phi, psi
    use mp, only: proc0
    use time, only: tt
    use grid, only: nkx, nky, nky_local, nkz
    use mpiio, only: mpiio_write_var
    use shearing_box, only: tsc
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_3D
    implicit none
    integer, dimension(3) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_3D)

    sizes(1) = nkz
    sizes(2) = nky
    sizes(3) = nkx
    subsizes(1) = nkz
    subsizes(2) = nky_local
    subsizes(3) = nkx
    starts(1) = 0
    starts(2) = 0
    starts(3) = 0

    !$acc update host(phi, psi)
    
    call mpiio_write_var(fh_phi, disp_phi, sizes, subsizes, starts, phi)
    call mpiio_write_var(fh_psi, disp_psi, sizes, subsizes, starts, psi)

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
!-----------------------------------------------!
  subroutine save_restart
    use fields, only: phi, omg, psi
    use fields, only: g
    use mp, only: proc0, iproc_fft, proc0_m, proc0_s, comm_fft, comm_fm
    use grid, only: nkx, nky, nky_local, nkz
    use grid, only: nm, nm_local, m_offset
    use time, only: tt, dt
    use params, only: restart_dir
    use file, only: open_output_file, close_file
    use mpiio, only: mpiio_write_one
    use shearing_box, only: tsc
    use time_stamp, only: put_time_stamp, timer_save_restart
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts
    integer, dimension(4) :: sizes4, subsizes4, starts4

    if (proc0) call put_time_stamp(timer_save_restart)

    sizes(1) = nkz
    sizes(2) = nky
    sizes(3) = nkx
    subsizes(1) = nkz
    subsizes(2) = nky_local
    subsizes(3) = nkx
    starts(1) = 0
    starts(2) = nky_local*iproc_fft
    starts(3) = 0

    !$acc update host(phi, omg, psi)
    !$acc update host(g)

    ! Only comm_fft group 0 writes the checkpoint; every group holds an identical
    ! copy, so writing from all groups would corrupt the same files.
    if (proc0_m .and. proc0_s) then
      call mpiio_write_one(phi, sizes, subsizes, starts, trim(restart_dir)//'phi.dat', comm_fft)
      call mpiio_write_one(omg, sizes, subsizes, starts, trim(restart_dir)//'omg.dat', comm_fft)
      call mpiio_write_one(psi, sizes, subsizes, starts, trim(restart_dir)//'psi.dat', comm_fft)
    endif

    ! g is distributed over the (comm_fft x comm_m) plane = comm_fm and is
    ! redundant across comm_s, so only the iproc_s==0 group writes it (collective
    ! over comm_fm), covering the full nky x nm decomposition in one file.
    sizes4    = (/ nkz, nky      , nkx, nm       /)
    subsizes4 = (/ nkz, nky_local, nkx, nm_local /)
    starts4   = (/ 0  , nky_local*iproc_fft, 0, m_offset /)
    if (proc0_s) then
      call mpiio_write_one(g(:,:,:,1:nm_local), sizes4, subsizes4, starts4, &
                           trim(restart_dir)//'g.dat', comm_fm)
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
    integer :: ierr

    ! Only comm_fft group 0 opened the MPI-IO handles (see init_io_decomp_*).
    if (proc0_m .and. proc0_s) then
    !3D
    call MPI_FILE_CLOSE(fh_phi,ierr)
    call MPI_FILE_CLOSE(fh_psi,ierr)


    !2D
    call MPI_FILE_CLOSE(fh_phi_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_phi_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_phi_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_psi_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_psi_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_psi_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_omg_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_omg_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_omg_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_jpa_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_jpa_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_jpa_r_y0,ierr)
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


