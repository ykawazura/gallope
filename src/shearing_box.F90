!-----------------------------------------------!
!> @author  YK
!! @brief   Shearing box setting
!-----------------------------------------------!
module shearing_box
  implicit none

  public init_shearing_time, remap_fld, get_imp_terms_tintg_with_shear
  public tsc ! time in shearing coordinate. must be 0 < tsc < tremap
  public tremap
  public nremap
  public k2t, k2t_inv, kxt, kxt_old1
  public timer_remap 
  public to_non_shearing_coordinate

  integer :: shear_flg = 0 ! 0 for wo shear, 1 for w shear
  real(8) :: tsc = 0.d0, tremap = 0.d0
  integer :: nremap = 0

  real(8) :: timer_remap(2) = 0.d0

  real(8), allocatable ::  k2t(:, :, :), k2t_inv(:, :, :), kxt(:, :), kxt_old1(:, :)

  !$acc declare create(shear_flg)
  !$acc declare create(tsc, tremap)
  !$acc declare create(nremap)
  !$acc declare create(k2t, k2t_inv, kxt, kxt_old1)
contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of time parameters,
!!          followed by input file reading
!-----------------------------------------------!
  subroutine init_shearing_time
    use mp, only: proc0
    use params, only: shear, q
    use grid, only: lx, ly, kx, ky, kz, k2, k2inv
    use grid, only: nkx, nky_local, nkz
    implicit none
    integer :: i, j, k

    allocate(k2t    , source=k2)
    allocate(k2t_inv, source=k2inv)

    allocate(kxt     (1:nkx, 1:nky_local))
    allocate(kxt_old1(1:nkx, 1:nky_local))
    do j = 1, nky_local
      do i = 1, nkx
        kxt(i, j) = kx(i)
      enddo
    enddo
    kxt_old1 = kxt

    if(shear) then
      shear_flg = 1
      tremap = ly/(lx*q*shear_flg)
      if(proc0) print '("shear is on; tremap = ", f5.3)', tremap

      do j = 1, nky_local
        do i = 1, nkx
          kxt(i, j) = kx(i) + q*shear_flg*tsc*ky(j)
          do k = 1, nkz
            k2t(k, j, i) = kxt(i, j)**2 + ky(j)**2 + kz(k)**2
            if(k2t(k, j, i) == 0.d0) then
              k2t_inv(k, j, i) = 0.d0
            else
              k2t_inv(k, j, i) = 1.0d0/k2t(k, j, i)
            endif
          enddo
        enddo
      enddo
      kxt_old1 = kxt
    endif

    !$acc update device(shear_flg)
    !$acc update device(tsc, tremap)
    !$acc update device(nremap)
    !$acc update device(k2t, k2t_inv, kxt, kxt_old1)

  end subroutine init_shearing_time



!-----------------------------------------------!
!> @author  YK
!! @brief   Remap fields
!-----------------------------------------------!
  subroutine remap_fld(fld)
    use params, only: pi, q
    use grid, only: nkx, nky_local, nkz
    use grid, only: ikx, iky
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: fld
    complex(8), dimension(:,:,:), allocatable :: tmp
    integer :: i, j, k, inew
    integer :: ikx_min, ikx_max, ikxnew
    integer, allocatable :: ikx_to_index(:)  ! look-up table

    ikx_min = minval(ikx)
    ikx_max = maxval(ikx)
    
    ! ikx -> index
    allocate(ikx_to_index(ikx_min:ikx_max))
    ikx_to_index = 0
    do i = 1, nkx
      ikx_to_index(ikx(i)) = i
    enddo

    allocate(tmp(nkz, nky_local, nkx), source=(0.d0, 0.d0))

    !$acc enter data copyin(ikx_to_index) create(tmp)
    !$acc data present(fld, ikx, iky, ikx_to_index, tmp)

    !$acc parallel loop collapse(3) private(ikxnew, inew)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ikxnew = ikx(i) + iky(j)
          if(ikxnew <= ikx_max .and. ikxnew >= ikx_min) then
            inew = ikx_to_index(ikxnew)
            if(inew > 0) then
              tmp(k, j, inew) = fld(k, j, i)
            endif
          endif
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fld(k, j, i) = tmp(k, j, i)
        enddo
      enddo
    enddo

    !$acc end data
    !$acc exit data delete(ikx_to_index, tmp)

    deallocate(tmp, ikx_to_index)
  end subroutine remap_fld


!-----------------------------------------------!
!> @author  YK
!! @brief   Transform to non-shearing coordinate
!-----------------------------------------------!
  subroutine to_non_shearing_coordinate(fld)
    use grid, only: xx, dly, nly
    use grid, only: nlx_local, nly, nlz_padded
    use params, only: q
    implicit none
    real(8), dimension(:,:,:), intent(inout) :: fld
    real(8), dimension(:,:,:), allocatable :: tmp
    integer :: i, j, k, dj, j_nsc
    real(8) :: tsc_local
    integer :: shear_flg_local

    !$acc update host(shear_flg, tsc)
    shear_flg_local = shear_flg
    tsc_local = tsc

    allocate(tmp(nlz_padded, nly, nlx_local), source=0.d0)

    !$acc data present(fld, xx) create(tmp)

    !$acc parallel loop gang private(dj) &
    !$acc& copyin(shear_flg_local, tsc_local)
    do i = 1, nlx_local
      dj = floor(q * shear_flg_local * tsc_local * xx(i) / dly)
      !$acc loop vector collapse(2) private(j_nsc)
      do j = 1, nly
        do k = 1, nlz_padded
          j_nsc = j + dj
          if (j_nsc < 1  ) j_nsc = nly + j_nsc
          if (j_nsc > nly) j_nsc = j_nsc - nly
          tmp(k, j, i) = fld(k, j_nsc, i)
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          fld(k, j, i) = tmp(k, j, i)
        enddo
      enddo
    enddo

    !$acc end data

    deallocate(tmp)
  end subroutine to_non_shearing_coordinate

!-----------------------------------------------!
!> @author  YK
!! @brief   Time integral of hyperdissipation
!!          -\int \mathrm{d}t\, q[k_x(t)^2 
!!                 + k_y^2 + k_z^2]^{2n}
!-----------------------------------------------!
  subroutine get_imp_terms_tintg_with_shear(imp_terms_tintg, t, kx, ky, kz, coeff, coeff_h, nexp)
    !$acc routine seq
    use params, only: q
    use grid, only: k2_max
    implicit none
    real(8), intent(out) :: imp_terms_tintg
    real(8), intent(in) :: t, kx, ky, kz, coeff, coeff_h
    integer, intent(in) :: nexp
    real(8) :: regular, hyper
    real(8) :: kxt, k2yz

    if(ky == 0.d0) then
      regular = (kx**2 + ky**2 + kz**2)*t
      hyper   = (kx**2 + ky**2 + kz**2)**nexp*t
    else
      kxt  = kx + q*t*ky
      k2yz = ky**2 + kz**2

      regular = (k2yz*kxt + kxt**3/3.d0)/(q*ky)
      select case(nexp)
      case (1) 
        hyper = (k2yz*kxt + kxt**3/3.d0)/(q*ky)
      case (2) 
        hyper = (k2yz**2*kxt + 2.d0/3.d0*k2yz*kxt**3 + kxt**5/5.d0)/(q*ky)
      case (3) 
        hyper = (k2yz**3*kxt + k2yz**2*kxt**3 + 3.d0/5.d0*k2yz*kxt**5 & 
             + kxt**7/7.d0)/(q*ky)
      case (4) 
        hyper = (k2yz**4*kxt + 4.d0/3.d0*k2yz**3*kxt**3 + 6.d0/5.d0*k2yz**2*kxt**5 & 
             + 4.d0/7.d0*k2yz*kxt**7 + kxt**9/9.d0)/(q*ky)
      end select
    endif
    imp_terms_tintg = -coeff*regular/k2_max - coeff_h*hyper/k2_max**nexp
  end subroutine get_imp_terms_tintg_with_shear

end module shearing_box

