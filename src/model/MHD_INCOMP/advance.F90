!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    25 Sep 2019
!! @brief   Time stepping for MHD_INCOMP
!-----------------------------------------------!
module advance
  use fields, only: nfields
  use fields, only: iux, iuy, iuz
  use fields, only: ibx, iby, ibz
  use mp, only: proc0, max_allreduce, sum_allreduce
  implicit none

  public solve, is_allocated, allocate_advance, deallocate_advance

  logical :: is_allocated = .false.

  integer :: counter = 0
  complex(8), allocatable, dimension(:,:,:)   :: ux_new, uy_new, uz_new
  complex(8), allocatable, dimension(:,:,:)   :: bx_new, by_new, bz_new
  complex(8), allocatable, dimension(:,:,:,:) :: w
  real   (8), allocatable, dimension(:,:,:,:) :: w_r, flx_r 
  complex(8), allocatable, dimension(:,:,:,:) :: flx, exp_terms

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(8), allocatable, dimension(:,:,:)   :: ux_tmp, uy_tmp, uz_tmp
  complex(8), allocatable, dimension(:,:,:)   :: bx_tmp, by_tmp, bz_tmp
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(8), allocatable, dimension(:,:,:)   :: ux_old2, uy_old2, uz_old2
  complex(8), allocatable, dimension(:,:,:)   :: bx_old2, by_old2, bz_old2
  complex(8), allocatable, dimension(:,:,:)   :: fux_old2, fuy_old2, fuz_old2
  complex(8), allocatable, dimension(:,:,:)   :: fbx_old2, fby_old2, fbz_old2
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2
  real   (8), allocatable, dimension(:,:)     :: kxt_old, kxt_old2

  integer :: cfl_unit, rms_unit

  ! Forward FFT variables
  integer, parameter :: nftran = 9
  integer, parameter :: iflx_uxx = 1, iflx_uxy = 2, iflx_uxz = 3 ! uu - bb 
  integer, parameter ::               iflx_uyy = 4, iflx_uyz = 5 ! uu - bb 
  integer, parameter ::                             iflx_uzz = 6 ! uu - bb 
  integer, parameter :: iflx_bx  = 7, iflx_by  = 8, iflx_bz  = 9 ! b x u
  
  !$acc declare copyin(counter)
  !$acc declare copyin(nftran)
  !$acc declare copyin(iflx_uxx, iflx_uxy, iflx_uxz)
  !$acc declare copyin(          iflx_uyy, iflx_uyz)
  !$acc declare copyin(                    iflx_uzz)
  !$acc declare copyin(iflx_bx , iflx_by , iflx_bz )

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use fields, only: ux_old, uy_old, uz_old
    use fields, only: bx_old, by_old, bz_old
    use grid, only: k2_max
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: nu, nu_h, eta_h, nu_h_exp, eta, eta_h_exp, zi, nonlinear, shear, q
    use dealias, only: filter
    use shearing_box, only: shear_flg, kxt, k2t, k2t_inv, tsc, tremap
    use force, only: driven, elsasser, update_force, get_force, normalize_force, normalize_force_els
    use force, only: fux, fuz, fuy, fux_old, fuy_old, fuz_old
    use force, only: fbx, fbz, fby, fbx_old, fby_old, fbz_old
    use force, only: fzpx, fzpy, fzpz, fzmx, fzmy, fzmz
    use advance_common, only: gear1, gear2, gear3
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use diagnostics_common, only: series_output
    use params, only: time_step_scheme
    implicit none
    real(8) :: imp_terms_tintg0(nfields), imp_terms_tintg1(nfields)
    real(8) :: imp_terms_tintg2(nfields), imp_terms_tintg3(nfields)  
    integer :: i, j, k, l

    if (proc0) call put_time_stamp(timer_advance)

    ! initialize tmp fields
    if(counter == 0) then
      call init_work_fields
      counter = 1
      !$acc update device (counter)
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      ! Calculate force terms
      if (driven) then
        ! at n
        if(elsasser) then
          call get_force('zpx', fzpx)
          call get_force('zpy', fzpy)
          call get_force('zpz', fzpz)
          call get_force('zmx', fzmx)
          call get_force('zmy', fzmy)
          call get_force('zmz', fzmz)
          call div_free_force(fzpx, fzpy, fzpz)
          call div_free_force(fzmx, fzmy, fzmz)
          call normalize_force_els(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          !$acc data present(fux , fuy , fuz , fbx , fby , fbz, &
          !$acc              fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          !$acc parallel loop collapse(3)
          do i = 1, nkx
            do j = 1, nky_local
              do k = 1, nkz
                fux(k,j,i) = fzpx(k,j,i) + fzmx(k,j,i)
                fuy(k,j,i) = fzpy(k,j,i) + fzmy(k,j,i)
                fuz(k,j,i) = fzpz(k,j,i) + fzmz(k,j,i)
                fbx(k,j,i) = fzpx(k,j,i) - fzmx(k,j,i)
                fby(k,j,i) = fzpy(k,j,i) - fzmy(k,j,i)
                fbz(k,j,i) = fzpz(k,j,i) - fzmz(k,j,i)
              enddo
            enddo
          enddo
          !$acc end data
        else
          call get_force('ux', fux)
          call get_force('uy', fuy)
          call get_force('uz', fuz)
          call get_force('bx', fbx)
          call get_force('by', fby)
          call get_force('bz', fbz)
          call div_free_force(fux, fuy, fuz)
          call div_free_force(fbx, fby, fbz)
          call normalize_force(fux, fuy, fuz, fbx, fby, fbz)
        endif
      endif

      !---------------  RK 1st step  ---------------
      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(ux, uy, uz, bx, by, bz, .true.)

      !$acc data present(exp_terms, flx, kxt, ky, kz, k2t_inv, &
      !$acc              ux    , uy    , uz    , bx    , by    , bz    , &
      !$acc              ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, &
      !$acc              fux   , fuy   , fuz   , fbx   , fby   , fbz   , &
      !$acc              kx, ky, kz)
      !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg1)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iux), exp_terms(k,j,i,iuy), exp_terms(k,j,i,iuz), &
                               exp_terms(k,j,i,ibx), exp_terms(k,j,i,iby), exp_terms(k,j,i,ibz), &
                               ux(k,j,i), uy(k,j,i), uz(k,j,i), bx(k,j,i), by(k,j,i), bz(k,j,i), &
                               flx(k,j,i,iflx_uxx), flx(k,j,i,iflx_uxy), flx(k,j,i,iflx_uxz), &
                                                    flx(k,j,i,iflx_uyy), flx(k,j,i,iflx_uyz), &
                                                                         flx(k,j,i,iflx_uzz), &
                               flx(k,j,i,iflx_bx ), flx(k,j,i,iflx_by ), flx(k,j,i,iflx_bz ), &
                               fux(k,j,i), fuy(k,j,i), fuz(k,j,i), &
                               fbx(k,j,i), fby(k,j,i), fbz(k,j,i), &
                               kxt(i,j), ky(j), kz(k), k2t_inv(k,j,i))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iux), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iux), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuy), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iuy), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuz), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iuz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ibx), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ibx), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(iby), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ibz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)

            ! update u
            call eSSPIFRK1(ux_tmp(k,j,i), ux(k,j,i), &
               exp_terms(k,j,i,iux), &
               imp_terms_tintg1(iux), imp_terms_tintg0(iux) &
            )
            call eSSPIFRK1(uy_tmp(k,j,i), uy(k,j,i), &
               exp_terms(k,j,i,iuy), &
               imp_terms_tintg1(iuy), imp_terms_tintg0(iuy) &
            )
            call eSSPIFRK1(uz_tmp(k,j,i), uz(k,j,i), &
               exp_terms(k,j,i,iuz), &
               imp_terms_tintg1(iuz), imp_terms_tintg0(iuz) &
            )

            ! update b
            call eSSPIFRK1(bx_tmp(k,j,i), bx(k,j,i), &
               exp_terms(k,j,i,ibx), &
               imp_terms_tintg1(ibx), imp_terms_tintg0(ibx) &
            )
            call eSSPIFRK1(by_tmp(k,j,i), by(k,j,i), &
               exp_terms(k,j,i,iby), &
               imp_terms_tintg1(iby), imp_terms_tintg0(iby) &
            )
            call eSSPIFRK1(bz_tmp(k,j,i), bz(k,j,i), &
               exp_terms(k,j,i,ibz), &
               imp_terms_tintg1(ibz), imp_terms_tintg0(ibz) &
            )

          enddo
        enddo
      enddo
      !$acc end data

      ! save explicit terms at the previous step
      !$acc data present(exp_terms, exp_terms0)
      !$acc parallel loop collapse(4)
      do l = 1, nfields
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              exp_terms0(k,j,i,l) = exp_terms(k,j,i,l)
            enddo
          enddo
        enddo
      enddo
      !$acc end data

      ! Dealiasing
      !$acc data present(ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, filter)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            ux_tmp(k,j,i) = ux_tmp(k,j,i)*filter(k,j,i)
            uy_tmp(k,j,i) = uy_tmp(k,j,i)*filter(k,j,i)
            uz_tmp(k,j,i) = uz_tmp(k,j,i)*filter(k,j,i)
            bx_tmp(k,j,i) = bx_tmp(k,j,i)*filter(k,j,i)
            by_tmp(k,j,i) = by_tmp(k,j,i)*filter(k,j,i)
            bz_tmp(k,j,i) = bz_tmp(k,j,i)*filter(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data

      !---------------  RK 2nd step  ---------------
      ! Calculate force terms
      if (driven) then
        ! at n + 2/3
        call update_force(2.d0/3.d0*dt)

        if(elsasser) then
          call get_force('zpx', fzpx)
          call get_force('zpy', fzpy)
          call get_force('zpz', fzpz)
          call get_force('zmx', fzmx)
          call get_force('zmy', fzmy)
          call get_force('zmz', fzmz)
          call div_free_force(fzpx, fzpy, fzpz)
          call div_free_force(fzmx, fzmy, fzmz)
          call normalize_force_els(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          !$acc data present(fux , fuy , fuz , fbx , fby , fbz, &
          !$acc              fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          !$acc parallel loop collapse(3)
          do i = 1, nkx
            do j = 1, nky_local
              do k = 1, nkz
                fux(k,j,i) = fzpx(k,j,i) + fzmx(k,j,i)
                fuy(k,j,i) = fzpy(k,j,i) + fzmy(k,j,i)
                fuz(k,j,i) = fzpz(k,j,i) + fzmz(k,j,i)
                fbx(k,j,i) = fzpx(k,j,i) - fzmx(k,j,i)
                fby(k,j,i) = fzpy(k,j,i) - fzmy(k,j,i)
                fbz(k,j,i) = fzpz(k,j,i) - fzmz(k,j,i)
              enddo
            enddo
          enddo
          !$acc end data
        else
          call get_force('ux', fux)
          call get_force('uy', fuy)
          call get_force('uz', fuz)
          call get_force('bx', fbx)
          call get_force('by', fby)
          call get_force('bz', fbz)
          call div_free_force(fux, fuy, fuz)
          call div_free_force(fbx, fby, fbz)
          call normalize_force(fux, fuy, fuz, fbx, fby, fbz)
        endif

        ! go to n + 1
        call update_force(1.d0/3.d0*dt)
      endif

      ! Calculate kxt at n + 2/3
      if(shear) then
        !$acc parallel loop collapse(2) private(k)
        do j = 1, nky_local
          do i = 1, nkx
            kxt(i,j) = kx(i) + q*shear_flg*(tsc + 2.d0/3.d0*dt)*ky(j)

            !$acc loop
            do k = 1, nkz
              k2t(k,j,i) = kxt(i,j)**2 + ky(j)**2 + kz(k)**2

              if(k2t(k,j,i) == 0.d0) then
                k2t_inv(k,j,i) = 0.d0
              else
                k2t_inv(k,j,i) = 1.d0/k2t(k,j,i)
              endif
            enddo
          enddo
        enddo
        !$acc end parallel
      endif

      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, .false.)

      !$acc data present(exp_terms, flx, kxt, ky, kz, k2t_inv, &
      !$acc              ux    , uy    , uz    , bx    , by    , bz    , &
      !$acc              ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, &
      !$acc              fux   , fuy   , fuz   , fbx   , fby   , fbz   , &
      !$acc              kx, ky, kz)
      !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg2)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iux), exp_terms(k,j,i,iuy), exp_terms(k,j,i,iuz), &
                               exp_terms(k,j,i,ibx), exp_terms(k,j,i,iby), exp_terms(k,j,i,ibz), &
                               ux_tmp(k,j,i), uy_tmp(k,j,i), uz_tmp(k,j,i), bx_tmp(k,j,i), by_tmp(k,j,i), bz_tmp(k,j,i), &
                               flx(k,j,i,iflx_uxx), flx(k,j,i,iflx_uxy), flx(k,j,i,iflx_uxz), &
                                                    flx(k,j,i,iflx_uyy), flx(k,j,i,iflx_uyz), &
                                                                         flx(k,j,i,iflx_uzz), &
                               flx(k,j,i,iflx_bx ), flx(k,j,i,iflx_by ), flx(k,j,i,iflx_bz ), &
                               fux(k,j,i), fuy(k,j,i), fuz(k,j,i), &
                               fbx(k,j,i), fby(k,j,i), fbz(k,j,i), &
                               kxt(i,j), ky(j), kz(k), k2t_inv(k,j,i))
            !
            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iux), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iux), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuy), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iuy), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuz), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iuz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ibx), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)

            ! update u
            call eSSPIFRK2(ux_tmp(k,j,i), ux_tmp(k,j,i), ux(k,j,i), &
               exp_terms(k,j,i,iux), &
               imp_terms_tintg2(iux), imp_terms_tintg0(iux) &
            )
            call eSSPIFRK2(uy_tmp(k,j,i), uy_tmp(k,j,i), uy(k,j,i), &
               exp_terms(k,j,i,iuy), &
               imp_terms_tintg2(iuy), imp_terms_tintg0(iuy) &
            )
            call eSSPIFRK2(uz_tmp(k,j,i), uz_tmp(k,j,i), uz(k,j,i), &
               exp_terms(k,j,i,iuz), &
               imp_terms_tintg2(iuz), imp_terms_tintg0(iuz) &
            )

            ! update b
            call eSSPIFRK2(bx_tmp(k,j,i), bx_tmp(k,j,i), bx(k,j,i), &
               exp_terms(k,j,i,ibx), &
               imp_terms_tintg2(ibx), imp_terms_tintg0(ibx) &
            )
            call eSSPIFRK2(by_tmp(k,j,i), by_tmp(k,j,i), by(k,j,i), &
               exp_terms(k,j,i,iby), &
               imp_terms_tintg2(iby), imp_terms_tintg0(iby) &
            )
            call eSSPIFRK2(bz_tmp(k,j,i), bz_tmp(k,j,i), bz(k,j,i), &
               exp_terms(k,j,i,ibz), &
               imp_terms_tintg2(ibz), imp_terms_tintg0(ibz) &
            )
          enddo
        enddo
      enddo
      !$acc end data

      ! Dealiasing
      !$acc data present(ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, filter)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            ux_tmp(k,j,i) = ux_tmp(k,j,i)*filter(k,j,i)
            uy_tmp(k,j,i) = uy_tmp(k,j,i)*filter(k,j,i)
            uz_tmp(k,j,i) = uz_tmp(k,j,i)*filter(k,j,i)
            bx_tmp(k,j,i) = bx_tmp(k,j,i)*filter(k,j,i)
            by_tmp(k,j,i) = by_tmp(k,j,i)*filter(k,j,i)
            bz_tmp(k,j,i) = bz_tmp(k,j,i)*filter(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data

      !---------------  RK 3rd step  ---------------
      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, .false.)

      !$acc data present(exp_terms, exp_terms0, flx, kxt, ky, kz, k2t_inv, &
      !$acc              ux    , uy    , uz    , bx    , by    , bz    , &
      !$acc              ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, &
      !$acc              ux_new, uy_new, uz_new, bx_new, by_new, bz_new, &
      !$acc              fux   , fuy   , fuz   , fbx   , fby   , fbz   , &
      !$acc              kx, ky, kz)
      !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg2, imp_terms_tintg3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iux), exp_terms(k,j,i,iuy), exp_terms(k,j,i,iuz), &
                               exp_terms(k,j,i,ibx), exp_terms(k,j,i,iby), exp_terms(k,j,i,ibz), &
                               ux_tmp(k,j,i), uy_tmp(k,j,i), uz_tmp(k,j,i), bx_tmp(k,j,i), by_tmp(k,j,i), bz_tmp(k,j,i), &
                               flx(k,j,i,iflx_uxx), flx(k,j,i,iflx_uxy), flx(k,j,i,iflx_uxz), &
                                                    flx(k,j,i,iflx_uyy), flx(k,j,i,iflx_uyz), &
                                                                         flx(k,j,i,iflx_uzz), &
                               flx(k,j,i,iflx_bx ), flx(k,j,i,iflx_by ), flx(k,j,i,iflx_bz ), &
                               fux(k,j,i), fuy(k,j,i), fuz(k,j,i), &
                               fbx(k,j,i), fby(k,j,i), fbz(k,j,i), &
                               kxt(i,j), ky(j), kz(k), k2t_inv(k,j,i))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iux), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iux), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iux), tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuy), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iuy), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iuy), tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuz), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iuz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iuz), tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ibx), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibx), tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(iby), tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibz), tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)

            ! update u
            call eSSPIFRK3(ux_new(k,j,i), ux_tmp(k,j,i), ux(k,j,i), &
               exp_terms(k,j,i,iux), exp_terms0(k,j,i,iux), &
               imp_terms_tintg3(iux), imp_terms_tintg2(iux), imp_terms_tintg0(iux) &
            )
            call eSSPIFRK3(uy_new(k,j,i), uy_tmp(k,j,i), uy(k,j,i), &
               exp_terms(k,j,i,iuy), exp_terms0(k,j,i,iuy), &
               imp_terms_tintg3(iuy), imp_terms_tintg2(iuy), imp_terms_tintg0(iuy) &
            )
            call eSSPIFRK3(uz_new(k,j,i), uz_tmp(k,j,i), uz(k,j,i), &
               exp_terms(k,j,i,iuz), exp_terms0(k,j,i,iuz), &
               imp_terms_tintg3(iuz), imp_terms_tintg2(iuz), imp_terms_tintg0(iuz) &
            )

            ! update b
            call eSSPIFRK3(bx_new(k,j,i), bx_tmp(k,j,i), bx(k,j,i), &
               exp_terms(k,j,i,ibx), exp_terms0(k,j,i,ibx), &
               imp_terms_tintg3(ibx), imp_terms_tintg2(ibx), imp_terms_tintg0(ibx) &
            )
            call eSSPIFRK3(by_new(k,j,i), by_tmp(k,j,i), by(k,j,i), &
               exp_terms (k,j,i,iby), exp_terms0(k,j,i,iby), &
               imp_terms_tintg3(iby), imp_terms_tintg2(iby), imp_terms_tintg0(iby) &
            )

            call eSSPIFRK3(bz_new(k,j,i), bz_tmp(k,j,i), bz(k,j,i), &
               exp_terms(k,j,i,ibz), exp_terms0(k,j,i,ibz), &
               imp_terms_tintg3(ibz), imp_terms_tintg2(ibz), imp_terms_tintg0(ibz) &
            )
          enddo
        enddo
      enddo
      !$acc end data

      ! save fields at the previous step
      !$acc data present(ux, uy, uz, bx, by, bz, ux_old, uy_old, uz_old, bx_old, by_old, bz_old)
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
      !$acc end data

      if (driven) then
        !$acc data present(fux, fuy, fuz, fbx, fby, fbz, fux_old, fuy_old, fuz_old, fbx_old, fby_old, fbz_old)
        !$acc parallel loop collapse(3)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              fux_old(k,j,i) = fux(k,j,i)
              fuy_old(k,j,i) = fuy(k,j,i)
              fuz_old(k,j,i) = fuz(k,j,i)
              fbx_old(k,j,i) = fbx(k,j,i)
              fby_old(k,j,i) = fby(k,j,i)
              fbz_old(k,j,i) = fbz(k,j,i)
            enddo
          enddo
        enddo
        !$acc end data
      endif

      ! Dealiasing
      !$acc data present(ux_new, uy_new, uz_new, bx_new, by_new, bz_new, filter)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            ux_new(k,j,i) = ux_new(k,j,i)*filter(k,j,i)
            uy_new(k,j,i) = uy_new(k,j,i)*filter(k,j,i)
            uz_new(k,j,i) = uz_new(k,j,i)*filter(k,j,i)
            bx_new(k,j,i) = bx_new(k,j,i)*filter(k,j,i)
            by_new(k,j,i) = by_new(k,j,i)*filter(k,j,i)
            bz_new(k,j,i) = bz_new(k,j,i)*filter(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      ! Calculate force terms
      if (driven) then
        call update_force(dt)

        if(elsasser) then
          call get_force('zpx', fzpx)
          call get_force('zpy', fzpy)
          call get_force('zpz', fzpz)
          call get_force('zmx', fzmx)
          call get_force('zmy', fzmy)
          call get_force('zmz', fzmz)
          call div_free_force(fzpx, fzpy, fzpz)
          call div_free_force(fzmx, fzmy, fzmz)
          call normalize_force_els(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          !$acc data present(fux , fuy , fuz , fbx , fby , fbz, &
          !$acc              fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          !$acc parallel loop collapse(3)
          do i = 1, nkx
            do j = 1, nky_local
              do k = 1, nkz
                fux(k,j,i) = fzpx(k,j,i) + fzmx(k,j,i)
                fuy(k,j,i) = fzpy(k,j,i) + fzmy(k,j,i)
                fuz(k,j,i) = fzpz(k,j,i) + fzmz(k,j,i)
                fbx(k,j,i) = fzpx(k,j,i) - fzmx(k,j,i)
                fby(k,j,i) = fzpy(k,j,i) - fzmy(k,j,i)
                fbz(k,j,i) = fzpz(k,j,i) - fzmz(k,j,i)
              enddo
            enddo
          enddo
          !$acc end data
        else
          call get_force('ux', fux)
          call get_force('uy', fuy)
          call get_force('uz', fuz)
          call get_force('bx', fbx)
          call get_force('by', fby)
          call get_force('bz', fbz)
          call div_free(fux, fuy, fuz)
          call div_free(fbx, fby, fbz)
          call normalize_force(fux, fuy, fuz, fbx, fby, fbz)
        endif
      endif

      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(ux, uy, uz, bx, by, bz, .true.)

      !$acc data present(exp_terms, exp_terms_old, exp_terms_old2, flx, kx, ky, kz, k2t, k2_max, k2t_inv, &
      !$acc              ux      , uy      , uz      , bx      , by      , bz      , &
      !$acc              ux_old  , uy_old  , uz_old  , bx_old  , by_old  , bz_old  , &
      !$acc              ux_old2 , uy_old2 , uz_old2 , bx_old2 , by_old2 , bz_old2 , &
      !$acc              ux_new  , uy_new  , uz_new  , bx_new  , by_new  , bz_new  , &
      !$acc              fux     , fuy     , fuz     , fbx     , fby     , fbz     , &
      !$acc              fux_old , fuy_old , fuz_old , fbx_old , fby_old , fbz_old , &
      !$acc              fux_old2, fuy_old2, fuz_old2, fbx_old2, fby_old2, fbz_old2, &
      !$acc              kxt, kxt_old, kxt_old2)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iux), exp_terms(k,j,i,iuy), exp_terms(k,j,i,iuz), &
                               exp_terms(k,j,i,ibx), exp_terms(k,j,i,iby), exp_terms(k,j,i,ibz), &
                               ux(k,j,i), uy(k,j,i), uz(k,j,i), bx(k,j,i), by(k,j,i), bz(k,j,i), &
                               flx(k,j,i,iflx_uxx), flx(k,j,i,iflx_uxy), flx(k,j,i,iflx_uxz), &
                                                    flx(k,j,i,iflx_uyy), flx(k,j,i,iflx_uyz), &
                                                                         flx(k,j,i,iflx_uzz), &
                               flx(k,j,i,iflx_bx ), flx(k,j,i,iflx_by ), flx(k,j,i,iflx_bz ), &
                               fux(k,j,i), fuy(k,j,i), fuz(k,j,i), &
                               fbx(k,j,i), fby(k,j,i), fbz(k,j,i), &
                               kxt(i,j), ky(j), kz(k), k2t_inv(k,j,i))

            ! 1st order 
            if(counter == 1) then
              ! update u
              call gear1(ux_new(k,j,i), ux(k,j,i), &
                 exp_terms(k,j,i,iux), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )
              call gear1(uy_new(k,j,i), uy(k,j,i), &
                 exp_terms(k,j,i,iuy), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )
              call gear1(uz_new(k,j,i), uz(k,j,i), &
                 exp_terms(k,j,i,iuz), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )

              ! update b
              call gear1(bx_new(k,j,i), bx(k,j,i), &
                 exp_terms(k,j,i,ibx), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )
              call gear1(by_new(k,j,i), by(k,j,i), &
                 exp_terms(k,j,i,iby), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )
              call gear1(bz_new(k,j,i), bz(k,j,i), &
                 exp_terms(k,j,i,ibz), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )

            ! 2nd order 
            elseif(counter == 2) then
              ! update u
              call gear2(ux_new(k,j,i), ux(k,j,i), ux_old(k,j,i), &
                 exp_terms    (k,j,i,iux), &
                 exp_terms_old(k,j,i,iux), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )
              call gear2(uy_new(k,j,i), uy(k,j,i), uy_old(k,j,i), &
                 exp_terms    (k,j,i,iuy), &
                 exp_terms_old(k,j,i,iuy), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )
              call gear2(uz_new(k,j,i), uz(k,j,i), uz_old(k,j,i), &
                 exp_terms    (k,j,i,iuz), &
                 exp_terms_old(k,j,i,iuz), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )

              ! update b
              call gear2(bx_new(k,j,i), bx(k,j,i), bx_old(k,j,i), &
                 exp_terms    (k,j,i,ibx), &
                 exp_terms_old(k,j,i,ibx), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )
              call gear2(by_new(k,j,i), by(k,j,i), by_old(k,j,i), &
                 exp_terms    (k,j,i,iby), &
                 exp_terms_old(k,j,i,iby), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )
              call gear2(bz_new(k,j,i), bz(k,j,i), bz_old(k,j,i), &
                 exp_terms    (k,j,i,ibz), &
                 exp_terms_old(k,j,i,ibz), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )

            ! 3rd order 
            else
              ! update u
              call gear3(ux_new(k,j,i), ux(k,j,i), ux_old(k,j,i), ux_old2(k,j,i), &
                 exp_terms     (k,j,i,iux), &
                 exp_terms_old (k,j,i,iux), &
                 exp_terms_old2(k,j,i,iux), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )
              call gear3(uy_new(k,j,i), uy(k,j,i), uy_old(k,j,i), uy_old2(k,j,i), &
                 exp_terms     (k,j,i,iuy), &
                 exp_terms_old (k,j,i,iuy), &
                 exp_terms_old2(k,j,i,iuy), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )
              call gear3(uz_new(k,j,i), uz(k,j,i), uz_old(k,j,i), uz_old2(k,j,i), &
                 exp_terms     (k,j,i,iuz), &
                 exp_terms_old (k,j,i,iuz), &
                 exp_terms_old2(k,j,i,iuz), &
                 nu*(k2t(k,j,i)/k2_max) + nu_h*(k2t(k,j,i)/k2_max)**nu_h_exp &
              )

              ! update b
              call gear3(bx_new(k,j,i), bx(k,j,i), bx_old(k,j,i), bx_old2(k,j,i), &
                 exp_terms     (k,j,i,ibx), &
                 exp_terms_old (k,j,i,ibx), &
                 exp_terms_old2(k,j,i,ibx), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )
              call gear3(by_new(k,j,i), by(k,j,i), by_old(k,j,i), by_old2(k,j,i), &
                 exp_terms     (k,j,i,iby), &
                 exp_terms_old (k,j,i,iby), &
                 exp_terms_old2(k,j,i,iby), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )
              call gear3(bz_new(k,j,i), bz(k,j,i), bz_old(k,j,i), bz_old2(k,j,i), &
                 exp_terms     (k,j,i,ibz), &
                 exp_terms_old (k,j,i,ibz), &
                 exp_terms_old2(k,j,i,ibz), &
                 eta*(k2t(k,j,i)/k2_max) + eta_h*(k2t(k,j,i)/k2_max)**eta_h_exp &
              )
            endif
          enddo
        enddo
      enddo

      if(counter <= 2) counter = counter + 1
      !$acc update device (counter)

      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            ux_old2(k,j,i) = ux_old(k,j,i) 
            ux_old(k,j,i)  = ux(k,j,i)

            uy_old2(k,j,i) = uy_old(k,j,i) 
            uy_old(k,j,i)  = uy(k,j,i)

            uz_old2(k,j,i) = uz_old(k,j,i) 
            uz_old(k,j,i)  = uz(k,j,i)

            bx_old2(k,j,i) = bx_old(k,j,i) 
            bx_old(k,j,i)  = bx(k,j,i)

            by_old2(k,j,i) = by_old(k,j,i) 
            by_old(k,j,i)  = by(k,j,i)

            bz_old2(k,j,i) = bz_old(k,j,i) 
            bz_old(k,j,i)  = bz(k,j,i)

            do l = 1, nfields
              exp_terms_old2(k,j,i,l) = exp_terms_old(k,j,i,l)
              exp_terms_old(k,j,i,l)  = exp_terms(k,j,i,l)
            enddo
          enddo
        enddo
      enddo

      if (driven) then
        !$acc parallel loop collapse(3)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              fux_old2(k,j,i) = fux_old(k,j,i) 
              fux_old(k,j,i)  = fux(k,j,i)

              fuy_old2(k,j,i) = fuy_old(k,j,i) 
              fuy_old(k,j,i)  = fuy(k,j,i)

              fuz_old2(k,j,i) = fuz_old(k,j,i) 
              fuz_old(k,j,i)  = fuz(k,j,i)

              fbx_old2(k,j,i) = fbx_old(k,j,i) 
              fbx_old(k,j,i)  = fbx(k,j,i)

              fby_old2(k,j,i) = fby_old(k,j,i) 
              fby_old(k,j,i)  = fby(k,j,i)

              fbz_old2(k,j,i) = fbz_old(k,j,i) 
              fbz_old(k,j,i)  = fbz(k,j,i)
            enddo
          enddo
        enddo
      endif

      if(shear) then
        !$acc parallel loop collapse(2)
        do j = 1, nky_local
          do i = 1, nkx
            kxt_old2(i,j) = kxt_old(i,j) 
            kxt_old(i,j)  = kxt(i,j)
          enddo
        enddo
      endif

      ! Dealiasing
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            ux_new(k,j,i) = ux_new(k,j,i)*filter(k,j,i)
            uy_new(k,j,i) = uy_new(k,j,i)*filter(k,j,i)
            uz_new(k,j,i) = uz_new(k,j,i)*filter(k,j,i)
            bx_new(k,j,i) = bx_new(k,j,i)*filter(k,j,i)
            by_new(k,j,i) = by_new(k,j,i)*filter(k,j,i)
            bz_new(k,j,i) = bz_new(k,j,i)*filter(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data
    endif

    !$acc data present(ux, uy, uz, bx, by, bz, ux_new, uy_new, uz_new, bx_new, by_new, bz_new)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux(k,j,i) = ux_new(k,j,i)
          uy(k,j,i) = uy_new(k,j,i)
          uz(k,j,i) = uz_new(k,j,i)
          bx(k,j,i) = bx_new(k,j,i)
          by(k,j,i) = by_new(k,j,i)
          bz(k,j,i) = bz_new(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

    tt  = tt  + dt
    tsc = tsc + dt
    !$acc update device (tt, tsc)

    if(shear .and. tsc > tremap) call remap

    if(shear) then
      !$acc parallel loop collapse(2) private(k)
      do j = 1, nky_local
        do i = 1, nkx
          kxt(i,j) = kx(i) + q*shear_flg*tsc*ky(j)

          !$acc loop
          do k = 1, nkz
            k2t(k,j,i) = kxt(i,j)**2 + ky(j)**2 + kz(k)**2

            if(k2t(k,j,i) == 0.d0) then
              k2t_inv(k,j,i) = 0.d0
            else
              k2t_inv(k,j,i) = 1.d0/k2t(k,j,i)
            endif
          enddo
        enddo
      enddo
      !$acc end parallel
    endif

    ! Div u & b cleaing
    call div_free(ux, uy, uz)
    call div_free(bx, by, bz)

    if(series_output) call output_series_modes

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Allocate fields used only here
!-----------------------------------------------!
  subroutine init_work_fields
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: q
    use shearing_box, only: shear_flg, tsc, kxt, k2t, k2t_inv
    use file, only: open_output_file
    use params, only: time_step_scheme
    use mp, only: proc0
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src
    integer :: i, j, k

    !$acc parallel loop collapse(2) private(k)
    do j = 1, nky_local
      do i = 1, nkx
        kxt(i, j) = kx(i) + q*shear_flg*tsc*ky(j)

        !$acc loop
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
    !$acc end parallel

    call allocate_advance

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      allocate(src(nkz, nky_local, nkx), source=(0.d0,0.d0))

      allocate( ux_old2, source=src)
      allocate( uy_old2, source=src)
      allocate( uz_old2, source=src)
      allocate( bx_old2, source=src)
      allocate( by_old2, source=src)
      allocate( bz_old2, source=src)

      allocate(fux_old2, source=src)
      allocate(fuy_old2, source=src)
      allocate(fuz_old2, source=src)
      allocate(fbx_old2, source=src)
      allocate(fby_old2, source=src)
      allocate(fbz_old2, source=src)

      allocate(exp_terms_old (nkz, nky_local, nkx, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(nkz, nky_local, nkx, nfields)); exp_terms_old2 = 0.d0

      allocate(kxt_old (1:nkx, 1:nky_local)); kxt_old  = kxt
      allocate(kxt_old2(1:nkx, 1:nky_local)); kxt_old2 = kxt

      deallocate(src)

      !$acc enter data copyin( ux_old2)
      !$acc enter data copyin( uy_old2)
      !$acc enter data copyin( uz_old2)
      !$acc enter data copyin( bx_old2)
      !$acc enter data copyin( by_old2)
      !$acc enter data copyin( bz_old2)

      !$acc enter data copyin(fux_old2)
      !$acc enter data copyin(fuy_old2)
      !$acc enter data copyin(fuz_old2)
      !$acc enter data copyin(fbx_old2)
      !$acc enter data copyin(fby_old2)
      !$acc enter data copyin(fbz_old2)
      
      !$acc enter data copyin(exp_terms_old )
      !$acc enter data copyin(exp_terms_old2)
      !$acc enter data copyin(kxt_old )
      !$acc enter data copyin(kxt_old2)
    endif


    if(proc0) call open_output_file (cfl_unit, 'cfl.dat')
    if(proc0) call open_output_file (rms_unit, 'rms.dat')

  end subroutine init_work_fields


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Calculate nonlinear terms via
!!          1. Calculate grad in Fourier space
!!          2. Inverse FFT
!!          3. Calculate nonlinear terms 
!!             in real space
!!          4. Forward FFT
!-----------------------------------------------!
  subroutine get_nonlinear_terms(ux, uy, uz, bx, by, bz, dt_reset)
    use grid, only: dlx, dly, dlz
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz
    use grid, only: ntot
    use mp, only: proc0, max_allreduce, sum_allreduce
    use params, only: zi
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms
    use advance_common, only: dt_adjust_while_running 
    use cuFFTmp, only: btran_c2r, ftran_r2c
    implicit none
    complex(8), dimension (:,:,:), intent(in) :: ux, uy, uz, bx, by, bz

    logical, intent(in) :: dt_reset

    integer :: i, j, k, l
    real   (8) :: ux_rms, uy_rms, uz_rms, bx_rms, by_rms, bz_rms
    real   (8) :: max_vel_x, max_vel_y, max_vel_z , dt_cfl, dt_digit
    !$acc declare create(ux_rms, uy_rms, uz_rms, bx_rms, by_rms, bz_rms)

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    !$acc data present(ux, uy, uz, bx, by, bz, w)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          w(k,j,i,iux) = ux(k,j,i)
          w(k,j,i,iuy) = uy(k,j,i)
          w(k,j,i,iuz) = uz(k,j,i)
          w(k,j,i,ibx) = bx(k,j,i)
          w(k,j,i,iby) = by(k,j,i)
          w(k,j,i,ibz) = bz(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

    ! 1. Inverse FFT
    do i = 1, nfields
      call btran_c2r(w(:,:,:,i), w_r(:,:,:,i))
    enddo


    ! write rms
    ux_rms = 0.d0; uy_rms = 0.d0; uz_rms = 0.d0
    bx_rms = 0.d0; by_rms = 0.d0; bz_rms = 0.d0
    !$acc data present(w_r)
    !$acc parallel loop collapse(3) reduction(+:ux_rms, uy_rms, uz_rms, bx_rms, by_rms, bz_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          ux_rms = ux_rms + sqrt(w_r(k,j,i,iux)**2)
          uy_rms = uy_rms + sqrt(w_r(k,j,i,iuy)**2)
          uz_rms = uz_rms + sqrt(w_r(k,j,i,iuz)**2)
          bx_rms = bx_rms + sqrt(w_r(k,j,i,ibx)**2)
          by_rms = by_rms + sqrt(w_r(k,j,i,iby)**2)
          bz_rms = bz_rms + sqrt(w_r(k,j,i,ibz)**2)
        end do
      end do
    end do
    !$acc end data

    !$acc update host (ux_rms, uy_rms, uz_rms, bx_rms, by_rms, bz_rms)

    call sum_allreduce(ux_rms); call sum_allreduce(uy_rms); call sum_allreduce(uz_rms)
    call sum_allreduce(bx_rms); call sum_allreduce(by_rms); call sum_allreduce(bz_rms)

    ux_rms = sqrt(ux_rms/ntot); uy_rms = sqrt(uy_rms/ntot); uz_rms = sqrt(uz_rms/ntot)
    bx_rms = sqrt(bx_rms/ntot); by_rms = sqrt(by_rms/ntot); bz_rms = sqrt(bz_rms/ntot)

    if(proc0) then
      write (unit=rms_unit, fmt="(100es30.21)") tt, ux_rms, uy_rms, uz_rms, bx_rms, by_rms, bz_rms
      call flush(rms_unit) 
    endif


    ! (get max_vel for dt reset)
    if(dt_reset) then
      max_vel_x = 0.d0; max_vel_y = 0.d0; max_vel_z = 0.d0
      !$acc data present(w_r)
      !$acc parallel loop collapse(3) reduction(max:max_vel_x, max_vel_y, max_vel_z)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz_padded
            max_vel_x = max( &
                          max_vel_x, &
                          abs(w_r(k,j,i,iux) + w_r(k,j,i,ibx)), &
                          abs(w_r(k,j,i,iux) - w_r(k,j,i,ibx))  &
                        )
            max_vel_y = max( &
                          max_vel_y, &
                          abs(w_r(k,j,i,iuy) + w_r(k,j,i,iby)), &
                          abs(w_r(k,j,i,iuy) - w_r(k,j,i,iby))  &
                        )
            max_vel_z = max( &
                          max_vel_z, &
                          abs(w_r(k,j,i,iuz) + w_r(k,j,i,ibz)), &
                          abs(w_r(k,j,i,iuz) - w_r(k,j,i,ibz))  &
                        )
          end do
        end do
      end do
      !$acc end data

      call max_allreduce(max_vel_x)
      call max_allreduce(max_vel_y)
      call max_allreduce(max_vel_z)
      dt_cfl = cfl*min(dlx/max_vel_x, dly/max_vel_y, dlz/max_vel_z)

      if(proc0) then
        write (unit=cfl_unit, fmt="(100es30.21)") tt, dt_cfl, max_vel_x, max_vel_y, max_vel_z
        call flush(cfl_unit) 
      endif

      if(dt_cfl < dt) then
        if(proc0) then
          print *
          write (*, '("dt is decreased from ", es12.4e3)', advance='no') dt
        endif

        dt_digit = (log10(dt_cfl)/abs(log10(dt_cfl)))*ceiling(abs(log10(dt_cfl)))
        ! dt = floor(dt_cfl*10.d0**(-dt_digit))*10.d0**dt_digit

        if (reset_method == 'multiply') then
          dt = 0.5d0*dt
        elseif (reset_method == 'decrement') then
          dt_digit = (log10(dt)/abs(log10(dt)))*ceiling(abs(log10(dt)))

          ! when dt = 0.0**01***
          if (dt*10.d0**(-dt_digit) - 1.0d0 < 1.0d0) then
            dt = 0.9d0*10.d0**dt_digit
          else
            dt = (dt*10.d0**(-dt_digit) - 1.0d0)*10.d0**dt_digit
          endif
        endif

        counter = 1

        if(proc0) then
          print '("  to ", es12.4e3)', dt
          print *
        endif

        !$acc update device (dt, counter)
      endif
      if(dt_cfl > 0.d0 .and. dt_cfl > increase_dt .and. dt < increase_dt) then
        if(proc0) then
          print *
          write (*, '("dt is increased from ", es12.4e3)', advance='no') dt
        endif

        dt = increase_dt

        counter = 1

        if(proc0) then
          print '("  to ", es12.4e3)', dt
          print *
        endif
      endif

      !$acc update device (dt, counter)
    endif

    ! When the file 'dt_adjust' including a float number is created,
    ! dt will be manually adjusted to that value while running.
    call dt_adjust_while_running() 

    ! 2. Calculate nonlinear terms in real space
    !$acc data present(w_r, flx_r)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          flx_r(k,j,i,iflx_uxx) = w_r(k,j,i,iux)*w_r(k,j,i,iux) - w_r(k,j,i,ibx)*w_r(k,j,i,ibx)
          flx_r(k,j,i,iflx_uxy) = w_r(k,j,i,iux)*w_r(k,j,i,iuy) - w_r(k,j,i,ibx)*w_r(k,j,i,iby)
          flx_r(k,j,i,iflx_uxz) = w_r(k,j,i,iux)*w_r(k,j,i,iuz) - w_r(k,j,i,ibx)*w_r(k,j,i,ibz)
          flx_r(k,j,i,iflx_uyy) = w_r(k,j,i,iuy)*w_r(k,j,i,iuy) - w_r(k,j,i,iby)*w_r(k,j,i,iby)
          flx_r(k,j,i,iflx_uyz) = w_r(k,j,i,iuy)*w_r(k,j,i,iuz) - w_r(k,j,i,iby)*w_r(k,j,i,ibz)
          flx_r(k,j,i,iflx_uzz) = w_r(k,j,i,iuz)*w_r(k,j,i,iuz) - w_r(k,j,i,ibz)*w_r(k,j,i,ibz)

          flx_r(k,j,i,iflx_bx ) = w_r(k,j,i,iby)*w_r(k,j,i,iuz) - w_r(k,j,i,ibz)*w_r(k,j,i,iuy)
          flx_r(k,j,i,iflx_by ) = w_r(k,j,i,ibz)*w_r(k,j,i,iux) - w_r(k,j,i,ibx)*w_r(k,j,i,iuz)
          flx_r(k,j,i,iflx_bz ) = w_r(k,j,i,ibx)*w_r(k,j,i,iuy) - w_r(k,j,i,iby)*w_r(k,j,i,iux)
        enddo
      enddo
    enddo
    !$acc end data

    ! 3. Forward FFT
    do i = 1, nftran
      call ftran_r2c(flx_r(:,:,:,i), flx(:,:,:,i))
    enddo

    !$acc data present(flx)
    !$acc parallel loop collapse(4)
    do l = 1, nftran
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            flx(k,j,i,l) = flx(k,j,i,l)/ntot
          enddo
        enddo
      enddo
    enddo
    !$acc end data

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms_ux, exp_terms_uy, exp_terms_uz, &
                           exp_terms_bx, exp_terms_by, exp_terms_bz, &
                           ux, uy, uz, bx, by, bz, &
                           flx_uxx, flx_uxy, flx_uxz, &
                                    flx_uyy, flx_uyz, &
                                             flx_uzz, &
                           flx_bx , flx_by , flx_bz , &
                           fux, fuy, fuz, &
                           fbx, fby, fbz, &
                           kxt, ky, kz, k2t_inv)
    !$acc routine seq
    use params, only: zi, q
    use shearing_box, only: shear_flg
    implicit none
    complex(8), intent(out) :: exp_terms_ux, exp_terms_uy, exp_terms_uz, &
                               exp_terms_bx, exp_terms_by, exp_terms_bz
    complex(8), intent(in ) :: ux, uy, uz, bx, by, bz
    complex(8), intent(in ) :: flx_uxx, flx_uxy, flx_uxz, &
                                        flx_uyy, flx_uyz, &
                                                 flx_uzz, &
                               flx_bx , flx_by , flx_bz 
    complex(8), intent(in ) :: fux, fuy, fuz
    complex(8), intent(in ) :: fbx, fby, fbz
    real(8)   , intent(in)  :: kxt, ky, kz, k2t_inv
    complex(8)              :: nl(nfields)
    complex(8)              :: p

    ! div (uu - bb)
    nl(iux) = -zi*( kxt*flx_uxx + ky*flx_uxy + kz*flx_uxz )
    nl(iuy) = -zi*( kxt*flx_uxy + ky*flx_uyy + kz*flx_uyz )
    nl(iuz) = -zi*( kxt*flx_uxz + ky*flx_uyz + kz*flx_uzz )

    ! curl (b x u)
    nl(ibx) = -zi*( ky *flx_bz - kz *flx_by )
    nl(iby) = -zi*( kz *flx_bx - kxt*flx_bz )
    nl(ibz) = -zi*( kxt*flx_by - ky *flx_bx )

    ! get pressure
    ! The constraint in shearing coordinates is kxt*ux + ky*uy + kz*uz = 0 with
    ! kxt = kx + q*tsc*ky, so d(kxt)/dt = q*ky and differentiating it gives
    ! kxt*dux/dt + ky*duy/dt + kz*duz/dt = -q*ky*ux.  Solving that for p:
    p = -zi*( kxt*nl(iux) + ky*nl(iuy) + kz*nl(iuz) &
              + 2.d0*shear_flg*kxt*uy - (2.d0 - 2.d0*q)*shear_flg*ky*ux )*k2t_inv

    exp_terms_ux = nl(iux) + fux - zi*kxt*p + 2.d0*shear_flg*uy
    exp_terms_uy = nl(iuy) + fuy - zi*ky *p - (2.d0 - q)*shear_flg*ux
    exp_terms_uz = nl(iuz) + fuz - zi*kz *p
    
    exp_terms_bx = nl(ibx) + fbx
    exp_terms_by = nl(iby) + fby - q*shear_flg*bx
    exp_terms_bz = nl(ibz) + fbz

  end subroutine get_ext_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Time integral of hyperdissipation
!-----------------------------------------------!
  subroutine get_imp_terms_tintg(imp_terms_tintg, t, kx, ky, kz, coeff, coeff_h, nexp)
    !$acc routine seq
    use grid, only: k2_max
    use params, only: shear
    use shearing_box, only: get_imp_terms_tintg_with_shear
    implicit none
    real(8), intent(out) :: imp_terms_tintg
    real(8), intent(in)  :: t, kx, ky, kz, coeff, coeff_h
    integer, intent(in)  :: nexp

    if(shear) then
      call get_imp_terms_tintg_with_shear(imp_terms_tintg, t, kx, ky, kz, coeff, coeff_h, nexp )
    else
      imp_terms_tintg = -(coeff*((kx**2 + ky**2 + kz**2)/k2_max) + coeff_h*((kx**2 + ky**2 + kz**2)/k2_max)**nexp)*t
    endif
  end subroutine get_imp_terms_tintg


!-----------------------------------------------!
!> @author  YK
!! @date    15 Jul 2021
!! @brief   Output series modes
!-----------------------------------------------!
  subroutine output_series_modes
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use mp, only: proc0, sum_reduce
    use time, only: tt
    use diagnostics_common, only: n_series_modes, series_modes
    use diagnostics_common, only: series_modes_unit
    implicit none
    complex(8), dimension(n_series_modes) :: ux_modes, uy_modes, uz_modes
    complex(8), dimension(n_series_modes) :: bx_modes, by_modes, bz_modes
    integer :: n, i, j, k

    ux_modes(:) = 0.d0
    uy_modes(:) = 0.d0
    uz_modes(:) = 0.d0
    bx_modes(:) = 0.d0
    by_modes(:) = 0.d0
    bz_modes(:) = 0.d0
    !$acc enter data copyin(ux_modes, uy_modes, uz_modes, bx_modes, by_modes, bz_modes)

    !$acc data present(ux, uy, uz, bx, by, bz)
    !$acc parallel loop private(i, j, k)
    do n = 1, n_series_modes
      i = series_modes(n, 1)
      j = series_modes(n, 2)
      k = series_modes(n, 3)

      if(       (i >= 1 .and. i <= nkx      ) &
          .and. (j >= 1 .and. j <= nky_local) &
          .and. (k >= 1 .and. k <= nkx      ) &
        ) then

        ux_modes(n) = ux(k, j, i)
        uy_modes(n) = uy(k, j, i)
        uz_modes(n) = uz(k, j, i)

        bx_modes(n) = bx(k, j, i)
        by_modes(n) = by(k, j, i)
        bz_modes(n) = bz(k, j, i)

      endif
    enddo
    !$acc end data
      
    !$acc update host (ux_modes, uy_modes, uz_modes, bx_modes, by_modes, bz_modes)

    !$acc exit data delete(ux_modes, uy_modes, uz_modes, bx_modes, by_modes, bz_modes)

    call sum_reduce(ux_modes, 0)
    call sum_reduce(uy_modes, 0)
    call sum_reduce(uz_modes, 0)
    call sum_reduce(bx_modes, 0)
    call sum_reduce(by_modes, 0)
    call sum_reduce(bz_modes, 0)

    do n = 1, n_series_modes
      if(proc0) then
        i = series_modes(n, 1)
        j = series_modes(n, 2)
        k = series_modes(n, 3)
999 format(es30.21, A6, 5es30.21e3)
        write (unit=series_modes_unit, fmt=999) tt, 'ux', kx(i), ky(j), kz(k), ux_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'uy', kx(i), ky(j), kz(k), uy_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'uz', kx(i), ky(j), kz(k), uz_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'bx', kx(i), ky(j), kz(k), bx_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'by', kx(i), ky(j), kz(k), by_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'bz', kx(i), ky(j), kz(k), bz_modes(n)
        call flush(series_modes_unit) 

      endif
    enddo

  end subroutine output_series_modes


!-----------------------------------------------!
!> @author  YK
!! @date    3 Mar 2021
!! @brief   Remap
!-----------------------------------------------!
  subroutine remap
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use mp, only: proc0
    use time_stamp, only: put_time_stamp
    use shearing_box, only: tsc, nremap, kxt, k2t, k2t_inv
    use shearing_box, only: timer_remap, remap_fld
    use dealias, only: filter
    implicit none
    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_remap)

    if(proc0) then
      print *
      print *, 'remapping...'
      print *
    endif

    call remap_fld(ux)
    call remap_fld(uy)
    call remap_fld(uz)
    call remap_fld(bx)
    call remap_fld(by)
    call remap_fld(bz)

    tsc = 0.d0
    nremap = nremap + 1
    counter = 1
    !$acc update device (tsc, counter)

    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          ux(k,j,i) = ux(k,j,i)*filter(k,j,i)
          uy(k,j,i) = uy(k,j,i)*filter(k,j,i)
          uz(k,j,i) = uz(k,j,i)*filter(k,j,i)
          bx(k,j,i) = bx(k,j,i)*filter(k,j,i)
          by(k,j,i) = by(k,j,i)*filter(k,j,i)
          bz(k,j,i) = bz(k,j,i)*filter(k,j,i)
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(2)
    do i = 1, nkx
      do j = 1, nky_local
        kxt(i,j) = kx(i)
      enddo
    enddo

    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          k2t(k, j, i) = kx(i)**2 + ky(j)**2 + kz(k)**2

          if(k2t(k, j, i) == 0.d0) then
            k2t_inv(k, j, i) = 0.d0
          else
            k2t_inv(k, j, i) = 1.d0/k2t(k, j, i)
          endif
        enddo
      enddo
    enddo

    if (proc0) call put_time_stamp(timer_remap)
  end subroutine remap


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Enforce div free for (wx, wy, wz)
!-----------------------------------------------!
  subroutine div_free(wx, wy, wz)
    use grid, only: ky, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: zi
    use shearing_box, only: kxt, k2t_inv
    implicit none
    complex(8), dimension (:,:,:), intent(inout) :: wx, wy, wz
    complex(8), allocatable, dimension(:,:,:)   :: nbl2inv_div_w ! nabla^-2 (div w)
    integer :: i, j, k

    allocate(nbl2inv_div_w(nkz, nky_local, nkx), source=(0.d0,0.d0))
    !$acc enter data create(nbl2inv_div_w)

    !$acc data present(wx, wy, wz, nbl2inv_div_w)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          nbl2inv_div_w(k,j,i) = -zi*(   kxt(i,j)*wx(k,j,i) &
                                       + ky(j)   *wy(k,j,i) &
                                       + kz(k)   *wz(k,j,i) )*k2t_inv(k,j,i)

          wx(k,j,i) = wx(k,j,i) - zi*kxt(i,j)*nbl2inv_div_w(k,j,i)
          wy(k,j,i) = wy(k,j,i) - zi*ky(j)   *nbl2inv_div_w(k,j,i)
          wz(k,j,i) = wz(k,j,i) - zi*kz(k)   *nbl2inv_div_w(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

    !$acc exit data delete(nbl2inv_div_w)
    deallocate(nbl2inv_div_w)

  end subroutine div_free


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Enforce div free for (fwx, fwy, fwz)
!-----------------------------------------------!
  subroutine div_free_force(fwx, fwy, fwz)
    use grid, only: ky, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: zi
    use shearing_box, only: kxt
    use mp, only: max_allreduce
    implicit none
    complex(8), dimension (:,:,:), intent(inout) :: fwx, fwy, fwz
    complex(8), allocatable, dimension(:,:,:)   :: nbl2inv_div_fw ! nabla^-2 (div fw)
    integer :: i, j, k
    real(8) :: fwx_max, fwy_max, fwz_max
    real(8) :: eps, sx, sy, sz, k2t, k2t_inv

    allocate(nbl2inv_div_fw(nkz, nky_local, nkx), source=(0.d0,0.d0))
    !$acc enter data create(nbl2inv_div_fw)

    fwx_max = 0.d0; fwy_max = 0.d0; fwz_max = 0.d0
    !$acc data present(fwx, fwy, fwz, nbl2inv_div_fw)
    !$acc parallel loop collapse(3) reduction(max:fwx_max, fwy_max, fwz_max)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fwx_max = max(fwx_max, abs(fwx(k,j,i)))
          fwy_max = max(fwy_max, abs(fwy(k,j,i)))
          fwz_max = max(fwz_max, abs(fwz(k,j,i)))
        end do
      end do
    end do

    call max_allreduce(fwx_max)
    call max_allreduce(fwy_max)
    call max_allreduce(fwz_max)

    eps = 1d-10
    if(fwx_max < eps) then
      sx = 0.d0
    else
      sx = 1.d0
    endif

    if(fwy_max < eps) then
      sy = 0.d0
    else
      sy = 1.d0
    endif

    if(fwz_max < eps) then
      sz = 0.d0
    else
      sz = 1.d0
    endif

    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          k2t = sx*kxt(i, j)**2 + sy*ky(j)**2 + sz*kz(k)**2
          if(k2t == 0.d0) then
            k2t_inv = 0.d0
          else
            k2t_inv = 1.d0/k2t
          endif
          nbl2inv_div_fw(k,j,i) = -zi*(  kxt(i,j)*fwx(k,j,i) &
                                       + ky(j)   *fwy(k,j,i) &
                                       + kz(k)   *fwz(k,j,i) )*k2t_inv

          fwx(k,j,i) = fwx(k,j,i) - sx*zi*kxt(i,j)*nbl2inv_div_fw(k,j,i)
          fwy(k,j,i) = fwy(k,j,i) - sy*zi*ky(j)   *nbl2inv_div_fw(k,j,i)
          fwz(k,j,i) = fwz(k,j,i) - sz*zi*kz(k)   *nbl2inv_div_fw(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

    !$acc exit data delete(nbl2inv_div_fw)
    deallocate(nbl2inv_div_fw)

  end subroutine div_free_force


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Allocates arrays for this module
!-----------------------------------------------!
  subroutine allocate_advance
    use grid, only: nkx, nky_local, nkz
    use grid, only: nlx_local, nly, nlz_padded
    use params, only: time_step_scheme
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src

    if (is_allocated) then
      return
    else
      is_allocated = .true.
    endif

    allocate(src(nkz, nky_local, nkx), source=(0.d0,0.d0))

    allocate(uy_new, source=src)
    allocate(ux_new, source=src)
    allocate(uz_new, source=src)
    allocate(bx_new, source=src)
    allocate(by_new, source=src)
    allocate(bz_new, source=src)

    allocate(w        (nkz       , nky_local, nkx, nfields)      , source=(0.d0, 0.d0))
    allocate(flx      (nkz       , nky_local, nkx, nftran )      , source=(0.d0, 0.d0))
    allocate(w_r      (nlz_padded, nly      , nlx_local, nfields), source=0.d0)
    allocate(flx_r    (nlz_padded, nly      , nlx_local, nftran ), source=0.d0)
    allocate(exp_terms(nkz, nky_local, nkx, nfields), source=(0.d0, 0.d0))

    !$acc enter data copyin(ux_new)
    !$acc enter data copyin(uy_new)
    !$acc enter data copyin(uz_new)
    !$acc enter data copyin(bx_new)
    !$acc enter data copyin(by_new)
    !$acc enter data copyin(bz_new)
    !$acc enter data copyin(w        )
    !$acc enter data copyin(flx      )
    !$acc enter data copyin(w_r      )
    !$acc enter data copyin(flx_r    )
    !$acc enter data copyin(exp_terms)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      allocate(ux_tmp, source=src)
      allocate(uy_tmp, source=src)
      allocate(uz_tmp, source=src)
      allocate(bx_tmp, source=src)
      allocate(by_tmp, source=src)
      allocate(bz_tmp, source=src)

      allocate(exp_terms0(nkz, nky_local, nkx, nfields)); exp_terms0  = 0.d0

      !$acc enter data copyin(ux_tmp)
      !$acc enter data copyin(uy_tmp)
      !$acc enter data copyin(uz_tmp)
      !$acc enter data copyin(bx_tmp)
      !$acc enter data copyin(by_tmp)
      !$acc enter data copyin(bz_tmp)
      !$acc enter data copyin(exp_terms0)
    endif

    deallocate(src)

  end subroutine allocate_advance


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Deallocates arrays for this module
!-----------------------------------------------!
  subroutine deallocate_advance
    use params, only: time_step_scheme
    implicit none

    if (is_allocated) then
      is_allocated = .false.
    else
      return
    endif

    !$acc exit data delete(ux_new)
    !$acc exit data delete(uy_new)
    !$acc exit data delete(uz_new)
    !$acc exit data delete(bx_new)
    !$acc exit data delete(by_new)
    !$acc exit data delete(bz_new)
    !$acc exit data delete(w        )
    !$acc exit data delete(flx      )
    !$acc exit data delete(w_r      )
    !$acc exit data delete(flx_r    )
    !$acc exit data delete(exp_terms)

    deallocate(ux_new)
    deallocate(uy_new)
    deallocate(uz_new)
    deallocate(bx_new)
    deallocate(by_new)
    deallocate(bz_new)

    deallocate(w        )
    deallocate(flx      )
    deallocate(w_r      )
    deallocate(flx_r    )
    deallocate(exp_terms)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then

      !$acc exit data delete(ux_tmp)
      !$acc exit data delete(uy_tmp)
      !$acc exit data delete(uz_tmp)
      !$acc exit data delete(bx_tmp)
      !$acc exit data delete(by_tmp)
      !$acc exit data delete(bz_tmp)
      !$acc exit data delete(exp_terms0)
      deallocate(ux_tmp)
      deallocate(uy_tmp)
      deallocate(uz_tmp)
      deallocate(bx_tmp)
      deallocate(by_tmp)
      deallocate(bz_tmp)

      deallocate(exp_terms0)
    endif

  end subroutine deallocate_advance


end module advance
