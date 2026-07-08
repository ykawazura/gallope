!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Time stepping for RMHD
!-----------------------------------------------!
module advance
  use fields, only: nfields
  use fields, only: nfields
  use fields, only: iomg, ipsi
  implicit none

  public solve, is_allocated, allocate_advance, deallocate_advance

  logical :: is_allocated = .false.

  integer :: counter = 0
  complex(8), allocatable, dimension(:,:,:)   :: phi_new
  complex(8), allocatable, dimension(:,:,:)   :: omg_new
  complex(8), allocatable, dimension(:,:,:)   :: psi_new
  complex(8), allocatable, dimension(:,:,:,:) :: w
  real   (8), allocatable, dimension(:,:,:,:) :: w_r, nonlin_r 
  complex(8), allocatable, dimension(:,:,:,:) :: nonlin, exp_terms

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(8), allocatable, dimension(:,:,:)   :: phi_tmp
  complex(8), allocatable, dimension(:,:,:)   :: omg_tmp
  complex(8), allocatable, dimension(:,:,:)   :: psi_tmp
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(8), allocatable, dimension(:,:,:)   :: phi_old2
  complex(8), allocatable, dimension(:,:,:)   :: omg_old2
  complex(8), allocatable, dimension(:,:,:)   :: psi_old2
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2
  complex(8), allocatable, dimension(:,:,:)   :: fphi_old2
  complex(8), allocatable, dimension(:,:,:)   :: fpsi_old2

  integer :: cfl_unit, rms_unit

  ! Forward FFT variables
  integer, parameter :: nbtran = 8
  integer, parameter :: idphi_dx = 1, idphi_dy = 2
  integer, parameter :: idomg_dx = 3, idomg_dy = 4
  integer, parameter :: idpsi_dx = 5, idpsi_dy = 6
  integer, parameter :: idjpa_dx = 7, idjpa_dy = 8
  
  !$acc declare copyin(counter)
  !$acc declare copyin(nbtran)
  !$acc declare copyin(idphi_dx, idphi_dy)
  !$acc declare copyin(idomg_dx, idomg_dy)
  !$acc declare copyin(idpsi_dx, idpsi_dy)
  !$acc declare copyin(idjpa_dx, idjpa_dy)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: phi, omg, psi
    use fields, only: phi_old, omg_old, psi_old
    use grid, only: kprp2, kprp2inv, kz2, kprp2_max, kz2_max
    use grid, only: kz
    use grid, only: nkx, nky_local, nkz
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp, &
                      zi, nonlinear
    use dealias, only: filter
    use force, only: fphi, fpsi, fphi_old, fpsi_old, driven, elsasser, update_force, get_force, normalize_force, get_force, normalize_force_els
    use force, only: fzppe, fzmpe
    use advance_common, only: gear1, gear2, gear3
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
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
          call get_force('zppe', fzppe)
          call get_force('zmpe', fzmpe)
          call normalize_force_els(fzppe, fzmpe)
          !$acc data present(fphi, fpsi, fzppe, fzmpe)
          !$acc parallel loop collapse(3)
          do i = 1, nkx
            do j = 1, nky_local
              do k = 1, nkz
                fphi(k,j,i) = fzppe(k,j,i) + fzmpe(k,j,i)
                fpsi(k,j,i) = fzppe(k,j,i) - fzmpe(k,j,i)
              enddo
            enddo
          enddo
          !$acc end data
        else
          call get_force('phi', fphi)
          call get_force('psi', fpsi)
          call normalize_force(fphi, fpsi)
        endif
      endif

      !---------------  RK 1st step  ---------------
      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, .true.)

      !$acc data present(exp_terms, nonlin, kz, kprp2, kz2, kprp2inv, &
      !$acc              phi    , psi    , omg    , &
      !$acc              phi_tmp, psi_tmp, omg_tmp, &
      !$acc              fphi   , fpsi)
      !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg1)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iomg), exp_terms(k,j,i,ipsi), &
                               phi(k,j,i), psi(k,j,i), &
                               nonlin(k,j,i,iomg), nonlin(k,j,i,ipsi), &
                               fphi(k,j,i), fpsi(k,j,i), &
                               kz(k), kprp2(k,j,i))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(k,j,i), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iomg), 2.d0/3.d0*dt, kprp2(k,j,i), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(k,j,i), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ipsi), 2.d0/3.d0*dt, kprp2(k,j,i), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)

            call eSSPIFRK1(omg_tmp(k,j,i), omg(k,j,i), &
               exp_terms(k,j,i,iomg), &
               imp_terms_tintg1(iomg), imp_terms_tintg0(iomg) &
            )

            call eSSPIFRK1(psi_tmp(k,j,i), psi(k,j,i), &
               exp_terms(k,j,i,ipsi), &
               imp_terms_tintg1(ipsi), imp_terms_tintg0(ipsi) &
            )

            phi_tmp(k,j,i) = omg_tmp(k,j,i)*(-kprp2inv(k,j,i))
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
      !$acc data present(phi_tmp, psi_tmp, omg_tmp, filter)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            phi_tmp(k,j,i) = phi_tmp(k,j,i)*filter(k,j,i)
            psi_tmp(k,j,i) = psi_tmp(k,j,i)*filter(k,j,i)
            omg_tmp(k,j,i) = omg_tmp(k,j,i)*filter(k,j,i)
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
          call get_force('zppe', fzppe)
          call get_force('zmpe', fzmpe)
          call normalize_force_els(fzppe, fzmpe)
          !$acc data present(fphi, fpsi, fzppe, fzmpe)
          !$acc parallel loop collapse(3)
          do i = 1, nkx
            do j = 1, nky_local
              do k = 1, nkz
                fphi(k,j,i) = fzppe(k,j,i) + fzmpe(k,j,i)
                fpsi(k,j,i) = fzppe(k,j,i) - fzmpe(k,j,i)
              enddo
            enddo
          enddo
          !$acc end data
        else
          call get_force('phi', fphi)
          call get_force('psi', fpsi)
          call normalize_force(fphi, fpsi)
        endif

        ! go to n + 1
        call update_force(1.d0/3.d0*dt)
      endif

      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, .false.)

      !$acc data present(exp_terms, nonlin, kz, kprp2, kz2, kprp2inv, &
      !$acc              phi    , psi    , omg    , &
      !$acc              phi_tmp, psi_tmp, omg_tmp, &
      !$acc              fphi   , fpsi)
      !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg2)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iomg), exp_terms(k,j,i,ipsi), &
                               phi_tmp(k,j,i), psi_tmp(k,j,i), &
                               nonlin(k,j,i,iomg), nonlin(k,j,i,ipsi), &
                               fphi(k,j,i), fpsi(k,j,i), &
                               kz(k), kprp2(k,j,i))
            !
            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(k,j,i), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iomg), 2.d0/3.d0*dt, kprp2(k,j,i), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(k,j,i), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ipsi), 2.d0/3.d0*dt, kprp2(k,j,i), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)

            call eSSPIFRK2(omg_tmp(k,j,i), omg_tmp(k,j,i), omg(k,j,i), &
               exp_terms(k,j,i,iomg), &
               imp_terms_tintg2(iomg), imp_terms_tintg0(iomg) &
            )

            call eSSPIFRK2(psi_tmp(k,j,i), psi_tmp(k,j,i), psi(k,j,i), &
               exp_terms(k,j,i,ipsi), &
               imp_terms_tintg2(ipsi), imp_terms_tintg0(ipsi) &
            )

            phi_tmp(k,j,i) = omg_tmp(k,j,i)*(-kprp2inv(k,j,i))
          enddo
        enddo
      enddo
      !$acc end data

      ! Dealiasing
      !$acc data present(phi_tmp, psi_tmp, omg_tmp, filter)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            phi_tmp(k,j,i) = phi_tmp(k,j,i)*filter(k,j,i)
            psi_tmp(k,j,i) = psi_tmp(k,j,i)*filter(k,j,i)
            omg_tmp(k,j,i) = omg_tmp(k,j,i)*filter(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data

      !---------------  RK 3rd step  ---------------
      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, .false.)

      !$acc data present(exp_terms, nonlin, kz, kprp2, kz2, kprp2inv, &
      !$acc              phi    , psi    , omg    , &
      !$acc              phi_tmp, psi_tmp, omg_tmp, &
      !$acc              phi_new, psi_new, omg_new, &
      !$acc              fphi   , fpsi)
      !$acc parallel loop collapse(3) private(imp_terms_tintg0, imp_terms_tintg2, imp_terms_tintg3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iomg), exp_terms(k,j,i,ipsi), &
                               phi_tmp(k,j,i), psi_tmp(k,j,i), &
                               nonlin(k,j,i,iomg), nonlin(k,j,i,ipsi), &
                               fphi(k,j,i), fpsi(k,j,i), &
                               kz(k), kprp2(k,j,i))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(k,j,i), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iomg), 2.d0/3.d0*dt, kprp2(k,j,i), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iomg),           dt, kprp2(k,j,i), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(k,j,i), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ipsi), 2.d0/3.d0*dt, kprp2(k,j,i), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ipsi),           dt, kprp2(k,j,i), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)

            call eSSPIFRK3(omg_new(k,j,i), omg_tmp(k,j,i), omg(k,j,i), &
               exp_terms (k,j,i,iomg), exp_terms0(k,j,i,iomg), &
               imp_terms_tintg3(iomg), imp_terms_tintg2(iomg), imp_terms_tintg0(iomg) &
            )
            call eSSPIFRK3(psi_new(k,j,i), psi_tmp(k,j,i), psi(k,j,i), &
               exp_terms (k,j,i,ipsi), exp_terms0(k,j,i,ipsi), &
               imp_terms_tintg3(ipsi), imp_terms_tintg2(ipsi), imp_terms_tintg0(ipsi) &
            )

            phi_new(k,j,i) = omg_new(k,j,i)*(-kprp2inv(k,j,i))
          enddo
        enddo
      enddo
      !$acc end data

      ! save fields at the previous step
      !$acc data present(phi, omg, psi, phi_old, omg_old, psi_old)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            phi_old(k,j,i) = phi(k,j,i)
            omg_old(k,j,i) = omg(k,j,i)
            psi_old(k,j,i) = psi(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data

      if (driven) then
        !$acc data present(fphi, fpsi, fphi_old, fpsi_old)
        !$acc parallel loop collapse(3)
        do i = 1, nkx
          do j = 1, nky_local
            do k = 1, nkz
              fphi_old(k,j,i) = fphi(k,j,i)
              fpsi_old(k,j,i) = fpsi(k,j,i)
            enddo
          enddo
        enddo
        !$acc end data
      endif

      ! Dealiasing
      !$acc data present(phi_new, psi_new, omg_new, filter)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            phi_new(k,j,i) = phi_new(k,j,i)*filter(k,j,i)
            omg_new(k,j,i) = omg_new(k,j,i)*filter(k,j,i)
            psi_new(k,j,i) = psi_new(k,j,i)*filter(k,j,i)
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
          call get_force('zppe', fzppe)
          call get_force('zmpe', fzmpe)
          call normalize_force_els(fzppe, fzmpe)
          !$acc data present(fphi, fpsi, fzppe, fzmpe)
          !$acc parallel loop collapse(3)
          do i = 1, nkx
            do j = 1, nky_local
              do k = 1, nkz
                fphi(k,j,i) = fzppe(k,j,i) + fzmpe(k,j,i)
                fpsi(k,j,i) = fzppe(k,j,i) - fzmpe(k,j,i)
              enddo
            enddo
          enddo
          !$acc end data
        else
          call get_force('phi', fphi)
          call get_force('psi', fpsi)
          call normalize_force(fphi, fpsi)
        endif
      endif

      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, .true.)

      !$acc data present(exp_terms, exp_terms_old, exp_terms_old2, nonlin, kz, kprp2, kz2, kprp2inv, &
      !$acc              kz2_max, kprp2_max, &
      !$acc              phi      , psi      , omg     , &
      !$acc              phi_old  , psi_old  , omg_old , &
      !$acc              phi_old2 , psi_old2 , omg_old2, &
      !$acc              phi_new  , psi_new  , omg_new , &
      !$acc              fphi     , fpsi     , &
      !$acc              fphi_old , fpsi_old , &
      !$acc              fphi_old2, fpsi_old2)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iomg), exp_terms(k,j,i,ipsi), &
                               phi(k,j,i), psi(k,j,i), &
                               nonlin(k,j,i,iomg), nonlin(k,j,i,ipsi), &
                               fphi(k,j,i), fpsi(k,j,i), &
                               kz(k), kprp2(k,j,i))

            ! 1st order 
            if(counter == 1) then
              call gear1(omg_new(k,j,i), omg(k,j,i), &
                 exp_terms(k,j,i,iomg), &
                 nupe_x*(kprp2(k,j,i)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
              )
              call gear1(psi_new(k,j,i), psi(k,j,i), &
                 exp_terms(k,j,i,ipsi), &
                 etape_x*(kprp2(k,j,i)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
              )

            ! 2nd order 
            elseif(counter == 2) then
              call gear2(omg_new(k,j,i), omg(k,j,i), omg_old(k,j,i), &
                 exp_terms    (k,j,i,iomg), &
                 exp_terms_old(k,j,i,iomg), &
                 nupe_x*(kprp2(k,j,i)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
              )
              call gear2(psi_new(k,j,i), psi(k,j,i), psi_old(k,j,i), &
                 exp_terms    (k,j,i,ipsi), &
                 exp_terms_old(k,j,i,ipsi), &
                 etape_x*(kprp2(k,j,i)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
              )

            ! 3rd order 
            else
              call gear3(omg_new(k,j,i), omg(k,j,i), omg_old(k,j,i), omg_old2(k,j,i), &
                 exp_terms     (k,j,i,iomg), &
                 exp_terms_old (k,j,i,iomg), &
                 exp_terms_old2(k,j,i,iomg), &
                 nupe_x*(kprp2(k,j,i)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
              )
              call gear3(psi_new(k,j,i), psi(k,j,i), psi_old(k,j,i), psi_old2(k,j,i), &
                 exp_terms     (k,j,i,ipsi), &
                 exp_terms_old (k,j,i,ipsi), &
                 exp_terms_old2(k,j,i,ipsi), &
                 etape_x*(kprp2(k,j,i)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
              )

            endif
            phi_new(k,j,i) = omg_new(k,j,i)*(-kprp2inv(k,j,i))
          enddo
        enddo
      enddo

      if(counter <= 2) counter = counter + 1
      !$acc update device (counter)

      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            phi_old2(k,j,i) = phi_old(k,j,i) 
            phi_old(k,j,i)  = phi(k,j,i)

            omg_old2(k,j,i) = omg_old(k,j,i) 
            omg_old(k,j,i)  = omg(k,j,i)

            psi_old2(k,j,i) = psi_old(k,j,i) 
            psi_old(k,j,i)  = psi(k,j,i)

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
              fphi_old2(k,j,i) = fphi_old(k,j,i) 
              fphi_old(k,j,i)  = fphi(k,j,i)

              fpsi_old2(k,j,i) = fpsi_old(k,j,i) 
              fpsi_old(k,j,i)  = fpsi(k,j,i)
            enddo
          enddo
        enddo
      endif

      ! Dealiasing
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            phi_new(k,j,i) = phi_new(k,j,i)*filter(k,j,i)
            psi_new(k,j,i) = psi_new(k,j,i)*filter(k,j,i)
            omg_new(k,j,i) = omg_new(k,j,i)*filter(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data
    endif

    !$acc data present(phi, psi, omg, phi_new, psi_new, omg_new)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi_new(k,j,i)
          psi(k,j,i) = psi_new(k,j,i)
          omg(k,j,i) = omg_new(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

    tt  = tt  + dt
    !$acc update device (tt)

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @brief   Allocate fields used only here
!-----------------------------------------------!
  subroutine init_work_fields
    use grid, only: kx, ky, kz
    use grid, only: nkx, nky_local, nkz
    use params, only: time_step_scheme
    use file, only: open_output_file
    use mp, only: proc0
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src


    call allocate_advance

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      allocate(src(nkz, nky_local, nkx), source=(0.d0,0.d0))

      allocate(phi_old2 , source=src)
      allocate(omg_old2 , source=src)
      allocate(psi_old2 , source=src)

      allocate(fphi_old2, source=src)
      allocate(fpsi_old2, source=src)

      allocate(exp_terms_old (nkz, nky_local, nkx, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(nkz, nky_local, nkx, nfields)); exp_terms_old2 = 0.d0

      deallocate(src)

      !$acc enter data copyin( phi_old2)
      !$acc enter data copyin( omg_old2)
      !$acc enter data copyin( psi_old2)

      !$acc enter data copyin(fphi_old2)
      !$acc enter data copyin(fpsi_old2)
      
      !$acc enter data copyin(exp_terms_old )
      !$acc enter data copyin(exp_terms_old2)
    endif


    if(proc0) call open_output_file (cfl_unit, 'cfl.dat')
    if(proc0) call open_output_file (rms_unit, 'rms.dat')

  end subroutine init_work_fields


!-----------------------------------------------!
!> @author  YK
!! @brief   Calculate nonlinear terms via
!!          1. Calculate grad in Fourier space
!!          2. Inverse FFT
!!          3. Calculate nonlinear terms 
!!             in real space
!!          4. Forward FFT
!-----------------------------------------------!
  subroutine get_nonlinear_terms(phi, psi, dt_reset)
    use grid, only: dlx, dly, dlz
    use grid, only: kprp2, kx, ky
    use grid, only: nlx, nlx_local, nly, nlz_padded
    use grid, only: nkx, nky_local, nkz, kprp2
    use grid, only: ntot
    use params, only: zi
    use mp, only: proc0, max_allreduce, sum_allreduce, comm_fft
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms
    use advance_common, only: dt_adjust_while_running 
    use cuFFTmp, only: btran_c2r, ftran_r2c
    implicit none
    complex(8), dimension (:,:,:), intent(in) :: phi, psi

    logical, intent(in) :: dt_reset

    integer :: i, j, k, l
    real   (8) :: ux_rms, uy_rms, bx_rms, by_rms
    real   (8) :: max_vel_x, max_vel_y , dt_cfl, dt_digit
    !$acc declare create(ux_rms, uy_rms, bx_rms, by_rms)

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    !$acc data present(phi, psi, w)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          w(k,j,i,idphi_dx) = zi*kx(i)                *phi(k,j,i)
          w(k,j,i,idphi_dy) = zi*ky(j)                *phi(k,j,i)
          w(k,j,i,idomg_dx) = zi*kx(i)*(-kprp2(k,j,i))*phi(k,j,i)
          w(k,j,i,idomg_dy) = zi*ky(j)*(-kprp2(k,j,i))*phi(k,j,i)
                               
          w(k,j,i,idpsi_dx) = zi*kx(i)                *psi(k,j,i)
          w(k,j,i,idpsi_dy) = zi*ky(j)                *psi(k,j,i)
          w(k,j,i,idjpa_dx) = zi*kx(i)*(-kprp2(k,j,i))*psi(k,j,i)
          w(k,j,i,idjpa_dy) = zi*ky(j)*(-kprp2(k,j,i))*psi(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

    ! 1. Inverse FFT
    do i = 1, nbtran
      call btran_c2r(w(:,:,:,i), w_r(:,:,:,i))
    enddo


    ! write rms
    ux_rms = 0.d0; uy_rms = 0.d0
    bx_rms = 0.d0; by_rms = 0.d0
    !$acc data present(w_r)
    !$acc parallel loop collapse(3) reduction(+:ux_rms, uy_rms, bx_rms, by_rms)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          ux_rms = ux_rms + sqrt(w_r(k,j,i,idphi_dy)**2)
          uy_rms = uy_rms + sqrt(w_r(k,j,i,idphi_dx)**2)
          bx_rms = bx_rms + sqrt(w_r(k,j,i,idpsi_dy)**2)
          by_rms = by_rms + sqrt(w_r(k,j,i,idpsi_dx)**2)
        end do
      end do
    end do
    !$acc end data

    !$acc update host (ux_rms, uy_rms, bx_rms, by_rms)

    ! rms is a grid-space partial sum normalized by the global ntot: reduce on
    ! comm_fft so redundant comm_m/comm_s groups do not scale it by P_m*P_s.
    call sum_allreduce(ux_rms, comm=comm_fft); call sum_allreduce(uy_rms, comm=comm_fft)
    call sum_allreduce(bx_rms, comm=comm_fft); call sum_allreduce(by_rms, comm=comm_fft)

    ux_rms = sqrt(ux_rms/ntot); uy_rms = sqrt(uy_rms/ntot)
    bx_rms = sqrt(bx_rms/ntot); by_rms = sqrt(by_rms/ntot)

    if(proc0) then
      write (unit=rms_unit, fmt="(100es30.21)") tt, ux_rms, uy_rms, bx_rms, by_rms
      call flush(rms_unit) 
    endif


    ! (get max_vel for dt reset)
    if(dt_reset) then
      max_vel_x = 0.d0; max_vel_y = 0.d0
      !$acc data present(w_r)
      !$acc parallel loop collapse(3) reduction(max:max_vel_x, max_vel_y)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz_padded
            max_vel_x = max( &
                          max_vel_x, &
                          abs(w_r(k,j,i,idphi_dy) + w_r(k,j,i,idpsi_dy)), &
                          abs(w_r(k,j,i,idphi_dy) - w_r(k,j,i,idpsi_dy))  &
                        )
            max_vel_y = max( &
                          max_vel_y, &
                          abs(w_r(k,j,i,idphi_dx) + w_r(k,j,i,idpsi_dx)), &
                          abs(w_r(k,j,i,idphi_dx) - w_r(k,j,i,idpsi_dx))  &
                        )
          end do
        end do
      end do
      !$acc end data

      ! CFL velocity max must span all ranks (world). max is duplication-
      ! invariant, so keeping it on world is both correct and lockstep-safe:
      ! every group computes the same dt, avoiding step-count divergence.
      call max_allreduce(max_vel_x)
      call max_allreduce(max_vel_y)
      dt_cfl = cfl*min(dlx/max_vel_x, dly/max_vel_y)

      if(proc0) then
        write (unit=cfl_unit, fmt="(100es30.21)") tt, dt_cfl, max_vel_x, max_vel_y
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
    !$acc data present(w_r, nonlin_r)
    !$acc parallel loop collapse(3)
    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          nonlin_r(k,j,i,iomg) = - w_r(k,j,i,idphi_dx)*w_r(k,j,i,idomg_dy) &
                                 + w_r(k,j,i,idphi_dy)*w_r(k,j,i,idomg_dx) &
                                 + w_r(k,j,i,idpsi_dx)*w_r(k,j,i,idjpa_dy) &
                                 - w_r(k,j,i,idpsi_dy)*w_r(k,j,i,idjpa_dx)
          nonlin_r(k,j,i,ipsi) = - w_r(k,j,i,idphi_dx)*w_r(k,j,i,idpsi_dy) &
                                 + w_r(k,j,i,idphi_dy)*w_r(k,j,i,idpsi_dx)
        enddo
      enddo
    enddo
    !$acc end data

    ! 3. Forward FFT
    do i = 1, nfields
      call ftran_r2c(nonlin_r(:,:,:,i), nonlin(:,:,:,i))
    enddo

    !$acc data present(nonlin)
    !$acc parallel loop collapse(4)
    do l = 1, nfields
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            nonlin(k,j,i,l) = nonlin(k,j,i,l)/ntot
          enddo
        enddo
      enddo
    enddo
    !$acc end data

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms_omg, exp_terms_psi, &
                           phi, psi, &
                           nonlin_omg, nonlin_psi, &
                           fphi, fpsi, &
                           kz, kprp2)
    !$acc routine seq
    use params, only: zi
    implicit none
    complex(8), intent(out) :: exp_terms_omg, exp_terms_psi
    complex(8), intent(in ) :: phi, psi
    complex(8), intent(in ) :: nonlin_omg, nonlin_psi
    complex(8), intent(in ) :: fphi, fpsi
    real(8)   , intent(in)  :: kz, kprp2

    exp_terms_omg = nonlin_omg - kprp2*fphi - zi*kz*kprp2*psi
    exp_terms_psi = nonlin_psi + fpsi + zi*kz*phi

  end subroutine get_ext_terms


!-----------------------------------------------!
!> @author  YK
!! @brief   Time integral of hyperdissipation
!-----------------------------------------------!
  subroutine get_imp_terms_tintg(imp_terms_tintg, t, kprp2, kz2, coeff_x, coeff_z, nexp_x, nexp_z)
    !$acc routine seq
    use grid, only: kprp2_max, kz2_max
    implicit none
    real(8), intent(out) :: imp_terms_tintg
    real(8), intent(in) :: t, kprp2, kz2, coeff_x, coeff_z
    integer, intent(in) :: nexp_x, nexp_z

    imp_terms_tintg = -(coeff_x*(kprp2/kprp2_max)**nexp_x + coeff_z*(kz2/kz2_max)**nexp_z)*t
  end subroutine get_imp_terms_tintg


!-----------------------------------------------!
!> @author  YK
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

    allocate(phi_new, source=src)
    allocate(omg_new, source=src)
    allocate(psi_new, source=src)

    allocate(w        (nkz       , nky_local, nkx      , nbtran ), source=(0.d0, 0.d0))
    allocate(nonlin   (nkz       , nky_local, nkx      , nfields), source=(0.d0, 0.d0))
    allocate(w_r      (nlz_padded, nly      , nlx_local, nbtran) , source=0.d0)
    allocate(nonlin_r (nlz_padded, nly      , nlx_local, nfields), source=0.d0)
    allocate(exp_terms(nkz       , nky_local, nkx      , nfields), source=(0.d0, 0.d0))

    !$acc enter data copyin(phi_new)
    !$acc enter data copyin(omg_new)
    !$acc enter data copyin(psi_new)
    !$acc enter data copyin(w        )
    !$acc enter data copyin(nonlin   )
    !$acc enter data copyin(w_r      )
    !$acc enter data copyin(nonlin_r )
    !$acc enter data copyin(exp_terms)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      allocate(phi_tmp, source=src)
      allocate(omg_tmp, source=src)
      allocate(psi_tmp, source=src)

      allocate(exp_terms0(nkz, nky_local, nkx, nfields)); exp_terms0  = 0.d0

      !$acc enter data copyin(phi_tmp)
      !$acc enter data copyin(omg_tmp)
      !$acc enter data copyin(psi_tmp)
      !$acc enter data copyin(exp_terms0)
    endif

    deallocate(src)

  end subroutine allocate_advance


!-----------------------------------------------!
!> @author  YK
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

    !$acc exit data delete(phi_new)
    !$acc exit data delete(omg_new)
    !$acc exit data delete(psi_new)
    !$acc exit data delete(w        )
    !$acc exit data delete(nonlin   )
    !$acc exit data delete(w_r      )
    !$acc exit data delete(nonlin_r )
    !$acc exit data delete(exp_terms)

    deallocate(phi_new)
    deallocate(omg_new)
    deallocate(psi_new)

    deallocate(w        )
    deallocate(nonlin   )
    deallocate(w_r      )
    deallocate(nonlin_r )
    deallocate(exp_terms)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then

      !$acc exit data delete(phi_tmp)
      !$acc exit data delete(omg_tmp)
      !$acc exit data delete(psi_tmp)
      !$acc exit data delete(exp_terms0)
      deallocate(phi_tmp)
      deallocate(omg_tmp)
      deallocate(psi_tmp)

      deallocate(exp_terms0)
    endif

  end subroutine deallocate_advance


end module advance
