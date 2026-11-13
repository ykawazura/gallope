!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Time stepping for HRMHD
!!          (RMHD + parallel flow upa and parallel
!!           magnetic field bpa, coupled by the Hall
!!           term rho and the slow-mode speed sgm)
!-----------------------------------------------!
module advance
  use fields, only: nfields
  use fields, only: iomg, ipsi, iupa, ibpa
  use model_specific, only: nbtran
  implicit none

  public solve, is_allocated, allocate_advance, deallocate_advance

  logical :: is_allocated = .false.

  integer :: counter = 0
  complex(8), allocatable, dimension(:,:,:)   :: phi_new
  complex(8), allocatable, dimension(:,:,:)   :: omg_new
  complex(8), allocatable, dimension(:,:,:)   :: psi_new
  complex(8), allocatable, dimension(:,:,:)   :: upa_new
  complex(8), allocatable, dimension(:,:,:)   :: bpa_new
  complex(8), allocatable, dimension(:,:,:,:) :: w
  real   (8), allocatable, dimension(:,:,:,:) :: w_r, nonlin_r
  complex(8), allocatable, dimension(:,:,:,:) :: nonlin, exp_terms

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(8), allocatable, dimension(:,:,:)   :: phi_tmp
  complex(8), allocatable, dimension(:,:,:)   :: omg_tmp
  complex(8), allocatable, dimension(:,:,:)   :: psi_tmp
  complex(8), allocatable, dimension(:,:,:)   :: upa_tmp
  complex(8), allocatable, dimension(:,:,:)   :: bpa_tmp
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  ! (the forcing history f*_old2 is never read, so it is not kept)
  complex(8), allocatable, dimension(:,:,:)   :: phi_old2
  complex(8), allocatable, dimension(:,:,:)   :: omg_old2
  complex(8), allocatable, dimension(:,:,:)   :: psi_old2
  complex(8), allocatable, dimension(:,:,:)   :: upa_old2
  complex(8), allocatable, dimension(:,:,:)   :: bpa_old2
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2

  integer :: cfl_unit, rms_unit

  ! Largest z phase speed of the linear waves: the KAW branch at kprp2_max.
  ! It does not depend on time, so it is computed once (init_work_fields).
  real(8) :: vmax_KAW = 0.d0

  ! The gradient count nbtran and the nonlinear terms are in model_specific
  ! (shared with the diagnostics).
  !$acc declare copyin(counter)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: phi, omg, psi, upa, bpa
    use fields, only: phi_old, omg_old, psi_old, upa_old, bpa_old
    use grid, only: kprp2, kprp2inv, kz2
    use grid, only: kz
    use grid, only: nkx, nky_local, nkz
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: nonlinear
    use force, only: fphi, fpsi, fupa, fbpa
    use force, only: fphi_old, fpsi_old, fupa_old, fbpa_old
    use force, only: driven, update_force
    use advance_common, only: gear1, gear2, gear3
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use params, only: time_step_scheme
    implicit none
    real(8) :: r_omg, r_psi, r_upa, r_bpa
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
      ! Calculate force terms at n
      if (driven) call set_force

      !---------------  RK 1st step  ---------------
      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, upa, bpa, .true.)

      !$acc data present(exp_terms, exp_terms0, nonlin, kz, kprp2, kz2, kprp2inv, &
      !$acc              phi    , psi    , omg    , upa    , bpa    , &
      !$acc              phi_tmp, psi_tmp, omg_tmp, upa_tmp, bpa_tmp, &
      !$acc              fphi   , fpsi   , fupa   , fbpa)
      !$acc parallel loop collapse(3) private(r_omg, r_psi, r_upa, r_bpa)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iomg), exp_terms(k,j,i,ipsi), &
                               exp_terms(k,j,i,iupa), exp_terms(k,j,i,ibpa), &
                               phi(k,j,i), psi(k,j,i), upa(k,j,i), bpa(k,j,i), &
                               nonlin(k,j,i,iomg), nonlin(k,j,i,ipsi), &
                               nonlin(k,j,i,iupa), nonlin(k,j,i,ibpa), &
                               fphi(k,j,i), fpsi(k,j,i), fupa(k,j,i), fbpa(k,j,i), &
                               kz(k), kprp2(k,j,i))

            ! Dissipation rates; their time integral from 0 to t is -rate*t
            call get_imp_rates(r_omg, r_psi, r_upa, r_bpa, kprp2(k,j,i), kz2(k))

            call eSSPIFRK1(omg_tmp(k,j,i), omg(k,j,i), exp_terms(k,j,i,iomg), -r_omg*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK1(psi_tmp(k,j,i), psi(k,j,i), exp_terms(k,j,i,ipsi), -r_psi*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK1(upa_tmp(k,j,i), upa(k,j,i), exp_terms(k,j,i,iupa), -r_upa*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK1(bpa_tmp(k,j,i), bpa(k,j,i), exp_terms(k,j,i,ibpa), -r_bpa*(2.d0/3.d0*dt), 0.d0)

            phi_tmp(k,j,i) = omg_tmp(k,j,i)*(-kprp2inv(k,j,i))

            ! save explicit terms at the previous step
            do l = 1, nfields
              exp_terms0(k,j,i,l) = exp_terms(k,j,i,l)
            enddo
          enddo
        enddo
      enddo
      !$acc end data

      call dealias_fields(phi_tmp, omg_tmp, psi_tmp, upa_tmp, bpa_tmp)
      call zero_unforced (phi_tmp, omg_tmp, psi_tmp, upa_tmp, bpa_tmp)

      !---------------  RK 2nd step  ---------------
      ! Calculate force terms
      if (driven) then
        ! at n + 2/3
        call update_force(2.d0/3.d0*dt)
        call set_force

        ! go to n + 1
        call update_force(1.d0/3.d0*dt)
      endif

      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, upa_tmp, bpa_tmp, .false.)

      !$acc data present(exp_terms, nonlin, kz, kprp2, kz2, kprp2inv, &
      !$acc              phi    , psi    , omg    , upa    , bpa    , &
      !$acc              phi_tmp, psi_tmp, omg_tmp, upa_tmp, bpa_tmp, &
      !$acc              fphi   , fpsi   , fupa   , fbpa)
      !$acc parallel loop collapse(3) private(r_omg, r_psi, r_upa, r_bpa)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iomg), exp_terms(k,j,i,ipsi), &
                               exp_terms(k,j,i,iupa), exp_terms(k,j,i,ibpa), &
                               phi_tmp(k,j,i), psi_tmp(k,j,i), upa_tmp(k,j,i), bpa_tmp(k,j,i), &
                               nonlin(k,j,i,iomg), nonlin(k,j,i,ipsi), &
                               nonlin(k,j,i,iupa), nonlin(k,j,i,ibpa), &
                               fphi(k,j,i), fpsi(k,j,i), fupa(k,j,i), fbpa(k,j,i), &
                               kz(k), kprp2(k,j,i))

            call get_imp_rates(r_omg, r_psi, r_upa, r_bpa, kprp2(k,j,i), kz2(k))

            call eSSPIFRK2(omg_tmp(k,j,i), omg_tmp(k,j,i), omg(k,j,i), exp_terms(k,j,i,iomg), &
                           -r_omg*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK2(psi_tmp(k,j,i), psi_tmp(k,j,i), psi(k,j,i), exp_terms(k,j,i,ipsi), &
                           -r_psi*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK2(upa_tmp(k,j,i), upa_tmp(k,j,i), upa(k,j,i), exp_terms(k,j,i,iupa), &
                           -r_upa*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK2(bpa_tmp(k,j,i), bpa_tmp(k,j,i), bpa(k,j,i), exp_terms(k,j,i,ibpa), &
                           -r_bpa*(2.d0/3.d0*dt), 0.d0)

            phi_tmp(k,j,i) = omg_tmp(k,j,i)*(-kprp2inv(k,j,i))
          enddo
        enddo
      enddo
      !$acc end data

      call dealias_fields(phi_tmp, omg_tmp, psi_tmp, upa_tmp, bpa_tmp)
      call zero_unforced (phi_tmp, omg_tmp, psi_tmp, upa_tmp, bpa_tmp)

      !---------------  RK 3rd step  ---------------
      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, upa_tmp, bpa_tmp, .false.)

      !$acc data present(exp_terms, exp_terms0, nonlin, kz, kprp2, kz2, kprp2inv, &
      !$acc              phi    , psi    , omg    , upa    , bpa    , &
      !$acc              phi_tmp, psi_tmp, omg_tmp, upa_tmp, bpa_tmp, &
      !$acc              phi_new, psi_new, omg_new, upa_new, bpa_new, &
      !$acc              fphi   , fpsi   , fupa   , fbpa)
      !$acc parallel loop collapse(3) private(r_omg, r_psi, r_upa, r_bpa)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iomg), exp_terms(k,j,i,ipsi), &
                               exp_terms(k,j,i,iupa), exp_terms(k,j,i,ibpa), &
                               phi_tmp(k,j,i), psi_tmp(k,j,i), upa_tmp(k,j,i), bpa_tmp(k,j,i), &
                               nonlin(k,j,i,iomg), nonlin(k,j,i,ipsi), &
                               nonlin(k,j,i,iupa), nonlin(k,j,i,ibpa), &
                               fphi(k,j,i), fpsi(k,j,i), fupa(k,j,i), fbpa(k,j,i), &
                               kz(k), kprp2(k,j,i))

            call get_imp_rates(r_omg, r_psi, r_upa, r_bpa, kprp2(k,j,i), kz2(k))

            call eSSPIFRK3(omg_new(k,j,i), omg_tmp(k,j,i), omg(k,j,i), &
                           exp_terms(k,j,i,iomg), exp_terms0(k,j,i,iomg), &
                           -r_omg*dt, -r_omg*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK3(psi_new(k,j,i), psi_tmp(k,j,i), psi(k,j,i), &
                           exp_terms(k,j,i,ipsi), exp_terms0(k,j,i,ipsi), &
                           -r_psi*dt, -r_psi*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK3(upa_new(k,j,i), upa_tmp(k,j,i), upa(k,j,i), &
                           exp_terms(k,j,i,iupa), exp_terms0(k,j,i,iupa), &
                           -r_upa*dt, -r_upa*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK3(bpa_new(k,j,i), bpa_tmp(k,j,i), bpa(k,j,i), &
                           exp_terms(k,j,i,ibpa), exp_terms0(k,j,i,ibpa), &
                           -r_bpa*dt, -r_bpa*(2.d0/3.d0*dt), 0.d0)

            phi_new(k,j,i) = omg_new(k,j,i)*(-kprp2inv(k,j,i))
          enddo
        enddo
      enddo
      !$acc end data

      ! save fields at the previous step
      !$acc data present(phi, omg, psi, upa, bpa, phi_old, omg_old, psi_old, upa_old, bpa_old)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            phi_old(k,j,i) = phi(k,j,i)
            omg_old(k,j,i) = omg(k,j,i)
            psi_old(k,j,i) = psi(k,j,i)
            upa_old(k,j,i) = upa(k,j,i)
            bpa_old(k,j,i) = bpa(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data

      if (driven) call save_force_old

      call dealias_fields(phi_new, omg_new, psi_new, upa_new, bpa_new)
      call zero_unforced (phi_new, omg_new, psi_new, upa_new, bpa_new)
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      ! Calculate force terms
      if (driven) then
        call update_force(dt)
        call set_force
      endif

      ! Calculate nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, upa, bpa, .true.)

      !$acc data present(exp_terms, exp_terms_old, exp_terms_old2, nonlin, kz, kprp2, kz2, kprp2inv, &
      !$acc              phi     , psi     , omg     , upa     , bpa     , &
      !$acc              phi_old , psi_old , omg_old , upa_old , bpa_old , &
      !$acc              phi_old2, psi_old2, omg_old2, upa_old2, bpa_old2, &
      !$acc              phi_new , psi_new , omg_new , upa_new , bpa_new , &
      !$acc              fphi    , fpsi    , fupa    , fbpa)
      !$acc parallel loop collapse(3) private(r_omg, r_psi, r_upa, r_bpa)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(k,j,i,iomg), exp_terms(k,j,i,ipsi), &
                               exp_terms(k,j,i,iupa), exp_terms(k,j,i,ibpa), &
                               phi(k,j,i), psi(k,j,i), upa(k,j,i), bpa(k,j,i), &
                               nonlin(k,j,i,iomg), nonlin(k,j,i,ipsi), &
                               nonlin(k,j,i,iupa), nonlin(k,j,i,ibpa), &
                               fphi(k,j,i), fpsi(k,j,i), fupa(k,j,i), fbpa(k,j,i), &
                               kz(k), kprp2(k,j,i))

            call get_imp_rates(r_omg, r_psi, r_upa, r_bpa, kprp2(k,j,i), kz2(k))

            ! 1st order
            if(counter == 1) then
              call gear1(omg_new(k,j,i), omg(k,j,i), exp_terms(k,j,i,iomg), r_omg)
              call gear1(psi_new(k,j,i), psi(k,j,i), exp_terms(k,j,i,ipsi), r_psi)
              call gear1(upa_new(k,j,i), upa(k,j,i), exp_terms(k,j,i,iupa), r_upa)
              call gear1(bpa_new(k,j,i), bpa(k,j,i), exp_terms(k,j,i,ibpa), r_bpa)

            ! 2nd order
            elseif(counter == 2) then
              call gear2(omg_new(k,j,i), omg(k,j,i), omg_old(k,j,i), &
                         exp_terms(k,j,i,iomg), exp_terms_old(k,j,i,iomg), r_omg)
              call gear2(psi_new(k,j,i), psi(k,j,i), psi_old(k,j,i), &
                         exp_terms(k,j,i,ipsi), exp_terms_old(k,j,i,ipsi), r_psi)
              call gear2(upa_new(k,j,i), upa(k,j,i), upa_old(k,j,i), &
                         exp_terms(k,j,i,iupa), exp_terms_old(k,j,i,iupa), r_upa)
              call gear2(bpa_new(k,j,i), bpa(k,j,i), bpa_old(k,j,i), &
                         exp_terms(k,j,i,ibpa), exp_terms_old(k,j,i,ibpa), r_bpa)

            ! 3rd order
            else
              call gear3(omg_new(k,j,i), omg(k,j,i), omg_old(k,j,i), omg_old2(k,j,i), &
                         exp_terms(k,j,i,iomg), exp_terms_old(k,j,i,iomg), exp_terms_old2(k,j,i,iomg), r_omg)
              call gear3(psi_new(k,j,i), psi(k,j,i), psi_old(k,j,i), psi_old2(k,j,i), &
                         exp_terms(k,j,i,ipsi), exp_terms_old(k,j,i,ipsi), exp_terms_old2(k,j,i,ipsi), r_psi)
              call gear3(upa_new(k,j,i), upa(k,j,i), upa_old(k,j,i), upa_old2(k,j,i), &
                         exp_terms(k,j,i,iupa), exp_terms_old(k,j,i,iupa), exp_terms_old2(k,j,i,iupa), r_upa)
              call gear3(bpa_new(k,j,i), bpa(k,j,i), bpa_old(k,j,i), bpa_old2(k,j,i), &
                         exp_terms(k,j,i,ibpa), exp_terms_old(k,j,i,ibpa), exp_terms_old2(k,j,i,ibpa), r_bpa)

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
            phi_old (k,j,i) = phi(k,j,i)

            omg_old2(k,j,i) = omg_old(k,j,i)
            omg_old (k,j,i) = omg(k,j,i)

            psi_old2(k,j,i) = psi_old(k,j,i)
            psi_old (k,j,i) = psi(k,j,i)

            upa_old2(k,j,i) = upa_old(k,j,i)
            upa_old (k,j,i) = upa(k,j,i)

            bpa_old2(k,j,i) = bpa_old(k,j,i)
            bpa_old (k,j,i) = bpa(k,j,i)

            do l = 1, nfields
              exp_terms_old2(k,j,i,l) = exp_terms_old(k,j,i,l)
              exp_terms_old (k,j,i,l) = exp_terms(k,j,i,l)
            enddo
          enddo
        enddo
      enddo
      !$acc end data

      if (driven) call save_force_old

      ! Calliope zeroes the time-n arrays here, which the copy-back below
      ! overwrites; the new fields are the ones to zero.
      call dealias_fields(phi_new, omg_new, psi_new, upa_new, bpa_new)
      call zero_unforced (phi_new, omg_new, psi_new, upa_new, bpa_new)
    endif

    !$acc data present(phi, psi, omg, upa, bpa, phi_new, psi_new, omg_new, upa_new, bpa_new)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi_new(k,j,i)
          psi(k,j,i) = psi_new(k,j,i)
          omg(k,j,i) = omg_new(k,j,i)
          upa(k,j,i) = upa_new(k,j,i)
          bpa(k,j,i) = bpa_new(k,j,i)
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
!! @brief   Set the forcing fphi, fpsi, fupa, fbpa
!!          at the current forcing time
!-----------------------------------------------!
  subroutine set_force
    use grid, only: nkx, nky_local, nkz
    use force, only: fphi, fpsi, fupa, fbpa, fzppe, fzmpe, elsasser
    use force, only: get_force, normalize_force, normalize_force_els
    implicit none
    integer :: i, j, k

    if(elsasser) then
      ! only the Alfvenic Elsasser fields are forced (fupa = fbpa = 0)
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
      call get_force('upa', fupa)
      call get_force('bpa', fbpa)
      call normalize_force(fphi, fpsi, fupa, fbpa)
    endif

  end subroutine set_force


!-----------------------------------------------!
!> @author  YK
!! @brief   f*_old = f* (read by the diagnostics)
!-----------------------------------------------!
  subroutine save_force_old
    use grid, only: nkx, nky_local, nkz
    use force, only: fphi, fpsi, fupa, fbpa
    use force, only: fphi_old, fpsi_old, fupa_old, fbpa_old
    implicit none
    integer :: i, j, k

    !$acc data present(fphi, fpsi, fupa, fbpa, fphi_old, fpsi_old, fupa_old, fbpa_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fphi_old(k,j,i) = fphi(k,j,i)
          fpsi_old(k,j,i) = fpsi(k,j,i)
          fupa_old(k,j,i) = fupa(k,j,i)
          fbpa_old(k,j,i) = fbpa(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine save_force_old


!-----------------------------------------------!
!> @author  YK
!! @brief   Apply the dealiasing filter
!-----------------------------------------------!
  subroutine dealias_fields(phi, omg, psi, upa, bpa)
    use grid, only: nkx, nky_local, nkz
    use dealias, only: filter
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: phi, omg, psi, upa, bpa
    integer :: i, j, k

    !$acc data present(phi, omg, psi, upa, bpa, filter)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          phi(k,j,i) = phi(k,j,i)*filter(k,j,i)
          omg(k,j,i) = omg(k,j,i)*filter(k,j,i)
          psi(k,j,i) = psi(k,j,i)*filter(k,j,i)
          upa(k,j,i) = upa(k,j,i)*filter(k,j,i)
          bpa(k,j,i) = bpa(k,j,i)*filter(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine dealias_fields


!-----------------------------------------------!
!> @author  YK
!! @brief   Keep the unforced fields zero at the
!!          forcing wavenumbers
!!
!!          omg is the evolved variable and phi is
!!          derived from it, so both are zeroed
!!          (Calliope zeroes phi only, which the
!!          next phi = -omg/kprp2 undoes).
!-----------------------------------------------!
  subroutine zero_unforced(phi, omg, psi, upa, bpa)
    use force, only: driven, elsasser, is_forced, zeroing_at_forcing_wavenumber
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: phi, omg, psi, upa, bpa

    if (.not. driven) return

    if(elsasser) then
      if ((.not. is_forced('zppa')) .and. (.not. is_forced('zmpa'))) then
        call zeroing_at_forcing_wavenumber(upa)
        call zeroing_at_forcing_wavenumber(bpa)
      endif
    else
      if (.not. is_forced('phi')) then
        call zeroing_at_forcing_wavenumber(phi)
        call zeroing_at_forcing_wavenumber(omg)
      endif
      if (.not. is_forced('psi')) call zeroing_at_forcing_wavenumber(psi)
      if (.not. is_forced('upa')) call zeroing_at_forcing_wavenumber(upa)
      if (.not. is_forced('bpa')) call zeroing_at_forcing_wavenumber(bpa)
    endif

  end subroutine zero_unforced


!-----------------------------------------------!
!> @author  YK
!! @brief   Allocate fields used only here
!-----------------------------------------------!
  subroutine init_work_fields
    use grid, only: nkx, nky_local, nkz
    use grid, only: kprp2_max
    use params, only: time_step_scheme, rho, sgm
    use file, only: open_output_file
    use mp, only: proc0
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src
    real(8) :: x


    call allocate_advance

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      allocate(src(nkz, nky_local, nkx), source=(0.d0,0.d0))

      allocate(phi_old2 , source=src)
      allocate(omg_old2 , source=src)
      allocate(psi_old2 , source=src)
      allocate(upa_old2 , source=src)
      allocate(bpa_old2 , source=src)

      allocate(exp_terms_old (nkz, nky_local, nkx, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(nkz, nky_local, nkx, nfields)); exp_terms_old2 = 0.d0

      deallocate(src)

      !$acc enter data copyin( phi_old2)
      !$acc enter data copyin( omg_old2)
      !$acc enter data copyin( psi_old2)
      !$acc enter data copyin( upa_old2)
      !$acc enter data copyin( bpa_old2)

      !$acc enter data copyin(exp_terms_old )
      !$acc enter data copyin(exp_terms_old2)
    endif

    ! Linear waves: lambda^4 - (1 + sgm^2 + rho^2 kprp2) lambda^2 + sgm^2 = 0,
    ! lambda = omega/kz. The KAW root grows with kprp2, so its maximum over the
    ! grid is at kprp2_max (max(1, sgm) at rho = 0).
    x = rho**2*kprp2_max
    vmax_KAW = sqrt(0.5d0*( 1.d0 + sgm**2 + x &
                           + sqrt( (x + (1.d0 + sgm)**2)*(x + (1.d0 - sgm)**2) ) ))
    if(proc0) print '(" Max. linear phase speed along z (KAW) = ", es12.4e3)', vmax_KAW

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
  subroutine get_nonlinear_terms(phi, psi, upa, bpa, dt_reset)
    use grid, only: dlx, dly, dlz
    use grid, only: nlx_local, nly, nlz_padded
    use grid, only: ntot
    use mp, only: proc0, max_allreduce, sum_allreduce, comm_fft
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms
    use advance_common, only: dt_adjust_while_running
    use model_specific, only: idphi_dx, idphi_dy, idpsi_dx, idpsi_dy
    use model_specific, only: gradients_to_real, nonlinear_to_spectral
    implicit none
    complex(8), dimension (:,:,:), intent(in) :: phi, psi, upa, bpa

    logical, intent(in) :: dt_reset

    integer :: i, j, k
    real   (8) :: ux_rms, uy_rms, bx_rms, by_rms
    real   (8) :: max_vel_x, max_vel_y , dt_cfl, dt_digit, dt_old
    !$acc declare create(ux_rms, uy_rms, bx_rms, by_rms)

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    ! 1. Gradients in real space
    call gradients_to_real(phi, psi, upa, bpa, w, w_r)


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
      ! The linear waves along z are explicit too: add their CFL limit.
      dt_cfl = cfl*min(dlx/max_vel_x, dly/max_vel_y, dlz/vmax_KAW)

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

      ! When the file 'dt_adjust' including a float number is created,
      ! dt will be manually adjusted to that value while running.
      ! This is done only at the first stage of a step, like the CFL reset
      ! above, so a step never mixes two dt values, and a change restarts
      ! the gear sequence.
      dt_old = dt
      call dt_adjust_while_running()
      if (dt /= dt_old) counter = 1

      ! dt and counter may have changed on the host
      !$acc update device (dt, counter)
    endif

    ! 2. Nonlinear terms in real space, back to Fourier space
    call nonlinear_to_spectral(w_r, nonlin_r, nonlin)

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @brief   Calculate explicit terms
!!
!!   d_t omg = -{phi,omg} + {psi,j} + d_z j - kprp2*fphi
!!   d_t psi =  {psi,phi + rho*b} + d_z(phi + rho*b) + fpsi
!!   d_t u   = -{phi,u} + sgm*({psi,b} + d_z b) + fupa
!!   d_t b   = -{phi,b} + {psi,sgm*u - rho*j} + d_z(sgm*u - rho*j) + fbpa
!!
!!   with {f,g} = f_x g_y - f_y g_x and j = -kprp2*psi
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms_omg, exp_terms_psi, exp_terms_upa, exp_terms_bpa, &
                           phi, psi, upa, bpa, &
                           nonlin_omg, nonlin_psi, nonlin_upa, nonlin_bpa, &
                           fphi, fpsi, fupa, fbpa, &
                           kz, kprp2)
    !$acc routine seq
    use params, only: zi, rho, sgm
    implicit none
    complex(8), intent(out) :: exp_terms_omg, exp_terms_psi, exp_terms_upa, exp_terms_bpa
    complex(8), intent(in ) :: phi, psi, upa, bpa
    complex(8), intent(in ) :: nonlin_omg, nonlin_psi, nonlin_upa, nonlin_bpa
    complex(8), intent(in ) :: fphi, fpsi, fupa, fbpa
    real(8)   , intent(in)  :: kz, kprp2

    exp_terms_omg = nonlin_omg - kprp2*fphi - zi*kz*kprp2*psi
    exp_terms_psi = nonlin_psi + fpsi + zi*kz*(phi + rho*bpa)
    exp_terms_upa = nonlin_upa + fupa + zi*kz*sgm*bpa
    exp_terms_bpa = nonlin_bpa + fbpa + zi*kz*(sgm*upa + rho*kprp2*psi)

  end subroutine get_ext_terms


!-----------------------------------------------!
!> @author  YK
!! @brief   Hyperdissipation rates
!!          (omg: nupe, psi: etape, upa: nupa,
!!           bpa: etapa)
!-----------------------------------------------!
  subroutine get_imp_rates(r_omg, r_psi, r_upa, r_bpa, kprp2, kz2)
    !$acc routine seq
    use grid, only: kprp2_max, kz2_max
    use params, only: nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp, &
                      nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                      etapa_x, etapa_x_exp, etapa_z, etapa_z_exp
    implicit none
    real(8), intent(out) :: r_omg, r_psi, r_upa, r_bpa
    real(8), intent(in)  :: kprp2, kz2

    r_omg = nupe_x *(kprp2/kprp2_max)**nupe_x_exp  + nupe_z *(kz2/kz2_max)**nupe_z_exp
    r_psi = etape_x*(kprp2/kprp2_max)**etape_x_exp + etape_z*(kz2/kz2_max)**etape_z_exp
    r_upa = nupa_x *(kprp2/kprp2_max)**nupa_x_exp  + nupa_z *(kz2/kz2_max)**nupa_z_exp
    r_bpa = etapa_x*(kprp2/kprp2_max)**etapa_x_exp + etapa_z*(kz2/kz2_max)**etapa_z_exp

  end subroutine get_imp_rates


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
    allocate(upa_new, source=src)
    allocate(bpa_new, source=src)

    allocate(w        (nkz       , nky_local, nkx      , nbtran ), source=(0.d0, 0.d0))
    allocate(nonlin   (nkz       , nky_local, nkx      , nfields), source=(0.d0, 0.d0))
    allocate(w_r      (nlz_padded, nly      , nlx_local, nbtran) , source=0.d0)
    allocate(nonlin_r (nlz_padded, nly      , nlx_local, nfields), source=0.d0)
    allocate(exp_terms(nkz       , nky_local, nkx      , nfields), source=(0.d0, 0.d0))

    !$acc enter data copyin(phi_new)
    !$acc enter data copyin(omg_new)
    !$acc enter data copyin(psi_new)
    !$acc enter data copyin(upa_new)
    !$acc enter data copyin(bpa_new)
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
      allocate(upa_tmp, source=src)
      allocate(bpa_tmp, source=src)

      allocate(exp_terms0(nkz, nky_local, nkx, nfields)); exp_terms0  = 0.d0

      !$acc enter data copyin(phi_tmp)
      !$acc enter data copyin(omg_tmp)
      !$acc enter data copyin(psi_tmp)
      !$acc enter data copyin(upa_tmp)
      !$acc enter data copyin(bpa_tmp)
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
    !$acc exit data delete(upa_new)
    !$acc exit data delete(bpa_new)
    !$acc exit data delete(w        )
    !$acc exit data delete(nonlin   )
    !$acc exit data delete(w_r      )
    !$acc exit data delete(nonlin_r )
    !$acc exit data delete(exp_terms)

    deallocate(phi_new)
    deallocate(omg_new)
    deallocate(psi_new)
    deallocate(upa_new)
    deallocate(bpa_new)

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
      !$acc exit data delete(upa_tmp)
      !$acc exit data delete(bpa_tmp)
      !$acc exit data delete(exp_terms0)
      deallocate(phi_tmp)
      deallocate(omg_tmp)
      deallocate(psi_tmp)
      deallocate(upa_tmp)
      deallocate(bpa_tmp)

      deallocate(exp_terms0)
    endif

  end subroutine deallocate_advance


end module advance
