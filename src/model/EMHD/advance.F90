!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @brief   Time stepping for EMHD
!!
!!   (1 + de^2 k^2) d_t b = curl(u x b*) + fb
!!                          - eta (k^2/k2_max)^eta_exp b
!!
!!   with b* = (1 + de^2 k^2) b and u = -i k x b.
!!   The guide field is the k = 0 mode of b; it
!!   enters only through b*, so the linear whistler
!!   propagation is part of the "nonlinear" term.
!!   The resistive term is integrated exactly
!!   (integrating factor / implicit), the rest
!!   explicitly. b is projected onto div b = 0
!!   after every step.
!-----------------------------------------------!
module advance
  use fields, only: nfields
  use fields, only: ibx, iby, ibz
  use model_specific, only: nbtran, nftran
  implicit none

  public solve, is_allocated, allocate_advance, deallocate_advance

  logical :: is_allocated = .false.

  integer :: counter = 0
  complex(8), allocatable, dimension(:,:,:)   :: bx_new, by_new, bz_new
  complex(8), allocatable, dimension(:,:,:,:) :: w
  real   (8), allocatable, dimension(:,:,:,:) :: w_r, nonlin_r
  complex(8), allocatable, dimension(:,:,:,:) :: nonlin, exp_terms

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(8), allocatable, dimension(:,:,:)   :: bx_tmp, by_tmp, bz_tmp
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(8), allocatable, dimension(:,:,:)   :: bx_old2, by_old2, bz_old2
  complex(8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2

  integer :: cfl_unit

  ! Whistler CFL: omega = (k.B)|k|/(1 + de^2 k^2) <= |B| Kw with
  !   Kw = max over the dealiased modes of k^2/(1 + de^2 k^2).
  ! Time independent, so computed once (init_work_fields).
  real(8) :: Kw = 0.d0

  !$acc declare copyin(counter)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: bx, by, bz
    use fields, only: bx_old, by_old, bz_old
    use fields, only: div_free
    use grid, only: k2
    use grid, only: nkx, nky_local, nkz
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: nonlinear, de2
    use force, only: fbx, fby, fbz
    use force, only: driven, update_force
    use advance_common, only: gear1, gear2, gear3
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use params, only: time_step_scheme
    implicit none
    real(8) :: r, d
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
      if(nonlinear) call get_nonlinear_terms(bx, by, bz, .true.)

      !$acc data present(exp_terms, exp_terms0, nonlin, k2, &
      !$acc              bx, by, bz, bx_tmp, by_tmp, bz_tmp, fbx, fby, fbz)
      !$acc parallel loop collapse(3) private(r, d)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            d = 1.d0 + de2*k2(k,j,i)
            exp_terms(k,j,i,ibx) = (nonlin(k,j,i,ibx) + fbx(k,j,i))/d
            exp_terms(k,j,i,iby) = (nonlin(k,j,i,iby) + fby(k,j,i))/d
            exp_terms(k,j,i,ibz) = (nonlin(k,j,i,ibz) + fbz(k,j,i))/d

            ! Resistive rate; its time integral from 0 to t is -r*t
            r = imp_rate(k2(k,j,i))

            call eSSPIFRK1(bx_tmp(k,j,i), bx(k,j,i), exp_terms(k,j,i,ibx), -r*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK1(by_tmp(k,j,i), by(k,j,i), exp_terms(k,j,i,iby), -r*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK1(bz_tmp(k,j,i), bz(k,j,i), exp_terms(k,j,i,ibz), -r*(2.d0/3.d0*dt), 0.d0)

            ! save explicit terms at the previous step
            do l = 1, nfields
              exp_terms0(k,j,i,l) = exp_terms(k,j,i,l)
            enddo
          enddo
        enddo
      enddo
      !$acc end data

      call dealias_fields(bx_tmp, by_tmp, bz_tmp)

      !---------------  RK 2nd step  ---------------
      if (driven) then
        ! at n + 2/3
        call update_force(2.d0/3.d0*dt)
        call set_force

        ! go to n + 1
        call update_force(1.d0/3.d0*dt)
      endif

      if(nonlinear) call get_nonlinear_terms(bx_tmp, by_tmp, bz_tmp, .false.)

      !$acc data present(exp_terms, nonlin, k2, &
      !$acc              bx, by, bz, bx_tmp, by_tmp, bz_tmp, fbx, fby, fbz)
      !$acc parallel loop collapse(3) private(r, d)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            d = 1.d0 + de2*k2(k,j,i)
            exp_terms(k,j,i,ibx) = (nonlin(k,j,i,ibx) + fbx(k,j,i))/d
            exp_terms(k,j,i,iby) = (nonlin(k,j,i,iby) + fby(k,j,i))/d
            exp_terms(k,j,i,ibz) = (nonlin(k,j,i,ibz) + fbz(k,j,i))/d

            r = imp_rate(k2(k,j,i))

            call eSSPIFRK2(bx_tmp(k,j,i), bx_tmp(k,j,i), bx(k,j,i), exp_terms(k,j,i,ibx), &
                           -r*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK2(by_tmp(k,j,i), by_tmp(k,j,i), by(k,j,i), exp_terms(k,j,i,iby), &
                           -r*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK2(bz_tmp(k,j,i), bz_tmp(k,j,i), bz(k,j,i), exp_terms(k,j,i,ibz), &
                           -r*(2.d0/3.d0*dt), 0.d0)
          enddo
        enddo
      enddo
      !$acc end data

      call dealias_fields(bx_tmp, by_tmp, bz_tmp)

      !---------------  RK 3rd step  ---------------
      if(nonlinear) call get_nonlinear_terms(bx_tmp, by_tmp, bz_tmp, .false.)

      !$acc data present(exp_terms, exp_terms0, nonlin, k2, &
      !$acc              bx, by, bz, bx_tmp, by_tmp, bz_tmp, bx_new, by_new, bz_new, &
      !$acc              fbx, fby, fbz)
      !$acc parallel loop collapse(3) private(r, d)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            d = 1.d0 + de2*k2(k,j,i)
            exp_terms(k,j,i,ibx) = (nonlin(k,j,i,ibx) + fbx(k,j,i))/d
            exp_terms(k,j,i,iby) = (nonlin(k,j,i,iby) + fby(k,j,i))/d
            exp_terms(k,j,i,ibz) = (nonlin(k,j,i,ibz) + fbz(k,j,i))/d

            r = imp_rate(k2(k,j,i))

            call eSSPIFRK3(bx_new(k,j,i), bx_tmp(k,j,i), bx(k,j,i), &
                           exp_terms(k,j,i,ibx), exp_terms0(k,j,i,ibx), &
                           -r*dt, -r*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK3(by_new(k,j,i), by_tmp(k,j,i), by(k,j,i), &
                           exp_terms(k,j,i,iby), exp_terms0(k,j,i,iby), &
                           -r*dt, -r*(2.d0/3.d0*dt), 0.d0)
            call eSSPIFRK3(bz_new(k,j,i), bz_tmp(k,j,i), bz(k,j,i), &
                           exp_terms(k,j,i,ibz), exp_terms0(k,j,i,ibz), &
                           -r*dt, -r*(2.d0/3.d0*dt), 0.d0)
          enddo
        enddo
      enddo
      !$acc end data

      ! save fields at the previous step
      !$acc data present(bx, by, bz, bx_old, by_old, bz_old)
      !$acc parallel loop collapse(3)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            bx_old(k,j,i) = bx(k,j,i)
            by_old(k,j,i) = by(k,j,i)
            bz_old(k,j,i) = bz(k,j,i)
          enddo
        enddo
      enddo
      !$acc end data

      if (driven) call save_force_old

      call dealias_fields(bx_new, by_new, bz_new)
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

      if(nonlinear) call get_nonlinear_terms(bx, by, bz, .true.)

      !$acc data present(exp_terms, exp_terms_old, exp_terms_old2, nonlin, k2, &
      !$acc              bx     , by     , bz     , &
      !$acc              bx_old , by_old , bz_old , &
      !$acc              bx_old2, by_old2, bz_old2, &
      !$acc              bx_new , by_new , bz_new , &
      !$acc              fbx    , fby    , fbz)
      !$acc parallel loop collapse(3) private(r, d)
      do i = 1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            d = 1.d0 + de2*k2(k,j,i)
            exp_terms(k,j,i,ibx) = (nonlin(k,j,i,ibx) + fbx(k,j,i))/d
            exp_terms(k,j,i,iby) = (nonlin(k,j,i,iby) + fby(k,j,i))/d
            exp_terms(k,j,i,ibz) = (nonlin(k,j,i,ibz) + fbz(k,j,i))/d

            r = imp_rate(k2(k,j,i))

            ! 1st order
            if(counter == 1) then
              call gear1(bx_new(k,j,i), bx(k,j,i), exp_terms(k,j,i,ibx), r)
              call gear1(by_new(k,j,i), by(k,j,i), exp_terms(k,j,i,iby), r)
              call gear1(bz_new(k,j,i), bz(k,j,i), exp_terms(k,j,i,ibz), r)

            ! 2nd order
            elseif(counter == 2) then
              call gear2(bx_new(k,j,i), bx(k,j,i), bx_old(k,j,i), &
                         exp_terms(k,j,i,ibx), exp_terms_old(k,j,i,ibx), r)
              call gear2(by_new(k,j,i), by(k,j,i), by_old(k,j,i), &
                         exp_terms(k,j,i,iby), exp_terms_old(k,j,i,iby), r)
              call gear2(bz_new(k,j,i), bz(k,j,i), bz_old(k,j,i), &
                         exp_terms(k,j,i,ibz), exp_terms_old(k,j,i,ibz), r)

            ! 3rd order
            else
              call gear3(bx_new(k,j,i), bx(k,j,i), bx_old(k,j,i), bx_old2(k,j,i), &
                         exp_terms(k,j,i,ibx), exp_terms_old(k,j,i,ibx), exp_terms_old2(k,j,i,ibx), r)
              call gear3(by_new(k,j,i), by(k,j,i), by_old(k,j,i), by_old2(k,j,i), &
                         exp_terms(k,j,i,iby), exp_terms_old(k,j,i,iby), exp_terms_old2(k,j,i,iby), r)
              call gear3(bz_new(k,j,i), bz(k,j,i), bz_old(k,j,i), bz_old2(k,j,i), &
                         exp_terms(k,j,i,ibz), exp_terms_old(k,j,i,ibz), exp_terms_old2(k,j,i,ibz), r)
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
            bx_old2(k,j,i) = bx_old(k,j,i)
            bx_old (k,j,i) = bx(k,j,i)

            by_old2(k,j,i) = by_old(k,j,i)
            by_old (k,j,i) = by(k,j,i)

            bz_old2(k,j,i) = bz_old(k,j,i)
            bz_old (k,j,i) = bz(k,j,i)

            do l = 1, nfields
              exp_terms_old2(k,j,i,l) = exp_terms_old(k,j,i,l)
              exp_terms_old (k,j,i,l) = exp_terms(k,j,i,l)
            enddo
          enddo
        enddo
      enddo
      !$acc end data

      if (driven) call save_force_old

      call dealias_fields(bx_new, by_new, bz_new)
    endif

    ! Remove the round-off divergence; the k = 0 mode (guide field) is
    ! left untouched since k2inv = 0 there.
    call div_free(bx_new, by_new, bz_new)

    !$acc data present(bx, by, bz, bx_new, by_new, bz_new)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          bx(k,j,i) = bx_new(k,j,i)
          by(k,j,i) = by_new(k,j,i)
          bz(k,j,i) = bz_new(k,j,i)
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
!! @brief   Set the forcing fbx, fby, fbz at the
!!          current forcing time
!!
!!          get_force writes only the stirred modes
!!          of the forced components; the projection
!!          then mixes the components, so all three
!!          are cleared first. Otherwise an unforced
!!          component would keep the projection of
!!          the previous call and accumulate.
!-----------------------------------------------!
  subroutine set_force
    use grid, only: nkx, nky_local, nkz
    use fields, only: div_free
    use force, only: fbx, fby, fbz
    use force, only: get_force, normalize_force
    implicit none
    integer :: i, j, k

    !$acc parallel loop collapse(3) present(fbx, fby, fbz)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fbx(k,j,i) = (0.d0, 0.d0)
          fby(k,j,i) = (0.d0, 0.d0)
          fbz(k,j,i) = (0.d0, 0.d0)
        enddo
      enddo
    enddo

    call get_force('bx', fbx)
    call get_force('by', fby)
    call get_force('bz', fbz)

    ! The forcing must keep div b = 0
    call div_free(fbx, fby, fbz)

    call normalize_force(fbx, fby, fbz)

  end subroutine set_force


!-----------------------------------------------!
!> @author  YK
!! @brief   fb*_old = fb* (read by the diagnostics)
!-----------------------------------------------!
  subroutine save_force_old
    use grid, only: nkx, nky_local, nkz
    use force, only: fbx, fby, fbz
    use force, only: fbx_old, fby_old, fbz_old
    implicit none
    integer :: i, j, k

    !$acc data present(fbx, fby, fbz, fbx_old, fby_old, fbz_old)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          fbx_old(k,j,i) = fbx(k,j,i)
          fby_old(k,j,i) = fby(k,j,i)
          fbz_old(k,j,i) = fbz(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine save_force_old


!-----------------------------------------------!
!> @author  YK
!! @brief   Apply the dealiasing filter
!-----------------------------------------------!
  subroutine dealias_fields(bx, by, bz)
    use grid, only: nkx, nky_local, nkz
    use dealias, only: filter
    implicit none
    complex(8), dimension(:,:,:), intent(inout) :: bx, by, bz
    integer :: i, j, k

    !$acc data present(bx, by, bz, filter)
    !$acc parallel loop collapse(3)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          bx(k,j,i) = bx(k,j,i)*filter(k,j,i)
          by(k,j,i) = by(k,j,i)*filter(k,j,i)
          bz(k,j,i) = bz(k,j,i)*filter(k,j,i)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine dealias_fields


!-----------------------------------------------!
!> @author  YK
!! @brief   Allocate fields used only here
!-----------------------------------------------!
  subroutine init_work_fields
    use grid, only: nkx, nky_local, nkz
    use grid, only: k2
    use dealias, only: filter
    use params, only: time_step_scheme, de2
    use file, only: open_output_file
    use mp, only: proc0, max_allreduce, comm_fft
    implicit none
    complex(8), allocatable, dimension(:,:,:) :: src
    integer :: i, j, k


    call allocate_advance

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    ! Persistent across the deallocate/allocate_advance cycle around the output
    if(time_step_scheme == 'gear3') then
      allocate(src(nkz, nky_local, nkx), source=(0.d0,0.d0))

      allocate(bx_old2, source=src)
      allocate(by_old2, source=src)
      allocate(bz_old2, source=src)

      allocate(exp_terms_old (nkz, nky_local, nkx, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(nkz, nky_local, nkx, nfields)); exp_terms_old2 = 0.d0

      deallocate(src)

      !$acc enter data copyin(bx_old2, by_old2, bz_old2)
      !$acc enter data copyin(exp_terms_old )
      !$acc enter data copyin(exp_terms_old2)
    endif

    ! k^2/(1 + de^2 k^2) increases with k^2, so the maximum sits at the
    ! largest retained wavenumber; the reduction is over the local slab.
    Kw = 0.d0
    !$acc parallel loop collapse(3) present(k2, filter) reduction(max:Kw)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          if (filter(k,j,i) > 0.d0) Kw = max(Kw, k2(k,j,i)/(1.d0 + de2*k2(k,j,i)))
        enddo
      enddo
    enddo
    call max_allreduce(Kw, comm=comm_fft)
    if(proc0) print '(" Max. k^2/(1 + de^2 k^2) (whistler CFL) = ", es12.4e3)', Kw

    if(proc0) call open_output_file (cfl_unit, 'cfl.dat')

  end subroutine init_work_fields


!-----------------------------------------------!
!> @author  YK
!! @brief   Calculate nonlinear terms via
!!          1. b* and u in Fourier space
!!          2. Inverse FFT
!!          3. b* x u in real space
!!          4. Forward FFT and curl
!!
!!          CFL (dt_reset): advection by u and the
!!          whistler frequency, bounded by
!!            omega <= max|b*| Kw,
!!          against the imaginary-axis stability
!!          limit C_s of the explicit scheme
!!          (sqrt(3) for 3-stage RK3; 0.63 for the
!!          3rd-order extrapolated Gear scheme).
!!          max|b*| includes the guide field, so the
!!          limit is finite even for zero fluctuations.
!-----------------------------------------------!
  subroutine get_nonlinear_terms(bx, by, bz, dt_reset)
    use grid, only: dlx, dly, dlz
    use grid, only: nlx_local, nly, nlz
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms
    use advance_common, only: dt_adjust_while_running
    use params, only: time_step_scheme
    use model_specific, only: ibsx, ibsy, ibsz, iux, iuy, iuz
    use model_specific, only: fields_to_real, nonlinear_to_spectral
    implicit none
    complex(8), dimension (:,:,:), intent(in) :: bx, by, bz

    logical, intent(in) :: dt_reset

    integer :: i, j, k
    real   (8) :: vx, vy, vz, bm, c_s, rate
    real   (8) :: dt_cfl, dt_digit, dt_old

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    ! 1. b* and u in real space
    call fields_to_real(bx, by, bz, w, w_r)

    ! (get max speeds for dt reset)
    if(dt_reset) then
      vx = 0.d0; vy = 0.d0; vz = 0.d0; bm = 0.d0
      ! k = nlz + 1 : nlz_padded is FFT padding, not physical points
      !$acc data present(w_r)
      !$acc parallel loop collapse(3) reduction(max:vx, vy, vz, bm)
      do i = 1, nlx_local
        do j = 1, nly
          do k = 1, nlz
            vx = max(vx, abs(w_r(k,j,i,iux)))
            vy = max(vy, abs(w_r(k,j,i,iuy)))
            vz = max(vz, abs(w_r(k,j,i,iuz)))
            bm = max(bm, sqrt(w_r(k,j,i,ibsx)**2 + w_r(k,j,i,ibsy)**2 + w_r(k,j,i,ibsz)**2))
          end do
        end do
      end do
      !$acc end data

      ! max is duplication-invariant: reduce on world so every rank gets
      ! the same dt and the step counts stay in lockstep.
      call max_allreduce(vx)
      call max_allreduce(vy)
      call max_allreduce(vz)
      call max_allreduce(bm)

      if (time_step_scheme == 'gear3') then
        c_s = 0.63d0
      else
        c_s = sqrt(3.d0)
      endif

      rate = max(vx/dlx, vy/dly, vz/dlz, bm*Kw/c_s)
      if (rate > 0.d0) then
        dt_cfl = cfl/rate
      else
        dt_cfl = huge(1.d0)
      endif

      if(proc0) then
        write (unit=cfl_unit, fmt="(100es30.21)") tt, dt_cfl, vx, vy, vz, bm*Kw
        call flush(cfl_unit)
      endif

      if(dt_cfl < dt) then
        if(proc0) then
          print *
          write (*, '("dt is decreased from ", es12.4e3)', advance='no') dt
        endif

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

    ! 2. b* x u in real space, back to Fourier space, curl
    call nonlinear_to_spectral(w_r, nonlin_r, nonlin)

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @brief   Resistive damping rate of b,
!!          eta (k^2/k2_max)^eta_exp/(1 + de^2 k^2)
!!
!!          Zero at k = 0 so that the guide field
!!          stays constant for any eta_exp.
!-----------------------------------------------!
  function imp_rate(k2) result(r)
    !$acc routine seq
    use grid, only: k2_max
    use params, only: eta, eta_exp, de2
    implicit none
    real(8), intent(in) :: k2
    real(8) :: r

    if (k2 == 0.d0) then
      r = 0.d0
    else
      r = eta*(k2/k2_max)**eta_exp/(1.d0 + de2*k2)
    endif

  end function imp_rate


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

    allocate(bx_new, source=src)
    allocate(by_new, source=src)
    allocate(bz_new, source=src)

    allocate(w        (nkz       , nky_local, nkx      , nbtran ), source=(0.d0, 0.d0))
    allocate(nonlin   (nkz       , nky_local, nkx      , nfields), source=(0.d0, 0.d0))
    allocate(w_r      (nlz_padded, nly      , nlx_local, nbtran ), source=0.d0)
    allocate(nonlin_r (nlz_padded, nly      , nlx_local, nftran ), source=0.d0)
    allocate(exp_terms(nkz       , nky_local, nkx      , nfields), source=(0.d0, 0.d0))

    !$acc enter data copyin(bx_new, by_new, bz_new)
    !$acc enter data copyin(w        )
    !$acc enter data copyin(nonlin   )
    !$acc enter data copyin(w_r      )
    !$acc enter data copyin(nonlin_r )
    !$acc enter data copyin(exp_terms)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      allocate(bx_tmp, source=src)
      allocate(by_tmp, source=src)
      allocate(bz_tmp, source=src)

      allocate(exp_terms0(nkz, nky_local, nkx, nfields)); exp_terms0  = 0.d0

      !$acc enter data copyin(bx_tmp, by_tmp, bz_tmp)
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

    !$acc exit data delete(bx_new, by_new, bz_new)
    !$acc exit data delete(w        )
    !$acc exit data delete(nonlin   )
    !$acc exit data delete(w_r      )
    !$acc exit data delete(nonlin_r )
    !$acc exit data delete(exp_terms)

    deallocate(bx_new, by_new, bz_new)
    deallocate(w        )
    deallocate(nonlin   )
    deallocate(w_r      )
    deallocate(nonlin_r )
    deallocate(exp_terms)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      !$acc exit data delete(bx_tmp, by_tmp, bz_tmp)
      !$acc exit data delete(exp_terms0)
      deallocate(bx_tmp, by_tmp, bz_tmp)
      deallocate(exp_terms0)
    endif

  end subroutine deallocate_advance


end module advance
