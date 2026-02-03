!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Model specific subroutines
!-----------------------------------------------!
module model_specific
  use force_common
  implicit none

  public init_model_specific

contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialize model specific subroutines
!-----------------------------------------------!
  subroutine init_model_specific
    use mp, only: proc0
    use grid, only: lx, ly, lz, kx, ky, kz
    use force, only :init_force, alloc_force
    use shearing_box, only: init_shearing_time
    use params, only: shear, pi
    use fields, only: bx, by, bz
    use grid, only: nkx, nky_local, nkz
    use mp, only: broadcast
    implicit none
    real(8) :: bx0, by0, bz0 ! mean magnetic field
    real(8) :: lmd_MRI
    integer :: i, j, k

    call init_shearing_time

    call init_force
    call alloc_force

    ! show minimum k in unit of omega0/vA
    if(shear) then
      !$acc parallel loop collapse(3) reduction(+:bx0, by0, bz0)
      do i =1, nkx
        do j = 1, nky_local
          do k = 1, nkz
            if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
              bx0 = abs(bx(k,j,i))
              by0 = abs(by(k,j,i))
              bz0 = abs(bz(k,j,i))
            endif
          end do
        end do
      end do
      call broadcast(bx0)
      call broadcast(by0)
      call broadcast(bz0)

      if(proc0) then
        write(*, "('--------------------------------------')", advance='no')
        write(*, "('--------------------------------------')")

        write(*, "('kx0*vA_x/Omega = ', f5.2, &
              & ',  ky0*vA_y/Omega = ', f5.2, &
              & ',  kz0*vA_z/Omega = ', f5.2)") &
              & bx0*kx(2), by0*ky(2), bz0*kz(2)
        write(*, "('  (vA is computed by the initial mean magnetic fields.)')")

        lmd_MRI = 2.d0*pi*bz0
        write(*, "(' Lx/lmd_MRI = ', f5.2, &
              & ',   Ly/lmd_MRI = ', f5.2, &
              & ',   Lz/lmd_MRI = ', f5.2)") &
              & lx/lmd_MRI, ly/lmd_MRI, lz/lmd_MRI
        write(*, "('  (lmd_MRI = 2pi*vA_z/Omega.)')")

        write(*, "('--------------------------------------')", advance='no')
        write(*, "('--------------------------------------')")
      endif
    endif

  end subroutine init_model_specific
end module model_specific

