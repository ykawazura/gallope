!-----------------------------------------------!
!> @author  YK
!! @brief   Dealiasing setting
!-----------------------------------------------!
module dealias
  implicit none

  public  init_filter, finish_filter

  ! de-aliasing filter; 2/3 or Hou-Li depending on input
  real(8), allocatable, dimension(:,:,:) :: filter 
  !$acc declare create(filter)

contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of filter
!-----------------------------------------------!
  subroutine init_filter
    use grid, only: nkx, nky_local, nkz
    use grid, only: kx, ky, kz
    use grid, only: kx_max, ky_max, kz_max
    implicit none
    integer :: i, j, k

    allocate(filter(nkz, nky_local, nkx))

    ! dealiasing filter
    filter = 1.d0

    do i = 1, nkx
      if(abs(kx(i)) >= kx_max*2.d0/3.d0) then
        filter(:, :, i) = 0.d0
      endif
    enddo
    do j = 1, nky_local
      if(abs(ky(j)) >= ky_max*2.d0/3.d0) then
        filter(:, j, :) = 0.d0
      endif
    enddo
    do k = 1, nkz
      if(abs(kz(k)) >= kz_max*2.d0/3.d0) then
        filter(k, :, :) = 0.d0
      endif
    enddo

    !$acc update device(filter)

  end subroutine init_filter


!-----------------------------------------------!
!> @author  YK
!! @brief   Initialization of filter
!-----------------------------------------------!
  subroutine finish_filter
    implicit none

    deallocate(filter)

  end subroutine finish_filter

end module dealias

