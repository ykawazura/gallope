!-----------------------------------------------!
!> @author  YK
!! @brief   Useful functions
!-----------------------------------------------!
module utils
  implicit none

  public :: cabs2
  public :: ranf
  public :: curl
  public :: check_floor

  private

  interface cabs2
    module procedure cabs2_scaler
    module procedure cabs2_array1
    module procedure cabs2_array2
    module procedure cabs2_array3
  end interface cabs2


contains


!-----------------------------------------------!
!> @author  YK
!! @brief   Get square of complex
!!          At Fugaku, sometimes abs(z) gives NaN,
!!          but z*conjg(z) is OK
!-----------------------------------------------!
  function cabs2_scaler(z) result(r)
    !$acc routine seq
    implicit none
    complex(8), intent(in) :: z
    real(8) :: r

    r = dble(z*conjg(z))

  end function cabs2_scaler

  function cabs2_array1(z) result(r)
    !$acc routine seq
    complex(8), intent(in) :: z(:)
    real(8) :: r(size(z))

    r = dble(z*conjg(z))
  end function cabs2_array1

  function cabs2_array2(z) result(r)
    !$acc routine seq
    complex(8), intent(in) :: z(:, :)
    real(8) :: r(size(z, 1), size(z, 2))

    r = dble(z*conjg(z))
  end function cabs2_array2

  function cabs2_array3(z) result(r)
    !$acc routine seq
    complex(8), intent(in) :: z(:, :, :)
    real(8) :: r(size(z, 1), size(z, 2), size(z, 3))

    r = dble(z*conjg(z))
  end function cabs2_array3


!-----------------------------------------------!
!> @author  YK
!! @brief   Get random number
!-----------------------------------------------!
  function ranf()
    use mp, only: proc_id
    implicit none
    integer :: seedsize, clock
    integer, allocatable :: seed(:)
    integer :: i, time(8)
    real(8) :: ranf
    
    ! Get the size of the seed array
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    
    ! Get current time components
    call date_and_time(values=time)
    
    ! Get system clock value
    call system_clock(count=clock)
    
    ! Create unique seed for each process
    do i = 1, seedsize
      ! Combine time components and clock with process ID
      seed(i) = clock + time(7) + time(8) + &
                (proc_id + 1) * 1000 + i * 100
    end do
    
    ! Set the random seed
    call random_seed(put=seed)
    
    ! Generate random number
    call random_number(ranf)
    
    deallocate(seed)

  end function ranf


!-----------------------------------------------!
!> @author  YK
!! @brief   Get curl
!-----------------------------------------------!
  subroutine curl(fx, fy, fz, curlfx, curlfy, curlfz)
    use grid, only: nkx, nky_local, nkz
    use grid, only: kx, ky, kz
    use params, only: zi
    implicit none
    integer :: i, j, k
    complex(8), dimension (:,:,:), intent(in ) :: fx, fy, fz
    complex(8), dimension (:,:,:), intent(out) :: curlfx, curlfy, curlfz

    !$acc parallel loop collapse(3) present(fx, fy, fz, curlfx, curlfy, curlfz)
    do i = 1, nkx
      do j = 1, nky_local
        do k = 1, nkz
          curlfx(k, j, i) = zi*(ky(j)*fz(k, j, i) - kz(k)*fy(k, j, i))
          curlfy(k, j, i) = zi*(kz(k)*fx(k, j, i) - kx(i)*fz(k, j, i))
          curlfz(k, j, i) = zi*(kx(i)*fy(k, j, i) - ky(j)*fx(k, j, i))
        enddo
      enddo
    enddo

  end subroutine curl


!-----------------------------------------------!
!> @author  YK
!! @brief   Check if the value is less than the floor value
!-----------------------------------------------!
  subroutine check_floor(u, u_min)
    use grid, only: nlx_local, nly, nlz_padded, ntot
    use mp, only: sum_allreduce
    use cuFFTmp, only: btran_c2r
    implicit none
    complex(8), dimension (:,:,:), intent(in) :: u
    real(8), intent(in) :: u_min
    real(8):: u_avg
    real(8), allocatable, dimension(:,:,:) :: u_r
    integer :: i, j, k

    allocate(u_r(nlz_padded, nly, nlx_local), source=0.d0)

    call btran_c2r(u, u_r)

    u_avg = sum(u_r)
    call sum_allreduce(u_avg)
    u_avg = u_avg/ntot

    do i = 1, nlx_local
      do j = 1, nly
        do k = 1, nlz_padded
          if(u_r(j, k, i) <= 0.d0) print *, '!!!  CAUTION negaive value  !!!'
          if(u_r(j, k, i) < u_min*u_avg) u_r(j, k, i) = u_min*u_avg
        enddo
      enddo
    enddo

    deallocate(u_r)

  end subroutine check_floor

end module utils

