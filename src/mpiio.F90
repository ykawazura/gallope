!=======================================================================
! This is based on the 2DECOMP&FFT library
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!=======================================================================

module mpiio
  use MPI
  implicit none

  private        ! Make everything private unless declared public

  public :: mpiio_write_one, mpiio_read_one, mpiio_write_var
  public :: mpiio_write_var_2d

  interface mpiio_write_one
    module procedure mpiio_write_one_complex
    module procedure mpiio_write_one_complex_4d
    module procedure mpiio_write_one_real
  end interface mpiio_write_one

  interface mpiio_read_one
     module procedure mpiio_read_one_complex
     module procedure mpiio_read_one_complex_4d
  end interface mpiio_read_one

  interface mpiio_write_var
     module procedure mpiio_write_var_complex
  end interface mpiio_write_var

  interface mpiio_write_var_2d
     module procedure mpiio_write_var_2d_real
  end interface mpiio_write_var_2d

  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 3D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_write_one_complex(var, sizes, subsizes, starts, filename, comm)

    implicit none

    complex(8), dimension(:,:,:), intent(IN) :: var
    integer, dimension(3), intent(IN) :: sizes, subsizes, starts
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: comm

    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer :: ierror, newtype, fh, data_type

    data_type = complex_type

#include "mpiio_write_one.F90"
    
    return
  end subroutine mpiio_write_one_complex

  subroutine mpiio_write_one_real(var, sizes, subsizes, starts, filename, comm)

    implicit none

    real(8), dimension(:,:,:), intent(IN) :: var
    integer, dimension(3), intent(IN) :: sizes, subsizes, starts
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: comm

    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer :: ierror, newtype, fh, data_type

    data_type = real_type

#include "mpiio_write_one.F90"

    return
  end subroutine mpiio_write_one_real

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 4D array to a file.
  ! Used for the Hermite-moment field g(nkz,nky_local,nkx,nm_local): the
  ! caller passes the interior slice g(:,:,:,1:nm_local) (contiguous, no ghosts).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_write_one_complex_4d(var, sizes, subsizes, starts, filename, comm)

    implicit none

    complex(8), dimension(:,:,:,:), intent(IN) :: var
    integer, dimension(4), intent(IN) :: sizes, subsizes, starts
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: comm

    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer :: ierror, newtype, fh, data_type

    data_type = complex_type

    call MPI_TYPE_CREATE_SUBARRAY(4, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype, ierror)
    call MPI_FILE_OPEN(comm, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh, filesize, ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh, disp, data_type, &
         newtype, 'native', MPI_INFO_NULL, ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3)*subsizes(4), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh, ierror)
    call MPI_TYPE_FREE(newtype, ierror)

    return
  end subroutine mpiio_write_one_complex_4d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to read from a file a single 3D array
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_read_one_complex(var, sizes, subsizes, starts, filename, comm)

    implicit none

    complex(8), dimension(:,:,:), intent(INOUT) :: var
    integer, dimension(3), intent(IN) :: sizes, subsizes, starts
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: comm

    ! integer(kind=MPI_OFFSET_KIND) :: disp
    integer :: ierror, newtype, fh, data_type
    
    data_type = complex_type

#include "mpiio_read_one.F90"

    return
  end subroutine mpiio_read_one_complex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to read from a file a single 4D array.
  ! Companion of mpiio_write_one_complex_4d (g restart read).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_read_one_complex_4d(var, sizes, subsizes, starts, filename, comm)

    implicit none

    complex(8), dimension(:,:,:,:), intent(INOUT) :: var
    integer, dimension(4), intent(IN) :: sizes, subsizes, starts
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: comm

    integer :: ierror, newtype, fh, data_type

    data_type = complex_type

    call MPI_TYPE_CREATE_SUBARRAY(4, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype, ierror)
    call MPI_FILE_OPEN(comm, filename, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, data_type, &
         newtype, 'native', MPI_INFO_NULL, ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3)*subsizes(4), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh, ierror)
    call MPI_TYPE_FREE(newtype, ierror)

    return
  end subroutine mpiio_read_one_complex_4d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the writing
  !  operation to prepare the writing of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_write_var_complex(fh, disp, sizes, subsizes, starts, var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, dimension(3), intent(IN) :: sizes, subsizes, starts
    complex(8), dimension(:,:,:), intent(IN) :: var

    integer :: ierror, newtype, data_type, bytes

    data_type = complex_type
    call MPI_TYPE_SIZE(data_type,bytes,ierror)

#include "mpiio_write_var.F90"

    return
  end subroutine mpiio_write_var_complex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the writing
  !  operation to prepare the writing of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_write_var_2d_real(fh, disp, sizes, subsizes, starts, var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, dimension(2), intent(IN) :: sizes, subsizes, starts
    real(8), dimension(:,:), intent(IN) :: var

    integer :: ierror, newtype, data_type, bytes

    data_type = real_type
    call MPI_TYPE_SIZE(data_type,bytes,ierror)

#include "mpiio_write_var_2d.F90"

    return
  end subroutine mpiio_write_var_2d_real

end module mpiio
