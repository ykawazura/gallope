! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2014 Dmitry Pekurovsky
!    Copyright (C) 2006-2014 University of California
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!
!----------------------------------------------------------------------------

! This sample program illustrates the
! use of P3DFFT library for highly scalable parallel 3D FFT.
!
! This program initializes a 3D array with random numbers, then
! performs forward transform, backward transform, and checks that
! the results are correct, namely the same as in the start except
! for a normalization factor. It can be used both as a correctness
! test and for timing the library functions.
!
! The program expects 'stdin' file in the working directory, with
! a single line of numbers : Nx,Ny,Nz,Ndim,Nrep. Here Nx,Ny,Nz
! are box dimensions, Ndim is the dimentionality of processor grid
! (1 or 2), and Nrep is the number of repititions. Optionally
! a file named 'dims' can also be provided to guide in the choice
! of processor geometry in case of 2D decomposition. It should contain
! two numbers in a line, with their product equal to the total number
! of tasks. Otherwise processor grid geometry is chosen automatically.
! For better performance, experiment with this setting, varying
! iproc and jproc. In many cases, minimizing iproc gives best results.
! Setting it to 1 corresponds to one-dimensional decomposition.
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu

program fft3d

   use p3dfft
   implicit none
   include 'mpif.h'

   integer i,n,nx,ny,nz
   integer m,x,y,z
   integer fstatus
   logical flg_inplace
   integer params(4), iounit
   logical have_params

   real(p3dfft_type), dimension(:,:,:),  allocatable :: BEG,FIN,CP
   complex(p3dfft_type), dimension(:,:,:),  allocatable :: AEND
   real(p3dfft_type) diff,cdiff,ccdiff,ans

   integer(i8) Ntot
   real(p3dfft_type) factor
   real(r8) rtime, rtime0,Nglob,prec
   integer ierr,nu,ndim,dims(2),nproc,proc_id
   integer istart(3),iend(3),isize(3)
   integer fstart(3),fend(3),fsize(3)
   integer iproc,jproc,nxc,nyc,nzc
   logical iex
   integer memsize(3)

   integer :: hostname_len
   character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
   character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: all_hostnames(:)
   integer :: unique_nodes
   logical :: is_unique
   integer :: j

   call MPI_INIT (ierr)
   call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
   call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

!vvvvv   Setting parameters   vvvvv
      nx   = 2048
      ny   = 2048
      nz   = 2048
      ndim = 2
      n    = 10
      ! A file "params.in" in the run directory, holding "nx ny nz nrepeat" on a
      ! single line, overrides the grid above (ndim stays 2, as in every run
      ! reported so far).  Same convention as the two cuFFTMp benchmarks, so the
      ! weak-scaling series does not need one binary per shape.  Absent the file
      ! the run is bit-for-bit the original one.
      if (proc_id == 0) then
         inquire(file='params.in', exist=have_params)
         if (have_params) then
            open(newunit=iounit, file='params.in', status='old', action='read')
            read(iounit, *) nx, ny, nz, n
            close(iounit)
         endif
      endif
      params = (/ nx, ny, nz, n /)
      call MPI_BCAST(params, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      nx = params(1)
      ny = params(2)
      nz = params(3)
      n  = params(4)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   timers = 0

   if (proc_id.eq.0) then
      if(p3dfft_type .eq. 4) then
         print *,'Single precision version'
      else if(p3dfft_type .eq. 8) then
         print *,'Double precision version'
      endif
   endif


   ! Get number of nodes
   if (proc_id.eq.0) then
      print *, 'P3DFFT test, random input'
      write(*, '(" nx = ", I5, ", ny = ", I5, ", nz = ", I5, ", nproc = ", I6, ", ndim = ", I5, ", nrepeat = ", I5)') &
            nx, ny, nz, nproc, ndim, n
   endif

   if (proc_id == 0) then
      allocate(all_hostnames(nproc))
   endif
   
   call MPI_GET_PROCESSOR_NAME(hostname, hostname_len, ierr)
   call MPI_Gather(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                   all_hostnames, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                   0, MPI_COMM_WORLD, ierr)
   
   if (proc_id == 0) then
      ! Count unique node number
      unique_nodes = 1
      do i = 2, nproc
         is_unique = .true.
         do j = 1, i-1
            if (trim(all_hostnames(i)) == trim(all_hostnames(j))) then
               is_unique = .false.
               exit
            endif
         enddo
         if (is_unique) unique_nodes = unique_nodes + 1
      enddo
      
      write(*, '(" ==========================================")') 
      write(*, '(" Total MPI processes: ", I6)') nproc
      write(*, '(" Unique nodes:        ", I6)') unique_nodes
      
      ! Show list of hostnames
      !do i = 1, nproc
      !   write(*, '(" Rank = ", I6, " on ", A)') i - 1, trim(all_hostnames(i))
      !enddo
      write(*, '(" ==========================================")') 
      
      deallocate(all_hostnames)
   endif

   !    nproc is devided into a iproc x jproc stencle
   !

   if(ndim .eq. 1) then
      dims(1) = 1
      dims(2) = nproc
   else if(ndim .eq. 2) then
      inquire(file='dims',exist=iex)
      if (iex) then
         if (proc_id.eq.0) print *, 'Reading proc. grid from file dims'
         open (999,file='dims')
         read (999,*) dims(1), dims(2)
         close (999)
         if(dims(1) * dims(2) .ne. nproc) then
            dims(2) = nproc / dims(1)
         endif
      else
         if (proc_id.eq.0) print *, 'Creating proc. grid with mpi_dims_create'
         dims(1) = 0
         dims(2) = 0
         call MPI_Dims_create(nproc,2,dims,ierr)
         if(dims(1) .gt. dims(2)) then
            dims(1) = dims(2)
            dims(2) = nproc / dims(1)
         endif
      endif
   endif

   iproc = dims(1)
   jproc = dims(2)

   if(proc_id .eq. 0) then
      print *,'Using processor grid ',iproc,' x ',jproc
   endif

   nxc = nx
   nyc = ny
   nzc = nz

   call p3dfft_setup (dims,nx,ny,nz,MPI_COMM_WORLD,nxc,nyc,nzc,.true.)
   call p3dfft_get_dims(istart,iend,isize,1)
   call p3dfft_get_dims(fstart,fend,fsize,2)

   !      print *,'Allocating BEG (',isize,istart,iend
   allocate (BEG(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), stat=ierr)
   if(ierr .ne. 0) then
      print *,'Error ',ierr,' allocating array BEG'
   endif
   allocate (CP(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), stat=ierr)
   if(ierr .ne. 0) then
      print *,'Error ',ierr,' allocating array CP'
   endif
   !      print *,'Allocating AEND (',fsize,fstart,fend
   allocate (AEND(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)), stat=ierr)
   if(ierr .ne. 0) then
      print *,'Error ',ierr,' allocating array AEND'
   endif
   allocate (FIN(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), stat=ierr)
   if(ierr .ne. 0) then
      print *,'Error ',ierr,' allocating array FIN'
   endif

   !
   ! initialize
   !

   ! start with x-z slabs in physical space
   !
   do z=istart(3),iend(3)
      do y=istart(2),iend(2)
         do x=istart(1),iend(1)
            call random_number(BEG(x,y,z))
            !		BEG(x,y,z) = (x-1)+2*(y-1)+3*(z-1)
         enddo
      enddo
   enddo

   CP = BEG

   !
   ! transform from physical space to wavenumber space
   ! (XgYiZj to XiYjZg)
   ! then transform back to physical space
   ! (XiYjZg to XgYiZj)
   !



   Nglob = nx * ny
   Nglob = Nglob * nz
   factor = 1.0d0/Nglob
   Ntot = fsize(1)*fsize(2)*fsize(3)
   !      if(proc_id .eq. 0) then
   !         print *,'Initial data: '
   !         call print_all_real(BEG,Ntot,proc_id,Nglob)
   !      endif

   ! Warming up
   do  m=1, 3
      ! Forward transform
      call p3dfft_ftran_r2c (BEG,AEND,'fft')
      ! Normalize
      call mult_array(AEND, Ntot,factor)
      ! Backward transform
      call p3dfft_btran_c2r (AEND,FIN,'tff')
   end do

   rtime  = 0.0
   rtime0 = 0.0
   do  m=1,n
      ! Barrier for correct timing
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      ! Forward transform
      rtime = rtime - MPI_wtime()
      call p3dfft_ftran_r2c (BEG,AEND,'fft')
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      rtime = rtime + MPI_wtime()

      ! Normalize
      call mult_array(AEND, Ntot,factor)
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      ! Backward transform
      rtime = rtime - MPI_wtime()
      call p3dfft_btran_c2r (AEND,FIN,'tff')
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      rtime = rtime + MPI_wtime()

      if(proc_id .eq. 0) then
         print *,'Iteration = ',m, ', time = ', rtime - rtime0
      endif
      rtime0 = rtime

   end do

   if (proc_id.eq.0) write(6,*)'proc_id = ', proc_id, ', cpu time per loop = ', rtime/dble(n)

   ! Free work space

   call p3dfft_clean

   ! Check results

   cdiff=0.0d0
   ccdiff = 0.0d0
   do 20 z=istart(3),iend(3)
      do 20 y=istart(2),iend(2)
         do 20 x=istart(1),iend(1)
         if(cdiff .lt. abs(CP(x,y,z)-FIN(x,y,z))) then
            cdiff = abs(CP(x,y,z)-FIN(x,y,z))
            !               print *,'x,y,z,cdiff=',x,y,z,cdiff
         endif
20 continue
   call MPI_Reduce(cdiff,ccdiff,1,p3dfft_mpireal,MPI_MAX,0, &
         MPI_COMM_WORLD,ierr)

   if(proc_id .eq. 0) then
      if(p3dfft_type .eq. 8) then
         prec = 1e-14
      else
         prec = 1e-5
      endif
      if(ccdiff .gt. prec * Nglob*0.25) then
         print *,'Results are incorrect'
      else
         print *,'Results are correct'
      endif
      write (6,*) 'max diff =',ccdiff
   endif

   call MPI_FINALIZE (ierr)

contains
   !=========================================================

   subroutine mult_array(X,nar,f)

      use p3dfft

      integer(i8) nar,i
      complex(p3dfft_type) X(nar)
      real(p3dfft_type) f

      do i=1,nar
         X(i) = X(i) * f
      enddo

      return
   end subroutine

   !=========================================================
   ! Translate one-dimensional index into three dimensions,
   !    print out significantly non-zero values
   !
   subroutine print_all(Ar,Nar,proc_id,Nglob)

      use p3dfft

      integer x,y,z,proc_id
      integer(i8) i,Nar
      real(r8) Nglob
      complex(p3dfft_type) Ar(1,1,*)
      integer Fstart(3),Fend(3),Fsize(3)

      call p3dfft_get_dims(Fstart,Fend,Fsize,2)
      if(proc_id .eq. 0) then
         print *,'Dimensions:', Fsize
      endif
      do i=1,Nar
         if(abs(Ar(1,1,i)) .gt. 1.25e-3 * Nglob) then
            z = (i-1)/(Fsize(1)*Fsize(2))
            y = (i-1 - z * Fsize(1)*Fsize(2))/Fsize(1)
            x = i-1-z*Fsize(1)*Fsize(2) - y*Fsize(1)
            print *,proc_id,': (',x+Fstart(1),y+Fstart(2),z+Fstart(3),') ',Ar(1,1,i)
         endif
      enddo

      return
   end subroutine

   !=========================================================
   ! Translate one-dimensional index into three dimensions,
   !    print out significantly non-zero values
   !
   subroutine print_all_real(Ar,Nar,proc_id,Nglob)

      use p3dfft

      integer x,y,z,proc_id
      integer(i8) i,Nar
      real(r8) Nglob
      real(p3dfft_type) Ar(1,1,*)
      integer Fstart(3),Fend(3),Fsize(3)

      call p3dfft_get_dims(Fstart,Fend,Fsize,1)
      do i=1,Nar
         if(abs(Ar(1,1,i)) .gt. Nglob *1.25e-8) then
            z = (i-1)/(Fsize(1)*Fsize(2))
            y = (i-1 - z * Fsize(1)*Fsize(2))/Fsize(1)
            x = i-1-z*Fsize(1)*Fsize(2) - y*Fsize(1)
            print *,'(',x+Fstart(1),y+Fstart(2),z+Fstart(3),') ',Ar(1,1,i)
         endif
      enddo

      return
   end subroutine

end
