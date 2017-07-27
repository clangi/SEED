!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!                                                                          !
!    Copyright (C) 2014, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger                !
!                                                                          !
!    Website: http://sourceforge.net/projects/campari/                     !
!                                                                          !
!    This tool is free software: you can redistribute it and/or modify     !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    The tool is distributed in the hope that it will be useful,           !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You can find a copy of the license at <http://www.gnu.org/licenses/>. !
!--------------------------------------------------------------------------!
! AUTHORSHIP INFO:                                                         !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN AUTHOR:   Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
module interfaces
!
  interface
    function handwrapped_getpid() bind(C, name='getpid')
      use, INTRINSIC:: ISO_C_BINDING
      implicit none
      integer(c_int) :: handwrapped_getpid
    end function handwrapped_getpid
  end interface
!
  interface
    function handwrapped_gethostname(hnam,hnaml) bind(C, name='gethostname')
      use, INTRINSIC:: ISO_C_BINDING
      implicit none
      character(len=1,kind=c_char), dimension(*), intent(out) :: hnam
      integer(c_size_t):: hnaml
!      type(C_PTR):: hnam
      integer(c_int) :: handwrapped_gethostname
    end function handwrapped_gethostname
  end interface
!
end module interfaces
!
!-------------------------------------------------------------------------
!
subroutine strlims(str,f,l)
!
  implicit none
!
  integer f,l,last,i,leng
  character(*) str
  character nix
!
  l = 0
  nix = char(0)
  leng = len(str)
  last = leng
  f = last+1
  do i=1,leng
    if (str(i:i).gt.' ') then
      f = i
      exit
    end if
  end do
  do i=1,leng
    if (str(i:i).eq.nix) then
      last = i - 1
      exit
    end if
  end do
!
  do i=last,1,-1
    if (str(i:i).gt.' ') then
      l = i
      exit
    end if
  end do
!
  if (f.gt.i) then
    f = 1
    i = 1
  end if
!
end subroutine strlims
!
!-------------------------------------------------------------------------------------
!
subroutine int2str(ii,string,strsz)
!
  implicit none
!
  integer mx(0:16),leng,i,ii,strsz,minsz,cpy,sumi
  character it(0:9)
  character(*) string
  logical just_r
  data it /'0','1','2','3','4','5','6','7','8','9'/
!
  if (strsz.le.0) then
    just_r = .true.
    strsz = 1
  else
    just_r = .false.
  end if
  minsz = strsz
  leng = len(string)
!
! invert sign if necessary (technically this routines only
! processes unsigned integers)
  if (ii .lt. 0) then
    i = -ii
  end if
!
! this is standard: the mx(...) are integers of course
!
  sumi = 0
  do i=16,0,-1
    mx(i) = (ii-sumi)/(10.0**i)
    if (i.gt.0) then
      sumi = sumi + (10.0**i) * mx(i)
    end if
  end do
!
! now we have the digits of our string in integer form (0 through 9)
!
! check for integer-size
!
  if ((mx(10).gt.0).OR.(mx(11).gt.0).OR.(mx(12).gt.0).OR.&
 &    (mx(13).gt.0).OR.&
 &    (mx(14).gt.0).OR.(mx(15).gt.0).OR.(mx(16).gt.0)) then
    write(*,*) 'Warning. Possibly bad result from int2str(...) (inte&
 &ger overflow).'
  end if
!
! find the final string length
!
  strsz = 1
  do i=16,1,-1
    if (mx(i).ne.0) then
      strsz = i+1
      exit
    end if
  end do
!
! correct if too large or small
  strsz = min(strsz,leng)
  strsz = max(strsz,minsz)
!
  if (strsz.gt.17) then
    write(*,*) 'Warning. Possibly bad result from int2str(...) (digi&
 &t overflow).'
  end if
!
! convert individual digits to a string of numeric characters
!
  do i=1,strsz
    string(i:i) = it(mx(strsz-i))
  end do
!
! left/right-justification
!
  if (just_r.EQV..true.) then
    do i=strsz,1,-1
      cpy = leng-strsz+i
      string(cpy:cpy) = string(i:i)
    end do
    do i=1,leng-strsz
      string(i:i) = ' '
    end do
  else
    do i = strsz+1,leng
      string(i:i) = ' '
    end do
  end if
!
end subroutine int2str
!
!-------------------------------------------------------------------------
!
program MPI_scriptwrapper
! program MPI_scriptwrapper()
!
  use mpi
  use interfaces
  use, INTRINSIC:: ISO_C_BINDING
!
  implicit none
!
  integer masterrank,myrank,mpi_nodes,ierr,j,k,t1,t2,t3,t4,t5,t6,ck,mstatus(MPI_STATUS_SIZE),i,endtag,argstat
  integer std_f_out,std_f_in
  parameter (std_f_out=6)
  parameter (std_f_in=5)
  character(len=1,kind=c_char), dimension(1001):: honam(0:1000)
  character(1000) hnam
  character(2000) sname
  integer(c_int) hnerr
  integer(c_size_t) hnml
!
  integer MAXARGLEN
  parameter (MAXARGLEN=150)
!
  type t_args
    character(MAXARGLEN):: it
    logical listit
    integer leng
  end type t_args
! 
  integer nargs
  type(t_args), ALLOCATABLE:: args(:)
!
  character(9) nostr
!
! some initialization
  masterrank = 0
!
! let's start up the universe
  call MPI_INIT(ierr)
!
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nodes,ierr)
  hnml = 1000
  honam(:) = ' '
  hnerr = handwrapped_gethostname(honam(0:(999)),hnml)
  do i=1,1000
    hnam(i:i) = honam(i-1)
  end do
  call strlims(hnam,t1,t2)
  if (t2.gt.t1) then
    if (hnam(t2:t2).eq.C_NULL_CHAR) t2 = t2 - 1
  end if
!
  write(std_f_out,*) 'Process ',myrank+1,' of ',mpi_nodes,' is on ',hnam(t1:t2),'.'
!
! the FORTRAN intrinsics can handle this well in modern Fortran
! note the zero-index argument is the program call itself and is
! not counted in COMMAND_ARGUMENT_COUNT
  nargs = COMMAND_ARGUMENT_COUNT()
!
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
  if (nargs.eq.0) then
    if (myrank.eq.0) then
      k = std_f_out
      write(k,*) 'USAGE: This script requires a mandatory argument, which is the base name (full path) of the script to be &
 &executed.'
      write(k,*) 'These names will be appended with the process number.'
      write(k,*)
      write(k,*) 'EXAMPLE: ./MPI_scriptwrapper.exe /users/yourname/tester'
      write(k,*)
      write(k,*) 'This will attempt to execute scripts /users/yourname/tester1.sh, /users/yourname/tester2.sh, ... in the current &
 &working directory (note suffix).'
      write(k,*) 'These scripts have to exist, be executable, and contain all relevant commands (including change of directories).'
    end if
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    stop 98
  end if
!
  allocate(args(nargs))
!
  do i=1,nargs
    call GET_COMMAND_ARGUMENT(i,args(i)%it,args(i)%leng,argstat)
  end do
!
! construct script name
  k = 0 
  call int2str(myrank+1,nostr,k)
  call strlims(nostr,t3,t4)
  sname = args(1)%it(1:args(1)%leng)//nostr(t3:t4)//'.sh'
  call strlims(sname,t5,t6)
  write(std_f_out,*) 'Process ',myrank+1,' of ',mpi_nodes,' will try to execute ',sname(t5:t6),' on ',hnam(t1:t2),'.'
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
  call EXECUTE_COMMAND_LINE(sname(t5:t6),WAIT=.true.,EXITSTAT=k,CMDSTAT=ck)
!
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (ck.ne.0) then
    write(std_f_out,*) 'Process ',myrank+1,' of ',mpi_nodes,' failed to execute ',sname(t5:t6),' on ',hnam(t1:t2),'.'
  else
    write(std_f_out,*) 'Process ',myrank+1,' of ',mpi_nodes,' finished executing ',sname(t5:t6),' on ',hnam(t1:t2),' with a status &
 &flag of ',k,' (this applies to the last command in the script).'
  end if
!
! close shop
  endtag = 777
!
  if (myrank.eq.masterrank) then
!
    do i=1,mpi_nodes-1
!      write(*,*) 'sending end to ',i+1
      call MPI_Send(j,1,MPI_INTEGER,i,endtag,MPI_COMM_WORLD,ierr)
    end do
!
  else
!
    call MPI_RECV(j,1,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.endtag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from master. Expected end signal (',endtag,').'
      stop 98
    end if
!
  end if
!
  call MPI_FINALIZE(ierr)
!
end

