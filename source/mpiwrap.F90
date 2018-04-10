module mpiWrap

  implicit none

# ifdef USE_MPI
#   include "mpif.h"
# else
    integer :: mpi_COMM_WORLD
# endif

  type mpi_info
     integer :: rank        ! global rank
     integer :: csize        ! global size  
     integer :: comm
     integer :: root
     integer, pointer :: work1(:)
     integer, pointer :: work2(:)
  end type mpi_info

  type(mpi_info)  :: sim

  interface mpiwrap_bcast
     module procedure mpiWrap_bcastI, mpiWrap_bcast1C, mpiWrap_bcast2DP,  &
       mpiWrap_bcastL, mpiwrap_bcast2I, mpiwrap_bcast1I, mpiwrap_bcast1DP,&
       mpiWrap_bcastDP, mpiWrap_bcastC, mpiwrap_bcast1L
  end interface

  interface mpiwrap_sum
     module procedure mpiWrap_sumDP, mpiWrap_sum2DP, mpiWrap_sum3DP, mpiWrap_sum1DPC, mpiWrap_sum1DP, mpiWrap_sumDPC
  end interface

  interface mpiwrap_collect
     module procedure mpiWrap_collectDP, mpiWrap_collect2DP
  end interface

  interface mpiwrap_max
     module procedure mpiWrap_maxDP, mpiWrap_maxI
  end interface

  interface mpiwrap_send
     module procedure mpiWrap_send_oneI, mpiWrap_send_array2DP, &
       mpiWrap_send_array2ZC, mpiwrap_send_arrayDP
  end interface

  interface mpiwrap_recv
     module procedure mpiWrap_recv_oneI, mpiWrap_recv_array2DP, &
       mpiWrap_recv_array2ZC, mpiwrap_send_arrayDP
  end interface

contains

  subroutine mpiWrap_Init()
    !
    ! Purpose 
    ! =======
    !   This subroutine initializes mpi procedure.
    !
    ! ... Local variables ...
#   ifdef USE_MPI
    integer :: ierr, rc
#   endif 

    ! ... Executable ... 
    
#   ifdef USE_MPI
       call mpi_INIT(ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error starting mpi program. Terminating.'
          call mpi_ABORT(mpi_COMM_WORLD, rc, ierr)
       end if
#   endif 

       ! Get mpi parameters
#   ifdef USE_MPI
       call mpi_COMM_RANK(mpi_COMM_WORLD, sim%rank, ierr)
       call mpi_COMM_SIZE(mpi_COMM_WORLD, sim%csize, ierr)
       if (sim%rank .eq. 0) then
          allocate(sim%work1(0:sim%csize))
          allocate(sim%work2(0:sim%csize))
       else
          allocate(sim%work1(1))
          allocate(sim%work2(1))
       endif
#   else
       sim%rank = 0
       sim%csize = 1
       allocate(sim%work1(0:1))
       allocate(sim%work2(0:1))
#   endif 

    sim%root = 0
    sim%comm = mpi_COMM_WORLD

  end subroutine mpiWrap_Init

  ! =========================================================

  subroutine mpiWrap_bcastL(value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    logical, target, intent(inout) :: value
    
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_bcast(value,1, MPI_LOGICAL, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcastL

  ! =========================================================

  subroutine mpiWrap_bcastI(value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, target, intent(inout) :: value
    
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_bcast(value,1, MPI_INTEGER, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcastI

  ! =========================================================

  subroutine mpiWrap_bcastDP(value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    real(kind=8), target, intent(inout) :: value
    
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_bcast(value,1, MPI_DOUBLE_PRECISION, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcastDP

  ! =========================================================

  subroutine mpiWrap_bcastC(strlen, value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: strlen
    character(len=strlen), target, intent(inout) :: value

#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_bcast(value, strlen, MPI_CHARACTER, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcastC

  ! =========================================================

  subroutine mpiWrap_bcast1L(n, value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n
    logical, target, intent(inout) :: value(n)

#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_bcast(value,n, MPI_LOGICAL, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcast1L

  ! =========================================================

  subroutine mpiWrap_bcast1C(n, strlen, value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n, strlen
    character(len=strlen), target, intent(inout) :: value(n)

#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc, p

       ! ... Executable ... 
       p = strlen * n
       call mpi_bcast(value,p, MPI_CHARACTER, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcast1C

  ! =========================================================

  subroutine mpiWrap_bcast1I(n,value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n
    integer, target, intent(inout) :: value(n)

#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_bcast(value,n, MPI_INTEGER, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcast1I

  ! =========================================================

  subroutine mpiWrap_bcast1DP(n,value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n
    real(kind=8), target, intent(inout) :: value(n)

#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_bcast(value,n, MPI_DOUBLE_PRECISION, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcast1DP

  ! =========================================================

  ! =========================================================

  subroutine mpiWrap_bcast2I(n,m,value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n, m
    integer, target, intent(inout) :: value(n,m)

#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc, p

       ! ... Executable ... 
       p = m*n
       call mpi_bcast(value,p, MPI_INTEGER, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcast2I

  ! =========================================================

  subroutine mpiWrap_bcast2DP(n,m,value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n, m
    real(kind=8), target, intent(inout) :: value(n,m)

#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc, p

       ! ... Executable ... 
       p = m*n
       call mpi_bcast(value,p, MPI_DOUBLE_PRECISION, sim%root, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   else

#   endif

  end subroutine mpiWrap_bcast2DP

  ! =========================================================

  subroutine mpiWrap_sum3DP(n, m, l, array)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n, m, l
    real(kind=8), target, intent(inout) :: array(n, m, l)
    
    real(kind=8), pointer :: tmp(:,:,:)
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       allocate(tmp(n,m,l))
       call mpi_ALLREDUCE(array, tmp, n*m*l, mpi_DOUBLE_PRECISION, mpi_SUM, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
       array=tmp
       deallocate(tmp)
#   else
       tmp => array
#   endif

  end subroutine mpiWrap_sum3DP

  ! =========================================================

  subroutine mpiWrap_sum2DP(n, m, array)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n, m
    real(kind=8), target :: array(n, m)
    
    real(kind=8), pointer :: tmp(:,:)
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       allocate(tmp(n,m))
       call mpi_ALLREDUCE(array(1,1), tmp(1,1), n*m, mpi_DOUBLE_PRECISION, mpi_SUM, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
       array=tmp
       deallocate(tmp)
#   else
       ! Prevents compiler from complaining
       tmp =>  array
#   endif

  end subroutine mpiWrap_sum2DP

  ! =========================================================

  subroutine mpiWrap_sum1DP(n, array)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n 
    real(kind=8), target :: array(n)
    
    real(kind=8), pointer :: tmp(:)
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       allocate(tmp(n))
       call mpi_ALLREDUCE(array, tmp, n, mpi_DOUBLE_PRECISION, mpi_SUM, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
       array=tmp
       deallocate(tmp)
#   else
       ! Prevents compiler from complaining
       tmp => array
#   endif

  end subroutine mpiWrap_sum1DP

  ! =========================================================

  subroutine mpiWrap_sum1DPC(n, array)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n 
    complex(kind=8), target :: array(n)
    complex(kind=8), pointer :: tmp(:)
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       allocate(tmp(n))
       call mpi_ALLREDUCE(array, tmp, n, mpi_DOUBLE_COMPLEX, mpi_SUM, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
       array=tmp
       deallocate(tmp)
#   else
       ! Prevents compiler from complaining
       tmp => array
#   endif

  end subroutine mpiWrap_sum1DPC

  ! =========================================================

  subroutine mpiWrap_sumDPC(val)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    complex(kind=8), intent(inout) :: val
    
    complex(kind=8) :: tmp
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_ALLREDUCE(val, tmp, 1, mpi_DOUBLE_COMPLEX, mpi_SUM, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
       val=tmp
#   else
       ! Prevents compiler from complaining
       tmp = val
#   endif

  end subroutine mpiWrap_sumDPC

  ! =========================================================

  subroutine mpiWrap_sumDP(val)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    real(kind=8), intent(inout) :: val

    real(kind=8) :: tmp
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_ALLREDUCE(val, tmp, 1, mpi_DOUBLE_PRECISION, mpi_SUM, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
       val=tmp
#   else
       ! Prevents compiler from complaining
       tmp = val
#   endif

  end subroutine mpiWrap_sumDP

  ! =========================================================

  subroutine mpiWrap_collect2DP(array1, array2, n, m1, m2)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n, m1, m2
    real(kind=8), target, intent(inout) :: array1(n, m1)
    real(kind=8), target, intent(inout) :: array2(n, m2)
    
#   ifdef USE_MPI
       ! ... Local variables ...
       !integer :: ierr, m, ip, mp, istat, nc
       integer :: ierr, ip, m, mv(sim%csize), displs(sim%csize)

       ! ... Executable ... 

       m = m2*n
       call mpi_gather(m, 1, MPI_INTEGER, mv, 1, MPI_INTEGER, sim%root, sim%comm, ierr)
       displs(1) = 0
       do ip = 2, sim%csize
          displs(ip) = displs(ip-1) + mv(ip-1)
       enddo
       call mpi_gatherv(array2, m, MPI_DOUBLE_PRECISION, array1, mv, displs, &
          MPI_DOUBLE_PRECISION, sim%root, sim%comm, ierr)

       !m  = 1
       !mp = m2
       !if (sim%rank .eq. sim%root)then
       !   array1(:,1:mp) = array2(:,1:mp)
       !   do ip = 1, sim%csize-1
       !      m = m + mp
       !      array => array1(1:n,m:m1)
       !      call MPI_RECV(mp,1,MPI_INTEGER,ip,2*ip-1,sim%comm,istat,ierr)
       !      nc = mp * n 
       !      !call MPI_RECV(array,nc,MPI_DOUBLE_PRECISION,ip,2*ip,sim%comm,istat,ierr)
       !      call MPI_RECV(array1(1,m),nc,MPI_DOUBLE_PRECISION,ip,2*ip,sim%comm,istat,ierr)
       !   enddo
       !else
       !   do ip = 1, sim%csize-1
       !      if (ip .eq. sim%rank) then
       !         call MPI_SEND(mp,1,MPI_INTEGER,sim%root,2*ip-1,sim%comm,ierr)
       !         nc = mp * n 
       !         call MPI_SEND(array2(1,1),mp*n,MPI_DOUBLE_PRECISION,sim%root,2*ip,sim%comm,ierr)
       !      endif
       !   enddo
       !endif
#   else
       real(kind=8), pointer :: array(:,:)
       array => array1
       array => array2
#   endif

  end subroutine mpiWrap_collect2DP

  ! =========================================================

  subroutine mpiWrap_collectDP(array1, array2, m1, m2)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: m1, m2
    real(kind=8), target, intent(inout) :: array1(m1)
    real(kind=8), target, intent(inout) :: array2(m2)
    
#   ifdef USE_MPI
       ! ... Local variables ...
       !integer :: ierr, m, ip, mp, istat, nc
       integer :: ierr, ip, mv(sim%csize), displs(sim%csize)

       ! ... Executable ... 

       call mpi_gather(m2, 1, MPI_INTEGER, mv, 1, MPI_INTEGER, sim%root, sim%comm, ierr)
       displs(1) = 0
       do ip = 2, sim%csize
          displs(ip) = displs(ip-1) + mv(ip-1)
       enddo
       call mpi_gatherv(array2, m2, MPI_DOUBLE_PRECISION, array1, mv, displs, &
          MPI_DOUBLE_PRECISION, sim%root, sim%comm, ierr)

       !if (sim%rank .eq. sim%root)then
       !   array1(1:m2) = array2(1:m2)
       !   m = m2 + 1
       !   do ip = 1, sim%csize-1
       !      call MPI_RECV(mp,1,MPI_INTEGER,ip,2*ip-1,sim%comm,istat,ierr)
       !      call MPI_RECV(array1(m),mp,MPI_INTEGER,ip,2*ip,sim%comm,istat,ierr)
       !      m = m + mp
       !   enddo
       !else
       !   ip = sim%rank
       !   call MPI_SEND(m2,1,MPI_INTEGER,sim%root,2*ip-1,sim%comm,ierr)
       !   call MPI_SEND(array2(1),m2,MPI_INTEGER,sim%root,2*ip,sim%comm,ierr)
       !endif
       !call mpi_barrier()
#   else
       real(kind=8), pointer :: array(:)
       array => array1
       array => array2
      
#   endif

  end subroutine mpiWrap_collectDP

  ! =========================================================

  subroutine mpiWrap_maxI(value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, target :: value

    integer :: tmp
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc
 
       ! ... Executable ... 
       call mpi_ALLREDUCE(value, tmp, 1, mpi_INTEGER, mpi_MAX, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
       value=tmp
#   else
       ! Prevents compiler from complaining
       tmp = value
#   endif

  end subroutine mpiWrap_maxI

  ! =========================================================

  subroutine mpiWrap_maxDP(value)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    real(kind=8), target :: value
    
    real(kind=8) :: tmp
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc
 
       ! ... Executable ... 
       call mpi_ALLREDUCE(value, tmp, 1, mpi_DOUBLE_PRECISION, mpi_MAX, sim%comm, ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in reducing data. Stop'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
       value=tmp
#   else
       ! Prevents compiler from complaining
       tmp = value
#   endif

  end subroutine mpiWrap_maxDP

  ! =========================================================

  integer function communicator()
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    
    communicator = sim%comm

  end function communicator

  ! =========================================================

  integer function communicator_size()
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    
    communicator_size = sim%csize

  end function communicator_size

  ! =========================================================

  logical function I_am_Root()
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    
    I_am_Root = .false.
    if (sim%root .eq. sim%rank) I_am_Root = .true.

  end function I_am_Root

  ! =========================================================

  logical function my_rank_is(n)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n
    
    my_rank_is = .false.
    if (n .eq. sim%rank) my_rank_is = .true.

  end function my_rank_is

  ! =========================================================

  subroutine mpiWrap_Barrier()
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer :: ierr
    
#   ifdef USE_MPI
       ! ... Executable ... 
       call mpi_barrier(sim%comm, ierr)
#   endif

  end subroutine mpiWrap_Barrier

  ! =========================================================

  subroutine mpiWrap_send_oneI(n, dest)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in) :: n, dest

    integer :: ierr
    
#   ifdef USE_MPI
       ! ... Executable ... 
       call mpi_send(n,1, MPI_INTEGER, dest, 1, sim%comm, ierr)
#   endif

  end subroutine mpiWrap_send_oneI

  ! =========================================================

  subroutine mpiWrap_recv_oneI(n, sender)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(out) :: n 
    integer, intent(in)  :: sender

    integer :: ierr
    
#   ifdef USE_MPI
    integer :: mpistatus(MPI_STATUS_SIZE)
       ! ... Executable ... 
       call mpi_recv(n,1, MPI_INTEGER, sender, 1, sim%comm, mpistatus, ierr)
#   endif

  end subroutine mpiWrap_recv_oneI

  ! =========================================================

  subroutine mpiWrap_send_arrayDP(a,n,dest)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in)      :: n, dest 
    real(kind=8), intent(in) :: a(n)
    integer :: ierr
    
#   ifdef USE_MPI
       ! ... Executable ... 
       call mpi_send(a, n, MPI_DOUBLE_PRECISION, dest, 1, sim%comm, ierr)
#   endif

  end subroutine mpiWrap_send_arrayDP

  ! =========================================================

  subroutine mpiWrap_send_array2DP(a,n,m,dest)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in)      :: n, m, dest 
    real(kind=8), intent(in) :: a(n,m)
    integer :: ierr
    
#   ifdef USE_MPI
       ! ... Executable ... 
       call mpi_send(a, n*m, MPI_DOUBLE_PRECISION, dest, 1, sim%comm, ierr)
#   endif

  end subroutine mpiWrap_send_array2DP

  ! =========================================================

  subroutine mpiWrap_send_array2ZC(a,n,m,dest)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in)      :: n, m, dest 
    complex(kind=8), intent(in) :: a(n,m)
    integer :: ierr
    
#   ifdef USE_MPI
       ! ... Executable ... 
       call mpi_send(a, n*m, MPI_DOUBLE_COMPLEX, dest, 1, sim%comm, ierr)
#   endif

  end subroutine mpiWrap_send_array2ZC

  ! =========================================================

  subroutine mpiWrap_recv_arrayDP(a,n,sender)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in)       :: n, sender 
    real(kind=8), intent(out) :: a(n)

    integer :: ierr
    
#   ifdef USE_MPI
    integer :: mpistatus(MPI_STATUS_SIZE)
       ! ... Executable ... 
       call mpi_recv(a, n, MPI_DOUBLE_PRECISION, sender,1, sim%comm, mpistatus, ierr)
#   endif

  end subroutine mpiWrap_recv_arrayDP

  ! =========================================================

  subroutine mpiWrap_recv_array2DP(a,n,m,sender)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in)       :: n, m, sender 
    real(kind=8), intent(out) :: a(n,m)

    integer :: ierr
    
#   ifdef USE_MPI
    integer :: mpistatus(MPI_STATUS_SIZE)
       ! ... Executable ... 
       call mpi_recv(a, n*m, MPI_DOUBLE_PRECISION, sender,1, sim%comm, mpistatus, ierr)
#   endif

  end subroutine mpiWrap_recv_array2DP

  ! =========================================================

  subroutine mpiWrap_recv_array2ZC(a,n,m,sender)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    integer, intent(in)       :: n, m, sender 
    complex(kind=8), intent(out) :: a(n,m)

    integer :: ierr
    
#   ifdef USE_MPI
    integer :: mpistatus(MPI_STATUS_SIZE)
       ! ... Executable ... 
       call mpi_recv(a, n*m, MPI_DOUBLE_COMPLEX, sender,1, sim%comm, mpistatus, ierr)
#   endif

  end subroutine mpiWrap_recv_array2ZC

  ! =========================================================

  subroutine mpiWrap_Final()
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes mpi procedure.
    !
    ! Argument
    ! ========
    
#   ifdef USE_MPI
       ! ... Local variables ...
       integer :: ierr, rc

       ! ... Executable ... 
       call mpi_FINALIZE(ierr)
       if (ierr .ne. mpi_SUCCESS) then
          print *,'Error in finalize mpi program. Terminating.'
          call mpi_ABORT(sim%comm, rc, ierr)
       end if
#   endif

  end subroutine mpiWrap_Final

  ! =========================================================

end module mpiWrap
