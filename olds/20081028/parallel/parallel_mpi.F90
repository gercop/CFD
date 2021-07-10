!-----------------------------------------------------------------------------!
! module for the MPI-parallelization                                          !
!                                                                             !
! subroutines: - initmpi (initializes the MPI Prallelization)                 !
!              - time_start (starts the clock for speed-up testing)           ! 
!	       - splitup (determines the best combination of how the domain   !
!		 should be split up in order to have the minimum comunication !
!                time                                                         ! 	 	
!              - decompose_domain (decomposes the domain)                     !
!              - exchange_ghostpoints (exchanges points between the domains)  !
!              - group_for_output (groups the flow field on one processor for ! 
!                                  the output)                                !
!              - time_end (computes the execution time of the code -- for     !
!                          speed-up testing)                                  !
!                                                                             !
!-----------------------------------------------------------------------------!
module parallel_mpi
use params

#include "mpif.h" 
!include "mpif.h"

!   implicit none

   save
real,dimension(:,:),allocatable :: local_uexact
real, dimension(ny,nx) :: uexact

!--- general MPI variables
   integer,public			::my_rank,p
   integer,public			::status(MPI_STATUS_SIZE)
   integer,public			::ierror
   
!--- local processor variables
   integer,public			::local_nx,local_ny,local_bx,local_by
   integer,public			::source,tag
   integer,public			::local_x1,local_x2,local_y1,local_y2
   integer,public			::global_x1,global_x2
   integer,public			::global_y1,global_y2
   integer,dimension(4),public		::local_sd
   real,public				::local_x_0,local_y_0   
   integer,dimension(2),public		::p_yx

!--- ghost points maybe better defined in different module  
   integer				::ghost_points
   parameter(ghost_points=1)


!--- variables for timing   
   integer,public	 		:: clock_start,clock_end,clock_rate
   real*8,public		 	:: elapsed_time

!---temporal
  integer,public	 		:: clock_s,clock_e
  real*8,public			 	:: commu_time
   
contains


!-----------------------------------------------------------------------------!
!  1st: This subroutine initializes the MPI Parallelization                   !
!-----------------------------------------------------------------------------!
subroutine initmpi
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,p,ierror)
  write(*,*) 'Greetings from processor ',my_rank
  if (my_rank==0) then
    write(*,*) 'There are ',p,' processors.'
  end if
  commu_time=0.
end subroutine initmpi



!-----------------------------------------------------------------------------!
!  2nd: This subroutine takes the time at the beginning of the program        !
!       execution							      !	
!-----------------------------------------------------------------------------!
subroutine time_start
!  if(my_rank==0) then
    call SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  if(my_rank==0) then
    call SYSTEM_CLOCK(COUNT=clock_start) ! Start timing
  end if 
end subroutine time_start


!-----------------------------------------------------------------------------!
!  3rd: This subroutine takes determines the best possible combination of     !
!       the processors in order to minimize the communication time.	      !
!       Note that this subroutine has been adopted from Andi Gross	      !	
!-----------------------------------------------------------------------------!
subroutine splitup(icomm_dims)
  implicit none
  integer,dimension(p,2)	:: idivisor
  integer,dimension(2)		:: icomm_dims,ncpu
  integer 			:: np,idivisors,divisor
  integer			:: i
  integer			:: icomm,icomb,ic,icc,ib,idigit,id,icommnew

!---calculate divisors of the processors
    idivisor(i,0)=0
    if (p.ne.1) then    
      divisor=2
      np=p
      i=1
      do while (np.ne.1)
        if (mod(np,divisor)==0) then
          idivisor(i,1)=divisor
          np=np/divisor
          idivisors=i
          i=i+1
        else
          divisor=divisor+1
        end if  
      end do
    else
      idivisors=1
      idivisor(1,1)=1
    end if   
!    do i = 1,imax
!      write(*,*) i,'divisor',divisors(i)
!    end do

!--- try all possible combinations for minimum communication time
  if (idivisors.gt.11) then
     write(*,*) 'number of possible combinations too large'
     write(*,*) 'rewrite splitup.f90 routine'
     stop 'mpi'
  endif
  icomm=0
  icomb=2**idivisors
  do ic=1,icomb
     icc=ic-1
     do id=idivisors-1,1,-1
        ib=2**id
        idigit=icc/ib
        icc=mod(icc,ib)
        idivisor(id+1,2)=idigit+1
     end do
     idigit=icc
     idivisor(1,2)=idigit+1
     do id=1,2
        ncpu(id)=1
        do i=1,idivisors
           if (idivisor(i,2).eq.id) ncpu(id)=ncpu(id)*idivisor(i,1)
        end do
     end do
     icommnew=ncpu(1)*(nx)+ncpu(2)*(ny)
     if ((ic.eq.1).or.(icommnew.lt.icomm)) then
        icomm=icommnew
        do id=1,2
           icomm_dims(id)=ncpu(id)
        end do
     end if
  end do
  write(*,*) 'Split up of computational domains:'
  write(*,*) icomm_dims(2),'domains in x-direction'
  write(*,*) icomm_dims(1),'domains in y-direction'
end subroutine splitup


!-----------------------------------------------------------------------------!
!  4th: This subroutine decomposes the domain				      !
!-----------------------------------------------------------------------------!
subroutine decompose_domain
  implicit none


!---subroutine splitup determines the best combination for the minimum
!---communication time
     call splitup(p_yx) 

!---local variables for the x-direction
  if (p_yx(2)==1) then
!---no domain decomposition since only one processor is used
    local_nx = nx
    local_x_0 = x_0
    global_x1=1
    global_x2=nx
    local_bx=1
  else
!---domain decomposition according to postion of domain
    if (mod(nx,p_yx(2))==0) then
      local_nx=nx/p_yx(2)
      local_x_0=x_0+mod(my_rank,p_yx(2))*local_nx*dx
      global_x1=mod(my_rank,p_yx(2))*local_nx+1
      global_x2=(mod(my_rank,p_yx(2))+1)*local_nx
    else 
      if (mod((my_rank+1),p_yx(2)).ne.0) then
        local_nx=(nx/p_yx(2))+1
        local_x_0=x_0+mod(my_rank,p_yx(2))*local_nx*dx
        global_x1=mod(my_rank,p_yx(2))*local_nx+1
        global_x2=(mod(my_rank,p_yx(2))+1)*local_nx
      else
        local_nx=nx-((p_yx(2)-1)*((nx/p_yx(2))+1))
        local_x_0=x_0+mod(my_rank,p_yx(2))*((nx/p_yx(2))+1)*dx
        global_x1=nx-local_nx+1
        global_x2=nx
      end if
    end if
!---add ghost-points in x-direction to domain
    if (mod(my_rank,p_yx(2))==0) then
      local_bx=1
      local_nx=local_nx+ghost_points
    else if (mod((my_rank+1),p_yx(2))==0) then
      local_bx=ghost_points+1
      local_x_0=local_x_0-dx
      local_nx=local_nx+ghost_points
    else
      local_bx=ghost_points+1
      local_nx=local_nx+2*ghost_points
      local_x_0=local_x_0-dx
    end if
  end if

!---local variables for the y-direction
  if (p_yx(1)==1) then
!---no domain decomposition since only one processor is used
    local_ny = ny
    local_y_0 = y_0
    global_y1 = 1
    global_y2 =ny 
    local_by=1
  else
!---domain decomposition according to postion of domain
    if (mod(ny,p_yx(1))==0) then
      local_ny=ny/p_yx(1)
      local_y_0=y_0+(my_rank/p_yx(2))*local_ny*dy
      global_y1=(my_rank/p_yx(2))*local_ny+1
      global_y2=((my_rank/p_yx(2))+1)*local_ny
    else 
      if (my_rank.lt.(p-p_yx(2))) then
        local_ny=(ny/p_yx(1))+1
        local_y_0=y_0+(my_rank/p_yx(2))*local_ny*dy
        global_y1=(my_rank/p_yx(2))*local_ny+1
        global_y2=((my_rank/p_yx(2))+1)*local_ny
      else
        local_ny=ny-((p_yx(1)-1)*((ny/p_yx(1))+1))
        local_y_0=y_0+(my_rank/p_yx(2))*((ny/p_yx(1))+1)*dy
        global_y1=ny-local_ny+1
        global_y2=ny
      end if
    end if
!---add ghost points in x-direction to domain
    if (my_rank.lt.p_yx(2)) then
      local_by=1
      local_ny=local_ny+ghost_points
    else if (my_rank.ge.(p-p_yx(2))) then
      local_by=ghost_points+1
      local_y_0=local_y_0-dy
      local_ny=local_ny+ghost_points
    else
      local_by=ghost_points+1
      local_ny=local_ny+2*ghost_points
      local_y_0=local_y_0-dy
    end if
  end if

!---local variables for source/destination for the exchange of points
  local_sd(:)=-1.
  if (mod(my_rank,p_yx(2)).ne.0) then
    local_sd(1)=my_rank-1
  end if
  if (my_rank.ge.p_yx(2)) then
    local_sd(2)=my_rank-p_yx(2)
  end if
  if (mod((my_rank+1),p_yx(2)).ne.0) then
    local_sd(3)=my_rank+1
  end if    
  if (my_rank.lt.(p-p_yx(2))) then
    local_sd(4)=my_rank+p_yx(2)
  end if

  
!---output for checking the domain decomposition
!  write(*,*) my_rank,'local_nx',local_nx
!  write(*,*) my_rank,'local_ny',local_ny
!  write(*,*) my_rank,'local_sd(1)',local_sd(1)
!  write(*,*) my_rank,'local_sd(2)',local_sd(2)
!  write(*,*) my_rank,'local_sd(3)',local_sd(3)
!  write(*,*) my_rank,'local_sd(4)',local_sd(4)

      write (*,*) my_rank,'local_ny :',local_ny
      write (*,*) my_rank,'local_by :',local_by    
      write (*,*) my_rank,'global_y1:',global_y1 
      write (*,*) my_rank,'global_y2:',global_y2 
  
end subroutine decompose_domain




!-----------------------------------------------------------------------------!
!  5th: This subroutine exchanges the ghost points			      !
!-----------------------------------------------------------------------------!
subroutine exchange_ghostpoints(uvel_work)
  real, dimension(local_ny,local_nx) :: uvel_work
  
!---exchange points in x-direction
  if (mod((my_rank+1),p_yx(2)).ne.0) then
    tag=1
    call MPI_SEND(uvel_work(:,local_nx-2*ghost_points+1:local_nx-ghost_points),&
                  local_ny*ghost_points,MPI_DOUBLE_PRECISION,&
		  local_sd(3),tag,MPI_COMM_WORLD,ierror)
  end if
  if (mod((my_rank),p_yx(2)).ne.0) then
    tag=2
    call MPI_SEND(uvel_work(:,1+ghost_points:2*ghost_points),&
    		  local_ny*ghost_points,&
                  MPI_DOUBLE_PRECISION,local_sd(1),tag,MPI_COMM_WORLD,ierror)
  end if
  if (mod((my_rank+1),p_yx(2)).ne.0) then
    tag=2
    call MPI_RECV(uvel_work(:,local_nx-ghost_points+1:local_nx),&
    	          local_ny*ghost_points,MPI_DOUBLE_PRECISION,local_sd(3),tag,&
		  MPI_COMM_WORLD,status,ierror)  
  end if
  if (mod((my_rank),p_yx(2)).ne.0) then
    tag=1
    call MPI_RECV(uvel_work(:,1:ghost_points),local_ny*ghost_points,&
                  MPI_DOUBLE_PRECISION,local_sd(1),tag,MPI_COMM_WORLD,&
		  status,ierror)  
  end if
  
!---exchange points in y-direction
  if (my_rank.lt.(p-p_yx(2))) then
    tag=3
    call MPI_SEND(uvel_work(local_ny-2*ghost_points+1:local_ny-ghost_points,:),&
                  local_nx*ghost_points,MPI_DOUBLE_PRECISION,&
		  local_sd(4),tag,MPI_COMM_WORLD,ierror)
  end if
  if (my_rank.ge.p_yx(2)) then
    tag=4
    call MPI_SEND(uvel_work(1+ghost_points:2*ghost_points,:),&
    		  local_nx*ghost_points,&
                  MPI_DOUBLE_PRECISION,local_sd(2),tag,MPI_COMM_WORLD,ierror)
  end if
  if (my_rank.lt.(p-p_yx(2))) then
    tag=4
    call MPI_RECV(uvel_work(local_ny-ghost_points+1:local_ny,:),&
    	          local_nx*ghost_points,MPI_DOUBLE_PRECISION,local_sd(4),tag,&
		  MPI_COMM_WORLD,status,ierror)  
  end if
  if (my_rank.ge.p_yx(2)) then
    tag=3
    call MPI_RECV(uvel_work(1:ghost_points,:),local_nx*ghost_points,&
                  MPI_DOUBLE_PRECISION,local_sd(2),tag,MPI_COMM_WORLD,&
		  status,ierror)  
  end if
  
end subroutine exchange_ghostpoints




!-----------------------------------------------------------------------------!
!  6th: This subroutine groups the whole flow field on one processor          !
!       Should be replaced by MPI-IO later				      !
!-----------------------------------------------------------------------------!
subroutine group_for_output(loc_u,glob_u)
  real,dimension(local_ny,local_nx)		::loc_u
  real,dimension(ny,nx)				::glob_u
  if (my_rank.ne.0) then
    tag=1
    call MPI_SEND(global_x1,1,MPI_INTEGER,0,tag,&
                  MPI_COMM_WORLD,ierror)
    tag=2
    call MPI_SEND(global_x2,1,MPI_INTEGER,0,tag,&
                  MPI_COMM_WORLD,ierror)
    tag=3
    call MPI_SEND(global_y1,1,MPI_INTEGER,0,tag,&
                  MPI_COMM_WORLD,ierror)
    tag=4
    call MPI_SEND(global_y2,1,MPI_INTEGER,0,tag,&
                  MPI_COMM_WORLD,ierror)
    tag=5
    call MPI_SEND(loc_u(local_by:(global_y2-global_y1+local_by),&
                  local_bx:(global_x2-global_x1+local_bx)),&
                  (global_y2-global_y1+1)*(global_x2-global_x1+1),&
                  MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierror)
  end if
 
  if (my_rank==0) then
    glob_u(:,:) = 0.
    glob_u(1:local_ny,1:local_nx)=loc_u(:,:)
    do source=1,p-1
    tag=1
    call MPI_RECV(global_x1,1,MPI_INTEGER,source,tag,&
                  MPI_COMM_WORLD,status,ierror)  
    tag=2
    call MPI_RECV(global_x2,1,MPI_INTEGER,source,tag,&
                  MPI_COMM_WORLD,status,ierror)  
    tag=3
    call MPI_RECV(global_y1,1,MPI_INTEGER,source,tag,&
                  MPI_COMM_WORLD,status,ierror)  
    tag=4
    call MPI_RECV(global_y2,1,MPI_INTEGER,source,tag,&
                  MPI_COMM_WORLD,status,ierror)  
    tag=5
    call MPI_RECV(glob_u(global_y1:global_y2,global_x1:global_x2),&
    		  (global_y2-global_y1+1)*(global_x2-global_x1+1),&
		  MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,ierror)  
    end do
  end if  
end subroutine group_for_output



!-----------------------------------------------------------------------------!
!  7th:  This subroutine computes the execution time of the code              !
!-----------------------------------------------------------------------------!
subroutine time_end
  if (my_rank==0) then
    call SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing
    ! Calculate the elapsed time in seconds:
    elapsed_time= dble(clock_end-clock_start)/dble(clock_rate)
    write(*,*) 'clock_start',clock_start
    write(*,*) 'clock_end',clock_end
    write(*,*) 'clock_rate',clock_rate
    write(*,*) 'time',elapsed_time
  end if
end subroutine time_end


end module parallel_mpi
