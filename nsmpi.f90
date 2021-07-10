ccccc ------------------------------------------------------------------------------------------------
ccccc module for the MPI-parallelization                                          
ccccc                                                                        
ccccc subroutines: - initmpi (initializes the MPI Prallelization)                 
ccccc  - time_start (starts the clock for speed-up testing)            
ccccc  - splitup (determines the best combination of how the domain   
ccccc should be split up in order to have the minimum comunication time  

ccccc  - decompose_domain (decomposes the domain)                     
ccccc  - exchange_ghostpoints (exchanges points between the domains)  
ccccc  - group_for_output (groups the flow field on one processor for the output)
ccccc  - time_end (computes the execution time of the code - for speed-up testing)
ccccc ------------------------------------------------------------------------------------------------      
      MODULE snmpi

      IMPLICIT NONE
      INCLUDE 'nspar.f90'
      !MS$IF (parallel.EQ.1)
        INCLUDE 'mpif.h'
      !MS$ENDIF

      save
      !General MPI variables
        integer, parameter              :: ghost_points=3
        integer, public		 	:: my_rank,pro
        integer, public			:: ierror  
        !MS$IF (parallel.EQ.1)
          integer, private              :: status(MPI_STATUS_SIZE)
        !MS$ENDIF
      !Local processor variables
        integer, public			:: local_imax,local_bx
        integer, public			:: global_x1,global_x2
        integer, dimension(2),public	:: local_sd
        integer, private	 	:: clock_start,clock_end,clock_rate  
        real*8,  private                :: elapsed_time
        real*8,  private                :: cpu_start,cpu_end

      contains
ccccc ************************************************************************************************
ccccc initmpi: Initialize MPI.
ccccc ************************************************************************************************
      SUBROUTINE initmpi

      !MS$IF (parallel.EQ.1)        
        call MPI_INIT(ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,pro,ierror)        
        !MS$IF (show_msg.EQ.1)
          if (my_rank.EQ.0) then
            write(*,*) 'THERE ARE ',pro,' PROCESSORS.'
            write(*,*)
          end if
          write(*,*) 'GREETINGS FROM PROCESSOR: ',my_rank
        !MS$ENDIF
      !MS$ELSE
        !MS$IF (show_msg.EQ.1)
          write(*,*) 'GREETINGS FROM PROCESSOR: 1'
        !MS$ENDIF
        my_rank = 0
        pro     = 1
      !MS$ENDIF

      RETURN
      END SUBROUTINE initmpi

ccccc ************************************************************************************************
ccccc time_start: This subroutine takes the time at the beginning of the program execution.	
ccccc ************************************************************************************************
      SUBROUTINE time_start

      call CPU_TIME(cpu_start)
      call SYSTEM_CLOCK(COUNT_RATE=clock_rate) !Find the rate
      if (my_rank.EQ.0) then
        call SYSTEM_CLOCK(COUNT=clock_start)   !Start timing
      end if 

      RETURN 
      END SUBROUTINE time_start

ccccc **********************************************************************
ccccc decompose_domain: This subroutine decomposes the domain.				      
ccccc **********************************************************************
      SUBROUTINE decompose(nx)

      integer nx

      !Local variables for the x-direction
      if (pro.EQ.1) then
        !No domain decomposition since only one processor is used
        local_imax = nx
        local_bx   = 1
        global_x1  = 1
        global_x2  = nx
      else
        !Domain decomposition according to postion of domain
        if (mod(nx,pro).EQ.0) then
          local_imax = nx/pro
          global_x1  = mod(my_rank,pro)*local_imax+1
          global_x2  = (mod(my_rank,pro)+1)*local_imax
        else 
          if (mod((my_rank+1),pro).NE.0) then
            local_imax = (nx/pro)+1
            global_x1  = mod(my_rank,pro)*local_imax+1
            global_x2  = (mod(my_rank,pro)+1)*local_imax
          else
            local_imax = nx-((pro-1)*((nx/pro)+1))
            global_x1  = nx-local_imax+1
            global_x2  = nx
          end if
        end if

        !Add ghost-points in x-direction to domain
        if (mod(my_rank,pro).EQ.0) then
          local_imax = local_imax+ghost_points
          local_bx   = 1
        else if (mod(my_rank+1,pro).EQ.0) then
          local_imax = local_imax+ghost_points
          local_bx   = ghost_points+1
        else
          local_imax = local_imax+2*ghost_points
          local_bx   = ghost_points+1
        end if
      end if

      !Local variables for source/destination for the exchange of points
      local_sd(1)=-1
      local_sd(2)=-1
      if (mod(my_rank,pro).NE.0) then
        local_sd(1) = my_rank-1
      end if
      if (mod((my_rank+1),pro).NE.0) then
        local_sd(2) = my_rank+1
      end if    
 
      !MS$IF (parallel.EQ.1).AND.(show_msg.EQ.1)
        write (*,*) my_rank,'LOCAL_IMAX :',local_imax
        write (*,*) my_rank,'LOCAL_BX   :',local_bx 
        write (*,*) my_rank,'GLOBAL_X1  :',global_x1 
        write (*,*) my_rank,'GLOBAL_X2  :',global_x2  
      !MS$ENDIF

      RETURN 
      END SUBROUTINE decompose

ccccc **********************************************************************
ccccc exchange_ghostpoints: This subroutine exchanges the ghost points			      	      
ccccc **********************************************************************
      SUBROUTINE exchange_ghostpoints(jmax,kmax,tag,E)

      integer k,jmax,kmax,tag,IB,IN,SZ
      real*8 E(local_imax,jmax,kmax)

      !MS$IF (parallel.EQ.1)
        do k=1,kmax
          if (mod((my_rank+1),pro).NE.0) then
            IB = local_imax-2*ghost_points+1
            IN = local_imax-ghost_points
            SZ = jmax*ghost_points  
            call MPI_SEND(E(IB:IN,1:jmax,k),SZ,
     &        MPI_DOUBLE_PRECISION,local_sd(2),
     &        (tag*100000)+k,MPI_COMM_WORLD,ierror)
          end if

          if (mod((my_rank),pro).NE.0) then
            IB = 1+ghost_points
            IN = 2*ghost_points
            SZ = jmax*ghost_points  
            call MPI_SEND(E(IB:IN,1:jmax,k),SZ,
     &        MPI_DOUBLE_PRECISION,local_sd(1),
     &        ((tag+1)*100000)+k,MPI_COMM_WORLD,ierror)
          end if

          if (mod((my_rank),pro).NE.0) then
            IB = 1
            IN = ghost_points
            SZ = jmax*ghost_points
            call MPI_RECV(E(IB:IN,1:jmax,k),SZ,
     &        MPI_DOUBLE_PRECISION,local_sd(1),
     &        (tag*100000)+k,MPI_COMM_WORLD,status,ierror)
          end if

          if (mod((my_rank+1),pro).NE.0) then
            IB = local_imax-ghost_points+1
            IN = local_imax
            SZ = jmax*ghost_points  
            call MPI_RECV(E(IB:IN,1:jmax,k),SZ,
     &        MPI_DOUBLE_PRECISION,local_sd(2),
     &        ((tag+1)*100000)+k,MPI_COMM_WORLD,status,ierror)  
          end if
        end do  
      !MS$ENDIF     
 
      RETURN 
      END SUBROUTINE exchange_ghostpoints

ccccc ************************************************************************************************
ccccc group_for_output: This subroutine groups the whole flow field on one processor.
ccccc ************************************************************************************************
      SUBROUTINE group_for_output(imax,jmax,kmax,
     &  local_U1,local_U2,local_U3,local_U4,local_U5,
     &  global_U1,global_U2,global_U3,global_U4,global_U5)

      integer i,j,k,imax,jmax,kmax,source,IB,IN,SZ,tag
      real*8 
     &  local_U1(local_imax,jmax,kmax),local_U2(local_imax,jmax,kmax),
     &  local_U3(local_imax,jmax,kmax),local_U4(local_imax,jmax,kmax),
     &  local_U5(local_imax,jmax,kmax),
     &  global_U1(imax,jmax,kmax),global_U2(imax,jmax,kmax),
     &  global_U3(imax,jmax,kmax),global_U4(imax,jmax,kmax),
     &  global_U5(imax,jmax,kmax)   

      !MS$IF (parallel.EQ.1)
        if (my_rank.NE.0) then

          IB = local_bx
          IN = (global_x2-global_x1+local_bx)
          SZ = (global_x2-global_x1+1)*jmax*kmax       

          tag = 1
          call MPI_SEND(global_x1,1,MPI_INTEGER,0,tag,
     &      MPI_COMM_WORLD,ierror)
          tag = 2
          call MPI_SEND(global_x2,1,MPI_INTEGER,0,tag,
     &      MPI_COMM_WORLD,ierror)
          tag = 3
          call MPI_SEND(local_U1(IB:IN,1:jmax,1:kmax),SZ,
     &      MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierror)
          tag = 4
          call MPI_SEND(local_U2(IB:IN,1:jmax,1:kmax),SZ,
     &      MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierror)
          tag = 5
          call MPI_SEND(local_U3(IB:IN,1:jmax,1:kmax),SZ,
     &      MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierror)
          tag = 6
          call MPI_SEND(local_U4(IB:IN,1:jmax,1:kmax),SZ,
     &      MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierror)
          tag = 7
          call MPI_SEND(local_U5(IB:IN,1:jmax,1:kmax),SZ,
     &      MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierror)
        end if
      !MS$ENDIF 

      if (my_rank.EQ.0) then

        do k=1,kmax
          do j=1,jmax
            do i=1,local_imax
              global_U1(i,j,k) = local_U1(i,j,k)
              global_U2(i,j,k) = local_U2(i,j,k)
              global_U3(i,j,k) = local_U3(i,j,k)
              global_U4(i,j,k) = local_U4(i,j,k)
              global_U5(i,j,k) = local_U5(i,j,k)
            end do
          end do
        end do 

        !MS$IF (parallel.EQ.1)
          do source=1,pro-1
            tag = 1
            call MPI_RECV(global_x1,1,MPI_INTEGER,source,tag,
     &        MPI_COMM_WORLD,status,ierror)  
            tag = 2
            call MPI_RECV(global_x2,1,MPI_INTEGER,source,tag,
     &        MPI_COMM_WORLD,status,ierror)  

            IB = global_x1
            IN = global_x2
            SZ = (global_x2-global_x1+1)*jmax*kmax

            tag = 3
            call MPI_RECV(global_U1(IB:IN,1:jmax,1:kmax),SZ,
     &        MPI_DOUBLE_PRECISION,source,tag,
     &        MPI_COMM_WORLD,status,ierror)  
            tag = 4
            call MPI_RECV(global_U2(IB:IN,1:jmax,1:kmax),SZ,
     &        MPI_DOUBLE_PRECISION,source,tag,
     &        MPI_COMM_WORLD,status,ierror)  
            tag = 5
            call MPI_RECV(global_U3(IB:IN,1:jmax,1:kmax),SZ,
     &        MPI_DOUBLE_PRECISION,source,tag,
     &        MPI_COMM_WORLD,status,ierror)  
            tag = 6
            call MPI_RECV(global_U4(IB:IN,1:jmax,1:kmax),SZ,
     &        MPI_DOUBLE_PRECISION,source,tag,
     &        MPI_COMM_WORLD,status,ierror)  
            tag = 7
            call MPI_RECV(global_U5(IB:IN,1:jmax,1:kmax),SZ,
     &        MPI_DOUBLE_PRECISION,source,tag,
     &        MPI_COMM_WORLD,status,ierror)  
          end do
        !MS$ENDIF
      end if  

      RETURN 
      END SUBROUTINE group_for_output

ccccc **********************************************************************
ccccc Protostr: Convert from process number to string values.
ccccc **********************************************************************
      SUBROUTINE protostr(np,fname,ext,fname_aux)

      IMPLICIT NONE

      character*10 fname
      character*16 fname_aux
      character*04 ext
      integer np

      if (np.EQ.0) then
        fname_aux=fname//'01'//ext
      else if (np.EQ.1) then
        fname_aux=fname//'02'//ext
      else if (np.EQ.2) then
        fname_aux=fname//'03'//ext
      else if (np.EQ.3) then
        fname_aux=fname//'04'//ext
      else if (np.EQ.4) then
        fname_aux=fname//'05'//ext
      else if (np.EQ.5) then
        fname_aux=fname//'06'//ext
      else if (np.EQ.6) then
        fname_aux=fname//'07'//ext
      else if (np.EQ.7) then
        fname_aux=fname//'08'//ext
      else if (np.EQ.8) then
        fname_aux=fname//'09'//ext
      else if (np.EQ.9) then
        fname_aux=fname//'10'//ext
      else if (np.EQ.10) then
        fname_aux=fname//'11'//ext
      else if (np.EQ.11) then
        fname_aux=fname//'12'//ext
      else if (np.EQ.12) then
        fname_aux=fname//'13'//ext
      end if

      RETURN
      END SUBROUTINE protostr

ccccc ************************************************************************************************
ccccc time_end: This subroutine computes the execution time of the code.              
ccccc ************************************************************************************************
      SUBROUTINE time_end

      if (my_rank.EQ.0) then
        call CPU_TIME(cpu_end)
        call SYSTEM_CLOCK(COUNT=clock_end) 
        elapsed_time=dble(clock_end-clock_start)/dble(clock_rate)

        write(*,*) '***************************************************'
        write(*,*) 'Clock_start.....:',clock_start
        write(*,*) 'Clock_end.......:',clock_end
        write(*,*) 'Clock_rate......:',clock_rate
        write(*,*) 'Time............:',elapsed_time
        write(*,*) 'CPU/Time........:',cpu_end-cpu_start
        write(*,*) '***************************************************'
      end if

      !MS$IF (parallel.EQ.1)
        call MPI_FINALIZE(ierror)
      !MS$ENDIF

      RETURN 
      END SUBROUTINE time_end

      END MODULE snmpi

