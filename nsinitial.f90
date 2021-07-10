      MODULE sninitial

      IMPLICIT NONE       
      INCLUDE 'nspar.f90'

      save
        real*8,  dimension(:,:,:), allocatable :: rho_0
        real*8,  dimension(:,:,:), allocatable :: u_0
        real*8,  dimension(:,:,:), allocatable :: v_0
        real*8,  dimension(:,:,:), allocatable :: tp_0

        real*8,  dimension(:,:,:), allocatable :: U1
        real*8,  dimension(:,:,:), allocatable :: U2
        real*8,  dimension(:,:,:), allocatable :: U3
        real*8,  dimension(:,:,:), allocatable :: U4
        real*8,  dimension(:,:,:), allocatable :: U5

      contains
ccccc **********************************************************************
ccccc init_cond: Calculate the initial condition for all problems.
ccccc **********************************************************************
      SUBROUTINE init_cond(local_imax,x,y,z)

      USE snmpi, only: my_rank

      integer i,j,k,local_imax
      real*8 x(local_imax),y(jmax),z(kmax)

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
        call ini_cond_ac(x,y,z)
      !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5)         
        call ini_cond_sl(local_imax,x,y,z)
      !MS$ELSEIF (tp_s.EQ.6)                
        call ini_cond_bl(local_imax,x,y,z)
      !MS$ENDIF    

      if (tmin.NE.0) then
        call ini_cond_file(tmin)
      end if

      if (my_rank.EQ.0) then
        call sv_profile(y,rho_0,u_0,v_0,tp_0)
      end if

      RETURN
      END SUBROUTINE init_cond

ccccc **********************************************************************
ccccc ini_cond_file: This routine calculate the initial condition from
ccccc                the last step simulated.
ccccc **********************************************************************
      SUBROUTINE ini_cond_file(tt)
      
      USE snmpi, only: my_rank,pro,global_x1,global_x2,ghost_points

      character*05 fname
      character*20 fname_aux
      integer error,tt,i,i_aux,j,k
      real*8 local_x1,local_x2,
     &       U1_full(imax,jmax,kmax),U2_full(imax,jmax,kmax),
     &       U3_full(imax,jmax,kmax),U4_full(imax,jmax,kmax),
     &       U5_full(imax,jmax,kmax)

      call get_filename(tt,fnamegraf,'.dat',fname_aux)
      
      !MS$IF (tp_dim.EQ.1)
        call read_BIN_1D(dirgraf,fname_aux,U1_full,U2_full,
     &    U5_full,error)
      !MS$ELSEIF (tp_dim.EQ.2)
        call read_BIN_2D(dirgraf,fname_aux,U1_full,U2_full,
     &    U3_full,U5_full,error)
      !MS$ELSEIF (tp_dim.EQ.3)
        call read_BIN_3D(dirgraf,fname_aux,U1_full,U2_full,
     &    U3_full,U4_full,U5_full,error)
      !MS$ENDIF

      if (error.NE.0) then
        write(*,*) 'INITIAL FILE ERROR!!!'       
        STOP 
      end if 

      !MS$IF (parallel.EQ.0)
        local_x1=1
        local_x2=imax
      !MS$ELSEIF (parallel.EQ.1)
        if (my_rank.EQ.0) then 
          local_x1 = global_x1
          local_x2 = global_x2+ghost_points
        else if (my_rank+1.EQ.pro) then 
          local_x1 = global_x1-ghost_points
          local_x2 = global_x2
        else
          local_x1 = global_x1-ghost_points
          local_x2 = global_x2+ghost_points
        end if
      !MS$ENDIF

      i_aux=0
      do i=local_x1,local_x2
        i_aux=i_aux+1
        do j=1,jmax
          do k=1,kmax
            U1(i_aux,j,k) = U1_full(i,j,k)
            U2(i_aux,j,k) = U2_full(i,j,k)
            U3(i_aux,j,k) = U3_full(i,j,k)
            U4(i_aux,j,k) = U4_full(i,j,k)
            U5(i_aux,j,k) = U5_full(i,j,k)
          end do
        end do
      end do

      RETURN
      END SUBROUTINE ini_cond_file

ccccc **********************************************************************
ccccc ini_cond_ac: This routine calculate the initial condition for  
ccccc              acoustic problem.
ccccc **********************************************************************
      SUBROUTINE ini_cond_ac(x,y,z)

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
        USE snexactly, only: exact_acoustic,
     &    se_rho,se_u,se_v,se_w,se_p,se_tp,se_Et
      !MS$ENDIF 

      integer i,j,k
      real*8 x(imax),y(jmax),z(kmax)

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
        allocate(se_rho(imax,jmax,kmax))
        allocate(se_u  (imax,jmax,kmax))
        allocate(se_v  (imax,jmax,kmax))
        allocate(se_w  (imax,jmax,kmax))
        allocate(se_p  (imax,jmax,kmax))
        allocate(se_tp (imax,jmax,kmax))
        allocate(se_Et (imax,jmax,kmax))

        call exact_acoustic(x,y,z,0)
 
        do k=1,kmax 
          do j=1,jmax
            do i=1,imax
              U1(i,j,k) = se_rho(i,j,k)
              U2(i,j,k) = se_rho(i,j,k)*se_u(i,j,k)
              U3(i,j,k) = se_rho(i,j,k)*se_v(i,j,k)
              U4(i,j,k) = se_rho(i,j,k)*se_w(i,j,k)
              U5(i,j,k) = se_rho(i,j,k)*( 
     &          +se_tp(i,j,k)/((gamma-1.d0)*gamma*Ma**2.d0)
     &          +(se_u(i,j,k)**2.d0+se_v(i,j,k)**2.d0+se_w(i,j,k)**2.d0)/2.d0) 
            end do
	  end do
        end do

        deallocate(se_rho)
        deallocate(se_u  )
        deallocate(se_v  )
        deallocate(se_w  )
        deallocate(se_p  )
        deallocate(se_tp )
        deallocate(se_Et )
      !MS$ENDIF 
	
      RETURN
      END SUBROUTINE ini_cond_ac

ccccc **********************************************************************
ccccc ini_cond_sl: This routine calculate the initial condition for  
ccccc              mixing layer problem.
ccccc ----------------------------------------------------------------------
ccccc Velocity Profile:
ccccc   u_0* = U1*+U2*/2 + U1*-U2*/2 tanh(2y*/delta_w)
ccccc
ccccc Velocity Profile Hyphotesis:
ccccc   u_0 = u_0* /(Ma c*) 
ccccc   M1  = Ma; M2 = -Ma; Ma = U_(1,2)* / c*
ccccc ----------------------------------------------------------------------
ccccc Temperature Profile:
cccc    T_0* = 1/2cp(-u_0*^2-U1*U2*+u_0*(U1*+U2*)) 
cccc           + (T1*-T2*)(u_0*/(U1*-U2*)) 
cccc           + (T2*U1*-T1*U2*)/(U1*-U2*)
cccc
ccccc Temperature Profile Hyphotesis:
ccccc   T_0  = T_0* / T_1*
ccccc   T_1* = T_0*;  T_2* = -T_0*; T_1* = 288.16d0;
ccccc **********************************************************************
      SUBROUTINE ini_cond_sl(local_imax,x,y,z)

      integer i,j,k,local_imax
      real*8 xx,yy,zz,u_1,v_1,w_1,x(local_imax),y(jmax),z(kmax),
     &  u(local_imax,jmax,kmax),
     &  v(local_imax,jmax,kmax), 
     &  w(local_imax,jmax,kmax)

      do k=1,kmax 
        zz = z(k) 
        do j=1,jmax
          yy = y(j) 
          do i=1,local_imax
            xx = x(i)             

            u_0  (1,j,k) = (M1+M2)/(2.d0*Ma)+(M1-M2)/(2.d0*Ma)*tanh(2.d0*yy) 
            v_0  (1,j,k) = 0.d0 
            rho_0(1,j,k) = 1.d0
            !MS$IF (tp_proc.EQ.1)
              tp_0(1,j,k) = 1.d0-Ma**2.d0*(gamma-1.d0)/2.d0*
     &          (u_0(1,j,k)**2.d0+M2/M1-u_0(1,j,k)*(1.d0+M2/M1))
            !MS$ELSEIF (tp_proc.EQ.2)
              tp_0(1,j,k) = 1.d0
            !MS$ENDIF

            !MS$IF (tp_dim.EQ.2)      !Two-dimensional problem                   
	      u_1 = 
     &          -02.d0*A_1*sigma*yy*dexp(-sigma*yy**2.d0)/alpha*dsin(alpha*xx/1.d0) !Fundamental wave
     &          -04.d0*A_2*sigma*yy*dexp(-sigma*yy**2.d0)/alpha*dsin(alpha*xx/2.d0) !Subharmonic wave
     &          -08.d0*A_3*sigma*yy*dexp(-sigma*yy**2.d0)/alpha*dsin(alpha*xx/4.d0) 
     &          -16.d0*A_4*sigma*yy*dexp(-sigma*yy**2.d0)/alpha*dsin(alpha*xx/8.d0) 
              v_1 = 
     &          +A_1*dcos(alpha*xx/1.d0)*dexp(-sigma*yy**2.d0)                      !Fundamental wave 
     &          +A_2*dcos(alpha*xx/2.d0)*dexp(-sigma*yy**2.d0)			    !Subharmonic wave	
     &          +A_3*dcos(alpha*xx/4.d0)*dexp(-sigma*yy**2.d0)
     &          +A_4*dcos(alpha*xx/8.d0)*dexp(-sigma*yy**2.d0)
              w_1 = 0.d0
            !MS$ELSEIF (tp_dim.EQ.3)  !Three-dimensional problem                           
              u_1 = 
     &          +2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*(
     &          +1.d0*A_1/alpha*dsin(alpha*xx/1.d0+phi)
     &          +2.d0*A_2/alpha*dsin(alpha*xx/2.d0+phi)
     &          +1.d0*A_3/alpha*dsin(alpha*xx/1.d0+beta*zz)+1.d0*A_3/alpha*dsin(alpha*xx/1.d0-beta*zz)
     &          +2.d0*A_4/alpha*dsin(alpha*xx/2.d0+beta*zz)+2.d0*A_4/alpha*dsin(alpha*xx/2.d0-beta*zz)
     &          )
     &          -dexp(-sigma*yy**2.d0)*(
     &          +1.d0*A_3*beta*dcos(alpha*xx/1.d0+beta*zz)/alpha-1.d0*A_3*beta*dcos(alpha*xx/1.d0-beta*zz)/alpha
     &          +2.d0*A_4*beta*dcos(alpha*xx/2.d0+beta*zz)/alpha-2.d0*A_4*beta*dcos(alpha*xx/2.d0-beta*zz)/alpha
     &          )

              v_1 = (
     &          +A_1*dcos(alpha*xx/1.d0+phi)
     &          +A_2*dcos(alpha*xx/2.d0+phi)
     &          +A_3*dcos(alpha*xx/1.d0+beta*zz)+A_3*dcos(alpha*xx/1.d0-beta*zz)
     &          +A_4*dcos(alpha*xx/2.d0+beta*zz)+A_4*dcos(alpha*xx/2.d0-beta*zz)
     &          )*dexp(-sigma*yy**2.d0)

              w_1 = v_1
            !MS$ENDIF

            !MS$IF (tp_s.EQ.4)      !Velocity profile fort temporal development
              u  (i,j,k) = u_0(1,j,k) + u_1  
              v  (i,j,k) = v_0(1,j,k) + v_1  
              w  (i,j,k) = 0.d0       + w_1  
            !MS$ELSEIF (tp_s.EQ.5) !Velocity profile for spatial development
              u  (i,j,k) = u_0(1,j,k)
              v  (i,j,k) = v_0(1,j,k)
              w  (i,j,k) = 0.d0
            !MS$ENDIF         

            U1(i,j,k) = rho_0(1,j,k)
            U2(i,j,k) = rho_0(1,j,k)*u(i,j,k)
            U3(i,j,k) = rho_0(1,j,k)*v(i,j,k)          
            U4(i,j,k) = rho_0(1,j,k)*w(i,j,k)          
            !MS$IF (tp_proc.EQ.1)
              U5(i,j,k) = rho_0(1,j,k)*(tp_0(1,j,k)/((gamma-1.d0)*gamma*Ma**2.d0)
     &          +(u(i,j,k)**2.d0+v(i,j,k)**2.d0+w(i,j,k)**2.d0)/2.d0) 
            !MS$ELSEIF (tp_proc.EQ.2)
              U5(i,j,k) = 0.d0
            !MS$ENDIF         
          end do
	end do
      end do

      RETURN
      END SUBROUTINE ini_cond_sl

ccccc **********************************************************************
ccccc read_profile: Reads and calculates initial condition from Profkom.
ccccc **********************************************************************
      SUBROUTINE read_profile(u_ini,v_ini,rho_ini,tp_ini)

      integer nmax
      parameter (nmax=1000) !Size of the vector in program profcom
      integer i
      real*8 u_ini(nmax),v_ini(nmax),rho_ini(nmax),tp_ini(nmax)
       
      open(100,file='profile/profout.dat',status='unknown',form='formatted')
      do i=1,nmax
        read(100,'(f30.16,f30.16,f30.16,f30.16)') 
     &    u_ini(i),tp_ini(i),rho_ini(i),v_ini(i)
      end do
      close(100)

      RETURN
      END SUBROUTINE read_profile

ccccc **********************************************************************
ccccc profile: Reads and calculates initial condition from Profkom.
ccccc **********************************************************************
      SUBROUTINE profile(local_imax,x,y,z)

      integer local_imax
      integer i,j,k,switch,st_hw
      integer nmax,nx
      parameter (nmax=1000) !Size of the vector in program profcom
      real*8 x(local_imax),y(jmax),z(kmax),eta(nmax),peta,poly_6,
     &  u_aux  (local_imax,jmax),v_aux (local_imax,jmax),
     &  rho_aux(local_imax,jmax),tp_aux(local_imax,jmax),
     &  u_ini(nmax),v_ini(nmax),rho_ini(nmax),tp_ini(nmax)
       
      call read_profile(u_ini,v_ini,rho_ini,tp_ini)
	
      nx=nmax
      do j=1,nx
        eta(j) = dble(j-1)*0.1d0 !0.1 Tamanho da malha no programa ProfCom
      end do        
     
      !Check when u-velocity switches from 1.0 to 0.0 in 
      !the free-stream (Problem of Profkom).      
      switch = nx
      st_hw  = 3 !This paramenter is the stencil of the 
                 !polynomial minus 1 divided by 2.
      do j=2,nx
        if ((u_ini(j)-u_ini(j-1)).LT.-0.9d0) then
          switch = j-st_hw 
          exit
        end if
      end do

      !Interpolation of initial data onto the physical grid 
      !of computation using polynomial interpolation.
      do i=1,local_imax
        do j=1,jmax
          !Get eta coordinate for the xy-grid point and check 
          !at what location k, eta(k) >= peta.
          peta = y(j)*dsqrt(Re/x(i))
          do k=1,nx
            if (eta(k).GE.peta) then
              exit
            end if
          end do

          if (k.LT.switch) then
            if (k.LE.1+st_hw) then
  	      u_aux(i,j) = 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,1)*u_ini(1) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,2)*u_ini(2) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,3)*u_ini(3) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,4)*u_ini(4) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,5)*u_ini(5) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,6)*u_ini(6) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,7)*u_ini(7)
              v_aux(i,j) = 
     &          (poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,1)*v_ini(1) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,2)*v_ini(2) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,3)*v_ini(3) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,4)*v_ini(4) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,5)*v_ini(5) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,6)*v_ini(6) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,7)*v_ini(7)
     &          )/dsqrt(Re*x(i))
              tp_aux(i,j) = 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,1)*tp_ini(1) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,2)*tp_ini(2) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,3)*tp_ini(3) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,4)*tp_ini(4) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,5)*tp_ini(5) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,6)*tp_ini(6) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,7)*tp_ini(7)
              rho_aux(i,j) = 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,1)*rho_ini(1) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,2)*rho_ini(2) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,3)*rho_ini(3) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,4)*rho_ini(4) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,5)*rho_ini(5) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,6)*rho_ini(6) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,7)*rho_ini(7)
            else if (k.GE.nx-st_hw) then
              u_aux(i,j) = 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,1)*u_ini(nx-6) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,2)*u_ini(nx-5) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,3)*u_ini(nx-4) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,4)*u_ini(nx-3) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,5)*u_ini(nx-2) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,6)*u_ini(nx-1) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,7)*u_ini(nx  )
              v_aux(i,j) =
     &          (poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,1)*v_ini(nx-6) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,2)*v_ini(nx-5) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,3)*v_ini(nx-4) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,4)*v_ini(nx-3) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,5)*v_ini(nx-2) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,6)*v_ini(nx-1) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,7)*v_ini(nx  )
     &          )/dsqrt(Re*x(i))
              tp_aux(i,j) = 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,1)*tp_ini(nx-6) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,2)*tp_ini(nx-5) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,3)*tp_ini(nx-4) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,4)*tp_ini(nx-3) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,5)*tp_ini(nx-2) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,6)*tp_ini(nx-1) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,7)*tp_ini(nx  )
              rho_aux(i,j) = 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,1)*rho_ini(nx-6) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,2)*rho_ini(nx-5) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,3)*rho_ini(nx-4) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,4)*rho_ini(nx-3) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,5)*rho_ini(nx-2) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,6)*rho_ini(nx-1) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),eta(nx-2),eta(nx-1),eta(nx),peta,7)*rho_ini(nx  )
            else
              u_aux(i,j) = 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,1)*u_ini(k-3) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,2)*u_ini(k-2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,3)*u_ini(k-1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,4)*u_ini(k  ) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,5)*u_ini(k+1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,6)*u_ini(k+2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,7)*u_ini(k+3)
              v_aux(i,j) =
     &          (poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,1)*v_ini(k-3) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,2)*v_ini(k-2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,3)*v_ini(k-1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,4)*v_ini(k  ) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,5)*v_ini(k+1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,6)*v_ini(k+2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,7)*v_ini(k+3)
     &          )/dsqrt(Re*x(i))
              tp_aux(i,j) = 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,1)*tp_ini(k-3) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,2)*tp_ini(k-2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,3)*tp_ini(k-1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,4)*tp_ini(k  ) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,5)*tp_ini(k+1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,6)*tp_ini(k+2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,7)*tp_ini(k+3)
              rho_aux(i,j) = 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,1)*rho_ini(k-3) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,2)*rho_ini(k-2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,3)*rho_ini(k-1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,4)*rho_ini(k  ) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,5)*rho_ini(k+1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,6)*rho_ini(k+2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,7)*rho_ini(k+3)
            end if
          else
            u_aux  (i,j) = u_aux  (i,j-1)
            v_aux  (i,j) = v_aux  (i,j-1)
            tp_aux (i,j) = tp_aux (i,j-1)
            rho_aux(i,j) = rho_aux(i,j-1)
          end if        
        end do
      end do

      do k=1,kmax         
        do j=1,jmax
          do i=1,local_imax
            rho_0(i,j,k) = rho_aux(i,j)
            u_0  (i,j,k) = u_aux  (i,j)
            v_0  (i,j,k) = v_aux  (i,j)
            tp_0 (i,j,k) = tp_aux (i,j)
          end do
        end do
      end do

      RETURN
      END SUBROUTINE profile

ccccc ************************************************************************
ccccc poly_6: This routine to get the coefficients for a 6th order polynomial.
ccccc ************************************************************************
      FUNCTION poly_6(x1,x2,x3,x4,x5,x6,x7,x,which)

      integer which
      real*8 x1,x2,x3,x4,x5,x6,x7,x,poly_6

      select case (which)
      case(1)
        poly_6=((x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/ 
     &         ((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)*(x1-x7))
      case(2)
        poly_6=((x-x1)*(x-x3)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/ 
     &         ((-x1+x2)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)*(x2-x7))
      case(3)
        poly_6=((x-x1)*(x-x2)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/ 
     &         ((-x1+x3)*(-x2+x3)*(x3-x4)*(x3-x5)*(x3-x6)*(x3-x7))
      case(4)
        poly_6=((x-x1)*(x-x2)*(x-x3)*(x-x5)*(x-x6)*(x-x7))/ 
     &         ((-x1+x4)*(-x2+x4)*(-x3+x4)*(x4-x5)*(x4-x6)*(x4-x7))
      case(5)
        poly_6=((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x6)*(x-x7))/ 
     &         ((-x1+x5)*(-x2+x5)*(-x3+x5)*(-x4+x5)*(x5-x6)*(x5-x7))
      case(6)
        poly_6=((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x7))/ 
     &         ((-x1+x6)*(-x2+x6)*(-x3+x6)*(-x4+x6)*(-x5+x6)*(x6-x7))
      case(7)
        poly_6=((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x6))/ 
     &         ((-x1+x7)*(-x2+x7)*(-x3+x7)*(-x4+x7)*(-x5+x7)*(-x6+x7))
      end select

      RETURN
      END FUNCTION poly_6

ccccc **********************************************************************
ccccc ini_cond_bl: This routine calculate the initial condition for  
ccccc              boundary layer problem.
ccccc **********************************************************************
      SUBROUTINE ini_cond_bl(local_imax,x,y,z)

      integer i,j,k,local_imax
      real*8 x(local_imax),y(jmax),z(kmax)

      call profile(local_imax,x,y,z)

      do k=1,kmax 
        do j=1,jmax
          do i=1,local_imax
            U1(i,j,k) = rho_0(i,j,k)
            U2(i,j,k) = rho_0(i,j,k)*u_0(i,j,k)
            U3(i,j,k) = rho_0(i,j,k)*v_0(i,j,k)
            U4(i,j,k) = 0.d0            
            U5(i,j,k) = rho_0(i,j,k)*( 
     &        +tp_0(i,j,k)/((gamma-1.d0)*gamma*Ma**2.d0)
     &        +(u_0(i,j,k)**2.d0+v_0(i,j,k)**2.d0)/2.d0) 
          end do
	end do
      end do

      RETURN
      END SUBROUTINE ini_cond_bl

      END MODULE sninitial

