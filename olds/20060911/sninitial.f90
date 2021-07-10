      MODULE sninitial

      save
        real*8,  dimension(:,:,:), allocatable :: rho_0
        real*8,  dimension(:,:,:), allocatable :: u_0
        real*8,  dimension(:,:,:), allocatable :: v_0
        real*8,  dimension(:,:,:), allocatable :: tp_0 

      contains
ccccc **********************************************************************
ccccc init_cond: Calculate the initial condition for all problems.
ccccc **********************************************************************
      SUBROUTINE init_cond(local_imax,x,y,z,E1,E2,E3,E4,E5)

      IMPLICIT NONE       
      INCLUDE 'snparv35.f90'

      integer i,j,k,local_imax
      real*8 x(local_imax),y(jmax),z(kmax),
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)         
        call ini_cond_slayer(local_imax,x,y,z,E1,E2,E3,E4,E5)
      !MS$ELSEIF (tp_s.EQ.6)                
        call ini_cond_blayer(local_imax,x,y,z,E1,E2,E3,E4,E5)
      !MS$ENDIF

      RETURN
      END SUBROUTINE init_cond

ccccc **********************************************************************
ccccc ini_cond_slayer: This routine calculate the initial condition for  
ccccc                  shear layer problem.
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
      SUBROUTINE ini_cond_slayer(local_imax,x,y,z,E1,E2,E3,E4,E5)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer i,j,k,local_imax
      real*8 x(local_imax),y(jmax),z(kmax),xx,yy,zz,
     &  Ma_c,u_1,v_1,w_1,A_1,A_2,Phi,alpha_1,alpha_2,
     &  E1 (local_imax,jmax,kmax),E2 (local_imax,jmax,kmax),
     &  E3 (local_imax,jmax,kmax),E4 (local_imax,jmax,kmax),
     &  E5 (local_imax,jmax,kmax),rho(local_imax,jmax,kmax),
     &  u  (local_imax,jmax,kmax),v  (local_imax,jmax,kmax),
     &  w  (local_imax,jmax,kmax),Et (local_imax,jmax,kmax),
     &  p  (local_imax,jmax,kmax),tp (local_imax,jmax,kmax)

      do k=1,kmax 
        zz = z(k) 
        do j=1,jmax
          yy = y(j) 
          do i=1,local_imax
            xx = x(i) 
            
            Ma_c = (M1+M2)/(2.d0*Ma)        

            u_0  (i,j,k) = Ma_c+(M1-M2)/(2.d0*Ma)*tanh(2.d0*yy)         
            v_0  (i,j,k) = 0.d0 
            rho_0(i,j,k) = 1.d0
            !MS$IF (tp_proc.EQ.1)
              tp_0(i,j,k) = 1.d0-Ma**2.d0*(gamma-1.d0)/2.d0*
     &          (u_0(i,j,k)**2.d0+M2/M1-u_0(i,j,k)*(1.d0+M2/M1))
            !MS$ELSEIF (tp_proc.EQ.2)
              tp_0(i,j,k) = 1.d0
            !MS$ENDIF

            if ((jmax.NE.1).AND.(kmax.EQ.1)) then !Two-Dimension Problem
              A_1     = A
              A_2     = A/2.d0
              alpha_1 = alpha
              alpha_2 = alpha/2.d0
	      u_1     = -2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*(
     &          A_1*dsin(alpha_1*xx)*alpha_2+
     &          A_2*dsin(alpha_2*xx)*alpha_1)/(alpha_1*alpha_2)
              v_1 = (A_1*dcos(alpha_1*xx)+A_2*dcos(alpha_2*xx))*
     &          dexp(-sigma*yy**2.d0) 
              w_1 = 0.d0
            else if ((jmax.NE.1).AND.(kmax.NE.1)) then !Three-Dimension Problem              
              A_1 = A
              A_2 = A/2.d0
              phi = 1.d0/2.d0*pi                             
              u_1 = -dexp(-sigma*yy**2.d0)*(
     &                +A_2*beta*dcos(alpha*xx+beta*zz)/alpha
     &                -A_2*beta*dcos(alpha*xx-beta*zz)/alpha)
     &              +2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*(
     &                +A_1*dsin(alpha*xx+phi)/alpha+ 
     &                +A_2*dsin(alpha*xx+beta*zz)/alpha
     &                +A_2*dsin(alpha*xx-beta*zz)/alpha)
              v_1 = (A_1*dcos(alpha*xx+phi)+
     &               A_2*dcos(alpha*xx+beta*zz)+
     &               A_2*dcos(alpha*xx-beta*zz))*
     &              dexp(-sigma*yy**2.d0)
              w_1 = (A_1*dcos(alpha*xx+phi)+
     &               A_2*dcos(alpha*xx+beta*zz)+
     &               A_2*dcos(alpha*xx-beta*zz))*
     &              dexp(-sigma*yy**2.d0)
            end if

            !MS$IF (tp_s.EQ.3)                     !Velocity profile for Temporal development
              u  (i,j,k) = u_0(i,j,k) + 1.d0/(Ma*c)*u_1  
              v  (i,j,k) = v_0(i,j,k) + 1.d0/(Ma*c)*v_1  
              w  (i,j,k) = 0.d0       + 1.d0/(Ma*c)*w_1  
            !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5) !Velocity profile for Spatial development
              u  (i,j,k) = u_0(i,j,k)
              v  (i,j,k) = v_0(i,j,k)
              w  (i,j,k) = 0.d0
            !MS$ENDIF         
            rho(i,j,k) = rho_0(i,j,k) !rho = rho_0 = rho_0/rho_oo=1.d0
            tp (i,j,k) = tp_0 (i,j,k) !tp  = tp_0  = tp_0/tp_oo        
            Et (i,j,k) = rho(i,j,k)*( 
     &        +tp(i,j,k)/((gamma-1.d0)*gamma*Ma**2.d0)
     &        +(u(i,j,k)**2.d0+v(i,j,k)**2.d0+w(i,j,k)**2.d0)/2.d0) 

            E1(i,j,k) = rho(i,j,k)
            E2(i,j,k) = rho(i,j,k)*u(i,j,k)
            E3(i,j,k) = rho(i,j,k)*v(i,j,k)          
            E4(i,j,k) = rho(i,j,k)*w(i,j,k)          
            E5(i,j,k) = Et(i,j,k)          
          end do
	end do
      end do

      RETURN
      END SUBROUTINE ini_cond_slayer

ccccc **********************************************************************
ccccc read_profile: Reads and calculates initial condition from Profkom.
ccccc **********************************************************************
      SUBROUTINE read_profile(local_imax,x,y,z)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*16 fbase
      character*18 fname
      character*72 inf(100)
      integer nmax,local_imax
      parameter (nmax=1000) !Size of the vector in program profcom
      integer i,j,k,ioerr,switch,ninf,imach,
     &  nt,nz,ny,nx,nparam,ntimes(512),st_hw
      real*8 x(local_imax),y(jmax),z(kmax),eta(nmax),peta,poly_6,
     &  u_aux  (local_imax,jmax),v_aux (local_imax,jmax),
     &  rho_aux(local_imax,jmax),tp_aux(local_imax,jmax),
     &  u_ini(nmax),v_ini(nmax),rho_ini(nmax),tp_ini(nmax)
       
      fbase = 'profile/profout'
      fname = 'profile/profout.s8'
      imach = 2
      !MS$IF (tp_s.EQ.6)
        call mkfname(fname,fbase,imach,ioerr)
        call readcd(50,50,fbase,imach,nt,nz,ny,nx,
     &    nparam,ntimes,inf,ninf)
        call readd(50,u_ini,  1,1,ioerr)
        call readd(50,tp_ini, 1,2,ioerr)
        call readd(50,rho_ini,1,3,ioerr)
        call readd(50,v_ini,  1,4,ioerr) 
        close(50)
      !MS$ENDIF 

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
          do k=1,nmax
            if (eta(k).GE.peta) then
              exit
            end if
          end do

          if (k.LT.switch) then
            if (k.LE.1+st_hw) then
  	      u_aux(i,j) = 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,1)*u_ini(1) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,2)*u_ini(2) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,3)*u_ini(3) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,4)*u_ini(4) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,5)*u_ini(5) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,6)*u_ini(6) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,7)*u_ini(7)
              v_aux(i,j) = 
     &          (poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,1)*v_ini(1) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,2)*v_ini(2) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,3)*v_ini(3) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,4)*v_ini(4) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,5)*v_ini(5) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,6)*v_ini(6) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,7)*v_ini(7)
     &          )/dsqrt(Re*x(i))
              tp_aux(i,j) = 
     &          poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,1)*tp_ini(1) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,2)*tp_ini(2) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,3)*tp_ini(3) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,4)*tp_ini(4) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,5)*tp_ini(5) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,6)*tp_ini(6) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,7)*tp_ini(7)
              rho_aux(i,j) = 
     &          poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,1)*rho_ini(1) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,2)*rho_ini(2) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,3)*rho_ini(3) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,4)*rho_ini(4) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,5)*rho_ini(5) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,6)*rho_ini(6) 
     &          +poly_6(eta(1),eta(2),eta(3),eta(4),eta(5),
     &          eta(6),eta(7),peta,7)*rho_ini(7)
            else if (k.GE.nx-st_hw) then
              u_aux(i,j) = 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,1)*u_ini(nx-6) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,2)*u_ini(nx-5) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,3)*u_ini(nx-4) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,4)*u_ini(nx-3) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,5)*u_ini(nx-2) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,6)*u_ini(nx-1) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,7)*u_ini(nx  )
              v_aux(i,j) =
     &          (poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,1)*v_ini(nx-6) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,2)*v_ini(nx-5) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,3)*v_ini(nx-4) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,4)*v_ini(nx-3) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,5)*v_ini(nx-2) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,6)*v_ini(nx-1) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,7)*v_ini(nx  )
     &          )/dsqrt(Re*x(i))
              tp_aux(i,j) = 
     &          poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,1)*tp_ini(nx-6) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,2)*tp_ini(nx-5) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,3)*tp_ini(nx-4) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,4)*tp_ini(nx-3) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,5)*tp_ini(nx-2) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,6)*tp_ini(nx-1) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,7)*tp_ini(nx  )
              rho_aux(i,j) = 
     &          poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,1)*rho_ini(nx-6) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,2)*rho_ini(nx-5) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,3)*rho_ini(nx-4) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,4)*rho_ini(nx-3) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,5)*rho_ini(nx-2) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,6)*rho_ini(nx-1) 
     &          +poly_6(eta(nx-6),eta(nx-5),eta(nx-4),eta(nx-3),
     &          eta(nx-2),eta(nx-1),eta(nx),peta,7)*rho_ini(nx  )
            else
              u_aux(i,j) = 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,1)*u_ini(k-3) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,2)*u_ini(k-2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,3)*u_ini(k-1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,4)*u_ini(k  ) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,5)*u_ini(k+1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,6)*u_ini(k+2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,7)*u_ini(k+3)
              v_aux(i,j) =
     &          (poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,1)*v_ini(k-3) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,2)*v_ini(k-2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,3)*v_ini(k-1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,4)*v_ini(k  ) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,5)*v_ini(k+1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,6)*v_ini(k+2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,7)*v_ini(k+3)
     &          )/dsqrt(Re*x(i))
              tp_aux(i,j) = 
     &          poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,1)*tp_ini(k-3) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,2)*tp_ini(k-2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,3)*tp_ini(k-1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,4)*tp_ini(k  ) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,5)*tp_ini(k+1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,6)*tp_ini(k+2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,7)*tp_ini(k+3)
              rho_aux(i,j) = 
     &          poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,1)*rho_ini(k-3) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,2)*rho_ini(k-2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,3)*rho_ini(k-1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,4)*rho_ini(k  ) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,5)*rho_ini(k+1) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,6)*rho_ini(k+2) 
     &          +poly_6(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),
     &          eta(k+2),eta(k+3),peta,7)*rho_ini(k+3)
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
      END SUBROUTINE read_profile

ccccc ************************************************************************
ccccc poly_6: This routine to get the coefficients for a 6th order polynomial.
ccccc ************************************************************************
      FUNCTION poly_6(x1,x2,x3,x4,x5,x6,x7,x,which)

      IMPLICIT NONE

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
ccccc ini_cond_blayer: This routine calculate the initial condition for  
ccccc                  boundary layer problem.
ccccc **********************************************************************
      SUBROUTINE ini_cond_blayer(local_imax,x,y,z,E1,E2,E3,E4,E5)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer i,j,k,local_imax
      real*8 x(local_imax),y(jmax),z(kmax),
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax)

      call read_profile(local_imax,x,y,z)

      do k=1,kmax 
        do j=1,jmax
          do i=1,local_imax
            E1(i,j,k) = rho_0(i,j,k)
            E2(i,j,k) = rho_0(i,j,k)*u_0(i,j,k)
            E3(i,j,k) = rho_0(i,j,k)*v_0(i,j,k)
            E4(i,j,k) = 0.d0            
            E5(i,j,k) = rho_0(i,j,k)*( 
     &        +tp_0(i,j,k)/((gamma-1.d0)*gamma*Ma**2.d0)
     &        +(u_0(i,j,k)**2.d0+v_0(i,j,k)**2.d0)/2.d0) 
          end do
	end do
      end do

      RETURN
      END SUBROUTINE ini_cond_blayer

      END MODULE sninitial
