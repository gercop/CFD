      MODULE snconf

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      real*8, dimension(:,:,:), allocatable :: 
     &  U1_t,U2_t,U3_t,U4_t,U5_t,
     &  U1_RK1,U2_RK1,U3_RK1,U4_RK1,U5_RK1,
     &  U1_RK2,U2_RK2,U3_Rk2,U4_RK2,U5_RK2

      real*8, dimension(:,:,:), allocatable :: 
     &  u,v,w,tp,p,mu,du_dx,du_dy,du_dz,dv_dx,
     &  dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,dtp_dx,dtp_dy,dtp_dz

      contains 
ccccc **********************************************************************
ccccc RK4: The 4th order temporary integration method for conservative 
ccccc      form of Navier-Stokes equations.
ccccc **********************************************************************
      SUBROUTINE RK4(local_imax,x,y,z,tt)
          
      USE sninitial, only: U1,U2,U3,U4,U5
      USE snfilter,  only: filter_x,filter_y,filter_z

      integer i,j,k,tt,local_imax
      real*8 t,dt2,dt6,x(local_imax),y(jmax),z(kmax)

      dt2=dt*0.5d0
      dt6=dt/6.d0			      

      allocate(U1_t  (local_imax,jmax,kmax),U2_t  (local_imax,jmax,kmax),
     &         U3_t  (local_imax,jmax,kmax),U4_t  (local_imax,jmax,kmax),
     &         U5_t  (local_imax,jmax,kmax),U1_RK1(local_imax,jmax,kmax),
     &         U2_RK1(local_imax,jmax,kmax),U3_RK1(local_imax,jmax,kmax),
     &         U4_RK1(local_imax,jmax,kmax),U5_RK1(local_imax,jmax,kmax),
     &         U1_RK2(local_imax,jmax,kmax),U2_RK2(local_imax,jmax,kmax),
     &         U3_RK2(local_imax,jmax,kmax),U4_RK2(local_imax,jmax,kmax),
     &         U5_RK2(local_imax,jmax,kmax))

      t = dble(tt)*dt
      call derivs(local_imax,x,y,z,t,U1,U2,U3,U4,U5)  
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            U1_t(i,j,k) = U1(i,j,k) + dt2*U1_RK1(i,j,k) 
            U2_t(i,j,k) = U2(i,j,k) + dt2*U2_RK1(i,j,k) 
            U3_t(i,j,k) = U3(i,j,k) + dt2*U3_RK1(i,j,k) 
            U4_t(i,j,k) = U4(i,j,k) + dt2*U4_RK1(i,j,k) 
            U5_t(i,j,k) = U5(i,j,k) + dt2*U5_RK1(i,j,k) 

            U1_RK2(i,j,k) = U1_RK1(i,j,k) 
            U2_RK2(i,j,k) = U2_RK1(i,j,k) 
            U3_RK2(i,j,k) = U3_RK1(i,j,k) 
            U4_RK2(i,j,k) = U4_RK1(i,j,k) 
            U5_RK2(i,j,k) = U5_RK1(i,j,k) 
          end do
        end do        
      end do

      t = dble(tt)*(dt+dt2)
      call derivs(local_imax,x,y,z,t,U1_t,U2_t,U3_t,U4_t,U5_t)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            U1_t(i,j,k) = U1(i,j,k) + dt2*U1_RK1(i,j,k) 
            U2_t(i,j,k) = U2(i,j,k) + dt2*U2_RK1(i,j,k) 
            U3_t(i,j,k) = U3(i,j,k) + dt2*U3_RK1(i,j,k)
            U4_t(i,j,k) = U4(i,j,k) + dt2*U4_RK1(i,j,k)
            U5_t(i,j,k) = U5(i,j,k) + dt2*U5_RK1(i,j,k)

            U1_RK2(i,j,k) = U1_RK2(i,j,k) + 2.d0*U1_RK1(i,j,k)
            U2_RK2(i,j,k) = U2_RK2(i,j,k) + 2.d0*U2_RK1(i,j,k)
            U3_RK2(i,j,k) = U3_RK2(i,j,k) + 2.d0*U3_RK1(i,j,k)			  
            U4_RK2(i,j,k) = U4_RK2(i,j,k) + 2.d0*U4_RK1(i,j,k)			  
            U5_RK2(i,j,k) = U5_RK2(i,j,k) + 2.d0*U5_RK1(i,j,k)			  
          end do
        end do
      end do

      t = dble(tt)*(dt+dt2)
      call derivs(local_imax,x,y,z,t,U1_t,U2_t,U3_t,U4_t,U5_t)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            U1_t(i,j,k) = U1(i,j,k) + dt*U1_RK1(i,j,k)
            U2_t(i,j,k) = U2(i,j,k) + dt*U2_RK1(i,j,k)
            U3_t(i,j,k) = U3(i,j,k) + dt*U3_RK1(i,j,k)
            U4_t(i,j,k) = U4(i,j,k) + dt*U4_RK1(i,j,k)
            U5_t(i,j,k) = U5(i,j,k) + dt*U5_RK1(i,j,k)
		
            U1_RK2(i,j,k) = U1_RK2(i,j,k) + 2.d0*U1_RK1(i,j,k)
            U2_RK2(i,j,k) = U2_RK2(i,j,k) + 2.d0*U2_RK1(i,j,k)
            U3_RK2(i,j,k) = U3_RK2(i,j,k) + 2.d0*U3_RK1(i,j,k)			  
            U4_RK2(i,j,k) = U4_RK2(i,j,k) + 2.d0*U4_RK1(i,j,k)			  
            U5_RK2(i,j,k) = U5_RK2(i,j,k) + 2.d0*U5_RK1(i,j,k)			  
          end do
        end do
      end do

      t = dble(tt)*(2.d0*dt)
      call derivs(local_imax,x,y,z,t,U1_t,U2_t,U3_t,U4_t,U5_t)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            U1(i,j,k) = U1(i,j,k) + dt6*(U1_RK1(i,j,k)+U1_RK2(i,j,k))	    
            U2(i,j,k) = U2(i,j,k) + dt6*(U2_RK1(i,j,k)+U2_RK2(i,j,k))	    
            U3(i,j,k) = U3(i,j,k) + dt6*(U3_RK1(i,j,k)+U3_RK2(i,j,k))        			
            U4(i,j,k) = U4(i,j,k) + dt6*(U4_RK1(i,j,k)+U4_RK2(i,j,k))        			
            U5(i,j,k) = U5(i,j,k) + dt6*(U5_RK1(i,j,k)+U5_RK2(i,j,k))        			
          end do
        end do
      end do

      deallocate(U1_t,U2_t,U3_t,U4_t,U5_t,
     &  U1_RK1,U2_RK1,U3_RK1,U4_RK1,U5_RK1,
     &  U1_RK2,U2_RK2,U3_RK2,U4_RK2,U5_RK2)
 
      call filter_x(local_imax,U1,U2,U3,U4,U5)
      call filter_y(local_imax,U1,U2,U3,U4,U5)
      call filter_z(local_imax,U1,U2,U3,U4,U5)

      RETURN
      END SUBROUTINE RK4

ccccc **************************************************************************
ccccc derivs: Calculate the spatial derivatives of conservation equation. 
ccccc **************************************************************************
      SUBROUTINE derivs(local_imax,x,y,z,t,U1,U2,U3,U4,U5)

      USE snmpi,    only: my_rank,pro,exchange_ghostpoints
      USE snmethod, only: der_x_DD,der_y_WD,der_y_WN,der_z_PP
      
      integer i,j,k,local_imax
      real*8 t,x(local_imax),y(jmax),z(kmax),lambda,C2,
     &  U1(local_imax,jmax,kmax),U2(local_imax,jmax,kmax),
     &  U3(local_imax,jmax,kmax),U4(local_imax,jmax,kmax),
     &  U5(local_imax,jmax,kmax),F (local_imax,jmax,kmax),
     &  G (local_imax,jmax,kmax),H (local_imax,jmax,kmax)
      parameter (lambda = -1.d0/((gamma-1.d0)*Ma**2.d0*Pr*Re))
      parameter (C2=110.4d0/288.15d0) !T_oo=288.15 for sea level. 

      !MS$IF (parallel.EQ.1)
        call exchange_ghostpoints(jmax,kmax,1,U1)
        call exchange_ghostpoints(jmax,kmax,3,U2)
        call exchange_ghostpoints(jmax,kmax,5,U3)
        call exchange_ghostpoints(jmax,kmax,7,U4)
        call exchange_ghostpoints(jmax,kmax,9,U5)
      !MS$ENDIF

      !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5).OR.(tp_s.EQ.6)                  
        call boundcond(local_imax,x,y,z,t,U1,U2,U3,U4,U5)
      !MS$ENDIF 

      allocate(u     (local_imax,jmax,kmax),v     (local_imax,jmax,kmax),
     &         w     (local_imax,jmax,kmax),p     (local_imax,jmax,kmax),
     &         tp    (local_imax,jmax,kmax),mu    (local_imax,jmax,kmax),
     &         du_dx (local_imax,jmax,kmax),du_dy (local_imax,jmax,kmax),
     &         du_dz (local_imax,jmax,kmax),dv_dx (local_imax,jmax,kmax),
     &         dv_dy (local_imax,jmax,kmax),dv_dz (local_imax,jmax,kmax),
     &         dw_dx (local_imax,jmax,kmax),dw_dy (local_imax,jmax,kmax),
     &         dw_dz (local_imax,jmax,kmax),dtp_dx(local_imax,jmax,kmax),
     &         dtp_dy(local_imax,jmax,kmax),dtp_dz(local_imax,jmax,kmax))

ccccc ----------------------------------------------------
ccccc CALCULATING THE VELOCITY COMPONENTS
ccccc ----------------------------------------------------
      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            u (i,j,k) = U2(i,j,k)/U1(i,j,k)
            v (i,j,k) = U3(i,j,k)/U1(i,j,k)
            w (i,j,k) = U4(i,j,k)/U1(i,j,k)
            !MS$IF (tp_proc.EQ.1)
              p (i,j,k) = (gamma-1.d0)*(U5(i,j,k)-1.d0/(2.d0*U1(i,j,k))*
     &          (U2(i,j,k)**2.d0+U3(i,j,k)**2.d0+U4(i,j,k)**2.d0))
            !MS$ELSEIF (tp_proc.EQ.2)
              if (U1(i,j,k).LT.0) then
                write(*,*) 
                write(*,*) 'ERROR..: Negative Density in Point (',i,j,k,').'	       
                write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.'	                    
                STOP          
              end if 
              p (i,j,k) = 1.d0/(gamma*Ma**2.d0)*U1(i,j,k)**gamma
            !MS$ENDIF
            tp(i,j,k) = gamma*Ma**2.d0*p(i,j,k)/U1(i,j,k)

            !MS$IF (tp_viscous.EQ.1)    
              mu(i,j,k) = tp(i,j,k)**(3.d0/2.d0)*((1.d0+C2)/(tp(i,j,k)+C2))            
            !MS$ELSEIF (tp_viscous.EQ.2)
              mu(i,j,k) = 1.d0
            !MS$ENDIF 
          end do
        end do
      end do    

      call der_x_DD(my_rank,pro,local_imax,u, du_dx ) 
      call der_x_DD(my_rank,pro,local_imax,v, dv_dx ) 
      call der_x_DD(my_rank,pro,local_imax,w, dw_dx ) 
      call der_x_DD(my_rank,pro,local_imax,tp,dtp_dx) 

      call der_y_WD(local_imax,u, du_dy )                          
      call der_y_WD(local_imax,v, dv_dy )                          
      call der_y_WD(local_imax,w, dw_dy )                          
      call der_y_WD(local_imax,tp,dtp_dy)                          

      call der_z_PP(local_imax,u, du_dz )                          
      call der_z_PP(local_imax,v, dv_dz )                          
      call der_z_PP(local_imax,w, dw_dz )                          
      call der_z_PP(local_imax,tp,dtp_dz) 

ccccc -------------------------------
ccccc CALCULATING THE FLUX VARIABLES
ccccc -------------------------------     
      call calculating_flux(local_imax,U2,U3,U4,U1_RK1) !Otimização: U2=F1; U3=G1; U4=H1;

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            F(i,j,k) = U2(i,j,k)**2.d0/U1(i,j,k)+p(i,j,k)-2.d0/3.d0*mu(i,j,k)/Re*
     &        (2.d0*du_dx(i,j,k)-dv_dy(i,j,k)-dw_dz(i,j,k))
            G(i,j,k) = U2(i,j,k)*U3(i,j,k)/U1(i,j,k)-mu(i,j,k)/Re*(du_dy(i,j,k)+dv_dx(i,j,k))
            H(i,j,k) = U2(i,j,k)*U4(i,j,k)/U1(i,j,k)-mu(i,j,k)/Re*(du_dz(i,j,k)+dw_dx(i,j,k))
          end do
        end do
      end do    
      call calculating_flux(local_imax,F,G,H,U2_RK1) 

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            F(i,j,k) = U2(i,j,k)*U3(i,j,k)/U1(i,j,k)-mu(i,j,k)/Re*(du_dy(i,j,k)+dv_dx(i,j,k))
            G(i,j,k) = U3(i,j,k)**2.d0/U1(i,j,k)+p(i,j,k)-2.d0/3.d0*mu(i,j,k)/Re*(
     &        2.d0*dv_dy(i,j,k)-du_dx(i,j,k)-dw_dz(i,j,k)) 
            H(i,j,k) = U3(i,j,k)*U4(i,j,k)/U1(i,j,k)-mu(i,j,k)/Re*(dv_dz(i,j,k)+dw_dy(i,j,k))
          end do
        end do
      end do    
      call calculating_flux(local_imax,F,G,H,U3_RK1) 

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            F(i,j,k) = U2(i,j,k)*U4(i,j,k)/U1(i,j,k)-mu(i,j,k)/Re*(du_dz(i,j,k)+dw_dx(i,j,k))
            G(i,j,k) = U3(i,j,k)*U4(i,j,k)/U1(i,j,k)-mu(i,j,k)/Re*(dv_dz(i,j,k)+dw_dy(i,j,k))
            H(i,j,k) = U4(i,j,k)**2.d0/U1(i,j,k)+p(i,j,k)-2.d0/3.d0*mu(i,j,k)/Re*(
     &        2.d0*dw_dz(i,j,k)-du_dx(i,j,k)-dv_dy(i,j,k))
          end do
        end do
      end do    
      call calculating_flux(local_imax,F,G,H,U4_RK1) 

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            F(i,j,k) = (U5(i,j,k)+p(i,j,k))*U2(i,j,k)/U1(i,j,k)
     &        -U2(i,j,k)/U1(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*
     &        (2.d0*du_dx(i,j,k)-dv_dy(i,j,k)-dw_dz(i,j,k))
     &        -U3(i,j,k)/U1(i,j,k)*mu(i,j,k)/Re*(du_dy(i,j,k)+dv_dx(i,j,k))
     &        -U4(i,j,k)/U1(i,j,k)*mu(i,j,k)/Re*(du_dz(i,j,k)+dw_dx(i,j,k))
     &        +lambda*dtp_dx(i,j,k)
            G(i,j,k) = (U5(i,j,k)+p(i,j,k))*U3(i,j,k)/U1(i,j,k)
     &        -U2(i,j,k)/U1(i,j,k)*mu(i,j,k)/Re*(du_dy(i,j,k)+dv_dx(i,j,k))
     &        -U3(i,j,k)/U1(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*
     &        (2.d0*dv_dy(i,j,k)-du_dx(i,j,k)-dw_dz(i,j,k)) 
     &        -U4(i,j,k)/U1(i,j,k)*mu(i,j,k)/Re*(dv_dz(i,j,k)+dw_dy(i,j,k))
     &        +lambda*dtp_dy(i,j,k)
            H(i,j,k) = (U5(i,j,k)+p(i,j,k))*U4(i,j,k)/U1(i,j,k)
     &        -U2(i,j,k)/U1(i,j,k)*mu(i,j,k)/Re*(du_dz(i,j,k)+dw_dx(i,j,k))
     &        -U3(i,j,k)/U1(i,j,k)*mu(i,j,k)/Re*(dv_dz(i,j,k)+dw_dy(i,j,k))
     &        -U4(i,j,k)/U1(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*
     &        (2.d0*dw_dz(i,j,k)-du_dx(i,j,k)-dv_dy(i,j,k))
     &        +lambda*dtp_dz(i,j,k)
          end do
        end do
      end do    
      !MS$IF (tp_proc.EQ.1)
        call calculating_flux(local_imax,F,G,H,U5_RK1) 
      !MS$ELSEIF (tp_proc.EQ.2)
        U5_RK1(:,:,:) = 0.d0
      !MS$ENDIF          

      deallocate(u,v,w,p,tp,mu,du_dx,du_dy,du_dz,dv_dx,dv_dy,
     &  dv_dz,dw_dx,dw_dy,dw_dz,dtp_dx,dtp_dy,dtp_dz)

      RETURN
      END SUBROUTINE derivs

ccccc **************************************************************************
ccccc calculating_flux: Calculate the spatial derivatives. 
ccccc **************************************************************************
      SUBROUTINE calculating_flux(local_imax,F,G,H,U_RK)

      USE snmpi,    only: my_rank,pro
      USE snmethod, only: der_x_DD,der_y_WN,der_z_PP

      integer i,j,k,local_imax
      real*8 U_RK(local_imax,jmax,kmax),
     &  F    (local_imax,jmax,kmax),G    (local_imax,jmax,kmax),
     &  H    (local_imax,jmax,kmax),dF_dx(local_imax,jmax,kmax),
     &  dG_dy(local_imax,jmax,kmax),dH_dz(local_imax,jmax,kmax)

      call der_x_DD(my_rank,pro,local_imax,F,dF_dx)
      call der_y_WN(local_imax,G,dG_dy)
      call der_z_PP(local_imax,H,dH_dz)

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            U_RK(i,j,k) = -dF_dx(i,j,k)-dG_dy(i,j,k)-dH_dz(i,j,k)
          end do
        end do
      end do	

      RETURN
      END SUBROUTINE calculating_flux

ccccc **********************************************************************
ccccc extrapolated_inlet: This routine extrapolated the variables at 
ccccc                     inflow boundaries.
ccccc **********************************************************************
      SUBROUTINE extrapolated_inlet(local_imax,j,k,fc)

      USE snmethod, only: stx_d1_ex_D

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer j,k,local_imax 
      real*8 fc(local_imax,jmax,kmax)

      !Extrapolation of sixth order
      fc(1,j,k) = 1.d0/stx_d1_ex_D(1,1)*(
     &  -stx_d1_ex_D(1,2)*fc(2,j,k)+
     &  -stx_d1_ex_D(1,3)*fc(3,j,k)+
     &  -stx_d1_ex_D(1,4)*fc(4,j,k)+
     &  -stx_d1_ex_D(1,5)*fc(5,j,k)+
     &  -stx_d1_ex_D(1,6)*fc(6,j,k)+
     &  -stx_d1_ex_D(1,7)*fc(7,j,k))

      RETURN
      END SUBROUTINE extrapolated_inlet

ccccc **********************************************************************
ccccc extrapolated_outlet: This routine extrapolated the variables at
ccccc                      outflow boundaries.
ccccc **********************************************************************
      SUBROUTINE extrapolated_outlet(local_imax,j,k,fc)

      USE snmethod, only: stx_d2_ex_D

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer j,k,local_imax
      real*8 fc(local_imax,jmax,kmax)

      !Extrapolation of fifth order
      fc(local_imax,j,k) = 1.d0/stx_d2_ex_D(local_imax,7)*(
     &  -stx_d2_ex_D(local_imax,1)*fc(local_imax-6,j,k)+
     &  -stx_d2_ex_D(local_imax,2)*fc(local_imax-5,j,k)+
     &  -stx_d2_ex_D(local_imax,3)*fc(local_imax-4,j,k)+
     &  -stx_d2_ex_D(local_imax,4)*fc(local_imax-3,j,k)+
     &  -stx_d2_ex_D(local_imax,5)*fc(local_imax-2,j,k)+
     &  -stx_d2_ex_D(local_imax,6)*fc(local_imax-1,j,k))

      RETURN
      END SUBROUTINE extrapolated_outlet

ccccc **********************************************************************
ccccc boundcond: This routine calculated the boundary conditions for spatial 
ccccc            development of shear layer flow.
ccccc **********************************************************************
      SUBROUTINE boundcond(local_imax,x,y,z,t,U1,U2,U3,U4,U5)

      USE snmpi,     only: my_rank,pro
      USE sninitial, only: rho_0,u_0,v_0,tp_0

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer i,j,k,local_imax
      real*8 t,yy,zz,u_1,v_1,w_1,
     &  x(local_imax),y(jmax),z(kmax),
     &  U1(local_imax,jmax,kmax),U2(local_imax,jmax,kmax),
     &  U3(local_imax,jmax,kmax),U4(local_imax,jmax,kmax),
     &  U5(local_imax,jmax,kmax),p (local_imax,jmax,kmax) 
                             
      !MS$IF (tp_s.EQ.5).OR.(tp_s.EQ.6)
        if (my_rank.EQ.0) then  
          do k=1,kmax 
            zz=z(k) 
            do j=1,jmax
              yy=y(j)          

              !MS$IF (tp_s.EQ.5)
                !MS$IF (tp_dim.EQ.2)      !Two-dimensional Problem                   
                  u_1 = 
     &              +02.d0*A_1*sigma/alpha*yy*dsin(omega*t/1.d0)*dexp(-sigma*yy**2.d0)
     &              +04.d0*A_2*sigma/alpha*yy*dsin(omega*t/2.d0)*dexp(-sigma*yy**2.d0)
     &              +08.d0*A_3*sigma/alpha*yy*dsin(omega*t/4.d0)*dexp(-sigma*yy**2.d0)
     &              +16.d0*A_4*sigma/alpha*yy*dsin(omega*t/8.d0)*dexp(-sigma*yy**2.d0)
                  v_1 = 
     &              +A_1*dcos(omega/1.d0*t)*dexp(-sigma*yy**2.d0)
     &              +A_2*dcos(omega/2.d0*t)*dexp(-sigma*yy**2.d0)
     &              +A_3*dcos(omega/4.d0*t)*dexp(-sigma*yy**2.d0)
     &              +A_4*dcos(omega/8.d0*t)*dexp(-sigma*yy**2.d0)
                  w_1 = 0.d0                  
                !MS$ELSEIF (tp_dim.EQ.3)  !Three-dimensional Problem                  
                  u_1 = 
     &              +2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*(
     &              +1.d0*A_1/alpha*dsin(-omega*t+phi)
     &              +1.d0*A_2/alpha*dsin(beta/1.d0*zz-omega/1.d0*t)+1.d0*A_2/alpha*dsin(-beta/1.d0*zz-omega/1.d0*t)
     &              +2.d0*A_3/alpha*dsin(beta/2.d0*zz-omega/2.d0*t)+2.d0*A_3/alpha*dsin(-beta/2.d0*zz-omega/2.d0*t)
     &              +4.d0*A_4/alpha*dsin(beta/4.d0*zz-omega/4.d0*t)+4.d0*A_4/alpha*dsin(-beta/4.d0*zz-omega/4.d0*t)) 
     &              -dexp(-sigma*yy**2.d0)*(
     &              +A_2*beta/alpha*dcos(-beta/1.d0*zz+omega/1.d0*t)+A_2*beta/alpha*dcos(beta/1.d0*zz+omega/1.d0*t)
     &              +A_3*beta/alpha*dcos(-beta/2.d0*zz+omega/2.d0*t)+A_3*beta/alpha*dcos(beta/2.d0*zz+omega/2.d0*t)
     &              +A_4*beta/alpha*dcos(-beta/4.d0*zz+omega/4.d0*t)+A_4*beta/alpha*dcos(beta/4.d0*zz+omega/4.d0*t))
                   v_1 = 
     &               (A_1*dcos(-omega*t+phi)
     &               +A_2*dcos(beta/1.d0*zz-omega/1.d0*t)+A_2*dcos(beta/1.d0*zz+omega/1.d0*t)
     &               +A_3*dcos(beta/2.d0*zz-omega/2.d0*t)+A_3*dcos(beta/2.d0*zz+omega/2.d0*t)
     &               +A_4*dcos(beta/4.d0*zz-omega/4.d0*t)+A_4*dcos(beta/4.d0*zz+omega/4.d0*t))*dexp(-sigma*yy**2)
                   w_1 = 
     &               (A_1*dcos(-omega*t+phi)
     &               +A_2*dcos(beta/1.d0*zz-omega/1.d0*t)+A_2*dcos(beta/1.d0*zz+omega/1.d0*t)
     &               +A_3*dcos(beta/2.d0*zz-omega/2.d0*t)+A_3*dcos(beta/2.d0*zz+omega/2.d0*t)
     &               +A_4*dcos(beta/4.d0*zz-omega/4.d0*t)+A_4*dcos(beta/4.d0*zz+omega/4.d0*t))*dexp(-sigma*yy**2)
                !MS$ENDIF
              !MS$ELSEIF (tp_s.EQ.6)
                u_1 = 0.d0
                v_1 = 0.d0
                w_1 = 0.d0
              !MS$ENDIF               

              do i=1,7
                !MS$IF (tp_proc.EQ.1)
                  p (i,j,k) = (gamma-1.d0)*(U5(i,j,k)-1.d0/(2.d0*U1(i,j,k))*(
     &              U2(i,j,k)**2.d0+U3(i,j,k)**2.d0+U4(i,j,k)**2.d0))
                !MS$ELSEIF (tp_proc.EQ.2)
                  p (i,j,k) = 1.d0/(gamma*Ma**2.d0)*U1(i,j,k)**gamma
                !MS$ENDIF
              end do
              call extrapolated_inlet(local_imax,j,k,p)   !Pressure is extrapolated to satisfy dp_dx=0

              U1(1,j,k) = rho_0(1,j,k)           
              U2(1,j,k) = rho_0(1,j,k)*(u_0(1,j,k)+1.d0/(Ma*c)*u_1) 
              U3(1,j,k) = rho_0(1,j,k)*(v_0(1,j,k)+1.d0/(Ma*c)*v_1)
              U4(1,j,k) = rho_0(1,j,k)*(          +1.d0/(Ma*c)*w_1)
              U5(1,j,k) = p(1,j,k)/(gamma-1.d0)+1.d0/(2.d0*U1(1,j,k))*(
     &          U2(1,j,k)**2.d0+U3(1,j,k)**2.d0+U4(1,j,k)**2.d0)
            end do
          end do
        end if
  
        if (my_rank+1.EQ.pro) then
          do k=1,kmax 
            do j=1,jmax
              p (local_imax,j,k) = 1.d0/(gamma*Ma**2.d0)  !Pressure is fixed in x(imax)

              call extrapolated_outlet(local_imax,j,k,U1) !Density is extrapolated to satisfy d2rho_dx2=0
              call extrapolated_outlet(local_imax,j,k,U2) !Velocity is extrapolated to satisfy d2u_dx2=0
              call extrapolated_outlet(local_imax,j,k,U3) !Velocity is extrapolated to satisfy d2v_dx2=0
              call extrapolated_outlet(local_imax,j,k,U4) !Velocity is extrapolated to satisfy d2w_dx2=0

              U5(local_imax,j,k) = p(local_imax,j,k)/(gamma-1.d0)
     &          +1.d0/(2.d0*U1(local_imax,j,k))*(U2(local_imax,j,k)**2.d0
     &          +U3(local_imax,j,k)**2.d0+U4(local_imax,j,k)**2.d0)          
            end do
          end do 
        end if
      !MS$ENDIF       

      !Free-slip Boundary Condition - Impermeability
      !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        do k=1,kmax
          do i=1,local_imax
            U3(i,1,   k) = 0.d0
            U3(i,jmax,k) = 0.d0            
          end do
        end do
      !MS$ENDIF       

      !Wall-slip boundary condition
      !MS$IF (tp_s.EQ.6)
        do k=1,kmax
          do i=1,local_imax
            U2(i,1,k) = 0.d0
            U3(i,1,k) = 0.d0
            U4(i,1,k) = 0.d0
          end do
        end do
      !MS$ENDIF       

      RETURN
      END SUBROUTINE boundcond

      END MODULE snconf 
