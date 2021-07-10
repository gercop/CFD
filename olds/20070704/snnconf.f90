      MODULE snconf

      IMPLICIT NONE
      INCLUDE 'snpar.f90'

      real*8, dimension(:,:,:), allocatable :: 
     &  E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_RK1,E2_RK1,E3_RK1,E4_RK1,E5_RK1,
     &  E1_RK2,E2_RK2,E3_Rk2,E4_RK2,E5_RK2

      contains 
ccccc **********************************************************************
ccccc RK4: The 4th order temporary integration method for conservative 
ccccc      form of Navier-Stokes equations.
ccccc **********************************************************************
      SUBROUTINE RK4(local_imax,x,y,z,tt)
          
      USE sninitial, only: E1,E2,E3,E4,E5
      USE snfilter,  only: filter_x,filter_y,filter_z

      integer i,j,k,tt,local_imax
      real*8 t,dt2,dt6,x(local_imax),y(jmax),z(kmax)

      dt2=dt*0.5d0
      dt6=dt/6.d0			      

      allocate(E1_t  (local_imax,jmax,kmax),E2_t  (local_imax,jmax,kmax),
     &         E3_t  (local_imax,jmax,kmax),E4_t  (local_imax,jmax,kmax),
     &         E5_t  (local_imax,jmax,kmax),E1_RK1(local_imax,jmax,kmax),
     &         E2_RK1(local_imax,jmax,kmax),E3_RK1(local_imax,jmax,kmax),
     &         E4_RK1(local_imax,jmax,kmax),E5_RK1(local_imax,jmax,kmax),
     &         E1_RK2(local_imax,jmax,kmax),E2_RK2(local_imax,jmax,kmax),
     &         E3_RK2(local_imax,jmax,kmax),E4_RK2(local_imax,jmax,kmax),
     &         E5_RK2(local_imax,jmax,kmax))

      t = dble(tt)*dt
      call derivs(local_imax,x,y,z,t,E1,E2,E3,E4,E5)  
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            E1_t(i,j,k) = E1(i,j,k) + dt2*E1_RK1(i,j,k) 
            E2_t(i,j,k) = E2(i,j,k) + dt2*E2_RK1(i,j,k) 
            E3_t(i,j,k) = E3(i,j,k) + dt2*E3_RK1(i,j,k) 
            E4_t(i,j,k) = E4(i,j,k) + dt2*E4_RK1(i,j,k) 
            E5_t(i,j,k) = E5(i,j,k) + dt2*E5_RK1(i,j,k) 

            E1_RK2(i,j,k) = E1_RK1(i,j,k) 
            E2_RK2(i,j,k) = E2_RK1(i,j,k) 
            E3_RK2(i,j,k) = E3_RK1(i,j,k) 
            E4_RK2(i,j,k) = E4_RK1(i,j,k) 
            E5_RK2(i,j,k) = E5_RK1(i,j,k) 
          end do
        end do        
      end do

      t = dble(tt)*(dt+dt2)
      call derivs(local_imax,x,y,z,t,E1_t,E2_t,E3_t,E4_t,E5_t)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            E1_t(i,j,k) = E1(i,j,k) + dt2*E1_RK1(i,j,k) 
            E2_t(i,j,k) = E2(i,j,k) + dt2*E2_RK1(i,j,k) 
            E3_t(i,j,k) = E3(i,j,k) + dt2*E3_RK1(i,j,k)
            E4_t(i,j,k) = E4(i,j,k) + dt2*E4_RK1(i,j,k)
            E5_t(i,j,k) = E5(i,j,k) + dt2*E5_RK1(i,j,k)

            E1_RK2(i,j,k) = E1_RK2(i,j,k) + 2.d0*E1_RK1(i,j,k)
            E2_RK2(i,j,k) = E2_RK2(i,j,k) + 2.d0*E2_RK1(i,j,k)
            E3_RK2(i,j,k) = E3_RK2(i,j,k) + 2.d0*E3_RK1(i,j,k)			  
            E4_RK2(i,j,k) = E4_RK2(i,j,k) + 2.d0*E4_RK1(i,j,k)			  
            E5_RK2(i,j,k) = E5_RK2(i,j,k) + 2.d0*E5_RK1(i,j,k)			  
          end do
        end do
      end do

      t = dble(tt)*(dt+dt2)
      call derivs(local_imax,x,y,z,t,E1_t,E2_t,E3_t,E4_t,E5_t)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            E1_t(i,j,k) = E1(i,j,k) + dt*E1_RK1(i,j,k)
            E2_t(i,j,k) = E2(i,j,k) + dt*E2_RK1(i,j,k)
            E3_t(i,j,k) = E3(i,j,k) + dt*E3_RK1(i,j,k)
            E4_t(i,j,k) = E4(i,j,k) + dt*E4_RK1(i,j,k)
            E5_t(i,j,k) = E5(i,j,k) + dt*E5_RK1(i,j,k)
		
            E1_RK2(i,j,k) = E1_RK2(i,j,k) + 2.d0*E1_RK1(i,j,k)
            E2_RK2(i,j,k) = E2_RK2(i,j,k) + 2.d0*E2_RK1(i,j,k)
            E3_RK2(i,j,k) = E3_RK2(i,j,k) + 2.d0*E3_RK1(i,j,k)			  
            E4_RK2(i,j,k) = E4_RK2(i,j,k) + 2.d0*E4_RK1(i,j,k)			  
            E5_RK2(i,j,k) = E5_RK2(i,j,k) + 2.d0*E5_RK1(i,j,k)			  
          end do
        end do
      end do

      t = dble(tt)*(2.d0*dt)
      call derivs(local_imax,x,y,z,t,E1_t,E2_t,E3_t,E4_t,E5_t)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            E1(i,j,k) = E1(i,j,k) + dt6*(E1_RK1(i,j,k)+E1_RK2(i,j,k))	    
            E2(i,j,k) = E2(i,j,k) + dt6*(E2_RK1(i,j,k)+E2_RK2(i,j,k))	    
            E3(i,j,k) = E3(i,j,k) + dt6*(E3_RK1(i,j,k)+E3_RK2(i,j,k))        			
            E4(i,j,k) = E4(i,j,k) + dt6*(E4_RK1(i,j,k)+E4_RK2(i,j,k))        			
            E5(i,j,k) = E5(i,j,k) + dt6*(E5_RK1(i,j,k)+E5_RK2(i,j,k))        			
          end do
        end do
      end do

      deallocate(E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_RK1,E2_RK1,E3_RK1,E4_RK1,E5_RK1,
     &  E1_RK2,E2_RK2,E3_RK2,E4_RK2,E5_RK2)
 
      call filter_x(local_imax,E1,E2,E3,E4,E5)
      call filter_y(local_imax,E1,E2,E3,E4,E5)
      call filter_z(local_imax,E1,E2,E3,E4,E5)

      RETURN
      END SUBROUTINE RK4

ccccc **************************************************************************
ccccc derivs: Calculate the spatial derivatives of conservation equation. 
ccccc **************************************************************************
      SUBROUTINE derivs(local_imax,x,y,z,t,E1,E2,E3,E4,E5)

      USE snmpi,    only: exchange_ghostpoints
      USE snmethod, only: der_x_DD,der_y_WD,der_y_WN,der_z_PP
      
      integer i,j,k,local_imax
      real*8 t,x(local_imax),y(jmax),z(kmax),lambda,C2,
     &  E1    (local_imax,jmax,kmax),E2    (local_imax,jmax,kmax),
     &  E3    (local_imax,jmax,kmax),E4    (local_imax,jmax,kmax),
     &  E5    (local_imax,jmax,kmax),F     (local_imax,jmax,kmax),
     &  G     (local_imax,jmax,kmax),H     (local_imax,jmax,kmax),
     &  u     (local_imax,jmax,kmax),v     (local_imax,jmax,kmax),
     &  w     (local_imax,jmax,kmax),tp    (local_imax,jmax,kmax),
     &  p     (local_imax,jmax,kmax),mu    (local_imax,jmax,kmax),
     &  du_dx (local_imax,jmax,kmax),du_dy (local_imax,jmax,kmax),
     &  du_dz (local_imax,jmax,kmax),dv_dx (local_imax,jmax,kmax),
     &  dv_dy (local_imax,jmax,kmax),dv_dz (local_imax,jmax,kmax),
     &  dw_dx (local_imax,jmax,kmax),dw_dy (local_imax,jmax,kmax),
     &  dw_dz (local_imax,jmax,kmax),dtp_dx(local_imax,jmax,kmax),
     &  dtp_dy(local_imax,jmax,kmax),dtp_dz(local_imax,jmax,kmax),
     &  dE1_dx(local_imax,jmax,kmax),dE2_dx(local_imax,jmax,kmax),
     &  dE1_dy(local_imax,jmax,kmax),dE3_dy(local_imax,jmax,kmax),
     &  dE1_dz(local_imax,jmax,kmax),dE4_dz(local_imax,jmax,kmax)
      parameter (lambda = -1.d0/((gamma-1.d0)*Ma**2.d0*Pr*Re))
      parameter (C2=110.4d0/288.15d0) !T_oo=288.15 for sea level. 

      !MS$IF (parallel.EQ.1)
        call exchange_ghostpoints(jmax,kmax,1,E1)
        call exchange_ghostpoints(jmax,kmax,3,E2)
        call exchange_ghostpoints(jmax,kmax,5,E3)
        call exchange_ghostpoints(jmax,kmax,7,E4)
        call exchange_ghostpoints(jmax,kmax,9,E5)
      !MS$ENDIF

ccccc ----------------------------------------------------
ccccc CALCULATING THE VELOCITY COMPONENTS
ccccc ----------------------------------------------------
      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            u (i,j,k) = E2(i,j,k)/E1(i,j,k)
            v (i,j,k) = E3(i,j,k)/E1(i,j,k)
            w (i,j,k) = E4(i,j,k)/E1(i,j,k)
            !MS$IF (tp_proc.EQ.1)
              p (i,j,k) = (gamma-1.d0)*(E5(i,j,k)
     &          -1.d0/(2.d0*E1(i,j,k))*(E2(i,j,k)**2.d0
     &          +E3(i,j,k)**2.d0+E4(i,j,k)**2.d0))
            !MS$ELSEIF (tp_proc.EQ.2)
              if (E1(i,j,k).LT.0) then
                write(*,*) 
                write(*,*) 'ERROR..: Negative Density in Point (',i,j,k,').'	       
                write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.'	                    
                STOP          
              end if 
              p (i,j,k) = 1.d0/(gamma*Ma**2.d0)*E1(i,j,k)**gamma
            !MS$ENDIF
            tp(i,j,k) = gamma*Ma**2.d0*p(i,j,k)/E1(i,j,k)

            !MS$IF (tp_viscous.EQ.1)    
              mu(i,j,k) = tp(i,j,k)**(3.d0/2.d0)*((1.d0+C2)/(tp(i,j,k)+C2))            
            !MS$ELSEIF (tp_viscous.EQ.2)
              mu(i,j,k) = 1.d0
            !MS$ENDIF 
          end do
        end do
      end do    

      !MS$IF (tp_s.EQ.5).OR.(tp_s.EQ.6)                  
        call boundcond(local_imax,x,y,z,t,E1,E2,E3,E4,E5,p,tp,mu)
      !MS$ENDIF 

      call der_x_DD(local_imax,u, du_dx ) 
      call der_x_DD(local_imax,v, dv_dx ) 
      call der_x_DD(local_imax,w, dw_dx ) 
      call der_x_DD(local_imax,tp,dtp_dx) 
      call der_x_DD(local_imax,E1,dE1_dx) 
      call der_x_DD(local_imax,E2,dE2_dx) 

      call der_y_WN(local_imax,u, du_dy )                          
      call der_y_WN(local_imax,v, dv_dy )                          
      call der_y_WN(local_imax,w, dw_dy )                          
      call der_y_WN(local_imax,tp,dtp_dy)                          
      call der_x_DD(local_imax,E1,dE1_dy) 
      call der_x_DD(local_imax,E3,dE3_dy) 

      call der_z_PP(local_imax,u, du_dz )                          
      call der_z_PP(local_imax,v, dv_dz )                          
      call der_z_PP(local_imax,w, dw_dz )                          
      call der_z_PP(local_imax,tp,dtp_dz) 
      call der_x_DD(local_imax,E1,dE1_dz) 
      call der_x_DD(local_imax,E4,dE4_dz) 

ccccc -------------------------------
ccccc CALCULATING THE FLUX VARIABLES
ccccc -------------------------------     
      call calculating_flux(local_imax,E2,E3,E4,E1_RK1) !Otimização: E2=F1; E3=G1; E4=H1;

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            F(i,j,k) = E2(i,j,k)**2.d0/E1(i,j,k)+p(i,j,k)-2.d0/3.d0*mu(i,j,k)/Re*
     &        (2.d0*(1.d0/E1(i,j,k)*dE2_dx(i,j,k)-E2(i,j,k)/(E1(i,j,k)**2.d0)*dE1_dx(i,j,k))
     &        -(1.d0/E1(i,j,k)*dE3_dy(i,j,k)-E3(i,j,k)/(E1(i,j,k)**2.d0)*dE1_dy(i,j,k))
     &        -(1.d0/E1(i,j,k)*dE4_dz(i,j,k)-E4(i,j,k)/(E1(i,j,k)**2.d0)*dE1_dz(i,j,k)) )
            G(i,j,k) = E2(i,j,k)*E3(i,j,k)/E1(i,j,k)-mu(i,j,k)/Re*(
     &        +(1.d0/E1(i,j,k)*dE2_dy(i,j,k)-E2(i,j,k)/(E1(i,j,k)**2.d0)*dE1_dy(i,j,k))
     &        +dv_dx(i,j,k))
            H(i,j,k) = E2(i,j,k)*E4(i,j,k)/E1(i,j,k)-mu(i,j,k)/Re*(du_dz(i,j,k)+dw_dx(i,j,k))
          end do
        end do
      end do    
      call calculating_flux(local_imax,F,G,H,E2_RK1) 

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            F(i,j,k) = E2(i,j,k)*E3(i,j,k)/E1(i,j,k)-mu(i,j,k)/Re*(du_dy(i,j,k)+dv_dx(i,j,k))
            G(i,j,k) = E3(i,j,k)**2.d0/E1(i,j,k)+p(i,j,k)-2.d0/3.d0*mu(i,j,k)/Re*(
     &        2.d0*dv_dy(i,j,k)-du_dx(i,j,k)-dw_dz(i,j,k)) 
            H(i,j,k) = E3(i,j,k)*E4(i,j,k)/E1(i,j,k)-mu(i,j,k)/Re*(dv_dz(i,j,k)+dw_dy(i,j,k))
          end do
        end do
      end do    
      call calculating_flux(local_imax,F,G,H,E3_RK1) 

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            F(i,j,k) = E2(i,j,k)*E4(i,j,k)/E1(i,j,k)-mu(i,j,k)/Re*(du_dz(i,j,k)+dw_dx(i,j,k))
            G(i,j,k) = E3(i,j,k)*E4(i,j,k)/E1(i,j,k)-mu(i,j,k)/Re*(dv_dz(i,j,k)+dw_dy(i,j,k))
            H(i,j,k) = E4(i,j,k)**2.d0/E1(i,j,k)+p(i,j,k)-2.d0/3.d0*mu(i,j,k)/Re*(
     &        2.d0*dw_dz(i,j,k)-du_dx(i,j,k)-dv_dy(i,j,k))
          end do
        end do
      end do    
      call calculating_flux(local_imax,F,G,H,E4_RK1) 

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            F(i,j,k) = (E5(i,j,k)+p(i,j,k))*E2(i,j,k)/E1(i,j,k)
     &        -E2(i,j,k)/E1(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*
     &        (2.d0*du_dx(i,j,k)-dv_dy(i,j,k)-dw_dz(i,j,k))
     &        -E3(i,j,k)/E1(i,j,k)*mu(i,j,k)/Re*(du_dy(i,j,k)+dv_dx(i,j,k))
     &        -E4(i,j,k)/E1(i,j,k)*mu(i,j,k)/Re*(du_dz(i,j,k)+dw_dx(i,j,k))
     &        +lambda*dtp_dx(i,j,k)
            G(i,j,k) = (E5(i,j,k)+p(i,j,k))*E3(i,j,k)/E1(i,j,k)
     &        -E2(i,j,k)/E1(i,j,k)*mu(i,j,k)/Re*(du_dy(i,j,k)+dv_dx(i,j,k))
     &        -E3(i,j,k)/E1(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*
     &        (2.d0*dv_dy(i,j,k)-du_dx(i,j,k)-dw_dz(i,j,k)) 
     &        -E4(i,j,k)/E1(i,j,k)*mu(i,j,k)/Re*(dv_dz(i,j,k)+dw_dy(i,j,k))
     &        +lambda*dtp_dy(i,j,k)
            H(i,j,k) = (E5(i,j,k)+p(i,j,k))*E4(i,j,k)/E1(i,j,k)
     &        -E2(i,j,k)/E1(i,j,k)*mu(i,j,k)/Re*(du_dz(i,j,k)+dw_dx(i,j,k))
     &        -E3(i,j,k)/E1(i,j,k)*mu(i,j,k)/Re*(dv_dz(i,j,k)+dw_dy(i,j,k))
     &        -E4(i,j,k)/E1(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*
     &        (2.d0*dw_dz(i,j,k)-du_dx(i,j,k)-dv_dy(i,j,k))
     &        +lambda*dtp_dz(i,j,k)
          end do
        end do
      end do    
      !MS$IF (tp_proc.EQ.1)
        call calculating_flux(local_imax,F,G,H,E5_RK1) 
      !MS$ELSEIF (tp_proc.EQ.2)
        E5_RK1(:,:,:) = 0.d0
      !MS$ENDIF          

      RETURN
      END SUBROUTINE derivs

ccccc **************************************************************************
ccccc calculating_flux: Calculate the spatial derivatives of conservation equation. 
ccccc **************************************************************************
      SUBROUTINE calculating_flux(local_imax,F,G,H,E_RK)

      USE snmethod, only: der_x_DD,der_y_WN,der_z_PP

      integer i,j,k,local_imax
      real*8 E_RK(local_imax,jmax,kmax),
     &  F    (local_imax,jmax,kmax),G    (local_imax,jmax,kmax),
     &  H    (local_imax,jmax,kmax),dF_dx(local_imax,jmax,kmax),
     &  dG_dy(local_imax,jmax,kmax),dH_dz(local_imax,jmax,kmax)

      call der_x_DD(local_imax,F,dF_dx)
      call der_y_WN(local_imax,G,dG_dy)
      call der_z_PP(local_imax,H,dH_dz)

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax
            E_RK(i,j,k) = -dF_dx(i,j,k)-dG_dy(i,j,k)-dH_dz(i,j,k)
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
      INCLUDE 'snpar.f90'

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
      INCLUDE 'snpar.f90'

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
ccccc boundcond: This routine calculated the boundary conditions for 
ccccc            spatial development of shear layer flow.
ccccc **********************************************************************
      SUBROUTINE boundcond(local_imax,x,y,z,t,E1,E2,E3,E4,E5,p,tp,mu)

      USE snmpi,     only: my_rank,pro
      USE sninitial, only: rho_0,u_0,v_0,tp_0

      IMPLICIT NONE
      INCLUDE 'snpar.f90'

      integer i,j,k,local_imax
      real*8 t,yy,zz,u_1,v_1,w_1,phi,C2,
     &  x(local_imax),y(jmax),z(kmax),
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax),p (local_imax,jmax,kmax),
     &  tp(local_imax,jmax,kmax),mu(local_imax,jmax,kmax)
      parameter (C2=110.4d0/288.15d0) !T_oo=288.15 for sea level.
     
      !MS$IF (tp_s.EQ.5).OR.(tp_s.EQ.6)
        do k=1,kmax 
          zz=z(k) 
          do j=1,jmax
            yy=y(j)          

            if (my_rank.EQ.0) then  
              !MS$IF (tp_s.EQ.5)
                if ((jmax.NE.1).AND.(kmax.EQ.1)) then      !Two-Dimensional Problem
                  u_1 = 
     &              -2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*
     &                A_1*dsin(-omega/1.d0*t)/alpha_1
     &              -2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*
     &                A_2*dsin(-omega/2.d0*t)/alpha_1
     &              -2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*
     &                A_3*dsin(-omega/4.d0*t)/alpha_1
                  v_1 = 
     &              +A_1*dcos(-omega/1.d0*t)*dexp(-sigma*yy**2.d0)
     &              +A_2*dcos(-omega/2.d0*t)*dexp(-sigma*yy**2.d0)
     &              +A_3*dcos(-omega/4.d0*t)*dexp(-sigma*yy**2.d0)
                  w_1 = 0.d0                  
                else if ((jmax.NE.1).AND.(kmax.NE.1)) then !Three-Dimensional Problem
                  phi = 1.d0/2.d0*pi                             
                  u_1 = 
     &              -dexp(-sigma*yy**2.d0)*(
     &                +A_2*beta*dcos(beta*zz-omega*t)/alpha_1
     &                -A_2*beta*dcos(beta*zz-omega*t)/alpha_1)
     &              +2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*(
     &                +A_1*dsin(-omega*t+phi)/alpha_1+ 
     &                +A_2*dsin(+beta*zz-omega*t)/alpha_1
     &                +A_2*dsin(-beta*zz-omega*t)/alpha_1)
                  v_1 = 
     &              (A_1*dcos(-omega*t+phi)+
     &               A_2*dcos(+beta*zz-omega*t)+
     &               A_2*dcos(-beta*zz-omega*t))*
     &              dexp(-sigma*yy**2.d0)
                  w_1 = 
     &              (A_1*dcos(-omega*t+phi)+
     &               A_2*dcos(+beta*zz-omega*t)+
     &               A_2*dcos(-beta*zz-omega*t))*
     &              dexp(-sigma*yy**2.d0)
                end if
              !MS$ELSEIF (tp_s.EQ.6)
                u_1 = 0.d0
                v_1 = 0.d0
                w_1 = 0.d0
              !MS$ENDIF               

              E1(1,j,k) = rho_0(1,j,k)           
              E2(1,j,k) = rho_0(1,j,k)*(u_0(1,j,k)+1.d0/(M1*c)*u_1) 
              E3(1,j,k) = rho_0(1,j,k)*(v_0(1,j,k)+1.d0/(M1*c)*v_1)
              E4(1,j,k) = rho_0(1,j,k)*(          +1.d0/(M1*c)*w_1)

              tp(1,j,k) = tp_0(1,j,k)     

              call extrapolated_inlet(local_imax,j,k,p)   !Pressure is extrapolated to satisfy dp_dx=0

              E5 (1,j,k) = E1(1,j,k)*(p(1,j,k)/(E1(1,j,k)*(gamma-1.d0))+1.d0/2.d0*(
     &          (E2(1,j,k)/E1(1,j,k))**2.d0+(E3(1,j,k)/E1(1,j,k))**2.d0+
     &          (E4(1,j,k)/E1(1,j,k))**2.d0))
              !MS$IF (tp_viscous.EQ.1)
                mu(1,j,k) = tp(1,j,k)**(3.d0/2.d0)*((1.d0+C2)
     &            /(tp(1,j,k)+C2))
              !MS$ELSEIF (tp_viscous.EQ.2)
                mu(1,j,k) = 1.d0
              !MS$ENDIF 
            end if         

            if (my_rank+1.EQ.pro) then
              call extrapolated_outlet(local_imax,j,k,E1) !Density is extrapolated to satisfy d2rho_dx2=0
              call extrapolated_outlet(local_imax,j,k,E2) !Velocity is extrapolated to satisfy d2u_dx2=0
              call extrapolated_outlet(local_imax,j,k,E3) !Velocity is extrapolated to satisfy d2v_dx2=0
              call extrapolated_outlet(local_imax,j,k,E4) !Velocity is extrapolated to satisfy d2w_dx2=0
              call extrapolated_outlet(local_imax,j,k,tp) !Temperature is extrapolated to satisfy d2tp_dx2=0
           
              p(local_imax,j,k) = 1.d0/(gamma*Ma**2.d0)   !Pressure is fixed in x(imax)

              E5(local_imax,j,k) = E1(local_imax,j,k)*(
     &          p(local_imax,j,k)/(E1(local_imax,j,k)*(gamma-1.d0))
     &          +0.5d0*((E2(local_imax,j,k)/E2(local_imax,j,k))**2.d0
     &          +(E3(local_imax,j,k)/E2(local_imax,j,k))**2.d0
     &          +(E4(local_imax,j,k)/E2(local_imax,j,k))**2.d0)) 
              !MS$IF (tp_viscous.EQ.1)    
                mu(local_imax,j,k) = tp(local_imax,j,k)**(3.d0/2.d0)*
     &            ((1.d0+C2)/(tp(local_imax,j,k)+C2))
              !MS$ELSEIF (tp_viscous.EQ.2)
                mu(local_imax,j,k) = 1.d0
              !MS$ENDIF 
            end if
          end do
        end do 
      !MS$ENDIF       

      !Free-slip Boundary Condition - Impermeability
      !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        do k=1,kmax
          do i=1,local_imax
            E3(i,1,   k) = 0.d0
            E3(i,jmax,k) = 0.d0            
          end do
        end do
      !MS$ENDIF       

      !Wall-slip boundary condition
      !MS$IF (tp_s.EQ.6)
        do k=1,kmax
          do i=1,local_imax
            E2(i,1,k) = 0.d0
            E3(i,1,k) = 0.d0
            E4(i,1,k) = 0.d0
          end do
        end do
      !MS$ENDIF       

      RETURN
      END SUBROUTINE boundcond

      END MODULE snconf 
