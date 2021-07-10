      MODULE snconf

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      real*8, dimension(:,:,:), allocatable :: 
     &  E1_t,E2_t,E3_t,E4_t,E5_t,E1_rk1,E1_rk2,E1_rk3,E1_rk4,
     &  E2_rk1,E2_rk2,E2_rk3,E2_rk4,E3_rk1,E3_rk2,E3_rk3,E3_rk4, 
     &  E4_rk1,E4_rk2,E4_rk3,E4_rk4,E5_rk1,E5_rk2,E5_rk3,E5_rk4,
     &  du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,dp_dx,
     &  dp_dy,dp_dz,drho_dx,drho_dy,drho_dz,dEt_dx,dEt_dy,dEt_dz,
     &  dtp_dx,dtp_dy,dtp_dz,d2u_dx2,d2u_dy2,d2u_dz2,d2v_dx2,d2v_dy2,
     &  d2v_dz2,d2w_dx2,d2w_dy2,d2w_dz2,d2tp_dx2,d2tp_dy2,d2tp_dz2,
     &  d2u_dxy,d2u_dyx,d2u_dxz,d2u_dzx,d2v_dxy,d2v_dyx,d2v_dyz,d2v_dzy,
     &  d2w_dxz,d2w_dzx,d2w_dyz,d2w_dzy,du_0_dy,d2u_0_dy2 

      contains 
ccccc **********************************************************************
ccccc flux_to_primitive: Convert primitive parameters to flux parameters.
ccccc **********************************************************************
      SUBROUTINE flux_to_primitive(local_imax,E1,E2,E3,E4,E5,
     &  rho,u,v,w,Et,p,tp,mu)

      integer i,j,k,local_imax   
      real*8 C2,
     &  E1 (local_imax,jmax,kmax),E2 (local_imax,jmax,kmax),
     &  E3 (local_imax,jmax,kmax),E4 (local_imax,jmax,kmax),
     &  E5 (local_imax,jmax,kmax),rho(local_imax,jmax,kmax),
     &  u  (local_imax,jmax,kmax),v  (local_imax,jmax,kmax),
     &  w  (local_imax,jmax,kmax),Et (local_imax,jmax,kmax),
     &  p  (local_imax,jmax,kmax),tp (local_imax,jmax,kmax),
     &  mu (local_imax,jmax,kmax),e  (local_imax,jmax,kmax)
      parameter (C2=110.4d0/288.15d0) !T_oo=288.15 for sea level.

      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax 
            rho(i,j,k) = E1(i,j,k)
            u  (i,j,k) = E2(i,j,k)/E1(i,j,k) 
            v  (i,j,k) = E3(i,j,k)/E1(i,j,k)
            w  (i,j,k) = E4(i,j,k)/E1(i,j,k)            
            Et (i,j,k) = E5(i,j,k)

            !MS$IF (tp_proc.EQ.1)
              e (i,j,k) = Et(i,j,k)/rho(i,j,k)
     &          -(u(i,j,k)**2.d0+v(i,j,k)**2.d0+w(i,j,k)**2.d0)/2.d0
              p (i,j,k) = (gamma-1.d0)*rho(i,j,k)*e(i,j,k)
              tp(i,j,k) = gamma*Ma**2.d0*p(i,j,k)/rho(i,j,k)
            !MS$ELSEIF (tp_proc.EQ.2)
              if (rho(i,j,k).LT.0) then
                write(*,*) 
                write(*,*) 'ERROR..: Negative Density in Point (',i,j,k,').'	       
                write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.'	                    
                STOP          
              end if 
              p (i,j,k)  = 1.d0/(gamma*Ma**2.d0)*rho(i,j,k)**gamma
              tp(i,j,k) = gamma*Ma**2.d0*p(i,j,k)/rho(i,j,k)
              e (i,j,k)  = p(i,j,k)/(rho(i,j,k)*(gamma-1.d0))
              Et(i,j,k) = rho(i,j,k)*(e(i,j,k)
     &          +(u(i,j,k)**2.d0+v(i,j,k)**2.d0+w(i,j,k)**2.d0)/2.d0) 
            !MS$ENDIF

            !MS$IF (tp_viscous.EQ.1)    
              mu(i,j,k) = tp(i,j,k)**(3.d0/2.d0)*((1.d0+C2)/(tp(i,j,k)+C2))            
            !MS$ELSEIF (tp_viscous.EQ.2)
              mu(i,j,k) = 1.d0
            !MS$ENDIF 
          end do
        end do
      end do
      
      RETURN
      END SUBROUTINE flux_to_primitive

ccccc **********************************************************************
ccccc RK4: The 4th order temporary integration method for conservative 
ccccc      form of Navier-Stokes equations.
ccccc **********************************************************************
      SUBROUTINE RK4(local_imax,x,y,z,tt)
          
      USE sninitial, only: E1,E2,E3,E4,E5
      USE snfilter,  only: filter_x,filter_y,filter_z

      integer i,j,k,tt,local_imax
      real*8 dt2,dt6,x(local_imax),y(jmax),z(kmax),t

      dt2=dt*0.5d0
      dt6=dt/6.d0			      

      t = dble(tt)*dt
      call derivs(local_imax,x,y,z,t,E1,E2,E3,E4,E5,
     &  E1_RK1,E2_RK1,E3_RK1,E4_RK1,E5_RK1)  
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            E1_t(i,j,k) = E1(i,j,k) + dt2*E1_rk1(i,j,k) 
            E2_t(i,j,k) = E2(i,j,k) + dt2*E2_rk1(i,j,k) 
            E3_t(i,j,k) = E3(i,j,k) + dt2*E3_rk1(i,j,k) 
            E4_t(i,j,k) = E4(i,j,k) + dt2*E4_rk1(i,j,k) 
            E5_t(i,j,k) = E5(i,j,k) + dt2*E5_rk1(i,j,k) 
          end do
        end do        
      end do
cc    call filter_x(local_imax,E1_t,E2_t,E3_t,E4_t,E5_t)

      t = dble(tt)*(dt+dt2)
      call derivs(local_imax,x,y,z,t,E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_RK2,E2_RK2,E3_RK2,E4_RK2,E5_RK2)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            E1_t(i,j,k) = E1(i,j,k) + dt2*E1_rk2(i,j,k) 
            E2_t(i,j,k) = E2(i,j,k) + dt2*E2_rk2(i,j,k) 
            E3_t(i,j,k) = E3(i,j,k) + dt2*E3_rk2(i,j,k)
            E4_t(i,j,k) = E4(i,j,k) + dt2*E4_rk2(i,j,k)
            E5_t(i,j,k) = E5(i,j,k) + dt2*E5_rk2(i,j,k)
          end do
        end do
      end do
cc    call filter_x(local_imax,E1_t,E2_t,E3_t,E4_t,E5_t)

      t = dble(tt)*(dt+dt2)
      call derivs(local_imax,x,y,z,t,E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_RK3,E2_RK3,E3_RK3,E4_RK3,E5_RK3)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            E1_t(i,j,k) = E1(i,j,k) + dt*E1_rk3(i,j,k)
            E2_t(i,j,k) = E2(i,j,k) + dt*E2_rk3(i,j,k)
            E3_t(i,j,k) = E3(i,j,k) + dt*E3_rk3(i,j,k)
            E4_t(i,j,k) = E4(i,j,k) + dt*E4_rk3(i,j,k)
            E5_t(i,j,k) = E5(i,j,k) + dt*E5_rk3(i,j,k)
		
            E1_rk3(i,j,k) = E1_rk2(i,j,k) + E1_rk3(i,j,k)
            E2_rk3(i,j,k) = E2_rk2(i,j,k) + E2_rk3(i,j,k)
            E3_rk3(i,j,k) = E3_rk2(i,j,k) + E3_rk3(i,j,k)			  
            E4_rk3(i,j,k) = E4_rk2(i,j,k) + E4_rk3(i,j,k)			  
            E5_rk3(i,j,k) = E5_rk2(i,j,k) + E5_rk3(i,j,k)			  
          end do
        end do
      end do
cc    call filter_x(local_imax,E1_t,E2_t,E3_t,E4_t,E5_t)

      t = dble(tt)*(2.d0*dt)
      call derivs(local_imax,x,y,z,t,E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_RK4,E2_RK4,E3_RK4,E4_RK4,E5_RK4)
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            E1(i,j,k)=E1(i,j,k)+dt6*(E1_rk1(i,j,k)+
     &        2.d0*E1_rk3(i,j,k)+E1_rk4(i,j,k))	    
            E2(i,j,k)=E2(i,j,k)+dt6*(E2_rk1(i,j,k)+
     &        2.d0*E2_rk3(i,j,k)+E2_rk4(i,j,k))	    
            E3(i,j,k)=E3(i,j,k)+dt6*(E3_rk1(i,j,k)+
     &        2.d0*E3_rk3(i,j,k)+E3_rk4(i,j,k))        			
            E4(i,j,k)=E4(i,j,k)+dt6*(E4_rk1(i,j,k)+
     &        2.d0*E4_rk3(i,j,k)+E4_rk4(i,j,k))        			
            E5(i,j,k)=E5(i,j,k)+dt6*(E5_rk1(i,j,k)+
     &        2.d0*E5_rk3(i,j,k)+E5_rk4(i,j,k))        			
          end do
        end do
      end do

      call filter_x(local_imax,E1,E2,E3,E4,E5)
      call filter_y(local_imax,E1,E2,E3,E4,E5)
      call filter_z(local_imax,E1,E2,E3,E4,E5)

      RETURN
      END SUBROUTINE RK4

ccccc **************************************************************************
ccccc derivs: Calculate the spatial derivatives of conservation equation. 
ccccc **************************************************************************
      SUBROUTINE derivs(local_imax,x,y,z,t,E1,E2,E3,E4,E5,
     &  E1_RK,E2_RK,E3_RK,E4_RK,E5_RK)

      USE snmpi,    only: my_rank,pro,exchange_ghostpoints
      USE sninitial,only: rho,u,v,w,Et,p,tp,mu
      USE snbound,  only: boundflow
      USE snmethod, only: stx_ex_6th_D,stxx_ex_6th_D
      
      integer i,j,k,local_imax
      real*8 lambda,t,x(local_imax),y(jmax),z(kmax),C2,
     &  E1     (local_imax,jmax,kmax),E2   (local_imax,jmax,kmax),
     &  E3     (local_imax,jmax,kmax),E4   (local_imax,jmax,kmax),
     &  E5     (local_imax,jmax,kmax),E1_RK(local_imax,jmax,kmax),
     &  E2_RK  (local_imax,jmax,kmax),E3_RK(local_imax,jmax,kmax),
     &  E4_RK  (local_imax,jmax,kmax),E5_RK(local_imax,jmax,kmax),
     &  dmu_dtp(local_imax,jmax,kmax)
      parameter (C2=110.4d0/288.15d0) !T_oo=288.15 for sea level. 

      if (pro.NE.1) then
        call exchange_ghostpoints(jmax,kmax,1,E1)
        call exchange_ghostpoints(jmax,kmax,3,E2)
        call exchange_ghostpoints(jmax,kmax,5,E3)
        call exchange_ghostpoints(jmax,kmax,7,E4)
        call exchange_ghostpoints(jmax,kmax,9,E5)
      end if

      call flux_to_primitive(local_imax,E1,E2,E3,E4,E5,rho,u,v,w,Et,p,tp,mu)
      call boundflow(my_rank,pro,local_imax,x,y,z,t,
     &  stx_ex_6th_D,stxx_ex_6th_D,rho,u,v,w,Et,p,tp,mu)

ccccc ----------------------------------
ccccc CALCULATING THE CONVECTIVE TERM
ccccc ----------------------------------      
      call der_x_D(my_rank,pro,local_imax,stx_ex_6th_D,u,  du_dx  )
      call der_x_D(my_rank,pro,local_imax,stx_ex_6th_D,v,  dv_dx  )    
      call der_x_D(my_rank,pro,local_imax,stx_ex_6th_D,w,  dw_dx  ) 
      call der_x_D(my_rank,pro,local_imax,stx_ex_6th_D,rho,drho_dx)
      call der_x_D(my_rank,pro,local_imax,stx_ex_6th_D,p,  dp_dx  )
      call der_x_D(my_rank,pro,local_imax,stx_ex_6th_D,Et, dEt_dx )
      call der_x_D(my_rank,pro,local_imax,stx_ex_6th_D,tp, dtp_dx )

      !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)                  
        call der_y_WN(local_imax,u,  du_dy  )
        call der_y_WD(local_imax,v,  dv_dy  )
        call der_y_WN(local_imax,w,  dw_dy  )
        call der_y_WN(local_imax,rho,drho_dy) 
        call der_y_WN(local_imax,p,  dp_dy  )        
        call der_y_WN(local_imax,Et, dEt_dy )        
        call der_y_NN(local_imax,tp, dtp_dy )
      !MS$ELSEIF (tp_s.EQ.6)
        call der_y_WN(local_imax,u,  du_dy  )
        call der_y_WN(local_imax,v,  dv_dy  )
        call der_y_WN(local_imax,w,  dw_dy  )
        call der_y_WD(local_imax,rho,drho_dy)
        call der_y_WN(local_imax,p,  dp_dy  )        
        call der_y_WD(local_imax,Et, dEt_dy )        
        call der_y_NN(local_imax,tp, dtp_dy )
      !MS$ENDIF 

      call der_z_P(local_imax,u,  du_dz  )    
      call der_z_P(local_imax,v,  dv_dz  )
      call der_z_P(local_imax,w,  dw_dz  )    
      call der_z_P(local_imax,rho,drho_dz)
      call der_z_P(local_imax,p,  dp_dz  )
      call der_z_P(local_imax,Et, dEt_dz )
      call der_z_P(local_imax,tp, dtp_dz )

ccccc ----------------------------------
ccccc CALCULATING THE VISCOUS TERM
ccccc ----------------------------------
      call der_xx_D(local_imax,u, d2u_dx2 )
      call der_xx_D(local_imax,v, d2v_dx2 )
      call der_xx_D(local_imax,w, d2w_dx2 )
      call der_xx_D(local_imax,tp,d2tp_dx2)

      !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        call der_yy_WN(local_imax,u, du_dy, d2u_dy2 )
        call der_yy_WD(local_imax,v, dv_dy, d2v_dy2 )
        call der_yy_WN(local_imax,w, dw_dy, d2w_dy2 )
        call der_yy_NN(local_imax,tp,dtp_dy,d2tp_dy2)
      !MS$ELSEIF (tp_s.EQ.6)
        call der_yy_WN(local_imax,u, du_dy, d2u_dy2 )
        call der_yy_WN(local_imax,v, dv_dy, d2v_dy2 )
        call der_yy_WN(local_imax,w, dw_dy, d2w_dy2 )
        call der_yy_NN(local_imax,tp,dtp_dy,d2tp_dy2)
      !MS$ENDIF
 
      call der_zz_P(local_imax,u, d2u_dz2 )
      call der_zz_P(local_imax,v, d2v_dz2 )
      call der_zz_P(local_imax,w, d2w_dz2 )
      call der_zz_P(local_imax,tp,d2tp_dz2)

ccccc ----------------------------------
ccccc CALCULATING THE CROSS DERIVATIVES
ccccc ----------------------------------
      call der_xy_DD(local_imax,u,d2u_dxy)
      call der_xy_DD(local_imax,v,d2v_dxy)
      call der_yx_ND(local_imax,u,d2u_dyx)
      call der_yx_DD(local_imax,v,d2v_dyx) 
      call der_xz_DP(local_imax,u,d2u_dxz)
      call der_xz_DP(local_imax,w,d2w_dxz)
      call der_yz_DP(local_imax,v,d2v_dyz) 
      call der_yz_NP(local_imax,w,d2w_dyz) 
      call der_zx_PD(local_imax,u,d2u_dzx) 
      call der_zx_PD(local_imax,w,d2w_dzx) 
      call der_zy_PD(local_imax,v,d2v_dzy)
      call der_zy_PD(local_imax,w,d2w_dzy)

ccccc ---------------------------------
ccccc CALCULATING NON-ALONGAMENT
ccccc ---------------------------------
      !For shear layer problem the viscous 
      !alongment in normal direction was canceled.
      !MS$IF (cancel_alongvis.EQ.0)
        d2u_0_dy2(:,:,:) = 0.d0
      !MS$ELSEIF (cancel_alongvis.EQ.1)
        call der_y_WN (local_imax,u_0,du_0_dy)
        call der_yy_WN(local_imax,u_0,du_0_dy,d2u_0_dy2)
      !MS$ENDIF

      lambda = -1.d0/((gamma-1.d0)*Ma**2.d0*Pr*Re)

      do k=1,kmax 
        do j=1,jmax 
          do i=1,local_imax

            !MS$IF (tp_viscous.EQ.1)    
              dmu_dtp(i,j,k) = (3.d0/2.d0)*dsqrt(tp(i,j,k))*(1.d0+C2)
     &          /(tp(i,j,k)+C2)-tp(i,j,k)**(3.d0/2.d0)*(1.d0+C2)
     &          /((tp(i,j,k)+C2)**2.d0)
            !MS$ELSEIF (tp_viscous.EQ.2)
              dmu_dtp(i,j,k) = 0.d0
            !MS$ENDIF 

            !MS$IF (tp_s.EQ.6)
              dp_dy(i,1,k) = 
     &          -2.d0*rho(i,1,k)*v(i,1,k)*dv_dy(i,1,k)
     &          -v(i,1,k)**2.d0*drho_dy(i,1,k)
     &          +mu(i,1,k)/Re*(
     &            +d2u_dyx(i,1,k) 
     &            +d2v_dx2(i,1,k)
     &            +4.d0/3.d0*d2v_dy2(i,1,k)
     &            -2.d0/3.d0*d2u_dxy(i,1,k) 
     &            -2.d0/3.d0*d2w_dzy(i,1,k) 
     &            +d2v_dz2(i,1,k)
     &            +d2w_dyz(i,1,k)) 
     &          +1.d0/Re*(
     &            +dmu_dtp(i,1,k)*dtp_dx(i,1,k)*du_dy(i,1,k)
     &            +dmu_dtp(i,1,k)*dtp_dx(i,1,k)*dv_dx(i,1,k)
     &            +4.d0/3.d0*dmu_dtp(i,1,k)*dtp_dy(i,1,k)*dv_dy(i,1,k)
     &            -2.d0/3.d0*dmu_dtp(i,1,k)*dtp_dy(i,1,k)*du_dx(i,1,k)
     &            -2.d0/3.d0*dmu_dtp(i,1,k)*dtp_dy(i,1,k)*dw_dz(i,1,k)
     &            +dmu_dtp(i,1,k)*dtp_dz(i,1,k)*dv_dz(i,1,k)
     &            +dmu_dtp(i,1,k)*dtp_dz(i,1,k)*dw_dy(i,1,k))
            !MS$ENDIF 

            E1_RK(i,j,k) =
     &        -rho(i,j,k)*du_dx(i,j,k)
     &        -u(i,j,k)*drho_dx(i,j,k)
     &        -rho(i,j,k)*dv_dy(i,j,k)
     &        -v(i,j,k)*drho_dy(i,j,k)
     &        -rho(i,j,k)*dw_dz(i,j,k)
     &        -w(i,j,k)*drho_dz(i,j,k)
            E2_RK(i,j,k) =
     &        -2.d0*rho(i,j,k)*u(i,j,k)*du_dx(i,j,k)
     &        -u(i,j,k)**2.d0*drho_dx(i,j,k)
     &        -dp_dx(i,j,k)
     &        -rho(i,j,k)*u(i,j,k)*dv_dy(i,j,k)
     &        -rho(i,j,k)*v(i,j,k)*du_dy(i,j,k) 
     &        -u(i,j,k)*v(i,j,k)*drho_dy(i,j,k)
     &        -rho(i,j,k)*u(i,j,k)*dw_dz(i,j,k)  
     &        -rho(i,j,k)*w(i,j,k)*du_dz(i,j,k)  
     &        -u(i,j,k)*w(i,j,k)*drho_dz(i,j,k)  
     &        +mu(i,j,k)/Re*(
     &          +4.d0/3.d0*d2u_dx2(i,j,k)
     &          -2.d0/3.d0*d2v_dyx(i,j,k) 
     &          -2.d0/3.d0*d2w_dzx(i,j,k) 
     &          +d2u_dy2(i,j,k) -0.d0*d2u_0_dy2(i,j,k)
     &          +d2v_dxy(i,j,k)      
     &          +d2u_dz2(i,j,k)
     &          +d2w_dxz(i,j,k))      
     &        +1.d0/Re*(
     &          +4.d0/3.d0*dmu_dtp(i,j,k)*dtp_dx(i,j,k)*du_dx(i,j,k)
     &          -2.d0/3.d0*dmu_dtp(i,j,k)*dtp_dx(i,j,k)*dv_dy(i,j,k)
     &          -2.d0/3.d0*dmu_dtp(i,j,k)*dtp_dx(i,j,k)*dw_dz(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dy(i,j,k)*du_dy(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dy(i,j,k)*dv_dx(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dz(i,j,k)*du_dz(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dz(i,j,k)*dw_dx(i,j,k)
     &        )
            E3_RK(i,j,k) = 
     &        -rho(i,j,k)*u(i,j,k)*dv_dx(i,j,k)
     &        -rho(i,j,k)*v(i,j,k)*du_dx(i,j,k)
     &        -u(i,j,k)*v(i,j,k)*drho_dx(i,j,k)
     &        -2.d0*rho(i,j,k)*v(i,j,k)*dv_dy(i,j,k)
     &        -v(i,j,k)**2.d0*drho_dy(i,j,k)
     &        -dp_dy(i,j,k)
     &        -rho(i,j,k)*v(i,j,k)*dw_dz(i,j,k)
     &        -rho(i,j,k)*w(i,j,k)*dv_dz(i,j,k)
     &        -v(i,j,k)*w(i,j,k)*drho_dz(i,j,k)
     &        +mu(i,j,k)/Re*(
     &          +d2u_dyx(i,j,k) 
     &          +d2v_dx2(i,j,k)
     &          +4.d0/3.d0*d2v_dy2(i,j,k)
     &          -2.d0/3.d0*d2u_dxy(i,j,k) 
     &          -2.d0/3.d0*d2w_dzy(i,j,k) 
     &          +d2v_dz2(i,j,k)
     &          +d2w_dyz(i,j,k)) 
     &        +1.d0/Re*(
     &          +dmu_dtp(i,j,k)*dtp_dx(i,j,k)*du_dy(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dx(i,j,k)*dv_dx(i,j,k)
     &          +4.d0/3.d0*dmu_dtp(i,j,k)*dtp_dy(i,j,k)*dv_dy(i,j,k)
     &          -2.d0/3.d0*dmu_dtp(i,j,k)*dtp_dy(i,j,k)*du_dx(i,j,k)
     &          -2.d0/3.d0*dmu_dtp(i,j,k)*dtp_dy(i,j,k)*dw_dz(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dz(i,j,k)*dv_dz(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dz(i,j,k)*dw_dy(i,j,k)
     &        )
            E4_RK(i,j,k) = 
     &        -rho(i,j,k)*u(i,j,k)*dw_dx(i,j,k)
     &        -rho(i,j,k)*w(i,j,k)*du_dx(i,j,k)
     &        -u(i,j,k)*w(i,j,k)*drho_dx(i,j,k)
     &        -rho(i,j,k)*v(i,j,k)*dw_dy(i,j,k)
     &        -rho(i,j,k)*w(i,j,k)*dv_dy(i,j,k)
     &        -v(i,j,k)*w(i,j,k)*drho_dy(i,j,k)
     &        -2.d0*rho(i,j,k)*w(i,j,k)*dw_dz(i,j,k)
     &        -w(i,j,k)**2.d0*drho_dz(i,j,k)
     &        -dp_dz(i,j,k)
     &        +mu(i,j,k)/Re*(
     &          +d2w_dx2(i,j,k)
     &          +d2u_dzx(i,j,k)            
     &          +d2v_dzy(i,j,k)
     &          +d2w_dy2(i,j,k)
     &          +4.d0/3.d0*d2w_dz2(i,j,k)
     &          -2.d0/3.d0*d2u_dxz(i,j,k) 
     &          -2.d0/3.d0*d2v_dyz(i,j,k)) 
     &        +1.d0/Re*(
     &          +dmu_dtp(i,j,k)*dtp_dx(i,j,k)*dw_dx(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dx(i,j,k)*du_dz(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dy(i,j,k)*dv_dz(i,j,k)
     &          +dmu_dtp(i,j,k)*dtp_dy(i,j,k)*dw_dy(i,j,k)
     &          +4.d0/3.d0*dmu_dtp(i,j,k)*dtp_dz(i,j,k)*dw_dz(i,j,k)
     &          -2.d0/3.d0*dmu_dtp(i,j,k)*dtp_dz(i,j,k)*du_dx(i,j,k)
     &          -2.d0/3.d0*dmu_dtp(i,j,k)*dtp_dz(i,j,k)*dv_dy(i,j,k)
     &        )
            !MS$IF (tp_proc.EQ.1)
              E5_RK(i,j,k) =
     &          -(Et(i,j,k)+p(i,j,k))*
     &          (du_dx(i,j,k)+dv_dy(i,j,k)+dw_dz(i,j,k))
     &          -u(i,j,k)*(dEt_dx(i,j,k)+dp_dx(i,j,k))
     &          -v(i,j,k)*(dEt_dy(i,j,k)+dp_dy(i,j,k))
     &          -w(i,j,k)*(dEt_dz(i,j,k)+dp_dz(i,j,k))
     &          -mu(i,j,k)*lambda*d2tp_dx2(i,j,k)
     &          -mu(i,j,k)*lambda*d2tp_dy2(i,j,k)
     &          -mu(i,j,k)*lambda*d2tp_dz2(i,j,k)
     &          -lambda*dmu_dtp(i,j,k)*dtp_dx(i,j,k)*dtp_dx(i,j,k)
     &          -lambda*dmu_dtp(i,j,k)*dtp_dy(i,j,k)*dtp_dy(i,j,k)
     &          -lambda*dmu_dtp(i,j,k)*dtp_dz(i,j,k)*dtp_dz(i,j,k)

     &          +u(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*(
     &            +2.d0*d2u_dx2(i,j,k)-d2v_dyx(i,j,k)-d2w_dzx(i,j,k))
     &          +u(i,j,k)*2.d0/3.d0*1.d0/Re*dmu_dtp(i,j,k)*dtp_dx(i,j,k)*(
     &            +2.d0*du_dx(i,j,k)-dv_dy(i,j,k)-dw_dz(i,j,k))
     &          +2.d0/3.d0*mu(i,j,k)/Re*du_dx(i,j,k)*(
     &            +2.d0*du_dx(i,j,k)-dv_dy(i,j,k)-dw_dz(i,j,k))
     &          +v(i,j,k)*mu(i,j,k)/Re*(d2u_dyx(i,j,k)+d2v_dx2(i,j,k))
     &          +v(i,j,k)/Re*dmu_dtp(i,j,k)*dtp_dx(i,j,k)*(du_dy(i,j,k)+dv_dx(i,j,k))
     &          +mu(i,j,k)/Re*dv_dx(i,j,k)*(du_dy(i,j,k)+dv_dx(i,j,k))
     &          +w(i,j,k)*mu(i,j,k)/Re*(d2w_dx2(i,j,k)+d2u_dzx(i,j,k))
     &          +w(i,j,k)/Re*dmu_dtp(i,j,k)*dtp_dx(i,j,k)*(dw_dx(i,j,k)+du_dz(i,j,k))
     &          +mu(i,j,k)/Re*dw_dx(i,j,k)*(dw_dx(i,j,k)+du_dz(i,j,k))

     &          +u(i,j,k)*mu(i,j,k)/Re*(d2u_dy2(i,j,k)+d2v_dxy(i,j,k))
     &          +u(i,j,k)/Re*dmu_dtp(i,j,k)*dtp_dy(i,j,k)*(du_dy(i,j,k)+dv_dx(i,j,k))
     &          +mu(i,j,k)/Re*du_dy(i,j,k)*(du_dy(i,j,k)+dv_dx(i,j,k))
     &          +v(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*(
     &            +2.d0*d2v_dy2(i,j,k)-d2u_dxy(i,j,k)-d2w_dzy(i,j,k))
     &          +v(i,j,k)*2.d0/3.d0*1.d0/Re*dmu_dtp(i,j,k)*dtp_dy(i,j,k)*(
     &            +2.d0*dv_dy(i,j,k)-du_dx(i,j,k)-dw_dz(i,j,k))
     &          +2.d0/3.d0*mu(i,j,k)/Re*dv_dy(i,j,k)*(
     &            +2.d0*dv_dy(i,j,k)-du_dx(i,j,k)-dw_dz(i,j,k))
     &          +w(i,j,k)*mu(i,j,k)/Re*(d2v_dzy(i,j,k)+d2w_dy2(i,j,k))
     &          +w(i,j,k)/Re*dmu_dtp(i,j,k)*dtp_dy(i,j,k)*(dv_dz(i,j,k)+dw_dy(i,j,k))
     &          +mu(i,j,k)/Re*dw_dy(i,j,k)*(dv_dz(i,j,k)+dw_dy(i,j,k))

     &          +u(i,j,k)*mu(i,j,k)/Re*(d2w_dxz(i,j,k)+d2u_dz2(i,j,k))
     &          +u(i,j,k)/Re*dmu_dtp(i,j,k)*dtp_dz(i,j,k)*(dw_dx(i,j,k)+du_dz(i,j,k))
     &          +mu(i,j,k)/Re*du_dz(i,j,k)*(dw_dx(i,j,k)+du_dz(i,j,k))
     &          +v(i,j,k)*mu(i,j,k)/Re*(d2v_dz2(i,j,k)+d2w_dyz(i,j,k))
     &          +v(i,j,k)/Re*dmu_dtp(i,j,k)*dtp_dz(i,j,k)*(dv_dz(i,j,k)+dw_dy(i,j,k))
     &          +mu(i,j,k)/Re*dv_dz(i,j,k)*(dv_dz(i,j,k)+dw_dy(i,j,k))
     &          +w(i,j,k)*2.d0/3.d0*mu(i,j,k)/Re*(
     &            +2.d0*d2w_dz2(i,j,k)-d2u_dxz(i,j,k)-d2v_dyz(i,j,k))
     &          +w(i,j,k)*2.d0/3.d0*1.d0/Re*dmu_dtp(i,j,k)*dtp_dz(i,j,k)*(
     &            +2.d0*dw_dz(i,j,k)-du_dx(i,j,k)-dv_dy(i,j,k))
     &          +2.d0/3.d0*mu(i,j,k)/Re*dw_dz(i,j,k)*(
     &            +2.d0*dw_dz(i,j,k)-du_dx(i,j,k)-dv_dy(i,j,k))
          !MS$ELSEIF (tp_proc.EQ.2)
            E5_RK(i,j,k)=0.d0
          !MS$ENDIF          
          end do
        end do
      end do	

      RETURN
      END SUBROUTINE derivs

ccccc **********************************************************************
ccccc der_x_D: This routine calculated the first derivatives in the
ccccc          x direction. Periodic and Dirichlet boundary condition 
ccccc          was adopted in according to the problem simulated.
ccccc **********************************************************************
      SUBROUTINE der_x_D(my_rank,pro,local_imax,stx,fc,dfc_dx)

      USE snmethod, only: der_x_ex_6th_D

      integer local_imax,my_rank,pro
      real*8 stx(local_imax,7),fc(local_imax,jmax,kmax),
     &  dfc_dx(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2).OR.(tp_s.EQ.3)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,fc,dfc_dx)
      !MS$ELSEIF (tp_s.EQ.4)        
        call der_x_cp_6th_D(dx,local_imax,jmax,kmax,fc,dfc_dx)
      !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx,fc,dfc_dx)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_x_D

ccccc **********************************************************************
ccccc der_y_WD: This routine calculated the first derivatives in the  
ccccc           y direction. 
ccccc **********************************************************************
      SUBROUTINE der_y_WD(local_imax,fc,dfc_dy)

      USE snmethod, only: der_y_cp_62th_D,der_y_cp_65th_WD,dcy_dpy

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
      !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_D(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_D(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF         
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)
          call der_y_cp_65th_WD(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WD(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF         
      !MS$ENDIF  

      RETURN
      END SUBROUTINE der_y_WD

ccccc **********************************************************************
ccccc der_y_WN: This routine calculated the first derivatives in the 
ccccc           y-direction. 
ccccc **********************************************************************
      SUBROUTINE der_y_WN(local_imax,fc,dfc_dy)

      USE snmethod, only: der_y_cp_62th_N,der_y_cp_65th_WN,dcy_dpy

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
      !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_N(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_N(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF 
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_65th_WN(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WN(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF 
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_y_WN

ccccc **********************************************************************
ccccc der_y_NN: This routine calculated the first derivatives in the  
ccccc           y-direction. Neumann boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_y_NN(local_imax,fc,dfc_dy)

      USE snmethod, only: der_y_cp_62th_N,der_y_cp_65th_NN,dcy_dpy

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
      !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_N(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_N(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF         
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_65th_NN(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_NN(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF         
      !MS$ENDIF  

      RETURN
      END SUBROUTINE der_y_NN

ccccc **********************************************************************
ccccc der_z_P: This routine calculated the first derivatives in the
ccccc          z direction. Periodic boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_z_P(local_imax,fc,dfc_dz)

      USE snmethod, only: der_z_cp_6th_P

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dz(local_imax,jmax,kmax)

      call der_z_cp_6th_P(dz,local_imax,jmax,kmax,fc,dfc_dz)

      RETURN
      END SUBROUTINE der_z_P

ccccc **********************************************************************
ccccc der_xx_D: This routine calculated the second derivatives in the
ccccc           x direction. Periodic and Neumann boundary condition 
ccccc           was adopted in according to the problem simulated.
ccccc **********************************************************************
      SUBROUTINE der_xx_D(local_imax,fc,d2fc_dx2)

      USE snmpi,    only: my_rank,pro
      USE snmethod, only: der_xx_ex_6th_D,stxx_ex_6th_D

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),d2fc_dx2(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2).OR.(tp_s.EQ.3)
        call der_xx_cp_6th_P(dx,local_imax,jmax,kmax,fc,d2fc_dx2)
      !MS$ELSEIF (tp_s.EQ.4)
        call der_xx_cp_6th_D(dx,local_imax,jmax,kmax,fc,d2fc_dx2)
      !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)  
        call der_xx_ex_6th_D(local_imax,jmax,kmax,stxx_ex_6th_D,fc,d2fc_dx2)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_xx_D

ccccc **********************************************************************
ccccc der_yy_WD: This routine calculated the second derivatives in the
ccccc            y-direction. Dirichlet boundary condition was used.
ccccc **********************************************************************
      SUBROUTINE der_yy_WD(local_imax,fc,dfc_dy,d2fc_dy2)

      USE snmethod, only: der_yy_cp_62th_D,der_yy_cp_65th_WD,
     &  dcy_dpy,d2cy_dpy2

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax),
     &  d2fc_dy2(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_yy_cp_6th_P(dy,local_imax,jmax,kmax,fc,d2fc_dy2)
      !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_yy_cp_62th_D(dy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_yy_cp_62th_D(dcy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ENDIF 
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)
          call der_yy_cp_65th_WD(dy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ELSEIF (stretching_in_y.EQ.1)
          call der_yy_cp_65th_WD(dcy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ENDIF
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_yy_WD

ccccc **********************************************************************
ccccc der_yy_WN: This routine calculated the second derivatives in the
ccccc            y-direction. Neumann boundary condition was used.
ccccc **********************************************************************
      SUBROUTINE der_yy_WN(local_imax,fc,dfc_dy,d2fc_dy2)

      USE snmethod, only: der_yy_cp_62th_N,der_yy_cp_65th_WN,
     &  dcy_dpy,d2cy_dpy2

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax),
     &  d2fc_dy2(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_yy_cp_6th_P(dy,local_imax,jmax,kmax,fc,d2fc_dy2)
      !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0) 
          call der_yy_cp_62th_N(dy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ELSEIF (stretching_in_y.EQ.1)     
          call der_yy_cp_62th_N(dcy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ENDIF  
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)
          call der_yy_cp_65th_WN(dy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ELSEIF (stretching_in_y.EQ.1)
          call der_yy_cp_65th_WN(dcy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ENDIF
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_yy_WN

ccccc **********************************************************************
ccccc der_yy_NN: This routine calculated the second derivatives in the
ccccc            y-direction. Neumann boundary condition was used.
ccccc **********************************************************************
      SUBROUTINE der_yy_NN(local_imax,fc,dfc_dy,d2fc_dy2)

      USE snmethod, only: der_yy_cp_62th_N,der_yy_cp_65th_NN,
     &  dcy_dpy,d2cy_dpy2

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax),
     &  d2fc_dy2(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_yy_cp_6th_P(dy,local_imax,jmax,kmax,fc,d2fc_dy2)
      !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0) 
          call der_yy_cp_62th_N(dy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ELSEIF (stretching_in_y.EQ.1)     
          call der_yy_cp_62th_N(dcy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ENDIF  
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0) 
          call der_yy_cp_65th_NN(dy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ELSEIF (stretching_in_y.EQ.1)     
          call der_yy_cp_65th_NN(dcy,local_imax,jmax,kmax,
     &      dcy_dpy,d2cy_dpy2,fc,dfc_dy,d2fc_dy2)
        !MS$ENDIF  
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_yy_NN

ccccc **********************************************************************
ccccc der_zz_P: This routine calculated the second derivatives in the
ccccc           z direction. Periodic boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_zz_P(local_imax,fc,d2fc_dz2)

      USE snmethod, only: der_zz_cp_6th_P

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),d2fc_dz2(local_imax,jmax,kmax)

      call der_zz_cp_6th_P(dz,local_imax,jmax,kmax,fc,d2fc_dz2)

      RETURN
      END SUBROUTINE der_zz_P

ccccc **********************************************************************
ccccc der_xy_DD: This routine calculated the second derivatives in the  
ccccc            x and y direction. Periodic boundary condition was 
ccccc            used for first derivatives and Dirichlet boundary    
ccccc            condition was used for the second derivatives.
ccccc **********************************************************************
      SUBROUTINE der_xy_DD(local_imax,fc,d2fc_dxy)

      USE snmpi,    only: my_rank,pro
      USE snmethod, only: der_x_ex_6th_D,der_y_cp_62th_D,
     &  der_y_cp_65th_WN,dcy_dpy,stx_ex_6th_D

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dx(local_imax,jmax,kmax),
     &  d2fc_dxy(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,fc,dfc_dx)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,dfc_dx,d2fc_dxy)
      !MS$ELSEIF (tp_s.EQ.3)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,fc,dfc_dx)
        call der_y_cp_62th_D(dcy_dpy,dfc_dx,d2fc_dxy)
      !MS$ELSEIF (tp_s.EQ.4)
        call der_x_cp_6th_D (dx,local_imax,jmax,kmax,fc,dfc_dx)
        call der_y_cp_62th_D(dcy_dpy,dfc_dx,d2fc_dxy)
      !MS$ELSEIF (tp_s.EQ.5)
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx_ex_6th_D,fc,dfc_dx)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_D(dy,local_imax,jmax,kmax,dcy_dpy,dfc_dx,d2fc_dxy)
        !MS$ELSEIF (stretching_in_y.EQ.1)  
          call der_y_cp_62th_D(dcy,local_imax,jmax,kmax,dcy_dpy,dfc_dx,d2fc_dxy)         
        !MS$ENDIF 
      !MS$ELSEIF (tp_s.EQ.6)  
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx_ex_6th_D,fc,dfc_dx)
        !MS$IF (stretching_in_y.EQ.0)
          call der_y_cp_65th_WN(dy,local_imax,jmax,kmax,dcy_dpy,dfc_dx,d2fc_dxy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WN(dcy,local_imax,jmax,kmax,dcy_dpy,dfc_dx,d2fc_dxy)
        !MS$ENDIF
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_xy_DD

ccccc **********************************************************************
ccccc der_xz_DP: This routine calculated the second derivatives in the 
ccccc            x and z direction. Periodic boundary condition was 
ccccc            used for first and second derivatives.
ccccc **********************************************************************
      SUBROUTINE der_xz_DP(local_imax,fc,d2fc_dxz)

      USE snmpi,    only: my_rank,pro
      USE snmethod, only: der_x_ex_6th_D,der_z_cp_6th_P,stx_ex_6th_D

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dx(local_imax,jmax,kmax),
     &  d2fc_dxz(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2).OR.(tp_s.EQ.3)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,fc,dfc_dx)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dx,d2fc_dxz)
      !MS$ELSEIF (tp_s.EQ.4)
        call der_x_cp_6th_D(dx,local_imax,jmax,kmax,fc,dfc_dx)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dx,d2fc_dxz)
      !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx_ex_6th_D,fc,dfc_dx)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dx,d2fc_dxz)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_xz_DP

ccccc **********************************************************************
ccccc der_zy_PD: This routine calculated the second derivatives in the 
ccccc            z and y direction. Periodic boundary condition was 
ccccc            used for first derivatives and Dirichlet boundary  
ccccc            condition was used for the second derivatives.
ccccc **********************************************************************
      SUBROUTINE der_zy_PD(local_imax,fc,d2fc_dzy)

      USE snmethod, only: der_z_cp_6th_P,der_y_cp_62th_D,
     &  der_y_cp_65th_WD,dcy_dpy

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dz(local_imax,jmax,kmax),
     &  d2fc_dzy(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,fc,dfc_dz)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,dfc_dz,d2fc_dzy)
      !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,fc,dfc_dz)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_D(dy,local_imax,jmax,kmax,dcy_dpy,dfc_dz,d2fc_dzy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_D(dcy,local_imax,jmax,kmax,dcy_dpy,dfc_dz,d2fc_dzy)
        !MS$ENDIF 
      !MS$ELSEIF (tp_s.EQ.6)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,fc,dfc_dz)
        !MS$IF (stretching_in_y.EQ.0)
          call der_y_cp_65th_WD(dy,local_imax,jmax,kmax,dcy_dpy,dfc_dz,d2fc_dzy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WD(dcy,local_imax,jmax,kmax,dcy_dpy,dfc_dz,d2fc_dzy)
        !MS$ENDIF         
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_zy_PD

ccccc **********************************************************************
ccccc der_yx_DD: This routine calculate the second derivatives in the 
ccccc            y and x direction. Dirichlet boundary condition 
ccccc            was used for first derivatives and periodic boundary 
ccccc            condition was used for the second derivatives.
ccccc **********************************************************************
      SUBROUTINE der_yx_DD(local_imax,fc,d2fc_dyx) 

      USE snmpi,    only: my_rank,pro
      USE snmethod, only: der_y_cp_62th_D,der_x_ex_6th_D,
     &  der_y_cp_65th_WN,dcy_dpy,d2cy_dpy2,stx_ex_6th_D

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax),
     &  d2fc_dyx(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)        
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,dfc_dy,d2fc_dyx)
      !MS$ELSEIF (tp_s.EQ.3)
        call der_y_cp_62th_D(dcy_dpy,fc,dfc_dy)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,dfc_dy,d2fc_dyx )
      !MS$ELSEIF (tp_s.EQ.4)
        call der_y_cp_62th_D(dcy_dpy,fc,dfc_dy)
        call der_x_cp_6th_D (dx,local_imax,jmax,kmax,dfc_dy,d2fc_dyx )
      !MS$ELSEIF (tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_D(dy,local_imax,jmax,kmax,
     &      dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_D(dcy,local_imax,jmax,kmax,
     &      dcy_dpy,fc,dfc_dy)
        !MS$ENDIF 
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx_ex_6th_D,dfc_dy,d2fc_dyx)
      !MS$ELSEIF (tp_s.EQ.6)  
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_65th_WN(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WN(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF 
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx_ex_6th_D,dfc_dy,d2fc_dyx)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_yx_DD

ccccc **********************************************************************
ccccc der_yx_ND: This routine calculate the second derivatives in the 
ccccc            y and x direction. Neumann boundary condition was 
ccccc            used for first derivatives and periodic boundary 
ccccc            condition was used for the second derivatives.
ccccc **********************************************************************
      SUBROUTINE der_yx_ND(local_imax,fc,d2fc_dyx) 

      USE snmpi,    only: my_rank,pro
      USE snmethod, only: der_y_cp_62th_N,der_x_ex_6th_D,
     &  der_y_cp_65th_WN,dcy_dpy,d2cy_dpy2,stx_ex_6th_D

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax),
     &  d2fc_dyx(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,dfc_dy,d2fc_dyx)
      !MS$ELSEIF (tp_s.EQ.3)
        call der_y_cp_62th_N(dcy_dpy,fc,dfc_dy)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,dfc_dy,d2fc_dyx)
      !MS$ELSEIF (tp_s.EQ.4)
        call der_y_cp_62th_N(dcy_dpy,fc,dfc_dy)
        call der_x_cp_6th_D (dx,local_imax,jmax,kmax,dfc_dy,d2fc_dyx)
      !MS$ELSEIF (tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_N(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_N(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx_ex_6th_D,dfc_dy,d2fc_dyx)
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_65th_WN(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WN(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF 
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx_ex_6th_D,dfc_dy,d2fc_dyx)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_yx_ND

ccccc **********************************************************************
ccccc der_zx_PD: This routine calculate the second derivatives in the 
ccccc            z and x direction. Periodic boundary condition was 
ccccc            used for first and second derivatives.
ccccc **********************************************************************
      SUBROUTINE der_zx_PD(local_imax,fc,d2fc_dzx) 

      USE snmpi,    only: my_rank,pro
      USE snmethod, only: der_z_cp_6th_P,der_x_ex_6th_D,stx_ex_6th_D

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dz(local_imax,jmax,kmax),
     &  d2fc_dzx(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2).OR.(tp_s.EQ.3)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,fc,dfc_dz)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,dfc_dz,d2fc_dzx)
      !MS$ELSEIF (tp_s.EQ.4)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,fc,dfc_dz)
        call der_x_cp_6th_D(dx,local_imax,jmax,kmax,dfc_dz,d2fc_dzx)
      !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,fc,dfc_dz)
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,stx_ex_6th_D,dfc_dz,d2fc_dzx)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_zx_PD

ccccc **********************************************************************
ccccc der_yz_DP: This routine calculate the second derivatives in the 
ccccc            y and z direction. Dirichlet boundary condition 
ccccc            was used for first derivatives and periodic boundary 
ccccc            condition was used for the second derivatives.
ccccc **********************************************************************
      SUBROUTINE der_yz_DP(local_imax,fc,d2fc_dyz) 

      USE snmethod, only: der_y_cp_62th_D,der_z_cp_6th_P,
     &  der_y_cp_65th_WN,dcy_dpy

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax),
     &  d2fc_dyz(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dy,d2fc_dyz)
      !MS$ELSEIF (tp_s.EQ.3)
        call der_y_cp_62th_D(dcy_dpy,fc,dfc_dy)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dy,d2fc_dyz)
      !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_D(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_D(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF 
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dy,d2fc_dyz)
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_65th_WN(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WN(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF 
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dy,d2fc_dyz)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_yz_DP

ccccc **********************************************************************
ccccc der_yz_NP: This routine calculate the second derivatives in the 
ccccc            y and z direction. Neumann boundary condition was 
ccccc            used for first derivatives and periodic boundary 
ccccc            condition was used for the second derivatives.
ccccc **********************************************************************
      SUBROUTINE der_yz_NP(local_imax,fc,d2fc_dyz)

      USE snmethod, only: der_y_cp_62th_N,der_z_cp_6th_P,
     &  der_y_cp_65th_WN,dcy_dpy

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax),
     &  d2fc_dyz(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dy,d2fc_dyz)
      !MS$ELSEIF (tp_s.EQ.3)
        call der_y_cp_62th_N(dcy_dpy,fc,dfc_dy)
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dy,d2fc_dyz)
      !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_N(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_N(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF         
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dy,d2fc_dyz)
      !MS$ELSEIF (tp_s.EQ.6)      
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_65th_WN(dy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WN(dcy,local_imax,jmax,kmax,dcy_dpy,fc,dfc_dy)
        !MS$ENDIF 
        call der_z_cp_6th_P(dz,local_imax,jmax,kmax,dfc_dy,d2fc_dyz)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_yz_NP

      END MODULE snconf 
