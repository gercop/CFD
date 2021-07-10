      MODULE snconf

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

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

      allocate(E1_t  (local_imax,jmax,kmax),E2_t  (local_imax,jmax,kmax),
     &         E3_t  (local_imax,jmax,kmax),E4_t  (local_imax,jmax,kmax),
     &         E5_t  (local_imax,jmax,kmax),E1_rk1(local_imax,jmax,kmax),
     &         E1_rk2(local_imax,jmax,kmax),E1_rk3(local_imax,jmax,kmax),
     &         E1_rk4(local_imax,jmax,kmax),E2_rk1(local_imax,jmax,kmax),
     &         E2_rk2(local_imax,jmax,kmax),E2_rk3(local_imax,jmax,kmax),
     &         E2_rk4(local_imax,jmax,kmax),E3_rk1(local_imax,jmax,kmax),
     &         E3_rk2(local_imax,jmax,kmax),E3_rk3(local_imax,jmax,kmax),
     &         E3_rk4(local_imax,jmax,kmax),E4_rk1(local_imax,jmax,kmax),
     &         E4_rk2(local_imax,jmax,kmax),E4_rk3(local_imax,jmax,kmax),
     &         E4_rk4(local_imax,jmax,kmax),E5_rk1(local_imax,jmax,kmax),
     &         E5_rk2(local_imax,jmax,kmax),E5_rk3(local_imax,jmax,kmax),
     &         E5_rk4(local_imax,jmax,kmax))

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

      deallocate(E1_t,E2_t,E3_t,E4_t,E5_t,E1_rk1,E1_rk2,E1_rk3,E1_rk4,
     &           E2_rk1,E2_rk2,E2_rk3,E2_rk4,E3_rk1,E3_rk2,E3_rk3,E3_rk4,
     &           E4_rk1,E4_rk2,E4_rk3,E4_rk4,E5_rk1,E5_rk2,E5_rk3,E5_rk4)
 
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
      USE snmethod, only: stx_ex_6th_D,stxx_ex_6th_D,
     &  der_x_D,der_y_NN,der_y_WN,der_y_WD,der_z_P,der_xx_D,der_yy_NN,
     &  der_yy_WN,der_yy_WD,der_zz_P,der_xy_DD,der_yx_ND,der_yx_DD,
     &  der_xz_DP,der_yz_DP,der_yz_NP,der_zx_PD,der_zy_PD 
      
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

      allocate(
     &  rho(local_imax,jmax,kmax),u  (local_imax,jmax,kmax),
     &  v  (local_imax,jmax,kmax),w  (local_imax,jmax,kmax),
     &  Et (local_imax,jmax,kmax),p  (local_imax,jmax,kmax),
     &  tp (local_imax,jmax,kmax),mu (local_imax,jmax,kmax))
      allocate(
     &  du_dx  (local_imax,jmax,kmax),du_dy  (local_imax,jmax,kmax),
     &  du_dz  (local_imax,jmax,kmax),dv_dx  (local_imax,jmax,kmax), 
     &  dv_dy  (local_imax,jmax,kmax),dv_dz  (local_imax,jmax,kmax),
     &  dw_dx  (local_imax,jmax,kmax),dw_dy  (local_imax,jmax,kmax),
     &  dw_dz  (local_imax,jmax,kmax),drho_dx(local_imax,jmax,kmax),
     &  drho_dy(local_imax,jmax,kmax),drho_dz(local_imax,jmax,kmax),
     &  dEt_dx (local_imax,jmax,kmax),dEt_dy (local_imax,jmax,kmax),
     &  dEt_dz (local_imax,jmax,kmax),dp_dx  (local_imax,jmax,kmax),
     &  dp_dy  (local_imax,jmax,kmax),dp_dz  (local_imax,jmax,kmax),
     &  dtp_dx (local_imax,jmax,kmax),dtp_dy (local_imax,jmax,kmax),
     &  dtp_dz (local_imax,jmax,kmax))
      allocate(
     &  d2u_dx2  (local_imax,jmax,kmax),d2u_dy2  (local_imax,jmax,kmax),
     &  d2u_dz2  (local_imax,jmax,kmax),d2v_dx2  (local_imax,jmax,kmax), 
     &  d2v_dy2  (local_imax,jmax,kmax),d2v_dz2  (local_imax,jmax,kmax),
     &  d2w_dx2  (local_imax,jmax,kmax),d2w_dy2  (local_imax,jmax,kmax),
     &  d2w_dz2  (local_imax,jmax,kmax),d2tp_dx2 (local_imax,jmax,kmax),
     &  d2tp_dy2 (local_imax,jmax,kmax),d2tp_dz2 (local_imax,jmax,kmax))
      allocate(
     &  d2u_dxy(local_imax,jmax,kmax),d2u_dyx(local_imax,jmax,kmax),
     &  d2u_dxz(local_imax,jmax,kmax),d2u_dzx(local_imax,jmax,kmax), 
     &  d2v_dxy(local_imax,jmax,kmax),d2v_dyx(local_imax,jmax,kmax),
     &  d2v_dyz(local_imax,jmax,kmax),d2v_dzy(local_imax,jmax,kmax),
     &  d2w_dxz(local_imax,jmax,kmax),d2w_dzx(local_imax,jmax,kmax),
     &  d2w_dyz(local_imax,jmax,kmax),d2w_dzy(local_imax,jmax,kmax))
      allocate(
     &  du_0_dy(local_imax,jmax,kmax),d2u_0_dy2(local_imax,jmax,kmax)) 

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

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)                  
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
      call der_xx_D(my_rank,pro,local_imax,u, d2u_dx2 )
      call der_xx_D(my_rank,pro,local_imax,v, d2v_dx2 )
      call der_xx_D(my_rank,pro,local_imax,w, d2w_dx2 )
      call der_xx_D(my_rank,pro,local_imax,tp,d2tp_dx2)

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
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
      call der_xy_DD(my_rank,pro,local_imax,u,d2u_dxy)
      call der_xy_DD(my_rank,pro,local_imax,v,d2v_dxy)
      call der_yx_ND(my_rank,pro,local_imax,u,d2u_dyx)
      call der_yx_DD(my_rank,pro,local_imax,v,d2v_dyx) 
      call der_xz_DP(my_rank,pro,local_imax,u,d2u_dxz)
      call der_xz_DP(my_rank,pro,local_imax,w,d2w_dxz)
      call der_yz_DP(local_imax,v,d2v_dyz) 
      call der_yz_NP(local_imax,w,d2w_dyz) 
      call der_zx_PD(my_rank,pro,local_imax,u,d2u_dzx) 
      call der_zx_PD(my_rank,pro,local_imax,w,d2w_dzx) 
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
     &          +d2u_dy2(i,j,k)-d2u_0_dy2(i,j,k)
     &          +d2v_dxy(i,j,k)      
     &          +d2u_dz2(i,j,k)
     &          +d2w_dxz(i,j,k)
     &         )      
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

      deallocate(rho,u,v,w,Et,p,tp,mu,du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,
     &  dw_dx,dw_dy,dw_dz,drho_dx,drho_dy,drho_dz,dEt_dx,dEt_dy,
     &  dEt_dz,dp_dx,dp_dy,dp_dz,dtp_dx,dtp_dy,dtp_dz,d2u_dx2,
     &  d2u_dy2,d2u_dz2,d2v_dx2,d2v_dy2,d2v_dz2,d2w_dx2,d2w_dy2,
     &  d2w_dz2,d2tp_dx2,d2tp_dy2,d2tp_dz2,d2u_dxy,d2u_dyx,d2u_dxz,
     &  d2u_dzx,d2v_dxy,d2v_dyx,d2v_dyz,d2v_dzy,d2w_dxz,d2w_dzx,
     &  d2w_dyz,d2w_dzy,du_0_dy,d2u_0_dy2)

      RETURN
      END SUBROUTINE derivs

ccccc **********************************************************************
ccccc extrapolated_inlet: This routine extrapolated the variables at 
ccccc                     inflow boundaries.
ccccc **********************************************************************
      SUBROUTINE extrapolated_inlet(local_imax,j,k,stx,fc)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer j,k,local_imax 
      real*8 fc(local_imax,jmax,kmax),stx(local_imax,7)

      !Extrapolation of sixth order
      fc(1,j,k) = 1.d0/stx(1,1)*(
     &  -stx(1,2)*fc(2,j,k)+
     &  -stx(1,3)*fc(3,j,k)+
     &  -stx(1,4)*fc(4,j,k)+
     &  -stx(1,5)*fc(5,j,k)+
     &  -stx(1,6)*fc(6,j,k)+
     &  -stx(1,7)*fc(7,j,k))

      RETURN
      END SUBROUTINE extrapolated_inlet

ccccc **********************************************************************
ccccc extrapolated_outlet: This routine extrapolated the variables at
ccccc                      outflow boundaries.
ccccc **********************************************************************
      SUBROUTINE extrapolated_outlet(local_imax,j,k,stxx,fc)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer j,k,local_imax
      real*8 fc(local_imax,jmax,kmax),stxx(local_imax,7)

      !Extrapolation of fifth order
      fc(local_imax,j,k) = 1.d0/stxx(local_imax,7)*(
     &  -stxx(local_imax,1)*fc(local_imax-6,j,k)+
     &  -stxx(local_imax,2)*fc(local_imax-5,j,k)+
     &  -stxx(local_imax,3)*fc(local_imax-4,j,k)+
     &  -stxx(local_imax,4)*fc(local_imax-3,j,k)+
     &  -stxx(local_imax,5)*fc(local_imax-2,j,k)+
     &  -stxx(local_imax,6)*fc(local_imax-1,j,k))

      RETURN
      END SUBROUTINE extrapolated_outlet

ccccc **********************************************************************
ccccc boundflow: This routine calculated the boundary conditions for 
ccccc            spatial development of shear layer flow.
ccccc **********************************************************************
      SUBROUTINE boundflow(my_rank,pro,local_imax,x,y,z,t,
     &  stx,stxx,rho,u,v,w,Et,p,tp,mu)

      USE sninitial, only: rho_0,u_0,v_0,tp_0

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer i,j,k,local_imax,my_rank,pro
      real*8 t,yy,zz,u_1,v_1,w_1,phi,C2,x(local_imax),
     &  y(jmax),z(kmax),stx(local_imax,7),stxx(local_imax,7),
     &  rho(local_imax,jmax,kmax),u (local_imax,jmax,kmax),
     &  v  (local_imax,jmax,kmax),w (local_imax,jmax,kmax),
     &  Et (local_imax,jmax,kmax),p (local_imax,jmax,kmax),
     &  tp (local_imax,jmax,kmax),mu(local_imax,jmax,kmax),
     &  e  (local_imax,jmax,kmax)
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
           
              u  (1,j,k)  = u_0  (1,j,k) + 1.d0/(M1*c)*u_1 
              v  (1,j,k)  = v_0  (1,j,k) + 1.d0/(M1*c)*v_1
              w  (1,j,k)  =              + 1.d0/(M1*c)*w_1
              rho(1,j,k)  = rho_0(1,j,k)
              tp (1,j,k)  = tp_0 (1,j,k)     

              call extrapolated_inlet(local_imax,j,k,stx,p)     !Pressure is extrapolated to satisfy dp_dx=0

              e  (1,j,k) = p(1,j,k)/(rho(1,j,k)*(gamma-1.d0)) 
              Et (1,j,k) = rho(1,j,k)*(e(1,j,k)+1.d0/2.d0*(
     &          u(1,j,k)**2.d0+v(1,j,k)**2.d0+w(1,j,k)**2.d0))
              !MS$IF (tp_viscous.EQ.1)
                mu(1,j,k) = tp(1,j,k)**(3.d0/2.d0)*((1.d0+C2)
     &            /(tp(1,j,k)+C2))
              !MS$ELSEIF (tp_viscous.EQ.2)
                mu(1,j,k) = 1.d0
              !MS$ENDIF 
            end if         

            if (my_rank+1.EQ.pro) then
              call extrapolated_outlet(local_imax,j,k,stxx,u  ) !Velocity is extrapolated to satisfy d2u_dx2=0
              call extrapolated_outlet(local_imax,j,k,stxx,v  ) !Velocity is extrapolated to satisfy d2v_dx2=0
              call extrapolated_outlet(local_imax,j,k,stxx,w  ) !Velocity is extrapolated to satisfy d2w_dx2=0
              call extrapolated_outlet(local_imax,j,k,stxx,rho) !Density is extrapolated to satisfy d2rho_dx2=0
              call extrapolated_outlet(local_imax,j,k,stxx,tp ) !Temperature is extrapolated to satisfy d2tp_dx2=0
           
              p(local_imax,j,k) = 1.d0/(gamma*Ma**2.d0)         !Pressure is fixed in x(imax)

              e  (local_imax,j,k) = p(local_imax,j,k)/
     &          (rho(local_imax,j,k)*(gamma-1.d0)) 
              Et (local_imax,j,k) = rho(local_imax,j,k)*(
     &          e(local_imax,j,k)+1.d0/2.d0*(u(local_imax,j,k)**2.d0
     &          +v(local_imax,j,k)**2.d0+w(local_imax,j,k)**2.d0)) 
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

      !Free-slip boundary condition - impermeability
      !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        do k=1,kmax
          do i=1,local_imax
            v  (i,1,   k) = 0.d0
            v  (i,jmax,k) = 0.d0            
          end do
        end do
      !MS$ENDIF       

      !Wall-slip boundary condition
      !MS$IF (tp_s.EQ.6)
        do k=1,kmax
          do i=1,local_imax
            u(i,1,k) = 0.d0
            v(i,1,k) = 0.d0
            w(i,1,k) = 0.d0
          end do
        end do
      !MS$ENDIF       

      RETURN
      END SUBROUTINE boundflow

      END MODULE snconf 
