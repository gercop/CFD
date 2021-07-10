      PROGRAM snttcl

      USE snmpi, only: my_rank,local_imax
      USE snmpi, only: initmpi,time_start,decompose,time_end

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      call initmpi
      call time_start
      call decompose(imax)
      call main(my_rank,local_imax)
      call time_end

      END

ccccc **********************************************************************
ccccc main: Main subroutine
ccccc **********************************************************************
      SUBROUTINE main(my_rank,local_imax)
    
      USE sninitial, only: init_cond
      USE snconf,    only: RK4

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer tt,my_rank,local_imax
      real*8 x(local_imax),y(jmax),z(kmax),
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax)

      call ini_iosfile

      call init_data(x,y,z,E1,E2,E3,E4,E5) 
      call init_cond(local_imax,x,y,z,E1,E2,E3,E4,E5)
      call veri_data(local_imax,x,y,z,E1,E2,E3,E4,E5)     
      call view_data(my_rank,local_imax,x,y,z)  

      call sv_all(x,y,z,0,E1,E2,E3,E4,E5)
      if (my_rank.EQ.0) then
        write(*,*)
        write(*,'(A36)') 'INITIAL CONDICIONAL........: tt = 0'
        write(*,*)
      end if

      do tt=1,tmax
        if (my_rank.EQ.0) write(*,'(A29,A1,I8,A1,f17.10,A1,f15.10)')
     &    'NUMERICAL SOLUTION (tt,t) = ','(',tt,',',dble(tt)*dt,')'

        call RK4(local_imax,x,y,z,tt,E1,E2,E3,E4,E5)

        if (((mod(tt,qtimes).EQ.0)).AND.(tt.GT.0)) then
          call sv_all(x,y,z,tt,E1,E2,E3,E4,E5)
        end if
      end do

      call end_data
      call end_iosfile

      RETURN
      END

ccccc **********************************************************************
ccccc init_data: Initialize data.
ccccc **********************************************************************
      SUBROUTINE init_data(x,y,z,E1,E2,E3,E4,E5)

      USE snmpi,     only: my_rank,pro,local_imax,global_x1
      USE snmethod,  only: stx_ex_6th_D,stxx_ex_6th_D,dcy_dpy,d2cy_dpy2
      USE sninitial, only: rho_0,u_0,v_0,tp_0
      USE snconf,    only: E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_rk1,E1_rk2,E1_rk3,E1_rk4,E2_rk1,E2_rk2,E2_rk3,E2_rk4, 
     &  E3_rk1,E3_rk2,E3_rk3,E3_rk4,E4_rk1,E4_rk2,E4_rk3,E4_rk4,
     &  E5_rk1,E5_rk2,E5_rk3,E5_rk4,
     &  du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,dp_dx,
     &  dp_dy,dp_dz,drho_dx,drho_dy,drho_dz,dEt_dx,dEt_dy,dEt_dz,
     &  dtp_dx,dtp_dy,dtp_dz,dmu_dx,dmu_dy,dmu_dz,d2u_dx2,d2u_dy2,
     &  d2u_dz2,d2v_dx2,d2v_dy2,d2v_dz2,d2w_dx2,d2w_dy2,d2w_dz2,
     &  d2tp_dx2,d2tp_dy2,d2tp_dz2,d2u_dxy,d2u_dyx,d2u_dxz,d2u_dzx, 
     &  d2v_dxy,d2v_dyx,d2v_dyz,d2v_dzy,d2w_dxz,d2w_dzx,d2w_dyz,d2w_dzy 

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer i,j,k,imax_spo,local_i
      real*8 x(local_imax),y(jmax),z(kmax),t,Lx_spo,cyy,
     &  py_h2,sp,aa,bb,B,B1,B2,du_0_dy(jmax),du_0_dy_max,delta_w,
     &  py_aux(jmax),dcy_dpy_aux(jmax),d2cy_dpy2_aux(jmax),
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax)

      allocate(rho_0(local_imax,jmax,kmax))
      allocate(u_0  (local_imax,jmax,kmax))
      allocate(v_0  (local_imax,jmax,kmax))
      allocate(tp_0 (local_imax,jmax,kmax))

      allocate(stx_ex_6th_D (local_imax,7))
      allocate(stxx_ex_6th_D(local_imax,7))

      allocate(dcy_dpy  (jmax))
      allocate(d2cy_dpy2(jmax))

      allocate(E1_t(local_imax,jmax,kmax))
      allocate(E2_t(local_imax,jmax,kmax))
      allocate(E3_t(local_imax,jmax,kmax))
      allocate(E4_t(local_imax,jmax,kmax))
      allocate(E5_t(local_imax,jmax,kmax))

      allocate(E1_rk1(local_imax,jmax,kmax))
      allocate(E1_rk2(local_imax,jmax,kmax))
      allocate(E1_rk3(local_imax,jmax,kmax))
      allocate(E1_rk4(local_imax,jmax,kmax))

      allocate(E2_rk1(local_imax,jmax,kmax))
      allocate(E2_rk2(local_imax,jmax,kmax))
      allocate(E2_rk3(local_imax,jmax,kmax))
      allocate(E2_rk4(local_imax,jmax,kmax))

      allocate(E3_rk1(local_imax,jmax,kmax))
      allocate(E3_rk2(local_imax,jmax,kmax))
      allocate(E3_rk3(local_imax,jmax,kmax))
      allocate(E3_rk4(local_imax,jmax,kmax))

      allocate(E4_rk1(local_imax,jmax,kmax))
      allocate(E4_rk2(local_imax,jmax,kmax))
      allocate(E4_rk3(local_imax,jmax,kmax))
      allocate(E4_rk4(local_imax,jmax,kmax))

      allocate(E5_rk1(local_imax,jmax,kmax))
      allocate(E5_rk2(local_imax,jmax,kmax))
      allocate(E5_rk3(local_imax,jmax,kmax))
      allocate(E5_rk4(local_imax,jmax,kmax))

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
     &  dtp_dz (local_imax,jmax,kmax),dmu_dx (local_imax,jmax,kmax),
     &  dmu_dy (local_imax,jmax,kmax),dmu_dz (local_imax,jmax,kmax))

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

      do k=1,kmax
        do j=1,jmax         
          do i=1,local_imax
            E1   (i,j,k) = 0.d0
            E2   (i,j,k) = 0.d0
            E3   (i,j,k) = 0.d0
            E4   (i,j,k) = 0.d0
            E5   (i,j,k) = 0.d0

            rho_0(i,j,k) = 0.d0
            u_0  (i,j,k) = 0.d0
            v_0  (i,j,k) = 0.d0
            tp_0 (i,j,k) = 0.d0
          end do
        end do
      end do

cccccc ************************************
cccccc ************************************
cccccc DOMAIN IN STREAMWISE DIRECTION
cccccc ************************************
cccccc ************************************
      !MS$IF(tp_s.EQ.1).OR.(tp_s.EQ.2).OR.(tp_s.EQ.3).OR.(tp_s.EQ.4)             
        do i=1,local_imax
          x(i) = dble(i-1)*dx
        end do
      !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)        
        !MS$IF (stretching_in_x.EQ.0)
          do i=1,local_imax
            if (my_rank.EQ.0) then 
              local_i = global_x1-1+i-1
            else if ((my_rank+1).EQ.pro) then              
              local_i = global_x1-1-3+i-1
            else
              local_i = global_x1-1-3+i-1
            end if
            x(i) = x_0 + dble(local_i)*dx
          end do         
        !MS$ELSEIF (stretching_in_x.EQ.1)                      
          imax_spo = imax-imax_phy    !Sponge domain - Points from imax_phy..(imax_phy+imax_spo)
          Lx_spo   = 700.d0           !Sponge domain extend to Lx_spo
          sp       = 004.d0           !Stretching parameter

          do i=1,imax_phy+1
            x(i) = x_0 + dble(i-1)*Lx/dble(imax_phy-1)
          end do 

          aa = ((Lx_spo-x(imax_phy))-dble(imax_spo)*
     &      (x(imax_phy+1)-x(imax_phy)))/
     &      (dble(imax_spo)**sp-dble(imax_spo))
          bb = (x(imax_phy+1)-x(imax_phy))-aa 

          do i=imax_phy,imax
            x(i) = x(imax_phy)+
     &        (aa*dble(i-imax_phy)**sp+bb*dble(i-imax_phy))
          end do
        !MS$ENDIF  
      !MS$ENDIF 

cccccc ************************************
cccccc ************************************
cccccc DOMAIN IN NORMALWISE DIRECTION
cccccc ************************************
cccccc ************************************
      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
        do j=1,jmax 
          y(j)         = dble(j-1)*dy
          dcy_dpy(j)   = 1.d0
          d2cy_dpy2(j) = 0.d0
        end do
      !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        do j=1,jmax 
          y(j)       = dble(j-1)*dy 
          du_0_dy(j) = c*(M1-M2)*(1.d0-tanh(2.d0*y(j))**2.d0)  
        end do
        call maximo_1D(jmax,du_0_dy,du_0_dy_max)
        delta_w  = c*(M1-M2)/du_0_dy_max !delta_w = (M1-M2)*c/du0_dy_m) = 1.0
        !MS$IF (stretching_in_y.EQ.0)
          do j=1,jmax 
            y(j)         = (y(j)-y(jmax)/2.d0)/delta_w 
            dcy_dpy(j)   = 1.d0
            d2cy_dpy2(j) = 0.d0
          end do    
        !MS$ELSEIF (stretching_in_y.EQ.1)  
          sp    = 06.d0
          py_h2 = y(jmax/2.d0)
          B     = 1.d0/(2.d0*sp)*dlog( 
     &      (1.d0+(dexp(+sp)-1.d0)*(py_h2/Ly))/  
     &      (1.d0+(dexp(-sp)-1.d0)*(py_h2/Ly)))
 
          do j=1,jmax 
            cyy       = dble(j-1)*dcy
            py_aux(j) = py_h2*(1.d0+(dsinh(sp*(cyy-B))/dsinh(sp*B)))

            dcy_dpy_aux(j)   = 
     &        dsinh(sp*B)/((sp*py_h2)*
     &        dsqrt(1.d0+(py_aux(j)/py_h2-1.d0)**2.d0*
     &        dsinh(sp*B)**2.d0))
            d2cy_dpy2_aux(j) = 
     &        (-dsinh(sp*B)**3.d0*(py_aux(j)/py_h2-1.d0))
     &        /(sp*py_h2**2.d0*( 1.d0+(py_aux(j)/py_h2-1.d0)**2.d0
     &        *dsinh(sp*B)**2.d0)**(3.d0/2.d0) ) 
          end do    

          do j=1,jmax 
            y        (j) = (py_aux       (j) - py_aux       (jmax)/2.d0)
            dcy_dpy  (j) = (dcy_dpy_aux  (j) - dcy_dpy_aux  (jmax)/2.d0)
            d2cy_dpy2(j) = (d2cy_dpy2_aux(j) - d2cy_dpy2_aux(jmax)/2.d0)
          end do
        !MS$ENDIF 
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)
          do j=1,jmax 
            y(j)         = dble(j-1)*dy 
            dcy_dpy(j)   = 1.d0
            d2cy_dpy2(j) = 0.d0
          end do
        !MS$ELSEIF (stretching_in_y.EQ.1)            
          B  = 1.04d0 
          B1 = ((B+1.d0)/(B-1.d0))
          B2 = dlog(B1)

          do j=1,jmax 
            cyy       = dble(j-1)*dcy
            py_aux(j) = Ly*( (B+1.d0)-(B-1.d0)*(B1**(1.d0-cyy)) )
     &        /(B1**(1.d0-cyy)+1.d0)
                       
            dcy_dpy_aux(j)   = 2.d0*Ly*B/
     &        (B*Ly-Ly+py_aux(j))/(B*Ly+Ly-py_aux(j))/B2
            d2cy_dpy2_aux(j) = 4.d0*B*Ly*(py_aux(j)-Ly)/
     &        (B*Ly-Ly+py_aux(j))**2/(B*Ly+Ly-py_aux(j))**2.d0/B2
          end do    

          do j=1,jmax 
            y        (j) = py_aux(j)
            dcy_dpy  (j) = dcy_dpy_aux(j)
            d2cy_dpy2(j) = d2cy_dpy2_aux(j)
          end do
        !MS$ENDIF 
      !MS$ENDIF 


ccccc ************************************
ccccc ************************************
ccccc DOMAIN IN SPANWISE DIRECTION
ccccc ************************************
ccccc ************************************
      do k=1,kmax
        z(k)=DBLE(k-1)*dz
      end do

ccccc ************************************
ccccc ************************************
      call init_grid(my_rank,pro,local_imax,x)

      !MS$IF (sv_IOS.EQ.1)
        call write_grid(x,y,z)
      !MS$ENDIF 

      !MS$IF (sv_TEC_BIN.EQ.1).OR.(sv_TEC_ASC.EQ.1)
        open (1,file=dirgraf//'x-grid-'//char(my_rank+48)//'.dat',
     &    access='sequential')
        write(1,*) 'VARIABLES= n,df'
        write(1,*) 'ZONE I=',local_imax,',F=POINT'
        do i=1,local_imax-1
          if (my_rank.EQ.0) then 
            local_i = global_x1-1+i-1
          else if ((my_rank+1).EQ.pro) then              
            local_i = global_x1-1-3+i-1
          else
            local_i = global_x1-1-3+i-1
          end if
          write(1,*) DBLE(local_i),abs(x(i+1)-x(i))
        end do
        write(1,*) DBLE(local_imax),abs(x(local_imax)-x(local_imax-1))
        close (unit=1)
        if (my_rank.EQ.0) then
          open (1,file=dirgraf//'y-grid.dat',access='sequential')
          write(1,*) 'VARIABLES= n,df'
          write(1,*) 'ZONE I=',jmax+1,',F=POINT'
          write(1,*) DBLE(0),abs(y(0)-y(1))
          do j=1,jmax
            write(1,*) DBLE(j),abs(y(j)-y(j-1))
          end do
          close (unit=1)
        end if   
      !MS$ENDIF 

      RETURN 
      END

ccccc **********************************************************************
ccccc init_grid: This routine calculates the stencils for all derivatives.
ccccc **********************************************************************
      SUBROUTINE init_grid(my_rank,pro,local_imax,x)

      USE snmethod, only: stencil_x_ex_6th_D,stencil_xx_ex_6th_D,
     &  stx_ex_6th_D,stxx_ex_6th_D

      IMPLICIT NONE

      integer my_rank,pro,local_imax
      real*8 x(local_imax)

      call stencil_x_ex_6th_D (my_rank,pro,local_imax,x,stx_ex_6th_D )
      call stencil_xx_ex_6th_D(my_rank,pro,local_imax,x,stxx_ex_6th_D)

      RETURN
      END 

ccccc **********************************************************************
ccccc cond_CFL: Determine the CFL condition.
ccccc **********************************************************************
      SUBROUTINE cond_CFL(local_imax,x,y,z,E1,E2,E3,E4,E5)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer i,j,k,local_imax
      real*8 x(local_imax),y(jmax),z(kmax),
     &  dt_max_1,dt_max_2,u_m,v_m,w_m,dx_aux,dy_aux,dz_aux,
     &  u (local_imax,jmax,kmax),v (local_imax,jmax,kmax),
     &  w (local_imax,jmax,kmax),E1(local_imax,jmax,kmax),
     &  E2(local_imax,jmax,kmax),E3(local_imax,jmax,kmax),
     &  E4(local_imax,jmax,kmax),E5(local_imax,jmax,kmax)
      
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            u(i,j,k) = E2(i,j,k)/E1(i,j,k)
            v(i,j,k) = E3(i,j,k)/E1(i,j,k)
            w(i,j,k) = E4(i,j,k)/E1(i,j,k)
          end do
        end do
      end do

      call maximo_3D(local_imax,jmax,kmax,u,u_m) 
      call maximo_3D(local_imax,jmax,kmax,v,v_m) 
      call maximo_3D(local_imax,jmax,kmax,w,w_m)   

      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            if (i.LT.local_imax) then 
              dx_aux=abs(x(i+1)-x(i))
            end if
            if (j.LT.jmax) then 
              dy_aux=abs(y(j+1)-y(j))
            end if
            if (k.LT.kmax) then 
              dz_aux=abs(z(k+1)-z(k))
            end if
           
            if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
              dt_max_1 = CFL/((u_m+1.d0/Ma)/dx_aux)
              dt_max_2 = D/(1.d0/Re*(1.d0/dx_aux**2.d0))
            else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
              dt_max_1 = CFL/(
     &          (u_m+1.d0/Ma)/dx_aux+
     &          (v_m+1.d0/Ma)/dy_aux)
              dt_max_2 = D/(1.d0/Re*(
     &          1.d0/dx_aux**2.d0+
     &          1.d0/dy_aux**2.d0))
            else
              dt_max_1 = CFL/(
     &          (u_m+1.d0/Ma)/dx_aux+
     &          (v_m+1.d0/Ma)/dy_aux+
     &          (w_m+1.d0/Ma)/dz_aux)
              dt_max_2 = D/(1.d0/Re*(
     &          1.d0/dx_aux**2.d0+
     &          1.d0/dy_aux**2.d0+
     &          1.d0/dz_aux**2.d0))
            end if            
  
            if (dt.GT.dt_max_1) then
              write(*,*)
              write(*,*) '------------------------------------------'
              write(*,*) 'ERROR: CFL Condition not satisfied.       '
              write(*,*) '------------------------------------------'
              write(*,*) 'DETAILS: DT < DT_MAX'
              write(*,*) 'DT      = ',dt
              write(*,*) 'DT_MAX  = ',dt_max_1
              write(*,'(a14,i3,i3,i3,a1)') '(I J K) =   (',i,j,k,')'
              write(*,*) '------------------------------------------'
              write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.  '
              write(*,*) '------------------------------------------'
              STOP
            end if

            if (dt.GT.dt_max_2) then
              write(*,*)
              write(*,*) '------------------------------------------'
	      write(*,*) 'ERROR: Viscous Condition not satisfied.   '
              write(*,*) '------------------------------------------'
              write(*,*) 'DETAILS: DT < DT_MAX'
              write(*,*) 'DT     = ',dt
              write(*,*) 'DT_MAX = ',dt_max_2
              write(*,'(a14,i3,i3,i3,a1)') '(I J K) =   (',i,j,k,')'
              write(*,*) '------------------------------------------'
              write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.  '
              write(*,*) '------------------------------------------' 
              STOP
            end if
          end do
        end do
      end do

      RETURN
      END

ccccc **********************************************************************
ccccc veri_data: This routine verifies the initial parameters.
ccccc **********************************************************************
      SUBROUTINE veri_data(local_imax,x,y,z,E1,E2,E3,E4,E5)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer local_imax
      real*8 x(local_imax),y(jmax),z(kmax),aux1,aux2,
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax)

      if ((jmax.EQ.1).AND.(kmax.EQ.1).AND.
     &  (kapa.NE.0).AND.(beta.NE.0)) then
        write(*,*)
        write(*,*) '------------------------------------------------'
        write(*,*) 'PARAMETERS ERRORS...                            '
        write(*,*) '------------------------------------------------'
	write(*,*) 'EXPLANATION: For 1D problem we most have kapa=0 '
	write(*,*) 'and beta=0.                                     '
        write(*,*) '------------------------------------------------'
        write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE!        '
        write(*,*) '------------------------------------------------'
        STOP
      end if

      if ((jmax.NE.1).AND.(kmax.NE.1).AND.(beta.EQ.0)) then
        write(*,*)
        write(*,*) '------------------------------------------------'
        write(*,*) 'PARAMETERS ERRORS...                            '
        write(*,*) '------------------------------------------------'
	write(*,*) 'EXPLANATION: For 3D problem we most have beta<>0'
        write(*,*) '------------------------------------------------'
        write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE!        '
        write(*,*) '------------------------------------------------'
        STOP
      end if

      if (((sv_u_m==1).OR.(sv_v_m==1).OR.
     &    (sv_rho_m==1).OR.(sv_p_m==1)).AND.(ftm.EQ.0)) then
          write(*,*)
          write(*,*) '-------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                       '
          write(*,*) '-------------------------------------------'
          write(*,*) 'EXPLANATION: Para gerar a amplitude maxima '
          write(*,*) 'da solução numérica em função do tempo é   '
	  write(*,*) 'necessario assumir que o codigo salve todos'
	  write(*,*) 'os tempos, portanto o parametro FTM        '
	  write(*,*) '(FullTime) tem que ser igual a 1.          ' 
          write(*,*) '-------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.   '
          write(*,*) '-------------------------------------------'
	  read(*,*)
          STOP
      end if

      if (((sv_e_u_m==1).OR.(sv_e_v_m==1).OR.
     &    (sv_e_rho_m==1).OR.(sv_e_p_m==1)).AND.(ftm.EQ.0)) then
          write(*,*)
          write(*,*) '-------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                       '
          write(*,*) '-------------------------------------------'
          write(*,*) 'EXPLANATION: Para gerar o erro numerico    '
          write(*,*) 'maximo e necessario assumir que o codigo   '
	  write(*,*) 'salve todos os tempos, portanto o parametro'
	  write(*,*) 'FTM (FullTime) deve ser igual a 1.         '
          write(*,*) '-------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.   '
          write(*,*) '-------------------------------------------'
	  read(*,*)
          STOP
      end if

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
      !MS$IF (tp_proc.EQ.2)
        if ((jmax.NE.1).AND.(kmax.NE.1)) then
          if ((alpha.NE.kapa).OR.(kapa.NE.beta).OR.(beta.NE.alpha)) then
            write(*,*)
            write(*,*) '-------------------------------------------'
            write(*,*) 'PARAMETERS ERRORS...                       '
            write(*,*) '-------------------------------------------'
            write(*,*) 'EXPLANATION: A solução analítica para o    '
            write(*,*) 'problema de ondas acústicas tridimensionais'
            write(*,*) 'com escoamento insetrópico foi gerada para ' 
            write(*,*) 'os seguintes parâmetros:                   '
            write(*,*) 'alpha=kapa=beta;                           '
            write(*,*) '-------------------------------------------'
            write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.   '
            write(*,*) '-------------------------------------------'
	    read(*,*)
            STOP
          end if  
        end if
      !MS$ENDIF
      !MS$ENDIF

      !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
      !MS$IF (tp_proc.EQ.1)
        if (((Re.NE.50.d0).AND.(Re.NE.100.d0).AND.(Re.NE.500.d0)).OR.
     &      (Ma.NE.1.d0).OR.(alpha.NE.2.d0*Pi).OR.(kapa.NE.2.d0*Pi).OR.
     &      (beta.NE.2.d0*Pi).OR.(Pr.NE.1.d0)) then
          write(*,*)
          write(*,*) '-------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                       '
          write(*,*) '-------------------------------------------'
          write(*,*) 'EXPLANATION: A solução analítica para o    '
          write(*,*) 'problema de ondas acústicas com a equação  '
          write(*,*) 'de energia na forma completa foi gerada    '
          write(*,*) 'para os seguintes parâmetros:              '
          write(*,*) 'Ma=1.0;                                    '
	  write(*,*) 'Re=50, 100 e 500;                          '
	  write(*,*) 'Pr=1;                                      '
	  write(*,*) 'alpha=kapa=beta=2 Pi;                      '
	  write(*,*) 'gamma=1.4;                                 '
          write(*,*) '-------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.   '
          write(*,*) '-------------------------------------------'
	  read(*,*)
          STOP
        end if
      !MS$ENDIF
      !MS$ENDIF

      !MS$IF (tp_s.EQ.3).AND.(tp_sp.EQ.1)
        if (aux1.GE.aux2) then
          write(*,*)
          write(*,*) '----------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                          '
          write(*,*) '----------------------------------------------'
          write(*,*) 'EXPLANATION: You have not to used a stretching'
          write(*,*) 'in shear zone that exceed 1.8% of dy_MAX.     '
          write(*,*) 
          write(*,*) 'Stretching:                                   '    
          write(*,*) '   At the boundaries: ', aux1
          write(*,*) '   In the shear zone: ', aux2
          write(*,*) '----------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.      '
          write(*,*) '----------------------------------------------'
	  read(*,*)
        end if  
      !MS$ENDIF

      call cond_CFL(local_imax,x,y,z,E1,E2,E3,E4,E5)

      RETURN
      END

ccccc **********************************************************************
ccccc end_data: Finalize data.
ccccc **********************************************************************
      SUBROUTINE end_data

      USE sninitial, only: rho_0,u_0,v_0,tp_0
      USE snmethod,  only: stx_ex_6th_D,stxx_ex_6th_D,
     &  global_stx_ex_6th_D,dcy_dpy,d2cy_dpy2    
      USE snconf,    only: E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_rk1,E1_rk2,E1_rk3,E1_rk4,E2_rk1,E2_rk2,E2_rk3,E2_rk4, 
     &  E3_rk1,E3_rk2,E3_rk3,E3_rk4,E4_rk1,E4_rk2,E4_rk3,E4_rk4,
     &  E5_rk1,E5_rk2,E5_rk3,E5_rk4,
     &  du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,dp_dx,
     &  dp_dy,dp_dz,drho_dx,drho_dy,drho_dz,dEt_dx,dEt_dy,dEt_dz,
     &  dtp_dx,dtp_dy,dtp_dz,dmu_dx,dmu_dy,dmu_dz,d2u_dx2,d2u_dy2,
     &  d2u_dz2,d2v_dx2,d2v_dy2,d2v_dz2,d2w_dx2,d2w_dy2,d2w_dz2,
     &  d2tp_dx2,d2tp_dy2,d2tp_dz2,d2u_dxy,d2u_dyx,d2u_dxz,d2u_dzx, 
     &  d2v_dxy,d2v_dyx,d2v_dyz,d2v_dzy,d2w_dxz,d2w_dzx,d2w_dyz,d2w_dzy 

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      deallocate(rho_0)
      deallocate(u_0  )
      deallocate(v_0  )
      deallocate(tp_0 )

      deallocate(stx_ex_6th_D )
      deallocate(stxx_ex_6th_D)

      deallocate(global_stx_ex_6th_D )

      deallocate(dcy_dpy  )
      deallocate(d2cy_dpy2)

      deallocate(E1_t)
      deallocate(E2_t)
      deallocate(E3_t)
      deallocate(E4_t)
      deallocate(E5_t)

      deallocate(E1_rk1)
      deallocate(E1_rk2)
      deallocate(E1_rk3)
      deallocate(E1_rk4)

      deallocate(E2_rk1)
      deallocate(E2_rk2)
      deallocate(E2_rk3)
      deallocate(E2_rk4)

      deallocate(E3_rk1)
      deallocate(E3_rk2)
      deallocate(E3_rk3)
      deallocate(E3_rk4)

      deallocate(E4_rk1)
      deallocate(E4_rk2)
      deallocate(E4_rk3)
      deallocate(E4_rk4)

      deallocate(E5_rk1)
      deallocate(E5_rk2)
      deallocate(E5_rk3)
      deallocate(E5_rk4)

      deallocate(du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,
     &  drho_dx,drho_dy,drho_dz,dEt_dx,dEt_dy,dEt_dz,dp_dx,dp_dy,dp_dz,
     &  dtp_dx,dtp_dy,dtp_dz,dmu_dx,dmu_dy,dmu_dz)

      deallocate(d2u_dx2,d2u_dy2,d2u_dz2,d2v_dx2,d2v_dy2,d2v_dz2,
     &  d2w_dx2,d2w_dy2,d2w_dz2,d2tp_dx2,d2tp_dy2,d2tp_dz2)

      deallocate(d2u_dxy,d2u_dyx,d2u_dxz,d2u_dzx,d2v_dxy,d2v_dyx,
     &  d2v_dyz,d2v_dzy,d2w_dxz,d2w_dzx,d2w_dyz,d2w_dzy)

      RETURN 
      END

