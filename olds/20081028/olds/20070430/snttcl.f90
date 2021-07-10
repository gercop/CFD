      PROGRAM snttcl

      USE snmpi, only: my_rank,pro,local_imax
      USE snmpi, only: initmpi,time_start,decompose,time_end

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      call initmpi
      call time_start
      call decompose(imax)
      call main(my_rank,pro,local_imax)
      call time_end

      END

ccccc **********************************************************************
ccccc main: Main subroutine
ccccc **********************************************************************
      SUBROUTINE main(my_rank,pro,local_imax)
    
      USE sninitial, only: init_cond,E1,E2,E3,E4,E5
      USE snconf,    only: RK4

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      integer tt,my_rank,pro,local_imax
      real*8 x(local_imax),y(jmax),z(kmax)

      call init_data(x,y,z) 
      call init_cond(my_rank,local_imax,x,y,z)    
      call veri_data(local_imax,x,y,z)     
      call view_data(my_rank,pro,local_imax,x,y,z)  
   
      call sv_all(x,y,z,tmin)
      if (my_rank.EQ.0) then
        write(*,*)
        write(*,'(A35,I5)') 'INITIAL CONDICIONAL........: tt = ',tmin
        write(*,*)
      end if

      do tt=tmin+1,tmax
        if (my_rank.EQ.0) write(*,'(A29,A1,I8,A1,f17.10,A1,f15.10)')
     &    'NUMERICAL SOLUTION (tt,t) = ','(',tt,',',dble(tt)*dt,')'

        call RK4(local_imax,x,y,z,tt)

        if (((mod(tt,qtimes).EQ.0)).AND.(tt.GT.0)) then
          call sv_all(x,y,z,tt)
        end if
      end do

      call end_data

      RETURN
      END

ccccc **********************************************************************
ccccc init_data: Initialize data.
ccccc **********************************************************************
      SUBROUTINE init_data(x,y,z)

      USE snmpi,     only: my_rank,pro,local_imax,global_x1,ghost_points
      USE snmethod,  only: dcy_dpy,d2cy_dpy2,stx_ex_6th_D,stxx_ex_6th_D
      USE sninitial, only: rho_0,u_0,v_0,tp_0,rho,u,v,w,Et,p,tp,mu,
     &  E1,E2,E3,E4,E5
      USE snconf,    only: E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_rk1,E1_rk2,E1_rk3,E1_rk4,E2_rk1,E2_rk2,E2_rk3,E2_rk4, 
     &  E3_rk1,E3_rk2,E3_rk3,E3_rk4,E4_rk1,E4_rk2,E4_rk3,E4_rk4,
     &  E5_rk1,E5_rk2,E5_rk3,E5_rk4,
     &  du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,dp_dx,
     &  dp_dy,dp_dz,drho_dx,drho_dy,drho_dz,dEt_dx,dEt_dy,dEt_dz,
     &  dtp_dx,dtp_dy,dtp_dz,d2u_dx2,d2u_dy2,d2u_dz2,d2v_dx2,d2v_dy2,
     &  d2v_dz2,d2w_dx2,d2w_dy2,d2w_dz2,d2tp_dx2,d2tp_dy2,d2tp_dz2,
     &  d2u_dxy,d2u_dyx,d2u_dxz,d2u_dzx,d2v_dxy,d2v_dyx,d2v_dyz,d2v_dzy,
     &  d2w_dxz,d2w_dzx,d2w_dyz,d2w_dzy,du_0_dy,d2u_0_dy2 

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      integer i,j,k
      real*8 x(local_imax),y(jmax),z(kmax),t

      allocate(rho_0(local_imax,jmax,kmax))
      allocate(u_0  (local_imax,jmax,kmax))
      allocate(v_0  (local_imax,jmax,kmax))
      allocate(tp_0 (local_imax,jmax,kmax))

      allocate(E1(local_imax,jmax,kmax))
      allocate(E2(local_imax,jmax,kmax))
      allocate(E3(local_imax,jmax,kmax))
      allocate(E4(local_imax,jmax,kmax))
      allocate(E5(local_imax,jmax,kmax))

      allocate(rho(local_imax,jmax,kmax))
      allocate(u  (local_imax,jmax,kmax))
      allocate(v  (local_imax,jmax,kmax))
      allocate(w  (local_imax,jmax,kmax))
      allocate(Et (local_imax,jmax,kmax))
      allocate(p  (local_imax,jmax,kmax))
      allocate(tp (local_imax,jmax,kmax))
      allocate(mu (local_imax,jmax,kmax))

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

      allocate(du_0_dy  (local_imax,jmax,kmax),
     &         d2u_0_dy2(local_imax,jmax,kmax)) 

      do k=1,kmax
        do j=1,jmax         
          do i=1,local_imax
            x    (i)     = 0.d0
            y    (j)     = 0.d0
            z    (k)     = 0.d0

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

      open (1,file=dirgraf//'f_grid_x_P'//char(my_rank+48)//'.cx',access='sequential')
      do i=1,local_imax
        read(1,'(f30.16)') x(i)
      end do
      close (unit=1)

      open (1,file=dirgraf//'f_grid_y.cy',access='sequential')
      do j=1,jmax
        read(1,'(f30.16,f30.16,f30.16)') y(j),dcy_dpy(j),d2cy_dpy2(j)
      end do
      close (unit=1)

      open (1,file=dirgraf//'f_grid_z.cz',access='sequential')
      do k=1,kmax
        read(1,'(f30.16)') z(k)
      end do
      close (unit=1)

      call init_grid(my_rank,pro,local_imax,x)

      RETURN 
      END

ccccc **********************************************************************
ccccc init_grid: This routine calculates the stencils for all derivatives.
ccccc **********************************************************************
      SUBROUTINE init_grid(my_rank,pro,local_imax,x)

      USE snmethod, only: stx_ex_6th_D,stxx_ex_6th_D,
     &  stencil_x_ex_6th_D,stencil_xx_ex_6th_D     

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      integer my_rank,pro,local_imax,tp_s_aux
      real*8 x(local_imax)

      !MS$IF (tp_s.EQ.6)
        tp_s_aux = 6
      !MS$ELSE
        tp_s_aux = 5
      !MS$ENDIF 

      call stencil_x_ex_6th_D (my_rank,pro,local_imax,tp_s_aux,x,stx_ex_6th_D )
      call stencil_xx_ex_6th_D(my_rank,pro,local_imax,x,stxx_ex_6th_D)

      RETURN
      END 

ccccc **********************************************************************
ccccc cond_CFL: Determine the CFL condition.
ccccc **********************************************************************
      SUBROUTINE cond_CFL(local_imax,x,y,z)

      USE sninitial, only: E1,E2,E3,E4,E5

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      integer i,j,k,local_imax
      real*8 x(local_imax),y(jmax),z(kmax),
     &  dt_max_1,dt_max_2,u_m,v_m,w_m,dx_aux,dy_aux,dz_aux,
     &  u(local_imax,jmax,kmax),v(local_imax,jmax,kmax),
     &  w(local_imax,jmax,kmax)
      
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
      SUBROUTINE veri_data(local_imax,x,y,z)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      integer local_imax
      real*8 x(local_imax),y(jmax),z(kmax),aux1,aux2

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
          STOP
        end if  
      !MS$ENDIF

      !MS$IF (tp_s.NE.5).AND.(tp_sp.NE.6).AND.(parallel.EQ.1)
        if (aux1.GE.aux2) then
          write(*,*)
          write(*,*) '----------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                          '
          write(*,*) '----------------------------------------------'
          write(*,*) 'EXPLANATION: The code was parallelizing just  '
          write(*,*) 'for the shear and boundary layer problems.    '
          write(*,*) '----------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.      '
          write(*,*) '----------------------------------------------'
	  read(*,*)
          STOP
        end if  
      !MS$ENDIF

      call cond_CFL(local_imax,x,y,z)

      RETURN
      END

ccccc **********************************************************************
ccccc end_data: Finalize data.
ccccc **********************************************************************
      SUBROUTINE end_data

      USE sninitial, only: rho,u,v,w,Et,p,tp,mu,rho_0,u_0,v_0,tp_0,
     &  E1,E2,E3,E4,E5
      USE snmethod,  only: stx_ex_6th_D,stxx_ex_6th_D,
     &  global_stx_ex_6th_D,global_stxx_ex_6th_D,dcy_dpy,d2cy_dpy2    
      USE snconf,    only: E1_t,E2_t,E3_t,E4_t,E5_t,
     &  E1_rk1,E1_rk2,E1_rk3,E1_rk4,E2_rk1,E2_rk2,E2_rk3,E2_rk4, 
     &  E3_rk1,E3_rk2,E3_rk3,E3_rk4,E4_rk1,E4_rk2,E4_rk3,E4_rk4,
     &  E5_rk1,E5_rk2,E5_rk3,E5_rk4,
     &  du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,dp_dx,
     &  dp_dy,dp_dz,drho_dx,drho_dy,drho_dz,dEt_dx,dEt_dy,dEt_dz,
     &  dtp_dx,dtp_dy,dtp_dz,d2u_dx2,d2u_dy2,d2u_dz2,d2v_dx2,d2v_dy2,
     &  d2v_dz2,d2w_dx2,d2w_dy2,d2w_dz2,d2tp_dx2,d2tp_dy2,d2tp_dz2,
     &  d2u_dxy,d2u_dyx,d2u_dxz,d2u_dzx,d2v_dxy,d2v_dyx,d2v_dyz,d2v_dzy,
     &  d2w_dxz,d2w_dzx,d2w_dyz,d2w_dzy,du_0_dy,d2u_0_dy2

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      deallocate(rho_0)
      deallocate(u_0  )
      deallocate(v_0  )
      deallocate(tp_0 )

      deallocate(E1)
      deallocate(E2)
      deallocate(E3)
      deallocate(E4)
      deallocate(E5)

      deallocate(rho)
      deallocate(u  )
      deallocate(v  )
      deallocate(w  )
      deallocate(Et )
      deallocate(p  )
      deallocate(tp )
      deallocate(mu )

      deallocate(stx_ex_6th_D )
      deallocate(stxx_ex_6th_D)

      deallocate(global_stx_ex_6th_D )
      deallocate(global_stxx_ex_6th_D)

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
     &  dtp_dx,dtp_dy,dtp_dz)

      deallocate(d2u_dx2,d2u_dy2,d2u_dz2,d2v_dx2,d2v_dy2,d2v_dz2,
     &  d2w_dx2,d2w_dy2,d2w_dz2,d2tp_dx2,d2tp_dy2,d2tp_dz2)

      deallocate(d2u_dxy,d2u_dyx,d2u_dxz,d2u_dzx,d2v_dxy,d2v_dyx,
     &  d2v_dyz,d2v_dzy,d2w_dxz,d2w_dzx,d2w_dyz,d2w_dzy)

      deallocate(du_0_dy,d2u_0_dy2)

      RETURN 
      END

