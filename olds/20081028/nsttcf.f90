      PROGRAM snttcl

      USE snmpi, only: my_rank,pro,local_imax
      USE snmpi, only: initmpi,time_start,decompose,time_end

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      call initmpi
      call time_start
      call decompose(imax)

      call main(local_imax)

      call time_end

      END

ccccc **********************************************************************
ccccc main: Main subroutine
ccccc **********************************************************************
      SUBROUTINE main(local_imax)
    
      USE snmpi,     only: my_rank,pro
      USE sninitial, only: init_cond,U1,U2,U3,U4,U5
      USE snconf,    only: RK4

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer tt,local_imax
      real*8 x(local_imax),y(jmax),z(kmax)

      call init_data(x,y,z) 
      call init_cond(local_imax,x,y,z)    
      call veri_data(local_imax,x,y,z)     
      call view_data(local_imax,x,y,z,my_rank,pro)  
      
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

      USE snmpi,     only: my_rank,pro,local_imax,protostr
      USE sninitial, only: rho_0,u_0,v_0,tp_0,U1,U2,U3,U4,U5
      USE snmethod,  only: stx_d1_ex_D,stx_d2_ex_D,
     &  sty_a_cp_D,sty_b_cp_D,sty_c_cp_D,
     &  sty_a_cp_N,sty_b_cp_N,sty_c_cp_N,
     &  stz_a_cp_P,stz_b_cp_P,stz_c_cp_P,
     &  d1cy_dpy1,d2cy_dpy2

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*16 fname_aux
      character*04 ext
      integer i,j,k
      real*8 x(local_imax),y(jmax),z(kmax),t

      allocate(rho_0(local_imax,jmax,kmax))
      allocate(u_0  (local_imax,jmax,kmax))
      allocate(v_0  (local_imax,jmax,kmax))
      allocate(tp_0 (local_imax,jmax,kmax))

      allocate(U1(local_imax,jmax,kmax))
      allocate(U2(local_imax,jmax,kmax))
      allocate(U3(local_imax,jmax,kmax))
      allocate(U4(local_imax,jmax,kmax))
      allocate(U5(local_imax,jmax,kmax))

      do k=1,kmax
        do j=1,jmax         
          do i=1,local_imax
            x    (i)     = 0.d0
            y    (j)     = 0.d0
            z    (k)     = 0.d0

            U1   (i,j,k) = 0.d0
            U2   (i,j,k) = 0.d0
            U3   (i,j,k) = 0.d0
            U4   (i,j,k) = 0.d0
            U5   (i,j,k) = 0.d0

            rho_0(i,j,k) = 0.d0
            u_0  (i,j,k) = 0.d0
            v_0  (i,j,k) = 0.d0
            tp_0 (i,j,k) = 0.d0
          end do
        end do
      end do

      !MS$IF (parallel.EQ.0)
        open (1,file=dirgraf//'f_grid_x.cx',access='sequential')
      !MS$ELSEIF (parallel.EQ.1)
        call protostr(my_rank,'f_grid_x_P','.cx',fname_aux)
        open (1,file=dirgraf//fname_aux,access='sequential')
      !MS$ENDIF
      do i=1,local_imax
        read (1,'(e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,
     &      e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)') x(i),
     &      stx_d1_ex_D(i,1),stx_d1_ex_D(i,2),stx_d1_ex_D(i,3),stx_d1_ex_D(i,4),
     &      stx_d1_ex_D(i,5),stx_d1_ex_D(i,6),stx_d1_ex_D(i,7),
     &      stx_d2_ex_D(i,1),stx_d2_ex_D(i,2),stx_d2_ex_D(i,3),stx_d2_ex_D(i,4),
     &      stx_d2_ex_D(i,5),stx_d2_ex_D(i,6),stx_d2_ex_D(i,7)
      end do
      close (unit=1)

      open (1,file=dirgraf//'f_grid_y.cy',access='sequential')
      do j=1,jmax
        read(1,'(e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)') y(j),
     &    sty_a_cp_D(j),sty_b_cp_D(j),sty_c_cp_D(j),sty_a_cp_N(j),sty_b_cp_N(j),sty_c_cp_N(j),
     &    d1cy_dpy1(j),d2cy_dpy2(j)    
      end do
      close (unit=1)

      open (1,file=dirgraf//'f_grid_z.cz',access='sequential')
      do k=1,kmax
        read(1,'(e30.16,e30.16,e30.16,e30.16)') z(k),stz_a_cp_P(k),
     &    stz_b_cp_P(k),stz_c_cp_P(k)  
      end do
      close (unit=1) 

      RETURN 
      END

ccccc **********************************************************************
ccccc cond_CFL: Determine the CFL condition.
ccccc **********************************************************************
      SUBROUTINE cond_CFL(local_imax,x,y,z)

      USE sninitial, only: U1,U2,U3,U4,U5

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer i,j,k,local_imax
      real*8 x(local_imax),y(jmax),z(kmax),
     &  dt_max_1,dt_max_2,u_m,v_m,w_m,dx_aux,dy_aux,dz_aux,
     &  u(local_imax,jmax,kmax),v(local_imax,jmax,kmax),
     &  w(local_imax,jmax,kmax)
      
      do k=1,kmax
        do j=1,jmax
          do i=1,local_imax
            u(i,j,k) = U2(i,j,k)/U1(i,j,k)
            v(i,j,k) = U3(i,j,k)/U1(i,j,k)
            w(i,j,k) = U4(i,j,k)/U1(i,j,k)
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
      INCLUDE 'nspar.f90'

      integer local_imax
      real*8 x(local_imax),y(jmax),z(kmax),aux1,aux2

      !MS$IF (tp_s.EQ.4)
        !MS$IF (parallel.NE.00)
          write(*,*)
          write(*,*) '------------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                            '
          write(*,*) '------------------------------------------------'
          write(*,*) 'EXPLANATION: For temporal shear layer problem   '
 	  write(*,*) 'the code has to work sequential computer        '
          write(*,*) 'proccess (parallel=0).                           '
          write(*,*) '------------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE!        '
          write(*,*) '------------------------------------------------'
          STOP
        !MS$ENDIF
        !MS$IF (stretching_in_x.NE.00)
          write(*,*)
          write(*,*) '------------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                            '
          write(*,*) '------------------------------------------------'
          write(*,*) 'EXPLANATION: For temporal shear layer problem   '
 	  write(*,*) 'the code has to work with equidistant grid in   '
          write(*,*) 'the longitudinal direction (stretching_in_x=0). '
          write(*,*) '------------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE!        '
          write(*,*) '------------------------------------------------'
          STOP
        !MS$ENDIF
      !MS$ENDIF

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

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
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

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
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

      !MS$IF (filter_in_y.EQ.1)
        if ((imax.EQ.1).OR.(jmax.EQ.1)) then
          write(*,*)
          write(*,*) '----------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                          '
          write(*,*) '----------------------------------------------'
          write(*,*) 'EXPLANATION: The filter in the normalwise     '
          write(*,*) 'direction have only been used for 2D          '
          write(*,*) 'simulations.                                  '
          write(*,*) '----------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.      '
          write(*,*) '----------------------------------------------'
	  read(*,*)
          STOP 
        end if
      !MS$ENDIF
      !MS$IF (filter_in_z.EQ.1)
        if ((imax.EQ.1).OR.(jmax.EQ.1).OR.(kmax.EQ.1)) then
          write(*,*)
          write(*,*) '----------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                          '
          write(*,*) '----------------------------------------------'
          write(*,*) 'EXPLANATION: The filter in the spanwise       '
          write(*,*) 'direction have only been used for 3D          '
          write(*,*) 'simulations.                                  '
          write(*,*) '----------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.      '
          write(*,*) '----------------------------------------------'
	  read(*,*)
          STOP 
        end if
      !MS$ENDIF

      !MS$IF (stretching_in_y.EQ.1)
        if ((imax.EQ.1).OR.(jmax.EQ.1)) then
          write(*,*)
          write(*,*) '----------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                          '
          write(*,*) '----------------------------------------------'
          write(*,*) 'EXPLANATION: The stretching in the normalwise '
          write(*,*) 'direction have only been used for 2D          '
          write(*,*) 'simulations.                                  '
          write(*,*) '----------------------------------------------'
          write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.      '
          write(*,*) '----------------------------------------------'
	  read(*,*)
          STOP 
        end if
      !MS$ENDIF
      !MS$IF (stretching_in_z.EQ.1)
        if ((imax.EQ.1).OR.(jmax.EQ.1).OR.(kmax.EQ.1)) then
          write(*,*)
          write(*,*) '----------------------------------------------'
          write(*,*) 'PARAMETERS ERRORS...                          '
          write(*,*) '----------------------------------------------'
          write(*,*) 'EXPLANATION: The stretching in the spanwise   '
          write(*,*) 'direction have only been used for 3D          '
          write(*,*) 'simulations.                                  '
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
ccccc sv_all: Save all variables.
ccccc **********************************************************************
      SUBROUTINE sv_all(x,y,z,tt)
        
      USE snmpi,     only: my_rank,group_for_output
      USE sninitial, only: U1,U2,U3,U4,U5

      IMPLICIT NONE   
      INCLUDE 'nspar.f90'

      integer tt
      real*8 x(imax),y(jmax),z(kmax),
     &  global_U1(imax,jmax,kmax),global_U2(imax,jmax,kmax),
     &  global_U3(imax,jmax,kmax),global_U4(imax,jmax,kmax),
     &  global_U5(imax,jmax,kmax)

      call group_for_output(imax,jmax,kmax,U1,U2,U3,U4,U5,
     &  global_U1,global_U2,global_U3,global_U4,global_U5)

      if ( (my_rank.NE.0).OR.((tt.EQ.tmin).AND.(tmin.GT.0)) ) RETURN
      call sv_SN_BIN('SN__S',x,y,z,tt,global_U1,global_U2,
     &  global_U3,global_U4,global_U5) 

      RETURN
      END

ccccc **********************************************************************
ccccc end_data: Finalize data.
ccccc **********************************************************************
      SUBROUTINE end_data

      USE sninitial, only: rho_0,u_0,v_0,tp_0,U1,U2,U3,U4,U5

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      deallocate(rho_0)
      deallocate(u_0  )
      deallocate(v_0  )
      deallocate(tp_0 )

      deallocate(U1)
      deallocate(U2)
      deallocate(U3)
      deallocate(U4)
      deallocate(U5)

      RETURN 
      END

