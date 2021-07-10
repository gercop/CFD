ccccc **********************************************************************
ccccc maximo_1D: Calculates maximum value.
ccccc **********************************************************************
      SUBROUTINE maximo_1D(nmax,var,var_nmax)

      IMPLICIT NONE

      integer n,nmax
      real*8 var(nmax),var_nmax

      var_nmax = var(1)
      do n=1,nmax	 
        if (var(n)>=var_nmax) then
          var_nmax=var(n)
         end if        
      end do  			

      RETURN
      END

ccccc **********************************************************************
ccccc minimo_1D: Calculates minimum value.
ccccc **********************************************************************
      SUBROUTINE minimo_1D(nmax,var,var_nmin)

      IMPLICIT NONE

      integer n,nmax
      real*8 var(nmax),var_nmin

      var_nmin = var(1)
      do n=1,nmax	 
        if (var(n)<=var_nmin) then
          var_nmin=var(n)
         end if        
      end do  			

      RETURN
      END

ccccc **********************************************************************
ccccc maximo_3D: Calculates maximum value.
ccccc **********************************************************************
      SUBROUTINE maximo_3D(imax,jmax,kmax,var,var_max)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 var(imax,jmax,kmax),var_max

      var_max = var(1,1,1)
      do k=1,kmax	 
        do j=1,jmax	 
          do i=1,imax	 
            if (var(i,j,k).GE.var_max) then
              var_max = var(i,j,k)
            end if
          end do
        end do
      end do  			

      RETURN
      END

ccccc **********************************************************************
ccccc view_data: View data
ccccc **********************************************************************
      SUBROUTINE view_data(my_rank,x,y,z)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer n,i,j,k,wlabel,my_rank
      real*8 x(imax),y(jmax),z(kmax),
     &  dx_i(imax-1),dy_j(jmax-1),dz_k(kmax-1),
     &  dx_max,dy_max,dz_max,dx_min,dy_min,dz_min

      if (my_rank.NE.0) RETURN

      open(1,file=dirgraf//'INIT_DATA.TXT',status='unknown')
      do n=0,1
        write (n,*)
        write (n,*)'***************************************************'
        if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 1D V35 *********'
        else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 2D V35 *********'
        else if ((jmax.NE.1).AND.(kmax.NE.1)) then
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 3D V35 *********'
        end if
        write (n,*)'---------------------------------------------------'

        !MS$IF (tp_s.EQ.1)
          write(n,*)'SIMULATION WITH PROGRESSIVE ACOUSTIC WAVE'
        !MS$ELSEIF (tp_s.EQ.2)
          write(n,*)'SIMULATION WITH STANDING ACOUSTIC WAVE'
        !MS$ELSEIF (tp_s.EQ.3)
          write(n,*)'FREE SHEAR LAYER SIMULATION - TEMPORAL DEVELOPMENT'
        !MS$ELSEIF (tp_s.EQ.4)
          write(n,*)'FREE SHEAR LAYER SIMULATION - SPATIAL DEVELOPMENT'
        !MS$ELSEIF (tp_s.EQ.5)
          write(n,*)'SHEAR/BOUNDARY SIMULATION - SPATIAL DEVELOPMENT'
        !MS$ELSEIF (tp_s.EQ.6)
          write(n,*)'BOUNDARY LAYER SIMULATION - SPATIAL DEVELOPMENT'
        !MS$ENDIF

        do i=1,imax-1
          dx_i(i)=abs(x(i+1)-x(i))
        end do
        do j=1,jmax-1
          dy_j(j)=abs(y(j+1)-y(j))
        end do
        do k=1,kmax-1
          dz_k(k)=abs(z(k+1)-z(k))
        end do

        call maximo_1D(imax-1,dx_i,dx_max)      
        call maximo_1D(jmax-1,dy_j,dy_max)
        call maximo_1D(kmax-1,dz_k,dz_max)

        call minimo_1D(imax-1,dx_i,dx_min)
        call minimo_1D(jmax-1,dy_j,dy_min)
        call minimo_1D(kmax-1,dz_k,dz_min)
 
        write (n,*)'--------------------------------------------------'
        write (n,*)'COMPUTATIONAL MESH'
        write (n,*)'--------------------------------------------------'
        write (n,'(A19,A2,F20.6 )')'TEMPORAL DOMAIN.: ' , '  ', 
     &    DBLE(tmax)*dt

        !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2).OR.(tp_s.EQ.3)
          write (n,'(A19,F7.2,F7.2,F7.2)') 'SPATIAL DOMAIN..: ', 
     &      DBLE(imax)*dx,DBLE(jmax)*dy,DBLE(kmax)*dz
          write (n,'(A19,A2,F20.15)')'DX..............: ','  ',dx
          write (n,'(A19,A2,F20.15)')'DY..............: ','  ',dy
          write (n,'(A19,A2,F20.15)')'DZ..............: ','  ',dz
        !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5).OR.(tp_s.EQ.6)
          !MS$IF (stretching_in_x.EQ.0)
            write (n,'(A19,A2,F20.15)')'DX..............: ','  ',dx  
          !MS$ELSEIF (stretching_in_x.EQ.1)
            write (n,'(A19,A2,F20.15)')'DX_MAX..........: ','  ',dx_max
            write (n,'(A19,A2,F20.15)')'DX_MIN..........: ','  ',dx_min
          !MS$ENDIF
          !MS$IF (stretching_in_y.EQ.0)
            write (n,'(A19,A2,F20.15)')'DY..............: ','  ',dy
          !MS$ELSEIF (stretching_in_y.EQ.1)
            write (n,'(A19,A2,F20.15)')'DY_MAX..........: ','  ',dy_max
            write (n,'(A19,A2,F20.15)')'DY_MIN..........: ','  ',dy_min
          !MS$ENDIF
          !MS$IF (stretching_in_z.EQ.0)
            write (n,'(A19,A2,F20.15)')'DZ..............: ','  ',dz
          !MS$ELSEIF (stretching_in_z.EQ.1)         
            write (n,'(A19,A2,F20.15)')'DZ_MAX..........: ','  ',dz_max
            write (n,'(A19,A2,F20.15)')'DZ_MIN..........: ','  ',dz_min
          !MS$ENDIF
        !MS$ENDIF

        write (n,'(A19,A2,F20.15)')'DT..............: ','  ',dt
        write (n,'(A19,A2,F20.15)')'LX..............: ','  ',Lx
        write (n,'(A19,A2,F20.15)')'LY..............: ','  ',Ly
        write (n,'(A19,A2,F20.15)')'LZ..............: ','  ',Lz
        write (n,'(A19,A4,I6,I6,I6)')'IMAX/JMAX/KMAX..: ','    ',
     &    imax,jmax,kmax
        write (n,'(A19,A6,I16)'   )'TMAX............: ' ,'    ',tmax

        write (n,*)'--------------------------------------------------'
        write (n,*)'PHYSICS PROPERTIES'
        write (n,*)'--------------------------------------------------'

        !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
          write (n,*) 'PERFIL..........:   TANGENTE HIPERBÓLICO' 
        !MS$ELSEIF (tp_s.EQ.6)
          write (n,*) 'PERFIL..........:        PROFCOM PROGRAM'
        !MS$ENDIF

        write (n,'(A19,f22.05)') 'MACH............: ' , Ma
        write (n,'(A19,f22.08)') 'REYNOLDS........: ' , Re
        write (n,'(A19,f22.08)') 'PRANDTL.........: ' , Pr

        !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
          write (n,'(A19,f22.05)') 'MACH_1 (Y>0)....: ' , M1
          write (n,'(A19,f22.05)') 'MACH_2 (Y<0)....: ' , M2
          write (n,'(A19,f22.05)') 'MACH CONVECTIVO.: ' , (M1-M2)/2.d0
          write (n,'(A19,f22.08)') 'SIGMA...........: ' , sigma
        !MS$ENDIF

        write (n,'(A19,f22.08)') 'SPEC. HEAT CV...: ' , cv
        write (n,'(A19,f22.08)') 'SPEC. HEAT CP...: ' , cp
        write (n,'(A19,f22.08)') 'GAMMA...........: ' , gamma      
        write (n,*)'--------------------------------------------------'
        !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2)
          write (n,*) 'WAVES PROPERTIES'        
        !MS$ELSEIF (tp_s.EQ.3).OR.(tp_s.EQ.4)
          write (n,*) 'DISTURBANCE PROPERTIES'
        !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)
          write (n,*) 'DISTURBANCE PROPERTIES'
        !MS$ENDIF
        write (n,*)'--------------------------------------------------'
        write (n,'(A19,f22.08)') 'AMPLITUDE.......: ' , A
        write (n,'(A19,f22.02)') 'WAVE VELOCITY...: ' , c
        write (n,'(A19,f22.06)') 'ALPHA...........: ' , alpha
        write (n,'(A19,f22.06)') 'KAPA............: ' , kapa
        write (n,'(A19,f22.06)') 'BETA............: ' , beta
        if (kmax.NE.1) then
          write (n,'(A19,f22.06)') 'THETA...........: ' , 
     &      atan(beta/alpha)*180.d0/pi
        end if 
        !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5)
          write (n,'(A19,f22.06)') 'FREQUENCY.......: ' , omega
        !MS$ENDIF

        write (n,*) '**************************************************'
        write (n,*) 'SPATIAL METHOD..: 6th Order Center Compact Schemes'  
        write (n,*) 'TEMPORAL METHOD.: 4th Order Runge-Kutta Schemes   '  
        !MS$IF (tp_form.EQ.1)
          write (n,*) 'FORMULATION.....: Non-Conservative              '
        !MS$ELSEIF (tp_form.EQ.2)
          write (n,*) 'FORMULATION.....: Conservative                  '
        !MS$ENDIF

        !MS$IF (tp_proc.EQ.1)
          write (n,*) 'HYPHOTESIS......: Energy equation in full form  '
        !MS$ELSEIF (tp_proc.EQ.2)
          write (n,*) 'HYPHOTESIS......: Isentropic flow               '
        !MS$ENDIF

        !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
          !MS$IF (cancel_alongvis.EQ.0)
            write (n,*) '                  Viscous alongment introduced'
          !MS$ELSEIF (cancel_alongvis.EQ.1)
            write (n,*) '                  Viscous alongment canceled  '
          !MS$ENDIF
        !MS$ENDIF

        wlabel=0  
        !MS$IF (filter_in_x.EQ.1)
          write (n,*) 'FILTER..........: In x-direction;               '
          wlabel=1
        !MS$ENDIF
        !MS$IF (filter_in_y.EQ.1)          
          if (wlabel.EQ.0) then
            if (jmax.NE.1) write(n,*)'FILTER..........: In y-direction;'
            wlabel=1
          else          
            if (jmax.NE.1) write(n,*)'                  In y-direction;'
          end if
        !MS$ENDIF
        !MS$IF (filter_in_z.EQ.1)
          if (wlabel.EQ.0) then
            if (kmax.NE.1) write(n,*)'FILTER..........: In z-direction;'
            wlabel=1
          else 
            if (kmax.NE.1) write(n,*)'                  In z-direction;'  
          end if 
        !MS$ENDIF

        wlabel=0 
        !MS$IF (stretching_in_x.EQ.1)
          write (n,*) 'STRETCHING......: In x-direction;               '
          wlabel=1 
        !MS$ENDIF
        !MS$IF (stretching_in_y.EQ.1)
          if (wlabel.EQ.0) then
            if (jmax.NE.1) write(n,*)'STRETCHING......: In y-direction;'
            wlabel=1 
          else
            if (jmax.NE.1) write(n,*)'                  In y-direction;'
          end if
        !MS$ENDIF
        !MS$IF (stretching_in_z.EQ.1)
          if (wlabel.EQ.0) then
            if (kmax.NE.1) write(n,*)'STRETCHING......: In z-direction;'
            wlabel=1 
          else
            if (kmax.NE.1) write(n,*)'                  In z-direction;'
          end if
        !MS$ENDIF

        write (n,*) 'CODE PRECISION..: Double Precision                '
        write (n,*) '**************************************************'  
      end do
      close (unit=1)	  	 

      read(*,*) 

      !MS$IF (sv_MAT.EQ.1)
        open(1,file=dirgraf//'f_inf.m',status='unknown')
        if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
          write (1,*) 'PRG_D=1;'
        else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
          write (1,*) 'PRG_D=2;'
        else if ((jmax.NE.1).AND.(kmax.NE.1)) then
          write (1,*) 'PRG_D=3;'
        end if
        write (1,*) 'Ma=',Ma,';'
        write (1,*) 'Re=',Re,';'
        write (1,*) 'c=',c,';'
        write (1,*) 'A=',A,';'
        write (1,*) 'Lx=',Lx,';'
        write (1,*) 'Ly=',Ly,';'  
        write (1,*) 'Lz=',Lz,';'  
        write (1,*) 'alpha=',alpha,';'
        write (1,*) 'kapa=',kapa,';'  
        write (1,*) 'beta=',beta,';'  
        write (1,*) 'theta=',atan(beta/alpha)*180.d0/pi,';'  
        !MS$IF (tp_s.EQ.1)
          write (1,*) 'tp_s=1;'
        !MS$ELSEIF (tp_s.EQ.2)
          write (1,*) 'tp_s=2;'
        !MS$ELSEIF (tp_s.EQ.3)
          write (1,*) 'tp_s=3;'
        !MS$ELSEIF (tp_s.EQ.4)
          write (1,*) 'tp_s=4;'
        !MS$ELSEIF (tp_s.EQ.5)
          write (1,*) 'tp_s=5;'
        !MS$ELSEIF (tp_s.EQ.6)
          write (1,*) 'tp_s=6;'
        !MS$ENDIF
        !MS$IF (tp_form.EQ.1)
          write (1,*) 'tp_f=1;'
        !MS$ELSEIF (tp_form.EQ.2)
          write (1,*) 'tp_f=2;'
        !MS$ENDIF
        write (1,*) 'ftm=',ftm,';'
        close (unit=1)	  	 
      !MS$ENDIF
  
      RETURN
      END

ccccc **********************************************************************
ccccc equation_vorticity: This routine calculates the vorticity of flow.
ccccc **********************************************************************
      SUBROUTINE equation_vorticity(u,v,w,wx,wy,wz)

      USE snconf, only: der_x_D,der_y_WD,der_z_P

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer i,j,k
      real*8 
     &  u    (imax,jmax,kmax),v    (imax,jmax,kmax),
     &  w    (imax,jmax,kmax),wx   (imax,jmax,kmax),
     &  wy   (imax,jmax,kmax),wz   (imax,jmax,kmax),
     &  dv_dz(imax,jmax,kmax),dw_dy(imax,jmax,kmax),
     &  dw_dx(imax,jmax,kmax),du_dz(imax,jmax,kmax),
     &  du_dy(imax,jmax,kmax),dv_dx(imax,jmax,kmax)

      call der_x_D (imax,v,dv_dx)
      call der_x_D (imax,w,dw_dx)
      call der_y_WD(imax,u,du_dy)
      call der_y_WD(imax,w,dw_dy)
      call der_z_P (imax,v,dv_dz)
      call der_z_P (imax,u,du_dz)

      do k=1,kmax
        do j=1,jmax         
          do i=1,imax
            wx(i,j,k) = dv_dz(i,j,k) - dw_dy(i,j,k)
            wy(i,j,k) = dw_dx(i,j,k) - du_dz(i,j,k)
            wz(i,j,k) = du_dy(i,j,k) - dv_dx(i,j,k)
          end do 
        end do
      end do

      RETURN
      END

ccccc **********************************************************************
ccccc calculate_iso: Calculate figures with isosuperficie.
ccccc **********************************************************************
      SUBROUTINE calculate_iso(u,v,w,Qi)

      USE snconf, only: der_x_D,der_y_WD,der_y_WN,der_z_P

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer i,j,k
      real*8 
     &  u     (imax,jmax,kmax),v     (imax,jmax,kmax),
     &  w     (imax,jmax,kmax),du_dx (imax,jmax,kmax),
     &  dv_dx (imax,jmax,kmax),dw_dx (imax,jmax,kmax),
     &  du_dy (imax,jmax,kmax),dv_dy (imax,jmax,kmax),
     &  dw_dy (imax,jmax,kmax),du_dz (imax,jmax,kmax),
     &  dv_dz (imax,jmax,kmax),dw_dz (imax,jmax,kmax),
     &  Qi    (imax,jmax,kmax)

      call der_x_D (imax,u,du_dx) 
      call der_x_D (imax,v,dv_dx) 
      call der_x_D (imax,w,dw_dx) 

      call der_y_WD(imax,u,du_dy) 
      call der_y_WN(imax,v,dv_dy) 
      call der_y_WD(imax,w,dw_dy) 

      call der_z_P (imax,u,du_dz) 
      call der_z_P (imax,v,dv_dz) 
      call der_z_P (imax,w,dw_dz) 

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            Qi(i,j,k) = 
     &        -du_dx(i,j,k)*dv_dy(i,j,k)
     &        -du_dx(i,j,k)*dw_dz(i,j,k)  
     &        -dv_dy(i,j,k)*dw_dz(i,j,k)  
     &        +du_dy(i,j,k)*dv_dx(i,j,k)  
     &        +du_dz(i,j,k)*dw_dx(i,j,k)  
     &        +dv_dz(i,j,k)*dw_dy(i,j,k)  
          end do
        end do
      end do   

      RETURN
      END      

ccccc **********************************************************************
ccccc get_filename: Get filename.
ccccc **********************************************************************
      SUBROUTINE get_filename(tt,fname,ext,fname_aux)

      IMPLICIT NONE

      character*05 fname
      character*26 fname_aux
      character*04 ext
      character*01 c1
      character*02 c2
      character*03 c3
      character*04 c4
      character*05 c5
      character*06 c6
      character*07 c7
      character*08 c8
      integer tt

      if (tt.LT.10) then
        write (c1,'(I1)'),int(tt)
        fname_aux=fname//'_000000000'//c1//ext
      else if (tt.LT.100) then
        write (c2,'(I2)'),int(tt)
        fname_aux=fname//'_00000000'//c2//ext
       else if (tt.LT.1000) then
        write (c3,'(I3)'),int(tt)
        fname_aux=fname//'_0000000'//c3//ext
       else if (tt.LT.10000) then
        write (c4,'(I4)'),int(tt)
        fname_aux=fname//'_000000'//c4//ext
       else if (tt.LT.100000) then
        write (c5,'(I5)'),int(tt)
        fname_aux=fname//'_00000'//c5//ext
       else if (tt.LT.1000000) then
        write (c6,'(I6)'),int(tt)
        fname_aux=fname//'_0000'//c6//ext
       else if (tt.LT.10000000) then
        write (c7,'(I7)'),int(tt)
        fname_aux=fname//'_000'//c7//ext
       else if (tt.LT.100000000) then
        write (c8,'(I8)'),int(tt)
        fname_aux=fname//'_00'//c8//ext
      end if

      RETURN
      END

ccccc **********************************************************************
cccccc sv_TCP_ALL_PAR_SASC_1D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_TCP_ALL_PAR_SASC_1D(fname,x,E1,E2,E3,E4,E5)

      USE snconf, only: flux_to_primitive
 
      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*26 fname
      integer i
      real*8 x(imax),
     &  E1 (imax,jmax,kmax),E2 (imax,jmax,kmax),
     &  E3 (imax,jmax,kmax),E4 (imax,jmax,kmax),
     &  E5 (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  w  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  mu (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  Qi (imax,jmax,kmax)

      call flux_to_primitive(imax,E1,E2,E3,E4,E5,rho,u,v,w,Et,p,tp,mu)
      call equation_vorticity(u,v,w,wx,wy,wz)

      open (1,file=dirgraf//'f_'//fname,access='sequential')
      write(1,*) 'VARIABLES= x,u,rho,p,T,Et'
      write(1,*) 'ZONE I=',imax,',F=POINT'
      do i=1,imax
        write(1,10) x(i),u(i,1,1),rho(i,1,1),
     &    p(i,1,1),tp(i,1,1),Et(i,1,1)
      end do
      close (unit=1)
      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if

   10 format(f20.15,e30.16,e30.16,e30.16,e30.16,e30.16)

      RETURN
      END

ccccc **********************************************************************
cccccc sv_TCP_ALL_PAR_SASC_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_TCP_ALL_PAR_SASC_2D(fname,x,y,E1,E2,E3,E4,E5)

      USE snconf, only: flux_to_primitive

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*26 fname
      integer i,j
      real*8 x(imax),y(jmax),     
     &  E1 (imax,jmax,kmax),E2 (imax,jmax,kmax),
     &  E3 (imax,jmax,kmax),E4 (imax,jmax,kmax),
     &  E5 (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  w  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  mu (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  Qi (imax,jmax,kmax)

      call flux_to_primitive(imax,E1,E2,E3,E4,E5,rho,u,v,w,Et,p,tp,mu)
      call equation_vorticity(u,v,w,wx,wy,wz)

      open (1,file=dirgraf//'f_'//fname,access='sequential')
      write(1,*) 'VARIABLES=x,y,u,v,wz,rho,p,T,Et'
      write(1,*) 'ZONE I=',imax,',J=',jmax,',F=POINT'
      do j=1,jmax
        do i=1,imax          
          write(1,10) x(i),y(j),u(i,j,1),v(i,j,1),wz(i,j,1),
     &      rho(i,j,1),p(i,j,1),tp(i,j,1),Et(i,j,1)
        end do
      end do
      close (unit=1)
      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if

   10 format(f20.15,f20.15,e30.16,e30.16,
     &  e30.16,e30.16,e30.16,e30.16,e30.16)

      RETURN
      END

ccccc **********************************************************************
cccccc sv_TCP_ALL_PAR_SASC_3D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_TCP_ALL_PAR_SASC_3D(fname,x,y,z,E1,E2,E3,E4,E5)

      USE snconf, only: flux_to_primitive

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*26 fname
      integer i,j,k
      real*8 x(imax),y(jmax),z(kmax),
     &  E1 (imax,jmax,kmax),E2 (imax,jmax,kmax),
     &  E3 (imax,jmax,kmax),E4 (imax,jmax,kmax),
     &  E5 (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  w  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  mu (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  Qi (imax,jmax,kmax)

      call flux_to_primitive(imax,E1,E2,E3,E4,E5,rho,u,v,w,Et,p,tp,mu)
      call equation_vorticity(u,v,w,wx,wy,wz)
      call calculate_iso(u,v,w,Qi)

      open (1,file=dirgraf//'f_'//fname,access='sequential')
      write(1,*) 'VARIABLES= x,y,z,u,v,w,wx,wy,wz,Q,rho,p,T,Et'
      write(1,*) 'ZONE I=',imax,',J=',jmax,',K=',kmax,',F=POINT'
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            write(1,10) x(i),y(j),z(k),u(i,j,k),v(i,j,k),w(i,j,k),
     &        wx(i,j,k),wy(i,j,k),wz(i,j,k),Qi(i,j,k),
     &        rho(i,j,k),p(i,j,k),tp(i,j,k),Et(i,j,k)

          end do
        end do
      end do
      close (unit=1)
      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if

   10 format(f20.15,f20.15,f20.15,e30.16,e30.16,e30.16,e30.16,e30.16,
     &  e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)

      RETURN
      END

ccccc **********************************************************************
cccccc sv_TCP_ALL_PAR_SBIN_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_TCP_ALL_PAR_SBIN_2D(fname,x,y,E1,E2,E3,E4,E5)

      USE snconf, only: flux_to_primitive

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,j,Debug,IORE,xy,NPts,NElm,
     &  TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd,VIsDouble
      real*8 x(imax),y(jmax),
     &  xx (imax,jmax,kmax),yy (imax,jmax,kmax),
     &  zz (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  w  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  mu (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  E1 (imax,jmax,kmax),E2 (imax,jmax,kmax),
     &  E3 (imax,jmax,kmax),E4 (imax,jmax,kmax),
     &  E5 (imax,jmax,kmax),Qi (imax,jmax,kmax)

      !MS$IF (sv_TEC_BIN.EQ.1)
        call flux_to_primitive(imax,E1,E2,E3,E4,E5,rho,u,v,w,Et,p,tp,mu)
        call equation_vorticity(u,v,w,wx,wy,wz)

        NULLCHR   = CHAR(0)
        Debug     = 0
        VIsDouble = 1
    
        IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &    'x y u v wz rho p T Et'//NULLCHR,
     &    dirgraf//'f_'//fname//NULLCHR,
     &    '.'//NULLCHR,Debug,VIsDouble)

        IORE = TecZne('Simple Zone'//NULLCHR,
     &    imax,jmax,kmax,'BLOCK'//NULLCHR,NULLCHR)

        xy = dble((imax)*(jmax))

        do i=1,imax
          do j=1,jmax
            xx(i,j,1) = x(i)
            yy(i,j,1) = y(j)
          end do
        end do

        IORE   = TecDat(xy,xx, VIsDouble)
        IORE   = TecDat(xy,yy, VIsDouble)
        IORE   = TecDat(xy,u,  VIsDouble)
        IORE   = TecDat(xy,v,  VIsDouble)
        IORE   = TecDat(xy,wz, VIsDouble)
        IORE   = TecDat(xy,rho,VIsDouble)
        IORE   = TecDat(xy,p,  VIsDouble)
        IORE   = TecDat(xy,tp, VIsDouble)    
        IORE   = TecDat(xy,Et, VIsDouble)
        IORE   = TecEnd()

        if (QTimes.GT.1) then
          write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
        end if
      !MS$ENDIF

      RETURN
      END

ccccc **********************************************************************
cccccc sv_TCP_ALL_PAR_SBIN_3D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_TCP_ALL_PAR_SBIN_3D(fname,x,y,z,E1,E2,E3,E4,E5)

      USE snconf, only: flux_to_primitive

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,j,k,Debug,IORE,xyz,NPts,NElm,
     &  TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd,VIsDouble
      real*8 x(imax),y(jmax),z(kmax),
     &  xx (imax,jmax,kmax),yy (imax,jmax,kmax),
     &  zz (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  w  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  mu (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  E1 (imax,jmax,kmax),E2 (imax,jmax,kmax),
     &  E3 (imax,jmax,kmax),E4 (imax,jmax,kmax),
     &  E5 (imax,jmax,kmax),Qi (imax,jmax,kmax)

      !MS$IF (sv_TEC_BIN.EQ.1)        
        call flux_to_primitive(imax,E1,E2,E3,E4,E5,rho,u,v,w,Et,p,tp,mu)
        call equation_vorticity(u,v,w,wx,wy,wz)
        call calculate_iso(u,v,w,Qi)
  
        NULLCHR   = CHAR(0)
        Debug     = 0  
        VIsDouble = 1
    
        IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &    'x y z u v w wx wy wz Q rho p T Et'//NULLCHR,
     &    dirgraf//'f_'//fname//NULLCHR,
     &    '.'//NULLCHR,Debug,VIsDouble)

        IORE = TecZne('Simple Zone'//NULLCHR,imax,jmax,kmax,
     &    'BLOCK'//NULLCHR,NULLCHR)

        xyz = dble((imax)*(jmax)*(kmax))

        do i=1,imax
          do j=1,jmax
            do k=1,kmax
              xx(i,j,k) = x(i)
              yy(i,j,k) = y(j)
              zz(i,j,k) = z(k)
            end do
          end do
        end do
  
        IORE   = TecDat(xyz,xx, VIsDouble)
        IORE   = TecDat(xyz,yy, VIsDouble)
        IORE   = TecDat(xyz,zz, VIsDouble)
        IORE   = TecDat(xyz,u,  VIsDouble)
        IORE   = TecDat(xyz,v,  VIsDouble)
        IORE   = TecDat(xyz,w,  VIsDouble)
        IORE   = TecDat(xyz,wx, VIsDouble)
        IORE   = TecDat(xyz,wy, VIsDouble)
        IORE   = TecDat(xyz,wz, VIsDouble)
        IORE   = TecDat(xyz,Qi, VIsDouble)
        IORE   = TecDat(xyz,rho,VIsDouble)
        IORE   = TecDat(xyz,p,  VIsDouble)
        IORE   = TecDat(xyz,tp, VIsDouble)    
        IORE   = TecDat(xyz,Et, VIsDouble)
        IORE   = TecEnd()

        if (QTimes.GT.1) then
          write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
        end if
      !MS$ENDIF

      RETURN
      END

ccccc **********************************************************************
ccccc sv_TCP_ALL_PAR_SN_S: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_TCP_ALL_SN_SASC(enbd,fname,x,y,z,tt,E1,E2,E3,E4,E5)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*05 fname
      character*26 fname_aux
      integer enbd,tt
      real*8 x(imax),y(jmax),z(kmax),
     &  E1(imax,jmax,kmax),E2(imax,jmax,kmax),
     &  E3(imax,jmax,kmax),E4(imax,jmax,kmax),
     &  E5(imax,jmax,kmax)

      if (enbd.EQ.0) RETURN
      
      call get_filename(tt,fname,'.dat',fname_aux)

      if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
        call sv_TCP_ALL_PAR_SASC_1D(fname_aux,x,E1,E2,E3,E4,E5)
      else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
        call sv_TCP_ALL_PAR_SASC_2D(fname_aux,x,y,E1,E2,E3,E4,E5)
      else if ((jmax.NE.1).AND.(kmax.NE.1)) then
        call sv_TCP_ALL_PAR_SASC_3D(fname_aux,x,y,z,E1,E2,E3,E4,E5)
      end if
       
      RETURN
      END

ccccc **********************************************************************
cccccc sv_TCP_ALL_PAR_SN_SBIN: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_TCP_ALL_SN_SBIN(enbd,fname,x,y,z,tt,E1,E2,E3,E4,E5)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*05 fname
      character*20 fname_aux
      integer enbd,tt
      real*8 x(imax),y(jmax),z(kmax),
     &  E1(imax,jmax,kmax),E2(imax,jmax,kmax),
     &  E3(imax,jmax,kmax),E4(imax,jmax,kmax),
     &  E5(imax,jmax,kmax)

      if (enbd.EQ.0) RETURN
      
      call get_filename(tt,fname,'.plt',fname_aux)

      if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
        call sv_TCP_ALL_PAR_SN_1D(fname_aux,x,E1,E2,E3,E4,E5)
      else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
        call sv_TCP_ALL_PAR_SBIN_2D(fname_aux,x,y,E1,E2,E3,E4,E5)
      else if ((jmax.NE.1).AND.(kmax.NE.1)) then
        call sv_TCP_ALL_PAR_SBIN_3D(fname_aux,x,y,z,E1,E2,E3,E4,E5)
      end if
       
      RETURN
      END

ccccc **********************************************************************
cccccc sv_TCP_ALL_SN_T: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_TCP_ALL_SN_T(enbd,fname,tt,u_m,v_m,w_m,rho_m,
     &  p_m,tp_m,Et_m)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*05 fname
      integer tt,enbd
      real*8 u_m,v_m,w_m,rho_m,p_m,tp_m,s_m,Et_m

      if (enbd.EQ.0) RETURN
      if ((tt.EQ.0).OR.(ftm.EQ.0)) then
        open (1,file=dirgraf//'f_'//fname//'.dat',access='sequential')
        write(1,*) 'VARIABLES= t,u_m,v_m,w_m,rho_m,p_m,T_m,Et_m,
     &    e_u_m,e_v_m,e_w_m,e_rho_m,e_p_m,e_T_m,e_Et_m'
      else
        open (1,file=dirgraf//'f_'//fname//'.dat',access='append')
      end if
      write(1,10) dble(tt)*dt,u_m,v_m,w_m,rho_m,p_m,tp_m,Et_m
      close (unit=1)      
      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname//'.dat'
      end if

   10 format(f20.15,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)

      RETURN
      END

ccccc **********************************************************************
ccccc write_grid: Write the coordinates in IOS.
ccccc **********************************************************************
      SUBROUTINE write_grid(x,y,z)

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      character*16 fbase
      character*72 inf(100)
      integer i,j,k,mt,mm1,mm2,mm3,mp,minf,imach,itape,iunit,itimes(512)
      real*4 x_aux(imax,jmax,kmax),y_aux(imax,jmax,kmax),z_aux(imax,jmax,kmax)
      real*8 x(imax),y(jmax),z(kmax)

      !MS$IF (sv_IOS.EQ.1)
        !Initialize parameters for output
        fbase     = dirgraf//'grid'
        itape     = 10
        iunit     = 10
        imach     = 1
        mm1       = imax
        mm2       = jmax
        mm3       = kmax
        mt        = 1
        mp        = 3
        minf      = 1
        inf(1)    = 'x-grid'
        inf(2)    = 'y-grid'
        inf(3)    = 'z-grid'
        inf(4)    = 'grid file for IOS'
        itimes(1) = 1     

        call writecd(itape,iunit,fbase,imach,mt,
     &    mm3,mm2,mm1,mp,itimes,inf,minf)

        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              x_aux(i,j,k) = x(i) 
              y_aux(i,j,k) = y(j)
              z_aux(i,j,k) = z(k)
            end do
          end do
        end do  
  
        call writed(itape,x_aux)
        call writed(itape,y_aux)
        call writed(itape,z_aux)
 
        close(itape)
      !MS$ENDIF
	
      RETURN
      END      

ccccc **********************************************************************
ccccc ini_iosfile_2D: Initialize the IOS.
ccccc **********************************************************************
      SUBROUTINE ini_iosfile_2D

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer mm1,mm2,mm3,mmt,mmp,minf,imach,itimes(512)
      character*16 fbase
      character*72 inf(100)

      !MS$IF (sv_IOS.EQ.1)
        !Initialize parameters for output
        fbase     = dirgraf//'DATAIOS'
        imach     = 1
        mm1       = imax
        mm2       = jmax
        mm3       = 1
        mmt       = tmax/qtimes+1
        mmp       = 6
        minf      = 1
        inf(1)    = 'u'
        inf(2)    = 'v'
        inf(3)    = 'wz'
        inf(4)    = 'rho'
        inf(5)    = 'p'
        inf(6)    = 'T'
        inf(7)    = 'DATAIOS'
        itimes(1) = 1

        !Open file for output  
        call writecd(25,25,fbase,imach,mmt,
     &    mm3,mm2,mm1,mmp,itimes,inf,minf)
      !MS$ENDIF

      RETURN
      END

ccccc **********************************************************************
ccccc ini_iosfile_3D: Initialize the IOS.
ccccc **********************************************************************
      SUBROUTINE ini_iosfile_3D

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer mm1,mm2,mm3,mmt,mmp,minf,imach,itimes(512)
      character*16 fbase
      character*72 inf(100)

      !MS$IF (sv_IOS.EQ.1)
        !Initialize parameters for output
        fbase     = dirgraf//'DATAIOS'
        imach     = 1
        mm1       = imax
        mm2       = jmax
        mm3       = kmax
        mmt       = tmax/qtimes+1
        mmp       = 11
        minf      = 1
        inf(1)    = 'u'
        inf(2)    = 'v'
        inf(3)    = 'w'
        inf(4)    = 'wx'
        inf(5)    = 'wy'
        inf(6)    = 'wz'
        inf(7)    = 'Q'
        inf(8)    = 'rho'
        inf(9)    = 'p'
        inf(10)   = 'T'
        inf(11)   = 'E'
        inf(12)   = 'DATAIOS'
        itimes(1) = 1
  
        !Open file for output  
        call writecd(25,25,fbase,imach,mmt,
     &    mm3,mm2,mm1,mmp,itimes,inf,minf)      
      !MS$ENDIF

      RETURN
      END

ccccc **********************************************************************
ccccc ini_iosfile: Initialize the IOS.
ccccc **********************************************************************
      SUBROUTINE ini_iosfile

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      !MS$IF (sv_IOS.EQ.1)
        if ((jmax.NE.1).AND.(kmax.EQ.1)) then
          call ini_iosfile_2D
        else if ((jmax.NE.1).AND.(kmax.NE.1)) then
          call ini_iosfile_3D
        end if
      !MS$ENDIF
      
      RETURN
      END

ccccc **********************************************************************
ccccc write_ios_2D: Write data in IOS.
ccccc **********************************************************************
      SUBROUTINE write_ios_2D

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      real*8 wx(imax,jmax,kmax),wy(imax,jmax,kmax),wz(imax,jmax,kmax)

      !MS$IF (sv_IOS.EQ.1)
        call equation_vorticity(u,v,w,wx,wy,wz)

        call writed(25,real(u)  )
        call writed(25,real(v)  )
        call writed(25,real(wz) )
        call writed(25,real(rho))
        call writed(25,real(p)  )
        call writed(25,real(tp) )
      !MS$ENDIF

      RETURN
      END

ccccc **********************************************************************
ccccc write_ios_3D: Write data in IOS.
ccccc **********************************************************************
      SUBROUTINE write_ios_3D

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      real*8 
     &  wx(imax,jmax,kmax),wy(imax,jmax,kmax),
     &  wz(imax,jmax,kmax),Qi(imax,jmax,kmax)

      !MS$IF (sv_IOS.EQ.1)
        call equation_vorticity(u,v,w,wx,wy,wz)
        call calculate_iso(u,v,w,Qi)

        call writed(25,real(u))
        call writed(25,real(v))
        call writed(25,real(w))
        call writed(25,real(wx))
        call writed(25,real(wy))
        call writed(25,real(wz))
        call writed(25,real(Qi))
        call writed(25,real(rho))
        call writed(25,real(p))
        call writed(25,real(tp))
        call writed(25,real(Et))
      !MS$ENDIF

      RETURN
      END

ccccc **********************************************************************
ccccc write_ios: Write data in IOS.
ccccc **********************************************************************
      SUBROUTINE write_ios

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      !MS$IF (sv_IOS.EQ.1)
        if ((jmax.NE.1).AND.(kmax.EQ.1)) then
          call write_ios_2D
        else if ((jmax.NE.1).AND.(kmax.NE.1)) then
          call write_ios_3D
        end if
      !MS$ENDIF

      RETURN
      END

ccccc **********************************************************************
ccccc end_iosfile: Finalize the IOS.
ccccc **********************************************************************
      SUBROUTINE end_iosfile

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      !MS$IF (sv_IOS.EQ.1)
        close(25)
      !MS$ENDIF

      RETURN
      END


ccccc **********************************************************************
ccccc sv_all: Save all variables.
ccccc **********************************************************************
      SUBROUTINE sv_all(x,y,z,tt,E1,E2,E3,E4,E5)
      
      USE snmpi, only: my_rank, group_for_output

      IMPLICIT NONE   
      INCLUDE 'snparv35.f90'

      integer i,j,k,tt
      real*8 x(imax),y(jmax),z(kmax),global_x(imax),
     &  E1       (imax,jmax,kmax),E2       (imax,jmax,kmax),
     &  E3       (imax,jmax,kmax),E4       (imax,jmax,kmax),
     &  E5       (imax,jmax,kmax),
     &  global_E1(imax,jmax,kmax),global_E2(imax,jmax,kmax),    
     &  global_E3(imax,jmax,kmax),global_E4(imax,jmax,kmax),    
     &  global_E5(imax,jmax,kmax),u_m,v_m,w_m,rho_m,p_m,tp_m,Et_m

      call group_for_output(imax,jmax,kmax,x,global_x,E1,E2,E3,E4,E5,
     &  global_E1,global_E2,global_E3,global_E4,global_E5)

      if (my_rank.NE.0) RETURN
ccccc *****************************************
ccccc Rotinas para salvar dados para o TECPLOT
ccccc *****************************************
      !MS$IF (sv_TEC_BIN.EQ.1)
        call sv_TCP_ALL_SN_SBIN(sv_TCP_SNS,'SN__S',x,y,z,tt,
     &    global_E1,global_E2,global_E3,global_E4,global_E5) 
      !MS$ENDIF

      !MS$IF (sv_TEC_ASC.EQ.1)
        call sv_TCP_ALL_SN_SASC(sv_TCP_SNS,'SN__S',x,y,z,tt,
     &    global_E1,global_E2,global_E3,global_E4,global_E5)
      !MS$ENDIF

      !MS$IF (sv_TEC_BIN.EQ.1).OR.(sv_TEC_ASC.EQ.1)
        call sv_TCP_ALL_SN_T(sv_TCP_SNT,'SN__T',tt,
     &    u_m,v_m,w_m,rho_m,p_m,tp_m,Et_m)
      !MS$ENDIF

ccccc *****************************************
ccccc Rotina para salvar dados em IOS
ccccc *****************************************
      !MS$IF (sv_IOS.EQ.1)
        call write_IOS
      !MS$ENDIF 

      RETURN
      END

