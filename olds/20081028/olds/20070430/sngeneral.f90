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
      SUBROUTINE view_data(my_rank,pro,local_imax,x,y,z)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      integer n,i,j,k,wlabel,my_rank,pro,local_imax
      real*8 x(local_imax),y(jmax),z(kmax),
     &  dx_i(local_imax-1),dy_j(jmax-1),dz_k(kmax-1),
     &  dx_max,dy_max,dz_max,dx_min,dy_min,dz_min

      if (my_rank.NE.0) RETURN

      open(1,file=dirgraf//'INIT_DATA.TXT',status='unknown')
      do n=0,1
        write (n,*)
        write (n,*)'***************************************************'
        if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 1D V36 *********'
        else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 2D V36 *********'
        else if ((jmax.NE.1).AND.(kmax.NE.1)) then
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 3D V36 *********'
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

        do i=1,local_imax-1
          dx_i(i)=abs(x(i+1)-x(i))
        end do
        do j=1,jmax-1
          dy_j(j)=abs(y(j+1)-y(j))
        end do
        do k=1,kmax-1
          dz_k(k)=abs(z(k+1)-z(k))
        end do

        call maximo_1D(local_imax-1,dx_i,dx_max)      
        call maximo_1D(jmax-1,dy_j,dy_max)
        call maximo_1D(kmax-1,dz_k,dz_max)

        call minimo_1D(local_imax-1,dx_i,dx_min)
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
        write (n,*) 'FORMULATION.....: Conservative                  '

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

      open(1,file=dirgraf//'INIT_DATA.DAT',status='unknown',form='formatted')
      write(1,*) Ma
      write(1,*) Re
      write(1,*) pro
      close(1)

      read(*,*)

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
cccccc read_BIN_1D: Read data in binary file.
ccccc **********************************************************************
      SUBROUTINE read_BIN_1D(fname,global_E1,global_E2,global_E5,error)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      integer error
      real*8 
     &  global_E1(imax,jmax,kmax),global_E2(imax,jmax,kmax),
     &  global_E5(imax,jmax,kmax)

      error=1
      open(unit=1,file=dirgraf//'f_'//fname,status='old',
     &  form='formatted',err=10)
      read(1,*) global_E1
      read(1,*) global_E2
      read(1,*) global_E5
      close(unit=1)

      write(*,*) 
      write(*,*) 'READ FILE.: ',dirgraf//'f_'//fname                
      error=0

   10 continue
      if (error.ne.0) then
        write(*,*) 'ARQUIVO ', dirgraf//'f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
cccccc read_BIN_2D: Read data in binary file.
ccccc **********************************************************************
      SUBROUTINE read_BIN_2D(fname,global_E1,global_E2,global_E3,
     &  global_E5,error)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      integer error
      real*8 
     &  global_E1(imax,jmax,kmax),global_E2(imax,jmax,kmax),
     &  global_E3(imax,jmax,kmax),global_E5(imax,jmax,kmax)

      error=1
      open(unit=1,file=dirgraf//'f_'//fname,status='old',
     &  form='formatted',err=10)
      read(1,*) global_E1
      read(1,*) global_E2
      read(1,*) global_E3
      read(1,*) global_E5
      close(unit=1)

      write(*,*) 
      write(*,*) 'READ FILE.: ',dirgraf//'f_'//fname                
      error=0

   10 continue
      if (error.ne.0) then
        write(*,*) 'ARQUIVO ', dirgraf//'f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
cccccc read_BIN_3D: Read data in binary file.
ccccc **********************************************************************
      SUBROUTINE read_BIN_3D(fname,global_E1,global_E2,global_E3,
     &  global_E4,global_E5,error)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      integer error
      real*8 
     &  global_E1(imax,jmax,kmax),global_E2(imax,jmax,kmax),
     &  global_E3(imax,jmax,kmax),global_E4(imax,jmax,kmax),
     &  global_E5(imax,jmax,kmax)

      error=1
      open(unit=1,file=dirgraf//'f_'//fname,status='old',
     &  form='formatted',err=10)
      read(1,*) global_E1
      read(1,*) global_E2
      read(1,*) global_E3
      read(1,*) global_E4
      read(1,*) global_E5
      close(unit=1)

      write(*,*) 
      write(*,*) 'READ FILE.: ',dirgraf//'f_'//fname                
      error=0

   10 continue
      if (error.ne.0) then
        write(*,*) 'ARQUIVO ', dirgraf//'f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
cccccc sv_profile: Write profile in ASC file.
ccccc **********************************************************************
      SUBROUTINE sv_profile(y,rho_0,u_0,v_0,tp_0)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      integer j
      real*8 y(jmax),
     &  rho_0(imax,jmax,kmax),u_0 (imax,jmax,kmax),
     &  v_0  (imax,jmax,kmax),tp_0(imax,jmax,kmax)

      open (1,file=dirgraf//'f_profile.dat',access='sequential')
      write(1,*) 'VARIABLES= y,rho_0,u_0,v_0,T_0'
      write(1,*) 'ZONE I=',jmax,',','F=POINT'
      do j=1,jmax
        write(1,10) y(j),rho_0(1,j,1),u_0(1,j,1),v_0(1,j,1),tp_0(1,j,1)
      end do
      close (unit=1)

      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_profile.dat'
      end if      

   10 format(f20.15,e30.16,e30.16,e30.16,e30.16)

      RETURN
      END

ccccc **********************************************************************
cccccc sv_BIN_1D: Write data in binary file.
ccccc **********************************************************************
      SUBROUTINE sv_BIN_1D(fname)

      USE sninitial, only: global_E1,global_E2,global_E5

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      open(1,file=dirgraf//'f_'//fname,status='new',form='formatted')
      write(1,*) global_E1
      write(1,*) global_E2
      write(1,*) global_E5
      close(1)

      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if      

      RETURN
      END

ccccc **********************************************************************
cccccc sv_BIN_2D: Write data in binary file.
ccccc **********************************************************************
      SUBROUTINE sv_BIN_2D(fname)

      USE sninitial, only: global_E1,global_E2,global_E3,global_E5

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      open(1,file=dirgraf//'f_'//fname,status='new',form='formatted')
      write(1,*) global_E1
      write(1,*) global_E2
      write(1,*) global_E3
      write(1,*) global_E5
      close(1)

      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if      

      RETURN
      END

ccccc **********************************************************************
cccccc sv_BIN_3D: Write data in binary file.
ccccc **********************************************************************
      SUBROUTINE sv_BIN_3D(fname)

      USE sninitial, only: global_E1,global_E2,global_E3,global_E4,global_E5

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname

      open(1,file=dirgraf//'f_'//fname,status='new',form='formatted')
      write(1,*) global_E1
      write(1,*) global_E2
      write(1,*) global_E3
      write(1,*) global_E4
      write(1,*) global_E5
      close(1)

      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if      

      RETURN
      END

ccccc **********************************************************************
cccccc sv_SN_BIN: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_SN_BIN(fname,tt)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*05 fname
      character*20 fname_aux
      integer tt
    
      call get_filename(tt,fname,'.dat',fname_aux)

      if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
        call sv_BIN_1D(fname_aux)
      else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
        call sv_BIN_2D(fname_aux)
      else if ((jmax.NE.1).AND.(kmax.NE.1)) then
        call sv_BIN_3D(fname_aux)
      end if
       
      RETURN
      END

ccccc **********************************************************************
ccccc sv_all: Save all variables.
ccccc **********************************************************************
      SUBROUTINE sv_all(x,y,z,tt)
        
      USE snmpi,     only: my_rank,group_for_output
      USE sninitial, only: E1,E2,E3,E4,E5,global_x,
     &  global_E1,global_E2,global_E3,global_E4,global_E5

      IMPLICIT NONE   
      INCLUDE 'snparv36.f90'

      integer tt
      real*8 x(local_imax),y(jmax),z(kmax)

      call group_for_output(imax,jmax,kmax,x,global_x,E1,E2,E3,E4,E5,
     &  global_E1,global_E2,global_E3,global_E4,global_E5)

      if ( (my_rank.NE.0).OR.((tt.EQ.tmin).AND.(tmin.GT.0)) ) RETURN
      call sv_SN_BIN('SN__S',tt) 

      RETURN
      END


