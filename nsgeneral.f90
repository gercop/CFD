ccccc **********************************************************************
ccccc Numtostr: Convert number to string values.
ccccc **********************************************************************
      SUBROUTINE numtostr(num,str)

      IMPLICIT NONE

      character*02 str
      integer num

      if (num.EQ.01) then
        str='01'
      else if (num.EQ.02) then
        str='02'
      else if (num.EQ.03) then
        str='03'
      else if (num.EQ.04) then
        str='04'
      else if (num.EQ.05) then
        str='05'
      else if (num.EQ.06) then
        str='06'
      else if (num.EQ.07) then
        str='07'
      else if (num.EQ.08) then
        str='08'
      else if (num.EQ.09) then
        str='09'
      else if (num.EQ.10) then
        str='10'
      end if

      RETURN
      END

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
ccccc maximo_abs_3D: Calculates maximum value.
ccccc **********************************************************************
      SUBROUTINE maximo_abs_3D(imax,jmax,kmax,var,var_max)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 var(imax,jmax,kmax),var_max

      var_max = abs(var(1,1,1))
      do k=1,kmax	 
        do j=1,jmax	 
          do i=1,imax	 
            if (abs(var(i,j,k)).GE.var_max) then
              var_max = abs(var(i,j,k))
            end if
          end do
        end do
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
ccccc view_data: View data
ccccc **********************************************************************
      SUBROUTINE view_data(local_imax,x,y,z,my_rank,pro)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

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
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 1D V37 *********'
        else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 2D V37 *********'
        else if ((jmax.NE.1).AND.(kmax.NE.1)) then
        write (n,*)'********* EXACT/MUMERICAL SOLUTION 3D V37 *********'
        end if
        write (n,*)'---------------------------------------------------'

        !MS$IF (tp_s.EQ.2)
          write(n,*)'PROGRESSIVE ACOUSTIC WAVE SIMULATION '
        !MS$ELSEIF (tp_s.EQ.3)
          write(n,*)'STANDING ACOUSTIC WAVE SIMULATION'
        !MS$ELSEIF (tp_s.EQ.4)
          write(n,*)'SHEAR LAYER SIMULATION - TEMPORAL DEVELOPMENT  '
        !MS$ELSEIF (tp_s.EQ.5)
          write(n,*)'SHEAR LAYER SIMULATION - SPATIAL DEVELOPMENT   '
        !MS$ELSEIF (tp_s.EQ.6)
          write(n,*)'BOUNDARY LAYER SIMULATION - SPATIAL DEVELOPMENT'
        !MS$ENDIF

        !MS$IF (parallel.EQ.0)
          write(n,*) 'SEQUENTIAL CODE SIMULATION'
        !MS$ELSEIF (parallel.EQ.1)
          write(n,'(a26,i6,a9)') 'PARALLEL CODE SIMULATION:',pro,' PROCCESS'
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

        !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
          write (n,'(A19,F7.2,F7.2,F7.2)') 'SPATIAL DOMAIN..: ',Lx,Ly,Lz
          write (n,'(A19,A2,F20.15)')'DX..............: ','  ',dx
          !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
            write (n,'(A19,A2,F20.15)')'DY..............: ','  ',dy
          !MS$ENDIF
          !MS$IF (tp_dim.EQ.3)
            write (n,'(A19,A2,F20.15)')'DZ..............: ','  ',dz
          !MS$ENDIF
        !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5).OR.(tp_s.EQ.6)
          !MS$IF (stretching_in_x.EQ.0)
            write (n,'(A19,A2,F20.15)')'DX..............: ','  ',dx  
          !MS$ELSEIF (stretching_in_x.EQ.1)
            write (n,'(A19,A2,F20.15)')'DX_MAX..........: ','  ',dx_max
            write (n,'(A19,A2,F20.15)')'DX_MIN..........: ','  ',dx_min
            write (n,'(A19,A2,F20.15)')'STR. PAR. IN X..: ','  ',sp_x
          !MS$ENDIF
          !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
            !MS$IF (stretching_in_y.EQ.0)
              write (n,'(A19,A2,F20.15)')'DY..............: ','  ',dy
            !MS$ELSEIF (stretching_in_y.EQ.1)
              write (n,'(A19,A2,F20.15)')'DY_MAX..........: ','  ',dy_max
              write (n,'(A19,A2,F20.15)')'DY_MIN..........: ','  ',dy_min
              write (n,'(A19,A2,F20.15)')'STR. PAR. IN Y..: ','  ',sp_y
            !MS$ENDIF
          !MS$ENDIF
          !MS$IF (tp_dim.EQ.3)
            !MS$IF (stretching_in_z.EQ.0)
              write (n,'(A19,A2,F20.15)')'DZ..............: ','  ',dz
            !MS$ELSEIF (stretching_in_z.EQ.1)         
              write (n,'(A19,A2,F20.15)')'DZ_MAX..........: ','  ',dz_max
              write (n,'(A19,A2,F20.15)')'DZ_MIN..........: ','  ',dz_min
            !MS$ENDIF
          !MS$ENDIF
        !MS$ENDIF

        write (n,'(A19,A2,F20.15)')'DT..............: ','  ',dt
        write (n,'(A19,A2,F20.15)')'LX..............: ','  ',Lx
        !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
          write (n,'(A19,A2,F20.15)')'LY..............: ','  ',Ly
        !MS$ENDIF
        !MS$IF (tp_dim.EQ.3)
          write (n,'(A19,A2,F20.15)')'LZ..............: ','  ',Lz
        !MS$ENDIF
        write (n,'(A19,A4,I6,I6,I6)')'IMAX/JMAX/KMAX..: ','    ',
     &    imax,jmax,kmax
        write (n,'(A19,A6,I16)'   )'TMAX............: ' ,'    ',tmax

        write (n,*)'--------------------------------------------------'
        write (n,*)'PHYSICS PROPERTIES'
        write (n,*)'--------------------------------------------------'

        !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5)
          write (n,*) 'PERFIL..........:     HYPERBOLIC TANGENT' 
        !MS$ELSEIF (tp_s.EQ.6)
          write (n,*) 'PERFIL..........:        PROFCOM PROGRAM'
        !MS$ENDIF

        write (n,'(A19,f22.05)') 'MACH............: ' , Ma
        write (n,'(A19,f22.08)') 'REYNOLDS........: ' , Re
        write (n,'(A19,f22.08)') 'PRANDTL.........: ' , Pr

        !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5)
          write (n,'(A19,f22.05)') 'MACH_1 (Y>0)....: ' , M1
          write (n,'(A19,f22.05)') 'MACH_2 (Y<0)....: ' , M2
          write (n,'(A19,f22.05)') 'MACH CONVECTIVO.: ' , (M1-M2)/2.d0
          write (n,'(A19,f22.08)') 'SIGMA...........: ' , sigma
        !MS$ENDIF

        write (n,*)'--------------------------------------------------'
        !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
          write (n,*) 'WAVES PROPERTIES'        
        !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5).OR.(tp_s.EQ.6)
          write (n,*) 'DISTURBANCE PROPERTIES'
        !MS$ENDIF
        write (n,*)'--------------------------------------------------'
        write (n,'(A19,f22.08)') 'AMPLITUDE 1.....: ' , A_1
        write (n,'(A19,f22.08)') 'AMPLITUDE 2.....: ' , A_2
        write (n,'(A19,f22.08)') 'AMPLITUDE 3.....: ' , A_3
        write (n,'(A19,f22.08)') 'AMPLITUDE 4.....: ' , A_4
        write (n,'(A19,f22.06)') 'ALPHA...........: ' , alpha
        !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
          write (n,'(A19,f22.06)') 'KAPA............: ' , kapa
        !MS$ENDIF
        !MS$IF (tp_dim.EQ.3)
          write (n,'(A19,f22.06)') 'BETA............: ' , beta
        !MS$ENDIF
        !MS$IF (tp_s.EQ.5).OR.(tp_s.EQ.6)
          write (n,'(A19,f22.06)') 'FREQ. (FUND)....: ' , omega
        !MS$ENDIF
        !MS$IF (tp_dim.EQ.3) 
          write (n,'(A19,f22.06)') 'PHI.............: ' , phi
        !MS$ENDIF
        !MS$IF (tp_dim.EQ.3) 
          write (n,'(A19,f22.06)') 'THETA...........: ' , 
     &      atan(beta/alpha)*180.d0/pi
        !MS$ENDIF

        write (n,*) '**************************************************'
        write (n,*) 'SPATIAL METHOD..: 6th Order Center Compact Schemes'  
        write (n,*) 'TEMPORAL METHOD.: 4th Order Runge-Kutta Schemes   '  
        write (n,*) 'FORMULATION.....: Conservative                  '

        !MS$IF (tp_proc.EQ.1)
          write (n,*) 'HYPHOTESIS......: Full Energy Equation          '
        !MS$ELSEIF (tp_proc.EQ.2)
          write (n,*) 'HYPHOTESIS......: Isentropic Flow               '
        !MS$ENDIF

        !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5)
          !MS$IF (cancel_alongvis.EQ.0)
            write (n,*) '                  Viscous Alongment Introduced'
          !MS$ELSEIF (cancel_alongvis.EQ.1)
            write (n,*) '                  Viscous Alongment Canceled  '
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
ccccc int2str: Convert integer for string data.
ccccc **********************************************************************
      SUBROUTINE int2str(tt,str)

      IMPLICIT NONE

      character*10 str
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
        str='000000000'//c1
      else if (tt.LT.100) then
        write (c2,'(I2)'),int(tt)
        str='00000000'//c2
       else if (tt.LT.1000) then
        write (c3,'(I3)'),int(tt)
        str='0000000'//c3
       else if (tt.LT.10000) then
        write (c4,'(I4)'),int(tt)
        str='000000'//c4
       else if (tt.LT.100000) then
        write (c5,'(I5)'),int(tt)
        str='00000'//c5
       else if (tt.LT.1000000) then
        write (c6,'(I6)'),int(tt)
        str='0000'//c6
       else if (tt.LT.10000000) then
        write (c7,'(I7)'),int(tt)
        str='000'//c7
       else if (tt.LT.100000000) then
        write (c8,'(I8)'),int(tt)
        str='00'//c8
      end if

      RETURN
      END

ccccc **********************************************************************
ccccc get_filename: Get filename.
ccccc **********************************************************************
      SUBROUTINE get_filename(tt,fname,ext,fname_aux)

      IMPLICIT NONE

      character*10 str
      character*05 fname
      character*20 fname_aux
      character*04 ext
      integer tt

      call int2str(tt,str)
      fname_aux = fname//'_'//str//ext

      RETURN
      END

ccccc **********************************************************************
cccccc read_BIN_1D: Read data in binary file.
ccccc **********************************************************************
      SUBROUTINE read_BIN_1D(dir,fname,U1,U2,U5,error)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*07 dir
      character*20 fname
      integer error
      real*8 
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),U5(imax,jmax,kmax)

      error=1
      open(unit=1,file=dir//'f_'//fname,status='old',
     &  form='unformatted',err=10)
      read(1) U1
      read(1) U2
      read(1) U5
      close(unit=1)

      write(*,*) 
      write(*,*) ' READ DATA IN FILE.: ',dir//'f_'//fname                
      error=0

   10 continue
      if (error.ne.0) then
        write(*,*) ' ARQUIVO ', dir//'f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
cccccc read_BIN_2D: Read data in binary file.
ccccc **********************************************************************
      SUBROUTINE read_BIN_2D(dir,fname,U1,U2,U3,U5,error)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*07 dir
      character*20 fname
      integer error
      real*8 
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U5(imax,jmax,kmax)

      error=1
      open(unit=1,file=dir//'f_'//fname,status='old',
     &  form='unformatted',err=10)
      read(1) U1
      read(1) U2
      read(1) U3
      read(1) U5
      close(unit=1)

      write(*,*) ' READ DATA IN FILE.: ',dir//'f_'//fname                
      error=0

   10 continue
      if (error.ne.0) then
        write(*,*) ' ARQUIVO ', dir//'f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
cccccc read_BIN_3D: Read data in binary file.
ccccc **********************************************************************
      SUBROUTINE read_BIN_3D(dir,fname,U1,U2,U3,U4,U5,error)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*07 dir
      character*20 fname
      integer error
      real*8 
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U4(imax,jmax,kmax),
     &  U5(imax,jmax,kmax)

      error=1
      open(unit=1,file=dir//'f_'//fname,status='old',
     &  form='unformatted',err=10)
      read(1) U1
      read(1) U2
      read(1) U3
      read(1) U4
      read(1) U5
      close(unit=1)

      write(*,*) ' READ DATA IN FILE.: ',dir//'f_'//fname                
      error=0

   10 continue
      if (error.ne.0) then
        write(*,*) ' ARQUIVO ', dir//'f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
cccccc sv_profile: Write profile in ASC file.
ccccc **********************************************************************
      SUBROUTINE sv_profile(y,rho_0,u_0,v_0,tp_0)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer j
      real*8 y(jmax),
     &  rho_0(1,jmax,kmax),u_0 (1,jmax,kmax),
     &  v_0  (1,jmax,kmax),tp_0(1,jmax,kmax)

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
      SUBROUTINE sv_BIN_1D(fname,U1,U2,U5)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      real*8
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U5(imax,jmax,kmax)

      open(1,file=dirgraf//'f_'//fname,form='unformatted',status='new')
      write(1) U1
      write(1) U2
      write(1) U5
      close(1)

      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if      

      RETURN
      END

ccccc **********************************************************************
cccccc sv_BIN_2D: Write data in binary file.
ccccc **********************************************************************
      SUBROUTINE sv_BIN_2D(fname,U1,U2,U3,U5)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      real*8
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U5(imax,jmax,kmax)

      open(1,file=dirgraf//'f_'//fname,form='unformatted',status='new')
      write(1) U1
      write(1) U2
      write(1) U3
      write(1) U5
      close(1)

      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if      

      RETURN
      END

ccccc **********************************************************************
cccccc sv_BIN_3D: Write data in binary file.
ccccc **********************************************************************
      SUBROUTINE sv_BIN_3D(fname,U1,U2,U3,U4,U5)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      real*8
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U4(imax,jmax,kmax),
     &  U5(imax,jmax,kmax)

      open(1,file=dirgraf//'f_'//fname,form='unformatted',status='new')
      write(1) U1
      write(1) U2
      write(1) U3
      write(1) U4
      write(1) U5
      close(1)

      if (QTimes.GT.1) then
        write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname
      end if      

      RETURN
      END

ccccc **********************************************************************
cccccc sv_TCP_ALL_PAR_SBIN_3D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
cc      SUBROUTINE sv_TCP_ALL_PAR_SBIN_3D(fname_aux,x,y,z,U1,U2,U3,U4,U5)

cc      IMPLICIT NONE
cc      INCLUDE 'nspar.f90'

cc      character*20 fname_aux
cc      character*1 NULLCHR
cc      integer i,j,k,Debug,IORE,xyz,NPts,NElm,
cc     &  TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd,VIsDouble
cc      real*8 Ma_aux,x(imax),y(jmax),z(kmax),
cc     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
cc     &  U3(imax,jmax,kmax),U4(imax,jmax,kmax),
cc     &  U5(imax,jmax,kmax),
cc     &  xx  (imax,jmax,kmax),yy  (imax,jmax,kmax),
cc     &  zz  (imax,jmax,kmax), 
cc     &  u   (imax,jmax,kmax),v   (imax,jmax,kmax),
cc     &  w   (imax,jmax,kmax),rho (imax,jmax,kmax),
cc     &  p   (imax,jmax,kmax),tp  (imax,jmax,kmax),
cc     &  Et  (imax,jmax,kmax),wx  (imax,jmax,kmax),
cc     &  wy  (imax,jmax,kmax),wz  (imax,jmax,kmax),
cc     &  Qi  (imax,jmax,kmax)       

cc      do k=1,kmax 
cc        do j=1,jmax 
cc          do i=1,imax
cc            u (i,j,k) = U2(i,j,k)/U1(i,j,k)
cc            v (i,j,k) = U3(i,j,k)/U1(i,j,k)
cc            w (i,j,k) = U4(i,j,k)/U1(i,j,k)
cc            !MS$IF (tp_proc.EQ.1)
cc              p (i,j,k) = (gamma-1.d0)*(U5(i,j,k)-1.d0/(2.d0*U1(i,j,k))*
cc     &          (U2(i,j,k)**2.d0+U3(i,j,k)**2.d0+U4(i,j,k)**2.d0))
cc            !MS$ELSEIF (tp_proc.EQ.2)
cc              if (U1(i,j,k).LT.0) then
cc                write(*,*) 
cc                write(*,*) 'ERROR..: Negative Density in Point (',i,j,k,').'	       
cc                write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.'	                    
cc                STOP          
cc              end if 
cc              p (i,j,k) = 1.d0/(gamma*Ma**2.d0)*U1(i,j,k)**gamma
cc            !MS$ENDIF
cc            tp(i,j,k) = gamma*Ma**2.d0*p(i,j,k)/U1(i,j,k)
cc          end do
cc        end do
cc      end do    

cc        NULLCHR   = CHAR(0)
cc        Debug     = 0  
cc        VIsDouble = 1
    
cc        IORE = TecIni('SIMPLE DATASET'//NULLCHR,
cc     &    'x y z u v w wx wy wz Q rho p T Et'//NULLCHR,
cc     &    dirgraf//'f_'//fname_aux//NULLCHR,
cc     &    '.'//NULLCHR,Debug,VIsDouble)

cc        IORE = TecZne('Simple Zone'//NULLCHR,imax,jmax,kmax,
cc     &    'BLOCK'//NULLCHR,NULLCHR)

cc        xyz = dble((imax)*(jmax)*(kmax))

cc        do i=1,imax
cc          do j=1,jmax
cc            do k=1,kmax
cc              xx(i,j,k) = x(i)
cc              yy(i,j,k) = y(j)
cc              zz(i,j,k) = z(k)
cc            end do
cc          end do
cc        end do
  
cc        call equation_vorticity(u,v,w,wx,wy,wz)
cc        call calculate_iso(u,v,w,Qi)

cc        IORE   = TecDat(xyz,xx, VIsDouble)
cc        IORE   = TecDat(xyz,yy, VIsDouble)
cc        IORE   = TecDat(xyz,zz, VIsDouble)
cc        IORE   = TecDat(xyz,u,  VIsDouble)
cc        IORE   = TecDat(xyz,v,  VIsDouble)
cc        IORE   = TecDat(xyz,w,  VIsDouble)
cc        IORE   = TecDat(xyz,wx, VIsDouble)
cc        IORE   = TecDat(xyz,wy, VIsDouble)
cc        IORE   = TecDat(xyz,wz, VIsDouble)
cc        IORE   = TecDat(xyz,Qi, VIsDouble)
cc        IORE   = TecDat(xyz,rho,VIsDouble)
cc        IORE   = TecDat(xyz,p,  VIsDouble)
cc        IORE   = TecDat(xyz,tp, VIsDouble)    
cc        IORE   = TecDat(xyz,Et, VIsDouble)
cc        IORE   = TecEnd()

cc        if (QTimes.GT.1) then
cc          write(*,*) 'SAVED FILE AS: ',dirgraf//'f_'//fname_aux
cc        end if

cc      RETURN
cc      END

ccccc **********************************************************************
cccccc sv_SN_BIN: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE sv_SN_BIN(fname,x,y,z,tt,U1,U2,U3,U4,U5)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*05 fname
      character*20 fname_aux
      integer tt           
      real*8 x(imax),y(jmax),z(kmax),
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U4(imax,jmax,kmax),
     &  U5(imax,jmax,kmax)

      call get_filename(tt,fname,'.dat',fname_aux)

      if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
        call sv_BIN_1D(fname_aux,U1,U2,U5)
      else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
        call sv_BIN_2D(fname_aux,U1,U2,U3,U5)
      else if ((jmax.NE.1).AND.(kmax.NE.1)) then
        call sv_BIN_3D(fname_aux,U1,U2,U3,U4,U5)        
      end if

cc      call get_filename(tt,fname,'.plt',fname_aux)
cc      call sv_TCP_ALL_PAR_SBIN_3D(fname_aux,x,y,z,U1,U2,U3,U4,U5)
       
      RETURN
      END


