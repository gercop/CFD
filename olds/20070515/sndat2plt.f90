cccccc *******************************************************************
cccccc Doctorate Project USP/EESC: Numerical simulation of the transition
cccccc                             for turbulence in a compressible 
cccccc                             boundary layer on a plane plate.
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos
cccccc Programmer\Guiding: Marcello A. Faraco de Medeiros
cccccc Version : 1.0
cccccc *******************************************************************
cccccc This program convert binary data for binary data to TECPLOT.
cccccc Comando para execução:
cccccc   ifort dat2plt.f -o dat2plt /usr/local/tecplot8/lib/tecio.a 
cccccc   /usr/local/tecplot8/lib/ctype.o
cccccc *******************************************************************
      PROGRAM sndat2plt

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      call main

      END

ccccc **********************************************************************
ccccc main: Main subroutine
ccccc **********************************************************************
      SUBROUTINE main

      USE snmethod, only: dcy_dpy,d2cy_dpy2

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname_aux
      character*05 fname
      character*04 ext
      integer tt,q_times,t_ini,t_end,error
      real*8 x(imax),y(jmax),z(kmax),
     &  stx_ex_6th_D(imax,7),
     &  rho(imax,jmax,kmax),u (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wx(imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz(imax,jmax,kmax),
     &  Qi (imax,jmax,kmax)

      fname   = 'SN__S'
      q_times = 0     
      t_ini   = 0
      t_end   = 0

      call init_grid(x,y,z,stx_ex_6th_D)

      write(*,*) 'PROGRAMA PARA CONVERTER DE BINARIO P/ TEPLOT.'  
      write(*,'(a23,\)') 'TEMPO INICIAL........:'
      read (*,*) t_ini
      write(*,'(a23,\)') 'TEMPO FINAL..........:'
      read (*,*) t_end
      write(*,'(a23,\)') 'QUANTIDADE DE TEMPOS.:'
      read (*,*) q_times       

      do tt=t_ini,t_end
        if (mod(tt,q_times).EQ.0) then
          call get_filename(tt,fname,'.dat',fname_aux)           
          if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
            call read_BIN_1D_PP(fname_aux,stx_ex_6th_D,u,rho,p,tp,Et,error)
            if (error.EQ.0) then
              call get_filename(tt,fname,'.plt',fname_aux)
              call write_TCP_1D(fname_aux,x,u,rho,p,tp,Et)
            end if
          else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
            call read_BIN_2D_PP(fname_aux,stx_ex_6th_D,u,v,wz,rho,p,tp,Et,error)
            if (error.EQ.0) then
              call get_filename(tt,fname,'.plt',fname_aux)
              call write_TCP_2D(fname_aux,x,y,u,v,wz,rho,p,tp,Et)
            end if
          else if ((jmax.NE.1).AND.(kmax.NE.1)) then
            call read_BIN_3D_PP(fname_aux,stx_ex_6th_D,u,v,w,wx,wy,wz,Qi,rho,p,tp,Et,error)
            if (error.EQ.0) then
              call get_filename(tt,fname,'.plt',fname_aux)
              call write_TCP_3D(fname_aux,x,y,z,u,v,w,wx,wy,wz,Qi,rho,p,tp,Et)
            end if
          end if
        end if
      end do

      deallocate(dcy_dpy  )
      deallocate(d2cy_dpy2)

      RETURN
      END

ccccc **********************************************************************
ccccc init_grid: This routine calculates the stencils for all derivatives.
ccccc **********************************************************************
      SUBROUTINE init_grid(x,y,z,stx_ex_6th_D)

      USE snmethod, only: dcy_dpy,d2cy_dpy2

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      integer i,j,k,pro_aux	
      real*8 x(imax),y(jmax),z(kmax),Ma_aux,Re_aux,
     &  stx_ex_6th_D(imax,7),stxx_ex_6th_D(imax,7)

      allocate(dcy_dpy  (jmax))
      allocate(d2cy_dpy2(jmax))

      open(1,file=dirgraf//'INIT_DATA.DAT',status='old',form='formatted')
      read(1,*) Ma_aux
      read(1,*) Re_aux
      read(1,*) pro_aux
      close(1)

      open(1,file=dirgraf//'f_grid_x.cx',access='sequential')
      do i=1,imax
        read (1,'(e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,
     &    e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)') x(i),
     &    stx_ex_6th_D (i,1),stx_ex_6th_D (i,2),stx_ex_6th_D (i,3),
     &    stx_ex_6th_D (i,4),stx_ex_6th_D (i,5),stx_ex_6th_D (i,6),
     &    stx_ex_6th_D (i,7),stxx_ex_6th_D(i,1),stxx_ex_6th_D(i,2),
     &    stxx_ex_6th_D(i,3),stxx_ex_6th_D(i,4),stxx_ex_6th_D(i,5),
     &    stxx_ex_6th_D(i,6),stxx_ex_6th_D(i,7)
      end do
      close (unit=1)

      open (1,file=dirgraf//'f_grid_y.cy',access='sequential')
      do j=1,jmax
        read(1,'(e30.16,e30.16,e30.16)') y(j),dcy_dpy(j),d2cy_dpy2(j)
      end do
      close (unit=1)

      open (1,file=dirgraf//'f_grid_z.cz',access='sequential')
      do k=1,kmax
        read(1,'(e30.16)') z(k)
      end do
      close (unit=1)

      RETURN
      END 

ccccc **********************************************************************
ccccc calculate_iso: Calculate figures with isosuperficie.
ccccc **********************************************************************
      SUBROUTINE calculate_iso(u,v,w,Qi)

      USE snmethod, only: stx_ex_6th_D,der_x_D,der_y_WD,der_y_WN,der_z_P

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      integer i,j,k
      real*8 dcy_dpy(jmax),
     &  u    (imax,jmax,kmax),v    (imax,jmax,kmax),
     &  w    (imax,jmax,kmax),du_dx(imax,jmax,kmax),
     &  dv_dx(imax,jmax,kmax),dw_dx(imax,jmax,kmax),
     &  du_dy(imax,jmax,kmax),dv_dy(imax,jmax,kmax),
     &  dw_dy(imax,jmax,kmax),du_dz(imax,jmax,kmax),
     &  dv_dz(imax,jmax,kmax),dw_dz(imax,jmax,kmax),
     &  Qi   (imax,jmax,kmax)

      call der_x_D (0,1,imax,stx_ex_6th_D,u,du_dx) 
      call der_x_D (0,1,imax,stx_ex_6th_D,v,dv_dx) 
      call der_x_D (0,1,imax,stx_ex_6th_D,w,dw_dx) 

      !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)                  
        call der_y_WN(imax,u,du_dy) 
        call der_y_WD(imax,v,dv_dy) 
        call der_y_WN(imax,w,dw_dy) 
      !MS$ELSEIF (tp_s.EQ.6)
        call der_y_WN(imax,u,du_dy) 
        call der_y_WN(imax,v,dv_dy) 
        call der_y_WN(imax,w,dw_dy) 
      !MS$ENDIF 
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
ccccc read_BIN_1D_PP: Read data in binary file for the post process.
ccccc **********************************************************************
      SUBROUTINE read_BIN_1D_PP(fname,stx_ex_6th_D,u,rho,p,tp,Et,error)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      integer i,j,k,error
      real*8 stx_ex_6th_D(imax,7),
     &  rho(imax,jmax,kmax),u (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),mu(imax,jmax,kmax),
     &  wx (imax,jmax,kmax),wy(imax,jmax,kmax),
     &  wz (imax,jmax,kmax),
     &  global_E1(imax,jmax,kmax),global_E2(imax,jmax,kmax),
     &  global_E3(imax,jmax,kmax),global_E4(imax,jmax,kmax),
     &  global_E5(imax,jmax,kmax)

      error=1
   
      open(unit=1,file=dirgraf//'f_'//fname,status='old',form='formatted',err=10)
      read(1,*) global_E1
      read(1,*) global_E2
      read(1,*) global_E5
      close(unit=1)

      write(*,*) 'READ FILE.: ','f_'//fname                
      error=0
   
      call flux_to_primitive(imax,global_E1,global_E2,global_E3,
     &  global_E4,global_E5,rho,u,v,w,Et,p,tp,mu)
      call equation_vorticity(stx_ex_6th_D,u,v,w,wx,wy,wz)

   10 continue
      if (error.ne.0) then
        write(*,*) 'ARQUIVO', ' f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
ccccc read_BIN_2D_PP: Read data in binary file for the post process.
ccccc **********************************************************************
      SUBROUTINE read_BIN_2D_PP(fname,stx_ex_6th_D,u,v,wz,rho,p,tp,Et,error)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      integer i,j,k,error
      real*8 stx_ex_6th_D(imax,7),
     &  rho(imax,jmax,kmax),u (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),mu(imax,jmax,kmax),
     &  wx (imax,jmax,kmax),wy(imax,jmax,kmax),
     &  wz (imax,jmax,kmax),
     &  global_E1(imax,jmax,kmax),global_E2(imax,jmax,kmax),
     &  global_E3(imax,jmax,kmax),global_E4(imax,jmax,kmax),
     &  global_E5(imax,jmax,kmax)

      error=1
  
      open(unit=1,file=dirgraf//'f_'//fname,status='old',form='formatted',err=10)
      read(1,*) global_E1
      read(1,*) global_E2
      read(1,*) global_E3
      read(1,*) global_E5
      close(unit=1)

      write(*,*) 'READ FILE.: ','f_'//fname                
      error=0
	
      call flux_to_primitive(imax,global_E1,global_E2,global_E3,
     &  global_E4,global_E5,rho,u,v,w,Et,p,tp,mu)
      call equation_vorticity(stx_ex_6th_D,u,v,w,wx,wy,wz)
 
   10 continue
      if (error.ne.0) then
        write(*,*) 'ARQUIVO', ' f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
ccccc read_BIN_3D_PP: Read data in binary file for the post process.
ccccc **********************************************************************
      SUBROUTINE read_BIN_3D_PP(fname,stx_ex_6th_D,u,v,w,
     &  wx,wy,wz,Qi,rho,p,tp,Et,error)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      integer i,j,k,error
      real*8 stx_ex_6th_D(imax,7),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w  (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p  (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  Qi (imax,jmax,kmax),mu (imax,jmax,kmax),
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

      write(*,*) 'READ FILE.: ','f_'//fname                
      error=0

      call flux_to_primitive(imax,global_E1,global_E2,global_E3,
     &  global_E4,global_E5,rho,u,v,w,Et,p,tp,mu)
      call equation_vorticity(stx_ex_6th_D,u,v,w,wx,wy,wz)
      call calculate_iso(u,v,w,Qi)

   10 continue
      if (error.ne.0) then
        write(*,*) 'ARQUIVO', ' f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
ccccc write_TCP_1D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_TCP_1D(fname,x,u,rho,p,tp,Et)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,Debug,IORE,xy,NPts,NElm,TecIni,
     &  TecDat,TecZne,TecNod,TecFil,TecEnd,VIsDouble
      real*8 x(imax),xx(imax,jmax,kmax),
     &  rho(imax,jmax,kmax),u(imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p(imax,jmax,kmax),
     &  tp (imax,jmax,kmax) 

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsDouble = 1
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y u v wz rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsDouble)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  imax,jmax,kmax,'BLOCK'//NULLCHR,NULLCHR)

      xy = dble(imax)

      do i=1,imax
        xx(i,1,1) = x(i)
      end do

      IORE = TecDat(xy,xx, VIsDouble)
      IORE = TecDat(xy,u,  VIsDouble)
      IORE = TecDat(xy,rho,VIsDouble)
      IORE = TecDat(xy,p,  VIsDouble)
      IORE = TecDat(xy,tp, VIsDouble)    
      IORE = TecDat(xy,Et, VIsDouble)
      IORE = TecEnd()

      if (QTimes.GT.1) then
        write(*,*) 'WRITE FILE: ','f_'//fname
      end if
      
      RETURN
      END

ccccc **********************************************************************
ccccc write_TCP_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_TCP_2D(fname,x,y,u,v,wz,rho,p,tp,Et)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,j,Debug,IORE,xy,NPts,NElm,TecIni,
     &  TecDat,TecZne,TecNod,TecFil,TecEnd,VIsDouble
      real*8 x(imax),y(jmax),
     &  xx (imax,jmax,kmax),yy (imax,jmax,kmax),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  wz (imax,jmax,kmax)

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsDouble = 1
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y u v wz rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsDouble)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  imax,jmax,kmax,'BLOCK'//NULLCHR,NULLCHR)

      xy = dble((imax)*(jmax))

      do i=1,imax
        do j=1,jmax
          xx(i,j,1) = x(i)
          yy(i,j,1) = y(j)
        end do
      end do

      IORE = TecDat(xy,xx, VIsDouble)
      IORE = TecDat(xy,yy, VIsDouble)
      IORE = TecDat(xy,u,  VIsDouble)
      IORE = TecDat(xy,v,  VIsDouble)
      IORE = TecDat(xy,wz, VIsDouble)
      IORE = TecDat(xy,rho,VIsDouble)
      IORE = TecDat(xy,p,  VIsDouble)
      IORE = TecDat(xy,tp, VIsDouble)    
      IORE = TecDat(xy,Et, VIsDouble)
      IORE = TecEnd()

      if (QTimes.GT.1) then
        write(*,*) 'WRITE FILE: ','f_'//fname
      end if
      
      RETURN
      END

ccccc **********************************************************************
ccccc write_TCP_3D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_TCP_3D(fname,x,y,z,u,v,w,wx,wy,wz,Qi,rho,p,tp,Et)

      IMPLICIT NONE
      INCLUDE 'snparv36.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,j,k,Debug,IORE,xy,NPts,NElm,TecIni,
     &  TecDat,TecZne,TecNod,TecFil,TecEnd,VIsDouble
      real*8 x(imax),y(jmax),z(kmax),
     &  xx (imax,jmax,kmax),yy (imax,jmax,kmax),
     &  zz (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  w  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  wx (imax,jmax,kmax),wy (imax,jmax,kmax),
     &  wz (imax,jmax,kmax),Qi (imax,jmax,kmax)  

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsDouble = 1
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y z u v w wx wy wz Qi rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsDouble)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  imax,jmax,kmax,'BLOCK'//NULLCHR,NULLCHR)

      xy = dble((imax)*(jmax)*(kmax))

      do i=1,imax
        do j=1,jmax
          do k=1,kmax
            xx(i,j,k) = x(i)
            yy(i,j,k) = y(j)
            zz(i,j,k) = y(k)
          end do
        end do
      end do

      IORE = TecDat(xy,xx, VIsDouble)
      IORE = TecDat(xy,yy, VIsDouble)
      IORE = TecDat(xy,zz, VIsDouble)
      IORE = TecDat(xy,u,  VIsDouble)
      IORE = TecDat(xy,v,  VIsDouble)
      IORE = TecDat(xy,w,  VIsDouble)
      IORE = TecDat(xy,wx, VIsDouble)
      IORE = TecDat(xy,wy, VIsDouble)
      IORE = TecDat(xy,wz, VIsDouble)
      IORE = TecDat(xy,Qi, VIsDouble)
      IORE = TecDat(xy,rho,VIsDouble)
      IORE = TecDat(xy,p,  VIsDouble)
      IORE = TecDat(xy,tp, VIsDouble)    
      IORE = TecDat(xy,Et, VIsDouble)
      IORE = TecEnd()

      if (QTimes.GT.1) then
        write(*,*) 'WRITE FILE: ','f_'//fname
      end if
      
      RETURN
      END


