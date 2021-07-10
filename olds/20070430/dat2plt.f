cccccc *******************************************************************
cccccc Doctorate Project USP/EESC: Numerical simulation of the transition
cccccc                             for turbulence in a compressible 
cccccc                             boundary layer on a plane plate.
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos
cccccc Programmer\Guiding: Marcello A. Faraco de Medeiros
cccccc Version : 3.6
cccccc *******************************************************************
cccccc This program convert binary data for binary data to TECPLOT.
cccccc Comando para execução:
cccccc   ifort dat2plt.f -o dat2plt /usr/local/tecplot8/lib/tecio.a 
cccccc   /usr/local/tecplot8/lib/ctype.o
cccccc *******************************************************************
      PROGRAM dat2plt

      IMPLICIT NONE
      INCLUDE '../snparv36.f90'

      character*20 fname_aux
      character*05 fname
      character*04 ext
      integer tt,q_times,t_ini,t_end,error
      real*8 x(imax),y(jmax),z(kmax),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w  (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p  (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  Qi (imax,jmax,kmax)

      fname   = 'SN__S'
      q_times = 0     
      t_ini   = 0
      t_end   = 0

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

          else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
            call read_BIN_2D(fname_aux,x,y,u,v,wz,rho,p,tp,Et,error)

            if (error.EQ.0) then
              call get_filename(tt,fname,'.plt',fname_aux)
              call write_TCP_2D(fname_aux,x,y,u,v,wz,rho,p,tp,Et)
            end if
          else if ((jmax.NE.1).AND.(kmax.NE.1)) then
            call read_BIN_3D(fname_aux,x,y,z,u,v,w,wx,wy,wz,
     &        Qi,rho,p,tp,Et,error)

            if (error.EQ.0) then
              call get_filename(tt,fname,'.plt',fname_aux)
              call write_TCP_3D(fname_aux,x,y,z,u,v,w,
     &          wx,wy,wz,Qi,rho,p,tp,Et)
            end if
          end if
        end if
      end do

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
ccccc read_BIN_2D: Read data in binary file.
ccccc **********************************************************************
      SUBROUTINE read_BIN_2D(fname,x,y,u,v,wz,rho,p,tp,Et,error)

      IMPLICIT NONE
      INCLUDE '../snparv36.f90'

      character*20 fname
      integer i,j,k,error
      real*8 x(imax),y(jmax),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w  (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p  (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wz (imax,jmax,kmax)

      error=1
   
      open(unit=1,file='f_'//fname,status='old',
     &  form='formatted',err=10)

      read(1,*) x
      read(1,*) y
      read(1,*) u
      read(1,*) v
      read(1,*) wz
      read(1,*) rho
      read(1,*) p
      read(1,*) tp
      read(1,*) Et
      close(unit=1)

      write(*,*) 'READ FILE.: ','f_'//fname                
      error=0

   10 continue
      if (error.ne.0) then
        write(*,*) 'ARQUIVO', ' f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
ccccc read_BIN_3D: Read data in binary file.
ccccc **********************************************************************
      SUBROUTINE read_BIN_3D(fname,x,y,z,u,v,w,wx,wy,wz,
     &  Qi,rho,p,tp,Et,error)

      IMPLICIT NONE
      INCLUDE '../snparv36.f90'

      character*20 fname
      integer i,j,k,error
      real*8 x(imax),y(jmax),z(kmax),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w  (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p  (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  Qi (imax,jmax,kmax)

      error=1
   
      open(unit=1,file='f_'//fname,status='old',
     &  form='formatted',err=10)

      read(1,*) x
      read(1,*) y
      read(1,*) z
      read(1,*) u
      read(1,*) v
      read(1,*) w
      read(1,*) wx
      read(1,*) wy
      read(1,*) wz
      read(1,*) Qi
      read(1,*) rho
      read(1,*) p
      read(1,*) tp
      read(1,*) Et
      close(unit=1)

      write(*,*) 'READ FILE.: ','f_'//fname                
      error=0

   10 continue
      if (error.ne.0) then
        write(*,*) 'ARQUIVO', ' f_'//fname, ' INEXISTENTE!!!' 
      end if
 
      RETURN
      END

ccccc **********************************************************************
ccccc write_TCP_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_TCP_2D(fname,x,y,u,v,wz,rho,p,tp,Et)

      IMPLICIT NONE
      INCLUDE '../snparv36.f90'

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
     &  'f_'//fname//NULLCHR,
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
      INCLUDE '../snparv36.f90'

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
     &  'f_'//fname//NULLCHR,
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


