cccccc *******************************************************************
cccccc Doctorate Project USP/EESC: Numerical simulation of the transition
cccccc                             for turbulence in a compressible 
cccccc                             boundary layer on a flat plate.
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos
cccccc Programmer\Guiding: Marcello A. Faraco de Medeiros
cccccc Version : 1.1
cccccc *******************************************************************
cccccc This program convert binary data for TECPLOT file.
cccccc Comando para execução:
cccccc   ifort dat2plt.f -o dat2plt /usr/local/tecplot8/lib/tecio.a 
cccccc   /usr/local/tecplot8/lib/ctype.o
cccccc *******************************************************************
      MODULE sndat2plt

      contains
ccccc **********************************************************************
ccccc dat2plt: Main subroutine.
ccccc **********************************************************************
      SUBROUTINE dat2plt(t_ini,t_end)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*02 numstr
      character*05 fname
      character*07 dir
      character*20 fname_aux
      integer i,j,k,tt,t_ini,t_end,error,extra_yy
      real*8 x(imax),y(jmax),z(kmax),
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),U3(imax,jmax,kmax),
     &  U4(imax,jmax,kmax),U5(imax,jmax,kmax),
     &  y_ini,x_full(imax),y_full(extra_y*jmax),z_full(kmax),
     &  U1_full(imax,extra_y*jmax,kmax),U2_full(imax,extra_y*jmax,kmax),
     &  U3_full(imax,extra_y*jmax,kmax),U4_full(imax,extra_y*jmax,kmax),
     &  U5_full(imax,extra_y*jmax,kmax),
     &  rho(imax,extra_y*jmax,kmax),u  (imax,extra_y*jmax,kmax),
     &  v  (imax,extra_y*jmax,kmax),w  (imax,extra_y*jmax,kmax),
     &  Et (imax,extra_y*jmax,kmax),p  (imax,extra_y*jmax,kmax),
     &  tp (imax,extra_y*jmax,kmax),wx (imax,extra_y*jmax,kmax),
     &  wy (imax,extra_y*jmax,kmax),wz (imax,extra_y*jmax,kmax),
     &  Qi (imax,extra_y*jmax,kmax)
      common/par/u,v,w,wx,wy,wz,Qi,rho,p,tp,Et

      fname   = 'SN__S'

      write(*,*) ' _______________________________________________'
      write(*,*) 
      write(*,*) ' CONVERTENDO ARQUIVOS DE FORMATO BINÁRIO PARA   '   
      write(*,*) ' O FORMATO DE LEITURA PADRÃO TECPLOT...         '  
      write(*,*)      

      call init_grid(x,y,z)

      !MS$IF (tp_dim.EQ.2)
        call makelayout_2D(x,y,z,t_ini,t_end)
      !MS$ELSEIF (tp_dim.EQ.3)
        call makelayout_3D(x,y,z,t_ini,t_end)
      !MS$ENDIF

      do tt=t_ini,t_end
        if (mod(tt,qtimes).EQ.0) then
          call get_filename(tt,fname,'.dat',fname_aux)           

          y_ini = 0.d0 
          do extra_yy=1,extra_y

            call numtostr(extra_yy,numstr)
            dir = 'vis-'//numstr//'/'
 
            !MS$IF (tp_dim.EQ.1)
              call read_BIN_1D(fname_aux,U1,U2,U5,error)
            !MS$ELSEIF (tp_dim.EQ.2)
              call read_BIN_2D(fname_aux,U1,U2,U3,U5,error)
            !MS$ELSEIF (tp_dim.EQ.3)
              call read_BIN_3D(dir,fname_aux,U1,U2,U3,U4,U5,error)
            !MS$ENDIF
            
            do i=1,imax
              do j=1,jmax
                do k=1,kmax
                  x_full(i         ) = x(i)
                  y_full(extra_yy*j) = y_ini+y(j)
                  z_full(k         ) = z(k)
                  U1_full(i,extra_yy*j,k) = U1(i,j,k)
                  U2_full(i,extra_yy*j,k) = U2(i,j,k)
                  U3_full(i,extra_yy*j,k) = U3(i,j,k)
                  U4_full(i,extra_yy*j,k) = U4(i,j,k)
                  U5_full(i,extra_yy*j,k) = U5(i,j,k)
                end do
              end do
            end do      

            y_ini = y_full(extra_yy*jmax)
          end do 
         

          call flux_to_primitive(imax,extra_y*jmax,kmax,Ma,gamma,
     &      U1_full,U2_full,U3_full,U4_full,U5_full,rho,u,v,w,Et,p,tp)
cc        call equation_vorticity(u,v,w,wx,wy,wz)

          if ((jmax.EQ.1).AND.(kmax.EQ.1)) then
            if (error.EQ.0) then
              call get_filename(tt,fname,'.plt',fname_aux)
              call write_PLT_TCP_1D(fname_aux,x,u,rho,p,tp,Et)
            end if
          else if ((jmax.NE.1).AND.(kmax.EQ.1)) then
            if (error.EQ.0) then
              call get_filename(tt,fname,'.plt',fname_aux)
              call write_PLT_TCP_2D(fname_aux,x,y,u,v,wz,rho,p,tp,Et)
            end if
          else if ((jmax.NE.1).AND.(kmax.NE.1)) then
            if (error.EQ.0) then
cc              call calculate_iso(u,v,w,Qi)
              call get_filename(tt,fname,'.plt',fname_aux)
              call write_PLT_TCP_3D(fname_aux,x_full,y_full,z_full)
            end if
          end if
        end if
      end do

      RETURN
      END SUBROUTINE dat2plt

ccccc **********************************************************************
ccccc Numtostr: Convert number to string values.
ccccc **********************************************************************
      SUBROUTINE numtostr(num,str)

      IMPLICIT NONE

      character*02 str
      integer num

      if (num.EQ.1) then
        str='01'
      else if (num.EQ.2) then
        str='02'
      else if (num.EQ.3) then
        str='03'
      else if (num.EQ.4) then
        str='04'
      else if (num.EQ.5) then
        str='05'
      else if (num.EQ.6) then
        str='06'
      else if (num.EQ.7) then
        str='07'
      else if (num.EQ.8) then
        str='08'
      else if (num.EQ.9) then
        str='09'
      else if (num.EQ.10) then
        str='10'
      end if

      RETURN
      END SUBROUTINE numtostr

ccccc **********************************************************************
ccccc write_PLT_TCP_1D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_PLT_TCP_1D(fname,x,u,rho,p,tp,Et)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,Debug,IORE,xy,NPts,NElm,TecIni,
     &  TecDat,TecZne,TecNod,TecFil,TecEnd,VIsPrecis
      real*8 x(imax),
     &  u  (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),
     &  Et (imax,jmax,kmax) 
      real*4 xx(imax,jmax,kmax),
     &  u_4  (imax,jmax,kmax),rho_4(imax,jmax,kmax),
     &  p_4  (imax,jmax,kmax),tp_4 (imax,jmax,kmax),  
     &  Et_4 (imax,jmax,kmax)

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsPrecis = 0
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y u v wz rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  imax,jmax,kmax,'BLOCK'//NULLCHR,NULLCHR)

      xy = dble(imax)

      do i=1,imax
        xx   (i,1,1) = x  (i)
        u_4  (i,1,1) = u  (i,1,1)
        rho_4(i,1,1) = rho(i,1,1)
        p_4  (i,1,1) = p  (i,1,1)
        tp_4 (i,1,1) = tp (i,1,1)
        Et_4 (i,1,1) = Et (i,1,1)
      end do

      IORE = TecDat(xy,xx,   VIsPrecis)
      IORE = TecDat(xy,u_4,  VIsPrecis)
      IORE = TecDat(xy,rho_4,VIsPrecis)
      IORE = TecDat(xy,p_4,  VIsPrecis)
      IORE = TecDat(xy,tp_4, VIsPrecis)    
      IORE = TecDat(xy,Et_4, VIsPrecis)
      IORE = TecEnd()

      write(*,*) ' WRITE FILE: ',dirgraf//'f_'//fname
      
      RETURN
      END SUBROUTINE write_PLT_TCP_1D

ccccc **********************************************************************
ccccc write_PLT_TCP_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_PLT_TCP_2D(fname,x,y,u,v,wz,rho,p,tp,Et)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,j,Debug,IORE,xy,NPts,NElm,TecIni,
     &  TecDat,TecZne,TecNod,TecFil,TecEnd,VIsPrecis
      real*8 x(imax),y(jmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  wz (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  Et (imax,jmax,kmax)
      real*4 xx(imax,jmax,kmax),yy(imax,jmax,kmax),
     &  u_4  (imax,jmax,kmax),v_4  (imax,jmax,kmax),
     &  wz_4 (imax,jmax,kmax),rho_4(imax,jmax,kmax),
     &  p_4  (imax,jmax,kmax),tp_4 (imax,jmax,kmax),  
     &  Et_4 (imax,jmax,kmax)

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsPrecis = 0
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y u v wz rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  imax,jmax,kmax,'BLOCK'//NULLCHR,NULLCHR)

      xy = dble((imax)*(jmax))

      do i=1,imax
        do j=1,jmax
          xx   (i,j,1) = x  (i)
          yy   (i,j,1) = y  (j)
          u_4  (i,j,1) = u  (i,j,1)
          v_4  (i,j,1) = v  (i,j,1)
          wz_4 (i,j,1) = wz (i,j,1)
          rho_4(i,j,1) = rho(i,j,1)
          p_4  (i,j,1) = p  (i,j,1)
          tp_4 (i,j,1) = tp (i,j,1)
          Et_4 (i,j,1) = Et (i,j,1)
        end do
      end do

      IORE = TecDat(xy,xx   ,VIsPrecis)
      IORE = TecDat(xy,yy   ,VIsPrecis)
      IORE = TecDat(xy,u_4  ,VIsPrecis)
      IORE = TecDat(xy,v_4  ,VIsPrecis)
      IORE = TecDat(xy,wz_4 ,VIsPrecis)
      IORE = TecDat(xy,rho_4,VIsPrecis)
      IORE = TecDat(xy,p_4  ,VIsPrecis)
      IORE = TecDat(xy,tp_4 ,VIsPrecis)    
      IORE = TecDat(xy,Et_4 ,VIsPrecis)
      IORE = TecEnd()

      write(*,*) ' WRITE FILE: ',dirgraf//'f_'//fname
      
      RETURN
      END SUBROUTINE write_PLT_TCP_2D

ccccc **********************************************************************
ccccc write_PLT_TCP_3D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_PLT_TCP_3D(fname,x,y,z)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,j,k,ii,jj,kk,iimax,jjmax,kkmax,Debug,IORE,xyz,NPts,NElm,
     &  TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd,VIsPrecis
      real*8 x(imax),y(extra_y*jmax),z(kmax),
     &  u  (imax,extra_y*jmax,kmax),v  (imax,extra_y*jmax,kmax),
     &  w  (imax,extra_y*jmax,kmax),wx (imax,extra_y*jmax,kmax),
     &  wy (imax,extra_y*jmax,kmax),wz (imax,extra_y*jmax,kmax),
     &  Qi (imax,extra_y*jmax,kmax),rho(imax,extra_y*jmax,kmax),
     &  p  (imax,extra_y*jmax,kmax),tp (imax,extra_y*jmax,kmax),
     &  Et (imax,extra_y*jmax,kmax)  
      real*4 
     &  xx   (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  yy   (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  zz   (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  u_4  (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  v_4  (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  w_4  (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  wx_4 (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  wy_4 (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  wz_4 (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  Qi_4 (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  rho_4(imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  p_4  (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  tp_4 (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z),
     &  Et_4 (imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z)  
      common/par/u,v,w,wx,wy,wz,Qi,rho,p,tp,Et

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsPrecis = 0
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y z u v w wx wy wz Qi rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  imax/qtimes_x,extra_y*(jmax/qtimes_y),kmax/qtimes_z,'BLOCK'//NULLCHR,NULLCHR)

      xyz = dble((imax/qtimes_x)*extra_y*(jmax/qtimes_y)*(kmax/qtimes_z))

      if (qtimes_x.EQ.1) then
        iimax = imax
      else
        iimax = imax-1
      end if

      if (qtimes_y.EQ.1) then
        jjmax = extra_y*jmax
      else
        jjmax = jmax-1
      end if

      if (qtimes_z.EQ.1) then
        kkmax = kmax
      else
        kkmax = kmax-1
      end if

      do i=1,iimax,qtimes_x
        do j=1,jjmax,qtimes_y
          do k=1,kkmax,qtimes_z
            ii = (i+qtimes_x-1)/qtimes_x
            jj = (j+qtimes_y-1)/qtimes_y
            kk = (k+qtimes_z-1)/qtimes_z

            xx   (ii,jj,kk) = x  (i)
            yy   (ii,jj,kk) = y  (j)
            zz   (ii,jj,kk) = z  (k)
            u_4  (ii,jj,kk) = u  (i,j,k)
            v_4  (ii,jj,kk) = v  (i,j,k)
            w_4  (ii,jj,kk) = w  (i,j,k)
            wx_4 (ii,jj,kk) = wx (i,j,k)
            wy_4 (ii,jj,kk) = wy (i,j,k)
            wz_4 (ii,jj,kk) = wz (i,j,k)
            Qi_4 (ii,jj,kk) = Qi (i,j,k)
            rho_4(ii,jj,kk) = rho(i,j,k)
            p_4  (ii,jj,kk) = p  (i,j,k)
            tp_4 (ii,jj,kk) = tp (i,j,k)
            Et_4 (ii,jj,kk) = Et (i,j,k)
          end do
        end do
      end do

      IORE = TecDat(xyz,xx,   VIsPrecis)
cc      IORE = TecDat(xyz,yy,   VIsPrecis)
cc      IORE = TecDat(xyz,zz,   VIsPrecis)
cc      IORE = TecDat(xyz,u_4,  VIsPrecis)
cc      IORE = TecDat(xyz,v_4,  VIsPrecis)
cc      IORE = TecDat(xyz,w_4,  VIsPrecis)
cc      IORE = TecDat(xyz,wx_4, VIsPrecis)
cc      IORE = TecDat(xyz,wy_4, VIsPrecis)
cc      IORE = TecDat(xyz,wz_4, VIsPrecis)
cc      IORE = TecDat(xyz,Qi_4, VIsPrecis)
cc      IORE = TecDat(xyz,rho_4,VIsPrecis)
cc      IORE = TecDat(xyz,p_4,  VIsPrecis)
cc      IORE = TecDat(xyz,tp_4, VIsPrecis)    
cc      IORE = TecDat(xyz,Et_4, VIsPrecis)
      IORE = TecEnd()

      write(*,*) ' WRITE FILE: ',dirgraf//'f_'//fname
      
      RETURN
      END SUBROUTINE write_PLT_TCP_3D

ccccc **********************************************************************
ccccc makelayout_2D: Save TECPLOT layout for movies.
ccccc **********************************************************************
      SUBROUTINE makelayout_2D(x,y,z,t_ini,t_end)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname_aux
      integer tt,t_ini,t_end
      real*8 x(imax),y(jmax),z(kmax)

      open (1,file=dirgraf//'lay_anime.lay',access='sequential')
      write(1,'(A08)     ') '#!MC 800'
      write(1,'(A21,A1,\)') '$!VarSet |LFDSFN1| = ',char(39)

      do tt=t_ini,t_end
        call get_filename(tt,'SN__S','.plt',fname_aux)           

        if (((mod(tt,qtimes).EQ.0)).AND.(tt.GT.0)) then
          write(1,'(A25,\)') '"f_'//fname_aux//'" '
        end if
      end do

      write(1,'(A1)') char(39)

      write(1,*) '$!VarSet |LFDSVL1| ='//char(39)
     &  //'"x" "y" "u" "v" "wz" "rho" "p" "T" "Et"'//char(39)
      write(1,*) '$!SETSTYLEBASE FACTORY'
      write(1,*) '$!PAPER'
      write(1,*) '  BACKGROUNDCOLOR = WHITE'
      write(1,*) '  ISTRANSPARENT = YES'
      write(1,*) '  ORIENTPORTRAIT = NO'
      write(1,*) '  SHOWGRID = YES'
      write(1,*) '  SHOWRULER = YES'
      write(1,*) '  SHOWPAPER = YES'
      write(1,*) '  PAPERSIZE = LETTER'
      write(1,*) '  PAPERSIZEINFO'
      write(1,*) '    {'
      write(1,*) '    LETTER'
      write(1,*) '      {'
      write(1,*) '      WIDTH = 8.5'
      write(1,*) '      HEIGHT = 11'
      write(1,*) '      LEFTHARDCLIPOFFSET = 0.125'
      write(1,*) '      RIGHTHARDCLIPOFFSET = 0.125'
      write(1,*) '      TOPHARDCLIPOFFSET = 0.125'
      write(1,*) '      BOTTOMHARDCLIPOFFSET = 0.125'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '  RULERSPACING = ONEINCH'
      write(1,*) '  PAPERGRIDSPACING = HALFINCH'
      write(1,*) '  REGIONINWORKAREA'
      write(1,*) '    {'
      write(1,*) '    X1 = -0.05'
      write(1,*) '    Y1 = -0.05'
      write(1,*) '    X2 = 11.05'
      write(1,*) '    Y2 = 8.55'
      write(1,*) '    }'
      write(1,*) '$!COLORSPECTRUM' 
      write(1,*) '  CONTOURCOLORMAP = SMRAINBOW'
      write(1,*) '  SURFACERGBMIN'
      write(1,*) '    {'
      write(1,*) '    R = 0'
      write(1,*) '    G = 0'
      write(1,*) '    B = 0'
      write(1,*) '    }'
      write(1,*) '  SURFACERGBMAX'
      write(1,*) '    {'
      write(1,*) '    R = 255'
      write(1,*) '    G = 255'
      write(1,*) '    B = 255'
      write(1,*) '    }'
      write(1,*) '$!COLORMAPCONTROL RESETTOFACTORY'
      write(1,*) '### Frame Number 1 ###'
      write(1,*) '$!READDATASET  '//char(39)//'|LFDSFN1|'//char(39) 
      write(1,*) '  INCLUDETEXT = NO'
      write(1,*) '  INCLUDEGEOM = NO'
      write(1,*) '  VARLOADMODE = BYNAME'
      write(1,*) '  VARNAMELIST ='//char(39)//'|LFDSVL1|'//char(39) 
      write(1,*) '$!REMOVEVAR |LFDSVL1|'
      write(1,*) '$!FRAMELAYOUT'
      write(1,*) '  SHOWHEADER = NO'
      write(1,*) '  HEADERCOLOR = RED'
      write(1,*) '  XYPOS'
      write(1,*) '    {'
      write(1,*) '    X = 0.039247'
      write(1,*) '    Y = 0.082739'
      write(1,*) '    }'
      write(1,*) '  WIDTH = 10.88'
      write(1,*) '  HEIGHT = 8.3763'
      write(1,*) '$!FRAMEMODE  = TWOD'
      write(1,*) '$!FRAMENAME  = '//char(39)//'frame 001'//char(39) 
      write(1,*) '$!ACTIVEFIELDZONES  =  [1]'
      write(1,*) '$!GLOBALCONTOUR'
      write(1,*) '  VAR = 5'
      write(1,*) '$!CONTOURLEVELS NEW'
      write(1,*) '  RAWDATA'
      write(1,*) '30'
      write(1,*) '-0.0256113478781'
      write(1,*) '-0.00845669437187'
      write(1,*) '0.00869795913439'
      write(1,*) '0.0258526126406'
      write(1,*) '0.0430072661469'
      write(1,*) '0.0601619196532'
      write(1,*) '0.0773165731594'
      write(1,*) '0.0944712266657'
      write(1,*) '0.111625880172'
      write(1,*) '0.128780533678'
      write(1,*) '0.145935187184'
      write(1,*) '0.163089840691'
      write(1,*) '0.180244494197'
      write(1,*) '0.197399147703'
      write(1,*) '0.214553801209'
      write(1,*) '0.231708454716'
      write(1,*) '0.248863108222'
      write(1,*) '0.266017761728'
      write(1,*) '0.283172415234'
      write(1,*) '0.300327068741'
      write(1,*) '0.317481722247'
      write(1,*) '0.334636375753'
      write(1,*) '0.35179102926'
      write(1,*) '0.368945682766'
      write(1,*) '0.386100336272'
      write(1,*) '0.403254989778'
      write(1,*) '0.420409643285'
      write(1,*) '0.437564296791'
      write(1,*) '0.454718950297'
      write(1,*) '0.471873603803'
      do tt=1,(t_end-t_ini)/Qtimes
        write(1,'(a10,i6,a1)') '$!FIELD  [',tt,']'
        write(1,*) '  MESH'
        write(1,*) '    {'
        write(1,*) '    COLOR = RED'
        write(1,*) '    }'
        write(1,*) '  CONTOUR'
        write(1,*) '    {'
        write(1,*) '    CONTOURTYPE = BOTHLINESANDFLOOD'
        write(1,*) '    COLOR = RED'
        write(1,*) '    }'
        write(1,*) '  VECTOR'
        write(1,*) '    {'
        write(1,*) '    COLOR = RED'
        write(1,*) '    }'
        write(1,*) '  SCATTER'
        write(1,*) '    {'
        write(1,*) '    COLOR = RED'
        write(1,*) '    }'
        write(1,*) '  SHADE'
        write(1,*) '    {'
        write(1,*) '    COLOR = RED'
        write(1,*) '    }'
        write(1,*) '  BOUNDARY'
        write(1,*) '    {'
        write(1,*) '    COLOR = RED'
        write(1,*) '    }'
      end do
      write(1,*) '$!TWODAXIS '
      write(1,*) '  XVAR = 1'
      write(1,*) '  YVAR = 2'
      write(1,*) '$!VIEW FIT'
      write(1,*) '$!TWODAXIS'
      write(1,*) '  GRIDAREA'
      write(1,*) '    {'
      write(1,*) '    EXTENTS'
      write(1,*) '      {'
      write(1,*) '      X1 = 7.4965'
      write(1,*) '      Y1 = 7.672'
      write(1,*) '      X2 = 96.447'
      write(1,*) '      Y2 = 96.487'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '  DEPXTOYRATIO = 1'
      write(1,*) '  AXISMODE = INDEPENDENT'
      write(1,*) '  XDETAIL'
      write(1,*) '    {'
      write(1,*) '    AXISPOSITION = 7.672'
      write(1,*) '    RANGEMIN = ',x(1)
      write(1,*) '    RANGEMAX = ',x(imax)
      write(1,*) '    GRSPACING = 50'
      write(1,*) '    }'
      write(1,*) '  YDETAIL'
      write(1,*) '    {'
      write(1,*) '    AXISPOSITION = 7.497'
      write(1,*) '    RANGEMIN = ',y(1)
      write(1,*) '    RANGEMAX = ',y(jmax)
      write(1,*) '    GRSPACING = 2'
      write(1,*) '    TITLE'
      write(1,*) '      {'
      write(1,*) '      OFFSET = 7.024'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '$!FIELDLAYERS'
      write(1,*) '  SHOWMESH = NO'
      write(1,*) '  SHOWCONTOUR = YES'
      write(1,*) '$!REMOVEVAR |LFDSFN1|'
      write(1,*) '$!SETSTYLEBASE CONFIG'

      close (unit=1) 

      RETURN
      END SUBROUTINE makelayout_2D

ccccc **********************************************************************
ccccc makelayout_3D: Save TECPLOT layout for movies.
ccccc **********************************************************************
      SUBROUTINE makelayout_3D(x,y,z,t_ini,t_end)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname_aux
      integer tt,t_ini,t_end
      real*8 x(imax),y(jmax),z(kmax)

      open (1,file=dirgraf//'lay_anime.lay',access='sequential')

      write(1,'(A9)') '#!MC 1000'
      write(1,*) '$!VarSet |LFBD| = '//char(39)//'/mnt/dados/usp/temp/version37/visual'//char(39)
      write(1,'(A22,A1,\)') '$!VarSet |LFDSFN1| = ',char(39)
      do tt=t_ini,t_end
        call get_filename(tt,'SN__S','.plt',fname_aux)           

        if (((mod(tt,qtimes).EQ.0)).AND.(tt.GT.0)) then
          write(1,'(A33,\)') '"|LFBD|/f_'//fname_aux//'"'
        end if
      end do
      write(1,'(A1)') char(39)

      write(1,*) '$!VarSet |LFDSVL1| = '//char(39)//'"x""y""z""u""v""w""wx""wy""wz""Qi""rho""p""T""Et"'//char(39)
      write(1,*) '$!SETSTYLEBASE FACTORY'
      write(1,*) '$!PAPER'
      write(1,*) '  BACKGROUNDCOLOR = WHITE'
      write(1,*) '  ISTRANSPARENT = YES'
      write(1,*) '  ORIENTPORTRAIT = NO'
      write(1,*) '  SHOWGRID = YES'
      write(1,*) '  SHOWRULER = YES'
      write(1,*) '  SHOWPAPER = YES'
      write(1,*) '  PAPERSIZE = LETTER'
      write(1,*) '  PAPERSIZEINFO'
      write(1,*) '    {'
      write(1,*) '    LETTER'
      write(1,*) '      {'
      write(1,*) '      WIDTH = 8.5'
      write(1,*) '      HEIGHT = 11'
      write(1,*) '      LEFTHARDCLIPOFFSET = 0.125'
      write(1,*) '      RIGHTHARDCLIPOFFSET = 0.125'
      write(1,*) '      TOPHARDCLIPOFFSET = 0.125'
      write(1,*) '      BOTTOMHARDCLIPOFFSET = 0.125'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '  RULERSPACING = ONEINCH'
      write(1,*) '  PAPERGRIDSPACING = HALFINCH'
      write(1,*) '  REGIONINWORKAREA'
      write(1,*) '    {'
      write(1,*) '    X1 = -0.05'
      write(1,*) '    Y1 = -0.05'
      write(1,*) '    X2 = 11.05'
      write(1,*) '    Y2 = 8.55'
      write(1,*) '    }'
      write(1,*) '$!COLORMAP'
      write(1,*) '  CONTOURCOLORMAP = SMRAINBOW'
      write(1,*) '$!COLORMAPCONTROL RESETTOFACTORY'
      write(1,*) '### Frame Number 1 ###'
      write(1,*) '$!READDATASET  ',char(39),'|LFDSFN1|',char(39) 
      write(1,*) '  INITIALPLOTTYPE = CARTESIAN3D'
      write(1,*) '  INCLUDETEXT = NO'
      write(1,*) '  INCLUDEGEOM = NO'
      write(1,*) '  VARLOADMODE = BYNAME'
      write(1,*) '  VARNAMELIST = ',char(39),'|LFDSVL1|',char(39) 
      write(1,*) '$!REMOVEVAR |LFDSVL1|'
      write(1,*) '$!REMOVEVAR |LFDSFN1|'
      write(1,*) '$!FRAMELAYOUT '
      write(1,*) '  SHOWHEADER = NO'
      write(1,*) '  HEADERCOLOR = RED'
      write(1,*) '  XYPOS'
      write(1,*) '    {'
      write(1,*) '    X = 0.077156'
      write(1,*) '    Y = 0.034265'
      write(1,*) '    }'
      write(1,*) '  WIDTH = 10.826'
      write(1,*) '  HEIGHT = 8.4021'
      write(1,*) '$!PLOTTYPE  = CARTESIAN3D'
      write(1,*) '$!FRAMENAME  = ',char(39),'Frame 001',char(39) 
      write(1,*) '$!ACTIVEFIELDZONES  =  [1]'
      write(1,*) '$!GLOBALRGB'
      write(1,*) '  RANGEMIN = 0'
      write(1,*) '  RANGEMAX = 1'
      write(1,*) '$!GLOBALCONTOUR  1'
      write(1,*) '  VAR = 12'
      write(1,*) '  DEFNUMLEVELS = 2'
      write(1,*) '  LEGEND'
      write(1,*) '    {'
      write(1,*) '    XYPOS'
      write(1,*) '      {'
      write(1,*) '      X = 95'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '  COLORCUTOFF'
      write(1,*) '    {'
      write(1,*) '    RANGEMIN = 4.46388232355'
      write(1,*) '    RANGEMAX = 4.46468895257'
      write(1,*) '    }'
      write(1,*) '  COLORMAPFILTER'
      write(1,*) '    {'
      write(1,*) '    CONTINUOUSCOLOR'
      write(1,*) '      {'
      write(1,*) '      CMIN = 4.46347900903'
      write(1,*) '      CMAX = 4.46509226709'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '$!CONTOURLEVELS NEW'
      write(1,*) '  CONTOURGROUP = 1'
      write(1,*) '  RAWDATA'
      write(1,*) '1'
      write(1,*) '4.5'
      write(1,*) '$!GLOBALCONTOUR  2'
      write(1,*) '  VAR = 10'
      write(1,*) '  LEGEND'
      write(1,*) '    {'
      write(1,*) '    XYPOS'
      write(1,*) '      {'
      write(1,*) '      X = 95'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '  COLORCUTOFF'
      write(1,*) '    {'
      write(1,*) '    RANGEMIN = -0.680775590241'
      write(1,*) '    RANGEMAX = -0.0118814334273'
      write(1,*) '    }'
      write(1,*) '  COLORMAPFILTER'
      write(1,*) '    {'
      write(1,*) '    CONTINUOUSCOLOR'
      write(1,*) '      {'
      write(1,*) '      CMIN = -1.01522266865'
      write(1,*) '      CMAX = 0.322565644979'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '$!GLOBALCONTOUR  3'
      write(1,*) '  VAR = 11'
      write(1,*) '  LEGEND'
      write(1,*) '    {'
      write(1,*) '    XYPOS'
      write(1,*) '      {'
      write(1,*) '      X = 95'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '  COLORCUTOFF'
      write(1,*) '    {'
      write(1,*) '    RANGEMIN = 0.946419730783'
      write(1,*) '    RANGEMAX = 1.01876227558'
      write(1,*) '    }'
      write(1,*) '  COLORMAPFILTER'
      write(1,*) '    {'
      write(1,*) '    CONTINUOUSCOLOR'
      write(1,*) '      {'
      write(1,*) '      CMIN = 0.910248458385'
      write(1,*) '      CMAX = 1.05493354797'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '$!GLOBALCONTOUR  4'
      write(1,*) '  VAR = 13'
      write(1,*) '  LEGEND'
      write(1,*) '    {'
      write(1,*) '    XYPOS'
      write(1,*) '      {'
      write(1,*) '      X = 95'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '  COLORCUTOFF'
      write(1,*) '    {'
      write(1,*) '    RANGEMIN = 0.977718070149'
      write(1,*) '    RANGEMAX = 1.00698707998'
      write(1,*) '    }'
      write(1,*) '  COLORMAPFILTER'
      write(1,*) '    {'
      write(1,*) '    CONTINUOUSCOLOR'
      write(1,*) '      {'
      write(1,*) '      CMIN = 0.963083565235'
      write(1,*) '      CMAX = 1.02162158489'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '$!GLOBALSCATTER'
      write(1,*) '  LEGEND'
      write(1,*) '    {'
      write(1,*) '    XYPOS'
      write(1,*) '      {'
      write(1,*) '      X = 95'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '  REFSCATSYMBOL'
      write(1,*) '    {'
      write(1,*) '    COLOR = RED'
      write(1,*) '    FILLCOLOR = RED'
      write(1,*) '    }'
      write(1,*) '$!FIELD  [1]'
      write(1,*) '  MESH'
      write(1,*) '    {'
      write(1,*) '    COLOR = RED'
      write(1,*) '    }'
      write(1,*) '  CONTOUR'
      write(1,*) '    {'
      write(1,*) '    CONTOURTYPE = FLOOD'
      write(1,*) '    COLOR = RED'
      write(1,*) '    USELIGHTINGEFFECT = YES'
      write(1,*) '    }'
      write(1,*) '  VECTOR'
      write(1,*) '    {'
      write(1,*) '    COLOR = RED'
      write(1,*) '    }'
      write(1,*) '  SCATTER'
      write(1,*) '    {'
      write(1,*) '    COLOR = RED'
      write(1,*) '    }'
      write(1,*) '  SHADE'
      write(1,*) '    {'
      write(1,*) '    COLOR = WHITE'
      write(1,*) '    }'
      write(1,*) '  BOUNDARY'
      write(1,*) '    {'
      write(1,*) '    SHOW = YES'
      write(1,*) '    COLOR = RED'
      write(1,*) '    }'
      write(1,*) '  POINTS'
      write(1,*) '    {'
      write(1,*) '    POINTSTOPLOT = SURFACENODES'
      write(1,*) '    }'
      write(1,*) '  SURFACES'
      write(1,*) '    {'
      write(1,*) '    SURFACESTOPLOT = EXPOSEDCELLFACES'
      write(1,*) '    }'
      write(1,*) '$!THREEDAXIS'
      write(1,*) '  XDETAIL'
      write(1,*) '    {'
      write(1,*) '    VARNUM = 1'
      write(1,*) '    }'
      write(1,*) '  YDETAIL'
      write(1,*) '    {'
      write(1,*) '    VARNUM = 2'
      write(1,*) '    }'
      write(1,*) '  ZDETAIL'
      write(1,*) '    {'
      write(1,*) '    VARNUM = 3'
      write(1,*) '    }'
      write(1,*) '$!VIEW FIT'
      write(1,*) '$!THREEDAXIS'
      write(1,*) '  AXISMODE = XYZDEPENDENT'
      write(1,*) '  XYDEPXTOYRATIO = 1'
      write(1,*) '  DEPXTOYRATIO = 1'
      write(1,*) '  DEPXTOZRATIO = 1'
      write(1,*) '$!THREEDAXIS'
      write(1,*) '  XDETAIL'
      write(1,*) '    {'
      write(1,*) '    SHOWAXIS = NO'
      write(1,*) '    RANGEMIN = -0.39270000000000004903'
      write(1,*) '    RANGEMAX = 8.2467000000000005855'
      write(1,*) '    GRSPACING = 2'
      write(1,*) '    AXISLINE'
      write(1,*) '      {'
      write(1,*) '      EDGE = 3'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '$!THREEDAXIS'
      write(1,*) '  YDETAIL'
      write(1,*) '    {'
      write(1,*) '    SHOWAXIS = NO'
      write(1,*) '    RANGEMIN = -5.5'
      write(1,*) '    RANGEMAX = 5.5'
      write(1,*) '    GRSPACING = 2'
      write(1,*) '    AXISLINE'
      write(1,*) '      {'
      write(1,*) '      EDGE = 2'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '$!THREEDAXIS' 
      write(1,*) '  ZDETAIL'
      write(1,*) '    {'
      write(1,*) '    SHOWAXIS = NO'
      write(1,*) '    RANGEMIN = -0.39270000000000004903'
      write(1,*) '    RANGEMAX = 8.2467000000000005855'
      write(1,*) '    GRSPACING = 2'
      write(1,*) '    AXISLINE'
      write(1,*) '      {'
      write(1,*) '      EDGE = 2'
      write(1,*) '      }'
      write(1,*) '    }'
      write(1,*) '$!GLOBALISOSURFACE' 
      write(1,*) '  SHOW = YES'
      write(1,*) '  ISOSURFACESELECTION = ALLCONTOURLEVELS'
      write(1,*) '  ISOVALUE1 = 4.46388232355'
      write(1,*) '  ISOVALUE2 = 4.46428563806'
      write(1,*) '  ISOVALUE3 = 4.46468895257'
      write(1,*) '  MARCHINGCUBEALGORITHM = CLASSICPLUS'
      write(1,*) '$!GLOBALSLICE'
      write(1,*) '  BOUNDARY'
      write(1,*) '    {'
      write(1,*) '    SHOW = NO'
      write(1,*) '    }'
      write(1,*) '$!GLOBALTHREED'
      write(1,*) '  AXISSCALEFACT'
      write(1,*) '    {'
      write(1,*) '    X = 1'
      write(1,*) '    Y = 1'
      write(1,*) '    Z = 1'
      write(1,*) '    }'
      write(1,*) '  ROTATEORIGIN'
      write(1,*) '    {'
      write(1,*) '    X = 3.927'
      write(1,*) '    Y = 0'
      write(1,*) '    Z = 3.927'
      write(1,*) '    }'
      write(1,*) '  LIGHTSOURCE'
      write(1,*) '    {'
      write(1,*) '    INTENSITY = 75'
      write(1,*) '    BACKGROUNDLIGHT = 30'
      write(1,*) '    }'
      write(1,*) '  LINELIFTFRACTION = 0.2'
      write(1,*) '  SYMBOLLIFTFRACTION = 0.6'
      write(1,*) '  VECTORLIFTFRACTION = 0.7'
      write(1,*) '$!THREEDVIEW'
      write(1,*) '  PSIANGLE = 42.8984'
      write(1,*) '  THETAANGLE = -155.267'
      write(1,*) '  ALPHAANGLE = 155.684'
      write(1,*) '  VIEWERPOSITION'
      write(1,*) '    {'
      write(1,*) '    X = 92.461889751'
      write(1,*) '    Y = 165.027886925'
      write(1,*) '    Z = 205.401564678'
      write(1,*) '    }'
      write(1,*) '  VIEWWIDTH = 39.3367'
      write(1,*) '$!FIELDLAYERS'
      write(1,*) '  SHOWMESH = NO'
      write(1,*) '  SHOWBOUNDARY = NO'
      write(1,*) '$!REMOVEVAR |LFBD|'
      write(1,*) '$!SETSTYLEBASE CONFIG'

      close (unit=1) 

      RETURN
      END SUBROUTINE makelayout_3D

      END MODULE sndat2plt

