cccccc *******************************************************************
cccccc Doctorate Project USP/EESC: Numerical simulation of the transition
cccccc                             for turbulence in a compressible 
cccccc                             boundary layer on a flat plate.
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos
cccccc Programmer\Guiding: Marcello A. Faraco de Medeiros
cccccc Version : 1.2
cccccc Update  : 18/09/2007
cccccc *******************************************************************
cccccc This program convert binary data for TECPLOT file.
cccccc Comando para execução:
cccccc   ifort dat2plt.f -o dat2plt /usr/local/tecplot8/lib/tecio.a 
cccccc   /usr/local/tecplot8/lib/ctype.o
cccccc *******************************************************************
      MODULE nsdat2plt

      contains
ccccc **********************************************************************
ccccc dat2plt: Main subroutine.
ccccc **********************************************************************
      SUBROUTINE dat2plt(t_ini,t_end)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*02 numstr
      character*07 dir
      character*20 fname_aux
      integer i,j,k,tt,t_ini,t_end,error,extra_xx,extra_yy,extra_zz,cont
      real*8 x(imax),y(jmax),z(kmax),x_adj,y_adj,z_adj,
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),U3(imax,jmax,kmax),
     &  U4(imax,jmax,kmax),U5(imax,jmax,kmax),

     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w  (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p  (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  Qi (imax,jmax,kmax),

     &  x_full(extra_x*imax),y_full(extra_y*jmax),z_full(extra_z*kmax),

     &  rho_full(extra_x*imax,extra_y*jmax,extra_z*kmax),u_full  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  v_full  (extra_x*imax,extra_y*jmax,extra_z*kmax),w_full  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  Et_full (extra_x*imax,extra_y*jmax,extra_z*kmax),p_full  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  tp_full (extra_x*imax,extra_y*jmax,extra_z*kmax),wx_full (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  wy_full (extra_x*imax,extra_y*jmax,extra_z*kmax),wz_full (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  Qi_full (extra_x*imax,extra_y*jmax,extra_z*kmax)

      common/par/u_full,v_full,w_full,wx_full,wy_full,wz_full,
     &  Qi_full,rho_full,p_full,tp_full,Et_full

      write(*,*) ' _______________________________________________'
      write(*,*) 
      write(*,*) ' CONVERTENDO ARQUIVOS DE FORMATO BINÁRIO PARA   '   
      write(*,*) ' O FORMATO DE LEITURA PADRÃO TECPLOT...         '  
      write(*,*)      

      do tt=t_ini,t_end
        if (mod(tt,qtimes).EQ.0) then
          call get_filename(tt,fnamegraf,'.dat',fname_aux)           

          x_adj = 0
          y_adj = y(1)
          z_adj = 0
          cont  = 0
          do extra_zz=1,extra_z
            do extra_yy=1,extra_y    
              do extra_xx=1,extra_x              
                if ((extra_x.EQ.1).AND.(extra_y.EQ.1).AND.(extra_z.EQ.1)) then
                  dir = dirgraf 
                else 
                  cont = cont+1
                  call numtostr(cont,numstr)
                  dir = 'vis-'//numstr//'/'
                end if
 
                call init_grid(dir,x,y,z)

                !MS$IF (tp_dim.EQ.1)
                  call read_BIN_1D(dir,fname_aux,U1,U2,U5,error)
                !MS$ELSEIF (tp_dim.EQ.2)
                  call read_BIN_2D(dir,fname_aux,U1,U2,U3,U5,error)
                !MS$ELSEIF (tp_dim.EQ.3)
                  call read_BIN_3D(dir,fname_aux,U1,U2,U3,U4,U5,error)
                !MS$ENDIF

                if (error.EQ.0) then
                  call flux_to_primitive(U1,U2,U3,U4,U5,rho,u,v,w,Et,p,tp)
                  call equation_vorticity(u,v,w,wx,wy,wz)

                  do i=1,imax
                    do j=1,jmax
                      do k=1,kmax
                        x_full  ((extra_xx-1)*imax+i) = x(i)-x(1)+x_adj
                        y_full  ((extra_yy-1)*jmax+j) = y(j)-y(1)+y_adj
                        z_full  ((extra_zz-1)*kmax+k) = z(k)-z(1)+z_adj
                        rho_full((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = rho(i,j,k)
                        u_full  ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = u  (i,j,k)
                        v_full  ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = v  (i,j,k)
                        w_full  ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = w  (i,j,k)
                        Et_full ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = Et (i,j,k)
                        p_full  ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = p  (i,j,k)
                        tp_full ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = tp (i,j,k)
                        wx_full ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = wx (i,j,k)
                        wy_full ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = wy (i,j,k)
                        wz_full ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = wz (i,j,k)
                        Qi_full ((extra_xx-1)*imax+i,(extra_yy-1)*jmax+j,(extra_zz-1)*kmax+k) = Qi (i,j,k)
                      end do
                    end do
                  end do         
 
                  x_adj = x_full(imax)+0.d0
                  y_adj = y_full(jmax)+0.d0
                  z_adj = z_full(kmax)+0.d0
                end if
              end do 
            end do
          end do  

          if (error.EQ.0) then
            call get_filename(tt,fnamegraf,'.plt',fname_aux)
            !MS$IF (tp_dim.EQ.1)    
              call write_PLT_TCP_1D(fname_aux,x_full)
            !MS$ELSEIF (tp_dim.EQ.2)
              call write_PLT_TCP_2D(fname_aux,x_full,y_full)
            !MS$ELSEIF (tp_dim.EQ.3)
              call calculate_iso(u,v,w,Qi)
              call write_PLT_TCP_3D(fname_aux,x_full,y_full,z_full)
            !MS$ENDIF
          end if 
        end if
      end do

      call makelayout(t_ini,t_end)

      RETURN
      END SUBROUTINE dat2plt

ccccc **********************************************************************
ccccc write_PLT_TCP_1D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_PLT_TCP_1D(fname,x)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,Debug,IORE,xy,NPts,NElm,TecIni,
     &  TecDat,TecZne,TecNod,TecFil,TecEnd,VIsPrecis
      real*8 x(extra_x*imax),
     &  u  (extra_x*imax,extra_y*jmax,extra_z*kmax),v  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  w  (extra_x*imax,extra_y*jmax,extra_z*kmax),wx (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  wy (extra_x*imax,extra_y*jmax,extra_z*kmax),wz (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  Qi (extra_x*imax,extra_y*jmax,extra_z*kmax),rho(extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  p  (extra_x*imax,extra_y*jmax,extra_z*kmax),tp (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  Et (extra_x*imax,extra_y*jmax,extra_z*kmax)  
      !MS$IF (precision_plot.EQ.1)
        real*4 
     &    xx   (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    u_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    rho_4(extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    p_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    tp_4 (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    Et_4 (extra_x*imax,extra_y*jmax,extra_z*kmax)  
      !MS$ELSEIF (precision_plot.EQ.2)
        real*8 
     &    xx   (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    u_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    rho_4(extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    p_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    tp_4 (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    Et_4 (extra_x*imax,extra_y*jmax,extra_z*kmax)  
      !MS$ENDIF 
      common/par/u,v,w,wx,wy,wz,Qi,rho,p,tp,Et

      NULLCHR   = CHAR(0)
      Debug     = 0
      !MS$IF (precision_plot.EQ.1)
        VIsPrecis = 0
      !MS$ELSEIF (precision_plot.EQ.2)
        VIsPrecis = 1
      !MS$ENDIF 
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y u v wz rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  extra_x*imax,jmax,kmax,'BLOCK'//NULLCHR,NULLCHR)

      xy = dble(extra_x*imax)

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

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
      
      RETURN
      END SUBROUTINE write_PLT_TCP_1D

ccccc **********************************************************************
ccccc write_PLT_TCP_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_PLT_TCP_2D(fname,x,y)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,j,Debug,IORE,xy,NPts,NElm,TecIni,
     &  TecDat,TecZne,TecNod,TecFil,TecEnd,VIsPrecis
      real*8 x(extra_x*imax),y(extra_y*jmax),
     &  u  (extra_x*imax,extra_y*jmax,extra_z*kmax),v  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  w  (extra_x*imax,extra_y*jmax,extra_z*kmax),wx (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  wy (extra_x*imax,extra_y*jmax,extra_z*kmax),wz (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  Qi (extra_x*imax,extra_y*jmax,extra_z*kmax),rho(extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  p  (extra_x*imax,extra_y*jmax,extra_z*kmax),tp (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  Et (extra_x*imax,extra_y*jmax,extra_z*kmax)  
      !MS$IF (precision_plot.EQ.1)
        real*4 
     &    xx   (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    yy   (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    u_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    v_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    wz_4 (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    rho_4(extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    p_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    tp_4 (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    Et_4 (extra_x*imax,extra_y*jmax,extra_z*kmax)  
      !MS$ELSEIF (precision_plot.EQ.2)
        real*8 
     &    xx   (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    yy   (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    u_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    v_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    wz_4 (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    rho_4(extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    p_4  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    tp_4 (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &    Et_4 (extra_x*imax,extra_y*jmax,extra_z*kmax)  
      !MS$ENDIF 
      common/par/u,v,w,wx,wy,wz,Qi,rho,p,tp,Et

      NULLCHR   = CHAR(0)
      Debug     = 0
      !MS$IF (precision_plot.EQ.1)
        VIsPrecis = 0
      !MS$ELSEIF (precision_plot.EQ.2)
        VIsPrecis = 1
      !MS$ENDIF 
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y u v wz rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  extra_x*imax,extra_y*jmax,kmax,'BLOCK'//NULLCHR,NULLCHR)

      xy = dble(extra_x*imax*extra_y*jmax)

      do i=1,imax
        do j=1,extra_y*jmax
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

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
      
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
      real*8 x(extra_x*imax),y(extra_y*jmax),z(extra_z*kmax),
     &  u  (extra_x*imax,extra_y*jmax,extra_z*kmax),v  (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  w  (extra_x*imax,extra_y*jmax,extra_z*kmax),wx (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  wy (extra_x*imax,extra_y*jmax,extra_z*kmax),wz (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  Qi (extra_x*imax,extra_y*jmax,extra_z*kmax),rho(extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  p  (extra_x*imax,extra_y*jmax,extra_z*kmax),tp (extra_x*imax,extra_y*jmax,extra_z*kmax),
     &  Et (extra_x*imax,extra_y*jmax,extra_z*kmax)  
      !MS$IF (precision_plot.EQ.1)
        real*4 
     &    xx   (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    yy   (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    zz   (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    u_4  (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    v_4  (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    w_4  (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    wx_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    wy_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    wz_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    Qi_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    rho_4(extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    p_4  (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    tp_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    Et_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z))  
      !MS$ELSEIF (precision_plot.EQ.2)
        real*8 
     &    xx   (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    yy   (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    zz   (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    u_4  (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    v_4  (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    w_4  (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    wx_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    wy_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    wz_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    Qi_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    rho_4(extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    p_4  (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    tp_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z)),
     &    Et_4 (extra_x*(imax/qtimes_x),extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z))  
      !MS$ENDIF 
      common/par/u,v,w,wx,wy,wz,Qi,rho,p,tp,Et

      NULLCHR   = CHAR(0) 
      Debug     = 0
      !MS$IF (precision_plot.EQ.1)
        VIsPrecis = 0
      !MS$ELSEIF (precision_plot.EQ.2)
        VIsPrecis = 1
      !MS$ENDIF 
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'x y z u v w wx wy wz Qi rho p T Et'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,extra_x*(imax/qtimes_x),
     &  extra_y*(jmax/qtimes_y),extra_z*(kmax/qtimes_z),
     &  'BLOCK'//NULLCHR,NULLCHR)

      xyz = dble(extra_x*(imax/qtimes_x)*extra_y*(jmax/qtimes_y)*extra_z*(kmax/qtimes_z))

      if (qtimes_x.EQ.1) then
        iimax = extra_x*imax
      else
        iimax = imax-1
      end if

      if (qtimes_y.EQ.1) then
        jjmax = extra_y*jmax
      else
        jjmax = jmax-1
      end if

      if (qtimes_z.EQ.1) then
        kkmax = extra_z*kmax
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
      IORE = TecDat(xyz,yy,   VIsPrecis)
      IORE = TecDat(xyz,zz,   VIsPrecis)
      IORE = TecDat(xyz,u_4,  VIsPrecis)
      IORE = TecDat(xyz,v_4,  VIsPrecis)
      IORE = TecDat(xyz,w_4,  VIsPrecis)
      IORE = TecDat(xyz,wx_4, VIsPrecis)
      IORE = TecDat(xyz,wy_4, VIsPrecis)
      IORE = TecDat(xyz,wz_4, VIsPrecis)
      IORE = TecDat(xyz,Qi_4, VIsPrecis)
      IORE = TecDat(xyz,rho_4,VIsPrecis)
      IORE = TecDat(xyz,p_4,  VIsPrecis)
      IORE = TecDat(xyz,tp_4, VIsPrecis)    
      IORE = TecDat(xyz,Et_4, VIsPrecis)
      IORE = TecEnd()

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
      
      RETURN
      END SUBROUTINE write_PLT_TCP_3D

ccccc **********************************************************************
ccccc file_exists: This routine verify if the exists.
ccccc **********************************************************************
      FUNCTION file_exists(name)

      IMPLICIT NONE

      logical file_exists
      character*22 name
        
      inquire(file=name,exist=file_exists)
 
      RETURN
      END FUNCTION file_exists
	
ccccc **********************************************************************
ccccc makelayout: Save TECPLOT layout for movies.
ccccc **********************************************************************
      SUBROUTINE makelayout(t_ini,t_end)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname_aux
      character*120 linha_tecplot      
      integer tt,t_ini,t_end,dim

      open (1,file=dirgraf//'lay_anime.lay',access='sequential')
      !MS$IF (tp_dim.EQ.1)

      !MS$ELSEIF (tp_dim.EQ.2)       
        open (2,file=dirgraf//'tec_2d_sl.lay',access='sequential')
      !MS$ELSEIF (tp_dim.EQ.3)
        open (2,file=dirgraf//'tec_3d_sl.lay',access='sequential')
      !MS$ENDIF

      read (2,'(a080)') linha_tecplot
      write(1,'(a080)') linha_tecplot         
      read (2,'(a120)') linha_tecplot
      write(1,'(a120)') linha_tecplot         

      write(1,'(A22,A1,\)') '$!VarSet |LFDSFN1| = ',char(39)
      do tt=t_ini,t_end
        call get_filename(tt,'SN__S','.plt',fname_aux)           

        if (((mod(tt,qtimes).EQ.0)).AND.(tt.GT.0)) then
          write(1,'(A33,\)') '"|LFBD|/f_'//fname_aux//'"'
        end if
      end do
      write(1,'(A1)') char(39)

      read (2,'(a120)') linha_tecplot
      do while (.not. EOF(2))
        read (2,'(a120)') linha_tecplot
        write(1,'(a120)') linha_tecplot         
      end do     

      close (unit=1) 
      close (unit=2) 

      RETURN
      END SUBROUTINE makelayout

      END MODULE nsdat2plt

