      MODULE sndat2mat

      contains
ccccc **********************************************************************
ccccc dat2plt: Main subroutine.
ccccc **********************************************************************
      SUBROUTINE dat2mat(t_ini,t_end)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname_aux
      character*05 fname
      integer tt,t_ini,t_end,error
      real*8 x(imax),y(jmax),z(kmax),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w  (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p  (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  Qi (imax,jmax,kmax),U1 (imax,jmax,kmax),
     &  U2 (imax,jmax,kmax),U3 (imax,jmax,kmax),
     &  U4 (imax,jmax,kmax),U5 (imax,jmax,kmax)

      call init_grid(dirgraf,x,y,z)

      write(*,*) ' _______________________________________________'
      write(*,*) 
      write(*,*) ' CONVERTENDO ARQUIVOS DE FORMATO BINÁRIO PARA   '   
      write(*,*) ' O FORMATO DE LEITURA PADRÃO MATLAB...          '  
      write(*,*) 

      do tt=t_ini,t_end
        if (mod(tt,qtimes).EQ.0) then
          call get_filename(tt,fnamegraf,'.dat',fname_aux)           

          !MS$IF (tp_dim.EQ.1)  
            call read_BIN_1D(dirgraf,fname_aux,U1,U2,U5,error)
          !MS$ELSEIF (tp_dim.EQ.2)
            call read_BIN_2D(dirgraf,fname_aux,U1,U2,U3,U5,error)
          !MS$ELSEIF (tp_dim.EQ.3)
            call read_BIN_3D(dirgraf,fname_aux,U1,U2,U3,U4,U5,error)
          !MS$ENDIF    

          if (error.EQ.0) then
            call flux_to_primitive(U1,U2,U3,U4,U5,rho,u,v,w,Et,p,tp)
            call equation_vorticity(u,v,w,wx,wy,wz)

            call get_filename(tt,fnamegraf,'.m',fname_aux) 

            !MS$IF (tp_dim.EQ.1)               
              call write_MAT_1D(tt,fname_aux,x,u,rho,p,tp,Et)            
            !MS$ELSEIF (tp_dim.EQ.2)
              call write_MAT_2D(tt,fname_aux,x,y,u,v,wz,rho,p,tp,Et)            
            !MS$ELSEIF (tp_dim.EQ.3)
              call calculate_iso(u,v,w,Qi)
              call write_MAT_3D(tt,fname_aux,x,y,z,u,v,w,wx,wy,wz,rho,p,tp,Et)
            !MS$ENDIF
          end if
        end if
      end do

      open (1,file=dirgraf//'mat_layout.m',status='unknown',access='sequential') 
      write(1,*) 'figure(1);'
      write(1,*) 'colormap(gray);'
      write(1,*) 'contourf(y,x,wz,400);'
      write(1,*) 'axis([ -5 +5 1 16]);'
      write(1,*) 'view(90,270);'
      write(1,*) 'shading flat;'
      write(1,*) 'axis tight;'
      write(1,*) 'grid on;'   
      write(1,*) 'axis off;' 
      close(1)    

      open (1,file=dirgraf//'mat_showall.m',status='unknown',form='formatted') 
      do tt=t_ini,t_end
        if (mod(tt,qtimes).EQ.0) then
          call get_filename(tt,'SN_WZ',';   ',fname_aux) 
          write(1,*) 'f_'//fname_aux
        end if
      end do
      close(1)

      RETURN
      END SUBROUTINE dat2mat

ccccc **********************************************************************
ccccc write_MAT_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_MAT_2D(tt,fname,x,y,u,v,wz,rho,p,tp,Et)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*20 fname_aux
      integer i,j,k,tt
      real*8 x(imax),y(jmax),
     &  xx (imax,jmax,kmax),yy (imax,jmax,kmax),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),  
     &  wz (imax,jmax,kmax)

      open (1,file=dirgraf//'f_SN__x.m',status='unknown',access='sequential') 
      do i=1,imax 
        write(1,'(a3,i6,a2,e30.16,a1)') 'x(',i,')=',x(i),';'
      end do
      close(1)

      open (1,file=dirgraf//'f_SN__y.m',status='unknown',access='sequential') 
      do j=1,jmax 
        write(1,'(a3,i6,a2,e30.16,a1)') 'y(',j,')=',y(j),';'
      end do
      close(1)

      open (1,file=dirgraf//'f_'//fname,status='unknown',access='sequential') 
      write(1,*) 'clear all' 
      write(1,*) 'close all;'
      do j=1,jmax
        do i=1,imax 
          write(1,'(a3,i6,a1,i6,a2,e30.16,a1)')
     &      'wz(',i,',',j,')=',wz(i,j,1),';'
        end do
      end do
      write(1,*) 'f_SN__x;'
      write(1,*) 'f_SN__y;'
      write(1,*) 'mat_layout;'
      call get_filename(tt,'re500','.eps',fname_aux)
      write(1,'(a6,a1,a7,a1,a1,a1,a5,a1,a1,a1,a28,a20,a1,a2)') 
     &  'print(',"'",'-depsc2',"'",',',"'",'-r300',"'",',',
     &  "'",'figs/wz_isen_twomode_mach04_',fname_aux,"'",');'     
      close(1)

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
      
      RETURN
      END SUBROUTINE write_MAT_2D

ccccc **********************************************************************
ccccc write_MAT_3D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_MAT_3D(tt,fname,x,y,z,u,v,w,wx,wy,wz,rho,p,tp,Et)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*20 fname_aux
      integer i,j,k,tt
      real*8 x(imax),y(jmax),z(kmax),
     &  xx (imax,jmax,kmax),yy (imax,jmax,kmax),
     &  zz (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  w  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),
     &  wx (imax,jmax,kmax),wy (imax,jmax,kmax),
     &  wz (imax,jmax,kmax)

      open (1,file=dirgraf//'f_SN__x.m',status='unknown',access='sequential') 
      do i=1,imax 
        write(1,'(a3,i6,a2,e30.16,a1)') 'x(',i,')=',x(i),';'
      end do
      close(1)

      open (1,file=dirgraf//'f_SN__y.m',status='unknown',access='sequential') 
      do j=1,jmax 
        write(1,'(a3,i6,a2,e30.16,a1)') 'y(',j,')=',y(j),';'
      end do
      close(1)

      open (1,file=dirgraf//'f_SN__z.m',status='unknown',access='sequential') 
      do k=1,kmax 
        write(1,'(a3,i6,a2,e30.16,a1)') 'z(',k,')=',z(k),';'
      end do
      close(1)

      open (1,file=dirgraf//'f_'//fname,status='unknown',access='sequential') 
      write(1,*) 'clear all' 
      write(1,*) 'close all;'
      do k=1,kmax
        do j=1,jmax
          do i=1,imax 
            write(1,'(a3,i6,a1,i6,a2,i6,a1,e30.16,a1)')
     &        ' u(',i,',',j,',',k,')=',u(i,j,k),';'
            write(1,'(a3,i6,a1,i6,a2,i6,a1,e30.16,a1)')
     &        ' v(',i,',',j,',',k,')=',v(i,j,k),';'
            write(1,'(a3,i6,a1,i6,a2,i6,a1,e30.16,a1)')
     &        ' w(',i,',',j,',',k,')=',w(i,j,k),';'
            write(1,'(a3,i6,a1,i6,a2,i6,a1,e30.16,a1)')
     &        'wz(',i,',',j,',',k,')=',wz(i,j,k),';'
          end do
        end do
      end do

      write(1,*) 'f_SN__x;'
      write(1,*) 'f_SN__y;'
      write(1,*) 'f_SN__z;'
      write(1,*) 'mat_layout;'
      call get_filename(tt,'4_re500','.eps',fname_aux)
      write(2,'(a6,a1,a7,a1,a1,a1,a5,a1,a1,a1,a33,a22,a1,a2)') 
     &  'print(',"'",'-depsc2',"'",',',"'",'-r300',"'",',',
     &  "'",'/extra/tese/figs/wz_twomode_mach0',fname_aux,"'",');'     
      close(1)

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
      
      RETURN
      END SUBROUTINE write_MAT_3D

      END MODULE sndat2mat

