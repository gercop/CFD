      MODULE nsdat2ASC

      contains
ccccc **********************************************************************
ccccc dat2asc: Main subroutine.
ccccc **********************************************************************
      SUBROUTINE dat2ASC(t_ini,t_end)

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

      write(*,*) ' _______________________________________________________'
      write(*,*) 
      write(*,*) ' CONVERTING FILES FROM BINARY FORMAT TO ASC PATTERN...  '   
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

            call get_filename(tt,fnamegraf,'.dem',fname_aux) 

            !MS$IF (tp_dim.EQ.1)               
            !MS$ELSEIF (tp_dim.EQ.2)
              call write_ASC_2D(tt,fname_aux,x,y,u,v,wz,rho,p,tp,Et)            
            !MS$ELSEIF (tp_dim.EQ.3)
              call calculate_iso(u,v,w,Qi)
            !MS$ENDIF
          end if
        end if
      end do

      RETURN
      END SUBROUTINE dat2ASC

ccccc **********************************************************************
ccccc write_ASC_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_ASC_2D(tt,fname,x,y,u,v,wz,rho,p,tp,Et)

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

      open (1,file=dirgraf//'f_'//fname,status='unknown',access='sequential') 
      do j=1,jmax
        do i=1,imax        
          write(1,*) x(i),y(j),u(i,j,1)
        end do
        write(1,*)
      end do
      close(1)

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
      
      RETURN
      END SUBROUTINE write_ASC_2D

      END MODULE nsdat2ASC

