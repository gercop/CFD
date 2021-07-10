      PROGRAM nsprocess

      USE nsdat2ASC, only: dat2ASC
      USE nsdat2plt, only: dat2plt
      USE nsdat2mat, only: dat2mat
      USE nsdat2fft, only: dat2fft
      
      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer qtimes_local,t_ini,t_end

      t_ini = 0
      t_end = 0

      write(*,*) 
      write(*,*) ' _______________________________________________________'
      write(*,*) 
      write(*,*) ' POST PROCESSING APPLICATION                            '
      write(*,*) ' _______________________________________________________'
      write(*,*) 
      write(*,'(a21,\)') 'INICIAL TIME......:'
      read (*,*) t_ini
      write(*,*) 
      write(*,'(a21,\)') 'END TIME..........:'
      read (*,'(i6,\)') t_end

      !MS$IF (tp_dat2ASC.EQ.1)
        call dat2ASC(t_ini,t_end)
      !MS$ENDIF
      !MS$IF (tp_dat2plt.EQ.1)
        call dat2plt(t_ini,t_end)
      !MS$ENDIF
      !MS$IF (tp_dat2mat.EQ.1)
        call dat2mat(t_ini,t_end)
      !MS$ENDIF
      !MS$IF (tp_dat2fft.EQ.1)
        qtimes_local = ((t_end-t_ini)/qtimes+1)
        !MS$IF (tp_dim.EQ.2)
          call dat2fft(t_ini,t_end,qtimes_local,imax,kmax,'SN_F1')
        !MS$ELSEIF (tp_dim.EQ.3)
          call dat2fft(t_ini,t_end,qtimes_local,imax,kmax,'SN_F2')
        !MS$ENDIF
      !MS$ENDIF

      write(*,*) ' _______________________________________________________'
      write(*,*) 
      write(*,*) ' END PROCESSING                                         '
      write(*,*) ' _______________________________________________________'
      write(*,*) 

      END

ccccc **********************************************************************
ccccc init_grid: This routine calculates the stencils for all derivatives.
ccccc **********************************************************************
      SUBROUTINE init_grid(dir,x,y,z)

      USE snmethod, only: stx_d1_ex_D,stx_d2_ex_D,
     &  sty_a_cp_D,sty_b_cp_D,sty_c_cp_D,
     &  sty_a_cp_N,sty_b_cp_N,sty_c_cp_N,
     &  stz_a_cp_P,stz_b_cp_P,stz_c_cp_P,
     &  d1cy_dpy1,d2cy_dpy2   

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*07 dir
      integer i,j,k,pro_aux	
      real*8 x(imax),y(jmax),z(kmax),Ma_aux,Re_aux

      open(1,file=dir//'INIT_DATA.DAT',status='old',form='formatted')
      read(1,*) Ma_aux
      read(1,*) Re_aux
      read(1,*) pro_aux
      close(1)

      open(1,file=dir//'f_grid_x.cx',access='sequential')
      do i=1,imax
        read (1,'(e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,
     &      e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)') x(i),
     &      stx_d1_ex_D(i,1),stx_d1_ex_D(i,2),stx_d1_ex_D(i,3),stx_d1_ex_D(i,4),
     &      stx_d1_ex_D(i,5),stx_d1_ex_D(i,6),stx_d1_ex_D(i,7),
     &      stx_d2_ex_D(i,1),stx_d2_ex_D(i,2),stx_d2_ex_D(i,3),stx_d2_ex_D(i,4),
     &      stx_d2_ex_D(i,5),stx_d2_ex_D(i,6),stx_d2_ex_D(i,7)
      end do
      close (unit=1)

      open (1,file=dir//'f_grid_y.cy',access='sequential')
      do j=1,jmax
        read(1,'(e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)') y(j),
     &    sty_a_cp_D(j),sty_b_cp_D(j),sty_c_cp_D(j),sty_a_cp_N(j),sty_b_cp_N(j),sty_c_cp_N(j),
     &    d1cy_dpy1(j),d2cy_dpy2(j)    
      end do
      close (unit=1)

      open (1,file=dir//'f_grid_z.cz',access='sequential')
      do k=1,kmax
        read(1,'(e30.16,e30.16,e30.16,e30.16)') z(k),stz_a_cp_P(k),
     &    stz_b_cp_P(k),stz_c_cp_P(k)  
      end do
      close (unit=1)

      RETURN
      END 

ccccc **********************************************************************
ccccc flux_to_primitive: Convert primitive parameters to flux parameters.
ccccc **********************************************************************
      SUBROUTINE flux_to_primitive(U1,U2,U3,U4,U5,rho,u,v,w,Et,p,tp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer i,j,k
      real*8 C2,
     &  U1 (imax,jmax,kmax),U2 (imax,jmax,kmax),
     &  U3 (imax,jmax,kmax),U4 (imax,jmax,kmax),
     &  U5 (imax,jmax,kmax),rho(imax,jmax,kmax),
     &  u  (imax,jmax,kmax),v  (imax,jmax,kmax),
     &  w  (imax,jmax,kmax),Et (imax,jmax,kmax),
     &  p  (imax,jmax,kmax),tp (imax,jmax,kmax),
     &  e  (imax,jmax,kmax)
      parameter (C2=110.4d0/288.15d0) !T_oo=288.15 for sea level.

      do k=1,kmax
        do j=1,jmax
          do i=1,imax              
            rho(i,j,k) = U1(i,j,k)
            u  (i,j,k) = U2(i,j,k)/U1(i,j,k) 
            v  (i,j,k) = U3(i,j,k)/U1(i,j,k)
            w  (i,j,k) = U4(i,j,k)/U1(i,j,k)                        

            !MS$IF (tp_proc.EQ.1)
              Et(i,j,k) = U5(i,j,k)
              e (i,j,k) = Et(i,j,k)/rho(i,j,k)
     &          -(u(i,j,k)**2.d0+v(i,j,k)**2.d0+w(i,j,k)**2.d0)/2.d0
              p (i,j,k) = (gamma-1.d0)*rho(i,j,k)*e(i,j,k)
              tp(i,j,k) = gamma*Ma**2.d0*p(i,j,k)/rho(i,j,k) 
            !MS$ELSEIF (tp_proc.EQ.2)
              if (rho(i,j,k).LT.0) then
                write(*,*) 
                write(*,*) 'ERROR..: Negative Density in Point (',i,j,k,').'	       
                write(*,*) 'WARNING: THE PROGRAMM WILL BE TERMINATE.'	                    
                STOP          
              end if 
              p (i,j,k) = 1.d0/(gamma*Ma**2.d0)*rho(i,j,k)**gamma
              tp(i,j,k) = gamma*Ma**2.d0*p(i,j,k)/rho(i,j,k)
              e (i,j,k) = p(i,j,k)/(rho(i,j,k)*(gamma-1.d0))
              Et(i,j,k) = rho(i,j,k)*(e(i,j,k)
     &          +(u(i,j,k)**2.d0+v(i,j,k)**2.d0+w(i,j,k)**2.d0)/2.d0) 
            !MS$ENDIF
          end do
        end do
      end do
      
      RETURN
      END

ccccc **********************************************************************
ccccc equation_vorticity: This routine calculates the vorticity of flow.
ccccc **********************************************************************
      SUBROUTINE equation_vorticity(u,v,w,wx,wy,wz)

      USE snmethod, only: der_x_DD,der_y_WN,der_z_PP

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer i,j,k
      real*8 
     &  u    (imax,jmax,kmax),v    (imax,jmax,kmax),
     &  w    (imax,jmax,kmax),wx   (imax,jmax,kmax),
     &  wy   (imax,jmax,kmax),wz   (imax,jmax,kmax),
     &  dv_dz(imax,jmax,kmax),dw_dy(imax,jmax,kmax),
     &  dw_dx(imax,jmax,kmax),du_dz(imax,jmax,kmax),
     &  du_dy(imax,jmax,kmax),dv_dx(imax,jmax,kmax)
    
      call der_x_DD(0,1,imax,v,dv_dx)
      call der_x_DD(0,1,imax,w,dw_dx)
      call der_y_WN(imax,u,du_dy)
      call der_y_WN(imax,w,dw_dy)
      call der_z_PP(imax,v,dv_dz)
      call der_z_PP(imax,u,du_dz)

      do k=1,kmax
        do j=1,jmax         
          do i=1,imax
            wx(i,j,k) = dv_dz(i,j,k)-dw_dy(i,j,k)
            wy(i,j,k) = dw_dx(i,j,k)-du_dz(i,j,k)
            wz(i,j,k) = du_dy(i,j,k)-dv_dx(i,j,k)
          end do 
        end do
      end do

      RETURN
      END

ccccc **********************************************************************
ccccc calculate_iso: Calculate figures with isosuperficie.
ccccc **********************************************************************
      SUBROUTINE calculate_iso(u,v,w,Qi)

      USE snmethod, only: der_x_DD,der_y_WD,der_y_WN,der_z_PP

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer i,j,k
      real*8 d1cy_dpy1(jmax),
     &  u    (imax,jmax,kmax),v    (imax,jmax,kmax),
     &  w    (imax,jmax,kmax),du_dx(imax,jmax,kmax),
     &  dv_dx(imax,jmax,kmax),dw_dx(imax,jmax,kmax),
     &  du_dy(imax,jmax,kmax),dv_dy(imax,jmax,kmax),
     &  dw_dy(imax,jmax,kmax),du_dz(imax,jmax,kmax),
     &  dv_dz(imax,jmax,kmax),dw_dz(imax,jmax,kmax),
     &  Qi   (imax,jmax,kmax)

      call der_x_DD(0,1,imax,u,du_dx) 
      call der_x_DD(0,1,imax,v,dv_dx) 
      call der_x_DD(0,1,imax,w,dw_dx) 
      call der_y_WD(imax,u,du_dy) 
      call der_y_WD(imax,v,dv_dy) 
      call der_y_WD(imax,w,dw_dy) 
      call der_z_PP(imax,u,du_dz) 
      call der_z_PP(imax,v,dv_dz) 
      call der_z_PP(imax,w,dw_dz) 

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


