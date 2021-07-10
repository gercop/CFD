      PROGRAM sngrid

      USE snmpi, only: initmpi,time_start,decompose,time_end,ierror

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      real*8 x(imax),y(jmax),z(kmax)

      call initmpi
      call time_start
      call decompose(imax)

      call main(x,y,z)
      call view_data(x,y,z)
      call save_grid(x,y,z)      

      !MS$IF (parallel.EQ.1)
        call MPI_FINALIZE(ierror)
      !MS$ENDIF
  
      END

ccccc **********************************************************************
ccccc main: Main subroutine
ccccc **********************************************************************
      SUBROUTINE main(x,y,z)
    
      USE snmethod, only: stx_d1_ex_D,stx_d2_ex_D,
     &  sty_a_cp_D,sty_b_cp_D,sty_c_cp_D,
     &  sty_a_cp_N,sty_b_cp_N,sty_c_cp_N,
     &  stz_a_cp_P,stz_b_cp_P,stz_c_cp_P,
     &  d1cy_dpy1,d2cy_dpy2,stx_d1_ex_6th_D,stx_d2_ex_6th_D,
     &  lhs_d1_62th_D,lhs_d1_62th_N,lhs_d1_6th_P

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer i,j,k,imax_spo
      real*8 x(imax),y(jmax),z(kmax),t,cyy,py_h2,aa,bb,B,B1,B2,
     &  du_0_dy_aux(jmax),delta_w,py_aux(jmax),
     &  d1cy_dpy1_aux(jmax),d2cy_dpy2_aux(jmax)

      do k=1,kmax
        do j=1,jmax         
          do i=1,imax
            x(i) = 0.d0
            y(j) = 0.d0
            z(k) = 0.d0
          end do
        end do
      end do

cccccc ************************************
cccccc ************************************
cccccc DOMAIN IN STREAMWISE DIRECTION
cccccc ************************************
cccccc ************************************
      !MS$IF(tp_s.EQ.2).OR.(tp_s.EQ.3).OR.(tp_s.EQ.4)             
        do i=1,imax
          x(i) = dble(i-1)*dx
        end do
      !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)        
        !MS$IF (stretching_in_x.EQ.0)
          do i=1,imax
            x(i) = x_0 + dble(i-1)*dx
          end do         
        !MS$ELSEIF (stretching_in_x.EQ.1)                      
          imax_spo = imax-imax_phy    !Sponge domain - Points from imax_phy..(imax_phy+imax_spo)          

          do i=1,imax_phy+1
            x(i) = x_0 + dble(i-1)*Lx/dble(imax_phy-1)
          end do 

          aa = ((Lx_spo-x(imax_phy))-dble(imax_spo)*
     &      (x(imax_phy+1)-x(imax_phy)))/
     &      (dble(imax_spo)**sp_x-dble(imax_spo))
          bb = (x(imax_phy+1)-x(imax_phy))-aa 

          do i=imax_phy,imax
            x(i) = x(imax_phy)+
     &        (aa*dble(i-imax_phy)**sp_x+bb*dble(i-imax_phy))
          end do
        !MS$ENDIF  
      !MS$ENDIF 

cccccc ************************************
cccccc ************************************
cccccc DOMAIN IN NORMALWISE DIRECTION
cccccc ************************************
cccccc ************************************
      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
        do j=1,jmax 
          y(j)         = dble(j-1)*dy
          d1cy_dpy1(j) = 1.d0
          d2cy_dpy2(j) = 0.d0
        end do
      !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        do j=1,jmax 
          y(j)           = dble(j-1)*dy 
          du_0_dy_aux(j) = 1.d0-tanh(2.d0*y(j))**2.d0  
        end do
        call maximo_1D(jmax,du_0_dy_aux,delta_w) !delta_w = (M1-M2)*c/du0_dy_max) = 1.0

        !MS$IF (stretching_in_y.EQ.0)
          do j=1,jmax 
            y(j)         = (y(j)-y(jmax)/2.d0)/delta_w 
            d1cy_dpy1(j) = 1.d0
            d2cy_dpy2(j) = 0.d0
          end do    
        !MS$ELSEIF (stretching_in_y.EQ.1)  
          py_h2 = y((jmax+1)/2) 
          B     = 1.d0/(2.d0*sp_y)*dlog( 
     &      (1.d0+(dexp(+sp_y)-1.d0)*(py_h2/Ly))/  
     &      (1.d0+(dexp(-sp_y)-1.d0)*(py_h2/Ly)))
 
          do j=1,jmax 
            cyy       = dble(j-1)*dcy
            py_aux(j) = py_h2*(1.d0+(dsinh(sp_y*(cyy-B))/dsinh(sp_y*B)))

            d1cy_dpy1_aux(j)   = 
     &        dsinh(sp_y*B)/((sp_y*py_h2)*
     &        dsqrt(1.d0+(py_aux(j)/py_h2-1.d0)**2.d0*
     &        dsinh(sp_y*B)**2.d0))
            d2cy_dpy2_aux(j) = 
     &        (-dsinh(sp_y*B)**3.d0*(py_aux(j)/py_h2-1.d0))
     &        /(sp_y*py_h2**2.d0*( 1.d0+(py_aux(j)/py_h2-1.d0)**2.d0
     &        *dsinh(sp_y*B)**2.d0)**(3.d0/2.d0) ) 
          end do    

          do j=1,jmax 
            y        (j) = (py_aux       (j) - py_aux       (jmax)/2.d0)
            d1cy_dpy1(j) = (d1cy_dpy1_aux(j) - d1cy_dpy1_aux(jmax)/2.d0)
            d2cy_dpy2(j) = (d2cy_dpy2_aux(j) - d2cy_dpy2_aux(jmax)/2.d0)
          end do
        !MS$ENDIF 
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)
          do j=1,jmax 
            y(j)         = dble(j-1)*dy 
            d1cy_dpy1(j) = 1.d0
            d2cy_dpy2(j) = 0.d0
          end do
        !MS$ELSEIF (stretching_in_y.EQ.1)            
          B1 = ((sp_y+1.d0)/(sp_y-1.d0))
          B2 = dlog(B1)

          do j=1,jmax 
            cyy       = dble(j-1)*dcy
            py_aux(j) = Ly*( (B+1.d0)-(B-1.d0)*(B1**(1.d0-cyy)) )
     &        /(B1**(1.d0-cyy)+1.d0)
                       
            d1cy_dpy1_aux(j) = 2.d0*Ly*sp_y/
     &        (sp_y*Ly-Ly+py_aux(j))/(sp_y*Ly+Ly-py_aux(j))/B2
            d2cy_dpy2_aux(j) = 4.d0*sp_y*Ly*(py_aux(j)-Ly)/
     &        (sp_y*Ly-Ly+py_aux(j))**2/(sp_y*Ly+Ly-py_aux(j))**2.d0/B2
          end do    

          do j=1,jmax 
            y        (j) = py_aux(j)
            d1cy_dpy1(j) = d1cy_dpy1_aux(j)
            d2cy_dpy2(j) = d2cy_dpy2_aux(j)
          end do
        !MS$ENDIF 
      !MS$ENDIF 

ccccc ************************************
ccccc ************************************
ccccc DOMAIN IN SPANWISE DIRECTION
ccccc ************************************
ccccc ************************************
      do k=1,kmax
        z(k)=DBLE(k-1)*dz
      end do

ccccc ************************************
ccccc CALCULATING STENCILS
ccccc ************************************
      call stx_d1_ex_6th_D(0,1,imax,6,x,stx_d1_ex_D )
      call stx_d2_ex_6th_D(0,1,imax,6,x,stx_d2_ex_D )
      if (jmax.NE.1) then 
        call lhs_d1_62th_D(sty_a_cp_D,sty_b_cp_D,sty_c_cp_D,d1cy_dpy1,jmax)
        call lhs_d1_62th_N(sty_a_cp_N,sty_b_cp_N,sty_c_cp_N,d1cy_dpy1,jmax)
      end if
      if (kmax.NE.1) then 
        call lhs_d1_6th_P(stz_a_cp_P,stz_b_cp_P,stz_c_cp_P,kmax)
      end if

      RETURN
      END

ccccc **********************************************************************
ccccc save_grid: Save Grid.
ccccc **********************************************************************
      SUBROUTINE save_grid(x,y,z)
    
      USE snmpi,    only: my_rank,pro,global_x1,global_x2,
     &  ghost_points,protostr
      USE snmethod, only: stx_d1_ex_D,stx_d2_ex_D,
     &  sty_a_cp_D,sty_b_cp_D,sty_c_cp_D,
     &  sty_a_cp_N,sty_b_cp_N,sty_c_cp_N,
     &  stz_a_cp_P,stz_b_cp_P,stz_c_cp_P,
     &  d1cy_dpy1,d2cy_dpy2

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*16 fname_aux
      character*04 ext
      integer i,j,k
      real*8 x(imax),y(jmax),z(kmax),local_x1,local_x2,global_x(imax)

      !MS$IF (parallel.EQ.1)
        call protostr(my_rank,'f_grid_x_P','.cx ',fname_aux)
        open (1,file=dirgraf//fname_aux,access='sequential')   
        if (my_rank.EQ.0) then 
          local_x1 = global_x1
          local_x2 = global_x2+ghost_points
        else if (my_rank+1.EQ.pro) then 
          local_x1 = global_x1-ghost_points
          local_x2 = global_x2
        else
          local_x1 = global_x1-ghost_points
          local_x2 = global_x2+ghost_points
        end if

        do i=local_x1,local_x2
          write(1,'(e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,
     &      e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)') x(i),
     &      stx_d1_ex_D(i,1),stx_d1_ex_D(i,2),stx_d1_ex_D(i,3),stx_d1_ex_D(i,4),
     &      stx_d1_ex_D(i,5),stx_d1_ex_D(i,6),stx_d1_ex_D(i,7),
     &      stx_d2_ex_D(i,1),stx_d2_ex_D(i,2),stx_d2_ex_D(i,3),stx_d2_ex_D(i,4),
     &      stx_d2_ex_D(i,5),stx_d2_ex_D(i,6),stx_d2_ex_D(i,7)
        end do
        close (unit=1)   
      !MS$ENDIF

      if (my_rank.NE.0) RETURN

      open (1,file=dirgraf//'f_grid_x.cx',access='sequential')   
      do i=1,imax
          write(1,'(e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,
     &      e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)') x(i),
     &      stx_d1_ex_D(i,1),stx_d1_ex_D(i,2),stx_d1_ex_D(i,3),stx_d1_ex_D(i,4),
     &      stx_d1_ex_D(i,5),stx_d1_ex_D(i,6),stx_d1_ex_D(i,7),
     &      stx_d2_ex_D(i,1),stx_d2_ex_D(i,2),stx_d2_ex_D(i,3),stx_d2_ex_D(i,4),
     &      stx_d2_ex_D(i,5),stx_d2_ex_D(i,6),stx_d2_ex_D(i,7)
      end do
      close (unit=1)

      open (1,file=dirgraf//'f_grid_y.cy',access='sequential')   
      do j=1,jmax
        write(1,'(e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16,e30.16)') y(j),
     &    sty_a_cp_D(j),sty_b_cp_D(j),sty_c_cp_D(j),sty_a_cp_N(j),sty_b_cp_N(j),sty_c_cp_N(j),
     &    d1cy_dpy1(j),d2cy_dpy2(j)    
      end do
      close (unit=1)

      open (1,file=dirgraf//'f_grid_z.cz',access='sequential')   
      do k=1,kmax
        write(1,'(e30.16,e30.16,e30.16,e30.16)') z(k),stz_a_cp_P(k),
     &    stz_b_cp_P(k),stz_c_cp_P(k)  
      end do
      close (unit=1)

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
      SUBROUTINE view_data(x,y,z)

      USE snmpi, only: my_rank

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer n,i,j,k
      real*8 x(imax),y(jmax),z(kmax),
     &  dx_i(imax-1),dy_j(jmax-1),dz_k(kmax-1),
     &  dx_max,dy_max,dz_max,dx_min,dy_min,dz_min

      if (my_rank.NE.0) RETURN

      open(1,file=dirgraf//'INIT_GRID.TXT',status='unknown')
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
          write(n,*)'SIMULATION WITH PROGRESSIVE ACOUSTIC WAVE'
        !MS$ELSEIF (tp_s.EQ.3)
          write(n,*)'SIMULATION WITH STANDING ACOUSTIC WAVE'
        !MS$ELSEIF (tp_s.EQ.4)
          write(n,*)'SHEAR LAYER SIMULATION - TEMPORAL DEVELOPMENT    '
        !MS$ELSEIF (tp_s.EQ.5)
          write(n,*)'SHEAR LAYER SIMULATION - SPATIAL DEVELOPMENT     '
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
        write (n,*)'SPATIAL COMPUTATIONAL MESH'
        write (n,*)'--------------------------------------------------'
        write (n,'(A19,F7.2,F7.2,F7.2)') 'SPATIAL DOMAIN..: ',Lx,Ly,Lz
        !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
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

        write (n,'(A19,A2,F20.15)')'LX..............: ','  ',Lx
        !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
          write (n,'(A19,A2,F20.15)')'LY..............: ','  ',Ly
        !MS$ENDIF
        !MS$IF (tp_dim.EQ.3)
          write (n,'(A19,A2,F20.15)')'LZ..............: ','  ',Lz
        !MS$ENDIF
        write (n,'(A19,A4,I6,I6,I6)')'IMAX/JMAX/KMAX..: ','    ',
     &    imax,jmax,kmax
        write (n,*)'***************************************************'
        write (n,*) '**************************************************'  
        write (n,*)
      end do
      close (unit=1)	  	 

      RETURN
      END
