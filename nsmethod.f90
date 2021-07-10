      MODULE snmethod

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      save
        real*8,  dimension(imax,7), public :: stx_d1_ex_D
        real*8,  dimension(imax,7), public :: stx_d2_ex_D

        real*8,  dimension(jmax),   public :: sty_a_cp_D
        real*8,  dimension(jmax),   public :: sty_b_cp_D
        real*8,  dimension(jmax),   public :: sty_c_cp_D
        real*8,  dimension(jmax),   public :: sty_a_cp_N
        real*8,  dimension(jmax),   public :: sty_b_cp_N
        real*8,  dimension(jmax),   public :: sty_c_cp_N

        real*8,  dimension(kmax),   public :: stz_a_cp_P
        real*8,  dimension(kmax),   public :: stz_b_cp_P
        real*8,  dimension(kmax),   public :: stz_c_cp_P

        real*8,  dimension(jmax),   public :: d1cy_dpy1
        real*8,  dimension(jmax),   public :: d2cy_dpy2

      contains
ccccc **********************************************************************
ccccc der_x_ex_6th_D: This routine calculated the first derivatives in 
ccccc                 the x-direction using the 6th order centered finite 
ccccc                 difference. Dirichlet boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_x_ex_6th_D(my_rank,pro,imax,jmax,kmax,fc,dfc_dx)

      integer i,j,k,imax,jmax,kmax,my_rank,pro
      real*8 fc(imax,jmax,kmax),dfc_dx(imax,jmax,kmax)

      if (my_rank.EQ.0) then 
        do k=1,kmax
          do j=1,jmax
            do i=1,3 
              dfc_dx(i,j,k) = 
     &          stx_d1_ex_D(i,1)*fc(1,j,k)+
     &          stx_d1_ex_D(i,2)*fc(2,j,k)+
     &          stx_d1_ex_D(i,3)*fc(3,j,k)+
     &          stx_d1_ex_D(i,4)*fc(4,j,k)+
     &          stx_d1_ex_D(i,5)*fc(5,j,k)+
     &          stx_d1_ex_D(i,6)*fc(6,j,k)+
     &          stx_d1_ex_D(i,7)*fc(7,j,k)
            end do 
        
            do i=4,imax-3 
              dfc_dx(i,j,k) = 
     &          stx_d1_ex_D(i,1)*fc(i-3,j,k)+
     &          stx_d1_ex_D(i,2)*fc(i-2,j,k)+
     &          stx_d1_ex_D(i,3)*fc(i-1,j,k)+
     &          stx_d1_ex_D(i,4)*fc(i  ,j,k)+
     &          stx_d1_ex_D(i,5)*fc(i+1,j,k)+
     &          stx_d1_ex_D(i,6)*fc(i+2,j,k)+
     &          stx_d1_ex_D(i,7)*fc(i+3,j,k)
            end do
          end do
        end do
      end if

      if ((my_rank.NE.0).AND.(my_rank+1.NE.pro)) then 
        do k=1,kmax
          do j=1,jmax
            do i=4,imax-3 
              dfc_dx(i,j,k) = 
     &          stx_d1_ex_D(i,1)*fc(i-3,j,k)+
     &          stx_d1_ex_D(i,2)*fc(i-2,j,k)+
     &          stx_d1_ex_D(i,3)*fc(i-1,j,k)+
     &          stx_d1_ex_D(i,4)*fc(i  ,j,k)+
     &          stx_d1_ex_D(i,5)*fc(i+1,j,k)+
     &          stx_d1_ex_D(i,6)*fc(i+2,j,k)+
     &          stx_d1_ex_D(i,7)*fc(i+3,j,k)
            end do
          end do
        end do
      end if

      if (my_rank+1.EQ.pro) then 
        do k=1,kmax
          do j=1,jmax
            do i=4,imax-3 
              dfc_dx(i,j,k) = 
     &          stx_d1_ex_D(i,1)*fc(i-3,j,k)+
     &          stx_d1_ex_D(i,2)*fc(i-2,j,k)+
     &          stx_d1_ex_D(i,3)*fc(i-1,j,k)+
     &          stx_d1_ex_D(i,4)*fc(i  ,j,k)+
     &          stx_d1_ex_D(i,5)*fc(i+1,j,k)+
     &          stx_d1_ex_D(i,6)*fc(i+2,j,k)+
     &          stx_d1_ex_D(i,7)*fc(i+3,j,k)
            end do

            do i=imax-2,imax 
              dfc_dx(i,j,k) = 
     &          stx_d1_ex_D(i,1)*fc(imax-6,j,k)+
     &          stx_d1_ex_D(i,2)*fc(imax-5,j,k)+
     &          stx_d1_ex_D(i,3)*fc(imax-4,j,k)+
     &          stx_d1_ex_D(i,4)*fc(imax-3,j,k)+
     &          stx_d1_ex_D(i,5)*fc(imax-2,j,k)+
     &          stx_d1_ex_D(i,6)*fc(imax-1,j,k)+
     &          stx_d1_ex_D(i,7)*fc(imax  ,j,k)
            end do  
          end do
        end do
      end if  

      RETURN
      END SUBROUTINE der_x_ex_6th_D

ccccc **********************************************************************
ccccc der_x_cp_6th_P: This routine calculated the first derivatives 
ccccc                 in the x-direction using the 6th order centered 
ccccc                 compact finite difference. Periodic boundary 
ccccc                 condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_x_cp_6th_P(dx,imax,jmax,kmax,fc,dfc_dx)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dx,fc(imax,jmax,kmax),dfc_dx(imax,jmax,kmax),
     &  dfc_dx_aux(imax),aa(imax),bb(imax),cc(imax),rhs(imax)  

      do k=1,kmax
        do j=1,jmax
          call rhs_x_6th_P(dx,imax,jmax,kmax,j,k,fc,rhs)
          call lhs_d1_6th_P(aa,bb,cc,imax)
          call trid_cyclic(aa,bb,cc,1.d0,1.d0,imax,rhs,dfc_dx_aux)

          do i=1,imax-1
            dfc_dx(i,j,k) = dfc_dx_aux(i)
	  end do
          dfc_dx(imax,j,k) = dfc_dx(1,j,k)          
        end do
      end do

      RETURN
      END SUBROUTINE der_x_cp_6th_P

ccccc **********************************************************************
ccccc der_y_cp_6th_P: This routine calculated the first derivatives 
ccccc                 in the y-direction using the 6th order centered 
ccccc                 compact finite-difference. Periodic boundary 
ccccc                 condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_6th_P(dy,imax,jmax,kmax,fc,dfc_dy)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,aa(jmax),bb(jmax),cc(jmax),rhs(jmax),
     &  fc(imax,jmax,kmax),dfc_dy_aux(jmax),dfc_dy(imax,jmax,kmax)

      !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
        do k=1,kmax
          do i=1,imax
            call rhs_y_6th_P(dy,imax,jmax,kmax,i,k,fc,rhs)
            call lhs_d1_6th_P(aa,bb,cc,jmax)
            call trid_cyclic(aa,bb,cc,1.d0,1.d0,jmax,rhs,dfc_dy_aux)

            do j=1,jmax-1
              dfc_dy(i,j,k) = dfc_dy_aux(j)
            end do
            dfc_dy(i,jmax,k) = dfc_dy(i,1,k)
          end do
        end do 
      !MS$ELSE
        dfc_dy(:,:,:) = 0.d0
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_y_cp_6th_P

ccccc **********************************************************************
ccccc der_y_cp_62th_D: This routine calculated the first derivatives 
ccccc                  in the y-direction using the 2nd and 6th order 
ccccc                  compact finite-difference schemes. Dirichlet 
ccccc                  boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_62th_D(dy,imax,jmax,kmax,mtrc,fc,dfc_dy)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc(jmax),fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),
     &  dfc_dy_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
        do k=1,kmax        
          do i=1,imax  
            call rhs_y_62th_D(dy,imax,jmax,kmax,i,k,fc,rhs)
            call trid(0,sty_a_cp_D,sty_b_cp_D,sty_c_cp_D,jmax,rhs,dfc_dy_aux)   
            do j=1,jmax
              dfc_dy(i,j,k) = dfc_dy_aux(j)
            end do	  
          end do
        end do
      !MS$ELSE
        dfc_dy(:,:,:) = 0.d0
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_y_cp_62th_D

ccccc **********************************************************************
ccccc der_y_cp_62th_N: This routine calculated the first derivatives 
ccccc                  in the y-direction using the 2nd and 6th order 
ccccc                  compact finite-difference schemes. Neumann 
ccccc                  boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_62th_N(dy,imax,jmax,kmax,mtrc,fc,dfc_dy)      

      IMPLICIT NONE
 
      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc(jmax),fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),
     &  dfc_dy_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
        do k=1,kmax
          do i=1,imax
            call rhs_y_62th_N(dy,imax,jmax,kmax,i,k,fc,rhs)
            call trid(0,sty_a_cp_N,sty_b_cp_N,sty_c_cp_N,jmax,rhs,dfc_dy_aux)   
            do j=1,jmax
              dfc_dy(i,j,k) = dfc_dy_aux(j)
            end do	  
          end do
        end do
      !MS$ELSE
        dfc_dy(:,:,:) = 0.d0
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_y_cp_62th_N

ccccc **********************************************************************
ccccc der_y_cp_65th_WD: This routine calculated the first derivatives 
ccccc                   in the y-direction using the 6th and 5th order 
ccccc                   compact finite difference. Dirichlet and Neumann 
ccccc                   boundary condition was used to simulated the wall.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_65th_WD(dy,imax,jmax,kmax,mtrc,fc,dfc_dy)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc(jmax),fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),
     &  dfc_dy_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
        do k=1,kmax        
          do i=1,imax
            call rhs_y_65th_WD(dy,imax,jmax,kmax,i,k,fc,rhs)
            call lhs_d1_65th_W(aa,bb,cc,mtrc,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,dfc_dy_aux)   

            do j=1,jmax
              dfc_dy(i,j,k) = dfc_dy_aux(j)
            end do	  
          end do
        end do
      !MS$ELSE
        dfc_dy(:,:,:) = 0.d0
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_y_cp_65th_WD						

ccccc **********************************************************************
ccccc der_y_cp_65th_WN: This routine calculated the first derivatives 
ccccc                   in the y-direction using the 6th and 5th order 
ccccc                   compact finite difference. Dirichlet and Neumann 
ccccc                   boundary condition was used to simulated the wall.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_65th_WN(dy,imax,jmax,kmax,mtrc,fc,dfc_dy)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc(jmax),fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),
     &  dfc_dy_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
        do k=1,kmax        
          do i=1,imax
            call rhs_y_65th_WN(dy,imax,jmax,kmax,i,k,fc,rhs)
            call lhs_d1_65th_W(aa,bb,cc,mtrc,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,dfc_dy_aux)   

            do j=1,jmax
              dfc_dy(i,j,k) = dfc_dy_aux(j)
            end do	  
          end do
        end do
      !MS$ELSE
        dfc_dy(:,:,:) = 0.d0
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_y_cp_65th_WN						

ccccc **********************************************************************
ccccc der_y_cp_65th_NN: This routine calculated the first derivatives 
ccccc                   in the y-direction using the 6th and 5th order 
ccccc                   compact finite difference. Neumann boundary 
ccccc                   condition was used to simulated the wall.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_65th_NN(dy,imax,jmax,kmax,mtrc,fc,dfc_dy)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc(jmax),fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),
     &  dfc_dy_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      !MS$IF (tp_dim.EQ.2).OR.(tp_dim.EQ.3)
        do k=1,kmax        
          do i=1,imax
            call rhs_y_65th_NN(dy,imax,jmax,kmax,i,k,fc,rhs)
            call lhs_d1_65th_N(aa,bb,cc,mtrc,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,dfc_dy_aux)   

            do j=1,jmax
              dfc_dy(i,j,k) = dfc_dy_aux(j)
            end do	  
          end do
        end do
      !MS$ELSE
        dfc_dy(:,:,:) = 0.d0
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_y_cp_65th_NN									

ccccc **********************************************************************
ccccc der_z_cp_6th_P: This routine calculated the first derivatives 
ccccc                 in the z-direction using the 6th order centered 
ccccc                 compact finite difference. Periodic boundary 
ccccc                 condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_z_cp_6th_P(dz,imax,jmax,kmax,fc,dfc_dz)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dz,fc(imax,jmax,kmax),dfc_dz(imax,jmax,kmax),
     &  dfc_dz_aux(kmax),aa(kmax),bb(kmax),cc(kmax),rhs(kmax)

      !MS$IF (tp_dim.EQ.3)
        do j=1,jmax
          do i=1,imax
            call rhs_z_6th_P(dz,imax,jmax,kmax,i,j,fc,rhs)
            call lhs_d1_6th_P(aa,bb,cc,kmax)
            call trid_cyclic(aa,bb,cc,1.d0,1.d0,kmax,rhs,dfc_dz_aux)
 
            do k=1,kmax-1
              dfc_dz(i,j,k) = dfc_dz_aux(k)
            end do
            dfc_dz(i,j,kmax) = dfc_dz(i,j,1)
          end do
        end do
      !MS$ELSE
        dfc_dz(:,:,:) = 0.d0
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_z_cp_6th_P

ccccc **********************************************************************
ccccc rhs_x_6th_P: This routine storing the right hand side of linear system 
ccccc              of first derivatives in the x-direction. The 6th order 
ccccc              compact finite difference schemes was used. Periodic 
ccccc              boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_x_6th_P(dx,imax,jmax,kmax,j,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,imax_1,jmax,kmax
      real*8 dx,fc(imax,jmax,kmax),r(imax)

      imax_1=imax-1
      r(1)=(-fc(imax_1-1,j,k)-28.d0*fc(imax_1,j,k)+
     &  28.d0*fc(2,j,k)+fc(3,j,k))/(12.d0*dx)
      r(2)=(-fc(imax_1,j,k)-28.d0*fc(1,j,k)+
     &  28.d0*fc(3,j,k)+fc(4,j,k))/(12.d0*dx)
      do i=3,imax_1-2
        r(i)=(-fc(i-2,j,k)-28.d0*fc(i-1,j,k)+
     &    28.d0*fc(i+1,j,k)+fc(i+2,j,k))/(12.d0*dx)
      end do
      r(imax_1-1)=(-fc(imax_1-3,j,k)-28.d0*fc(imax_1-2,j,k)+
     &  28.d0*fc(imax_1,j,k)+fc(1,j,k))/(12.d0*dx)
      r(imax_1)=(-fc(imax_1-2,j,k)-28.d0*fc(imax_1-1,j,k)+
     &  28.d0*fc(1,j,k)+fc(2,j,k))/(12.d0*dx)

      !The routine trid and trid_cyclic was modified to despise 
      !the last point of computational domain for periodic 
      !boundary conditions. Therefore the last point of vector  
      !r is equal zero as follow below.
      r(imax)=0.d0

      RETURN
      END SUBROUTINE rhs_x_6th_P

ccccc **********************************************************************
ccccc rhs_y_6th_P: This routine storing the right hand side of linear system 
ccccc              of first derivatives in the y-direction. The 6th order 
ccccc              compact finite difference schemes was used. Periodic
ccccc              boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_y_6th_P(dy,imax,jmax,kmax,i,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,jmax_1,kmax
      real*8 dy,fc(imax,jmax,kmax),r(jmax)
      
      jmax_1=jmax-1
      r(1)=(-fc(i,jmax_1-1,k)-28.d0*fc(i,jmax_1,k)+
     &  28.d0*fc(i,2,k)+fc(i,3,k))/(12.d0*dy)
      r(2)=(-fc(i,jmax_1,k)-28.d0*fc(i,1,k)+
     &  28.d0*fc(i,3,k)+fc(i,4,k))/(12.d0*dy)
      do j=3,jmax_1-2
        r(j)=(-fc(i,j-2,k)-28.d0*fc(i,j-1,k)+
     &    28.d0*fc(i,j+1,k)+fc(i,j+2,k))/(12.d0*dy)
      end do
      r(jmax_1-1)=(-fc(i,jmax_1-3,k)-28.d0*fc(i,jmax_1-2,k)+
     &  28.d0*fc(i,jmax_1,k)+fc(i,1,k))/(12.d0*dy)
      r(jmax_1)=(-fc(i,jmax_1-2,k)-28.d0*fc(i,jmax_1-1,k)+
     &  28.d0*fc(i,1,k)+fc(i,2,k))/(12.d0*dy)

      !The routine trid and trid_cyclic was modified to despise 
      !the last point of computational domain for periodic 
      !boundary conditions. Therefore the last point of vector  
      !r is equal zero as follow below.
      r(jmax)=0.d0
 
      RETURN
      END SUBROUTINE rhs_y_6th_P

ccccc **********************************************************************
ccccc rhs_y_62th_D: This routine storing the right hand side of linear 
ccccc               system of first derivatives in the y-direction. The 
ccccc               6th and 2nd order compact finite difference schemes 
ccccc               was used. Dirichlet boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_y_62th_D(dy,imax,jmax,kmax,i,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,fc(imax,jmax,kmax),r(jmax) 
	            
      r(1)=(-3.d0*fc(i,1,k)+4.d0*fc(i,2,k)-fc(i,3,k))/(2.d0*dy)
      r(2)=(-fc(i,1,k)+fc(i,3,k))/(2.d0*dy)      
      do j=3,jmax-2 
        r(j)=(-fc(i,j-2,k)-28.d0*fc(i,j-1,k)
     &    +28.d0*fc(i,j+1,k)+fc(i,j+2,k))/(12.d0*dy)
      end do	
      r(jmax-1)=(-fc(i,jmax-2,k)+fc(i,jmax,k))/(2.d0*dy)      
      r(jmax)=(fc(i,jmax-2,k)-4.d0*fc(i,jmax-1,k)
     &  +3.d0*fc(i,jmax,k))/(2.d0*dy)

      RETURN
      END SUBROUTINE rhs_y_62th_D

ccccc **********************************************************************
ccccc rhs_y_65th_WD: This routine storing the right hand side of linear 
ccccc                system of first derivatives in the y-direction. The 
ccccc                6th and 5th order compact finite difference schemes 
ccccc                was used. Dirichlet and Neumann boundary condition 
ccccc                was adopted to simulated the wall.
ccccc **********************************************************************
      SUBROUTINE rhs_y_65th_WD(dy,imax,jmax,kmax,i,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,fc(imax,jmax,kmax),r(jmax) 
	            
      r(1)=(-74.d0*fc(i,1,k)+16.d0*fc(i,2,k)+72.d0*fc(i,3,k)
     &  -16.d0*fc(i,4,k)+2.d0*fc(i,5,k))/(24.d0*dy)
      r(2)=(-406.d0*fc(i,1,k)-300.d0*fc(i,2,k)
     &  +760.d0*fc(i,3,k)-80.d0*fc(i,4,k)
     &  +30.d0*fc(i,5,k)-4.d0*fc(i,6,k))/(120.d0*dy)
      do j=3,jmax-2 
        r(j)=(-fc(i,j-2,k)-28.d0*fc(i,j-1,k)
     &    +28.d0*fc(i,j+1,k)+fc(i,j+2,k))/(12.d0*dy)
      end do	
      r(jmax-1)=(-fc(i,jmax-2,k)+fc(i,jmax,k))/(2.d0*dy)      
      r(jmax)=(1.d0*fc(i,jmax-2,k)-4.d0*fc(i,jmax-1,k)
     &  +3.d0*fc(i,jmax,k))/(2.d0*dy)

      RETURN
      END SUBROUTINE rhs_y_65th_WD

ccccc **********************************************************************
ccccc rhs_y_65th_WN: This routine storing the right hand side of linear 
ccccc                system of first derivatives in the y-direction. The 
ccccc                6th and 5th order compact finite difference schemes 
ccccc                was used. Dirichlet and Neumann boundary condition 
ccccc                was adopted to simulated the wall.
ccccc **********************************************************************
      SUBROUTINE rhs_y_65th_WN(dy,imax,jmax,kmax,i,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,fc(imax,jmax,kmax),r(jmax) 
	            
      r(1)=(-74.d0*fc(i,1,k)+16.d0*fc(i,2,k)+72.d0*fc(i,3,k)
     &  -16.d0*fc(i,4,k)+2.d0*fc(i,5,k))/(24.d0*dy)
      r(2)=(-406.d0*fc(i,1,k)-300.d0*fc(i,2,k)
     &  +760.d0*fc(i,3,k)-80.d0*fc(i,4,k)
     &  +30.d0*fc(i,5,k)-4.d0*fc(i,6,k))/(120.d0*dy)
      do j=3,jmax-2 
        r(j)=(-fc(i,j-2,k)-28.d0*fc(i,j-1,k)
     &    +28.d0*fc(i,j+1,k)+fc(i,j+2,k))/(12.d0*dy)
      end do	
      r(jmax-1)=(-fc(i,jmax-2,k)+fc(i,jmax,k))/(2.d0*dy)      
      r(jmax)=0.d0

      RETURN
      END SUBROUTINE rhs_y_65th_WN

ccccc **********************************************************************
ccccc rhs_y_62th_N: This routine storing the right hand side of linear 
ccccc               system of first derivatives in the y-direction. The 
ccccc               6th and 2nd order compact finite difference schemes 
ccccc               was used. Neumann boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_y_62th_N(dy,imax,jmax,kmax,i,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,fc(imax,jmax,kmax),r(jmax)
  
      r(1)=0.d0
      r(2)=(-fc(i,1,k)+fc(i,3,k))/(2.d0*dy)      
      do j=3,jmax-2 
        r(j)=(-fc(i,j-2,k)-28.d0*fc(i,j-1,k)
     &    +28.d0*fc(i,j+1,k)+fc(i,j+2,k))/(12.d0*dy)
      end do	
      r(jmax-1)=(-fc(i,jmax-2,k)+fc(i,jmax,k))/(2.d0*dy)      
      r(jmax)=0.d0

      RETURN
      END SUBROUTINE rhs_y_62th_N

ccccc **********************************************************************
ccccc rhs_y_65th_NN: This routine storing the right hand side of linear 
ccccc                system of first derivatives in the y-direction. The 
ccccc                6th and 5th order compact finite difference schemes 
ccccc                was used. Neumann boundary condition 
ccccc                was adopted to simulated the wall.
ccccc **********************************************************************
      SUBROUTINE rhs_y_65th_NN(dy,imax,jmax,kmax,i,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,fc(imax,jmax,kmax),r(jmax) 
	            
      r(1)=0.d0
      r(2)=(-132.d0*fc(i,1,k)-900.d0*fc(i,2,k)
     &  +1360.d0*fc(i,3,k)-480.d0*fc(i,4,k)
     &  +180.d0*fc(i,5,k)-28.d0*fc(i,6,k))/(120.d0*dy)
      do j=3,jmax-2 
        r(j)=(-fc(i,j-2,k)-28.d0*fc(i,j-1,k)
     &    +28.d0*fc(i,j+1,k)+fc(i,j+2,k))/(12.d0*dy)
      end do	
      r(jmax-1)=(-fc(i,jmax-2,k)+fc(i,jmax,k))/(2.d0*dy)      
      r(jmax)=0.d0

      RETURN
      END SUBROUTINE rhs_y_65th_NN

ccccc **********************************************************************
ccccc rhs_z_6th_P: This routine storing the right hand side of linear system 
ccccc              of first derivatives in the z-direction. The 6th order 
ccccc              compact finite difference schemes was used. Periodic 
ccccc              boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_z_6th_P(dz,imax,jmax,kmax,i,j,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax,kmax_1
      real*8 dz,fc(imax,jmax,kmax),r(kmax)

      kmax_1=kmax-1
      r(1)=(-fc(i,j,kmax_1-1)-28.d0*fc(i,j,kmax_1)+
     &  28.d0*fc(i,j,2)+fc(i,j,3))/(12.d0*dz)
      r(2)=(-fc(i,j,kmax_1)-28.d0*fc(i,j,1)+
     &  28.d0*fc(i,j,3)+fc(i,j,4))/(12.d0*dz)
      do k=3,kmax_1-2
        r(k)=(-fc(i,j,k-2)-28.d0*fc(i,j,k-1)+
     &    28.d0*fc(i,j,k+1)+fc(i,j,k+2))/(12.d0*dz)
      end do
      r(kmax_1-1)=(-fc(i,j,kmax_1-3)-28.d0*fc(i,j,kmax_1-2)+
     &  28.d0*fc(i,j,kmax_1)+fc(i,j,1))/(12.d0*dz)
      r(kmax_1)=(-fc(i,j,kmax_1-2)-28.d0*fc(i,j,kmax_1-1)+
     &  28.d0*fc(i,j,1)+fc(i,j,2))/(12.d0*dz)
      !The routine trid and trid_cyclic was modified to despise 
      !the last point of computational domain for periodic 
      !boundary conditions. Therefore the last point of vector  
      !r is equal zero as follow below.
      r(kmax)=0.d0

      RETURN
      END SUBROUTINE rhs_z_6th_P

ccccc **********************************************************************
ccccc lhs_d1_6th_P: This routine storing the left hand side of linear system 
ccccc               of first derivatives. The 6th order compact finite 
ccccc               difference schemes was used. Periodic boundary 
ccccc               conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE lhs_d1_6th_P(aa,bb,cc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax,lmax_1
      real*8 aa(lmax),bb(lmax),cc(lmax)
  
      lmax_1=lmax-1

      aa(1)=0.d0 
      bb(1)=3.d0
      cc(1)=1.d0
  
      do l=2,lmax_1-1 
        aa(l)=1.d0
        bb(l)=3.d0
        cc(l)=1.d0
      end do
    
      aa(lmax_1)=1.d0 
      bb(lmax_1)=3.d0
      cc(lmax_1)=0.d0
    
      !The routine trid and trid_cyclic was modified to despise 
      !the last point of computational domain for periodic 
      !boundary conditions. Therefore the last point of vector  
      !aa, bb and cc is equal zero as follow below.
      aa(lmax)=0.d0 
      bb(lmax)=0.d0
      cc(lmax)=0.d0
    		
      RETURN
      END SUBROUTINE lhs_d1_6th_P

ccccc **********************************************************************
ccccc lhs_d1_62th_D: This routine storing the left hand side of linear 
ccccc                system of first derivatives. The 6th and 2nd order  
ccccc                compact finite-difference schemes was used. Dirichlet  
ccccc                boundary conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE lhs_d1_62th_D(aa,bb,cc,mtrc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc(lmax)

      aa(1)=0.d0
      bb(1)=1.d0/mtrc(1)
      cc(1)=0.d0

      aa(2)=0.d0
      bb(2)=1.d0/mtrc(2)
      cc(2)=0.d0

      do l=3,lmax-2
        aa(l)=1.d0/mtrc(l-1)*1.d0
        bb(l)=1.d0/mtrc(l  )*3.d0
        cc(l)=1.d0/mtrc(l+1)*1.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0/mtrc(lmax-1)
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0/mtrc(lmax)
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d1_62th_D

ccccc **********************************************************************
ccccc lhs_d1_62th_N: This routine storing the left hand side of linear 
ccccc                system of first derivatives. The 6th order compact 
ccccc                finite-difference schemes was used. The Neumann 
ccccc                boundary conditions was adopted. A stretching is 
ccccc                used in this schemes.
ccccc **********************************************************************
      SUBROUTINE lhs_d1_62th_N(aa,bb,cc,mtrc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc(lmax)

      aa(1)=0.d0
      bb(1)=1.d0/mtrc(1)
      cc(1)=0.d0

      aa(2)=0.d0 
      bb(2)=1.d0/mtrc(2)
      cc(2)=0.d0

      do l=3,lmax-2
        aa(l)=1.d0/mtrc(l-1)*1.d0
        bb(l)=1.d0/mtrc(l  )*3.d0
        cc(l)=1.d0/mtrc(l+1)*1.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0/mtrc(lmax-1)
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0/mtrc(lmax)
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d1_62th_N

ccccc **********************************************************************
ccccc lhs_d1_65th_W: This routine storing the left hand side of linear 
ccccc                system of first derivatives. The 6th and 5th order  
ccccc                compact finite difference schemes was used. Dirichlet  
ccccc                and Neumann boundary conditions was adopted to 
ccccc                simulate the wall.
ccccc **********************************************************************
      SUBROUTINE lhs_d1_65th_W(aa,bb,cc,mtrc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc(lmax)

      aa(1)=0.d0
      bb(1)=1.d0/mtrc(1)*1.d0
      cc(1)=1.d0/mtrc(2)*4.d0

      aa(2)=1.d0/mtrc(1)*1.d0
      bb(2)=1.d0/mtrc(2)*6.d0
      cc(2)=1.d0/mtrc(3)*2.d0

      do l=3,lmax-2
        aa(l)=1.d0/mtrc(l-1)*1.d0
        bb(l)=1.d0/mtrc(l  )*3.d0
        cc(l)=1.d0/mtrc(l+1)*1.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0/mtrc(lmax-1)*1.d0
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0/mtrc(lmax)*1.d0
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d1_65th_W

ccccc **********************************************************************
ccccc lhs_d1_65th_N: This routine storing the left hand side of linear 
ccccc                system of first derivatives. The 6th and 5th order  
ccccc                compact finite difference schemes was used. Neumann 
ccccc                and Dirichlet boundary conditions was adopted to 
ccccc                simulate the wall.
ccccc **********************************************************************
      SUBROUTINE lhs_d1_65th_N(aa,bb,cc,mtrc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc(lmax)

      aa(1)=0.d0
      bb(1)=1.d0/mtrc(1)*1.d0
      cc(1)=0.d0

      aa(2)=0.d0                       
      bb(2)=1.d0/mtrc(2)*6.d0
      cc(2)=1.d0/mtrc(3)*2.d0

      do l=3,lmax-2
        aa(l)=1.d0/mtrc(l-1)*1.d0
        bb(l)=1.d0/mtrc(l  )*3.d0
        cc(l)=1.d0/mtrc(l+1)*1.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0/mtrc(lmax-1)*1.d0
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0/mtrc(lmax)*1.d0
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d1_65th_N

ccccc **********************************************************************
ccccc coeffs_d1_ex_5th: This routine gets the coefficients for the 5th  
ccccc                   order stencils for first derivatives.
ccccc **********************************************************************
      SUBROUTINE coeffs_d1_ex_5th(x1,x2,x3,x4,x5,x6,x,which,coeffs)

      IMPLICIT NONE

      integer which
      real*8 x1,x2,x3,x4,x5,x6,x,coeffs

      if (which.EQ.1) then
        coeffs = 1.d0
     &    /(x1-x2)*(x-x3)/(x1-x3)*(x-x4)/(x1-x4)*(x-x5)/(x1-x5)*(x-x6)/(x1-x6)
     &    +(x-x2)/(x1-x2)/(x1-x3)*(x-x4)/(x1-x4)*(x-x5)/(x1-x5)*(x-x6)/(x1-x6)
     &    +(x-x2)/(x1-x2)*(x-x3)/(x1-x3)/(x1-x4)*(x-x5)/(x1-x5)*(x-x6)/(x1-x6)
     &    +(x-x2)/(x1-x2)*(x-x3)/(x1-x3)*(x-x4)/(x1-x4)/(x1-x5)*(x-x6)/(x1-x6)
     &    +(x-x2)/(x1-x2)*(x-x3)/(x1-x3)*(x-x4)/(x1-x4)*(x-x5)/(x1-x5)/(x1-x6)
      else if (which.EQ.2) then
        coeffs = 1.d0
     &    /(x2-x1)*(x-x3)/(x2-x3)*(x-x4)/(x2-x4)*(x-x5)/(x2-x5)*(x-x6)/(x2-x6)
     &    +(x-x1)/(x2-x1)/(x2-x3)*(x-x4)/(x2-x4)*(x-x5)/(x2-x5)*(x-x6)/(x2-x6)
     &    +(x-x1)/(x2-x1)*(x-x3)/(x2-x3)/(x2-x4)*(x-x5)/(x2-x5)*(x-x6)/(x2-x6)
     &    +(x-x1)/(x2-x1)*(x-x3)/(x2-x3)*(x-x4)/(x2-x4)/(x2-x5)*(x-x6)/(x2-x6)
     &    +(x-x1)/(x2-x1)*(x-x3)/(x2-x3)*(x-x4)/(x2-x4)*(x-x5)/(x2-x5)/(x2-x6)
      else if (which.EQ.3) then
        coeffs = 1.d0
     &    /(x3-x1)*(x-x2)/(x3-x2)*(x-x4)/(x3-x4)*(x-x5)/(x3-x5)*(x-x6)/(x3-x6)
     &    +(x-x1)/(x3-x1)/(x3-x2)*(x-x4)/(x3-x4)*(x-x5)/(x3-x5)*(x-x6)/(x3-x6)
     &    +(x-x1)/(x3-x1)*(x-x2)/(x3-x2)/(x3-x4)*(x-x5)/(x3-x5)*(x-x6)/(x3-x6)
     &    +(x-x1)/(x3-x1)*(x-x2)/(x3-x2)*(x-x4)/(x3-x4)/(x3-x5)*(x-x6)/(x3-x6)
     &    +(x-x1)/(x3-x1)*(x-x2)/(x3-x2)*(x-x4)/(x3-x4)*(x-x5)/(x3-x5)/(x3-x6)
      else if (which.EQ.4) then
        coeffs = 1.d0
     &    /(x4-x1)*(x-x2)/(x4-x2)*(x-x3)/(x4-x3)*(x-x5)/(x4-x5)*(x-x6)/(x4-x6)
     &    +(x-x1)/(x4-x1)/(x4-x2)*(x-x3)/(x4-x3)*(x-x5)/(x4-x5)*(x-x6)/(x4-x6)
     &    +(x-x1)/(x4-x1)*(x-x2)/(x4-x2)/(x4-x3)*(x-x5)/(x4-x5)*(x-x6)/(x4-x6)
     &    +(x-x1)/(x4-x1)*(x-x2)/(x4-x2)*(x-x3)/(x4-x3)/(x4-x5)*(x-x6)/(x4-x6)
     &    +(x-x1)/(x4-x1)*(x-x2)/(x4-x2)*(x-x3)/(x4-x3)*(x-x5)/(x4-x5)/(x4-x6)
      else if(which.EQ.5) then
        coeffs = 1.d0
     &    /(x5-x1)*(x-x2)/(x5-x2)*(x-x3)/(x5-x3)*(x-x4)/(x5-x4)*(x-x6)/(x5-x6)
     &    +(x-x1)/(x5-x1)/(x5-x2)*(x-x3)/(x5-x3)*(x-x4)/(x5-x4)*(x-x6)/(x5-x6)
     &    +(x-x1)/(x5-x1)*(x-x2)/(x5-x2)/(x5-x3)*(x-x4)/(x5-x4)*(x-x6)/(x5-x6)
     &    +(x-x1)/(x5-x1)*(x-x2)/(x5-x2)*(x-x3)/(x5-x3)/(x5-x4)*(x-x6)/(x5-x6)
     &    +(x-x1)/(x5-x1)*(x-x2)/(x5-x2)*(x-x3)/(x5-x3)*(x-x4)/(x5-x4)/(x5-x6)
      else if(which.EQ.6) then
        coeffs = 1.d0
     &    /(x6-x1)*(x-x2)/(x6-x2)*(x-x3)/(x6-x3)*(x-x4)/(x6-x4)*(x-x5)/(x6-x5)
     &    +(x-x1)/(x6-x1)/(x6-x2)*(x-x3)/(x6-x3)*(x-x4)/(x6-x4)*(x-x5)/(x6-x5)
     &    +(x-x1)/(x6-x1)*(x-x2)/(x6-x2)/(x6-x3)*(x-x4)/(x6-x4)*(x-x5)/(x6-x5)
     &    +(x-x1)/(x6-x1)*(x-x2)/(x6-x2)*(x-x3)/(x6-x3)/(x6-x4)*(x-x5)/(x6-x5)
     &    +(x-x1)/(x6-x1)*(x-x2)/(x6-x2)*(x-x3)/(x6-x3)*(x-x4)/(x6-x4)/(x6-x5)
      end if

      RETURN
      END SUBROUTINE coeffs_d1_ex_5th

ccccc **********************************************************************
ccccc coeffs_d1_ex_6th: This routine gets the coefficients for the 6th  
ccccc                   order stencils for first derivatives.
ccccc **********************************************************************
      SUBROUTINE coeffs_d1_ex_6th(x1,x2,x3,x4,x5,x6,x7,x,which,coeffs)

      IMPLICIT NONE

      integer which
      real*8 x1,x2,x3,x4,x5,x6,x7,x,coeffs

      if (which.EQ.1) then
        coeffs = 
     &   ((x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x6))/             
     &    ((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)*(x1-x7))+ 
     &   ((x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x7))/                  
     &    ((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)*(x1-x7))+
     &   ((x-x2)*(x-x3)*(x-x4)*(x-x6)*(x-x7))/                  
     &    ((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)*(x1-x7))+
     &   ((x-x2)*(x-x3)*(x-x5)*(x- x6)*(x- x7))/                  
     &    ((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)*(x1-x7))+
     &   ((x-x2)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/                  
     &    ((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)*(x1-x7))+
     &   ((x-x3)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/                  
     &    ((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)*(x1-x7))
      else if (which.EQ.2) then
        coeffs =
     &   ((x-x1)*(x-x3)*(x-x4)*(x-x5)*(x-x6))/                   
     &    ((-x1+x2)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)*(x2-x7))+
     &   ((x-x1)*(x-x3)*(x-x4)*(x-x5)*(x-x7))/                   
     &    ((-x1+x2)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)*(x2-x7))+
     &   ((x-x1)*(x-x3)*(x-x4)*(x-x6)*(x-x7))/                   
     &    ((-x1+x2)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)*(x2-x7))+
     &   ((x-x1)*(x-x3)*(x-x5)*(x-x6)*(x-x7))/                   
     &    ((-x1+x2)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)*(x2-x7))+
     &   ((x-x1)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/                   
     &    ((-x1+x2)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)*(x2-x7))+
     &   ((x-x3)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/                   
     &    ((-x1+x2)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)*(x2-x7))
      else if (which.EQ.3) then
        coeffs =
     &   ((x-x1)*(x-x2)*(x-x4)*(x-x5)*(x-x6))/                    
     &    ((-x1+x3)*(-x2+x3)*(x3-x4)*(x3-x5)*(x3-x6)*(x3-x7))+
     &   ((x-x1)*(x-x2)*(x-x4)*(x-x5)*(x-x7))/                    
     &    ((-x1+x3)*(-x2+x3)*(x3-x4)*(x3-x5)*(x3-x6)*(x3-x7))+
     &   ((x-x1)*(x-x2)*(x-x4)*(x-x6)*(x-x7))/                    
     &    ((-x1+x3)*(-x2+x3)*(x3-x4)*(x3-x5)*(x3-x6)*(x3-x7))+
     &   ((x-x1)*(x-x2)*(x-x5)*(x-x6)*(x-x7))/                    
     &    ((-x1+x3)*(-x2+x3)*(x3-x4)*(x3-x5)*(x3-x6)*(x3-x7))+
     &   ((x-x1)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/                    
     &    ((-x1+x3)*(-x2+x3)*(x3-x4)*(x3-x5)*(x3-x6)*(x3-x7))+
     &   ((x-x2)*(x-x4)*(x-x5)*(x-x6)*(x-x7))/                    
     &    ((-x1+x3)*(-x2+x3)*(x3-x4)*(x3-x5)*(x3-x6)*(x3-x7))
      else if (which.EQ.4) then
        coeffs =
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x5)*(x-x6))/                     
     &    ((-x1+x4)*(-x2+x4)*(-x3+x4)*(x4-x5)*(x4-x6)*(x4-x7))+
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x5)*(x-x7))/                     
     &    ((-x1+x4)*(-x2+x4)*(-x3+x4)*(x4-x5)*(x4-x6)*(x4-x7))+
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x6)*(x-x7))/                     
     &    ((-x1+x4)*(-x2+x4)*(-x3+x4)*(x4-x5)*(x4-x6)*(x4-x7))+
     &   ((x-x1)*(x-x2)*(x-x5)*(x-x6)*(x-x7))/                     
     &    ((-x1+x4)*(-x2+x4)*(-x3+x4)*(x4-x5)*(x4-x6)*(x4-x7))+
     &   ((x-x1)*(x-x3)*(x-x5)*(x-x6)*(x-x7))/                     
     &    ((-x1+x4)*(-x2+x4)*(-x3+x4)*(x4-x5)*(x4-x6)*(x4-x7))+
     &   ((x-x2)*(x-x3)*(x-x5)*(x-x6)*(x-x7))/                     
     &    ((-x1+x4)*(-x2+x4)*(-x3+x4)*(x4-x5)*(x4-x6)*(x4-x7))
      else if(which.EQ.5) then
        coeffs =
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x6))/                      
     &    ((-x1+x5)*(-x2+x5)*(-x3+x5)*(-x4+x5)*(x5-x6)*(x5-x7))+
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x7))/                      
     &    ((-x1+x5)*(-x2+x5)*(-x3+x5)*(-x4+x5)*(x5-x6)*(x5-x7))+
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x6)*(x-x7))/                      
     &    ((-x1+x5)*(-x2+x5)*(-x3+x5)*(-x4+x5)*(x5-x6)*(x5-x7))+
     &   ((x-x1)*(x-x2)*(x-x4)*(x-x6)*(x-x7))/                      
     &    ((-x1+x5)*(-x2+x5)*(-x3+x5)*(-x4+x5)*(x5-x6)*(x5-x7))+
     &   ((x-x1)*(x-x3)*(x-x4)*(x-x6)*(x-x7))/                      
     &    ((-x1+x5)*(-x2+x5)*(-x3+x5)*(-x4+x5)*(x5-x6)*(x5-x7))+
     &   ((x-x2)*(x-x3)*(x-x4)*(x-x6)*(x-x7))/                      
     &    ((-x1+x5)*(-x2+x5)*(-x3+x5)*(-x4+x5)*(x5-x6)*(x5-x7)) 
      else if(which.EQ.6) then
        coeffs =
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x5))/                       
     &    ((-x1+x6)*(-x2+x6)*(-x3+x6)*(-x4+x6)*(-x5+x6)*(x6-x7))+
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x7))/                       
     &    ((-x1+x6)*(-x2+x6)*(-x3+x6)*(-x4+x6)*(-x5+x6)*(x6-x7))+
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x5)*(x-x7))/                       
     &    ((-x1+x6)*(-x2+x6)*(-x3+x6)*(-x4+x6)*(-x5+x6)*(x6-x7))+
     &   ((x-x1)*(x-x2)*(x-x4)*(x-x5)*(x-x7))/                       
     &    ((-x1+x6)*(-x2+x6)*(-x3+x6)*(-x4+x6)*(-x5+x6)*(x6-x7))+
     &   ((x-x1)*(x-x3)*(x-x4)*(x-x5)*(x-x7))/                       
     &    ((-x1+x6)*(-x2+x6)*(-x3+x6)*(-x4+x6)*(-x5+x6)*(x6-x7))+
     &   ((x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x7))/                       
     &    ((-x1+x6)*(-x2+x6)*(-x3+x6)*(-x4+x6)*(-x5+x6)*(x6-x7))
      else if(which.EQ.7) then
        coeffs =
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x5))/                        
     &    ((-x1+x7)*(-x2+x7)*(-x3+x7)*(-x4+x7)*(-x5+x7)*(-x6+x7))+
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x6))/                        
     &    ((-x1+x7)*(-x2+x7)*(-x3+x7)*(-x4+x7)*(-x5+x7)*(-x6+x7))+
     &   ((x-x1)*(x-x2)*(x-x3)*(x-x5)*(x-x6))/                        
     &    ((-x1+x7)*(-x2+x7)*(-x3+x7)*(-x4+x7)*(-x5+x7)*(-x6+x7))+
     &   ((x-x1)*(x-x2)*(x-x4)*(x-x5)*(x-x6))/                        
     &    ((-x1+x7)*(-x2+x7)*(-x3+x7)*(-x4+x7)*(-x5+x7)*(-x6+x7))+
     &   ((x-x1)*(x-x3)*(x-x4)*(x-x5)*(x-x6))/                        
     &    ((-x1+x7)*(-x2+x7)*(-x3+x7)*(-x4+x7)*(-x5+x7)*(-x6+x7))+
     &   ((x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x6))/                        
     &    ((-x1+x7)*(-x2+x7)*(-x3+x7)*(-x4+x7)*(-x5+x7)*(-x6+x7))
      end if

      RETURN
      END SUBROUTINE coeffs_d1_ex_6th


ccccc **********************************************************************
ccccc coeffs_d2_ex_6th: This routine gets the coefficients for the 6th  
ccccc                   order stencils for the second derivatives.
ccccc **********************************************************************
      SUBROUTINE coeffs_d2_ex_6th(x1,x2,x3,x4,x5,x6,x7,x,which,coeffs)

      IMPLICIT NONE

      integer which
      real*8 x1,x2,x3,x4,x5,x6,x7,x,coeffs,s1,s2

      if (which.EQ.1) then
        coeffs = 
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x5)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x6)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x4)*(x-x5)*(x-x6)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x3)*(x-x4)*(x-x5)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x3)*(x-x4)*(x-x5)*(x-x6)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x5)*(x-x6)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x4)*(x-x6)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x3)*(x-x4)*(x-x6)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x4)*(x-x5)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x6)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x5)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x5)*(x-x6)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x4)*(x-x5)*(x-x6)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
     &    +2.d0*(x-x3)*(x-x5)*(x-x6)*(x-x7)
     &      /(x1-x2)/(x1-x3)/(x1-x4)/(x1-x5)/(x1-x6)/(x1-x7)
      else if (which.EQ.2) then
        coeffs = 
     &    2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x5)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x6)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x6)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x5)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x5)*(x-x6)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x4)*(x-x5)*(x-x6)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x3)*(x-x5)*(x-x6)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x3)*(x-x4)*(x-x6)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x3)*(x-x4)*(x-x5)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x3)*(x-x4)*(x-x5)*(x-x6)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x5)*(x-x6)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x4)*(x-x6)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x4)*(x-x5)*(x-x7)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
     &    +2.d0*(x-x1)*(x-x4)*(x-x5)*(x-x6)
     &      /(x2-x1)/(x2-x3)/(x2-x4)/(x2-x5)/(x2-x6)/(x2-x7)
      else if (which.EQ.3) then
        coeffs = 
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x6)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x4)*(x-x5)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x4)*(x-x5)*(x-x6)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x6)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x5)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x5)*(x-x6)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x4)*(x-x5)*(x-x6)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x2)*(x-x5)*(x-x6)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x2)*(x-x4)*(x-x6)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x2)*(x-x4)*(x-x5)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x2)*(x-x4)*(x-x5)*(x-x6)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x5)*(x-x6)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x4)*(x-x6)*(x-x7)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x5)
     &      /(x3-x1)/(x3-x2)/(x3-x4)/(x3-x5)/(x3-x6)/(x3-x7)
      else if (which.EQ.4) then
        coeffs = 
     &    +2.d0*(x-x1)*(x-x2)*(x-x6)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x5)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x5)*(x-x6)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x6)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x5)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x5)*(x-x6)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x5)*(x-x6)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x6)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x5)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x3)*(x-x5)*(x-x6)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x2)*(x-x5)*(x-x6)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x6)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x5)*(x-x7)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x5)*(x-x6)
     &      /(x4-x1)/(x4-x2)/(x4-x3)/(x4-x5)/(x4-x6)/(x4-x7)
      else if (which.EQ.5) then
        coeffs = 
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x7)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x6)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x4)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x3)*(x-x4)*(x-x6)*(x-x7)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x2)*(x-x4)*(x-x6)*(x-x7)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x6)*(x-x7)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x7) 
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x6)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x4)*(x-x6)*(x-x7) 
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x6)*(x-x7)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x7)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x6)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x6)*(x-x7)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x7)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x6)
     &      /(x5-x1)/(x5-x2)/(x5-x3)/(x5-x4)/(x5-x6)/(x5-x7)
      else if (which.EQ.6) then
        coeffs = 
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x4)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x5)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7) 
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x5)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x1)*(x-x4)*(x-x5)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x5)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x5)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x5)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x3)*(x-x4)*(x-x5)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x2)*(x-x4)*(x-x5)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x5)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x7)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x5)
     &      /(x6-x1)/(x6-x2)/(x6-x3)/(x6-x4)/(x6-x5)/(x6-x7)
      else if (which.EQ.7) then
        coeffs = 
     &    +2.d0*(x-x1)*(x-x2)*(x-x5)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x2)*(x-x4)*(x-x5)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x3)*(x-x4)*(x-x5)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x2)*(x-x4)*(x-x5)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x2)*(x-x3)*(x-x5)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x2)*(x-x3)*(x-x4)*(x-x5)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x4)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x2)*(x-x3)*(x-x5)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x4)*(x-x5)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x3)*(x-x5)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x6)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
     &    +2.d0*(x-x1)*(x-x3)*(x-x4)*(x-x5)
     &      /(x7-x1)/(x7-x2)/(x7-x3)/(x7-x4)/(x7-x5)/(x7-x6)
      end if 

      RETURN
      END SUBROUTINE coeffs_d2_ex_6th

ccccc **********************************************************************
ccccc stx_d1_ex_6th_D: This routine calculates the stencils for 1st 
ccccc                  derivatives in x direction with a 6th order 
ccccc                  explicit schemes. The non-periodic boundary 
ccccc                  conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE stx_d1_ex_6th_D(my_rank,pro,imax,tp_s,x,stx)

      IMPLICIT NONE

      integer i,imax,my_rank,pro,tp_s
      real*8 x(imax),stx(imax,7)

      if (my_rank.EQ.0) then 
        do i=1,3
          if ((tp_s.EQ.4).OR.(tp_s.EQ.5)) then   
            call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),1,stx(i,1))
            call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),2,stx(i,2))
            call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),3,stx(i,3))
            call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),4,stx(i,4))
            call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),5,stx(i,5))
            call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),6,stx(i,6))
            call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),7,stx(i,7))
          else if (tp_s.EQ.6) then
            call coeffs_d1_ex_5th(x(1),x(2),x(3),x(4),x(5),x(6),x(i),1,stx(i,1))
            call coeffs_d1_ex_5th(x(1),x(2),x(3),x(4),x(5),x(6),x(i),2,stx(i,2))
            call coeffs_d1_ex_5th(x(1),x(2),x(3),x(4),x(5),x(6),x(i),3,stx(i,3))
            call coeffs_d1_ex_5th(x(1),x(2),x(3),x(4),x(5),x(6),x(i),4,stx(i,4))
            call coeffs_d1_ex_5th(x(1),x(2),x(3),x(4),x(5),x(6),x(i),5,stx(i,5))
            call coeffs_d1_ex_5th(x(1),x(2),x(3),x(4),x(5),x(6),x(i),6,stx(i,6))
            stx(i,7)=0.d0
          end if   
        end do
      end if 

      do i=4,imax-3
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),1,stx(i,1))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),2,stx(i,2))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),3,stx(i,3))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),4,stx(i,4))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),5,stx(i,5))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),6,stx(i,6))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),7,stx(i,7))   
      end do

      if (my_rank+1.EQ.pro) then
        do i=imax-2,imax
          if ((tp_s.EQ.4).OR.(tp_s.EQ.5)) then
            call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),1,stx(i,1))
            call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),2,stx(i,2))
            call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),3,stx(i,3))
            call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),4,stx(i,4))
            call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),5,stx(i,5))
            call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),6,stx(i,6))
            call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),7,stx(i,7))
          else if (tp_s.EQ.6) then
            stx(i,1) = 0.d0
            call coeffs_d1_ex_5th(x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),1,stx(i,2))
            call coeffs_d1_ex_5th(x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),2,stx(i,3))
            call coeffs_d1_ex_5th(x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),3,stx(i,4))
            call coeffs_d1_ex_5th(x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),4,stx(i,5))
            call coeffs_d1_ex_5th(x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),5,stx(i,6))
            call coeffs_d1_ex_5th(x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),6,stx(i,7))
          end if
        end do
      end if

      RETURN
      END SUBROUTINE stx_d1_ex_6th_D 

ccccc **********************************************************************
ccccc stx_d2_ex_6th_D: This routine calculates the stencils for 2nd 
ccccc                  derivatives in x direction with a 6th order 
ccccc                  explicit schemes. The non-periodic boundary 
ccccc                  conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE stx_d2_ex_6th_D(my_rank,pro,imax,tp_s,x,stxx)

      IMPLICIT NONE

      integer i,imax,my_rank,pro,tp_s
      real*8 x(imax),stxx(imax,7)

      if (my_rank.EQ.0) then
        do i=1,3
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),1,stxx(i,1))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),2,stxx(i,2))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),3,stxx(i,3))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),4,stxx(i,4))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),5,stxx(i,5))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),6,stxx(i,6))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(i),7,stxx(i,7))   
        end do
      end if

      do i=4,imax-3
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),1,stxx(i,1))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),2,stxx(i,2))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),3,stxx(i,3))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),4,stxx(i,4))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),5,stxx(i,5))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),6,stxx(i,6))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),x(i+3),x(i),7,stxx(i,7))   
      end do

      if (my_rank+1.EQ.pro) then
        do i=imax-2,imax
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),1,stxx(i,1))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),2,stxx(i,2))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),3,stxx(i,3))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),4,stxx(i,4))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),5,stxx(i,5))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),6,stxx(i,6))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),x(imax-2),x(imax-1),x(imax),x(i),7,stxx(i,7))
        end do
      end if 

      RETURN
      END SUBROUTINE stx_d2_ex_6th_D

ccccc **********************************************************************
ccccc der_x_DD: This routine calculated the first derivatives in the
ccccc           x direction. Periodic and Dirichlet boundary condition 
ccccc           was adopted in according to the problem simulated.
ccccc **********************************************************************
      SUBROUTINE der_x_DD(my_rank,pro,local_imax,fc,dfc_dx)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer local_imax,my_rank,pro
      real*8 fc(local_imax,jmax,kmax),dfc_dx(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3).OR.(tp_s.EQ.4)
        call der_x_cp_6th_P(dx,local_imax,jmax,kmax,fc,dfc_dx)
      !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)
        call der_x_ex_6th_D(my_rank,pro,local_imax,jmax,kmax,fc,dfc_dx)
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_x_DD

ccccc **********************************************************************
ccccc der_y_WD: This routine calculated the first derivatives in the  
ccccc           y direction. 
ccccc **********************************************************************
      SUBROUTINE der_y_WD(local_imax,fc,dfc_dy)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(jmax)

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
      !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_D(dy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_D(dcy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ENDIF         
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)
          call der_y_cp_65th_WD(dy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WD(dcy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ENDIF         
      !MS$ENDIF  

      RETURN
      END SUBROUTINE der_y_WD

ccccc **********************************************************************
ccccc der_y_WN: This routine calculated the first derivatives in the 
ccccc           y-direction. 
ccccc **********************************************************************
      SUBROUTINE der_y_WN(local_imax,fc,dfc_dy)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),
     &  dfc_dy(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
      !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_N(dy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_N(dcy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ENDIF 
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_65th_WN(dy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_WN(dcy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ENDIF 
      !MS$ENDIF

      RETURN
      END SUBROUTINE der_y_WN

ccccc **********************************************************************
ccccc der_y_NN: This routine calculated the first derivatives in the  
ccccc           y-direction. Neumann boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_y_NN(local_imax,fc,dfc_dy)

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dy(local_imax,jmax,kmax)

      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3)
        call der_y_cp_6th_P(dy,local_imax,jmax,kmax,fc,dfc_dy)
      !MS$ELSEIF (tp_s.EQ.4).OR.(tp_s.EQ.5)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_62th_N(dy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_62th_N(dcy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ENDIF         
      !MS$ELSEIF (tp_s.EQ.6)
        !MS$IF (stretching_in_y.EQ.0)     
          call der_y_cp_65th_NN(dy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ELSEIF (stretching_in_y.EQ.1) 
          call der_y_cp_65th_NN(dcy,local_imax,jmax,kmax,d1cy_dpy1,fc,dfc_dy)
        !MS$ENDIF         
      !MS$ENDIF  

      RETURN
      END SUBROUTINE der_y_NN

ccccc **********************************************************************
ccccc der_z_PP: This routine calculated the first derivatives in the
ccccc           z direction. Periodic boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_z_PP(local_imax,fc,dfc_dz)

      integer local_imax
      real*8 fc(local_imax,jmax,kmax),dfc_dz(local_imax,jmax,kmax)

      call der_z_cp_6th_P(dz,local_imax,jmax,kmax,fc,dfc_dz)

      RETURN
      END SUBROUTINE der_z_PP

ccccc **********************************************************************
ccccc trid: Solve the linear system (Revised on March, 7).
ccccc fct.: Factor (0 - Non Periodic; 1 - Periodic) 
ccccc **********************************************************************
      SUBROUTINE trid(fct,a,b,c,lmax,r,x)   
    
      integer i,lmax,fct
      real*8 bet,gam(lmax),a(lmax),b(lmax),c(lmax),
     &  u(lmax),r(lmax),x(lmax)

      bet=b(1)
      u(1)=r(1)/bet

      do i=2,lmax-fct
        gam(i)=c(i-1)/bet
        bet=b(i)-a(i)*gam(i)
        if (bet.EQ.0.) pause 'tridag failed'
     	u(i)=(r(i)-a(i)*u(i-1))/bet
      end do

      do i=lmax-1-fct,1,-1
        u(i)=u(i)-gam(i+1)*u(i+1)
      end do

      do i=1,lmax-fct
        x(i)=u(i)
      end do

      RETURN
      END SUBROUTINE trid

ccccc **********************************************************************
ccccc trid_cyclic: Solve the cyclic linear system.
ccccc **********************************************************************
      SUBROUTINE trid_cyclic(a,b,c,alpha,beta,lmax,r,x)

      integer l,lmax,fct
      real*8 alpha,beta,a(lmax),b(lmax),c(lmax),
     &  r(lmax),x(lmax),fact,gamma,bb(lmax),u(lmax),z(lmax)

      fct=1 !Factor (0 - Non Periodic; 1 - Periodic) 
      gamma=-b(1)
      bb(1)=b(1)-gamma 
      bb(lmax-fct)=b(lmax-fct)-alpha*beta/gamma
      do l=2,lmax-fct-1
        bb(l)=b(l)
      end do 
      call trid(fct,a,bb,c,lmax,r,x) 
      u(1)=gamma
      u(lmax-fct)=alpha
      do l=2,lmax-fct-1
        u(l)=0.d0
      end do
      call trid(fct,a,bb,c,lmax,u,z)
      fact=(x(1)+beta*x(lmax-fct)/gamma)/
     &  (1.d0+z(1)+beta*z(lmax-fct)/gamma)
      do l=1,lmax-fct
        x(l)=x(l)-fact*z(l)
      end do

      RETURN
      END SUBROUTINE trid_cyclic

      END MODULE snmethod
