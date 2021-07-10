      MODULE snmethod

      save
        real*8,  dimension(:,:), allocatable :: stx_ex_6th_D
        real*8,  dimension(:,:), allocatable :: stxx_ex_6th_D

        real*8,  dimension(:,:), allocatable :: global_stx_ex_6th_D

        real*8,  dimension(:),   allocatable :: dcy_dpy
        real*8,  dimension(:),   allocatable :: d2cy_dpy2

      contains
ccccc **********************************************************************
ccccc der_x_ex_6th_D: This routine calculated the first derivatives in 
ccccc                 the x-direction using the 6th order centered finite 
ccccc                 difference. Dirichlet boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_x_ex_6th_D(imax,jmax,kmax,stx,fc,dfc_dx)

      USE snmpi, only: my_rank,pro

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 fc(imax,jmax,kmax),dfc_dx(imax,jmax,kmax),stx(imax,7)

      do k=1,kmax
        do j=1,jmax

          if (my_rank.EQ.0) then 
            do i=1,3 
              dfc_dx(i,j,k) = 
     &          stx(i,1)*fc(1,j,k)+
     &          stx(i,2)*fc(2,j,k)+
     &          stx(i,3)*fc(3,j,k)+
     &          stx(i,4)*fc(4,j,k)+
     &          stx(i,5)*fc(5,j,k)+
     &          stx(i,6)*fc(6,j,k)+
     &          stx(i,7)*fc(7,j,k)     
            end do 
          end if  

          do i=4,imax-3 
            dfc_dx(i,j,k) = 
     &        stx(i,1)*fc(i-3,j,k)+
     &        stx(i,2)*fc(i-2,j,k)+
     &        stx(i,3)*fc(i-1,j,k)+
     &        stx(i,4)*fc(i  ,j,k)+
     &        stx(i,5)*fc(i+1,j,k)+
     &        stx(i,6)*fc(i+2,j,k)+
     &        stx(i,7)*fc(i+3,j,k)
          end do

          if (my_rank+1.EQ.pro) then 
            do i=imax-2,imax 
              dfc_dx(i,j,k) = 
     &          stx(i,1)*fc(imax-6,j,k)+
     &          stx(i,2)*fc(imax-5,j,k)+
     &          stx(i,3)*fc(imax-4,j,k)+
     &          stx(i,4)*fc(imax-3,j,k)+
     &          stx(i,5)*fc(imax-2,j,k)+
     &          stx(i,6)*fc(imax-1,j,k)+
     &          stx(i,7)*fc(imax,j,k)
            end do  
          end if  
        end do
      end do

      RETURN
      END SUBROUTINE der_x_ex_6th_D

ccccc **********************************************************************
ccccc der_xx_ex_6th_D: This routine calculated the second derivatives 
ccccc                  in the x-direction using the 6th order finite 
ccccc                  difference. Dirichlet boundary condition was 
ccccc                  adopted.
ccccc **********************************************************************
      SUBROUTINE der_xx_ex_6th_D(imax,jmax,kmax,stxx,fc,d2fc_dx2)

      USE snmpi, only: my_rank,pro

      IMPLICIT NONE
 
      integer i,j,k,imax,jmax,kmax
      real*8 fc(imax,jmax,kmax),d2fc_dx2(imax,jmax,kmax),
     &  stxx(imax,7)

      do k=1,kmax
        do j=1,jmax

          if (my_rank.EQ.0) then 
            do i=1,3 
              d2fc_dx2(i,j,k) = 
     &          stxx(i,1)*fc(1,j,k)+
     &          stxx(i,2)*fc(2,j,k)+
     &          stxx(i,3)*fc(3,j,k)+
     &          stxx(i,4)*fc(4,j,k)+
     &          stxx(i,5)*fc(5,j,k)+
     &          stxx(i,6)*fc(6,j,k)+
     &          stxx(i,7)*fc(7,j,k)
            end do
          end if

          do i=4,imax-3 
            d2fc_dx2(i,j,k) = 
     &        stxx(i,1)*fc(i-3,j,k)+
     &        stxx(i,2)*fc(i-2,j,k)+
     &        stxx(i,3)*fc(i-1,j,k)+
     &        stxx(i,4)*fc(i  ,j,k)+
     &        stxx(i,5)*fc(i+1,j,k)+
     &        stxx(i,6)*fc(i+2,j,k)+
     &        stxx(i,7)*fc(i+3,j,k)
          end do

          if (my_rank+1.EQ.pro) then 
            do i=imax-2,imax 
              d2fc_dx2(i,j,k) = 
     &          stxx(i,1)*fc(imax-6,j,k)+
     &          stxx(i,2)*fc(imax-5,j,k)+
     &          stxx(i,3)*fc(imax-4,j,k)+
     &          stxx(i,4)*fc(imax-3,j,k)+
     &          stxx(i,5)*fc(imax-2,j,k)+
     &          stxx(i,6)*fc(imax-1,j,k)+
     &          stxx(i,7)*fc(imax-0,j,k)
            end do  
          end if
        end do
      end do

      RETURN
      END SUBROUTINE der_xx_ex_6th_D

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
ccccc der_x_cp_6th_D: This routine calculated the first derivatives 
ccccc                 in the x-direction using the 6th and 5th order 
ccccc                 compact finite difference. Dirichlet boundary
ccccc                 condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_x_cp_6th_D(dx,imax,jmax,kmax,fc,dfc_dx)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dx,fc(imax,jmax,kmax),dfc_dx(imax,jmax,kmax),
     &  dfc_dx_aux(imax),aa(imax),bb(imax),cc(imax),rhs(imax) 

      do k=1,kmax        
        do j=1,jmax
          call rhs_x_6th_D(dx,imax,jmax,kmax,j,k,fc,rhs) 
          call lhs_d1_6th_D(aa,bb,cc,imax)	
          call trid(0,aa,bb,cc,imax,rhs,dfc_dx_aux)   

          do i=1,imax
            dfc_dx(i,j,k) = dfc_dx_aux(i)
          end do	  
        end do
      end do

      RETURN
      END SUBROUTINE der_x_cp_6th_D

ccccc **********************************************************************
ccccc der_xx_cp_6th_P: This routine calculated the first derivatives 
ccccc                  in the x-direction using the 6th order centered 
ccccc                  compact finite difference. Periodic boundary 
ccccc                  condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_xx_cp_6th_P(dx,imax,jmax,kmax,fc,d2fc_dx2)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dx,fc(imax,jmax,kmax),d2fc_dx2(imax,jmax,kmax),
     &  d2fc_dx2_aux(imax),aa(imax),bb(imax),cc(imax),rhs(imax) 

      do k=1,kmax
        do j=1,jmax
          call rhs_xx_6th_P(dx,imax,jmax,kmax,j,k,fc,rhs)
          call lhs_d2_6th_P(aa,bb,cc,imax)
          call trid_cyclic(aa,bb,cc,2.d0,2.d0,imax,rhs,d2fc_dx2_aux)

          do i=1,imax-1
            d2fc_dx2(i,j,k) = d2fc_dx2_aux(i)
	  end do
          d2fc_dx2(imax,j,k) = d2fc_dx2(1,j,k)          
        end do
      end do

      RETURN
      END SUBROUTINE der_xx_cp_6th_P

ccccc **********************************************************************
ccccc der_xx_cp_6th_D: This routine calculated the second derivatives 
ccccc                  in the x-direction using the 6th and 5th order 
ccccc                  compact finite difference. Dirichlet boundary
ccccc                  condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_xx_cp_6th_D(dx,imax,jmax,kmax,fc,d2fc_dx2)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dx,fc(imax,jmax,kmax),d2fc_dx2(imax,jmax,kmax),
     &  d2fc_dx2_aux(imax),aa(imax),bb(imax),cc(imax),rhs(imax)

      do k=1,kmax        
        do j=1,jmax
          call rhs_xx_6th_D(dx,imax,jmax,kmax,j,k,fc,rhs)
          call lhs_d2_6th_D(aa,bb,cc,imax)	
          call trid(0,aa,bb,cc,imax,rhs,d2fc_dx2_aux)   

          do i=1,imax
            d2fc_dx2(i,j,k) = d2fc_dx2_aux(i)
          end do	  
         end do
      end do

      RETURN
      END SUBROUTINE der_xx_cp_6th_D

ccccc **********************************************************************
ccccc der_y_cp_6th_P: This routine calculated the first derivatives 
ccccc                 in the y-direction using the 6th order centered 
ccccc                 compact finite difference. Periodic boundary 
ccccc                 condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_6th_P(dy,imax,jmax,kmax,fc,dfc_dy)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,aa(jmax),bb(jmax),cc(jmax),rhs(jmax),
     &  fc(imax,jmax,kmax),dfc_dy_aux(jmax),dfc_dy(imax,jmax,kmax)

      if (jmax.NE.1) then
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
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax  
              dfc_dy(i,j,k) = 0.d0
            end do
          end do
        end do
      end if

      RETURN
      END SUBROUTINE der_y_cp_6th_P

ccccc **********************************************************************
ccccc der_y_cp_62th_D: This routine calculated the first derivatives 
ccccc                  in the y-direction using the 6th and 2nd order 
ccccc                  compact finite difference. Dirichlet boundary
ccccc                  condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_62th_D(dy,imax,jmax,kmax,mtrc,fc,dfc_dy)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc(jmax),fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),
     &  dfc_dy_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      if (jmax.NE.1) then
        do k=1,kmax        
          do i=1,imax  
            call rhs_y_62th_D(dy,imax,jmax,kmax,i,k,fc,rhs)
            call lhs_d1_62th_D(aa,bb,cc,mtrc,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,dfc_dy_aux)   

            do j=1,jmax
              dfc_dy(i,j,k) = dfc_dy_aux(j)
            end do	  
          end do
        end do
      else 
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              dfc_dy(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

      RETURN
      END SUBROUTINE der_y_cp_62th_D

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

      if (jmax.NE.1) then
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
      else 
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              dfc_dy(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

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

      if (jmax.NE.1) then
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
      else 
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              dfc_dy(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

      RETURN
      END SUBROUTINE der_y_cp_65th_WN						

ccccc **********************************************************************
ccccc der_y_cp_62th_N: This routine calculated the first derivatives 
ccccc                  in the y-direction using the 6th and 2nd order 
ccccc                  compact finite difference. The Neumann boundary
ccccc                  condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_y_cp_62th_N(dy,imax,jmax,kmax,mtrc,fc,dfc_dy)      

      IMPLICIT NONE
 
      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc(jmax),fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),
     &  dfc_dy_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      if (jmax.NE.1) then
        do k=1,kmax
          do i=1,imax
            call rhs_y_62th_N(dy,imax,jmax,kmax,i,k,fc,rhs)
            call lhs_d1_6th_N(aa,bb,cc,mtrc,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,dfc_dy_aux)   

            do j=1,jmax
              dfc_dy(i,j,k) = dfc_dy_aux(j)
            end do	  
          end do
        end do
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              dfc_dy(i,j,k) = 0.d0
            end do
          end do
        end do  
      end if

      RETURN
      END SUBROUTINE der_y_cp_62th_N

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

      if (jmax.NE.1) then
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
      else 
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              dfc_dy(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

      RETURN
      END SUBROUTINE der_y_cp_65th_NN									

ccccc **********************************************************************
ccccc der_yy_cp_6th_P: This routine calculated the second derivatives 
ccccc                  in the y-direction using the 6th order centered 
ccccc                  compact finite difference. Periodic boundary 
ccccc                  condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_yy_cp_6th_P(dy,imax,jmax,kmax,fc,d2fc_dy2)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,fc(imax,jmax,kmax),d2fc_dy2(imax,jmax,kmax),
     &  d2fc_dy2_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      if (jmax.NE.1) then
        do k=1,kmax
          do i=1,imax
            call rhs_yy_6th_P(dy,imax,jmax,kmax,i,k,fc,rhs)
            call lhs_d2_6th_P(aa,bb,cc,jmax)
            call trid_cyclic(aa,bb,cc,2.d0,2.d0,jmax,rhs,d2fc_dy2_aux)

            do j=1,jmax-1
              d2fc_dy2(i,j,k) = d2fc_dy2_aux(j)
            end do
            d2fc_dy2(i,jmax,k) = d2fc_dy2(i,1,k)
          end do
        end do 
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax  
              d2fc_dy2(i,j,k) = 0.d0
            end do
          end do
        end do
      end if

      RETURN
      END SUBROUTINE der_yy_cp_6th_P

ccccc **********************************************************************
ccccc der_yy_cp_62th_D: This routine calculated the second derivatives 
ccccc                   in the y-direction using the 6th and 2nd order 
ccccc                   compact finite difference. Dirichlet boundary
ccccc                   condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_yy_cp_62th_D(dy,imax,jmax,kmax,
     &  mtrc_d1,mtrc_d2,fc,dfc_dy,d2fc_dy2)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),fc(imax,jmax,kmax),
     &  dfc_dy(imax,jmax,kmax),d2fc_dy2(imax,jmax,kmax),
     &  d2fc_dy2_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      if (jmax.NE.1) then
        do k=1,kmax        
          do i=1,imax
            call rhs_yy_62th_D(dy,imax,jmax,kmax,i,k,
     &        mtrc_d1,mtrc_d2,fc,dfc_dy,rhs)
            call lhs_d2_62th_D(aa,bb,cc,mtrc_d1,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,d2fc_dy2_aux)   

            do j=1,jmax
              d2fc_dy2(i,j,k)=d2fc_dy2_aux(j)
            end do	  
          end do
        end do
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              d2fc_dy2(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

      RETURN
      END SUBROUTINE der_yy_cp_62th_D

ccccc **********************************************************************
ccccc der_yy_cp_65th_WD: This routine calculated the second derivatives 
ccccc                    in the y-direction using the 6th and 5th order 
ccccc                    compact finite difference. Dirichlet boundary 
ccccc                    condition was used to simulated the wall.
ccccc **********************************************************************
      SUBROUTINE der_yy_cp_65th_WD(dy,imax,jmax,kmax,mtrc_d1,mtrc_d2,fc,dfc_dy,d2fc_dy2)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),fc(imax,jmax,kmax),
     &  dfc_dy(imax,jmax,kmax),d2fc_dy2(imax,jmax,kmax),
     &  d2fc_dy2_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      if (jmax.NE.1) then
        do k=1,kmax        
          do i=1,imax
            call rhs_yy_65th_WD(dy,imax,jmax,kmax,i,k,
     &        mtrc_d1,mtrc_d2,fc,dfc_dy,rhs)
            call lhs_d2_65th_W(aa,bb,cc,mtrc_d1,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,d2fc_dy2_aux)   

            do j=1,jmax
              d2fc_dy2(i,j,k)=d2fc_dy2_aux(j)
            end do	  
          end do
        end do
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              d2fc_dy2(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

      RETURN
      END SUBROUTINE der_yy_cp_65th_WD						

ccccc **********************************************************************
ccccc der_yy_cp_65th_WN: This routine calculated the second derivatives 
ccccc                    in the y-direction using the 6th and 5th order 
ccccc                    compact finite difference. Dirichlet and Neumann
ccccc                    boundary condition was used to simulated the wall.
ccccc **********************************************************************
      SUBROUTINE der_yy_cp_65th_WN(dy,imax,jmax,kmax,mtrc_d1,mtrc_d2,fc,dfc_dy,d2fc_dy2)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),fc(imax,jmax,kmax),
     &  dfc_dy(imax,jmax,kmax),d2fc_dy2(imax,jmax,kmax),
     &  d2fc_dy2_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      if (jmax.NE.1) then
        do k=1,kmax        
          do i=1,imax
            call rhs_yy_65th_WN(dy,imax,jmax,kmax,i,k,
     &        mtrc_d1,mtrc_d2,fc,dfc_dy,rhs)
            call lhs_d2_65th_W(aa,bb,cc,mtrc_d1,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,d2fc_dy2_aux)   

            do j=1,jmax
              d2fc_dy2(i,j,k)=d2fc_dy2_aux(j)
            end do	  
          end do
        end do
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              d2fc_dy2(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

      RETURN
      END SUBROUTINE der_yy_cp_65th_WN						

ccccc **********************************************************************
ccccc der_yy_cp_62th_N: This routine calculated the second derivatives 
ccccc                   in the y-direction using the 6th and 2nd order 
ccccc                   compact finite difference. Neumann boundary
ccccc                   condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_yy_cp_62th_N(dy,imax,jmax,kmax,
     &  mtrc_d1,mtrc_d2,fc,dfc_dy,d2fc_dy2)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),fc(imax,jmax,kmax),
     &  dfc_dy(imax,jmax,kmax),d2fc_dy2(imax,jmax,kmax),
     &  d2fc_dy2_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      if (jmax.NE.1) then
        do k=1,kmax
          do i=1,imax
            call rhs_yy_62th_N(dy,imax,jmax,kmax,i,k,
     &        mtrc_d1,mtrc_d2,fc,dfc_dy,rhs)
            call lhs_d2_6th_N(aa,bb,cc,mtrc_d1,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,d2fc_dy2_aux)   

            do j=1,jmax
              d2fc_dy2(i,j,k)=d2fc_dy2_aux(j)
            end do	  
          end do
        end do
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              d2fc_dy2(i,j,k) = 0.d0
            end do
          end do
        end do  
      end if

      RETURN
      END SUBROUTINE der_yy_cp_62th_N

ccccc **********************************************************************
ccccc der_yy_cp_65th_NN: This routine calculated the second derivatives 
ccccc                    in the y-direction using the 6th and 5th order 
ccccc                    compact finite difference. Neumann boundary
ccccc                    condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_yy_cp_65th_NN(dy,imax,jmax,kmax,mtrc_d1,mtrc_d2,fc,dfc_dy,d2fc_dy2)      

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),fc(imax,jmax,kmax),
     &  dfc_dy(imax,jmax,kmax),d2fc_dy2(imax,jmax,kmax),
     &  d2fc_dy2_aux(jmax),aa(jmax),bb(jmax),cc(jmax),rhs(jmax)

      if (jmax.NE.1) then
        do k=1,kmax
          do i=1,imax
            call rhs_yy_65th_NN(dy,imax,jmax,kmax,i,k,
     &        mtrc_d1,mtrc_d2,fc,dfc_dy,rhs)
            call lhs_d2_65th_NN(aa,bb,cc,mtrc_d1,jmax)	
            call trid(0,aa,bb,cc,jmax,rhs,d2fc_dy2_aux)   

            do j=1,jmax
              d2fc_dy2(i,j,k)=d2fc_dy2_aux(j)
            end do	  
          end do
        end do
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              d2fc_dy2(i,j,k) = 0.d0
            end do
          end do
        end do  
      end if

      RETURN
      END SUBROUTINE der_yy_cp_65th_NN

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

      if (kmax.NE.1) then
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
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              dfc_dz(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

      RETURN
      END SUBROUTINE der_z_cp_6th_P

ccccc **********************************************************************
ccccc der_zz_cp_6th_P: This routine calculated the second derivatives 
ccccc                  in the z-direction using the 6th order centered 
ccccc                  compact finite difference. Periodic boundary 
ccccc                  condition was adopted.
ccccc **********************************************************************
      SUBROUTINE der_zz_cp_6th_P(dz,imax,jmax,kmax,fc,d2fc_dz2)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dz,fc(imax,jmax,kmax),d2fc_dz2(imax,jmax,kmax),
     &  d2fc_dz2_aux(kmax),aa(kmax),bb(kmax),cc(kmax),rhs(kmax)

      if (kmax.NE.1) then
        do j=1,jmax
          do i=1,imax
            call rhs_zz_6th_P(dz,imax,jmax,kmax,i,j,fc,rhs)
            call lhs_d2_6th_P(aa,bb,cc,kmax)
            call trid_cyclic(aa,bb,cc,2.d0,2.d0,kmax,rhs,d2fc_dz2_aux)
 
            do k=1,kmax-1
              d2fc_dz2(i,j,k) = d2fc_dz2_aux(k)
            end do
            d2fc_dz2(i,j,kmax) = d2fc_dz2(i,j,1)
          end do
        end do
      else
        do k=1,kmax
          do j=1,jmax
            do i=1,imax      
              d2fc_dz2(i,j,k) = 0.d0
            end do
          end do  
        end do
      end if

      RETURN
      END SUBROUTINE der_zz_cp_6th_P

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
ccccc rhs_x_6th_D: This routine storing the right hand side of linear 
ccccc              system of first derivatives in the x-direction. The 
ccccc              6th order compact finite difference schemes was used.
ccccc              Dirichlet boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_x_6th_D(dx,imax,jmax,kmax,j,k,fc,r) 

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dx,fc(imax,jmax,kmax),r(imax)

      r(1)=(-1764.d0*fc(1,j,k)+4320.d0*fc(2,j,k)-5400.d0*fc(3,j,k)        
     &  +4800.d0*fc(4,j,k)-2700.d0*fc(5,j,k)+864.d0*fc(6,j,k)
     &  -120.d0*fc(7,j,k))/(720.d0*dx)                                   
      r(2)=(-120.d0*fc(1,j,k)-924.d0*fc(2,j,k)+1800.d0*fc(3,j,k)       
     &  -1200.d0*fc(4,j,k)+600.d0*fc(5,j,k)-180.d0*fc(6,j,k)
     &  +24.d0*fc(7,j,k))/(720.d0*dx)                                    
      do i=3,imax-2 
        r(i)=(-fc(i-2,j,k)-28.d0*fc(i-1,j,k)                       
     &    +28.d0*fc(i+1,j,k)+fc(i+2,j,k))/(12.d0*dx)
      end do	
      r(imax-1)=(-24.d0*fc(imax-6,j,k)+180.d0*fc(imax-5,j,k)           
     &  -600.d0*fc(imax-4,j,k)+1200.d0*fc(imax-3,j,k)
     &  -1800.d0*fc(imax-2,j,k)+924.d0*fc(imax-1,j,k)
     &  +120.d0*fc(imax,j,k))/(720.d0*dx)                                
      r(imax)=(120.d0*fc(imax-6,j,k)-864.d0*fc(imax-5,j,k)           
     &  +2700.d0*fc(imax-4,j,k)-4800.d0*fc(imax-3,j,k)
     &  +5400.d0*fc(imax-2,j,k)-4320.d0*fc(imax-1,j,k)
     &  +1764.d0*fc(imax,j,k))/(720.d0*dx)                              

      RETURN
      END SUBROUTINE rhs_x_6th_D

ccccc **********************************************************************
ccccc rhs_xx_6th_P: This routine storing the right hand side of linear  
ccccc               system of second derivatives in the x-direction. The 
ccccc               6th order compact finite difference schemes was used. 
ccccc               Periodic boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_xx_6th_P(dx,imax,jmax,kmax,j,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,imax_1,jmax,kmax
      real*8 dx,fc(imax,jmax,kmax),r(imax)

      imax_1=imax-1

      r(1)=(3.d0*fc(imax_1-1,j,k)+48.d0*fc(imax_1,j,k)-102.d0*fc(1,j,k)
     &  +48.d0*fc(2,j,k)+3.d0*fc(3,j,k))/(4.d0*dx*dx)
      r(2)=(3.d0*fc(imax_1,j,k)+48.d0*fc(1,j,k)-102.d0*fc(2,j,k)
     &  +48.d0*fc(3,j,k)+3.d0*fc(4,j,k))/(4.d0*dx*dx)
      do i=3,imax_1-2
        r(i)=(3.d0*fc(i-2,j,k)+48.d0*fc(i-1,j,k)-102.d0*fc(i,j,k)
     &    +48.d0*fc(i+1,j,k)+3.d0*fc(i+2,j,k))/(4.d0*dx*dx)
      end do
      r(imax_1-1)=(3.d0*fc(imax_1-3,j,k)+48.d0*fc(imax_1-2,j,k)
     &  -102.d0*fc(imax_1-1,j,k)+48.d0*fc(imax_1,j,k)
     &  +3.d0*fc(1,j,k))/(4.d0*dx*dx)
      r(imax_1)=(3.d0*fc(imax_1-2,j,k)+48.d0*fc(imax_1-1,j,k)
     &  -102.d0*fc(imax_1,j,k)+48.d0*fc(1,j,k)
     &  +3.d0*fc(2,j,k))/(4.d0*dx*dx)

      !The routine trid and trid_cyclic was modified to despise 
      !the last point of computational domain for periodic 
      !boundary conditions. Therefore the last point of vector  
      !r is equal zero as follow below.
      r(imax)=0.d0

      RETURN
      END SUBROUTINE rhs_xx_6th_P

ccccc **********************************************************************
ccccc rhs_xx_6th_D: This routine storing the right hand side of linear  
ccccc               system of second derivatives in the x-direction. The  
ccccc               6th order compact finite difference schemes  was used.    
ccccc               Dirichlet boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_xx_6th_D(dx,imax,jmax,kmax,j,k,fc,r) 

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax 
      real*8 dx,fc(imax,jmax,kmax),r(imax)
 
      r(1)=(13132.d0*fc(1,j,k)-56196.d0*fc(2,j,k)+110754.d0*fc(3,j,k)                        
     &  -132860.d0*fc(4,j,k)+103320.d0*fc(5,j,k)-50652.d0*fc(6,j,k)
     &  +14266.d0*fc(7,j,k)-1764.d0*fc(8,j,k))/(2520.d0*dx*dx)         
      r(2)=(1764.d0*fc(1,j,k)-980.d0*fc(2,j,k)-6804.d0*fc(3,j,k) 
     &  +11970.d0*fc(4,j,k)-9380.d0*fc(5,j,k)+4536.d0*fc(6,j,k)
     &  -1260.d0*fc(7,j,k)+154.d0*fc(8,j,k))/(2520.d0*dx*dx)             
      do i=3,imax-2 
        r(i)=(3.d0*fc(i-2,j,k)+48.d0*fc(i-1,j,k)-102.d0*fc(i,j,k)  
     &    +48.d0*fc(i+1,j,k)+3.d0*fc(i+2,j,k))/(4.d0*dx*dx)
      end do	
      r(imax-1)=(154.d0*fc(imax-7,j,k)-1260.d0*fc(imax-6,j,k)
     &  +4536.d0*fc(imax-5,j,k)-9380.d0*fc(imax-4,j,k)
     &  +11970.d0*fc(imax-3,j,k)-6804.d0*fc(imax-2,j,k) 
     &  -980.d0*fc(imax-1,j,k)+1764.d0*fc(imax,j,k))
     &  /(2520.d0*dx*dx)                                                 
      r(imax)=(-1764.d0*fc(imax-7,j,k)+14266.d0*fc(imax-6,j,k)
     &  -50652.d0*fc(imax-5,j,k)+103320.d0*fc(imax-4,j,k)
     &  -132860.d0*fc(imax-3,j,k)+110754.d0*fc(imax-2,j,k)
     &  -56196.d0*fc(imax-1,j,k)+13132.d0*fc(imax,j,k))/(2520.d0*dx*dx)  

      RETURN
      END SUBROUTINE rhs_xx_6th_D

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
	            
      r(1)=(-3.d0*fc(i,1,k)+4.d0*fc(i,2,k)-1.d0*fc(i,3,k))/(2.d0*dy)
      r(2)=(-fc(i,1,k)+fc(i,3,k))/(2.d0*dy)      
      do j=3,jmax-2 
        r(j)=(-fc(i,j-2,k)-28.d0*fc(i,j-1,k)
     &    +28.d0*fc(i,j+1,k)+fc(i,j+2,k))/(12.d0*dy)
      end do	
      r(jmax-1)=(-fc(i,jmax-2,k)+fc(i,jmax,k))/(2.d0*dy)      
      r(jmax)=(1.d0*fc(i,jmax-2,k)-4.d0*fc(i,jmax-1,k)
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
ccccc rhs_yy_6th_P: This routine storing the right hand side of linear  
ccccc               system of second derivatives in the y-direction. The 
ccccc               6th order compact finite difference schemes was used. 
ccccc               Periodic boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_yy_6th_P(dy,imax,jmax,kmax,i,k,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,jmax_1,kmax
      real*8 dy,fc(imax,jmax,kmax),r(imax)

      jmax_1=jmax-1

      r(1)=(3.d0*fc(i,jmax_1-1,k)+48.d0*fc(i,jmax_1,k)-102.d0*fc(i,1,k)
     &  +48.d0*fc(i,2,k)+3.d0*fc(i,3,k))/(4.d0*dy*dy)
      r(2)=(3.d0*fc(i,jmax_1,k)+48.d0*fc(i,1,k)-102.d0*fc(i,2,k)
     &  +48.d0*fc(i,3,k)+3.d0*fc(i,4,k))/(4.d0*dy*dy)
      do j=3,jmax_1-2
        r(j)=(3.d0*fc(i,j-2,k)+48.d0*fc(i,j-1,k)-102.d0*fc(i,j,k)
     &    +48.d0*fc(i,j+1,k)+3.d0*fc(i,j+2,k))/(4.d0*dy*dy)
      end do
      r(jmax_1-1)=(3.d0*fc(i,jmax_1-3,k)+48.d0*fc(i,jmax_1-2,k)
     &  -102.d0*fc(i,jmax_1-1,k)+48.d0*fc(i,jmax_1,k)
     &  +3.d0*fc(i,1,k))/(4.d0*dy*dy)
      r(jmax_1)=(3.d0*fc(i,jmax_1-2,k)+48.d0*fc(i,jmax_1-1,k)
     &  -102.d0*fc(i,jmax_1,k)+48.d0*fc(i,1,k)
     &  +3.d0*fc(i,2,k))/(4.d0*dy*dy)

      !The routine trid and trid_cyclic was modified to despise 
      !the last point of computational domain for periodic 
      !boundary conditions. Therefore the last point of vector  
      !r is equal zero as follow below.
      r(jmax)=0.d0

      RETURN
      END SUBROUTINE rhs_yy_6th_P

ccccc **********************************************************************
ccccc rhs_yy_62th_D: This routine storing the right hand side of linear  
ccccc                system of second derivatives in the y-direction. The  
ccccc                6th and 2nd order compact finite difference schemes 
ccccc                was used. Dirichlet boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_yy_62th_D(dy,imax,jmax,kmax,i,k,
     &  mtrc_d1,mtrc_d2,fc,dfc_dy,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),
     &  fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),r(jmax)

      r(1)=(fc(i,1,k)-2.d0*fc(i,2,k)+fc(i,3,k))/(dy*dy)
     &  +mtrc_d2(1)/(mtrc_d1(1)**3.d0)*dfc_dy(i,1,k)
      r(2)=(fc(i,1,k)-2.d0*fc(i,2,k)+fc(i,3,k))/(dy*dy)
     &  +mtrc_d2(2)/(mtrc_d1(2)**3.d0)*dfc_dy(i,2,k)
      do j=3,jmax-2 
        r(j)=(3.d0*fc(i,j-2,k)+48.d0*fc(i,j-1,k)-102.d0*fc(i,j,k)
     &    +48.d0*fc(i,j+1,k)+3.d0*fc(i,j+2,k))/(4.d0*dy*dy)
     &    +02.d0*mtrc_d2(j-1)/(mtrc_d1(j-1)**3.d0)*dfc_dy(i,j-1,k) 
     &    +11.d0*mtrc_d2(j  )/(mtrc_d1(j  )**3.d0)*dfc_dy(i,j  ,k) 
     &    +02.d0*mtrc_d2(j+1)/(mtrc_d1(j+1)**3.d0)*dfc_dy(i,j+1,k) 
      end do	
      r(jmax-1)=(fc(i,jmax-2,k)-2.d0*fc(i,jmax-1,k)+fc(i,jmax,k))/(dy*dy)
     &  +mtrc_d2(jmax-1)/(mtrc_d1(jmax-1)**3.d0)*dfc_dy(i,jmax-1,k)
      r(jmax)=(fc(i,jmax-2,k)-2.d0*fc(i,jmax-1,k)+fc(i,jmax,k))/(dy*dy)
     &  +mtrc_d2(jmax)/(mtrc_d1(jmax)**3.d0)*dfc_dy(i,jmax,k) 

      RETURN
      END SUBROUTINE rhs_yy_62th_D

ccccc **********************************************************************
ccccc rhs_yy_65th_WD: This routine storing the right hand side of linear  
ccccc                 system of second derivatives in the y-direction. The  
ccccc                 6th and 5th order compact finite difference schemes 
ccccc                 was used. Dirichlet boundary condition was 
ccccc                 adopted to simulate the wall.
ccccc **********************************************************************
      SUBROUTINE rhs_yy_65th_WD(dy,imax,jmax,kmax,i,k,
     &  mtrc_d1,mtrc_d2,fc,dfc_dy,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),
     &  fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),r(jmax)

      r(1)=(9775.d0*fc(i,1,k)-20285.d0*fc(i,2,k)
     &  +11170.d0*fc(i,3,k)-550.d0*fc(i,4,k)-145.d0*fc(i,5,k)
     &  +35.d0*fc(i,6,k))/(60.d0*dy*dy)
     &  +013.d0*mtrc_d2(1)/(mtrc_d1(1)**3.d0)*dfc_dy(i,1,k)
     &  +137.d0*mtrc_d2(2)/(mtrc_d1(2)**3.d0)*dfc_dy(i,2,k)
      r(2)=(4834.d0*fc(i,1,k)-8424.d0*fc(i,2,k)+1890.d0*fc(i,3,k)
     &  +2320.d0*fc(i,4,k)-810.d0*fc(i,5,k)+216.d0*fc(i,6,k)
     &  -26.d0*fc(i,7,k))/(360.d0*dy*dy)
     &  +01.d0*mtrc_d2(1)/(mtrc_d1(1)**3.d0)*dfc_dy(i,1,k)
     &  +12.d0*mtrc_d2(2)/(mtrc_d1(2)**3.d0)*dfc_dy(i,2,k)
     &  +03.d0*mtrc_d2(3)/(mtrc_d1(3)**3.d0)*dfc_dy(i,3,k)
      do j=3,jmax-2 
        r(j)=(3.d0*fc(i,j-2,k)+48.d0*fc(i,j-1,k)-102.d0*fc(i,j,k)
     &    +48.d0*fc(i,j+1,k)+3.d0*fc(i,j+2,k))/(4.d0*dy*dy)
     &    +02.d0*mtrc_d2(j-1)/(mtrc_d1(j-1)**3.d0)*dfc_dy(i,j-1,k) 
     &    +11.d0*mtrc_d2(j  )/(mtrc_d1(j  )**3.d0)*dfc_dy(i,j  ,k) 
     &    +02.d0*mtrc_d2(j+1)/(mtrc_d1(j+1)**3.d0)*dfc_dy(i,j+1,k) 
      end do	
      r(jmax-1)=(fc(i,jmax-2,k)-2.d0*fc(i,jmax-1,k)+fc(i,jmax,k))/(dy*dy)
     &  +mtrc_d2(jmax-1)/(mtrc_d1(jmax-1)**3.d0)*dfc_dy(i,jmax-1,k)
      r(jmax)=(fc(i,jmax-2,k)-2.d0*fc(i,jmax-1,k)+fc(i,jmax,k))/(dy*dy)
     &  +mtrc_d2(jmax)/(mtrc_d1(jmax)**3.d0)*dfc_dy(i,jmax,k)  

      RETURN
      END SUBROUTINE rhs_yy_65th_WD

ccccc **********************************************************************
ccccc rhs_yy_65th_WN: This routine storing the right hand side of linear  
ccccc                 system of second derivatives in the y-direction. The  
ccccc                 6th and 5th order compact finite difference schemes 
ccccc                 was used. Dirichlet and Neumann boundary condition  
ccccc                 was adopted to simulate the wall.
ccccc **********************************************************************
      SUBROUTINE rhs_yy_65th_WN(dy,imax,jmax,kmax,i,k,
     &  mtrc_d1,mtrc_d2,fc,dfc_dy,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),
     &  fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),r(jmax)

      r(1)=(9775.d0*fc(i,1,k)-20285.d0*fc(i,2,k)
     &  +11170.d0*fc(i,3,k)-550.d0*fc(i,4,k)-145.d0*fc(i,5,k)
     &  +35.d0*fc(i,6,k))/(60.d0*dy*dy)
     &  +013.d0*mtrc_d2(1)/(mtrc_d1(1)**3.d0)*dfc_dy(i,1,k)
     &  +137.d0*mtrc_d2(2)/(mtrc_d1(2)**3.d0)*dfc_dy(i,2,k)
      r(2)=(4834.d0*fc(i,1,k)-8424.d0*fc(i,2,k)+1890.d0*fc(i,3,k)
     &  +2320.d0*fc(i,4,k)-810.d0*fc(i,5,k)+216.d0*fc(i,6,k)
     &  -26.d0*fc(i,7,k))/(360.d0*dy*dy)
     &  +01.d0*mtrc_d2(1)/(mtrc_d1(1)**3.d0)*dfc_dy(i,1,k)
     &  +12.d0*mtrc_d2(2)/(mtrc_d1(2)**3.d0)*dfc_dy(i,2,k)
     &  +03.d0*mtrc_d2(3)/(mtrc_d1(3)**3.d0)*dfc_dy(i,3,k)
      do j=3,jmax-2 
        r(j)=(3.d0*fc(i,j-2,k)+48.d0*fc(i,j-1,k)-102.d0*fc(i,j,k)
     &    +48.d0*fc(i,j+1,k)+3.d0*fc(i,j+2,k))/(4.d0*dy*dy)
     &    +02.d0*mtrc_d2(j-1)/(mtrc_d1(j-1)**3.d0)*dfc_dy(i,j-1,k) 
     &    +11.d0*mtrc_d2(j  )/(mtrc_d1(j  )**3.d0)*dfc_dy(i,j  ,k) 
     &    +02.d0*mtrc_d2(j+1)/(mtrc_d1(j+1)**3.d0)*dfc_dy(i,j+1,k) 
      end do	
      r(jmax-1)=(fc(i,jmax-2,k)-2.d0*fc(i,jmax-1,k)+fc(i,jmax,k))/(dy*dy)
     &  +mtrc_d2(jmax-1)/(mtrc_d1(jmax-1)**3.d0)*dfc_dy(i,jmax-1,k)
      r(jmax)=(-7.d0*fc(i,jmax,k)+8.d0*fc(i,jmax-1,k)-fc(i,jmax-2,k))
     &  /(2.d0*dy*dy)
     &  +mtrc_d2(jmax)/(mtrc_d1(jmax)**3.d0)*dfc_dy(i,jmax,k)
      
      RETURN
      END SUBROUTINE rhs_yy_65th_WN

ccccc **********************************************************************
ccccc rhs_yy_62th_N: This routine storing the right hand side of linear  
ccccc                system of second derivatives in the y-direction. The  
ccccc                6th and 2nd order compact finite difference schemes   
ccccc                was used. Neumann boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_yy_62th_N(dy,imax,jmax,kmax,i,k,
     &  mtrc_d1,mtrc_d2,fc,dfc_dy,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),
     &  fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),r(jmax)

      r(1)=(-7.d0*fc(i,1,k)+8.d0*fc(i,2,k)-fc(i,3,k))/(2.d0*dy*dy)
     &  +mtrc_d2(1)/(mtrc_d1(1)**3.d0)*dfc_dy(i,1,k)
      r(2)=(fc(i,1,k)-2.d0*fc(i,2,k)+fc(i,3,k))/(dy*dy)
     &  +mtrc_d2(2)/(mtrc_d1(2)**3.d0)*dfc_dy(i,2,k)
      do j=3,jmax-2 
        r(j)=(3.d0*fc(i,j-2,k)+48.d0*fc(i,j-1,k)-102.d0*fc(i,j,k)
     &    +48.d0*fc(i,j+1,k)+3.d0*fc(i,j+2,k))/(4.d0*dy*dy)
     &    +02.d0*mtrc_d2(j-1)/(mtrc_d1(j-1)**3.d0)*dfc_dy(i,j-1,k) 
     &    +11.d0*mtrc_d2(j  )/(mtrc_d1(j  )**3.d0)*dfc_dy(i,j  ,k) 
     &    +02.d0*mtrc_d2(j+1)/(mtrc_d1(j+1)**3.d0)*dfc_dy(i,j+1,k) 
      end do	
      r(jmax-1)=(fc(i,jmax-2,k)-2.d0*fc(i,jmax-1,k)+fc(i,jmax,k))/(dy*dy)
     &  +mtrc_d2(jmax-1)/(mtrc_d1(jmax-1)**3.d0)*dfc_dy(i,jmax-1,k)
      r(jmax)=(-7.d0*fc(i,jmax,k)+8.d0*fc(i,jmax-1,k)-fc(i,jmax-2,k))
     &  /(2.d0*dy*dy)
     &  +mtrc_d2(jmax)/(mtrc_d1(jmax)**3.d0)*dfc_dy(i,jmax,k)

      RETURN
      END SUBROUTINE rhs_yy_62th_N

ccccc **********************************************************************
ccccc rhs_yy_65th_NN: This routine storing the right hand side of linear  
ccccc                 system of second derivatives in the y-direction. The  
ccccc                 6th and 5th order compact finite difference schemes   
ccccc                 was used. Neumann boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_yy_65th_NN(dy,imax,jmax,kmax,i,k,
     &  mtrc_d1,mtrc_d2,fc,dfc_dy,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 dy,mtrc_d1(jmax),mtrc_d2(jmax),
     &  fc(imax,jmax,kmax),dfc_dy(imax,jmax,kmax),r(jmax)

      r(1)=(9775.d0*fc(i,1,k)-20285.d0*fc(i,2,k)
     &  +11170.d0*fc(i,3,k)-550.d0*fc(i,4,k)-145.d0*fc(i,5,k)
     &  +35.d0*fc(i,6,k))/(60.d0*dy*dy)
     &  -13.d0*(-50.d0*fc(i,1,k)+96.d0*fc(i,2,k)-72.d0*fc(i,3,k)
     &  +32.d0*fc(i,4,k)-6.d0*fc(i,5,k))/(60.d0*dy*dy)
     &  +013.d0*mtrc_d2(1)/(mtrc_d1(1)**3.d0)*dfc_dy(i,1,k)
     &  +137.d0*mtrc_d2(2)/(mtrc_d1(2)**3.d0)*dfc_dy(i,2,k)
      r(2)=(4834.d0*fc(i,1,k)-8424.d0*fc(i,2,k)+1890.d0*fc(i,3,k)
     &  +2320.d0*fc(i,4,k)-810.d0*fc(i,5,k)+216.d0*fc(i,6,k)
     &  -26.d0*fc(i,7,k))/(360.d0*dy*dy)
     &  -13.d0*(-50.d0*fc(i,1,k)+96.d0*fc(i,2,k)-72.d0*fc(i,3,k)
     &  +32.d0*fc(i,4,k)-6.d0*fc(i,5,k))/(60.d0*dy*dy)
     &  +01.d0*mtrc_d2(1)/(mtrc_d1(1)**3.d0)*dfc_dy(i,1,k)
     &  +12.d0*mtrc_d2(2)/(mtrc_d1(2)**3.d0)*dfc_dy(i,2,k)
     &  +03.d0*mtrc_d2(3)/(mtrc_d1(3)**3.d0)*dfc_dy(i,3,k)
      do j=3,jmax-2 
        r(j)=(3.d0*fc(i,j-2,k)+48.d0*fc(i,j-1,k)-102.d0*fc(i,j,k)
     &    +48.d0*fc(i,j+1,k)+3.d0*fc(i,j+2,k))/(4.d0*dy*dy)
     &    +02.d0*mtrc_d2(j-1)/(mtrc_d1(j-1)**3.d0)*dfc_dy(i,j-1,k) 
     &    +11.d0*mtrc_d2(j  )/(mtrc_d1(j  )**3.d0)*dfc_dy(i,j  ,k)
     &    +02.d0*mtrc_d2(j+1)/(mtrc_d1(j+1)**3.d0)*dfc_dy(i,j+1,k) 
      end do	
      r(jmax-1)=(fc(i,jmax-2,k)-2.d0*fc(i,jmax-1,k)+fc(i,jmax,k))/(dy*dy)
     &  +mtrc_d2(jmax-1)/(mtrc_d1(jmax-1)**3.d0)*dfc_dy(i,jmax-1,k)
      r(jmax)=(-7.d0*fc(i,jmax,k)+8.d0*fc(i,jmax-1,k)-fc(i,jmax-2,k))
     &  /(2.d0*dy*dy)
     &  +mtrc_d2(jmax)/(mtrc_d1(jmax)**3.d0)*dfc_dy(i,jmax,k)

      RETURN
      END SUBROUTINE rhs_yy_65th_NN

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
ccccc rhs_zz_6th_P: This routine storing the right hand side of linear  
ccccc               system of second derivatives in the z-direction. The 
ccccc               6th order compact finite difference schemes was used. 
ccccc               Periodic boundary condition was adopted.
ccccc **********************************************************************
      SUBROUTINE rhs_zz_6th_P(dz,imax,jmax,kmax,i,j,fc,r)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax,kmax_1
      real*8 dz,fc(imax,jmax,kmax),r(imax)

      kmax_1=kmax-1

      r(1)=(3.d0*fc(i,j,kmax_1-1)+48.d0*fc(i,j,kmax_1)-102.d0*fc(i,j,1)
     &  +48.d0*fc(i,j,2)+3.d0*fc(i,j,3))/(4.d0*dz*dz)
      r(2)=(3.d0*fc(i,j,kmax_1)+48.d0*fc(i,j,1)-102.d0*fc(i,j,2)
     &  +48.d0*fc(i,j,3)+3.d0*fc(i,j,4))/(4.d0*dz*dz)
      do k=3,kmax_1-2
        r(k)=(3.d0*fc(i,j,k-2)+48.d0*fc(i,j,k-1)-102.d0*fc(i,j,k)
     &    +48.d0*fc(i,j,k+1)+3.d0*fc(i,j,k+2))/(4.d0*dz*dz)
      end do
      r(kmax_1-1)=(3.d0*fc(i,j,kmax_1-3)+48.d0*fc(i,j,kmax_1-2)
     &  -102.d0*fc(i,j,kmax_1-1)+48.d0*fc(i,j,kmax_1)
     &  +3.d0*fc(i,j,1))/(4.d0*dz*dz)
      r(kmax_1)=(3.d0*fc(i,j,kmax_1-2)+48.d0*fc(i,j,kmax_1-1)
     &  -102.d0*fc(i,j,kmax_1)+48.d0*fc(i,j,1)
     &  +3.d0*fc(i,j,2))/(4.d0*dz*dz)

      !The routine trid and trid_cyclic was modified to despise 
      !the last point of computational domain for periodic 
      !boundary conditions. Therefore the last point of vector  
      !r is equal zero as follow below.
      r(kmax)=0.d0

      RETURN
      END SUBROUTINE rhs_zz_6th_P

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
ccccc lhs_d1_6th_D: This routine storing the left hand side of linear 
ccccc               system of first derivatives. The 6th order compact 
ccccc               finite difference schemes was used. Dirichlet 
ccccc               boundary conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE lhs_d1_6th_D(aa,bb,cc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax)
  
      aa(1)=0.d0
      bb(1)=1.d0
      cc(1)=0.d0

      aa(2)=0.d0
      bb(2)=1.d0
      cc(2)=0.d0

      do l=3,lmax-2
        aa(l)=1.d0
        bb(l)=3.d0
        cc(l)=1.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d1_6th_D

ccccc **********************************************************************
ccccc lhs_d1_62th_D: This routine storing the left hand side of linear 
ccccc                system of first derivatives. The 6th and 2nd order  
ccccc                compact finite difference schemes was used. Dirichlet  
ccccc                boundary conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE lhs_d1_62th_D(aa,bb,cc,mtrc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc(lmax)

      aa(1)=0.d0
      bb(1)=1.d0/mtrc(1)*1.d0
      cc(1)=0.d0

      aa(2)=0.d0
      bb(2)=1.d0/mtrc(2)*1.d0
      cc(2)=0.d0

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
      END SUBROUTINE lhs_d1_62th_D

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
ccccc lhs_d1_6th_N: This routine storing the left hand side of linear 
ccccc               system of first derivatives. The 6th order compact 
ccccc               finite difference schemes was used. The Neumann 
ccccc               boundary conditions was adopted. A stretching is 
ccccc               used in this schemes.
ccccc **********************************************************************
      SUBROUTINE lhs_d1_6th_N(aa,bb,cc,mtrc,lmax)
                             
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
      END SUBROUTINE lhs_d1_6th_N

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
ccccc lhs_d2_6th_P: This routine storing the left hand side of linear 
ccccc               system of second derivatives. The 6th order compact 
ccccc               finite difference schemes was used. Periodic 
ccccc               boundary conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE lhs_d2_6th_P(aa,bb,cc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax,lmax_1
      real*8 aa(lmax),bb(lmax),cc(lmax)
  
      lmax_1=lmax-1

      aa(1)=00.d0 
      bb(1)=11.d0
      cc(1)=02.d0
  
      do l=2,lmax_1-1 
        aa(l)=02.d0
        bb(l)=11.d0
        cc(l)=02.d0
      end do
    
      aa(lmax_1)=02.d0 
      bb(lmax_1)=11.d0
      cc(lmax_1)=00.d0
    
      !The routine trid and trid_cyclic was modified to despise 
      !the last point of computational domain for periodic 
      !boundary conditions. Therefore the last point of vector  
      !aa, bb and cc is equal zero as follow below.
      aa(lmax)=0.d0 
      bb(lmax)=0.d0
      cc(lmax)=0.d0
    		
      RETURN
      END SUBROUTINE lhs_d2_6th_P

ccccc **********************************************************************
ccccc lhs_d2_6th_D: This routine storing the left hand side of linear 
ccccc                system of second derivatives. The 6th order compact  
ccccc                finite difference schemes was used for Dirichlet 
ccccc                boundary condition.
ccccc **********************************************************************
      SUBROUTINE lhs_d2_6th_D(aa,bb,cc,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax)
  
      aa(1)=0.d0
      bb(1)=1.d0
      cc(1)=0.d0

      aa(2)=0.d0 
      bb(2)=1.d0
      cc(2)=0.d0

      do l=3,lmax-2
        aa(l)=02.d0
        bb(l)=11.d0
        cc(l)=02.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d2_6th_D

ccccc **********************************************************************
ccccc lhs_d2_62th_D: This routine storing the left hand side of linear 
ccccc                system of second derivatives. The 6th and 2nd order 
ccccc                compact finite difference schemes was used for 
ccccc                Dirichlet boundary condition.
ccccc **********************************************************************
      SUBROUTINE lhs_d2_62th_D(aa,bb,cc,mtrc_d1,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc_d1(lmax)
  
      aa(1)=0.d0
      bb(1)=1.d0/(mtrc_d1(1)**2.d0)*1.d0
      cc(1)=0.d0

      aa(2)=0.d0
      bb(2)=1.d0/(mtrc_d1(2)**2.d0)*1.d0
      cc(2)=0.d0

      do l=3,lmax-2
        aa(l)=1.d0/(mtrc_d1(l-1)**2.d0)*02.d0
        bb(l)=1.d0/(mtrc_d1(l  )**2.d0)*11.d0
        cc(l)=1.d0/(mtrc_d1(l+1)**2.d0)*02.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0/(mtrc_d1(lmax-1)**2.d0)*1.d0
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0/(mtrc_d1(lmax)**2.d0)*1.d0
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d2_62th_D

ccccc **********************************************************************
ccccc lhs_d2_65th_W: This routine storing the left hand side of linear 
ccccc                system of second derivatives. The 6th and 5th order 
ccccc                compact finite difference schemes was used for 
ccccc                Dirichlet and Neumann boundary condition to simulate
ccccc                the wall.
ccccc **********************************************************************
      SUBROUTINE lhs_d2_65th_W(aa,bb,cc,mtrc_d1,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc_d1(lmax)
  
      aa(1)=0.d0
      bb(1)=1.d0/(mtrc_d1(1)**2.d0)*013.d0
      cc(1)=1.d0/(mtrc_d1(2)**2.d0)*137.d0

      aa(2)=1.d0/(mtrc_d1(1)**2.d0)*01.d0 
      bb(2)=1.d0/(mtrc_d1(2)**2.d0)*12.d0
      cc(2)=1.d0/(mtrc_d1(3)**2.d0)*03.d0

      do l=3,lmax-2
        aa(l)=1.d0/(mtrc_d1(l-1)**2.d0)*02.d0
        bb(l)=1.d0/(mtrc_d1(l  )**2.d0)*11.d0
        cc(l)=1.d0/(mtrc_d1(l+1)**2.d0)*02.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0/(mtrc_d1(lmax-1)**2.d0)*1.d0
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0/(mtrc_d1(lmax)**2.d0)*1.d0
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d2_65th_W

ccccc **********************************************************************
ccccc lhs_d2_6th_N: This routine storing the left hand side of linear 
ccccc               system of second derivatives. The 6th order compact  
ccccc               finite difference schemes was used for Neumann 
ccccc               boundary condition.
ccccc **********************************************************************
      SUBROUTINE lhs_d2_6th_N(aa,bb,cc,mtrc_d1,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc_d1(lmax)
  
      aa(1)=0.d0
      bb(1)=1.d0/(mtrc_d1(1)**2.d0)*1.d0
      cc(1)=0.d0

      aa(2)=0.d0 
      bb(2)=1.d0/(mtrc_d1(2)**2.d0)*1.d0
      cc(2)=0.d0

      do l=3,lmax-2
        aa(l)=1.d0/(mtrc_d1(l-1)**2.d0)*02.d0
        bb(l)=1.d0/(mtrc_d1(l  )**2.d0)*11.d0
        cc(l)=1.d0/(mtrc_d1(l+1)**2.d0)*02.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0/(mtrc_d1(lmax-1)**2.d0)*1.d0
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0/(mtrc_d1(lmax)**2.d0)*1.d0
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d2_6th_N

ccccc **********************************************************************
ccccc lhs_d2_65th_NN: This routine storing the left hand side of linear 
ccccc                 system of second derivatives. The 6th and 5th order 
ccccc                 compact finite difference schemes was used for 
ccccc                 Neumann boundary condition.
ccccc **********************************************************************
      SUBROUTINE lhs_d2_65th_NN(aa,bb,cc,mtrc_d1,lmax)
                             
      IMPLICIT NONE
      
      integer l,lmax
      real*8 aa(lmax),bb(lmax),cc(lmax),mtrc_d1(lmax)

      aa(1)=0.d0                                  
      bb(1)=1.d0/(mtrc_d1(1)**2.d0)*013.d0         
      cc(1)=1.d0/(mtrc_d1(2)**2.d0)*137.d0        

      aa(2)=1.d0/(mtrc_d1(1)**2.d0)*01.d0          
      bb(2)=1.d0/(mtrc_d1(2)**2.d0)*12.d0         
      cc(2)=1.d0/(mtrc_d1(3)**2.d0)*03.d0          

      do l=3,lmax-2
        aa(l)=1.d0/(mtrc_d1(l-1)**2.d0)*02.d0
        bb(l)=1.d0/(mtrc_d1(l  )**2.d0)*11.d0
        cc(l)=1.d0/(mtrc_d1(l+1)**2.d0)*02.d0
      end do

      aa(lmax-1)=0.d0
      bb(lmax-1)=1.d0/(mtrc_d1(lmax-1)**2.d0)*1.d0
      cc(lmax-1)=0.d0

      aa(lmax)=0.d0
      bb(lmax)=1.d0/(mtrc_d1(lmax)**2.d0)*1.d0
      cc(lmax)=0.d0

      RETURN
      END SUBROUTINE lhs_d2_65th_NN

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
ccccc stencil_x_ex_6th_D: This routine calculates the stencils for 1st 
ccccc                     derivatives in x direction with a 6th order 
ccccc                     explicit schemes. The non-periodic boundary 
ccccc                     conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE stencil_x_ex_6th_D(my_rank,pro,imax,x,stx)

      IMPLICIT NONE

      integer i,imax,my_rank,pro
      real*8 x(imax),stx(imax,7)

      if (my_rank.EQ.0) then 
        do i=1,3
          call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),
     &      x(i),1,stx(i,1))
          call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),
     &      x(i),2,stx(i,2))
          call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),
     &      x(i),3,stx(i,3))
          call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),
     &      x(i),4,stx(i,4))
          call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),
     &      x(i),5,stx(i,5))
          call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),
     &      x(i),6,stx(i,6))
          call coeffs_d1_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),x(7),
     &      x(i),7,stx(i,7))
        end do
      end if 

      do i=4,imax-3
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),1,stx(i,1))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),2,stx(i,2))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),3,stx(i,3))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),4,stx(i,4))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),5,stx(i,5))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),6,stx(i,6))
        call coeffs_d1_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),7,stx(i,7))   
      end do

      if (my_rank+1.EQ.pro) then
        do i=imax-2,imax
          call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),1,stx(i,1))
          call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),2,stx(i,2))
          call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),3,stx(i,3))
          call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),4,stx(i,4))
          call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),5,stx(i,5))
          call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),6,stx(i,6))
          call coeffs_d1_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),7,stx(i,7))
        end do
      end if

      RETURN
      END SUBROUTINE stencil_x_ex_6th_D 

ccccc **********************************************************************
ccccc stencil_xx_ex_6th_D: This routine calculates the stencils for 2nd 
ccccc                      derivatives in x direction with a 6th order 
ccccc                      explicit schemes. The non-periodic boundary 
ccccc                      conditions was adopted.
ccccc **********************************************************************
      SUBROUTINE stencil_xx_ex_6th_D(my_rank,pro,imax,x,stxx)

      IMPLICIT NONE

      integer i,imax,my_rank,pro
      real*8 x(imax),stxx(imax,7)

      if (my_rank.EQ.0) then
        do i=1,3
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),
     &      x(7),x(i),1,stxx(i,1))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),
     &      x(7),x(i),2,stxx(i,2))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),
     &      x(7),x(i),3,stxx(i,3))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),
     &      x(7),x(i),4,stxx(i,4))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),
     &      x(7),x(i),5,stxx(i,5))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),
     &      x(7),x(i),6,stxx(i,6))
          call coeffs_d2_ex_6th(x(1),x(2),x(3),x(4),x(5),x(6),
     &      x(7),x(i),7,stxx(i,7))   
        end do
      end if

      do i=4,imax-3
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),1,stxx(i,1))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),2,stxx(i,2))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),3,stxx(i,3))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),4,stxx(i,4))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),5,stxx(i,5))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),6,stxx(i,6))
        call coeffs_d2_ex_6th(x(i-3),x(i-2),x(i-1),x(i),x(i+1),x(i+2),
     &    x(i+3),x(i),7,stxx(i,7))   
      end do

      if (my_rank+1.EQ.pro) then
        do i=imax-2,imax
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),1,stxx(i,1))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),2,stxx(i,2))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),3,stxx(i,3))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),4,stxx(i,4))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),5,stxx(i,5))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),6,stxx(i,6))
          call coeffs_d2_ex_6th(x(imax-6),x(imax-5),x(imax-4),x(imax-3),
     &      x(imax-2),x(imax-1),x(imax),x(i),7,stxx(i,7))
        end do
      end if 

      RETURN
      END SUBROUTINE stencil_xx_ex_6th_D

ccccc **********************************************************************
ccccc trid: Solve the linear system.
ccccc ----------------------------------------------------------------------
ccccc fct.....: Factor (0 - Non Periodic; 1 - Periodic) 
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

      do i=lmax-fct-1,1,-1
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

      IMPLICIT NONE

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
