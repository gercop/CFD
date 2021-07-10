!-----------------------------------------------------------------------------------------------------------------------!
! module which defines and calculates grid information.
!-----------------------------------------------------------------------------------------------------------------------!
module grid

#include "setup.h"
! include "setup.h"

! modules
   use params
   use parallel_mpi
   implicit none

! variables
   real, dimension(nx), public, save :: pgrid_x, cgrid_x
   real, dimension(ny), public, save :: pgrid_y, cgrid_y
#if (SPACE_INT == 2)
   real, public, save :: scale_x, scale_y
   real, dimension(nx,7), public, save :: st_coef_x
   real, dimension(ny,7), public, save :: st_coef_y
#endif ! SPACE_INT

! functions
#if (SPACE_INT == 2)
   private :: fd1
#endif ! SPACE_INT

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: Calculates grids.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine grid_init

!----- Variables -------------------------------------------------------------------------------------------------------!
! local variables
   integer :: i


#if (SPACE_INT == 2)
!----- Calculates scale factors for derivatives ------------------------------------------------------------------------!
! x-direction
   if (x_L.le.x_0) then
      print *,'ERROR->GRID: Outflow has to be downstream of inflow (x_L > x_0).'
      stop
   end if
   scale_x = 2.0/(x_L - x_0)

! y-direction
   if (y_H.le.y_0) then
      print *,'ERROR->GRID: Wall has to be below free stream (y_H > y_0).'
      stop
   end if
   scale_y = 2.0/(y_H - y_0)
#endif ! SPACE_INT

!----- Calculates computational grid -----------------------------------------------------------------------------------!
! x-direction
   do i=1,local_nx
#if (SPACE_INT == 2)
      cgrid_x(i) = asin(-alpha*dcos(pi*dble(i-1)/dble(nx-1)))/asin(alpha)
#elif (SPACE_INT == 1)
      cgrid_x(i) = (i-1)*dx
#endif ! SPACE_INT
   end do

! y-direction
   do i=1,local_ny
#if (SPACE_INT == 2)
      cgrid_y(i) = asin(-alpha*dcos(pi*dble(i-1)/dble(ny-1)))/asin(alpha)
#elif (SPACE_INT == 1)
      cgrid_y(i) = (i-1)*dy
#endif ! SPACE_INT
   end do

!----- Calculates physical grid ----------------------------------------------------------------------------------------!
! x-direction
   do i=1,local_nx
#if (SPACE_INT == 2)
      pgrid_x(i) = (x_L - x_0)/2.0*cgrid_x(i)+(x_L + x_0)/2.0
#elif (SPACE_INT == 1)
      pgrid_x(i) = local_x_0 + cgrid_x(i)
#endif ! SPACE_INT
   end do

! y-direction
   do i=1,local_ny
#if (SPACE_INT == 2)
      pgrid_y(i) = (y_H - y_0)/2.0*cgrid_y(i)+(y_H + y_0)/2.0
#elif (SPACE_INT == 1)
      pgrid_y(i) = local_y_0 + cgrid_y(i)
#endif ! SPACE_INT
   end do

#if (SPACE_INT == 2)
!----- Calculates 6th order stencil coefficients -----------------------------------------------------------------------!
! x-direction
   do i=1,4
      st_coef_x(i,1) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),1)
      st_coef_x(i,2) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),2)
      st_coef_x(i,3) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),3)
      st_coef_x(i,4) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),4)
      st_coef_x(i,5) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),5)
      st_coef_x(i,6) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),6)
      st_coef_x(i,7) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),7)
   end do

   do i=5,nx-4
      st_coef_x(i,1) = fd1(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                           cgrid_x(i),1)
      st_coef_x(i,2) = fd1(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                           cgrid_x(i),2)
      st_coef_x(i,3) = fd1(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                           cgrid_x(i),3)
      st_coef_x(i,4) = fd1(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                           cgrid_x(i),4)
      st_coef_x(i,5) = fd1(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                           cgrid_x(i),5)
      st_coef_x(i,6) = fd1(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                           cgrid_x(i),6)
      st_coef_x(i,7) = fd1(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                           cgrid_x(i),7)   
   end do

   do i=nx-3,nx
      st_coef_x(i,1) = fd1(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                           cgrid_x(nx),cgrid_x(i),1)
      st_coef_x(i,2) = fd1(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                           cgrid_x(nx),cgrid_x(i),2)
      st_coef_x(i,3) = fd1(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                           cgrid_x(nx),cgrid_x(i),3)
      st_coef_x(i,4) = fd1(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                           cgrid_x(nx),cgrid_x(i),4)
      st_coef_x(i,5) = fd1(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                           cgrid_x(nx),cgrid_x(i),5)
      st_coef_x(i,6) = fd1(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                           cgrid_x(nx),cgrid_x(i),6)
      st_coef_x(i,7) = fd1(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                           cgrid_x(nx),cgrid_x(i),7)
   end do

! y-direction
   do i=1,4
      st_coef_y(i,1) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),1)
      st_coef_y(i,2) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),2)
      st_coef_y(i,3) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),3)
      st_coef_y(i,4) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),4)
      st_coef_y(i,5) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),5)
      st_coef_y(i,6) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),6)
      st_coef_y(i,7) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),7)
   end do

   do i=5,ny-4
      st_coef_y(i,1) = fd1(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                           cgrid_y(i),1)
      st_coef_y(i,2) = fd1(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                           cgrid_y(i),2)
      st_coef_y(i,3) = fd1(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                           cgrid_y(i),3)
      st_coef_y(i,4) = fd1(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                           cgrid_y(i),4)
      st_coef_y(i,5) = fd1(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                           cgrid_y(i),5)
      st_coef_y(i,6) = fd1(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                           cgrid_y(i),6)
      st_coef_y(i,7) = fd1(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                           cgrid_y(i),7)   
   end do

   do i=ny-3,ny
      st_coef_y(i,1) = fd1(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                           cgrid_y(ny),cgrid_y(i),1)
      st_coef_y(i,2) = fd1(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                           cgrid_y(ny),cgrid_y(i),2)
      st_coef_y(i,3) = fd1(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                           cgrid_y(ny),cgrid_y(i),3)
      st_coef_y(i,4) = fd1(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                           cgrid_y(ny),cgrid_y(i),4)
      st_coef_y(i,5) = fd1(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                           cgrid_y(ny),cgrid_y(i),5)
      st_coef_y(i,6) = fd1(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                           cgrid_y(ny),cgrid_y(i),6)
      st_coef_y(i,7) = fd1(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                           cgrid_y(ny),cgrid_y(i),7)
   end do
#endif ! SPACE_INT

   end subroutine grid_init

#if (SPACE_INT == 2)
!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: Function to get the coefficients for the 6th order stencils.
!-----------------------------------------------------------------------------------------------------------------------!
   function fd1(x1,x2,x3,x4,x5,x6,x7,x,which)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input variables
   integer, intent(in) :: which
   real, intent(in) :: x1,x2,x3,x4,x5,x6,x7,x

! output variable
   real :: fd1

   if(which.eq.1) then
   fd1 =((x - x2)*(x - x3)*(x - x4)*(x - x5)*(x - x6))/                  &
         ((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        ((x - x2)*(x - x3)*(x - x4)*(x - x5)*(x - x7))/                  &
         ((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        ((x - x2)*(x - x3)*(x - x4)*(x - x6)*(x - x7))/                  &
         ((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        ((x - x2)*(x - x3)*(x - x5)*(x - x6)*(x - x7))/                  &
         ((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        ((x - x2)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/                  &
         ((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        ((x - x3)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/                  &
         ((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7))
   else if(which.eq.2) then
   fd1 =((x - x1)*(x - x3)*(x - x4)*(x - x5)*(x - x6))/                   &
         ((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        ((x - x1)*(x - x3)*(x - x4)*(x - x5)*(x - x7))/                   &
         ((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        ((x - x1)*(x - x3)*(x - x4)*(x - x6)*(x - x7))/                   &
         ((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        ((x - x1)*(x - x3)*(x - x5)*(x - x6)*(x - x7))/                   &
         ((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        ((x - x1)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/                   &
         ((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        ((x - x3)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/                   &
         ((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7))
   else if(which.eq.3) then
   fd1 =((x - x1)*(x - x2)*(x - x4)*(x - x5)*(x - x6))/                    &
         ((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        ((x - x1)*(x - x2)*(x - x4)*(x - x5)*(x - x7))/                    &
         ((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        ((x - x1)*(x - x2)*(x - x4)*(x - x6)*(x - x7))/                    &
         ((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        ((x - x1)*(x - x2)*(x - x5)*(x - x6)*(x - x7))/                    &
         ((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        ((x - x1)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/                    &
         ((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        ((x - x2)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/                    &
         ((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7))
   else if(which.eq.4) then
   fd1 =((x - x1)*(x - x2)*(x - x3)*(x - x5)*(x - x6))/                     &
         ((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        ((x - x1)*(x - x2)*(x - x3)*(x - x5)*(x - x7))/                     &
         ((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        ((x - x1)*(x - x2)*(x - x3)*(x - x6)*(x - x7))/                     &
         ((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        ((x - x1)*(x - x2)*(x - x5)*(x - x6)*(x - x7))/                     &
         ((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        ((x - x1)*(x - x3)*(x - x5)*(x - x6)*(x - x7))/                     &
         ((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        ((x - x2)*(x - x3)*(x - x5)*(x - x6)*(x - x7))/                     &
         ((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7))
   else if(which.eq.5) then
   fd1 =((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x6))/                      &
         ((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        ((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x7))/                      &
         ((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        ((x - x1)*(x - x2)*(x - x3)*(x - x6)*(x - x7))/                      &
         ((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        ((x - x1)*(x - x2)*(x - x4)*(x - x6)*(x - x7))/                      &
         ((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        ((x - x1)*(x - x3)*(x - x4)*(x - x6)*(x - x7))/                      &
         ((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        ((x - x2)*(x - x3)*(x - x4)*(x - x6)*(x - x7))/                      &
         ((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) 
   else if(which.eq.6) then
   fd1 =((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x5))/                       &
         ((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        ((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x7))/                       &
         ((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        ((x - x1)*(x - x2)*(x - x3)*(x - x5)*(x - x7))/                       &
         ((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        ((x - x1)*(x - x2)*(x - x4)*(x - x5)*(x - x7))/                       &
         ((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        ((x - x1)*(x - x3)*(x - x4)*(x - x5)*(x - x7))/                       &
         ((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        ((x - x2)*(x - x3)*(x - x4)*(x - x5)*(x - x7))/                       &
         ((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7))
   else if(which.eq.7) then
   fd1 =((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x5))/                        &
         ((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        ((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x6))/                        &
         ((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        ((x - x1)*(x - x2)*(x - x3)*(x - x5)*(x - x6))/                        &
         ((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        ((x - x1)*(x - x2)*(x - x4)*(x - x5)*(x - x6))/                        &
         ((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        ((x - x1)*(x - x3)*(x - x4)*(x - x5)*(x - x6))/                        &
         ((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        ((x - x2)*(x - x3)*(x - x4)*(x - x5)*(x - x6))/                        &
         ((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7))
   end if

   end function fd1    
#endif ! SPACE_INT

end module grid
