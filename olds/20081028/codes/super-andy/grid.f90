!----------------------------------------------------------------------------------------------------------------------------------!
! module which defines and calculates grid information.
!----------------------------------------------------------------------------------------------------------------------------------!
module grid

#include 'setup.h'

!----- Modules --------------------------------------------------------------------------------------------------------------------!
#if defined(DEBUG)
   use debug
#endif ! DEBUG
   use params

   implicit none

!----- Constants ------------------------------------------------------------------------------------------------------------------!
! stencil half width around central point for interpolation of initial condition.
#if (SPACE_INT==1)
   integer, parameter, public                     :: st_w = 3
#else
   integer, parameter, public                     :: st_w = 7
#endif ! SPACE_INT

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! stencil width and stencil mid point
   integer, public          		          :: st_hw, st_mid

   real, dimension(nx), public, save              :: pgrid_x, cgrid_x
   real, dimension(ny), public, save              :: pgrid_y, cgrid_y
   real, public, save                             :: x_scale_1st, y_scale_1st, xoy_scale_1st, x_scale_2nd, y_scale_2nd
   real, dimension(nx,st_w), public, save         :: x_coef_1st, x_coef_2nd
   real, dimension(ny,st_w), public, save         :: y_coef_1st, y_coef_2nd

   real, dimension(nx,st_w), public, save         :: x_coef_1st_1, x_coef_1st_2, x_coef_1st_3
   real, dimension(ny,st_w), public, save         :: y_coef_1st_1, y_coef_1st_2, y_coef_1st_3

!----- Functions ------------------------------------------------------------------------------------------------------------------!
   private                                        :: poly_2, poly_6, fd1, fd2

   interface poly
     module procedure poly_2
     module procedure poly_6
   end interface poly

contains

!----------------------------------------------------------------------------------------------------------------------------------!
! first subroutine: Calculates grids.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine grid_init

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! local variables
   integer                                        :: i, j

!----- Calculates scale factors for derivatives -----------------------------------------------------------------------------------!
! x-direction
   if (x_L.le.x_0) then
      print *,'ERROR->GRID: Outflow has to be downstream of inflow (x_L > x_0).'
      stop
   end if
   x_scale_1st =  2.0d0/(x_L - x_0)
   x_scale_2nd = x_scale_1st**2

! y-direction
   if (y_H.le.y_0) then
      print *,'ERROR->GRID: Wall has to be below free stream (y_H > y_0).'
      stop
   end if
   y_scale_1st =  2.0d0/(y_H - y_0)
   y_scale_2nd = y_scale_1st**2

! x over y
   xoy_scale_1st = x_scale_1st/y_scale_1st

!----- Calculates computational grid ----------------------------------------------------------------------------------------------!
! x-direction
   do i=1,nx
      cgrid_x(i) = dasin(-alpha*dcos(pi*dble(i-1)/dble(nx-1)))/dasin(alpha)
   end do

! y-direction
   do i=1,ny
      cgrid_y(i) = dasin(-alpha*dcos(pi*dble(i-1)/dble(ny-1)))/dasin(alpha)
   end do

!----- Calculates physical grid ---------------------------------------------------------------------------------------------------!

open(12,file='grid.dat',status='unknown')
! x-direction
   do i=1,nx
      pgrid_x(i) = (x_L - x_0)/2.0d0*cgrid_x(i)+(x_L + x_0)/2.0d0
   end do

! y-direction
   do i=1,ny
      pgrid_y(i) = (y_H - y_0)/2.0d0*cgrid_y(i)+(y_H + y_0)/2.0d0
   end do

open(12,file='grid.dat',status='unknown')
   do i=2,nx
      write(12,*) pgrid_x(i)-pgrid_x(i-1)
   end do
   write(12,*) '&'
   do i=2,ny
      write(12,*) pgrid_y(i)-pgrid_y(i-1)
   end do

close(12)

!----- Calculates some stencil infos ----------------------------------------------------------------------------------------------!
   st_hw  = INT((st_w-1)/2)
   st_mid = st_hw+1

!----- Calculates stencil coefficients --------------------------------------------------------------------------------------------!

#if (SPACE_INT==1)
   do i = 1,nx-2
     do j = 1,st_w  
       x_coef_1st_1(i,j) = fd1(cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i),j)
     end do
   end do
   do i = 2,nx-1
     do j = 1,st_w  
       x_coef_1st_2(i,j) = fd1(cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i),j)
     end do
   end do
   do i = 3,nx
     do j = 1,st_w  
       x_coef_1st_3(i,j) = fd1(cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i),j)
     end do
   end do

   do i = 1,ny-2
     do j = 1,st_w  
       y_coef_1st_1(i,j) = fd1(cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i),j)
     end do
   end do
   do i = 2,ny-1
     do j = 1,st_w  
       y_coef_1st_2(i,j) = fd1(cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i),j)
     end do
   end do
   do i = 3,ny
     do j = 1,st_w  
       y_coef_1st_3(i,j) = fd1(cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i),j)
     end do
   end do
#endif !SPACE_INT


! x-direction
   do i=1,1+st_hw
! first derivative
      do j=1,st_w
#if (SPACE_INT==1)
         x_coef_1st(i,j) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(i),j)
#else
         x_coef_1st(i,j) = fd1(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),j)
#endif ! SPACE_INT
      end do
! second derivative
      do j=1,st_w
#if (SPACE_INT==1)
         x_coef_2nd(i,j) = fd2(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(i),j)
#else
         x_coef_2nd(i,j) = fd2(cgrid_x(1),cgrid_x(2),cgrid_x(3),cgrid_x(4),cgrid_x(5),cgrid_x(6),cgrid_x(7),cgrid_x(i),j)
#endif ! SPACE_INT
      end do
   end do

   do i=2+st_hw,nx-st_hw-1
! first derivative
      do j=1,st_w
#if (SPACE_INT==1)
         x_coef_1st(i,j) = fd1(cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i),j)
#else
         x_coef_1st(i,j) = fd1(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                               cgrid_x(i),j)
#endif ! SPACE_INT
      end do
! second derivative
      do j=1,st_w
#if (SPACE_INT==1)
         x_coef_2nd(i,j) = fd2(cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i),j)
#else
         x_coef_2nd(i,j) = fd2(cgrid_x(i-3),cgrid_x(i-2),cgrid_x(i-1),cgrid_x(i),cgrid_x(i+1),cgrid_x(i+2),cgrid_x(i+3),&
                               cgrid_x(i),j)
#endif ! SPACE_INT
      end do
   end do

   do i=nx-st_hw,nx
! first derivative
      do j=1,st_w
#if (SPACE_INT==1)
         x_coef_1st(i,j) = fd1(cgrid_x(nx-2),cgrid_x(nx-1),cgrid_x(nx),cgrid_x(i),j)
#else
         x_coef_1st(i,j) = fd1(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                               cgrid_x(nx),cgrid_x(i),j)
#endif ! SPACE_INT
      end do
! second derivative
      do j=1,st_w
#if (SPACE_INT==1)
         x_coef_2nd(i,j) = fd2(cgrid_x(nx-2),cgrid_x(nx-1),cgrid_x(nx),cgrid_x(i),j)
#else
         x_coef_2nd(i,j) = fd2(cgrid_x(nx-6),cgrid_x(nx-5),cgrid_x(nx-4),cgrid_x(nx-3),cgrid_x(nx-2),cgrid_x(nx-1),&
                               cgrid_x(nx),cgrid_x(i),j)
#endif ! SPACE_INT
      end do
   end do

! y-direction
   do i=1,1+st_hw
! first derivative
      do j=1,st_w
#if (SPACE_INT==1)
         y_coef_1st(i,j) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(i),j)
#else
         y_coef_1st(i,j) = fd1(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),j)
#endif ! SPACE_INT
      end do
! second derivative
      do j=1,st_w
#if (SPACE_INT==1)
         y_coef_2nd(i,j) = fd2(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(i),j)
#else
         y_coef_2nd(i,j) = fd2(cgrid_y(1),cgrid_y(2),cgrid_y(3),cgrid_y(4),cgrid_y(5),cgrid_y(6),cgrid_y(7),cgrid_y(i),j)
#endif ! SPACE_INT
      end do
   end do

   do i=2+st_hw,ny-st_hw-1
! first derivative
      do j=1,st_w
#if (SPACE_INT==1)
         y_coef_1st(i,j) = fd1(cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i),j)
#else
         y_coef_1st(i,j) = fd1(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                               cgrid_y(i),j)
#endif ! SPACE_INT
      end do
! second derivative
      do j=1,st_w
#if (SPACE_INT==1)
         y_coef_2nd(i,j) = fd2(cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i),j)
#else
         y_coef_2nd(i,j) = fd2(cgrid_y(i-3),cgrid_y(i-2),cgrid_y(i-1),cgrid_y(i),cgrid_y(i+1),cgrid_y(i+2),cgrid_y(i+3),&
                               cgrid_y(i),j)
#endif ! SPACE_INT
      end do
   end do

   do i=ny-st_hw,ny
! first derivative
      do j=1,st_w
#if (SPACE_INT==1)
         y_coef_1st(i,j) = fd1(cgrid_y(ny-2),cgrid_y(ny-1),cgrid_y(ny),cgrid_y(i),j)
#else
         y_coef_1st(i,j) = fd1(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                               cgrid_y(ny),cgrid_y(i),j)
#endif ! SPACE_INT
      end do
! second derivative
      do j=1,st_w
#if (SPACE_INT==1)
         y_coef_2nd(i,j) = fd2(cgrid_y(ny-2),cgrid_y(ny-1),cgrid_y(ny),cgrid_y(i),j)
#else
         y_coef_2nd(i,j) = fd2(cgrid_y(ny-6),cgrid_y(ny-5),cgrid_y(ny-4),cgrid_y(ny-3),cgrid_y(ny-2),cgrid_y(ny-1),&
                               cgrid_y(ny),cgrid_y(i),j)
#endif ! SPACE_INT
      end do
   end do

   end subroutine grid_init

!----------------------------------------------------------------------------------------------------------------------------------!
! second subroutine: Function to get the coefficients for a 2nd order polynomial.
!----------------------------------------------------------------------------------------------------------------------------------!
   function poly_2(x1,x2,x3,x,which)

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! input variables
   integer, intent(in)                            :: which
   real, intent(in)                               :: x1,x2,x3,x

! output variable
   real                                           :: poly_2

!----- Main -----------------------------------------------------------------------------------------------------------------------!
   select case (which)
   case(1)
      poly_2 =((x - x2)*(x - x3))/(( x1 - x2)*( x1 - x3))
   case(2)
      poly_2 =((x - x1)*(x - x3))/((-x1 + x2)*( x2 - x3))
   case(3)
      poly_2 =((x - x1)*(x - x2))/((-x1 + x3)*(-x2 + x3))
   end select

   end function poly_2

!----------------------------------------------------------------------------------------------------------------------------------!
! third subroutine: Function to get the coefficients for a 6th order polynomial.
!----------------------------------------------------------------------------------------------------------------------------------!
   function poly_6(x1,x2,x3,x4,x5,x6,x7,x,which)

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! input variables
   integer, intent(in)                            :: which
   real, intent(in)                               :: x1,x2,x3,x4,x5,x6,x7,x

! output variable
   real                                           :: poly_6

!----- Main -----------------------------------------------------------------------------------------------------------------------!
   select case (which)
   case(1)
      poly_6 =((x - x2)*(x - x3)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/ &
              ((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7))
   case(2)
      poly_6 =((x - x1)*(x - x3)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/ &
              ((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7))
   case(3)
      poly_6 =((x - x1)*(x - x2)*(x - x4)*(x - x5)*(x - x6)*(x - x7))/ &
              ((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7))
   case(4)
      poly_6 =((x - x1)*(x - x2)*(x - x3)*(x - x5)*(x - x6)*(x - x7))/ &
              ((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7))
   case(5)
      poly_6 =((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x6)*(x - x7))/ &
              ((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7))
   case(6)
      poly_6 =((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x5)*(x - x7))/ &
              ((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7))
   case(7)
      poly_6 =((x - x1)*(x - x2)*(x - x3)*(x - x4)*(x - x5)*(x - x6))/ &
              ((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7))
   end select

   end function poly_6
#if (SPACE_INT==1)
!----------------------------------------------------------------------------------------------------------------------------------!
! fourth subroutine: Function to get the coefficients for the 2nd order stencils of a first derivative.
!----------------------------------------------------------------------------------------------------------------------------------!
   function fd1(x1,x2,x3,x,which)

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! input variables
   integer, intent(in)                            :: which
   real, intent(in)                               :: x1,x2,x3,x

! output variable
   real                                           :: fd1

!----- Main -----------------------------------------------------------------------------------------------------------------------!
   select case (which)
   case(1)
      fd1 = (x - x2)/(( x1 - x2)*( x1 - x3)) + (x - x3)/(( x1 - x2)*( x1 - x3))
   case(2)
      fd1 = (x - x1)/((-x1 + x2)*( x2 - x3)) + (x - x3)/((-x1 + x2)*( x2 - x3))
   case(3)
      fd1 = (x - x1)/((-x1 + x3)*(-x2 + x3)) + (x - x2)/((-x1 + x3)*(-x2 + x3))
   end select

   end function fd1
#else
!----------------------------------------------------------------------------------------------------------------------------------!
! fifth subroutine: Function to get the coefficients for the 6th order stencils of a first derivative.
!----------------------------------------------------------------------------------------------------------------------------------!
   function fd1(x1,x2,x3,x4,x5,x6,x7,x,which)

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! input variables
   integer, intent(in)                            :: which
   real, intent(in)                               :: x1,x2,x3,x4,x5,x6,x7,x

! output variable
   real                                           :: fd1

!----- Main -----------------------------------------------------------------------------------------------------------------------!
   select case (which)
   case(1)
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
   case(2)
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
   case(3)
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
   case(4)
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
   case(5)
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
   case(6)
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
   case(7)
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
   end select

   end function fd1
#endif !SPACE_INT
!----------------------------------------------------------------------------------------------------------------------------------!
! sixth subroutine: Function to get the coefficients for the 6th order stencils of a second derivative.
!----------------------------------------------------------------------------------------------------------------------------------!
#if (SPACE_INT==1)
   function fd2(x1,x2,x3,x,which)

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! input variables
   integer, intent(in)                            :: which
   real, intent(in)                               :: x1,x2,x3,x

! output variable
   real                                           :: fd2

!----- Main -----------------------------------------------------------------------------------------------------------------------!
   select case (which)
   case(1)
      fd2 = 2.0/(( x1 - x2)*( x1 - x3))
   case(2)
      fd2 = 2.0/((-x1 + x2)*( x2 - x3))
   case(3)
      fd2 = 2.0/((-x1 + x3)*(-x2 + x3))
   end select

   end function fd2


#else
   function fd2(x1,x2,x3,x4,x5,x6,x7,x,which)

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! input variables
   integer, intent(in)                            :: which
   real, intent(in)                               :: x1,x2,x3,x4,x5,x6,x7,x

! output variable
   real                                           :: fd2

!----- Main -----------------------------------------------------------------------------------------------------------------------!
   select case (which)
   case(1)
   fd2 =(2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x5))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x6))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x5)*(x - x6))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x5)*(x - x6))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x5)*(x - x6))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x5)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x5)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x5)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x6)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x6)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x6)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x2)*(x - x5)*(x - x6)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x3)*(x - x5)*(x - x6)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7)) + &
        (2.0d0*(x - x4)*(x - x5)*(x - x6)*(x - x7))/((x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)*(x1 - x6)*(x1 - x7))
   case(2)
   fd2 =(2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x5))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x6))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x5)*(x - x6))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x5)*(x - x6))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x5)*(x - x6))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x5)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x5)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x5)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x6)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x6)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x6)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x1)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x3)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7)) + &
        (2.0d0*(x - x4)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x2)*(x2 - x3)*(x2 - x4)*(x2 - x5)*(x2 - x6)*(x2 - x7))
   case(3)
   fd2 =(2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x5))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x6))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x5)*(x - x6))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x5)*(x - x6))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x5)*(x - x6))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x5)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x5)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x5)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x6)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x6)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x6)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x1)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x2)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7)) + &
        (2.0d0*(x - x4)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x3)*(-x2 + x3)*(x3 - x4)*(x3 - x5)*(x3 - x6)*(x3 - x7))
   case(4)
   fd2 =(2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x5))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x6))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x5)*(x - x6))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x5)*(x - x6))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x5)*(x - x6))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x5)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x5)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x5)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x6)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x6)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x6)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x1)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x2)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7)) + &
        (2.0d0*(x - x3)*(x - x5)*(x - x6)*(x - x7))/((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*(x4 - x5)*(x4 - x6)*(x4 - x7))
   case(5)
   fd2 =(2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x4))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x6))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x6))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x6))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x6))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x6)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x6)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x6)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x6)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x6)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x6)*(x - x7))/((-x1 + x5)*(-x2 + x5)*(-x3 + x5)*(-x4 + x5)*(x5 - x6)*(x5 - x7)) 
   case(6)
   fd2 =(2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x4))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x5))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x5))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x5))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x5))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x5)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x5)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x5)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x5)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x5)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x5)*(x - x7))/((-x1 + x6)*(-x2 + x6)*(-x3 + x6)*(-x4 + x6)*(-x5 + x6)*(x6 - x7))
   case(7)
   fd2 =(2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x4))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x5))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x5))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x5))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x5))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x3)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x4)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x4)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x4)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x2)*(x - x5)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x3)*(x - x5)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x2)*(x - x3)*(x - x5)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x1)*(x - x4)*(x - x5)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x2)*(x - x4)*(x - x5)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7)) + &
        (2.0d0*(x - x3)*(x - x4)*(x - x5)*(x - x6))/((-x1 + x7)*(-x2 + x7)*(-x3 + x7)*(-x4 + x7)*(-x5 + x7)*(-x6 + x7))
   end select

   end function fd2
#endif !SPACE_INT

end module grid
