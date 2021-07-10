module space_integration

# include "setup.h"

! modules
   use params, only: nx, ny
   use grid
   use parallel_mpi
   
   implicit none

contains

#if (SPACE_INT == 1)
!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: first order b/f differences + second order central differences.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine bfc_diff(input, idh, ify, iby, jfx, jbx, output)

! backward finite differences:
!      in x: ify = 0, jfx = 0, iby = 0, jbx = 1
!      in y: ify = 0, jfx = 0, iby = 1, jbx = 0
!
! forward finite differences:
!      in x: ify = 0, jfx = 1, iby = 0, jbx = 0
!      in y: ify = 1, jfx = 0, iby = 0, jbx = 0
!
! central finite differences:
!      in x: ify = 0, jfx = 1, iby = 0, jbx = 1, dh has to be 2dx
!      in y: ify = 1, jfx = 0, iby = 1, jbx = 0, dh has to be 2dy
!
! There is no error check for jfx, ify, jbx, and iby since computer time has to be saved. Be sure what you are doing!

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   integer, intent(in) :: ify, iby, jfx, jbx
   real, intent(in) :: idh
   real, dimension(local_ny,local_nx), intent(in) :: input

! output
   real, dimension(local_ny,local_nx), intent(out) :: output

! local variables x->j, y->i
   integer :: i, j

!----- Derivative ------------------------------------------------------------------------------------------------------!
   do i=(1+iby),(local_ny-ify)
      do j=(1+jbx),(local_nx-jfx)
	 output(i,j) = idh*(input(i+ify,j+jfx)-input(i-iby,j-jbx))
      end do
   end do

   end subroutine bfc_diff

#elif (SPACE_INT == 2)
!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: 6th order finite differences (a la Zhong) in x-direction.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine zhong_diff_x(input, output)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, dimension(ny,nx), intent(in) :: input

! output
   real, dimension(ny,nx), intent(out) :: output

! local variables x->j, y->i
   integer :: i, j

!----- Derivative ------------------------------------------------------------------------------------------------------!
   do i=1,ny
      do j=1,4
         output(i,j) = scale_x*( st_coef_x(j,1)*input(i,1)+st_coef_x(j,2)*input(i,2)+st_coef_x(j,3)*input(i,3)+ &
                                 st_coef_x(j,4)*input(i,4)+st_coef_x(j,5)*input(i,5)+st_coef_x(j,6)*input(i,6)+ &
                                 st_coef_x(j,7)*input(i,7) )
      end do
      do j=5,nx-4
         output(i,j) = scale_x*( st_coef_x(j,1)*input(i,j-3)+st_coef_x(j,2)*input(i,j-2)+st_coef_x(j,3)*input(i,j-1)+ &
                                 st_coef_x(j,4)*input(i,j)  +st_coef_x(j,5)*input(i,j+1)+st_coef_x(j,6)*input(i,j+2)+ &
                                 st_coef_x(j,7)*input(i,j+3) )
      end do
      do j=nx-3,nx
         output(i,j) = scale_x*( st_coef_x(j,1)*input(i,nx-6)+st_coef_x(j,2)*input(i,nx-5)+st_coef_x(j,3)*input(i,nx-4)+ &
                                 st_coef_x(j,4)*input(i,nx-3)+st_coef_x(j,5)*input(i,nx-2)+st_coef_x(j,6)*input(i,nx-1)+ &
                                 st_coef_x(j,7)*input(i,nx) )
      end do
   end do

   end subroutine zhong_diff_x

!-----------------------------------------------------------------------------------------------------------------------!
! third subroutine: 6th order finite differences (a la Zhong) in y-direction.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine zhong_diff_y(input, output)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, dimension(ny,nx), intent(in) :: input

! output
   real, dimension(ny,nx), intent(out) :: output

! local variables x->j, y->i
   integer :: i, j

!----- Derivative ------------------------------------------------------------------------------------------------------!
   do j=1,nx
      do i=1,4
         output(i,j) = scale_y*( st_coef_y(i,1)*input(1,j)+st_coef_y(i,2)*input(2,j)+st_coef_y(i,3)*input(3,j)+ &
                                 st_coef_y(i,4)*input(4,j)+st_coef_y(i,5)*input(5,j)+st_coef_y(i,6)*input(6,j)+ &
                                 st_coef_y(i,7)*input(7,j) )
      end do
      do i=5,ny-4
         output(i,j) = scale_y*( st_coef_y(i,1)*input(i-3,j)+st_coef_y(i,2)*input(i-2,j)+st_coef_y(i,3)*input(i-1,j)+ &
                                 st_coef_y(i,4)*input(i,j)  +st_coef_y(i,5)*input(i+1,j)+st_coef_y(i,6)*input(i+2,j)+ &
                                 st_coef_y(i,7)*input(i+3,j) )
      end do
      do i=ny-3,ny
         output(i,j) = scale_y*( st_coef_y(i,1)*input(ny-6,j)+st_coef_y(i,2)*input(ny-5,j)+st_coef_y(i,3)*input(ny-4,j)+ &
                                 st_coef_y(i,4)*input(ny-3,j)+st_coef_y(i,5)*input(ny-2,j)+st_coef_y(i,6)*input(ny-1,j)+ &
                                 st_coef_y(i,7)*input(ny,j) )
      end do
   end do

   end subroutine zhong_diff_y
#endif ! SPACE_INT

end module space_integration
