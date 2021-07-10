module space_integration

! modules
   use params
   use grid

   implicit none

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: 6th order finite differences (a la Zhong) for a first derivative.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine zhong_1st_diff(input, output)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, dimension(nx), intent(in) :: input

! output
   real, dimension(nx), intent(out) :: output

! local variables
   integer :: j

!----- Derivative ------------------------------------------------------------------------------------------------------!
   do j=1,4
      output(j) = x_scale_1st*( x_coef_1st(j,1)*input(1)+x_coef_1st(j,2)*input(2)+x_coef_1st(j,3)*input(3)+ &
                                x_coef_1st(j,4)*input(4)+x_coef_1st(j,5)*input(5)+x_coef_1st(j,6)*input(6)+ &
                                x_coef_1st(j,7)*input(7) )
   end do
   do j=5,nx-4
      output(j) = x_scale_1st*( x_coef_1st(j,1)*input(j-3)+x_coef_1st(j,2)*input(j-2)+x_coef_1st(j,3)*input(j-1)+ &
                                x_coef_1st(j,4)*input(j)  +x_coef_1st(j,5)*input(j+1)+x_coef_1st(j,6)*input(j+2)+ &
                                x_coef_1st(j,7)*input(j+3) )
   end do
   do j=nx-3,nx
      output(j) = x_scale_1st*( x_coef_1st(j,1)*input(nx-6)+x_coef_1st(j,2)*input(nx-5)+x_coef_1st(j,3)*input(nx-4)+ &
                                x_coef_1st(j,4)*input(nx-3)+x_coef_1st(j,5)*input(nx-2)+x_coef_1st(j,6)*input(nx-1)+ &
                                x_coef_1st(j,7)*input(nx) )
   end do

   end subroutine zhong_1st_diff

!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: 6th order finite differences (a la Zhong) for a second derivative.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine zhong_2nd_diff(input, output)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, dimension(nx), intent(in) :: input

! output
   real, dimension(nx), intent(out) :: output

! local variables
   integer :: j

!----- Derivative ------------------------------------------------------------------------------------------------------!
   do j=1,4
      output(j) = x_scale_2nd*( x_coef_2nd(j,1)*input(1)+x_coef_2nd(j,2)*input(2)+x_coef_2nd(j,3)*input(3)+ &
                                x_coef_2nd(j,4)*input(4)+x_coef_2nd(j,5)*input(5)+x_coef_2nd(j,6)*input(6)+ &
                                x_coef_2nd(j,7)*input(7) )
   end do
   do j=5,nx-4
      output(j) = x_scale_2nd*( x_coef_2nd(j,1)*input(j-3)+x_coef_2nd(j,2)*input(j-2)+x_coef_2nd(j,3)*input(j-1)+ &
                                x_coef_2nd(j,4)*input(j)  +x_coef_2nd(j,5)*input(j+1)+x_coef_2nd(j,6)*input(j+2)+ &
                                x_coef_2nd(j,7)*input(j+3) )
   end do
   do j=nx-3,nx
      output(j) = x_scale_2nd*( x_coef_2nd(j,1)*input(nx-6)+x_coef_2nd(j,2)*input(nx-5)+x_coef_2nd(j,3)*input(nx-4)+ &
                                x_coef_2nd(j,4)*input(nx-3)+x_coef_2nd(j,5)*input(nx-2)+x_coef_2nd(j,6)*input(nx-1)+ &
                                x_coef_2nd(j,7)*input(nx) )
   end do

   end subroutine zhong_2nd_diff

end module space_integration
