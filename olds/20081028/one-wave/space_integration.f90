module space_integration

! modules
   use params
   use grid

   implicit none

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: 6th order finite differences (a la Zhong).
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine zhong_diff(input, output)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, dimension(nx), intent(in) :: input

! output
   real, dimension(nx), intent(out) :: output

! local variables x->j
   integer :: j

!----- Derivative ------------------------------------------------------------------------------------------------------!
   do j=1,4
      output(j) = scale_x*( st_coef_x(j,1)*input(1)+st_coef_x(j,2)*input(2)+st_coef_x(j,3)*input(3)+ &
                            st_coef_x(j,4)*input(4)+st_coef_x(j,5)*input(5)+st_coef_x(j,6)*input(6)+ &
                            st_coef_x(j,7)*input(7) )
   end do
   do j=5,nx-4
      output(j) = scale_x*( st_coef_x(j,1)*input(j-3)+st_coef_x(j,2)*input(j-2)+st_coef_x(j,3)*input(j-1)+ &
                            st_coef_x(j,4)*input(j)  +st_coef_x(j,5)*input(j+1)+st_coef_x(j,6)*input(j+2)+ &
                            st_coef_x(j,7)*input(j+3) )
   end do
   do j=nx-3,nx
      output(j) = scale_x*( st_coef_x(j,1)*input(nx-6)+st_coef_x(j,2)*input(nx-5)+st_coef_x(j,3)*input(nx-4)+ &
                            st_coef_x(j,4)*input(nx-3)+st_coef_x(j,5)*input(nx-2)+st_coef_x(j,6)*input(nx-1)+ &
                            st_coef_x(j,7)*input(nx) )
   end do

   end subroutine zhong_diff

end module space_integration
