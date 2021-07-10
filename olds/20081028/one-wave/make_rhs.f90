!-----------------------------------------------------------------------------------------------------------------------!
! module which contains all routines important for the integration.
!-----------------------------------------------------------------------------------------------------------------------!
module make_rhs

! modules
   use params
   use space_integration

   implicit none

! variables
   real, dimension(nx), private :: rhs_x
   real, dimension(nx), public, save :: deriv_x

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: allocation of stress arrays.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine rhs_init

   deriv_x(:) = 0.0

   end subroutine rhs_init

!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: Calculates rhs for x- and y- derivatives.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine rhsxmo(du_in)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, dimension(nx), intent(in) :: du_in

! local variables x->j
   integer :: j

!----- Make RHS --------------------------------------------------------------------------------------------------------!
   do j=1,nx
      rhs_x(j) = c*du_in(j)
   end do

   call zhong_diff(rhs_x, deriv_x)

   end subroutine rhsxmo

end module make_rhs
