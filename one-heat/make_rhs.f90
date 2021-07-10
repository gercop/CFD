!-----------------------------------------------------------------------------------------------------------------------!
! module which contains all routines important for the integration.
!-----------------------------------------------------------------------------------------------------------------------!
module make_rhs

#include 'setup.h'

! modules
   use params
   use space_integration

   implicit none

! variables
   real, dimension(nx), private :: rhs_x
   real, dimension(nx), public, save :: deriv_x

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: Initialization.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine rhs_init

   deriv_x(:) = 0.0

   end subroutine rhs_init

!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: Calculates rhs.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine rhsxmo(du_in)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, dimension(nx), intent(in) :: du_in

! local variables
   integer :: j

!----- Make RHS --------------------------------------------------------------------------------------------------------!
#if (SEC_DER == 1)
   call zhong_1st_diff(du_in, rhs_x)
#endif

   do j=1,nx
      rhs_x(j) = a*rhs_x(j)
   end do

#if (SEC_DER == 1)
   call zhong_1st_diff(rhs_x, deriv_x)
#else
   call zhong_2nd_diff(rhs_x, deriv_x)
#endif

   end subroutine rhsxmo

end module make_rhs
