!-----------------------------------------------------------------------------------------------------------------------!
! module which calculates the time integration.
!-----------------------------------------------------------------------------------------------------------------------!
module time_integration

#include 'setup.h'

! modules
   use make_rhs

   implicit none

! variables
   private :: deriv_x

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: integration.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine timei_calc(input, vdt, output)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, intent(in) :: vdt
   real, dimension(nx), intent(in)  :: input

! output
   real, dimension(nx), intent(out) :: output

! local variables x->j
   integer :: j

!----- Integration -----------------------------------------------------------------------------------------------------!
#if (INFLOW_BC == 1) || (INFLOW_BC == 3)
   do j=2,nx
      output(j) = input(j)-vdt*deriv_x(j)
   end do
#else ! INFLOW_BC
   do j=1,nx
      output(j) = input(j)-vdt*deriv_x(j)
   end do
#endif ! INFLOW_BC

   end subroutine timei_calc

end module time_integration
