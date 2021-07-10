!-----------------------------------------------------------------------------------------------------------------------!
! module which calculates the time integration.
!-----------------------------------------------------------------------------------------------------------------------!
module time_integration

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
   do j=2,nx-1
      output(j) = input(j)+vdt*deriv_x(j)
   end do

   end subroutine timei_calc

end module time_integration
