!-----------------------------------------------------------------------------------------------------------------------!
! module which calculates the time integration.
!-----------------------------------------------------------------------------------------------------------------------!
module time_integration

! modules
   use make_rhs
   use parallel_mpi

   implicit none

   private :: deriv_x, deriv_y

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: integration.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine timei_calc(input, vdt, output)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, intent(in) :: vdt
   real, dimension(local_ny,local_nx), intent(in)  :: input

! output
   real, dimension(local_ny,local_nx), intent(out) :: output

! local variables x->j, y->i
   integer :: i, j

!----- Integration -----------------------------------------------------------------------------------------------------!
   Y: do i=2,(local_ny-1)
      X: do j=2,(local_nx-1)
         output(i,j) = input(i,j)-vdt*( deriv_x(i,j) + deriv_y(i,j) )
      end do X
   end do Y

   end subroutine timei_calc

end module time_integration
