subroutine bc_in(input)

use grid
use sutherland

implicit none

! Inflow boundary.

!----- Variables -------------------------------------------------------------------------------------------------------!
! in and output
real, dimension(local_ny,local_nx), intent(inout) :: input

! local variables x->j, y->i
integer :: i, j

!----- Inflow ----------------------------------------------------------------------------------------------------------!
! primitive variables
j=1
Y1: do i=2,local_ny-1

   input(i,j) = (1.0-exp((pgrid_y(i)-1.0)*d/dy_visc))/(1.0-exp(-d/dy_visc))

end do Y1 

end subroutine bc_in
