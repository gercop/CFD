subroutine bc_wall(input)

use grid
use sutherland
use parallel_mpi

implicit none

! Wall boundary.

!----- Variables -------------------------------------------------------------------------------------------------------!
! in and output
real, dimension(local_ny,local_nx), intent(inout) :: input

! local variables x->j, y->i
integer :: i, j

!----- Wall ------------------------------------------------------------------------------------------------------------!
! primitive variables
i=1
X1: do j=1,local_nx

   input(i,j) = (1.0-exp((pgrid_x(j)-1.0)*c/dy_visc))/(1.0-exp(-c/dy_visc))

end do X1

end subroutine bc_wall
