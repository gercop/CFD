subroutine bc_out(input)

use params	
use parallel_mpi

implicit none

! Outflow boundary.

!----- Variables -------------------------------------------------------------------------------------------------------!
! in and output
real, dimension(local_ny,local_nx), intent(inout) :: input

! local variables x->j, y->i
integer :: i, j

!----- Outflow ---------------------------------------------------------------------------------------------------------!
! primitive variables
j=local_nx
Y1: do i=1,local_ny

   input(i,j) = 0.0

end do Y1

end subroutine bc_out
