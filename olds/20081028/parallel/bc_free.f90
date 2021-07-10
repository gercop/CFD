subroutine bc_free(input)

use params
use parallel_mpi
implicit none

! Free-stream boundary.

!----- Variables -------------------------------------------------------------------------------------------------------!
! in and output
real, dimension(local_ny,local_nx), intent(inout) :: input

! local variables x->j, y->i
integer :: i, j

!----- Free-stream -----------------------------------------------------------------------------------------------------!
i=local_ny
X1: do j=1,local_nx

   input(i,j) = 0.0

end do X1

end subroutine bc_free
