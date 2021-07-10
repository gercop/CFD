!-----------------------------------------------------------------------------------------------------------------------!
! module which contains all routines important for the integration.
!-----------------------------------------------------------------------------------------------------------------------!
module make_rhs

#include "setup.h"
! include "setup.h"

! modules
   use interfaces
   use params
   use space_integration
   use stresses

   implicit none

! variables
   real, dimension(:,:),allocatable, public :: rhs_x, rhs_y
   real, dimension(:,:),allocatable, public, save :: deriv_x, deriv_y

contains

  subroutine rhs_alloc
   allocate(rhs_x(local_ny,local_nx),rhs_y(local_ny,local_nx))
   allocate(deriv_x(local_ny,local_nx),deriv_y(local_ny,local_nx))
  end subroutine rhs_alloc


!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: allocation of stress arrays.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine rhs_init

   
   deriv_x(:,:) = 0.0
   deriv_y(:,:) = 0.0

   end subroutine rhs_init

!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: Calculates rhs for x- and y- derivatives (x-momentum).
!-----------------------------------------------------------------------------------------------------------------------!
#if (SPACE_INT == 1)
   subroutine rhsxmo(du_in,ify, iby, jfx, jbx )
#else
   subroutine rhsxmo(du_in)
#endif

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
#if (SPACE_INT == 1)
   integer, intent(in) :: ify, iby, jfx, jbx
#endif ! SPACE_INT
   real, dimension(local_ny,local_nx), intent(in) :: du_in

! local variables x->j, y->i
   integer :: i, j

!----- Make RHS --------------------------------------------------------------------------------------------------------!
   Y: do i=1,local_ny
      X: do j=1,local_nx

         rhs_x(i,j) = c*du_in(i,j)-txx(i,j)
         rhs_y(i,j) = d*du_in(i,j)-tyy(i,j)

      end do X
   end do Y

! x derivative
#if (SPACE_INT == 1)
   call bfc_diff(rhs_x, idx, 0, 0, jfx, jbx, deriv_x)
#elif (SPACE_INT == 2)
   call zhong_diff_x(rhs_x, deriv_x)
#endif ! SPACE_INT

! y derivative
#if (SPACE_INT == 1)
   call bfc_diff(rhs_y, idy, ify, iby, 0, 0, deriv_y)
#elif (SPACE_INT == 2)
   call zhong_diff_y(rhs_y, deriv_y)
#endif ! SPACE_INT

   end subroutine rhsxmo

end module make_rhs
