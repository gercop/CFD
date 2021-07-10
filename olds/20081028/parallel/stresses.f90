!-----------------------------------------------------------------------------------------------------------------------!
! module which calculates the stresses and makes them globally available.
!-----------------------------------------------------------------------------------------------------------------------!
module stresses

#include "setup.h"

! modules
   use interfaces
   use space_integration
   use sutherland
   use parallel_mpi

   implicit none

! variables
   real, dimension(:,:),allocatable, public,save :: deriv
   real, dimension(:,:),allocatable, public, save :: txx, tyy

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: initialization
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine stress_alloc
     allocate(deriv(local_ny,local_nx),txx(local_ny,local_nx),&
              tyy(local_ny,local_nx))
   end subroutine stress_alloc
   
   subroutine stress_init


!----- Initialization --------------------------------------------------------------------------------------------------!
   txx(:,:)  =0.0
   tyy(:,:)  =0.0

   end subroutine stress_init

!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: Calculates stresses.
!-----------------------------------------------------------------------------------------------------------------------!
#if (SPACE_INT ==1)
   subroutine stress_calc(input,ify, iby, jfx, jbx) 
#else
   subroutine stress_calc(input)
#endif   

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
#if (SPACE_INT == 1)
   integer, intent(in) :: ify, iby, jfx, jbx
#endif ! SPACE_INT
   real, dimension(local_ny,local_nx), intent(in) :: input

! local variables x->j, y->i
   integer :: i, j

!----- Calculate stresses ----------------------------------------------------------------------------------------------!      
! x derivative
#if (SPACE_INT == 1)
   call bfc_diff(input, idx, 0, 0, jfx, jbx, deriv)
#elif (SPACE_INT == 2)
   call zhong_diff_x(input, deriv)
#endif ! SPACE_INT

   Y1: do i=1,local_ny
      X1: do j=1+jbx,local_nx-jfx
	 txx(i,j)   = visc(i,j)*deriv(i,j)
      end do X1
   end do Y1

! y derivative
#if (SPACE_INT == 1)
   call bfc_diff(input, idy, ify, iby, 0, 0, deriv)
#elif (SPACE_INT == 2)
   call zhong_diff_y(input, deriv)
#endif ! SPACE_INT

   Y2: do i=1+iby,local_ny-ify
      X2: do j=1,local_nx
         tyy(i,j)   = visc(i,j)*deriv(i,j)
      end do X2
   end do Y2

   end subroutine stress_calc

end module stresses
