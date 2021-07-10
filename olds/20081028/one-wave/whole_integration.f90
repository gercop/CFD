module whole_integration

#include 'setup.h'

! modules
   use time_integration

   implicit none

! variables
   real, private :: t_intern
   real, dimension(nx), private :: suvel, nuvel

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: rk2 [Heun] time integration.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine rk2i_calc(time, uvel)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, intent(in) :: time
   real, dimension(nx), intent(inout) :: uvel

!----- Integration -----------------------------------------------------------------------------------------------------!
   t_intern = time

!----- First substep ---------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(uvel)
   call timei_calc(uvel, dt, suvel)
   call timei_calc(uvel, dti2, nuvel)

! Boundary conditions:
   t_intern = time + dt
   suvel(1) = dsin(w*pi*t_intern)

!----- Second substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(nuvel, dti2, uvel)

! Boundary conditions:
   uvel(1) = dsin(w*pi*t_intern)

   end subroutine rk2i_calc

!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: SSP-rk3 with 4 stages time integration.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine rk3i_calc(time, uvel)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, intent(in) :: time
   real, dimension(nx), intent(inout) :: uvel

! local
   integer :: i

!----- Integration -----------------------------------------------------------------------------------------------------!
   t_intern = time

!----- First substep ---------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(uvel)
   call timei_calc(uvel, dti2, suvel)

! Boundary conditions:
   t_intern = time + dti2
   suvel(1) = dsin(w*pi*t_intern)

!----- Second substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(suvel, dti2, suvel)

! Boundary conditions:
   t_intern = time + dt
   suvel(1) = dsin(w*pi*t_intern)

!----- Third substep ---------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   do i=1,nx
      suvel(i) = n2i3*uvel(i)+n1i3*suvel(i)
   end do
   call timei_calc(suvel, dti6, suvel)

! Boundary conditions:
   t_intern = time + dti2
   suvel(1) = dsin(w*pi*t_intern)

!----- Fourth substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(suvel, dti2, uvel)

! Boundary conditions:
   t_intern = time + dt
   uvel(1) = dsin(w*pi*t_intern)

   end subroutine rk3i_calc

!-----------------------------------------------------------------------------------------------------------------------!
! third subroutine: standard rk4 time integration.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine rk4i_calc(time, uvel)

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
   real, intent(in) :: time
   real, dimension(nx), intent(inout) :: uvel

!----- Integration -----------------------------------------------------------------------------------------------------!
   t_intern = time

!----- First substep ---------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(uvel)
   call timei_calc(uvel, dti2, suvel)
   call timei_calc(uvel, dti6, nuvel)

! Boundary conditions:
#if (INFLOW_BC == 1)
   t_intern = time + dti2
   suvel(1) = dsin(w*pi*t_intern)
#elif (INFLOW_BC == 3)
   suvel(1) = dsin(w*pi*t_intern)+dti2*w*pi*dcos(w*pi*t_intern)
#endif ! INFLOW_BC

!----- Second substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(uvel, dti2, suvel)
   call timei_calc(nuvel, dti3, nuvel)

! Boundary conditions:
#if (INFLOW_BC == 1)
   suvel(1) = dsin(w*pi*t_intern)
#elif (INFLOW_BC == 3)
   suvel(1) = dsin(w*pi*t_intern)+dti2*w*pi*dcos(w*pi*t_intern)-(dti2*w*pi)**2*dsin(w*pi*t_intern)
#endif ! INFLOW_BC

!----- Third substep ---------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(uvel, dt, suvel)
   call timei_calc(nuvel, dti3, nuvel)

! Boundary conditions:
#if (INFLOW_BC == 1)
   t_intern = time + dt
   suvel(1) = dsin(w*pi*t_intern)
#elif (INFLOW_BC == 3)
   suvel(1) = dsin(w*pi*t_intern)+dt*w*pi*dcos(w*pi*t_intern)-(dt*w*pi)**2/2.0d0*dsin(w*pi*t_intern) &
              -(dt*w*pi)**3/4.0d0*dcos(w*pi*t_intern)
#endif ! INFLOW_BC

!----- Fourth substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(nuvel, dti6, uvel)

! Boundary conditions:
   t_intern = time + dt
   uvel(1) = dsin(w*pi*t_intern)

   end subroutine rk4i_calc

end module whole_integration
