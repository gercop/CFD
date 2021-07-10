module whole_integration

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
   suvel(1) = 0.0d0
   suvel(nx)= 0.0d0

!----- Second substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(nuvel, dti2, uvel)

! Boundary conditions:
   uvel(1) = 0.0d0
   uvel(nx)= 0.0d0

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
   suvel(1) = 0.0d0
   suvel(nx)= 0.0d0

!----- Second substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(suvel, dti2, suvel)

! Boundary conditions:
   t_intern = time + dt
   suvel(1) = 0.0d0
   suvel(nx)= 0.0d0

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
   suvel(1) = 0.0d0
   suvel(nx)= 0.0d0

!----- Fourth substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(suvel, dti2, uvel)

! Boundary conditions:
   t_intern = time + dt
   uvel(1)  = 0.0d0
   uvel(nx) = 0.0d0

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
   t_intern = time + dti2
   suvel(1) = 0.0d0
   suvel(nx)= 0.0d0

!----- Second substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(uvel, dti2, suvel)
   call timei_calc(nuvel, dti3, nuvel)

! Boundary conditions:
   suvel(1) = 0.0d0
   suvel(nx)= 0.0d0

!----- Third substep ---------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(uvel, dt, suvel)
   call timei_calc(nuvel, dti3, nuvel)

! Boundary conditions:
   t_intern = time + dt
   suvel(1) = 0.0d0
   suvel(nx)= 0.0d0

!----- Fourth substep --------------------------------------------------------------------------------------------------!

! momentum:	    	 
! x-direction
   call rhs_init
   call rhsxmo(suvel)
   call timei_calc(nuvel, dti6, uvel)

! Boundary conditions:
   t_intern = time + dt
   uvel(1) = 0.0d0
   uvel(nx)= 0.0d0

   end subroutine rk4i_calc

end module whole_integration
