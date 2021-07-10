module whole_integration

#include "setup.h"
! include "setup.h"

! modules
   use parallel_mpi
   use time_integration
   use output_ios
   
   implicit none

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: rk2 [Heun] time integration.
!-----------------------------------------------------------------------------------------------------------------------!
#if (SPACE_INT ==1) 
   subroutine rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#else
   subroutine rk2i_calc(local_uvel)
#endif   
 integer::i,j
!----- Variables -------------------------------------------------------------------------------------------------------!
! input
#if (SPACE_INT == 1)
   integer, intent(in) :: ify, iby, jfx, jbx
#endif ! SPACE_INT
   real, dimension(ny,nx) :: uvel
   real,dimension(local_ny,local_nx),intent(inout) ::local_uvel
   
! local variables
#if (SPACE_INT == 1)
   integer :: lify, liby, ljfx, ljbx
#endif ! SPACE_INT
   real, dimension(local_ny,local_nx) :: suvel, nuvel

!----- Integration -----------------------------------------------------------------------------------------------------!
! switch integration direction for Pre inner derivatives and Cor outer derivatives.
#if (SPACE_INT == 1)
   lify=iby; liby=ify; ljfx=jbx; ljbx=jfx
#endif ! SPACE_INT

! ------------------------------ Predictor -----------------------------      
! Boundary conditions:
! inflow
if (mod(my_rank,p_yx(2))==0) then   
   call bc_in(local_uvel)
end if

! outflow
if (mod((my_rank+1),p_yx(2))==0) then
   call bc_out(local_uvel)
end if

! wall
if (my_rank.lt.p_yx(2)) then
   call bc_wall(local_uvel)
end if

! free-stream
if (my_rank.ge.(p-p_yx(2)))then
   call bc_free(local_uvel)
end if

! Sutherland + stresses:	    
   call suther_calc
   call stress_init
#if (SPACE_INT == 1)
   call stress_calc(local_uvel, lify, liby, ljfx, ljbx )
#else
   call stress_calc(local_uvel)
#endif

! momentum:	    	 
! x-direction
   call rhs_init
#if (SPACE_INT == 1)
   call rhsxmo(local_uvel, ify, iby, jfx, jbx)
#else
   call rhsxmo(local_uvel)
#endif

   call timei_calc(local_uvel, dt, suvel)
   call timei_calc(local_uvel, dti2, nuvel)

!---- excange the ghostpoints between the different processors
    call SYSTEM_CLOCK(COUNT=clock_s) ! Start timing

   if (p.ne.1) then
     call exchange_ghostpoints(suvel)
   end if 
!   call exchange_ghostpoints(nuvel)
    call SYSTEM_CLOCK(COUNT=clock_e) ! Stop timing
    ! Calculate the elapsed time in seconds:
    commu_time=commu_time+dble(clock_e-clock_s)/dble(clock_rate)


! ------------------------------ Corrector -----------------------------
! Boundary conditions:
! inflow
  if (mod(my_rank,p_yx(2))==0) then   
    call bc_in(suvel)
  end if
! outflow
  if (mod((my_rank+1),p_yx(2))==0) then
    call bc_out(suvel)
  end if
! wall
  if (my_rank.lt.p_yx(2)) then
    call bc_wall(suvel)
  end if
! free-stream
  if (my_rank.ge.(p-p_yx(2)))then
    call bc_free(suvel)
  end if 





! Sutherland + stresses + heat fluxes:	    
   call suther_calc
   call stress_init
#if (SPACE_INT == 1)
   call stress_calc(suvel, ify, iby, jfx, jbx)
#else
   call stress_calc(suvel)
#endif


! momentum:	    	 
! x-direction
   call rhs_init
#if (SPACE_INT == 1)
   call rhsxmo(suvel, lify, liby, ljfx, ljbx)
#else
   call rhsxmo(suvel)
#endif

   call timei_calc(nuvel, dti2, local_uvel)

!----excange the ghostpoints between the different processors
    call SYSTEM_CLOCK(COUNT=clock_s) ! Start timing

  if (p.ne.1) then
    call exchange_ghostpoints(local_uvel)
  end if
    call SYSTEM_CLOCK(COUNT=clock_e) ! Stop timing
    ! Calculate the elapsed time in seconds:
    commu_time=commu_time+dble(clock_e-clock_s)/dble(clock_rate)
    

!  call group_for_output(local_uvel,uvel)
!  call group_for_output(local_uexact,uexact)
!  if (my_rank==0) then
!    call output_init
!    call writo(uvel, uexact)
!    call output_final
!   stop
!  end if



   end subroutine rk2i_calc






!-----------------------------------------------------------------------------------------------------------------------!
! second subroutine: standard rk4 time integration.
!-----------------------------------------------------------------------------------------------------------------------!
#if (SPACE_INT == 1)
   subroutine rk4i_calc(uvel, ify, iby, jfx, jbx)
#else
   subroutine rk4i_calc(uvel)
#endif

!----- Variables -------------------------------------------------------------------------------------------------------!
! input
#if (SPACE_INT == 1)
   integer, intent(in) :: ify, iby, jfx, jbx
#endif ! SPACE_INT
   real, dimension(ny,nx), intent(inout) :: uvel

! local variables
#if (SPACE_INT == 1)
   integer :: lify, liby, ljfx, ljbx
#endif ! SPACE_INT
   real, dimension(ny,nx) :: suvel, nuvel

!----- Integration -----------------------------------------------------------------------------------------------------!
! switch integration direction for Pre inner derivatives and Cor outer derivatives.
#if (SPACE_INT == 1)
   lify=iby; liby=ify; ljfx=jbx; ljbx=jfx
#endif ! SPACE_INT

!----- First substep ---------------------------------------------------------------------------------------------------!
! Boundary conditions:
! inflow
   call bc_in(uvel)

! outflow
   call bc_out(uvel)

! wall
   call bc_wall(uvel)

! free-stream
   call bc_free(uvel)

! Sutherland + stresses:	    
   call suther_calc
   call stress_init
#if (SPACE_INT == 1)
   call stress_calc(uvel, lify, liby, ljfx, ljbx)
#else
   call stress_calc(uvel)
#endif


! momentum:	    	 
! x-direction
   call rhs_init
#if (SPACE_INT == 1)
   call rhsxmo(uvel, ify, iby, jfx, jbx)
#else
   call rhsxmo(uvel)
#endif

   call timei_calc(uvel, dti2, suvel)
   call timei_calc(uvel, dti6, nuvel)

!----- Second substep --------------------------------------------------------------------------------------------------!
! Boundary conditions:
! inflow
   call bc_in(suvel)

! outflow
   call bc_out(suvel)

! wall
   call bc_wall(suvel)

! free-stream
   call bc_free(suvel)

! Sutherland + stresses + heat fluxes:	    
   call suther_calc
   call stress_init
#if (SPACE_INT == 1)
   call stress_calc(suvel, ify, iby, jfx, jbx)
#else
   call stress_calc(suvel)
#endif


! momentum:	    	 
! x-direction
   call rhs_init
#if (SPACE_INT == 1)
   call rhsxmo(suvel, lify, liby, ljfx, ljbx)
#else
   call rhsxmo(suvel)
#endif

   call timei_calc(uvel, dti2, suvel)
   call timei_calc(nuvel, dti3, nuvel)

!----- Third substep ---------------------------------------------------------------------------------------------------!
! Boundary conditions:
! inflow
   call bc_in(suvel)

! outflow
   call bc_out(suvel)

! wall
   call bc_wall(suvel)

! free-stream
   call bc_free(suvel)

! Sutherland + stresses + heat fluxes:	    
   call suther_calc
   call stress_init
#if (SPACE_INT == 1)
   call stress_calc(suvel, lify, liby, ljfx, ljbx)
#else
   call stress_calc(suvel)
#endif


! momentum:	    	 
! x-direction
   call rhs_init
#if (SPACE_INT == 1)
   call rhsxmo(suvel, ify, iby, jfx, jbx)
#else
   call rhsxmo(suvel)
#endif

   call timei_calc(uvel, dt, suvel)
   call timei_calc(nuvel, dti3, nuvel)

!----- Fourth substep --------------------------------------------------------------------------------------------------!
! Boundary conditions:
! inflow
   call bc_in(suvel)

! outflow
   call bc_out(suvel)

! wall
   call bc_wall(suvel)

! free-stream
   call bc_free(suvel)

! Sutherland + stresses + heat fluxes:	    
   call suther_calc
   call stress_init
#if (SPACE_INT ==1) 
   call stress_calc(suvel, ify, iby, jfx, jbx)
#else
   call stress_calc(suvel)
#endif


! momentum:	    	 
! x-direction
   call rhs_init
#if (SPACE_INT == 1)
   call rhsxmo(suvel, lify, liby, ljfx, ljbx)
#else
   call rhsxmo(suvel)
#endif

   call timei_calc(nuvel, dti6, uvel)

   end subroutine rk4i_calc

end module whole_integration
