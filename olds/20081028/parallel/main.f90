program onss

! Main program
!
! Uses MacCormack method in Tannehill pg.230-231

!----- Modules ---------------------------------------------------------------------------------------------------------!
use whole_integration
use parallel_mpi
use stresses
use output_ios
 
#include 'setup.h'

!---superandi
!#define TIME_INT 1
!#define SPACE_INT 1

!implicit none

!----- Local variables -------------------------------------------------------------------------------------------------!
integer :: i, ify, iby, j, jfx, jbx, t_step

real :: time, x
!real, dimension(ny,nx) :: uvel, uexact
!real,dimension(:,:),allocatable :: local_uvel, local_uexact
real, dimension(ny,nx) :: uvel
real,dimension(:,:),allocatable :: local_uvel
 
!---initialize MPI parallelization--------------------------------------------------------------------------------------!
call initmpi
call time_start

!----- Initialization of parameters ------------------------------------------------------------------------------------!
call params_init

!----- Domain decomposition in x-direction only (so far)----------------------------------------------------------------!
call decompose_domain
allocate(local_uvel(local_ny,local_nx),local_uexact(local_ny,local_nx))

!----- Initialization of grid ------------------------------------------------------------------------------------------!
call grid_init  
  
!---allocate arrays for several subroutines---------------------------------------------------------------------------!
call suther_init
call stress_alloc
call rhs_alloc

!----- Initialization of parameters ------------------------------------------------------------------------------------!
time       = 0.0
local_uvel (:,:) = 0.0
do i=1,local_ny
  do j=1,local_nx
    local_uexact(i,j)=(1.0-exp((pgrid_x(j)-1.0)*c/dy_visc))/ &
                      (1.0-exp(-c/dy_visc))* &
                      (1.0-exp((pgrid_y(i)-1.0)*d/dy_visc))/ &
                      (1.0-exp(-d/dy_visc))
  end do
end do
local_uvel (:,:) = local_uexact(:,:)
local_uvel (:,:) = 0.0

!----- Integration -----------------------------------------------------------------------------------------------------!
if (my_rank==0) then
print *,'TIME INTEGRATION METHOD:'
#if (TIME_INT == 1)
print *,' -> Runge-Kutta second order'
#elif (TIME_INT == 2)
print *,' -> Runge-Kutta fourth order'
#else ! TIME_INT
print *,'ERROR: Integration method (',TIME_INT,') not known!'
stop
#endif ! TIME_INT

print *,'SPACE INTEGRATION METHOD:'
#if (SPACE_INT == 1)
print *,' -> Second order finite split differences'
#elif (SPACE_INT == 2)
print *,' -> Sixth order finite differences (a la Zhong)'
#else
print *,'ERROR: Integration method (',SPACE_INT,') not know!'
stop
#endif
print *,''

end if

T_INT: do t_step=t_be,t_en

! time step
   time = t_step*dt

#if (SPACE_INT == 1)

   select case (mod(t_step,8))

!----- CASE 1: OUTER Pre and INNER Cor: F, F ---------------------------------------------------------------------------!
   case (1)
      jfx = 1
      jbx = 0
      ify = 1
      iby = 0
#if (TIME_INT == 1)
      call rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel, ify, iby, jfx, jbx)
#endif ! TIME_INT

!----- CASE 2: OUTER Pre and INNER Cor: B, B ---------------------------------------------------------------------------!
   case (2)
      jfx = 0
      jbx = 1
      ify = 0
      iby = 1
#if (TIME_INT == 1)
      call rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel, ify, iby, jfx, jbx)
#endif ! TIME_INT

!----- CASE 3: OUTER Pre and INNER Cor: F, F ---------------------------------------------------------------------------!
   case (3)
      jfx = 1
      jbx = 0
      ify = 1
      iby = 0
#if defined(RK2)
      call rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel, ify, iby, jfx, jbx)
#endif ! TIME_INT

!----- CASE 4: OUTER Pre and INNER Cor: B, F ---------------------------------------------------------------------------!
   case (4)
      jfx = 0
      jbx = 1
      ify = 0
      iby = 1
#if (TIME_INT == 1)
      call rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel, ify, iby, jfx, jbx)
#endif ! TIME_INT

!----- CASE 5: OUTER Pre and INNER Cor: F, B ---------------------------------------------------------------------------!
   case (5)
      jfx = 0
      jbx = 1
      ify = 1
      iby = 0
#if (TIME_INT == 1)
      call rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel, ify, iby, jfx, jbx)
#endif ! TIME_INT

!----- CASE 6: OUTER Pre and INNER Cor: B, F ---------------------------------------------------------------------------!
   case (6)
      jfx = 1
      jbx = 0
      ify = 0
      iby = 1
#if (TIME_INT == 1)
      call rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel, ify, iby, jfx, jbx)
#endif ! TIME_INT

!----- CASE 7: OUTER Pre and INNER Cor: F, B ---------------------------------------------------------------------------!
   case (7)
      jfx = 1
      jbx = 0
      ify = 0
      iby = 1
#if (TIME_INT == 1)
      call rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel, ify, iby, jfx, jbx)
#endif ! TIME_INT

!----- CASE 8: OUTER Pre and INNER Cor: B, B ---------------------------------------------------------------------------!
   case (0)
      jfx = 0
      jbx = 1
      ify = 1
      iby = 0
#if (TIME_INT == 1)
      call rk2i_calc(local_uvel, ify, iby, jfx, jbx)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel, ify, iby, jfx, jbx)
#endif ! TIME_INT

   end select

#elif (SPACE_INT == 2)

#if (TIME_INT == 1)
      call rk2i_calc(local_uvel)
#elif (TIME_INT == 2)
      call rk4i_calc(uvel)
#endif ! TIME_INT

#endif ! SPACE_INT

!if (my_rank==0) then
!  write(*,*) 'timestep ',t_step,' done...'
!end if

end do T_INT 

!----- Output ----------------------------------------------------------------------------------------------------------!
! write data


!---stop time for speed-up testing and finalize MPI parallelization
call time_end

if (p==1) then
  uvel(:,:) = local_uvel(:,:)
  uexact (:,:) = local_uexact(:,:)
else
  call group_for_output(local_uvel,uvel)
  call group_for_output(local_uexact,uexact)
!  if (my_rank==0) then
!    call output_init
!    call writo(uvel, uexact)
!    call output_final
!  end if
end if
!---postprocessing tools for error check
if (my_rank==0) then
  call tools(uvel,uexact)
end if
write(*,*) my_rank,'commu_time',commu_time  

call MPI_FINALIZE(ierror)

! end program      
end program onss
