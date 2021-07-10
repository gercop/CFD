program onss

! Main program

! heat eqn pp. 126 and following (Tannehill).

#include 'setup.h'

!----- Modules ---------------------------------------------------------------------------------------------------------!
use whole_integration

implicit none

!----- Local variables -------------------------------------------------------------------------------------------------!
integer, parameter :: filenr = 21
integer :: i, ioerr, j, jfx, jbx, t_step

real :: time, x
real, dimension(nx) :: uvel, uexact

!----- Initialization of parameters ------------------------------------------------------------------------------------!
call params_init
call grid_init

!----- Initialization of parameters ------------------------------------------------------------------------------------!
time = 0.0d0

do j=1,nx
   uvel(j) = 5.0d0*dsin(2.0d0*pi*pgrid_x(j))
end do

!----- Integration -----------------------------------------------------------------------------------------------------!
print *,'TIME INTEGRATION METHOD:'
#if (TIME_INT == 1)
print *,' -> Runge-Kutta second order'
#elif (TIME_INT == 2)
print *,' -> SSP-Runge-Kutta third order'
#elif (TIME_INT == 3)
print *,' -> Runge-Kutta fourth order'
#else ! TIME_INT
print *,'ERROR: Integration method (',TIME_INT,') not known!'
stop
#endif ! TIME_INT

print *,'SPACE DISCRETIZATION:'
#if (SEC_DER == 1)
print *,' -> Twice first derivative'
#elif (SEC_DER == 2)
print *,' -> One second derivative'
#else ! SEC_DER
print *,'ERROR: Discretization method (',SEC_DER,') not known!'
#endif ! SEC_DER

T_INT: do t_step=t_be,t_en

! time step
   time = dble(t_step-1)*dt

#if (TIME_INT == 1)
   call rk2i_calc(time, uvel)
#elif (TIME_INT == 2)
   call rk3i_calc(time, uvel)
#elif (TIME_INT == 3)
   call rk4i_calc(time, uvel)
#endif ! TIME_INT

end do T_INT 

do j=1,nx
   uexact(j)=5.0d0*dsin(2.0d0*pi*pgrid_x(j))*dexp(-a*(2.0d0*pi)**2*dble(t_step-1)*dt)
end do

write(*,*) 'REACHED LAST STEP:',(t_step-1),' TIME:',dble(t_step-1)*dt

!----- Output ----------------------------------------------------------------------------------------------------------!
! write results into result files.
open(filenr, iostat=ioerr, file=wfname1, status='unknown')
if (ioerr.ne.0) then
   print *,'ERROR: Could not open the result file', wfname1
   stop
end if

do j=1,nx
   write(filenr,'(2F20.16)', iostat=ioerr) pgrid_x(j), (uvel(j)-uexact(j))
   if (ioerr.ne.0) then
      print *,'ERROR: Could not write into the result file'
      stop
   end if
end do

close(filenr)

open(filenr, iostat=ioerr, file=wfname2, status='unknown')
if (ioerr.ne.0) then
   print *,'ERROR: Could not open the result file', wfname2
   stop
end if

do j=1,nx
   write(filenr,'(2F20.16)', iostat=ioerr) pgrid_x(j), uexact(j)
   if (ioerr.ne.0) then
      print *,'ERROR: Could not write into the result file'
      stop
   end if
end do

write(filenr,*) ''
do j=1,nx
   write(filenr,'(2F20.16)', iostat=ioerr) pgrid_x(j), uvel(j)
   if (ioerr.ne.0) then
      print *,'ERROR: Could not write into the result file'
      stop
   end if
end do

close(filenr)

! end program      
end program onss
