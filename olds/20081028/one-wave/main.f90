program onss

! Main program

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
time    = 0.0
uvel(:) = 0.0

do j=1,nx
   uexact(j)=sin(w*pi*(2.2d0+(1.0d0-pgrid_x(j))))
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


T_INT: do t_step=t_be,t_en

! time step
   time = (t_step-1)*dt

#if (TIME_INT == 1)
      call rk2i_calc(time, uvel)
#elif (TIME_INT == 2)
      call rk3i_calc(time, uvel)
#elif (TIME_INT == 3)
      call rk4i_calc(time, uvel)
#endif ! TIME_INT

end do T_INT 

write(*,*) 'REACHED LAST STEP:',(t_step-1),' TIME:',(t_step-1)*dt

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
