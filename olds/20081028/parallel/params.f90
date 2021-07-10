!-----------------------------------------------------------------------------------------------------------------------!
! module which defines input parameters and calculates some other parameters which are used by other subroutines.
!-----------------------------------------------------------------------------------------------------------------------!
module params

   implicit none

   save

!----- Constants -------------------------------------------------------------------------------------------------------!
! IC and BC
   real, parameter, public :: dy_visc=2.0, c=3.0, d=3.0

! grid
   integer, parameter, public :: nx=200, ny=200
   real, public		      :: dx, dy
   real, parameter, public :: x_0=0.0, x_L=1.0, y_0=0.0, y_H=1.0
   real, parameter, public :: alpha = 0.8

! time
   integer, parameter, public :: t_be=1, t_en=10000
   real, parameter, public :: dt=1D-8

! final files
!   character(LEN=*), parameter, public :: wfname = 'rk4_solution_41x41.eas'
   character*11, parameter, public :: fbase = 'burger'


! ----- Variables ------------------------------------------------------------------------------------------------------!
   real, public :: pi
   real, public :: dti2, dti3, dti6, idx, idy

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: Calculates some parameters for later.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine params_init

   dx=(x_l-x_0)/(nx-1)
   dy=(y_H-y_0)/(ny-1)
   
   dti2= dt/2.0
   dti3= dt/3.0
   dti6= dt/6.0
   idx = 1.0/dx
   idy = 1.0/dy

   pi = 4.0 * atan(1.0)

   end subroutine params_init

end module params
