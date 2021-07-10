!-----------------------------------------------------------------------------------------------------------------------!
! module which defines input parameters and calculates some other parameters which are used by other subroutines.
!-----------------------------------------------------------------------------------------------------------------------!
module params

   implicit none

   save

!----- Constants -------------------------------------------------------------------------------------------------------!
   real, parameter, public :: a=0.1d0

! grid
   integer, parameter, public :: nx=41
   real, parameter, public :: x_0=0.0d0, x_L=1.0d0
   real, parameter, public :: alpha =0.91 ! 0.91

! time
   integer, parameter, public :: t_be=1, t_en=1000
   real, parameter, public :: dt=2d-4

! final files
   character(LEN=*), parameter, public :: wfname1 = 'error_41.dat', wfname2 = 'solution_41.dat'

! ----- Variables ------------------------------------------------------------------------------------------------------!
   real, public :: pi
   real, public :: n1i3, n2i3, dti2, dti3, dti6

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: Calculates some parameters for later.
!-----------------------------------------------------------------------------------------------------------------------!
   subroutine params_init

   n1i3 = 1.0d0/3.0d0
   n2i3 = 2.0d0/3.0d0

   dti2 = dt   /2.0d0
   dti3 = dt   /3.0d0
   dti6 = dt   /6.0d0

   pi = 4.0d0 * datan(1.0d0)

   end subroutine params_init

end module params
