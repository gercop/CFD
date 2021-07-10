!-----------------------------------------------------------------------------------------------------------------------!
! module which calculates sutherland's law and makes the viscosity globally available.
!-----------------------------------------------------------------------------------------------------------------------!
module sutherland

! modules
   use params
   use parallel_mpi

   implicit none

! variables
   real, dimension(:,:),allocatable, public, save  :: visc

contains

!-----------------------------------------------------------------------------------------------------------------------!
! first subroutine: Calculation of viscosity.
!-----------------------------------------------------------------------------------------------------------------------!
  subroutine suther_init
    allocate(visc(local_ny,local_nx))
  end subroutine
   
   subroutine suther_calc
   

! local variables x->j, y->i
   integer :: i, j
 
!----- Sutherland ------------------------------------------------------------------------------------------------------!
   visc(:,:) = dy_visc

   end subroutine suther_calc

end module sutherland
