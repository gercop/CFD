module interfaces

interface

   subroutine bc_in(input)
   use params
   real, dimension(ny,nx), intent(inout) :: input
   end subroutine bc_in

   subroutine bc_out(input)
   use params
   real, dimension(ny,nx), intent(inout) :: input
   end subroutine bc_out

   subroutine bc_wall(input)
   use params
   real, dimension(ny,nx), intent(inout) :: input
   end subroutine bc_wall

   subroutine bc_free(input)
   use params
   real, dimension(ny,nx), intent(inout) :: input
   end subroutine bc_free

!   subroutine writo(inp1, inp2)
!   use params
!   real, dimension(ny,nx), intent(in) :: inp1, inp2
!   end subroutine writo

end interface

end module interfaces
