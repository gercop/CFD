subroutine tools(inp1, inp2)

!--- modules
   use params
   use output_ios
   
!-- input
  real, dimension(ny,nx), intent(in) 	:: inp1, inp2 
  integer				:: i,j
  real*8				:: error,error_max
  real, dimension(ny,nx)	 	:: error1
  error_max=0.
  do i = 2,ny-1
    do j = 2,nx-1
      error = abs(inp1(i,j)-inp2(i,j))
      if (error>=error_max) then
        error_max=error
      end if
    end do
  end do  
  write(*,*) 'Maximum error: ',error_max
  
  error1=0.
  do i = 2,ny-1
    do j = 2,nx-1
      error1(i,j) = abs(inp1(i,j)-inp2(i,j))
    end do
  end do  
  call output_init
  call writo(error1,inp1)
  call output_final
  write(*,*) 'output done'
end subroutine tools

