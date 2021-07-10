!-----------------------------------------------------------------------------!
!module which handles the output for the ios format			      !
!-----------------------------------------------------------------------------!
module output_ios

!--- modules
  use params

  implicit none

!----- Variables -------------------------------------------------------------!

!--- variables for ios file
  integer 				:: m1,m2,m3,mt,mp,minf
  integer 				:: itape,imach,iunit
  integer,dimension(512)		:: itimes
!  character*11 				:: fbase
  character*72,dimension(100)		:: inf

!--- additional variables
  real,dimension(ny,nx,1)		:: array1,array2

contains

subroutine output_init
!--- initialize parameters for output
!  fbase = 'burger'
  itape = 25
  iunit = 25
  imach = 2
  m1=ny
  m2=nx
  m3=1
  mt    = 1
  mp    = 2
  minf  = 1
  inf(1)= 'u-numerical'
  inf(2)= 'u-exact'
  inf(3)= 'Burgers equation'
  itimes(1)= 1
!--- open file for output  
  call writecd(itape,iunit,fbase,imach,mt,m3,m2,m1,mp,itimes,inf,minf)
end subroutine output_init


subroutine writo(inp1, inp2)
!-- input
  real, dimension(ny,nx), intent(in) 	:: inp1, inp2
  integer				:: i,j

!--- convert inp1,inp2 to ios arrays
  do i = 1,ny
    do j = 1,nx
      array1(i,j,1)=inp1(i,j)
      array2(i,j,1)=inp2(i,j)
    end do
  end do  
call writed(itape,array1)
call writed(itape,array2)
end subroutine writo

subroutine output_final
  close(itape)
end subroutine output_final

end module output_ios
