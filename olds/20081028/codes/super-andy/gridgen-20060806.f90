program write_grid

!***********************************************************!
!							    !
!This program writes a grid_file for x3dp.		    !
!							    !
!							    !
!***********************************************************!


implicit none

integer 				:: m1,m2,m3,mt,mp,minf
integer 				:: itape,imach,iunit
integer,dimension(512)			:: itimes
integer					:: i,k,j
      
real 					:: deltax,deltay,deltaz,dummi
real*8, dimension(:,:,:),allocatable	:: array1,array2,array3		
      
character*11 				:: fbase
character*72,dimension(100)		:: inf
      
	

!-------------------------------------------------------------!
!------------------------INPUT--------------------------------!
!-------------------------------------------------------------!


write (*,*) 'Input deltax'
read (*,*) deltax

write (*,*) 'Input deltay'
read (*,*) deltay

write (*,*) 'Input deltaz'
read (*,*) deltaz

write(*,*) 'Input points in x-direction'
read (*,*) m1

write(*,*) 'Input points in y-direction'
read (*,*) m2

write(*,*) 'Input points in z-direction'
read (*,*) m3

allocate(array1(m1,m2,m3),array2(m1,m2,m3),array3(m1,m2,m3))

 
do j=1,m2
 	do k=1,m3
 		do i=1, m1
 			array1(i,j,k) = i*deltax/m1
			array2(i,j,k) = j*deltay/m2
 			array3(i,j,k) = k*deltaz/m3
 		end do
 	end do
end do
 
!-------------------------------------------------------------!
!-------------------OUTPUT------------------------------------!
!-------------------------------------------------------------!  
  
 
fbase = 'grid'
itape = 10
iunit = 10
imach = 1
mt    = 1
mp    = 3
minf  = 1
inf(1)= 'x-grid'
inf(2)= 'y-grid'
inf(3)= 'z-grid'
inf(4)= 'grid file for x3dp'
itimes(1)= 1
     
call writecd(itape,iunit,fbase,imach,mt,m3,m2,m1,mp,itimes,inf,minf)

call writed(itape,array1)
call writed(itape,array2)
call writed(itape,array3)

close(itape)
	
write (*,*) 'Grid generation sucessful: ',m1,'x',m2,'x',m3,' points'
      
end program write_grid

