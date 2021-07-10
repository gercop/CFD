!----------------------------------------------------------------------------------------------------------------------------------!
! module which defines and calculates boundary conditions.
!----------------------------------------------------------------------------------------------------------------------------------!
module filter

#include 'setup.h'

!----- Modules --------------------------------------------------------------------------------------------------------------------!
#if defined(DEBUG)
   use debug
#endif ! DEBUG
   use params
   use grid
   use flow_prop
   use inp_outp

   implicit none

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! local (x->i, y->j)
   integer, private                               :: i, j, k
   integer,parameter,private                      :: st_f = st_w+2
   integer,parameter,private                      :: st_mid_f = (st_w+2)/2+1
   real,parameter, private                        :: factor = 1.0d0
   real*8,dimension(st_f,st_f),private            :: a
   real*8,dimension(st_f),private                 :: b
   real*8,dimension(nx,st_f),private              :: coef_filter_x
   real*8,dimension(ny,st_f),private              :: coef_filter_y
   real, dimension(nx,ny),private                 :: d_in_1, du_in_1, dv_in_1, te_in_1
contains


!----------------------------------------------------------------------------------------------------------------------------------!
! first subroutine: initialize filter coefficients
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine filter_init
   integer                                        :: fact

!---x-direction   
   
!---1st four points
   do i = 1,4
!---assigne a
     j=1
     do k = 1,st_f
        a(j,k) = 1.0d0 
     end do 
     do j = 2,st_f-1
       do k = 1,st_f
         a(j,k) = (pgrid_x(k)-pgrid_x(i))**dble(j-1)
       end do
     end do  
     j = st_f
     do k = 1,st_f
       a(j,k) = (-1.0d0)**(i+k)
     end do
!---assigne b
     b(1)=1.0d0
     do j = 2,st_f
       b(j) = 0.0d0
     end do
!---solve matrix
     call linsolve( a, b, st_f, st_f)  
!---assigne coeff
     do k = 1,st_f
       coef_filter_x(i,k) = b(k)
       write(*,*) 'i=',i,'k=,',k,',coef_filter =',coef_filter_x(i,k)
     end do
   end do
   
!---central points
   do i = 5,nx-4
!---assigne a
     j=1
     do k = 1,st_f
        a(j,k) = 1.0d0 
     end do 
     do j = 2,st_f-1
       do k = 1,st_f
         a(j,k) = (pgrid_x(i-st_mid_f+k)-pgrid_x(i))**dble(j-1)
       end do
     end do  
     j = st_f
     do k = 1,st_f
       a(j,k) = (-1.0d0)**((st_f/2)-st_f+k)
     end do
!---assigne b
     b(1)=1.0d0
     do j = 2,st_f
       b(j) = 0.0d0
     end do
!---solve matrix
     call linsolve( a, b, st_f, st_f)  
!---assigne coeff
     do k = 1,st_f
       coef_filter_x(i,k) = b(k)
       write(*,*) 'i=',i,'k=,',k,',coef_filter =',coef_filter_x(i,k)
     end do
   end do
   
!---last 4 points
   do i = nx-3,nx
!---assigne a
     j=1
     do k = 1,st_f
        a(j,k) = 1.0d0 
     end do 
     do j = 2,st_f-1
       do k = 1,st_f
         a(j,k) = (pgrid_x(nx-9+k)-pgrid_x(i))**dble(j-1)
       end do
     end do  
     j = st_f
     do k = 1,st_f
       a(j,k) = (-1.0d0)**(i+k)
     end do
!---assigne b
     b(1)=1.0d0
     do j = 2,st_f
       b(j) = 0.0d0
     end do
!---solve matrix
     call linsolve( a, b, st_f, st_f)  
!---assigne coeff
     do k = 1,st_f
       coef_filter_x(i,k) = b(k)
       write(*,*) 'i=',i,'k=,',k,',coef_filter =',coef_filter_x(i,k)
     end do
   end do



!---y-direction

!---1st four points
   do i = 1,4
!---assigne a
     j=1
     do k = 1,st_f
        a(j,k) = 1.0d0 
     end do 
     do j = 2,st_f-1
       do k = 1,st_f
         a(j,k) = (pgrid_y(k)-pgrid_y(i))**dble(j-1)
       end do
     end do  
     j = st_f
     do k = 1,st_f
       a(j,k) = (-1.0d0)**(i+k)
     end do
!---assigne b
     b(1)=1.0d0
     do j = 2,st_f
       b(j) = 0.0d0
     end do
!---solve matrix
     call linsolve( a, b, st_f, st_f)  
!---assigne coeff
     do k = 1,st_f
       coef_filter_y(i,k) = b(k)
       write(*,*) 'i=',i,'k=,',k,',coef_filter =',coef_filter_y(i,k)
     end do
   end do
   
!---central points
   do i = 5,ny-4
!---assigne a
     j=1
     do k = 1,st_f
        a(j,k) = 1.0d0 
     end do 
     do j = 2,st_f-1
       do k = 1,st_f
         a(j,k) = (pgrid_y(i-st_mid_f+k)-pgrid_y(i))**dble(j-1)
       end do
     end do  
     j = st_f
     do k = 1,st_f
       a(j,k) = (-1.0d0)**((st_f/2)-st_f+k)
     end do
!---assigne b
     b(1)=1.0d0
     do j = 2,st_f
       b(j) = 0.0d0
     end do
!---solve matrix
     call linsolve( a, b, st_f, st_f)  
!---assigne coeff
     do k = 1,st_f
       coef_filter_y(i,k) = b(k)
       write(*,*) 'i=',i,'k=,',k,',coef_filter =',coef_filter_y(i,k)
     end do
   end do
   
!---last 4 points
   do i = ny-3,ny
!---assigne a
     j=1
     do k = 1,st_f
        a(j,k) = 1.0d0 
     end do 
     do j = 2,st_f-1
       do k = 1,st_f
         a(j,k) = (pgrid_y(ny-9+k)-pgrid_y(i))**dble(j-1)
       end do
     end do  
     j = st_f
     do k = 1,st_f
       a(j,k) = (-1.0d0)**(i+k)
     end do
!---assigne b
     b(1)=1.0d0
     do j = 2,st_f
       b(j) = 0.0d0
     end do
!---solve matrix
     call linsolve( a, b, st_f, st_f)  
!---assigne coeff
     do k = 1,st_f
       coef_filter_y(i,k) = b(k)
       write(*,*) 'i=',i,'k=,',k,',coef_filter =',coef_filter_y(i,k)
     end do
   end do

   end subroutine filter_init

!----------------------------------------------------------------------------------------------------------------------------------!
! second subroutine: Filter in x 
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine filter_x(d_in, du_in, dv_in, te_in)

!----- Variables ------------------------------------------------------------------------------------------------------------------!

! in and output
   real, dimension(nx,ny), intent(inout)          :: d_in, du_in, dv_in, te_in
!   real, dimension(nx,ny)                         :: d_in_1, du_in_1, dv_in_1, te_in_1

   du_in_1=du_in
   dv_in_1=dv_in
   d_in_1=d_in
   te_in_1=te_in
!   open(11,file='du_test.dat',status='unknown')
!   do j = 1,ny
!     write(11,*) du_in(3,j)-dens_ana(3,j)*uvel_ana(3,j)
!   end do
!   close(11)
   
   
#if (SPACE_INT == 1)
 
!       i=2
!   do j=2,ny-1
!       du_in(i,j) = (-15.0d0*du_in(i,j)+4.0d0*du_in(i+1,j)-6.0d0*du_in(i+2,j)+4.0d0*du_in(i+3,j)-du_in(i+4,j))/(16.0d0)
!       dv_in(i,j) = (-15.0d0*dv_in(i,j)+4.0d0*dv_in(i+1,j)-6.0d0*dv_in(i+2,j)+4.0d0*dv_in(i+3,j)-dv_in(i+4,j))/(16.0d0)
!       d_in(i,j) = (-15.0d0*d_in(i,j)+4.0d0*d_in(i+1,j)-6.0d0*d_in(i+2,j)+4.0d0*d_in(i+3,j)-d_in(i+4,j))/(16.0d0)
!       te_in(i,j) = (-15.0d0*te_in(i,j)+4.0d0*te_in(i+1,j)-6.0d0*te_in(i+2,j)+4.0d0*te_in(i+3,j)-te_in(i+4,j))/(16.0d0)
!   end do    

    i=2
    do j = 2,ny-1
       du_in(i,j) = (1.0d0*du_in(i-1,j)+12.0d0*du_in(i,j)+6.0d0*du_in(i+1,j)-4.0d0*du_in(i+2,j)+du_in(i+3,j))/(16.0d0)
       dv_in(i,j) = (1.0d0*dv_in(i-1,j)+12.0d0*dv_in(i,j)+6.0d0*dv_in(i+1,j)-4.0d0*dv_in(i+2,j)+dv_in(i+3,j))/(16.0d0)
       d_in(i,j) = (1.0d0*d_in(i-1,j)+12.0d0*d_in(i,j)+6.0d0*d_in(i+1,j)-4.0d0*d_in(i+2,j)+d_in(i+3,j))/(16.0d0)
       te_in(i,j) = (1.0d0*te_in(i-1,j)+12.0d0*te_in(i,j)+6.0d0*te_in(i+1,j)-4.0d0*te_in(i+2,j)+te_in(i+3,j))/(16.0d0)
    end do  
   

   do i = 3,nx-2
     do j = 2,ny-1
       du_in(i,j) = (-du_in(i-2,j)+4.0d0*du_in(i-1,j)+10.0d0*factor*du_in(i,j)+4.0d0*du_in(i+1,j)-du_in(i+2,j))/(6.0d0+10.0d0*factor)
       dv_in(i,j) = (-dv_in(i-2,j)+4.0d0*dv_in(i-1,j)+10.0d0*factor*dv_in(i,j)+4.0d0*dv_in(i+1,j)-dv_in(i+2,j))/(6.0d0+10.0d0*factor)
       d_in(i,j) = (-d_in(i-2,j)+4.0d0*d_in(i-1,j)+10.0d0*factor*d_in(i,j)+4.0d0*d_in(i+1,j)-d_in(i+2,j))/(6.0d0+10.0d0*factor)
       te_in(i,j) = (-te_in(i-2,j)+4.0d0*te_in(i-1,j)+10.0d0*factor*te_in(i,j)+4.0d0*te_in(i+1,j)-te_in(i+2,j))/(6.0d0+10.0d0*factor)
     end do
   end do    

    j=2
   do i = 2,nx-1
      du_in(i,j) = (1.0d0*du_in(i,j-1)+12.0d0*du_in(i,j)+6.0d0*du_in(i,j+1)-4.0d0*du_in(i,j+2)+du_in(i,j+3))/(16.0d0)
      dv_in(i,j) = (1.0d0*dv_in(i,j-1)+12.0d0*dv_in(i,j)+6.0d0*dv_in(i,j+1)-4.0d0*dv_in(i,j+2)+dv_in(i,j+3))/(16.0d0)
      d_in(i,j) = (1.0d0*d_in(i,j-1)+12.0d0*d_in(i,j)+6.0d0*d_in(i,j+1)-4.0d0*d_in(i,j+2)+d_in(i,j+3))/(16.0d0)
      te_in(i,j) = (1.0d0*te_in(i,j-1)+12.0d0*te_in(i,j)+6.0d0*te_in(i,j+1)-4.0d0*te_in(i,j+2)+te_in(i,j+3))/(16.0d0)
   end do    
  

   do j = 3,ny-2
     do i = 2,nx-1
       du_in(i,j) = (-du_in(i,j-2)+4.0d0*du_in(i,j-1)+10.0d0*factor*du_in(i,j)+4.0d0*du_in(i,j+1)-du_in(i,j+2))/(6.0d0+10.0d0*factor)
       dv_in(i,j) = (-dv_in(i,j-2)+4.0d0*dv_in(i,j-1)+10.0d0*factor*dv_in(i,j)+4.0d0*dv_in(i,j+1)-dv_in(i,j+2))/(6.0d0+10.0d0*factor)
       d_in(i,j) = (-d_in(i,j-2)+4.0d0*d_in(i,j-1)+10.0d0*factor*d_in(i,j)+4.0d0*d_in(i,j+1)-d_in(i,j+2))/(6.0d0+10.0d0*factor)
       te_in(i,j) = (-te_in(i,j-2)+4.0d0*te_in(i,j-1)+10.0d0*factor*te_in(i,j)+4.0d0*te_in(i,j+1)-te_in(i,j+2))/(6.0d0+10.0d0*factor)
     end do
   end do    
#else
 do i = 2,4
   do j=2,ny-1
      du_in_1(i,j) = coef_filter_x(i,1)*du_in(1,j)+coef_filter_x(i,2)*du_in(2,j)+coef_filter_x(i,3)*du_in(3,j)&
		  +coef_filter_x(i,4)*du_in(4,j)+coef_filter_x(i,5)*du_in(5,j)+coef_filter_x(i,6)*du_in(6,j)&
		  +coef_filter_x(i,7)*du_in(7,j)+coef_filter_x(i,8)*du_in(8,j)+coef_filter_x(i,9)*du_in(9,j)
       dv_in_1(i,j) = coef_filter_x(i,1)*dv_in(1,j)+coef_filter_x(i,2)*dv_in(2,j)+coef_filter_x(i,3)*dv_in(3,j)&
		  +coef_filter_x(i,4)*dv_in(4,j)+coef_filter_x(i,5)*dv_in(5,j)+coef_filter_x(i,6)*dv_in(6,j)&
		  +coef_filter_x(i,7)*dv_in(7,j)+coef_filter_x(i,8)*dv_in(8,j)+coef_filter_x(i,9)*dv_in(9,j)
       d_in_1(i,j) = coef_filter_x(i,1)*d_in(1,j)+coef_filter_x(i,2)*d_in(2,j)+coef_filter_x(i,3)*d_in(3,j)&
		  +coef_filter_x(i,4)*d_in(4,j)+coef_filter_x(i,5)*d_in(5,j)+coef_filter_x(i,6)*d_in(6,j)&
		  +coef_filter_x(i,7)*d_in(7,j)+coef_filter_x(i,8)*d_in(8,j)+coef_filter_x(i,9)*d_in(9,j)
       te_in_1(i,j) = coef_filter_x(i,1)*te_in(1,j)+coef_filter_x(i,2)*te_in(2,j)+coef_filter_x(i,3)*te_in(3,j)&
		  +coef_filter_x(i,4)*te_in(4,j)+coef_filter_x(i,5)*te_in(5,j)+coef_filter_x(i,6)*te_in(6,j)&
		  +coef_filter_x(i,7)*te_in(7,j)+coef_filter_x(i,8)*te_in(8,j)+coef_filter_x(i,9)*te_in(9,j)
     end do
  end do
  
  do i = 5,nx-4
    do j=2,ny-1
       du_in_1(i,j) = coef_filter_x(i,1)*du_in(i-4,j)+coef_filter_x(i,2)*du_in(i-3,j)+coef_filter_x(i,3)*du_in(i-2,j)&
		   +coef_filter_x(i,4)*du_in(i-1,j)+coef_filter_x(i,5)*du_in(i,j)+coef_filter_x(i,6)*du_in(i+1,j)&
		   +coef_filter_x(i,7)*du_in(i+2,j)+coef_filter_x(i,8)*du_in(i+3,j)+coef_filter_x(i,9)*du_in(i+4,j)
       dv_in_1(i,j) = coef_filter_x(i,1)*dv_in(i-4,j)+coef_filter_x(i,2)*dv_in(i-3,j)+coef_filter_x(i,3)*dv_in(i-2,j)&
		   +coef_filter_x(i,4)*dv_in(i-1,j)+coef_filter_x(i,5)*dv_in(i,j)+coef_filter_x(i,6)*dv_in(i+1,j)&
		   +coef_filter_x(i,7)*dv_in(i+2,j)+coef_filter_x(i,8)*dv_in(i+3,j)+coef_filter_x(i,9)*dv_in(i+4,j)
       d_in_1(i,j) = coef_filter_x(i,1)*d_in(i-4,j)+coef_filter_x(i,2)*d_in(i-3,j)+coef_filter_x(i,3)*d_in(i-2,j)&
		   +coef_filter_x(i,4)*d_in(i-1,j)+coef_filter_x(i,5)*d_in(i,j)+coef_filter_x(i,6)*d_in(i+1,j)&
		   +coef_filter_x(i,7)*d_in(i+2,j)+coef_filter_x(i,8)*d_in(i+3,j)+coef_filter_x(i,9)*d_in(i+4,j)
       te_in_1(i,j) = coef_filter_x(i,1)*te_in(i-4,j)+coef_filter_x(i,2)*te_in(i-3,j)+coef_filter_x(i,3)*te_in(i-2,j)&
		   +coef_filter_x(i,4)*te_in(i-1,j)+coef_filter_x(i,5)*te_in(i,j)+coef_filter_x(i,6)*te_in(i+1,j)&
		   +coef_filter_x(i,7)*te_in(i+2,j)+coef_filter_x(i,8)*te_in(i+3,j)+coef_filter_x(i,9)*te_in(i+4,j)

    end do
  end do
  do i =nx-3,nx
    do j=2,ny-1
       du_in_1(i,j) = coef_filter_x(i,1)*du_in(nx-8,j)+coef_filter_x(i,2)*du_in(nx-7,j)+coef_filter_x(i,3)*du_in(nx-6,j)&
		       +coef_filter_x(i,4)*du_in(nx-5,j)+coef_filter_x(i,5)*du_in(nx-4,j)+coef_filter_x(i,6)*du_in(nx-3,j)&
		       +coef_filter_x(i,7)*du_in(nx-2,j)+coef_filter_x(i,8)*du_in(nx-1,j)+coef_filter_x(i,9)*du_in(nx,j)
       dv_in_1(i,j) = coef_filter_x(i,1)*dv_in(nx-8,j)+coef_filter_x(i,2)*dv_in(nx-7,j)+coef_filter_x(i,3)*dv_in(nx-6,j)&
		       +coef_filter_x(i,4)*dv_in(nx-5,j)+coef_filter_x(i,5)*dv_in(nx-4,j)+coef_filter_x(i,6)*dv_in(nx-3,j)&
		       +coef_filter_x(i,7)*dv_in(nx-2,j)+coef_filter_x(i,8)*dv_in(nx-1,j)+coef_filter_x(i,9)*dv_in(nx,j)
       d_in_1(i,j) = coef_filter_x(i,1)*d_in(nx-8,j)+coef_filter_x(i,2)*d_in(nx-7,j)+coef_filter_x(i,3)*d_in(nx-6,j)&
			+coef_filter_x(i,4)*d_in(nx-5,j)+coef_filter_x(i,5)*d_in(nx-4,j)+coef_filter_x(i,6)*d_in(nx-3,j)&
			+coef_filter_x(i,7)*d_in(nx-2,j)+coef_filter_x(i,8)*d_in(nx-1,j)+coef_filter_x(i,9)*d_in(nx,j)
       te_in_1(i,j) = coef_filter_x(i,1)*te_in(nx-8,j)+coef_filter_x(i,2)*te_in(nx-7,j)+coef_filter_x(i,3)*te_in(nx-6,j)&
			+coef_filter_x(i,4)*te_in(nx-5,j)+coef_filter_x(i,5)*te_in(nx-4,j)+coef_filter_x(i,6)*te_in(nx-3,j)&
			+coef_filter_x(i,7)*te_in(nx-2,j)+coef_filter_x(i,8)*te_in(nx-1,j)+coef_filter_x(i,9)*te_in(nx,j)
    end do
  end do
   du_in(:,:)=du_in_1(:,:)
   dv_in(:,:)=dv_in_1(:,:)
   d_in(:,:)=d_in_1(:,:)
   te_in(:,:)=te_in_1(:,:)


#endif
   end subroutine filter_x
   
!----------------------------------------------------------------------------------------------------------------------------------!
! third subroutine: Filter in y 
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine filter_y(d_in, du_in, dv_in, te_in)

!----- Variables ------------------------------------------------------------------------------------------------------------------!

! in and output
   real, dimension(nx,ny), intent(inout)          :: d_in, du_in, dv_in, te_in

!   du_in = du_in -dens_ana*uvel_ana
!   dv_in = dv_in -dens_ana*vvel_ana
!   d_in = d_in -dens_ana
!   te_in = te_in -dens_ana*(temp_ana*igg1m2+0.5*( uvel_ana*uvel_ana+vvel_ana*vvel_ana ) )
   du_in_1=du_in 
   dv_in_1=dv_in
   d_in_1=d_in
   te_in_1=te_in
!   open(11,file='du_test.dat',status='unknown')
!   do j = 1,ny
!     write(11,*) du_in(3,j)-dens_ana(3,j)*uvel_ana(3,j)
!   end do
!   close(11)
   
   
#if (SPACE_INT == 1)
 
!       i=2
!   do j=2,ny-1
!       du_in(i,j) = (-15.0d0*du_in(i,j)+4.0d0*du_in(i+1,j)-6.0d0*du_in(i+2,j)+4.0d0*du_in(i+3,j)-du_in(i+4,j))/(16.0d0)
!       dv_in(i,j) = (-15.0d0*dv_in(i,j)+4.0d0*dv_in(i+1,j)-6.0d0*dv_in(i+2,j)+4.0d0*dv_in(i+3,j)-dv_in(i+4,j))/(16.0d0)
!       d_in(i,j) = (-15.0d0*d_in(i,j)+4.0d0*d_in(i+1,j)-6.0d0*d_in(i+2,j)+4.0d0*d_in(i+3,j)-d_in(i+4,j))/(16.0d0)
!       te_in(i,j) = (-15.0d0*te_in(i,j)+4.0d0*te_in(i+1,j)-6.0d0*te_in(i+2,j)+4.0d0*te_in(i+3,j)-te_in(i+4,j))/(16.0d0)
!   end do    

!    i=2
!    do j = 2,ny-1
!       du_in(i,j) = (1.0d0*du_in(i-1,j)+12.0d0*du_in(i,j)+6.0d0*du_in(i+1,j)-4.0d0*du_in(i+2,j)+du_in(i+3,j))/(16.0d0)
!       dv_in(i,j) = (1.0d0*dv_in(i-1,j)+12.0d0*dv_in(i,j)+6.0d0*dv_in(i+1,j)-4.0d0*dv_in(i+2,j)+dv_in(i+3,j))/(16.0d0)
!       d_in(i,j) = (1.0d0*d_in(i-1,j)+12.0d0*d_in(i,j)+6.0d0*d_in(i+1,j)-4.0d0*d_in(i+2,j)+d_in(i+3,j))/(16.0d0)
!       te_in(i,j) = (1.0d0*te_in(i-1,j)+12.0d0*te_in(i,j)+6.0d0*te_in(i+1,j)-4.0d0*te_in(i+2,j)+te_in(i+3,j))/(16.0d0)
!    end do  
   

!   do i = 3,nx-2
!     do j = 2,ny-1
!       du_in(i,j) = (-du_in(i-2,j)+4.0d0*du_in(i-1,j)+10.0d0*factor*du_in(i,j)+4.0d0*du_in(i+1,j)-du_in(i+2,j))/(6.0d0+10.0d0*factor)
!       dv_in(i,j) = (-dv_in(i-2,j)+4.0d0*dv_in(i-1,j)+10.0d0*factor*dv_in(i,j)+4.0d0*dv_in(i+1,j)-dv_in(i+2,j))/(6.0d0+10.0d0*factor)
!       d_in(i,j) = (-d_in(i-2,j)+4.0d0*d_in(i-1,j)+10.0d0*factor*d_in(i,j)+4.0d0*d_in(i+1,j)-d_in(i+2,j))/(6.0d0+10.0d0*factor)
!       te_in(i,j) = (-te_in(i-2,j)+4.0d0*te_in(i-1,j)+10.0d0*factor*te_in(i,j)+4.0d0*te_in(i+1,j)-te_in(i+2,j))/(6.0d0+10.0d0*factor)
!     end do
!   end do    

    j=2
   do i = 2,nx-1
      du_in(i,j) = (1.0d0*du_in(i,j-1)+12.0d0*du_in(i,j)+6.0d0*du_in(i,j+1)-4.0d0*du_in(i,j+2)+du_in(i,j+3))/(16.0d0)
      dv_in(i,j) = (1.0d0*dv_in(i,j-1)+12.0d0*dv_in(i,j)+6.0d0*dv_in(i,j+1)-4.0d0*dv_in(i,j+2)+dv_in(i,j+3))/(16.0d0)
      d_in(i,j) = (1.0d0*d_in(i,j-1)+12.0d0*d_in(i,j)+6.0d0*d_in(i,j+1)-4.0d0*d_in(i,j+2)+d_in(i,j+3))/(16.0d0)
      te_in(i,j) = (1.0d0*te_in(i,j-1)+12.0d0*te_in(i,j)+6.0d0*te_in(i,j+1)-4.0d0*te_in(i,j+2)+te_in(i,j+3))/(16.0d0)
   end do    
  

   do j = 3,ny-2
     do i = 2,nx-1
       du_in(i,j) = (-du_in(i,j-2)+4.0d0*du_in(i,j-1)+10.0d0*factor*du_in(i,j)+4.0d0*du_in(i,j+1)-du_in(i,j+2))/(6.0d0+10.0d0*factor)
       dv_in(i,j) = (-dv_in(i,j-2)+4.0d0*dv_in(i,j-1)+10.0d0*factor*dv_in(i,j)+4.0d0*dv_in(i,j+1)-dv_in(i,j+2))/(6.0d0+10.0d0*factor)
       d_in(i,j) = (-d_in(i,j-2)+4.0d0*d_in(i,j-1)+10.0d0*factor*d_in(i,j)+4.0d0*d_in(i,j+1)-d_in(i,j+2))/(6.0d0+10.0d0*factor)
       te_in(i,j) = (-te_in(i,j-2)+4.0d0*te_in(i,j-1)+10.0d0*factor*te_in(i,j)+4.0d0*te_in(i,j+1)-te_in(i,j+2))/(6.0d0+10.0d0*factor)
     end do
   end do    
!  open(12,file='testing.dat',status='unknown')
!  do j = 1,ny
!    write(12,*) pgrid_y(j),du_in_1(3,j)-dens_ana(3,j)*uvel_ana(3,j)
!  end do
!  write(12,*) '&'
!  do j = 1,ny
!    write(12,*) pgrid_y(j),du_in(3,j)-dens_ana(3,j)*uvel_ana(3,j)
!  end do
!  close(12)  
!  stop 'ssaannddii'
#else

!---and now in y

  do j = 2,4
    do i=2,nx
       du_in_1(i,j) = coef_filter_y(j,1)*du_in(i,1)+coef_filter_y(j,2)*du_in(i,2)+coef_filter_y(j,3)*du_in(i,3)&
		   +coef_filter_y(j,4)*du_in(i,4)+coef_filter_y(j,5)*du_in(i,5)+coef_filter_y(j,6)*du_in(i,6)&
		   +coef_filter_y(j,7)*du_in(i,7)+coef_filter_y(j,8)*du_in(i,8)+coef_filter_y(j,9)*du_in(i,9)
       dv_in_1(i,j) = coef_filter_y(j,1)*dv_in(i,1)+coef_filter_y(j,2)*dv_in(i,2)+coef_filter_y(j,3)*dv_in(i,3)&
		   +coef_filter_y(j,4)*dv_in(i,4)+coef_filter_y(j,5)*dv_in(i,5)+coef_filter_y(j,6)*dv_in(i,6)&
		   +coef_filter_y(j,7)*dv_in(i,7)+coef_filter_y(j,8)*dv_in(i,8)+coef_filter_y(j,9)*dv_in(i,9)
       d_in_1(i,j) = coef_filter_y(j,1)*d_in(i,1)+coef_filter_y(j,2)*d_in(i,2)+coef_filter_y(j,3)*d_in(i,3)&
		   +coef_filter_y(j,4)*d_in(i,4)+coef_filter_y(j,5)*d_in(i,5)+coef_filter_y(j,6)*d_in(i,6)&
		   +coef_filter_y(j,7)*d_in(i,7)+coef_filter_y(j,8)*d_in(i,8)+coef_filter_y(j,9)*d_in(i,9)
       te_in_1(i,j) = coef_filter_y(j,1)*te_in(i,1)+coef_filter_y(j,2)*te_in(i,2)+coef_filter_y(j,3)*te_in(i,3)&
		   +coef_filter_x(j,4)*te_in(i,4)+coef_filter_x(j,5)*te_in(i,5)+coef_filter_y(j,6)*te_in(i,6)&
		   +coef_filter_x(j,7)*te_in(i,7)+coef_filter_x(j,8)*te_in(i,8)+coef_filter_y(j,9)*te_in(i,9)
!	 uvel(i,j) = coef_filter_x(j,1)*uvel(i,1)+coef_filter_x(j,2)*uvel(i,2)+coef_filter_x(j,3)*uvel(i,3)&
!		    +coef_filter_x(j,4)*uvel(i,4)+coef_filter_x(j,5)*uvel(i,5)+coef_filter_x(j,6)*uvel(i,6)&
!		    +coef_filter_x(j,7)*uvel(i,7)+coef_filter_x(j,8)*uvel(i,8)+coef_filter_x(j,9)*uvel(i,9)
!	dv_in_1(i,j) = coef_filter_x(j,1)*vvel(i,1)+coef_filter_x(j,2)*vvel(i,2)+coef_filter_x(j,3)*vvel(i,3)&
!		    +coef_filter_x(j,4)*vvel(i,4)+coef_filter_x(j,5)*vvel(i,5)+coef_filter_x(j,6)*vvel(i,6)&
!		    +coef_filter_x(j,7)*vvel(i,7)+coef_filter_x(j,8)*vvel(i,8)+coef_filter_x(j,9)*vvel(i,9)
!	d_in_1(i,j) = coef_filter_x(j,1)*d_in(i,1)+coef_filter_x(j,2)*d_in(i,2)+coef_filter_x(j,3)*d_in(i,3)&
!		    +coef_filter_x(j,4)*d_in(i,4)+coef_filter_x(j,5)*d_in(i,5)+coef_filter_x(j,6)*d_in(i,6)&
!		    +coef_filter_x(j,7)*d_in(i,7)+coef_filter_x(j,8)*d_in(i,8)+coef_filter_x(j,9)*d_in(i,9)
!	te_in_1(i,j) = coef_filter_x(j,1)*temp(i,1)+coef_filter_x(j,2)*temp(i,2)+coef_filter_x(j,3)*temp(i,3)&
!		    +coef_filter_x(j,4)*temp(i,4)+coef_filter_x(j,5)*temp(i,5)+coef_filter_x(j,6)*temp(i,6)&
!		    +coef_filter_x(j,7)*temp(i,7)+coef_filter_x(j,8)*temp(i,8)+coef_filter_x(j,9)*temp(i,9)
     end do
  end do
  
  do j = 5,ny-4
    do i=2,nx
       du_in_1(i,j) = coef_filter_y(j,1)*du_in(i,j-4)+coef_filter_y(j,2)*du_in(i,j-3)+coef_filter_y(j,3)*du_in(i,j-2)&
		   +coef_filter_y(j,4)*du_in(i,j-1)+coef_filter_y(j,5)*du_in(i,j)+coef_filter_y(j,6)*du_in(i,j+1)&
		   +coef_filter_y(j,7)*du_in(i,j+2)+coef_filter_y(j,8)*du_in(i,j+3)+coef_filter_y(j,9)*du_in(i,j+4)
       dv_in_1(i,j) = coef_filter_y(j,1)*dv_in(i,j-4)+coef_filter_y(j,2)*dv_in(i,j-3)+coef_filter_y(j,3)*dv_in(i,j-2)&
		   +coef_filter_y(j,4)*dv_in(i,j-1)+coef_filter_y(j,5)*dv_in(i,j)+coef_filter_y(j,6)*dv_in(i,j+1)&
		   +coef_filter_y(j,7)*dv_in(i,j+2)+coef_filter_y(j,8)*dv_in(i,j+3)+coef_filter_y(j,9)*dv_in(i,j+4)
       d_in_1(i,j) = coef_filter_y(j,1)*d_in(i,j-4)+coef_filter_y(j,2)*d_in(i,j-3)+coef_filter_y(j,3)*d_in(i,j-2)&
		   +coef_filter_y(j,4)*d_in(i,j-1)+coef_filter_y(j,5)*d_in(i,j)+coef_filter_y(j,6)*d_in(i,j+1)&
		   +coef_filter_y(j,7)*d_in(i,j+2)+coef_filter_y(j,8)*d_in(i,j+3)+coef_filter_y(j,9)*d_in(i,j+4)
       te_in_1(i,j) = coef_filter_y(j,1)*te_in(i,j-4)+coef_filter_y(j,2)*te_in(i,j-3)+coef_filter_y(j,3)*te_in(i,j-2)&
		   +coef_filter_y(j,4)*te_in(i,j-1)+coef_filter_y(j,5)*te_in(i,j)+coef_filter_y(j,6)*te_in(i,j+1)&
		   +coef_filter_y(j,7)*te_in(i,j+2)+coef_filter_y(j,8)*te_in(i,j+3)+coef_filter_y(j,9)*te_in(i,j+4)
!	uvel(i,j) = coef_filter_x(st_mid_f,1)*uvel(i,j-4)+coef_filter_x(st_mid_f,2)*uvel(i,j-3)+coef_filter_x(st_mid_f,3)*uvel(i,j-2)&
!		    +coef_filter_x(st_mid_f,4)*uvel(i,j-1)+coef_filter_x(st_mid_f,5)*uvel(i,j)+coef_filter_x(st_mid_f,6)*uvel(i,j+1)&
!		    +coef_filter_x(st_mid_f,7)*uvel(i,j+2)+coef_filter_x(st_mid_f,8)*uvel(i,j+3)+coef_filter_x(st_mid_f,9)*uvel(i,j+4)
!	vvel(i,j) = coef_filter_x(st_mid_f,1)*vvel(i,j-4)+coef_filter_x(st_mid_f,2)*vvel(i,j-3)+coef_filter_x(st_mid_f,3)*vvel(i,j-2)&
!		    +coef_filter_x(st_mid_f,4)*vvel(i,j-1)+coef_filter_x(st_mid_f,5)*vvel(i,j)+coef_filter_x(st_mid_f,6)*vvel(i,j+1)&
!		    +coef_filter_x(st_mid_f,7)*vvel(i,j+2)+coef_filter_x(st_mid_f,8)*vvel(i,j+3)+coef_filter_x(st_mid_f,9)*vvel(i,j+4)
!	d_in(i,j) = coef_filter_x(st_mid_f,1)*d_in(i,j-4)+coef_filter_x(st_mid_f,2)*d_in(i,j-3)+coef_filter_x(st_mid_f,3)*d_in(i,j-2)&
!		    +coef_filter_x(st_mid_f,4)*d_in(i,j-1)+coef_filter_x(st_mid_f,5)*d_in(i,j)+coef_filter_x(st_mid_f,6)*d_in(i,j+1)&
!		    +coef_filter_x(st_mid_f,7)*d_in(i,j+2)+coef_filter_x(st_mid_f,8)*d_in(i,j+3)+coef_filter_x(st_mid_f,9)*d_in(i,j+4)
!	temp(i,j) = coef_filter_x(st_mid_f,1)*temp(i,j-4)+coef_filter_x(st_mid_f,2)*temp(i,j-3)+coef_filter_x(st_mid_f,3)*temp(i,j-2)&
!		    +coef_filter_x(st_mid_f,4)*temp(i,j-1)+coef_filter_x(st_mid_f,5)*temp(i,j)+coef_filter_x(st_mid_f,6)*temp(i,j+1)&
!		    +coef_filter_x(st_mid_f,7)*temp(i,j+2)+coef_filter_x(st_mid_f,8)*temp(i,j+3)+coef_filter_x(st_mid_f,9)*temp(i,j+4)
    end do
  end do
 ! stop
  do j =ny-3,ny
    do i=2,ny
       du_in_1(i,ny-9+j) = coef_filter_y(j,1)*du_in(i,ny-8)+coef_filter_y(j,2)*du_in(i,ny-7)+coef_filter_y(j,3)*du_in(i,ny-6)&
			+coef_filter_y(j,4)*du_in(i,ny-5)+coef_filter_y(j,5)*du_in(i,ny-4)+coef_filter_y(j,6)*du_in(i,ny-3)&
			+coef_filter_y(j,7)*du_in(i,ny-2)+coef_filter_y(j,8)*du_in(i,ny-1)+coef_filter_y(j,9)*du_in(i,ny)
       dv_in_1(i,ny-9+j) = coef_filter_y(j,1)*dv_in(i,ny-8)+coef_filter_y(j,2)*dv_in(i,ny-7)+coef_filter_y(j,3)*dv_in(i,ny-6)&
			+coef_filter_y(j,4)*dv_in(i,ny-5)+coef_filter_y(j,5)*dv_in(i,ny-4)+coef_filter_y(j,6)*dv_in(i,ny-3)&
			+coef_filter_y(j,7)*dv_in(i,ny-2)+coef_filter_y(j,8)*dv_in(i,ny-1)+coef_filter_y(j,9)*dv_in(i,ny)
       d_in_1(i,ny-9+j) = coef_filter_y(j,1)*d_in(i,ny-8)+coef_filter_y(j,2)*d_in(i,ny-7)+coef_filter_y(j,3)*d_in(i,ny-6)&
			+coef_filter_y(j,4)*d_in(i,ny-5)+coef_filter_y(j,5)*d_in(i,ny-4)+coef_filter_y(j,6)*d_in(i,ny-3)&
			+coef_filter_y(j,7)*d_in(i,ny-2)+coef_filter_y(j,8)*d_in(i,ny-1)+coef_filter_y(j,9)*d_in(i,ny)
       te_in_1(i,ny-9+j) = coef_filter_y(j,1)*te_in(i,ny-8)+coef_filter_y(j,2)*te_in(i,ny-7)+coef_filter_y(j,3)*te_in(i,ny-6)&
			+coef_filter_y(j,4)*te_in(i,ny-5)+coef_filter_y(j,5)*te_in(i,ny-4)+coef_filter_y(j,6)*te_in(i,ny-3)&
			+coef_filter_y(j,7)*te_in(i,ny-2)+coef_filter_y(j,8)*te_in(i,ny-1)+coef_filter_y(j,9)*te_in(i,ny)
!	uvel(i,ny-9+j) = coef_filter_x(j,1)*uvel(i,ny-8)+coef_filter_x(j,2)*uvel(i,ny-7)+coef_filter_x(j,3)*uvel(i,ny-6)&
!			 +coef_filter_x(j,4)*uvel(i,ny-5)+coef_filter_x(j,5)*uvel(i,ny-4)+coef_filter_x(j,6)*uvel(i,ny-3)&
!			 +coef_filter_x(j,7)*uvel(i,ny-2)+coef_filter_x(j,8)*uvel(i,ny-1)+coef_filter_x(j,9)*uvel(i,ny)
!	vvel(i,ny-9+j) = coef_filter_x(j,1)*vvel(i,ny-8)+coef_filter_x(j,2)*vvel(i,ny-7)+coef_filter_x(j,3)*vvel(i,ny-6)&
!			 +coef_filter_x(j,4)*vvel(i,ny-5)+coef_filter_x(j,5)*vvel(i,ny-4)+coef_filter_x(j,6)*vvel(i,ny-3)&
!			 +coef_filter_x(j,7)*vvel(i,ny-2)+coef_filter_x(j,8)*vvel(i,ny-1)+coef_filter_x(j,9)*vvel(i,ny)
!	d_in(i,ny-9+j) = coef_filter_x(j,1)*d_in(i,ny-8)+coef_filter_x(j,2)*d_in(i,ny-7)+coef_filter_x(j,3)*d_in(i,ny-6)&
!			 +coef_filter_x(j,4)*d_in(i,ny-5)+coef_filter_x(j,5)*d_in(i,ny-4)+coef_filter_x(j,6)*d_in(i,ny-3)&
!			 +coef_filter_x(j,7)*d_in(i,ny-2)+coef_filter_x(j,8)*d_in(i,ny-1)+coef_filter_x(j,9)*d_in(i,ny)
!	temp(i,ny-9+j) = coef_filter_x(j,1)*temp(i,ny-8)+coef_filter_x(j,2)*temp(i,ny-7)+coef_filter_x(j,3)*temp(i,ny-6)&
!			 +coef_filter_x(j,4)*temp(i,ny-5)+coef_filter_x(j,5)*temp(i,ny-4)+coef_filter_x(j,6)*temp(i,ny-3)&
!			 +coef_filter_x(j,7)*temp(i,ny-2)+coef_filter_x(j,8)*temp(i,ny-1)+coef_filter_x(j,9)*temp(i,ny)
    end do
  end do
!  open(12,file='testing.dat',status='unknown')
!  do j = 1,ny
!    write(12,*) pgrid_y(j),du_in(3,j)-dens_ana(3,j)*uvel_ana(3,j)
!  end do
!  write(12,*) '&'
!  do j = 1,ny
!    write(12,*) pgrid_y(j),du_in_1(3,j)-dens_ana(3,j)*uvel_ana(3,j)
!  end do
!  close(12)  
!  stop 'ssaannddii'
  
   du_in(:,:)=du_in_1(:,:)!+dens_ana*uvel_ana
   dv_in(:,:)=dv_in_1(:,:)!+dens_ana*vvel_ana
   d_in(:,:)=d_in_1(:,:)!+dens_ana
   te_in(:,:)=te_in_1(:,:)!+dens_ana*(temp_ana*igg1m2+0.5*( uvel_ana*uvel_ana+vvel_ana*vvel_ana ) )
!   do i = 2,nx
!   call conserv_x(d_in, du_in, dv_in, te_in, i, 1, ny)
!   end do
#endif
   end subroutine filter_y
   
!----------------------------------------------------------------------------------------------------------------------------------!
! first function : computes the factorial
!----------------------------------------------------------------------------------------------------------------------------------!

function fact(n)
   integer                     :: n,fact,i
   fact=1
   do i = 1,n
     fact=fact*i
   end do 
end function

end module filter






!
!
!
!
!
!
!
 

!
!
!
!
!
!
!
!
