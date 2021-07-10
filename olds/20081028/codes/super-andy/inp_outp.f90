!----------------------------------------------------------------------------------------------------------------------------------!
! module which regulates input output operations.
!----------------------------------------------------------------------------------------------------------------------------------!
module inp_outp

#include 'setup.h'

!----- Modules --------------------------------------------------------------------------------------------------------------------!
#if defined(DEBUG)
   use debug
#endif ! DEBUG
   use flow_prop
   use params
   use grid

   implicit none

!----- Variables ------------------------------------------------------------------------------------------------------------------!
! IOS
   integer(KIND=4), parameter,private             :: ftape=11, funit=11, err=8
   integer(KIND=4), private                       :: imach, mt, mx, my, mz, nparam, ninf
   integer(KIND=4), dimension(512), private       :: times

   character(LEN=124), private                    :: fbase
   character(LEN=72), dimension(128), private     :: info

! superandi's test variables (couette)
   real, dimension(nx,ny), public, save           :: dens_ana, temp_ana, uvel_ana
   real, dimension(nx,ny), public, save           :: vvel_ana, pres_ana

contains

#if (READ_IC == 1)
!----------------------------------------------------------------------------------------------------------------------------------!
! first subroutine: Reads and calculates initial condition from Profkom.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine readi

! local variables
   integer                                        :: aerr, i, ioerr, j, k, switch
   real                                           :: peta
   real, dimension(:), allocatable                :: eta, u_init, v_init, d_init, t_init

!----- Read file (IOS) ------------------------------------------------------------------------------------------------------------!
! read file info
   call mkfname(RFNAME, fbase, imach, err)
   call readcd(ftape, funit, fbase, imach, mt, mz, my, mx, nparam, times, info, ninf)

! allocate memory
   allocate(eta(mx),u_init(mx),v_init(mx),d_init(mx),t_init(mx),stat=aerr)
   if (aerr.ne.0) then
      print *,'ERROR->INP_OUTP: Allocation error!'
      stop
   end if

! read data
   call readd(ftape, u_init, 1, 1, err)
   call readd(ftape, t_init, 1, 2, err)
   call readd(ftape, d_init, 1, 3, err)
   call readd(ftape, v_init, 1, 4, err)

   close(ftape)

! eta
   do i=1,mx
      eta(i) = dble(i-1)*deta
   end do  

! check when u-velocity switches from 1.0 to 0.0 in the free-stream (Problem of Profkom).
   switch = mx
   do i=2,mx
      if ((u_init(i)-u_init(i-1)).lt.-0.9d0) then
         switch = i-st_hw !-st_hw is used because of the long stencils for the interpolation.
         exit
      end if
   end do

! interpolation of initial data onto the physical grid of computation using polynomial interpolation.
   do i=1,nx
      do j=1,ny

! get eta coordinate for the xy-grid point and check at what location k, eta(k) >= peta.
         peta = pgrid_y(j)*sqrt(re/pgrid_x(i));
         do k=1,mx
            if (eta(k).ge.peta) then
               exit
            end if
         end do

         if (k.lt.switch) then

! the interpolation
         if (k.le.1+st_hw) then
! u-velocity
#if (SPACE_INT==1)
	    uvel(i,j) = poly(eta(1),eta(2),eta(3),peta,1)*u_init(1) &
                       +poly(eta(1),eta(2),eta(3),peta,2)*u_init(2) &
                       +poly(eta(1),eta(2),eta(3),peta,3)*u_init(3)
#else
	    uvel(i,j) = poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,1)*u_init(1) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,2)*u_init(2) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,3)*u_init(3) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,4)*u_init(4) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,5)*u_init(5) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,6)*u_init(6) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,7)*u_init(7)
#endif ! SPACE_INT
! v-velocity
#if (SPACE_INT==1)
            vvel(i,j) =(poly(eta(1),eta(2),eta(3),peta,1)*v_init(1) &
                       +poly(eta(1),eta(2),eta(3),peta,2)*v_init(2) &
                       +poly(eta(1),eta(2),eta(3),peta,3)*v_init(3))&
                       /sqrt(re*pgrid_x(i))
#else
            vvel(i,j) =(poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,1)*v_init(1) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,2)*v_init(2) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,3)*v_init(3) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,4)*v_init(4) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,5)*v_init(5) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,6)*v_init(6) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,7)*v_init(7))&
                       /sqrt(re*pgrid_x(i))
#endif ! SPACE_INT
! temperature
#if (SPACE_INT==1)
            temp(i,j) = poly(eta(1),eta(2),eta(3),peta,1)*t_init(1) &
                       +poly(eta(1),eta(2),eta(3),peta,2)*t_init(2) &
                       +poly(eta(1),eta(2),eta(3),peta,3)*t_init(3)
#else
            temp(i,j) = poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,1)*t_init(1) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,2)*t_init(2) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,3)*t_init(3) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,4)*t_init(4) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,5)*t_init(5) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,6)*t_init(6) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,7)*t_init(7)
#endif ! SPACE_INT
! density
#if (SPACE_INT==1)
            dens(i,j) = poly(eta(1),eta(2),eta(3),peta,1)*d_init(1) &
                       +poly(eta(1),eta(2),eta(3),peta,2)*d_init(2) &
                       +poly(eta(1),eta(2),eta(3),peta,3)*d_init(3)
#else
            dens(i,j) = poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,1)*d_init(1) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,2)*d_init(2) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,3)*d_init(3) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,4)*d_init(4) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,5)*d_init(5) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,6)*d_init(6) &
                       +poly(eta(1),eta(2),eta(3),eta(4),eta(5),eta(6),eta(7),peta,7)*d_init(7)
#endif ! SPACE_INT
         else if (k.ge.mx-st_hw) then
! u-velocity
#if (SPACE_INT==1)
            uvel(i,j) = poly(eta(mx-2),eta(mx-1),eta(mx),peta,1)*u_init(mx-2) &
                       +poly(eta(mx-2),eta(mx-1),eta(mx),peta,2)*u_init(mx-1) &
                       +poly(eta(mx-2),eta(mx-1),eta(mx),peta,3)*u_init(mx  )
#else
            uvel(i,j) = poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,1)*u_init(mx-6) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,2)*u_init(mx-5) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,3)*u_init(mx-4) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,4)*u_init(mx-3) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,5)*u_init(mx-2) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,6)*u_init(mx-1) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,7)*u_init(mx  )
#endif ! SPACE_INT
! v-velocity
#if (SPACE_INT==1)
            vvel(i,j) =(poly(eta(mx-2),eta(mx-1),eta(mx),peta,1)*v_init(mx-2) &
                       +poly(eta(mx-2),eta(mx-1),eta(mx),peta,2)*v_init(mx-1) &
                       +poly(eta(mx-2),eta(mx-1),eta(mx),peta,3)*v_init(mx  ))&
                       /sqrt(re*pgrid_x(i))
#else
            vvel(i,j) =(poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,1)*v_init(mx-6) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,2)*v_init(mx-5) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,3)*v_init(mx-4) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,4)*v_init(mx-3) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,5)*v_init(mx-2) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,6)*v_init(mx-1) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,7)*v_init(mx  ))&
                       /sqrt(re*pgrid_x(i))
#endif ! SPACE_INT
! temperature
#if (SPACE_INT==1)
            temp(i,j) = poly(eta(mx-2),eta(mx-1),eta(mx),peta,1)*t_init(mx-2) &
                       +poly(eta(mx-2),eta(mx-1),eta(mx),peta,2)*t_init(mx-1) &
                       +poly(eta(mx-2),eta(mx-1),eta(mx),peta,3)*t_init(mx  )
#else
            temp(i,j) = poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,1)*t_init(mx-6) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,2)*t_init(mx-5) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,3)*t_init(mx-4) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,4)*t_init(mx-3) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,5)*t_init(mx-2) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,6)*t_init(mx-1) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,7)*t_init(mx  )
#endif ! SPACE_INT
! density
#if (SPACE_INT==1)
            dens(i,j) = poly(eta(mx-2),eta(mx-1),eta(mx),peta,1)*d_init(mx-2) &
                       +poly(eta(mx-2),eta(mx-1),eta(mx),peta,2)*d_init(mx-1) &
                       +poly(eta(mx-2),eta(mx-1),eta(mx),peta,3)*d_init(mx  )
#else
            dens(i,j) = poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,1)*d_init(mx-6) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,2)*d_init(mx-5) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,3)*d_init(mx-4) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,4)*d_init(mx-3) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,5)*d_init(mx-2) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,6)*d_init(mx-1) &
                       +poly(eta(mx-6),eta(mx-5),eta(mx-4),eta(mx-3),eta(mx-2),eta(mx-1),eta(mx),peta,7)*d_init(mx  )
#endif ! SPACE_INT
         else
! u-velocity
#if (SPACE_INT==1)
            uvel(i,j) = poly(eta(k-1),eta(k),eta(k+1),peta,1)*u_init(k-1) &
                       +poly(eta(k-1),eta(k),eta(k+1),peta,2)*u_init(k  ) &
                       +poly(eta(k-1),eta(k),eta(k+1),peta,3)*u_init(k+1)
#else
            uvel(i,j) = poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,1)*u_init(k-3) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,2)*u_init(k-2) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,3)*u_init(k-1) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,4)*u_init(k  ) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,5)*u_init(k+1) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,6)*u_init(k+2) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,7)*u_init(k+3)
#endif ! SPACE_INT
! v-velocity
#if (SPACE_INT==1)
            vvel(i,j) =(poly(eta(k-1),eta(k),eta(k+1),peta,1)*v_init(k-1) &
                       +poly(eta(k-1),eta(k),eta(k+1),peta,2)*v_init(k  ) &
                       +poly(eta(k-1),eta(k),eta(k+1),peta,3)*v_init(k+1))&
                       /sqrt(re*pgrid_x(i))
#else
            vvel(i,j) =(poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,1)*v_init(k-3) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,2)*v_init(k-2) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,3)*v_init(k-1) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,4)*v_init(k  ) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,5)*v_init(k+1) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,6)*v_init(k+2) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,7)*v_init(k+3))&
                       /sqrt(re*pgrid_x(i))
#endif ! SPACE_INT
! temperature
#if (SPACE_INT==1)
            temp(i,j) = poly(eta(k-1),eta(k),eta(k+1),peta,1)*t_init(k-1) &
                       +poly(eta(k-1),eta(k),eta(k+1),peta,2)*t_init(k  ) &
                       +poly(eta(k-1),eta(k),eta(k+1),peta,3)*t_init(k+1)
#else
            temp(i,j) = poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,1)*t_init(k-3) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,2)*t_init(k-2) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,3)*t_init(k-1) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,4)*t_init(k  ) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,5)*t_init(k+1) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,6)*t_init(k+2) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,7)*t_init(k+3)
#endif ! SPACE_INT
! density
#if (SPACE_INT==1)
            dens(i,j) = poly(eta(k-1),eta(k),eta(k+1),peta,1)*d_init(k-1) &
                       +poly(eta(k-1),eta(k),eta(k+1),peta,2)*d_init(k  ) &
                       +poly(eta(k-1),eta(k),eta(k+1),peta,3)*d_init(k+1)
#else
            dens(i,j) = poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,1)*d_init(k-3) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,2)*d_init(k-2) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,3)*d_init(k-1) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,4)*d_init(k  ) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,5)*d_init(k+1) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,6)*d_init(k+2) &
                       +poly(eta(k-3),eta(k-2),eta(k-1),eta(k),eta(k+1),eta(k+2),eta(k+3),peta,7)*d_init(k+3)
#endif ! SPACE_INT
         end if

         else
! u-velocity
            uvel(i,j) = uvel(i,j-1)
! v-velocity
            vvel(i,j) = vvel(i,j-1)
! temperature
            temp(i,j) = temp(i,j-1)
! density
            dens(i,j) = dens(i,j-1)
         end if

! calculate pressure
!         pres(i,j) = igm2*dens(i,j)*temp(i,j)
         pres(i,j) = igm2

      end do
   end do

! calculation of conservative variables
   call conserv

! fill up inflow arrays
   do j=1,ny
      u_infl(j) = uvel(1,j)
      v_infl(j) = vvel(1,j)
      t_infl(j) = temp(1,j)
      d_infl(j) = dens(1,j)
   end do
   uvel_ana(:,:) = uvel(:,:)
   vvel_ana(:,:) = vvel(:,:)
   temp_ana(:,:) = temp(:,:)
   dens_ana(:,:) = dens(:,:)
   pres_ana(:,:) = pres(:,:)

! deallocation of memory
   deallocate(eta,u_init,v_init,d_init,t_init)

   end subroutine readi

#else ! READ_IC

!----------------------------------------------------------------------------------------------------------------------------------!
! second subroutine: Reads initial condition from continuation file.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine readi

! local variables
   integer                                        :: i, j, k

!----- Read file (IOS) ------------------------------------------------------------------------------------------------------------!
! read file info
   call mkfname(CFNAME, fbase, imach, err)
   call readcd(ftape, funit, fbase, imach, mt, mz, my, mx, nparam, times, info, ninf)

! read data
   call readd(ftape, uvel, 1, 1, err)
   call readd(ftape, vvel, 1, 2, err)
   call readd(ftape, temp, 1, 3, err)
   call readd(ftape, dens, 1, 4, err)
   call readd(ftape, pres, 1, 5, err)

   close(ftape)


! calculation of conservative variables
   call conserv

! fill up inflow arrays
   do j=1,ny
      u_infl(j) = uvel(1,j)
      v_infl(j) = vvel(1,j)
      t_infl(j) = temp(1,j)
      d_infl(j) = dens(1,j)
   end do
   uvel_ana(:,:) = uvel(:,:)
   vvel_ana(:,:) = vvel(:,:)
   temp_ana(:,:) = temp(:,:)
   dens_ana(:,:) = dens(:,:)
   pres_ana(:,:) = pres(:,:)


   end subroutine readi
#endif ! READ_IC

!----------------------------------------------------------------------------------------------------------------------------------!
! third subroutine: Calculates the initial condition for a compressible couette flow.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine couettei_calc

! local variables
   integer                                        :: i, j

!----- Calc couette flow ----------------------------------------------------------------------------------------------------------!
   do i=1,nx
      do j=1,ny
         uvel_ana(i,j) = pgrid_y(j)
         vvel_ana(i,j) = 0.0d0
         temp_ana(i,j) = 1.0d0+0.5d0*g1m2p*(1.0d0-pgrid_y(j)*pgrid_y(j))
         dens_ana(i,j) = 1.0d0/temp_ana(i,j)
         pres_ana(i,j) = igm2
      end do
   end do

   uvel(:,:) = uvel_ana(:,:)
   vvel(:,:) = vvel_ana(:,:)
   temp(:,:) = temp_ana(:,:)
   dens(:,:) = dens_ana(:,:)
   pres(:,:) = pres_ana(:,:)

!   uvel(:,:) = 0.0d0
!   vvel(:,:) = 0.0d0
!   temp(:,:) = 1.0d0
!   dens(:,:) = 1.0d0
!   pres(:,:) = igm2

! calculation of conservative variables
   call conserv

   end subroutine couettei_calc

!----------------------------------------------------------------------------------------------------------------------------------!
! fourth subroutine: Writes out data from last timestep.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine writo

!----- Write file (IOS) -----------------------------------------------------------------------------------------------------------!
! configure cd file.
   fbase   = WFNAME
   imach   = 2
   mt      = 1
   mx      = INT(nx, KIND=4)
   my      = INT(ny, KIND=4)
   mz      = 1
   nparam  = 5
   ninf    = 0
   info(1) = 'u-velocity'
   info(2) = 'v-velocity'
   info(3) = 'temperature'
   info(4) = 'density'
   info(5) = 'pressure'
!   info(6) = 'debug1'
!   info(7) = 'debug2'
!   info(8) = 'debug3'
!   info(9) = 'debug4'
!   info(10)= 'debug5'
!   info(11)= 'debug6'
   times(1)= 1

! write cd file
   call writecd(ftape, funit, fbase, imach, mt, mz, my, mx, nparam, times, info, ninf)

! write data
   call writed(ftape, uvel)
   call writed(ftape, vvel)
   call writed(ftape, temp)
   call writed(ftape, dens)
   call writed(ftape, pres)
!   call writed(ftape, debug1)
!   call writed(ftape, debug2)
!   call writed(ftape, debug3)
!   call writed(ftape, debug4)
!   call writed(ftape, debug5)
!   call writed(ftape, debug6)

! close the file
   close(ftape)

   end subroutine writo

!----------------------------------------------------------------------------------------------------------------------------------!
! fifth subroutine: Initializes output for intermediate steps.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine init_inter

! local variables
   integer                                        :: i

! configure cd file.
   fbase   = IFNAME
   imach   = 2
   mt      = INT(t_en/t_inter+1, KIND=4)
   mx      = INT(nx, KIND=4)
   my      = INT(ny, KIND=4)
   mz      = 1
   nparam  = 5
   ninf    = 0
   info(1) = 'u-velocity'
   info(2) = 'v-velocity'
   info(3) = 'density'
   info(4) = 'temperature'
   info(5) = 'pressure'
!   info(6) = 'debug1'
!   info(7) = 'debug2'
!   info(8) = 'debug3'
   do i=1,mt
     times(i)= INT(i*t_inter, KIND=4)
   end do

! write cd file
   call writecd(ftape, funit, fbase, imach, mt, mz, my, mx, nparam, times, info, ninf)

   end subroutine init_inter

!----------------------------------------------------------------------------------------------------------------------------------!
! sixth subroutine: Writes out intermediate steps.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine write_inter

! write data
   call writed(ftape, dens*uvel-dens_ana*uvel_ana)
   call writed(ftape, vvel-vvel_ana)
   call writed(ftape, dens-dens_ana)
   call writed(ftape, temp-temp_ana)
   call writed(ftape, pres-pres_ana)
!   call writed(ftape, superandi1-dens_ana)
!   call writed(ftape, superandi2-dens_ana)
!    call writed(ftape, superandi1+superandi2)
  
!   call writed(ftape, dens*uvel)
!   call writed(ftape, vvel)
!   call writed(ftape, dens)
!   call writed(ftape, temp)
!   call writed(ftape, pres)

   end subroutine write_inter

!----------------------------------------------------------------------------------------------------------------------------------!
! seventh subroutine: Closes file for intermediate steps.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine close_inter

! close the file
   close(ftape)

   end subroutine close_inter

!----------------------------------------------------------------------------------------------------------------------------------!
! eighth subroutine: Writes grid.
!----------------------------------------------------------------------------------------------------------------------------------!
   subroutine write_grid

! local variables   
   real, dimension(nx,ny,1)                       :: array_grid_x, array_grid_y, array_grid_z
   integer                                        :: i, j

! configure cd file.
   fbase   = GFNAME
   imach   = 2
   mt      = 1
   mx      = INT(nx, KIND=4)
   my      = INT(ny, KIND=4)
   mz      = 1
   nparam  = 3
   ninf    = 1
   info(1) = 'x-grid'
   info(2) = 'y-grid'
   info(3) = 'z-grid'
   info(4) = 'grid file for x3dp'
   times(1)= 1

! create arrays for grid file
   do i=1,nx
     do j=1,ny
      array_grid_x(i,j,1) = pgrid_x(i)
      array_grid_y(i,j,1) = pgrid_y(j)
      array_grid_z(i,j,1) = 1.0d0
     end do
   end do

! write cd file
   call writecd(ftape, funit, fbase, imach, mt, mz, my, mx, nparam, times, info, ninf)

! write out grid
   call writed(ftape,array_grid_x)
   call writed(ftape,array_grid_y)
   call writed(ftape,array_grid_z)
   
! close grid file   
   close(ftape)
   
   end subroutine write_grid

end module inp_outp




