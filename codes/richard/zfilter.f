c
c Version 2.0, 19 June 2003
c
c--------1---------2---------3---------4---------5---------6---------7--
c
      subroutine zfilter_init
c       Z filter routines
c       can't use the compact filter for the set-up of the domains
c       or domain-decomposition
c       use higher order explicit filter
c
c--------1---------2---------3---------4---------5---------6---------7--

      implicit none
      include 'param.h'
      include 'const.h'
      include 'filter.h'
      include 'knobs.h'
      include 'coord.h'
      
      integer i,j,k,ifil,nzo1,nzo2,nn1,nn2
      real*8 AA(7,7),bb(-3:3),h1(4:nzp-3,-3:3),h2(4:nz2tot-3,-3:3)
      real*8 AA1(7,7),bb1(1:7),hb(1:6,1:7)
      
      if(I_BL.ne.1) then
      nzo1=9  !total number of points that are filtered
      nzo2=7  !length of ramp
      nn1=nz2tot-nzo1
      nn2=nn1+nzo2

      
      do j=1,nn1-1
         fout_z(j) = 1.d0     
      end do
      
      do j=nn1,nn2
         fout_z(j) = 0.5d0*(1.d0+(cos(pi*(z2(j)-z2(nn1))
     &	             /(z2(nn2)-z2(nn1)))))     
      end do
      
      do j=nn2+1,nz2tot
         fout_z(j) = 0.d0
      end do
      
      open(89,file='z_out_filter.dat',status='unknown')
          do j=1,nz2tot
	     write(89,*) z2(j),fout_z(j)	 
          end do
      close(89)	  
      
crs   set up ramping of higher modes at outflow
      ifil=30

      do j=1,nz2tot-(ifil+1)
         fout_mod(j) = 1.d0     
      end do
      
      do j=nz2tot-ifil,nz2tot
         fout_mod(j) = (0.5d0*(1.d0+(dcos(pi*(dble(j-(nz2tot-ifil)))
     &	              /dble(ifil)))))**(2.d-1)
      end do
      

crs  set up ramping of higher modes in  spectral direction 
crs  if underresolved in theta

crs      do k=1,kh
crs      kscale(k) = (0.5d0*(1.d0+dcos(dble(k-1)*pi/dble(kh+1))))**dkn2
crs      write(*,*) 'KSCALE', kscale(k)
crs      end do

crs      if(isw15.eq.0) then
crs         do k=1,kh
crs            kscale(k) = 1.0d0	    
crs            write(*,*) 'KSCALE', kscale(k)
crs         end do
crs      else
crs         do k=1,isw15
crs            kscale(k) = 1.0d0	    
crs            write(*,*) 'KSCALE', kscale(k)
crs         end do
crs	 
crs	 do k=isw15,kh
crs            kscale(k) = (0.5d0*(1.d0+
crs     &	    dcos(dble(k-isw15)*pi/dble(kh-isw15+1))))	    
crs            write(*,*) 'KSCALE', kscale(k)
crs         end do
crs      end if	 
	 
      end if
      
      
      open(45,file='zfilter_weights.dat',status='unknown')

crs   ******************************************************************
crs   start with upper domain (for BL-comp, length nzp=nz1tot, 
crs   else nzp=nz2tot)
crs   ******************************************************************      

crs    Boundary points (inflow)
         j=1
	   do k=0,6
	      hb(j,k+1) = -(z1(j)-z1(j+k))	  	  
	   end do
	 j=2
	   hb(j,1) = -(z1(j)-z1(j-1))
	   do k=0,5
	       hb(j,k+2) = -(z1(j)-z1(j+k))	    
   	   end do
	 j=3
	   hb(j,1) = -(z1(j)-z1(j-2))
	   hb(j,2) = -(z1(j)-z1(j-1))
	   do k=0,4
	       hb(j,k+3) = -(z1(j)-z1(j+k))	    
   	   end do 

crscrs   ------------------ boundary points -> j=1 ------------------------
crs   set up rhs
      bb1(1)=1.d0
      bb1(2)=0.d0
      bb1(3)=0.d0
      bb1(4)=0.d0
      bb1(5)=0.d0
      bb1(6)=0.d0
      bb1(7)=0.d0
      
      j=1
crs   set up matrix
      AA1(1,1) = 1.d0
      AA1(1,2) = 1.d0
      AA1(1,3) = 1.d0
      AA1(1,4) = 1.d0
      AA1(1,5) = 1.d0
      AA1(1,6) = 1.d0
      AA1(1,7) = 1.d0
      
      AA1(2,1) = 0.d0
      AA1(2,2) = hb(j,2)
      AA1(2,3) = hb(j,3)
      AA1(2,4) = hb(j,4)
      AA1(2,5) = hb(j,5)
      AA1(2,6) = hb(j,6)
      AA1(2,7) = hb(j,7)
      
      AA1(3,1) = 0.d0
      AA1(3,2) = hb(j,2)**2
      AA1(3,3) = hb(j,3)**2
      AA1(3,4) = hb(j,4)**2
      AA1(3,5) = hb(j,5)**2
      AA1(3,6) = hb(j,6)**2
      AA1(3,7) = hb(j,7)**2
                  
      AA1(4,1) = 0.d0
      AA1(4,2) = hb(j,2)**3
      AA1(4,3) = hb(j,3)**3
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = hb(j,5)**3
      AA1(4,6) = hb(j,6)**3
      AA1(4,7) = hb(j,7)**3
      
      AA1(5,1) = 0.d0
      AA1(5,2) = hb(j,2)**4
      AA1(5,3) = hb(j,3)**4
      AA1(5,4) = hb(j,4)**4
      AA1(5,5) = hb(j,5)**4
      AA1(5,6) = hb(j,6)**4
      AA1(5,7) = hb(j,7)**4
      
      AA1(6,1) = 0.d0
      AA1(6,2) = hb(j,2)**5
      AA1(6,3) = hb(j,3)**5
      AA1(6,4) = hb(j,4)**5
      AA1(6,5) = hb(j,5)**5
      AA1(6,6) = hb(j,6)**5
      AA1(6,7) = hb(j,7)**5
      
      AA1(7,1) = 1.d0
      AA1(7,2) = -1.d0
      AA1(7,3) = 1.d0
      AA1(7,4) = -1.d0
      AA1(7,5) = 1.d0
      AA1(7,6) = -1.d0
      AA1(7,7) = 1.d0
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      do k=1,7
	   zf1b(j,k) = bb1(k)
      end do
      
      write(45,*) 'coefficients for j=1 are (domain 1)'
      do i=1,7
         write(45,*) zf1b(j,i)
      end do
crs
crs
crscrs   ------------------ boundary points -> j=2 ------------------------
crs   set up rhs
      bb1(1)=0.d0
      bb1(2)=1.d0
      bb1(3)=0.d0
      bb1(4)=0.d0
      bb1(5)=0.d0
      bb1(6)=0.d0
      bb1(7)=0.d0
      
      j=2
crs   set up matrix
      AA1(1,1) = hb(j,1)
      AA1(1,2) = 0.0d0
      AA1(1,3) = hb(j,3)
      AA1(1,4) = hb(j,4)
      AA1(1,5) = hb(j,5)
      AA1(1,6) = hb(j,6)
      AA1(1,7) = hb(j,7)
      
      AA1(2,1) = 1.d0
      AA1(2,2) = 1.d0
      AA1(2,3) = 1.d0
      AA1(2,4) = 1.d0
      AA1(2,5) = 1.d0
      AA1(2,6) = 1.d0
      AA1(2,7) = 1.d0
      
      AA1(3,1) = hb(j,1)**2
      AA1(3,2) = 0.0d0
      AA1(3,3) = hb(j,3)**2
      AA1(3,4) = hb(j,4)**2
      AA1(3,5) = hb(j,5)**2
      AA1(3,6) = hb(j,6)**2
      AA1(3,7) = hb(j,7)**2
                  
      AA1(4,1) = hb(j,1)**3
      AA1(4,2) = 0.0d0
      AA1(4,3) = hb(j,3)**3
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = hb(j,5)**3
      AA1(4,6) = hb(j,6)**3
      AA1(4,7) = hb(j,7)**3
      
      AA1(5,1) = hb(j,1)**4
      AA1(5,2) = 0.0d0
      AA1(5,3) = hb(j,3)**4
      AA1(5,4) = hb(j,4)**4
      AA1(5,5) = hb(j,5)**4
      AA1(5,6) = hb(j,6)**4
      AA1(5,7) = hb(j,7)**4
      
      AA1(6,1) = hb(j,1)**5
      AA1(6,2) = 0.0d0
      AA1(6,3) = hb(j,3)**5
      AA1(6,4) = hb(j,4)**5
      AA1(6,5) = hb(j,5)**5
      AA1(6,6) = hb(j,6)**5
      AA1(6,7) = hb(j,7)**5
      
      AA1(7,1) = -1.d0
      AA1(7,2) = 1.d0
      AA1(7,3) = -1.d0
      AA1(7,4) = 1.d0
      AA1(7,5) = -1.d0
      AA1(7,6) = 1.d0
      AA1(7,7) = -1.d0
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      do k=1,7
	   zf1b(j,k) = bb1(k)
      end do
      
      write(45,*) 'coefficients for j=2 are (domain 1)'
      do i=1,7
         write(45,*) zf1b(j,i)
      end do

crs
crscrs   ------------------ boundary points -> j=3 ------------------------
crs   set up rhs
      bb1(1)=0.d0
      bb1(2)=0.d0
      bb1(3)=1.d0
      bb1(4)=0.d0
      bb1(5)=0.d0
      bb1(6)=0.d0
      bb1(7)=0.d0
      
      j=3
crs   set up matrix
      AA1(1,1) = hb(j,1)**2
      AA1(1,2) = hb(j,2)**2
      AA1(1,3) = 0.d0
      AA1(1,4) = hb(j,4)**2
      AA1(1,5) = hb(j,5)**2
      AA1(1,6) = hb(j,6)**2
      AA1(1,7) = hb(j,7)**2
                  
      AA1(2,1) = hb(j,1)
      AA1(2,2) = hb(j,2)
      AA1(2,3) = 0.d0
      AA1(2,4) = hb(j,4)
      AA1(2,5) = hb(j,5)
      AA1(2,6) = hb(j,6)
      AA1(2,7) = hb(j,7)
      
      AA1(3,1) = 1.d0
      AA1(3,2) = 1.d0
      AA1(3,3) = 1.d0
      AA1(3,4) = 1.d0
      AA1(3,5) = 1.d0
      AA1(3,6) = 1.d0
      AA1(3,7) = 1.d0
      
      AA1(4,1) = hb(j,1)**3
      AA1(4,2) = hb(j,2)**3
      AA1(4,3) = 0.d0
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = hb(j,5)**3
      AA1(4,6) = hb(j,6)**3
      AA1(4,7) = hb(j,7)**3
      
      AA1(5,1) = hb(j,1)**4
      AA1(5,2) = hb(j,2)**4
      AA1(5,3) = 0.d0
      AA1(5,4) = hb(j,4)**4
      AA1(5,5) = hb(j,5)**4
      AA1(5,6) = hb(j,6)**4
      AA1(5,7) = hb(j,7)**4
      
      AA1(6,1) = hb(j,1)**5
      AA1(6,2) = hb(j,2)**5
      AA1(6,3) = 0.d0
      AA1(6,4) = hb(j,4)**5
      AA1(6,5) = hb(j,5)**5
      AA1(6,6) = hb(j,6)**5
      AA1(6,7) = hb(j,7)**5
      
      AA1(7,1) = 1.d0
      AA1(7,2) = -1.d0
      AA1(7,3) = 1.d0
      AA1(7,4) = -1.d0
      AA1(7,5) = 1.d0
      AA1(7,6) = -1.d0
      AA1(7,7) = 1.d0
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      do k=1,7
	   zf1b(j,k) = bb1(k)
      end do
      
      write(45,*) 'coefficients for j=3 are (domain 1)'
      do i=1,7
         write(45,*) zf1b(j,i)
      end do
      
crs   ------------------------ interior --------------------------------

        do j=4,nzp-3	     
	   do k=-3,3
	      h1(j,k) = -(zall(j)-zall(j+k))		
	   end do
	end do


      do j=4,nzp-3
      
crs   set up rhs
         bb(-3)= 0.d0
         bb(-2)= 0.d0
         bb(-1)= 0.d0
         bb(0) = 1.d0
         bb(1) = 0.d0
         bb(2) = 0.d0
         bb(3) = 0.d0

         
crs   set up matrix
         AA(1,1) = h1(j,-3)**5
         AA(1,2) = h1(j,-2)**5
         AA(1,3) = h1(j,-1)**5
         AA(1,4) = 0.d0
         AA(1,5) = h1(j,1)**5
         AA(1,6) = h1(j,2)**5
         AA(1,7) = h1(j,3)**5
	 
         AA(2,1) = h1(j,-3)**4
         AA(2,2) = h1(j,-2)**4
         AA(2,3) = h1(j,-1)**4
         AA(2,4) = 0.d0
         AA(2,5) = h1(j,1)**4
         AA(2,6) = h1(j,2)**4
         AA(2,7) = h1(j,3)**4
	 
         AA(3,1) = h1(j,-3)**3
         AA(3,2) = h1(j,-2)**3
         AA(3,3) = h1(j,-1)**3
         AA(3,4) = 0.0d0
         AA(3,5) = h1(j,1)**3
         AA(3,6) = h1(j,2)**3
         AA(3,7) = h1(j,3)**3
      
         AA(4,1) = 1.d0
         AA(4,2) = 1.d0
         AA(4,3) = 1.d0
         AA(4,4) = 1.d0
         AA(4,5) = 1.d0
         AA(4,6) = 1.d0
         AA(4,7) = 1.d0
      
         AA(5,1) = h1(j,-3)
         AA(5,2) = h1(j,-2)
         AA(5,3) = h1(j,-1)
         AA(5,4) = 0.0d0
         AA(5,5) = h1(j,1)
         AA(5,6) = h1(j,2)
         AA(5,7) = h1(j,3)
      
         AA(6,1) = h1(j,-3)**2
         AA(6,2) = h1(j,-2)**2
         AA(6,3) = h1(j,-1)**2
         AA(6,4) = 0.0d0
         AA(6,5) = h1(j,1)**2
         AA(6,6) = h1(j,2)**2
         AA(6,7) = h1(j,3)**2
            
         AA(7,1) = -1.d0
         AA(7,2) = 1.d0
         AA(7,3) = -1.d0
         AA(7,4) = 1.d0
         AA(7,5) = -1.d0
         AA(7,6) = 1.d0
         AA(7,7) = -1.d0
      
           
         call LinSolve( AA, bb, 7, 7)
	 
	 do k=-3,3
	    zf1(j,k) = bb(k)
         end do
	 
         write(45,*) 'coefficients for domain 1 are'
         do k=-3,3
            write(45,*) j,k,zall(j),zf1(j,k)
         end do
      
      end do

crs   outflow (either B-L or wake)
      
      	 j=4
	   hb(j,1) = -(zall(nzp-2)-zall(nzp-6))
	   hb(j,2) = -(zall(nzp-2)-zall(nzp-5))
	   hb(j,3) = -(zall(nzp-2)-zall(nzp-4))
	   hb(j,4) = -(zall(nzp-2)-zall(nzp-3))
	   hb(j,6) = -(zall(nzp-2)-zall(nzp-1))
	   hb(j,7) = -(zall(nzp-2)-zall(nzp))
	   
	   
	 j=5
	   hb(j,1) = -(zall(nzp-1)-zall(nzp-6))
	   hb(j,2) = -(zall(nzp-1)-zall(nzp-5))
	   hb(j,3) = -(zall(nzp-1)-zall(nzp-4))
	   hb(j,4) = -(zall(nzp-1)-zall(nzp-3))
	   hb(j,5) = -(zall(nzp-1)-zall(nzp-2))
	   hb(j,7) = -(zall(nzp-1)-zall(nzp))

         j=6
	   hb(j,1) = -(zall(nzp)-zall(nzp-6))
	   hb(j,2) = -(zall(nzp)-zall(nzp-5))
	   hb(j,3) = -(zall(nzp)-zall(nzp-4))
	   hb(j,4) = -(zall(nzp)-zall(nzp-3))
	   hb(j,5) = -(zall(nzp)-zall(nzp-2))
	   hb(j,6) = -(zall(nzp)-zall(nzp-1)) 


crscrs   ------------------ boundary points -> j=nz2tot-2 --------------------
crscrs   set up rhs
crs   set up rhs
      bb1(1)=0.d0
      bb1(2)=0.d0
      bb1(3)=0.d0
      bb1(4)=0.d0
      bb1(5)=1.d0
      bb1(6)=0.d0
      bb1(7)=0.d0
      
      j=4
crs   set up matrix
      AA1(1,1) = 1.d0
      AA1(1,2) = -1.d0
      AA1(1,3) = 1.d0
      AA1(1,4) = -1.d0
      AA1(1,5) = 1.d0
      AA1(1,6) = -1.d0
      AA1(1,7) = 1.d0
            
      AA1(2,1) = hb(j,1)**5
      AA1(2,2) = hb(j,2)**5
      AA1(2,3) = hb(j,3)**5
      AA1(2,4) = hb(j,4)**5
      AA1(2,5) = 0.d0
      AA1(2,6) = hb(j,6)**5
      AA1(2,7) = hb(j,7)**5
      
      AA1(3,1) = hb(j,1)**4
      AA1(3,2) = hb(j,2)**4
      AA1(3,3) = hb(j,3)**4
      AA1(3,4) = hb(j,4)**4
      AA1(3,5) = 0.d0
      AA1(3,6) = hb(j,6)**4
      AA1(3,7) = hb(j,7)**4
      
      AA1(4,1) = hb(j,1)**3
      AA1(4,2) = hb(j,2)**3
      AA1(4,3) = hb(j,3)**3
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = 0.d0
      AA1(4,6) = hb(j,6)**3
      AA1(4,7) = hb(j,7)**3
      
      AA1(5,1) = 1.d0
      AA1(5,2) = 1.d0
      AA1(5,3) = 1.d0
      AA1(5,4) = 1.d0
      AA1(5,5) = 1.d0
      AA1(5,6) = 1.d0
      AA1(5,7) = 1.d0
      
      AA1(6,1) = hb(j,1)
      AA1(6,2) = hb(j,2)
      AA1(6,3) = hb(j,3)
      AA1(6,4) = hb(j,4)
      AA1(6,5) = 0.d0
      AA1(6,6) = hb(j,6)
      AA1(6,7) = hb(j,7)
                  
      AA1(7,1) = hb(j,1)**2
      AA1(7,2) = hb(j,2)**2
      AA1(7,3) = hb(j,3)**2
      AA1(7,4) = hb(j,4)**2
      AA1(7,5) = 0.0d0
      AA1(7,6) = hb(j,6)**2
      AA1(7,7) = hb(j,7)**2
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      if (I_BL.ne.1) then
         do k=1,7
	     zf2b(j,k) = bb1(k)
         end do
      else
         do k=1,7
	     zf1b(j,k) = bb1(k)
         end do
      endif	 	 
      
      write(45,*) 'coefficients for outflow-2'
      if (I_BL.ne.1) then
         do i=1,7
            write(45,*) zf2b(j,i)
         end do
      else
         do i=1,7
            write(45,*) zf1b(j,i)
         end do
      end if

crscrs   ------------------ boundary points -> j=nz2tot-1 --------------------
crs   set up rhs
      bb1(1)=0.d0
      bb1(2)=0.d0
      bb1(3)=0.d0      
      bb1(4)=0.d0
      bb1(5)=0.d0
      bb1(6)=1.d0
      bb1(7)=0.d0
      
      j=5
crs   set up matrix
      AA1(1,1) = -1.d0
      AA1(1,2) = 1.d0
      AA1(1,3) = -1.d0
      AA1(1,4) = 1.d0
      AA1(1,5) = -1.d0
      AA1(1,6) = 1.d0
      AA1(1,7) = -1.d0
      
      AA1(2,1) = hb(j,1)**5
      AA1(2,2) = hb(j,2)**5
      AA1(2,3) = hb(j,3)**5
      AA1(2,4) = hb(j,4)**5
      AA1(2,5) = hb(j,5)**5
      AA1(2,6) = 0.0d0
      AA1(2,7) = hb(j,7)**5
            
      AA1(3,1) = hb(j,1)**4
      AA1(3,2) = hb(j,2)**4
      AA1(3,3) = hb(j,3)**4
      AA1(3,4) = hb(j,4)**4
      AA1(3,5) = hb(j,5)**4
      AA1(3,6) = 0.0d0
      AA1(3,7) = hb(j,7)**4
      
      AA1(4,1) = hb(j,1)**3
      AA1(4,2) = hb(j,2)**3
      AA1(4,3) = hb(j,3)**3
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = hb(j,5)**3
      AA1(4,6) = 0.0d0
      AA1(4,7) = hb(j,7)**3
      
      AA1(5,1) = hb(j,1)**2
      AA1(5,2) = hb(j,2)**2
      AA1(5,3) = hb(j,3)**2
      AA1(5,4) = hb(j,4)**2
      AA1(5,5) = hb(j,5)**2
      AA1(5,6) = 0.0d0
      AA1(5,7) = hb(j,7)**2
                  
      AA1(6,1) = 1.d0
      AA1(6,2) = 1.d0
      AA1(6,3) = 1.d0
      AA1(6,4) = 1.d0
      AA1(6,5) = 1.d0
      AA1(6,6) = 1.d0
      AA1(6,7) = 1.d0
      
      AA1(7,1) = hb(j,1)
      AA1(7,2) = hb(j,2)
      AA1(7,3) = hb(j,3)
      AA1(7,4) = hb(j,4)
      AA1(7,5) = hb(j,5)
      AA1(7,6) = 0.0d0
      AA1(7,7) = hb(j,7)
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      if (I_BL.ne.1) then
         do k=1,7
	     zf2b(j,k) = bb1(k)
         end do
      else
         do k=1,7
	     zf1b(j,k) = bb1(k)
         end do
      endif	 	 
      
      write(45,*) 'coefficients for outflow-1'
      if (I_BL.ne.1) then
         do i=1,7
            write(45,*) zf2b(j,i)
         end do
      else
         do i=1,7
            write(45,*) zf1b(j,i)
         end do
      end if


crscrs   ------------------ boundary points -> j=nz2tot --------------------
crs   set up rhs
      bb1(1)=0.d0
      bb1(2)=0.d0
      bb1(3)=0.d0
      bb1(4)=0.d0
      bb1(5)=0.d0
      bb1(6)=0.d0
      bb1(7)=1.d0
      
      j=6
crs   set up matrix
      AA1(1,1) = 1.d0
      AA1(1,2) = -1.d0
      AA1(1,3) = 1.d0
      AA1(1,4) = -1.d0
      AA1(1,5) = 1.d0
      AA1(1,6) = -1.d0
      AA1(1,7) = 1.d0
      
      AA1(2,1) = hb(j,1)**5
      AA1(2,2) = hb(j,2)**5
      AA1(2,3) = hb(j,3)**5
      AA1(2,4) = hb(j,4)**5
      AA1(2,5) = hb(j,5)**5
      AA1(2,6) = hb(j,6)**5
      AA1(2,7) = 0.0d0
      
      AA1(3,1) = hb(j,1)**4
      AA1(3,2) = hb(j,2)**4
      AA1(3,3) = hb(j,3)**4
      AA1(3,4) = hb(j,4)**4
      AA1(3,5) = hb(j,5)**4
      AA1(3,6) = hb(j,6)**4
      AA1(3,7) = 0.0d0
      
      AA1(4,1) = hb(j,1)**3
      AA1(4,2) = hb(j,2)**3
      AA1(4,3) = hb(j,3)**3
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = hb(j,5)**3
      AA1(4,6) = hb(j,6)**3
      AA1(4,7) = 0.0d0
      
      AA1(5,1) = hb(j,1)**2
      AA1(5,2) = hb(j,2)**2
      AA1(5,3) = hb(j,3)**2
      AA1(5,4) = hb(j,4)**2
      AA1(5,5) = hb(j,5)**2
      AA1(5,6) = hb(j,6)**2
      AA1(5,7) = 0.0d0
                  
      AA1(6,1) = hb(j,1)
      AA1(6,2) = hb(j,2)
      AA1(6,3) = hb(j,3)
      AA1(6,4) = hb(j,4)
      AA1(6,5) = hb(j,5)
      AA1(6,6) = hb(j,6)
      AA1(6,7) = 0.0d0
      
      AA1(7,1) = 1.d0
      AA1(7,2) = 1.d0
      AA1(7,3) = 1.d0
      AA1(7,4) = 1.d0
      AA1(7,5) = 1.d0
      AA1(7,6) = 1.d0
      AA1(7,7) = 1.d0
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      if (I_BL.ne.1) then
         do k=1,7
	     zf2b(j,k) = bb1(k)
         end do
      else
         do k=1,7
	     zf1b(j,k) = bb1(k)
         end do
      endif	 	 
      
      write(45,*) 'coefficients for outflow'
      if (I_BL.ne.1) then
         do i=1,7
            write(45,*) zf2b(j,i)
         end do
      else
         do i=1,7
            write(45,*) zf1b(j,i)
         end do
      end if


      
      if(I_BL.eq.0) then 
crs   ******************************************************************
crs    WAKE REGION
crs   ******************************************************************      

crs    Boundary points (base)
         j=1
	   do k=0,6
	      hb(j,k+1) = -(z2(j)-z2(j+k))	  	  
	   end do
	 j=2
	   hb(j,1) = -(z2(j)-z2(j-1))
	   do k=0,5
	       hb(j,k+2) = -(z2(j)-z2(j+k))	    
   	   end do
	 j=3
	   hb(j,1) = -(z2(j)-z2(j-2))
	   hb(j,2) = -(z2(j)-z2(j-1))
	   do k=0,4
	       hb(j,k+3) = -(z2(j)-z2(j+k))	    
   	   end do 

crscrs   ------------------ boundary points -> j=1 ------------------------
crs   set up rhs
      bb1(1)=1.d0
      bb1(2)=0.d0
      bb1(3)=0.d0
      bb1(4)=0.d0
      bb1(5)=0.d0
      bb1(6)=0.d0
      bb1(7)=0.d0
      
      j=1
crs   set up matrix
      AA1(1,1) = 1.d0
      AA1(1,2) = 1.d0
      AA1(1,3) = 1.d0
      AA1(1,4) = 1.d0
      AA1(1,5) = 1.d0
      AA1(1,6) = 1.d0
      AA1(1,7) = 1.d0
      
      AA1(2,1) = 0.d0
      AA1(2,2) = hb(j,2)
      AA1(2,3) = hb(j,3)
      AA1(2,4) = hb(j,4)
      AA1(2,5) = hb(j,5)
      AA1(2,6) = hb(j,6)
      AA1(2,7) = hb(j,7)
      
      AA1(3,1) = 0.d0
      AA1(3,2) = hb(j,2)**2
      AA1(3,3) = hb(j,3)**2
      AA1(3,4) = hb(j,4)**2
      AA1(3,5) = hb(j,5)**2
      AA1(3,6) = hb(j,6)**2
      AA1(3,7) = hb(j,7)**2
                  
      AA1(4,1) = 0.d0
      AA1(4,2) = hb(j,2)**3
      AA1(4,3) = hb(j,3)**3
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = hb(j,5)**3
      AA1(4,6) = hb(j,6)**3
      AA1(4,7) = hb(j,7)**3
      
      AA1(5,1) = 0.d0
      AA1(5,2) = hb(j,2)**4
      AA1(5,3) = hb(j,3)**4
      AA1(5,4) = hb(j,4)**4
      AA1(5,5) = hb(j,5)**4
      AA1(5,6) = hb(j,6)**4
      AA1(5,7) = hb(j,7)**4
      
      AA1(6,1) = 0.d0
      AA1(6,2) = hb(j,2)**5
      AA1(6,3) = hb(j,3)**5
      AA1(6,4) = hb(j,4)**5
      AA1(6,5) = hb(j,5)**5
      AA1(6,6) = hb(j,6)**5
      AA1(6,7) = hb(j,7)**5
      
      AA1(7,1) = 1.d0
      AA1(7,2) = -1.d0
      AA1(7,3) = 1.d0
      AA1(7,4) = -1.d0
      AA1(7,5) = 1.d0
      AA1(7,6) = -1.d0
      AA1(7,7) = 1.d0
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      do k=1,7
	   zf2b(j,k) = bb1(k)
      end do
      
      write(45,*) 'coefficients for j=1 are (domain 2)'
      do i=1,7
         write(45,*) zf2b(j,i)
      end do
crs
crs
crscrs   ------------------ boundary points -> j=2 ------------------------
crs   set up rhs
      bb1(1)=0.d0
      bb1(2)=1.d0
      bb1(3)=0.d0
      bb1(4)=0.d0
      bb1(5)=0.d0
      bb1(6)=0.d0
      bb1(7)=0.d0
      
      j=2
crs   set up matrix
      AA1(1,1) = hb(j,1)
      AA1(1,2) = 0.0d0
      AA1(1,3) = hb(j,3)
      AA1(1,4) = hb(j,4)
      AA1(1,5) = hb(j,5)
      AA1(1,6) = hb(j,6)
      AA1(1,7) = hb(j,7)
      
      AA1(2,1) = 1.d0
      AA1(2,2) = 1.d0
      AA1(2,3) = 1.d0
      AA1(2,4) = 1.d0
      AA1(2,5) = 1.d0
      AA1(2,6) = 1.d0
      AA1(2,7) = 1.d0
      
      AA1(3,1) = hb(j,1)**2
      AA1(3,2) = 0.0d0
      AA1(3,3) = hb(j,3)**2
      AA1(3,4) = hb(j,4)**2
      AA1(3,5) = hb(j,5)**2
      AA1(3,6) = hb(j,6)**2
      AA1(3,7) = hb(j,7)**2
                  
      AA1(4,1) = hb(j,1)**3
      AA1(4,2) = 0.0d0
      AA1(4,3) = hb(j,3)**3
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = hb(j,5)**3
      AA1(4,6) = hb(j,6)**3
      AA1(4,7) = hb(j,7)**3
      
      AA1(5,1) = hb(j,1)**4
      AA1(5,2) = 0.0d0
      AA1(5,3) = hb(j,3)**4
      AA1(5,4) = hb(j,4)**4
      AA1(5,5) = hb(j,5)**4
      AA1(5,6) = hb(j,6)**4
      AA1(5,7) = hb(j,7)**4
      
      AA1(6,1) = hb(j,1)**5
      AA1(6,2) = 0.0d0
      AA1(6,3) = hb(j,3)**5
      AA1(6,4) = hb(j,4)**5
      AA1(6,5) = hb(j,5)**5
      AA1(6,6) = hb(j,6)**5
      AA1(6,7) = hb(j,7)**5
      
      AA1(7,1) = -1.d0
      AA1(7,2) = 1.d0
      AA1(7,3) = -1.d0
      AA1(7,4) = 1.d0
      AA1(7,5) = -1.d0
      AA1(7,6) = 1.d0
      AA1(7,7) = -1.d0
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      do k=1,7
	   zf2b(j,k) = bb1(k)
      end do
      
      write(45,*) 'coefficients for j=2 are (domain 2)'
      do i=1,7
         write(45,*) zf2b(j,i)
      end do

crs
crscrs   ------------------ boundary points -> j=3 ------------------------
crs   set up rhs
      bb1(1)=0.d0
      bb1(2)=0.d0
      bb1(3)=1.d0
      bb1(4)=0.d0
      bb1(5)=0.d0
      bb1(6)=0.d0
      bb1(7)=0.d0
      
      j=3
crs   set up matrix
      AA1(1,1) = hb(j,1)**2
      AA1(1,2) = hb(j,2)**2
      AA1(1,3) = 0.d0
      AA1(1,4) = hb(j,4)**2
      AA1(1,5) = hb(j,5)**2
      AA1(1,6) = hb(j,6)**2
      AA1(1,7) = hb(j,7)**2
                  
      AA1(2,1) = hb(j,1)
      AA1(2,2) = hb(j,2)
      AA1(2,3) = 0.d0
      AA1(2,4) = hb(j,4)
      AA1(2,5) = hb(j,5)
      AA1(2,6) = hb(j,6)
      AA1(2,7) = hb(j,7)
      
      AA1(3,1) = 1.d0
      AA1(3,2) = 1.d0
      AA1(3,3) = 1.d0
      AA1(3,4) = 1.d0
      AA1(3,5) = 1.d0
      AA1(3,6) = 1.d0
      AA1(3,7) = 1.d0
      
      AA1(4,1) = hb(j,1)**3
      AA1(4,2) = hb(j,2)**3
      AA1(4,3) = 0.d0
      AA1(4,4) = hb(j,4)**3
      AA1(4,5) = hb(j,5)**3
      AA1(4,6) = hb(j,6)**3
      AA1(4,7) = hb(j,7)**3
      
      AA1(5,1) = hb(j,1)**4
      AA1(5,2) = hb(j,2)**4
      AA1(5,3) = 0.d0
      AA1(5,4) = hb(j,4)**4
      AA1(5,5) = hb(j,5)**4
      AA1(5,6) = hb(j,6)**4
      AA1(5,7) = hb(j,7)**4
      
      AA1(6,1) = hb(j,1)**5
      AA1(6,2) = hb(j,2)**5
      AA1(6,3) = 0.d0
      AA1(6,4) = hb(j,4)**5
      AA1(6,5) = hb(j,5)**5
      AA1(6,6) = hb(j,6)**5
      AA1(6,7) = hb(j,7)**5
      
      AA1(7,1) = 1.d0
      AA1(7,2) = -1.d0
      AA1(7,3) = 1.d0
      AA1(7,4) = -1.d0
      AA1(7,5) = 1.d0
      AA1(7,6) = -1.d0
      AA1(7,7) = 1.d0
      
           
      call LinSolve( AA1, bb1, 7, 7)
      
      do k=1,7
	   zf2b(j,k) = bb1(k)
      end do
      
      write(45,*) 'coefficients for j=3 are (domain 2)'
      do i=1,7
         write(45,*) zf2b(j,i)
      end do
      
         
      end if ! if(I_BL.eq.0)


c *** Report Setup ***

      write(*,*) '*** zfilter_init setup report ***********************'
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)


      return
      end

 
c--------1---------2---------3---------4---------5---------6---------7--
c
      subroutine zfilter_exp1(nd)
c
c     Filter in z in place for one r/z plane of wave space data.
c
c--------1---------2---------3---------4---------5---------6---------7--


      implicit none
      include 'param.h'
#      include "common.h"
      include 'dgrp.h'
      include 'filter.h'


c --- Arguments ---
       
      integer i,j,k,l,nd,jj
      
      if(ndom1.eq.1) then
      
      do l=1,neqn
           do k=0,kh
	     
	       j=1
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(1,1)*fitmp1(i,j,k,l,nd)
     &	                        + zf1b(1,2)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(1,3)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(1,4)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(1,5)*fitmp1(i,j+4,k,l,nd) 
     &	                        + zf1b(1,6)*fitmp1(i,j+5,k,l,nd) 
     &	                        + zf1b(1,7)*fitmp1(i,j+6,k,l,nd)
               end do
	     
	       j=2
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(2,1)*fitmp1(i,j-1,k,l,nd)
     &	                        + zf1b(2,2)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(2,3)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(2,4)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(2,5)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(2,6)*fitmp1(i,j+4,k,l,nd) 
     &	                        + zf1b(2,7)*fitmp1(i,j+5,k,l,nd)
               end do
	     
	       j=3
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(3,1)*fitmp1(i,j-2,k,l,nd)
     &	                        + zf1b(3,2)*fitmp1(i,j-1,k,l,nd) 
     &	                        + zf1b(3,3)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(3,4)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(3,5)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(3,6)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(3,7)*fitmp1(i,j+4,k,l,nd)
               end do
	       
	     
	     do j=4,nz1-3	         
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(j,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
	     
	     j=nz1-2
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(j,3)*fitmp2(i+mrd,1,k,l,1)
      
	     end do
	     
	     j=nz1-1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp2(i+mrd,1,k,l,1)
     &	                       + zf1(j,3)*fitmp2(i+mrd,2,k,l,1)
      
	     end do
	     
	     j=nz1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp2(i+mrd,1,k,l,1) 
     &	                       + zf1(j,2)*fitmp2(i+mrd,2,k,l,1)
     &	                       + zf1(j,3)*fitmp2(i+mrd,3,k,l,1)
      
	     end do
	     

	   end do
	   end do
      

      else !if ndom1.ne.1


             
       if(nd.eq.1) then !inflow
           do l=1,neqn
           do k=0,kh
	     
	       j=1
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(1,1)*fitmp1(i,j,k,l,nd)
     &	                        + zf1b(1,2)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(1,3)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(1,4)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(1,5)*fitmp1(i,j+4,k,l,nd) 
     &	                        + zf1b(1,6)*fitmp1(i,j+5,k,l,nd) 
     &	                        + zf1b(1,7)*fitmp1(i,j+6,k,l,nd)
               end do
	     
	       j=2
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(2,1)*fitmp1(i,j-1,k,l,nd)
     &	                        + zf1b(2,2)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(2,3)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(2,4)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(2,5)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(2,6)*fitmp1(i,j+4,k,l,nd) 
     &	                        + zf1b(2,7)*fitmp1(i,j+5,k,l,nd)
               end do
	     
	       j=3
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(3,1)*fitmp1(i,j-2,k,l,nd)
     &	                        + zf1b(3,2)*fitmp1(i,j-1,k,l,nd) 
     &	                        + zf1b(3,3)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(3,4)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(3,5)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(3,6)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(3,7)*fitmp1(i,j+4,k,l,nd)
               end do
	      	       
	     
	     do j=4,nz1-3	         
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(j,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
crs	     call writeit(3,1,mr1,nz1,1,cu12(1,1,0,1,nd)
crs     &	     ,fitmp1(1,1,0,1,nd),fitmp1(1,1,0,1,nd))
crs             stop 'z-filter'

	     
	     j=nz1-2
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(j,3)*fitmp1(i,1,k,l,nd+1)
      
	     end do
	     
	     j=nz1-1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,1,k,l,nd+1)
     &	                       + zf1(j,3)*fitmp1(i,2,k,l,nd+1)
      
	     end do
	     
	     j=nz1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,1,k,l,nd+1) 
     &	                       + zf1(j,2)*fitmp1(i,2,k,l,nd+1)
     &	                       + zf1(j,3)*fitmp1(i,3,k,l,nd+1)
      
	     end do
	     

	   end do
	   end do
       
       elseif(nd.eq.ndom1) then !last domain of region 1
                     
	   do l=1,neqn
           do k=0,kh
	     
	     j=1
	     jj=(nd-1)*nz1+1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1-2,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,nz1-1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     
	     j=2
	     jj=(nd-1)*nz1+2
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1-1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     
	     j=3
	     jj=(nd-1)*nz1+3
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd)
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd)
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	   
	     do j=4,nz1-3
	     jj=(nd-1)*nz1+j	         
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
	     j=nz1-2
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i+mrd,1,k,l,1)
      
	     end do
	     
	     j=nz1-1
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i+mrd,1,k,l,1)
     &	                       + zf1(jj,3)*fitmp2(i+mrd,2,k,l,1)
      
	     end do
	     
	     j=nz1
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i+mrd,1,k,l,1) 
     &	                       + zf1(jj,2)*fitmp2(i+mrd,2,k,l,1)
     &	                       + zf1(jj,3)*fitmp2(i+mrd,3,k,l,1)
      
	     end do
	     

	   end do
	   end do
	   
	     
	   
       else ! any interior domain
       
           do l=1,neqn
           do k=0,kh
	     
	     j=1
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1-2,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,nz1-1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     
	     j=2
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1-1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     
	     j=3
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	   
	     do j=4,nz1-3
	     jj=(nd-1)*nz1+j	         
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
	     j=nz1-2
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,1,k,l,nd+1)
      
	     end do
	     
	     j=nz1-1
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,1,k,l,nd+1)
     &	                       + zf1(jj,3)*fitmp1(i,2,k,l,nd+1)
      
	     end do
	     
	     j=nz1
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,1,k,l,nd+1) 
     &	                       + zf1(jj,2)*fitmp1(i,2,k,l,nd+1)
     &	                       + zf1(jj,3)*fitmp1(i,3,k,l,nd+1)
      
	     end do
	     

	   end do
	   end do
	   
       end if
       
       end if !if ndom1.eq.1	   	   
       
       return
       
       end	       	        


c--------1---------2---------3---------4---------5---------6---------7--
c
      subroutine zfilter_exp2(nd)
c
c     Filter in z in place for one r/z plane of wave space data.
c
c--------1---------2---------3---------4---------5---------6---------7--


      implicit none
      include 'param.h'
#      include "common.h"
      include 'filter.h'


c --- Arguments ---
       
      integer i,j,k,l,nd,jj
      
       
       if(nd.eq.1) then !interface with region 1
           do l=1,neqn
           do k=0,kh
	   
	     j=1
	     jj=nz1tot+j
	     do i=mrd1,mr2
	     cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i-mrd,nz1-2,k,l,ndom1) 
     &	                      + zf1(jj,-2)*fitmp1(i-mrd,nz1-1,k,l,ndom1) 
     &	                      + zf1(jj,-1)*fitmp1(i-mrd,nz1,k,l,ndom1) 
     &	                      + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                      + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                      + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                      + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     
	     j=2
	     jj=nz1tot+j
	     do i=mrd1,mr2
	     cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i-mrd,nz1-1,k,l,ndom1) 
     &	                      + zf1(jj,-2)*fitmp1(i-mrd,nz1,k,l,ndom1) 
     &	                      + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                      + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                      + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                      + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                      + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     
	     j=3
	     jj=nz1tot+j
	     do i=mrd1,mr2
	     cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i-mrd,nz1,k,l,ndom1) 
     &	                      + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                      + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                      + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                      + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                      + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                      + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     
crs      base -> one-sided stencils	     
crs	       j=1
crs	       do i=1,mrd	     
crs	       cu22(i,j,k,l,nd) = zf2b(1,1)*fitmp2(i,j,k,l,nd)
crs     &	                        + zf2b(1,2)*fitmp2(i,j+1,k,l,nd) 
crs     &	                        + zf2b(1,3)*fitmp2(i,j+2,k,l,nd) 
crs     &	                        + zf2b(1,4)*fitmp2(i,j+3,k,l,nd) 
crs     &	                        + zf2b(1,5)*fitmp2(i,j+4,k,l,nd) 
crs     &	                        + zf2b(1,6)*fitmp2(i,j+5,k,l,nd) 
crs     &	                        + zf2b(1,7)*fitmp2(i,j+6,k,l,nd)
crs               end do
	     
	       j=2
	       do i=1,mrd	     
	       cu22(i,j,k,l,nd) = zf2b(2,1)*fitmp2(i,j-1,k,l,nd)
     &	                        + zf2b(2,2)*fitmp2(i,j,k,l,nd) 
     &	                        + zf2b(2,3)*fitmp2(i,j+1,k,l,nd) 
     &	                        + zf2b(2,4)*fitmp2(i,j+2,k,l,nd) 
     &	                        + zf2b(2,5)*fitmp2(i,j+3,k,l,nd) 
     &	                        + zf2b(2,6)*fitmp2(i,j+4,k,l,nd) 
     &	                        + zf2b(2,7)*fitmp2(i,j+5,k,l,nd)
               end do
	     
	       j=3
	       do i=1,mrd	     
	       cu22(i,j,k,l,nd) = zf2b(3,1)*fitmp2(i,j-2,k,l,nd)
     &	                        + zf2b(3,2)*fitmp2(i,j-1,k,l,nd) 
     &	                        + zf2b(3,3)*fitmp2(i,j,k,l,nd) 
     &	                        + zf2b(3,4)*fitmp2(i,j+1,k,l,nd) 
     &	                        + zf2b(3,5)*fitmp2(i,j+2,k,l,nd)
     &	                        + zf2b(3,6)*fitmp2(i,j+3,k,l,nd) 
     &	                        + zf2b(3,7)*fitmp2(i,j+4,k,l,nd)
               end do

	   
	     do j=4,nz2-3
	     jj=nz1tot+j	         
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
	     j=nz2-2
	     jj=nz1tot+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,1,k,l,nd+1)
      
	     end do
	     
	     j=nz2-1
	     jj=nz1tot+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,1,k,l,nd+1)
     &	                       + zf1(jj,3)*fitmp2(i,2,k,l,nd+1)
      
	     end do
	     
	     j=nz2
	     jj=nz1tot+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,1,k,l,nd+1) 
     &	                       + zf1(jj,2)*fitmp2(i,2,k,l,nd+1)
     &	                       + zf1(jj,3)*fitmp2(i,3,k,l,nd+1)
      
	     end do
	     

	   end do
	   end do
       
       elseif(nd.eq.ndom2) then !last domain of region 2
           
           do l=1,neqn
           do k=0,kh
	     
	     j=1
	     jj=nz1tot+(nd-1)*nz2+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,nz2-2,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp2(i,nz2-1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp2(i,nz2,k,l,nd-1) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     
	     j=2
	     jj=nz1tot+(nd-1)*nz2+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,nz2-1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp2(i,nz2,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     
	     j=3
	     jj=nz1tot+(nd-1)*nz2+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,nz2,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd)
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd)
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	   
	     do j=4,nz2-3
	     jj=nz1tot+(nd-1)*nz2+j	         
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
crs	     call writeit(3,1,mr2,nz2,1,cu22(1,1,0,1,nd)
crs     &	     ,fitmp2(1,1,0,1,nd),fitmp2(1,1,0,1,nd))
crs             stop
	     
	     j=nz2-2
	     jj=4
	     do i=1,mr2
	       cu22(i,j,k,l,nd) = zf2b(jj,1)*fitmp2(i,j-4,k,l,nd)
     &	                        + zf2b(jj,2)*fitmp2(i,j-3,k,l,nd) 
     &	                        + zf2b(jj,3)*fitmp2(i,j-2,k,l,nd) 
     &	                        + zf2b(jj,4)*fitmp2(i,j-1,k,l,nd) 
     &	                        + zf2b(jj,5)*fitmp2(i,j,k,l,nd) 
     &	                        + zf2b(jj,6)*fitmp2(i,j+1,k,l,nd) 
     &	                        + zf2b(jj,7)*fitmp2(i,j+2,k,l,nd)
             end do	     
             
	     j=nz2-1
	     jj=5
	     do i=1,mr2
	       cu22(i,j,k,l,nd) = zf2b(jj,1)*fitmp2(i,j-5,k,l,nd) 
     &	                        + zf2b(jj,2)*fitmp2(i,j-4,k,l,nd) 
     &	                        + zf2b(jj,3)*fitmp2(i,j-3,k,l,nd) 
     &	                        + zf2b(jj,4)*fitmp2(i,j-2,k,l,nd) 
     &	                        + zf2b(jj,5)*fitmp2(i,j-1,k,l,nd) 
     &	                        + zf2b(jj,6)*fitmp2(i,j,k,l,nd) 
     &	                        + zf2b(jj,7)*fitmp2(i,j+1,k,l,nd) 
             end do
	     
	     j=nz2
	     jj=6
	     do i=1,mr2
	       cu22(i,j,k,l,nd) = zf2b(jj,1)*fitmp2(i,j-6,k,l,nd) 
     &  		        + zf2b(jj,2)*fitmp2(i,j-5,k,l,nd) 
     &  		        + zf2b(jj,3)*fitmp2(i,j-4,k,l,nd) 
     &  		        + zf2b(jj,4)*fitmp2(i,j-3,k,l,nd) 
     &  		        + zf2b(jj,5)*fitmp2(i,j-2,k,l,nd) 
     &  		        + zf2b(jj,6)*fitmp2(i,j-1,k,l,nd) 
     &  		        + zf2b(jj,7)*fitmp2(i,j,k,l,nd) 
	     end do

	     
	   end do
	   end do
	   
       else ! any interior domain
       
           do l=1,neqn
           do k=0,kh
	     
	     j=1
	     jj=nz1tot+(nd-1)*nz2+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,nz2-2,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp2(i,nz2-1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp2(i,nz2,k,l,nd-1) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     
	     j=2
	     jj=nz1tot+(nd-1)*nz2+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,nz2-1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp2(i,nz2,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     
	     j=3
	     jj=nz1tot+(nd-1)*nz2+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,nz2,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	   
	     do j=4,nz2-3
	     jj=nz1tot+(nd-1)*nz2+j	         
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
	     j=nz2-2
	     jj=nz1tot+(nd-1)*nz2+j	
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp2(i,1,k,l,nd+1)
      
	     end do
	     
	     j=nz2-1
	     jj=nz1tot+(nd-1)*nz2+j	
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp2(i,1,k,l,nd+1)
     &	                       + zf1(jj,3)*fitmp2(i,2,k,l,nd+1)
      
	     end do
	     
	     j=nz2
	     jj=nz1tot+(nd-1)*nz2+j
	     do i=1,mr2
	      cu22(i,j,k,l,nd) = zf1(jj,-3)*fitmp2(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp2(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp2(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp2(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp2(i,1,k,l,nd+1) 
     &	                       + zf1(jj,2)*fitmp2(i,2,k,l,nd+1)
     &	                       + zf1(jj,3)*fitmp2(i,3,k,l,nd+1)
      
	     end do
	     

	   end do
	   end do
	   
       end if	   	   
       
       return
       
       end
       
       
c--------1---------2---------3---------4---------5---------6---------7--
c
      subroutine zfilter_exp_bl(nd)
c
c     Filter in z in place for one r/z plane of wave space data.
c
c--------1---------2---------3---------4---------5---------6---------7--


      implicit none
      include 'param.h'
#      include "common.h"
      include 'dgrp.h'
      include 'filter.h'


c --- Arguments ---
       
      integer i,j,k,l,nd,jj
      
      if(ndom1.eq.1) then
      
      do l=1,neqn
           do k=0,kh
	     
	       j=1
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(1,1)*fitmp1(i,j,k,l,nd)
     &	                        + zf1b(1,2)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(1,3)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(1,4)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(1,5)*fitmp1(i,j+4,k,l,nd) 
     &	                        + zf1b(1,6)*fitmp1(i,j+5,k,l,nd) 
     &	                        + zf1b(1,7)*fitmp1(i,j+6,k,l,nd)
               end do
	     
	       j=2
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(2,1)*fitmp1(i,j-1,k,l,nd)
     &	                        + zf1b(2,2)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(2,3)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(2,4)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(2,5)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(2,6)*fitmp1(i,j+4,k,l,nd) 
     &	                        + zf1b(2,7)*fitmp1(i,j+5,k,l,nd)
               end do
	     
	       j=3
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(3,1)*fitmp1(i,j-2,k,l,nd)
     &	                        + zf1b(3,2)*fitmp1(i,j-1,k,l,nd) 
     &	                        + zf1b(3,3)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(3,4)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(3,5)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(3,6)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(3,7)*fitmp1(i,j+4,k,l,nd)
               end do
	       
	     
	     do j=4,nz1-3	         
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(j,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
	     
	     j=nz1-2
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(j,3)*fitmp2(i+mrd,1,k,l,1)
      
	     end do
	     
	     j=nz1-1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp2(i+mrd,1,k,l,1)
     &	                       + zf1(j,3)*fitmp2(i+mrd,2,k,l,1)
      
	     end do
	     
	     j=nz1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp2(i+mrd,1,k,l,1) 
     &	                       + zf1(j,2)*fitmp2(i+mrd,2,k,l,1)
     &	                       + zf1(j,3)*fitmp2(i+mrd,3,k,l,1)
      
	     end do
	     
	     

	   end do
	   end do
      

      else !if ndom1.ne.1


             
       if(nd.eq.1) then !inflow
           do l=1,neqn
           do k=0,kh
	     
	       j=1
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(1,1)*fitmp1(i,j,k,l,nd)
     &	                        + zf1b(1,2)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(1,3)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(1,4)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(1,5)*fitmp1(i,j+4,k,l,nd) 
     &	                        + zf1b(1,6)*fitmp1(i,j+5,k,l,nd) 
     &	                        + zf1b(1,7)*fitmp1(i,j+6,k,l,nd)
               end do
	     
	       j=2
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(2,1)*fitmp1(i,j-1,k,l,nd)
     &	                        + zf1b(2,2)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(2,3)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(2,4)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(2,5)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(2,6)*fitmp1(i,j+4,k,l,nd) 
     &	                        + zf1b(2,7)*fitmp1(i,j+5,k,l,nd)
               end do
	     
	       j=3
	       do i=1,mr1	     
	       cu12(i,j,k,l,nd) = zf1b(3,1)*fitmp1(i,j-2,k,l,nd)
     &	                        + zf1b(3,2)*fitmp1(i,j-1,k,l,nd) 
     &	                        + zf1b(3,3)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(3,4)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(3,5)*fitmp1(i,j+2,k,l,nd) 
     &	                        + zf1b(3,6)*fitmp1(i,j+3,k,l,nd) 
     &	                        + zf1b(3,7)*fitmp1(i,j+4,k,l,nd)
               end do
	      	       
	     
	     do j=4,nz1-3	         
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(j,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
crs	     call writeit(3,1,mr1,nz1,1,cu12(1,1,0,1,nd)
crs     &	     ,fitmp1(1,1,0,1,nd),fitmp1(1,1,0,1,nd))
crs             stop 'z-filter'

	     
	     j=nz1-2
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(j,3)*fitmp1(i,1,k,l,nd+1)
      
	     end do
	     
	     j=nz1-1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(j,2)*fitmp1(i,1,k,l,nd+1)
     &	                       + zf1(j,3)*fitmp1(i,2,k,l,nd+1)
      
	     end do
	     
	     j=nz1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(j,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(j,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(j,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(j,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(j,1)*fitmp1(i,1,k,l,nd+1) 
     &	                       + zf1(j,2)*fitmp1(i,2,k,l,nd+1)
     &	                       + zf1(j,3)*fitmp1(i,3,k,l,nd+1)
      
	     end do
	     

	   end do
	   end do
       
       elseif(nd.eq.ndom1) then !last domain of region 1
           
           do l=1,neqn
           do k=0,kh
	     
	     j=1
	     jj=(nd-1)*nz1+1
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1-2,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,nz1-1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     
	     j=2
	     jj=(nd-1)*nz1+2
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1-1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     
	     j=3
	     jj=(nd-1)*nz1+3
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd)
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd)
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	   
	     do j=4,nz1-3
	     jj=(nd-1)*nz1+j	         
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
	     j=nz1-2
	     jj=4
	     do i=1,mr1
	       cu12(i,j,k,l,nd) = zf1b(jj,1)*fitmp1(i,j-4,k,l,nd) 
     &	                        + zf1b(jj,2)*fitmp1(i,j-3,k,l,nd) 
     &	                        + zf1b(jj,3)*fitmp1(i,j-2,k,l,nd) 
     &	                        + zf1b(jj,4)*fitmp1(i,j-1,k,l,nd) 
     &	                        + zf1b(jj,5)*fitmp1(i,j,k,l,nd)
     &	                        + zf1b(jj,6)*fitmp1(i,j+1,k,l,nd) 
     &	                        + zf1b(jj,7)*fitmp1(i,j+2,k,l,nd)
             end do	     
             
	     j=nz1-1
	     jj=5
	     do i=1,mr1
	       cu12(i,j,k,l,nd) = zf1b(jj,1)*fitmp1(i,j-5,k,l,nd) 
     &	                        + zf1b(jj,2)*fitmp1(i,j-4,k,l,nd) 
     &	                        + zf1b(jj,3)*fitmp1(i,j-3,k,l,nd) 
     &	                        + zf1b(jj,4)*fitmp1(i,j-2,k,l,nd) 
     &	                        + zf1b(jj,5)*fitmp1(i,j-1,k,l,nd)
     &	                        + zf1b(jj,6)*fitmp1(i,j,k,l,nd) 
     &	                        + zf1b(jj,7)*fitmp1(i,j+1,k,l,nd) 
             end do
	     
	     j=nz1
	     jj=6
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1b(jj,1)*fitmp1(i,j-6,k,l,nd) 
     &	                       + zf1b(jj,2)*fitmp1(i,j-5,k,l,nd) 
     &	                       + zf1b(jj,3)*fitmp1(i,j-4,k,l,nd) 
     &	                       + zf1b(jj,4)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1b(jj,5)*fitmp1(i,j-2  ,k,l,nd)
     &	                       + zf1b(jj,6)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1b(jj,7)*fitmp1(i,j  ,k,l,nd)
      
	     end do
	          	     
	   end do
	   end do
	   
       else ! any interior domain
       
           do l=1,neqn
           do k=0,kh
	     
	     j=1
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1-2,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,nz1-1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     
	     j=2
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1-1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     
	     j=3
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,nz1,k,l,nd-1) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	   
	     do j=4,nz1-3
	     jj=(nd-1)*nz1+j	         
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,j+3,k,l,nd)
      
	     end do
	     end do
	     
	     j=nz1-2
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,j+2,k,l,nd) 
     &	                       + zf1(jj,3)*fitmp1(i,1,k,l,nd+1)
      
	     end do
	     
	     j=nz1-1
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,j+1,k,l,nd) 
     &	                       + zf1(jj,2)*fitmp1(i,1,k,l,nd+1)
     &	                       + zf1(jj,3)*fitmp1(i,2,k,l,nd+1)
      
	     end do
	     
	     j=nz1
	     jj=(nd-1)*nz1+j
	     do i=1,mr1
	      cu12(i,j,k,l,nd) = zf1(jj,-3)*fitmp1(i,j-3,k,l,nd) 
     &	                       + zf1(jj,-2)*fitmp1(i,j-2,k,l,nd) 
     &	                       + zf1(jj,-1)*fitmp1(i,j-1,k,l,nd) 
     &	                       + zf1(jj,0)*fitmp1(i,j,k,l,nd) 
     &	                       + zf1(jj,1)*fitmp1(i,1,k,l,nd+1) 
     &	                       + zf1(jj,2)*fitmp1(i,2,k,l,nd+1)
     &	                       + zf1(jj,3)*fitmp1(i,3,k,l,nd+1)
      
	     end do
	     

	   end do
	   end do
	   
       end if
       
       end if !if ndom1.eq.1	   	   
       
       return
       
       end




c--------1---------2---------3---------4---------5---------6---------7--
c
      subroutine filter_temp1(nd)
c
c     Filter in z in place for one r/z plane of wave space data.
c
c--------1---------2---------3---------4---------5---------6---------7--


      implicit none
      include 'param.h'
#      include "common.h"
      include 'filter.h'


      integer i,j,k,l,nd
      
crs   copy array into temp first
      
      do l=1,neqn
         do k=0,kh
	    do j=1,nz1
	       do i=1,mr1
	          fitmp1(i,j,k,l,nd) = cu12(i,j,k,l,nd)
	       end do
	    end do
	 end do
       end do
       
       return
       end



c--------1---------2---------3---------4---------5---------6---------7--
c
      subroutine filter_temp2(nd)
c
c     Filter in z in place for one r/z plane of wave space data.
c
c--------1---------2---------3---------4---------5---------6---------7--


      implicit none
      include 'param.h'
#      include "common.h"
      include 'filter.h'

      integer i,j,k,l,nd
      
crs   copy array into temp first
      
      do l=1,neqn
         do k=0,kh
	    do j=1,nz2
	       do i=1,mr2
	          fitmp2(i,j,k,l,nd) = cu22(i,j,k,l,nd)
	       end do
	    end do
	 end do
       end do
       
       return
       end
       
       
c--------1---------2---------3---------4---------5---------6---------7--
c
      subroutine azfilter1(f1)
c
c     Filter in z in place for one r/z plane of wave space data.
c
c--------1---------2---------3---------4---------5---------6---------7--


      implicit none
      include 'param.h'
#      include "common.h"
      include 'filter.h'


      integer i,j,k
      real*8 f1(mr1,nz1,0:kh)
      
crs   scale higher modes
      

         do k=1,kh
	    do j=1,nz1
	       do i=1,mr1
	          f1(i,j,k) = kscale(k)*f1(i,j,k)
	       end do
	    end do
	 end do

       
       return
       end

c--------1---------2---------3---------4---------5---------6---------7--
c
      subroutine azfilter2(f2)
c
c     Filter in z in place for one r/z plane of wave space data.
c
c--------1---------2---------3---------4---------5---------6---------7--


      implicit none
      include 'param.h'
#      include "common.h"
      include 'filter.h'


      integer i,j,k
      real*8 f2(mr2,nz2,0:kh)
      
crs   scale higher modes
      

         do k=1,kh
	    do j=1,nz2
	       do i=1,mr2
	          f2(i,j,k) = kscale(k)*f2(i,j,k)
	       end do
	    end do
	 end do

       
       return
       end

       

