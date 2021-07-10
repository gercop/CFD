c ifc -r8 -o grid.o -c grid.f
c ifc -r8 -o grid grid.o -lios

      integer imax,jmax,i,j
      double precision Lx,Ly,dcx,dcy,x_0,B,B1,B2
      double precision dx,dy,cyy

      parameter (x_0 = 0.2d0)
      parameter (Lx=5.d0)
      parameter (Ly=0.12d0)

      parameter (imax=200)
      parameter (jmax=200)

      double precision x(0:imax),y(0:imax),py_aux(0:jmax)

      double precision geo(0:imax,0:jmax,2)

      dx  = Lx/dble(imax)       !Step size in the x-direction
      dy  = Ly/dble(jmax)       !Step size in the y-direction
      dcx = 1.d0/dble(imax)     !Step size in the x-direction
      dcy = 1.d0/dble(jmax)     !Step size in the y-direction

!In x-direction
      do i=0,imax
         x(i) = x_0 + dble(i)*dx
c         print '(f18.14)',x(i)
      end do

!In y-direction
      B  = 1.01d0
      B1 = ((B+1.d0)/(B-1.d0))
      B2 = dlog(B1)
      do j=0,jmax
         cyy       = dble(j)*dcy
         py_aux(j) = Ly*((B+1.d0)-(B-1.d0)*
     &        (B1**(1.d0-cyy)))/(B1**(1.d0-cyy)+1.d0)
      end do
      do j=0,jmax
         y(j) = py_aux(j)
         print '(e20.14)',max(y(j),0.)
      end do

      do i=0,imax
         do j=0,jmax
            geo(i,j,1)=x(i)
            geo(i,j,2)=y(j)
         enddo
      enddo

      call gridtoios(geo,2,imax,jmax,0)

      end






      subroutine gridtoios(geo,idim,ix,jx,kx)
* write grid to ios file 'grid'
      dimension geo(0:ix,0:jx,0:kx,idim)
      character filebase*128,info(128)*72
      integer times(512)

      itape=10
      iunit=11
      filebase='grid'
c      imach=10 ! ascii real*8
      imach=2 ! unformatted real*8
      ntime=1
      m3=kx+1
      m2=jx+1
      m1=ix+1
      nparam=idim
      do n=1,ntime
         times(n)=1
      enddo
      info(1)='x-grid'
      info(2)='y-grid'
      if (idim.eq.3) info(3)='z-grid'
      info(idim+1)='grid'
      ncomment=1
      call writecd(itape,iunit,filebase,imach,ntime,m3,m2,m1,
     &     nparam,times,info,ncomment)
      do id=1,idim
         call writed(itape,geo(0,0,0,id))
      enddo
      close(itape)

      end
