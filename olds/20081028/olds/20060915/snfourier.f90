ccccc **********************************************************************
ccccc sv_Fourier: This routine save the Fourier coefficients.
ccccc **********************************************************************
      SUBROUTINE sv_Fourier

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer i,j,k
      real*8 data(imax*2),dx2,f,
     &  u    (0:imax,0:jmax,0:kmax),v     (0:imax,0:jmax,0:kmax),
     &  w    (0:imax,0:jmax,0:kmax),rho   (0:imax,0:jmax,0:kmax),
     &  tp   (0:imax,0:jmax,0:kmax),Et    (0:imax,0:jmax,0:kmax),
     &  p    (0:imax,0:jmax,0:kmax)
      common/par1/u,v,w
      common/par2/rho,p,tp

      do k=0,kmax
        do j=0,jmax

          do i=0,imax
cc            data((i+1)*2-1) = v(i,j,k)
cc            data((i+1)*2  ) = 0.d0
          end do	

          do i=1,imax
            data(i*2-1) = cos(2.d0*pi*1.d0*dble(i-1)*dx)*
     &        tanh(2.d0*dble(j)*dy) 
            data(i*2  ) = 0.d0
          end do	

          call FFT(data,imax,1)

          open (1,file='f_f2c.dat',access='sequential')      
          do i=1,imax/2
            f = dble(i)/(imax*dx)
            write(1,*) data((i)*2-1),f
          end do					   
          close (unit=1)
          write(*,*) 'primeira t'
          read(*,*)	
        end do	
      end do	

      RETURN
      END


ccccc **********************************************************************
ccccc FFT: This routine carried out the Fourier Transform.
ccccc **********************************************************************
      SUBROUTINE FFT(data,nn,isign)

      integer isign,nn
      real*8 data(2*nn)
        ! Replaces data(1:2*nn) by its discrete Fourier transform, if isign 
        ! is input as 1; or replaces data(1:2*nn) by nn times its inverse 
        ! discrete Fourier transform, if isign is input as -1. Data is a 
        ! complex array of length nn or, equivalently, a real array of 
        ! length 2*nn. nn MUST be an integer power of 2 (this is not checked for!).
      integer i,istep,j,m,mmax,n
      real*8 tempi,tempr      
      real*8 theta,wi,wpi,wpr,wr,wtemp 
        !Double precision for the trigonometric recurrences. 

      n = 2*nn
      j = 1
      do i=1,n,2 !This is the bit-reversal section of the routine.
        if(j.gt.i)then
          tempr     = data(j) !Exchange the two complex numbers.
          tempi     = data(j+1)
          data(j)   = data(i)
          data(j+1) = data(i+1)
          data(i)   = tempr
          data(i+1) = tempi
        end if
        m = nn
    1   if ((m.ge.2).and.(j.gt.m)) then
          j = j-m
          m = m/2
          goto 1
        endif
        j = j+m
      end do 
      mmax = 2 !Here begins the Danielson-Lanczos section of the routine.
    2 if (n.gt.mmax) then !then Outer loop executed log2 nn times.
        istep = 2*mmax
        theta = 6.28318530717959d0/(isign*mmax) !Initialize for the trigonometric recurrence.
        wpr   = -2.d0*sin(0.5d0*theta)**2.d0
        wpi   = sin(theta)
        wr    = 1.d0
        wi    = 0.d0
        do m=1,mmax,2 !Here are the two nested inner loops.
          do i=m,n,istep
            j = i+mmax !This is the Danielson-Lanczos formula:
            tempr     = sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi     = sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)   = data(i)-tempr
            data(j+1) = data(i+1)-tempi
            data(i)   = data(i)+tempr
            data(i+1) = data(i+1)+tempi
          end do 
          wtemp = wr 
          wr    = wr*wpr-wi*wpi+wr
          wi    = wi*wpr+wtemp*wpi+wi
        end do 
        mmax = istep
        goto 2 
      end if 
      
      RETURN
      END

