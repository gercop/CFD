ccccc **********************************************************************
ccccc sv_Fourier: This routine save the Fourier coefficients.
ccccc **********************************************************************
      SUBROUTINE sv_Fourier(tt,field)

      IMPLICIT NONE
      INCLUDE 'snparv36.f'

      integer i,j,k
      integer tt,field
      real*8 data(imax*2),dx2,f,
     &  u    (0:imax,0:jmax,0:kmax),v     (0:imax,0:jmax,0:kmax),
     &  w    (0:imax,0:jmax,0:kmax),rho   (0:imax,0:jmax,0:kmax),
     &  tp   (0:imax,0:jmax,0:kmax),Et    (0:imax,0:jmax,0:kmax),
     &  p    (0:imax,0:jmax,0:kmax),freq_x  (0:imax-1,1:500),
     &  rea(0:imax-1,1:500), im(0:imax-1,1:500), time(0:imax-1,1:500),
     &  pwrsp(0:imax-1,1:500)
      common/par1/u,v,w
      common/par2/rho,p,tp
      common/parsfft/rea,im,pwrsp,time,freq_x

      do i=0,imax-1
        data((i+1)*2-1) = v(i,jmax/2,0)
        data((i+1)*2  ) = 0.d0
	time(i,field)=dt*dble(tt)
      end do	

      call FFT(data,imax,1)

      do i=0,imax/2
        f = dble(i)/(imax*dx)
        freq_x(i,field)=f
	rea(i,field)= data((i+1)*2-1)
	im(i,field)= data((i+1)*2  )
	pwrsp(i,field)= (rea(i,field)**2.d0+im(i,field)**2.d0)**0.5d0
      end do
      do i=imax/2+1,imax-1
        f = dble(i)/(imax*dx)
        freq_x(i,field)=f
	rea(i,field)= data((i+1)*2-1)
	im(i,field)= data((i+1)*2  )
	pwrsp(i,field)= (rea(i,field)**2.d0+im(i,field)**2.d0)**0.5d0 !power spectrum
      end do

      RETURN
      END


ccccc **********************************************************************
ccccc sv_Fourier: This routine save the Fourier coefficients.
ccccc **********************************************************************
      SUBROUTINE sv_Fourier_temp(x)

      IMPLICIT NONE
      INCLUDE 'snparv36.f'

      integer i,j,k,fieldnbr,fournbr
      character*20 fname_aux
      character*01 c1
      character*02 c2
      character*03 c3
      character*04 c4
      character*05 c5
      character*06 c6
      character*07 c7
      character*08 c8
      integer tt
      character*1 NULLCHR
      integer Debug,IORE,xy,NPts,NElm,
     &  TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd,VIsDouble
      real*8 x(0:imax),y(0:jmax),z(0:kmax),data(fournbr*2-2),dx2,f,
     &  xx  (0:imax,0:fournbr-1,0:kmax),v (0:imax,0:fournbr-1),
     &  freq_t (0:imax,0:fournbr-1,0:kmax),Rea(0:imax,0:fournbr-1),
     &  Im(0:imax,0:fournbr-1), time(0:fournbr-1)
      character*23 filename(0:fournbr-1)
      common/parlay/fieldnbr,fournbr

      open (1,file='SN_V0.log',access='sequential')
      rewind(unit=1)
      read(1,*)(filename(i-1),i=1,fournbr)
      close(unit=1)
      do j=0,fournbr-1
        open (1,file=dirgraf//filename(j),access='sequential')
        read(1,*)time(j)
        do i=0,imax
          read(1,*)v(i,j)
        enddo
        close(unit=1)
      enddo
      do i=0,imax
        do j=0,fournbr-1
          data((j+1)*2-1) = v(i,j)*((dtanh((6.d0*time(j)/
     &  (time(fournbr-1)-time(0))-1.d0)*5.d0)+dtanh((-6.d0*time(j)/
     &  (time(fournbr-1)-time(0))+5.d0)*5.d0))/2.d0)
          data((j+1)*2  ) = 0.d0
        end do
        call FFT(data,fournbr,1)
        do j=0,fournbr/2
          Rea(i,j)= data((j+1)*2-1)
          Im(i,j)= data((j+1)*2  )
        end do
        do j=fournbr/2+1,fournbr-1
          Rea(i,j)= data((j+1)*2-1)
          Im(i,j)= data((j+1)*2  )
        end do
      end do
      do j=0,fournbr/2
        f = dble(j)/(time(fournbr-1))
        freq_t(0,j,0)=f
      end do
      do j=fournbr/2+1,fournbr-1
        f = dble(j)/time(fournbr-1)
        freq_t(0,j,0)=f
      end do

      fname_aux='SN_V0_0000000000.plt'

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsDouble = 1
    
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'X Freq_t Re_F_F_T Im_F_F_T'//NULLCHR,
     &  dirgraf//'f_'//fname_aux//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsDouble)

      IORE = TecZne('Simple Zone'//NULLCHR,imax+1,fournbr,1,
     &  'BLOCK'//NULLCHR,NULLCHR)

      xy = dble((imax+1)*(fournbr))
      do i=0,imax
        do j=0,fournbr-1
          xx(i,j,0) = x(i)
          freq_t(i,j,0)=freq_t(0,j,0)
        end do
      end do

      IORE   = TecDat(xy,xx, VIsDouble)
      IORE   = TecDat(xy,freq_t, VIsDouble)
      IORE   = TecDat(xy,Rea,  VIsDouble)
      IORE   = TecDat(xy,Im, VIsDouble)

      IORE   = TecEnd()

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

