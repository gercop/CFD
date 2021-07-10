cccccc *******************************************************************
cccccc Doctorate Project USP/EESC: Numerical simulation of the transition
cccccc                             for turbulence in a compressible 
cccccc                             boundary layer on a plane plate.
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos
cccccc Programmer\Guiding: Marcello A. Faraco de Medeiros
cccccc Version : 1.1
cccccc *******************************************************************
cccccc This program carried out Fourier Spectral Analysis.
cccccc *******************************************************************
      MODULE sndat2fft

      contains
ccccc **********************************************************************
ccccc dat2fft: Main subroutine.
ccccc **********************************************************************
      SUBROUTINE dat2fft

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname_aux
      integer i,tt,field,t_ini,t_end,error
      real*8 x(imax),y(jmax),z(kmax),
     &  par(imax,jmax,kmax),freq_x(imax/2),
     &  par_r(imax/2),par_i(imax/2),pwr_sp(imax/2),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w  (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p  (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  U1 (imax,jmax,kmax),U2 (imax,jmax,kmax),
     &  U3 (imax,jmax,kmax),U4 (imax,jmax,kmax),
     &  U5 (imax,jmax,kmax)

      t_ini = 0
      t_end = 0

      write(*,*) 
      write(*,*) 'PROGRAMA PARA REALIZAR A ANÁLISE ESPECTRAL:'  
      write(*,'(a23,\)') 'TEMPO INICIAL........:'
      read (*,*) t_ini
      write(*,'(a23,\)') 'TEMPO FINAL..........:'
      read (*,*) t_end

      field=1     

      do tt=t_ini,t_end
        if (mod(tt,qtimes).EQ.0) then
          call get_filename(tt,'SN__S','.dat',fname_aux)           

          if ((jmax.NE.1).AND.(kmax.EQ.1)) then
            call read_BIN_2D(fname_aux,U1,U2,U3,U4,U5,error)
          end if

          call flux_to_primitive(imax,U1,U2,U3,U4,U5,rho,u,v,w,Et,p,tp)

          !MS$IF (tp_s.EQ.5).OR.(tp_s.EQ.6)
            call windowing(x,v)
          !MS$ENDIF 
          call fourier(v,freq_x,par_r,par_i,pwr_sp) 

          if ((jmax.NE.1).AND.(kmax.EQ.1)) then
            if (error.EQ.0) then
              call get_filename(tt,'SN_FF','.plt',fname_aux)
              call write_FFT_TCP_2D(fname_aux,freq_x,par_r,par_i,pwr_sp)
            end if
          end if    
        end if
      end do

      END SUBROUTINE dat2fft

ccccc **********************************************************************
ccccc windowing: 
ccccc **********************************************************************
      SUBROUTINE windowing(x,par)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer i,j,k
      real*8 x(imax),par(imax,jmax,kmax),fc,wp,x0,x1
      parameter (wp=1.d0)
 
      x0 = x(1+5)
      x1 = x(imax-5)

      do i=1,imax
        do j=1,jmax
          do k=1,kmax
            fc = ( tanh(wp*(x(i)-x0))+tanh(wp*(x1-x(i))) )/2.d0
            par(i,j,k) = par(i,j,k)*fc
          end do	
        end do	
      end do	


      RETURN
      END SUBROUTINE windowing

ccccc **********************************************************************
ccccc fourier: This routine calculate the Fourier coefficients.
ccccc **********************************************************************
      SUBROUTINE fourier(par,freq_x,par_r,par_i,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer tt,i
      real*8 data(imax*2),par(imax,jmax,kmax),
     &  freq_x(imax/2),par_r(imax/2),par_i(imax/2),pwr_sp(imax/2)

      do i=1,imax
        !Test case: data(i*2-1) = dcos(2.d0*pi*1.d0*dble(i-1)*dx)
        data(i*2-1) = par(i,jmax/2+1,1)        
        data(i*2  ) = 0.d0
      end do	

      call DFT(data,imax,1)

      do i=1,imax/2
        freq_x(i) = dble(i-1)/(imax*dx)
	par_r (i) = data(i*2-1)
	par_i (i) = data(i*2  )
	pwr_sp(i) = (par_r(i)**2.d0+par_i(i)**2.d0)**0.5d0
      end do

      RETURN
      END SUBROUTINE fourier

ccccc **********************************************************************
ccccc write_FFT_TCP_2D: Save data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_FFT_TCP_2D(fname,freq_x,par_r,par_i,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,size,Debug,IORE,NPts,NElm,TecIni,
     &  TecDat,TecZne,TecNod,TecFil,TecEnd,VIsPrecis
      real*8 
     &  freq_x(imax/2),par_r (imax/2),
     &  par_i (imax/2),pwr_sp(imax/2)
      real*4 
     &  freq_x_4(imax/2),par_r_4 (imax/2),
     &  par_i_4 (imax/2),pwr_sp_4(imax/2)

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsPrecis = 0
  
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'f v_r v_i pwr_sp'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  imax/2,1,1,'BLOCK'//NULLCHR,NULLCHR)

      size = dble(imax/2)

      do i=1,imax/2
        freq_x_4(i) = freq_x(i)
        par_r_4 (i) = par_r (i)
        par_i_4 (i) = par_i (i)
        pwr_sp_4(i) = pwr_sp(i)
      end do

      IORE = TecDat(size,freq_x_4,VIsPrecis)
      IORE = TecDat(size,par_r_4, VIsPrecis)
      IORE = TecDat(size,par_i_4, VIsPrecis)
      IORE = TecDat(size,pwr_sp_4,VIsPrecis)
      IORE = TecEnd()

      if (QTimes.GT.1) then
        write(*,*) 'WRITE FILE: ',dirgraf//'f_'//fname
      end if
      
      RETURN
      END SUBROUTINE write_FFT_TCP_2D

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
      END SUBROUTINE FFT

ccccc **********************************************************************
ccccc DFT: This routine carried out the Fourier Transform.
ccccc **********************************************************************
      SUBROUTINE DFT(data,nn,isign)

      integer isign,nn
      real*8 data(2*nn),dataaux(2*nn)
        ! Replaces data(1:2*nn) by its discrete Fourier transform, if isign
        ! is input as 1; or replaces data(1:2*nn) by nn times its inverse
        ! discrete Fourier transform, if isign is input as -1. Data is a
        ! complex array of length nn or, equivalently, a real array of
        ! length 2*nn.
      integer i,n
      real*8 tempi,tempr     
      real*8 wi,wr
      complex*16 w
      real*8 Pi

      Pi= 2.d0*dasin(1.d0)
      if (isign.eq.1) then ! Discrete fourier transform
      w = cmplx(dcos(2.d0*Pi/DBLE(nn)),dsin(2.d0*Pi/DBLE(nn)))
      do n=-nn/2,nn/2-1,1
        tempr = 0.d0
        tempi = 0.d0
        do k=0,nn-1
          wr = dreal(w**(DBLE(n)*DBLE(k)))
          wi = dimag(w**(DBLE(n)*DBLE(k)))
          tempr=tempr+data((k+1)*2-1)*wr-data((k+1)*2)*wi
          tempi=tempi+data((k+1)*2-1)*wi+data((k+1)*2)*wr
        enddo
        dataaux(2*(n+nn/2)+1)=tempr
        dataaux(2*(n+nn/2)+2)=tempi
      enddo
      n=0
      data(1)=dataaux(2*(n+nn/2)+1)
      data(2)=dataaux(2*(n+nn/2)+2)
      do n=1,nn/2-1,1
        data(2*n+1)=dataaux(2*(n+nn/2)+1)
        data(2*n+2)=dataaux(2*(n+nn/2)+2)
      enddo
      do n=-nn/2,-1,1
        data(2*(n+nn/2)+nn+1)=dataaux(2*(n+nn/2)+1)
        data(2*(n+nn/2)+nn+2)=dataaux(2*(n+nn/2)+2)
      enddo
      elseif (isign.eq.-1) then ! Inverse discrete fourier transform
        write(*,*)"The procedure for inverse Discrete Fourier",
     &            "Transform is not yet implemented!"
        write(*,*)"Press <ENTER>"
        read(*,*)
        STOP
      end if
           
      RETURN
      END SUBROUTINE DFT

      END MODULE sndat2fft
