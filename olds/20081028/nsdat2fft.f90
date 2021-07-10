      MODULE sndat2fft

      contains
ccccc **********************************************************************
ccccc dat2fft: Main subroutine.
ccccc **********************************************************************
      SUBROUTINE dat2fft(t_ini,t_end,qtimes_local,M,N,fname)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*05 fname
      character*20 fname_aux
      character*17 fname_time
      integer i,k,tt,tt_local,qtimes_local,t_ini,t_end,M,N,error
      real*8 x(imax),y(jmax),z(kmax),
     &  rho(imax,jmax,kmax),u  (imax,jmax,kmax),
     &  v  (imax,jmax,kmax),w  (imax,jmax,kmax),
     &  Et (imax,jmax,kmax),p  (imax,jmax,kmax),
     &  tp (imax,jmax,kmax),wx (imax,jmax,kmax),
     &  wy (imax,jmax,kmax),wz (imax,jmax,kmax),
     &  U1 (imax,jmax,kmax),U2 (imax,jmax,kmax),
     &  U3 (imax,jmax,kmax),U4 (imax,jmax,kmax),
     &  U5 (imax,jmax,kmax)
      !MS$IF (tp_dim.EQ.2)
        real*8 pwr_sp((M-1)/2,qtimes_local)
      !MS$ELSEIF (tp_dim.EQ.3)
        real*8 pwr_sp((M-1)/2,(N-1)/2,qtimes_local)
      !MS$ENDIF

      call init_grid(dirgraf,x,y,z)
      if (y(jmax/2+1).GT.0.d0) then
        write(*,*) 'ERROR: IN Y DIRECTION IT NEEDS TO CHOOSE THE POINT'
        write(*,*) 'IN THE MIDLLE OF THE COMPUTATIONAL DOMAIN!        '
        write(*,*) 'THE PROGRAM WILL STOP.                            '
        stop
      end if

      write(*,*) ' _______________________________________________'
      write(*,*) 
      write(*,*) ' EXECUTANDO ANÁLISE DE FOURIER...               '   
      write(*,*)      
      	
      tt_local=0
      do tt=t_ini,t_end
        if (mod(tt,qtimes).EQ.0) then
          call get_filename(tt,fnamegraf,'.dat',fname_aux)           
  
          !MS$IF     (tp_dim.EQ.1)
            call read_BIN_1D(dirgraf,fname_aux,U1,U2,U5,error)
          !MS$ELSEIF (tp_dim.EQ.2)
            call read_BIN_2D(dirgraf,fname_aux,U1,U2,U3,U5,error)
          !MS$ELSEIF (tp_dim.EQ.3)
            call read_BIN_3D(dirgraf,fname_aux,U1,U2,U3,U4,U5,error)
          !MS$ENDIF          
  
          tt_local=tt_local+1
          if (error.EQ.0) then            
            call flux_to_primitive(U1,U2,U3,U4,U5,rho,u,v,w,Et,p,tp)
            call get_filename(tt,fname,'.plt',fname_aux)
            !MS$IF (tp_s.EQ.5).OR.(tp_s.EQ.6)
              call windowing(x,v)
            !MS$ENDIF

            !MS$IF     (tp_dim.EQ.1)                    
            !MS$ELSEIF (tp_dim.EQ.2)
              call DFT_1D(fname_aux,tt_local,qtimes_local,M,v,pwr_sp) 
            !MS$ELSEIF (tp_dim.EQ.3)
              call DFT_2D(fname_aux,tt_local,qtimes_local,M,N,v,pwr_sp) 
            !MS$ENDIF
          end if
        end if
      end do

      fname_time = fname//'_IN_TIME.plt'
      !MS$IF     (tp_dim.EQ.1)
      !MS$ELSEIF (tp_dim.EQ.2)
        call write_intime_2D(t_ini,t_end,qtimes_local,M,fname_time,pwr_sp)
      !MS$ELSEIF (tp_dim.EQ.3)
        call write_intime_3D(t_ini,t_end,qtimes_local,M,N,fname_time,pwr_sp)
      !MS$ENDIF

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
            fc         = (tanh(wp*(x(i)-x0))+tanh(wp*(x1-x(i))))/2.d0
            par(i,j,k) = par(i,j,k)*fc
          end do	
        end do	
      end do	

      RETURN
      END SUBROUTINE windowing

ccccc **********************************************************************
ccccc DFT_1D: This routine calculate the Fourier coefficients.
ccccc **********************************************************************
      SUBROUTINE DFT_1D(fname_aux,tt_local,qtimes_local,M,par,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'
	
      character*20 fname_aux
      integer i,tt_local,qtimes_local,M            
      real*8 par(imax,jmax,kmax),
     &  data(2*M),freq_x((M-1)/2),
     &  par_r ((M-1)/2,qtimes_local),
     &  par_i ((M-1)/2,qtimes_local),
     &  pwr_sp((M-1)/2,qtimes_local)
       
      do i=1,M-1
        data(i*2-1)=par(i,jmax/2+1,1)
        data(i*2  )=0.d0
      end do	

      call DFT(data,M-1,1)
      
      do i=1,(M-1)/2       
        freq_x(i         ) = dble(i-1  )
        par_r (i,tt_local) = data(i*2-1)
	par_i (i,tt_local) = data(i*2  )
	pwr_sp(i,tt_local) = 2.d0/dble(M-1)*
     &    dsqrt(par_r(i,tt_local)**2.d0+par_i(i,tt_local)**2.d0)        
      end do
      
      call write_2D(tt_local,fname_aux,M,freq_x,par_r,par_i,pwr_sp)

      RETURN
      END SUBROUTINE DFT_1D

ccccc **********************************************************************
ccccc FFT_2D: This routine calculate the Fourier coefficients.
ccccc **********************************************************************
      SUBROUTINE FFT_2D(fname_aux,tt_local,qtimes_local,M,N,par,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname_aux	
      integer i,k,tt_local,qtimes_local,M,N,dims(2)            
      real*8 par(imax,jmax,kmax),data(2*(M-1)*(N-1)),
     &  freq_x((M-1)/2),freq_z((N-1)/2),
     &  par_r ((M-1)/2,(N-1)/2,qtimes_local),
     &  par_i ((M-1)/2,(N-1)/2,qtimes_local),
     &  pwr_sp((M-1)/2,(N-1)/2,qtimes_local)

      do k=1,N-1
        do i=1,M-1
          data(i*2-1+(k-1)*2*(M-1))=par(i,jmax/2+1,k)
          data(i*2  +(k-1)*2*(M-1))=0.d0
        end do
      end do	

      dims(1)=M-1
      dims(2)=N-1
      !call FFT_ND(data,dims,2,1)
      
      do i=1,(M-1)/2
        do k=1,(N-1)/2 
          freq_x(i           ) = dble(i-1)
          freq_z(k           ) = dble(k-1)
	  par_r (i,k,tt_local) = data(i*2-1+(k-1)*2*(M-1))
	  par_i (i,k,tt_local) = data(i*2  +(k-1)*2*(M-1))
	  pwr_sp(i,k,tt_local) = 2.d0/dble((M-1)*(N-1))*
     &      dsqrt(par_r(i,k,tt_local)**2.d0+par_i(i,k,tt_local)**2.d0)
        end do
      end do
          	        
      call write_3D(tt_local,fname_aux,M,N,freq_x,freq_z,par_r,par_i,pwr_sp)              
  
      RETURN
      END SUBROUTINE FFT_2D

ccccc **********************************************************************
ccccc DFT_2D: This routine carried out the DFT for n dimension.
ccccc **********************************************************************
      SUBROUTINE DFT_2D(fname_aux,tt_local,qtimes_local,M,N,par,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname_aux      
      integer i,k,M,N,tt_local,qtimes_local
      real*8 par(imax,jmax,kmax),data_dft(2*(M-1),2*(N-1)),
     &  data_x(2*(M-1)),data_z(2*(N-1)),
     &  freq_x((M-1)/2),freq_z((N-1)/2),
     &  par_r ((M-1)/2,(N-1)/2,qtimes_local),
     &  par_i ((M-1)/2,(N-1)/2,qtimes_local),
     &  pwr_sp((M-1)/2,(N-1)/2,qtimes_local)

      !DFT in the streamwise direction
      do k=1,N-1
        do i=1,M-1
          data_x(i*2-1)=par(i,jmax/2+1,k)
          data_x(i*2  )=0.d0
        end do

        call DFT(data_x,M-1,1)
        do i=1,2*(M-1)
          data_dft(i,k)=data_x(i)
        end do
      end do

      !DFT in the spanwise direction
      do i=1,M-1
        do k=1,N-1
          data_z(k*2-1)=data_dft(i*2-1,k)
          data_z(k*2  )=data_dft(i*2  ,k)
        end do

        call DFT(data_z,N-1,1)
        do k=1,2*(N-1)
          data_dft(i,k)=data_z(k)
        end do
      end do

      do i=1,(M-1)/2
        do k=1,(N-1)/2       
          freq_x(i           ) = dble(i-1)
          freq_z(k           ) = dble(k-1)
	  par_r (i,k,tt_local) = data_dft(i,k*2-1)
	  par_i (i,k,tt_local) = data_dft(i,k*2  )
	  pwr_sp(i,k,tt_local) = 2.d0/dble((M-1)*(N-1))*
     &      dsqrt(par_r(i,k,tt_local)**2.d0+par_i(i,k,tt_local)**2.d0)
        end do
      end do

      call write_3D(tt_local,fname_aux,M,N,freq_x,freq_z,par_r,par_i,pwr_sp)              

      RETURN
      END SUBROUTINE DFT_2D

ccccc **********************************************************************
ccccc write_2D: Write data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_2D(tt,fname,M,freq_x,par_r,par_i,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,tt,M,size,Debug,IORE,TecIni,
     &  TecDat,TecZne,TecEnd,VIsPrecis
      real*8 freq_x(M/2),
     &  par_r   (M/2,1,qtimes),
     &  par_i   (M/2,1,qtimes),
     &  pwr_sp  (M/2,1,qtimes)
      real*4 freq_x_4(M/2),
     &  par_r_4 (M/2),
     &  par_i_4 (M/2),
     &  pwr_sp_4(M/2)

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsPrecis = 0
  
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'f v_r v_i pwr_sp'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  M/2,1,1,'BLOCK'//NULLCHR,NULLCHR)

      size = dble(M/2)

      do i=1,M/2
        freq_x_4(i) = freq_x(i   )
        par_r_4 (i) = par_r (i,1,tt)
        par_i_4 (i) = par_i (i,1,tt)
        pwr_sp_4(i) = pwr_sp(i,1,tt)
      end do

      IORE = TecDat(size,freq_x_4,VIsPrecis)
      IORE = TecDat(size,par_r_4 ,VIsPrecis)
      IORE = TecDat(size,par_i_4 ,VIsPrecis)
      IORE = TecDat(size,pwr_sp_4,VIsPrecis)
      IORE = TecEnd()

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
        
      RETURN
      END SUBROUTINE write_2D

ccccc **********************************************************************
ccccc write_FFT_TCP_3D: Write 3D data in file *.dat for TECPLOT.
ccccc **********************************************************************
      SUBROUTINE write_3D(tt,fname,M,N,freq_x,freq_z,par_r,par_i,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*20 fname
      character*1 NULLCHR
      integer i,k,tt,M,N,size,Debug,IORE,TecIni,
     &  TecDat,TecZne,TecEnd,VIsPrecis
      real*8 freq_x((M-1)/2),freq_z((N-1)/2),
     &  par_r   ((M-1)/2,(N-1)/2,qtimes),
     &  par_i   ((M-1)/2,(N-1)/2,qtimes),
     &  pwr_sp  ((M-1)/2,(N-1)/2,qtimes)
      real*4 freq_x_4((M-1)/2,(N-1)/2),freq_z_4((M-1)/2,(N-1)/2),
     &  par_r_4 ((M-1)/2,(N-1)/2),
     &  par_i_4 ((M-1)/2,(N-1)/2),
     &  pwr_sp_4((M-1)/2,(N-1)/2)

      NULLCHR   = CHAR(0)
      Debug     = 0
      VIsPrecis = 0
  
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'fx fz v_r v_i pwr_sp'//NULLCHR,
     &  dirgraf//'f_'//fname//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsPrecis)

      IORE = TecZne('Simple Zone'//NULLCHR,
     &  (M-1)/2,(N-1)/2,1,'BLOCK'//NULLCHR,NULLCHR)

      size = dble(((M-1)/2)*((N-1)/2))

      do i=1,(M-1)/2
        do k=1,(N-1)/2
          freq_x_4(i,k) = freq_x(i     )
          freq_z_4(i,k) = freq_z(k     )
          par_r_4 (i,k) = par_r (i,k,tt)
          par_i_4 (i,k) = par_i (i,k,tt)
          pwr_sp_4(i,k) = pwr_sp(i,k,tt)
        end do	
      end do

      IORE = TecDat(size,freq_x_4,VIsPrecis)
      IORE = TecDat(size,freq_z_4,VIsPrecis)
      IORE = TecDat(size,par_r_4 ,VIsPrecis)
      IORE = TecDat(size,par_i_4 ,VIsPrecis)
      IORE = TecDat(size,pwr_sp_4,VIsPrecis)
      IORE = TecEnd()

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
        
      RETURN
      END SUBROUTINE write_3D

ccccc **********************************************************************
ccccc write_intime_2D: Write data for FFT modes in time.
ccccc **********************************************************************
      SUBROUTINE write_intime_2D(t_ini,t_end,qtimes_local,M,fname,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*17 fname
      character*1 NULLCHR
      integer Ai,tt,tt_local,qtimes_local,M,t_ini,t_end,size,
     &  Debug,IORE,TecIni,TecZne,TecDat,TecEnd,VIsPrecis
      real*8 pwr_sp((M-1)/2,qtimes_local)
      real*4 time_4(qtimes_local),pwr_sp_4(qtimes_local)

      NULLCHR      = CHAR(0)
      Debug        = 0
      VIsPrecis    = 0  
      
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,'TIME M1 M2 M3 M4 M5 M6',
     &  dirgraf//'f_'//fname//NULLCHR,'.'//NULLCHR,Debug,VIsPrecis)
      IORE = TecZne('Simple Zone'//NULLCHR,qtimes_local,1,1,'BLOCK'//NULLCHR,NULLCHR)

      size = dble(qtimes_local)

      tt_local = 0
      do tt=t_ini,t_end           
        if (mod(tt,qtimes).EQ.0) then
          tt_local = tt_local+1
          time_4(tt_local) = dble(tt)*dt
        end if
      end do 
      IORE = TecDat(size,time_4,VIsPrecis)

      do Ai=2,7
        tt_local = 0
        do tt=t_ini,t_end           
          if (mod(tt,qtimes).EQ.0) then
            tt_local = tt_local+1
            pwr_sp_4(tt_local) = pwr_sp(Ai,tt_local)
          end if
        end do 
        IORE = TecDat(size,pwr_sp_4,VIsPrecis)
      end do

      IORE = TecEnd()

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
       
      RETURN
      END SUBROUTINE write_intime_2D

ccccc **********************************************************************
ccccc write_intime_3D: Write data for FFT modes in time.
ccccc **********************************************************************
      SUBROUTINE write_intime_3D(t_ini,t_end,qtimes_local,M,N,fname,pwr_sp)

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      character*17 fname
      character*1 NULLCHR
      integer Ai,Bi,tt,tt_local,qtimes_local,M,N,t_ini,t_end,size,
     &  Debug,IORE,TecIni,TecZne,TecDat,TecEnd,VIsPrecis
      real*8 pwr_sp((M-1)/2,(N-1)/2,qtimes_local)
      real*4 time_4 (qtimes_local),pwr_sp_4(qtimes_local)

      NULLCHR      = CHAR(0)
      Debug        = 0
      VIsPrecis    = 0  
      
      IORE = TecIni('SIMPLE DATASET'//NULLCHR,'TIME M1_2D M2_2D M3_2D M4_2D 
     &  M0_00 M0_01 M0_02 M0_03 M0_04 M0_05 M0_06 M0_07 M0_08 M0_09 M0_10 M0_11 M0_12
     &  M1_00 M1_01 M1_02 M1_03 M1_04 M1_05 M1_06 M1_07 M1_08 M1_09 M1_10 M1_11 M1_12
     &  M2_00 M2_01 M2_02 M2_03 M2_04 M2_05 M2_06 M2_07 M2_08 M2_09 M2_10 M2_11 M2_12',
     &  dirgraf//'f_'//fname//NULLCHR,'.'//NULLCHR,Debug,VIsPrecis)
      IORE = TecZne('Simple Zone'//NULLCHR,qtimes_local,1,1,'BLOCK'//NULLCHR,NULLCHR)

      size = dble(qtimes_local)

      tt_local = 0
      do tt=t_ini,t_end           
        if (mod(tt,qtimes).EQ.0) then
          tt_local = tt_local+1
          time_4(tt_local) = dble(tt)*dt
        end if
      end do 
      IORE = TecDat(size,time_4,VIsPrecis)

      do Ai=2,5
        tt_local = 0
        do tt=t_ini,t_end           
          if (mod(tt,qtimes).EQ.0) then
            tt_local = tt_local+1
            pwr_sp_4(tt_local) = pwr_sp(Ai,1,tt_local)
          end if
        end do 
        IORE = TecDat(size,pwr_sp_4,VIsPrecis)
      end do

      do Ai=1,3
        do Bi=1,13
          tt_local = 0
          do tt=t_ini,t_end           
            if (mod(tt,qtimes).EQ.0) then
              tt_local = tt_local+1
              pwr_sp_4(tt_local) = pwr_sp(Ai,Bi,tt_local)
            end if
          end do 
          IORE = TecDat(size,pwr_sp_4,VIsPrecis)
        end do
      end do


      IORE = TecEnd()

      write(*,*) ' WRITE DATA IN FILE: ',dirgraf//'f_'//fname
      write(*,*)
       
      RETURN
      END SUBROUTINE write_intime_3D

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
        wpr   = -2.d0*dsin(0.5d0*theta)**2.d0
        wpi   = dsin(theta)
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

ccccc **********************************************************************
ccccc FFT_ND: This routine carried out the FTT for n dimension.
ccccc **********************************************************************
      SUBROUTINE FFT_ND(data,nn,ndim,isign)

      integer isign,ndim,nn(ndim)
      integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
      real data(*)
        !Replaces data by its ndim-dimensional discrete Fourier transform, if isign is input as 1.
        !The variable nn(1:ndim) is an integer array containing the lengths of each dimension (number 
        !of complex values), which MUST all be powers of 2. data is a real array of length twice the
        !product of these lengths, in which the data are stored as in a multidimensional complex
        !FORTRAN array. If isign is input as −1, data is replaced by its inverse transform times
        !the product of the lengths of all dimensions.
      real tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp !Double precision for trigonometric recurrences.
      
      ntot=1
      do idim=1,ndim !Compute total number of complex values.
        ntot=ntot*nn(idim)
      end do
      nprev=1
      	
      do idim=1,ndim !Main loop over the dimensions.
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do i2=1,ip2,ip1 !This is the bit-reversal section of the routine.
          if (i2.lt.i2rev) then
            do i1=i2,i2+ip1-2,2
              do i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
              end do
            end do
          end if
          	
          ibit=ip2/2
   01     if ((ibit.GE.ip1).AND.(i2rev.GT.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
            goto 01
	  end if          
          i2rev=i2rev+ibit
        end do
                
        ifp1=ip1 !Here begins the Danielson-Lanczos section of the routine.
   02   if (ifp1.lt.ip2) then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1) !Initialize for the trig. recurrence.
          wpr=-2.d0*dsin(0.5d0*theta)**2.d0 
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do i3=1,ifp1,ip1
            do i1=i3,i3+ip1-2,2
              do i2=i1,ip3,ifp2
                k1=i2 !Danielson-Lanczos formula:
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi                
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
              end do 
            end do 
            wtemp=wr !Trigonometric recurrence.
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
          end do 
          ifp1=ifp2
          goto 02
        end if
        nprev=n*nprev
      enddo 

      RETURN
      END SUBROUTINE FFT_ND

      END MODULE sndat2fft
