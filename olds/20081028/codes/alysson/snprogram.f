      PROGRAM snttcl

      IMPLICIT NONE
      INCLUDE 'snparv36.f'

      integer tt,i,j,k,fieldnbr,fournbr,ttmax,nmore
      real*8 x(0:imax),y(0:jmax),z(0:kmax),
     &  u    (0:imax,0:jmax,0:kmax),v     (0:imax,0:jmax,0:kmax),
     &  w    (0:imax,0:jmax,0:kmax),rho   (0:imax,0:jmax,0:kmax),
     &  p    (0:imax,0:jmax,0:kmax),tp    (0:imax,0:jmax,0:kmax),
     &  s    (0:imax,0:jmax,0:kmax),Et    (0:imax,0:jmax,0:kmax),
     &  E1   (0:imax,0:jmax,0:kmax),E2    (0:imax,0:jmax,0:kmax),
     &  E3   (0:imax,0:jmax,0:kmax),E4    (0:imax,0:jmax,0:kmax),
     &  E5   (0:imax,0:jmax,0:kmax),   
     &  wx   (0:imax,0:jmax,0:kmax),wy    (0:imax,0:jmax,0:kmax),
     &  wz   (0:imax,0:jmax,0:kmax),
     &  u_0  (0:imax,0:jmax,0:kmax),tp_0  (0:imax,0:jmax,0:kmax),
     &  se_u (0:imax,0:jmax,0:kmax),se_v  (0:imax,0:jmax,0:kmax), 
     &  se_w (0:imax,0:jmax,0:kmax),se_rho(0:imax,0:jmax,0:kmax),
     &  se_p (0:imax,0:jmax,0:kmax),se_tp (0:imax,0:jmax,0:kmax), 
     &  se_s (0:imax,0:jmax,0:kmax),se_Et (0:imax,0:jmax,0:kmax),
     &  e_u  (0:imax,0:jmax,0:kmax),e_v   (0:imax,0:jmax,0:kmax),
     &  e_w  (0:imax,0:jmax,0:kmax),e_rho (0:imax,0:jmax,0:kmax),
     &  e_p  (0:imax,0:jmax,0:kmax),e_tp  (0:imax,0:jmax,0:kmax),
     &  e_s  (0:imax,0:jmax,0:kmax),e_Et  (0:imax,0:jmax,0:kmax),
     &  u_m,v_m,w_m,rho_m,p_m,tp_m,s_m,Et_m,
     &  se_u_m,se_v_m,se_w_m,se_rho_m,se_p_m,se_tp_m,se_s_m,se_Et_m,
     &  e_u_m,e_v_m,e_w_m,e_rho_m,e_p_m,e_tp_m,e_s_m,e_Et_m,
     &  stx_ex_6th_D(0:imax,0:6),stxx_ex_6th_D(0:imax,0:6),
     &  dpy_dcy(0:jmax),d2py_dcy2(0:jmax),freq_x  (0:imax-1,1:500),
     &  rea(0:imax-1,1:500), im(0:imax-1,1:500), time(0:imax-1,1:500),
     &  pwrsp(0:imax-1,1:500),reast(1:500), imst(1:500), timest(1:500),
     &  pwrspst(1:500),freq_xst(1:500),f
      character*1 NULLCHR
      integer Debug,IORE,xy,NPts,NElm,
     &  TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd,VIsDouble
      character*20 fname_aux
      common/par1/u,v,w
      common/par2/rho,p,tp
      common/par3/u_0,tp_0
      common/par4/s
      common/par5/Et
      common/par7/wx,wy,wz
      common/par8/se_u,se_v,se_w,se_rho,se_p,se_tp,se_s,se_Et
      common/par9/e_u,e_v,e_w,e_rho,e_p,e_tp,e_s,e_Et
      common/par50/stx_ex_6th_D
      common/par51/stxx_ex_6th_D
      common/parlay/fieldnbr,fournbr
      common/parsfft/rea,im,pwrsp,time,freq_x

      call ini_iosfile
      call ini_data(x,y,z)     
      call ini_cond(x,y,z)
      call verify_data(x,y,z)     
      call view_data

      call sv_all(x,y,z,0)
      call sv_Fourier(0,1)
      !MS$IF (tp_form.EQ.2)
        call primitive_to_flux(E1,E2,E3,E4,E5) 
      !MS$ENDIF
      fieldnbr=1
      tt=1
      ttmax=tmax
      do while (tt.le.ttmax)
        write(*,'(A29,A1,I8,A1,f17.10,A1,f15.10)')
     &    'NUMERICAL SOLUTION (tt,t) = ','(',tt,',',dble(tt)*dt,')'
        !MS$IF (tp_ci.EQ.1).OR.(tp_ci.EQ.2)
          call exact_acoustic(x,y,z,tt)
          call numerical_error
        !MS$ENDIF
        !MS$IF (tp_form.EQ.1)
          call RK4_FNCON(x,y,z,tt)
        !MS$ELSEIF (tp_form.EQ.2)
          call RK4(x,y,z,tt,E1,E2,E3,E4,E5)
        !MS$ENDIF
        
	if (((mod(tt,qtimes).EQ.0)).AND.(tt.GT.0)) then
          fieldnbr=fieldnbr+1
          !MS$IF (tp_form.EQ.2)
            call flux_to_primitive(E1,E2,E3,E4,E5) 
          !MS$ENDIF
          call sv_all(x,y,z,tt)
	  if (fieldnbr.le.500) then 
	    call sv_Fourier(tt,fieldnbr)
	  else
	    write(*,*)"******    Field Number Trepassed 500!!    ******"
	    write(*,*)"******        no more spatial fft         ******"
	    write(*,*)"******    computation will be stored      ******"
	  endif
        end if
	if (tt.eq.ttmax) then
	  write(*,*)"******Simulation time already reached",ttmax,
     &              "******"
          write(*,*)"**** Insert the number of time steps to run ****"
	  read(*,*) nmore
	  ttmax=ttmax+nmore
	  nmore=0
	endif
	tt=tt+1
      end do
      if (fieldnbr.le.500) then
        NULLCHR   = CHAR(0)
        Debug     = 0
        VIsDouble = 1
        fname_aux='SNFFT_0FULLTIME0.plt'
        IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'Freq_x t R_F_F_T I_F_F_T Pwr_F_F_T'//NULLCHR,
     &  dirgraf//'f_'//fname_aux//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsDouble)

        IORE = TecZne('Simple Zone'//NULLCHR,imax,fieldnbr,1,
     &  'BLOCK'//NULLCHR,NULLCHR)

        xy = dble(imax)*dble(fieldnbr)

        IORE   = TecDat(xy,freq_x, VIsDouble)
        IORE   = TecDat(xy,time, VIsDouble)
        IORE   = TecDat(xy,rea,  VIsDouble)
        IORE   = TecDat(xy,im, VIsDouble)
        IORE   = TecDat(xy,pwrsp, VIsDouble)

        IORE   = TecEnd()

        do i=1,fieldnbr
          f = dble(1)*2.d0*Pi/(imax*dx)
          freq_xst(i)=f
	  timest(i) = time (1,i)
	  reast(i)= rea(1,i)
	  imst(i) = im(1,i)
	  pwrspst(i)= (reast(i)**2.d0+
     &                  imst(i)**2.d0)**0.5d0
        end do
        NULLCHR   = CHAR(0)
        Debug     = 0
        VIsDouble = 1
        fname_aux='SNFFT_1STMODEFT0.plt'
        
	IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'Freq_x t R_F_F_T I_F_F_T Pwr_F_F_T'//NULLCHR,
     &  dirgraf//'f_'//fname_aux//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsDouble)
        IORE = TecZne('Simple Zone'//NULLCHR,1,fieldnbr,1,
     &  'BLOCK'//NULLCHR,NULLCHR)

        xy = dble(fieldnbr)

        IORE   = TecDat(xy,freq_xst, VIsDouble)
        IORE   = TecDat(xy,timest, VIsDouble)
        IORE   = TecDat(xy,reast,  VIsDouble)
        IORE   = TecDat(xy,imst, VIsDouble)
        IORE   = TecDat(xy,pwrspst, VIsDouble)

        IORE   = TecEnd()
      else
        NULLCHR   = CHAR(0)
        Debug     = 0
        VIsDouble = 1
        fname_aux='SNFFT_0FULLTIME0.plt'
        
	IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'Freq_x t R_F_F_T I_F_F_T Pwr_F_F_T'//NULLCHR,
     &  dirgraf//'f_'//fname_aux//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsDouble)

        IORE = TecZne('Simple Zone'//NULLCHR,imax,500,1,
     &  'BLOCK'//NULLCHR,NULLCHR)

        xy = dble(imax)*dble(500)

        IORE   = TecDat(xy,freq_x, VIsDouble)
        IORE   = TecDat(xy,time, VIsDouble)
        IORE   = TecDat(xy,rea,  VIsDouble)
        IORE   = TecDat(xy,im, VIsDouble)
        IORE   = TecDat(xy,pwrsp, VIsDouble)

        IORE   = TecEnd()

        do i=1,500
          f = dble(1)*2.d0*Pi/(imax*dx)
          freq_xst(i)=f
	  timest(i) = time (1,i)
	  reast(i)= rea(1,i)
	  imst(i) = im(1,i)
	  pwrspst(i)= (reast(i)**2.d0+
     &                  imst(i)**2.d0)**0.5d0
        end do
        NULLCHR   = CHAR(0)
        Debug     = 0
        VIsDouble = 1
        fname_aux='SNFFT_1STMODEFT0.plt'
        
	IORE = TecIni('SIMPLE DATASET'//NULLCHR,
     &  'Freq_x t R_F_F_T I_F_F_T Pwr_F_F_T'//NULLCHR,
     &  dirgraf//'f_'//fname_aux//NULLCHR,
     &  '.'//NULLCHR,Debug,VIsDouble)
        
	IORE = TecZne('Simple Zone'//NULLCHR,1,500,1,
     &  'BLOCK'//NULLCHR,NULLCHR)

        xy = dble(500)

        IORE   = TecDat(xy,freq_xst, VIsDouble)
        IORE   = TecDat(xy,timest, VIsDouble)
        IORE   = TecDat(xy,reast,  VIsDouble)
        IORE   = TecDat(xy,imst, VIsDouble)
        IORE   = TecDat(xy,pwrspst, VIsDouble)

        IORE   = TecEnd()
      endif
      
      call end_iosfile

!***************Tecplot Layout file printing (in "snlay.f")********************!

!      call snlayheader

      END

cccccc **********************************************************************
cccccc ini_data: Initialize data.
cccccc **********************************************************************
      SUBROUTINE ini_data(x,y,z)

      IMPLICIT NONE
      INCLUDE 'snparv36.f'

      integer i,j,k,id,ii
      real*8 x(0:imax),y(0:jmax),z(0:kmax),cxx,t,
     &  px_h,B,tau,sp,xs,xz,aa,bb,
     &  u      (0:imax,0:jmax,0:kmax),v     (0:imax,0:jmax,0:kmax),
     &  w      (0:imax,0:jmax,0:kmax),rho   (0:imax,0:jmax,0:kmax),
     &  p      (0:imax,0:jmax,0:kmax),tp    (0:imax,0:jmax,0:kmax),
     &  s      (0:imax,0:jmax,0:kmax),Et    (0:imax,0:jmax,0:kmax),
     &  wx     (0:imax,0:jmax,0:kmax),wy    (0:imax,0:jmax,0:kmax),
     &  wz     (0:imax,0:jmax,0:kmax),yc    (0:jmax),
     &  se_u   (0:imax,0:jmax,0:kmax),se_v  (0:imax,0:jmax,0:kmax), 
     &  se_w   (0:imax,0:jmax,0:kmax),se_rho(0:imax,0:jmax,0:kmax),
     &  se_p   (0:imax,0:jmax,0:kmax),se_tp (0:imax,0:jmax,0:kmax),  
     &  se_s   (0:imax,0:jmax,0:kmax),se_Et (0:imax,0:jmax,0:kmax),
     &  u_0    (0:imax,0:jmax,0:kmax),tp_0  (0:imax,0:jmax,0:kmax),
     &  yp     (0:imax,0:jmax,0:kmax),dypdyc(0:imax,0:jmax,0:kmax),
     &  d2ypdyc2(0:imax,0:jmax,0:kmax),xx,yy
      common/par1/u,v,w
      common/par2/rho,p,tp
      common/par3/u_0,tp_0
      common/par4/s
      common/par5/Et
      common/par7/wx,wy,wz
      common/par8/se_u,se_v,se_w,se_rho,se_p,se_tp,se_s,se_Et
      common/parmetric/dypdyc,d2ypdyc2
      
      do k=0,kmax
        do j=0,jmax
          do i=0,imax
            se_u  (i,j,k) = 0.d0
            se_v  (i,j,k) = 0.d0
            se_w  (i,j,k) = 0.d0
            se_rho(i,j,k) = 0.d0
            se_p  (i,j,k) = 0.d0
            se_tp (i,j,k) = 0.d0
            se_s  (i,j,k) = 0.d0
            se_Et (i,j,k) = 0.d0

            u     (i,j,k) = 0.d0
            v     (i,j,k) = 0.d0
            w     (i,j,k) = 0.d0
            rho   (i,j,k) = 0.d0
            p     (i,j,k) = 0.d0
            tp    (i,j,k) = 0.d0
            s     (i,j,k) = 0.d0
            Et    (i,j,k) = 0.d0

            wx    (i,j,k) = 0.d0
            wy    (i,j,k) = 0.d0                
            wz    (i,j,k) = 0.d0

            u_0   (i,j,k) = 0.d0
            tp_0  (i,j,k) = 0.d0
          end do   
    	end do  
      end do  

      !MS$IF (tp_ci.EQ.1).OR.(tp_ci.EQ.2).OR.(tp_ci.EQ.3)
        do i=0,imax
          x(i) = dble(i)*dx
        end do
      !MS$ELSEIF (tp_ci.EQ.4)
        !MS$IF (stretching_in_x.EQ.0)
          do i=0,imax
            x(i) = dble(i)*dx
          end do
        !MS$ELSEIF (stretching_in_x.EQ.1)  
          sp      = 1.5d0
          ii      = 100
          id      = imax-ii
          x(imax) = 100.d0
          do i=0,ii+1
            x(i) = dble(i)*dx
          end do 
          aa = ((x(imax)-x(ii))-dble(id)*(x(ii+1)-x(ii)))/
     &      (dble(id)**sp-dble(id))
          bb = (x(ii+1)-x(ii))-aa 
          do i=ii,imax
            x(i) = x(ii)+(aa*dble(i-ii)**sp+bb*dble(i-ii))
          end do
          
          sp = 0.91d0
          do i=0,imax
            cxx  = dble(i)/dble(imax) 
            x(i) = Lx*(0.5d0+dasin(-sp*cos(Pi*cxx))/(2.d0*asin(sp)))
          end do
        !MS$ENDIF  
      !MS$ENDIF 

      !MS$IF (tp_ci.EQ.1).OR.(tp_ci.EQ.2)
        do j=0,jmax 
          y(j)=dble(j)*dy
        end do
      !MS$ELSEIF (tp_ci.EQ.3).OR.(tp_ci.EQ.4)
        do j=0,jmax
          yc(j)=DBLE(j)/jmax
        end do
        !MS$IF (stretching_in_y.EQ.0)
          do j=0,jmax
            y(j)=DBLE(j)*dy
            y(j)=y(j)-Ly/2.d0
          end do
        !MS$ELSEIF (stretching_in_y.EQ.1)  
          tau=1.8d0
	  B=1.d0/(2.d0*tau)*dlog((1.d0+(dexp(tau)-1.d0)*1.d0/2.d0)/
     &                       (1.d0+(dexp(-tau)-1.d0)*1.d0/2.d0))
          do j=0,jmax
            y(j)=Ly/2.d0*(1.d0+dsinh(tau*(yc(j)-B))/dsinh(tau*B))
            y(j)=y(j)-Ly/2.d0
          end do
        !MS$ENDIF 
        do i=0,imax
          do j=0,jmax
            do k=0,kmax
              yp(i,j,k)=y(j)
            enddo
          enddo
        enddo
      !MS$ENDIF 
      call dnfdxin_comp(yp,dypdyc,'678symmd1                ',
     &         imax,jmax,kmax,'y','d','d',0,imax,0,jmax,0,kmax,dcy)
      call dnfdxin_comp(yp,d2ypdyc2,'678symmd2                ',
     &         imax,jmax,kmax,'y','d','d',0,imax,0,jmax,0,kmax,dcy)
     
      do k=0,kmax
        z(k)=DBLE(k)*dz
      end do

      if (sv_MAT.EQ.1) then
        call M_S(1,'x.m     ','x     ',x,imax,0,0,-1)
        call M_S(1,'y.m     ','y     ',y,jmax,0,0,-1)
        call M_S(1,'z.m     ','z     ',z,kmax,0,0,-1)

        call M_S(1,'dx.m    ','dx    ',dx,0,0,0,-1)
        call M_S(1,'dy.m    ','dy    ',dy,0,0,0,-1)
        call M_S(1,'dz.m    ','dz    ',dz,0,0,0,-1)
        call M_S(1,'dt.m    ','dt    ',dt,0,0,0,-1)

        call M_S(1,'imax.m  ','imax  ',DBLE(imax),0,0,0,-1)
        call M_S(1,'jmax.m  ','jmax  ',DBLE(jmax),0,0,0,-1)
        call M_S(1,'kmax.m  ','kmax  ',DBLE(kmax),0,0,0,-1)
        call M_S(1,'tmax.m  ','tmax  ',DBLE(tmax),0,0,0,-1)    
      end if

      RETURN 
      END
