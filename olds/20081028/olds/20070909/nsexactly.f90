      MODULE snexactly

      IMPLICIT NONE       
      INCLUDE 'nspar.f90'

      save
        real*8, dimension(:,:,:), allocatable :: se_rho
        real*8, dimension(:,:,:), allocatable :: se_u
        real*8, dimension(:,:,:), allocatable :: se_v
        real*8, dimension(:,:,:), allocatable :: se_w
        real*8, dimension(:,:,:), allocatable :: se_p
        real*8, dimension(:,:,:), allocatable :: se_tp
        real*8, dimension(:,:,:), allocatable :: se_Et

      contains
cccccc **********************************************************************
cccccc exact_acoustic: Exact solution for acoutic wave problem.
cccccc **********************************************************************
      SUBROUTINE exact_acoustic(x,y,z,tt)

      IMPLICIT NONE       
      INCLUDE 'nspar.f90'

      integer i,j,k,tt
      complex*16 img,w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,w_c_1,w_c_2,
     &  rho_c_1,rho_c_2,tp_c_1,tp_c_2
      real*8 x(imax),y(jmax),z(kmax),
     &  xx,yy,zz,t,u_1,v_1,w_1,rho_1,p_1,tp_1

      img = cmplx(0.d0,1.d0)

      !MS$IF (tp_proc.EQ.1)
      !MS$ELSEIF (tp_proc.EQ.2)
        if ((jmax.EQ.0).AND.(kmax.EQ.0)) then               !One-dimensional      
          call coef_1D_adia(w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,
     &      w_c_1,w_c_2,rho_c_1,rho_c_2,tp_c_1,tp_c_2)
        else if ((jmax.NE.0).AND.(kmax.EQ.0)) then          !Two-dimensional      
          call coef_2D_adia(w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,
     &      w_c_1,w_c_2,rho_c_1,rho_c_2,tp_c_1,tp_c_2)
        else if ((jmax.NE.0).AND.(kmax.NE.0)) then          !Three-Dimensional      
          call coef_3D_adia(w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,
     &      w_c_1,w_c_2,rho_c_1,rho_c_2,tp_c_1,tp_c_2)
        end if 
      !MS$ENDIF

      t=DBLE(tt)*dt

      do k=0,kmax
        zz = z(k)
        do j=0,jmax
          yy = y(j)
          do i=0,imax
            xx = x(i)

            u_1   = dreal(
     &        u_c_1*exp(img*alpha*xx+img*kapa*yy+img*beta*zz-w1*t)+
     &        u_c_2*exp(img*alpha*xx+img*kapa*yy+img*beta*zz-w2*t))
            v_1   = dreal(
     &        v_c_1*exp(img*alpha*xx+img*kapa*yy+img*beta*zz-w1*t)+
     &        v_c_2*exp(img*alpha*xx+img*kapa*yy+img*beta*zz-w2*t))
            w_1   = dreal(
     &        w_c_1*exp(img*alpha*xx+img*kapa*yy+img*beta*zz-w1*t)+
     &        w_c_2*exp(img*alpha*xx+img*kapa*yy+img*beta*zz-w2*t))
            rho_1 = dreal(
     &        rho_c_1*exp(img*alpha*xx+img*kapa*yy+img*beta*zz-w1*t)+
     &        rho_c_2*exp(img*alpha*xx+img*kapa*yy+img*beta*zz-w2*t))

            !MS$IF (tp_proc.EQ.1)
              tp_1 = dreal(
     &          tp_c_1*dexp(img*alpha*xx+img*kapa*yy+img*beta*zz-w1*t)+
     &          tp_c_2*dexp(img*alpha*xx+img*kapa*yy+img*beta*zz-w2*t))
              p_1  = 1.d0/(gamma*Ma**2.d0)*(rho_1+tp_1+tp_1*rho_1)
            !MS$ELSEIF (tp_proc.EQ.2)
              p_1  = 1.d0/(gamma*Ma**2.d0)*((1.d0+rho_1)**gamma-1.d0)
              tp_1 = (1.d0+gamma*Ma**2.d0*p_1)/(1.d0+rho_1)-1.d0
            !MS$ENDIF

            se_u  (i,j,k) = (                       + u_1   ) 
            se_v  (i,j,k) = (                       + v_1   ) 
            se_w  (i,j,k) = (                       + w_1   ) 
            se_rho(i,j,k) = ( 1.d0                  + rho_1 ) 
            se_tp (i,j,k) = ( 1.d0                  + tp_1  )
            se_p  (i,j,k) = ( 1.d0/(gamma*Ma**2.d0) + p_1   )          
            se_Et (i,j,k) = se_rho(i,j,k)*(se_p(i,j,k)/(se_rho(i,j,k)*
     &        (gamma-1.d0))+(se_u(i,j,k)**2.d0+se_v(i,j,k)**2.d0
     &        +se_w(i,j,k)**2.d0)/2.d0)
          end do 
	end do  	
      end do	 

      RETURN 
      END SUBROUTINE exact_acoustic

cccccc **********************************************************************
cccccc coef_1D_adia: Fourier Coeficients (adiabatic).
cccccc **********************************************************************
      SUBROUTINE coef_1D_adia(w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,
     &  w_c_1,w_c_2,rho_c_1,rho_c_2,tp_c_1,tp_c_2)

      IMPLICIT NONE       
      INCLUDE 'nspar.f90'

      complex*16 img,w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,w_c_1,w_c_2,
     &  rho_c_1,rho_c_2,tp_c_1,tp_c_2,cte1,cte2,cte3

      img = cmplx(0.d0,1.d0) 

      cte1 = img*sqrt(9.d0*Re**2.d0-4.d0*(Ma*alpha)**2.d0)
      cte2 = 2.d0*alpha*Ma
      cte3 = 3.d0*img*Re*A_1*Ma

      w1      = alpha/(3.d0*Re*Ma)*(cte2+cte1)
      u_c_1   = A_1
      v_c_1   = 0.d0
      w_c_1   = 0.d0
      rho_c_1 = cte3/(cte2+cte1)
      tp_c_1  = 0.d0

      !MS$IF (tp_s.EQ.2) !Progressive Waves
        w2      = 0.d0
        u_c_2   = 0.d0      
        v_c_2   = 0.d0
        w_c_2   = 0.d0
        rho_c_2 = 0.d0
        tp_c_2  = 0.d0
      !MS$ELSEIF (tp_s.EQ.3) !Standing Waves
        w2      = alpha/(3.d0*Re*Ma)*(cte2-cte1)
        u_c_2   = A_1
        v_c_2   = 0.d0
        w_c_2   = 0.d0
        rho_c_2 = -cte3/(-cte2+cte1)
        tp_c_2  = 0.d0
      !MS$ENDIF

      RETURN 
      END SUBROUTINE coef_1D_adia

cccccc **********************************************************************
cccccc coef_2D_adia: Fourier Coeficients (adiabatic).
cccccc **********************************************************************
      SUBROUTINE coef_2D_adia(w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,
     &    w_c_1,w_c_2,rho_c_1,rho_c_2,tp_c_1,tp_c_2)

      IMPLICIT NONE       
      INCLUDE 'nspar.f90'

      complex*16 img,w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,w_c_1,w_c_2,
     &  rho_c_1,rho_c_2,tp_c_1,tp_c_2,cte1,cte2,cte3,cte4 

      img = cmplx(0.d0,1.d0) 

      cte1 = 2.d0*Ma*alpha**2.d0+2.d0*Ma*kapa**2.d0
      cte2 = img*sqrt((alpha**2.d0+kapa**2.d0)*
     &       (-4.d0*Ma**2.d0*alpha**2.d0
     &        -4.d0*Ma**2.d0*kapa**2.d0
     &        +9.d0*Re**2.d0))
      cte3 = (alpha**2.d0+kapa**2.d0)      
      cte4 = 3.d0*img*Ma*Re

      w1      = 1.d0/(3.d0*Re*Ma)*(cte1+cte2)
      u_c_1   = A_1
      v_c_1   = kapa/alpha*u_c_1
      w_c_1   = 0.d0
      rho_c_1 = cte4*u_c_1*cte3/(alpha*(cte1+cte2))
      tp_c_1  = 0.d0

      !MS$IF (tp_s.EQ.2) !Progressive Waves
        w2      = 0.d0
        u_c_2   = 0.d0
        v_c_2   = 0.d0
        w_c_2   = 0.d0
        rho_c_2 = 0.d0
        tp_c_2  = 0.d0
      !MS$ELSEIF (tp_s.EQ.3) !Standing Waves
        w2      = 1.d0/(3.d0*Re*Ma)*(cte1-cte2)
        u_c_2   = A_1
        v_c_2   = kapa/alpha*u_c_2
        rho_c_2 = cte4*u_c_2*cte3/(alpha*(cte1-cte2))
        w_c_2   = 0.d0
        tp_c_2  = 0.d0
      !MS$ENDIF

      RETURN 
      END SUBROUTINE coef_2D_adia

cccccc **********************************************************************
cccccc coef_3D_adia: Fourier Coeficients (adiabatic).
cccccc **********************************************************************
      SUBROUTINE coef_3D_adia(w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,
     &    w_c_1,w_c_2,rho_c_1,rho_c_2,tp_c_1,tp_c_2)

      IMPLICIT NONE       
      INCLUDE 'nspar.f90'

      complex*16 img,w1,w2,u_c_1,u_c_2,v_c_1,v_c_2,
     &  w_c_1,w_c_2,rho_c_1,rho_c_2,tp_c_1,tp_c_2,cte1,cte2

      img = cmplx(0.d0,1.d0) 

      cte1    = 4.d0*Ma**2.d0*alpha**2.d0-3.d0*Re**2.d0 
      cte2    = img*sqrt(-cte1)

      w1      = (2.d0*alpha**2.d0*Ma+alpha*cte2)/Re/Ma
      u_c_1   = A_1
      v_c_1   = -1.d0*(81.d0*Re**2*cte2**5-87.d0*Ma**3*alpha**3*
     &cte1**2+738.d0*Ma**5*alpha**5*cte1-240.d0*Ma**5*alpha**5*Re**2
     &-325.d0*Ma**4*alpha**4*cte2**3-576.d0*Ma**6*alpha**6*cte2-27.d0*
     &Ma*alpha*cte1**3+117.d0*Ma**2*alpha**2*cte2**5-513.d0*Re**2*Ma*
     &alpha*cte1**2-1599.d0*Ma**3*alpha**3*Re**2*cte1+1287.d0*Ma**2*
     &alpha**2*Re**2*cte2**3+160.d0*Ma**7*alpha**7+984.d0*Ma**4*
     &alpha**4*Re**2*cte2)*A_1/(40.d0*Ma**4*alpha**4-74.d0*Ma**3*
     &alpha**3*cte2+25.d0*Ma**2*alpha**2*cte1+18.d0*Ma*alpha*cte2**3
     &-9.d0*cte1**2-24.d0*Re**2*Ma**2*alpha**2+42.d0*Re**2*Ma*
     &alpha*cte2-18.d0*Re**2*cte1)/(5.d0*Ma**2*alpha**2-8.d0*Ma*
     &alpha*cte2+3.d0*cte1)/(4.d0*Ma*alpha-3.d0*cte2)
      w_c_1   = -1.d0*A_1*(201.d0*Re**2*cte2*Ma**3*alpha**3-249.d0*
     &Ma**2*alpha**2*Re**2*cte1+135.d0*Re**2*cte2**3*Ma*alpha-27.d0*
     &Re**2*cte1**2+40.d0*Ma**6*alpha**6-114.d0*Ma**5*alpha**5*cte2+99
     &.d0*Ma**4*alpha**4*cte1-7.d0*Ma**3*alpha**3*cte2**3-27.d0*Ma
     &**2*alpha**2*cte1**2-60.d0*Ma**4*alpha**4*Re**2+9.d0*Ma*alpha*
     &cte2**5)/(40.d0*Ma**4*alpha**4-74.d0*Ma**3*alpha**3*cte2+2
     &5.d0*Ma**2*alpha**2*cte1+18.d0*Ma*alpha*cte2**3-9.d0*cte1**2
     &-24.d0*Re**2*Ma**2*alpha**2+42.d0*Re**2*Ma*alpha*cte2-18.d0*
     &Re**2*cte1)/(5.d0*Ma**2*alpha**2-8.d0*Ma*alpha*cte2+3.d0*cte1)
      rho_c_1 = cmplx(0.d0,3.d0)*A_1*Re*Ma*(4.d0*Ma**3*alpha**3-11.d0*
     &Ma**2*alpha**2*cte2+10.d0*Ma*alpha*cte1-3.d0*cte2**3)/(40.d0*
     &Ma**4*alpha**4-74.d0*Ma**3*alpha**3*cte2+25.d0*Ma**2*alpha**2*
     &cte1+18.d0*Ma*alpha*cte2**3-9.d0*cte1**2-24.d0*Re**2*Ma**2*
     &alpha**2+42.d0*Re**2*Ma*alpha*cte2-18.d0*Re**2*cte1)
      tp_c_1  = 0.d0

      !MS$IF (tp_s.EQ.2) !Progressive Waves
        w2      = cmplx(0.d0,0.d0)
        u_c_2   = cmplx(0.d0,0.d0)
        v_c_2   = cmplx(0.d0,0.d0)
        w_c_2   = cmplx(0.d0,0.d0)
        rho_c_2 = cmplx(0.d0,0.d0)
        tp_c_2  = cmplx(0.d0,0.d0)
      !MS$ELSEIF (tp_s.EQ.3) !Standing Waves
        w2      = cmplx(dreal(w1),dimag(w1*(-1.d0)))
        u_c_2   = u_c_1
        v_c_2   = cmplx(dreal(v_c_1),dimag(v_c_1*(-1.d0)))
        w_c_2   = cmplx(dreal(w_c_1),dimag(w_c_1*(-1.d0)))
        rho_c_2 = cmplx(dreal(rho_c_1*(-1.d0)),dimag(rho_c_1))
        tp_c_2  = cmplx(0.d0,0.d0)
      !MS$ENDIF

      RETURN 
      END SUBROUTINE coef_3D_adia

      END MODULE snexactly
