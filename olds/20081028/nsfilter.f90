      MODULE snfilter

      contains
ccccc **********************************************************************
ccccc filter_x: Numerical filter in x-direction.
ccccc **********************************************************************
      SUBROUTINE filter_x(local_imax,U1,U2,U3,U4,U5)

      USE snmpi, only: my_rank,pro,exchange_ghostpoints
           
      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer local_imax
      real*8 
     &  U1(local_imax,jmax,kmax),U2(local_imax,jmax,kmax),
     &  U3(local_imax,jmax,kmax),U4(local_imax,jmax,kmax),
     &  U5(local_imax,jmax,kmax)

      !MS$IF (filter_in_x.EQ.1)
        !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3).OR.(tp_s.EQ.4)
          call filter_x_cp_4th_P(imax,jmax,kmax,U1,U2,U3,U4,U5)        
        !MS$ENDIF
        !MS$IF (tp_s.EQ.5).OR.(tp_s.EQ.6)
          if (pro.NE.1) then
            call exchange_ghostpoints(jmax,kmax,1,U1)
            call exchange_ghostpoints(jmax,kmax,3,U2)
            call exchange_ghostpoints(jmax,kmax,5,U3)
            call exchange_ghostpoints(jmax,kmax,7,U4)
            call exchange_ghostpoints(jmax,kmax,9,U5)
          end if
          call filter_x_ex_6th(local_imax,jmax,kmax,U1,U2,U3,U4,U5)                             
        !MS$ENDIF
      !MS$ENDIF

      RETURN
      END SUBROUTINE filter_x

ccccc **********************************************************************
ccccc filter_y: Numerical filter in y-direction.
ccccc **********************************************************************
      SUBROUTINE filter_y(local_imax,U1,U2,U3,U4,U5)
           
      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer local_imax
      real*8 
     &  U1(local_imax,jmax,kmax),U2(local_imax,jmax,kmax),
     &  U3(local_imax,jmax,kmax),U4(local_imax,jmax,kmax),
     &  U5(local_imax,jmax,kmax)

      !MS$IF (filter_in_y.EQ.1)
        !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5).OR.(tp_s.EQ.6)
          call filter_y_cp_4th(local_imax,jmax,kmax,U1,U2,U3,U4,U5)
        !MS$ENDIF
      !MS$ENDIF

      RETURN
      END SUBROUTINE filter_y

ccccc **********************************************************************
ccccc filter_z: Numerical filter in z-direction.
ccccc **********************************************************************
      SUBROUTINE filter_z(local_imax,U1,U2,U3,U4,U5)
           
      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer local_imax
      real*8 
     &  U1(local_imax,jmax,kmax),U2(local_imax,jmax,kmax),
     &  U3(local_imax,jmax,kmax),U4(local_imax,jmax,kmax),
     &  U5(local_imax,jmax,kmax)

      !MS$IF (filter_in_z.EQ.1)
        call filter_z_cp_4th_P(local_imax,jmax,kmax,U1,U2,U3,U4,U5)
      !MS$ENDIF

      RETURN
      END SUBROUTINE filter_z

cccccc**********************************************************************
cccccc filter_x_cp_4th_P: Compact filter of fourth order in x-direction for 
cccccc                    periodic boundary conditions.
cccccc **********************************************************************
      SUBROUTINE filter_x_cp_4th_P(imax,jmax,kmax,U1,U2,U3,U4,U5)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U4(imax,jmax,kmax),U5(imax,jmax,kmax),
     &  AA(imax,imax),LL(imax,imax),UU(imax,imax),  
     &  U1_f(imax),U2_f(imax),U3_f(imax),U4_f(imax),U5_f(imax),
     &  rhs1(imax),rhs2(imax),rhs3(imax),rhs4(imax),rhs5(imax)

      call lhsf_P(imax,AA)
      call ludecomp(1,imax,AA,LL,UU)

      do k=1,kmax
        do j=1,jmax
          call rhsf_x_P(imax,jmax,kmax,j,k,U1,rhs1)
          call solve(1,imax,LL,UU,rhs1,U1_f)
          call rhsf_x_P(imax,jmax,kmax,j,k,U2,rhs2)
          call solve(1,imax,LL,UU,rhs2,U2_f)
          call rhsf_x_P(imax,jmax,kmax,j,k,U3,rhs3)
          call solve(1,imax,LL,UU,rhs3,U3_f)
          call rhsf_x_P(imax,jmax,kmax,j,k,U4,rhs4)
          call solve(1,imax,LL,UU,rhs4,U4_f)
          call rhsf_x_P(imax,jmax,kmax,j,k,U5,rhs5)
          call solve(1,imax,LL,UU,rhs5,U5_f)

          U1_f(imax)=U1_f(1)
          U2_f(imax)=U2_f(1)
          U3_f(imax)=U3_f(1)
          U4_f(imax)=U4_f(1)
          U5_f(imax)=U5_f(1)

          do i=1,imax
            U1(i,j,k)=U1_f(i)
            U2(i,j,k)=U2_f(i)
            U3(i,j,k)=U3_f(i)
            U4(i,j,k)=U4_f(i)
            U5(i,j,k)=U5_f(i)
          end do
        end do
      end do   

      RETURN
      END SUBROUTINE filter_x_cp_4th_P

cccccc **********************************************************************
cccccc filter_x_cp_4th: Compact filter of fourth order in x-direction for 
cccccc                  non-periodic boundary conditions.
cccccc **********************************************************************
      SUBROUTINE filter_x_cp_4th(imax,jmax,kmax,U1,U2,U3,U4,U5)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U4(imax,jmax,kmax),U5(imax,jmax,kmax),
     &  AA(imax,imax),LL(imax,imax),UU(imax,imax),
     &  U1_f(imax),U2_f(imax),U3_f(imax),U4_f(imax),U5_f(imax),
     &  rhs1(imax),rhs2(imax),rhs3(imax),rhs4(imax),rhs5(imax)

      call lhsf_L(imax,AA)
      call ludecomp(0,imax,AA,LL,UU) 

      do k=1,kmax
        do j=1,jmax
          call rhsf_x_L(imax,jmax,kmax,j,k,U1,rhs1)
          call solve(0,imax,LL,UU,rhs1,U1_f)
          call rhsf_x_L(imax,jmax,kmax,j,k,U2,rhs2)
          call solve(0,imax,LL,UU,rhs2,U2_f)
          call rhsf_x_L(imax,jmax,kmax,j,k,U3,rhs3)
          call solve(0,imax,LL,UU,rhs3,U3_f)
          call rhsf_x_L(imax,jmax,kmax,j,k,U4,rhs4)
          call solve(0,imax,LL,UU,rhs4,U4_f)
          call rhsf_x_L(imax,jmax,kmax,j,k,U5,rhs5)
          call solve(0,imax,LL,UU,rhs5,U5_f)
         
          do i=1,imax
            U1(i,j,k)=U1_f(i)           
            U2(i,j,k)=U2_f(i)
            U3(i,j,k)=U3_f(i)
            U4(i,j,k)=U4_f(i)
            U5(i,j,k)=U5_f(i)
          end do
        end do
      end do   

      RETURN
      END SUBROUTINE filter_x_cp_4th

cccccc **********************************************************************
cccccc filter_x_ex_6th: Explicit filter of sixth order in x-direction for 
cccccc                  non-periodic boundary conditions.
cccccc **********************************************************************
      SUBROUTINE filter_x_ex_6th(imax,jmax,kmax,U1,U2,U3,U4,U5)

      USE snmpi, only: my_rank,pro

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 aa,bb,cc,dd,
     &  U1  (imax,jmax,kmax),U2  (imax,jmax,kmax),
     &  U3  (imax,jmax,kmax),U4  (imax,jmax,kmax),
     &  U5  (imax,jmax,kmax),U1_f(imax,jmax,kmax),
     &  U2_f(imax,jmax,kmax),U3_f(imax,jmax,kmax),
     &  U4_f(imax,jmax,kmax),U5_f(imax,jmax,kmax)

      aa = +11.d0/16.d0
      bb = +15.d0/32.d0
      cc = -03.d0/16.d0
      dd = +01.d0/32.d0

      if (my_rank.EQ.0) then
        do k=1,kmax
          do j=1,jmax            
            U1_f(1,j,k) = 15.d0/16.d0*U1(1,j,k)+1.d0/16.d0*(4.d0*U1(2,j,k)
     &        -6.d0*U1(3,j,k)+4.d0*U1(4,j,k)-U1(5,j,k))
            U2_f(1,j,k) = 15.d0/16.d0*U2(1,j,k)+1.d0/16.d0*(4.d0*U2(2,j,k)
     &        -6.d0*U2(3,j,k)+4.d0*U2(4,j,k)-U2(5,j,k))
            U3_f(1,j,k) = 15.d0/16.d0*U3(1,j,k)+1.d0/16.d0*(4.d0*U3(2,j,k)
     &        -6.d0*U3(3,j,k)+4.d0*U3(4,j,k)-U3(5,j,k))
            U4_f(1,j,k) = 15.d0/16.d0*U4(1,j,k)+1.d0/16.d0*(4.d0*U4(2,j,k)
     &        -6.d0*U4(3,j,k)+4.d0*U4(4,j,k)-U4(5,j,k))
            U5_f(1,j,k) = 15.d0/16.d0*U5(1,j,k)+1.d0/16.d0*(4.d0*U5(2,j,k)
     &        -6.d0*U5(3,j,k)+4.d0*U5(4,j,k)-U5(5,j,k))

            U1_f(2,j,k) = 3.d0/4.d0*U1(2,j,k)+1.d0/16.d0*(U1(1,j,k)
     &        +6.d0*U1(3,j,k)-4.d0*U1(4,j,k)+U1(5,j,k))
            U2_f(2,j,k) = 3.d0/4.d0*U2(2,j,k)+1.d0/16.d0*(U2(1,j,k)
     &        +6.d0*U2(3,j,k)-4.d0*U2(4,j,k)+U2(5,j,k))
            U3_f(2,j,k) = 3.d0/4.d0*U3(2,j,k)+1.d0/16.d0*(U3(1,j,k)
     &        +6.d0*U3(3,j,k)-4.d0*U3(4,j,k)+U3(5,j,k))
            U4_f(2,j,k) = 3.d0/4.d0*U4(2,j,k)+1.d0/16.d0*(U4(1,j,k)
     &        +6.d0*U4(3,j,k)-4.d0*U4(4,j,k)+U4(5,j,k))
            U5_f(2,j,k) = 3.d0/4.d0*U5(2,j,k)+1.d0/16.d0*(U5(1,j,k)
     &        +6.d0*U5(3,j,k)-4.d0*U5(4,j,k)+U5(5,j,k))

            U1_f(3,j,k) = 5.d0/8.d0*U1(3,j,k)+1.d0/16.d0*(-U1(1,j,k)
     &        +4.d0*U1(2,j,k)+4.d0*U1(4,j,k)-U1(5,j,k))
            U2_f(3,j,k) = 5.d0/8.d0*U2(3,j,k)+1.d0/16.d0*(-U2(1,j,k)
     &        +4.d0*U2(2,j,k)+4.d0*U2(4,j,k)-U2(5,j,k))
            U3_f(3,j,k) = 5.d0/8.d0*U3(3,j,k)+1.d0/16.d0*(-U3(1,j,k)
     &        +4.d0*U3(2,j,k)+4.d0*U3(4,j,k)-U3(5,j,k))
            U4_f(3,j,k) = 5.d0/8.d0*U4(3,j,k)+1.d0/16.d0*(-U4(1,j,k)
     &        +4.d0*U4(2,j,k)+4.d0*U4(4,j,k)-U4(5,j,k))
            U5_f(3,j,k) = 5.d0/8.d0*U5(3,j,k)+1.d0/16.d0*(-U5(1,j,k)
     &        +4.d0*U5(2,j,k)+4.d0*U5(4,j,k)-U5(5,j,k))

          do i=4,imax-3
              U1_f(i,j,k) = 
     &          +dd/2.d0*U1(i-3,j,k)+cc/2.d0*U1(i-2,j,k)+bb/2.d0*U1(i-1,j,k)+aa/1.d0*U1(i,j,k)
     &          +bb/2.d0*U1(i+1,j,k)+cc/2.d0*U1(i+2,j,k)+dd/2.d0*U1(i+3,j,k)
              U2_f(i,j,k) = 
     &          +dd/2.d0*U2(i-3,j,k)+cc/2.d0*U2(i-2,j,k)+bb/2.d0*U2(i-1,j,k)+aa/1.d0*U2(i,j,k)
     &          +bb/2.d0*U2(i+1,j,k)+cc/2.d0*U2(i+2,j,k)+dd/2.d0*U2(i+3,j,k)
              U3_f(i,j,k) = 
     &          +dd/2.d0*U3(i-3,j,k)+cc/2.d0*U3(i-2,j,k)+bb/2.d0*U3(i-1,j,k)+aa/1.d0*U3(i,j,k)
     &          +bb/2.d0*U3(i+1,j,k)+cc/2.d0*U3(i+2,j,k)+dd/2.d0*U3(i+3,j,k)
              U4_f(i,j,k) = 
     &          +dd/2.d0*U4(i-3,j,k)+cc/2.d0*U4(i-2,j,k)+bb/2.d0*U4(i-1,j,k)+aa/1.d0*U4(i,j,k)
     &          +bb/2.d0*U4(i+1,j,k)+cc/2.d0*U4(i+2,j,k)+dd/2.d0*U4(i+3,j,k)
              U5_f(i,j,k) = 
     &          +dd/2.d0*U5(i-3,j,k)+cc/2.d0*U5(i-2,j,k)+bb/2.d0*U5(i-1,j,k)+aa/1.d0*U5(i,j,k)
     &          +bb/2.d0*U5(i+1,j,k)+cc/2.d0*U5(i+2,j,k)+dd/2.d0*U5(i+3,j,k)
            end do
          end do
        end do
      end if

      if ((my_rank.NE.0).AND.(my_rank+1.NE.pro)) then
        do k=1,kmax
          do j=1,jmax            
            do i=4,imax-3
              U1_f(i,j,k) = 
     &          +dd/2.d0*U1(i-3,j,k)+cc/2.d0*U1(i-2,j,k)+bb/2.d0*U1(i-1,j,k)+aa/1.d0*U1(i,j,k)
     &          +bb/2.d0*U1(i+1,j,k)+cc/2.d0*U1(i+2,j,k)+dd/2.d0*U1(i+3,j,k)
              U2_f(i,j,k) = 
     &          +dd/2.d0*U2(i-3,j,k)+cc/2.d0*U2(i-2,j,k)+bb/2.d0*U2(i-1,j,k)+aa/1.d0*U2(i,j,k)
     &          +bb/2.d0*U2(i+1,j,k)+cc/2.d0*U2(i+2,j,k)+dd/2.d0*U2(i+3,j,k)
              U3_f(i,j,k) = 
     &          +dd/2.d0*U3(i-3,j,k)+cc/2.d0*U3(i-2,j,k)+bb/2.d0*U3(i-1,j,k)+aa/1.d0*U3(i,j,k)
     &          +bb/2.d0*U3(i+1,j,k)+cc/2.d0*U3(i+2,j,k)+dd/2.d0*U3(i+3,j,k)
              U4_f(i,j,k) = 
     &          +dd/2.d0*U4(i-3,j,k)+cc/2.d0*U4(i-2,j,k)+bb/2.d0*U4(i-1,j,k)+aa/1.d0*U4(i,j,k)
     &          +bb/2.d0*U4(i+1,j,k)+cc/2.d0*U4(i+2,j,k)+dd/2.d0*U4(i+3,j,k)
              U5_f(i,j,k) = 
     &          +dd/2.d0*U5(i-3,j,k)+cc/2.d0*U5(i-2,j,k)+bb/2.d0*U5(i-1,j,k)+aa/1.d0*U5(i,j,k)
     &          +bb/2.d0*U5(i+1,j,k)+cc/2.d0*U5(i+2,j,k)+dd/2.d0*U5(i+3,j,k)
            end do
          end do
        end do
      end if

      if (my_rank+1.EQ.pro) then
        do k=1,kmax
          do j=1,jmax            
            do i=4,imax-3
              U1_f(i,j,k) = 
     &          +dd/2.d0*U1(i-3,j,k)+cc/2.d0*U1(i-2,j,k)+bb/2.d0*U1(i-1,j,k)+aa/1.d0*U1(i,j,k)
     &          +bb/2.d0*U1(i+1,j,k)+cc/2.d0*U1(i+2,j,k)+dd/2.d0*U1(i+3,j,k)
              U2_f(i,j,k) = 
     &          +dd/2.d0*U2(i-3,j,k)+cc/2.d0*U2(i-2,j,k)+bb/2.d0*U2(i-1,j,k)+aa/1.d0*U2(i,j,k)
     &          +bb/2.d0*U2(i+1,j,k)+cc/2.d0*U2(i+2,j,k)+dd/2.d0*U2(i+3,j,k)
              U3_f(i,j,k) = 
     &          +dd/2.d0*U3(i-3,j,k)+cc/2.d0*U3(i-2,j,k)+bb/2.d0*U3(i-1,j,k)+aa/1.d0*U3(i,j,k)
     &          +bb/2.d0*U3(i+1,j,k)+cc/2.d0*U3(i+2,j,k)+dd/2.d0*U3(i+3,j,k)
              U4_f(i,j,k) = 
     &          +dd/2.d0*U4(i-3,j,k)+cc/2.d0*U4(i-2,j,k)+bb/2.d0*U4(i-1,j,k)+aa/1.d0*U4(i,j,k)
     &          +bb/2.d0*U4(i+1,j,k)+cc/2.d0*U4(i+2,j,k)+dd/2.d0*U4(i+3,j,k)
              U5_f(i,j,k) = 
     &          +dd/2.d0*U5(i-3,j,k)+cc/2.d0*U5(i-2,j,k)+bb/2.d0*U5(i-1,j,k)+aa/1.d0*U5(i,j,k)
     &          +bb/2.d0*U5(i+1,j,k)+cc/2.d0*U5(i+2,j,k)+dd/2.d0*U5(i+3,j,k)
            end do

            U1_f(imax-2,j,k) = 5.d0/8.d0*U1(imax-2,j,k)+1.d0/16.d0*(-U1(imax,j,k)
     &        +4.d0*U1(imax-1,j,k)+4.d0*U1(imax-3,j,k)-U1(imax-4,j,k))
            U2_f(imax-2,j,k) = 5.d0/8.d0*U2(imax-2,j,k)+1.d0/16.d0*(-U2(imax,j,k)
     &        +4.d0*U2(imax-1,j,k)+4.d0*U2(imax-3,j,k)-U2(imax-4,j,k))
            U3_f(imax-2,j,k) = 5.d0/8.d0*U3(imax-2,j,k)+1.d0/16.d0*(-U3(imax,j,k)
     &        +4.d0*U3(imax-1,j,k)+4.d0*U3(imax-3,j,k)-U3(imax-4,j,k))
            U4_f(imax-2,j,k) = 5.d0/8.d0*U4(imax-2,j,k)+1.d0/16.d0*(-U4(imax,j,k)
     &        +4.d0*U4(imax-1,j,k)+4.d0*U4(imax-3,j,k)-U4(imax-4,j,k))
            U5_f(imax-2,j,k) = 5.d0/8.d0*U5(imax-2,j,k)+1.d0/16.d0*(-U5(imax,j,k)
     &        +4.d0*U5(imax-1,j,k)+4.d0*U5(imax-3,j,k)-U5(imax-4,j,k))

            U1_f(imax-1,j,k) = 3.d0/4.d0*U1(imax-1,j,k)+1.d0/16.d0*(U1(imax,j,k)
     &        +6.d0*U1(imax-2,j,k)-4.d0*U1(imax-3,j,k)+U1(imax-4,j,k))
            U2_f(imax-1,j,k) = 3.d0/4.d0*U2(imax-1,j,k)+1.d0/16.d0*(U2(imax,j,k)
     &        +6.d0*U2(imax-2,j,k)-4.d0*U2(imax-3,j,k)+U2(imax-4,j,k))
            U3_f(imax-1,j,k) = 3.d0/4.d0*U3(imax-1,j,k)+1.d0/16.d0*(U3(imax,j,k)
     &        +6.d0*U3(imax-2,j,k)-4.d0*U3(imax-3,j,k)+U3(imax-4,j,k))
            U4_f(imax-1,j,k) = 3.d0/4.d0*U4(imax-1,j,k)+1.d0/16.d0*(U4(imax,j,k)
     &        +6.d0*U4(imax-2,j,k)-4.d0*U4(imax-3,j,k)+U4(imax-4,j,k))
            U5_f(imax-1,j,k) = 3.d0/4.d0*U5(imax-1,j,k)+1.d0/16.d0*(U5(imax,j,k)
     &        +6.d0*U5(imax-2,j,k)-4.d0*U5(imax-3,j,k)+U5(imax-4,j,k))

            U1_f(imax,j,k) = 15.d0/16.d0*U1(imax,j,k)+1.d0/16.d0*(4.d0*U1(imax-1,j,k)
     &        -6.d0*U1(imax-2,j,k)+4.d0*U1(imax-3,j,k)-U1(imax-4,j,k))
            U2_f(imax,j,k) = 15.d0/16.d0*U2(imax,j,k)+1.d0/16.d0*(4.d0*U2(imax-1,j,k)
     &        -6.d0*U2(imax-2,j,k)+4.d0*U2(imax-3,j,k)-U2(imax-4,j,k))
            U3_f(imax,j,k) = 15.d0/16.d0*U3(imax,j,k)+1.d0/16.d0*(4.d0*U3(imax-1,j,k)
     &        -6.d0*U3(imax-2,j,k)+4.d0*U3(imax-3,j,k)-U3(imax-4,j,k))
            U4_f(imax,j,k) = 15.d0/16.d0*U4(imax,j,k)+1.d0/16.d0*(4.d0*U4(imax-1,j,k)
     &        -6.d0*U4(imax-2,j,k)+4.d0*U4(imax-3,j,k)-U4(imax-4,j,k))
            U5_f(imax,j,k) = 15.d0/16.d0*U5(imax,j,k)+1.d0/16.d0*(4.d0*U5(imax-1,j,k)
     &        -6.d0*U5(imax-2,j,k)+4.d0*U5(imax-3,j,k)-U5(imax-4,j,k))
          end do
        end do
      end if

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            U1(i,j,k)=U1_f(i,j,k)
            U2(i,j,k)=U2_f(i,j,k)
            U3(i,j,k)=U3_f(i,j,k)
            U4(i,j,k)=U4_f(i,j,k)
            U5(i,j,k)=U5_f(i,j,k)
          end do
        end do 
      end do

      RETURN
      END SUBROUTINE filter_x_ex_6th

cccccc **********************************************************************
cccccc filter_y_cp_4th: Compact filter of fourth order in y-direction for 
cccccc                  non-periodic boundary conditions.
cccccc **********************************************************************
      SUBROUTINE filter_y_cp_4th(imax,jmax,kmax,U1,U2,U3,U4,U5)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U4(imax,jmax,kmax),U5(imax,jmax,kmax),
     &  AA(jmax,jmax),LL(jmax,jmax),UU(jmax,jmax),
     &  U1_f(jmax),U2_f(jmax),U3_f(jmax),U4_f(jmax),U5_f(jmax),
     &  rhs1(jmax),rhs2(jmax),rhs3(jmax),rhs4(jmax),rhs5(jmax)

      call lhsf_L(jmax,AA)
      call ludecomp(0,jmax,AA,LL,UU)   

      do k=1,kmax
        do i=1,imax
          call rhsf_y_L(imax,jmax,kmax,i,k,U1,rhs1)
          call solve(0,jmax,LL,UU,rhs1,U1_f)
          call rhsf_y_L(imax,jmax,kmax,i,k,U2,rhs2)
          call solve(0,jmax,LL,UU,rhs2,U2_f)
          call rhsf_y_L(imax,jmax,kmax,i,k,U3,rhs3)
          call solve(0,jmax,LL,UU,rhs3,U3_f)
          call rhsf_y_L(imax,jmax,kmax,i,k,U4,rhs4)
          call solve(0,jmax,LL,UU,rhs4,U4_f)
          call rhsf_y_L(imax,jmax,kmax,i,k,U5,rhs5)
          call solve(0,jmax,LL,UU,rhs5,U5_f)
          
          do j=1,jmax
            U1(i,j,k)=U1_f(j)           
            U2(i,j,k)=U2_f(j)
            U3(i,j,k)=U3_f(j)
            U4(i,j,k)=U4_f(j)
            U5(i,j,k)=U5_f(j)
          end do
        end do
      end do   

      RETURN
      END SUBROUTINE filter_y_cp_4th

cccccc **********************************************************************
cccccc filter_z_cp_4th_P: Compact filter of fourth order in z-direction for
cccccc                    periodic boundary conditions.
cccccc **********************************************************************
      SUBROUTINE filter_z_cp_4th_P(imax,jmax,kmax,U1,U2,U3,U4,U5)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 
     &  U1(imax,jmax,kmax),U2(imax,jmax,kmax),
     &  U3(imax,jmax,kmax),U4(imax,jmax,kmax),U5(imax,jmax,kmax),
     &  AA(kmax,kmax),LL(kmax,kmax),UU(kmax,kmax),
     &  U1_f(kmax),U2_f(kmax),U3_f(kmax),U4_f(kmax),U5_f(kmax),
     &  rhs1(kmax),rhs2(kmax),rhs3(kmax),rhs4(kmax),rhs5(kmax)

      call lhsf_P(kmax,AA)
      call ludecomp(1,kmax,AA,LL,UU)

      do j=1,jmax
        do i=1,imax
          call rhsf_z_P(imax,jmax,kmax,i,j,U1,rhs1)
          call solve(1,kmax,LL,UU,rhs1,U1_f)
          call rhsf_z_P(imax,jmax,kmax,i,j,U2,rhs2)
          call solve(1,kmax,LL,UU,rhs2,U2_f)
          call rhsf_z_P(imax,jmax,kmax,i,j,U3,rhs3)
          call solve(1,kmax,LL,UU,rhs3,U3_f)
          call rhsf_z_P(imax,jmax,kmax,i,j,U4,rhs4)
          call solve(1,kmax,LL,UU,rhs4,U4_f)
          call rhsf_z_P(imax,jmax,kmax,i,j,U5,rhs5)
          call solve(1,kmax,LL,UU,rhs5,U5_f)

          U1_f(kmax)=U1_f(1)
          U2_f(kmax)=U2_f(1)
          U3_f(kmax)=U3_f(1)
          U4_f(kmax)=U4_f(1)
          U5_f(kmax)=U5_f(1)

          do k=1,kmax
            U1(i,j,k)=U1_f(k)
            U2(i,j,k)=U2_f(k)
            U3(i,j,k)=U3_f(k)
            U4(i,j,k)=U4_f(k)
            U5(i,j,k)=U5_f(k)
          end do
        end do
      end do   

      RETURN
      END SUBROUTINE filter_z_cp_4th_P

cccccc **********************************************************************
cccccc lhsf_P: Left hand side of numerical filter in x-direction
cccccc         for periodic boundary condition.
cccccc **********************************************************************
      SUBROUTINE lhsf_P(lmax,aa)

      IMPLICIT NONE

      integer i,j,lmax,lmax_1
      real*8 aa(lmax,lmax),alpha,beta

      lmax_1=lmax-1
      alpha = 0.6522474d0 
      beta  = 0.1702929d0 

      do i=1,lmax
        do j=1,lmax
          aa(i,j)=0.d0
        end do
      end do  

      aa(1,lmax_1-1) = beta
      aa(1,lmax_1)   = alpha
      aa(1,1)        = 1.d0
      aa(1,2)        = alpha
      aa(1,3)        = beta

      aa(2,lmax_1)   = beta
      aa(2,1)        = alpha
      aa(2,2)        = 1.d0
      aa(2,3)        = alpha
      aa(2,4)        = beta

      do i=3,lmax_1-2
        aa(i,i-2) = beta
        aa(i,i-1) = alpha
        aa(i,i-0) = 1.d0
        aa(i,i+1) = alpha
        aa(i,i+2) = beta
      end do

      aa(lmax_1-1,lmax_1-3) = beta
      aa(lmax_1-1,lmax_1-2) = alpha
      aa(lmax_1-1,lmax_1-1) = 1.d0
      aa(lmax_1-1,lmax_1-0) = alpha
      aa(lmax_1-1,1)        = beta

      aa(lmax_1,lmax_1-2)  = beta
      aa(lmax_1,lmax_1-1)  = alpha
      aa(lmax_1,lmax_1-0)  = 1.d0
      aa(lmax_1,1)         = alpha
      aa(lmax_1,2)         = beta

      RETURN
      END SUBROUTINE lhsf_P

cccccc **********************************************************************
cccccc lhsf_L: Left hand side of numerical filter in y-direction
cccccc         for non-periodic boundary condition.
cccccc **********************************************************************
      SUBROUTINE lhsf_L(lmax,aa)

      IMPLICIT NONE

      integer i,j,l,lmax
      real*8 aa(lmax,lmax),alpha,beta

      alpha = 0.6522474d0
      beta  = 0.1702929d0

      do i=1,lmax
        do j=1,lmax
          aa(i,j)=0.d0
        end do
      end do  

      aa(1,1) = 1.d0
      aa(2,2) = 1.d0
      aa(3,3) = 1.d0

      do l=4,lmax-3
        aa(l,l-2) = beta
        aa(l,l-1) = alpha
        aa(l,l  ) = 1.d0
        aa(l,l+1) = alpha
        aa(l,l+2) = beta
      end do

      aa(lmax-2,lmax-2) = 1.d0
      aa(lmax-1,lmax-1) = 1.d0
      aa(lmax-0,lmax-0) = 1.d0

      RETURN
      END SUBROUTINE lhsf_L

cccccc **********************************************************************
cccccc rhsf_x_P: Right hand side of numerical filter in the x-direction
cccccc           for periodic boundary.
cccccc **********************************************************************
      SUBROUTINE rhsf_x_P(imax,jmax,kmax,j,k,E,rhs)

      IMPLICIT NONE

      integer i,j,k,imax,imax_1,jmax,kmax
      real*8 rhs(imax),E(imax,jmax,kmax),a,b,c,d

      a = 0.98918560d0 
      b = 1.32118000d0 
      c = 0.33335480d0
      d = 0.00135985d0 

      imax_1=imax-1
      rhs(1)=a*(E(1,j,k))+
     &  b/2.d0*(E(2,j,k)+E(imax_1,j,k))+
     &  c/2.d0*(E(3,j,k)+E(imax_1-1,j,k))+
     &  d/2.d0*(E(4,j,k)+E(imax_1-2,j,k))
      rhs(2)=a*(E(2,j,k))+
     &  b/2.d0*(E(3,j,k)+E(1,j,k))+
     &  c/2.d0*(E(4,j,k)+E(imax_1,j,k))+
     &  d/2.d0*(E(5,j,k)+E(imax_1-1,j,k))
      rhs(3)=a*(E(3,j,k))+
     &  b/2.d0*(E(4,j,k)+E(2,j,k))+
     &  c/2.d0*(E(5,j,k)+E(1,j,k))+
     &  d/2.d0*(E(6,j,k)+E(imax_1,j,k))
      do i=4,imax_1-3
        rhs(i)=a*(E(i,j,k))+
     &    b/2.d0*(E(i+1,j,k)+E(i-1,j,k))+
     &    c/2.d0*(E(i+2,j,k)+E(i-2,j,k))+
     &    d/2.d0*(E(i+3,j,k)+E(i-3,j,k))
      end do
      rhs(imax_1-2)=a*(E(imax_1-2,j,k))+
     &  b/2.d0*(E(imax_1-1,j,k)+E(imax_1-3,j,k))+
     &  c/2.d0*(E(imax_1,j,k)  +E(imax_1-4,j,k))+
     &  d/2.d0*(E(1,j,k)       +E(imax_1-5,j,k))
      rhs(imax_1-1)=a*(E(imax_1-1,j,k))+
     &  b/2.d0*(E(imax_1,j,k)+E(imax_1-2,j,k))+
     &  c/2.d0*(E(1,j,k)     +E(imax_1-3,j,k))+
     &  d/2.d0*(E(2,j,k)     +E(imax_1-4,j,k))
      rhs(imax_1)=a*(E(imax_1,j,k))+
     &  b/2.d0*(E(1,j,k)     +E(imax_1-1,j,k))+
     &  c/2.d0*(E(2,j,k)     +E(imax_1-2,j,k))+
     &  d/2.d0*(E(3,j,k)     +E(imax_1-3,j,k))
      rhs(imax)=0.d0

      RETURN
      END SUBROUTINE rhsf_x_P

cccccc **********************************************************************
cccccc rhsf_x_L: Right hand side of numerical filter in the x-direction
cccccc           for non-periodic boundary condition.
cccccc **********************************************************************
      SUBROUTINE rhsf_x_L(imax,jmax,kmax,j,k,E,rhs)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 rhs(imax),E(imax,jmax,kmax),a,b,c,d,alpha,beta

      alpha = 0.6522474d0
      beta  = 0.1702929d0

      d     = 0.00135985d0
      c     = -1.d0/8.d0*(1.d0-2.d0*alpha-14.d0*beta+16.d0*d)
      b     = +1.d0/2.d0*(1.d0+2.d0*alpha+2.d0*beta-2.d0*d  )
      a     = +1.d0/8.d0*(5.d0+6.d0*alpha-6.d0*beta+16.d0*d )

      rhs(1)=15.d0/16.d0*E(1,j,k)+1.d0/16.d0*(4.d0*E(2,j,k)
     &  -6.d0*E(3,j,k)+4.d0*E(4,j,k)-E(5,j,k))
      rhs(2)=3.d0/4.d0*E(2,j,k)+1.d0/16.d0*(E(1,j,k)
     &  +6.d0*E(3,j,k)-4.d0*E(4,j,k)+E(5,j,k))
      rhs(3)=5.d0/8.d0*E(3,j,k)+1.d0/16.d0*(-E(1,j,k)
     &  +4.d0*E(2,j,k)+4.d0*E(4,j,k)-E(5,j,k))
      do i=4,imax-3
        rhs(i)=
     &    a*(E(i,j,k))+
     &    b/2.d0*(E(i+1,j,k)+E(i-1,j,k))+
     &    c/2.d0*(E(i+2,j,k)+E(i-2,j,k))+
     &    d/2.d0*(E(i+3,j,k)+E(i-3,j,k))
      end do
      rhs(imax-2)=5.d0/8.d0*E(imax-2,j,k)+1.d0/16.d0*(-E(imax,j,k)
     &  +4.d0*E(imax-1,j,k)+4.d0*E(imax-3,j,k)-E(imax-4,j,k) )
      rhs(imax-1)=3.d0/4.d0*E(imax-1,j,k)+1.d0/16.d0*(E(imax,j,k)
     &  +6.d0*E(imax-2,j,k)-4.d0*E(imax-3,j,k)+E(imax-4,j,k) )
      rhs(imax)  =15.d0/16.d0*E(imax,j,k)+1.d0/16.d0*(4.d0*E(imax-1,j,k)
     &  -6.d0*E(imax-2,j,k)+4.d0*E(imax-3,j,k)-E(imax-4,j,k) )

      RETURN
      END SUBROUTINE rhsf_x_L

cccccc **********************************************************************
cccccc rhsf_y_L: Right hand side of numerical filter in the y-direction
cccccc           for non-periodic boundary condition.
cccccc **********************************************************************
      SUBROUTINE rhsf_y_L(imax,jmax,kmax,i,k,E,rhs)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 rhs(jmax),E(imax,jmax,kmax),a,b,c,d,alpha,beta

      alpha = 0.6522474d0
      beta  = 0.1702929d0

      d     = 0.00135985d0
      c     = -1.d0/8.d0*(1.d0-2.d0*alpha-14.d0*beta+16.d0*d)
      b     = +1.d0/2.d0*(1.d0+2.d0*alpha+02.d0*beta-02.d0*d)
      a     = +1.d0/8.d0*(5.d0+6.d0*alpha-06.d0*beta+16.d0*d)

      rhs(1)=15.d0/16.d0*E(i,1,k)+1.d0/16.d0*(4.d0*E(i,2,k)
     &  -6.d0*E(i,3,k)+4.d0*E(i,4,k)-E(i,5,k) )
      rhs(2)=3.d0/4.d0*E(i,2,k)+1.d0/16.d0*(E(i,1,k)
     &  +6.d0*E(i,3,k)-4.d0*E(i,4,k)+E(i,5,k) )
      rhs(3)=5.d0/8.d0*E(i,3,k)+1.d0/16.d0*(-E(i,1,k)
     &  +4.d0*E(i,2,k)+4.d0*E(i,4,k)-E(i,5,k) )
      do j=4,jmax-3
        rhs(j)=
     &    a*(E(i,j,k))+
     &    b/2.d0*(E(i,j+1,k)+E(i,j-1,k))+
     &    c/2.d0*(E(i,j+2,k)+E(i,j-2,k))+
     &    d/2.d0*(E(i,j+3,k)+E(i,j-3,k))
      end do
      rhs(jmax-2)=5.d0/8.d0*E(i,jmax-2,k)+1.d0/16.d0*(-E(i,jmax,k)
     &  +4.d0*E(i,jmax-1,k)+4.d0*E(i,jmax-3,k)-E(i,jmax-4,k))
      rhs(jmax-1)=3.d0/4.d0*E(i,jmax-1,k)+1.d0/16.d0*(E(i,jmax,k)
     &  +6.d0*E(i,jmax-2,k)-4.d0*E(i,jmax-3,k)+E(i,jmax-4,k))
      rhs(jmax)  =15.d0/16.d0*E(i,jmax,k)+1.d0/16.d0*(4.d0*E(i,jmax-1,k)
     &  -6.d0*E(i,jmax-2,k)+4.d0*E(i,jmax-3,k)-E(i,jmax-4,k))

      RETURN
      END SUBROUTINE rhsf_y_L

cccccc **********************************************************************
cccccc rhsf_z_P: Right hand side of numerical filter in the z-direction.
cccccc **********************************************************************
      SUBROUTINE rhsf_z_P(imax,jmax,kmax,i,j,E,rhs)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax,kmax_1
      real*8 rhs(kmax),E(imax,jmax,kmax),a,b,c,d

      a = 0.98918560d0 
      b = 1.32118000d0 
      c = 0.33335480d0
      d = 0.00135985d0 

      kmax_1=kmax-1
      rhs(1)=a*(E(i,j,1))+
     &  b/2.d0*(E(i,j,2)+E(i,j,kmax_1))+
     &  c/2.d0*(E(i,j,3)+E(i,j,kmax_1-1))+
     &  d/2.d0*(E(i,j,4)+E(i,j,kmax_1-2))
      rhs(2)=a*(E(i,j,2))+
     &  b/2.d0*(E(i,j,3)+E(i,j,1))+
     &  c/2.d0*(E(i,j,4)+E(i,j,kmax_1))+
     &  d/2.d0*(E(i,j,5)+E(i,j,kmax_1-1))
      rhs(3)=a*(E(i,j,3))+
     &  b/2.d0*(E(i,j,4)+E(i,j,2))+
     &  c/2.d0*(E(i,j,5)+E(i,j,1))+
     &  d/2.d0*(E(i,j,6)+E(i,j,kmax_1))
      do k=4,kmax_1-3
        rhs(k)=a*(E(i,j,k))+
     &    b/2.d0*(E(i,j,k+1)+E(i,j,k-1))+
     &    c/2.d0*(E(i,j,k+2)+E(i,j,k-2))+
     &    d/2.d0*(E(i,j,k+3)+E(i,j,k-3))
      end do
      rhs(kmax_1-2)=a*(E(i,j,kmax_1-2))+
     &  b/2.d0*(E(i,j,kmax_1-1)+E(i,j,kmax_1-3))+
     &  c/2.d0*(E(i,j,kmax_1)  +E(i,j,kmax_1-4))+
     &  d/2.d0*(E(i,j,1)       +E(i,j,kmax_1-5))
      rhs(kmax_1-1)=a*(E(i,j,kmax_1-1))+
     &  b/2.d0*(E(i,j,kmax_1)  +E(i,j,kmax_1-2))+
     &  c/2.d0*(E(i,j,1)       +E(i,j,kmax_1-3))+
     &  d/2.d0*(E(i,j,2)       +E(i,j,kmax_1-4))
      rhs(kmax_1)=a*(E(i,j,kmax_1))+
     &  b/2.d0*(E(i,j,1)       +E(i,j,kmax_1-1))+
     &  c/2.d0*(E(i,j,2)       +E(i,j,kmax_1-2))+
     &  d/2.d0*(E(i,j,3)       +E(i,j,kmax_1-3))
      rhs(kmax)=0.d0

      RETURN
      END SUBROUTINE rhsf_z_P

cccccc **********************************************************************
cccccc ludecomp: This routine apply LU decomposition.
cccccc **********************************************************************
      SUBROUTINE ludecomp(fct,lmax,a,LL,UU)

      IMPLICIT NONE

      integer i,j,k,it,lmax,fct
      real*8 a(lmax,lmax),LL(lmax,lmax),UU(lmax,lmax),soma

      do i=1,lmax
        do j=1,lmax
          LL(i,j)=0.d0 !zera a matriz L
          UU(i,j)=0.d0 !zera a matriz U
        end do
      end do

      !Convert from A to LU
      do it=1,lmax-fct
        i=it
        do j=it,lmax-fct
          soma=0.d0
          do k=1,i-1
           soma=soma+LL(i,k)*UU(k,j)
          end do
          UU(i,j)=a(i,j)-soma
        end do
        j=it
        do i=it,lmax-fct
          soma=0.d0
          do k=1,j-1
           soma=soma+LL(i,k)*UU(k,j)
          end do
          LL(i,j)=(a(i,j)-soma)/UU(j,j)
        end do
      end do
      do i=1,lmax-fct
        LL(i,i)=1.d0
      end do

      RETURN
      END SUBROUTINE ludecomp

cccccc **********************************************************************
cccccc solve: This routine solve the linear system with LU decomposition.
cccccc **********************************************************************
      SUBROUTINE solve(fct,lmax,LL,UU,b,x)

      IMPLICIT NONE

      integer i,j,lmax,fct
      real*8 LL(lmax,lmax),UU(lmax,lmax),
     &       b(lmax),x(lmax),y(lmax),soma
 
      !Resolve a triangular inferior
      y(1)=b(1)/LL(1,1)
      do i=2,lmax-fct
        soma=0.d0
        do j=1,i-1
          soma=soma+LL(i,j)*y(j)
        end do
        y(i)=(b(i)-soma)/LL(i,i)
      end do

      !Resolve a triangular superior
      x(lmax-fct)=y(lmax-fct)/UU(lmax-fct,lmax-fct)
      do i=lmax-fct-1,1,-1
        soma=0.d0
        do j=i+1,lmax-fct
          soma=soma+UU(i,j)*x(j)
        end do
        x(i)=(y(i)-soma)/UU(i,i)
      end do

      RETURN
      END SUBROUTINE solve

      END MODULE snfilter    
