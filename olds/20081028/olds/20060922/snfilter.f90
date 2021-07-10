      MODULE snfilter
 
      contains
ccccc **********************************************************************
ccccc filter_x: Numerical filter.
ccccc **********************************************************************
      SUBROUTINE filter_x(local_imax,E1,E2,E3,E4,E5)
           
      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer local_imax
      real*8 
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax)

      !MS$IF (filter_in_x.EQ.1)
        !MS$IF (tp_s.EQ.1).OR.(tp_s.EQ.2).OR.(tp_s.EQ.3)
          call filter_x_P_4th(imax,jmax,kmax,E1,E2,E3,E4,E5)        
        !MS$ENDIF
        !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5).OR.(tp_s.EQ.6)
          call filter_x_L_6th(local_imax,jmax,kmax,E1,E2,E3,E4,E5)                             
        !MS$ENDIF
      !MS$ENDIF

      RETURN
      END SUBROUTINE filter_x

ccccc **********************************************************************
ccccc filter_y: Numerical filter.
ccccc **********************************************************************
      SUBROUTINE filter_y(local_imax,E1,E2,E3,E4,E5)
           
      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer local_imax
      real*8 
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax)

      !MS$IF (filter_in_y.EQ.1)
        !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5).OR.(tp_s.EQ.6)
          call filter_y_L_4th(local_imax,jmax,kmax,E1,E2,E3,E4,E5)
        !MS$ENDIF
      !MS$ENDIF

      RETURN
      END SUBROUTINE filter_y

ccccc **********************************************************************
ccccc filter_z: Numerical filter.
ccccc **********************************************************************
      SUBROUTINE filter_z(local_imax,E1,E2,E3,E4,E5)
           
      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer local_imax
      real*8 
     &  E1(local_imax,jmax,kmax),E2(local_imax,jmax,kmax),
     &  E3(local_imax,jmax,kmax),E4(local_imax,jmax,kmax),
     &  E5(local_imax,jmax,kmax)

      !MS$IF (filter_in_z.EQ.1)
        call filter_z_P_4th(local_imax,jmax,kmax,E1,E2,E3,E4,E5)
      !MS$ENDIF

      RETURN
      END SUBROUTINE filter_z

cccccc**********************************************************************
cccccc filter_x_P_4th: Numerical filter in x-direction for periodic 
cccccc                 boundary condition.
cccccc **********************************************************************
      SUBROUTINE filter_x_P_4th(imax,jmax,kmax,E1,E2,E3,E4,E5)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 
     &  E1(imax,jmax,kmax),E2(imax,jmax,kmax),
     &  E3(imax,jmax,kmax),E4(imax,jmax,kmax),E5(imax,jmax,kmax),
     &  AA(imax,imax),LL(imax,imax),UU(imax,imax),  
     &  E1_f(imax),E2_f(imax),E3_f(imax),E4_f(imax),E5_f(imax),
     &  rhs1(imax),rhs2(imax),rhs3(imax),rhs4(imax),rhs5(imax)

      call lhsf_P(imax,AA)
      call ludecomp(1,imax,AA,LL,UU)

      do k=1,kmax
        do j=1,jmax
          call rhsf_x_P(imax,jmax,kmax,j,k,E1,rhs1)
          call solve(1,imax,LL,UU,rhs1,E1_f)
          call rhsf_x_P(imax,jmax,kmax,j,k,E2,rhs2)
          call solve(1,imax,LL,UU,rhs2,E2_f)
          call rhsf_x_P(imax,jmax,kmax,j,k,E3,rhs3)
          call solve(1,imax,LL,UU,rhs3,E3_f)
          call rhsf_x_P(imax,jmax,kmax,j,k,E4,rhs4)
          call solve(1,imax,LL,UU,rhs4,E4_f)
          call rhsf_x_P(imax,jmax,kmax,j,k,E5,rhs5)
          call solve(1,imax,LL,UU,rhs5,E5_f)

          E1_f(imax)=E1_f(1)
          E2_f(imax)=E2_f(1)
          E3_f(imax)=E3_f(1)
          E4_f(imax)=E4_f(1)
          E5_f(imax)=E5_f(1)

          do i=1,imax
            E1(i,j,k)=E1_f(i)
            E2(i,j,k)=E2_f(i)
            E3(i,j,k)=E3_f(i)
            E4(i,j,k)=E4_f(i)
            E5(i,j,k)=E5_f(i)
          end do
        end do
      end do   

      RETURN
      END SUBROUTINE filter_x_P_4th

cccccc **********************************************************************
cccccc filter_x_L_4th: Numerical filter in x-direction for non-periodic 
cccccc                 boundary condition.
cccccc **********************************************************************
      SUBROUTINE filter_x_L_4th(imax,jmax,kmax,E1,E2,E3,E4,E5)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 
     &  E1(imax,jmax,kmax),E2(imax,jmax,kmax),
     &  E3(imax,jmax,kmax),E4(imax,jmax,kmax),E5(imax,jmax,kmax),
     &  AA(imax,imax),LL(imax,imax),UU(imax,imax),
     &  E1_f(imax),E2_f(imax),E3_f(imax),E4_f(imax),E5_f(imax),
     &  rhs1(imax),rhs2(imax),rhs3(imax),rhs4(imax),rhs5(imax)

      call lhsf_L(imax,AA)
      call ludecomp(0,imax,AA,LL,UU) 

      do k=1,kmax
        do j=1,jmax
          call rhsf_x_L(imax,jmax,kmax,j,k,E1,rhs1)
          call solve(0,imax,LL,UU,rhs1,E1_f)
          call rhsf_x_L(imax,jmax,kmax,j,k,E2,rhs2)
          call solve(0,imax,LL,UU,rhs2,E2_f)
          call rhsf_x_L(imax,jmax,kmax,j,k,E3,rhs3)
          call solve(0,imax,LL,UU,rhs3,E3_f)
          call rhsf_x_L(imax,jmax,kmax,j,k,E4,rhs4)
          call solve(0,imax,LL,UU,rhs4,E4_f)
          call rhsf_x_L(imax,jmax,kmax,j,k,E5,rhs5)
          call solve(0,imax,LL,UU,rhs5,E5_f)
         
          do i=1,imax
            E1(i,j,k)=E1_f(i)           
            E2(i,j,k)=E2_f(i)
            E3(i,j,k)=E3_f(i)
            E4(i,j,k)=E4_f(i)
            E5(i,j,k)=E5_f(i)
          end do
        end do
      end do   

      RETURN
      END SUBROUTINE filter_x_L_4th

cccccc **********************************************************************
cccccc filter_x_L_6th: Filter in x-direction for stretching grid.
cccccc **********************************************************************
      SUBROUTINE filter_x_L_6th(imax,jmax,kmax,E1,E2,E3,E4,E5)

      USE snmpi, only: my_rank,pro

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 aa,bb,cc,dd,
     &  E1  (imax,jmax,kmax),E2  (imax,jmax,kmax),
     &  E3  (imax,jmax,kmax),E4  (imax,jmax,kmax),
     &  E5  (imax,jmax,kmax),E1_f(imax,jmax,kmax),
     &  E2_f(imax,jmax,kmax),E3_f(imax,jmax,kmax),
     &  E4_f(imax,jmax,kmax),E5_f(imax,jmax,kmax)

      aa = +11.d0/16.d0
      bb = +15.d0/32.d0
      cc = -03.d0/16.d0
      dd = +01.d0/32.d0

      do k=1,kmax
        do j=1,jmax            

          if (my_rank.EQ.0) then
            E1_f(1,j,k) = 15.d0/16.d0*E1(1,j,k)+1.d0/16.d0*(4.d0*E1(2,j,k)
     &        -6.d0*E1(3,j,k)+4.d0*E1(4,j,k)-E1(5,j,k))
            E2_f(1,j,k) = 15.d0/16.d0*E2(1,j,k)+1.d0/16.d0*(4.d0*E2(2,j,k)
     &        -6.d0*E2(3,j,k)+4.d0*E2(4,j,k)-E2(5,j,k))
            E3_f(1,j,k) = 15.d0/16.d0*E3(1,j,k)+1.d0/16.d0*(4.d0*E3(2,j,k)
     &        -6.d0*E3(3,j,k)+4.d0*E3(4,j,k)-E3(5,j,k))
            E4_f(1,j,k) = 15.d0/16.d0*E4(1,j,k)+1.d0/16.d0*(4.d0*E4(2,j,k)
     &        -6.d0*E4(3,j,k)+4.d0*E4(4,j,k)-E4(5,j,k))
            E5_f(1,j,k) = 15.d0/16.d0*E5(1,j,k)+1.d0/16.d0*(4.d0*E5(2,j,k)
     &        -6.d0*E5(3,j,k)+4.d0*E5(4,j,k)-E5(5,j,k))

            E1_f(2,j,k) = 3.d0/4.d0*E1(2,j,k)+1.d0/16.d0*(E1(1,j,k)
     &        +6.d0*E1(3,j,k)-4.d0*E1(4,j,k)+E1(5,j,k))
            E2_f(2,j,k) = 3.d0/4.d0*E2(2,j,k)+1.d0/16.d0*(E2(1,j,k)
     &        +6.d0*E2(3,j,k)-4.d0*E2(4,j,k)+E2(5,j,k))
            E3_f(2,j,k) = 3.d0/4.d0*E3(2,j,k)+1.d0/16.d0*(E3(1,j,k)
     &        +6.d0*E3(3,j,k)-4.d0*E3(4,j,k)+E3(5,j,k))
            E4_f(2,j,k) = 3.d0/4.d0*E4(2,j,k)+1.d0/16.d0*(E4(1,j,k)
     &        +6.d0*E4(3,j,k)-4.d0*E4(4,j,k)+E4(5,j,k))
            E5_f(2,j,k) = 3.d0/4.d0*E5(2,j,k)+1.d0/16.d0*(E5(1,j,k)
     &        +6.d0*E5(3,j,k)-4.d0*E5(4,j,k)+E5(5,j,k))

            E1_f(3,j,k) = 5.d0/8.d0*E1(3,j,k)+1.d0/16.d0*(-E1(1,j,k)
     &        +4.d0*E1(2,j,k)+4.d0*E1(4,j,k)-E1(5,j,k))
            E2_f(3,j,k) = 5.d0/8.d0*E2(3,j,k)+1.d0/16.d0*(-E2(1,j,k)
     &        +4.d0*E2(2,j,k)+4.d0*E2(4,j,k)-E2(5,j,k))
            E3_f(3,j,k) = 5.d0/8.d0*E3(3,j,k)+1.d0/16.d0*(-E3(1,j,k)
     &        +4.d0*E3(2,j,k)+4.d0*E3(4,j,k)-E3(5,j,k))
            E4_f(3,j,k) = 5.d0/8.d0*E4(3,j,k)+1.d0/16.d0*(-E4(1,j,k)
     &        +4.d0*E4(2,j,k)+4.d0*E4(4,j,k)-E4(5,j,k))
            E5_f(3,j,k) = 5.d0/8.d0*E5(3,j,k)+1.d0/16.d0*(-E5(1,j,k)
     &        +4.d0*E5(2,j,k)+4.d0*E5(4,j,k)-E5(5,j,k))
          end if

          do i=4,imax-3
            E1_f(i,j,k) = 
     &        +dd/2.d0*E1(i-3,j,k)
     &        +cc/2.d0*E1(i-2,j,k)
     &        +bb/2.d0*E1(i-1,j,k)
     &        +aa/1.d0*E1(i  ,j,k)
     &        +bb/2.d0*E1(i+1,j,k)
     &        +cc/2.d0*E1(i+2,j,k)
     &        +dd/2.d0*E1(i+3,j,k)

            E2_f(i,j,k) = 
     &        +dd/2.d0*E2(i-3,j,k)
     &        +cc/2.d0*E2(i-2,j,k)
     &        +bb/2.d0*E2(i-1,j,k)
     &        +aa/1.d0*E2(i  ,j,k)
     &        +bb/2.d0*E2(i+1,j,k)
     &        +cc/2.d0*E2(i+2,j,k)
     &        +dd/2.d0*E2(i+3,j,k)

            E3_f(i,j,k) = 
     &        +dd/2.d0*E3(i-3,j,k)
     &        +cc/2.d0*E3(i-2,j,k)
     &        +bb/2.d0*E3(i-1,j,k)
     &        +aa/1.d0*E3(i  ,j,k)
     &        +bb/2.d0*E3(i+1,j,k)
     &        +cc/2.d0*E3(i+2,j,k)
     &        +dd/2.d0*E3(i+3,j,k)

            E4_f(i,j,k) = 
     &        +dd/2.d0*E4(i-3,j,k)
     &        +cc/2.d0*E4(i-2,j,k)
     &        +bb/2.d0*E4(i-1,j,k)
     &        +aa/1.d0*E4(i  ,j,k)
     &        +bb/2.d0*E4(i+1,j,k)
     &        +cc/2.d0*E4(i+2,j,k)
     &        +dd/2.d0*E4(i+3,j,k)

            E5_f(i,j,k) = 
     &        +dd/2.d0*E5(i-3,j,k)
     &        +cc/2.d0*E5(i-2,j,k)
     &        +bb/2.d0*E5(i-1,j,k)
     &        +aa/1.d0*E5(i  ,j,k)
     &        +bb/2.d0*E5(i+1,j,k)
     &        +cc/2.d0*E5(i+2,j,k)
     &        +dd/2.d0*E5(i+3,j,k)
          end do

          if (my_rank+1.EQ.pro) then
            E1_f(imax-2,j,k) = 5.d0/8.d0*E1(imax-2,j,k)+1.d0/16.d0*
     &        (-E1(imax,j,k)+4.d0*E1(imax-1,j,k)+4.d0*E1(imax-3,j,k)
     &        -E1(imax-4,j,k))
            E2_f(imax-2,j,k) = 5.d0/8.d0*E2(imax-2,j,k)+1.d0/16.d0*
     &        (-E2(imax,j,k)+4.d0*E2(imax-1,j,k)+4.d0*E2(imax-3,j,k)
     &        -E2(imax-4,j,k))
            E3_f(imax-2,j,k) = 5.d0/8.d0*E3(imax-2,j,k)+1.d0/16.d0*
     &        (-E3(imax,j,k)+4.d0*E3(imax-1,j,k)+4.d0*E3(imax-3,j,k)
     &        -E3(imax-4,j,k))
            E4_f(imax-2,j,k) = 5.d0/8.d0*E4(imax-2,j,k)+1.d0/16.d0*
     &        (-E4(imax,j,k)+4.d0*E4(imax-1,j,k)+4.d0*E4(imax-3,j,k)
     &        -E4(imax-4,j,k))
            E5_f(imax-2,j,k) = 5.d0/8.d0*E5(imax-2,j,k)+1.d0/16.d0*
     &        (-E5(imax,j,k)+4.d0*E5(imax-1,j,k)+4.d0*E5(imax-3,j,k)
     &        -E5(imax-4,j,k))

            E1_f(imax-1,j,k) = 3.d0/4.d0*E1(imax-1,j,k)+1.d0/16.d0*
     &        (E1(imax,j,k)+6.d0*E1(imax-2,j,k)-4.d0*E1(imax-3,j,k)
     &        +E1(imax-4,j,k) )
            E2_f(imax-1,j,k) = 3.d0/4.d0*E2(imax-1,j,k)+1.d0/16.d0*
     &        (E2(imax,j,k)+6.d0*E2(imax-2,j,k)-4.d0*E2(imax-3,j,k)
     &        +E2(imax-4,j,k) )
            E3_f(imax-1,j,k) = 3.d0/4.d0*E3(imax-1,j,k)+1.d0/16.d0*
     &        (E3(imax,j,k)+6.d0*E3(imax-2,j,k)-4.d0*E3(imax-3,j,k)
     &        +E3(imax-4,j,k) )
            E4_f(imax-1,j,k) = 3.d0/4.d0*E4(imax-1,j,k)+1.d0/16.d0*
     &        (E4(imax,j,k)+6.d0*E4(imax-2,j,k)-4.d0*E4(imax-3,j,k)
     &        +E4(imax-4,j,k) )
            E5_f(imax-1,j,k) = 3.d0/4.d0*E5(imax-1,j,k)+1.d0/16.d0*
     &        (E5(imax,j,k)+6.d0*E5(imax-2,j,k)-4.d0*E5(imax-3,j,k)
     &        +E5(imax-4,j,k) )

            E1_f(imax,j,k) = 15.d0/16.d0*E1(imax,j,k)+1.d0/16.d0*
     &        (4.d0*E1(imax-1,j,k)-6.d0*E1(imax-2,j,k)
     &        +4.d0*E1(imax-3,j,k)-E1(imax-4,j,k) )
            E2_f(imax,j,k) = 15.d0/16.d0*E2(imax,j,k)+1.d0/16.d0*
     &        (4.d0*E2(imax-1,j,k)-6.d0*E2(imax-2,j,k)
     &        +4.d0*E2(imax-3,j,k)-E2(imax-4,j,k) )
            E3_f(imax,j,k) = 15.d0/16.d0*E3(imax,j,k)+1.d0/16.d0*
     &        (4.d0*E3(imax-1,j,k)-6.d0*E3(imax-2,j,k)
     &        +4.d0*E3(imax-3,j,k)-E3(imax-4,j,k) )
            E4_f(imax,j,k) = 15.d0/16.d0*E4(imax,j,k)+1.d0/16.d0*
     &        (4.d0*E4(imax-1,j,k)-6.d0*E4(imax-2,j,k)
     &        +4.d0*E4(imax-3,j,k)-E4(imax-4,j,k) )
            E5_f(imax,j,k) = 15.d0/16.d0*E5(imax,j,k)+1.d0/16.d0*
     &        (4.d0*E5(imax-1,j,k)-6.d0*E5(imax-2,j,k)
     &        +4.d0*E5(imax-3,j,k)-E5(imax-4,j,k) )
          end if
        end do
      end do

      if (my_rank.EQ.0) then 
        do k=1,kmax
          do j=1,jmax
            do i=1,imax-6
              E1(i,j,k)=E1_f(i,j,k)
              E2(i,j,k)=E2_f(i,j,k)
              E3(i,j,k)=E3_f(i,j,k)
              E4(i,j,k)=E4_f(i,j,k)
              E5(i,j,k)=E5_f(i,j,k)
            end do
          end do 
        end do
      end if

      if ((my_rank.NE.0).AND.(my_rank+1.NE.pro)) then
        do k=1,kmax
          do j=1,jmax
            do i=7,imax-6
              E1(i,j,k)=E1_f(i,j,k)
              E2(i,j,k)=E2_f(i,j,k)
              E3(i,j,k)=E3_f(i,j,k)
              E4(i,j,k)=E4_f(i,j,k)
              E5(i,j,k)=E5_f(i,j,k)
            end do
          end do 
        end do
      end if 

      if (my_rank+1.EQ.pro) then
        do k=1,kmax
          do j=1,jmax
            do i=7,imax
              E1(i,j,k)=E1_f(i,j,k)
              E2(i,j,k)=E2_f(i,j,k)
              E3(i,j,k)=E3_f(i,j,k)
              E4(i,j,k)=E4_f(i,j,k)
              E5(i,j,k)=E5_f(i,j,k)
            end do
          end do 
        end do
      end if
 
      RETURN
      END SUBROUTINE filter_x_L_6th

cccccc **********************************************************************
cccccc filter_y_L_4th: Numerical filter in y-direction for non-periodic
cccccc                 boundary condition.
cccccc **********************************************************************
      SUBROUTINE filter_y_L_4th(imax,jmax,kmax,E1,E2,E3,E4,E5)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 
     &  E1(imax,jmax,kmax),E2(imax,jmax,kmax),
     &  E3(imax,jmax,kmax),E4(imax,jmax,kmax),E5(imax,jmax,kmax),
     &  AA(jmax,jmax),LL(jmax,jmax),UU(jmax,jmax),
     &  E1_f(jmax),E2_f(jmax),E3_f(jmax),E4_f(jmax),E5_f(jmax),
     &  rhs1(jmax),rhs2(jmax),rhs3(jmax),rhs4(jmax),rhs5(jmax)

      call lhsf_L(jmax,AA)
      call ludecomp(0,jmax,AA,LL,UU)   

      do k=1,kmax
        do i=1,imax
          call rhsf_y_L(imax,jmax,kmax,i,k,E1,rhs1)
          call solve(0,jmax,LL,UU,rhs1,E1_f)
          call rhsf_y_L(imax,jmax,kmax,i,k,E2,rhs2)
          call solve(0,jmax,LL,UU,rhs2,E2_f)
          call rhsf_y_L(imax,jmax,kmax,i,k,E3,rhs3)
          call solve(0,jmax,LL,UU,rhs3,E3_f)
          call rhsf_y_L(imax,jmax,kmax,i,k,E4,rhs4)
          call solve(0,jmax,LL,UU,rhs4,E4_f)
          call rhsf_y_L(imax,jmax,kmax,i,k,E5,rhs5)
          call solve(0,jmax,LL,UU,rhs5,E5_f)
          
          do j=1,jmax
            E1(i,j,k)=E1_f(j)           
            E2(i,j,k)=E2_f(j)
            E3(i,j,k)=E3_f(j)
            E4(i,j,k)=E4_f(j)
            E5(i,j,k)=E5_f(j)
          end do
        end do
      end do   

      RETURN
      END SUBROUTINE filter_y_L_4th

cccccc **********************************************************************
cccccc filter_z_P_4th: Numerical filter in z-direction.
cccccc **********************************************************************
      SUBROUTINE filter_z_P_4th(imax,jmax,kmax,E1,E2,E3,E4,E5)

      IMPLICIT NONE

      integer i,j,k,imax,jmax,kmax
      real*8 
     &  E1(imax,jmax,kmax),E2(imax,jmax,kmax),
     &  E3(imax,jmax,kmax),E4(imax,jmax,kmax),E5(imax,jmax,kmax),
     &  AA(kmax,kmax),LL(kmax,kmax),UU(kmax,kmax),
     &  E1_f(kmax),E2_f(kmax),E3_f(kmax),E4_f(kmax),E5_f(kmax),
     &  rhs1(kmax),rhs2(kmax),rhs3(kmax),rhs4(kmax),rhs5(kmax)

      call lhsf_P(kmax,AA)
      call ludecomp(1,kmax,AA,LL,UU)

      do j=1,jmax
        do i=1,imax
          call rhsf_z_P(imax,jmax,kmax,i,j,E1,rhs1)
          call solve(1,kmax,LL,UU,rhs1,E1_f)
          call rhsf_z_P(imax,jmax,kmax,i,j,E2,rhs2)
          call solve(1,kmax,LL,UU,rhs2,E2_f)
          call rhsf_z_P(imax,jmax,kmax,i,j,E3,rhs3)
          call solve(1,kmax,LL,UU,rhs3,E3_f)
          call rhsf_z_P(imax,jmax,kmax,i,j,E4,rhs4)
          call solve(1,kmax,LL,UU,rhs4,E4_f)
          call rhsf_z_P(imax,jmax,kmax,i,j,E5,rhs5)
          call solve(1,kmax,LL,UU,rhs5,E5_f)

          E1_f(kmax)=E1_f(1)
          E2_f(kmax)=E2_f(1)
          E3_f(kmax)=E3_f(1)
          E4_f(kmax)=E4_f(1)
          E5_f(kmax)=E5_f(1)

          do k=1,kmax
            E1(i,j,k)=E1_f(k)
            E2(i,j,k)=E2_f(k)
            E3(i,j,k)=E3_f(k)
            E4(i,j,k)=E4_f(k)
            E5(i,j,k)=E5_f(k)
          end do
        end do
      end do   

      RETURN
      END SUBROUTINE filter_z_P_4th

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
