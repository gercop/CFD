      MODULE snbound
      
      contains
ccccc **********************************************************************
ccccc extrapolated_inlet: This routine extrapolated the variables at 
ccccc                     inflow boundaries.
ccccc **********************************************************************
      SUBROUTINE extrapolated_inlet(local_imax,j,k,fc)

      USE snmethod

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer j,k,local_imax 
      real*8 fc(local_imax,jmax,kmax)

      !Extrapolation of sixth order
      fc(1,j,k) = 1.d0/stx_ex_6th_D(1,1)*(
     &  -stx_ex_6th_D(1,2)*fc(2,j,k)+
     &  -stx_ex_6th_D(1,3)*fc(3,j,k)+
     &  -stx_ex_6th_D(1,4)*fc(4,j,k)+
     &  -stx_ex_6th_D(1,5)*fc(5,j,k)+
     &  -stx_ex_6th_D(1,6)*fc(6,j,k)+
     &  -stx_ex_6th_D(1,7)*fc(7,j,k))

      RETURN
      END SUBROUTINE extrapolated_inlet

ccccc **********************************************************************
ccccc extrapolated_outlet: This routine extrapolated the variables at
ccccc                      outflow boundaries.
ccccc **********************************************************************
      SUBROUTINE extrapolated_outlet(local_imax,j,k,fc)

      USE snmethod

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer j,k,local_imax
      real*8 fc(local_imax,jmax,kmax)

      !Extrapolation of fifth order
      fc(local_imax,j,k) = 1.d0/stxx_ex_6th_D(local_imax,7)*(
     &  -stxx_ex_6th_D(local_imax,1)*fc(local_imax-6,j,k)+
     &  -stxx_ex_6th_D(local_imax,2)*fc(local_imax-5,j,k)+
     &  -stxx_ex_6th_D(local_imax,3)*fc(local_imax-4,j,k)+
     &  -stxx_ex_6th_D(local_imax,4)*fc(local_imax-3,j,k)+
     &  -stxx_ex_6th_D(local_imax,5)*fc(local_imax-2,j,k)+
     &  -stxx_ex_6th_D(local_imax,6)*fc(local_imax-1,j,k))

      RETURN
      END SUBROUTINE extrapolated_outlet

ccccc **********************************************************************
ccccc boundflow: This routine calculated the boundary conditions for 
ccccc            spatial development of shear layer flow.
ccccc **********************************************************************
      SUBROUTINE boundflow(local_imax,x,y,z,t,rho,u,v,w,Et,p,tp,mu)

      USE snmpi,     only: my_rank,pro
      USE sninitial, only: rho_0,u_0,v_0,tp_0

      IMPLICIT NONE
      INCLUDE 'snparv35.f90'

      integer i,j,k,local_imax
      real*8 t,yy,zz,u_1,v_1,w_1,A_1,A_2,phi,
     &  x(local_imax),y(jmax),z(kmax),C2,
     &  rho(local_imax,jmax,kmax),u  (local_imax,jmax,kmax),
     &  v  (local_imax,jmax,kmax),w  (local_imax,jmax,kmax),
     &  Et (local_imax,jmax,kmax),p  (local_imax,jmax,kmax),
     &  tp (local_imax,jmax,kmax),mu (local_imax,jmax,kmax),
     &  e  (local_imax,jmax,kmax)
      parameter (C2=110.4d0/288.15d0) !T_oo=288.15 for sea level.
     
      !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5).OR.(tp_s.EQ.6)
        do k=1,kmax 
          zz=z(k) 
          do j=1,jmax
            yy=y(j)          
 
            if (my_rank.EQ.0) then
              !MS$IF (tp_s.EQ.4).OR.(tp_s.EQ.5)
                if ((jmax.NE.1).AND.(kmax.EQ.1)) then      !Two-Dimensional Problem
                  u_1 = 
     &              -2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*
     &                A*dsin(-omega*t)/alpha
     &              -2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*
     &                A/10.d0*dsin(-omega/2.d0*t)/alpha
     &              -2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*
     &                A/20.d0*dsin(-omega/4.d0*t)/alpha
                  v_1 = 
     &              +A*dcos(-omega*t)*dexp(-sigma*yy**2.d0)
     &              +A/10.d0*dcos(-omega/2.d0*t)*dexp(-sigma*yy**2.d0)
     &              +A/20.d0*dcos(-omega/4.d0*t)*dexp(-sigma*yy**2.d0)
                  w_1 = 0.d0
                else if ((jmax.NE.1).AND.(kmax.NE.1)) then !Three-Dimensional Problem
                  A_1 = A
                  A_2 = A/2.d0
                  phi = 1.d0/2.d0*pi                             
                  u_1 = 
     &              -dexp(-sigma*yy**2.d0)*(
     &                +A_2*beta*dcos(beta*zz-omega*t)/alpha
     &                -A_2*beta*dcos(beta*zz-omega*t)/alpha)
     &              +2.d0*sigma*yy*dexp(-sigma*yy**2.d0)*(
     &                +A_1*dsin(-omega*t+phi)/alpha+ 
     &                +A_2*dsin(+beta*zz-omega*t)/alpha
     &                +A_2*dsin(-beta*zz-omega*t)/alpha)
                  v_1 = 
     &              (A_1*dcos(-omega*t+phi)+
     &               A_2*dcos(+beta*zz-omega*t)+
     &               A_2*dcos(-beta*zz-omega*t))*
     &              dexp(-sigma*yy**2.d0)
                  w_1 = 
     &              (A_1*dcos(-omega*t+phi)+
     &               A_2*dcos(+beta*zz-omega*t)+
     &               A_2*dcos(-beta*zz-omega*t))*
     &              dexp(-sigma*yy**2.d0)
                end if
              !MS$ELSEIF (tp_s.EQ.6)
                u_1 = 0.d0
                v_1 = 0.d0
                w_1 = 0.d0
              !MS$ENDIF               
          
              u  (1,j,k)  = u_0  (1,j,k) + 1.d0/(M1*c)*u_1 
              v  (1,j,k)  = v_0  (1,j,k) + 1.d0/(M1*c)*v_1
              w  (1,j,k)  =              + 1.d0/(M1*c)*w_1
              rho(1,j,k)  = rho_0(1,j,k)
              tp (1,j,k)  = tp_0 (1,j,k)     

              call extrapolated_inlet(local_imax,j,k,p)        !Pressure is extrapolated to satisfy dp_dx=0

              !MS$IF (tp_proc.EQ.1)
cc                rho(1,j,k) = (gamma*Ma**2.d0)*p(1,j,k)/tp(1,j,k)
              !MS$ELSEIF (tp_proc.EQ.2)
cc                rho(1,j,k) = (p(1,j,k)*gamma*Ma**2.d0)**(1.d0/gamma)
              !MS$ENDIF         
              e  (1,j,k) = p(1,j,k)/(rho(1,j,k)*(gamma-1.d0)) 
              Et (1,j,k) = rho(1,j,k)*(e(1,j,k)+1.d0/2.d0*(
     &          u(1,j,k)**2.d0+v(1,j,k)**2.d0+w(1,j,k)**2.d0))
              mu (1,j,k) = tp(1,j,k)**(3.d0/2.d0)*((1.d0+C2)/(tp(1,j,k)+C2))
            end if

            if (my_rank+1.EQ.pro) then
              call extrapolated_outlet(local_imax,j,k,rho)     !Density is extrapolated to satisfy d2u_dx2=0
              call extrapolated_outlet(local_imax,j,k,u  )     !Velocity is extrapolated to satisfy d2u_dx2=0
              call extrapolated_outlet(local_imax,j,k,v  )     !Velocity is extrapolated to satisfy d2v_dx2=0
              call extrapolated_outlet(local_imax,j,k,w  )     !Velocity is extrapolated to satisfy d2w_dx2=0
              call extrapolated_outlet(local_imax,j,k,tp )     !Temperature is extrapolated to satisfy d2tp_dx2=0
           
              p(local_imax,j,k) = 1.d0/(gamma*Ma**2.d0)        !Pressure is fixed in x(imax)

              !MS$IF (tp_proc.EQ.1)
cc                rho(local_imax,j,k) = (gamma*Ma**2.d0)*p(local_imax,j,k)/tp(local_imax,j,k)
              !MS$ELSEIF (tp_proc.EQ.2)
cc                rho(local_imax,j,k) = (p(local_imax,j,k)*gamma*Ma**2.d0)**(1.d0/gamma)
              !MS$ENDIF         
              e  (local_imax,j,k) = p(local_imax,j,k)/(rho(local_imax,j,k)*(gamma-1.d0)) 
              Et (local_imax,j,k) = rho(local_imax,j,k)*(e(local_imax,j,k)+1.d0/2.d0*(
     &          u(local_imax,j,k)**2.d0+v(local_imax,j,k)**2.d0+w(local_imax,j,k)**2.d0)) 
              mu (local_imax,j,k) = tp(local_imax,j,k)**(3.d0/2.d0)*((1.d0+C2)/(tp(local_imax,j,k)+C2))
            end if
          end do
        end do 
      !MS$ENDIF       

      !Free-slip boundary condition - impermeability
      !MS$IF (tp_s.EQ.3).OR.(tp_s.EQ.4).OR.(tp_s.EQ.5)
        do k=1,kmax
          do i=1,local_imax
            v(i,1,k)    = 0.d0
            v(i,jmax,k) = 0.d0
          end do
        end do
      !MS$ENDIF       

      !Wall-slip boundary condition
      !MS$IF (tp_s.EQ.6)
        do k=1,kmax
          do i=1,local_imax
            u(i,1,k) = 0.d0
            v(i,1,k) = 0.d0
            w(i,1,k) = 0.d0

            !MS$IF (tp_proc.EQ.1)
cc              rho(i,jmax,k) = (gamma*Ma**2.d0)*p(i,jmax,k)/tp(i,jmax,k)
            !MS$ELSEIF (tp_proc.EQ.2)
cc              rho(i,jmax,k) = (p(i,jmax,k)*gamma*Ma**2.d0)**(1.d0/gamma)
            !MS$ENDIF         
cc            e (i,jmax,k) = p(i,jmax,k)/(rho(i,jmax,k)*(gamma-1.d0)) 
cc            Et(i,jmax,k) = rho(i,jmax,k)*(e(i,jmax,k)+1.d0/2.d0*(
cc     &        u(i,jmax,k)**2.d0+v(i,jmax,k)**2.d0+w(i,jmax,k)**2.d0))             
          end do
        end do
      !MS$ENDIF       

      RETURN
      END SUBROUTINE boundflow

      END MODULE snbound

