      subroutine getsimilarity2(f,fp,fpp,g,gp,vise,rgas,prandtl,
     &     rrhovise,te,pe,deta,ue2rhe,ny)
* integration subroutine for similarity solution getsimilarity.f
************************************************************************
*
      dimension f(0:ny),fp(0:ny),fpp(0:ny),g(0:ny),gp(0:ny)

      do n=1,ny
         t    = te*g(n-1)
         rho  = pe/(rgas*t)
         vis  = vise*sutherland(t/te,te,vise)
         c    = rho*vis*rrhovise
         rc   = 1./c
         if (n.eq.1) cold = c
         dcde = (c - cold)/deta
         cold = c
         rhs    = rc*(fpp(n-1)*dcde + f(n-1)*fpp(n-1))
         fpp(n) = fpp(n-1) - deta*rhs
         rhs    = rc*( gp(n-1)*dcde + prandtl*(f(n-1)*gp(n-1) + 
     &        c*ue2rhe*fpp(n-1)**2) )
         gp(n)  = gp(n-1)  - deta*rhs
         fp(n)  = fp(n-1) + deta*fpp(n)
         f (n)  = f (n-1) + deta*fp (n)
         g (n)  = g (n-1) + deta*gp (n)
      enddo

      end
