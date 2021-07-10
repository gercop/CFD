      subroutine getsimilarity(re_x,ue,te,pe,rgas,vise,gam,prandtl,
     &     simy,simu,simv,simr,simt,isimilarity,iprint)
* similarity solution of compressible laminar boundary layer
*
* everything is dimensional
*
************************************************************************
*
      dimension simy(0:isimilarity)
      dimension simu(0:isimilarity),simv(0:isimilarity)
      dimension simr(0:isimilarity),simt(0:isimilarity)

      parameter (ny=10000         ) ! integration steps
      parameter (irelax   = 100   ) ! max number iterations
      parameter (tol      = 1.e-10) ! iter. abortion criterion
      parameter (etamax   = 20.   ) ! compute up to this eta
      dimension f(0:ny),fp(0:ny),fpp(0:ny),g(0:ny),gp(0:ny)

      no=ny/isimilarity ! output every no profile point

      rhoe     = pe/(rgas*te)
      cp       = gam*rgas/(gam-1.)
      rhovise  = rhoe*vise
      rrhovise = 1./rhovise
      ue2rhe   = ue**2/(cp*te)
      reue     = rhoe*ue
c      he       = rgas/(gam-1.)*te+0.5*ue**2+pe/rhoe
      xi       = re_x*vise**2
      deta     = etamax/float(ny)
      sqr2xi   = sqrt(2.*xi)
      fac1     = sqr2xi/ue*deta
      fac2     = reue/vise
      fac3     = reue*vise/sqr2xi

* adiabatic wall flat plate BCs
      f (0)    = 0.
      fp(0)    = 0.
      gp(0)    = 0.
* initial guess (incompressible, laminar BL)
      fpp(0)   = 0.5
      g(0)     = 1+sqrt(prandtl)*ue**2/(gam*rgas*te)*0.5*(gam-1.)
* deltas
      dfpp     = fpp(0)*1.e-5
      dg       =   g(0)*1.e-5
      do i=1,irelax
         call getsimilarity2(f,fp,fpp,g,gp,vise,rgas,prandtl,
     &        rrhovise,te,pe,deta,ue2rhe,ny)
         fpe    = fp(ny)
         ge     =  g(ny)
         err=max(abs(fpe-1.),abs(ge-1.))
c         print '(a,i5,a,e12.5)','it ',i,' err= ',err
         if (err.lt.tol) goto 1
         fpp(0) = fpp(0)+dfpp ! delta in fpp direction
         call getsimilarity2(f,fp,fpp,g,gp,vise,rgas,prandtl,
     &        rrhovise,te,pe,deta,ue2rhe,ny)
         fpefp  = fp(ny)
         gefp   =  g(ny)
         fpp(0) = fpp(0)-dfpp
         g(0)   = g(0)+dg     ! delta in g direction
         call getsimilarity2(f,fp,fpp,g,gp,vise,rgas,prandtl,
     &        rrhovise,te,pe,deta,ue2rhe,ny)
         fpegp  = fp(ny)
         gegp   =  g(ny)
         g(0)   = g(0)-dg
* new guess
         dfpdg   = (fpegp-fpe)/dg
         dfpdfpp = (fpefp-fpe)/dfpp
         dgdg    = ( gegp- ge)/dg
         dgdfpp  = ( gefp- ge)/dfpp
         rjac    = 1./(dfpdg*dgdfpp-dfpdfpp*dgdg)
         ddg     = rjac*( dgdfpp*(1.-fpe)-dfpdfpp*(1.-ge))
         ddfpp   = rjac*(-dgdg  *(1.-fpe)+dfpdg  *(1.-ge))
         fpp(0)  = fpp(0)+ddfpp
         g(0)    = g(0)  +ddg
      enddo
      print '(a)','no convergence in simlarity'
      call stopmpi
 1    continue
c      print '(a)','analytical BL solution to similarity.dat'
c      open (99,file='similarity.dat',status='unknown',form='formatted')
c      print '(a)','output format: y,eta p,rho,T u v'
c      print '(a)','               1 2   3     5 6 7'
      eta=0.
      y=0.
      adelta1 = 0.
      atheta  = 0.
      adelta3 = 0.
c      addelta = 0.
      atw=g(0)*te
      arw=pe/(rgas*atw)
c      hw=rgas/(gam-1.)*atw+pe/arw
      do n=0,ny
         uu=fp(n)*ue
         tt=g(n)*te
         rho=pe/(rgas*tt)
         vv=-fac3/rho*(f(n)-eta*fp(n))
c         if (mod(n,no).eq.0) write(99,'(7(1x,e12.5))') 
c     &        y,eta,pe,rho,tt,uu,vv
         if (mod(n,no).eq.0) then
            isim=n/no
            simy(isim)=y
            simu(isim)=uu
            simv(isim)=vv
            simr(isim)=rho
            simt(isim)=tt
         endif
         eta=eta+deta
         dy=fac1/rho
         y=y+dy
c         hh=rgas/(gam-1.)*tt+0.5*uu**2+pe/rho
         rruu=rho*uu
         adelta1 = adelta1 +           (1.-rruu/reue    )*dy
         atheta  = atheta  + rruu/reue*(1.-  uu/  ue    )*dy
         adelta3 = adelta3 + rruu/reue*(1.- (uu/  ue)**2)*dy
c         addelta = addelta + rruu/reue*((hh-he)/(hw-he) )*dy
      enddo
c      close(99)
      readelta1 = fac2*adelta1
      reatheta  = fac2*atheta
      readelta3 = fac2*adelta3
c      readdelta = fac2*addelta
      avisw     = vise*sutherland(atw/te,te,vise)
      acf=sqrt(2.)*arw*avisw/(rhovise*sqrt(re_x))*fpp(0)
      if (iprint.eq.1) then
         print '(a)','similarity solution:'
         print '(a,e12.5)','T_w [K]   = ',atw
         print '(a,e12.5)','c_F       = ',acf
         print '(a,e12.5)','Re_delta1 = ',readelta1
         print '(a,e12.5)','Re_theta  = ',reatheta
         print '(a,e12.5)','Re_delta3 = ',readelta3
c         print '(a,e12.5)','Re_Delta  = ',readdelta
         print '(a,e12.5)','H         = ',readelta1/reatheta
         print '(a,e12.5)','H_32      = ',readelta3/reatheta
      endif

      end
