c**********************************************************************
c     Programa de solução da equação de Poisson discretização de
c     diferenças finitas de 2a. ordem de solução por lsor
      program lsor_paralel

c     declaração de variaveis
c     implicit none
      include 'mpif.h'
      include 'index.for'
      integer my_rank                      !number of process
      integer p                            !number of the process
      integer status(MPI_STATUS_SIZE)      !return status for receiver
      integer mp1, mp2, m
      integer numproc

c     Start up MPI
      call MPI_Init(ierr)

c     Find out process rank
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

c     Find out number of processes
      call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

c     We have p-1 nodes in the environment,
c     because the master is not processing!
c     So the number of nodes is numproc -->

      numproc = p - 1

c     calculates the value of points in the x direction of each node
      m = (m_orig + 2*(numproc - 1)) / numproc             !

C     Stop the program if the number of total points is not exactly
c     divided by nodes and if the  number of points in the x
c     direction can not be used in multigrid program and if the
c     number of points in the y direction can not be used in
c     multigrid program

c      if (((m_orig + dble(2 * numproc) - 2.d0) / dble(numproc)+
c     &   1.d0.ne.m).or.((dble(m - 1) / 8.d0).ne.((m - 1) / 8))) then
c       mp1 = 8 * ((m - 1) / 8) + 1;
c       mp2 = 8 * (((m - 1) / 8) + 1) + 1;
c       if (my_rank.eq.0) then
c         write(*,*)' Number of points in the x direction can
c     &               not be used in multigrid solver'
c         write(*,*)' The number of points in the x direction
c     &                should be',mp1*numproc-2*numproc+2,' or ',
c     &                mp2*numproc-2*numproc+2
c       end if
c       call MPI_Finalize(ierr)
c       end
c      end if
c      if ( (dble(n - 1) / 8.d0).ne.((n - 1) / 8) ) then
c        if (my_rank.eq.0) then
c          write(*,*)' Number of points in the y direction can
c     &                not be used in multigrid solver'
c          write(*,*)' The number of points in the y direction
c     &                should be ', 8 * ((n - 1) / 8) + 1,' or ',
c     &                8 * (((n - 1) / 8) + 1) + 1)
c        end if
c        call MPI_Finalize(ierr)
c        end
c      end if

      if (my_rank.eq.0) then !Master functions
        call master(numproc,my_rank,m)
      else !Nodes functions
        call nodes(numproc,my_rank,m)
      endif

c     shut down MPI
      call MPI_Finalize(ierr)

      end !main

c**********************************************************************
c     solves the tridiagonal matrix by LU decomposetion
      subroutine tridv(r,m,a,b,c)

c     Declaration of variables
      implicit none
      integer i,m
      real*8 a(m),b(m),c(m),r(m),u(m)
      real*8 bet,gam(m)

      bet=b(1)
      u(1)=r(1)/bet
      do i=2,m
        gam(i) = c(i-1)/bet
        bet = b(i) - a(i)*gam(i)
        u(i) = (r(i) - a(i)*u(i-1))/bet
      end do
      do i=m-1,1,-1
        u(i)=u(i) - gam(i+1)*u(i+1)
      end do
      do i=1,m
        r(i)=u(i)
      end do

      return
      end !tridv

c**********************************************************************
c     find the error in the poisson equation
      real*8 function erro(fc,sc,m)

      implicit none
      include 'index.for'
      integer i,j,m
      real*8 dxx, dyy, fc(m,n),sc(m,n)

      dxx = 1.d0 / dble(m_orig - 1)
      dyy = 1.d0 / dble(n - 1)
      dxx = dxx * dxx
      dyy = dyy * dyy

      erro=0.d0
      do j = 2, n-1
        do i = 2, m-1
          erro = max(erro, dabs(sc(i,j) -
     &          (fc(i + 1,j) - 2.d0*fc(i,j) + fc(i - 1,j))/dxx -
     &          (fc(i,j + 1) - 2.d0*fc(i,j) + fc(i,j - 1))/dyy  ))
        end do
      end do

      return
      end !erro

c**********************************************************************
c     writes the final matrix
      subroutine escreve(fc)

      implicit none
      include 'index.for'
      integer i,j
      real*8 fc(m_orig,n)

      do j = 1, n
        do i = 1, m_orig
          write(*,10) i,j,fc(i,j)
        end do
      end do

   10 format(2(1X,I5),4X,1D13.8)

      return
      end !escreve

c**********************************************************************
c     Subroutine that gives the inital values of sc in all
c     domain and fc at the boundaries
      subroutine initiate(fc,sc,rank,numproc,m)

      implicit none
      include 'index.for'
      integer i, j, m, rank, numproc, offset
      real*8 x, y, dx, dy, fc(m,n), sc(m,n)

      dx = 1.d0/dble(m_orig-1)
      dy = 1.d0/dble(n-1)

      offset = (rank - 1)*(m - 2)

c     gives the values of sc in all the domain and null
c     value for fc
      do i = 1, m
        x = dble(i-1 + offset) * dx
        do j = 1, n
          y = dble(j-1) * dy
          sc(i,j) = -2.d0*dcos(x)*dsin(y)
          fc(i,j) = 0.d0
        end do
      end do

c     gives the value of fc at bottom and top boundaries
      do i = 1, m
        fc(i,1) = dcos(dble(i-1 + offset)*dx)*dsin(0.d0)
        fc(i,n) = dcos(dble(i-1 + offset)*dx)*dsin(dble(n-1)*dy)
      end do

c     gives the value of fc at left boundary if rank = 1 and
c     right boundary
c     if rank = last
c     because 0 is the master and it is not processing...

      do j = 1,n
        if (rank .eq. 1) then
          fc(1,j) = dcos(0.d0)*dsin(dble(j-1)*dy)
        endif
        if (rank .eq. numproc) then
          fc(m,j) = dcos(dble(m_orig-1)*dx)*dsin(dble(j-1)*dy)
        endif
      end do

      return
      end !initiate

c**********************************************************************
c     solves poisson equation using Gauss-Seidel method
      subroutine gauss_seidel(fc,sc,m)

      implicit none
      include 'index.for'
      integer i, j,m
      real*8 dx, dy, dyy, aa,fc(m,n),sc(m,n)

      dx = 1.d0 / dble(m_orig-1)
      dy = 1.D0 / dble(n-1)
      dyy = dy * dy
      aa = dyy / (dx*dx)

      do j = 2, n-1
        do i = 2, m-1
          fc(i,j) = ( - sc(i,j) * dyy +
     &             aa*(fc(i + 1,j) + fc(i - 1,j))+
     &             fc(i,j + 1) + fc(i,j - 1) )/(2.d0+2.d0*aa)
        end do
      end do

      return
      end !gauss_seidel

c**********************************************************************
c     solves poisson equation using line sor solver
      subroutine lsor_only(fc,sc,m)

      implicit none
      include 'index.for'
      integer i, j, m
      real*8 aa, rf, dx, dy, dyy, fc(m,n), sc(m,n)
      real*8 a(n), b(n), c(n), rhs(n)

      dx=1.d0/dble(m_orig-1)
      dy=1.d0/dble(n-1)
      aa = dy*dy/dx/dx
      dyy = dy*dy
      rf = 1.9d0

c     LHS of the tridiagonal matrix
      a(1)=0.d0
      b(1)=1.d0
      c(1)=0.d0
      do j=2,n-1
        a(j)=1.d0
        b(j)=-2.d0*(1.d0+aa)
        c(j)=1.d0
      end do
      a(n)=0.d0
      b(n)=1.d0
      c(n)=0.d0

      do i=2,m-1
        rhs(1) = fc(i,1)
        do j=2,n-1
          rhs(j) = sc(i,j)*dyy - aa*( fc(i+1,j) + fc(i-1,j) )
        end do
        rhs(n) = fc(i,n)

        call tridv(rhs,n,a,b,c)

        do j=2,n-1
          fc(i,j) = fc(i,j) + rf*( rhs(j) - fc(i,j) )
        end do

      end do

      return
      end !lsor_only

C***********************************************************************
C     Sub program of the function of the nodes
      subroutine nodes(num_proc,my_rank,m)

      implicit none
      include 'mpif.h'
      include 'index.for'
      integer tag, m, my_rank,num_proc,i,j
      real*8 fc(m,n), sc(m,n), dx, dy
      real*8 erro, recv_col(io¾ËÒ•OQküÕ½ŒY=iaÒÆ|x]òÚëğõJä¦Dú6Ó&&şvÒ8‰ }x±ôë®\âÃÜ²ğıo^vşİ’p#}jJ11'ƒŸæë­/_i\‚yº–*ÆWî€+µÈàÓL^cyôó>(á(GZ_Y!/½9;zJËe\/	®‚WåŸXƒƒ˜Ûbç½HùGÊOÆnl	nóZæ»hdìÌà¯áçÎ¡¹ğâãgÈÙ×+:äQ~5&£¦°’kvãpÈáY”G†
ÖjÅea¸AŒdj)µÙ‚×‚¬Ô´	Öš5Ù7p¾98å›d:3¥MßĞEFğºò:a*Is‰+Šú›4Æ‘óû™ÈUÛ-xğóôšæ‡
ÂkW>óÌÎÏ¤xş3ıN7Ş¾
Ş®m„+ùÁc}Jğ“5y@ÚÌ\0öª3q~[‰Bf5°µXÉ%åòOx¯ ]L¬‘¶F´Š…F‰ƒšD`çVXz7g‘yâ>Õœo3Åãóeb¡MÍ¼6Œ¢ÏßD- L•ı¨ã­5fñQ¬<$"=¯SÊÊåŠİøš>—}/Ã«Ÿ…
Ï§ådş
foûÛ8÷÷°"ğ€Iz¯øŸïFîŠˆß‰GF´¼¿NØüc„·â×ët÷íûÑ¸auss_seidel(fc,sc,m)
C        call mult_grid(fc,sc,m)

c       Send the colun of the end of domain
        if (my_rank.lt.num_proc ) then

C         take the colun to comunicate
          do j=1, n
            send_col(j) = fc(m-1,j)
          end do

C         send the colun
C         if tag = my_rank*10 + 2, the colun of end of domain is send or recei
C         if tag = my_rank*10 + 1, the colun of begin of domain is send or rece
          tag = my_rank*10 + 2
          call MPI_Send(send_col(1),n,MPI_DOUBLE_PRECISION, 0,
     &                tag,MPI_COMM_WORLD, ierr)
        endif
        
C       send the colun from begin of domain and receive of the end
        if (my_rank.gt.1) then
C         take the colun to comunicate
          do j=1, n
            send_col(j) = fc(2,j)
          end do

C         send the colun
          tag = my_rank*10 + 1
          call MPI_Send(send_col(1),n,MPI_DOUBLE_PRECISION, 0,
     &        tag,MPI_COMM_WORLD, ierr)

        endif

C       Receive the colun from the master
        if (my_rank.lt.num_proc) then !colun from the beging of domain of next node
C         receive the colunm
          tag = (my_rank + 1)*10 + 1
          call MPI_Recv(recv_col,n,MPI_DOUBLE_PRECISION, 0,
     &                tag, MPI_COMM_WORLD, status, ierr)

C         incorporate the colunm on matrix
          do j = 1, n
            fc(m,j) = recv_col(j)
          end do
        endif
        if (my_rank.gt.1) then
C         receive the colunm
          tag = (my_rank - 1)*10 + 2
          call MPI_Recv(recv_col,n,MPI_DOUBLE_PRECISION, 0,
     &        tag, MPI_COMM_WORLD, status, ierr)

C         incorporate the colunm on matrix
          do j = 1, n
            fc(1,j) = recv_col(j)
          end do
        endif

c       calculates the error in the Poisson equation
        eerr = erro(fc,sc,m)

c       Send the error to the Master
        tag = ERR_MSG + my_rank
        call MPI_Send(eerr, 1, MPI_DOUBLE_PRECISION, 0,
     &      tag, MPI_COMM_WORLD, ierr)
        
c       Receive the message to end or continue the program
c       from the MASTER
        tag = ERR_MSG + my_rank
        call MPI_Recv(sum, 1, MPI_DOUBLE_PRECISION, 0,
     &      tag, MPI_COMM_WORLD, status, ierr)

      end do

c       Send the matrices to the Master to join in the final matrix
        call MPI_Send(fc(1,1), n * m, MPI_DOUBLE_PRECISION, 0,
     &      FINAL+my_rank, MPI_COMM_WORLD, ierr)

      return
      end !nodes

C***********************************************************************
C     Subroutine of master node
      subroutine master(numproc,my_rank,m)

      implicit none
      include 'mpif.h'
      include 'index.for'
      integer numproc, my_rank, m, pos
      real*8 buffercol(2*numproc-1,n)
      real*8 fc_final(m_orig,n), sc_final(m_orig,n)
      real*8 recv_col(n), send_col(n)
      real*8 sum, eerr
      integer i, j, it, ierr, tag
      integer ERR_MSG, INFO_MSG,FINAL
      parameter (ERR_MSG = 1000, INFO_MSG = 399, FINAL = 53)
      integer status(MPI_STATUS_SIZE)      !return status for receiver
      
      eerr = 1.d0
      sum = eerr
      do while (sum.gt.1.d-8)
c       Parallel processing of the Poisson equation
C       receive, from the nodes, the colun from the end of domain
        do i = 1, numproc - 1
          tag = i*10 + 2
          call MPI_Recv(recv_col, n, MPI_DOUBLE_PRECISION, i,tag,
     &        MPI_COMM_WORLD, status,ierr)
          do j = 1,n
            buffercol(i*2,j) = recv_col(j)
          end do
        end do

C       Receive, from the nodes, the coluns from the beging of domain
        do i = 2, numproc
          tag = i*10 + 1
          call MPI_Recv(recv_col, n, MPI_DOUBLE_PRECISION, i,tag,
     &        MPI_COMM_WORLD, status,ierr)
          do j = 1,n
            buffercol(i*2-1,j) = recv_col(j)
          end do
        end do

c       Send, to the nodes, the coluns from the beging of domain
        do i = 1, numproc-1
          do j = 1, n
            send_col(j) = buffercol(2*(i+1) - 1,j)
          end do
          tag = (i+1)*10 + 1
          call MPI_Send(send_col(1), n, MPI_DOUBLE_PRECISION, i, tag,
     &              MPI_COMM_WORLD, ierr)
        end do

c       Send, to the nodes, the coluns from the end of domains
        do i = 2, numproc
          do j = 1, n
            send_col(j) = buffercol(2*(i-1),j)
          end do
          tag = (i-1)*10 + 2
          call MPI_Send(send_col(1), n, MPI_DOUBLE_PRECISION, i, tag,
     &              MPI_COMM_WORLD, ierr)
        end do

c       Receive the errors from the nodes
        sum = 0.d0
        do i = 1, nu©Ì"Î7Šã1)9BXÄü¤2yÑQ jÅ.9İŠ4Š«…ô5Ïäb…Ó--Î*–ß›¢y¤tr¾(×‹.4»ä¦¼ğC>U8=Ä?÷³ğm'ôÜ6Ò`‡¤(ü$e¾£2~Ëô¶*=„Z¶|/UšPîÜÆB±u®2cª«zól“
ühIXÊkúœéXõ}>D®CYªÎ?P}Ùˆõ¿èSÉ°2xş”É)˜ı÷­¦3ìhú–Ûò-~ºê£×<ÏIùYèßÄv>q[€íü?|İCİÏáñ’
£U³G U6@ıŒÏÑé¨ VÁb$TÒï²Ê° Ê,w­

 ò]REV‰VA#Şß£ sB…iZ€ş£oÆÁ7°´nƒŸÁW»å—a‘TvåÒÒ
¾ëX”¼ƒê˜R[ª>J6{/…7P‚RstğÓ£ìKp–Äa5c›lô'ÍÒ·ÏËõï†Oü9cb®
y"‘¤s%Xœj~®[*FàF»N·Bì¬ôçÒ¢t[|½1Ë2×"æ&½Idj‘·üìáAVcèÛÿÚ‹çÖ—å²ç3loú@Œ«µÃıëá°#ûÌnÖÏ¶a9õçÑtÜ+ĞÏâ	~¹<ô!lÅĞöğ´¾6l¿Cûµ¼o´İ×‹“_PRECISION, 1,
     &    FINAL+1, MPI_COMM_WORLD, status, ierr) !receive from node 1
      do i = 2, numproc
        pos = i * (m - 2)
        call MPI_Recv(fc_final(pos,1), n * m, MPI_DOUBLE_PRECISION,
     &      i, FINAL+i, MPI_COMM_WORLD, status, ierr)
      end do

c     write the matrix
c     call escreve(fc_final)

      return
      end !master




ccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                 c
c         Begin of the Multigrid solver           c
c                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fas(fc,sc,imax,jmax,dx,dy)

c Full Multigrid Approximation      
      implicit none
      integer i,j,imax,jmax,m,n,k
      real*8 fc(imax,jmax),sc(imax,jmax),res(imax,jmax)
      real*8 fc2(imax,jmax),sc2(imax,jmax),fct2(imax,jmax)
      real*8 fc3(imax,jmax),sc3(imax,jmax),fct3(imax,jmax)
      real*8 fc4(imax,jmax),sc4(imax,jmax),fct4(imax,jmax)
      real*8 maxr,lp,erro,dx,dy

      lp=0.
      m=imax
      n=jmax

c main loop of V multigrid

   21 call lsormg(fc,sc,imax,jmax,dx,dy,2) ! grid 1
      call defect(res,fc,sc,imax,jmax,dx,dy)
      call restrict(fc2,fct2,sc2,fc,res,m,n,dx,dy)
      call newrhs(sc2,fc2,m,n,dx,dy)

      call lsormg(fc2,sc2,m,n,dx,dy,2) ! grid 2
      call defect(res,fc2,sc2,m,n,dx,dy)
      call restrict(fc3,fct3,sc3,fc2,res,m,n,dx,dy)
      call newrhs(sc3,fc3,m,n,dx,dy)

      call lsormg(fc3,sc3,m,n,dx,dy,2) ! grid 3
      call defect(res,fc3,sc3,m,n,dx,dy)
      call restrict(fc4,fct4,sc4,fc3,res,m,n,dx,dy)
      call newrhs(sc4,fc4,m,n,dx,dy)

      call lsormg(fc4,sc4,m,n,dx,dy,40) ! grid 4

      call correction(fc4,fct4,m,n)
      call intpol(fc3,fc4,m,n,dx,dy)
      call lsormg(fc3,sc3,m,n,dx,dy,1) ! grid 3

      call correction(fc3,fct3,m,n)
      call intpol(fc2,fc3,m,n,dx,dy)
      call lsormg(fc2,sc2,m,n,dx,dy,1) ! grid 2

      call correction(fc2,fct2,m,n)
      call intpol(fc,fc2,m,n,dx,dy)
      call lsormg(fc,sc,imax,jmax,dx,dy,1) ! grid 1

c     call maxerr(maxr,fc,sc,m,n,erro,dx,dy)
c     lp=lp+1
c     write(*,*)real(lp),real(maxr)
c     if(abs(maxr).gt.erro) goto 21

c     write(*,*)'loops=',real(lp),'    coresponding iterations=',
c    &          real(lp*(2.+2./4.+2./16.+50./64.+1./16.+1./4.+1.))


      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine lsormg(fc,sc,m,n,dx,dy,itt)

c solves the poisson equation using a 2nd order LSOR
      implicit none
      integer i,j,m,n,it,z,itt
      real*8 fc(m,n),sc(m,n)
      real*8 a(m),b(m),c(m),rhs(m)
      real*8 dx,dy,dyy,maxr,aa,rf

      rf=1.5
      if (itt.eq.1) rf=1.

      dyy=dy*dy
      aa=dyy/dx/dx

c LHS of the tridiagonal matrix
      a(1)=0.d0
      b(1)=1.d0
      c(1)=0.d0
      do j=2,n-1
        a(j)=1.d0
        b(j)=-2.d0*(1.d0+aa)
        c(j)=1.d0
      end do
      a(n)=0.d0
      b(n)=1.d0
      c(n)=0.d0
            
c 2nd order LSOR
      do it=1,itt
        do i=2,m-1
          rhs(1)=fc(i,1)
          do j=2,n-1
            rhs(j)=sc(i,j)*dyy-aa*(fc(i+1,j)+fc(i-1,j))
          end do
          rhs(n)=fc(i,n)
          call tridv2(rhs,n,a,b,c) ! solve the matrix
          do j=2,n-1
            fc(i,j)=fc(i,j)+rf*(rhs(j)-fc(i,j))
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine defect(res,fc,sc,m,n,dx,dy)

c return the defect of the solved equation on variable res
      implicit none
      integer i,j,m,n
      real*8 fc(m,n),sc(m,n),res(m,n)
      real*8 dx,dy,dyy,dxx

      dxx=dx*dx
      dyy=dy*dy

c 2nd order defect
      do j=2,n-1
        do i=2,m-1
          res(i,j)=sc(i,j)-
     &            (fc(i+1,j)-2.d0*fc(i,j)+fc(i-1,j))/dxx-
     &            (fc(i,j+1)-2.d0*fc(i,j)+fc(i,j-1))/dyy
        end do
      end do

      do i=1,m
        res(i,1)=0.d0
        res(i,n)=0.d0
      end do
      do j=1,n
        res(1,j)=0.d0
        res(m,j)=0.d0
      end do

      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine restrict(fct,fcc,scc,fc,sc,m,n,dx,dy)

c calculate the distribuition for font and source term on coarser grid
c fine in -> coarse out (fc,sc-> in, fcc,scc-> out)
      implicit none
      integer m,n,ic,jc,jf
      real*8 fc(m,n),fcc(m,n),fct(m,n)
      real*8 sc(m,n),scc(m,n),dx,dy
      
      m=m/2+1
      n=n/2+1
      dx=2.d0*dx
      dy=2.d0*dy

c internal points
      do jc=2,n-1
        jf=2*jc-1
        do ic=2,m-1
          fcc(ic,jc)=fc(2*ic-1,jf)
          fct(ic,jc)=fcc(ic,jc)
          scc(ic,jc)=6.25d-2*( 4.d0*sc(2*ic-1,jf)+
     &               2.d0*(sc(2*ic-1,jf-1)+sc(2*ic-1,jf+1)+
     &               sc(2*ic,jf)+sc(2*ic-2,jf))+
     &               sc(2*ic-3,jf)+sc(2*ic+1,jf)+
     &               sc(2*ic-1,jf+2)+sc(2*ic-1,jf-2) )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine newrhs(sc,fc,m,n,dx,dy)

c calculate the rhs for the next coarse grid
      implicit none
      integer i,j,m,n
      real*8 fc(m,n),sc(m,n)
      real*8 dxx,dyy,dx,dy
      
      dxx=dx*dx
      dyy=dy*dy
     
c 2nd order SOR
      do j=2,n-1
        do i=2,m-1
          sc(i,j)=sc(i,j)+
     &            (fc(i-1,j)+fc(i+1,j)-2.d0*fc(i,j))/dxx+
     &            (fc(i,j+1)+fc(i,j-1)-2.d0*fc(i,j))/dyy
        end do
      end do

      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine correction(fc,fct,m,n)

c calculate the correction to use on the next fine grid
      implicit none
      integer i,j,m,n
      real*8 fc(m,n),fct(m,n)
      
      do j=2,n-1
        do i=2,m-1
          fc(i,j)=fc(i,j)-fct(i,j)
        end do
      end do

      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine intpol(fc,corr,m,n,dx,dy)

c coarse in -> fine out (fcc-> in, fc-> out)
      implicit none
      integer i,j,jc,m,n
      real*8 fc(m,n),corr(m,n),tmp(m,n),dx,dy

c intepolacao bilinear
      do jc=1,n
        j=2*jc-1
        do i=1,m
          tmp(2*i-1,j)=corr(i,jc)
        end do
      end do
      m=m*2-1
      n=n*2-1
      dx=.5d0*dx
      dy=.5d0*dy
      do j=1,n,2
        do i=2,m-1,2
          tmp(i,j)=.5d0*(tmp(i+1,j)+tmp(i-1,j))
        end do
      end do
      do j=2,n-1,2
        do i=1,m
          tmp(i,j)=.5d0*(tmp(i,j+1)+tmp(i,j-1))
        end do
      end do

c adition to new lhs      
      do j=2,n-1
        do i=2,m-1
          fc(i,j)=fc(i,j)+tmp(i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine maxerr(maxr,fc,sc,m,n,erro,dx,dy)

c return the maximum error of the solved equation on maxr
      implicit none
      integer i,j,m,n
      real*8 fc(m,n),sc(m,n)
      real*8 dx,dy,dyy,dxx,erro,maxr

      dxx=dx*dx
      dyy=dy*dy
      maxr=0.d0

c 2nd order appx.
      do j=2,n-1
        do i=2,m-1
          maxr=abs(sc(i,j)-
     &         ((fc(i+1,j)-2.d0*fc(i,j)+fc(i-1,j))/dxx)-
     &         ((fc(i,j+1)-2.d0*fc(i,j)+fc(i,j-1))/dyy) )
          if(maxr.gt.erro) return
        end do
      end do

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridv2(r,lmax,a,b,c) ! retirado do Numerical Recipes
      
c resolve um sistema tridiagonal
      implicit none
      integer i,lmax
      real*8 a(lmax),b(lmax),c(lmax),r(lmax),u(lmax)
      real*8 bet,gam(lmax)
      
      bet=b(1)
      u(1)=r(1)/bet
      do i=2,lmax
        gam(i)=c(i-1)/bet
        bet=b(i)-a(i)*gam(i)
        u(i)=(r(i)-a(i)*u(i-1))/bet
      end do
      do i=lmax-1,1,-1
        u(i)=u(i)-gam(i+1)*u(i+1)
      end do
      do i=1,lmax
        r(i)=u(i)
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                 c
c         End of the Multigrid solver             c
c                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
