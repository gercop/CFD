cccccc ********************************************************************
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos.              *            
cccccc Programmer\Guiding: Marcello Augusto Faraco de Medeiros.           *
cccccc Date      :                                                        *
cccccc Version	 : 1.0                                                    *
cccccc File      : converte.f90                                           *
cccccc ********************************************************************
cccccc Commands: ifort -fixed converte.f90 -o converte ; converte
cccccc ********************************************************************
      PROGRAM convert

      call main

      END

ccccc **********************************************************************
ccccc main: Main subroutine
ccccc **********************************************************************
      SUBROUTINE main

      IMPLICIT NONE    

      character*26 fname_aux
      character*50 linha
      integer tt,i,j,imax,jmax,tmax,qtimes,qtimes_aux
      parameter (imax=81,jmax=121,tmax=134,qtimes=200) 
      real*8 x(imax),y(jmax),var(imax,jmax)

      do tt=0,tmax
        qtimes_aux = int(tt*qtimes)

        write(*,'(A29,A1,I8,A1,f17.10,A1)')
     &    'CONVERTENDO ARQUIVOS (tt) = ','(',tt,')'

        call get_filename(qtimes_aux,'f_sn_wz','.dat',fname_aux) 

        open (1,file=fname_aux,status='old',access='sequential')
        read(1,'(a50)') linha
        read(1,'(a50)') linha
        do j=1,jmax
          do i=1,imax
            read(1,'(f20.15,f20.15,e30.16)') x(i),y(j),var(i,j)
          end do
        end do
        close(1)

        call get_filename(qtimes_aux,'f_sn_wz','.m',fname_aux) 
        open (1,file=fname_aux,status='unknown',access='sequential') 
        write(1,*) 'clear all' 
        write(1,*) 'close all;'
        do j=1,jmax
          do i=1,imax 
            write(1,'(a3,i6,a1,i6,a2,e30.16,a1)')
     &        'wz(',i,',',j,')=',var(i,j),';'
          end do
        end do
        write(1,*) 'figure;' 
        write(1,*) 'surface(wz(:,:));'
        write(1,*) 'view(90,270);'
        write(1,*) 'shading interp;'
        write(1,*) 'axis tight;'
        write(1,*) 'grid on;'   
        write(1,*) 'axis off;'   
        close(1)
      end do

      RETURN
      END

ccccc **********************************************************************
ccccc get_filename: Get filename.
ccccc **********************************************************************
      SUBROUTINE get_filename(tt,fname,ext,fname_aux)

      IMPLICIT NONE

      character*07 fname
      character*26 fname_aux
      character*04 ext
      character*01 c1
      character*02 c2
      character*03 c3
      character*04 c4
      character*05 c5
      character*06 c6
      character*07 c7
      character*08 c8
      integer tt

      if (tt.LT.10) then
        write (c1,'(I1)'),int(tt)
        fname_aux=fname//'_000000000'//c1//ext
      else if (tt.LT.100) then
        write (c2,'(I2)'),int(tt)
        fname_aux=fname//'_00000000'//c2//ext
       else if (tt.LT.1000) then
        write (c3,'(I3)'),int(tt)
        fname_aux=fname//'_0000000'//c3//ext
       else if (tt.LT.10000) then
        write (c4,'(I4)'),int(tt)
        fname_aux=fname//'_000000'//c4//ext
       else if (tt.LT.100000) then
        write (c5,'(I5)'),int(tt)
        fname_aux=fname//'_00000'//c5//ext
       else if (tt.LT.1000000) then
        write (c6,'(I6)'),int(tt)
        fname_aux=fname//'_0000'//c6//ext
       else if (tt.LT.10000000) then
        write (c7,'(I7)'),int(tt)
        fname_aux=fname//'_000'//c7//ext
       else if (tt.LT.100000000) then
        write (c8,'(I8)'),int(tt)
        fname_aux=fname//'_00'//c8//ext
      end if

      RETURN
      END


