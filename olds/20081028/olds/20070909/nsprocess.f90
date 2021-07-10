      PROGRAM snafterpro

      USE sndat2plt, only: dat2plt
      USE sndat2mat, only: dat2mat
      USE sndat2fft, only: dat2fft

      IMPLICIT NONE
      INCLUDE 'nspar.f90'

      integer t_ini,t_end

      t_ini = 0
      t_end = 0

      write(*,*) 
      write(*,*) ' _______________________________________________'
      write(*,*) 
      write(*,*) ' PROGRAMA DE PÓS-PROCESSAMENTO                  '
      write(*,*) ' _______________________________________________'
      write(*,*) 
      write(*,'(a23,\)') 'TEMPO INICIAL.......:'
      read (*,*) t_ini
      write(*,*) 
      write(*,'(a23,\)') 'TEMPO FINAL.........:'
      read (*,'(i6,\)') t_end

      call dat2plt(t_ini,t_end)
cc      call dat2mat(t_ini,t_end)
cc      call dat2fft

      write(*,*) ' _______________________________________________'
      write(*,*) 
      write(*,*) ' FINALIZAÇÃO DO PROGRAMA DE PÓS-PROCESSAMENTO   '
      write(*,*) ' _______________________________________________'
      write(*,*) 

      END
