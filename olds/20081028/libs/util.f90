MODULE util
  
  USE kinds, ONLY : dbl
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: quicksort, my_getarg, file_exists, wstrlen
  
CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE my_getarg(m,string)
    IMPLICIT NONE
    INTEGER:: M
    CHARACTER(len=*) :: STRING

    RETURN
  END SUBROUTINE my_getarg
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  FUNCTION file_exists(name)
    IMPLICIT NONE
    LOGICAL:: file_exists
    CHARACTER(len=*) :: name

    INQUIRE(file=name,exist=file_exists)

    RETURN
  END FUNCTION file_exists
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  FUNCTION wstrlen(s)
    IMPLICIT NONE
    INTEGER:: wstrlen
    CHARACTER(len=*) :: s
    INTEGER:: n

    n=LEN(s)
!    DO WHILE ( n>0 .AND. ( s(n:n)==' ' .OR. s(n:n)=='	' ) )
    DO
       IF ( n < 1 ) EXIT
       IF ( ( s(n:n) /= ' ' .AND. s(n:n) /= '	' ) ) EXIT
       n=n-1
    END DO
    wstrlen=n

    RETURN
  END FUNCTION wstrlen
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Sorting routine for the reciprocal space vectors (g)
! KB07AD HANDLES DOUBLE PRECISION VARIABLES
! STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
! THE WORK-SPACE 'MARK' OF LENGTH 50 PERMITS UP TO 2**(50/2)
! NUMBERS TO BE SORTED. THIS IS MORE THAN THE IBM VIRTUAL
! MEMORY SPACE WILL HOLD .
! ----------------------------------------------------------------------
  SUBROUTINE quicksort(COUNT,N,INDEX)
    IMPLICIT NONE
!     Arguments
    INTEGER N
    REAL ( dbl ) :: COUNT(N)
    INTEGER:: INDEX(N)
!     Variables
    INTEGER:: MARK(50),I,K,K1,LNGTH,IP,IFK,INTEST,INT,IY,J,IS1,IFKA, &
         MLOOP,IF,IS,LA,M
    REAL ( dbl ) :: AV,X
!--------------------------------------------------------------==
!  SET INDEX ARRAY TO ORIGINAL ORDER .
!--------------------------------------------------------------==
    DO I=1,N
      INDEX(I)=I
    END DO
!--------------------------------------------------------------==
! CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.
!--------------------------------------------------------------==
    IF(N.EQ.1)GOTO 200
    IF(N.GE.1)GOTO 30
    WRITE(*,20)
20  FORMAT(/,20X, &
         '***KB07AD***NO NUMBERS TO BE SORTED ** RETURN TO CALLING PROGRAM' &
         )
    GOTO 200
!--------------------------------------------------------------==
! 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER
! THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.
!--------------------------------------------------------------==
30  M=12
!--------------------------------------------------------------==
! SET UP INITIAL VALUES.
!--------------------------------------------------------------==
    LA=2
    IS=1
    IF=N
    DO MLOOP=1,N
!--------------------------------------------------------------==
!  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE.
!--------------------------------------------------------------==
      IFKA=IF-IS
      IF((IFKA+1).GT.M)GOTO 70
!--------------------------------------------------------------==
! FINAL SORTING  ( A SIMPLE BUBBLE SORT )
!--------------------------------------------------------------==
      IS1=IS+1
      DO J=IS1,IF
        I=J
40      IF(COUNT(I-1).LT.COUNT(I))GOTO 60
        IF(COUNT(I-1).GT.COUNT(I))GOTO 50
        IF(INDEX(I-1).LT.INDEX(I))GOTO 60
50      AV=COUNT(I-1)
        COUNT(I-1)=COUNT(I)
        COUNT(I)=AV
        INT=INDEX(I-1)
        INDEX(I-1)=INDEX(I)
        INDEX(I)=INT
        I=I-1
        IF(I.GT.IS)GOTO  40
60    END DO
      LA=LA-2
      GOTO 170
!--------------------------------------------------------------==
!                *******  QUICKSORT  ********
! SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS
! THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENTS
! HIGHEST ADDRESS.
!--------------------------------------------------------------==
70    IY=(IS+IF)/2
      X=COUNT(IY)
      INTEST=INDEX(IY)
      COUNT(IY)=COUNT(IF)
      INDEX(IY)=INDEX(IF)
!--------------------------------------------------------------==
! THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END
! OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE
! OF X .
!--------------------------------------------------------------==
      K=1
      IFK=IF
!--------------------------------------------------------------==
! WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE
! INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS
! NECESSARY, UNTIL THEY MEET .
!--------------------------------------------------------------==
      DO I=IS,IF
        IF(X.GT.COUNT(I))GOTO 110
        IF(X.LT.COUNT(I))GOTO 80
        IF(INTEST.GT.INDEX(I))GOTO 110
80      IF(I.GE.IFK)GOTO 120
        COUNT(IFK)=COUNT(I)
        INDEX(IFK)=INDEX(I)
        K1=K
        DO K=K1,IFKA
          IFK=IF-K
          IF(COUNT(IFK).GT.X)GOTO 100
          IF(COUNT(IFK).LT.X)GOTO 90
          IF(INTEST.LE.INDEX(IFK))GOTO 100
90        IF(I.GE.IFK)GOTO 130
          COUNT(I)=COUNT(IFK)
          INDEX(I)=INDEX(IFK)
          GO TO 110
100     END DO
        GOTO 120
110   END DO
!--------------------------------------------------------------==
! RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER
! WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO
! 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL
! TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
! INDEPENDENTLY .
!--------------------------------------------------------------==
120   COUNT(IFK)=X
      INDEX(IFK)=INTEST
      IP=IFK
      GOTO 140
130   COUNT(I)=X
      INDEX(I)=INTEST
      IP=I
!--------------------------------------------------------------==
!  STORE THE LONGER SUBDIVISION IN WORKSPACE.
!--------------------------------------------------------------==
140   IF((IP-IS).GT.(IF-IP))GOTO 150
      MARK(LA)=IF
      MARK(LA-1)=IP+1
      IF=IP-1
      GOTO 160
150   MARK(LA)=IP-1
      MARK(LA-1)=IS
      IS=IP+1
!--------------------------------------------------------------==
! FIND THE LENGTH OF THE SHORTER SUBDIVISION.
!--------------------------------------------------------------==
160   LNGTH=IF-IS
      IF(LNGTH.LE.0)GOTO 180
!--------------------------------------------------------------==
! IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
!--------------------------------------------------------------==
      LA=LA+2
      GOTO 190
170   IF(LA.LE.0)GOTO 200
!--------------------------------------------------------------==
! OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
!--------------------------------------------------------------==
180   IF=MARK(LA)
      IS=MARK(LA-1)
190 END DO
!--------------------------------------------------------------==
200 RETURN
  END SUBROUTINE quicksort




END MODULE util
