C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE linsolve( A, B, M1, M2)
C ----------------------------------------------------------------------
C	
C	SOLVES 
C     
C     A*X=B 
C
C  FOR GENERAL MATRIX A 
C
C  NOTE:  A IS OVERWRITTEN, AND SOLUTION X IS 
C         RETURNED IN B.
C
C ----------------------------------------------------------------------

 
      IMPLICIT NONE
      INTEGER M1,M2
      REAL*8 A(M1,M2), B(M1)
      INTEGER N, INFO, IPIV(M1)

      N=M1

C LAPACK ROUTINES
	CALL DGETRF( N, N, A, N, IPIV, INFO)
	CALL DGETRS( 'NO TRANSPOSE', N, 1, A, N,
     &   		IPIV, B, N, INFO )
c      RETURN
      END 
      


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
C ----------------------------------------------------------------------
*
*  -- LAPACK ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETF2 COMPUTES AN LU FACTORIZATION OF A GENERAL M-BY-N MATRIX A
*  USING PARTIAL PIVOTING WITH ROW INTERCHANGES.
*
*  THE FACTORIZATION HAS THE FORM
*     A = P * L * U
*  WHERE P IS A PERMUTATION MATRIX, L IS LOWER TRIANGULAR WITH UNIT
*  DIAGONAL ELEMENTS (LOWER TRAPEZOIDAL IF M > N), AND U IS UPPER
*  TRIANGULAR (UPPER TRAPEZOIDAL IF M < N).
*
*  THIS IS THE RIGHT-LOOKING LEVEL 2 BLAS VERSION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE M BY N MATRIX TO BE FACTORED.
*          ON EXIT, THE FACTORS L AND U FROM THE FACTORIZATION
*          A = P*L*U; THE UNIT DIAGONAL ELEMENTS OF L ARE NOT STORED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (MIN(M,N))
*          THE PIVOT INDICES; FOR 1 <= I <= MIN(M,N), ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0: IF INFO = K, U(K,K) IS EXACTLY ZERO. THE FACTORIZATION
*               HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY
*               SINGULAR, AND DIVISION BY ZERO WILL OCCUR IF IT IS USED
*               TO SOLVE A SYSTEM OF EQUATIONS.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            J, JP
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        FIND PIVOT AND TEST FOR SINGULARITY.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           APPLY THE INTERCHANGE TO COLUMNS 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           COMPUTE ELEMENTS J+1:M OF J-TH COLUMN.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           UPDATE TRAILING SUBMATRIX.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     END OF DGETF2
*
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
C ----------------------------------------------------------------------
*
*  -- LAPACK ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETRF COMPUTES AN LU FACTORIZATION OF A GENERAL M-BY-N MATRIX A
*  USING PARTIAL PIVOTING WITH ROW INTERCHANGES.
*
*  THE FACTORIZATION HAS THE FORM
*     A = P * L * U
*  WHERE P IS A PERMUTATION MATRIX, L IS LOWER TRIANGULAR WITH UNIT
*  DIAGONAL ELEMENTS (LOWER TRAPEZOIDAL IF M > N), AND U IS UPPER
*  TRIANGULAR (UPPER TRAPEZOIDAL IF M < N).
*
*  THIS IS THE RIGHT-LOOKING LEVEL 3 BLAS VERSION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE M-BY-N MATRIX TO BE FACTORED.
*          ON EXIT, THE FACTORS L AND U FROM THE FACTORIZATION
*          A = P*L*U; THE UNIT DIAGONAL ELEMENTS OF L ARE NOT STORED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (MIN(M,N))
*          THE PIVOT INDICES; FOR 1 <= I <= MIN(M,N), ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  IF INFO = I, U(I,I) IS EXACTLY ZERO. THE FACTORIZATION
*                HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY
*                SINGULAR, AND DIVISION BY ZERO WILL OCCUR IF IT IS USED
*                TO SOLVE A SYSTEM OF EQUATIONS.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     DETERMINE THE BLOCK SIZE FOR THIS ENVIRONMENT.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        USE UNBLOCKED CODE.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        USE BLOCKED CODE.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           FACTOR DIAGONAL AND SUBDIAGONAL BLOCKS AND TEST FOR EXACT
*           SINGULARITY.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           ADJUST INFO AND THE PIVOT INDICES.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           APPLY INTERCHANGES TO COLUMNS 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              APPLY INTERCHANGES TO COLUMNS J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              COMPUTE BLOCK ROW OF U.
*
               CALL DTRSM( 'LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 UPDATE TRAILING SUBMATRIX.
*
                  CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     END OF DGETRF
*
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
C ----------------------------------------------------------------------
*
*  -- LAPACK ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETRS SOLVES A SYSTEM OF LINEAR EQUATIONS
*     A * X = B  OR  A' * X = B
*  WITH A GENERAL N-BY-N MATRIX A USING THE LU FACTORIZATION COMPUTED
*  BY DGETRF.
*
*  ARGUMENTS
*  =========
*
*  TRANS   (INPUT) CHARACTER*1
*          SPECIFIES THE FORM OF THE SYSTEM OF EQUATIONS:
*          = 'N':  A * X = B  (NO TRANSPOSE)
*          = 'T':  A'* X = B  (TRANSPOSE)
*          = 'C':  A'* X = B  (CONJUGATE TRANSPOSE = TRANSPOSE)
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  NRHS    (INPUT) INTEGER
*          THE NUMBER OF RIGHT HAND SIDES, I.E., THE NUMBER OF COLUMNS
*          OF THE MATRIX B.  NRHS >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          THE FACTORS L AND U FROM THE FACTORIZATION A = P*L*U
*          AS COMPUTED BY DGETRF.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  IPIV    (INPUT) INTEGER ARRAY, DIMENSION (N)
*          THE PIVOT INDICES FROM DGETRF; FOR 1<=I<=N, ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  B       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDB,NRHS)
*          ON ENTRY, THE RIGHT HAND SIDE MATRIX B.
*          ON EXIT, THE SOLUTION MATRIX X.
*
*  LDB     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY B.  LDB >= MAX(1,N).
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            NOTRAN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        SOLVE A * X = B.
*
*        APPLY ROW INTERCHANGES TO THE RIGHT HAND SIDES.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        SOLVE L*X = B, OVERWRITING B WITH X.
*
         CALL DTRSM( 'LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        SOLVE U*X = B, OVERWRITING B WITH X.
*
         CALL DTRSM( 'LEFT', 'UPPER', 'NO TRANSPOSE', 'NON-UNIT', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        SOLVE A' * X = B.
*
*        SOLVE U'*X = B, OVERWRITING B WITH X.
*
         CALL DTRSM( 'LEFT', 'UPPER', 'TRANSPOSE', 'NON-UNIT', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        SOLVE L'*X = B, OVERWRITING B WITH X.
*
         CALL DTRSM( 'LEFT', 'LOWER', 'TRANSPOSE', 'UNIT', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        APPLY ROW INTERCHANGES TO THE SOLUTION VECTORS.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     END OF DGETRS
*
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
C ----------------------------------------------------------------------
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1999
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLASWP PERFORMS A SERIES OF ROW INTERCHANGES ON THE MATRIX A.
*  ONE ROW INTERCHANGE IS INITIATED FOR EACH OF ROWS K1 THROUGH K2 OF A.
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE MATRIX OF COLUMN DIMENSION N TO WHICH THE ROW
*          INTERCHANGES WILL BE APPLIED.
*          ON EXIT, THE PERMUTED MATRIX.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.
*
*  K1      (INPUT) INTEGER
*          THE FIRST ELEMENT OF IPIV FOR WHICH A ROW INTERCHANGE WILL
*          BE DONE.
*
*  K2      (INPUT) INTEGER
*          THE LAST ELEMENT OF IPIV FOR WHICH A ROW INTERCHANGE WILL
*          BE DONE.
*
*  IPIV    (INPUT) INTEGER ARRAY, DIMENSION (M*ABS(INCX))
*          THE VECTOR OF PIVOT INDICES.  ONLY THE ELEMENTS IN POSITIONS
*          K1 THROUGH K2 OF IPIV ARE ACCESSED.
*          IPIV(K) = L IMPLIES ROWS K AND L ARE TO BE INTERCHANGED.
*
*  INCX    (INPUT) INTEGER
*          THE INCREMENT BETWEEN SUCCESSIVE VALUES OF IPIV.  IF IPIV
*          IS NEGATIVE, THE PIVOTS ARE APPLIED IN REVERSE ORDER.
*
*  FURTHER DETAILS
*  ===============
*
*  MODIFIED BY
*   R. C. WHALEY, COMPUTER SCIENCE DEPT., UNIV. OF TENN., KNOXVILLE, USA
*
* =====================================================================
*
*     .. LOCAL SCALARS ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     INTERCHANGE ROW I WITH ROW IPIV(I) FOR EACH OF ROWS K1 THROUGH K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     END OF DLASWP
*
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      INTEGER FUNCTION IEEECK( ISPEC, ZERO, ONE )
C ----------------------------------------------------------------------
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1998
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  PURPOSE
*  =======
*
*  IEEECK IS CALLED FROM THE ILAENV TO VERIFY THAT INFINITY AND
*  POSSIBLY NAN ARITHMETIC IS SAFE (I.E. WILL NOT TRAP).
*
*  ARGUMENTS
*  =========
*
*  ISPEC   (INPUT) INTEGER
*          SPECIFIES WHETHER TO TEST JUST FOR INIFINITY ARITHMETIC
*          OR WHETHER TO TEST FOR INFINITY AND NAN ARITHMETIC.
*          = 0: VERIFY INFINITY ARITHMETIC ONLY.
*          = 1: VERIFY INFINITY AND NAN ARITHMETIC.
*
*  ZERO    (INPUT) REAL
*          MUST CONTAIN THE VALUE 0.0
*          THIS IS PASSED TO PREVENT THE COMPILER FROM OPTIMIZING
*          AWAY THIS CODE.
*
*  ONE     (INPUT) REAL
*          MUST CONTAIN THE VALUE 1.0
*          THIS IS PASSED TO PREVENT THE COMPILER FROM OPTIMIZING
*          AWAY THIS CODE.
*
*  RETURN VALUE:  INTEGER
*          = 0:  ARITHMETIC FAILED TO PRODUCE THE CORRECT ANSWERS
*          = 1:  ARITHMETIC PRODUCED THE CORRECT ANSWERS
*
*     .. LOCAL SCALARS ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. EXECUTABLE STATEMENTS ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     RETURN IF WE WERE ONLY ASKED TO CHECK INFINITY ARITHMETIC
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
C ----------------------------------------------------------------------
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1999
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  PURPOSE
*  =======
*
*  ILAENV IS CALLED FROM THE LAPACK ROUTINES TO CHOOSE PROBLEM-DEPENDENT
*  PARAMETERS FOR THE LOCAL ENVIRONMENT.  SEE ISPEC FOR A DESCRIPTION OF
*  THE PARAMETERS.
*
*  THIS VERSION PROVIDES A SET OF PARAMETERS WHICH SHOULD GIVE GOOD,
*  BUT NOT OPTIMAL, PERFORMANCE ON MANY OF THE CURRENTLY AVAILABLE
*  COMPUTERS.  USERS ARE ENCOURAGED TO MODIFY THIS SUBROUTINE TO SET
*  THE TUNING PARAMETERS FOR THEIR PARTICULAR MACHINE USING THE OPTION
*  AND PROBLEM SIZE INFORMATION IN THE ARGUMENTS.
*
*  THIS ROUTINE WILL NOT FUNCTION CORRECTLY IF IT IS CONVERTED TO ALL
*  LOWER CASE.  CONVERTING IT TO ALL UPPER CASE IS ALLOWED.
*
*  ARGUMENTS
*  =========
*
*  ISPEC   (INPUT) INTEGER
*          SPECIFIES THE PARAMETER TO BE RETURNED AS THE VALUE OF
*          ILAENV.
*          = 1: THE OPTIMAL BLOCKSIZE; IF THIS VALUE IS 1, AN UNBLOCKED
*               ALGORITHM WILL GIVE THE BEST PERFORMANCE.
*          = 2: THE MINIMUM BLOCK SIZE FOR WHICH THE BLOCK ROUTINE
*               SHOULD BE USED; IF THE USABLE BLOCK SIZE IS LESS THAN
*               THIS VALUE, AN UNBLOCKED ROUTINE SHOULD BE USED.
*          = 3: THE CROSSOVER POINT (IN A BLOCK ROUTINE, FOR N LESS
*               THAN THIS VALUE, AN UNBLOCKED ROUTINE SHOULD BE USED)
*          = 4: THE NUMBER OF SHIFTS, USED IN THE NONSYMMETRIC
*               EIGENVALUE ROUTINES
*          = 5: THE MINIMUM COLUMN DIMENSION FOR BLOCKING TO BE USED;
*               RECTANGULAR BLOCKS MUST HAVE DIMENSION AT LEAST K BY M,
*               WHERE K IS GIVEN BY ILAENV(2,...) AND M BY ILAENV(5,...)
*          = 6: THE CROSSOVER POINT FOR THE SVD (WHEN REDUCING AN M BY N
*               MATRIX TO BIDIAGONAL FORM, IF MAX(M,N)/MIN(M,N) EXCEEDS
*               THIS VALUE, A QR FACTORIZATION IS USED FIRST TO REDUCE
*               THE MATRIX TO A TRIANGULAR FORM.)
*          = 7: THE NUMBER OF PROCESSORS
*          = 8: THE CROSSOVER POINT FOR THE MULTISHIFT QR AND QZ METHODS
*               FOR NONSYMMETRIC EIGENVALUE PROBLEMS.
*          = 9: MAXIMUM SIZE OF THE SUBPROBLEMS AT THE BOTTOM OF THE
*               COMPUTATION TREE IN THE DIVIDE-AND-CONQUER ALGORITHM
*               (USED BY XGELSD AND XGESDD)
*          =10: IEEE NAN ARITHMETIC CAN BE TRUSTED NOT TO TRAP
*          =11: INFINITY ARITHMETIC CAN BE TRUSTED NOT TO TRAP
*
*  NAME    (INPUT) CHARACTER*(*)
*          THE NAME OF THE CALLING SUBROUTINE, IN EITHER UPPER CASE OR
*          LOWER CASE.
*
*  OPTS    (INPUT) CHARACTER*(*)
*          THE CHARACTER OPTIONS TO THE SUBROUTINE NAME, CONCATENATED
*          INTO A SINGLE CHARACTER STRING.  FOR EXAMPLE, UPLO = 'U',
*          TRANS = 'T', AND DIAG = 'N' FOR A TRIANGULAR ROUTINE WOULD
*          BE SPECIFIED AS OPTS = 'UTN'.
*
*  N1      (INPUT) INTEGER
*  N2      (INPUT) INTEGER
*  N3      (INPUT) INTEGER
*  N4      (INPUT) INTEGER
*          PROBLEM DIMENSIONS FOR THE SUBROUTINE NAME; THESE MAY NOT ALL
*          BE REQUIRED.
*
* (ILAENV) (OUTPUT) INTEGER
*          >= 0: THE VALUE OF THE PARAMETER SPECIFIED BY ISPEC
*          < 0:  IF ILAENV = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE.
*
*  FURTHER DETAILS
*  ===============
*
*  THE FOLLOWING CONVENTIONS HAVE BEEN USED WHEN CALLING ILAENV FROM THE
*  LAPACK ROUTINES:
*  1)  OPTS IS A CONCATENATION OF ALL OF THE CHARACTER OPTIONS TO
*      SUBROUTINE NAME, IN THE SAME ORDER THAT THEY APPEAR IN THE
*      ARGUMENT LIST FOR NAME, EVEN IF THEY ARE NOT USED IN DETERMINING
*      THE VALUE OF THE PARAMETER SPECIFIED BY ISPEC.
*  2)  THE PROBLEM DIMENSIONS N1, N2, N3, N4 ARE SPECIFIED IN THE ORDER
*      THAT THEY APPEAR IN THE ARGUMENT LIST FOR NAME.  N1 IS USED
*      FIRST, N2 SECOND, AND SO ON, AND UNUSED PROBLEM DIMENSIONS ARE
*      PASSED A VALUE OF -1.
*  3)  THE PARAMETER VALUE RETURNED BY ILAENV IS CHECKED FOR VALIDITY IN
*      THE CALLING SUBROUTINE.  FOR EXAMPLE, ILAENV IS USED TO RETRIEVE
*      THE OPTIMAL BLOCKSIZE FOR STRTRI AS FOLLOWS:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. LOCAL SCALARS ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
*
*     INVALID VALUE FOR ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     CONVERT NAME TO UPPER CASE IF THE FIRST CHARACTER IS LOWER CASE.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII CHARACTER SET
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC CHARACTER SET
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        PRIME MACHINES:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  BLOCK SIZE
*
*     IN THESE EXAMPLES, SEPARATE CODE IS PROVIDED FOR SETTING NB FOR
*     REAL AND COMPLEX.  WE ASSUME THAT NB WILL TAKE THE SAME VALUE IN
*     SINGLE OR DOUBLE PRECISION.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  MINIMUM BLOCK SIZE
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
      END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  CROSSOVER POINT
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  NUMBER OF SHIFTS (USED BY XHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  MINIMUM COLUMN DIMENSION (NOT USED)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  CROSSOVER POINT FOR SVD (USED BY XGELSS AND XGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  NUMBER OF PROCESSORS (NOT USED)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  CROSSOVER POINT FOR MULTISHIFT (USED BY XHSEQR)
*
      ILAENV = 50
      RETURN
*
  900 CONTINUE
*
*     ISPEC = 9:  MAXIMUM SIZE OF THE SUBPROBLEMS AT THE BOTTOM OF THE
*                 COMPUTATION TREE IN THE DIVIDE-AND-CONQUER ALGORITHM
*                 (USED BY XGELSD AND XGESDD)
*
      ILAENV = 25
      RETURN
*
 1000 CONTINUE
*
*     ISPEC = 10: IEEE NAN ARITHMETIC CAN BE TRUSTED NOT TO TRAP
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 ) 
      END IF
      RETURN
*
 1100 CONTINUE
*
*     ISPEC = 11: INFINITY ARITHMETIC CAN BE TRUSTED NOT TO TRAP
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 ) 
      END IF
      RETURN
*
*     END OF ILAENV
*
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      LOGICAL FUNCTION LSAME( CA, CB )
C ----------------------------------------------------------------------
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     SEPTEMBER 30, 1994
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          CA, CB
*     ..
*
*  PURPOSE
*  =======
*
*  LSAME RETURNS .TRUE. IF CA IS THE SAME LETTER AS CB REGARDLESS OF
*  CASE.
*
*  ARGUMENTS
*  =========
*
*  CA      (INPUT) CHARACTER*1
*  CB      (INPUT) CHARACTER*1
*          CA AND CB SPECIFY THE SINGLE CHARACTERS TO BE COMPARED.
*
* =====================================================================
*
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ICHAR
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST IF THE CHARACTERS ARE EQUAL
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     NOW TEST FOR EQUIVALENCE IF BOTH CHARACTERS ARE ALPHABETIC.
*
      ZCODE = ICHAR( 'Z' )
*
*     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
*     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
*     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
*     ICHAR('A') ON AN EBCDIC MACHINE.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
*        UPPER CASE 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
*        UPPER CASE 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
*        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     END OF LSAME
*
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------      
      SUBROUTINE XERBLA( SRNAME, INFO )
C ----------------------------------------------------------------------
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     SEPTEMBER 30, 1994
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  PURPOSE
*  =======
*
*  XERBLA  IS AN ERROR HANDLER FOR THE LAPACK ROUTINES.
*  IT IS CALLED BY AN LAPACK ROUTINE IF AN INPUT PARAMETER HAS AN
*  INVALID VALUE.  A MESSAGE IS PRINTED AND EXECUTION STOPS.
*
*  INSTALLERS MAY CONSIDER MODIFYING THE STOP STATEMENT IN ORDER TO
*  CALL SYSTEM-SPECIFIC EXCEPTION-HANDLING FACILITIES.
*
*  ARGUMENTS
*  =========
*
*  SRNAME  (INPUT) CHARACTER*6
*          THE NAME OF THE ROUTINE WHICH CALLED XERBLA.
*
*  INFO    (INPUT) INTEGER
*          THE POSITION OF THE INVALID PARAMETER IN THE PARAMETER LIST
*          OF THE CALLING ROUTINE.
*
* =====================================================================
*
*     .. EXECUTABLE STATEMENTS ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** ON ENTRY TO ', A6, ' PARAMETER NUMBER ', I2, ' HAD ',
     $      'AN ILLEGAL VALUE' )
*
*     END OF XERBLA
*
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
C ----------------------------------------------------------------------
*     .. SCALAR ARGUMENTS ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
*
*     C := ALPHA*OP( A )*OP( B ) + BETA*C,
*
*  WHERE  OP( X ) IS ONE OF
*
*     OP( X ) = X   OR   OP( X ) = X',
*
*  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
*  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
*
*  PARAMETERS
*  ==========
*
*  TRANSA - CHARACTER*1.
*           ON ENTRY, TRANSA SPECIFIES THE FORM OF OP( A ) TO BE USED IN
*           THE MATRIX MULTIPLICATION AS FOLLOWS:
*
*              TRANSA = 'N' OR 'N',  OP( A ) = A.
*
*              TRANSA = 'T' OR 'T',  OP( A ) = A'.
*
*              TRANSA = 'C' OR 'C',  OP( A ) = A'.
*
*           UNCHANGED ON EXIT.
*
*  TRANSB - CHARACTER*1.
*           ON ENTRY, TRANSB SPECIFIES THE FORM OF OP( B ) TO BE USED IN
*           THE MATRIX MULTIPLICATION AS FOLLOWS:
*
*              TRANSB = 'N' OR 'N',  OP( B ) = B.
*
*              TRANSB = 'T' OR 'T',  OP( B ) = B'.
*
*              TRANSB = 'C' OR 'C',  OP( B ) = B'.
*
*           UNCHANGED ON EXIT.
*
*  M      - INTEGER.
*           ON ENTRY,  M  SPECIFIES  THE NUMBER  OF ROWS  OF THE  MATRIX
*           OP( A )  AND OF THE  MATRIX  C.  M  MUST  BE AT LEAST  ZERO.
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY,  N  SPECIFIES THE NUMBER  OF COLUMNS OF THE MATRIX
*           OP( B ) AND THE NUMBER OF COLUMNS OF THE MATRIX C. N MUST BE
*           AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  K      - INTEGER.
*           ON ENTRY,  K  SPECIFIES  THE NUMBER OF COLUMNS OF THE MATRIX
*           OP( A ) AND THE NUMBER OF ROWS OF THE MATRIX OP( B ). K MUST
*           BE AT LEAST  ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, KA ), WHERE KA IS
*           K  WHEN  TRANSA = 'N' OR 'N',  AND IS  M  OTHERWISE.
*           BEFORE ENTRY WITH  TRANSA = 'N' OR 'N',  THE LEADING  M BY K
*           PART OF THE ARRAY  A  MUST CONTAIN THE MATRIX  A,  OTHERWISE
*           THE LEADING  K BY M  PART OF THE ARRAY  A  MUST CONTAIN  THE
*           MATRIX A.
*           UNCHANGED ON EXIT.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. WHEN  TRANSA = 'N' OR 'N' THEN
*           LDA MUST BE AT LEAST  MAX( 1, M ), OTHERWISE  LDA MUST BE AT
*           LEAST  MAX( 1, K ).
*           UNCHANGED ON EXIT.
*
*  B      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDB, KB ), WHERE KB IS
*           N  WHEN  TRANSB = 'N' OR 'N',  AND IS  K  OTHERWISE.
*           BEFORE ENTRY WITH  TRANSB = 'N' OR 'N',  THE LEADING  K BY N
*           PART OF THE ARRAY  B  MUST CONTAIN THE MATRIX  B,  OTHERWISE
*           THE LEADING  N BY K  PART OF THE ARRAY  B  MUST CONTAIN  THE
*           MATRIX B.
*           UNCHANGED ON EXIT.
*
*  LDB    - INTEGER.
*           ON ENTRY, LDB SPECIFIES THE FIRST DIMENSION OF B AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. WHEN  TRANSB = 'N' OR 'N' THEN
*           LDB MUST BE AT LEAST  MAX( 1, K ), OTHERWISE  LDB MUST BE AT
*           LEAST  MAX( 1, N ).
*           UNCHANGED ON EXIT.
*
*  BETA   - DOUBLE PRECISION.
*           ON ENTRY,  BETA  SPECIFIES THE SCALAR  BETA.  WHEN  BETA  IS
*           SUPPLIED AS ZERO THEN C NEED NOT BE SET ON INPUT.
*           UNCHANGED ON EXIT.
*
*  C      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDC, N ).
*           BEFORE ENTRY, THE LEADING  M BY N  PART OF THE ARRAY  C MUST
*           CONTAIN THE MATRIX  C,  EXCEPT WHEN  BETA  IS ZERO, IN WHICH
*           CASE C NEED NOT BE SET ON ENTRY.
*           ON EXIT, THE ARRAY  C  IS OVERWRITTEN BY THE  M BY N  MATRIX
*           ( ALPHA*OP( A )*OP( B ) + BETA*C ).
*
*  LDC    - INTEGER.
*           ON ENTRY, LDC SPECIFIES THE FIRST DIMENSION OF C AS DECLARED
*           IN  THE  CALLING  (SUB)  PROGRAM.   LDC  MUST  BE  AT  LEAST
*           MAX( 1, M ).
*           UNCHANGED ON EXIT.
*
*
*  LEVEL 3 BLAS ROUTINE.
*
*  -- WRITTEN ON 8-FEBRUARY-1989.
*     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
*     IAIN DUFF, AERE HARWELL.
*     JEREMY DU CROZ, NUMERICAL ALGORITHMS GROUP LTD.
*     SVEN HAMMARLING, NUMERICAL ALGORITHMS GROUP LTD.
*
*
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     .. LOCAL SCALARS ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     SET  NOTA  AND  NOTB  AS  TRUE IF  A  AND  B  RESPECTIVELY ARE NOT
*     TRANSPOSED AND SET  NROWA, NCOLA AND  NROWB  AS THE NUMBER OF ROWS
*     AND  COLUMNS OF  A  AND THE  NUMBER OF  ROWS  OF  B  RESPECTIVELY.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     AND IF  ALPHA.EQ.ZERO.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     START THE OPERATIONS.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           FORM  C := ALPHA*A*B + BETA*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           FORM  C := ALPHA*A'*B + BETA*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           FORM  C := ALPHA*A*B' + BETA*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           FORM  C := ALPHA*A'*B' + BETA*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     END OF DGEMM .
*
      END
 
      
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C ----------------------------------------------------------------------
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGER   PERFORMS THE RANK 1 OPERATION
*
*     A := ALPHA*X*Y' + A,
*
*  WHERE ALPHA IS A SCALAR, X IS AN M ELEMENT VECTOR, Y IS AN N ELEMENT
*  VECTOR AND A IS AN M BY N MATRIX.
*
*  PARAMETERS
*  ==========
*
*  M      - INTEGER.
*           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF THE MATRIX A.
*           M MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF THE MATRIX A.
*           N MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
*           UNCHANGED ON EXIT.
*
*  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( M - 1 )*ABS( INCX ) ).
*           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE M
*           ELEMENT VECTOR X.
*           UNCHANGED ON EXIT.
*
*  INCX   - INTEGER.
*           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           X. INCX MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCY ) ).
*           BEFORE ENTRY, THE INCREMENTED ARRAY Y MUST CONTAIN THE N
*           ELEMENT VECTOR Y.
*           UNCHANGED ON EXIT.
*
*  INCY   - INTEGER.
*           ON ENTRY, INCY SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           Y. INCY MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
*           BEFORE ENTRY, THE LEADING M BY N PART OF THE ARRAY A MUST
*           CONTAIN THE MATRIX OF COEFFICIENTS. ON EXIT, A IS
*           OVERWRITTEN BY THE UPDATED MATRIX.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
*           MAX( 1, M ).
*           UNCHANGED ON EXIT.
*
*
*  LEVEL 2 BLAS ROUTINE.
*
*  -- WRITTEN ON 22-OCTOBER-1986.
*     JACK DONGARRA, ARGONNE NATIONAL LAB.
*     JEREMY DU CROZ, NAG CENTRAL OFFICE.
*     SVEN HAMMARLING, NAG CENTRAL OFFICE.
*     RICHARD HANSON, SANDIA NATIONAL LABS.
*
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
*     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     END OF DGER  .
*
      END
      

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE  DSCAL(N,DA,DX,INCX)
C ----------------------------------------------------------------------
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C     MODIFIED 3/93 TO RETURN IF INCX .LE. 0.
C     MODIFIED 12/3/93, ARRAY(1) DECLARATIONS CHANGED TO ARRAY(*)
C
      DOUBLE PRECISION DA,DX(*)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF( N.LE.0 .OR. INCX.LE.0 )RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
C ----------------------------------------------------------------------
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C     MODIFIED 12/3/93, ARRAY(1) DECLARATIONS CHANGED TO ARRAY(*)
C
      DOUBLE PRECISION DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END
      

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
C ----------------------------------------------------------------------
*     .. SCALAR ARGUMENTS ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DTRSM  SOLVES ONE OF THE MATRIX EQUATIONS
*
*     OP( A )*X = ALPHA*B,   OR   X*OP( A ) = ALPHA*B,
*
*  WHERE ALPHA IS A SCALAR, X AND B ARE M BY N MATRICES, A IS A UNIT, OR
*  NON-UNIT,  UPPER OR LOWER TRIANGULAR MATRIX  AND  OP( A )  IS ONE  OF
*
*     OP( A ) = A   OR   OP( A ) = A'.
*
*  THE MATRIX X IS OVERWRITTEN ON B.
*
*  PARAMETERS
*  ==========
*
*  SIDE   - CHARACTER*1.
*           ON ENTRY, SIDE SPECIFIES WHETHER OP( A ) APPEARS ON THE LEFT
*           OR RIGHT OF X AS FOLLOWS:
*
*              SIDE = 'L' OR 'L'   OP( A )*X = ALPHA*B.
*
*              SIDE = 'R' OR 'R'   X*OP( A ) = ALPHA*B.
*
*           UNCHANGED ON EXIT.
*
*  UPLO   - CHARACTER*1.
*           ON ENTRY, UPLO SPECIFIES WHETHER THE MATRIX A IS AN UPPER OR
*           LOWER TRIANGULAR MATRIX AS FOLLOWS:
*
*              UPLO = 'U' OR 'U'   A IS AN UPPER TRIANGULAR MATRIX.
*
*              UPLO = 'L' OR 'L'   A IS A LOWER TRIANGULAR MATRIX.
*
*           UNCHANGED ON EXIT.
*
*  TRANSA - CHARACTER*1.
*           ON ENTRY, TRANSA SPECIFIES THE FORM OF OP( A ) TO BE USED IN
*           THE MATRIX MULTIPLICATION AS FOLLOWS:
*
*              TRANSA = 'N' OR 'N'   OP( A ) = A.
*
*              TRANSA = 'T' OR 'T'   OP( A ) = A'.
*
*              TRANSA = 'C' OR 'C'   OP( A ) = A'.
*
*           UNCHANGED ON EXIT.
*
*  DIAG   - CHARACTER*1.
*           ON ENTRY, DIAG SPECIFIES WHETHER OR NOT A IS UNIT TRIANGULAR
*           AS FOLLOWS:
*
*              DIAG = 'U' OR 'U'   A IS ASSUMED TO BE UNIT TRIANGULAR.
*
*              DIAG = 'N' OR 'N'   A IS NOT ASSUMED TO BE UNIT
*                                  TRIANGULAR.
*
*           UNCHANGED ON EXIT.
*
*  M      - INTEGER.
*           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF B. M MUST BE AT
*           LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF B.  N MUST BE
*           AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY,  ALPHA SPECIFIES THE SCALAR  ALPHA. WHEN  ALPHA IS
*           ZERO THEN  A IS NOT REFERENCED AND  B NEED NOT BE SET BEFORE
*           ENTRY.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, K ), WHERE K IS M
*           WHEN  SIDE = 'L' OR 'L'  AND IS  N  WHEN  SIDE = 'R' OR 'R'.
*           BEFORE ENTRY  WITH  UPLO = 'U' OR 'U',  THE  LEADING  K BY K
*           UPPER TRIANGULAR PART OF THE ARRAY  A MUST CONTAIN THE UPPER
*           TRIANGULAR MATRIX  AND THE STRICTLY LOWER TRIANGULAR PART OF
*           A IS NOT REFERENCED.
*           BEFORE ENTRY  WITH  UPLO = 'L' OR 'L',  THE  LEADING  K BY K
*           LOWER TRIANGULAR PART OF THE ARRAY  A MUST CONTAIN THE LOWER
*           TRIANGULAR MATRIX  AND THE STRICTLY UPPER TRIANGULAR PART OF
*           A IS NOT REFERENCED.
*           NOTE THAT WHEN  DIAG = 'U' OR 'U',  THE DIAGONAL ELEMENTS OF
*           A  ARE NOT REFERENCED EITHER,  BUT ARE ASSUMED TO BE  UNITY.
*           UNCHANGED ON EXIT.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM.  WHEN  SIDE = 'L' OR 'L'  THEN
*           LDA  MUST BE AT LEAST  MAX( 1, M ),  WHEN  SIDE = 'R' OR 'R'
*           THEN LDA MUST BE AT LEAST MAX( 1, N ).
*           UNCHANGED ON EXIT.
*
*  B      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDB, N ).
*           BEFORE ENTRY,  THE LEADING  M BY N PART OF THE ARRAY  B MUST
*           CONTAIN  THE  RIGHT-HAND  SIDE  MATRIX  B,  AND  ON EXIT  IS
*           OVERWRITTEN BY THE SOLUTION MATRIX  X.
*
*  LDB    - INTEGER.
*           ON ENTRY, LDB SPECIFIES THE FIRST DIMENSION OF B AS DECLARED
*           IN  THE  CALLING  (SUB)  PROGRAM.   LDB  MUST  BE  AT  LEAST
*           MAX( 1, M ).
*           UNCHANGED ON EXIT.
*
*
*  LEVEL 3 BLAS ROUTINE.
*
*
*  -- WRITTEN ON 8-FEBRUARY-1989.
*     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
*     IAIN DUFF, AERE HARWELL.
*     JEREMY DU CROZ, NUMERICAL ALGORITHMS GROUP LTD.
*     SVEN HAMMARLING, NUMERICAL ALGORITHMS GROUP LTD.
*
*
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     .. LOCAL SCALARS ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     AND WHEN  ALPHA.EQ.ZERO.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     START THE OPERATIONS.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           FORM  B := ALPHA*INV( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           FORM  B := ALPHA*INV( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           FORM  B := ALPHA*B*INV( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           FORM  B := ALPHA*B*INV( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     END OF DTRSM .
*
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C ----------------------------------------------------------------------
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C     MODIFIED 3/93 TO RETURN IF INCX .LE. 0.
C     MODIFIED 12/3/93, ARRAY(1) DECLARATIONS CHANGED TO ARRAY(*)
C
      DOUBLE PRECISION DX(*),DMAX
      INTEGER I,INCX,IX,N
C
      IDAMAX = 0
      IF( N.LT.1 .OR. INCX.LE.0 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
