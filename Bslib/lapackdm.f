      SUBROUTINE DLARFGDUM (N, ALPHA, X, INCX, TAU)
C
      INTEGER          N, INCX
C
      DOUBLE PRECISION ALPHA(*), X(*), TAU(*)
C
      CALL DLARFG_SM (N, ALPHA, X, INCX, TAU)
C
      RETURN
      END
C
C
      SUBROUTINE DSYEVDUM (ITYPE, N, A, LDA, W, WORK, LWORK, INFO)
C
      INTEGER          ITYPE, INFO, LDA, LWORK, N
C
      DOUBLE PRECISION A( LDA, * ), W( * ), WORK( * )
C
      IF (ITYPE .EQ. 0) THEN
         CALL DSYEV_SM ("V", "U", N, A, LDA, W, WORK, LWORK, INFO)
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE DGESVDDUM (ITYPE, M, N, A, LDA, 
     *                      SV, U, LDU, V, LDV, 
     *                      WORK, LWORK, INFO)
C
      INTEGER          ITYPE, M, N, LDA, LDU, LDV, LWORK, INFO
C
      DOUBLE PRECISION A(LDA,*), SV(*), U(LDU,*), V(LDV,*), WORK(*)
C
      IF (ITYPE .EQ. 0) THEN
         CALL DGESVD_SM ("N", "N", M, N, 
     *                A, LDA, 
     *                SV, U, LDU, V, LDV, 
     *                WORK, LWORK, INFO)
      ELSE IF (ITYPE .EQ. 1) THEN
         CALL DGESVD_SM ("A", "A", M, N, 
     *                A, LDA, 
     *                SV, U, LDU, V, LDV, 
     *                WORK, LWORK, INFO)
      ELSE IF (ITYPE .EQ. 2) THEN
         CALL DGESVD_SM ("O", "A", M, N, 
     *                A, LDA, 
     *                SV, U, LDU, V, LDV, 
     *                WORK, LWORK, INFO)
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE ZLARFGDUM (N, ALPHA, X, INCX, TAU)
C
      INTEGER          N, INCX
C
      DOUBLE COMPLEX   ALPHA(*), X(*), TAU(*)
C
      CALL ZLARFG_SM (N, ALPHA, X, INCX, TAU)
C
      RETURN
      END
C
CC
C      SUBROUTINE ZHEEVDUM (ITYPE, N, A, LDA, W, 
C     *                     WORK, LWORK, RWORK, INFO)
CC
C      INTEGER          ITYPE, INFO, LDA, LWORK, N
CC
C      DOUBLE COMPLEX   A( LDA, * ), WORK( * )
C      DOUBLE PRECISION W( * ), RWORK( * )
CC
C      IF (ITYPE .EQ. 0) THEN
C         CALL ZHEEV_SM ("V", "U", N, A, LDA, W, 
C     *                  WORK, LWORK, RWORK, INFO)
C      END IF
CC
C      RETURN
C      END
CC
CC
      SUBROUTINE ZGESVDDUM (ITYPE, M, N, A, LDA, 
     *                      SV, U, LDU, V, LDV, 
     *                      WORK, LWORK, RWORK, INFO)
C
      INTEGER          ITYPE, M, N, LDA, LDU, LDV, LWORK, INFO
C
      DOUBLE COMPLEX   A(LDA,*), U(LDU,*), V(LDV,*), WORK(*)
      DOUBLE PRECISION SV(*), RWORK(*)
C
      IF (ITYPE .EQ. 0) THEN
         CALL ZGESVD_SM ("N", "N", M, N, 
     *                A, LDA, 
     *                SV, U, LDU, V, LDV, 
     *                WORK, LWORK, RWORK, INFO)
      ELSE IF (ITYPE .EQ. 1) THEN
         CALL ZGESVD_SM ("A", "A", M, N, 
     *                A, LDA, 
     *                SV, U, LDU, V, LDV, 
     *                WORK, LWORK, RWORK, INFO)
      ELSE IF (ITYPE .EQ. 2) THEN
         CALL ZGESVD_SM ("O", "A", M, N, 
     *                A, LDA, 
     *                SV, U, LDU, V, LDV, 
     *                WORK, LWORK, RWORK, INFO)
      END IF
C
      RETURN
      END

