MODULE SVD
IMPLICIT NONE

    CONTAINS
    
    SUBROUTINE SVDCMP_ROUTINE(A,M,N,MP,NP,W,DUMMY)
       IMPLICIT NONE
       INTEGER, PARAMETER :: DBP = SELECTED_REAL_KIND (15,307)
   
       INTEGER :: M,N,MP,NP
       REAL(DBP) :: A(MP,NP), W(NP), DUMMY(MP,NP)
       
       DUMMY(:,:) = A(:,:)
       CALL SVDCMP(DUMMY,M,N,MP,NP,W,A)
       
    END SUBROUTINE
    
    SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
       IMPLICIT NONE
       INTEGER, PARAMETER :: DBP = SELECTED_REAL_KIND (15,307)
   
       INTEGER :: M,N,MP,NP, I,J,L,K, ITS, NM, JJ
       REAL(DBP) :: G, SSCALE, ANORM, S, F,H,X,Z,Y,C

       REAL(DBP), PARAMETER :: EPS = 3.0D-15

       INTEGER, PARAMETER :: NMAX = 1000        !MAXIMUM ANTICIPATED VALUE OF N
       REAL(DBP) :: A(MP,NP), W(NP), V(NP,NP), RV1(NMAX)


    !   PRINT *, 'ENTER THE TOLERANCE OR PRECISION REQUIRED'
    !   READ *, EPS
    
       !PRINT *, 'PRECISION CHOSEN AS ',EPS

       IF (M.LT.N) THEN
   	    PRINT *, 'YOU MUST AUGMENT A WITH EXTRA ZERO ROWS'
          CALL EXIT(10)
       ENDIF 

                !HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM
    !(SEE FORSYTHE,MALCOLM,MOLER, "COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS"
       G=0.0D0
       SSCALE = 0.0D0
       ANORM = 0.0D0
       DO I = 1,N
          L = I + 1
          RV1(I) = SSCALE*G
          G = 0.0D0
          S = 0.0D0
          SSCALE = 0.0D0
          IF (I.LE.M) THEN
             DO K = I,M
                SSCALE = SSCALE + DABS(A(K,I))
             END DO       ! K LOOP
    !         IF (SSCALE.NE.0.0D0) THEN
			    IF ( DABS(SSCALE-0.0D0).GT.EPS ) THEN
                DO K = I,M
                   A(K,I) = A(K,I) / SSCALE
                   S = S + A(K,I)*A(K,I)
                END DO    ! K LOOP
                F = A(I,I)
                G = - SIGN(SQRT(S),F)
                H = F*G - S
                A(I,I) = F - G
                IF (I.NE.N) THEN
                   DO J = L,N
                      S = 0.0D0
                      DO K = I,M
                         S = S + A(K,I)*A(K,J)
                      END DO      ! K LOOP
                      F = S / H
                      DO K = I, M 
                         A(K,J) = A(K,J) + F*A(K,I)
                      END DO   ! K LOOP
                   END DO      ! J LOOP
                END IF
                DO K = I, M 
                   A(K,I) = SSCALE * A(K,I)
                END DO         ! K LOOP
             END IF
          END IF

          W(I) = SSCALE * G
          G = 0.0D0
          S = 0.0D0
          SSCALE = 0.0D0
          IF ((I.LE.M).AND.(I.NE.N)) THEN
             DO K = L, N
                SSCALE = SSCALE + DABS(A(I,K))
             END DO         ! K LOOP
    !         IF (SSCALE.NE.0.0D0) THEN
			    IF ( DABS(SSCALE-0.0D0).GT.EPS ) THEN
                DO K = L, N
                   A(I,K) = A(I,K) /SSCALE
                   S = S + A(I,K) * A(I,K)
                END DO      ! K LOOP 
                F = A(I,L) 
                G = - SIGN(SQRT(S),F)
                H = F * G - S
                A(I,L) = F - G
                DO K = L, N
                   RV1(K) = A(I,K) / H
                END DO      ! K LOOP
                IF (I.NE.M) THEN
                   DO J = L, M 
                      S = 0.0D0
                      DO K = L, N 
                         S = S + A(J,K)*A(I,K)
                      END DO   ! K LOOP
                      DO K = L, N 
                         A(J,K) = A(J,K) + S*RV1(K)
                      END DO   ! K LOOP
                   END DO      ! J LOOP
                END IF
				    DO K = L, N
                   A(I,K) = SSCALE * A(I,K)
           	    END DO
             END IF
          END IF
          ANORM = MAX(ANORM, (DABS(W(I)) + DABS(RV1(I))))
       END DO

    ! ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS
       DO I = N, 1, -1 
          IF (I.LT.N) THEN
    !         IF (G.NE.0.0D0) THEN
			    IF ( DABS(G-0.0D0).GT.EPS ) THEN
                DO J = L, N       ! DOUBLE DIVISION TO AVOID POSSIBLE OVERFLOW
                   V(J,I) = (A(I,J) / A(I,L)) / G
                END DO      ! J LOOP
                DO J = L, N
                   S = 0.0D0
                   DO K = L, N
                      S = S + A(I,K)*V(K,J)
                   END DO   ! K LOOP
                   DO K = L, N
                      V(K,J) = V(K,J) + S * V(K,I)
                   END DO   ! K LOOP
           	    END DO      ! J LOOP
             END IF
             DO J = L, N 
                V(I,J) = 0.0D0
                V(J,I) = 0.0D0
             END DO
          END IF
          V(I,I) = 1.0D0
          G = RV1(I)
          L = I
       END DO

    ! ACCUMULATION OF LEFT-HAND TRANSFORMATIONS
       DO I = N, 1, -1
          L = 1 + I
          G = W(I)
          IF (I.LT.N) THEN
             DO J = L, N
                A(I,J) = 0.0D0
             END DO
          END IF
    !      IF (G.NE.0.0D0) THEN
          IF ( DABS(G-0.0D0).GT.EPS ) THEN
             G = 1.0D0 / G
             IF (I.NE.N) THEN
                DO J = L,N 
                   S = 0.0D0
                   DO K = L, M
                      S = S + A(K,I)*A(K,J)
                   END DO   ! K LOOP
                   F = (S/A(I,I)) * G
                   DO K = I, M 
                      A(K,J) = A(K,J) + F * A(K,I)
                   END DO   ! K LOOP
                END DO      ! J LOOP
             END IF
             DO J = I, M 
                A(J,I) = A(J,I) * G
             END DO         ! J LOOP
          ELSE
             DO J = I, M
                A(J,I) = 0.0D0
             END DO         ! J LOOP
          END IF
          A(I,I) = A(I,I) + 1.0D0
       END DO               ! I LOOP

    ! DIAGONALIZATION OF THE BIDIGONAL FORM
       DO K = N, 1, -1                  !LOOP OVER SINGULAR VALUES
          DO ITS = 1,30                 !LOOP OVER ALLOWED ITERATIONS
             DO L = K, 1, -1            !TEST FOR SPLITTING
                NM = L - 1              ! NOTE THAT RV1(1) IS ALWAYS ZERO
    !           IF ( (DABS(RV1(L))+ANORM) .EQ. ANORM ) GO TO 2
    !          	IF ( (DABS(W(NM))+ANORM) .EQ. ANORM ) GO TO 1
                IF ( DABS((DABS(RV1(L))+ANORM) - ANORM).LT.EPS ) GO TO 2
          	    IF ( DABS((DABS(W(NM))+ANORM) - ANORM).LT.EPS ) GO TO 1
             END DO      !  L LOOP

1            C = 0.0D0                  ! CANCELLATION OF RV1(L), IF L>1 :
             S = 1.0D0
             DO I = L, K
                F = S * RV1(I)
    !            IF ( (DABS(F)+ANORM) .NE. ANORM ) THEN
                IF ( DABS( (DABS(F)+ANORM) - ANORM) .GT. EPS ) THEN

                   G = W(I)
                   H = SQRT(F*F + G*G)
                   W(I) = H
                   H = 1.0D0 / H
                   C = G * H
                   S = -F * H
                   DO J = 1, M
                      Y = A(J,NM)
                      Z = A(J,I)
                      A(J,NM) = (Y*C) + (Z*S)
                      A(J,I) = -(Y*S) + (Z*C)
                   END DO   ! J LOOP
                END IF
             END DO         ! I LOOP
2            Z = W(K) 
             IF (L .EQ. K) THEN         ! CONVERGENCE
				    IF (Z .LT. 0.0D0) THEN  ! SINGULAR VALUE IS MADE NON-NEGATIVE
                   W(K) = -Z
                   DO J = 1,N
                      V(J,K) = -V(J,K)
                   END DO         ! J LOOP
	    		    END IF
                GO TO 3
             END IF
             IF (ITS.EQ.300) THEN
         	    PRINT*, 'NO CONVERGENCE IN 30 ITERATIONS'
                CALL EXIT(10)
             ENDIF
             X = W(L)          ! SHIFT FROM BOTTOM 2-BY-2 MINOR
             NM = K - 1
             Y = W(NM)
             G = RV1(NM)
             H = RV1(K)
             F = ( (Y-Z)*(Y+Z) + (G-H)*(G+H) ) / ( 2.0D0*H*Y)
             G = SQRT(F*F + 1.0D0)
             F = ( (X-Z)*(X+Z) + H*((Y/(F+SIGN(G,F))) - H) ) / X

    ! NEXT   QR TRANSFORMATION
             C = 1.0D0
             S = 1.0D0
             DO J = L, NM
          	    I = J + 1
                G = RV1(I)
                Y = W(I)
                H = S*G
                G = C*G
                Z = SQRT(F*F + H*H)
                RV1(J) = Z
                C = F/Z
                S = H/Z
                F = (X*C) + (G*S)
                G = -(X*S) + (G*C)
                H = Y*S
                Y = Y*C
                DO JJ = 1, N
                   X = V(JJ,J)
                   Z = V(JJ,I)
                   V(JJ,J) = (X*C) + (Z*S)
                   V(JJ,I) = -(X*S) + (Z*C)
           	    END DO
                Z = SQRT(F*F + H*H)
                W(J) = Z
    !            IF (Z.NE.0.0D0) THEN
                IF (  DABS(Z-0.0D0).GT.EPS  ) THEN
                   Z = 1.0D0 / Z
                   C = F*Z
                   S = H*Z
                END IF
                F = (G*C) + (Y*S)
                X = -(G*S) + (Y*C)
                DO JJ = 1, M
                   Y = A(JJ,J)
                   Z = A(JJ,I)
                   A(JJ,J) = (Y*C) + (Z*S)
                   A(JJ,I) = -(Y*S) + (Z*C)
                END DO
             END DO         ! J LOOP
             RV1(L) = 0.0D0
             RV1(K) = F
             W(K) = X
            END DO            ! ITS LOOP
3       CONTINUE
        END DO               ! K LOOP

        RETURN
    END SUBROUTINE SVDCMP
    
END MODULE SVD
