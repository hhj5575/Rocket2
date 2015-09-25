MODULE HIGHORDER_3D
    USE SURFACES_3D
    USE OPERATORS_3D
    
    IMPLICIT NONE
    
    TYPE SPLINE3
        INTEGER :: D, N
        REAL(8), ALLOCATABLE :: T(:) ! T(N)
        REAL(8), ALLOCATABLE :: X(:,:) ! X(D,N)
        REAL(8), ALLOCATABLE :: COEFFS(:,:,:) ! COEFFS(4,D,N-1)
    END TYPE
    
    CONTAINS
    
    !! SOLVERS
    
    SUBROUTINE INDEX_SPARSE(A_INDICES,I0,J0, J)
	IMPLICIT NONE
        INTEGER :: A_INDICES(:,:)
        INTEGER :: I0,J0
        
        INTEGER :: J
        
        J = 1
        DO WHILE(.TRUE.)
            IF(A_INDICES(I0,J)==J0) THEN
                EXIT
            END IF
            
            J = J+1
        END DO
        
    END SUBROUTINE
    
    SUBROUTINE MATMUL_SPARSE(A,A_INDICES,B,N, C)
	IMPLICIT NONE
        REAL(8) :: A(:,:)
        INTEGER :: A_INDICES(:,:)
        
        INTEGER :: N
        REAL(8) :: B(:)
        REAL(8) :: C(:)
        
        INTEGER :: I,J
        
        DO I=1,N
            C(I) = 0.
            
            DO J=1,N
                IF(A_INDICES(I,J)==0) THEN
                    EXIT
                END IF
                
                C(I) = C(I) + A(I,J) * B(A_INDICES(I,J))
            END DO
	END DO
    END SUBROUTINE
    
    SUBROUTINE CG_SOLVER_SPARSE(A, A_INDICES, B, N, X)
	IMPLICIT NONE
        REAL(8) :: A(:,:)
        INTEGER :: A_INDICES(:,:)
        REAL(8) :: B(N)
        REAL(8) :: X(N)
        REAL(8), ALLOCATABLE :: D(:)
        REAL(8), ALLOCATABLE :: AD(:)
        REAL(8), ALLOCATABLE :: R(:)
        REAL(8), ALLOCATABLE :: OLDR(:)
        INTEGER :: N
        INTEGER :: MAX_ITER
        REAL(8) :: DELTA, BETA, ALPHA
        INTEGER :: K
        
        ALLOCATE(R(N), OLDR(N), D(N), AD(N))
        MAX_ITER = 1000
        DELTA = 1E-10
        X(:) = 0.
        D(:) = 0.
        
        CALL MATMUL_SPARSE(A,A_INDICES,X,N,AD)
        
        R = B - AD
        D = R
        
        DO K = 1, MAX_ITER
            
            IF (DOT_PRODUCT(R,R)< DELTA) THEN
                RETURN
            END IF
            
            OLDR = R
            
            CALL MATMUL_SPARSE(A,A_INDICES,D,N,AD)
            ALPHA = DOT_PRODUCT(R,R)/DOT_PRODUCT(D,AD)
            X = X + ALPHA*D
            R = R - ALPHA*AD
            
            BETA = DOT_PRODUCT(R,R)/DOT_PRODUCT(OLDR,OLDR)
            D = R + BETA*D
            
        END DO
        
	DEALLOCATE(R, OLDR, D, AD)
    END SUBROUTINE CG_SOLVER_SPARSE
    
    SUBROUTINE CG_SOLVER(A, B, N, X)
	IMPLICIT NONE
        REAL(8) :: A(N,N)
        REAL(8) :: B(N)
        REAL(8) :: X(N)
        REAL(8), ALLOCATABLE :: D(:)
        REAL(8), ALLOCATABLE :: AD(:)
        REAL(8), ALLOCATABLE :: R(:)
        REAL(8), ALLOCATABLE :: OLDR(:)
        INTEGER :: N
        INTEGER :: MAX_ITER
        REAL(8) :: DELTA, BETA, ALPHA
        INTEGER :: K
        
        ALLOCATE(R(N), OLDR(N), D(N), AD(N))
        MAX_ITER = 1000
        DELTA = 1E-10
        X(:) = 0.
        D(:) = 0.
        
        R = B - MATMUL(A,X)
        D = R
        
        DO K = 1, MAX_ITER
            
            IF (DOT_PRODUCT(R,R)< DELTA) THEN
                RETURN
            END IF
            
            OLDR = R
            
            AD = MATMUL(A,D)
            ALPHA = DOT_PRODUCT(R,R)/DOT_PRODUCT(D,AD)
            X = X + ALPHA*D
            R = R - ALPHA*AD
            
            BETA = DOT_PRODUCT(R,R)/DOT_PRODUCT(OLDR,OLDR)
            D = R + BETA*D
            
        END DO
        
	DEALLOCATE(R, OLDR, D, AD)
    END SUBROUTINE CG_SOLVER
        
    SUBROUTINE PCG_SOLVER(A, B, N, X)
	IMPLICIT NONE
        REAL(8), DIMENSION(:,:) :: A
        REAL(8), DIMENSION(:) :: B
        REAL(8), DIMENSION(:) :: X
        REAL(8), ALLOCATABLE :: D(:)
        REAL(8), ALLOCATABLE :: AD(:)
        REAL(8), ALLOCATABLE :: Z(:)
        REAL(8), ALLOCATABLE :: R(:)
        REAL(8), ALLOCATABLE :: OLDZ(:)
        REAL(8), ALLOCATABLE :: OLDR(:)
        REAL(8), ALLOCATABLE :: CHOL(:,:)
        INTEGER :: N
        INTEGER :: MAX_ITER
        REAL(8) :: DELTA, BETA, ALPHA
        INTEGER :: K
        
        ALLOCATE(R(N), OLDR(N), D(N), Z(N), OLDZ(N), CHOL(N,N), AD(N))
        MAX_ITER = 1000
        DELTA = 1E-10
        X(:) = 0.
        D(:) = 0.
        R = B - MATMUL(A,X)
        CALL CHOLESKY(N,A,CHOL)
        CALL LU(N,R,CHOL,TRANSPOSE(CHOL), Z)
        D = Z
        
        DO K = 1, MAX_ITER
            
            IF (DOT_PRODUCT(R,R)< DELTA) THEN
                RETURN
            END IF
            
            OLDR = R
            OLDZ = Z
            
            AD = MATMUL(A,D)
            ALPHA = DOT_PRODUCT(R,Z)/DOT_PRODUCT(D,AD)
            X = X + ALPHA*D
            R = R - ALPHA*AD
            
            CALL LU(N,R,CHOL,TRANSPOSE(CHOL), Z)
            
            BETA = DOT_PRODUCT(R,Z)/DOT_PRODUCT(OLDR,OLDZ)
            D = Z + BETA*D
            
        END DO
        
	DEALLOCATE(R, OLDR, D, Z, OLDZ, CHOL, AD)
    END SUBROUTINE PCG_SOLVER
    
    SUBROUTINE CHOLESKY(N,A,R)
	IMPLICIT NONE
        INTEGER :: N
        REAL(8), DIMENSION(:,:) :: A
        REAL(8), DIMENSION(:,:) :: R
        INTEGER :: I,K,S
        REAL(8) :: SUM
        
        R = 0.
        
        DO K = 1, N
            SUM = A(K,K)
            DO S = 1, K-1
                SUM = SUM - R(K,S)*R(K,S)
            END DO
            R(K,K) = SQRT(MAX(0.,SUM))
            
            DO I = K+1,N
                !IF (A(I,K) .NE. 0.) THEN
                    SUM = A(I,K)
                    DO S = 1, K-1
                        SUM = SUM - R(I,S)*R(K,S)
                    END DO
                    R(I,K) = SUM/R(K,K)
                !END IF
            END DO
        END DO
    END SUBROUTINE CHOLESKY
    
    SUBROUTINE LU(N,B,L,U,X)
	IMPLICIT NONE
        INTEGER :: N
        REAL(8), DIMENSION(:) :: B
        REAL(8), DIMENSION(:,:) :: L
        REAL(8), DIMENSION(:,:) :: U
        REAL(8), DIMENSION(:) :: X
        REAL(8), ALLOCATABLE :: Z(:)
        INTEGER :: I, J
        REAL(8) :: SUM
        
        ALLOCATE(Z(N))
        Z = 0.
        Z(1) = B(1)/L(1,1)
        
        DO I = 2, N
            SUM = B(I)
            DO J = 1, I-1
                SUM = SUM - L(I,J)*Z(J)
            END DO
            
            Z(I) = SUM/L(I,I)
        END DO
        
        X = 0.
        X(N) = Z(N)/U(N,N)
        DO I = N-1,1,-1
            SUM = Z(I)
            DO J= I+1,N
                SUM = SUM - U(I,J)*X(J)
            END DO
            X(I) = SUM/U(I,I)
        END DO
        DEALLOCATE(Z)
    END SUBROUTINE LU
    
    SUBROUTINE CUMSUM(A, A_SIZE,    B)
        IMPLICIT NONE
        REAL(8) :: A(:), B(:)
        INTEGER :: A_SIZE
        INTEGER :: I
        
        B(1) = A(1)
        
        DO I=2,A_SIZE
            B(I) = B(I-1) + A(I)
        END DO
    END SUBROUTINE
    
    SUBROUTINE TRIDIAGONAL(C, B, N, X)
        IMPLICIT NONE
        INTEGER :: N
        REAL(8) :: C(3,N)
        REAL(8) :: B(N)
        REAL(8) :: X(N)
        INTEGER :: I
        REAL(8) :: XMULT
        REAL(8), ALLOCATABLE :: A_TEMP(:), D_TEMP(:), C_TEMP(:), B_TEMP(:)
        
        ALLOCATE(A_TEMP(N-1))
        ALLOCATE(D_TEMP(N))
        ALLOCATE(C_TEMP(N-1))
        ALLOCATE(B_TEMP(N))
        
        A_TEMP(:) = C(1,2:N)
        D_TEMP(:) = C(2,:)
        C_TEMP(:) = C(3,1:N-1)
        B_TEMP(:) = B(:)
        
        DO I=2,N
            XMULT = A_TEMP(I-1)/D_TEMP(I-1)
            D_TEMP(I) = D_TEMP(I) - XMULT * C_TEMP(I-1)
            B_TEMP(I) = B_TEMP(I) - XMULT * B_TEMP(I-1)
        END DO
        
        X(N) = B_TEMP(N)/D_TEMP(N)
        DO I=N-1,1,-1
            X(I) = (B_TEMP(I)-C_TEMP(I)*X(I+1))/D_TEMP(I)
        END DO
        
        DEALLOCATE(B_TEMP)
        DEALLOCATE(D_TEMP)
        
    END SUBROUTINE
    
    
    !! High order for boundaries
    
    SUBROUTINE NEW_SPLINE(D,N,X,CS)
        IMPLICIT NONE
        INTEGER :: D,N
        REAL(8) :: X(D,N)
        TYPE(SPLINE3) :: CS
        
        CS%D = D
        CS%N = N
        ALLOCATE(CS%T(N))
        ALLOCATE(CS%X(D,N))
        ALLOCATE(CS%COEFFS(4,D,N-1))
        CS%X(:,:) = X(:,:)
         
        CALL CSVN(CS%D,CS%N,CS%X,CS)
    END SUBROUTINE
    
    SUBROUTINE DELETE_SPLINE(CS)
        IMPLICIT NONE
        TYPE(SPLINE3) :: CS
        
        DEALLOCATE(CS%T)
        DEALLOCATE(CS%X)
        DEALLOCATE(CS%COEFFS)
    END SUBROUTINE
    
    SUBROUTINE SAVESPLINE(CS, FILE_NUM)
        IMPLICIT NONE
        TYPE(SPLINE3) :: CS
        INTEGER :: FILE_NUM
        INTEGER :: I, J, K
        CHARACTER(500) ::STR, TEMPSTR
        CHARACTER(50000) :: STR2
        REAL(8), ALLOCATABLE :: DT(:)
        
        WRITE(STR, *), FILE_NUM
        STR = './output/surface/spline/spline_3d_' // TRIM(ADJUSTL(STR)) // '.txt'
        WRITE(*,*) TRIM(STR)
    
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        
        ALLOCATE(DT(CS%N-1))
        
        DT(:) = CS%T(2:CS%N) - CS%T(1:(CS%N-1)) ! DT(N-1), COEFFS(4,D,N-1)
        
        STR2 = ''
        TEMPSTR = ''
        DO I=1,CS%N-1
            WRITE(TEMPSTR, *) DT(I)
            STR2 = TRIM(ADJUSTL(STR2)) // ' ' // TRIM(ADJUSTL(TEMPSTR))
            TEMPSTR = ''
        END DO
        
        WRITE(21,'(A)') STR2(1:LEN_TRIM(STR2))
        
        DO K=1,CS%D
            DO J=4,1,-1
                STR2 = ''
                TEMPSTR = ''
                DO I=1,CS%N-1
                    WRITE(TEMPSTR, *) CS%COEFFS(J,K,I)
                    STR2 = TRIM(ADJUSTL(STR2)) // ' ' // TRIM(ADJUSTL(TEMPSTR))
                    TEMPSTR = ''
                END DO
                
                WRITE(21,'(A)') STR2(1:LEN_TRIM(STR2))
            END DO
        END DO
        
        DEALLOCATE(DT)
        
    END SUBROUTINE
    
    SUBROUTINE SPLINE_COMPUTE(CS, T0, X0)
        IMPLICIT NONE
        TYPE(SPLINE3) :: CS
        REAL(8) :: T0, X0(:)
        INTEGER :: I, I0
        REAL(8) :: T1
        
        IF(T0<CS%T(1)) THEN
            T0 = CS%T(1)
        ELSE IF(T0>CS%T(CS%N)) THEN
            T0 = CS%T(CS%N)
        END IF
        
        DO I=1,CS%N-1
            IF(T0 .LE. CS%T(I+1)) THEN
                I0 = I
                EXIT
            END IF
        END DO
        
        T1 = T0 - CS%T(I0)
        
        DO I=1,CS%D
            X0(I) = CS%COEFFS(1,I,I0) + CS%COEFFS(2,I,I0) * T1 + CS%COEFFS(3,I,I0) * T1**2 + CS%COEFFS(4,I,I0) * T1**3
        END DO
    END SUBROUTINE
    
    SUBROUTINE ONE_SPLINE_COMPUTE(CS, I0, WEIGHT1, WEIGHT2, X0)
        IMPLICIT NONE
        TYPE(SPLINE3) :: CS
        INTEGER :: I0
        REAL(8) :: WEIGHT1, WEIGHT2
        
        REAL(8) :: X0(:)
        REAL(8) :: T1
        
        T1 = WEIGHT1 * CS%T(I0) + WEIGHT2 * CS%T(I0+1) - CS%T(I0)
        
        X0(:) = CS%COEFFS(1,:,I0) + CS%COEFFS(2,:,I0) * T1 + CS%COEFFS(3,:,I0) * T1**2 + CS%COEFFS(4,:,I0) * T1**3
    END SUBROUTINE
    
    SUBROUTINE ONE_SPLINE_EXACT_LOCATION(CS, I0, WEIGHT1, WEIGHT2, V_RESULT)
        IMPLICIT NONE
        TYPE(SPLINE3) :: CS
        INTEGER :: I0
        REAL(8) :: WEIGHT1, WEIGHT2
        
        REAL(8) :: V_RESULT(:)
        REAL(8) :: T_WEIGHT1, T_WEIGHT2, A, B, TEMPWEIGHT1, TEMPWEIGHT2, WEIGHT_THRESHOLD, R1, R2
        INTEGER :: ITER
        
        T_WEIGHT1 = WEIGHT1
        T_WEIGHT2 = WEIGHT2
        A = 0.
        B = 1.
        
        DO ITER = 1,1000
            CALL ONE_SPLINE_COMPUTE(CS, I0, T_WEIGHT1, T_WEIGHT2, V_RESULT)
            R1 = SQRT(DOT_PRODUCT(V_RESULT-CS%X(:,I0),V_RESULT-CS%X(:,I0)))
            R2 = SQRT(DOT_PRODUCT(V_RESULT-CS%X(:,I0+1),V_RESULT-CS%X(:,I0+1)))
            
            TEMPWEIGHT1 = R2/(R1+R2)
            TEMPWEIGHT2 = R1/(R1+R2)
            
            WEIGHT_THRESHOLD = 1./50.
            
            IF(TEMPWEIGHT1 < WEIGHT1 + WEIGHT_THRESHOLD .AND. TEMPWEIGHT1 > WEIGHT1 - WEIGHT_THRESHOLD) THEN
                EXIT
            END IF
            
            IF(TEMPWEIGHT1 < WEIGHT1) THEN
                A = T_WEIGHT1
                T_WEIGHT1 = (B+T_WEIGHT1)/2.
                
                T_WEIGHT2 = 1 - T_WEIGHT1
            ELSE
                B = T_WEIGHT1
                T_WEIGHT1 = (A+T_WEIGHT1)/2.
                
                T_WEIGHT2 = 1 - T_WEIGHT1
            END IF
        END DO
    END SUBROUTINE
    
    SUBROUTINE SPLINE_MIDPOINT(D, N, X, I0, X0)
        IMPLICIT NONE
        INTEGER :: D
        INTEGER :: N
        REAL(8) :: X(D,N)
        INTEGER :: I0
        REAL(8) :: X0(D)
        TYPE(SPLINE3) :: CS
        
        CALL NEW_SPLINE(D, N, X, CS)
        
        CALL SPLINE_COMPUTE(CS, (CS%T(I0) + CS%T(I0+1))/2., X0)
        
        CALL DELETE_SPLINE(CS)
        
    END SUBROUTINE
    
    SUBROUTINE CSVN(D, N, X, CS)
        IMPLICIT NONE
        INTEGER :: D
        INTEGER :: N
        REAL(8) :: X(D,N)
        TYPE(SPLINE3) :: CS
        
        INTEGER :: ENDCONDS
        REAL(8), ALLOCATABLE :: DT(:), DX(:,:), T(:)
        
        IF(SQRT(DOT_PRODUCT(X(:,1) - X(:,N), X(:,1) - X(:,N))) < 1E-16) THEN
            ENDCONDS = 0
        ELSE
            ENDCONDS = 1
        END IF
        
        ALLOCATE(DT(N))
        DT(:) = 0.
        
        !IF(N==1) THEN
        !    DT = 0
        !ELSE
            ALLOCATE(DX(D, N-1))
            
            DX(:,:) = X(:,2:N) - X(:,1:N-1)
            DX(:,:) = DX(:,:)**2
            
            DT(2:N) = SUM(DX(:,:), 1)
            
            DEALLOCATE(DX)
        !END IF
        
        DT(:) = SQRT(SQRT(DT(:)))
        
        ALLOCATE(T(N))
        
        CALL CUMSUM(DT, N, T)
        
        CALL CSAPE(D, N, T, X, ENDCONDS, CS)
        
        DEALLOCATE(DT)
        DEALLOCATE(T)
    END SUBROUTINE
    
    SUBROUTINE CSAPE(D, N, X, Y, CONDS,    CS)
        IMPLICIT NONE
        INTEGER :: D, N
        REAL(8) :: X(N), Y(D,N)
        INTEGER :: CONDS
        TYPE(SPLINE3) :: CS
        REAL(8), ALLOCATABLE :: DX(:), DIVDIF(:,:), C(:,:), B(:,:), S(:,:), C3(:,:), C4(:,:)
        INTEGER :: I
        
        ALLOCATE(DX(N-1))
        ALLOCATE(DIVDIF(N-1,D))
        
        DX(:) = X(2:N)-X(1:N-1)
        DO I=1,D
            DIVDIF(:,I) = Y(I,2:N)-Y(I,1:N-1)
            DIVDIF(:,I) = DIVDIF(:,I)/DX(:)
        END DO
        
        IF(CONDS==1) THEN
            ALLOCATE(C(3,N))
            C(:,:) = 0.
            C(1,2:N-1) = DX(2:N-1)
            C(2,2:N-1) = 2.*(DX(2:N-1) + DX(1:N-2))
            C(3,2:N-1) = DX(1:N-2)
        ELSE
            ALLOCATE(C(N,N))
            C(:,:) = 0.
            DO I=2,N-1
                C(I,I-1) = DX(I)
                C(I,I) = 2.*(DX(I) + DX(I-1))
                C(I,I+1) = DX(I-1)
            END DO
        END IF
        
        ALLOCATE(B(N,D))
        B(:,:) = 0.
        
        DO I=1,D
            B(2:N-1,I) = 3.*(DX(2:N-1)*DIVDIF(1:N-2,I) + DX(1:N-2)*DIVDIF(2:N-1,I))
        END DO
        
        IF(CONDS==1) THEN
            C(2,1) = 2.
            C(3,1) = 1.
            B(1,:) = 3.*DIVDIF(1,:)
            
            C(1,N) = 1.
            C(2,N) = 2.
            B(N,:) = 3.*DIVDIF(N-1,:);
        ELSE
            C(1,1) = 1.
            C(1,N) = -1.
            
            C(N,1) = 2.*DX(N-1)
            C(N,2) = DX(N-1)
            
            C(N,N-1) = C(N,N-1) + DX(1)
            C(N,N) = C(N,N) + 2.*DX(1)
            
            B(N,:) = 3.*(DX(N-1)*DIVDIF(1,:) + DX(1)*DIVDIF(N-1,:))
        END IF
        
        ALLOCATE(S(N,D))
        
        IF(CONDS==1) THEN
            DO I=1,D
                CALL TRIDIAGONAL(C, B(:,I), N, S(:,I))
            END DO
        ELSE
            DO I=1,D
                CALL CG_SOLVER(C,B(:,I),N,S(:,I))
                !CALL PCG_SOLVER(C,B(:,I),N,S(:,I))
            END DO
        END IF
        
        ALLOCATE(C4(D,N-1))
        DO I=1,D
            C4(I,:) = (S(1:N-1,I) + S(2:N,I) - 2.*DIVDIF(1:N-1,I)) / DX(:)
        END DO
        
        ALLOCATE(C3(D,N-1))
        DO I=1,D
            C3(I,:) = (DIVDIF(1:N-1,I) - S(1:N-1,I)) / DX(:) - C4(I,:)
        END DO
        
        
        !CS%D = D
        !CS%N = N
        
        CS%T(:) = X(:)
        
        CS%COEFFS(1,:,:) = Y(:,1:N-1)
        
        DO I=1,D
            CS%COEFFS(2,I,:) = S(1:N-1,I)
        END DO
        
        CS%COEFFS(3,:,:) = C3(:,:)
        
        DO I=1,D
            CS%COEFFS(4,I,:) = C4(I,:)/DX(:)
        END DO
        
        DEALLOCATE(C3)
        DEALLOCATE(C4)
        DEALLOCATE(S)
        DEALLOCATE(B)
        DEALLOCATE(C)
        DEALLOCATE(DIVDIF)
        DEALLOCATE(DX)
        
    END SUBROUTINE
    
    SUBROUTINE RECONSTRUCTING_HIGH_ORDER_BOUNDARY(V, D, N, X, I0, V_RESULT)
        IMPLICIT NONE
        INTEGER :: D,N
        REAL(8) :: V(D)
        REAL(8) :: X(D,N)
        INTEGER :: I0
        REAL(8) :: V_RESULT(D)
        
        REAL(8) :: WEIGHT1, WEIGHT2
        
        TYPE(SPLINE3) :: CS
        
        CALL NEW_SPLINE(D, N, X, CS)
        
        CALL COORDINATE_EDGE_POINT(V, CS%X(:,I0), CS%X(:,I0+1), WEIGHT1, WEIGHT2)
        CALL ONE_SPLINE_EXACT_LOCATION(CS, I0, WEIGHT1, WEIGHT2, V_RESULT)
        
        CALL DELETE_SPLINE(CS)
        
    END SUBROUTINE
    
    
    !! High order for regions
    
    SUBROUTINE RECONSTRUCTING_HIGH_ORDER_SURFACE(POINT_NUM, POINT, FACE, FACE_INDEX, V, POLY_ORDER, RING_NUM, CONNECTION_NUM, CONNECTION, POINT_TYPE, V_RESULT)
        IMPLICIT NONE
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(:,:)
        INTEGER :: FACE(:,:)
        INTEGER :: FACE_INDEX
        REAL(8) :: V(3)
        INTEGER :: POLY_ORDER
        INTEGER :: RING_NUM
        INTEGER :: CONNECTION_NUM(:)
        INTEGER :: CONNECTION(:,:)
        INTEGER :: POINT_TYPE(:)
        REAL(8) :: V_RESULT(3)
        
        INTEGER :: NEIGHBOR_NUM, NEIGHBOR(3)
        REAL(8) :: WEIGHT(3)
        INTEGER :: I,J,K, ITER, NEIGHBOR_POINT
        
        INTEGER, ALLOCATABLE :: POINT_INDEX_CLUSTER(:)
        REAL(8), ALLOCATABLE :: POINT_CLUSTER(:,:)
        INTEGER :: POINT_CLUSTER_NUM, BEFORE_POINT_CLUSTER_NUM
        
        LOGICAL, ALLOCATABLE :: USED_POINT(:)
        REAL(8), ALLOCATABLE :: NEWPOINT(:,:)
        REAL(8) :: FITTING_RESULT(3), FITTING, FITTING_POINT(3)
        REAL(8) :: V1(3), V2(3), V3(3), N(3), NEW_AXIS(3,3), NEW_SCALE, NEW_V(3), NEW_V_TEMP(3)
        
        NEIGHBOR_NUM = 0
        NEIGHBOR(:) = 0
        WEIGHT(:) = 0.
        DO K=1,3
            IF(POINT_TYPE(FACE(K,FACE_INDEX))==1) THEN
                NEIGHBOR_NUM = NEIGHBOR_NUM + 1
                NEIGHBOR(NEIGHBOR_NUM) = FACE(K,FACE_INDEX)
            END IF
        END DO
        
        IF(NEIGHBOR_NUM==0) THEN
            V_RESULT = V
            RETURN
        ELSE IF(NEIGHBOR_NUM==1) THEN
            WEIGHT(1) = 1.
        ELSE IF(NEIGHBOR_NUM==2) THEN
            CALL COORDINATE_EDGE_POINT(V,POINT(:,NEIGHBOR(1)),POINT(:,NEIGHBOR(2)),WEIGHT(1),WEIGHT(2))
        ELSE
            CALL COORDINATE_FACE_POINT(V,POINT(:,NEIGHBOR(1)),POINT(:,NEIGHBOR(2)),POINT(:,NEIGHBOR(3)), WEIGHT(1), WEIGHT(2), WEIGHT(3))
        END IF
        
        
        V1 = POINT(:,FACE(2,FACE_INDEX)) - POINT(:,FACE(1,FACE_INDEX))
        V2 = POINT(:,FACE(3,FACE_INDEX)) - POINT(:,FACE(1,FACE_INDEX))
        V3 = POINT(:,FACE(3,FACE_INDEX)) - POINT(:,FACE(2,FACE_INDEX))
        
        ALLOCATE(NEWPOINT(3,POINT_NUM))
        
        CALL VEC_CURL1(V1, V2, N)
        NEW_AXIS(:,3) = N/SQRT(DOT_PRODUCT(N,N))
        NEW_AXIS(:,1) = V1/SQRT(DOT_PRODUCT(V1,V1))
        CALL VEC_CURL1(NEW_AXIS(:,3), NEW_AXIS(:,1), NEW_AXIS(:,2))
        NEW_SCALE = (SQRT(DOT_PRODUCT(V1,V1)) + SQRT(DOT_PRODUCT(V2,V2)) + SQRT(DOT_PRODUCT(V3,V3)))/3.
        
        DO I=1,POINT_NUM
            NEWPOINT(1,I) = DOT_PRODUCT(POINT(:,I), NEW_AXIS(:,1))/NEW_SCALE
            NEWPOINT(2,I) = DOT_PRODUCT(POINT(:,I), NEW_AXIS(:,2))/NEW_SCALE
            NEWPOINT(3,I) = DOT_PRODUCT(POINT(:,I), NEW_AXIS(:,3))/NEW_SCALE
        END DO
        
        NEW_V(1) = DOT_PRODUCT(V, NEW_AXIS(:,1))/NEW_SCALE
        NEW_V(2) = DOT_PRODUCT(V, NEW_AXIS(:,2))/NEW_SCALE
        NEW_V(3) = DOT_PRODUCT(V, NEW_AXIS(:,3))/NEW_SCALE
        
        FITTING_RESULT(:) = 0.
        
        ALLOCATE(POINT_INDEX_CLUSTER(POINT_NUM), USED_POINT(POINT_NUM))
        
        DO ITER=1,NEIGHBOR_NUM
            POINT_INDEX_CLUSTER(:) = 0
            USED_POINT(:) = .FALSE.
            
            POINT_CLUSTER_NUM = 1
            POINT_INDEX_CLUSTER(1) = NEIGHBOR(ITER)
            USED_POINT(POINT_INDEX_CLUSTER(1)) = .TRUE.
            
            DO I=1,RING_NUM
                BEFORE_POINT_CLUSTER_NUM = POINT_CLUSTER_NUM
                
                DO J=1,BEFORE_POINT_CLUSTER_NUM
                    IF(POINT_TYPE(POINT_INDEX_CLUSTER(J))==1) THEN
                    
                    DO K=1,CONNECTION_NUM(POINT_INDEX_CLUSTER(J))
                        CALL POINT_NEIGHBOR_POINT(FACE, CONNECTION, POINT_INDEX_CLUSTER(J), K, NEIGHBOR_POINT)
                        IF(.NOT. USED_POINT(NEIGHBOR_POINT)) THEN
                            POINT_CLUSTER_NUM = POINT_CLUSTER_NUM + 1
                            POINT_INDEX_CLUSTER(POINT_CLUSTER_NUM) = NEIGHBOR_POINT
                            USED_POINT(NEIGHBOR_POINT) = .TRUE.
                        END IF
                    END DO
                    
                    END IF
                END DO
            END DO
            
            ALLOCATE(POINT_CLUSTER(3,POINT_CLUSTER_NUM))
            
            POINT_CLUSTER(:,:) = NEWPOINT(:,POINT_INDEX_CLUSTER(1:POINT_CLUSTER_NUM))
            
            DO I=1,POINT_CLUSTER_NUM
                POINT_CLUSTER(:,I) = POINT_CLUSTER(:,I) - NEWPOINT(:,NEIGHBOR(ITER))
            END DO
            NEW_V_TEMP = NEW_V - NEWPOINT(:,NEIGHBOR(ITER))
            
            CALL LOCAL_POLY_FITTING(POLY_ORDER, POINT_CLUSTER_NUM, POINT_CLUSTER, NEW_V_TEMP, FITTING)
            
            FITTING_POINT(1:2) = NEW_V_TEMP(1:2)
            FITTING_POINT(3) = FITTING
            
            FITTING_POINT(:) = FITTING_POINT(:) + NEWPOINT(:,NEIGHBOR(ITER))
            
            FITTING_RESULT = FITTING_RESULT + WEIGHT(ITER) * FITTING_POINT
            
            DEALLOCATE(POINT_CLUSTER)
        END DO
        
        V_RESULT(:) = NEW_SCALE*(FITTING_RESULT(1) * NEW_AXIS(:,1) + FITTING_RESULT(2) * NEW_AXIS(:,2) + FITTING_RESULT(3) * NEW_AXIS(:,3))
        
        DEALLOCATE(NEWPOINT)
        
        DEALLOCATE(POINT_INDEX_CLUSTER, USED_POINT)
        
    END SUBROUTINE
    
    
    SUBROUTINE LOCAL_POLY_FITTING(D,M, POINT, V, FITTING)
        IMPLICIT NONE
        INTEGER :: D
        INTEGER :: M
        REAL(8) :: POINT(3,M), V(3), FITTING
        
        REAL(8), ALLOCATABLE :: VANDE(:,:), F(:), WEIGHT(:,:), FACTOR(:), A(:,:), B(:), C(:,:), X(:)
        REAL(8) :: MESH_RESOL
        INTEGER :: I, J, K, L, N, NUM
        
        N = (D+1)*(D+2)/2
        
        ALLOCATE(VANDE(M-1, N-1), F(M-1), WEIGHT(M-1,M-1), FACTOR(D+1))
        
        WEIGHT(:,:) = 0.
        
        FACTOR(1) = 1
        DO I = 1,D
            FACTOR(I+1) = FACTOR(I) * I
        END DO
        
        DO I = 1,M-1
            F(I) = POINT(3,I+1)
            
            NUM = 0
            DO L = 1, D
                DO J = 0, L
                    K = L - J
                    
                    NUM = NUM+1
                    VANDE(I,NUM) = POINT(1,I+1)**(REAL(J)) * POINT(2,I+1)**(REAL(K)) / (FACTOR(J+1) * FACTOR(K+1))
                END DO
            END DO
        END DO
        
        MESH_RESOL = 1.
        
        DO I = 1,M-1
            WEIGHT(I,I) = EXP( -DOT_PRODUCT(POINT(:,I+1),POINT(:,I+1))/MESH_RESOL)
            !WEIGHT(I,I) = 1.
        END DO
        
        ALLOCATE(A(N-1,N-1), B(N-1), C(M-1,N-1), X(N-1))
        
        C = MATMUL(WEIGHT,VANDE)
        
        A = MATMUL(TRANSPOSE(C),C)
        B = MATMUL(TRANSPOSE(C),MATMUL(WEIGHT,F))
        
        CALL CG_SOLVER(A,B, N-1, X)
        !CALL PCG_SOLVER(A,B, N-1, X)
        
        FITTING = 0.
        NUM = 0
        DO L = 1, D
            DO J = 0, L
                K = L - J
                NUM = NUM + 1
                FITTING = FITTING + X(NUM) * V(1)**(REAL(J)) * V(2)**(REAL(K)) / (FACTOR(J+1) * FACTOR(K+1))
            END DO
        END DO
        
        DEALLOCATE(A, B, C, X)
        
        DEALLOCATE(VANDE, F, WEIGHT, FACTOR)
        
    END SUBROUTINE LOCAL_POLY_FITTING
    
    
    !! High order methods
    
    SUBROUTINE PROJECTION_RIDGE_POINT_HIGHORDER(V, FACE_INDEX, EDGE_INDEX, POINT_NUM, POINT, FACE_NUM, FACE, CONNECTION_NUM, CONNECTION, RIDGE_FLAG, ALL_RIDGE, ALL_RIDGE_NUM, DIVIDED_BOUNDARY)
        IMPLICIT NONE
        REAL(8) :: V(3)
        INTEGER :: FACE_INDEX
        INTEGER :: EDGE_INDEX
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(:,:)
        INTEGER :: FACE_NUM
        INTEGER :: FACE(:,:)
        INTEGER :: CONNECTION_NUM(:)
        INTEGER :: CONNECTION(:,:)
        INTEGER :: RIDGE_FLAG
        INTEGER, OPTIONAL :: ALL_RIDGE(:,:)
        INTEGER, OPTIONAL :: ALL_RIDGE_NUM(:)
        INTEGER, OPTIONAL :: DIVIDED_BOUNDARY(:,:)
        
        REAL(8) :: V_RESULT(3)
        
        INTEGER, ALLOCATABLE :: RIDGE(:)
        REAL(8), ALLOCATABLE :: RIDGEPOINTS(:,:)
        INTEGER :: RIDGE_NUM
        
        INTEGER :: I0, I, I1, I2
        REAL(8) :: TEMPCOORD
        
        I1 = FACE(EDGE_INDEX, FACE_INDEX)
        I2 = FACE(MOD(EDGE_INDEX,3)+1, FACE_INDEX)
        
        CALL PROJECTION_EDGE_POINT(V, POINT(:,I1), POINT(:,I2), TEMPCOORD)
        
        IF(PRESENT(ALL_RIDGE)) THEN
            RIDGE_NUM = ALL_RIDGE_NUM(RIDGE_FLAG)
            ALLOCATE(RIDGE(RIDGE_NUM))
            RIDGE(:) = ALL_RIDGE(RIDGE_FLAG, 1:RIDGE_NUM)
        ELSE
            RIDGE_NUM = 0
            ALLOCATE(RIDGE(POINT_NUM))
            CALL FIND_RIDGE_ARRAY(POINT_NUM, FACE_NUM, FACE, CONNECTION_NUM, CONNECTION, DIVIDED_BOUNDARY, RIDGE_FLAG, RIDGE, RIDGE_NUM)
        END IF
        
        IF(RIDGE_NUM > 1) THEN
        
        DO I=1,RIDGE_NUM-1
            IF((RIDGE(I)==I1 .AND. RIDGE(I+1)==I2) .OR. (RIDGE(I)==I2 .AND. RIDGE(I+1)==I1)) THEN
                I0 = I
            END IF
        END DO
        
        ALLOCATE(RIDGEPOINTS(3,RIDGE_NUM))
        
        DO I=1,RIDGE_NUM
            RIDGEPOINTS(:,I) = POINT(:,RIDGE(I))
        END DO
        
        CALL RECONSTRUCTING_HIGH_ORDER_BOUNDARY(V, 3, RIDGE_NUM, RIDGEPOINTS, I0, V_RESULT)
        
        V(:) = V_RESULT(:)
        
        DEALLOCATE(RIDGEPOINTS)
        
        END IF
        
        DEALLOCATE(RIDGE)
        
    END SUBROUTINE
    
    SUBROUTINE PROJECTION_FACE_POINT_HIGHORDER(V, FACE_INDEX, POINT_NUM, POINT, FACE, CONNECTION_NUM, CONNECTION, POINT_TYPE)
        IMPLICIT NONE
        REAL(8) :: V(3)
        INTEGER :: FACE_INDEX
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(:,:)
        INTEGER :: FACE(:,:)
        INTEGER :: CONNECTION_NUM(:)
        INTEGER :: CONNECTION(:,:)
        INTEGER :: POINT_TYPE(:)
        
        REAL(8) :: V_RESULT(3)
        
        CALL PROJECTION_FACE_POINT(V, POINT(:,FACE(1,FACE_INDEX)), POINT(:,FACE(2,FACE_INDEX)), POINT(:,FACE(3,FACE_INDEX)))
        CALL RECONSTRUCTING_HIGH_ORDER_SURFACE(POINT_NUM, POINT, FACE, FACE_INDEX, V, 3, 2, CONNECTION_NUM, CONNECTION, POINT_TYPE, V_RESULT)
        V(:) = V_RESULT(:)
        
    END SUBROUTINE
    
    SUBROUTINE DISTANCE_FACE_POINT_HIGHORDER(V,FACE_INDEX, POINT_NUM, POINT, FACE, CONNECTION_NUM, CONNECTION, POINT_TYPE, R)
        IMPLICIT NONE
        REAL(8) :: V(3)
        INTEGER :: FACE_INDEX
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(:,:)
        INTEGER :: FACE(:,:)
        INTEGER :: CONNECTION_NUM(:)
        INTEGER :: CONNECTION(:,:)
        INTEGER :: POINT_TYPE(:)
        REAL(8) :: R
        
        REAL(8) :: V1(3), V2(3)
        
        V1(:) = V(:)
        
        CALL PROJECTION_FACE_POINT(V1, POINT(:,FACE(1,FACE_INDEX)), POINT(:,FACE(2,FACE_INDEX)), POINT(:,FACE(3,FACE_INDEX)))
        CALL RECONSTRUCTING_HIGH_ORDER_SURFACE(POINT_NUM, POINT, FACE, FACE_INDEX, V1, 3, 2, CONNECTION_NUM, CONNECTION, POINT_TYPE, V2)
        R = SQRT(DOT_PRODUCT(V2-V,V2-V))
    END SUBROUTINE
    
    SUBROUTINE SIGNED_DISTANCE_FACE_POINT_HIGHORDER(V, L, FACE_INDEX, POINT_NUM, POINT, FACE, CONNECTION_NUM, CONNECTION, POINT_TYPE, R)
        IMPLICIT NONE
        REAL(8) :: V(3), L(3)
        INTEGER :: FACE_INDEX
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(:,:)
        INTEGER :: FACE(:,:)
        INTEGER :: CONNECTION_NUM(:)
        INTEGER :: CONNECTION(:,:)
        INTEGER :: POINT_TYPE(:)
        REAL(8) :: R
        
        REAL(8) :: V1(3), V2(3)
        
        V1(:) = V(:)
        
        CALL PROJECTION_FACE_POINT(V1, POINT(:,FACE(1,FACE_INDEX)), POINT(:,FACE(2,FACE_INDEX)), POINT(:,FACE(3,FACE_INDEX)))
        CALL RECONSTRUCTING_HIGH_ORDER_SURFACE(POINT_NUM, POINT, FACE, FACE_INDEX, V1, 3, 2, CONNECTION_NUM, CONNECTION, POINT_TYPE, V2)
        R = SQRT(DOT_PRODUCT(V2-V,V2-V))
        
        IF(DOT_PRODUCT(L,V2-V)<0) THEN
            R = -R
        END IF
    END SUBROUTINE
    
    SUBROUTINE NEWPOINT_HIGHORDER(FACE_INDEX, EDGE_INDEX, POINT_NUM, POINT, FACE_NUM, FACE, CONNECTION_NUM, CONNECTION, POINT_TYPE, DIVIDED_BOUNDARY, V)
        IMPLICIT NONE
        INTEGER :: FACE_INDEX
        INTEGER :: EDGE_INDEX
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(:,:)
        INTEGER :: FACE_NUM
        INTEGER :: FACE(:,:)
        INTEGER :: CONNECTION_NUM(:)
        INTEGER :: CONNECTION(:,:)
        INTEGER :: POINT_TYPE(:)
        INTEGER, OPTIONAL :: DIVIDED_BOUNDARY(:,:)
        
        INTEGER :: RIDGE_FLAG
        REAL(8) :: V(3)
        
        V = (POINT(:,FACE(EDGE_INDEX,FACE_INDEX)) + POINT(:,FACE(MOD(EDGE_INDEX,3)+1,FACE_INDEX)))/2.
        
        IF(DIVIDED_BOUNDARY(EDGE_INDEX,FACE_INDEX)==0) THEN
            CALL PROJECTION_FACE_POINT_HIGHORDER(V, FACE_INDEX, POINT_NUM, POINT, FACE, CONNECTION_NUM, CONNECTION, POINT_TYPE)
        ELSE
            RIDGE_FLAG = DIVIDED_BOUNDARY(EDGE_INDEX,FACE_INDEX)
            CALL PROJECTION_RIDGE_POINT_HIGHORDER(V, FACE_INDEX, EDGE_INDEX, POINT_NUM, POINT, FACE_NUM, FACE, CONNECTION_NUM, CONNECTION, RIDGE_FLAG, DIVIDED_BOUNDARY = DIVIDED_BOUNDARY)
        END IF
        
    END SUBROUTINE

!    SUBROUTINE HIGHORDER_BOUNDARY_SMOOTHING(TEMP_POINT_NUM, TEMPPOINT, TEMP_FACE_NUM, TEMPFACE, TEMP_CONNECTION_NUM, TEMP_CONNECTION, TEMP_NEWEDGELENGTH, DIVIDED_BOUNDARY_ARRAY, BOUNDARY_FLAG)
!    IMPLICIT NONE
!    INTEGER :: TEMP_POINT_NUM
!    REAL(8) :: TEMPPOINT(:,:)
!    INTEGER :: TEMP_FACE_NUM
!    INTEGER :: TEMPFACE(:,:)
!    INTEGER :: TEMP_CONNECTION(:,:)
!    INTEGER :: TEMP_CONNECTION_NUM(:)
!    REAL(8) :: TEMP_NEWEDGELENGTH(:,:)
!    INTEGER :: DIVIDED_BOUNDARY_ARRAY(:,:)
!    INTEGER :: BOUNDARY_FLAG
!    REAL(8), ALLOCATABLE :: DISPLACEMENT(:,:)
!
!    ALLOCATE(DISPLACEMENT(3,TEMP_POINT_NUM))
!    DISPLACEMENT(:,:) = 0.
!    
!    
!    
!    END SUBROUTINE
    
END MODULE
