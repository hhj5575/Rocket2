MODULE SURFACE_MODULE_2D
    
    IMPLICIT NONE
    
    REAL(8), PARAMETER :: MINERROR = 1E-14, MINERROR2 = 1E-11, PI = 3.141592653589793238, MINLOCERROR = 1E-2, MINNORMERROR = 5E-2, LARGE_NUMBER = 100000.
    
    REAL(8) :: DOMAIN_MAX(2), DOMAIN_MIN(2)

    INTEGER :: SURFACE_AREA_ITER
    INTEGER :: SURFACE_PRESSURE_ITER	
    REAL(8) :: SURFACE_TOTAL_TIME
    REAL(8) :: SURFACE_TIME_STEP
    REAL(8) :: SURFACE_AREA_ARRAY(100000)
    REAL(8) :: TOTAL_PRESSURE_ARRAY(100000,3)
    INTEGER :: SURFACE_FLAG_ARRAY(100000,10)
    INTEGER :: RESTART_FLAG
    INTEGER :: RESTART_ITERATION
    REAL(8) :: CHI_C_LOW, CHI_C_HIGH
    REAL(8) :: INTERFACE_THRESHOLD
    INTEGER :: SMALL_REGION_POINT_NUM
    REAL(8) :: THIN_REGION_EXISTENCE
    REAL(8) :: THIN_REGION_ATTACHMENT
    REAL(8) :: EDGE_SPLITTING
    REAL(8) :: EDGE_COLLAPSING

    TYPE SURFACE_TYPE
        
        REAL(8) :: MESH_SIZE
        REAL(8) :: MESH_SIZE_MAX
        REAL(8), ALLOCATABLE :: SURFACE_POINTS(:,:)
        INTEGER :: SURFACE_POINTS_NUM
        INTEGER, ALLOCATABLE :: SURFACE_EDGES(:,:)
        INTEGER :: SURFACE_EDGES_NUM
	INTEGER :: SURFACE_PATCHES_NUM
	INTEGER, ALLOCATABLE :: SURFACE_PATCHES_TOPCHANGE_TYP(:)
        
        REAL(8), ALLOCATABLE :: SURFACE_INITIAL_EDGE_LENGTH(:)
        
        REAL(8), ALLOCATABLE :: EDGE_B_RATE(:)
        REAL(8), ALLOCATABLE :: POINT_VELOCITY(:,:)
        REAL(8), ALLOCATABLE :: POINT_DISPLACEMENT(:,:)
        INTEGER, ALLOCATABLE :: POINT_EDGE_CONNECTION(:,:)
        INTEGER, ALLOCATABLE :: POINT_TYPE(:)
        
        INTEGER, ALLOCATABLE :: EDGE_LOCATION(:)
        INTEGER, ALLOCATABLE :: EDGE_ONINTERFACE(:)
        
        INTEGER, ALLOCATABLE :: POINT_RELATEDPT(:,:)
        INTEGER, ALLOCATABLE :: POINT_RELATEDEDGE(:,:)
        
        REAL(8), ALLOCATABLE :: EDGE_PRESSURE(:)
        REAL(8), ALLOCATABLE :: POINT_FORCE(:,:)
        
        INTEGER, ALLOCATABLE :: EDGE_IMPACT_ZONE(:,:)
	INTEGER, ALLOCATABLE :: EDGE_ABLATION_FLAG(:)
        
    END TYPE SURFACE_TYPE
    
    TYPE(SURFACE_TYPE), TARGET :: SURFACE_FLUID
    TYPE(SURFACE_TYPE), TARGET :: SURFACE_PROPEL
    TYPE(SURFACE_TYPE), TARGET :: SURFACE_CASE
    
    INTEGER :: INTERFACE_FLUID_POINTS_NUM
    REAL(8), ALLOCATABLE, TARGET :: INTERFACE_FLUID_POINTS(:,:)
    INTEGER, ALLOCATABLE, TARGET :: INTERFACE_FLUID_POINTS_LOC(:,:)
    
    INTEGER :: INTERFACE_STRUCT_POINTS_NUM
    REAL(8), ALLOCATABLE, TARGET :: INTERFACE_STRUCT_POINTS(:,:)
    INTEGER, ALLOCATABLE, TARGET :: INTERFACE_STRUCT_POINTS_LOC(:,:)
    
    
    CONTAINS
    
    SUBROUTINE DISTANCE_LINE_POINT(V,W1,W2,     D) 
        REAL(8) :: V(2),W1(2),W2(2),VEC(2)
        REAL(8) :: L2,INNER,D
        INTEGER :: I
        
        L2 = DOT_PRODUCT(W2-W1,W2-W1)
        INNER = DOT_PRODUCT(W1-V,W2-W1)
        !$OMP PARALLEL DO PRIVATE(I)
        DO I=1,2
           VEC(I) = W1(I)- V(I) - INNER/L2 * (W2(I)-W1(I))
        END DO
        !$OMP END PARALLEL DO
        
        D = SQRT(DOT_PRODUCT(VEC,VEC))
    
    END SUBROUTINE DISTANCE_LINE_POINT
    
    SUBROUTINE MINMOD(A,B,      T)
        REAL(8) :: A,B,T
        INTEGER :: S1,S2
    
        CALL SIGN1(A,S1)
        CALL SIGN1(B,S2)
        T = (S1+S2)/2 * MIN(ABS(A),ABS(B))
        
    END SUBROUTINE MINMOD
    
    SUBROUTINE SIGN1(A,       T)
        REAL(8) :: A 
        INTEGER :: T
    
        IF(A.LT.0.0) THEN 
            T=-1
            
        ELSE IF(A.EQ.0.0) THEN
            T=0
  
        ELSE 
            T=1
        END IF
  
    END SUBROUTINE SIGN1
    
    SUBROUTINE SIGN2(R,DX,       T) ! X
        REAL(8) :: R,DX
        REAL(8) :: T
        
        T = R/SQRT(R**2 + DX**2)
  
    END SUBROUTINE SIGN2
    
    SUBROUTINE VEC_CURL1(V,W,       R) ! DIM=2
        REAL(8) :: V(2),W(2),R
    
        R = V(1)*W(2) - V(2)*W(1)
    
    END SUBROUTINE VEC_CURL1
    
    SUBROUTINE VEC_CURL2(V1,V2,W1,W2,       R) ! DIM=2 
        REAL(8) :: V1(2),V2(2),W1(2),W2(2),R
    
       
        R = (V2(1)-V1(1))*(W2(2)-W1(2)) - (V2(2)-V1(2))*(W2(1)-W1(1))
    
    END SUBROUTINE VEC_CURL2
    
    SUBROUTINE INIT_RANDOM_SEED()

        IMPLICIT NONE
        
        INTEGER :: I, N, CLOCK
        INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
        CALL RANDOM_SEED(SIZE = N)
        ALLOCATE(SEED(N))
        CALL SYSTEM_CLOCK(COUNT=CLOCK)
        SEED = CLOCK + 37*(/(I-1, I=1, N)/)
        CALL RANDOM_SEED(PUT = SEED)
        DEALLOCATE(SEED)

    END SUBROUTINE INIT_RANDOM_SEED
    
    !TYPE AVL_TREE
    !    INTEGER :: KEYTYPE1
    !    REAL(8) :: KEYTYPE2
    !    INTEGER :: KEYTYPE3(2)
    !    
    !    INTEGER :: BALANCE
    !END TYPE AVL_TREE
    !
    !SUBROUTINE AVL_CW_ROTATE(T, NODE_INDEX, CHILD_INDEX
    !    TYPE 
    
    !! SOLVERS    
        
    SUBROUTINE PCG_SOLVER(A, B, N, X)
	IMPLICIT NONE
        REAL(8), DIMENSION(:,:) :: A
        REAL(8), DIMENSION(:) :: B
        REAL(8), DIMENSION(:) :: X
        REAL(8), ALLOCATABLE :: D(:)
        REAL(8), ALLOCATABLE :: Z(:)
        REAL(8), ALLOCATABLE :: R(:)
        REAL(8), ALLOCATABLE :: OLDZ(:)
        REAL(8), ALLOCATABLE :: OLDR(:)
        REAL(8), ALLOCATABLE :: CHOL(:,:)
        INTEGER :: N
        INTEGER :: MAX_ITER
        REAL(8) :: DELTA, BETA, ALPHA
        INTEGER :: K
        REAL(8), ALLOCATABLE :: tempR(:)
        
        ALLOCATE(R(N), OLDR(N), D(N), Z(N), OLDZ(N), CHOL(N,N), tempr(n))
        MAX_ITER = 1000
        DELTA = 10.0**(-10.0)
        X = 0.0
        R = B - MATMUL(A,X)
!write(*,*) 'before CHOLESKY'
        CALL CHOLESKY(N,A,CHOL)
!write(*,*) 'after CHOLESKY'
        DO K = 1, MAX_ITER
!write(*,*) 'before LU',k
            CALL LU(N,R,CHOL,TRANSPOSE(CHOL), Z)
!write(*,*) 'after LU',k
            tempr = matmul( matmul(chol, TRANSPOSE(CHOL)), z)
            IF (K .EQ. 1) THEN
                D = Z
            ELSE
!write(*,*) 'oldr,oldz' , DOT_PRODUCT(OLDR,OLDZ)
                BETA = DOT_PRODUCT(R,Z)/DOT_PRODUCT(OLDR,OLDZ)
                D = Z + BETA*D
            END IF
!write(*,*) 'DOT_PRODUCT(D,MATMUL(A,D))' , DOT_PRODUCT(D,MATMUL(A,D))
!write(*,*) 'DOT_PRODUCT(D,MATMUL(A,D))' , D
!write(*,*) 'DOT_PRODUCT(D,MATMUL(A,D))' , A
            ALPHA = DOT_PRODUCT(R,Z)/DOT_PRODUCT(D,MATMUL(A,D))
            X = X + ALPHA*D
            OLDR = R
            OLDZ = Z
            R = R - ALPHA* MATMUL(A,D)
            IF (DOT_PRODUCT(R,R)< DELTA) THEN
                RETURN
            END IF
        END DO
        
	DEALLOCATE(R, OLDR, D, Z, OLDZ, CHOL)
    END SUBROUTINE PCG_SOLVER
    
    SUBROUTINE CHOLESKY(N,A,R)
	IMPLICIT NONE
        INTEGER :: N
        REAL(8), DIMENSION(:,:) :: A
        REAL(8), DIMENSION(:,:) :: R
        INTEGER :: I,K,S
        REAL(8) :: SUM
        
        R = 0.0
        
        DO K = 1, N
            SUM = A(K,K)
            DO S = 1, K-1
                SUM = SUM - R(K,S)*R(K,S)
            END DO
            R(K,K) = SQRT(max(0.0,SUM))
            
            DO I = K+1,N
                IF (A(I,K) .NE. 0.0) THEN
                    SUM = A(I,K)
                    DO S = 1, K-1
                        SUM = SUM - R(I,S)*R(K,S)
                    END DO
                    R(I,K) = SUM/R(K,K)
                END IF
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
        Z = 0
        Z(1) = B(1)/L(1,1)
        
        DO I = 2, N
            SUM = B(I)
            DO J = 1, I-1
                SUM = SUM - L(I,J)*Z(J)
            END DO
!write(*,*) 'L(I,I)',i, L(I,I)
if(abs(L(I,I))< minerror) then
            Z(I) = SUM/(L(I,I)+minerror )
else
            Z(I) = SUM/L(I,I)
end if
        END DO
        X = 0
!write(*,*) 'u(I,I)',i, u(I,I)

if(abs(U(N,N))< minerror) then
            X(N) = Z(N)/(U(N,N)+minerror )
else
            X(N) = Z(N)/U(N,N)
end if
        !X(N) = Z(N)/U(N,N)
        DO I = N-1,1,-1
            SUM = Z(I)
            DO J= I+1,N
                SUM = SUM - U(I,J)*X(J)
            END DO
if(abs(U(N,N))< minerror) then
            X(I) = SUM/(U(I,I)+minerror )
else
            X(I) = SUM/U(I,I)
end if
            !X(I) = SUM/U(I,I)
        END DO
        DEALLOCATE(Z)
    END SUBROUTINE LU
    
END MODULE
