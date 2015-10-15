MODULE SURFACES_3D

    USE SURFACE_MODULE_3D
    
    CONTAINS

    SUBROUTINE UNSIGNED_DISTANCE_EDGE_POINT(V,W1,W2,DIST)
    
        REAL(8) :: V(3),W1(3),W2(3)
    
        REAL(8) :: DIST1, DIST2, POINTDIST, DIST
        REAL(8) :: A1(3),A2(3)
        REAL(8) :: INNER
        
        DIST1 = SQRT(DOT_PRODUCT(V-W1,V-W1))
        DIST2 = SQRT(DOT_PRODUCT(V-W2,V-W2))
        
        IF(SQRT(DOT_PRODUCT(W1-W2,W1-W2))<MINERROR) THEN
            DIST = DIST1
            RETURN
        END IF
            
        IF(DIST1 < DIST2) THEN
            A1 = W1
            A2 = W2
            
            POINTDIST = DIST1
        ELSE
            A1 = W2
            A2 = W1
            
            POINTDIST = DIST2
        END IF
        
        INNER = DOT_PRODUCT(A1-V,A2-A1)
    
        IF(INNER>0) THEN
            DIST = POINTDIST
        ELSE
            CALL DISTANCE_LINE_POINT(V,A1,A2,DIST)
        END IF
    
    END SUBROUTINE UNSIGNED_DISTANCE_EDGE_POINT

    SUBROUTINE PROJECTION_EDGE_POINT(V,W1,W2,COORD)
    
        REAL(8) :: V(3),W1(3),W2(3),L(3), R
    
        REAL(8) :: DIST1, DIST2, POINTDIST
        REAL(8) :: A1(3),A2(3)
        REAL(8) :: INNER, INNER2
	REAL(8) :: COORD

	LOGICAL :: FLAG
        
        DIST1 = SQRT(DOT_PRODUCT(V-W1,V-W1))
        DIST2 = SQRT(DOT_PRODUCT(V-W2,V-W2))
            
        IF(DIST1 < DIST2) THEN
            A1 = W1
            A2 = W2
            
            POINTDIST = DIST1
	    FLAG = .FALSE.
        ELSE
            A1 = W2
            A2 = W1
            
            POINTDIST = DIST2
	    FLAG = .TRUE.
        END IF
        
        INNER = DOT_PRODUCT(A1-V,A2-A1)
        
        L = A2-A1
        R = SQRT(DOT_PRODUCT(L,L))
        L = L/R
        INNER2 = DOT_PRODUCT(V-A1,L)
        
        V = A1 + INNER2 * L

	IF(FLAG) THEN
	    COORD = 1-INNER2/R
	ELSE
	    COORD = INNER2/R
	END IF
        
        IF(COORD<0.) THEN
            COORD = 0.
            V = W1
        ELSE IF(COORD>1.) THEN
            COORD = 1.
            V = W2
        END IF
    
    END SUBROUTINE PROJECTION_EDGE_POINT
    
    SUBROUTINE UNSIGNED_DISTANCE_FACE_POINT(V,W1,W2,W3,DIST)
        
        REAL(8) :: V(3) 
    
        INTEGER :: IMIN
        REAL(8) :: DIST, POINTDIST, DIST1, DIST2, DIST3
        REAL(8) :: A1(3),A2(3),A3(3), W1(3),W2(3),W3(3)
        REAL(8) :: N(3), C1(3), C2(3), C3(3)
        REAL(8) :: ROT1, ROT2, ROT3, INNER1, INNER2
        
        POINTDIST = SQRT(DOT_PRODUCT(V-W1, V-W1))
        IMIN = 1
        IF(POINTDIST > SQRT(DOT_PRODUCT(V-W2, V-W2))) THEN
            POINTDIST = SQRT(DOT_PRODUCT(V-W2, V-W2))
            IMIN = 2
        END IF
        IF(POINTDIST > SQRT(DOT_PRODUCT(V-W3, V-W3))) THEN
            POINTDIST = SQRT(DOT_PRODUCT(V-W3, V-W3))
            IMIN = 3
        END IF
            
        IF(IMIN == 1) THEN
            A1 = W1
            A2 = W2
            A3 = W3
        ELSE IF(IMIN == 2) THEN
            A1 = W2
            A2 = W3
            A3 = W1
        ELSE
            A1 = W3
            A2 = W1
            A3 = W2
        END IF
        
        CALL VEC_CURL2(A1,A2,A1,A3,N)
    
        CALL VEC_CURL2(V,A1,A2,A1,C1)
        CALL VEC_CURL2(V,A2,A3,A2,C2)
        CALL VEC_CURL2(V,A3,A1,A3,C3)
        
        ROT1 = DOT_PRODUCT(C1,N)
        ROT2 = DOT_PRODUCT(C2,N)
        ROT3 = DOT_PRODUCT(C3,N)
        
        INNER1 = DOT_PRODUCT(A1-V, A2-A1)
        INNER2 = DOT_PRODUCT(A1-V, A3-A1)
        
        IF(INNER1>0 .AND. INNER2>0) THEN
            DIST = POINTDIST
        ELSE IF(ROT1<0 .AND. ROT2<0 .AND. ROT3<0) THEN
            CALL DISTANCE_FACE_POINT(V,A1,A2,A3,DIST)
        ELSE
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A1,A2,DIST1)
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A1,A3,DIST2)
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A2,A3,DIST3)
	    DIST = MIN(DIST1,DIST2,DIST3)
        END IF
    
    END SUBROUTINE UNSIGNED_DISTANCE_FACE_POINT
    
    SUBROUTINE PROJECTION_FACE_POINT(V,V1,V2,V3)
        
        REAL(8) :: V(3)
        REAL(8) :: POINTDIST
        INTEGER :: IMIN, I
        
        REAL(8) :: A1(3),A2(3),A3(3),V1(3),V2(3),V3(3)
        REAL(8) :: N(3), C1(3), C2(3), C3(3)
        REAL(8) :: ROT1, ROT2, ROT3, INNER1, INNER2
        
        REAL(8) :: PTONSURF(3)
        REAL(8) :: L1(3),L2(3),L3(3)
        REAL(8) :: L1_NORM, L2_NORM, L3_NORM
        REAL(8) :: W1,W2,W3
        REAL(8) :: D1, D2
	REAL(8) :: DIST, TEMPDIST
        
        POINTDIST = SQRT(DOT_PRODUCT(V-V1, V-V1))
        IMIN = 1
        IF(POINTDIST > SQRT(DOT_PRODUCT(V-V2, V-V2))) THEN
            POINTDIST = SQRT(DOT_PRODUCT(V-V2, V-V2))
            IMIN = 2
        END IF
        IF(POINTDIST > SQRT(DOT_PRODUCT(V-V3, V-V3))) THEN
            POINTDIST = SQRT(DOT_PRODUCT(V-V3, V-V3))
            IMIN = 3
        END IF
            
        IF(IMIN == 1) THEN
            A1 = V1
            A2 = V2
            A3 = V3
        ELSE IF(IMIN == 2) THEN
            A1 = V2
            A2 = V3
            A3 = V1
        ELSE
            A1 = V3
            A2 = V1
            A3 = V2
        END IF
        
        CALL VEC_CURL2(A1,A2,A1,A3,N)
    
        CALL VEC_CURL2(V,A1,A2,A1,C1)
        CALL VEC_CURL2(V,A2,A3,A2,C2)
        CALL VEC_CURL2(V,A3,A1,A3,C3)
        
        ROT1 = DOT_PRODUCT(C1,N)
        ROT2 = DOT_PRODUCT(C2,N)
        ROT3 = DOT_PRODUCT(C3,N)
        
        INNER1 = DOT_PRODUCT(A1-V, A2-A1)
        INNER2 = DOT_PRODUCT(A1-V, A3-A1)
        
        IF(INNER1>0 .AND. INNER2>0) THEN
            V = A1
        ELSE IF(ROT1<0 .AND. ROT2<0 .AND. ROT3<0) THEN
            N = N/SQRT(DOT_PRODUCT(N,N))
            
            D1 = DOT_PRODUCT(V-A1,N)
            
            PTONSURF = V - D1*N
            
            CALL VEC_CURL2(PTONSURF,A2,PTONSURF,A3,  L1)
            CALL VEC_CURL2(PTONSURF,A3,PTONSURF,A1,  L2)
            CALL VEC_CURL2(PTONSURF,A1,PTONSURF,A2,  L3)
            
            L1_NORM = SQRT(DOT_PRODUCT(L1,L1))
            L2_NORM = SQRT(DOT_PRODUCT(L2,L2))
            L3_NORM = SQRT(DOT_PRODUCT(L3,L3))
            
            W1 = L1_NORM/(L1_NORM+L2_NORM+L3_NORM)
            W2 = L2_NORM/(L1_NORM+L2_NORM+L3_NORM)
            W3 = L3_NORM/(L1_NORM+L2_NORM+L3_NORM)
            
            V = W1*A1 + W2*A2 + W3*A3
        ELSE
	    DIST = MAX(DOMAIN_MAX(1)-DOMAIN_MIN(1),DOMAIN_MAX(2)-DOMAIN_MIN(2),DOMAIN_MAX(3)-DOMAIN_MIN(3))
	    DO I = 1,3
		IF(I==1) THEN
		    CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A1,A2,TEMPDIST) 
		ELSE IF(I==2) THEN
            	    CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A1,A3,TEMPDIST)
		ELSE
            	    CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A2,A3,TEMPDIST)
		END IF

		IF(TEMPDIST<DIST) THEN
		    IMIN = I
		    DIST = TEMPDIST
		END IF
	    END DO

	    IF(IMIN==1) THEN

		D1 = DOT_PRODUCT(A2-A1,A2-V)
		D2 = DOT_PRODUCT(A2-A1,A2-A1)
		W1 = D1/D2
		D1 = DOT_PRODUCT(A1-A2,A1-V)
		W2 = D1/D2
		V = W1*A1 + W2*A2
	
	    ELSE IF(IMIN==2) THEN

		D1 = DOT_PRODUCT(A3-A1,A3-V)
		D2 = DOT_PRODUCT(A3-A1,A3-A1)
		W1 = D1/D2
		D1 = DOT_PRODUCT(A1-A3,A1-V)
		W3 = D1/D2
		V = W1*A1 + W3*A3

	    ELSE

		D1 = DOT_PRODUCT(A3-A2,A3-V)
		D2 = DOT_PRODUCT(A3-A2,A3-A2)
		W2 = D1/D2
		D1 = DOT_PRODUCT(A2-A3,A2-V)
		W3 = D1/D2
		V = W2*A2 + W3*A3

	    END IF
        END IF
    
    END SUBROUTINE PROJECTION_FACE_POINT
    
    
    
    
    SUBROUTINE LINE_FACE_INTERSECTING(POINT_NUM, POINT, FACE_NUM, FACE, I0, L, V, SGN, T)
    
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(3,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: FACE_NUM
        INTEGER :: FACE(3,FACE_NUM) ! FACE(3,FACE_NUM)
        INTEGER :: I0
        REAL(8) :: L(3), V(3) 
        REAL(8) :: A(3),B(3),C(3)
        REAL(8) :: M(3)
        REAL(8) :: L1(3), L2(3)
        REAL(8) :: A1,A2,B1,B2,C1,C2
        REAL(8) :: D1,D2,D3
        REAL(8) :: BC
        INTEGER :: SGN
        REAL(8) :: R,S,T
        
        REAL(8) :: TEMP
        INTEGER :: I
        
        DO I=1,3
            A(I) = POINT(I,FACE(1,I0))
            B(I) = POINT(I,FACE(2,I0))
            C(I) = POINT(I,FACE(3,I0))
        END DO
    
        CALL RANDOM_NUMBER(M)
        CALL VEC_CURL1(L,M,L1)
        TEMP = SQRT(DOT_PRODUCT(L1,L1))
        IF(TEMP<MINERROR) THEN
            SGN = -2
            RETURN
        END IF
        
        L1 = L1 / SQRT(DOT_PRODUCT(L1,L1))
        CALL VEC_CURL1(L,L1,L2)
        
        A1 = DOT_PRODUCT(V-A,L1)
        A2 = DOT_PRODUCT(V-A,L2)
        B1 = DOT_PRODUCT(B-A,L1)
        B2 = DOT_PRODUCT(B-A,L2)
        C1 = DOT_PRODUCT(C-A,L1)
        C2 = DOT_PRODUCT(C-A,L2)
        
        D1 = DOT_PRODUCT(V-A,L)
        D2 = DOT_PRODUCT(B-A,L)
        D3 = DOT_PRODUCT(C-A,L)
    
        BC = B1*C2 - B2*C1 
    
        IF ( ABS(BC) < MINERROR) THEN
            SGN = -2
            RETURN
        END IF
    
        S = (A1*C2 - A2*C1)/BC
        R = (A1*B2 - A2*B1)/(-BC)
        T = S*D2 + R*D3 - D1
        
        IF(S<0 .OR. R<0 .OR. S+R>1) THEN
            SGN = -3
        ELSE IF(T >= 0 .AND. S > 0.AND. R > 0 .AND. S+R < 1) THEN
            SGN = 1
        ELSE IF(T >= 0 .AND. (S == 0 .OR. R == 0 .OR. S + R == 1)) THEN
            SGN = -1
        ELSE 
            SGN = 0
        END IF
    
    END SUBROUTINE LINE_FACE_INTERSECTING
    
    
    
    
    
    
    SUBROUTINE LINE_SURFACE_INTERSECTING(POINT_NUM, POINT, FACE_NUM, FACE,V,RETURNS)

        INTEGER :: POINT_NUM
        REAL(8) :: POINT(3,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: FACE_NUM
        INTEGER :: FACE(3,FACE_NUM) ! FACE(3,FACE_NUM)
        INTEGER :: NUM
        INTEGER :: SGN
        REAL(8) :: L(3), V(3), T
        REAL(8) :: RETURNS
        INTEGER :: I
        
        SGN = -1
        
        DO WHILE(SGN == -1 .OR. SGN == -2)

            NUM = 0
            CALL RANDOM_NUMBER(L)
            L = L/SQRT(DOT_PRODUCT(L,L))
            
            DO I=1,FACE_NUM
                CALL LINE_FACE_INTERSECTING(POINT_NUM, POINT, FACE_NUM, FACE,I,L,V,SGN,T)
                
                IF (SGN == -1 .OR. SGN == -2) THEN
                    EXIT
                END IF
                
                IF (SGN == 1) THEN
                    NUM = NUM+1
                END IF
            END DO
        END DO
        
        IF(MOD(NUM,2) == 0 ) THEN
            RETURNS = 1
        ELSE
            RETURNS = -1
        END IF
        
    END SUBROUTINE LINE_SURFACE_INTERSECTING
    
    
    !! MODIFIED
    SUBROUTINE LINE_SURFACE_INTERSECTING_DISTANCE(POINT_NUM, POINT, FACE_NUM, FACE, L, V, TMIN, IMIN)

        INTEGER :: POINT_NUM
        REAL(8) :: POINT(3,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: FACE_NUM
        INTEGER :: FACE(3,FACE_NUM) ! FACE(3,FACE_NUM)
        
        INTEGER :: NUM
        INTEGER :: SGN
        
        REAL(8) :: L(3), L_TEMP(3), V(3), T, TMIN
        
        INTEGER :: I, IMIN

        NUM = 0
        L_TEMP = L/SQRT(DOT_PRODUCT(L,L))

        DO I=1,FACE_NUM
            CALL LINE_FACE_INTERSECTING(POINT_NUM, POINT, FACE_NUM, FACE,I,L_TEMP,V,SGN,T)

            IF (SGN.NE.-2 .AND. SGN.NE.-3 .AND. T > -SURFACE_FLUID%MESH_SIZE/10.) THEN
                NUM = NUM+1
                IF(NUM==1) THEN
                    TMIN = T
                    IMIN = I
                ELSE IF(T < TMIN) THEN
                    TMIN = T
                    IMIN = I
                END IF
            END IF
        END DO
        
        IF(NUM==0) THEN
            IMIN = 0
        END IF
        
    END SUBROUTINE LINE_SURFACE_INTERSECTING_DISTANCE
    !! END MODIFIED
    
    
    SUBROUTINE DISTANCE_SURFACE_POINT(POINT_NUM, POINT, FACE_NUM, FACE, V, DIST, IDX)

        INTEGER :: POINT_NUM
        REAL(8) :: POINT(3,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: FACE_NUM
        INTEGER :: FACE(3,FACE_NUM) ! FACE(3,FACE_NUM)
        REAL(8) :: DIST, TEMPDIST
        REAL(8) :: V(3)
        INTEGER :: I, IDX
        !REAL(8) :: TEMP
        
        DO I=1,FACE_NUM
            CALL UNSIGNED_DISTANCE_FACE_POINT(V, POINT(:,FACE(1,I)), POINT(:,FACE(2,I)), POINT(:,FACE(3,I)), TEMPDIST)
            
            IF (I==1) THEN
                DIST = TEMPDIST
		IDX = I
            ELSEIF(TEMPDIST <= DIST) THEN
                DIST = TEMPDIST
		IDX = I
            END IF
        END DO 
  
        !IF (DIST < MINERROR) THEN
        !    DIST = 0.
        !    RETURN
        !END IF
        !
        !CALL LINE_SURFACE_INTERSECTING(POINT_NUM, POINT, FACE_NUM, FACE,V,TEMP)
        !DIST = TEMP * DIST
        
    END SUBROUTINE DISTANCE_SURFACE_POINT
    
    SUBROUTINE INTERPOLATION_ONE_FACE(V, V1, V2, V3, DATA1, DATA2, DATA3,       INTER_DATA)
        
        REAL(8) :: V(3)
        
        REAL(8) :: V1(3),V2(3),V3(3)
        REAL(8) :: DATA1, DATA2, DATA3
        
        REAL(8) :: INTER_DATA
        
        INTEGER :: IMIN
        
        REAL(8) :: A1(3),A2(3),A3(3), TEMPDIST, POINTDIST, DIST1, DIST2, DIST3
        REAL(8) :: B1, B2, B3
        REAL(8) :: N(3), C1(3), C2(3), C3(3)
        REAL(8) :: ROT1, ROT2, ROT3, INNER1, INNER2
        
        REAL(8) :: PTONSURF(3)
        REAL(8) :: L1(3),L2(3),L3(3)
        REAL(8) :: L1_NORM, L2_NORM, L3_NORM
        REAL(8) :: W1,W2,W3
        REAL(8) :: D1, D2
        
        IMIN = 1
        POINTDIST = SQRT(DOT_PRODUCT(V-V1, V-V1))
        TEMPDIST = SQRT(DOT_PRODUCT(V-V2, V-V2))
        IF(TEMPDIST < POINTDIST) THEN
            POINTDIST = TEMPDIST
            IMIN = 2
        END IF
        TEMPDIST = SQRT(DOT_PRODUCT(V-V3, V-V3))
        IF(TEMPDIST < POINTDIST) THEN
            POINTDIST = TEMPDIST
            IMIN = 3
        END IF
        
        IF(IMIN == 1) THEN
            A1(:) = V1
            A2(:) = V2
            A3(:) = V3
            
            B1 = DATA1
            B2 = DATA2
            B3 = DATA3
        ELSEIF(IMIN == 2) THEN
            A1(:) = V2
            A2(:) = V3
            A3(:) = V1
            
            B1 = DATA2
            B2 = DATA3
            B3 = DATA1
        ELSE
            A1(:) = V3
            A2(:) = V1
            A3(:) = V2
            
            B1 = DATA3
            B2 = DATA1
            B3 = DATA2
        END IF
        
        CALL VEC_CURL2(A1,A2,A1,A3,N)
    
        CALL VEC_CURL2(V,A1,A2,A1,C1)
        CALL VEC_CURL2(V,A2,A3,A2,C2)
        CALL VEC_CURL2(V,A3,A1,A3,C3)
        
        ROT1 = DOT_PRODUCT(C1,N)
        ROT2 = DOT_PRODUCT(C2,N)
        ROT3 = DOT_PRODUCT(C3,N)
        
        INNER1 = DOT_PRODUCT(A1-V, A2-A1)
        INNER2 = DOT_PRODUCT(A1-V, A3-A1)
        
        IF(INNER1>0 .AND. INNER2>0) THEN
            INTER_DATA = B1
        ELSE IF(ROT1<0 .AND. ROT2<0 .AND. ROT3<0) THEN
            N = N/SQRT(DOT_PRODUCT(N,N))
            
            D1 = DOT_PRODUCT(V-A1,N)
            
            PTONSURF = V - D1*N
            
            CALL VEC_CURL2(PTONSURF,A2,PTONSURF,A3,  L1)
            CALL VEC_CURL2(PTONSURF,A3,PTONSURF,A1,  L2)
            CALL VEC_CURL2(PTONSURF,A1,PTONSURF,A2,  L3)
            
            L1_NORM = SQRT(DOT_PRODUCT(L1,L1))
            L2_NORM = SQRT(DOT_PRODUCT(L2,L2))
            L3_NORM = SQRT(DOT_PRODUCT(L3,L3))
            
            W1 = L1_NORM/(L1_NORM+L2_NORM+L3_NORM)
            W2 = L2_NORM/(L1_NORM+L2_NORM+L3_NORM)
            W3 = L3_NORM/(L1_NORM+L2_NORM+L3_NORM)
            
            INTER_DATA = W1*B1 + W2*B2 + W3*B3
        ELSE
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A1,A2,DIST1)
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A1,A3,DIST2)
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A2,A3,DIST3)
	    TEMPDIST = MIN(DIST1,DIST2,DIST3)
            IF(TEMPDIST==DIST1) THEN
                D1 = DOT_PRODUCT(A2-A1,A2-V)
                D2 = DOT_PRODUCT(A2-A1,A2-A1)
                W1 = D1/D2
                D1 = DOT_PRODUCT(A1-A2,A1-V)
                W2 = D1/D2
                
                INTER_DATA = W1*B1 + W2*B2
            ELSE IF(TEMPDIST==DIST2) THEN
                D1 = DOT_PRODUCT(A3-A1,A3-V)
                D2 = DOT_PRODUCT(A3-A1,A3-A1)
                W1 = D1/D2
                D1 = DOT_PRODUCT(A1-A3,A1-V)
                W3 = D1/D2
                
                INTER_DATA = W1*B1 + W3*B3
            ELSE
                D1 = DOT_PRODUCT(A3-A2,A3-V)
                D2 = DOT_PRODUCT(A3-A2,A3-A2)
                W2 = D1/D2
                D1 = DOT_PRODUCT(A2-A3,A2-V)
                W3 = D1/D2
                
                INTER_DATA = W2*B2 + W3*B3
            END IF
        END IF
        
    END SUBROUTINE INTERPOLATION_ONE_FACE
    
    
    
    SUBROUTINE COORDINATE_EDGE_POINT(V,V1,V2,WEIGHT1,WEIGHT2)
    
        REAL(8) :: V(3),V1(3),V2(3),L(3), R
    
        REAL(8) :: DIST1, DIST2, POINTDIST
        REAL(8) :: A1(3),A2(3)
        REAL(8) :: INNER, INNER2
	REAL(8) :: WEIGHT1,WEIGHT2

	LOGICAL :: FLAG
        
        DIST1 = SQRT(DOT_PRODUCT(V-V1,V-V1))
        DIST2 = SQRT(DOT_PRODUCT(V-V2,V-V2))
            
        IF(DIST1 < DIST2) THEN
            A1 = V1
            A2 = V2
            
            POINTDIST = DIST1
	    FLAG = .FALSE.
        ELSE
            A1 = V2
            A2 = V1
            
            POINTDIST = DIST2
	    FLAG = .TRUE.
        END IF
        
        INNER = DOT_PRODUCT(A1-V,A2-A1)
        
        L = A2-A1
        R = SQRT(DOT_PRODUCT(L,L))
        L = L/R
        INNER2 = DOT_PRODUCT(V-A1,L)

	IF(FLAG) THEN
            WEIGHT1 = INNER2/R
	    WEIGHT2 = 1-INNER2/R
	ELSE
            WEIGHT1 = 1-INNER2/R
	    WEIGHT2 = INNER2/R
	END IF
    
    END SUBROUTINE COORDINATE_EDGE_POINT
    
    SUBROUTINE COORDINATE_FACE_POINT(V, V1, V2, V3, WEIGHT1, WEIGHT2, WEIGHT3)
        
        REAL(8) :: V(3)
        
        REAL(8) :: V1(3),V2(3),V3(3)
        
        INTEGER :: IMIN
        
        REAL(8) :: A1(3),A2(3),A3(3), TEMPDIST, POINTDIST, DIST1, DIST2, DIST3
        REAL(8) :: N(3), C1(3), C2(3), C3(3)
        REAL(8) :: ROT1, ROT2, ROT3, INNER1, INNER2
        
        REAL(8) :: PTONSURF(3)
        REAL(8) :: L1(3),L2(3),L3(3)
        REAL(8) :: L1_NORM, L2_NORM, L3_NORM
        REAL(8) :: W1,W2,W3, WEIGHT1,WEIGHT2,WEIGHT3
        REAL(8) :: D1, D2
        
        IMIN = 1
        POINTDIST = SQRT(DOT_PRODUCT(V-V1, V-V1))
        TEMPDIST = SQRT(DOT_PRODUCT(V-V2, V-V2))
        IF(TEMPDIST < POINTDIST) THEN
            POINTDIST = TEMPDIST
            IMIN = 2
        END IF
        TEMPDIST = SQRT(DOT_PRODUCT(V-V3, V-V3))
        IF(TEMPDIST < POINTDIST) THEN
            POINTDIST = TEMPDIST
            IMIN = 3
        END IF
        
        IF(IMIN == 1) THEN
            A1(:) = V1
            A2(:) = V2
            A3(:) = V3
        ELSE IF(IMIN == 2) THEN
            A1(:) = V2
            A2(:) = V3
            A3(:) = V1
        ELSE
            A1(:) = V3
            A2(:) = V1
            A3(:) = V2
        END IF
        
        CALL VEC_CURL2(A1,A2,A1,A3,N)
    
        CALL VEC_CURL2(V,A1,A2,A1,C1)
        CALL VEC_CURL2(V,A2,A3,A2,C2)
        CALL VEC_CURL2(V,A3,A1,A3,C3)
        
        ROT1 = DOT_PRODUCT(C1,N)
        ROT2 = DOT_PRODUCT(C2,N)
        ROT3 = DOT_PRODUCT(C3,N)
        
        INNER1 = DOT_PRODUCT(A1-V, A2-A1)
        INNER2 = DOT_PRODUCT(A1-V, A3-A1)
        
        IF(INNER1>0 .AND. INNER2>0) THEN
            W1 = 1.
        ELSE IF(ROT1<0 .AND. ROT2<0 .AND. ROT3<0) THEN
            N = N/SQRT(DOT_PRODUCT(N,N))
            
            D1 = DOT_PRODUCT(V-A1,N)
            
            PTONSURF = V - D1*N
            
            CALL VEC_CURL2(PTONSURF,A2,PTONSURF,A3,  L1)
            CALL VEC_CURL2(PTONSURF,A3,PTONSURF,A1,  L2)
            CALL VEC_CURL2(PTONSURF,A1,PTONSURF,A2,  L3)
            
            L1_NORM = SQRT(DOT_PRODUCT(L1,L1))
            L2_NORM = SQRT(DOT_PRODUCT(L2,L2))
            L3_NORM = SQRT(DOT_PRODUCT(L3,L3))
            
            W1 = L1_NORM/(L1_NORM+L2_NORM+L3_NORM)
            W2 = L2_NORM/(L1_NORM+L2_NORM+L3_NORM)
            W3 = L3_NORM/(L1_NORM+L2_NORM+L3_NORM)
        ELSE
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A1,A2,DIST1)
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A1,A3,DIST2)
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,A2,A3,DIST3)
	    TEMPDIST = MIN(DIST1,DIST2,DIST3)
            IF(TEMPDIST==DIST1) THEN
                D1 = DOT_PRODUCT(A2-A1,A2-V)
                D2 = DOT_PRODUCT(A2-A1,A2-A1)
                W1 = D1/D2
                D1 = DOT_PRODUCT(A1-A2,A1-V)
                W2 = D1/D2
                
                W3 = 0.
            ELSE IF(TEMPDIST==DIST2) THEN
                D1 = DOT_PRODUCT(A3-A1,A3-V)
                D2 = DOT_PRODUCT(A3-A1,A3-A1)
                W1 = D1/D2
                D1 = DOT_PRODUCT(A1-A3,A1-V)
                W3 = D1/D2
                
                W2 = 0.
            ELSE
                D1 = DOT_PRODUCT(A3-A2,A3-V)
                D2 = DOT_PRODUCT(A3-A2,A3-A2)
                W2 = D1/D2
                D1 = DOT_PRODUCT(A2-A3,A2-V)
                W3 = D1/D2
                
                W1 = 0.
            END IF
        END IF
        
        IF(IMIN == 1) THEN
            WEIGHT1 = W1
            WEIGHT2 = W2
            WEIGHT3 = W3
        ELSE IF(IMIN == 2) THEN
            WEIGHT1 = W3
            WEIGHT2 = W1
            WEIGHT3 = W2
        ELSE
            WEIGHT1 = W2
            WEIGHT2 = W3
            WEIGHT3 = W1
        END IF
    END SUBROUTINE COORDINATE_FACE_POINT
    
    SUBROUTINE MESH_QUALITY_SQUARE_ONE(POINT1, POINT2, POINT3, POINT4, INITIAL_AREA,        F_SIZE_SHAPE, F_SIZE, F_SHAPE)
        IMPLICIT NONE
        REAL(8) :: POINT1(3), POINT2(3), POINT3(3), POINT4(3)
        REAL(8) :: INITIAL_AREA
        
        REAL(8) :: A(2,2,4), LAMBDA(2,2,4), ALPHA(4)
        
        REAL(8) :: F_SHAPE, F_SIZE, F_SIZE_SHAPE
        
        INTEGER :: J,K,L
        
        REAL(8) :: POINT_3D(3,4)
        
        REAL(8) :: POINT_2D(2,4)
        
        REAL(8) :: V1(3), V2(3), V3(3), V4(3), V5(3), W1(3), W2(3), W3(3), W4(3), N(3), N1(3), N2(3), N3(3), N4(3), U(3), V(3), R1, R2, AREA
        
            POINT_3D(:,1) = POINT1
            POINT_3D(:,2) = POINT2
            POINT_3D(:,3) = POINT3
            POINT_3D(:,4) = POINT4
            
            W1 = POINT_3D(:,2)-POINT_3D(:,1)
            W2 = POINT_3D(:,4)-POINT_3D(:,1)
            
            W3 = POINT_3D(:,4)-POINT_3D(:,3)
            W4 = POINT_3D(:,2)-POINT_3D(:,3)

	    V1 = POINT_3D(:,3)-POINT_3D(:,2)
	    V2 = POINT_3D(:,1)-POINT_3D(:,2)

	    V3 = POINT_3D(:,1)-POINT_3D(:,4)
	    V4 = POINT_3D(:,3)-POINT_3D(:,4)
            
            V5 = POINT_3D(:,3)-POINT_3D(:,1)
            
            CALL VEC_CURL1(V1,V2,N1)
            
            R1 = SQRT(DOT_PRODUCT(N1,N1))
            
            !N1 = N1/R1
            
            CALL VEC_CURL1(V3,V4,N2)
            
            R2 = SQRT(DOT_PRODUCT(N2,N2))

            !N2 = N2/R2

            AREA = (R1+R2)/2.

            CALL VEC_CURL1(W1,W2,N3)
	    CALL VEC_CURL1(W3,W4,N4)
            
            IF(AREA .NE. 0.) THEN
            
            F_SIZE = MIN(AREA/INITIAL_AREA, INITIAL_AREA/AREA)
            
            N = N1/R1 + N2/R2
            
            N = N/SQRT(DOT_PRODUCT(N,N))
            
            U = V5/SQRT(DOT_PRODUCT(V5,V5))
            CALL VEC_CURL1(N,U,V)
            
            POINT_2D(1,1) = DOT_PRODUCT(POINT_3D(:,1), U)
            POINT_2D(2,1) = DOT_PRODUCT(POINT_3D(:,1), V)
            
            POINT_2D(1,2) = DOT_PRODUCT(POINT_3D(:,2), U)
            POINT_2D(2,2) = DOT_PRODUCT(POINT_3D(:,2), V)
            
            POINT_2D(1,3) = DOT_PRODUCT(POINT_3D(:,3), U)
            POINT_2D(2,3) = DOT_PRODUCT(POINT_3D(:,3), V)
            
            POINT_2D(1,4) = DOT_PRODUCT(POINT_3D(:,4), U)
            POINT_2D(2,4) = DOT_PRODUCT(POINT_3D(:,4), V)
            
            DO L = 1,4
                A(1,1,L) = POINT_2D(1,MOD(L,4)+1) - POINT_2D(1,L)
                A(1,2,L) = POINT_2D(1,MOD(L+2,4)+1) - POINT_2D(1,L)
                A(2,1,L) = POINT_2D(2,MOD(L,4)+1) - POINT_2D(2,L)
                A(2,2,L) = POINT_2D(2,MOD(L+2,4)+1) - POINT_2D(2,L)
            
                ALPHA(L) = A(1,1,L) * A(2,2,L) - A(1,2,L) * A(2,1,L)
                
                DO J = 1,2
                    DO K = 1,2
                        LAMBDA(J,K,L) = (A(1,J,L) * A(1,K,L) + A(2,J,L) * A(2,K,L))
                    END DO
                END DO
            END DO
            
            IF(ALPHA(1) <= 0. .OR. ALPHA(2) <= 0. .OR. ALPHA(3) <= 0. .OR. ALPHA(4) <= 0.) THEN
                F_SHAPE = 0.
            ELSE
                F_SHAPE = 8. / ((LAMBDA(1,1,1) + LAMBDA(2,2,1))/ALPHA(1) + (LAMBDA(1,1,2) + LAMBDA(2,2,2))/ALPHA(2) + (LAMBDA(1,1,3) + LAMBDA(2,2,3))/ALPHA(3) + (LAMBDA(1,1,4) + LAMBDA(2,2,4))/ALPHA(4))
            END IF
            
            ELSE
            
            F_SIZE = 0.
            F_SHAPE = 0.
            
            END IF
            
            F_SIZE_SHAPE = F_SIZE * F_SHAPE
            
            if(f_size_shape<0) then
                j = 1
            end if
    END SUBROUTINE MESH_QUALITY_SQUARE_ONE
    
    SUBROUTINE MESH_QUALITY_SQUARE(POINT, FACE_NUM, FACE, INITIAL_AREA, QUALITY_ARRAY, FILE_NUM, TYP)
        IMPLICIT NONE

        REAL(8) :: POINT(:,:)
        INTEGER :: FACE_NUM
        INTEGER :: FACE(:,:)
        REAL(8) :: INITIAL_AREA(:)
        REAL(8) :: QUALITY_ARRAY(:)
        
        REAL(8) :: F_SHAPE, F_SIZE, F_SIZE_SHAPE
        
        INTEGER, OPTIONAL :: FILE_NUM
        INTEGER, OPTIONAL :: TYP
        
        INTEGER :: I
        
        CHARACTER(100) :: STR, STR2
        
        IF(PRESENT(FILE_NUM)) THEN
            WRITE(STR, *), FILE_NUM
            WRITE(STR2, *), TYP
            STR = './output/surface/feature_quality_square_2d' // TRIM(ADJUSTL(STR)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
            WRITE(*,*) TRIM(STR)
            
            OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        END IF
        
        DO I = 1, FACE_NUM
            CALL MESH_QUALITY_SQUARE_ONE(POINT(:,FACE(1,I)), POINT(:,FACE(2,I)), POINT(:,FACE(3,I)), POINT(:,FACE(4,I)), INITIAL_AREA(I), F_SIZE_SHAPE, F_SIZE, F_SHAPE)
        
            QUALITY_ARRAY(I) = F_SIZE_SHAPE
            
            IF(PRESENT(FILE_NUM)) THEN
                WRITE(21,*) QUALITY_ARRAY(I)
            END IF
        END DO
        
        IF(PRESENT(FILE_NUM)) THEN
            CLOSE(21)
        END IF
    
    END SUBROUTINE MESH_QUALITY_SQUARE
    
    SUBROUTINE MESH_QUALITY_TRIANGLE_ONE(POINT1, POINT2, POINT3, INITIAL_AREA,        F_SIZE_SHAPE, F_SIZE, F_SHAPE)
        IMPLICIT NONE
        REAL(8) :: POINT1(3), POINT2(3), POINT3(3)
        REAL(8) :: INITIAL_AREA
        
        REAL(8) :: A(2,2), LAMBDA(2,2), ALPHA
        
        REAL(8) :: F_SHAPE, F_SIZE, F_SIZE_SHAPE
        
        INTEGER :: J,K
        
        REAL(8) :: POINT_3D(3,3)
        
        REAL(8) :: POINT_2D(2,3)
        
        REAL(8) :: V1(3), V2(3), N(3), U(3), V(3), R, AREA
        
        
            POINT_3D(:,1) = POINT1
            POINT_3D(:,2) = POINT2
            POINT_3D(:,3) = POINT3
            
            V1 = POINT_3D(:,2)-POINT_3D(:,1)
            V2 = POINT_3D(:,3)-POINT_3D(:,1)
            
            CALL VEC_CURL1(V1,V2,N)
            
            R = SQRT(DOT_PRODUCT(N,N))
            
            AREA = R/2.
            
            IF(AREA .NE. 0.) THEN
            
            F_SIZE = MIN(AREA/INITIAL_AREA, INITIAL_AREA/AREA)
            
            N = N/R
            
            U = V1/SQRT(DOT_PRODUCT(V1,V1))
            CALL VEC_CURL1(N,U,V)
            
            POINT_2D(1,1) = DOT_PRODUCT(POINT_3D(:,1), U)
            POINT_2D(2,1) = DOT_PRODUCT(POINT_3D(:,1), V)
            
            POINT_2D(1,2) = DOT_PRODUCT(POINT_3D(:,2), U)
            POINT_2D(2,2) = DOT_PRODUCT(POINT_3D(:,2), V)
            
            POINT_2D(1,3) = DOT_PRODUCT(POINT_3D(:,3), U)
            POINT_2D(2,3) = DOT_PRODUCT(POINT_3D(:,3), V)
            
            A(1,1) = POINT_2D(1,2) - POINT_2D(1,1)
            A(1,2) = POINT_2D(1,3) - POINT_2D(1,1)
            A(2,1) = POINT_2D(2,2) - POINT_2D(2,1)
            A(2,2) = POINT_2D(2,3) - POINT_2D(2,1)
            
            ALPHA = A(1,1) * A(2,2) - A(1,2) * A(2,1)
            
            DO J = 1,2
                DO K = 1,2
                    LAMBDA(J,K) = A(1,J) * A(1,K) + A(2,J) * A(2,K)
                END DO
            END DO
            
            F_SHAPE = SQRT(3.) * ALPHA / (LAMBDA(1,1) + LAMBDA(2,2) - LAMBDA(1,2))
            
            ELSE
                F_SIZE = 0.
                F_SHAPE = 0.
            END IF
            
            F_SIZE_SHAPE = F_SIZE * F_SHAPE
    
    END SUBROUTINE MESH_QUALITY_TRIANGLE_ONE
    
    SUBROUTINE MESH_QUALITY_TRIANGLE(POINT, FACE_NUM, FACE, INITIAL_AREA, QUALITY_ARRAY, FILE_NUM, TYP)
        IMPLICIT NONE

        REAL(8) :: POINT(:,:)
        INTEGER :: FACE_NUM
        INTEGER :: FACE(:,:)
        REAL(8) :: INITIAL_AREA(:)
        REAL(8) :: QUALITY_ARRAY(:)
        
        REAL(8) :: F_SHAPE, F_SIZE, F_SIZE_SHAPE
        
        INTEGER, OPTIONAL :: FILE_NUM
        INTEGER, OPTIONAL :: TYP
        
        INTEGER :: I
        
        CHARACTER(100) :: STR, STR2
        
        IF(PRESENT(FILE_NUM)) THEN
            WRITE(STR, *), FILE_NUM
            WRITE(STR2, *), TYP
            STR = './output/surface/feature_quality_triangle_2d' // TRIM(ADJUSTL(STR)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
            WRITE(*,*) TRIM(STR)
            
            OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        END IF
        
        DO I = 1, FACE_NUM
            CALL MESH_QUALITY_TRIANGLE_ONE(POINT(:,FACE(1,I)), POINT(:,FACE(2,I)), POINT(:,FACE(3,I)), INITIAL_AREA(I), F_SIZE_SHAPE, F_SIZE, F_SHAPE)
        
            QUALITY_ARRAY(I) = F_SHAPE
            
            IF(PRESENT(FILE_NUM)) THEN
                WRITE(21,*) QUALITY_ARRAY(I)
            END IF
        END DO
        
        IF(PRESENT(FILE_NUM)) THEN
            CLOSE(21)
        END IF
    
    END SUBROUTINE MESH_QUALITY_TRIANGLE
    
    
    SUBROUTINE DIHEDRAL_ANGLE(POINT_NUM, POINT, FACE_NUM, FACE, I, I1,        ANGLE)
    IMPLICIT NONE
    INTEGER :: POINT_NUM
    REAL(8) :: POINT(3,POINT_NUM)
    INTEGER :: FACE_NUM
    INTEGER :: FACE(3,FACE_NUM)
    INTEGER :: I, I1
    REAL(8) :: V1(3), V2(3), V3(3), V4(3), N(3), N1(3), R, ANGLE

    V1 = POINT(:,FACE(2,I)) - POINT(:,FACE(1,I))
    V2 = POINT(:,FACE(3,I)) - POINT(:,FACE(1,I))

    CALL VEC_CURL1(V1,V2,N)
    N = N/SQRT(DOT_PRODUCT(N,N))

    V3 = POINT(:,FACE(2,I1)) - POINT(:,FACE(1,I1))
    V4 = POINT(:,FACE(3,I1)) - POINT(:,FACE(1,I1))

    CALL VEC_CURL1(V3,V4,N1)
    N1 = N1/SQRT(DOT_PRODUCT(N1,N1))

    R = DOT_PRODUCT(N,N1)
    ANGLE = ACOS(MAX(-1., MIN(1., R )) )

    END SUBROUTINE DIHEDRAL_ANGLE
    
    
    
    SUBROUTINE FACE_NEIGHBOR_FACE(FACE, CONNECTION_NUM, CONNECTION, J0, DIR, J1)
    IMPLICIT NONE
    INTEGER :: FACE(:,:)
    INTEGER :: CONNECTION_NUM(:)
    INTEGER :: CONNECTION(:,:)
    INTEGER :: J0, J1
    INTEGER :: DIR
    INTEGER :: I0, NUM, I
    I0 = FACE(DIR,J0)
    NUM = CONNECTION_NUM(I0)

    DO I=1,NUM
        IF(J0 == CONNECTION(I,I0)) THEN
            J1 = CONNECTION(MOD(I+NUM-2, NUM) + 1,I0)
            EXIT
        END IF
    END DO

    END SUBROUTINE FACE_NEIGHBOR_FACE
    
    
    
    SUBROUTINE POINT_NEIGHBOR_POINT(FACE, CONNECTION, I0, DIR,      I1)
    IMPLICIT NONE
    INTEGER :: FACE(:,:)
    INTEGER :: CONNECTION(:,:)
    INTEGER :: I0, I1, I
    INTEGER :: DIR
    INTEGER :: J0

    J0 = CONNECTION(DIR,I0)

    DO I=1,3
        IF(I0 == FACE(I,J0)) THEN
            I1 = FACE(MOD(I, 3) + 1,J0)
            EXIT
        END IF
    END DO

    END SUBROUTINE POINT_NEIGHBOR_POINT
    
END MODULE SURFACES_3D
