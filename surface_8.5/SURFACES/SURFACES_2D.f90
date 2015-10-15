MODULE SURFACES_2D
    USE SURFACE_MODULE_2D
    
    CONTAINS
    
    SUBROUTINE UNSIGNED_DISTANCE_EDGE_POINT(V,W1,W2,DIST)
    
        REAL(8) :: V(2),W1(2),W2(2)
    
        REAL(8) :: DIST1, DIST2, POINTDIST, DIST
        REAL(8) :: A1(2),A2(2)
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
    
        REAL(8) :: V(2),W1(2),W2(2),L(2), R
    
        REAL(8) :: DIST1, DIST2, POINTDIST
        REAL(8) :: A1(2),A2(2)
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
    
    SUBROUTINE COORDINATE_EDGE_POINT(V,W1,W2,COORD)
    
        REAL(8) :: V(2),W1(2),W2(2),L(2), R
    
        REAL(8) :: DIST1, DIST2, POINTDIST
        REAL(8) :: A1(2),A2(2)
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
        
	IF(FLAG) THEN
	    COORD = 1-INNER2/R
	ELSE
	    COORD = INNER2/R
	END IF
    
    END SUBROUTINE COORDINATE_EDGE_POINT
    
    SUBROUTINE LINE_EDGE_INTERSECTING(POINT_NUM, POINT, EDGE_NUM, EDGE, I0, L, V, SGN, T)
    
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(2,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: EDGE_NUM
        INTEGER :: EDGE(2,EDGE_NUM) ! EDGE(2,EDGE_NUM)
        INTEGER :: I0
        REAL(8) :: L(2), V(2) 
        REAL(8) :: A(2),B(2), C(2)
        REAL(8) :: L1(2)
        REAL(8) :: A1,B1,D1,D2, D3,D4
        INTEGER :: SGN
        REAL(8) :: S,T
        
        A = POINT(:,EDGE(1,I0))
        B = POINT(:,EDGE(2,I0))
    
        L1(1) = -L(2)
        L1(2) = L(1)
        
        C(1) = -(B(2)-A(2))
        C(2) = B(1) - A(1)
        
        A1 = DOT_PRODUCT(V-A,L1)
        B1 = DOT_PRODUCT(B-A,L1)
        D1 = DOT_PRODUCT(V-A,L)
        D2 = DOT_PRODUCT(B-A,L)
        
        D3 = DOT_PRODUCT(A-V,C)
        D4 = DOT_PRODUCT(L,C)
        
        !T = D3/D4

        !write(*,*) 'temp'
        
        IF(ABS(B1)<MINERROR) THEN
            SGN = -2
	    !WRITE(*,*) 'LINE_EDGE_INTERSECTING : PARALLEL!'
            RETURN
        END IF
        
        S = A1/B1
        T = S*D2 - D1
        
        IF(S<0 .OR. S>1) THEN
            SGN = -3
	    !WRITE(*,*) 'LINE_EDGE_INTERSECTING : INTERSECTION OCCURED OUTSIDE OF THE EDGE!'
        ELSE IF( T>=0 .AND. S>0 .AND. S<1) THEN
            SGN = 1
        ELSE IF(T>=0 .AND. (S==0 .OR. S==1)) THEN
            SGN = -1
        ELSE
            SGN = 0
        END IF
    
    END SUBROUTINE LINE_EDGE_INTERSECTING
    
    SUBROUTINE LINE_SURFACE_INTERSECTING(POINT_NUM, POINT, EDGE_NUM, EDGE,V,RETURNS)

        INTEGER :: POINT_NUM
        REAL(8) :: POINT(2,POINT_NUM)
        INTEGER :: EDGE_NUM
        INTEGER :: EDGE(2,EDGE_NUM)
        INTEGER :: NUM
        INTEGER :: SGN
        REAL(8) :: L(2), V(2), T
        REAL(8) :: RETURNS
        
        INTEGER :: I
        
        SGN = -1
        
        DO WHILE(SGN == -1 .OR. SGN == -2)
            NUM = 0
            
            CALL RANDOM_NUMBER(L)
            L = L/SQRT(DOT_PRODUCT(L,L))
            
            DO I=1,EDGE_NUM
                CALL LINE_EDGE_INTERSECTING(POINT_NUM, POINT, EDGE_NUM, EDGE,I,L,V,SGN,T)
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
    SUBROUTINE LINE_SURFACE_INTERSECTING_DISTANCE(POINT_NUM, POINT, EDGE_NUM, EDGE, L, V, TMIN, IMIN)

        INTEGER :: POINT_NUM
        REAL(8) :: POINT(2,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: EDGE_NUM
        INTEGER :: EDGE(2,EDGE_NUM) ! EDGE(2,EDGE_NUM)
        
        INTEGER :: NUM
        INTEGER :: SGN
        
        REAL(8) :: L(2), L_TEMP(2), V(2), T, TMIN
        
        INTEGER :: I, IMIN

        NUM = 0
        L_TEMP = L/SQRT(DOT_PRODUCT(L,L))

        DO I=1,EDGE_NUM
            CALL LINE_EDGE_INTERSECTING(POINT_NUM, POINT, EDGE_NUM, EDGE,I,L_TEMP,V,SGN,T)

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
    
    SUBROUTINE PROJECTION_SURFACE_POINT(V, POINT_NUM, POINT, EDGE_NUM, EDGE, COORD)

        REAL(8) :: V(2)
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(2,POINT_NUM)
        INTEGER :: EDGE_NUM
        INTEGER :: EDGE(2,EDGE_NUM) 
        REAL(8) :: COORD
        REAL(8) :: DIST, TEMPDIST
        INTEGER :: I, MINI
  
        DO I=1,EDGE_NUM
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,POINT(:,EDGE(1,I)),POINT(:,EDGE(2,I)),TEMPDIST)
      
            IF (I==1) THEN
                MINI = I
                DIST = TEMPDIST
            ELSEIF(TEMPDIST <= DIST) THEN
                MINI = I
                DIST = TEMPDIST
            END IF
        END DO
        
        CALL PROJECTION_EDGE_POINT(V,POINT(:,EDGE(1,MINI)),POINT(:,EDGE(2,MINI)),COORD)
        
    END SUBROUTINE PROJECTION_SURFACE_POINT
    
    SUBROUTINE DISTANCE_SURFACE_POINT(POINT_NUM, POINT, EDGE_NUM, EDGE, V, DIST)

        INTEGER :: POINT_NUM
        REAL(8) :: POINT(2,POINT_NUM)
        INTEGER :: EDGE_NUM
        INTEGER :: EDGE(2,EDGE_NUM) 
        REAL(8) :: DIST, TEMPDIST
        REAL(8) :: V(2)
        !REAL(8) :: TEMP
        
        INTEGER :: I
  
        DO I=1,EDGE_NUM
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,POINT(:,EDGE(1,I)),POINT(:,EDGE(2,I)),TEMPDIST)
      
            IF (I==1) THEN
                DIST = TEMPDIST
            ELSEIF(TEMPDIST <= DIST) THEN
                DIST = TEMPDIST
            END IF
        END DO
  
        !IF (DIST < MINERROR) THEN
        !    DIST = 0.
        !    RETURN
        !END IF
        !
        !CALL LINE_SURFACE_INTERSECTING(POINT_NUM, POINT, EDGE_NUM, EDGE,V,TEMP)
        !DIST = TEMP * DIST
  
    END SUBROUTINE DISTANCE_SURFACE_POINT
    
    SUBROUTINE FIND_PERPENDICULAR_POINT(POINT_NUM, POINT, EDGE_NUM, EDGE, V,    VEC)
        
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(2,POINT_NUM)
        INTEGER :: EDGE_NUM
        INTEGER :: EDGE(2,EDGE_NUM)
        
        REAL(8) :: V(2), VEC(2)
        
        REAL(8) :: TEMPDIST, POINTDIST
        REAL(8) :: DIST
        INTEGER :: I0
        
        INTEGER :: I, IMIN
        
        REAL(8) :: A1(2),A2(2)
        INTEGER :: B1, B2
        REAL(8) :: INNER1, INNER2
        
        REAL(8) :: L(2)
               
        I0 = 1
    
        DO I=1,EDGE_NUM
            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,POINT(:,EDGE(1,I)),POINT(:,EDGE(2,I)),TEMPDIST)
      
            IF (I==1) THEN
                DIST = TEMPDIST
                I0 = I
            ELSEIF(TEMPDIST <= DIST) THEN
                DIST = TEMPDIST
                I0 = I
            END IF
        END DO
        
        
    
        IMIN = 1
        
        DO I=1,2
            TEMPDIST = SQRT(DOT_PRODUCT(V-POINT(:,EDGE(I,I0)),V-POINT(:,EDGE(I,I0))))
        
            IF (I==1) THEN
                POINTDIST = TEMPDIST
            ELSEIF(TEMPDIST < POINTDIST) THEN
                POINTDIST = TEMPDIST
                IMIN = I
            END IF
        END DO     
            
        IF(IMIN == 1) THEN
            A1 = POINT(:,EDGE(1,I0))
            A2 = POINT(:,EDGE(2,I0))
                                
            B1 = EDGE(1,I0)
            B2 = EDGE(2,I0)
        ELSE
            A1 = POINT(:,EDGE(2,I0))
            A2 = POINT(:,EDGE(1,I0))
                
            B1 = EDGE(2,I0)
            B2 = EDGE(1,I0)
        END IF
        
        INNER1 = DOT_PRODUCT(V-A1,A2-A1)
    
        IF(INNER1<0) THEN
            VEC = POINT(:,B1)
        ELSE 
            L = A2-A1
            L = L/SQRT(DOT_PRODUCT(L,L))
            INNER2 = DOT_PRODUCT(V-A1,L)
             
            VEC = A1 + INNER2 * L
        END IF
     
    END SUBROUTINE FIND_PERPENDICULAR_POINT
    
    SUBROUTINE FIND_B_RATE_EDGE(V, V1, V2, B_RATE1, B_RATE2, B_RATE_CENTER,         INTER_B_RATE)
        REAL(8) :: V(2)
        REAL(8) :: V1(2),V2(2)
        REAL(8) :: B_RATE1, B_RATE2, B_RATE_CENTER
        
        REAL(8) :: A1(2), A2(2), INNER1, INNER2, L(2), B1, B2, BC
        REAL(8) :: PTONEDGE(2), L1(2), L2(2), L1_NORM, L2_NORM, W1, W2
        
        REAL(8) :: INTER_B_RATE
        
        IF(DOT_PRODUCT(V-V1,V-V1) < DOT_PRODUCT(V-V2,V-V2)) THEN
            A1 = V1
            A2 = V2
            
            B1 = B_RATE1
            B2 = B_RATE2
            BC = B_RATE_CENTER
        ELSE
            A1 = V2
            A2 = V1
            
            B1 = B_RATE2
            B2 = B_RATE1
            BC = B_RATE_CENTER
        ENDIF
        
        INNER1 = DOT_PRODUCT(V-A1,A2-A1)
    
        IF(INNER1<0) THEN
            INTER_B_RATE = B1
        ELSE 
            L = A2-A1
            L = L/SQRT(DOT_PRODUCT(L,L))
            INNER2 = DOT_PRODUCT(V-A1,L)
             
            PTONEDGE = A1 + INNER2 * L
            L1 = PTONEDGE - A1
            L2 = PTONEDGE - A2
            
            L1_NORM = SQRT(DOT_PRODUCT(L1,L1))
            L2_NORM = SQRT(DOT_PRODUCT(L2,L2))
                 
            W1 = L1_NORM/(L1_NORM+L2_NORM)
            W2 = L2_NORM/(L1_NORM+L2_NORM)
             
            INTER_B_RATE = B1 * (1.-2.*W1)*W2 + B2 * W1*(1.-2.*W2) + BC * 4.*W1*W2
    
        END IF
        
        
    END SUBROUTINE FIND_B_RATE_EDGE
    
    SUBROUTINE FIND_VELOCITY_EDGE(V, V1, V2, VELOCITY1, VELOCITY2,         INTER_VELOCITY)
        REAL(8) :: V(2)
        REAL(8) :: V1(2), V2(2)
        REAL(8) :: VELOCITY1(2), VELOCITY2(2)
        
        REAL(8) :: A1(2), A2(2), INNER1, INNER2, L(2), VEL1(2), VEL2(2)
        REAL(8) :: PTONEDGE(2), L1(2), L2(2), L1_NORM, L2_NORM, W1, W2
        
        REAL(8) :: INTER_VELOCITY(2)
        
        IF(DOT_PRODUCT(V-V1,V-V1) < DOT_PRODUCT(V-V2,V-V2)) THEN
            A1 = V1
            A2 = V2
            
            VEL1 = VELOCITY1
            VEL2 = VELOCITY2
        ELSE
            A1 = V2
            A2 = V1
            
            VEL1 = VELOCITY2
            VEL2 = VELOCITY1
        ENDIF
        
        INNER1 = DOT_PRODUCT(V-A1,A2-A1)
    
        IF(INNER1<0) THEN
            INTER_VELOCITY = VEL1
        ELSE 
            L = A2-A1
            L = L/SQRT(DOT_PRODUCT(L,L))
            INNER2 = DOT_PRODUCT(V-A1,L)
             
            PTONEDGE = A1 + INNER2 * L
            L1 = PTONEDGE - A1
            L2 = PTONEDGE - A2
            
            L1_NORM = SQRT(DOT_PRODUCT(L1,L1))
            L2_NORM = SQRT(DOT_PRODUCT(L2,L2))
                 
            W1 = L1_NORM/(L1_NORM+L2_NORM)
            W2 = L2_NORM/(L1_NORM+L2_NORM)
             
            INTER_VELOCITY = VEL1 * W2 + VEL2 * W1
    
        END IF
    END SUBROUTINE FIND_VELOCITY_EDGE
    
    SUBROUTINE TWO_EDGE_INTERSECTING(V1, V2, V3, V4, TF, V)
        REAL(8) :: V1(2), V2(2), V3(2), V4(2)
        LOGICAL :: TF
        REAL(8) :: V(2)
        
        REAL(8) :: W1(2), W2(2), W3(2),W4(2), INN, INN2, N0, N1, N2, N3, N4, C, C2, INN3, INN4, T, S
        
        N0 = SQRT(DOT_PRODUCT(V1-V3,V1-V3))
        IF(N0 < MINERROR2) THEN
            TF = .TRUE.
            V = (V1+V3)/2.
            RETURN
        END IF
        
        W1 = V4-V3
        W2(1) = -W1(2)
        W2(2) = W1(1)
        
        W3 = V2-V1
        W4(1) = -W3(2)
        W4(2) = W3(1)
        
        INN = DOT_PRODUCT(V2-V1,W2)
        N1 = SQRT(DOT_PRODUCT(V2-V1,V2-V1))
        N2 = SQRT(DOT_PRODUCT(W2,W2))
        C = INN/(N1*N2)
        
        INN3 = DOT_PRODUCT(V4-V3,W4)
        N3 = SQRT(DOT_PRODUCT(V4-V3,V4-V3))
        N4 = SQRT(DOT_PRODUCT(W4,W4))
        C2 = INN3/(N3*N4)
        
        IF(ABS(C) < 0.01 .OR. ABS(C2) < 0.01) THEN
            TF = .FALSE.
            RETURN
        END IF
        
        INN2 = DOT_PRODUCT(V3-V1,W2)
        T = INN2 / INN
        
        INN4 = DOT_PRODUCT(V1-V3,W4)
        S = INN4 / INN3
        
        IF(T<-MINERROR .OR. T > 1+MINERROR .OR. S<-MINERROR .OR. S>1+MINERROR) THEN
            TF = .FALSE.
        ELSE
            TF = .TRUE.
        END IF
        V = V1 + T*(V2-V1)
    END SUBROUTINE TWO_EDGE_INTERSECTING
    
END MODULE
