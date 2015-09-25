MODULE EXTRAPOLATION_2D
    USE SURFACE_MODULE_2D
    USE SURFACES_2D
    USE OPERATORS_2D
    USE SRMS_DATATRANS_MOD_2D
    IMPLICIT NONE
    
    CONTAINS
    
    RECURSIVE SUBROUTINE CONNECTING_INTERFACE_CLUSTER(I0, DIR, POINTLOC_NUM, POINTLOC, POINT_USED, TEMPEDGE_NUM, TEMPEDGE)
        IMPLICIT NONE
        INTEGER :: I0
        INTEGER :: DIR
        INTEGER :: POINTLOC_NUM
        INTEGER :: POINTLOC(2,POINTLOC_NUM)
        INTEGER :: POINT_USED(:,:)
        INTEGER :: TEMPEDGE_NUM
        INTEGER :: TEMPEDGE(:,:)
        
        INTEGER :: I1, PROP_CASE
        INTEGER :: I, IMIN
        
        REAL(8) :: R,R2, RMIN, V(2), W1(2), W2(2)
        
        LOGICAL :: B
        
        IF(DIR<=0) THEN !BEFORE
            IF(POINTLOC(2,I0)==0) THEN
                I1 = SURFACE_FLUID%SURFACE_EDGES(1,SURFACE_FLUID%POINT_EDGE_CONNECTION(1,POINTLOC(1,I0)))
                IF(POINT_USED(1,I1) > 0) THEN
                    CALL CONNECTING_INTERFACE_CLUSTER(POINT_USED(1,I1),-1,POINTLOC_NUM,POINTLOC,POINT_USED,TEMPEDGE_NUM,TEMPEDGE)
                    
                    TEMPEDGE_NUM = TEMPEDGE_NUM + 1
                    TEMPEDGE(1,TEMPEDGE_NUM) = POINT_USED(1,I1)
                    TEMPEDGE(2,TEMPEDGE_NUM) = I0
                END IF
            ELSE
                PROP_CASE = POINTLOC(2,I0) ! 1 or 2
                
                IF(PROP_CASE==1) THEN
                    I1 = SURFACE_PROPEL%SURFACE_EDGES(2,SURFACE_PROPEL%POINT_EDGE_CONNECTION(2,POINTLOC(1,I0)))
                ELSE
                    I1 = SURFACE_CASE%SURFACE_EDGES(2,SURFACE_CASE%POINT_EDGE_CONNECTION(2,POINTLOC(1,I0)))
                END IF
                
                IF(POINT_USED(PROP_CASE+1,I1) > 0) THEN
                    CALL CONNECTING_INTERFACE_CLUSTER(POINT_USED(PROP_CASE+1,I1),-1,POINTLOC_NUM,POINTLOC,POINT_USED,TEMPEDGE_NUM,TEMPEDGE)
                    
                    TEMPEDGE_NUM = TEMPEDGE_NUM + 1
                    TEMPEDGE(1,TEMPEDGE_NUM) = POINT_USED(PROP_CASE+1,I1)
                    TEMPEDGE(2,TEMPEDGE_NUM) = I0
                ELSE
                    
                    B = .FALSE.
                    RMIN = MAX(DOMAIN_MAX(1)-DOMAIN_MIN(1),DOMAIN_MAX(2)-DOMAIN_MIN(2))
                    
                    !!$OMP PARALLEL DO PRIVATE(I,W1,W2,V,R,B,R2)
                    DO I=1,POINTLOC_NUM
                        IF(POINTLOC(2,I).NE.PROP_CASE) THEN
                            IF(PROP_CASE==1) THEN
                                W1 = SURFACE_PROPEL%SURFACE_POINTS(:,POINTLOC(1,I0))
                                W2 = SURFACE_PROPEL%SURFACE_POINTS(:,I1)
                                
                                V = SURFACE_CASE%SURFACE_POINTS(:,POINTLOC(1,I))
                            ELSE
                                W1 = SURFACE_CASE%SURFACE_POINTS(:,POINTLOC(1,I0))
                                W2 = SURFACE_CASE%SURFACE_POINTS(:,I1)
                                
                                V = SURFACE_PROPEL%SURFACE_POINTS(:,POINTLOC(1,I))
                            END IF
                            
                            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,W1,W2,     R)
                            
                            IF(R<SURFACE_PROPEL%MESH_SIZE * 5.) THEN
                                B = .TRUE.
                                R2 = SQRT(DOT_PRODUCT(V-W1,V-W1))
                                IF(R2<RMIN) THEN
                                    RMIN = R2
                                    IMIN = I
                                END IF
                            END IF
                        END IF
                    END DO
                    !!$OMP END PARALLEL DO
                    
                    IF(B) THEN
                        CALL CONNECTING_INTERFACE_CLUSTER(IMIN,-1,POINTLOC_NUM,POINTLOC,POINT_USED,TEMPEDGE_NUM,TEMPEDGE)
                                
                        TEMPEDGE_NUM = TEMPEDGE_NUM + 1
                        TEMPEDGE(1,TEMPEDGE_NUM) = IMIN
                        TEMPEDGE(2,TEMPEDGE_NUM) = I0
                    END IF
                    
                END IF
            END IF
        END IF
        
        IF(DIR>=0) THEN !AFTER
            IF(POINTLOC(2,I0)==0) THEN
                I1 = SURFACE_FLUID%SURFACE_EDGES(2,SURFACE_FLUID%POINT_EDGE_CONNECTION(2,POINTLOC(1,I0)))
                IF(POINT_USED(1,I1) > 0) THEN
                    TEMPEDGE_NUM = TEMPEDGE_NUM + 1
                    TEMPEDGE(1,TEMPEDGE_NUM) = I0
                    TEMPEDGE(2,TEMPEDGE_NUM) = POINT_USED(1,I1)
                    
                    CALL CONNECTING_INTERFACE_CLUSTER(POINT_USED(1,I1),1,POINTLOC_NUM,POINTLOC,POINT_USED,TEMPEDGE_NUM,TEMPEDGE)
                END IF
            ELSE
                PROP_CASE = POINTLOC(2,I0) ! 1 or 2
                
                IF(PROP_CASE==1) THEN
                    I1 = SURFACE_PROPEL%SURFACE_EDGES(1,SURFACE_PROPEL%POINT_EDGE_CONNECTION(1,POINTLOC(1,I0)))
                ELSE
                    I1 = SURFACE_CASE%SURFACE_EDGES(1,SURFACE_CASE%POINT_EDGE_CONNECTION(1,POINTLOC(1,I0)))
                END IF
                
                IF(POINT_USED(PROP_CASE+1,I1) > 0) THEN
                    TEMPEDGE_NUM = TEMPEDGE_NUM + 1
                    TEMPEDGE(1,TEMPEDGE_NUM) = I0
                    TEMPEDGE(2,TEMPEDGE_NUM) = POINT_USED(PROP_CASE+1,I1)
                    
                    CALL CONNECTING_INTERFACE_CLUSTER(POINT_USED(PROP_CASE+1,I1),1,POINTLOC_NUM,POINTLOC,POINT_USED,TEMPEDGE_NUM,TEMPEDGE)
                ELSE
                    
                    B = .FALSE.
                    RMIN = MAX(DOMAIN_MAX(1)-DOMAIN_MIN(1),DOMAIN_MAX(2)-DOMAIN_MIN(2))
                    
                    !!$OMP PARALLEL DO PRIVATE(I,W1,W2,V,R,B,R2)
                    DO I=1,POINTLOC_NUM
                        IF(POINTLOC(2,I).NE.PROP_CASE) THEN
                            IF(PROP_CASE==1) THEN
                                W1 = SURFACE_PROPEL%SURFACE_POINTS(:,POINTLOC(1,I0))
                                W2 = SURFACE_PROPEL%SURFACE_POINTS(:,I1)
                                
                                V = SURFACE_CASE%SURFACE_POINTS(:,POINTLOC(1,I))
                            ELSE
                                W1 = SURFACE_CASE%SURFACE_POINTS(:,POINTLOC(1,I0))
                                W2 = SURFACE_CASE%SURFACE_POINTS(:,I1)
                                
                                V = SURFACE_PROPEL%SURFACE_POINTS(:,POINTLOC(1,I))
                            END IF
                            
                            CALL UNSIGNED_DISTANCE_EDGE_POINT(V,W1,W2,     R)
                            
                            IF(R<SURFACE_PROPEL%MESH_SIZE * 5.) THEN
                                B = .TRUE.
                                R2 = SQRT(DOT_PRODUCT(V-W1,V-W1))
                                IF(R2<RMIN) THEN
                                    RMIN = R2
                                    IMIN = I
                                END IF
                            END IF
                        END IF
                    END DO
                    !!$OMP END PARALLEL DO
                    
                    IF(B) THEN
                        TEMPEDGE_NUM = TEMPEDGE_NUM + 1
                        TEMPEDGE(1,TEMPEDGE_NUM) = I0
                        TEMPEDGE(2,TEMPEDGE_NUM) = IMIN
                        
                        CALL CONNECTING_INTERFACE_CLUSTER(IMIN,1,POINTLOC_NUM,POINTLOC,POINT_USED,TEMPEDGE_NUM,TEMPEDGE)
                    END IF
                    
                END IF
            END IF
        END IF
    
    END SUBROUTINE CONNECTING_INTERFACE_CLUSTER
    
    
    SUBROUTINE FIND_INTERFACE_CLUSTER(TYP)
	IMPLICIT NONE
        INTEGER :: TYP
        INTEGER, ALLOCATABLE :: POINTLOC(:,:)
        INTEGER, ALLOCATABLE :: POINT_USED(:,:)
        
        INTEGER :: NEWPOINT_NUM
        REAL(8), ALLOCATABLE :: NEWPOINT(:,:)
        INTEGER :: NEWEDGE_NUM
        INTEGER, ALLOCATABLE :: NEWEDGE(:,:)
        
        INTEGER, ALLOCATABLE :: TEMPEDGE(:,:)
        INTEGER :: I, J, I1, I2
        INTEGER :: R
        REAL(8) :: R_THRESHOLD, V1(2), V2(2)
        
        INTEGER :: TYP_POINT_NUM
        REAL(8), ALLOCATABLE :: TYP_POINT(:,:)
        INTEGER, ALLOCATABLE :: TYP_POINTLOC(:,:)
        REAL(8) :: R1, R2, T1, T2
        
        R = MAX(MAX(SURFACE_FLUID%SURFACE_POINTS_NUM,SURFACE_PROPEL%SURFACE_POINTS_NUM),SURFACE_CASE%SURFACE_POINTS_NUM)
        
        ALLOCATE(POINT_USED(3,R))
            
        POINT_USED(:,:) = 0
        NEWPOINT_NUM = 0
    
        IF(TYP==0) THEN
            ALLOCATE(POINTLOC(2,SURFACE_FLUID%SURFACE_POINTS_NUM))
            
            DO I=1,SURFACE_FLUID%SURFACE_POINTS_NUM
                I1 = SURFACE_FLUID%POINT_EDGE_CONNECTION(1,I)
                I2 = SURFACE_FLUID%POINT_EDGE_CONNECTION(2,I)
                
                IF(SURFACE_FLUID%EDGE_ONINTERFACE(I1)==1 .OR. SURFACE_FLUID%EDGE_ONINTERFACE(I2)==1 .OR. & 
                    SURFACE_FLUID%EDGE_ONINTERFACE(I1)==2 .OR. SURFACE_FLUID%EDGE_ONINTERFACE(I2)==2) THEN
                    NEWPOINT_NUM = NEWPOINT_NUM + 1
                    POINTLOC(1,NEWPOINT_NUM) = I
                    POINTLOC(2,NEWPOINT_NUM) = 0
                    
                    POINT_USED(1,I) = NEWPOINT_NUM
                END IF
            END DO
        ELSE
            ALLOCATE(POINTLOC(2,SURFACE_PROPEL%SURFACE_POINTS_NUM + SURFACE_CASE%SURFACE_POINTS_NUM))
            
            DO I=1,SURFACE_PROPEL%SURFACE_POINTS_NUM
                I1 = SURFACE_PROPEL%POINT_EDGE_CONNECTION(1,I)
                I2 = SURFACE_PROPEL%POINT_EDGE_CONNECTION(2,I)
                
                !IF(SURFACE_PROPEL%EDGE_ONINTERFACE(I1)==0 .OR. SURFACE_PROPEL%EDGE_ONINTERFACE(I2)==0) THEN
                IF(SURFACE_PROPEL%EDGE_ONINTERFACE(I1)==0 .OR. SURFACE_PROPEL%EDGE_ONINTERFACE(I2)==0 .or. SURFACE_PROPEL%EDGE_ONINTERFACE(I1)==-1 .OR. SURFACE_PROPEL%EDGE_ONINTERFACE(I2)==-1) THEN
                    NEWPOINT_NUM = NEWPOINT_NUM + 1
                    POINTLOC(1,NEWPOINT_NUM) = I
                    POINTLOC(2,NEWPOINT_NUM) = 1
                    
                    POINT_USED(2,I) = NEWPOINT_NUM
                END IF
            END DO
            
            DO I=1,SURFACE_CASE%SURFACE_POINTS_NUM
                I1 = SURFACE_CASE%POINT_EDGE_CONNECTION(1,I)
                I2 = SURFACE_CASE%POINT_EDGE_CONNECTION(2,I)
                
                IF(SURFACE_CASE%EDGE_ONINTERFACE(I1)==0 .OR. SURFACE_CASE%EDGE_ONINTERFACE(I2)==0) THEN
                    NEWPOINT_NUM = NEWPOINT_NUM + 1
                    POINTLOC(1,NEWPOINT_NUM) = I
                    POINTLOC(2,NEWPOINT_NUM) = 2
                    
                    POINT_USED(3,I) = NEWPOINT_NUM
                END IF
            END DO
        END IF
        
        ALLOCATE(NEWPOINT(2,NEWPOINT_NUM))
        !$OMP PARALLEL DO PRIVATE(I)
        DO I=1,NEWPOINT_NUM
            IF(POINTLOC(2,I)==0) THEN
                NEWPOINT(:,I) = SURFACE_FLUID%SURFACE_POINTS(:,POINTLOC(1,I))
            ELSE IF(POINTLOC(2,I)==1) THEN
                NEWPOINT(:,I) = SURFACE_PROPEL%SURFACE_POINTS(:,POINTLOC(1,I))
            ELSE IF(POINTLOC(2,I)==2) THEN
                NEWPOINT(:,I) = SURFACE_CASE%SURFACE_POINTS(:,POINTLOC(1,I))
            END IF
        END DO
        !$OMP END PARALLEL DO
        
        ALLOCATE(TEMPEDGE(2,NEWPOINT_NUM))
        NEWEDGE_NUM = 0
        CALL CONNECTING_INTERFACE_CLUSTER(1, 0, NEWPOINT_NUM, POINTLOC,POINT_USED,NEWEDGE_NUM,TEMPEDGE)
        
        ALLOCATE(NEWEDGE(2,NEWEDGE_NUM))
        
	!$OMP PARALLEL DO PRIVATE(I)
        DO I=1,NEWEDGE_NUM
            NEWEDGE(:,I) = TEMPEDGE(:,I)
        END DO
        !$OMP END PARALLEL DO
        
        DEALLOCATE(TEMPEDGE)
        DEALLOCATE(POINT_USED)
        
        
        ALLOCATE(TYP_POINT(2,NEWPOINT_NUM))
        ALLOCATE(TYP_POINTLOC(2,NEWPOINT_NUM))
        
        TYP_POINT_NUM = NEWEDGE_NUM + 1
        
	!$OMP PARALLEL DO PRIVATE(I)
        DO I=1,NEWEDGE_NUM
            TYP_POINT(:,I) = NEWPOINT(:,NEWEDGE(1,I))
            TYP_POINTLOC(:,I) = POINTLOC(:,NEWEDGE(1,I))
        END DO
	!$OMP END PARALLEL DO
        
        TYP_POINT(:,TYP_POINT_NUM) = NEWPOINT(:,NEWEDGE(2,NEWEDGE_NUM))
        TYP_POINTLOC(:,TYP_POINT_NUM) = POINTLOC(:,NEWEDGE(2,NEWEDGE_NUM))
        
        I = 1
        DO WHILE(I<=TYP_POINT_NUM-2)
            IF(TYP_POINTLOC(2,I)==0) THEN
                V1 = SURFACE_FLUID%SURFACE_POINTS(:,SURFACE_FLUID%SURFACE_EDGES(1,SURFACE_FLUID%POINT_EDGE_CONNECTION(1,TYP_POINTLOC(1,I)) ) ) - SURFACE_FLUID%SURFACE_POINTS(:,TYP_POINTLOC(1,I))
                V2 = SURFACE_FLUID%SURFACE_POINTS(:,SURFACE_FLUID%SURFACE_EDGES(2,SURFACE_FLUID%POINT_EDGE_CONNECTION(2,TYP_POINTLOC(1,I)) ) ) - SURFACE_FLUID%SURFACE_POINTS(:,TYP_POINTLOC(1,I))
            ELSE IF(TYP_POINTLOC(2,I)==1) THEN
                V1 = SURFACE_PROPEL%SURFACE_POINTS(:,SURFACE_PROPEL%SURFACE_EDGES(1,SURFACE_PROPEL%POINT_EDGE_CONNECTION(1,TYP_POINTLOC(1,I)) ) ) - SURFACE_PROPEL%SURFACE_POINTS(:,TYP_POINTLOC(1,I))
                V2 = SURFACE_PROPEL%SURFACE_POINTS(:,SURFACE_PROPEL%SURFACE_EDGES(2,SURFACE_PROPEL%POINT_EDGE_CONNECTION(2,TYP_POINTLOC(1,I)) ) ) - SURFACE_PROPEL%SURFACE_POINTS(:,TYP_POINTLOC(1,I))
            ELSE IF(TYP_POINTLOC(2,I)==2) THEN
                V1 = SURFACE_CASE%SURFACE_POINTS(:,SURFACE_CASE%SURFACE_EDGES(1,SURFACE_CASE%POINT_EDGE_CONNECTION(1,TYP_POINTLOC(1,I)) ) ) - SURFACE_CASE%SURFACE_POINTS(:,TYP_POINTLOC(1,I))
                V2 = SURFACE_CASE%SURFACE_POINTS(:,SURFACE_CASE%SURFACE_EDGES(2,SURFACE_CASE%POINT_EDGE_CONNECTION(2,TYP_POINTLOC(1,I)) ) ) - SURFACE_CASE%SURFACE_POINTS(:,TYP_POINTLOC(1,I))
            END IF
            R_THRESHOLD = (SQRT(DOT_PRODUCT(V1,V1)) + SQRT(DOT_PRODUCT(V2,V2)))/2. / 10.
            
            IF(SQRT(DOT_PRODUCT(TYP_POINT(:,I) - TYP_POINT(:,I+1), TYP_POINT(:,I) - TYP_POINT(:,I+1))) < R_THRESHOLD) THEN
                
                IF(TYP_POINTLOC(2,I)==1) THEN
                    DO J=I+2,TYP_POINT_NUM
                        TYP_POINT(:,J-1) = TYP_POINT(:,J)
                        TYP_POINTLOC(:,J-1) = TYP_POINTLOC(:,J)
                    END DO
                ELSE
                    DO J=I+1,TYP_POINT_NUM
                        TYP_POINT(:,J-1) = TYP_POINT(:,J)
                        TYP_POINTLOC(:,J-1) = TYP_POINTLOC(:,J)
                    END DO
                END IF
                
                TYP_POINT_NUM = TYP_POINT_NUM - 1
            ELSE IF(I>=2) THEN
                CALL UNSIGNED_DISTANCE_EDGE_POINT(TYP_POINT(:,I-1),TYP_POINT(:,I),TYP_POINT(:,I+1),R1)
                CALL UNSIGNED_DISTANCE_EDGE_POINT(TYP_POINT(:,I+1),TYP_POINT(:,I-1),TYP_POINT(:,I),R2)
                CALL COORDINATE_EDGE_POINT(TYP_POINT(:,I-1),TYP_POINT(:,I),TYP_POINT(:,I+1),T1)
                CALL COORDINATE_EDGE_POINT(TYP_POINT(:,I+1),TYP_POINT(:,I-1),TYP_POINT(:,I),T2)
                
                IF((R1 < R_THRESHOLD .AND. T1 > 0. .AND. T1 < 1.) .OR. (R2 < R_THRESHOLD .AND. T2 > 0. .AND. T2 < 1.)) THEN
                    DO J=I+1,TYP_POINT_NUM
                        TYP_POINT(:,J-1) = TYP_POINT(:,J)
                        TYP_POINTLOC(:,J-1) = TYP_POINTLOC(:,J)
                    END DO
                
                    TYP_POINT_NUM = TYP_POINT_NUM - 1
                ELSE
                    I = I + 1
                END IF
            ELSE
                I = I + 1
            END IF
        END DO
        
        DEALLOCATE(NEWPOINT)
        DEALLOCATE(NEWEDGE)
        DEALLOCATE(POINTLOC)
        
        
        IF(TYP==0) THEN
            IF(ALLOCATED(INTERFACE_FLUID_POINTS)) THEN
                DEALLOCATE(INTERFACE_FLUID_POINTS)
            END IF
            IF(ALLOCATED(INTERFACE_FLUID_POINTS_LOC)) THEN
                DEALLOCATE(INTERFACE_FLUID_POINTS_LOC)
            END IF
            
            INTERFACE_FLUID_POINTS_NUM = TYP_POINT_NUM
            ALLOCATE(INTERFACE_FLUID_POINTS(2,INTERFACE_FLUID_POINTS_NUM))
            ALLOCATE(INTERFACE_FLUID_POINTS_LOC(2,INTERFACE_FLUID_POINTS_NUM))
	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1,INTERFACE_FLUID_POINTS_NUM
                INTERFACE_FLUID_POINTS(:,I) = TYP_POINT(:,I)
                INTERFACE_FLUID_POINTS_LOC(:,I) = TYP_POINTLOC(:,I)
            END DO
	    !$OMP END PARALLEL DO
        ELSE
            IF(ALLOCATED(INTERFACE_STRUCT_POINTS)) THEN
                DEALLOCATE(INTERFACE_STRUCT_POINTS)
            END IF
            IF(ALLOCATED(INTERFACE_STRUCT_POINTS_LOC)) THEN
                DEALLOCATE(INTERFACE_STRUCT_POINTS_LOC)
            END IF
            
            INTERFACE_STRUCT_POINTS_NUM = TYP_POINT_NUM
            ALLOCATE(INTERFACE_STRUCT_POINTS(2,INTERFACE_STRUCT_POINTS_NUM))
            ALLOCATE(INTERFACE_STRUCT_POINTS_LOC(2,INTERFACE_STRUCT_POINTS_NUM))
	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1,INTERFACE_STRUCT_POINTS_NUM
                INTERFACE_STRUCT_POINTS(:,I) = TYP_POINT(:,I)
                INTERFACE_STRUCT_POINTS_LOC(:,I) = TYP_POINTLOC(:,I)
            END DO
	    !$OMP END PARALLEL DO
        END IF
        
        DEALLOCATE(TYP_POINT)
        DEALLOCATE(TYP_POINTLOC)
    
    END SUBROUTINE FIND_INTERFACE_CLUSTER
    
    SUBROUTINE UPDATE_INTERFACE_CLUSTER(TYP)
        IMPLICIT NONE
        INTEGER :: TYP
        INTEGER :: I
        
        IF(TYP==0) THEN
	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1,INTERFACE_FLUID_POINTS_NUM
                INTERFACE_FLUID_POINTS(:,I) = SURFACE_FLUID%SURFACE_POINTS(:,INTERFACE_FLUID_POINTS_LOC(1,I))
            END DO
	    !$OMP END PARALLEL DO
        ELSE
	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1,INTERFACE_STRUCT_POINTS_NUM
                IF(INTERFACE_STRUCT_POINTS_LOC(2,I)==1) THEN
                    INTERFACE_STRUCT_POINTS(:,I) = SURFACE_PROPEL%SURFACE_POINTS(:,INTERFACE_STRUCT_POINTS_LOC(1,I))
                ELSE IF(INTERFACE_STRUCT_POINTS_LOC(2,I)==2) THEN
                    INTERFACE_STRUCT_POINTS(:,I) = SURFACE_CASE%SURFACE_POINTS(:,INTERFACE_STRUCT_POINTS_LOC(1,I))
                END IF
            END DO
	    !$OMP END PARALLEL DO
        END IF
        
    END SUBROUTINE UPDATE_INTERFACE_CLUSTER
    
    SUBROUTINE EXTRAPOLATION_COMMON_REFINEMENT(TYP1, TYP2, DATA1, DATA2, TARGET1, TARGET2)
	IMPLICIT NONE
        INTEGER :: TYP1, TYP2
    
        INTEGER :: I
        
        INTEGER :: TYP1_POINT_NUM
        REAL(8), POINTER, DIMENSION(:,:) :: TYP1_POINT
        INTEGER, POINTER, DIMENSION(:,:) :: TYP1_POINTLOC
        
        INTEGER :: TYP2_POINT_NUM
        REAL(8), POINTER, DIMENSION(:,:) :: TYP2_POINT
        INTEGER, POINTER, DIMENSION(:,:) :: TYP2_POINTLOC
        
        REAL(8) :: DATA1(:), DATA2(:), TARGET1(:), TARGET2(:)
        
        REAL(8), ALLOCATABLE :: CR_SOURCE(:), CR_TARGET(:)
        
        CALL UPDATE_INTERFACE_CLUSTER(0)
        CALL UPDATE_INTERFACE_CLUSTER(1)
        
        IF (TYP1==0) THEN
            TYP1_POINT_NUM = INTERFACE_FLUID_POINTS_NUM
            TYP1_POINT => INTERFACE_FLUID_POINTS
            TYP1_POINTLOC => INTERFACE_FLUID_POINTS_LOC
        ELSE
            TYP1_POINT_NUM = INTERFACE_STRUCT_POINTS_NUM
            TYP1_POINT => INTERFACE_STRUCT_POINTS
            TYP1_POINTLOC => INTERFACE_STRUCT_POINTS_LOC
        END IF
        
        IF (TYP2==0) THEN
            TYP2_POINT_NUM = INTERFACE_FLUID_POINTS_NUM
            TYP2_POINT => INTERFACE_FLUID_POINTS
            TYP2_POINTLOC => INTERFACE_FLUID_POINTS_LOC
        ELSE
            TYP2_POINT_NUM = INTERFACE_STRUCT_POINTS_NUM
            TYP2_POINT => INTERFACE_STRUCT_POINTS
            TYP2_POINTLOC => INTERFACE_STRUCT_POINTS_LOC
        END IF
        
        
        ALLOCATE(CR_SOURCE(TYP1_POINT_NUM))
        
	!$OMP PARALLEL DO PRIVATE(I)
        DO I=1,TYP1_POINT_NUM
            IF(TYP1_POINTLOC(2,I)==0 .OR. TYP1_POINTLOC(2,I)==1) THEN
                CR_SOURCE(I) = DATA1(TYP1_POINTLOC(1,I))
            ELSE
                CR_SOURCE(I) = DATA2(TYP1_POINTLOC(1,I))
            END IF
        END DO
        !$OMP END PARALLEL DO
        
        ALLOCATE(CR_TARGET(TYP2_POINT_NUM))
        
        
        CALL DATA_TRANSFER(TYP1_POINT_NUM, TYP1_POINT, CR_SOURCE, TYP2_POINT_NUM, TYP2_POINT, CR_TARGET)
        
        
        !$OMP PARALLEL DO PRIVATE(I)
        DO I=1,TYP2_POINT_NUM
            IF(TYP2_POINTLOC(2,I)==0 .OR. TYP2_POINTLOC(2,I)==1) THEN
                TARGET1(TYP2_POINTLOC(1,I)) = CR_TARGET(I)
            ELSE
                TARGET2(TYP2_POINTLOC(1,I)) = CR_TARGET(I)
            END IF
        END DO
        !$OMP END PARALLEL DO
        
        DEALLOCATE(CR_SOURCE)
        DEALLOCATE(CR_TARGET)
        
    END SUBROUTINE EXTRAPOLATION_COMMON_REFINEMENT
    
END MODULE
