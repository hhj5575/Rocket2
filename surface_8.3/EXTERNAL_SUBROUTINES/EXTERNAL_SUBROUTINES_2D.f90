MODULE SURFACE_EXTERNAL_COMMON_2D

    USE SURFACE_MODULE_2D
    USE PROPA_RECONST_REINITIAL_2D
    USE OPERATORS_2D
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE READINPUT

	! Local Variables
	INTEGER :: io
	CHARACTER(500) :: cDum
	CHARACTER(500) :: Fname
	LOGICAL :: Exist
!	TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT

!        IF (TYP==0) THEN
!            SURFACE_CURRENT => SURFACE_FLUID
!        ELSE IF (TYP==1) THEN
!            SURFACE_CURRENT => SURFACE_PROPEL
!        ELSE IF (TYP==2) THEN
!            SURFACE_CURRENT => SURFACE_CASE
!        END IF


	! Find Input file
	Fname = "./input/surface.inp"
	INQUIRE(FILE=Fname, EXIST=Exist)
	IF( .NOT. Exist ) THEN
	  WRITE(*,*) 'surface >> Error : Cannot find input file < surface.inp >'
	  STOP
	END IF

	! Read Input File
	OPEN(newunit=io,File=Fname)
	READ(io,*) 
	READ(io,*)
	READ(io,*) cDum, cDum, RESTART_FLAG
	READ(io,*) cDum, cDum, RESTART_ITERATION 
	READ(io,*)
    READ(io,*)
    READ(io,*)
    READ(io,*)
    READ(io,*)
	READ(io,*) cDum, cDum, CHI_C_LOW
	READ(io,*) cDum, cDum, CHI_C_HIGH
	READ(io,*) cDum, cDum, INTERFACE_THRESHOLD
	READ(io,*)
	READ(io,*) cDum, cDum, SMALL_REGION_POINT_NUM
	READ(io,*) cDum, cDum, THIN_REGION_EXISTENCE
	READ(io,*) cDum, cDum, THIN_REGION_ATTACHMENT
	READ(io,*) cDum, cDum, EDGE_SPLITTING
	READ(io,*) cDum, cDum, EDGE_COLLAPSING
	CLOSE(io)

	END SUBROUTINE

    
    SUBROUTINE SAVINGPOINT(FILE_NUM, TYP)
	IMPLICIT NONE
        INTEGER :: TYP
        
        INTEGER :: FILE_NUM
        INTEGER :: I
        CHARACTER(500) :: STR, STR2
        
        INTEGER :: POINT_NUM
        REAL(8), POINTER, DIMENSION(:,:) :: POINT
        
        IF (TYP==0) THEN
            POINT_NUM = SURFACE_FLUID%SURFACE_POINTS_NUM
            POINT => SURFACE_FLUID%SURFACE_POINTS
        END IF
        IF (TYP==1) THEN
            POINT_NUM = SURFACE_PROPEL%SURFACE_POINTS_NUM
            POINT => SURFACE_PROPEL%SURFACE_POINTS
        END IF
        IF (TYP==2) THEN
            POINT_NUM = SURFACE_CASE%SURFACE_POINTS_NUM
            POINT => SURFACE_CASE%SURFACE_POINTS
        END IF
        
        WRITE(STR, *), FILE_NUM
        WRITE(STR2, *), TYP
        STR = './output/surface/pointedges/point_2d' // TRIM(ADJUSTL(STR)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        WRITE(*,*) TRIM(STR)
    
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        
        !$OMP DO ORDERED 
        DO I = 1, POINT_NUM
            WRITE(21,*) POINT(1,I)
            WRITE(21,*) POINT(2,I)
        END DO
        !$OMP END DO
    
        CLOSE(21)
    END SUBROUTINE SAVINGPOINT
    
    SUBROUTINE SAVINGEDGE(FILE_NUM, TYP)
	IMPLICIT NONE
        INTEGER :: TYP
        
        INTEGER :: FILE_NUM
        INTEGER :: I
        CHARACTER(500) :: STR, STR2
        
        INTEGER :: EDGE_NUM
        INTEGER, POINTER, DIMENSION(:,:) :: EDGE
        
        IF (TYP==0) THEN
            EDGE_NUM = SURFACE_FLUID%SURFACE_EDGES_NUM
            EDGE => SURFACE_FLUID%SURFACE_EDGES
        END IF
        IF (TYP==1) THEN
            EDGE_NUM = SURFACE_PROPEL%SURFACE_EDGES_NUM
            EDGE => SURFACE_PROPEL%SURFACE_EDGES
        END IF
        IF (TYP==2) THEN
            EDGE_NUM = SURFACE_CASE%SURFACE_EDGES_NUM
            EDGE => SURFACE_CASE%SURFACE_EDGES
        END IF
        
        WRITE(STR, *), FILE_NUM
        WRITE(STR2, *), TYP
        STR = './output/surface/pointedges/edge_2d' // TRIM(ADJUSTL(STR)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        WRITE(*,*) TRIM(STR)
    
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        
        !$OMP DO ORDERED 
        DO I = 1, EDGE_NUM
            WRITE(21,*) EDGE(1,I)
            WRITE(21,*) EDGE(2,I)
        END DO
        !$OMP END DO
    
        CLOSE(21)
    END SUBROUTINE SAVINGEDGE

    SUBROUTINE SAVINGDATA_TECPLOT(FILE_NUM)
        IMPLICIT NONE
        INTEGER :: FILE_NUM
        INTEGER :: I
        CHARACTER(500) :: STR, STR2, STR3, STR4, STR5, STR6
        INTEGER :: L, L2, L3, L4, L5, L6

        REAL(8) :: R
    
        INTEGER :: J, J0
        LOGICAL :: B
        
        TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT
        
        STR = ''
        WRITE(STR, '(A)') './output/surface/tec_data_2d/surface_data_0000000.plt'
        L = LEN_TRIM(STR)
        
        STR2 = ''
        WRITE(STR2, '(I)') FILE_NUM
        STR2 = TRIM(ADJUSTL(STR2))
        L2 = LEN_TRIM(STR2)
        
        WRITE(STR(L-4-L2+1:L-4),'(A)') STR2(1:L2)
        
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        WRITE(*,*) TRIM(STR)
        
        WRITE(21,'(A)') 'TITLE="surface_data"'
        WRITE(21,'(A)') 'VARIABLES= "x", "y", "oninterface", "point_type", "force1", "force2", "b_rate", "fluid_dis", "propel_dis","case_dis", "related_fluid_edge", "propel_reledge", "ablation", "x2","y2", "initlength" '
        
        DO J=0,2
            
            IF(J==0) THEN
                SURFACE_CURRENT => SURFACE_FLUID
            ELSE IF(J==1) THEN
                SURFACE_CURRENT => SURFACE_PROPEL
            ELSE IF(J==2) THEN
                SURFACE_CURRENT => SURFACE_CASE
            END IF
            STR3 = ''
            WRITE(STR3, '(F)') SURFACE_TOTAL_TIME
            STR3 = TRIM(ADJUSTL(STR3))
            L3 = LEN_TRIM(STR3)
            
            STR4 = ''
            WRITE(STR4, '(I)') SURFACE_CURRENT%SURFACE_POINTS_NUM
            STR4 = TRIM(ADJUSTL(STR4))
            L4 = LEN_TRIM(STR4)
            
            STR5 = ''
            WRITE(STR5, '(I)') SURFACE_CURRENT%SURFACE_EDGES_NUM
            STR5 = TRIM(ADJUSTL(STR5))
            L5 = LEN_TRIM(STR5)
            
            STR6 = ', DATAPACKING = BLOCK, ZONETYPE = FELINESEG, VARLOCATION = ([1,2,4,5,6,11,12,14,15]=NODAL, [3,7,8,9,10,13,16]=CELLCENTERED)'
            STR6 = TRIM(ADJUSTL(STR6))
            L6 = LEN_TRIM(STR6)
            
            IF(J==0) THEN
                WRITE(21,'(A,A,A,A,A,A,A)') 'ZONE SolutionTime = ', STR3(1:L3),', T = "surface_fluid", NODES = ', STR4(1:L4),', ELEMENTS =  ', STR5(1:L5), STR6(1:L6)
            ELSE IF(J==1) THEN
                WRITE(21,'(A,A,A,A,A,A,A)') 'ZONE SolutionTime = ', STR3(1:L3),', T = "surface_propel", NODES = ', STR4(1:L4),', ELEMENTS =  ', STR5(1:L5), STR6(1:L6)
            ELSE IF(J==2) THEN
                WRITE(21,'(A,A,A,A,A,A,A)') 'ZONE SolutionTime = ', STR3(1:L3),', T = "surface_case", NODES = ', STR4(1:L4),', ELEMENTS =  ', STR5(1:L5), STR6(1:L6)
            END IF
            
            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
                WRITE(21,'(F)') SURFACE_CURRENT%SURFACE_POINTS(1,I) + SURFACE_CURRENT%POINT_DISPLACEMENT(1,I)
            END DO
            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
                WRITE(21,'(F)') SURFACE_CURRENT%SURFACE_POINTS(2,I) + SURFACE_CURRENT%POINT_DISPLACEMENT(2,I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
                WRITE(21,'(I)') SURFACE_CURRENT%EDGE_ONINTERFACE(I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
                WRITE(21,'(I)') SURFACE_CURRENT%POINT_TYPE(I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
                WRITE(21,'(F)') SURFACE_CURRENT%POINT_FORCE(1,I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
                WRITE(21,'(F)') SURFACE_CURRENT%POINT_FORCE(2,I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
                WRITE(21,'(F)') SURFACE_CURRENT%EDGE_B_RATE(I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
                J0 = SURFACE_CURRENT%EDGE_IMPACT_ZONE(0+1, I)
                IF(J0==0) THEN
                    WRITE(21,'(A)') '10000'
                ELSE
		    IF(J==0) THEN
                    CALL DISTANCE_EDGE_EDGE_TYPE(I,0,J0,0,0,  R,B)
		    ELSE
                    CALL DISTANCE_EDGE_EDGE_TYPE(I,J,J0,0,1,  R,B)
		    END IF
                    IF(.NOT. B) THEN
                       WRITE(21,'(A)') '10000'
                    ELSE
                        WRITE(21,'(F)') R
                    END IF
                END IF
            END DO

            DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
                J0 = SURFACE_CURRENT%EDGE_IMPACT_ZONE(1+1, I)
                IF(J0==0) THEN
                    WRITE(21,'(A)') '10000'
                ELSE
		    IF(J==1) THEN
                    CALL DISTANCE_EDGE_EDGE_TYPE(I,1,J0,1,0,  R,B)
		    ELSE
                    CALL DISTANCE_EDGE_EDGE_TYPE(I,J,J0,1,1,  R,B)
		    END IF

                    IF(.NOT. B) THEN
                       WRITE(21,'(A)') '10000'
                    ELSE
                        WRITE(21,'(F)') R
                    END IF
                END IF
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
                J0 = SURFACE_CURRENT%EDGE_IMPACT_ZONE(2+1, I)
                IF(J0==0) THEN
                    WRITE(21,'(A)') '10000'
                ELSE
		    IF(J==2) THEN
                    CALL DISTANCE_EDGE_EDGE_TYPE(I,2,J0,2,0,  R,B)
		    ELSE
                    CALL DISTANCE_EDGE_EDGE_TYPE(I,J,J0,2,1,  R,B)
		    END IF

                    IF(.NOT. B) THEN
                       WRITE(21,'(A)') '10000'
                    ELSE
                        WRITE(21,'(F)') R
                    END IF
                END IF
            END DO

            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		WRITE(21,'(I)') SURFACE_CURRENT%POINT_RELATEDEDGE(0+1,I)
            END DO

            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		WRITE(21,'(I)') SURFACE_CURRENT%POINT_RELATEDEDGE(1+1,I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
                WRITE(21,'(I)') SURFACE_CURRENT%EDGE_ABLATION_FLAG(I)
            END DO

            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		WRITE(21,'(F)') SURFACE_CURRENT%SURFACE_POINTS(1,I)
            END DO

            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		WRITE(21,'(F)') SURFACE_CURRENT%SURFACE_POINTS(2,I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
                WRITE(21,'(F)') SURFACE_CURRENT%SURFACE_INITIAL_EDGE_LENGTH(I)
            END DO
            
            DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
                WRITE(21,'(I,I)') SURFACE_CURRENT%SURFACE_EDGES(1,I), SURFACE_CURRENT%SURFACE_EDGES(2,I)
            END DO
            
        END DO
    
        CLOSE(21)
    
    END SUBROUTINE SAVINGDATA_TECPLOT

    SUBROUTINE SAVINGINTERFACE_TECPLOT(FILE_NUM)
        IMPLICIT NONE
        INTEGER :: FILE_NUM
        INTEGER :: I
        CHARACTER(500) :: STR, STR2, STR3, STR4, STR5, STR6
        INTEGER :: L, L2, L3, L4, L5, L6
    
        INTEGER :: J
        
        INTEGER :: POINT_NUM
        REAL(8), POINTER, DIMENSION(:,:) :: POINT
        INTEGER, POINTER, DIMENSION(:,:) :: POINT_LOC
        
        CALL UPDATE_INTERFACE_CLUSTER(0)
        CALL UPDATE_INTERFACE_CLUSTER(1)
        
        STR = ''
        WRITE(STR, '(A)') './output/surface/tec_interface_2d/surface_interface_0000000.plt'
        L = LEN_TRIM(STR)
        
        STR2 = ''
        WRITE(STR2, '(I)') FILE_NUM
        STR2 = TRIM(ADJUSTL(STR2))
        L2 = LEN_TRIM(STR2)
        
        WRITE(STR(L-4-L2+1:L-4),'(A)') STR2(1:L2)
        
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        WRITE(*,*) TRIM(STR)
        
        WRITE(21,'(A)') 'TITLE="surface_interface"'
        WRITE(21,'(A)') 'VARIABLES= "x", "y", "x2", "y2", "force1", "force2", "pressure", "brate"'
        
        DO J=0,1
        
            IF (J==0) THEN
                POINT_NUM = INTERFACE_FLUID_POINTS_NUM
                POINT => INTERFACE_FLUID_POINTS
                POINT_LOC => INTERFACE_FLUID_POINTS_LOC
            ELSE IF(J==1) THEN
                POINT_NUM = INTERFACE_STRUCT_POINTS_NUM
                POINT => INTERFACE_STRUCT_POINTS
                POINT_LOC => INTERFACE_STRUCT_POINTS_LOC
            END IF
            
            STR3 = ''
            WRITE(STR3, '(F)') SURFACE_TOTAL_TIME
            STR3 = TRIM(ADJUSTL(STR3))
            L3 = LEN_TRIM(STR3)
            
            STR4 = ''
            WRITE(STR4, '(I)') POINT_NUM
            STR4 = TRIM(ADJUSTL(STR4))
            L4 = LEN_TRIM(STR4)
            
            STR5 = ''
            WRITE(STR5, '(I)') POINT_NUM - 1
            STR5 = TRIM(ADJUSTL(STR5))
            L5 = LEN_TRIM(STR5)
            
            STR6 = ', DATAPACKING = BLOCK, ZONETYPE = FELINESEG, VARLOCATION = ([1,2,3,4,5,6]=NODAL, [7,8]=CELLCENTERED)'
            STR6 = TRIM(ADJUSTL(STR6))
            L6 = LEN_TRIM(STR6)
            
            IF(J==0) THEN
                WRITE(21,'(A,A,A,A,A,A,A)') 'ZONE SolutionTime = ', STR3(1:L3),', T = "surface_fluid", NODES = ', STR4(1:L4),', ELEMENTS =  ', STR5(1:L5), STR6(1:L6)
            ELSE IF(J==1) THEN
                WRITE(21,'(A,A,A,A,A,A,A)') 'ZONE SolutionTime = ', STR3(1:L3),', T = "surface_propel", NODES = ', STR4(1:L4),', ELEMENTS =  ', STR5(1:L5), STR6(1:L6)
            ELSE IF(J==2) THEN
                WRITE(21,'(A,A,A,A,A,A,A)') 'ZONE SolutionTime = ', STR3(1:L3),', T = "surface_case", NODES = ', STR4(1:L4),', ELEMENTS =  ', STR5(1:L5), STR6(1:L6)
            END IF
            
!            DO I = 1, POINT_NUM
!                WRITE(21,'(F)') POINT(1,I)
!            END DO
!            DO I = 1, POINT_NUM
!                WRITE(21,'(F)') POINT(2,I)
!            END DO

            DO I = 1, POINT_NUM
                IF(POINT_LOC(2,I)==0) THEN
                    WRITE(21,'(F)') POINT(1,I)+SURFACE_FLUID%POINT_DISPLACEMENT(1,POINT_LOC(1,I))
                ELSE IF(POINT_LOC(2,I)==1) THEN
                    WRITE(21,'(F)') POINT(1,I)+SURFACE_PROPEL%POINT_DISPLACEMENT(1,POINT_LOC(1,I))
                ELSE
                    WRITE(21,'(F)') POINT(1,I)+SURFACE_CASE%POINT_DISPLACEMENT(1,POINT_LOC(1,I))
                END IF
            END DO

            DO I = 1, POINT_NUM
                IF(POINT_LOC(2,I)==0) THEN
                    WRITE(21,'(F)') POINT(2,I)+SURFACE_FLUID%POINT_DISPLACEMENT(2,POINT_LOC(1,I))
                ELSE IF(POINT_LOC(2,I)==1) THEN
                    WRITE(21,'(F)') POINT(2,I)+SURFACE_PROPEL%POINT_DISPLACEMENT(2,POINT_LOC(1,I))
                ELSE
                    WRITE(21,'(F)') POINT(2,I)+SURFACE_CASE%POINT_DISPLACEMENT(2,POINT_LOC(1,I))
                END IF
            END DO

            DO I = 1, POINT_NUM
                WRITE(21,'(F)') POINT(1,I)
            END DO

            DO I = 1, POINT_NUM
                WRITE(21,'(F)') POINT(2,I)
            END DO
            
            DO I = 1, POINT_NUM
                IF(POINT_LOC(2,I)==0) THEN
                    WRITE(21,'(F)') SURFACE_FLUID%POINT_FORCE(1,POINT_LOC(1,I))
                ELSE IF(POINT_LOC(2,I)==1) THEN
                    WRITE(21,'(F)') SURFACE_PROPEL%POINT_FORCE(1,POINT_LOC(1,I))
                ELSE
                    WRITE(21,'(F)') SURFACE_CASE%POINT_FORCE(1,POINT_LOC(1,I))
                END IF
            END DO
            
            DO I = 1, POINT_NUM
                IF(POINT_LOC(2,I)==0) THEN
                    WRITE(21,'(F)') SURFACE_FLUID%POINT_FORCE(2,POINT_LOC(1,I))
                ELSE IF(POINT_LOC(2,I)==1) THEN
                    WRITE(21,'(F)') SURFACE_PROPEL%POINT_FORCE(2,POINT_LOC(1,I))
                ELSE
                    WRITE(21,'(F)') SURFACE_CASE%POINT_FORCE(2,POINT_LOC(1,I))
                END IF
            END DO
            
            DO I = 1, POINT_NUM-1
                IF(POINT_LOC(2,I)==0) THEN
                    WRITE(21,'(F)') SURFACE_FLUID%EDGE_PRESSURE(SURFACE_FLUID%POINT_EDGE_CONNECTION(2,POINT_LOC(1,I)))
                ELSE IF(POINT_LOC(2,I)==1) THEN
                    WRITE(21,'(F)') SURFACE_PROPEL%EDGE_PRESSURE(SURFACE_PROPEL%POINT_EDGE_CONNECTION(1,POINT_LOC(1,I)))
                ELSE
                    WRITE(21,'(F)') SURFACE_CASE%EDGE_PRESSURE(SURFACE_CASE%POINT_EDGE_CONNECTION(1,POINT_LOC(1,I)))
                END IF
            END DO
            
            DO I = 1, POINT_NUM-1
                IF(POINT_LOC(2,I)==0) THEN
                    WRITE(21,'(F)') SURFACE_FLUID%EDGE_B_RATE(SURFACE_FLUID%POINT_EDGE_CONNECTION(2,POINT_LOC(1,I)))
                ELSE IF(POINT_LOC(2,I)==1) THEN
                    WRITE(21,'(F)') SURFACE_PROPEL%EDGE_B_RATE(SURFACE_PROPEL%POINT_EDGE_CONNECTION(1,POINT_LOC(1,I)))
                ELSE
                    WRITE(21,'(F)') SURFACE_CASE%EDGE_B_RATE(SURFACE_CASE%POINT_EDGE_CONNECTION(1,POINT_LOC(1,I)))
                END IF
            END DO
            
            DO I = 1, POINT_NUM - 1
                WRITE(21,'(I,I)') I, I+1
            END DO
            
        END DO
    
        CLOSE(21)
    
    END SUBROUTINE SAVINGINTERFACE_TECPLOT
    
    SUBROUTINE SAVINGBRATE(FILE_NUM, TYP)
        IMPLICIT NONE
        INTEGER :: TYP
        
        INTEGER :: FILE_NUM
        INTEGER :: I
        CHARACTER(500) :: STR
        
        INTEGER :: EDGE_NUM
        REAL(8), POINTER, DIMENSION(:) :: B_RATE
        
        TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT
        
        IF (TYP==0) THEN
            SURFACE_CURRENT => SURFACE_FLUID
        ELSE IF (TYP==1) THEN
            SURFACE_CURRENT => SURFACE_PROPEL
        ELSE IF (TYP==2) THEN
            SURFACE_CURRENT => SURFACE_CASE
        END IF
        
        EDGE_NUM = SURFACE_CURRENT%SURFACE_EDGES_NUM
        B_RATE => SURFACE_CURRENT%EDGE_B_RATE
        
        WRITE(STR, *), FILE_NUM
        STR = './output/surface/brates/brate_2d' // TRIM(ADJUSTL(STR)) // '.txt'
        WRITE(*,*) TRIM(STR)
    
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        
        !$OMP DO ORDERED 
        DO I = 1, EDGE_NUM
            WRITE(21,*) B_RATE(I)
        END DO
    	!$OMP END DO 
    
        CLOSE(21)
    END SUBROUTINE SAVINGBRATE    

    SUBROUTINE SAVINGDISPLACEMENT(FILE_NUM, TYP)
	IMPLICIT NONE
        INTEGER :: TYP
        
        INTEGER :: FILE_NUM
        INTEGER :: I
        CHARACTER(500) :: STR, STR2
        
        INTEGER :: POINT_NUM
        REAL(8), POINTER, DIMENSION(:,:) :: DISPLACEMENT
        
        TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT
        
        IF (TYP==0) THEN
            SURFACE_CURRENT => SURFACE_FLUID
        ELSE IF (TYP==1) THEN
            SURFACE_CURRENT => SURFACE_PROPEL
        ELSE IF (TYP==2) THEN
            SURFACE_CURRENT => SURFACE_CASE
        END IF
        
            POINT_NUM = SURFACE_CURRENT%SURFACE_POINTS_NUM
            DISPLACEMENT => SURFACE_CURRENT%POINT_VELOCITY

        WRITE(STR, *), FILE_NUM
        WRITE(STR2, *), TYP
        STR = './output/surface/displacements/displacement_2d' // TRIM(ADJUSTL(STR)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        WRITE(*,*) TRIM(STR)
    
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        
        !$OMP DO ORDERED 
        DO I = 1, POINT_NUM
            WRITE(21,*) DISPLACEMENT(1,I)
            WRITE(21,*) DISPLACEMENT(2,I)
        END DO
    	!$OMP END DO
    
        CLOSE(21)
    END SUBROUTINE SAVINGDISPLACEMENT

    SUBROUTINE SAVINGPRESSURE(FILE_NUM, TYP)
	IMPLICIT NONE
        INTEGER :: TYP
        
        INTEGER :: FILE_NUM
        INTEGER :: I
        CHARACTER(500) :: STR, STR2
        
        TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT
        
        IF (TYP==0) THEN
            SURFACE_CURRENT => SURFACE_FLUID
        ELSE IF (TYP==1) THEN
            SURFACE_CURRENT => SURFACE_PROPEL
        ELSE IF (TYP==2) THEN
            SURFACE_CURRENT => SURFACE_CASE
        END IF

        WRITE(STR, *), FILE_NUM
        WRITE(STR2, *), TYP
        STR = './output/surface/pressures/pressure_2d' // TRIM(ADJUSTL(STR)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        WRITE(*,*) TRIM(STR)
    
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        
        IF(TYP==0) THEN
            !$OMP DO ORDERED 
            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
                WRITE(21,*) SURFACE_CURRENT%EDGE_PRESSURE(I)
            END DO
    	    !$OMP END DO
        ELSE
            !$OMP DO ORDERED 
            DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
                WRITE(21,'(F,F)') SURFACE_CURRENT%POINT_FORCE(1,I), SURFACE_CURRENT%POINT_FORCE(2,I)
            END DO
    	    !$OMP END DO
        END IF
    
        CLOSE(21)
    END SUBROUTINE SAVINGPRESSURE

    SUBROUTINE COMPUTE_SURFACE_AREA(TYP)
	IMPLICIT NONE
        INTEGER :: TYP
        
        INTEGER :: I, I1, I2
        REAL(8) :: S, V(2), R, Y1, Y2
        CHARACTER(500) :: STR, STR2
        
        TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT
        
        IF (TYP==0) THEN
            SURFACE_CURRENT => SURFACE_FLUID
        ELSE IF (TYP==1) THEN
            SURFACE_CURRENT => SURFACE_PROPEL
        ELSE IF (TYP==2) THEN
            SURFACE_CURRENT => SURFACE_CASE
        END IF
        
	SURFACE_AREA_ITER = SURFACE_AREA_ITER + 1
        S = 0.
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            IF(SURFACE_CURRENT%EDGE_ONINTERFACE(I)==1-TYP) THEN
                I1 = SURFACE_CURRENT%SURFACE_EDGES(1,I)
                I2 = SURFACE_CURRENT%SURFACE_EDGES(2,I)
                Y1 = SURFACE_CURRENT%SURFACE_POINTS(2,I1)
		Y2 = SURFACE_CURRENT%SURFACE_POINTS(2,I2)

                V = SURFACE_CURRENT%SURFACE_POINTS(:,I2) - SURFACE_CURRENT%SURFACE_POINTS(:,I1)
                
                R = SQRT(DOT_PRODUCT(V,V))
                
                S = S + PI*R*(Y1+Y2)
            END IF
        END DO
        SURFACE_AREA_ARRAY(SURFACE_AREA_ITER) = S
        
        WRITE(STR2, *), TYP
        STR = './output/surface/surfacearea2d_' // TRIM(ADJUSTL(STR2)) // '.txt'
        WRITE(*,*) TRIM(STR)
        
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        
        DO I = 1, SURFACE_AREA_ITER
            WRITE(21,'(F,I,I,I)') SURFACE_AREA_ARRAY(I), SURFACE_FLAG_ARRAY(I,1), SURFACE_FLAG_ARRAY(I,2), SURFACE_PROPEL%SURFACE_PATCHES_NUM
        END DO
    
        CLOSE(21)
        
    END SUBROUTINE COMPUTE_SURFACE_AREA

   SUBROUTINE SAVING_SURFACE(FILE_NUM, TYP)
	IMPLICIT NONE

        INTEGER :: TYP
        
        INTEGER :: FILE_NUM
        INTEGER :: I
        CHARACTER(500) :: STR, str1,STR2
        
	TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT
        
         
        IF (TYP==0) THEN
            SURFACE_CURRENT => SURFACE_FLUID
        END IF
        IF (TYP==1) THEN
            SURFACE_CURRENT => SURFACE_PROPEL
        END IF
        IF (TYP==2) THEN
            SURFACE_CURRENT => SURFACE_CASE
        END IF
        
        WRITE(STR1, *), FILE_NUM
        WRITE(STR2, *), TYP

        STR = './output/surface/restart_2d/mesh_num_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        WRITE(21,*) SURFACE_CURRENT%MESH_SIZE
        WRITE(21,*) SURFACE_CURRENT%MESH_SIZE_MAX
        WRITE(21,*) SURFACE_CURRENT%SURFACE_POINTS_NUM
        WRITE(21,*) SURFACE_CURRENT%SURFACE_EDGES_NUM
	WRITE(21,*) SURFACE_CURRENT%SURFACE_PATCHES_NUM
	WRITE(21,*) INTERFACE_FLUID_POINTS_NUM
	WRITE(21,*) INTERFACE_STRUCT_POINTS_NUM
        CLOSE(21)


        STR = './output/surface/restart_2d/point_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
            WRITE(21,*) SURFACE_CURRENT%SURFACE_POINTS(1,I)
            WRITE(21,*) SURFACE_CURRENT%SURFACE_POINTS(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/point_velocity_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
            WRITE(21,*) SURFACE_CURRENT%POINT_VELOCITY(1,I)
            WRITE(21,*) SURFACE_CURRENT%POINT_VELOCITY(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/point_displacement_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
            WRITE(21,*) SURFACE_CURRENT%POINT_DISPLACEMENT(1,I)
            WRITE(21,*) SURFACE_CURRENT%POINT_DISPLACEMENT(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/point_force_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
            WRITE(21,*) SURFACE_CURRENT%POINT_FORCE(1,I)
            WRITE(21,*) SURFACE_CURRENT%POINT_FORCE(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)


        STR = './output/surface/restart_2d/point_edge_connection_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
            WRITE(21,*) SURFACE_CURRENT%POINT_EDGE_CONNECTION(1,I)
            WRITE(21,*) SURFACE_CURRENT%POINT_EDGE_CONNECTION(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/point_type_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
            WRITE(21,*) SURFACE_CURRENT%POINT_TYPE(I)
        END DO
        !$OMP END DO
        CLOSE(21)


        STR = './output/surface/restart_2d/point_relatedpt_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
            WRITE(21,*) SURFACE_CURRENT%POINT_RELATEDPT(1,I)
            WRITE(21,*) SURFACE_CURRENT%POINT_RELATEDPT(2,I)
            WRITE(21,*) SURFACE_CURRENT%POINT_RELATEDPT(3,I)
        END DO
        !$OMP END DO
        CLOSE(21)


        STR = './output/surface/restart_2d/point_relatededge_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
            WRITE(21,*) SURFACE_CURRENT%POINT_RELATEDEDGE(1,I)
            WRITE(21,*) SURFACE_CURRENT%POINT_RELATEDEDGE(2,I)
            WRITE(21,*) SURFACE_CURRENT%POINT_RELATEDEDGE(3,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/edge_b_rate_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            WRITE(21,*) SURFACE_CURRENT%EDGE_B_RATE(I)
        END DO
        !$OMP END DO
        CLOSE(21)

         STR = './output/surface/restart_2d/initial_edge_length_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            WRITE(21,*) SURFACE_CURRENT%SURFACE_INITIAL_EDGE_LENGTH(I)
        END DO
        !$OMP END DO
        CLOSE(21)


        STR = './output/surface/restart_2d/edge_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            WRITE(21,*) SURFACE_CURRENT%SURFACE_EDGES(1,I)
            WRITE(21,*) SURFACE_CURRENT%SURFACE_EDGES(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/edge_location_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            WRITE(21,*) SURFACE_CURRENT%EDGE_LOCATION(I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/edge_oninterface_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            WRITE(21,*) SURFACE_CURRENT%EDGE_ONINTERFACE(I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/edge_pressure_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            WRITE(21,*) SURFACE_CURRENT%EDGE_PRESSURE(I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/edge_impact_zone_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            WRITE(21,*) SURFACE_CURRENT%EDGE_IMPACT_ZONE(1,I)
            WRITE(21,*) SURFACE_CURRENT%EDGE_IMPACT_ZONE(2,I)
            WRITE(21,*) SURFACE_CURRENT%EDGE_IMPACT_ZONE(3,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/edge_ablation_flag_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
            WRITE(21,*) SURFACE_CURRENT%EDGE_ABLATION_FLAG(I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/topchange_typ_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, SURFACE_CURRENT%SURFACE_PATCHES_NUM
            WRITE(21,*) SURFACE_CURRENT%SURFACE_PATCHES_TOPCHANGE_TYP(I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/interface_fluid_points_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, INTERFACE_FLUID_POINTS_NUM
            WRITE(21,*) INTERFACE_FLUID_POINTS(1,I)
    	    WRITE(21,*) INTERFACE_FLUID_POINTS(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/interface_fluid_points_loc_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, INTERFACE_FLUID_POINTS_NUM
            WRITE(21,*) INTERFACE_FLUID_POINTS_LOC(1,I)
    	    WRITE(21,*) INTERFACE_FLUID_POINTS_LOC(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)


        STR = './output/surface/restart_2d/interface_struct_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, INTERFACE_STRUCT_POINTS_NUM
            WRITE(21,*) INTERFACE_STRUCT_POINTS(1,I)
    	    WRITE(21,*) INTERFACE_STRUCT_POINTS(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)

        STR = './output/surface/restart_2d/interface_struct_points_loc_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
        OPEN(UNIT=21, FILE = STR, ACTION = "WRITE", STATUS = "REPLACE")
        !$OMP DO ORDERED 
        DO I = 1, INTERFACE_STRUCT_POINTS_NUM
            WRITE(21,*) INTERFACE_STRUCT_POINTS_LOC(1,I)
    	    WRITE(21,*) INTERFACE_STRUCT_POINTS_LOC(2,I)
        END DO
        !$OMP END DO
        CLOSE(21)

    END SUBROUTINE SAVING_SURFACE

    SUBROUTINE READ_RESTART(TYP)

    INTEGER :: TYP, IO, I
    CHARACTER(500) :: STR, STR1, STR2
     TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT

	WRITE(STR1,*) RESTART_ITERATION
        WRITE(STR2,*) TYP
	
	IF(TYP==0) THEN
	SURFACE_CURRENT => SURFACE_FLUID
	ELSEIF(TYP==1) THEN
	SURFACE_CURRENT => SURFACE_PROPEL
	ELSE
	SURFACE_CURRENT => SURFACE_CASE
	END IF	 

	WRITE(*,*) 'START_READING_SURFACE_INFORMATION'

        STR = './output/surface/restart_2d/mesh_num_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
        OPEN(NEWUNIT=IO,FILE= STR)
    	    READ(IO,*) SURFACE_CURRENT%MESH_SIZE
    	    READ(IO,*) SURFACE_CURRENT%MESH_SIZE_MAX
    	    READ(IO,*) SURFACE_CURRENT%SURFACE_POINTS_NUM
    	    READ(IO,*) SURFACE_CURRENT%SURFACE_EDGES_NUM
    	    READ(IO,*) SURFACE_CURRENT%SURFACE_PATCHES_NUM
    	    READ(IO,*) INTERFACE_FLUID_POINTS_NUM
    	    READ(IO,*) INTERFACE_STRUCT_POINTS_NUM		
        CLOSE(IO)

        STR = './output/surface/restart_2d/point_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%SURFACE_POINTS(2,SURFACE_CURRENT%SURFACE_POINTS_NUM))
        OPEN(NEWUNIT=IO,FILE= STR)
	    DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		READ(IO,*) SURFACE_CURRENT%SURFACE_POINTS(1,I)
	        READ(IO,*) SURFACE_CURRENT%SURFACE_POINTS(2,I)
	    END DO		
        CLOSE(IO)

        STR = './output/surface/restart_2d/point_velocity_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%POINT_VELOCITY(2,SURFACE_CURRENT%SURFACE_POINTS_NUM))
        OPEN(NEWUNIT=IO,FILE= STR)
	    DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		READ(IO,*) SURFACE_CURRENT%POINT_VELOCITY(1,I)
	        READ(IO,*) SURFACE_CURRENT%POINT_VELOCITY(2,I)
	    END DO		
        CLOSE(IO)

        STR = './output/surface/restart_2d/point_displacement_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%POINT_DISPLACEMENT(2,SURFACE_CURRENT%SURFACE_POINTS_NUM))
        OPEN(NEWUNIT=IO,FILE= STR)
	    DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		READ(IO,*) SURFACE_CURRENT%POINT_DISPLACEMENT(1,I)
	        READ(IO,*) SURFACE_CURRENT%POINT_DISPLACEMENT(2,I)
	    END DO		
        CLOSE(IO)

        STR = './output/surface/restart_2d/point_force_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%POINT_FORCE(2,SURFACE_CURRENT%SURFACE_POINTS_NUM))
        OPEN(NEWUNIT=IO,FILE= STR)
	    DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		READ(IO,*) SURFACE_CURRENT%POINT_FORCE(1,I)
	        READ(IO,*) SURFACE_CURRENT%POINT_FORCE(2,I)
	    END DO		
        CLOSE(IO)

        STR = './output/surface/restart_2d/point_edge_connection_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%POINT_EDGE_CONNECTION(2,SURFACE_CURRENT%SURFACE_POINTS_NUM))
        OPEN(NEWUNIT=IO,FILE= STR)
	    DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
		READ(IO,*) SURFACE_CURRENT%POINT_EDGE_CONNECTION(1,I)
	        READ(IO,*) SURFACE_CURRENT%POINT_EDGE_CONNECTION(2,I)
	    END DO		
        CLOSE(IO)

	STR = './output/surface/restart_2d/point_type_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%POINT_TYPE(SURFACE_CURRENT%SURFACE_POINTS_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
	    READ(21,*) SURFACE_CURRENT%POINT_TYPE(I)
	END DO
	!$OMP END DO
	CLOSE(21)


	STR = './output/surface/restart_2d/point_relatedpt_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%POINT_RELATEDPT(3,SURFACE_CURRENT%SURFACE_POINTS_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
	    READ(21,*) SURFACE_CURRENT%POINT_RELATEDPT(1,I)
	    READ(21,*) SURFACE_CURRENT%POINT_RELATEDPT(2,I)
	    READ(21,*) SURFACE_CURRENT%POINT_RELATEDPT(3,I)
	END DO
	!$OMP END DO
	CLOSE(21)


	STR = './output/surface/restart_2d/point_relatededge_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%POINT_RELATEDEDGE(3,SURFACE_CURRENT%SURFACE_POINTS_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_POINTS_NUM
	    READ(21,*) SURFACE_CURRENT%POINT_RELATEDEDGE(1,I)
	    READ(21,*) SURFACE_CURRENT%POINT_RELATEDEDGE(2,I)
	    READ(21,*) SURFACE_CURRENT%POINT_RELATEDEDGE(3,I)
	END DO
	!$OMP END DO
	CLOSE(21)

	STR = './output/surface/restart_2d/edge_b_rate_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%EDGE_B_RATE(SURFACE_CURRENT%SURFACE_EDGES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
	    READ(21,*) SURFACE_CURRENT%EDGE_B_RATE(I)
	END DO
	!$OMP END DO
	CLOSE(21)

	STR = './output/surface/restart_2d/initial_edge_length_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%SURFACE_INITIAL_EDGE_LENGTH(SURFACE_CURRENT%SURFACE_EDGES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
	    READ(21,*) SURFACE_CURRENT%SURFACE_INITIAL_EDGE_LENGTH(I)
	END DO
	!$OMP END DO
	CLOSE(21)


	STR = './output/surface/restart_2d/edge_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%SURFACE_EDGES(2,SURFACE_CURRENT%SURFACE_EDGES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
	    READ(21,*) SURFACE_CURRENT%SURFACE_EDGES(1,I)
	    READ(21,*) SURFACE_CURRENT%SURFACE_EDGES(2,I)
	END DO
	!$OMP END DO
	CLOSE(21)

	STR = './output/surface/restart_2d/edge_location_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%EDGE_LOCATION(SURFACE_CURRENT%SURFACE_EDGES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
	    READ(21,*) SURFACE_CURRENT%EDGE_LOCATION(I)
	END DO
	!$OMP END DO
	CLOSE(21)

	STR = './output/surface/restart_2d/edge_oninterface_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%EDGE_ONINTERFACE(SURFACE_CURRENT%SURFACE_EDGES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
	    READ(21,*) SURFACE_CURRENT%EDGE_ONINTERFACE(I)
	END DO
	!$OMP END DO
	CLOSE(21)

	STR = './output/surface/restart_2d/edge_pressure_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%EDGE_PRESSURE(SURFACE_CURRENT%SURFACE_EDGES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
	    READ(21,*) SURFACE_CURRENT%EDGE_PRESSURE(I)
	END DO
	!$OMP END DO
	CLOSE(21)

	STR = './output/surface/restart_2d/edge_impact_zone_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%EDGE_IMPACT_ZONE(3,SURFACE_CURRENT%SURFACE_EDGES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
	    READ(21,*) SURFACE_CURRENT%EDGE_IMPACT_ZONE(1,I)
	    READ(21,*) SURFACE_CURRENT%EDGE_IMPACT_ZONE(2,I)
	    READ(21,*) SURFACE_CURRENT%EDGE_IMPACT_ZONE(3,I)
	END DO
	!$OMP END DO
	CLOSE(21)

	STR = './output/surface/restart_2d/edge_ablation_flag_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%EDGE_ABLATION_FLAG(SURFACE_CURRENT%SURFACE_EDGES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_EDGES_NUM
	    READ(21,*) SURFACE_CURRENT%EDGE_ABLATION_FLAG(I)
	END DO
	!$OMP END DO
	CLOSE(21)

	STR = './output/surface/restart_2d/topchange_typ_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
	STR = trim(STR)
	ALLOCATE(SURFACE_CURRENT%SURFACE_PATCHES_TOPCHANGE_TYP(SURFACE_CURRENT%SURFACE_PATCHES_NUM))
	OPEN(UNIT=21, FILE = STR)
	!$OMP DO ORDERED 
	DO I = 1, SURFACE_CURRENT%SURFACE_PATCHES_NUM
	    READ(21,*) SURFACE_CURRENT%SURFACE_PATCHES_TOPCHANGE_TYP(I)
	END DO
	!$OMP END DO
	CLOSE(21)

	IF(TYP == 0) THEN
		STR = './output/surface/restart_2d/interface_fluid_points_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
		STR = trim(STR)
		ALLOCATE(INTERFACE_FLUID_POINTS(2,INTERFACE_FLUID_POINTS_NUM))
		OPEN(UNIT=21, FILE = STR)
		!$OMP DO ORDERED 
		DO I = 1, INTERFACE_FLUID_POINTS_NUM
		    READ(21,*) INTERFACE_FLUID_POINTS(1,I)
		    READ(21,*) INTERFACE_FLUID_POINTS(2,I)
		END DO
		!$OMP END DO
		CLOSE(21)

		STR = './output/surface/restart_2d/interface_fluid_points_loc_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
		STR = trim(STR)
		ALLOCATE(INTERFACE_FLUID_POINTS_LOC(2,INTERFACE_FLUID_POINTS_NUM))
		OPEN(UNIT=21, FILE = STR)
		!$OMP DO ORDERED 
		DO I = 1, INTERFACE_FLUID_POINTS_NUM
		    READ(21,*) INTERFACE_FLUID_POINTS_LOC(1,I)
		    READ(21,*) INTERFACE_FLUID_POINTS_LOC(2,I)
		END DO
		!$OMP END DO
		CLOSE(21)


		STR = './output/surface/restart_2d/interface_struct_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
		STR = trim(STR)
		ALLOCATE(INTERFACE_STRUCT_POINTS(2,INTERFACE_STRUCT_POINTS_NUM))
		OPEN(UNIT=21, FILE = STR)
		!$OMP DO ORDERED 
		DO I = 1, INTERFACE_STRUCT_POINTS_NUM
		    READ(21,*) INTERFACE_STRUCT_POINTS(1,I)
		    READ(21,*) INTERFACE_STRUCT_POINTS(2,I)
		END DO
		!$OMP END DO
		CLOSE(21)

		STR = './output/surface/restart_2d/interface_struct_points_loc_2d' // TRIM(ADJUSTL(STR1)) // '_' // TRIM(ADJUSTL(STR2)) // '.txt'
		STR = trim(STR)
		ALLOCATE(INTERFACE_STRUCT_POINTS_LOC(2,INTERFACE_STRUCT_POINTS_NUM))
		OPEN(UNIT=21, FILE = STR)
		!$OMP DO ORDERED 
		DO I = 1, INTERFACE_STRUCT_POINTS_NUM
		    READ(21,*) INTERFACE_STRUCT_POINTS_LOC(1,I)
		    READ(21,*) INTERFACE_STRUCT_POINTS_LOC(2,I)
		END DO
		!$OMP END DO
		CLOSE(21)
	END IF

	WRITE(*,*) 'END_READING_SURFACE_INFORMATION'

END SUBROUTINE READ_RESTART

END MODULE
    
    

!! THIS SUBROUTINE IS MAIN SUBROUTINE OF SURFACE_2D MODULE
    

SUBROUTINE SURFACE_2D(TYPE1, TYPE2, TIMESTEP, &
           F_POINT_NUM,F_POINT,F_EDGE_NUM,F_EDGE,F_LOC,F_B_RATE,F_PRESSURE,F_BCFLAG, &
           P_POINT_NUM,P_POINT,P_EDGE_NUM,P_EDGE,P_LOC,P_VELOCITY,P_DISPLACEMENT,P_FORCE,P_EDGE_LENGTH,P_CORNER_INDEX, &
           C_POINT_NUM,C_POINT,C_EDGE_NUM,C_EDGE,C_VELOCITY,C_DISPLACEMENT,C_FORCE, &
           FLAG, FLAG_ARRAY, FILENUM, PATCHNUM)

    USE SURFACE_MODULE_2D
    USE PROPA_RECONST_REINITIAL_2D
    USE SURFACE_EXTERNAL_COMMON_2D
    USE REMESHING_2D
    IMPLICIT NONE

    !TYPE(SURFACE_TYPE), POINTER :: SURFACE_CURRENT
    
    INTEGER, OPTIONAL :: TYPE1
    INTEGER, OPTIONAL :: TYPE2

    REAL(8), OPTIONAL :: TIMESTEP
        
    INTEGER, OPTIONAL :: F_POINT_NUM
    REAL(8), ALLOCATABLE, OPTIONAL :: F_POINT(:,:)
    INTEGER, OPTIONAL :: F_EDGE_NUM
    INTEGER, ALLOCATABLE, OPTIONAL :: F_EDGE(:,:)
    INTEGER, ALLOCATABLE, OPTIONAL :: F_LOC(:)
    REAL(8), OPTIONAL :: F_B_RATE(:)
    REAL(8), OPTIONAL :: F_PRESSURE(:)
    INTEGER, ALLOCATABLE, OPTIONAL :: F_BCFLAG(:)
        
    INTEGER, OPTIONAL :: P_POINT_NUM
    REAL(8), ALLOCATABLE, OPTIONAL :: P_POINT(:,:)
    INTEGER, OPTIONAL :: P_EDGE_NUM
    INTEGER, ALLOCATABLE, OPTIONAL :: P_EDGE(:,:)
    INTEGER, ALLOCATABLE, OPTIONAL :: P_LOC(:)
    REAL(8), OPTIONAL :: P_VELOCITY(:,:)
    REAL(8), OPTIONAL :: P_DISPLACEMENT(:,:)
    REAL(8), OPTIONAL :: P_FORCE(:,:)
    
    REAL(8), ALLOCATABLE, OPTIONAL :: P_EDGE_LENGTH(:)
    INTEGER, ALLOCATABLE, OPTIONAL :: P_CORNER_INDEX(:)
        
    INTEGER, OPTIONAL :: C_POINT_NUM
    REAL(8), ALLOCATABLE, OPTIONAL :: C_POINT(:,:)
    INTEGER, OPTIONAL :: C_EDGE_NUM
    INTEGER, ALLOCATABLE, OPTIONAL :: C_EDGE(:,:)
    REAL(8), OPTIONAL :: C_VELOCITY(:,:)
    REAL(8), OPTIONAL :: C_DISPLACEMENT(:,:)
    REAL(8), OPTIONAL :: C_FORCE(:,:)
        
    LOGICAL, OPTIONAL :: FLAG
    INTEGER, ALLOCATABLE, OPTIONAL :: FLAG_ARRAY(:)
        
    INTEGER, OPTIONAL :: FILENUM
    INTEGER, OPTIONAL :: PATCHNUM
        
    REAL(8) :: MINX, MAXX, MINY, MAXY
    
    INTEGER, ALLOCATABLE :: TEMP_ONINTERFACE(:)
    INTEGER, ALLOCATABLE :: TEMP_CONNECTION(:,:)
!    INTEGER, ALLOCATABLE :: NUM_PATCH(:)
    !LOGICAL :: TEMPFLAG
    
    INTEGER :: I
    INTEGER :: TYP1_POINT_NUM
    REAL(8), ALLOCATABLE :: TYP1_POINT(:,:)
    INTEGER, ALLOCATABLE :: TYP1_POINTLOC(:,:)
    REAL(8), ALLOCATABLE :: DATA_DISPLACEMENT(:,:)

    !REAL(8), ALLOCATABLE :: TEMPPOINT(:,:)
    !INTEGER, ALLOCATABLE :: TEMPEDGE(:,:)
    
    LOGICAL :: TEMPFLAG
        
    !REAL(8) :: A(2), B(2)
    IF(TYPE1==-2) THEN
        
        WRITE(*,*) 'SURFACE FINILAZATION STARTED'

	CALL FIND_PATCH_NUM(0)
        CALL FIND_PATCH_NUM(1)
	CALL FIND_PATCH_NUM(2)
        CALL RESET_SURFACE(0, .FALSE., 0, 0) !SURFACE_FLUID%SURFACE_PATCHES_NUM)
        CALL RESET_SURFACE(1, .FALSE., 0, 0) !SURFACE_PROPEL%SURFACE_PATCHES_NUM)
        CALL RESET_SURFACE(2, .FALSE., 0, 0) !SURFACE_CASE%SURFACE_PATCHES_NUM)
        
        IF(ALLOCATED(INTERFACE_FLUID_POINTS)) THEN
            DEALLOCATE(INTERFACE_FLUID_POINTS)
        END IF
        IF(ALLOCATED(INTERFACE_FLUID_POINTS_LOC)) THEN
            DEALLOCATE(INTERFACE_FLUID_POINTS_LOC)
        END IF
        
        IF(ALLOCATED(INTERFACE_STRUCT_POINTS)) THEN
            DEALLOCATE(INTERFACE_STRUCT_POINTS)
        END IF
        IF(ALLOCATED(INTERFACE_STRUCT_POINTS_LOC)) THEN
            DEALLOCATE(INTERFACE_STRUCT_POINTS_LOC)
        END IF
        
        WRITE(*,*) 'SURFACE FINILAZATION FINISHED'
    
    END IF
        
    IF(TYPE1==-1) THEN

	CALL READINPUT

	SELECT CASE(RESTART_FLAG)
	CASE(0)
	  WRITE(*,*) 'SURFACE >> Start from Restart File'
	CASE(1)
	  WRITE(*,*) 'SURFACE >> Start from Initial Condition'
	CASE DEFAULT
	  WRITE(*,*) 'SURFACE >> Invalid Restart Flag'
	  STOP
	END SELECT
        
	IF(RESTART_FLAG .EQ. 1) THEN

		WRITE(*,*) 'SURFACE INITIALIZATION STARTED'
		CALL INIT_RANDOM_SEED()
		
		SURFACE_TOTAL_TIME = 0.
		
		SURFACE_FLUID%SURFACE_POINTS_NUM = 0
		SURFACE_FLUID%SURFACE_EDGES_NUM = 0
		SURFACE_PROPEL%SURFACE_POINTS_NUM = 0
		SURFACE_PROPEL%SURFACE_EDGES_NUM = 0
		SURFACE_CASE%SURFACE_POINTS_NUM = 0
		SURFACE_CASE%SURFACE_EDGES_NUM = 0
		
		INTERFACE_FLUID_POINTS_NUM = 0
		INTERFACE_STRUCT_POINTS_NUM = 0
		
		MINX = F_POINT(1,1)
		MAXX = F_POINT(1,1)
		MINY = F_POINT(2,1)
		MAXY = F_POINT(2,1)
		!$OMP PARALLEL DO PRIVATE(I),REDUCTION(MAX:MAXX,MAXY), REDUCTION(MIN:MINX,MINY)
		DO I = 2, F_POINT_NUM
		    IF(F_POINT(1,I) < MINX) THEN
		        MINX = F_POINT(1,I)
		    ENDIF
		    IF(F_POINT(1,I) > MAXX) THEN
		        MAXX = F_POINT(1,I)
		    ENDIF
		    IF(F_POINT(2,I) < MINY) THEN
		        MINY = F_POINT(2,I)
		    ENDIF
		    IF(F_POINT(2,I) > MAXY) THEN
		        MAXY = F_POINT(2,I)
		    ENDIF
		END DO
		!$OMP END PARALLEL DO
		
		!$OMP PARALLEL DO PRIVATE(I), REDUCTION(MAX:MAXX,MAXY), REDUCTION(MIN:MINX,MINY)
		DO I = 2, P_POINT_NUM
		    IF(P_POINT(1,I) < MINX) THEN
		        MINX = P_POINT(1,I)
		    ENDIF
		    IF(P_POINT(1,I) > MAXX) THEN
		        MAXX = P_POINT(1,I)
		    ENDIF
		    IF(P_POINT(2,I) < MINY) THEN
		        MINY = P_POINT(2,I)
		    ENDIF
		    IF(P_POINT(2,I) > MAXY) THEN
		        MAXY = P_POINT(2,I)
		    ENDIF
		END DO
		!$OMP END PARALLEL DO
		
		!$OMP PARALLEL DO PRIVATE(I),REDUCTION(MAX:MAXX,MAXY), REDUCTION(MIN:MINX,MINY)
		DO I = 2, C_POINT_NUM
		    IF(C_POINT(1,I) < MINX) THEN
		        MINX = C_POINT(1,I)
		    ENDIF
		    IF(C_POINT(1,I) > MAXX) THEN
		        MAXX = C_POINT(1,I)
		    ENDIF
		    IF(C_POINT(2,I) < MINY) THEN
		        MINY = C_POINT(2,I)
		    ENDIF
		    IF(C_POINT(2,I) > MAXY) THEN
		        MAXY = C_POINT(2,I)
		    ENDIF
		END DO
		!$OMP END PARALLEL DO
	    
		DOMAIN_MAX(1) = MAXX + (MAXX-MINX)*0.1
		DOMAIN_MIN(1) = MINX - (MAXX-MINX)*0.1
		DOMAIN_MAX(2) = MAXY + (MAXY-MINY)*0.1
		DOMAIN_MIN(2) = MINY - (MAXY-MINY)*0.1
		IF(DOMAIN_MAX(1) - DOMAIN_MIN(1) > DOMAIN_MAX(2) - DOMAIN_MIN(2)) THEN
		    DOMAIN_MAX(2) = (DOMAIN_MAX(2) + DOMAIN_MIN(2))/2. + (DOMAIN_MAX(1) - DOMAIN_MIN(1))/2.
		    DOMAIN_MIN(2) = (DOMAIN_MAX(2) + DOMAIN_MIN(2))/2. - (DOMAIN_MAX(1) - DOMAIN_MIN(1))/2.
		ELSE
		    DOMAIN_MAX(1) = (DOMAIN_MAX(1) + DOMAIN_MIN(1))/2. + (DOMAIN_MAX(2) - DOMAIN_MIN(2))/2.
		    DOMAIN_MIN(1) = (DOMAIN_MAX(1) + DOMAIN_MIN(1))/2. - (DOMAIN_MAX(2) - DOMAIN_MIN(2))/2.
		END IF
		
		!CALL FIND_PATCH_NUM(0)
		!CALL FIND_PATCH_NUM(1)
		!CALL FIND_PATCH_NUM(2)
		CALL RESET_SURFACE(0, .TRUE., F_POINT_NUM, F_EDGE_NUM, SURFACE_POINTS = F_POINT, SURFACE_EDGES = F_EDGE, &
		EDGE_LOCATION = F_LOC, EDGE_ABLATION_FLAG = F_BCFLAG)! , SURFACE_FLUID%SURFACE_PATCHES_NUM)
		
		CALL RESET_SURFACE(1, .TRUE., P_POINT_NUM, P_EDGE_NUM, SURFACE_POINTS = P_POINT, SURFACE_EDGES = P_EDGE, &
		EDGE_LOCATION = P_LOC) !SURFACE_PROPEL%SURFACE_PATCHES_NUM)
		
		CALL RESET_SURFACE(2, .TRUE., C_POINT_NUM, C_EDGE_NUM, SURFACE_POINTS = C_POINT, SURFACE_EDGES = C_EDGE) !SURFACE_CASE%SURFACE_PATCHES_NUM)
		
		CALL FIND_INTERFACE(0)
		CALL FIND_INTERFACE(1)
		CALL FIND_INTERFACE(2)
		
		CALL FIND_POINT_TYPE(0)
		CALL FIND_POINT_TYPE(1)
		CALL FIND_POINT_TYPE(2)
		
		CALL FIND_RELATEDPT(1,0,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))
		CALL FIND_RELATEDPT(0,1,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))

		CALL FIND_RELATEDEDGE(0,2,.FALSE.)
		CALL FIND_RELATEDEDGE(1,0,.FALSE.)

		CALL FIND_RELATEDEDGE(1,2,.TRUE.)
		CALL FIND_RELATEDEDGE(2,0,.TRUE.)
		
		CALL FIND_INTERFACE_CLUSTER(0)
		CALL FIND_INTERFACE_CLUSTER(1)
		
		CALL FIND_IMPACT_ZONE(0,0)
		CALL FIND_IMPACT_ZONE(0,2)

		CALL FIND_IMPACT_ZONE(1,1)
	       
		CALL FIND_PATCH_NUM(0)
		CALL FIND_PATCH_NUM(1)
		CALL FIND_PATCH_NUM(2) 
	       
		WRITE(*,*) 'SURFACE INITIALIZATION FINISHED'

	ELSEIF(RESTART_FLAG .EQ. 0) THEN

		CALL READ_RESTART(0)
		CALL READ_RESTART(1)
		CALL READ_RESTART(2)
	END IF

    END IF
        
        
    
    IF(TYPE1==0) THEN

        IF(TYPE2==0) THEN
                
            WRITE(*,*) 'FLUID SURFACE MOVING STARTED'
     
	    !$OMP PARALLEL DO PRIVATE(I)
            DO I = 1, SURFACE_FLUID%SURFACE_EDGES_NUM
                IF(SURFACE_FLUID%EDGE_ONINTERFACE(I)==1 .OR. SURFACE_FLUID%EDGE_ABLATION_FLAG(I) == -2) THEN
                    IF(PRESENT(F_B_RATE)) THEN
                    SURFACE_FLUID%EDGE_B_RATE(I) = F_B_RATE(I)
                    ELSE
                    SURFACE_FLUID%EDGE_B_RATE(I) = 0.008
                    END IF
                ELSE 
                    SURFACE_FLUID%EDGE_B_RATE(I) = 0.
                END IF
            END DO
            !$OMP END PARALLEL DO

            CALL COMPUTE_SURFACE_AREA(0)
                
            CALL FLUID_MOVE(TIMESTEP)
            SURFACE_TOTAL_TIME = SURFACE_TOTAL_TIME + TIMESTEP
            
	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1, SURFACE_FLUID%SURFACE_POINTS_NUM
                F_POINT(1,I) = SURFACE_FLUID%SURFACE_POINTS(1,I) + SURFACE_FLUID%POINT_DISPLACEMENT(1,I)
                F_POINT(2,I) = SURFACE_FLUID%SURFACE_POINTS(2,I) + SURFACE_FLUID%POINT_DISPLACEMENT(2,I)
            END DO
            !$OMP END PARALLEL DO

    	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1, SURFACE_PROPEL%SURFACE_POINTS_NUM
                P_POINT(1,I) = SURFACE_PROPEL%SURFACE_POINTS(1,I)
                P_POINT(2,I) = SURFACE_PROPEL%SURFACE_POINTS(2,I)
            END DO
	    !$OMP END PARALLEL DO
 
            WRITE(*,*) 'FLUID SURFACE MOVING ENDED'
        
        ELSE IF(TYPE2==1) THEN
            
            WRITE(*,*) 'FLUID SURFACE REMESHING STARTED'

            CALL ATTACH_FLUID_CASE_IMPACT_ZONE()
            CALL ZIPPER_FLUID_IMPACT_ZONE(FLAG)
            
            
            ALLOCATE(TYP1_POINT(2,INTERFACE_FLUID_POINTS_NUM))
            ALLOCATE(TYP1_POINTLOC(2,INTERFACE_FLUID_POINTS_NUM))
            ALLOCATE(DATA_DISPLACEMENT(2,SURFACE_FLUID%SURFACE_POINTS_NUM))
            
            TYP1_POINT_NUM = INTERFACE_FLUID_POINTS_NUM
            TYP1_POINT = INTERFACE_FLUID_POINTS
            TYP1_POINTLOC = INTERFACE_FLUID_POINTS_LOC
            DATA_DISPLACEMENT = SURFACE_FLUID%POINT_DISPLACEMENT
            
            !FLAG = .FALSE.
            !CALL EDGE_SPLITTING(0, TEMPFLAG)
            !FLAG = FLAG .OR. TEMPFLAG
            !CALL EDGE_COLLAPSING(0, TEMPFLAG)
            !FLAG = FLAG .OR. TEMPFLAG
            CALL EDGE_SPLITTING_COLLAPSING(0, FLAG)
            IF(FLAG) THEN
                CALL FIND_RELATEDPT(1,0,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))  
		CALL FIND_RELATEDPT(0,1,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))
                CALL FIND_RELATEDEDGE(0,2,.FALSE.)		
		CALL FIND_RELATEDEDGE(1,0,.FALSE.)
		CALL FIND_RELATEDEDGE(2,0,.TRUE.)
                CALL FIND_IMPACT_ZONE(0,0)
                CALL FIND_IMPACT_ZONE(0,2)
                
                CALL FIND_INTERFACE_CLUSTER(0)
                
                F_POINT_NUM = SURFACE_FLUID%SURFACE_POINTS_NUM
                DEALLOCATE(F_POINT)
                ALLOCATE(F_POINT(2,F_POINT_NUM))
                F_EDGE_NUM = SURFACE_FLUID%SURFACE_EDGES_NUM
                DEALLOCATE(F_EDGE)
                ALLOCATE(F_EDGE(2,F_EDGE_NUM))
                
                CALL INTERPOLATE_FLUID_DISPLACEMENT(TYP1_POINT_NUM, TYP1_POINT, TYP1_POINTLOC, DATA_DISPLACEMENT)
                
                ! testtest
                ! interpolation CALL 
                ! testtest
                
		!$OMP PARALLEL DO PRIVATE(I)
                DO I=1, SURFACE_FLUID%SURFACE_POINTS_NUM
                    F_POINT(1,I) = SURFACE_FLUID%SURFACE_POINTS(1,I) + SURFACE_FLUID%POINT_DISPLACEMENT(1,I)
                    F_POINT(2,I) = SURFACE_FLUID%SURFACE_POINTS(2,I) + SURFACE_FLUID%POINT_DISPLACEMENT(2,I)
                END DO
		!$OMP END PARALLEL DO
                
		!$OMP PARALLEL DO PRIVATE(I)
                DO I=1, SURFACE_FLUID%SURFACE_EDGES_NUM
                    F_EDGE(1,I) = SURFACE_FLUID%SURFACE_EDGES(1,I)
                    F_EDGE(2,I) = SURFACE_FLUID%SURFACE_EDGES(2,I)
                END DO
        	!$OMP END PARALLEL DO
                
                DEALLOCATE(F_LOC)
                ALLOCATE(F_LOC(F_EDGE_NUM))
                
                !$OMP PARALLEL DO PRIVATE(I)
                DO I=1, SURFACE_FLUID%SURFACE_EDGES_NUM
                    F_LOC(I) = SURFACE_FLUID%EDGE_LOCATION(I)
                END DO
		!$OMP END PARALLEL DO
            END IF

            DEALLOCATE(TYP1_POINT)
            DEALLOCATE(TYP1_POINTLOC)
            DEALLOCATE(DATA_DISPLACEMENT)

            WRITE(*,*) 'FLUID SURFACE REMESHING ENDED'

        END IF
            
            
    ELSE IF(TYPE1==1) THEN

        IF(TYPE2==0) THEN

            WRITE(*,*) 'STRUCT SURFACE MOVING STARTED'


            !$OMP PARALLEL DO PRIVATE(I)
            DO I = 1, SURFACE_PROPEL%SURFACE_POINTS_NUM
                SURFACE_PROPEL%POINT_VELOCITY(:,I) = P_DISPLACEMENT(:,I) - SURFACE_PROPEL%POINT_DISPLACEMENT(:,I)
                SURFACE_PROPEL%POINT_DISPLACEMENT(:,I) = P_DISPLACEMENT(:,I)
            END DO
            !$OMP END PARALLEL DO

	    !$OMP PARALLEL DO PRIVATE(I)
            DO I = 1, SURFACE_CASE%SURFACE_POINTS_NUM
                SURFACE_CASE%POINT_VELOCITY(:,I) = C_DISPLACEMENT(:,I) - SURFACE_CASE%POINT_DISPLACEMENT(:,I)
                SURFACE_CASE%POINT_DISPLACEMENT(:,I) = C_DISPLACEMENT(:,I)
            END DO
	    !$OMP END PARALLEL DO
            
            

            CALL STRUCT_MOVE()

	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1, SURFACE_PROPEL%SURFACE_POINTS_NUM
                P_POINT(1,I) = SURFACE_PROPEL%SURFACE_POINTS(1,I)
                P_POINT(2,I) = SURFACE_PROPEL%SURFACE_POINTS(2,I)
            END DO
            !$OMP END PARALLEL DO
            
	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1, SURFACE_CASE%SURFACE_POINTS_NUM
                C_POINT(1,I) = SURFACE_CASE%SURFACE_POINTS(1,I)
                C_POINT(2,I) = SURFACE_CASE%SURFACE_POINTS(2,I)
            END DO
            !$OMP END PARALLEL DO
            
	    !$OMP PARALLEL DO PRIVATE(I)
            DO I=1, SURFACE_FLUID%SURFACE_POINTS_NUM
                F_POINT(1,I) = SURFACE_FLUID%SURFACE_POINTS(1,I) + SURFACE_FLUID%POINT_DISPLACEMENT(1,I)
                F_POINT(2,I) = SURFACE_FLUID%SURFACE_POINTS(2,I) + SURFACE_FLUID%POINT_DISPLACEMENT(2,I)
            END DO
            !$OMP END PARALLEL DO
           
            WRITE(*,*) 'STRUCT SURFACE MOVING ENDED'

        ELSE IF(TYPE2==1) THEN

            WRITE(*,*) 'STRUCT SURFACE REMESHING STARTED'

	    SURFACE_PRESSURE_ITER = SURFACE_PRESSURE_ITER + 1

	    !ALLOCATE(NUM_PATCH(SURFACE_PROPEL%SURFACE_PATCHES_NUM)) 
	    !NUM_PATCH = 0
            
    	    !CALL ATTACH_PROPEL_FLUID_IMPACT_ZONE()]
            CALL ZIPPER_PROPEL_IMPACT_ZONE(FLAG)
            TEMPFLAG = .FALSE.
            IF(.NOT. FLAG) THEN
	        DO I=1,SURFACE_PROPEL%SURFACE_PATCHES_NUM
		    IF(SURFACE_PROPEL%SURFACE_PATCHES_TOPCHANGE_TYP(I) .EQ. 11) THEN
		       TEMPFLAG = .TRUE.
		    END IF
	        END DO
            END IF
            
            IF(.NOT. TEMPFLAG) THEN
                CALL EDGE_SPLITTING_COLLAPSING(1, TEMPFLAG, FLAG_ARRAY = FLAG_ARRAY)
		IF(.NOT. FLAG) THEN
	              SURFACE_FLAG_ARRAY(SURFACE_PRESSURE_ITER,1:SURFACE_PROPEL%SURFACE_PATCHES_NUM) = FLAG_ARRAY
                END IF
                FLAG = FLAG .OR. TEMPFLAG
            END IF
		
            
            IF(FLAG) THEN

!		call REARRANGE_INDEX(1)

		DO I=1,SURFACE_PROPEL%SURFACE_PATCHES_NUM
		   IF(SURFACE_PROPEL%SURFACE_PATCHES_TOPCHANGE_TYP(I) .GE. 2) THEN
	               FLAG_ARRAY(I) = SURFACE_PROPEL%SURFACE_PATCHES_TOPCHANGE_TYP(I)
		   END IF
		       SURFACE_FLAG_ARRAY(SURFACE_PRESSURE_ITER,I) = FLAG_ARRAY(I)
		END DO

		!DO I=1,SURFACE_PROPEL%SURFACE_PATCHES_NUM
		!    DO J=1,SURFACE_PROPEL%SURFACE_EDGES_NUM
		!	IF(SURFACE_PROPEL%EDGE_LOCATION(J) .EQ. I) THEN
		!	NUM_PATCH(I) = NUM_PATCH(I)+1
		!	END IF 
		!    END DO
		!
		!    IF(MOD(NUM_PATCH(I),2) .EQ. 1 .AND. SURFACE_PROPEL%SURFACE_PATCHES_TOPCHANGE_TYP(I) .NE. 11) THEN
	       	!	call NEW_STR_BDRY(I)
		!    END IF
		!END DO
                
!CALL CLASSIFY_PATCH(1)


                CALL FIND_RELATEDPT(1,0,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))
		CALL FIND_RELATEDPT(0,1,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))
		CALL FIND_RELATEDEDGE(1,0,.FALSE.)
                CALL FIND_RELATEDEDGE(1,2,.TRUE.)
  
        !CALL SAVINGDATA_TECPLOT(10000)
                
                CALL FIND_IMPACT_ZONE(1,1)
                
                CALL FIND_CASE_INTERFACE(TEMPFLAG)
                CALL FIND_INTERFACE_CLUSTER(1)
                
                P_POINT_NUM = SURFACE_PROPEL%SURFACE_POINTS_NUM
                DEALLOCATE(P_POINT)
                ALLOCATE(P_POINT(2,P_POINT_NUM))
                P_EDGE_NUM = SURFACE_PROPEL%SURFACE_EDGES_NUM
                DEALLOCATE(P_EDGE)
                ALLOCATE(P_EDGE(2,P_EDGE_NUM))
                
		!$OMP PARALLEL DO PRIVATE(I)
                DO I=1, SURFACE_PROPEL%SURFACE_POINTS_NUM
                    P_POINT(1,I) = SURFACE_PROPEL%SURFACE_POINTS(1,I)
                    P_POINT(2,I) = SURFACE_PROPEL%SURFACE_POINTS(2,I)
                END DO
		!$OMP END PARALLEL DO
                
		!$OMP PARALLEL DO PRIVATE(I)
                DO I=1, SURFACE_PROPEL%SURFACE_EDGES_NUM
                    P_EDGE(1,I) = SURFACE_PROPEL%SURFACE_EDGES(1,I)
                    P_EDGE(2,I) = SURFACE_PROPEL%SURFACE_EDGES(2,I)
                END DO
		!$OMP END PARALLEL DO
                
                DEALLOCATE(P_LOC)
                ALLOCATE(P_LOC(P_EDGE_NUM))
                
		!$OMP PARALLEL DO PRIVATE(I)
                DO I=1, SURFACE_PROPEL%SURFACE_EDGES_NUM
                    P_LOC(I) = SURFACE_PROPEL%EDGE_LOCATION(I)
                END DO
		!$OMP END PARALLEL DO
                
                DEALLOCATE(P_EDGE_LENGTH)
                ALLOCATE(P_EDGE_LENGTH(P_EDGE_NUM))
                DEALLOCATE(P_CORNER_INDEX)
                ALLOCATE(P_CORNER_INDEX(P_POINT_NUM))
                
                CALL SAVING_COEFFS_SURFACE_STRUCT(P_EDGE_LENGTH, P_CORNER_INDEX)
            END IF

  write(*,*) 'returned FLAG_ARRAY IS :', flag_array
            WRITE(*,*) 'STRUCT SURFACE REMESHING ENDED'
              
        END IF
        
    ELSE IF(TYPE1 ==2) THEN

        WRITE(*,*) 'STRUCT SURFACE RESET STARTED'
        
        ALLOCATE(TEMP_ONINTERFACE(P_EDGE_NUM))
        ALLOCATE(TEMP_CONNECTION(2,P_POINT_NUM))
        CALL LOADING_COEFFS_SURFACE_STRUCT(P_POINT_NUM, P_POINT, P_EDGE_NUM, P_EDGE, P_LOC, TEMP_ONINTERFACE, TEMP_CONNECTION)
        CALL RESET_SURFACE(1, .TRUE., P_POINT_NUM, P_EDGE_NUM, SURFACE_POINTS = P_POINT, SURFACE_EDGES = P_EDGE, POINT_EDGE_CONNECTION = TEMP_CONNECTION, EDGE_LOCATION = P_LOC, EDGE_ONINTERFACE = TEMP_ONINTERFACE) 
!, & SURFACE_PROPEL%SURFACE_PATCHES_NUM)
        DEALLOCATE(TEMP_ONINTERFACE)
        DEALLOCATE(TEMP_CONNECTION)      
        CALL FIND_INTERFACE_CLUSTER(1)              
        CALL FIND_IMPACT_ZONE(1,1)       
        CALL FIND_RELATEDPT(1,0,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))
	CALL FIND_RELATEDPT(0,1,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))
        CALL FIND_RELATEDEDGE(1,2,.TRUE.)
	CALL FIND_RELATEDEDGE(1,0,.FALSE.)

        WRITE(*,*) 'STRUCT SURFACE RESET ENDED'
  
    ELSE IF(TYPE1==3) THEN

        WRITE(*,*) 'PRESSURE TRANSFER STARTED'
        

        !$OMP PARALLEL DO PRIVATE(I)
        DO I = 1, SURFACE_PROPEL%SURFACE_POINTS_NUM
           SURFACE_PROPEL%POINT_DISPLACEMENT(:,I) = P_DISPLACEMENT(:,I)
        END DO
        !$OMP END PARALLEL DO

	!$OMP PARALLEL DO PRIVATE(I)
        DO I = 1, SURFACE_CASE%SURFACE_POINTS_NUM
            SURFACE_CASE%POINT_DISPLACEMENT(:,I) = C_DISPLACEMENT(:,I)
        END DO
	!$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(I)
        DO I = 1, SURFACE_FLUID%SURFACE_EDGES_NUM
            SURFACE_FLUID%EDGE_PRESSURE(I) = F_PRESSURE(I)
        END DO
        !$END OMP PARALLEL DO
        
        CALL PRESSURE_TRANSFER()
        
	!$OMP PARALLEL DO PRIVATE(I)
        DO I = 1, SURFACE_PROPEL%SURFACE_POINTS_NUM
            P_FORCE(1,I) = SURFACE_PROPEL%POINT_FORCE(1,I)
            P_FORCE(2,I) = SURFACE_PROPEL%POINT_FORCE(2,I)
        END DO
        !$END OMP PARALLEL DO
        
	!$OMP PARALLEL DO PRIVATE(I)
        DO I = 1, SURFACE_CASE%SURFACE_POINTS_NUM
            C_FORCE(1,I) = SURFACE_CASE%POINT_FORCE(1,I)
            C_FORCE(2,I) = SURFACE_CASE%POINT_FORCE(2,I)
        END DO
        !$END OMP PARALLEL DO
        
        WRITE(*,*) 'PRESSURE TRANSFER ENDED'
  
    ELSE IF(TYPE1==4) THEN

        WRITE(*,*) 'PATCH NUMBER FINDING STARTED'
        
        CALL FIND_PATCH_NUM(1)
        PATCHNUM = SURFACE_PROPEL%SURFACE_PATCHES_NUM
        
        WRITE(*,*) 'PATCH NUMBER FINDING ENDED'

    ELSE IF(TYPE1==5) THEN

        WRITE(*,*) 'EDGE_REARRANGEMENT STARTED'

	!CALL REARRANGE_EDGE(1, P_POINT_NUM, P_POINT, P_EDGE_NUM, P_EDGE)
	CALL REARRANGE_PROPEL_INDEX(P_POINT_NUM, P_POINT, P_EDGE_NUM, P_EDGE, P_LOC)	
	CALL FIND_INTERFACE_CLUSTER(1)
        CALL FIND_IMPACT_ZONE(1,1)       
        CALL FIND_RELATEDPT(1,0,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))
	CALL FIND_RELATEDPT(0,1,MAX(DOMAIN_MAX(1) - DOMAIN_MIN(1), DOMAIN_MAX(2) - DOMAIN_MIN(2)))
        CALL FIND_RELATEDEDGE(1,2,.TRUE.)
	CALL FIND_RELATEDEDGE(1,0,.FALSE.)
      
!	call REARRANGE_INDEX(1) 
        
        WRITE(*,*) 'EDGE_REARRANGEMENT ENDED'

    ELSE IF(TYPE1==100) THEN

!        CALL SAVINGPOINT(FILENUM,0)
!        CALL SAVINGEDGE(FILENUM,0)
!        CALL SAVINGPOINT(FILENUM,1)
!        CALL SAVINGEDGE(FILENUM,1)
!        CALL SAVINGPOINT(FILENUM,2)
!        CALL SAVINGEDGE(FILENUM,2)
        
        CALL SAVINGDATA_TECPLOT(FILENUM)

    ELSE IF(TYPE1==101) THEN

        !CALL SAVINGDISPLACEMENT(FILENUM,0)
!        CALL SAVINGDISPLACEMENT(FILENUM,1)
!        CALL SAVINGDISPLACEMENT(FILENUM,2)

    ELSE IF(TYPE1==102) THEN

        !CALL SAVINGBRATE(FILENUM,0)
        !CALL SAVINGBRATE(FILENUM,1)
        !CALL SAVINGBRATE(FILENUM,2)

    ELSE IF(TYPE1==103) THEN

        !CALL SAVINGPRESSURE(FILENUM,0)
        !CALL SAVINGPRESSURE(FILENUM,1)
        !CALL SAVINGPRESSURE(FILENUM,2)
    
    ELSE IF(TYPE1==104) THEN
        CALL SAVINGINTERFACE_TECPLOT(FILENUM)
 
    ELSE IF(TYPE1==105) THEN
	!CALL SAVING_SURFACE(FILENUM,0)
	!CALL SAVING_SURFACE(FILENUM,1)
	!CALL SAVING_SURFACE(FILENUM,2)
    END IF
    
END SUBROUTINE SURFACE_2D

