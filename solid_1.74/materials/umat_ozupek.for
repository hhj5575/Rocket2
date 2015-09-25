      SUBROUTINE UMAT_01(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
      CHARACTER*8 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),HDDSDDE(NTENS,NTENS),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),DFGRD0(3,3),DFGRD1(3,3),
     3 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DSTRESS(6),DDSTRESS(6,6)

C
C    LOCAL ARRAYS
C ----------------------------------------------------------------
C    BBAR   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C    DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
C ----------------------------------------------------------------
C
C
C     NUMBER OF PRONY SERIES

      DIMENSION BBAR(NTENS),DISTGR(3,3)
      DIMENSION DELTA(3,3),BBBAR(3,3),RELE((NPROPS-44)/2)
      DIMENSION TAU((NPROPS-44)/2)
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)
      DIMENSION X1(3,3),X2(6),X3(6,6)
      DIMENSION HDDSDDE1(6,6),HDDSDDE2(6,6),HDDSDDE3(6,6)
      DIMENSION HN(6,((NPROPS-44)/2)),CBAR(3,3),PRESTR(6)   
      DIMENSION HISTR(((NPROPS-44)/2)*6),HISTR1(((NPROPS-44)/2)*6)  

C      write(*,*) 'UMAT*: DFGRD1', DFGRD1
C      write(*,*) 'UMAT*: STATEV', STATEV
C      write(*,*) 'UMAT*: TEMP', TEMP

      CYCEG = 0.0

      PRONY=PROPS(1)
      NPRONY=PRONY

C     INFINITE PROPERTY            
      EEQ=PROPS(2)

C     COEFFICIENTS OF PRONY SERIES
      DO I=1,NPRONY
         RELE(I)=PROPS(2+I)
      ENDDO
      
C     RELAXATION TIME
      DO I=1,NPRONY
         TAU(I)=PROPS(2+NPRONY+I)
      ENDDO
      
C     DEWETTING FUNCTION 1=> ACTIVE 0=>DEACTIVE      
      DEWETTING=PROPS(2+2*NPRONY+1)

C     CYC FUNCTION 1=> ACTIVE 0=>DEACTIVE          
      CYCFUNCTION=PROPS(2+2*NPRONY+2)

C     SCALE FACTOR TO CONVERT THE UNIT (DEFAULT : BAR)         
      SCALE_FACTOR=PROPS(2+2*NPRONY+3)
      
C     SCALE FACTOR OF TIME (DEFAULT : MIN)       
      SCALE_TIME=PROPS(2+2*NPRONY+4)

C     REFERENCE TEMPERATURE FOR THERMAL EXPANSION    
      TEMPREF=PROPS(2+2*NPRONY+5)
      
C     THERMALE EXPANSION COEFFICIENT     
      AL_THER=PROPS(2+2*NPRONY+6)
      
C     DILATATION PARAMETERS      
      W1=PROPS(2+2*NPRONY+7)
      W2=PROPS(2+2*NPRONY+8)
      W3=PROPS(2+2*NPRONY+9)
      WN=PROPS(2+2*NPRONY+10)    
      BULK=PROPS(2+2*NPRONY+11)*SCALE_FACTOR


C     DEWETTING PARAMETERS      
C     SIGC=(DP1*TEMP+DP2)/(TEMP+DP3)
      DP1=PROPS(2+2*NPRONY+12)
      DP2=PROPS(2+2*NPRONY+13)
      DP3=PROPS(2+2*NPRONY+14)
      
C     DAMAGE FUNCTION CONSTANT
      ALPHA_G=PROPS(2+2*NPRONY+15)

C     FREE ENERGY PARAMETERS         
      BETA1=PROPS(2+2*NPRONY+16)*SCALE_FACTOR
      BETA4=PROPS(2+2*NPRONY+17)*SCALE_FACTOR
      BETA5=PROPS(2+2*NPRONY+18)  
      
C     WLF CONSTANTS
      WLFC1=PROPS(2+2*NPRONY+19)
      WLFC2=PROPS(2+2*NPRONY+20)


C     CYCLIC LOAD PARAMETERS
      D_UN_UP_S= PROPS(2+2*NPRONY+21)
      D_UN_LB  = PROPS(2+2*NPRONY+22)
      D_UN_LB_S= PROPS(2+2*NPRONY+23)
      D_UN_P1  = PROPS(2+2*NPRONY+24)
      D_UN_P2  = PROPS(2+2*NPRONY+25)
      D_UN_P3  = PROPS(2+2*NPRONY+26)
      D_UN_P4  = PROPS(2+2*NPRONY+27)
      D_UN_P5  = PROPS(2+2*NPRONY+28)
      D_UN_P6  = PROPS(2+2*NPRONY+29)
      D_RE_UB  = PROPS(2+2*NPRONY+30)
      D_RE_UP_S= PROPS(2+2*NPRONY+31)
      D_RE_LB  = PROPS(2+2*NPRONY+32)
      D_RE_LB_S= PROPS(2+2*NPRONY+33)
      D_RE_P1  = PROPS(2+2*NPRONY+34)
      D_RE_P2  = PROPS(2+2*NPRONY+35)
      D_RE_P3  = PROPS(2+2*NPRONY+36)
      D_RE_P4  = PROPS(2+2*NPRONY+37)
      D_RE_P5  = PROPS(2+2*NPRONY+38)
      D_RE_P6  = PROPS(2+2*NPRONY+39)
      D_RE_P7  = PROPS(2+2*NPRONY+40)
      D_RE_P8  = PROPS(2+2*NPRONY+41)
      D_UPPER  = PROPS(2+2*NPRONY+42)                    
                 

C      AL_THER=0.D0

C ----------------------------------------------------------------
C     STATE VARIABLES
C ----------------------------------------------------------------


C     STATEV(1)-(6) : PREVIOUS STRESS      
C     STATEV(7)-(6+NPRONY*6) : PREVIOUS HN1 


C     6+NPRONY*6 + ADD1

C     STATEV(6+NPRONY*6+1) : OCTAHEDRAL STRESS 

C     STATEV(6+NPRONY*6+2) : PSEUDO HYDROPRESSURE
C     STATEV(6+NPRONY*6+3) : PREVIOUS PSEUDO HYDROPRESSURE   
C     STATEV(6+NPRONY*6+4) : VOID
C     STATEV(6+NPRONY*6+5) : PREVIOUS VOID
C     STATEV(6+NPRONY*6+6) : DAMAGED BULKMODULUS
C     STATEV(6+NPRONY*6+7) : PREVIOUS BULK

C     STATEV(6+NPRONY*6+8) : OCTAHEDRAL STRAIN AT PREVIOUS STEP
C     STATEV(6+NPRONY*6+9) : PREVIOUS DAMAGEG FUNCTION

C     STATEV(6+NPRONY*6+10) : NUMBER OF CYLCIC LOADING
C     STATEV(6+NPRONY*6+11) : 
C     STATEV(6+NPRONY*6+12) : XIR1, MAXIMUM DURING LOADING
C     STATEV(6+NPRONY*6+13) : MINIMUM DURING UNLOADING
C     STATEV(6+NPRONY*6+14) : PREVIOUS STATEV(6+NPRONY*6+13)
C     STATEV(6+NPRONY*6+15) : IF=1, RELOADING 
C     STATEV(6+NPRONY*6+16) : PREVIOUS DAMAGEF   


C
C    JACOBIAN AND DISTORTION TENSOR
C
      DET=DFGRD1(1, 1)*DFGRD1(2, 2)*DFGRD1(3, 3)
     1   -DFGRD1(1, 2)*DFGRD1(2, 1)*DFGRD1(3, 3)
      IF(NSHR.EQ.3) THEN
        DET=DET+DFGRD1(1, 2)*DFGRD1(2, 3)*DFGRD1(3, 1)
     1         +DFGRD1(1, 3)*DFGRD1(3, 2)*DFGRD1(2, 1)
     2         -DFGRD1(1, 3)*DFGRD1(3, 1)*DFGRD1(2, 2)
     3         -DFGRD1(2, 3)*DFGRD1(3, 2)*DFGRD1(1, 1)
      END IF
      SCALE=DET**(-ONE/THREE)
      DO K1=1, 3
        DO K2=1, 3
          DISTGR(K2, K1)=SCALE*DFGRD1(K2, K1)
        END DO
      END DO
C
C    CALCULATE LEFT CAUCHY-GREEN TENSOR
C
      BBAR(1)=DISTGR(1, 1)**2+DISTGR(1, 2)**2+DISTGR(1, 3)**2
      BBAR(2)=DISTGR(2, 1)**2+DISTGR(2, 2)**2+DISTGR(2, 3)**2
      BBAR(3)=DISTGR(3, 3)**2+DISTGR(3, 1)**2+DISTGR(3, 2)**2
      BBAR(4)=DISTGR(1, 1)*DISTGR(2, 1)+DISTGR(1, 2)*DISTGR(2, 2)
     1       +DISTGR(1, 3)*DISTGR(2, 3)
      IF(NSHR.EQ.3) THEN
        BBAR(5)=DISTGR(1, 1)*DISTGR(3, 1)+DISTGR(1, 2)*DISTGR(3, 2)
     1         +DISTGR(1, 3)*DISTGR(3, 3)
        BBAR(6)=DISTGR(2, 1)*DISTGR(3, 1)+DISTGR(2, 2)*DISTGR(3, 2)
     1         +DISTGR(2, 3)*DISTGR(3, 3)
      END IF      
C

      IF (KINC .EQ. 1) THEN
      STATEV(6+NPRONY*6+7)=BULK
      STATEV(6+NPRONY*6+9)=1.D0
      STATEV(6+NPRONY*6+10)=0.D0
C      STATEV(6+NPRONY*6+3)=0.D0
      ENDIF


      CALL CALINVARIANT(BBAR,BAR1,BAR2,BOCTA,NSHR)

C     VOID EVOLUTION
      STATEV(6+NPRONY*6+4)=STATEV(6+NPRONY*6+4)
     1 +W1*(BOCTA**WN-STATEV(6+NPRONY*6+8)**WN)
     2         *2.7183**(STATEV(6+NPRONY*6+2)/W2)
          
      IF (STATEV(6+NPRONY*6+4) .LE. 0.D0) THEN
      STATEV(6+NPRONY*6+4)=0.D0
      ENDIF

C     STATEV(6+NPRONY*6+6) : DAMAGED BULK MODULUS 
      IF (DEWETTING .EQ. 1.D0) THEN
      STATEV(6+NPRONY*6+6)=BULK/(1+W3*BULK*STATEV(6+NPRONY*6+4))
     1 *(1-STATEV(6+NPRONY*6+4)) 
      
      IF (STATEV(6+NPRONY*6+6) .GE. STATEV(6+NPRONY*6+7) ) THEN
      STATEV(6+NPRONY*6+6)=STATEV(6+NPRONY*6+7)
      ENDIF   
            
      ELSE
      STATEV(6+NPRONY*6+6)=BULK 
      ENDIF      
      
      IF (STATEV(6+NPRONY*6+6) .LE. 5.D0) THEN
      STATEV(6+NPRONY*6+6)=5.D0
      ENDIF
C
C     CYCLIC LOAD DAMAGE FUNCTION
C
      IF (CYCFUNCTION .EQ. 1.D0) THEN
      CALL LOADTYPE(BOCTA,STATEV(6+NPRONY*6+8),STATEV(6+NPRONY*6+9),
     1 STATEV(6+NPRONY*6+10),STATEV(6+NPRONY*6+12),
     2 KINC,AA1,XIR1,LOADCASE,STATEV(6+NPRONY*6+14),
     3 STATEV(6+NPRONY*6+15),STATEV(6+NPRONY*6+16),
     4 STATEV(6+NPRONY*6+13),XIRU)

 
C       write(*,*) 'UMAT******* 5: loadcase', LOADCASE

       CALL CALDAMAGEF(AA1,LOADCASE,STATEV(6+NPRONY*6+10),
     1 STATEV(6+NPRONY*6+15),
     2 XIR1,XIRU,BOCTA,DAMAGEF,DERDF,
     3 D_UN_UP_S, D_UN_LB, D_UN_LB_S, D_UN_P1, D_UN_P2,D_UN_P3,
     4 D_UN_P4, D_UN_P5, D_UN_P6, D_RE_UB, D_RE_UP_S,D_RE_LB,D_RE_LB_S,
     5 D_RE_P1,D_RE_P2,D_RE_P3,D_RE_P4,D_RE_P5,D_RE_P6,D_RE_P7,
     6 D_RE_P8,D_UPPER)     
     

      ELSE


      DAMAGEF=1.D0
      DERDF=0.D0
      ENDIF


C     CRITICAL STRESS (DAMAGE INITIATION)
      SIGC=(DP1*TEMP+DP2)/(TEMP+DP3)*SCALE_FACTOR

      IF (DEWETTING .EQ. 1.D0) THEN
      CALL CALDAMAGEG(STATEV(6+NPRONY*6+1),STATEV(6+NPRONY*6+4),
     1 STATEV(6+NPRONY*6+5),TEMP,DAMAGEG,
     2 STATEV(6+NPRONY*6+9),SIGC,ALPHA_G)
      ELSE
      DAMAGEG=1.D0
      STATEV(6+NPRONY*6+4)=0.D0
      STATEV(6+NPRONY*6+6)=BULK
      END IF

C
C     CALCULATE STRESS
C
C     PREVIOUS SETP STRESS
      DO I=1,NTENS
         PRESTR(I)=STATEV(I)
      ENDDO
      
C     HISTORY TERMS   
      DO I=1,NPRONY*NTENS
         HISTR(I)=STATEV(I+6)
      ENDDO


      AT=10**((-WLFC1*(TEMP-20))/(WLFC2+(TEMP-20)))
      DETA=DTIME/AT/SCALE_TIME


C     CURRENT STEP STRESS 
      CALL CALSTRESS(BETA1,BETA4,BETA5,FIRSTD,SECONDD,TRBBAR,BBAR,
     1 EG,EK,PR,DET,DAMAGEF,DAMAGEG,STATEV(6+NPRONY*6+6),
     2 STATEV(6+NPRONY*6+4),NDI,NSHR,TEMP,
     3 BAR1,DETA,AT,HISTR,PRESTR,EEQ,RELE,TAU,DSTRESS,HISTR1,
     4 STRESS,AL_THER,TEMPREF,NPRONY)


C     STATEV(6+NPRONY*6+2) : PEUDO HYDRO PRESSURE       
      STATEV(6+NPRONY*6+2)=PR   


C     0.3 IS UPPER AND LOWER BOUND FOR BETTER CONVERGENCE
C     RECOMMENDED VALUE 0~0.3
      IF (STATEV(6+NPRONY*6+2) .LE. -0.2D0) THEN
      STATEV(6+NPRONY*6+2)=STATEV(6+NPRONY*6+3)
      ENDIF
      IF (STATEV(6+NPRONY*6+2) .GE. 0.2D0 ) THEN
      STATEV(6+NPRONY*6+2)=STATEV(6+NPRONY*6+3)
      ENDIF

      
         
      DO K=1,NTENS
         STATEV(K)=DSTRESS(K)           
      ENDDO
      DO I=1,NPRONY*NTENS
         STATEV(I+6)=HISTR1(I)
      ENDDO  
  

      CALL CALVONMISES(STRESS,STATEV(6+NPRONY*6+1),NTENS)


      CALL HYPERJ(HDDSDDE1,HDDSDDE2,BBAR,EG,EK,SECONDD,DET,NTENS,
     1 BAR1,NSHR,TRBBAR,SEEG)

C     CYCLING LOADING TERMS TANAGENT STIFFNESS

      IF (STATEV(6+NPRONY*6+10) .GE. 1.D0) THEN      
      CALL CYCSTIFF(FIRSTD,DET,DERDF,XIR1,BBAR,BAR1,BAR2,HDDSDDE3,
     1 CYCEG)      
      ELSE     
      DO I = 1, 6
        DO J = 1, 6
          HDDSDDE3(I,J) =0.D0
        ENDDO
      ENDDO     
      ENDIF


C     PSEUOD ELASTIC TANGENT STIFFNESS
      DO I=1,NTENS
       DO J=1,NTENS
         HDDSDDE(I,J)=HDDSDDE1(I,J)
     1    +SEEG*DAMAGEG*DAMAGEF*HDDSDDE2(I,J)
     2    +CYCEG*DAMAGEG*HDDSDDE3(I,J)
       ENDDO
      ENDDO   



C     VISCOELASTIC MODULUS         
          ALPHA=EEQ
      DO I=1,NPRONY
          ALPHA=ALPHA+RELE(I)*TAU(I)*(1-2.7183**(-DETA/TAU(I)))/DETA
      END DO

C     VISCOELASTIC TANGENT STIFFNESS
      DO I=1,NTENS
      DO J=1,NTENS
        DDSDDE(I,J)=ALPHA*HDDSDDE(I,J)
      END DO
      END DO     
      
      
      
C     UPDATE STATEV VARIABLES      

      STATEV(6+NPRONY*6+8)=BOCTA   
      STATEV(6+NPRONY*6+9)=DAMAGEG                     
      STATEV(6+NPRONY*6+7)=STATEV(6+NPRONY*6+6)
      STATEV(6+NPRONY*6+3)=STATEV(6+NPRONY*6+2)
      STATEV(6+NPRONY*6+14)=STATEV(6+NPRONY*6+13)
      STATEV(6+NPRONY*6+16)=DAMAGEF


C      write(*,*) 'UMAT*: STRESS', STRESS
C      write(*,*) 'UMAT*: DDSDDE', DDSDDE
      RETURN
      END
      

      
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C  SUBROUTINES 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
          

 
C-----------------------------------------------------------------------
      SUBROUTINE CALVONMISES(STRESS,OUTPUT,NTENS)
C
C     CALCULATE VOM MISES STRESS
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION STRESS(NTENS), OUTPUT

      
      
      IF (NTENS .EQ. 4) THEN
      OUTPUT=(((STRESS(1)-STRESS(2))**2+(STRESS(2)-STRESS(3))**2
     1 +(STRESS(3)-STRESS(1))**2
     2 +6*(STRESS(4)**2))/2.D0)**(0.5D0)         

  
      ELSE
       OUTPUT=(((STRESS(1)-STRESS(2))**2+(STRESS(2)-STRESS(3))**2
     1 +(STRESS(3)-STRESS(1))**2
     2 +6*(STRESS(4)**2+STRESS(5)**2+STRESS(6)**2))/2.D0)**(0.5D0)       

      ENDIF

      RETURN
      END

      
C-----------------------------------------------------------------------
      SUBROUTINE CYCSTIFF(FIRSTD,DET,DERDF,XIR1,BBAR,BAR1,BAR2,HDDSDDE3,
     1 CYCEG)
C
C     CALCULATE JACOBIAN OF CYCLIC LOADING
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION X3(6,6),BBBAR(3,3),XTERM,X1(3,3),X2(6),
     1 HDDSDDE3(6,6),DET,DERDF,XIR1,BBAR(6),BAR1,BAR2,CYCEG,FIRSTD


      DO I=1,6
        DO J=1,6
        X3(I,J)=0.D0
        ENDDO
      ENDDO

      CYCEG=2.D0*FIRSTD/DET*DERDF*1.D0/XIR1

      BBBAR(1,1)=BBAR(1)
      BBBAR(1,2)=BBAR(4)
      BBBAR(1,3)=BBAR(5)
      
      BBBAR(2,1)=BBAR(4)
      BBBAR(2,2)=BBAR(2)
      BBBAR(2,3)=BBAR(6)
      
      BBBAR(3,1)=BBAR(5)
      BBBAR(3,2)=BBAR(6)
      BBBAR(3,3)=BBAR(3)       
      
      XTERM=(2.D0*BAR1**2-6.D0*BAR2)**(0.5D0)

      
      IF(XTERM .LE. 1.E-20) THEN
        XTERM = 1.E-10
      ENDIF


      DO I = 1, 3
        DO J = I, 3
          X1(I,J) = 0.D0
          DO K = 1, 3
            X1(I,J) = X1(I,J) + BBBAR(I,K)*BBBAR(K,J)
          ENDDO
        ENDDO
      ENDDO
      
      
      X2(1)=X1(1,1)
      X2(2)=X1(2,2)
      X2(3)=X1(3,3)
      X2(4)=X1(1,2)
      X2(5)=X1(1,3)
      X2(6)=X1(2,3)
      
      
      DO I=1,6
        DO J=1,6
            X3(I,J)=-1.D0/3.D0/XTERM*BAR1*BBAR(I)*BBAR(J)               
        ENDDO
      ENDDO
        
      DO I=1,6
        DO J=1,6
            X3(I,J)=X3(I,J)+1.D0/XTERM*BBAR(I)*X2(J)               
        ENDDO
      ENDDO          
        
        
      DO I=1,6
        DO J=1,3
            X3(I,J)=X3(I,J)-1.D0/9.D0*XTERM*BBAR(I)               
        ENDDO
      ENDDO      
      
      
      DO I=1,3
        DO J=1,6
            X3(I,J)=X3(I,J)+BAR1*BAR1/9.D0/XTERM*BBAR(J)               
        ENDDO
      ENDDO       
      
      DO I=1,3
        DO J=1,6
            X3(I,J)=X3(I,J)-BAR1/3.D0/XTERM*X2(J)               
        ENDDO
      ENDDO        
      
      DO I=1,3
        DO J=1,3
            X3(I,J)=X3(I,J)+BAR1/27.D0*XTERM               
        ENDDO
      ENDDO        
   
      
      DO I = 1, 6
        DO J = 1, 6
          HDDSDDE3(I,J) = (X3(I,J)+X3(J,I))/2.D0
        ENDDO
      ENDDO
     
C
      RETURN
      END      
      
 
      
C-----------------------------------------------------------------------
      SUBROUTINE HYPERJ(HDDSDDE1,HDDSDDE2,BBAR,EG,EK,SECONDD,DET,NTENS,
     1 BAR1,NSHR,TRBBAR,SEEG)
C
C     CALCULATE JACOBIAN OF FREE ENERGY
C
C-----------------------------------------------------------------------
      
      
      DOUBLE PRECISION HDDSDDE1(6,6), HDDSDDE2(6,6),BBAR(6),EG,EK,
     1 SECONDD,DET,BAR1,EG23,SEEG,TRBBAR,TWO,THREE
     
      INTEGER NTENS,NSHR  
      
      TWO=2.D0
      THREE=3.D0
      
      EG23=EG*TWO/THREE
      HDDSDDE1(1, 1)= EG23*(BBAR(1)+TRBBAR)+EK
      HDDSDDE1(2, 2)= EG23*(BBAR(2)+TRBBAR)+EK
      HDDSDDE1(3, 3)= EG23*(BBAR(3)+TRBBAR)+EK
      HDDSDDE1(1, 2)=-EG23*(BBAR(1)+BBAR(2)-TRBBAR)+EK
      HDDSDDE1(1, 3)=-EG23*(BBAR(1)+BBAR(3)-TRBBAR)+EK
      HDDSDDE1(2, 3)=-EG23*(BBAR(2)+BBAR(3)-TRBBAR)+EK
      HDDSDDE1(1, 4)= EG23*BBAR(4)/TWO
      HDDSDDE1(2, 4)= EG23*BBAR(4)/TWO
      HDDSDDE1(3, 4)=-EG23*BBAR(4)
      HDDSDDE1(4, 4)= EG*(BBAR(1)+BBAR(2))/TWO
      IF(NSHR.EQ.3) THEN
        HDDSDDE1(1, 5)= EG23*BBAR(5)/TWO
        HDDSDDE1(2, 5)=-EG23*BBAR(5)
        HDDSDDE1(3, 5)= EG23*BBAR(5)/TWO
        HDDSDDE1(1, 6)=-EG23*BBAR(6)
        HDDSDDE1(2, 6)= EG23*BBAR(6)/TWO
        HDDSDDE1(3, 6)= EG23*BBAR(6)/TWO
        HDDSDDE1(5, 5)= EG*(BBAR(1)+BBAR(3))/TWO
        HDDSDDE1(6, 6)= EG*(BBAR(2)+BBAR(3))/TWO
        HDDSDDE1(4, 5)= EG*BBAR(6)/TWO
        HDDSDDE1(4, 6)= EG*BBAR(5)/TWO
        HDDSDDE1(5, 6)= EG*BBAR(4)/TWO
      END IF
      
      DO K1=1, NTENS
        DO K2=1, K1-1
          HDDSDDE1(K1, K2)=HDDSDDE1(K2, K1)
        END DO
      END DO
      

C     SECOND DERIVATIVE TERMS
      
      SEEG=4*SECONDD/DET

      HDDSDDE2(1, 1)=BBAR(1)*BBAR(1)-2.D0/3.D0*BAR1*BBAR(1)
     1 +BAR1**2.D0/9.D0
      HDDSDDE2(2, 2)=BBAR(2)*BBAR(2)-2.D0/3.D0*BAR1*BBAR(2)
     1 +BAR1**2.D0/9.D0
      HDDSDDE2(3, 3)=BBAR(3)*BBAR(3)-2.D0/3.D0*BAR1*BBAR(3)
     1 +BAR1**2.D0/9.D0
      HDDSDDE2(1, 2)=BBAR(1)*BBAR(2)-1.D0/3.D0*BAR1*BBAR(1)
     1 -1.D0/3.D0*BAR1*BBAR(2)+BAR1**2.D0/9.D0
      HDDSDDE2(1, 3)=BBAR(1)*BBAR(3)-1.D0/3.D0*BAR1*BBAR(1)
     1 -1.D0/3.D0*BAR1*BBAR(3)+BAR1**2.D0/9.D0      
      HDDSDDE2(2, 3)=BBAR(2)*BBAR(3)-1.D0/3.D0*BAR1*BBAR(2)
     1 -1.D0/3.D0*BAR1*BBAR(3)+BAR1**2.D0/9.D0
      HDDSDDE2(1, 4)=BBAR(1)*BBAR(4)-1.D0/3.D0*BBAR(4)
      HDDSDDE2(2, 4)=BBAR(2)*BBAR(4)-1.D0/3.D0*BBAR(4)
      HDDSDDE2(3, 4)=BBAR(3)*BBAR(4)-1.D0/3.D0*BBAR(4)
      HDDSDDE2(4, 4)=BBAR(4)*BBAR(4)
      IF(NSHR.EQ.3) THEN
        HDDSDDE2(1, 5)=BBAR(1)*BBAR(5)-1.D0/3.D0*BBAR(5)
        HDDSDDE2(2, 5)=BBAR(2)*BBAR(5)-1.D0/3.D0*BBAR(5)
        HDDSDDE2(3, 5)=BBAR(3)*BBAR(5)-1.D0/3.D0*BBAR(5)
        HDDSDDE2(1, 6)=BBAR(1)*BBAR(6)-1.D0/3.D0*BBAR(6)
        HDDSDDE2(2, 6)=BBAR(2)*BBAR(6)-1.D0/3.D0*BBAR(6) 
        HDDSDDE2(3, 6)=BBAR(3)*BBAR(6)-1.D0/3.D0*BBAR(6) 
        HDDSDDE2(5, 5)=BBAR(5)*BBAR(5) 
        HDDSDDE2(6, 6)=BBAR(6)*BBAR(6) 
        HDDSDDE2(4, 5)=BBAR(4)*BBAR(5) 
        HDDSDDE2(4, 6)=BBAR(4)*BBAR(6) 
        HDDSDDE2(5, 6)=BBAR(5)*BBAR(6) 
      END IF
      
      DO K1=1, NTENS
        DO K2=1, K1-1
          HDDSDDE2(K1, K2)=HDDSDDE2(K2, K1)
        END DO
      END DO


C
      RETURN
      END      
      
      
      
 
C-----------------------------------------------------------------------
      SUBROUTINE CALSTRESS(BETA1,BETA4,BETA5,FIRSTD,SECONDD,TRBBAR,BBAR,
     1 EG,EK,PR,DET,DAMAGEF,DAMAGEG,DBULK,VOID,NDI,NSHR,TEMP,BAR1,DETA,
     2 AT, HISTR,PRESTR,EEQ,RELE,TAU,DSTRESS,HISTR1,STRESS,
     3 AL_THER,TEMPI,NPRONY)
C
C     CALCULATE VISCOELASTIC STRESS
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION BETA1,BETA4,BETA5,FIRSTD,SECONDD,TRBBAR,
     1 BBAR(6),EG,EK,PR,DET,DAMAGEF,DAMAGEG,DBULK,VOID,TEMP,
     2 PRESTR(6),DTIME,EEQ,RELE(NPRONY),DSTRESS(6),
     3 HN1(6,6*NPRONY),TAU(NPRONY),HISTR(6*NPRONY),
     4 STRESS(NDI+NSHR),HISTR1(6*NPRONY),AT,DETA,BAR1,DETC,DETTH,
     5 AL_THER,TEMPI
     
      INTEGER NDI,NSHR 
      
      
      ONE=1.D0
      TWO=2.D0
      THREE=3.D0
      
      FIRSTD=BETA1-BETA4*(2.7183**(-BETA5*(BAR1-3)))
      SECONDD=BETA4*BETA5*(2.7183**(-BETA5*(BAR1-3)))


C     DILATATION DUE TO VOID 
      DETC=1.D0+VOID 
      
C     THERMAL EXPANSION      
      DETTH=(1+AL_THER*(TEMP-TEMPI))**3      
   
      TRBBAR=(BBAR(1)+BBAR(2)+BBAR(3))/THREE
      EG=TWO*FIRSTD/DET*DAMAGEF*DAMAGEG
      EK=DBULK/((DETC*DETTH))*(TWO*DET/(DETC*DETTH)-ONE)
      PR=DBULK/((DETC*DETTH))*(DET/(DETC*DETTH)-ONE)      


      DO K1=1,NDI
        DSTRESS(K1)=EG*(BBAR(K1)-TRBBAR)+PR
      END DO
      DO K1=NDI+1,NDI+NSHR
        DSTRESS(K1)=EG*BBAR(K1)
      END DO    
      
      
      K=1
      DO I=1,NPRONY
        DO J=1,NDI+NSHR  
        HN1(J,I)=2.7183**(-DETA/TAU(I))*HISTR(K)
     1   +(DSTRESS(J)-PRESTR(J))/DETA*TAU(I)
     2   *(1-2.7183**(-DETA/TAU(I)))
        K=K+1
        END DO        
      ENDDO
     
       K=1
      DO I=1,NPRONY
        DO J=1,NDI+NSHR      
        HISTR1(K)=HN1(J,I)
        K=K+1
        END DO
      ENDDO     
      

      DO K1=1,NDI+NSHR
        STRESS(K1)=EEQ*DSTRESS(K1)
      ENDDO
              
        
      DO K1=1,NPRONY
        DO K2=1,NDI+NSHR  
        STRESS(K2)=STRESS(K2)+RELE(K1)*HN1(K2,K1) 
        ENDDO        
      ENDDO
   
      RETURN
      END


      
C-----------------------------------------------------------------------
      SUBROUTINE LOADTYPE(BOCTA,PREOCTA,PREDG,CYCN,BMAX,KINC,AA1,
     1 XIR1,LOADCASE,PBMIN,RELOAD,PDAMAGEF,BMIN,XIRU)
C
C     SORT LOADTYPE
C
C-----------------------------------------------------------------------


      DOUBLE PRECISION BOCTA,CYCN,BMAX,AA1,PREOCTA,PREDG,
     1 XIR1,PBMIN,RELOAD,PDAMAGEF,BMIN,XIRU
     
      INTEGER LOADCASE

      IF (BOCTA .LT. PREOCTA .AND. PREDG .LT. 1.D0) THEN
      CYCN=CYCN+1.D0 
      ENDIF


C     MAXIMUM 
      IF (CYCN .EQ. 1.D0) THEN
      BMAX= PREOCTA
      ENDIF
      XIR1=BMAX


C     STATE VARIABLE 2

      IF (KINC .GT. 1.D0 .AND. CYCN .GE. 1.D0) THEN
          
      AA1=BOCTA/XIR1
      ELSE
      AA1=1
      ENDIF      

C     LOADING

      IF (CYCN .LT. 1.D0) THEN
      LOADCASE=1  
      PBMIN=1.D0
      RELOAD=0.D0    
      PDAMAGEF=1.D0  
      ENDIF
      
      
C     UNLOADING
      IF (CYCN .GE. 1.D0 .AND. BOCTA .LT.  PREOCTA 
     1 .AND. RELOAD .NE. 1.D0) THEN 
      LOADCASE=2
      ENDIF
      

C     MINIMUM DURING UNLOADING

      IF ( CYCN .GE. 1.D0 .AND. LOADCASE .EQ. 2) THEN  
      IF ( BOCTA  .LE.  PREOCTA) THEN
      BMIN=BOCTA      
      ENDIF
  
      ENDIF
      

      XIRU=BMIN     
      
      IF (BMIN .GE. PBMIN ) THEN
      XIRU=PBMIN
      ENDIF

C     RELOADING
      IF (CYCN .GE. 1.D0 .AND. BOCTA .GE.  PREOCTA) THEN
      RELOAD=1.D0
      LOADCASE=3
      ENDIF

C     UNLOADING FROM RELOADING      
      IF (RELOAD .EQ. 1.D0 .AND. BOCTA .LT.  PREOCTA) THEN
      LOADCASE=4
      ENDIF      
      
      
      RETURN
      END
      
      
C-----------------------------------------------------------------------
      SUBROUTINE CALDAMAGEG(VONMISES,VOID,PVOID,TEMP,DAMAGEG,
     1 PDAMAGEG,SIGC,ALPHA_G)
C
C     CAL DAMAGEFUNCTION G DUE TO DEWETTING
C
C-----------------------------------------------------------------------

C     STATEV(6+NPRONY*6+1) : VONMISES
C     STATEV(6+NPRONY*6+4) : VOID
C     STATEV(6+NPRONY*6+5) : PVOID
C     STATEV(6+NPRONY*6+9) : PDAMAGEG

      DOUBLE PRECISION VONMISES,VOID,PVOID,TEMP,DAMAGEG,
     1 PDAMAGEG,SIGC,ALPHA_G 
     
      DAMAGEG=1.D0
 
C     OZUPEK 
C       IF (STATEV(6+NPRONY*6+4) .GE. 0.008) THEN
C        VAL=STATEV(6+NPRONY*6+4)/0.008
C        DAMAGEG=1-0.28*LOG10(VAL)           
C      ENDIF

      IF (VONMISES .GE. SIGC .AND. PVOID .LE. 10E-10) THEN
        
        PVOID=VOID
        
      ENDIF      
 
      IF (VONMISES .GE. SIGC) THEN

        VAL=VOID/PVOID
        
        IF (VAL .LE. 1.D0) THEN
        VAL=1.D0
        ENDIF
        
        
        DAMAGEG=1-ALPHA_G*LOG10(VAL)      
        
      ENDIF

      IF (DAMAGEG .GE. PDAMAGEG) THEN
      DAMAGEG=PDAMAGEG
      ENDIF      

      IF (DAMAGEG .LE. 0.2D0) THEN
        DAMAGEG=0.2D0
      ENDIF

      RETURN
      END          
      


C-----------------------------------------------------------------------
      SUBROUTINE CALINVARIANT(BBAR,BAR1,BAR2,BOCTA,NSHR)
C
C     CAL INVARIANT OF STRAIN TENSOR
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION BAR1, BAR2, BBAR(3+NSHR), BOCTA, AAA

C     BAR1 : FIRST INVARIANT
C     BAR2 : SECOND INVARIANT     

      IF(NSHR.EQ.3) THEN
       BAR1=BBAR(1)+BBAR(2)+BBAR(3)
       BAR2=BBAR(1)*BBAR(2)+BBAR(2)*BBAR(3)+BBAR(3)*BBAR(1)
     1      -BBAR(4)**2-BBAR(5)**2-BBAR(6)**2     
      ELSE
      BAR1=BBAR(1)+BBAR(2)+BBAR(3)
      BAR2=BBAR(1)*BBAR(2)+BBAR(2)*BBAR(3)+BBAR(3)*BBAR(1)
     1      -BBAR(4)**2     
      END IF
      
      AAA=2.D0*BAR1**2-6.D0*BAR2
         
      IF (AAA .LE. 1E-10)THEN
      BOCTA=0.D0
      ELSE
      
     
C     BOCTA : OCTEHEDRAL SHEAR STRAIN
      BOCTA=1.D0/6.D0*((2.D0*BAR1**2-6.D0*BAR2)**(0.5D0))
      
      END IF     
 
      RETURN
      END      
      
      

      
C-----------------------------------------------------------------------
      SUBROUTINE CALDAMAGEF(VAL,LOADCASE,CYCN,RELN,XIR1,XIRU,
     1 BOCTA,DAMAGEF,DERDF,
     2 D_UN_UP_S, D_UN_LB, D_UN_LB_S, D_UN_P1, D_UN_P2,D_UN_P3,
     3 D_UN_P4, D_UN_P5, D_UN_P6, D_RE_UB, D_RE_UP_S,D_RE_LB,D_RE_LB_S,
     4 D_RE_P1,D_RE_P2,D_RE_P3,D_RE_P4,D_RE_P5,D_RE_P6,D_RE_P7,
     5 D_RE_P8,D_UPPER)
C
C     CALCULATE CYCLIC LOADING DAMAGE FUNCTION F
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION VAL, DAMAGEF,DERDF, CYCN, RELN, DAMAAGEFU, DERFU,
     1 DAMAGEFR, DERFR, SVAR3, DSVAR3,XIR1,XIRU,BOCTA,
     2 D_UN_UP_S, D_UN_LB, D_UN_LB_S, D_UN_P1, D_UN_P2,D_UN_P3,
     3 D_UN_P4, D_UN_P5, D_UN_P6, D_RE_UB, D_RE_UP_S,D_RE_LB,D_RE_LB_S,
     4 D_RE_P1,D_RE_P2,D_RE_P3,D_RE_P4,D_RE_P5,D_RE_P6,D_RE_P7,
     5 D_RE_P8,D_UPPER
      INTEGER LOADCASE


      IF (VAL .GE. 1.D0) THEN
      DAMAGEFU=1.D0
      DERFU=D_UN_UP_S
      ELSEIF (VAL .LE. D_UN_LB) THEN
      DAMAGEFU=1.D0
      DERFU=D_UN_LB_S    
      ELSE
      DAMAGEFU=D_UN_P1*VAL**5+D_UN_P2*VAL**4+D_UN_P3*VAL**3
     1          +D_UN_P4*VAL**2+D_UN_P5*VAL+D_UN_P6
     
      DERFU=D_UN_P1*5*VAL**4+D_UN_P2*4*VAL**3+D_UN_P3*3*VAL**2
     1          +D_UN_P4*2*VAL+D_UN_P5 
      ENDIF
      

      IF (VAL .GE. D_RE_UB) THEN
      DAMAGEFR=1.D0
      DERFR=D_RE_UP_S
      ELSEIF (VAL .LE. D_RE_LB) THEN
      DAMAGEFR=1.D0
      DERFR=D_RE_LB_S     
      ELSE
      DAMAGEFR=D_RE_P1*VAL**7+D_RE_P2*VAL**6+D_RE_P3*VAL**5
     1              +D_RE_P4*VAL**4+D_RE_P5*VAL**3+D_RE_P6*VAL**2
     2              +D_RE_P7*VAL+D_RE_P8  
      DERFR=D_RE_P1*7*VAL**6+D_RE_P2*6*VAL**5+D_RE_P3*5*VAL**4
     1              +D_RE_P4*4*VAL**3+D_RE_P5*3*VAL**2+D_RE_P6*2*VAL
     2              +D_RE_P7
      ENDIF
      
       
      IF (LOADCASE .EQ. 3 .OR. LOADCASE .EQ. 4) THEN
      
      SVAR3=(BOCTA-XIRU)/(XIR1-XIRU)
      DSVAR3=(XIR1)/(XIR1-XIRU)

      ELSE
      SVAR3=0.D0
      DSVAR3=0.D0
      ENDIF
      

      IF (SVAR3 .LE. 0.D0) THEN
      SVAR3=0.D0
      ENDIF       
       
      IF (LOADCASE .EQ. 1) THEN
      DAMAGEF=1.D0     
      DERDF=0.D0     

      ELSEIF (LOADCASE .EQ. 2 ) THEN

      DAMAGEF=DAMAGEFU     
      DERDF=DERFU          
      
      ELSEIF (LOADCASE .EQ. 3 .OR. LOADCASE .EQ. 4 ) THEN

      DAMAGEF=DAMAGEFU+SVAR3*(DAMAGEFR-DAMAGEFU)      
      DERDF=DERFU+SVAR3*(DERFR-DERFU)
     1 +DSVAR3*(DAMAGEFR-DAMAGEFU)
      
      ELSE
      DAMAGEF=1.D0
      DERDF=0.D0
      ENDIF
  
      IF (VAL .GE. D_UPPER) THEN
      LOADCASE=1
      CYCN=0.D0
      RELN=0.D0
      DAMAGEF=1.D0
      ENDIF  
      
      IF (DAMAGEF .GE. 1.D0) THEN
      DAMAGEF=1.D0
      ENDIF
C
      RETURN
      END 
      