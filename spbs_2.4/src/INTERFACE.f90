MODULE SPBS_INTERFACE_MOD

USE SPBS_TYPE_MOD
USE SPBS_REGION_MOD, ONLY : Region
IMPLICIT NONE; PRIVATE
PUBLIC INTERFACESET
PUBLIC INTERFACEIN
PUBLIC INTERFACEOUT

CONTAINS

SUBROUTINE INTERFACESET(nDim, nBNode, nBFace, Bxyz, Bf2n)

! In or Out Variables
INTEGER, INTENT(in   ) :: nDim
INTEGER, INTENT(in   ) :: nBNode, nBFace
REAL(8), INTENT(in   ) :: Bxyz(nDim,nBNode)
INTEGER, INTENT(in   ) :: Bf2n(nDim,nBFace)
! Local Variables
TYPE(t_SGrid) , POINTER :: pSGrid
TYPE(t_Nozzle), POINTER :: pNozzle
INTEGER :: BNode

! Set Surface Grid Pointer
pSGrid => Region%SGrid

! Set Nozzle Pointer
pNozzle => Region%Nozzle

! Set Dimension, Number of Node && Face
pSGrid%nDim = nDim
pSGrid%nBNode = nBNode
pSGrid%nBFace = nBFace

! Allocate Coordinate && Grid Connectivity
ALLOCATE( pSGrid%Bxyz(nDim,nBNode) )
ALLOCATE( pSGrid%Bf2n(nDim,nBFace) )

! Set Coordinate && Grid Connectivity
pSGrid%Bxyz(:,:) = Bxyz(:,:)
pSGrid%Bf2n(:,:) = Bf2n(:,:)

! Set Reference Location
IF( pSGrid%nDim == 2 ) THEN
  pSGrid%AxisCom = 1
ELSE IF( pSGrid%nDim == 3 ) THEN
  pSGrid%AxisCom = 2
END IF
pSGrid%Ref_Location = pSGrid%Bxyz(pSGrid%AxisCom,1)
DO BNode = 2, pSGrid%nBNode
  IF( pSGrid%Bxyz(pSGrid%AxisCom,BNode) < pSGrid%Ref_Location ) pSGrid%Ref_Location = pSGrid%Bxyz(pSGrid%AxisCom,BNode)
END DO

END SUBROUTINE

SUBROUTINE INTERFACEIN(nDim, nBNode, nBFace, Bxyz, InVars, Pressure)

! In or Out Variables
INTEGER, INTENT(in   ) :: nDim
INTEGER, INTENT(in   ) :: nBNode, nBFace
REAL(8), INTENT(in   ) :: Bxyz(nDim,nBNode)
REAL(8), INTENT(in   ) :: InVars(3,nBFace)
REAL(8), INTENT(in   ) :: Pressure
! Local Variables
TYPE(t_Input) , POINTER :: pInput
TYPE(t_SGrid) , POINTER :: pSGrid
TYPE(t_Prop)  , POINTER :: pProp
TYPE(t_Burn)  , POINTER :: pBurn
TYPE(t_Nozzle), POINTER :: pNozzle
INTEGER :: BFace, SubBFace
INTEGER :: n1, n2
REAL(8) :: Cen(2), SubCen(2)
REAL(8) :: Sign_Loc1, Sign_Loc2

! Set Input Pointer
pInput => Region%Input

! Set Surface Grid Pointer
pSGrid => Region%SGrid

! Set Propellant Pointer
pProp => Region%Prop

! Set Nozzle Pointer
pNozzle => Region%Nozzle

! Coordinate Change
pSGrid%Bxyz(:,:) = Bxyz(:,:)

! Set Fluid Boundary Data
! Including Non Propellant Region
! (We Need Fluid Data on Axis)
DO BFace = 1, pSGrid%nBFace
  pBurn => Region%Burn(BFace)
  pBurn%P_fluid = InVars(1,BFace)
  pBurn%T_fluid = InVars(2,BFace)
  pBurn%G_fluid = InVars(3,BFace)
END DO

! Reset Mass Velocity
! Using Center Value
IF( nDim == 2 ) THEN
  DO BFace = 1, pSGrid%nBFace
    pBurn => Region%Burn(BFace)
    SELECT CASE(pBurn%Ignited)
    CASE(0,1,2)
      SELECT CASE(pInput%indk_ero)
      CASE(1) ! weighting with center line
        n1 = pSGrid%Bf2n(1,BFace)
        n2 = pSGrid%Bf2n(2,BFace)
        Cen(:) = 0.5d0*( pSGrid%Bxyz(:,n1) + pSGrid%Bxyz(:,n2) )
        DO SubBFace = 1, pSGrid%nBFace
          n1 = pSGrid%Bf2n(1,SubBFace)
          n2 = pSGrid%Bf2n(2,SubBFace)
          SubCen(:) = 0.5d0*( pSGrid%Bxyz(:,n1) + pSGrid%Bxyz(:,n2) )
          Sign_Loc1 = Cen(1) - pSGrid%Bxyz(1,n1)
          Sign_Loc2 = Cen(1) - pSGrid%Bxyz(1,n2)
          IF( SubCen(2) < 1d-12 .AND. Sign_Loc1*Sign_Loc2 <= 0d0 ) THEN
            pBurn%G_fluid = (1d0-pProp%k_ero)*InVars(3,BFace) + pProp%k_ero*InVars(3,SubBFace)
            EXIT
          END IF
        END DO
      CASE(2) ! simply weighting
        pBurn%G_fluid = pProp%k_ero*pBurn%G_fluid
      END SELECT
    END SELECT
  END DO
END IF

! Set Reference Pressure for Nozzle Ablation
pNozzle%P_Ref = Pressure

END SUBROUTINE

SUBROUTINE INTERFACEOUT(nBFace, bFlag, OutVars, Prop_Dens, IgniteAPN)

! In or Out Variables
INTEGER, INTENT(in   ) :: nBFace
INTEGER, INTENT(  out) :: bFlag(nBFace)
REAL(8), INTENT(  out) :: OutVars(2,nBFace)
REAL(8), INTENT(  out) :: Prop_Dens
LOGICAL, INTENT(  out) :: IgniteAPN
! Local Variables
TYPE(t_SGrid) , POINTER :: pSGrid
TYPE(t_Burn)  , POINTER :: pBurn
TYPE(t_Prop)  , POINTER :: pProp
TYPE(t_Nozzle), POINTER :: pNozzle
INTEGER :: BFace, iDim
REAL(8) :: Cen(Region%SGrid%nDim)
REAL(8) :: Diameter2, Diameter2_Ref

! Set Pointer
pSGrid  => Region%SGrid
pProp   => Region%Prop 
pNozzle => Region%Nozzle

! Set Nozzle Ablation Range and Maximum Speed Location
pNozzle%Location_Range(1) = 2d0**1023
pNozzle%Location_Range(2) = -2d0**1023
Diameter2_Ref = 2d0**1023
DO BFace = 1, pSGrid%nBFace
  pBurn => Region%Burn(BFace)
  IF( pBurn%Ignited == -2 ) THEN
    DO iDim = 1, pSGrid%nDim
      Cen(iDim)= SUM( pSGrid%Bxyz(iDim,pSGrid%Bf2n(:,BFace)) ) / DBLE(pSGrid%nDim)
    END DO
    pNozzle%Location_Range(1) = MIN( pNozzle%Location_Range(1), Cen(pSGrid%AxisCom) )
    pNozzle%Location_Range(2) = MAX( pNozzle%Location_Range(2), Cen(pSGrid%AxisCom) )
    Diameter2 = SUM( Cen(:)**2 ) - Cen(pSGrid%AxisCom)**2
    IF( Diameter2 < Diameter2_Ref ) THEN
      Diameter2_Ref = Diameter2
      pNozzle%Max_Loc = Cen(pSGrid%AxisCom)
    END IF
  END IF
END DO

! Set bRate / fTemp / bFlag / IgniteAPN
IgniteAPN = .true.
DO BFace = 1, pSGrid%nBFace
  pBurn => Region%Burn(BFace)
  ! Set OutVars
  IF( pBurn%Ignited == -2 ) THEN ! Nozzle
    bFlag(BFace) = -2
    OutVars(1,BFace) = ABLATION(pSGrid,pNozzle,BFace)
    OutVars(2,BFace) = 0d0
  ELSE IF( pBurn%Ignited == -1 ) THEN ! Wall
    bFlag(BFace) = -1
    OutVars(1,BFace) = 0d0
    OutVars(2,BFace) = 0d0
  ELSE IF( pBurn%Ignited == 0 ) THEN ! Propellant (not ignited)
    bFlag(BFace) = 0
    OutVars(1,BFace) = 0d0
    OutVars(2,BFace) = 0d0
    IgniteAPN = .false.
  ELSE IF( pBurn%Ignited == 1 ) THEN ! Propellant (APN model)
    bFlag(BFace) = 1
    OutVars(1,BFace) = pBurn%bRate
    OutVars(2,BFace) = pBurn%fTemp
  ELSE IF( pBurn%Ignited == 2 ) THEN ! Propellant (ZN model)
    bFlag(BFace) = 1
    OutVars(1,BFace) = pBurn%bRate
    OutVars(2,BFace) = pBurn%fTemp
    IgniteAPN = .false.
  END IF
END DO

! Set Propellant Density
Prop_Dens = pProp%rhoc

END SUBROUTINE

FUNCTION ABLATION(pSGrid,pNozzle,BFace)

! In or Out Variables
TYPE(t_SGrid) , INTENT(in   ) :: pSGrid
TYPE(t_Nozzle), INTENT(in   ) :: pNozzle
INTEGER       , INTENT(in   ) :: BFace
REAL(8) :: ABLATION
! Local Variables
REAL(8) :: Cen_Coord

Cen_Coord = SUM( pSGrid%Bxyz(pSGrid%AxisCom,pSGrid%Bf2n(:,BFace)) ) / DBLE(pSGrid%nDim)

IF( Cen_Coord >= pNozzle%Location_Range(1) - 1d-12 .AND. Cen_Coord <= pNozzle%Max_Loc ) THEN
  ABLATION = pNozzle%Ablation_Speed * ( Cen_Coord - pNozzle%Location_Range(1) ) / (pNozzle%Max_Loc - pNozzle%Location_Range(1) )
ELSE IF( Cen_Coord <= pNozzle%Location_Range(2) + 1d-12 .AND. Cen_Coord > pNozzle%Max_Loc ) THEN
  ABLATION = pNozzle%Ablation_Speed * ( Cen_Coord - pNozzle%Location_Range(2) ) / (pNozzle%Max_Loc - pNozzle%Location_Range(2) )
END IF

END FUNCTION

END MODULE SPBS_INTERFACE_MOD
