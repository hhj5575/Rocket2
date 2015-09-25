SUBROUTINE SPBS(nDim, OP, CurrentTime, DeltaTime, nBNode, nBFace, Bxyz, Bf2n, bFlag, InVars, OutVars, Prop_Dens, Ref_Pressure, IgniteAPN)
!--------------------------------------------!
! Description:                               !
!   Solid Propellant Burning Simulation      !
!--------------------------------------------!
! Burning Model:          Errosive Burning:  !
!   1 : APN Model           L&R Model        !
!   2 : ZN  Model           (Optional)       !
!--------------------------------------------!
!                           Written by: LCS  !
!--------------------------------------------!
! Main Operation Subroutines
USE SPBS_INITIAL_MOD, ONLY : INITIAL
USE SPBS_SOLVER_MOD , ONLY : SOLVER
USE SPBS_REMESH_MOD , ONLY : REMESH
USE SPBS_POST_MOD   , ONLY : POST
USE SPBS_REGION_MOD , ONLY : FINISH
! Integration Interface Subroutines
USE SPBS_INTERFACE_MOD, ONLY : INTERFACESET
USE SPBS_INTERFACE_MOD, ONLY : INTERFACEIN
USE SPBS_INTERFACE_MOD, ONLY : INTERFACEOUT
! In or Out Variables
INTEGER, INTENT(in   ) :: nDim, OP
REAL(8), INTENT(in   ) :: CurrentTime
REAL(8), INTENT(in   ) :: DeltaTime
INTEGER, INTENT(in   ) :: nBNode, nBFace
REAL(8), INTENT(in   ) :: Bxyz(nDim,nBNode)
INTEGER, INTENT(in   ) :: Bf2n(nDim,nBFace)
INTEGER, INTENT(inout) :: bFlag(nBFace)
REAL(8), INTENT(inout) :: InVars(3,nBFace) ! 1: Pressure, 2: Temperature, 3: Tagential mass velocity
REAL(8), INTENT(  out) :: OutVars(2,nBFace) ! 1: Burning rate, 2: Adiabatic flame temperature
REAL(8), INTENT(  out) :: Prop_Dens
REAL(8), INTENT(in   ) :: Ref_Pressure
LOGICAL, INTENT(  out) :: IgniteAPN

! SPBS Module Operation
SELECT CASE(OP)
CASE(1)
  ! Operation 1 : Initialize Region
  CALL INTERFACESET(nDim, nBNode, nBFace, Bxyz, Bf2n)
  CALL INITIAL(nDim, nBFace, bFlag, InVars)
  CALL INTERFACEOUT(nBFace, bFlag, OutVars, Prop_Dens, IgniteAPN)
CASE(2)
  ! Operation 2 : Solve Region
  CALL INTERFACEIN(nDim, nBNode, nBFace, Bxyz, InVars, Ref_Pressure)
  CALL SOLVER(CurrentTime, DeltaTime)
  CALL INTERFACEOUT(nBFace, bFlag, OutVars, Prop_Dens, IgniteAPN)
CASE(3)
  ! Operation 3 : Remesh Region
  CALL REMESH(nDim, nBNode, nBFace, Bxyz, Bf2n, bFlag)
  CALL INTERFACEOUT(nBFace, bFlag, OutVars, Prop_Dens, IgniteAPN)
CASE(4)
  ! Operation 4 : Post Region
  CALL POST(CurrentTime)
CASE(5)
  ! Operation 5 : Finish Region
  CALL FINISH
END SELECT

END SUBROUTINE SPBS
