MODULE SPBS_REMESH_MOD

USE SPBS_TYPE_MOD
USE SPBS_REGION_MOD, ONLY : Region
IMPLICIT NONE; PRIVATE
PUBLIC REMESH

CONTAINS

SUBROUTINE REMESH(nDim, nBNode, nBFace, Bxyz, Bf2n, bFlag)

! In or Out Variables
INTEGER, INTENT(in   ) :: nDim
INTEGER, INTENT(in   ) :: nBNode, nBFace
REAL(8), INTENT(in   ) :: Bxyz(nDim,nBNode)
INTEGER, INTENT(in   ) :: Bf2n(nDim,nBFace)
INTEGER, INTENT(in   ) :: bFlag(nBFace)
! Local Variables
TYPE(t_SGrid), POINTER :: pSGrid
TYPE(t_Prop) , POINTER :: pProp
TYPE(t_Burn) , POINTER :: pBurn
INTEGER :: BFace

! Set Pointer
pSGrid => Region%SGrid
pProp  => Region%Prop

! Check All Surface Ignited with APN Model
DO BFace = 1, pSGrid%nBFace
  pBurn => Region%Burn(BFace)
  SELECT CASE(pBurn%Ignited)
  CASE(0,2)
    WRITE(*,*) 'SPBS >> Error : For remeshing, all surface must be APN burning state'
    STOP
  END SELECT
END DO

! Deallocate Coordinate %% Grid Connectivity
DEALLOCATE( pSGrid%Bxyz )
DEALLOCATE( pSGrid%Bf2n )

! Deallocate Burning Surface
DO BFace = 2, pSGrid%nBFace
  pBurn => Region%Burn(BFace)
  IF( ALLOCATED( pBurn%T  ) ) DEALLOCATE( pBurn%T  )
  IF( ALLOCATED( pBurn%T0 ) ) DEALLOCATE( pBurn%T0 )
END DO
DEALLOCATE( Region%Burn )

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

! Allocate Burning Surface && Set Ignited Flag
pProp%nProp = 0
ALLOCATE( Region%Burn(pSGrid%nBFace) )
DO BFace = 1, pSGrid%nBFace
  pBurn => Region%Burn(BFace)
  IF( bFlag(BFace) == -1 ) THEN
    pBurn%Ignited = -1
  ELSE IF( bFlag(BFace) == -2 ) THEN
    pBurn%Ignited = -2
  ELSE IF( bFlag(BFace) >= 0 ) THEN
    pProp%nProp = pProp%nProp + 1
    pBurn%Ignited = 1
    ALLOCATE( pBurn%T(1) )
    pBurn%T(1) = pProp%T_ignite
  END IF
END DO
pProp%nIgnited_APN = pProp%nProp
pProp%nIgnited_ZN  = 0

! Do not Need Time Step Calculation 

WRITE(*,*) 'SPBS >> Burning Region Remeshing Process'
WRITE(*,*) 'SPBS >> Propellant Surface    :', pProp%nProp
WRITE(*,*) 'SPBS >> Ignited with APN(L&R) :', pProp%nIgnited_APN
WRITE(*,*) 'SPBS >> Ignited with ZN       :', pProp%nIgnited_ZN
WRITE(*,*)

END SUBROUTINE

END MODULE SPBS_REMESH_MOD
