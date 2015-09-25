MODULE SPBS_BURN_MOD

USE SPBS_TYPE_MOD
IMPLICIT NONE; PRIVATE
PUBLIC BURN_APN
PUBLIC BURN_ZN
PUBLIC BURN_LR
PUBLIC CHECK_BRATE

CONTAINS

SUBROUTINE BURN_APN(pProp, pBurn)

! In or Out Variables
TYPE(t_prop), INTENT(in   ) :: pProp
TYPE(t_Burn), INTENT(inout) :: pBurn
! Local Variables
REAL(8), PARAMETER :: Pref = 6894733.26d0

pBurn%bRate = pProp%a_APN * ( pBurn%P_fluid / Pref )**pProp%n_APN
pBurn%mFlux = pBurn%bRate * pProp%rhoc
pBurn%fTemp = pProp%T_star0

END SUBROUTINE

SUBROUTINE BURN_ZN(pLGrid, pProp, pBurn, dT)

! In or Out Variables
TYPE(t_LGrid), INTENT(in   ) :: pLGrid
TYPE(t_prop) , INTENT(in   ) :: pProp
TYPE(t_Burn) , INTENT(inout) :: pBurn
REAL(8)      , INTENT(in   ) :: dT
! Local Variables

END SUBROUTINE

SUBROUTINE BURN_LR(indLR, pSGrid, pProp, pBurn, BFace)

! In or Out Variables
INTEGER      , INTENT(in   ) :: indLR
TYPE(t_SGrid), INTENT(in   ) :: pSGrid
TYPE(t_prop) , INTENT(in   ) :: pProp
TYPE(t_Burn) , INTENT(inout) :: pBurn
INTEGER      , INTENT(in   ) :: BFace
! Local Variables
INTEGER :: iter, Com
INTEGER :: Node(pSGrid%nDim)
REAL(8) :: Cen(pSGrid%nDim)
REAL(8) :: L, Error
REAL(8) :: modified_alpha
REAL(8) :: bRate_temp, bRate_erosive
INTEGER, PARAMETER :: MaxIter = 500
REAL(8), PARAMETER :: tol = 1d-12

if( pBurn%G_fluid < 1d-15 ) return

Node(:) = pSGrid%Bf2n(:,BFace)
DO Com = 1, pSGrid%nDim
  Cen(Com) = SUM( pSGrid%Bxyz(Com,Node(:)) )/pSGrid%nDim
END DO

! Calculate L
SELECT CASE(indLR)
CASE(1)
  L = Cen(pSGrid%AxisCom) - pSGrid%Ref_Location
CASE(2)
  L = DSQRT( SUM( Cen(:)**2 ) - Cen(pSGrid%AxisCom)**2 )
CASE(3)
  L = DSQRT( SUM( Cen(:)**2 ) - Cen(pSGrid%AxisCom)**2 )
  L = 0.9d0+0.189d0*L*(1d0+0.043d0*L*(1d0+0.023d0*L))
CASE(4)
  L = DSQRT( SUM( Cen(:)**2 ) - Cen(pSGrid%AxisCom)**2 )
  L = 0.9d0+0.189d0*L*(1d0+0.043d0*L*(1d0+0.023d0*L))
  L = 0.5d0*L**5d0
END SELECT

! Modify alpha with step function
modified_alpha = pProp%alpha
IF( pBurn%G_fluid >= pProp%StepInd ) then
  modified_alpha = pProp%alpha * pProp%StepMag
END IF

! Erosive Burning Rate Initial Guess
bRate_erosive = 0d0

! Calculate Erosive Burning Rate
DO iter = 1, MaxIter
  bRate_temp = modified_alpha*(pBurn%G_fluid**0.8d0)*(L**(-0.2d0))*DEXP(-pProp%beta*pProp%rhoc*(pBurn%bRate+bRate_erosive)/pBurn%G_fluid)
  Error = DABS( bRate_temp - bRate_erosive )
  bRate_erosive = bRate_temp 
  IF( Error < tol ) EXIT
END DO
IF( iter == MaxIter+1 ) THEN
  WRITE(*,*) 'SPBS >> Erosive burning divergence'
  WRITE(*,*) 'SPBS >> Error =', Error
  STOP
END IF

! Update Burning Rate
pBurn%bRate = pBurn%bRate + bRate_erosive

END SUBROUTINE

FUNCTION CHECK_BRATE(pProp, pBurn)

! In or Out Variables
TYPE(t_prop), INTENT(in   ) :: pProp
TYPE(t_Burn), INTENT(inout) :: pBurn
LOGICAL :: CHECK_BRATE
! Local Variables
REAL(8) :: bRate
REAL(8), PARAMETER :: Pref = 6894733.26d0

! Check Dynamic Burning Model And APN Model Burning Rate
bRate = pProp%a_APN * ( pBurn%P_fluid / Pref )**pProp%n_APN
IF( DABS( pBurn%bRate - bRate ) <= 1d-12 ) THEN
  CHECK_BRATE = .TRUE.
ELSE
  CHECK_BRATE = .FALSE.
END IF

END FUNCTION

END MODULE SPBS_BURN_MOD
