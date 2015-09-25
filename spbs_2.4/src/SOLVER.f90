MODULE SPBS_SOLVER_MOD

USE SPBS_TYPE_MOD
USE SPBS_REGION_MOD, ONLY : Region
IMPLICIT NONE; PRIVATE
PUBLIC SOLVER

CONTAINS

SUBROUTINE SOLVER(CurrentTime, DeltaTime)

USE SPBS_BURN_MOD, ONLY : BURN_APN
USE SPBS_BURN_MOD, ONLY : BURN_ZN
USE SPBS_BURN_MOD, ONLY : BURN_LR
USE SPBS_BURN_MOD, ONLY : CHECK_BRATE
! In or Out Variables
REAL(8), INTENT(in   ) :: CurrentTime
REAL(8), INTENT(in   ) :: DeltaTime
! Local Variables
TYPE(t_Input) , POINTER :: pInput
TYPE(t_SGrid) , POINTER :: pSGrid
TYPE(t_LGrid) , POINTER :: pLGrid
TYPE(t_Prop)  , POINTER :: pProp
TYPE(t_Burn)  , POINTER :: pBurn
TYPE(t_Nozzle), POINTER :: pNozzle
INTEGER :: iter
INTEGER :: BFace
REAL(8) :: Time, dTmin
REAL(8) :: TargetTime

! Set Pointer
pInput => Region%Input
pSGrid => Region%SGrid
pLGrid => Region%LGrid
pProp  => Region%Prop
pNozzle => Region%Nozzle

! Boundary Loop
DO BFace = 1, pSGrid%nBFace

  ! Set Pointer
  pBurn => Region%Burn(BFace)

  ! Non Propellant Surface
  IF( pBurn%Ignited < 0 ) CYCLE
 
  ! Set Time && Target Time
  Time       = CurrentTime
  TargetTime = CurrentTime + DeltaTime

  ! Calculate Film Coefficient
  CALL CALCFILM(pInput, pSGrid, pProp, pBurn, BFace)

  DO iter = 1, pInput%IterMax

    IF( pBurn%Ignited == 0 ) THEN
      ! Calculate Time Step
      dTmin = pBurn%dTime
      IF( Time + dTmin >= TargetTime ) dTmin = TargetTime - Time
      ! 1-D Heat Equation
      CALL CONVECTION(pLGrid, pProp, pBurn, dTmin)
      ! Verify Propellant Ignition
      IF( pBurn%T(1) >= pProp%T_ignite ) THEN
        SELECT CASE(pInput%iBurn)
        CASE(1)
          pProp%nIgnited_APN = pProp%nIgnited_APN + 1
          pBurn%Ignited = 1
        CASE(2)
          pProp%nIgnited_ZN = pProp%nIgnited_ZN + 1
          pBurn%Ignited = 2
        END SELECT
      END IF
    END IF
      
    IF( pBurn%Ignited == 1 ) THEN
      ! Calculate Time Step
      dTmin = TargetTime - Time
      ! APN Model
      CALL BURN_APN(pProp, pBurn)
      ! L&R Model
      IF( pInput%indLR /= 0 ) CALL BURN_LR(pInput%indLR, pSGrid, pProp, pBurn, BFace)

    ELSE IF( pBurn%Ignited == 2 ) THEN
      ! Calculate Time Step
      dTmin = TargetTime - Time
      ! ZN Model
      CALL BURN_ZN(pLGrid, pProp, pBurn, dTmin)
      IF( CHECK_BRATE(pProp, pBurn) ) THEN
        pProp%nIgnited_ZN  = pProp%nIgnited_ZN  - 1
        pProp%nIgnited_APN = pProp%nIgnited_APN + 1
        pBurn%Ignited = 1
      END IF

    END IF

    ! Time Update
    Time = Time + dTmin
    IF( Time >= TargetTime-1d-15 ) EXIT

  END DO

END DO

WRITE(*,*)
WRITE(*,*) 'SPBS >> Burning Region Result Overview'
WRITE(*,*) 'SPBS >> Propellant Surface    :', pProp%nProp
IF( pInput%indLR == 0 ) THEN
  WRITE(*,*) 'SPBS >> Ignited with APN      :', pProp%nIgnited_APN
ELSE
  WRITE(*,*) 'SPBS >> Ignited with APN(L&R) :', pProp%nIgnited_APN
END IF
WRITE(*,*) 'SPBS >> Ignited with ZN       :', pProp%nIgnited_ZN
WRITE(*,*) 'SPBS >> Max Burning Rate      :', MAXVAL(Region%Burn(:)%bRate)

! Calculate Ablation Maximum Speed
IF( CurrentTime + DeltaTime >= pNozzle%Ablation_Time ) THEN
  pNozzle%Ablation_Speed = 2d0 * pNozzle%Ablation_Coeff * (pNozzle%P_ref/6894.73326d0)**0.8d0
ELSE
  pNozzle%Ablation_Speed = 0d0
END IF
WRITE(*,*) 'SPBS >> Max Ablation Speed    :', pNozzle%Ablation_Speed

END SUBROUTINE

SUBROUTINE CALCFILM(pInput, pSGrid, pProp, pBurn, BFace)

! In or Out Variables
TYPE(t_Input), INTENT(in   ) :: pInput
TYPE(t_SGrid), INTENT(in   ) :: pSGrid
TYPE(t_Prop) , INTENT(inout) :: pProp
TYPE(t_Burn) , INTENT(inout) :: pBurn
INTEGER      , INTENT(in   ) :: BFace
! Local Variable
REAL(8), PARAMETER :: inv3 = 1d0/3d0
REAL(8), PARAMETER :: inv45 = 4d0/5d0
INTEGER :: Node(pSGrid%nDim), Com
REAL(8) :: Cen(pSGrid%nDim), Dia, Length, Re

IF( pInput%iFilm == 1 ) THEN
  pBurn%film = pProp%film_const
  RETURN
END IF

Node(:) = pSGrid%Bf2n(:,BFace)
DO Com = 1, pSGrid%nDim
  Cen(Com) = SUM( pSGrid%Bxyz(Com,Node(:)) )/pSGrid%nDim
END DO

IF( pInput%iFilm == 2 ) THEN
  Dia = DSQRT( SUM( Cen(:)**2 ) - Cen(pSGrid%AxisCom)**2 )
  Re = pBurn%G_fluid*Dia/pProp%mu
  pBurn%film = (pProp%k/Dia)*0.023d0*pProp%Pr**0.4d0*Re**inv45
ELSE IF( pInput%iFilm == 3 ) THEN
  Length = Cen(pSGrid%AxisCom) - pSGrid%Ref_Location
  Re = pBurn%G_fluid*Dia/pProp%mu
  pBurn%film = (pProp%k/Length)*0.029d0*pProp%Pr**inv3*Re**inv45
END IF

END SUBROUTINE

SUBROUTINE CONVECTION(pLGrid,pProp,pBurn,dT)

! In or Out Variables
TYPE(t_LGrid), INTENT(in   ) :: pLGrid
TYPE(t_Prop) , INTENT(in   ) :: pProp
TYPE(t_Burn) , INTENT(inout) :: pBurn
REAL(8)      , INTENT(in   ) :: dT
! Local Variables
INTEGER :: Node
REAL(8) :: First
REAL(8) :: Second
REAL(8) :: fsprime
REAL(8) :: coe, rhs, add

! Set T0
pBurn%T0(:) = pBurn%T(:)

! Update Temperature Profile
DO Node = 2, pLGrid%nNode-1
  First = ( pBurn%T0(Node+1) - pBurn%T0(Node-1) ) / ( 2d0 * pLGrid%dz )
  Second = ( pBurn%T0(Node+1) - 2d0 * pBurn%T0(Node) + pBurn%T0(Node-1) ) * pLGrid%dzsqinv
  pBurn%T(Node) = pBurn%T0(Node) + dT * ( pProp%alfac * ( pLGrid%dzdx(Node)**2 * Second + pLGrid%dzdx2(Node) * First ) )
END DO

! Update BC
pBurn%T(pLGrid%nNode) = pProp%T_init
fsprime = - pBurn%film / pProp%lamc
! 2nd Order
!coe =  3d0 + fsprime * 2d0 * pLGrid%dx1
!rhs = -2d0 * pLGrid%dx1 * fsprime * T_fluid
!add =  4d0 * pBurn%T(2) - pBurn%T(3)
! 3rd Order
coe = 11d0 + fsprime * 6d0 * pLGrid%dx1
rhs = -6d0 * pLGrid%dx1 * fsprime * pBurn%T_fluid
add = 18d0 * pBurn%T(2) - 9d0 * pBurn%T(3) + 2d0 * pBurn%T(4)
pBurn%T(1) = ( add - rhs ) / coe

END SUBROUTINE

END MODULE SPBS_SOLVER_MOD
