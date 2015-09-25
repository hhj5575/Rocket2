MODULE SPBS_INITIAL_MOD

USE SPBS_TYPE_MOD
USE SPBS_REGION_MOD, ONLY : Region
IMPLICIT NONE; PRIVATE
PUBLIC INITIAL

CONTAINS

SUBROUTINE INITIAL(nDim, nBFace, bFlag, InVars)

USE SPBS_REGION_MOD, ONLY : ALLOC

! In or Out Variables
INTEGER, INTENT(in   ) :: nDim
INTEGER, INTENT(in   ) :: nBFace
INTEGER, INTENT(in   ) :: bFlag(nBFace)
REAL(8), INTENT(in   ) :: InVars(3,nBFace)
! Local Variables
TYPE(t_Input) , POINTER :: pInput
TYPE(t_SGrid) , POINTER :: pSGrid
TYPE(t_LGrid) , POINTER :: pLGrid
TYPE(t_Prop)  , POINTER :: pProp
TYPE(t_Burn)  , POINTER :: pBurn
TYPE(t_Nozzle), POINTER :: pNozzle
INTEGER :: BFace

! Set Pointer
pInput  => Region%Input
pSGrid  => Region%SGrid
pLGrid  => Region%LGrid
pProp   => Region%Prop
pNozzle => Region%Nozzle

! Read Input File
CALL READINPUT(pInput,pProp,pNozzle,nDim)

! Make 1-D Line Grid
CALL MAKELGRID(pLGrid)

! Initialize Region
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
    CALL ALLOC(pLGrid,pBurn)
    CALL INITIALIZE(pInput,pSGrid,pLGrid,pProp,pBurn,BFace,InVars(1,BFace))
  END IF
END DO

IF( pInput%indSB == 0 ) THEN
  pProp%nIgnited_APN = 0
  pProp%nIgnited_ZN  = 0
ELSE IF( pInput%indSB == 1 ) THEN
  pProp%nIgnited_APN = pProp%nProp
  pProp%nIgnited_ZN  = 0
END IF

WRITE(*,*) 'SPBS >> Propellant Surface    :', pProp%nProp
WRITE(*,*) 'SPBS >> Ignited with APN(L&R) :', pProp%nIgnited_APN
WRITE(*,*) 'SPBS >> Ignited with ZN       :', pProp%nIgnited_ZN
WRITE(*,*)

! Initialize Nozzle Ablation Speed
pNozzle%Ablation_Speed = 0d0

END SUBROUTINE

SUBROUTINE READINPUT(pInput, pProp, pNozzle, nDim)

! In or Out Variables
TYPE(t_Input) , INTENT(  out) :: pInput
TYPE(t_Prop)  , INTENT(  out) :: pProp
TYPE(t_Nozzle), INTENT(inout) :: pNozzle
INTEGER       , INTENT(in   ) :: nDim
! Local Variables
INTEGER :: io
CHARACTER(99) :: Dummy
CHARACTER(99) :: Fname
LOGICAL :: Exist

! Find Input File
Fname = "./input/burn.inp"
INQUIRE(File=Fname, Exist=Exist)
IF( .not. Exist ) THEN
  WRITE(*,*) 'SPBS >> Error : Cannot find input file < burn.inp >'
  STOP
END IF

! Read Input File
OPEN(newunit=io,File=Fname)
READ(io,*) ! SPBS Input File
READ(io,*) ! Control
READ(io,*) Dummy, Dummy, pInput%indSB
READ(io,*) Dummy, Dummy, pInput%iFilm
READ(io,*) Dummy, Dummy, pInput%iBurn
READ(io,*) Dummy, Dummy, pInput%indLR
READ(io,*) Dummy, Dummy, pInput%indk_ero
READ(io,*) ! Solid Propellant
READ(io,*) Dummy, Dummy, pProp%rhoc
READ(io,*) Dummy, Dummy, pProp%MW
READ(io,*) Dummy, Dummy, pProp%alfac
READ(io,*) Dummy, Dummy, pProp%Cp
READ(io,*) Dummy, Dummy, pProp%T_init
READ(io,*) Dummy, Dummy, pProp%T_ignite
READ(io,*) Dummy, Dummy, pProp%T_star0
READ(io,*) ! Convection(film) Coefficient
READ(io,*) Dummy, Dummy, pProp%film_const
READ(io,*) Dummy, Dummy, pProp%k
READ(io,*) Dummy, Dummy, pProp%Pr
READ(io,*) Dummy, Dummy, pProp%mu
READ(io,*) ! APN Model
READ(io,*) Dummy, Dummy, pProp%a_APN
READ(io,*) Dummy, Dummy, pProp%n_APN
READ(io,*) ! ZN Model
READ(io,*) ! L&R Model
READ(io,*) Dummy, Dummy, pProp%alpha
READ(io,*) Dummy, Dummy, pProp%beta
READ(io,*) Dummy, Dummy, pProp%k_ero
READ(io,*) Dummy, Dummy, pProp%StepInd
READ(io,*) Dummy, Dummy, pProp%StepMag
READ(io,*) ! Numerics
READ(io,*) Dummy, Dummy, pInput%IterMax
READ(io,*) ! Nozzle Ablation
READ(io,*) Dummy, Dummy, pNozzle%Ablation_Time
READ(io,*) Dummy, Dummy, pNozzle%Ablation_Coeff
CLOSE(io)

! Set lamc
pProp%lamc = pProp%alfac * pProp%rhoc * pProp%Cp

! Project Information
Write(*,*)
Write(*,*) 'SPBS >> ======================================================'
Write(*,*) 'SPBS >> SPBS Simulation Overview'
Write(*,*) 'SPBS >> ======================================================'
Write(*,*) 'SPBS >>'

! Surface Burning Information
SELECT CASE(pInput%indSB)
CASE(0)
  SELECT CASE(pInput%iFilm)
  CASE(1)
    WRITE(*,*) 'SPBS >> Film Coeff Model : Constant'
  CASE(2)
    WRITE(*,*) 'SPBS >> Film Coeff Model : Fully Developed Pipe Flow'
  CASE(3)
    WRITE(*,*) 'SPBS >> Film Coeff Model : Flat Plate Flow'
  CASE DEFAULT
    WRITE(*,*) 'SPBS >> Invalid Index For Film Coefficient'
    STOP
  END SELECT
CASE(1)
  WRITE(*,*) 'SPBS >> Whole Surface Burning'
  pInput%iBurn = 1 ! In whole surface burning case, APN model used
CASE DEFAULT
  WRITE(*,*) 'SPBS >> Invalid Index For Surface Burning'
  Stop
END SELECT

! Burning Model
SELECT CASE(pInput%iBurn)
CASE(1)
  WRITE(*,*) 'SPBS >> Burning  Model : APN Model'
CASE(2)
  WRITE(*,*) 'SPBS >> Burning  Model : ZN Model'
CASE DEFAULT
  WRITE(*,*) 'SPBS >> Invalid Burning Model'
  STOP
END SELECT

! Option for L&R Errosive Burning Model
SELECT CASE(pInput%indLR)
CASE(0)
  WRITE(*,*) 'SPBS >> Index LR Model : Not Using L&R Model'
CASE(1)
  WRITE(*,*) 'SPBS >> Index LR Model : Using L&R Model with x'
CASE(2)
  WRITE(*,*) 'SPBS >> Index LR Model : Using L&R Model with D'
CASE(3)
  WRITE(*,*) 'SPBS >> Index LR Model : Using L&R Model with f1(D)'
CASE(4)
  WRITE(*,*) 'SPBS >> Index LR Model : Using L&R Model with f2(D)'
CASE DEFAULT
  WRITE(*,*) 'SPBS >> Invalid Index for L&R Model'
  STOP
END SELECT

! Option for L&R Errosive Burning Weighting
SELECT CASE(pInput%indLR)
CASE(1:4)
  SELECT CASE(pInput%indk_ero)
  CASE(1)
    IF( nDim == 3 ) THEN
      WRITE(*,*) 'SPBS >> Index LR Model Weighting must be 2 in 3-D simulation'
      STOP
    END IF
    WRITE(*,*) 'SPBS >> Index LR Model Weighting : Weighting with Centerline'
  CASE(2)
    WRITE(*,*) 'SPBS >> Index LR Model Weighting : Weighting Itself'
  CASE DEFAULT
    WRITE(*,*) 'SPBS >> Invalid Index for L&R Model Weighting'
    STOP
  END SELECT
END SELECT

END SUBROUTINE

SUBROUTINE MAKELGRID(pLGrid)

! In or Out Variables
TYPE(t_LGrid) , INTENT(  out) :: pLGrid
! Local Variables
INTEGER :: Node
REAL(8) :: xMax
REAL(8) :: beta, bpm
REAL(8) :: bp1, bm1
REAL(8) :: cdzdx, cdzdx2
REAL(8) :: term1, term2

! Allocate LGrid Array
pLGrid%nNode = 101
ALLOCATE( pLGrid%x(pLGrid%nNode), pLGrid%z(pLGrid%nNode) )
ALLOCATE( pLGrid%dzdx(pLGrid%nNode), pLGrid%dzdx2(pLGrid%nNode) )

! Set Delta zeta
pLGrid%dz = 1d0 / DBLE(pLGrid%nNode - 1)

! Anderson, Tannehill and Pletcher
! Boundary Layer LGrid Control
xMax = -0.002d0 !(m)
beta = 1.01d0
bp1  = beta + 1d0
bm1  = beta - 1d0
bpm  = bp1 / bm1
cdzdx  = 2d0 * beta / DLOG(bpm) / xMax
cdzdx2 = -2d0 / xMax * cdzdx
DO Node = 1, pLGrid%nNode
  pLGrid%z(Node)     = DBLE(Node - 1) * pLGrid%dz
  Term1             = bpm**( 1d0 - pLGrid%z(Node) )
  pLGrid%x(Node)     = xMax * ( bp1 - bm1 * Term1 ) / ( Term1 + 1d0 )
  Term2             = beta**2 - ( 1d0 - pLGrid%x(Node) / xMax )**2
  pLGrid%dzdx(Node)  = cdzdx / Term2
  pLGrid%dzdx2(Node) = cdzdx2 * ( 1d0 - pLGrid%x(Node) / xMax ) / Term2**2
END DO

! Set dx1, dz2inv, dzsqinv
pLGrid%dx1     = pLGrid%dz / pLGrid%dzdx(1)
pLGrid%dz2inv  = 0.5d0 / pLGrid%dz
pLGrid%dzsqinv = 1d0 / pLGrid%dz**2

END SUBROUTINE

SUBROUTINE INITIALIZE(pInput,pSGrid,pLGrid,pProp,pBurn,BFace,P_fluid)

USE SPBS_BURN_MOD, ONLY : BURN_APN
USE SPBS_BURN_MOD, ONLY : BURN_LR
! In or Out Variables
TYPE(t_Input), INTENT(in   ) :: pInput
TYPE(t_SGrid), INTENT(in   ) :: pSGrid
TYPE(t_LGrid), INTENT(in   ) :: pLGrid
TYPE(t_Prop) , INTENT(in   ) :: pProp
TYPE(t_Burn) , INTENT(inout) :: pBurn
INTEGER      , INTENT(in   ) :: BFace
REAL(8)      , INTENT(in   ) :: P_fluid

! Maximum Allowable Time Step
pBurn%dTime = 0.5d0 * pLGrid%dz**2 / MAXVAL(DABS(pLGrid%dzdx(:)))**2 / pProp%alfac

! Initial Temperature
IF( pInput%indSB == 1 ) THEN
  ! Whole Surface Ignited with APN Model
  pBurn%Ignited = 1
  pBurn%T(:) = pProp%T_ignite
  pBurn%P_fluid = P_fluid
  CALL BURN_APN(pProp, pBurn)
  IF( pInput%indLR /= 0 ) CALL BURN_LR(pInput%indLR, pSGrid, pProp, pBurn, BFace)
ELSE
  pBurn%Ignited = 0
  pBurn%T(:) = pProp%T_init
END IF


! Steady State Temperature

! Restart Temperature Profile

END SUBROUTINE

END MODULE SPBS_INITIAL_MOD
