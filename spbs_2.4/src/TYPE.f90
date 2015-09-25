MODULE SPBS_TYPE_MOD

IMPLICIT NONE; PUBLIC

! Input Data Type
TYPE t_Input
  ! Basic Input Quantities
  INTEGER :: indSB ! whole surface burning
  INTEGER :: iFilm ! film coefficient method
  INTEGER :: iBurn ! burning model
  INTEGER :: indLR ! L&R model on/off
  INTEGER :: indk_ero ! L&R model weighting
  ! Iteration Info
  INTEGER :: IterMax
END TYPE t_Input

! Surface Grid Data Type
TYPE t_SGrid
  ! Basic Surface Grid Quantities
  INTEGER :: nDim
  INTEGER :: nBNode, nBFace
  REAL(8), ALLOCATABLE :: Bxyz(:,:)
  INTEGER, ALLOCATABLE :: Bf2n(:,:)
  ! Reference Location
  REAL(8) :: Ref_Location
  ! Axis Location Component
  INTEGER :: AxisCom
END TYPE t_SGrid

! Line Grid Data Type (for 1-D Heat Eq.)
TYPE t_LGrid
  ! Basic Line Grid Quantities
  INTEGER :: nNode
  REAL(8) :: dz, dx1
  REAL(8) :: dz2inv, dzsqinv
  REAL(8), ALLOCATABLE :: x(:), z(:)
  REAL(8), ALLOCATABLE :: dzdx(:), dzdx2(:)
END TYPE t_LGrid

! Propellant Data Type
TYPE t_Prop
  ! Total Number of Propellant Surface && Ignited Surface
  INTEGER :: nProp, nIgnited_APN, nIgnited_ZN
  ! Solid Propellant Quantities
  REAL(8) :: rhoc, MW, alfac, Cp, lamc
  REAL(8) :: T_init, T_ignite, T_star0
  ! Convection(film) Coefficient
  REAL(8) :: film_const ! film coefficient (Continuously Change)
  REAL(8) :: k, Pr, mu ! for fully developed pipe flow or flat plate flow
  ! APN Model Coeff
  REAL(8) :: a_APN, n_APN
  ! L&R Model Coeff
  REAL(8) :: alpha, beta, k_ero
  REAL(8) :: StepInd, StepMag
END TYPE t_Prop

! Burn Data Type
TYPE t_Burn
  ! Basic Burn Quantities
  INTEGER :: Ignited
  REAL(8) :: film
  REAL(8) :: P_fluid, T_fluid, G_fluid
  REAL(8) :: bRate, mFlux, fTemp
  REAL(8), ALLOCATABLE :: T(:), T0(:)
  ! Time Step
  REAL(8) :: dTime
END TYPE t_Burn

! Nozzle Data Type
TYPE t_Nozzle
  ! Nozzle Ablation Maximum Speed
  REAL(8) :: Ablation_Speed
  ! Ablation Time / Coefficient
  REAL(8) :: Ablation_Time
  REAL(8) :: Ablation_Coeff
  ! Nozzle Ablation Range / Maximum Speed Location
  REAL(8) :: Location_Range(2)
  REAL(8) :: Max_Loc
  ! Reference Pressure
  REAL(8) :: P_Ref
END TYPE

! Region Data Type
TYPE t_Region
  TYPE(t_Input) :: Input
  TYPE(t_SGrid) :: SGrid
  TYPE(t_LGrid) :: LGrid
  TYPE(t_Prop)  :: Prop
  TYPE(t_Burn), ALLOCATABLE :: Burn(:)
  TYPE(t_Nozzle) :: Nozzle
END TYPE t_Region

END MODULE SPBS_TYPE_MOD
