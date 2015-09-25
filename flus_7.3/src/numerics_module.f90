module numerics_module

use config_module
use grid_module
implicit none; private
public SetConv
public ComputeConvRoe
public ComputeConvRoeM2
public ComputeConvAUSMP
public ComputeConvAUSMPWP
public ComputeDiffAvgGrad
public ComputeTurbConvSca
public ComputeTurbDiffAvgGrad
public ComputeTurbSrcs
public ComputeAxiSrcs

! Global variables
integer :: nDim ! dimension
integer :: nDimP1 ! dimension+1
integer :: nDimP2 ! dimension+2
integer :: nDimP3 ! dimension+3
real(8) :: Gamma ! heat capacity ratio
real(8) :: Gamm1 ! gamma-1
real(8) :: GasConstant ! specific gas constant
real(8) :: Prandtl_Lam ! laminar prandtl number
real(8) :: Prandtl_Turb ! turbulent prandtl number
integer :: nCVar ! number of conservative variable
integer :: iLamMu ! laminar viscosity index in pv
integer :: iEdyMu ! eddy viscosity index in pv
integer :: iBlend ! blending function index in pv
real(8), parameter :: Eps = 1d-12 ! epsilon for weighting

contains

subroutine SetConv(Grid,Conf)
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  ! Set dimension
  nDim = Grid%nDim
  nDimP1 = Grid%nDim+1
  nDimP2 = Grid%nDim+2
  nDimP3 = Grid%nDim+3

  ! Set gamma / gamma-1 / gas constant
  ! Laminar / turbulent Prandtl number
  Gamma = Conf%Gamma
  Gamm1 = Conf%Gamma-1d0
  GasConstant = Conf%GasConstant
  Prandtl_Lam = Conf%Prandtl_Lam
  Prandtl_Turb = Conf%Prandtl_Turb

  ! Set number of variable
  nCVar = Grid%nDim+2

  ! Set index of the primitive variable
  select case(Conf%FlowModel)
  case(2) ! N-S equation
    iLamMu = Grid%nDim+2 ! laminar viscosity index in pv
    iEdyMu = Grid%nDim+3 ! eddy viscosity index in pv
  case(3) ! RANS equation
    iLamMu = Grid%nDim+4 ! laminar viscosity index in pv
    iEdyMu = Grid%nDim+5 ! eddy viscosity index in pv
    iBlend = Grid%nDim+6 ! blending function index in pv
  end select

end subroutine

subroutine ComputeConvRoe( pvl, pvr, UnitNormal, FaceVelocity, Flux )
  real(8), intent(in) :: pvl(0:)
  real(8), intent(in) :: pvr(0:)
  real(8), intent(in) :: UnitNormal(:)
  real(8), intent(in) :: FaceVelocity
  real(8), intent(out) :: Flux(0:)

  real(8) :: Temperature_l, Velocity_l(nDim), Pressure_l, Density_l, Energy_l, Enthalpy_l, ProjVelocity_l, RelVelocity_l
  real(8) :: Temperature_r, Velocity_r(nDim), Pressure_r, Density_r, Energy_r, Enthalpy_r, ProjVelocity_r, RelVelocity_r
  real(8) :: RoeRef1, RoeRef2, RoeDensity, RoeVelocity(nDim), RoeEnthalpy, RoeSoundSpeed, RoeProjVelocity, RoeRelVelocity
  real(8) :: Delta_Density, Delta_Velocity(nDim), Delta_Pressure, Delta_ProjVelocity
  real(8) :: Lambda(4), DeltaWave(4), EigenMatrix(4,0:nCVar-1), Flux_l(0:nCVar-1), Flux_r(0:nCVar-1)

  integer :: iCVar

  ! Set primitive variable
  Temperature_l = pvl(0)
  Temperature_r = pvr(0)
  Velocity_l(:) = pvl(1:nDim)
  Velocity_r(:) = pvr(1:nDim)
  Pressure_l = pvl(nDimP1)
  Pressure_r = pvr(nDimP1)
  Density_l = Pressure_l / ( GasConstant * Temperature_l )
  Density_r = Pressure_r / ( GasConstant * Temperature_r )
  Energy_l = GasConstant * Temperature_l / Gamm1 + 0.5d0 * sum( Velocity_l(:)**2 )
  Energy_r = GasConstant * Temperature_r / Gamm1 + 0.5d0 * sum( Velocity_r(:)**2 )
  Enthalpy_l = Energy_l + GasConstant * Temperature_l
  Enthalpy_r = Energy_r + GasConstant * Temperature_r
  ProjVelocity_l = sum( Velocity_l(:) * UnitNormal(:) )
  ProjVelocity_r = sum( Velocity_r(:) * UnitNormal(:) )
  RelVelocity_l = ProjVelocity_l - FaceVelocity
  RelVelocity_r = ProjVelocity_r - FaceVelocity

  ! Roe average variable
  RoeRef1 = dsqrt( Density_r / Density_l )
  RoeRef2 = 1d0 / ( RoeRef1 + 1d0 )
  RoeDensity = RoeRef1 * Density_l
  RoeVelocity(:) = ( RoeRef1 * Velocity_r(:) + Velocity_l(:) ) * RoeRef2
  RoeEnthalpy = ( RoeRef1 * Enthalpy_r + Enthalpy_l ) * RoeRef2
  RoeProjVelocity = ( RoeRef1 * ProjVelocity_r + ProjVelocity_l ) * RoeRef2
  RoeRelVelocity = RoeProjVelocity - FaceVelocity
  RoeSoundSpeed = dsqrt( Gamm1 * ( RoeEnthalpy - 0.5d0 * sum( RoeVelocity(:)**2 ) ) )

  ! Compute eigenvalue(wave speed)
  Lambda(1) = dabs( RoeRelVelocity - RoeSoundSpeed )
  Lambda(2) = dabs( RoeRelVelocity )
  Lambda(3) = dabs( RoeRelVelocity )
  Lambda(4) = dabs( RoeRelVelocity + RoeSoundSpeed )

  ! Compute delta primitive
  Delta_Density = Density_r - Density_l
  Delta_Velocity(:) = Velocity_r(:) - Velocity_l(:)
  Delta_Pressure = Pressure_r - Pressure_l
  Delta_ProjVelocity = ProjVelocity_r - ProjVelocity_l

  ! Compute delta wave strength
  DeltaWave(1) = ( Delta_Pressure - RoeDensity*RoeSoundSpeed*Delta_ProjVelocity )/( 2d0*RoeSoundSpeed**2 )
  DeltaWave(2) = Delta_Density - Delta_Pressure/RoeSoundSpeed**2
  DeltaWave(3) = RoeDensity
  DeltaWave(4) = ( Delta_Pressure + RoeDensity*RoeSoundSpeed*Delta_ProjVelocity )/( 2d0*RoeSoundSpeed**2 )

  ! Compute right eigenvector matrix
  EigenMatrix(1,0) = 1d0
  EigenMatrix(2,0) = 1d0
  EigenMatrix(3,0) = 0d0
  EigenMatrix(4,0) = 1d0
  EigenMatrix(1,1:nDim) = RoeVelocity(:) - RoeSoundSpeed*UnitNormal(:)
  EigenMatrix(2,1:nDim) = RoeVelocity(:)
  EigenMatrix(3,1:nDim) = Delta_Velocity(:) - Delta_ProjVelocity*UnitNormal(:)
  EigenMatrix(4,1:nDim) = RoeVelocity(:) + RoeSoundSpeed*UnitNormal(:)
  EigenMatrix(1,nDimP1) = RoeEnthalpy - RoeSoundSpeed*RoeProjVelocity
  EigenMatrix(2,nDimP1) = 0.5d0 * sum( RoeVelocity(:)**2 )
  EigenMatrix(3,nDimP1) = sum( RoeVelocity(1:nDim) * Delta_Velocity(1:nDim) ) - RoeProjVelocity*Delta_ProjVelocity
  EigenMatrix(4,nDimP1) = RoeEnthalpy + RoeSoundSpeed*RoeProjVelocity

  ! Compute left flux
  Flux_l(0) = Density_l * RelVelocity_l
  Flux_l(1:nDim) = Density_l * RelVelocity_l * Velocity_l(:) + Pressure_l * UnitNormal(:)
  Flux_l(nDimP1) = Density_l * RelVelocity_l * Enthalpy_l + Pressure_l * FaceVelocity

  ! Compute right flux
  Flux_r(0) = Density_r * RelVelocity_r
  Flux_r(1:nDim) = Density_r * RelVelocity_r * Velocity_r(:) + Pressure_r * UnitNormal(:)
  Flux_r(nDimP1) = Density_r * RelVelocity_r * Enthalpy_r + Pressure_r * FaceVelocity

  ! Compute convective flux at interface
  do iCVar = 0, nCVar-1
    Flux(iCVar) = 0.5d0*( Flux_l(iCVar) + Flux_r(iCVar) - sum( Lambda(:)*DeltaWave(:)*EigenMatrix(:,iCVar) ) )
  end do

end subroutine

subroutine ComputeConvRoeM2( pvl, pvr, UnitNormal, FaceVelocity, PRateMin, Flux )
  real(8), intent(in) :: pvl(0:)
  real(8), intent(in) :: pvr(0:)
  real(8), intent(in) :: UnitNormal(:)
  real(8), intent(in) :: FaceVelocity
  real(8), intent(in) :: PRateMin
  real(8), intent(out) :: Flux(0:)

  real(8) :: Temperature_l, Velocity_l(nDim), Pressure_l, Density_l, Energy_l, Enthalpy_l, ProjVelocity_l, RelVelocity_l
  real(8) :: Temperature_r, Velocity_r(nDim), Pressure_r, Density_r, Energy_r, Enthalpy_r, ProjVelocity_r, RelVelocity_r
  real(8) :: RoeRef1, RoeRef2, RoeDensity, RoeVelocity(nDim), RoeEnthalpy, RoeSoundSpeed, RoeProjVelocity, RoeRelVelocity, RoeRelMach
  real(8) :: Delta_Density, Delta_Velocity(nDim), Delta_Pressure, Delta_ProjVelocity, Delta_Enthalpy
  real(8) :: DelQ(0:nCVar-1), BDelQ(0:nCVar-1), Flux_l(0:nCVar-1), Flux_r(0:nCVar-1)
  real(8) :: ff, hh, gg, Coeff, b1, b2

  ! Set primitive variable
  Temperature_l = pvl(0)
  Temperature_r = pvr(0)
  Velocity_l(:) = pvl(1:nDim)
  Velocity_r(:) = pvr(1:nDim)
  Pressure_l = pvl(nDimP1)
  Pressure_r = pvr(nDimP1)
  Density_l = Pressure_l / ( GasConstant * Temperature_l )
  Density_r = Pressure_r / ( GasConstant * Temperature_r )
  Energy_l = GasConstant * Temperature_l / Gamm1 + 0.5d0 * sum( Velocity_l(:)**2 )
  Energy_r = GasConstant * Temperature_r / Gamm1 + 0.5d0 * sum( Velocity_r(:)**2 )
  Enthalpy_l = Energy_l + GasConstant * Temperature_l
  Enthalpy_r = Energy_r + GasConstant * Temperature_r
  ProjVelocity_l = sum( Velocity_l(:) * UnitNormal(:) )
  ProjVelocity_r = sum( Velocity_r(:) * UnitNormal(:) )
  RelVelocity_l = ProjVelocity_l - FaceVelocity
  RelVelocity_r = ProjVelocity_r - FaceVelocity

  ! Roe average variable
  RoeRef1 = dsqrt( Density_r / Density_l )
  RoeRef2 = 1d0 / ( RoeRef1 + 1d0 )
  RoeDensity = RoeRef1 * Density_l
  RoeVelocity(:) = ( RoeRef1 * Velocity_r(:) + Velocity_l(:) ) * RoeRef2
  RoeEnthalpy = ( RoeRef1 * Enthalpy_r + Enthalpy_l ) * RoeRef2
  RoeProjVelocity = ( RoeRef1 * ProjVelocity_r + ProjVelocity_l ) * RoeRef2
  RoeRelVelocity = RoeProjVelocity - FaceVelocity
  RoeSoundSpeed = dsqrt( Gamm1 * ( RoeEnthalpy - 0.5d0 * sum( RoeVelocity(:)**2 ) ) )
  RoeRelMach = RoeRelVelocity / RoeSoundSpeed

  ! RoeM function1 f
  if( sum( RoeVelocity(:)**2 ) < Eps ) then
    ff = 1d0
  else
    hh = 1d0 - PRateMin
    ff = dabs( RoeRelMach )**hh
  end if

  ! RoeM function2 g
  if( dabs( RoeRelMach ) == 0d0 ) then
    gg = 1d0
  else
    gg = dabs( RoeRelMach )**( 1d0 - dmin1( Pressure_r/Pressure_l, Pressure_l/Pressure_r ) )
  end if

  ! Compute delta primitive
  Delta_Density = Density_r - Density_l
  Delta_Velocity(:) = Velocity_r(:) - Velocity_l(:)
  Delta_Pressure = Pressure_r - Pressure_l
  Delta_ProjVelocity = ProjVelocity_r - ProjVelocity_l
  Delta_Enthalpy = Enthalpy_r - Enthalpy_l

  ! Compute delta conservative
  DelQ(0) = Delta_Density
  DelQ(1:nDim) = Density_r * Velocity_r(:) - Density_l * Velocity_l(:)
  DelQ(nDimP1) = Density_r * Enthalpy_r - Density_l * Enthalpy_l

  ! Compute B*DelQ
  Coeff = Delta_Density - ff * Delta_Pressure / RoeSoundSpeed**2
  BDelQ(0) = Coeff
  BDelQ(1:nDim) = Coeff * RoeVelocity(:) + RoeDensity * ( Delta_Velocity(:) - UnitNormal(:) * Delta_ProjVelocity )
  BDelQ(nDimP1) = Coeff * RoeEnthalpy + RoeDensity * Delta_Enthalpy

  ! Calculate b1, b2
  b1 = dmax1( 0d0, RoeRelVelocity + RoeSoundSpeed, RelVelocity_r + RoeSoundSpeed )
  b2 = dmin1( 0d0, RoeRelVelocity - RoeSoundSpeed, RelVelocity_l - RoeSoundSpeed )

  ! Compute left flux
  Flux_l(0) = Density_l * RelVelocity_l
  Flux_l(1:nDim) = Density_l * RelVelocity_l * Velocity_l(:) + Pressure_l * UnitNormal(:)
  Flux_l(nDimP1) = Density_l * RelVelocity_l * Enthalpy_l + Pressure_l * FaceVelocity

  ! Compute right flux
  Flux_r(0) = Density_r * RelVelocity_r
  Flux_r(1:nDim) = Density_r * RelVelocity_r * Velocity_r(:) + Pressure_r * UnitNormal(:)
  Flux_r(nDimP1) = Density_r * RelVelocity_r * Enthalpy_r + Pressure_r * FaceVelocity

  Flux(:) = ( b1 * Flux_l(:) - b2 * Flux_r(:) ) / ( b1 - b2 ) + ( b1 * b2 ) / ( b1 - b2 ) * ( DelQ(:) - gg / ( 1d0 + dabs(RoeRelMach) ) * BDelQ(:) )

end subroutine

subroutine ComputeConvAUSMP( pvl, pvr, UnitNormal, FaceVelocity, Flux )
  real(8), intent(in) :: pvl(0:)
  real(8), intent(in) :: pvr(0:)
  real(8), intent(in) :: UnitNormal(:)
  real(8), intent(in) :: FaceVelocity
  real(8), intent(out) :: Flux(0:)

  real(8) :: Temperature_l, Velocity_l(nDim), Pressure_l, Density_l, Energy_l, Enthalpy_l, ProjVelocity_l, RelVelocity_l
  real(8) :: Temperature_r, Velocity_r(nDim), Pressure_r, Density_r, Energy_r, Enthalpy_r, ProjVelocity_r, RelVelocity_r
  real(8) :: Critical_SoundSpeed_l, Modified_SoundSpeed_l
  real(8) :: Critical_SoundSpeed_r, Modified_SoundSpeed_r
  real(8) :: Mach_l, Mach_l_p, Mach_Avg_l, Pressure_l_p
  real(8) :: Mach_r, Mach_r_m, Mach_Avg_r, Pressure_r_m
  real(8) :: Mach_Avg, Pressure_Avg, SoundSpeed
  real(8), parameter :: alpha = 3d0 / 16d0
  real(8), parameter :: beta  = 1d0 / 8d0

  ! Set primitive variable
  Temperature_l = pvl(0)
  Temperature_r = pvr(0)
  Velocity_l(:) = pvl(1:nDim)
  Velocity_r(:) = pvr(1:nDim)
  Pressure_l = pvl(nDimP1)
  Pressure_r = pvr(nDimP1)
  Density_l = Pressure_l / ( GasConstant * Temperature_l )
  Density_r = Pressure_r / ( GasConstant * Temperature_r )
  Energy_l = GasConstant * Temperature_l / Gamm1 + 0.5d0 * sum( Velocity_l(:)**2 )
  Energy_r = GasConstant * Temperature_r / Gamm1 + 0.5d0 * sum( Velocity_r(:)**2 )
  Enthalpy_l = Energy_l + GasConstant * Temperature_l
  Enthalpy_r = Energy_r + GasConstant * Temperature_r
  ProjVelocity_l = sum( Velocity_l(:) * UnitNormal(:) )
  ProjVelocity_r = sum( Velocity_r(:) * UnitNormal(:) )
  RelVelocity_l = ProjVelocity_l - FaceVelocity
  RelVelocity_r = ProjVelocity_r - FaceVelocity

  ! Critical speed of sound is not changed in ALE method
  Critical_SoundSpeed_l = dsqrt( 2d0 * Gamm1 * Enthalpy_l / ( Gamma + 1d0 ) )
  Critical_SoundSpeed_r = dsqrt( 2d0 * Gamm1 * Enthalpy_r / ( Gamma + 1d0 ) )
  Modified_SoundSpeed_l = Critical_SoundSpeed_l**2 / dmax1( Critical_SoundSpeed_l, dabs(ProjVelocity_l) )
  Modified_SoundSpeed_r = Critical_SoundSpeed_r**2 / dmax1( Critical_SoundSpeed_r, dabs(ProjVelocity_r) )
  SoundSpeed = dmin1( Modified_SoundSpeed_l, Modified_SoundSpeed_r )

  ! Left / right Mach
  Mach_l = RelVelocity_l/SoundSpeed
  Mach_r = RelVelocity_r/SoundSpeed

  ! Mach number i plus
  if( dabs(Mach_l) >= 1d0 ) then
    Mach_l_p = 0.5d0 * ( Mach_l + dabs( Mach_l ) )
    Pressure_l_p = Mach_l_p / Mach_l
  else
    Mach_l_p = 0.25d0 * ( Mach_l + 1d0 )**2 + beta * ( Mach_l**2 - 1d0 )**2
    Pressure_l_p = 0.25d0 * ( Mach_l + 1d0 )**2 * ( 2d0 - Mach_l ) + alpha * Mach_l * ( Mach_l**2 - 1d0 )**2
  end if

  ! Mach number j minus
  if( dabs(Mach_r) >= 1d0 ) then
    Mach_r_m = 0.5d0*( Mach_r - dabs( Mach_r ) )
    Pressure_r_m = Mach_r_m / Mach_r
  else
    Mach_r_m = -0.25d0 * ( Mach_r - 1d0 )**2 - beta * ( Mach_r**2 - 1d0 )**2
    Pressure_r_m = 0.25d0*( Mach_r - 1d0 )**2 * ( 2d0 + Mach_r ) - alpha * Mach_r * ( Mach_r**2 - 1d0 )**2
  end if

  ! Average mach number i j / pressure
  Mach_Avg = Mach_l_p + Mach_r_m
  Mach_Avg_l = 0.5d0 * ( Mach_Avg + dabs( Mach_Avg ) )
  Mach_Avg_r = 0.5d0 * ( Mach_Avg - dabs( Mach_Avg ) )
  Pressure_Avg = Pressure_l_p * Pressure_l + Pressure_r_m * Pressure_r

  ! Compute convective flux at interface
  Flux(0) = SoundSpeed * ( Mach_Avg_l * Density_l + Mach_Avg_r * Density_r )
  Flux(1:nDim) = SoundSpeed * ( Mach_Avg_l * Density_l * Velocity_l(:) + Mach_Avg_r * Density_r * Velocity_r(:) ) + Pressure_Avg * UnitNormal(:)
  Flux(nDimP1) = SoundSpeed * ( Mach_Avg_l * Density_l * Enthalpy_l + Mach_Avg_r * Density_r * Enthalpy_r ) + Pressure_Avg * FaceVelocity

end subroutine

subroutine ComputeConvAUSMPWP( pvl, pvr, UnitNormal, FaceVelocity, PMinRate, Flux )
  real(8), intent(in) :: pvl(0:)
  real(8), intent(in) :: pvr(0:)
  real(8), intent(in) :: UnitNormal(:)
  real(8), intent(in) :: FaceVelocity
  real(8), intent(in) :: PMinRate
  real(8), intent(out) :: Flux(0:)

  real(8) :: Temperature_l, Velocity_l(nDim), Pressure_l, Density_l, Energy_l, Enthalpy_l, ProjVelocity_l, RelVelocity_l
  real(8) :: Temperature_r, Velocity_r(nDim), Pressure_r, Density_r, Energy_r, Enthalpy_r, ProjVelocity_r, RelVelocity_r
  real(8) :: Mach_l, Mach_l_p, Mach_l_p_Avg, Pressure_l_p, TanVelocity2_l
  real(8) :: Mach_r, Mach_r_m, Mach_r_m_Avg, Pressure_r_m, TanVelocity2_r
  real(8) :: Pressure_Avg, SoundSpeed, NormEnthalpy, fl, fr, ww
  real(8), parameter :: alpha = 3d0 / 16d0

  ! Set primitive variable
  Temperature_l = pvl(0)
  Temperature_r = pvr(0)
  Velocity_l(:) = pvl(1:nDim)
  Velocity_r(:) = pvr(1:nDim)
  Pressure_l = pvl(nDimP1)
  Pressure_r = pvr(nDimP1)
  Density_l = Pressure_l / ( GasConstant * Temperature_l )
  Density_r = Pressure_r / ( GasConstant * Temperature_r )
  Energy_l = GasConstant * Temperature_l / Gamm1 + 0.5d0 * sum( Velocity_l(:)**2 )
  Energy_r = GasConstant * Temperature_r / Gamm1 + 0.5d0 * sum( Velocity_r(:)**2 )
  Enthalpy_l = Energy_l + GasConstant * Temperature_l
  Enthalpy_r = Energy_r + GasConstant * Temperature_r
  ProjVelocity_l = sum( Velocity_l(:) * UnitNormal(:) )
  ProjVelocity_r = sum( Velocity_r(:) * UnitNormal(:) )
  RelVelocity_l = ProjVelocity_l - FaceVelocity
  RelVelocity_r = ProjVelocity_r - FaceVelocity

  ! Compute normal Enthalpy
  TanVelocity2_l = sum( Velocity_l(:)**2 ) - ProjVelocity_l**2
  TanVelocity2_r = sum( Velocity_r(:)**2 ) - ProjVelocity_r**2
  NormEnthalpy = 0.5d0 * ( Enthalpy_l + Enthalpy_r - 0.5d0 * ( TanVelocity2_l + TanVelocity2_r ) )

  ! Compute speed of sound
  SoundSpeed = dsqrt( 2d0 * Gamm1 / ( Gamma + 1d0 ) * NormEnthalpy )

  ! Modify speed of sound
  if( 0.5d0 * ( ProjVelocity_l + ProjVelocity_r ) >= 0d0 ) then
    SoundSpeed = SoundSpeed**2 / dmax1( dabs( ProjVelocity_l ), SoundSpeed )
  else
    SoundSpeed = SoundSpeed**2 / dmax1( dabs( ProjVelocity_r ), SoundSpeed )
  end if

  ! Left / right Mach
  Mach_l = RelVelocity_l / SoundSpeed
  Mach_r = RelVelocity_r / SoundSpeed

  ! Mach number i plus
  if( dabs(Mach_l) >= 1d0 ) then
    Mach_l_p = 0.5d0 * ( Mach_l + dabs( Mach_l ) )
    Pressure_l_p = Mach_l_p / Mach_l
  else
    Mach_l_p = 0.25d0 * ( Mach_l + 1d0 )**2
    Pressure_l_p = 0.25d0 * ( Mach_l + 1d0 )**2 * ( 2d0 - Mach_l ) + alpha * Mach_l * ( Mach_l**2 - 1d0 )**2
  end if

  ! Mach number j minus
  if( dabs(Mach_r) >= 1d0 ) then
    Mach_r_m = 0.5d0*( Mach_r - dabs( Mach_r ) )
    Pressure_r_m = Mach_r_m / Mach_r
  else
    Mach_r_m = -0.25d0 * ( Mach_r - 1d0 )**2
    Pressure_r_m = 0.25d0*( Mach_r - 1d0 )**2 * ( 2d0 + Mach_r ) - alpha * Mach_r * ( Mach_r**2 - 1d0 )**2
  end if

  ! Average pressure
  Pressure_Avg = Pressure_l_p * Pressure_l + Pressure_r_m * Pressure_r

  ! f calculation
  if( Pressure_Avg < Eps ) then
    fl = 0d0
    fr = 0d0
  else
    fl = ( Pressure_l / Pressure_Avg - 1d0 )*dmin1( 1d0, PMinRate )**2d0
    fr = ( Pressure_r / Pressure_Avg - 1d0 )*dmin1( 1d0, PMinRate )**2d0
  end if

  ! Weighting function
  ww = 1d0 - dmin1( Pressure_l/Pressure_r, Pressure_r/Pressure_l )**3

  IF ( Mach_l_p + Mach_r_m >= 0d0 ) then
    Mach_l_p_Avg = Mach_l_p + Mach_r_m * ( ( 1d0 - ww ) * ( 1d0 + fr ) - fl )
    Mach_r_m_Avg = Mach_r_m * ww * ( 1d0 + fr )
  else
    Mach_l_p_Avg = Mach_l_p * ww * ( 1d0 + fl )
    Mach_r_m_Avg = Mach_r_m + Mach_l_p * ( ( 1d0 - ww ) * ( 1d0 + fl ) - fr )
  end if

  Flux(0) = SoundSpeed * ( Mach_l_p_Avg * Density_l + Mach_r_m_Avg * Density_r )
  Flux(1:nDim) = SoundSpeed * ( Mach_l_p_Avg * Density_l * Velocity_l(:) + Mach_r_m_Avg * Density_r * Velocity_r(:) ) + Pressure_Avg * UnitNormal(:)
  Flux(nDimP1) = SoundSpeed * ( Mach_l_p_Avg * Density_l * Enthalpy_l + Mach_r_m_Avg * Density_r * Enthalpy_r ) + Pressure_Avg * FaceVelocity

end subroutine

subroutine ComputeDiffAvgGrad( pvl, pvr, gpvl, gpvr, limiterl, limiterr, cenl, cenr, UnitNormal, Flux )
  real(8), intent(in) :: pvl(0:)
  real(8), intent(in) :: pvr(0:)
  real(8), intent(in) :: gpvl(:,0:)
  real(8), intent(in) :: gpvr(:,0:)
  real(8), intent(in) :: limiterl(0:)
  real(8), intent(in) :: limiterr(0:)
  real(8), intent(in) :: cenl(:)
  real(8), intent(in) :: cenr(:)
  real(8), intent(in) :: UnitNormal(:)
  real(8), intent(out) :: Flux(0:)

  integer :: iDim, jDim, iVar
  real(8) :: Mean_Laminar_Viscosity
  real(8) :: Mean_Eddy_Viscosity
  real(8) :: Mean_Conductivity
  real(8) :: Total_Viscosity
  real(8) :: Vector(nDim), Dist2
  real(8) :: Mean_Grad(nDim,0:nDim)
  real(8) :: Proj_Mean_Grad(0:nDim)
  real(8) :: Mean_Vel(nDim), Div_Vel
  real(8) :: Flux_Tensor(nDim,0:nCVar-1)
  real(8), parameter :: Two3 = 2d0/3d0

  ! Compute vector going from left to right
  Vector(:) = cenr(:) - cenl(:)
  Dist2 = sum( Vector(:)**2 )

  ! Compute mean viscosities and conductivity
  Mean_Laminar_Viscosity = 0.5d0 * ( pvl(iLamMu) + pvr(iLamMu) )
  Mean_Eddy_Viscosity = 0.5d0 * ( pvl(iEdyMu) + pvr(iEdyMu) )
  Mean_Conductivity = Gamma / Gamm1 * GasConstant * ( Mean_Laminar_Viscosity / Prandtl_Lam + Mean_Eddy_Viscosity / Prandtl_Turb )

  ! Compute total viscosity
  Total_Viscosity = Mean_Laminar_Viscosity + Mean_Eddy_Viscosity

  ! Compute modified mean gradient velocity and temperature
  do iVar = 0, nDim
    Mean_Grad(:,iVar) = 0.5d0 * ( limiterl(iVar) * gpvl(:,iVar) + limiterr(iVar) * gpvr(:,iVar) )
    Proj_Mean_Grad(iVar) = sum( Mean_Grad(:,iVar) * Vector(:) )
    Mean_Grad(:,iVar) = Mean_Grad(:,iVar) - ( Proj_Mean_Grad(iVar) - ( pvr(iVar) - pvl(iVar) ) ) * Vector(:) / Dist2
  end do

  ! Compute mean velocity
  Mean_Vel(:) = 0.5d0 * ( pvl(1:nDim) + pvr(1:nDim) )

  ! Compute velocity divergence
  Div_vel = 0d0
  do iDim  = 1, nDim
    Div_Vel = Div_Vel + Mean_Grad(iDim,iDim)
  end do

  ! Compute flux tensor
  Flux_Tensor(:,0) = 0d0
  do jDim = 1, nDim
    do iDim = 1, nDim
      Flux_Tensor(iDim,jDim) = Total_Viscosity * ( Mean_Grad(jDim,iDim) + Mean_Grad(iDim,jDim) )
      if( iDim == jDim ) Flux_Tensor(jDim,iDim) = Flux_Tensor(jDim,iDim) - Two3 * Total_Viscosity * Div_Vel
    end do
  end do
  do iDim = 1, nDim
    Flux_Tensor(iDim,nDimP1) = sum( Flux_Tensor(:,iDim) * Mean_Vel(:) ) + Mean_Conductivity * Mean_Grad(iDim,0)
  end do

  ! Update Flux
  do iVar = 0, nCVar-1
    Flux(iVar) = -sum( Flux_Tensor(:,iVar) * UnitNormal(:) )
  end do

end subroutine

subroutine ComputeTurbConvSca(pvl,pvr,UnitNormal,FaceVelocity,Flux)
  real(8), intent(in) :: pvl(0:)
  real(8), intent(in) :: pvr(0:)
  real(8), intent(in) :: UnitNormal(:)
  real(8), intent(in) :: FaceVelocity
  real(8), intent(out) :: Flux(0:)

  real(8) :: Density_l, Kine_l, Omega_l, Velocity_l(nDim), ProjVelocity_l, RelVelocity_l
  real(8) :: Density_r, Kine_r, Omega_r, Velocity_r(nDim), ProjVelocity_r, RelVelocity_r
  real(8) :: RelVelocity_Avg, RelVelocity_Avg_l, RelVelocity_Avg_r

  ! Set primitive variable
  Density_l = pvl(nDimP1) / ( GasConstant * pvl(0) )
  Density_r = pvr(nDimP1) / ( GasConstant * pvr(0) )
  Kine_l = pvl(nDimP2)
  Kine_r = pvr(nDimP2)
  Omega_l = pvl(nDimP3)
  Omega_r = pvr(nDimP3)
  Velocity_l(:) = pvl(1:nDim)
  Velocity_r(:) = pvr(1:nDim)
  ProjVelocity_l = sum( Velocity_l(:) * UnitNormal(:) )
  ProjVelocity_r = sum( Velocity_r(:) * UnitNormal(:) )
  RelVelocity_l = ProjVelocity_l - FaceVelocity
  RelVelocity_r = ProjVelocity_r - FaceVelocity

  ! Set average relative velocity
  RelVelocity_Avg = 0.5d0 * ( RelVelocity_l + RelVelocity_r )
  RelVelocity_Avg_l = 0.5d0*(RelVelocity_Avg+dabs(RelVelocity_Avg))
  RelVelocity_Avg_r = 0.5d0*(RelVelocity_Avg-dabs(RelVelocity_Avg))

  ! Compute flux
  Flux(nDimP2) = RelVelocity_Avg_l*Density_l*Kine_l+RelVelocity_Avg_r*Density_r*Kine_r
  Flux(nDimP3) = RelVelocity_Avg_l*Density_l*Omega_l+RelVelocity_Avg_r*Density_r*Omega_r

end subroutine

subroutine ComputeTurbDiffAvgGrad(pvl,pvr,gpvl,gpvr,limiterl,limiterr,cenl,cenr,UnitNormal,Flux)
  real(8), intent(in) :: pvl(0:)
  real(8), intent(in) :: pvr(0:)
  real(8), intent(in) :: gpvl(:,0:)
  real(8), intent(in) :: gpvr(:,0:)
  real(8), intent(in) :: limiterl(0:)
  real(8), intent(in) :: limiterr(0:)
  real(8), intent(in) :: cenl(:)
  real(8), intent(in) :: cenr(:)
  real(8), intent(in) :: UnitNormal(:)
  real(8), intent(out) :: Flux(0:)

  integer :: iVar
  real(8) :: Vector(nDim), Dist2
  real(8) :: Sigma_Kine_l, Sigma_Omega_l
  real(8) :: Sigma_Kine_r, Sigma_Omega_r
  real(8) :: Diff_l_Kine, Diff_l_Omega
  real(8) :: Diff_r_Kine, Diff_r_Omega
  real(8) :: Diff_Kine, Diff_Omega
  real(8) :: Mean_Grad(nDim,nDimP2:nDimP3)
  real(8) :: Proj_Mean_Grad(nDimP2:nDimP3)

  ! Compute vector going from left to right
  Vector(:) = cenr(:) - cenl(:)
  Dist2 = sum( Vector(:)**2 )

  ! Compute blended constant
  Sigma_Kine_l  = pvl(iBlend)*0.85d0 + (1d0-pvl(iBlend))*1d0
  Sigma_Kine_r  = pvr(iBlend)*0.85d0 + (1d0-pvr(iBlend))*1d0
  Sigma_Omega_l = pvl(iBlend)*0.5d0  + (1d0-pvl(iBlend))*0.856d0
  Sigma_Omega_r = pvr(iBlend)*0.5d0  + (1d0-pvr(iBlend))*0.856d0

  ! Compute mean effective viscosity
  Diff_l_Kine  = pvl(iLamMu) + sigma_Kine_l  * pvl(iEdyMu)
  Diff_r_Kine  = pvr(iLamMu) + sigma_Kine_r  * pvr(iEdyMu)
  Diff_l_Omega = pvl(iLamMu) + sigma_Omega_l * pvl(iEdyMu)
  Diff_r_Omega = pvr(iLamMu) + sigma_Omega_r * pvr(iEdyMu)

  Diff_Kine  = 0.5d0*( Diff_l_Kine + Diff_r_Kine )
  Diff_Omega = 0.5d0*( Diff_l_Omega + Diff_r_Omega )

  ! Compute modified mean gradient
  do iVar = nDimP2, nDimP3
    Mean_Grad(:,iVar) = 0.5d0*( limiterl(iVar) * gpvl(:,iVar) + limiterr(iVar) * gpvr(:,iVar) )
    Proj_Mean_Grad(iVar) = sum( Mean_Grad(:,iVar) * Vector(:) )
    Mean_Grad(:,iVar) = Mean_Grad(:,iVar) - ( Proj_Mean_Grad(iVar) - ( pvr(iVar) - pvl(iVar) ) ) * Vector(:) / Dist2
  end do

  ! Compute flux
  Flux(nDimP2) = - Diff_Kine  * sum( Mean_Grad(:,nDimP2) * UnitNormal(:) )
  Flux(nDimP3) = - Diff_Omega * sum( Mean_Grad(:,nDimP3) * UnitNormal(:) )

end subroutine

subroutine ComputeTurbSrcs(pv,gpv,Srcs)
  real(8), intent(in) :: pv(0:)
  real(8), intent(in) :: gpv(:,0:)
  real(8), intent(out) :: Srcs(0:)

  integer :: iDim
  real(8) :: Density, Kine, Omega, EddyMu, StrainMag2
  real(8) :: alfa_Blend, beta_Blend, F1
  real(8) :: Div, Pk, Pw, Dk, Dw, CDkw
  real(8), parameter :: inv3 = 1d0/3d0

  ! Set primitive
  Density = pv(nDimP1) / ( GasConstant * pv(0) )
  Kine = pv(nDimP2)
  Omega = pv(nDimP3)
  EddyMu = pv(iEdyMu)
  F1 = pv(iBlend)

  ! Compute blended constants
  alfa_Blend = F1*5d0/9d0  + (1d0-F1)*0.44d0
  beta_Blend = F1*3d0/40d0 + (1d0-F1)*0.0828d0

  ! Compute divergence
  Div = 0d0
  do iDim = 1, nDim
    Div = Div + gpv(iDim,iDim)
  end do

  StrainMag2 = 0d0
  ! Compute shear strain magnitude (add diagonal part)
  do iDim = 1, nDim
    StrainMag2 = StrainMag2 + 2d0*( gpv(iDim,iDim) - inv3*Div )**2
  end do

  ! Compute shear strain magnitude (add off diagonals)
  StrainMag2 = StrainMag2 + ( gpv(1,2) + gpv(2,1) )**2
  if( nDim == 3 ) then
    StrainMag2 = StrainMag2 + ( gpv(1,3) + gpv(3,1) )**2
    StrainMag2 = StrainMag2 + ( gpv(2,3) + gpv(3,2) )**2
  end if

  ! Compute production term
  Pk = EddyMu * StrainMag2 - 2d0/3d0*Density*Kine*Div
  Pk = dmin1(pk,10d0*0.09d0*Density*Kine*Omega)
  Pk = dmax1(pk,0d0)
  Pw = alfa_Blend*Density*Pk/EddyMu

  ! Compute dissipation term
  Dk = 0.09d0*Density*Omega*Kine
  Dw = beta_Blend*Density*Omega**2

  ! Compute cross diffusion
  CDkw = sum( gpv(1:nDim,nDimP2)*gpv(1:nDim,nDimP3) )
  CDkw = CDkw*2d0*Density*0.856d0/Omega

  ! Compute residual
  Srcs(nDimP2) = - ( Pk - Dk )
  Srcs(nDimP3) = - ( Pw - Dw + (1d0-F1)*CDkw )

end subroutine

subroutine ComputeAxiSrcs(pv,cen,Srcs)
  real(8), intent(in) :: pv(0:)
  real(8), intent(in) :: cen(:)
  real(8), intent(out) :: Srcs(0:)

  real(8) :: Density

  Density = pv(nDimP1) / ( GasConstant * pv(0) )
  Srcs(0) = Density * pv(2)
  Srcs(1) = Density * pv(2) * pv(1)
  Srcs(2) = Density * pv(2) * pv(2)
  Srcs(3) = pv(2) * ( Gamma / Gamm1 * pv(nDimP1) + 0.5d0 * Density * sum( pv(1:2)**2 ) )
  Srcs(:) = Srcs(:) / dabs( cen(2) )

end subroutine

end module
