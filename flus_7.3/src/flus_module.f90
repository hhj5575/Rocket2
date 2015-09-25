module flus_module
!---------------------------------------------------------------------------------!
! Description:                                                                    !
!   FLuid Unstructured Simulation                                                 !
!---------------------------------------------------------------------------------!
! Time Scheme:                Flux Scheme:              Flux LimIter:             !
!   1 : Grid Moving Test        1 : Roe                   1 : 1st Order           !
!   2 : Euler Explicit          2 : ReoM2                 2 : Barth               !
!   3 : TVD-RK 3rd              3 : AUSM+                 3 : Venkatakrishnan     !
!   4 : ESDIRK64                4 : AUSMPW+               4 : MLP-u1              !
!                                                         5 : MLP-u2              !
!---------------------------------------------------------------------------------!
! Mesh Smooth:                Remesh:                   Data Transfer:            !
!   IDW Smoothing               2D : TRIGEN Library       Consistent              !
!                               3D : TETGEN Library       Interpolation Method    !
!---------------------------------------------------------------------------------!
! Post Format:                Boundary Condition:                                 !
!   FLUS Format                 1 : Inviscid Wall         2 : Viscous Wall        !
!   (In-house)                  3 : Pressure Inlet        4 : Injection Inlet     !
!                               5 : Pressure Outlet       6 : Pressure Farfield   !
!                               7 : Symmetric Wall        8 : Axisymmetric Wall   !
!                               9 : Membrane Nozzle      10 : Propellant Wall     !
!                              11 : Ignitor Model        12 : Pressure Sensor     !
!                              13 : Nozzle Throat        14 : APN model           !
!---------------------------------------------------------------------------------!
!                                                                Written by: LCS  !
!---------------------------------------------------------------------------------!
use config_module
use grid_module
use mixture_module
use process_module
use interface_module
implicit none; private
public flus

type(t_Conf) :: Conf
type(t_Grid) :: Grid
type(t_Mixt) :: Mixt

contains

subroutine flus(Operation,Iteration,DeltaTime,nBNode,nBFace,Bxyz,Bf2n,BPatch,BCFlag,InVars,OutVars,Ref_Pressure)
  integer, intent(in) :: Operation ! operation flag
  integer, intent(in) :: Iteration ! system Iteration
  real(8), intent(in) :: DeltaTime ! system delta time
  integer, intent(inout) :: nBNode ! number of boundary nodes
  integer, intent(inout) :: nBFace ! number of boundary faces
  real(8), intent(inout), allocatable :: Bxyz(:,:) ! boundary node coordinates
  integer, intent(inout), allocatable :: Bf2n(:,:) ! boundary face node connectivity
  integer, intent(inout), allocatable :: BPatch(:) ! boundary face patch info
  integer, intent(inout), allocatable :: BCFlag(:) ! boundary condition flag
  real(8), intent(inout), allocatable :: InVars(:,:) ! inward variables
  real(8), intent(inout), allocatable :: OutVars(:,:) ! outward variables
  real(8), intent(out) :: Ref_Pressure ! reference pressure

  ! select flus module operation
  select case(Operation)
  case(1) ! Operation 1 : initiale process
    call InitialProcess(Grid,Mixt,Conf)
    call InterfaceSet(Grid,Conf,nBNode,nBFace,Bxyz,Bf2n,BPatch,BCFlag,InVars,OutVars)
    call InterfaceOut(Grid,Mixt,Conf,nBFace,OutVars,Ref_Pressure)
  case(2) ! Operation 2 : solve process
    call InterfaceIn(Grid,Mixt,DeltaTime,Bxyz,BCFlag,InVars)
    call SolverProcess(Grid,Mixt,Conf,DeltaTime,Bxyz,BCFlag,InVars)
    call InterfaceOut(Grid,Mixt,Conf,nBFace,OutVars,Ref_Pressure)
  case(3) ! Operation 3 : remesh process
    call RemeshProcess(Grid,Mixt,Conf,nBNode,nBFace,Bxyz,Bf2n,BPatch)
    call InterfaceDel(Bxyz,Bf2n,BPatch,BCFlag,InVars,OutVars)
    call InterfaceSet(Grid,Conf,nBNode,nBFace,Bxyz,Bf2n,BPatch,BCFlag,InVars,OutVars)
    call InterfaceOut(Grid,Mixt,Conf,nBFace,OutVars,Ref_Pressure)
  case(4) ! Operation 4 : post process
    call PostProcess(Grid,Mixt,Conf,Iteration)
  case(5) ! Operation 5 : finish process
    call InterfaceDel(Bxyz,Bf2n,BPatch,BCFlag,InVars,OutVars)
    call FinishProcess(Grid,Mixt,Conf)
  end select

end subroutine

end module
