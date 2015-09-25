module time_module

use config_module
use grid_module
use mixture_module
implicit none; private
public GridMovingTest
public ExplicitEuler
public TVDRK3
public LUSGS

contains

subroutine GridMovingTest(Grid,Mixt)
  type(t_Grid), intent(inout) :: Grid
  type(t_Mixt), intent(in) :: Mixt

  ! Update grid
  call UpdateGrid(Grid,Mixt%TimeStepMin)
  ! Update metric
  call SetFaceMetric(Grid)
  call SetCellMetric(Grid)

end subroutine

subroutine ExplicitEuler(Grid,Mixt,Conf)
  type(t_Grid), intent(inout) :: Grid
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Conf), intent(in) :: Conf

  integer :: iCell

  ! Set face velocity for ALE
  if( Conf%MoveGrid ) then
    call SetFaceVelocity(Grid,Mixt%TimeStepMin)
  end if

  ! Compute residual
  call ComputeResidual(Mixt,Grid,Conf)

  ! Update conservative variables
  do iCell = 1, Grid%nCell
    Mixt%cv(:,iCell) = Mixt%cv(:,iCell) - Mixt%TimeStep(iCell) * Mixt%Residual(:,iCell)
  end do

  ! Update grid / metric
  if( Conf%MoveGrid ) then
    call UpdateGrid(Grid,Mixt%TimeStepMin)
    call SetFaceMetric(Grid)
    call SetCellMetric(Grid)
    if( Conf%FlowModel == 3 ) then
      call SetWallDistance(Grid,Conf)
    end if
  end if

  ! Postprocessing for mixture
  call Postprocessing(Mixt,Grid,Conf)

end subroutine

subroutine TVDRK3(Grid,Mixt,Conf)
  type(t_Grid), intent(inout) :: Grid
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Conf), intent(in) :: Conf

  integer :: iStage, iCell
  real(8), parameter :: inv3 = 1d0 / 3d0
  real(8), parameter :: Coeff(3) = (/1d0,-0.5d0,0.5d0/)

  ! Set old conservative variables
  Mixt%cv0(:,:) = Mixt%cv(:,:)

  do iStage = 1, 3

    ! Set face velocity for ALE
    if( Conf%MoveGrid  ) then
      call SetFaceVelocity(Grid,Coeff(iStage)*Mixt%TimeStepMin)
    end if

    ! Compute residual
    call ComputeResidual(Mixt,Grid,Conf)

    ! Update conservative variables
    do iCell = 1, Grid%nCell
      select case(iStage)
      case(1); Mixt%cv(:,iCell) = Mixt%cv(:,iCell) - Mixt%TimeStep(iCell) * Mixt%Residual(:,iCell)
      case(2); Mixt%cv(:,iCell) = ( 3d0 * Mixt%cv0(:,iCell) + Mixt%cv(:,iCell) - Mixt%TimeStep(iCell) * Mixt%Residual(:,iCell) ) * 0.25d0
      case(3); Mixt%cv(:,iCell) = ( Mixt%cv0(:,iCell) + 2d0 * Mixt%cv(:,iCell) - Mixt%TimeStep(iCell) * Mixt%Residual(:,iCell) * 2d0 ) * inv3
      end select
    end do

    ! Update grid / metric
    if( Conf%MoveGrid ) then
      call UpdateGrid(Grid,Coeff(iStage)*Mixt%TimeStepMin)
      call SetFaceMetric(Grid)
      call SetCellMetric(Grid)
      if( Conf%FlowModel == 3 ) then
        call SetWallDistance(Grid,Conf)
      end if
    end if

    ! Postprocessing for mixture
    call Postprocessing(Mixt,Grid,Conf)

  end do

end subroutine

subroutine LUSGS(Grid,Mixt,Conf)
  type(t_Grid), intent(in) :: Grid
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Conf), intent(in) :: Conf

  integer :: iLocal, sign
  integer :: iFace, iCell, iNeiCell
  real(8) :: Ldw(0:Mixt%nCVar-1)
  real(8) :: Udw(0:Mixt%nCVar-1)
  real(8) :: df(0:Mixt%nCVar-1)

  ! Compute Residual
  call ComputeResidual(Mixt,Grid,Conf)

  ! Compute diagonal term
  call ComputeDiagonalTerm(Mixt,Grid,Conf)

  ! Initialize dwst / dw
  Mixt%dwst(:,:) = 0d0
  Mixt%dw(:,:) = 0d0

  ! Lower Sweep of LUSGS
  do iCell = 1, Grid%nCell
    Ldw(:) = 0d0
    do iLocal = 1, Grid%nFaceCell
      iFace = Grid%c2f(iLocal,iCell)
      ! Set neighboring cell & sign
      if ( Grid%f2c(1,iFace) == iCell ) then
        sign =  1d0; iNeiCell = Grid%f2c(2,iFace)
      else
        sign = -1d0; iNeiCell = Grid%f2c(1,iFace)
      end if
      ! Compute delta flux and update lower triangle matrix
      if ( iNeiCell < iCell ) then
        call ComputeDeltaFlux(Mixt,Grid,iNeiCell,iFace,Mixt%dwst(:,iNeiCell),df(:))
        Ldw(:) = Ldw(:) + 0.5d0 * ( sign * df(:) - Mixt%ram(:,iFace) * Mixt%dwst(:,iNeiCell) )
      end if
    end do
    ! Update middle step delta conservative variables
    Mixt%dwst(:,iCell) = ( -Mixt%Residual(:,iCell) - Ldw(:) ) / Mixt%diag(:,iCell)
  end do

  ! Upper Sweep of LUSGS
  do iCell = Grid%nCell, 1, -1
    Udw(:) = 0d0
    do iLocal = 1, Grid%nFaceCell
      iFace = Grid%c2f(iLocal,iCell)
      ! Set neighboring cell & sign
      if ( Grid%f2c(1,iFace) == iCell ) then
        sign =  1d0; iNeiCell = Grid%f2c(2,iFace)
      else
        sign = -1d0; iNeiCell = Grid%f2c(1,iFace)
      end if
      ! Compute delta flux and update upper triangle matrix
      if ( iNeiCell > iCell ) then
        call ComputeDeltaFlux(Mixt,Grid,iNeiCell,iFace,Mixt%dw(:,iNeiCell),df(:))
        Udw(:) = Udw(:) + 0.5d0 * ( sign * df(:) - Mixt%ram(:,iFace) * Mixt%dw(:,iNeiCell) )
      end if
    end do
    ! Update final step delta conservative variables
    Mixt%dw(:,iCell) = Mixt%dwst(:,iCell) - Udw(:) / Mixt%diag(:,iCell)
  end do

  ! Update conservative variables
  do iCell = 1, Grid%nCell
    Mixt%cv(:,iCell) = Mixt%cv(:,iCell) + Mixt%dw(:,iCell) * Grid%Vol(iCell)
  end do

  ! Postprocessing for mixture
  call Postprocessing(Mixt,Grid,Conf)

end subroutine

end module
