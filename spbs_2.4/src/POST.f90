MODULE SPBS_POST_MOD

USE SPBS_TYPE_MOD
USE SPBS_REGION_MOD, ONLY : Region
IMPLICIT NONE; PRIVATE
INTEGER, SAVE :: PostNum = -1
REAL(8) :: PostTime
PUBLIC POST
CONTAINS

SUBROUTINE POST(InitialTime)

! In or Out Variables
REAL(8), INTENT(in   ) :: InitialTime
! Local Variables
TYPE(t_SGrid), POINTER :: pSGrid
TYPE(t_LGrid), POINTER :: pLGrid
TYPE(t_Burn) , POINTER :: pBurn(:)

! Set Post Time / Number
PostTime = InitialTime
PostNum = PostNum + 1

! Set Pointer
pSGrid => Region%SGrid
pLGrid => Region%LGrid
pBurn  => Region%Burn

CALL PLTPOST_ASCII(pSGrid,pBurn)
WRITE(*,*) 'SPBS >> Burning region is postprocessed at Time', PostTime

END SUBROUTINE

SUBROUTINE PLTPOST_ASCII(pSGrid,pBurn)

! In or Out Variables
TYPE(t_SGrid), INTENT(in   ) :: pSGrid
TYPE(t_Burn) , INTENT(in   ) :: pBurn(:)
! Local Variables
INTEGER :: io
INTEGER :: Count
INTEGER :: Node, BFace
CHARACTER(99) :: cNum, Fname

WRITE(cNum,"(i6.6)") PostNum
Fname = "./output/burn/burn_"//TRIM(cNum)//".plt"

OPEN(newunit=io,file=TRIM(Fname))

SELECT CASE(pSGrid%nDim)
CASE(2)
  WRITE(io,*) 'TITLE = Result'
  WRITE(io,*) 'variables ="x","y","Temp"'
  WRITE(io,*) 'ZONE t="Burn Region", NODES=',pSGrid%nBNode,'ELEMENTS=', pSGrid%nBFace
  WRITE(io,*) 'VARLOCATION=([3]=CELLCENTERED)'
  WRITE(io,*) 'ZONETYPE=FELINESEG,DATAPACKING=BLOCK,SOLUTIONTIME=',PostTime
CASE(3)
  WRITE(io,*) 'TITLE = Result'
  WRITE(io,*) 'variables ="x","y","z","Temp"'
  WRITE(io,*) 'ZONE t="Burn Region", NODES=',pSGrid%nBNode,'ELEMENTS=', pSGrid%nBFace
  WRITE(io,*) 'VARLOCATION=([4]=CELLCENTERED)'
  WRITE(io,*) 'ZONETYPE=FETRIANGLE,DATAPACKING=BLOCK,SOLUTIONTIME=',PostTime
END SELECT

DO Count = 1, pSGrid%nDim
  DO Node = 1, pSGrid%nBNode
    WRITE(io,*) pSGrid%Bxyz(Count,Node)
  END DO
END DO

DO BFace = 1, pSGrid%nBFace
  IF( pBurn(BFace)%Ignited < 0 ) THEN
    WRITE(io,*) 0d0
  ELSE
    WRITE(io,*) pBurn(BFace)%T(1)
  END IF
END DO

DO BFace = 1 , pSGrid%nBFace
  WRITE(io,*) pSGrid%Bf2n(:,BFace)
END DO

CLOSE(io)

END SUBROUTINE

!SUBROUTINE VTKPOST_ASCII(pGrid,pMixt)
!
!! In or Out Variables
!TYPE(t_Grid) , INTENT(in   ) :: pGrid
!TYPE(t_Mixt) , INTENT(in   ) :: pMixt
!! Local Variables
!INTEGER :: io
!INTEGER :: Node, Cell
!CHARACTER(99) :: cNum, Fname
!
!WRITE(cNum,"(i6.6)") PostNum
!Fname = "./output/fluid_"//TRIM(cNum)//".vtk"
!
!OPEN(newunit=io,file=TRIM(Fname))
!WRITE(io,'(a26)') '# vtk DataFile Version 3.0'
!WRITE(io,*) 'Result file'
!WRITE(io,'(a5)' ) 'ASCII'
!WRITE(io,'(a25)') 'DATASET UNSTRUCTURED_GRID'
!WRITE(io,*) 'FIELD FIELD_DATA 1'
!WRITE(io,*) 'TIME 1 1 DOUBLE'
!WRITE(io,*) PostTime
!WRITE(io,*) 'POINTS', pGrid%nNode, 'DOUBLE'
!DO Node = 1, pGrid%nNode
!  WRITE(io,*) pGrid%xyz(:,Node), 0d0
!END DO
!WRITE(io,*)
!WRITE(io,*) 'CELLS', pGrid%nCell, 4 * pGrid%nCell
!DO Cell = 1, pGrid%nCell
!  WRITE(io,*) 3, pGrid%c2n(:,Cell) - 1
!END DO
!WRITE(io,*) 'CELL_TYPES', pGrid%nCell
!DO Cell = 1, pGrid%nCell
!  WRITE(io,*) 5
!END DO
!WRITE(io,*) 'CELL_DATA', pGrid%nCell
!WRITE(io,*) 'SCALARS DENSITY DOUBLE 1'
!WRITE(io,*) 'LOOKUP_TABLE DEFAULT'
!DO Cell = 1, pGrid%nCell
!  WRITE(io,*) pMixt%pv(1,Cell)
!END DO
!WRITE(io,*) 'VECTORS VELOCITY DOUBLE'
!DO Cell = 1, pGrid%nCell
!  WRITE(io,*) pMixt%pv(2:3,Cell), 0d0
!END DO
!WRITE(io,*) 'SCALARS PRESSURE DOUBLE 1'
!WRITE(io,*) 'LOOKUP_TABLE DEFAULT'
!DO Cell = 1, pGrid%nCell
!  WRITE(io,*) pMixt%pv(4,Cell)
!END DO
!WRITE(io,*) 'SCALARS MACH DOUBLE 1'
!WRITE(io,*) 'LOOKUP_TABLE DEFAULT'
!DO Cell = 1, pGrid%nCell
!  WRITE(io,*) DSQRT( pMixt%pv(2,Cell)**2 + pMixt%pv(3,Cell)**2 ) / DSQRT( pMixt%gamma * pMixt%pv(4,Cell) / pMixt%pv(1,Cell) )
!END DO
!WRITE(io,*) 'SCALARS TEMPERATURE DOUBLE 1'
!WRITE(io,*) 'LOOKUP_TABLE DEFAULT'
!DO Cell = 1, pGrid%nCell
!  WRITE(io,*) pMixt%pv(4,Cell) / ( pMixt%pv(1,Cell) * pMixt%Rs )
!END DO
!CLOSE(io)
!
!END SUBROUTINE

END MODULE SPBS_POST_MOD
