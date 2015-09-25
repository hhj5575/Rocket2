MODULE SPBS_REGION_MOD

USE SPBS_TYPE_MOD
IMPLICIT NONE; PRIVATE
TYPE(t_Region), PUBLIC, TARGET :: Region
PUBLIC ALLOC
PUBLIC FINISH

CONTAINS

SUBROUTINE ALLOC(pLGrid,pBurn)

! In or Out Variables
TYPE(t_LGrid), INTENT(in   ) :: pLGrid
TYPE(t_Burn) , INTENT(  out) :: pBurn
! Local Variables
INTEGER :: nNode

nNode = pLGrid%nNode

ALLOCATE( pBurn%T(nNode) )
ALLOCATE( pBurn%T0(nNode) )

END SUBROUTINE

SUBROUTINE FINISH

! Local Variables
TYPE(t_SGrid), POINTER :: pSGrid
TYPE(t_LGrid), POINTER :: pLGrid

! Set Pointer
pSGrid => Region%SGrid
pLGrid => Region%LGrid

! Deallocate SGrid
DEALLOCATE( pSGrid%Bxyz )
DEALLOCATE( pSGrid%Bf2n )

! Deallocate LGrid
DEALLOCATE( pLGrid%x )
DEALLOCATE( pLGrid%z )
DEALLOCATE( pLGrid%dzdx )
DEALLOCATE( pLGrid%dzdx2 )

! Deallocate Burn
DEALLOCATE( Region%Burn )

END SUBROUTINE

END MODULE SPBS_REGION_MOD
