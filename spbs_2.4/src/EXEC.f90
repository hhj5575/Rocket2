PROGRAM EXEC

IMPLICIT NONE
INTEGER :: nDim
INTEGER :: nBNode, nBFace
INTEGER :: iter
REAL(8) :: Time, dTime, TargetTime
REAL(8), ALLOCATABLE :: Bxyz(:,:)
INTEGER, ALLOCATABLE :: Bf2n(:,:)
REAL(8), ALLOCATABLE :: InVars(:,:)
REAL(8), ALLOCATABLE :: OutVars(:,:)
INTEGER, ALLOCATABLE :: bFlag(:)
REAL(8) :: Prop_Dens
REAL(8) :: Ref_Pressure
LOGICAL :: IgniteAPN

nDim = 2
nBNode = 6
nBFace = 4

ALLOCATE( Bxyz(nDim,nBNode) )
ALLOCATE( Bf2n(nDim,nBFace) )
ALLOCATE( InVars(3,nBFace) )
ALLOCATE( OutVars(2,nBFace) )
ALLOCATE( bFlag(nBFace) )

Bxyz(1,1) = 0d0
Bxyz(2,1) = 0d0
Bxyz(1,2) = 1d0
Bxyz(2,2) = 1d0
Bxyz(1,3) = 2d0
Bxyz(2,3) = 1d0
Bxyz(1,4) = 3d0
Bxyz(2,4) = 2d0
Bxyz(1,5) = 2d0
Bxyz(2,5) = 0d0
Bxyz(1,6) = 3d0
Bxyz(2,6) = 0d0


Bf2n(1,1) = 1
Bf2n(2,1) = 2
Bf2n(1,2) = 2
Bf2n(2,2) = 3
Bf2n(1,3) = 3
Bf2n(2,3) = 4
Bf2n(1,4) = 5
Bf2n(2,4) = 6

bFlag(1) = -1
bFlag(2) = -1
bFlag(3) = 0
bFlag(4) = -1

Time       = 0d0
dTime      = 0.0001d0
TargetTime = 0.001d0

Ref_Pressure = 101325

! Initialize SPBS Module
CALL SPBS(nDim,1,Time,dTime,nBNode,nBFace,Bxyz,Bf2n,bFlag,InVars,OutVars,Prop_Dens,Ref_Pressure,IgniteAPN)

DO iter = 1, 200

  WRITE(*,*)
  WRITE(*,*) 'Outer Iteration =', iter
  IF( Time + dTime >= TargetTime ) dTime = TargetTime - Time
  ! Solve SPBS Module
  InVars(1,:) = 6894733.26d0
  InVars(2,:) = 1000d0
  InVars(3,:) = 50d0
  InVars(3,4) = 25d0

  IF( iter == 5 ) THEN
    nBNode = 4
    nBFace = 3
    CALL SPBS(nDim,3,Time,0d0,nBNode,nBFace,Bxyz(:,:nBNode),Bf2n(:,:nBFace),bFlag(:nBFace),InVars(:,:nBFace),OutVars(:,:nBFace),Prop_Dens,Ref_Pressure,IgniteAPN)
  END IF

  CALL SPBS(nDim,2,Time,dTime,nBNode,nBFace,Bxyz(:,:nBNode),Bf2n(:,:nBFace),bFlag(:nBFace),InVars(:,:nBFace),OutVars(:,:nBFace),Prop_Dens,Ref_Pressure,IgniteAPN)
  Time = Time + dTime

  CALL SPBS(nDim,4,Time,dTime,nBNode,nBFace,Bxyz(:,:nBNode),Bf2n(:,:nBFace),bFlag(:nBFace),InVars(:,:nBFace),OutVars(:,:nBFace),Prop_Dens,Ref_Pressure,IgniteAPN)

  IF( Time == TargetTime ) EXIT

END DO

CALL SPBS(nDim,5,Time,dTime,nBNode,nBFace,Bxyz,Bf2n,bFlag,InVars,OutVars,Prop_Dens,Ref_Pressure,IgniteAPN)

END PROGRAM
