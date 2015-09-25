module interface_module

use mpi
use config_module
use grid_module
use mixture_module
implicit none; private
public InterfaceDel
public InterfaceSet
public InterfaceIn
public InterfaceOut

! Number of interface variables
integer, parameter :: nInVar = 2 ! number of inward variables
integer, parameter :: nOutVar = 3 ! number of outward variables

contains

subroutine InterfaceDel(Bxyz,Bf2n,BPatch,BCFlag,InVars,OutVars)
  real(8), intent(inout), allocatable :: Bxyz(:,:)
  integer, intent(inout), allocatable :: Bf2n(:,:)
  integer, intent(inout), allocatable :: BPatch(:)
  integer, intent(inout), allocatable :: BCFlag(:)
  real(8), intent(inout), allocatable :: InVars(:,:)
  real(8), intent(inout), allocatable :: OutVars(:,:)

  if( allocated( Bxyz ) ) deallocate( Bxyz )
  if( allocated( Bf2n ) ) deallocate( Bf2n )
  if( allocated( BPatch ) ) deallocate( BPatch )
  if( allocated( BCFlag ) ) deallocate( BCFlag )
  if( allocated( InVars ) ) deallocate( InVars )
  if( allocated( OutVars ) ) deallocate( OutVars )

end subroutine

subroutine InterfaceSet(Grid,Conf,nBNode,nBFace,Bxyz,Bf2n,BPatch,BCFlag,InVars,OutVars)
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(out) :: nBNode
  integer, intent(out) :: nBFace
  real(8), intent(out), allocatable :: Bxyz(:,:)
  integer, intent(out), allocatable :: Bf2n(:,:)
  integer, intent(out), allocatable :: BPatch(:)
  integer, intent(out), allocatable :: BCFlag(:)
  real(8), intent(out), allocatable :: InVars(:,:)
  real(8), intent(out), allocatable :: OutVars(:,:)

  integer :: rank, size, dest, srcs, ier
  integer :: iDim, iLocal, iNode, iBuffer
  integer :: iBC, iBound ,iDomain, nDomain, GlobalID
  integer :: iNodeBound, iRecvNode, nRecvNode
  integer :: iFaceBound, iRecvFace, nRecvFace

  ! MPI comm array
  real(8), allocatable :: Buffer_Send_xyz(:)
  real(8), allocatable :: Buffer_Recv_xyz(:)
  integer, allocatable :: Buffer_Send_Patch(:)
  integer, allocatable :: Buffer_Recv_Patch(:)
  integer, allocatable :: Buffer_Send_BCFlag(:)
  integer, allocatable :: Buffer_Recv_BCFlag(:)
  integer, allocatable :: Buffer_Send_f2n(:)
  integer, allocatable :: Buffer_Recv_f2n(:)

  ! MPI request / status
  integer :: Send_Req(4), Send_Stat(MPI_STATUS_SIZE,4)
  integer :: Recv_Req(4), Recv_Stat(MPI_STATUS_SIZE,4)

  ! Get rank of each core / number of cores
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)
  call MPI_Comm_size(MPI_COMM_WORLD,size,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Creating global boundary interface'

  ! Set number of domains
  nDomain = size

  ! Allocate global boundary data
  if( rank == 0 ) then ! master core
    nBNode = sum( Grid%BoundRecv(:)%nNode )
    nBFace = sum( Grid%BoundRecv(:)%nFace )
  else ! the others
    nBNode = 0
    nBFace = 0
  end if
  allocate( Bxyz(Grid%nDim,nBNode) )
  allocate( Bf2n(Grid%nNodeFace,nBFace) )
  allocate( BPatch(nBFace) )
  allocate( BCFlag(nBFace) )
  allocate( InVars(nInVar,nBFace) )
  allocate( OutVars(nOutVar,nBFace) )

  ! Allocate buffer send data
  allocate( Buffer_Send_xyz(Grid%nNodeBound*Grid%nDim) )
  allocate( Buffer_Send_Patch(Grid%nFaceBound) )
  allocate( Buffer_Send_BCFlag(Grid%nFaceBound) )
  allocate( Buffer_Send_f2n(Grid%nFaceBound*Grid%nNodeFace) )
  Buffer_Send_xyz(:) = 0
  Buffer_Send_Patch(:) = 0
  Buffer_Send_BCFlag(:) = 0
  Buffer_Send_f2n(:) = 0

  ! Set buffer send node data
  do iNodeBound = 1, Grid%nNodeBound
    do iDim = 1, Grid%nDim
      iBuffer = Grid%nDim*(iNodeBound-1)+iDim
      Buffer_Send_xyz(iBuffer) = Grid%xyz(iDim,iNodeBound)
    end do
  end do

  ! Set buffer send face data
  do iFaceBound = 1, Grid%nFaceBound
    Buffer_Send_Patch(iFaceBound) = Grid%Patch(iFaceBound)
    do iBC = 1, Conf%nBC
      if( Conf%iPatch(iBC) == Grid%Patch(iFaceBound) ) then
        select case(Conf%BCtype(iBC))
        case(10) ! propellant surface (not ignited)
          Buffer_Send_BCFlag(iFaceBound) =  0
        case(13) ! nozzle throat
          Buffer_Send_BCFlag(iFaceBound) = -2
        case default ! the other
          Buffer_Send_BCFlag(iFaceBound) = -1
        end select
      end if
    end do
    do iLocal = 1, Grid%nNodeFace
      iNode = Grid%f2n(iLocal,iFaceBound)
      iBuffer = Grid%nNodeFace*(iFaceBound-1)+iLocal
      Buffer_Send_f2n(iBuffer) = Grid%NodeGID(iNode)
    end do
  end do

  ! Create global boundary interface
  do iBound = 1, Grid%nBound
    iDomain = Grid%BoundList(iBound)

    if( rank == iDomain-1 ) then

      if( rank /= 0 ) then ! send buffer array
        dest = 0 ! master core
        ! Non-blocking send comm
        call MPI_Isend(Buffer_Send_xyz,   Grid%nNodeBound*Grid%nDim,     MPI_REAL8,  dest,1,MPI_COMM_WORLD,Send_Req(1),ier)
        call MPI_Isend(Buffer_Send_Patch, Grid%nFaceBound,               MPI_INTEGER,dest,2,MPI_COMM_WORLD,Send_Req(2),ier)
        call MPI_Isend(Buffer_Send_BCFlag,Grid%nFaceBound,               MPI_INTEGER,dest,3,MPI_COMM_WORLD,Send_Req(3),ier)
        call MPI_Isend(Buffer_Send_f2n,   Grid%nFaceBound*Grid%nNodeFace,MPI_INTEGER,dest,4,MPI_COMM_WORLD,Send_Req(4),ier)
        ! Wait for the set of non-blocking comm
        call MPI_Waitall(4,Send_Req,Send_Stat,ier)
      end if

    end if

    if( rank == 0 ) then

      ! Set number of recv nodes / faces
      nRecvNode = Grid%BoundRecv(iBound)%nNode
      nRecvFace = Grid%BoundRecv(iBound)%nFace

      ! Allocate buffer recv data
      allocate( Buffer_Recv_xyz(nRecvNode*Grid%nDim) )
      allocate( Buffer_Recv_Patch(nRecvFace) )
      allocate( Buffer_Recv_BCFlag(nRecvFace) )
      allocate( Buffer_Recv_f2n(nRecvFace*Grid%nNodeFace) )

      if( iDomain-1 /= 0 ) then
        srcs = iDomain-1 ! from boundary domain
        ! Non-blocking send comm
        call MPI_Irecv(Buffer_Recv_xyz,   nRecvNode*Grid%nDim,     MPI_REAL8,  srcs,1,MPI_COMM_WORLD,Recv_Req(1),ier)
        call MPI_Irecv(Buffer_Recv_Patch, nRecvFace,               MPI_INTEGER,srcs,2,MPI_COMM_WORLD,Recv_Req(2),ier)
        call MPI_Irecv(Buffer_Recv_BCFlag,nRecvFace,               MPI_INTEGER,srcs,3,MPI_COMM_WORLD,Recv_Req(3),ier)
        call MPI_Irecv(Buffer_Recv_f2n,   nRecvFace*Grid%nNodeFace,MPI_INTEGER,srcs,4,MPI_COMM_WORLD,Recv_Req(4),ier)
        ! Wait for the set of non-blocking comm
        call MPI_Waitall(4,Recv_Req,Recv_Stat,ier)
      else ! master core, simply copy buffer
        Buffer_Recv_xyz(:)     = Buffer_Send_xyz(:)
        Buffer_Recv_Patch(:)   = Buffer_Send_Patch(:)
        Buffer_Recv_BCFlag(:)  = Buffer_Send_BCFlag(:)
        Buffer_Recv_f2n(:)     = Buffer_Send_f2n(:)
      end if

      ! Set boundary node data
      do iRecvNode = 1, nRecvNode
        GlobalID = Grid%BoundRecv(iBound)%Node(iRecvNode)
        do iDim = 1, Grid%nDim
          iBuffer = Grid%nDim*(iRecvNode-1)+iDim
          Bxyz(iDim,GlobalID) = Buffer_Recv_xyz(iBuffer)
        end do
      end do

      ! Set boundary face data
      do iRecvFace = 1, nRecvFace
        GlobalID = Grid%BoundRecv(iBound)%Face(iRecvFace)
        BPatch(GlobalID) = Buffer_Recv_Patch(iRecvFace)
        BCFlag(GlobalID) = Buffer_Recv_BCFlag(iRecvFace)
        do iLocal = 1, Grid%nNodeFace
          iBuffer = Grid%nNodeFace*(iRecvFace-1)+iLocal
          Bf2n(iLocal,GlobalID) = Buffer_Recv_f2n(iBuffer)
        end do
      end do

      ! Deallocate buffer recv data
      deallocate( Buffer_Recv_xyz )
      deallocate( Buffer_Recv_Patch )
      deallocate( Buffer_Recv_BCFlag )
      deallocate( Buffer_Recv_f2n )

    end if

  end do

  ! Deallocate buffer send data
  deallocate( Buffer_Send_xyz )
  deallocate( Buffer_Send_Patch )
  deallocate( Buffer_Send_BCFlag )
  deallocate( Buffer_Send_f2n )

  ! MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

subroutine InterfaceIn(Grid,Mixt,DeltaTime,Bxyz,BCFlag,InVars)
  type(t_Grid), intent(inout) :: Grid
  type(t_Mixt), intent(inout) :: Mixt
  real(8), intent(in) :: DeltaTime
  real(8), intent(in) :: Bxyz(:,:)
  integer, intent(in) :: BCFlag(:)
  real(8), intent(in) :: InVars(:,:)

  integer :: rank, dest, srcs, ier
  integer :: iDim, iLocal, iBuffer
  integer :: iBound, iDomain, GlobalID
  integer :: iNodeBound, iSendNode, nSendNode
  integer :: iFaceBound, iSendFace, nSendFace
  real(8) :: inv_DeltaTime

  ! MPI comm buffer array
  real(8), allocatable :: Buffer_Send_Bxyz(:)
  real(8), allocatable :: Buffer_Recv_Bxyz(:)
  integer, allocatable :: Buffer_Send_BCFlag(:)
  integer, allocatable :: Buffer_Recv_BCFlag(:)
  real(8), allocatable :: Buffer_Send_InVars(:)
  real(8), allocatable :: Buffer_Recv_InVars(:)

  ! MPI comm request / stat
  integer :: Send_Req(3), Send_Stat(MPI_STATUS_SIZE,3)
  integer :: Recv_Req(3), Recv_Stat(MPI_STATUS_SIZE,3)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Transfer interface data (global => local)'

  ! Allocate / initialize buffer recv data
  allocate( Buffer_Recv_Bxyz(Grid%nNodeBound*Grid%nDim) )
  allocate( Buffer_Recv_BCFlag(Grid%nFaceBound) )
  allocate( Buffer_Recv_InVars(Grid%nFaceBound*nInVar) )
  Buffer_Recv_Bxyz(:) = 0
  Buffer_Recv_BCFlag(:) = 0
  Buffer_Recv_InVars(:) = 0

  do iBound = 1, Grid%nBound
    iDomain = Grid%BoundList(iBound)

    if( rank == 0 ) then ! master core

      ! Set number of send nodes / faces
      nSendNode = Grid%BoundRecv(iBound)%nNode
      nSendFace = Grid%BoundRecv(iBound)%nFace

      ! Allocate / initialize buffer send data
      allocate( Buffer_Send_Bxyz(nSendNode*Grid%nDim) )
      allocate( Buffer_Send_BCFlag(nSendFace) )
      allocate( Buffer_Send_InVars(nSendFace*nInVar) )
      Buffer_Send_Bxyz(:)  = 0
      Buffer_Send_BCFlag(:) = 0
      Buffer_Send_InVars(:) = 0

      ! Set boundary node data
      do iSendNode = 1, nSendNode
        GlobalID = Grid%BoundRecv(iBound)%Node(iSendNode)
        do iDim = 1, Grid%nDim
          iBuffer = Grid%nDim*(iSendNode-1)+iDim
          Buffer_Send_Bxyz(iBuffer) = Bxyz(iDim,GlobalID)
        end do
      end do

      ! Set boundary face data
      do iSendFace = 1, nSendFace
        GlobalID = Grid%BoundRecv(iBound)%Face(iSendFace)
        Buffer_Send_BCFlag(iSendFace) = BCFlag(GlobalID)
        do iLocal = 1, nInVar
          iBuffer = 2*(iSendFace-1)+iLocal
          Buffer_Send_InVars(iBuffer) = InVars(iLocal,GlobalID)
        end do
      end do

      if( iDomain-1 /= 0 ) then ! sent buffer array
        dest = iDomain-1 ! to boundary domain
        call MPI_Isend(Buffer_Send_Bxyz,  nSendNode*Grid%nDim,MPI_REAL8,  dest,1,MPI_COMM_WORLD,Send_Req(1),ier)
        call MPI_Isend(Buffer_Send_BCFlag,nSendFace,          MPI_INTEGER,dest,2,MPI_COMM_WORLD,Send_Req(2),ier)
        call MPI_Isend(Buffer_Send_InVars,nSendFace*nInVar,   MPI_REAL8,  dest,3,MPI_COMM_WORLD,Send_Req(3),ier)
        ! Wait for the set of non-blocking comm
        call MPI_Waitall(3,Send_Req,Send_Stat,ier)
      else ! iDomain is master core. so, simply copy buffer send
        Buffer_Recv_Bxyz(:)   = Buffer_Send_Bxyz(:)
        Buffer_Recv_BCFlag(:) = Buffer_Send_BCFlag(:)
        Buffer_Recv_InVars(:) = Buffer_Send_InVars(:)
      end if

      ! Deallocate buffer send data
      deallocate( Buffer_Send_Bxyz )
      deallocate( Buffer_Send_BCFlag )
      deallocate( Buffer_Send_InVars )

    end if

    if( rank == iDomain-1 ) then

      if( rank /= 0 ) then ! receive buffers array
        srcs = 0 ! from master
        call MPI_Irecv(Buffer_Recv_Bxyz,  Grid%nNodeBound*Grid%nDim,MPI_REAL8,  srcs,1,MPI_COMM_WORLD,Recv_Req(1),ier)
        call MPI_Irecv(Buffer_Recv_BCFlag,Grid%nFaceBound,          MPI_INTEGER,srcs,2,MPI_COMM_WORLD,Recv_Req(2),ier)
        call MPI_Irecv(Buffer_Recv_InVars,Grid%nFaceBound*nInVar,   MPI_REAL8,  srcs,3,MPI_COMM_WORLD,Recv_Req(3),ier)
        ! Wait for the set of non-blocking comm
        call MPI_Waitall(3,Recv_Req,Recv_Stat,ier)
      end if

      ! Set boundary node data
      inv_DeltaTime = 1d0 / DeltaTime
      do iNodeBound = 1, Grid%nNodeBound
        do iDim = 1, Grid%nDim
          iBuffer = Grid%nDim*(iNodeBound-1)+iDim
          Grid%gv(iDim,iNodeBound) = ( Buffer_Recv_Bxyz(iBuffer) - Grid%xyz(iDim,iNodeBound) ) * inv_DeltaTime
        end do
      end do

      ! Set boundary face data
      do iFaceBound = 1, Grid%nFaceBound
        iBuffer = nInVar*(iFaceBound-1)
        Mixt%bFlag(iFaceBound) = Buffer_Recv_BCFlag(iFaceBound)
        Mixt%mFlux(iFaceBound) = Buffer_Recv_InVars(iBuffer+1)
        Mixt%fTemp(iFaceBound) = Buffer_Recv_InVars(iBuffer+2)
      end do

    end if

  end do

  ! Deallocate buffer recv data
  deallocate( Buffer_Recv_Bxyz )
  deallocate( Buffer_Recv_BCFlag )
  deallocate( Buffer_Recv_InVars )

  ! MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

subroutine InterfaceOut(Grid,Mixt,Conf,nBFace,OutVars,Ref_Pressure)
  type(t_Grid), intent(in) :: Grid
  type(t_Mixt), intent(in) :: Mixt
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: nBFace
  real(8), intent(out) :: OutVars(nOutVar,nBFace)
  real(8), intent(out) :: Ref_Pressure

  integer :: rank, dest, srcs, ier
  integer :: iCell, iFace, iBound, iDomain, iBC
  integer :: iBuffer, GlobalID, iList, List1, List2
  integer :: iFaceBound, iRecvFace, nRecvFace
  real(8) :: Temperature, Pressure, Density, NormVel, TangVel
  real(8) :: Area, p_sum, p_sum_total, a_sum, a_sum_total
  real(8), parameter :: pi = 4d0*datan(1d0)

  ! MPI comm variable
  real(8), allocatable :: Buffer_Send_OutVars(:)
  real(8), allocatable :: Buffer_Recv_OutVars(:)

  ! MPI request / status
  integer :: Send_Req, Send_Stat(MPI_STATUS_SIZE)
  integer :: Recv_Req, Recv_Stat(MPI_STATUS_SIZE)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Transfer interface data (local => global)'

  ! Allocate / initialize buffer send data
  allocate( Buffer_Send_OutVars(Grid%nFaceBound*nOutVar) )
  Buffer_Send_OutVars(:) = 0d0

  ! Set outward variables (Pressure/Temperature/Tangential mass velocity)
  do iFaceBound = 1, Grid%nFaceBound
    iCell = Grid%f2c(1,iFaceBound)
    Temperature = Mixt%pv(0,iCell)
    Pressure = Mixt%pv(Grid%nDim+1,iCell)
    Density = Pressure / ( Mixt%GasConstant * Temperature )
    NormVel = sum( Mixt%pv(1:Grid%nDim,iCell)*Grid%fn(:,iFaceBound) )
    TangVel = dsqrt( dmax1( sum( Mixt%pv(1:Grid%nDim,iCell)**2 ) - NormVel**2, 0d0 ) )
    iBuffer = nOutVar*(iFaceBound-1)
    Buffer_Send_OutVars(iBuffer+1) = Pressure
    Buffer_Send_OutVars(iBuffer+2) = Temperature
    Buffer_Send_OutVars(iBuffer+3) = Density * TangVel
  end do

  ! Collet to master core
  do iBound = 1, Grid%nBound
    iDomain = Grid%BoundList(iBound)

    if( rank == iDomain-1 ) then
      if( iDomain-1 /= 0 ) then
        dest = 0 ! to master
        ! Non-blocking send comm
        call MPI_Isend(Buffer_Send_OutVars,Grid%nFaceBound*nOutVar,MPI_REAL8,dest,1,MPI_COMM_WORLD,Send_Req,ier)
        ! Wait for the non-blocking comm
        call MPI_Wait(Send_Req,Send_Stat,ier)
      end if
    end if

    if( rank == 0 ) then

      ! Set number of recv nodes / faces
      nRecvFace = Grid%BoundRecv(iBound)%nFace

      ! Allocate buffer recv data
      allocate( Buffer_Recv_OutVars(nRecvFace*nOutVar) )

      if( iDomain-1 /= 0 ) then
        srcs = iDomain-1 ! from boundary domain
        ! Non-blocking recv comm
        call MPI_Irecv(Buffer_Recv_OutVars,nRecvFace*nOutVar,MPI_REAL8,srcs,1,MPI_COMM_WORLD,Recv_Req,ier)
        ! Wait for the non-blocking comm
        call MPI_Wait(Recv_Req,Recv_Stat,ier)
      else ! just copy buffer
        Buffer_Recv_OutVars(:) = Buffer_Send_OutVars(:)
      end if

      ! Set outward variables from buffer recv
      do iRecvFace = 1, nRecvFace
        GlobalID = Grid%BoundRecv(iBound)%Face(iRecvFace)
        iBuffer = nOutVar*(iRecvFace-1)
        OutVars(1,GlobalID) = Buffer_Recv_OutVars(iBuffer+1)
        OutVars(2,GlobalID) = Buffer_Recv_OutVars(iBuffer+2)
        OutVars(3,GlobalID) = Buffer_Recv_OutVars(iBuffer+3)
      end do

      ! Deallocate buffer recv data
      deallocate( Buffer_Recv_OutVars )

    end if

  end do

  ! Deallocate buffer send data
  deallocate( Buffer_Send_OutVars )

  ! Set reference pressure
  do iBC = 1, Conf%nBC
    if( Conf%BCType(iBC) == 12 .and. Conf%bv(1,iBC) > 0d0 ) then
      List1 = Grid%ListIndex(iBC)
      List2 = Grid%ListIndex(iBC+1)-1
      p_sum = 0d0; a_sum = 0d0; Ref_Pressure = 0
      do iList = List1, List2
        iFace = Grid%List(iList)
        iCell = Grid%f2c(1,iFace)
        if( Conf%Axisymmetric ) then
          Area = Grid%fa(iFace) * 2d0 * pi * Grid%fc(2,iFace)
        else ! not axisymmetric case
          Area = Grid%fa(iFace)
        end if
        p_sum = p_sum + Mixt%pv(Grid%nDim+1,iCell) * Area
        a_sum = a_sum + Area
      end do
      call MPI_Allreduce(p_sum,p_sum_total,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call MPI_Allreduce(a_sum,a_sum_total,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      if( a_sum_total > 0 ) Ref_Pressure = p_sum_total / a_sum_total
      if( rank == 0 ) write(*,*) 'FLUS >> Reference Pressure at Monitor :', Ref_Pressure
    end if
  end do

  ! MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

end module
