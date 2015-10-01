MODULE SRMS_INITIAL_MOD

USE SRMS_TYPE_MOD
IMPLICIT NONE; PRIVATE
PUBLIC INITIAL

CONTAINS

SUBROUTINE INITIAL

! Read Input File
CALL READINPUT

Control%Dim = 3
! Initialize Solid Region
CALL FEM_GXL(Control%Dim, 0, iter_solid, dTmin, n_remesh, remesh_num, &
               SProp%Patch, SProp%Remesh_flag, SProp%corner_index, SProp%edge_length, &
               SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, SProp%Disp, &
               SProp%ridge_group_num, SProp%ridge_num, SProp%ridge, &
               SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad, SCase%Disp)

IF(n_remesh == 2) THEN
  CALL FEM_GXL(Control%Dim, 5, iter_solid, dTmin, n_remesh, remesh_num, &
               SProp%Patch, SProp%Remesh_flag, SProp%corner_index, SProp%edge_length, &
               SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, SProp%Disp, &
               SProp%ridge_group_num, SProp%ridge_num, SProp%ridge, &
               SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad, SCase%Disp)
ELSE
  CALL FEM_GXL(Control%Dim, 1, iter_solid, dTmin, n_remesh, remesh_num, &
               SProp%Patch, SProp%Remesh_flag, SProp%corner_index, SProp%edge_length, &
               SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, SProp%Disp, &
               SProp%ridge_group_num, SProp%ridge_num, SProp%ridge, &
               SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad, SCase%Disp)
END IF

! Initialize Fluid Region
CALL FLUS(Control%Dim,1,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%Patch,Fluid%BcFlag,Fluid%InVars,Fluid%OutVars,Fluid%Ref_Pressure)

! Initialize Burn Region
CALL SPBS(Control%Dim,1,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%BcFlag,Fluid%OutVars,Fluid%InVars,Fluid%Prop_Dens,Fluid%Ref_Pressure)

! Initialize Surface Region
CALL SURFACE(Control%Dim, TYPE1 = -1, TIMESTEP = dTmin, &
            F_POINT_NUM = Fluid%nBNode, F_POINT = Fluid%xyz, F_EDGE_NUM = Fluid%nBFace, F_EDGE = Fluid%f2n, F_LOC = Fluid%Patch, F_BCFLAG = Fluid%BcFlag, &
            P_POINT_NUM = SProp%nBNode, P_POINT = SProp%xyz, P_EDGE_NUM = SProp%nBFace, P_EDGE = SProp%f2n, P_LOC = SProp%Patch, P_CORNER_INDEX = SProp%Corner_Index, &
	    P_RIDGE_GROUP_NUM = SProp%Ridge_Group_Num, P_RIDGE_NUM = SProp%Ridge_num, P_RIDGE = SProp%Ridge, &
            C_POINT_NUM = SCase%nBNode, C_POINT = SCase%xyz, C_EDGE_NUM = SCase%nBFace, C_EDGE = SCase%f2n)
CALL SURFACE(Control%Dim, TYPE1 = 104, FILENUM = 0)

END SUBROUTINE

SUBROUTINE READINPUT

INTEGER :: io
CHARACTER(99) :: dummy

! Read Input
OPEN(newunit=io, file='./input/test.inp')
READ(io,*)
READ(io,*)
READ(io,*) dummy, dummy, Control%Dim
READ(io,*)
READ(io,*) dummy, dummy, Control%InitialTime
READ(io,*) dummy, dummy, Control%TargetTime
READ(io,*) dummy, dummy, Control%SystemTimeStep
READ(io,*)
READ(io,*) dummy, dummy, Control%IterMax
READ(io,*) dummy, dummy, Control%OutputIteration
READ(io,*)
READ(io,*) dummy, dummy, regression_rate_x
READ(io,*) dummy, dummy, restart_iteration
CLOSE(io)

! Set Time Variables
iter  = 0
Time  = Control%InitialTime
dTmin = Control%SystemTimeStep

END SUBROUTINE

END MODULE SRMS_INITIAL_MOD
