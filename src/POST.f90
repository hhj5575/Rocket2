MODULE SRMS_POST_MOD

USE SRMS_TYPE_MOD
IMPLICIT NONE; PRIVATE
PUBLIC POST

CONTAINS

SUBROUTINE POST

! Post Initial Fluid / Burn / Surface Region
CALL FLUS(Control%Dim,4,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%Patch,Fluid%BcFlag,Fluid%InVars,Fluid%OutVars,Fluid%Ref_Pressure)
!CALL SPBS(Control%Dim,4,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%BcFlag,Fluid%OutVars,Fluid%InVars,Fluid%Prop_Dens,Fluid%Ref_Pressure)
CALL SURFACE(Control%Dim, TYPE1 = 100, FILENUM = iter/Control%OutputIteration)
CALL SURFACE(Control%Dim, TYPE1 = 105, FILENUM = iter)
END SUBROUTINE

END MODULE SRMS_POST_MOD
