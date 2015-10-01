PROGRAM SRMS

USE OMP_LIB
USE SRMS_TYPE_MOD
USE SRMS_INITIAL_MOD
USE SRMS_POST_MOD
IMPLICIT NONE
INTEGER :: tid_max
INTEGER :: I, SP
REAL(8) :: press, b_rate
character(len=50) :: fn

tid_max = OMP_GET_MAX_THREADS()

! Initialize
CALL INITIAL

! Initial Post
CALL POST

!dTmin = 0.001
!Control%IterMax = 10000
!Control%Dim = 3
! Solve Fluid Region
SP = restart_iteration
Time = SP*dTmin
DO iter = SP+1, Control%IterMax
  ! Time Step Check
  If( Time + dTmin >= Control%TargetTime ) dTmin = Control%TargetTime - Time
  ! Print Current System Step
  WRITE(*,*)
  WRITE(*,*) 'SRMS >> ======================================================'
  WRITE(*,*) 'SRMS >> Start System Iteration Step'
  WRITE(*,*) 'SRMS >> System Iteration   =', iter
  WRITE(*,*) 'SRMS >> System Time        =', Time
  WRITE(*,*) 'SRMS >> System Time Step   =', dTmin
  WRITE(*,*) 'SRMS >> ======================================================'
  WRITE(*,*)

CALL SURFACE(Control%Dim, TYPE1 = 50)

  ! Solid Surface Remesh
  CALL SURFACE(Control%Dim, TYPE1 = 1, TYPE2 = 1, &
               P_POINT_NUM = SProp%nBNode, P_POINT = SProp%xyz, P_EDGE_NUM = SProp%nBFace, P_EDGE = SProp%f2n, P_LOC = SProp%Patch, &
               P_EDGE_LENGTH = SProp%Edge_Length, P_CORNER_INDEX = SProp%Corner_Index, P_RIDGE_GROUP_NUM = SProp%Ridge_Group_Num, P_RIDGE_NUM = SProp%Ridge_num, &
	       P_RIDGE = SProp%Ridge, FLAG = SProp%n_remesh, FLAG_ARRAY = SProp%Remesh_flag)
  IF(SProp%n_remesh == .TRUE.) THEN
    CALL SURFACE(Control%Dim, TYPE1 = 100, FILENUM = 10000)
    n_remesh = 3
  END IF

  ! Solid Movement Process
  CALL FEM_GXL(Control%Dim, 3, iter_solid, dTmin, n_remesh, remesh_num, &
               SProp%Patch, SProp%Remesh_flag, SProp%corner_index, SProp%edge_length, &
               SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, SProp%Disp, &
               SProp%ridge_group_num, SProp%ridge_num, SProp%ridge, &
               SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad, SCase%Disp)
  ! Solve Eulerian Description Stage
  CALL FEM_GXL(Control%Dim, 4, iter_solid, dTmin, n_remesh, remesh_num, &
               SProp%Patch, SProp%Remesh_flag, SProp%corner_index, SProp%edge_length, &
               SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, SProp%Disp, &
               SProp%ridge_group_num, SProp%ridge_num, SProp%ridge, &
               SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad, SCase%Disp)
  !call check_bound_face(SProp%nBNode, SProp%xyz, SProp%nBFace, SProp%f2n)
  ! Remeshing Stage
  IF (n_remesh == 2) THEN
    patch_num = maxval(SProp%Patch)
	ALLOCATE(patch_remesh_flag(patch_num))
	patch_remesh_flag = SProp%Remesh_flag
	
    ! Solid Remeshing
    CALL FEM_GXL(Control%Dim, 5, iter_solid, dTmin, n_remesh, remesh_num, &
                 SProp%Patch, SProp%Remesh_flag, SProp%corner_index, SProp%edge_length, &
                 SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, SProp%Disp, &
                 SProp%ridge_group_num, SProp%ridge_num, SProp%ridge, &
                 SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad, SCase%Disp)
    ! Solid Surface Re-initialization
    TEMPFLAG = .FALSE.
    CALL SURFACE(Control%Dim, TYPE1 = 4, PATCHNUM = PATCHNUM)
    DO I=1,PATCHNUM
        IF(patch_remesh_flag(I)==1) THEN
            TEMPFLAG = .TRUE.
        END IF
    END DO
    WRITE (*,*) 'patch_remesh_flag:', patch_remesh_flag
    !IF(TEMPFLAG) THEN         CALL LOADING_COEFFS_SURFACE_STRUCT(1, P_POINT_NUM, P_POINT, 2*P_FACE_NUM, TEMP_FACE, P_CORNER_INDEX, P_RIDGE_GROUP_NUM, P_RIDGE_NUM, P_RIDGE)
      CALL SURFACE(Control%Dim, TYPE1 = 2, TYPE2 = 0, &
                   P_POINT_NUM = SProp%nBNode, P_POINT = SProp%xyz, P_EDGE_NUM = SProp%nBFace, P_EDGE = SProp%f2n, P_LOC = SProp%Patch, P_CORNER_INDEX = SProp%Corner_index, &
		   P_RIDGE_GROUP_NUM = SProp%ridge_group_num, P_RIDGE_NUM = SProp%ridge_num,P_RIDGE = SProp%ridge, &
                   C_POINT_NUM = SCase%nBNode, C_POINT = SCase%xyz, C_EDGE_NUM = SCase%nBFace, C_EDGE = SCase%f2n)
    !ELSE
    !  CALL SURFACE(Control%Dim, TYPE1 = 5, P_POINT_NUM = SProp%nBNode, P_POINT = SProp%xyz, P_EDGE_NUM = SProp%nBFace, P_EDGE = SProp%f2n, P_LOC = SProp%Patch)
    !END IF
    DEALLOCATE(patch_remesh_flag)
  END IF

  CALL FEM_GXL(Control%Dim, 100, iter_solid, dTmin, n_remesh, remesh_num, &
               SProp%Patch, SProp%Remesh_flag, SProp%corner_index, SProp%edge_length, &
               SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, SProp%Disp, &
               SProp%ridge_group_num, SProp%ridge_num, SProp%ridge, &
               SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad, SCase%Disp)
			   
  ! Fluid -> Solid
  !Fluid%OutVars(1,:) = Fluid%OutVars(1,:) - 101325d0
  if (iter <= 200) then
	Fluid%OutVars(1,:) = (Time+dTmin)*8000.0*6894.73326
  else
    Fluid%OutVars(1,:) = 1600.0*6894.73326
  endif
  write (*,*) 'pressure:', Fluid%OutVars(1,1)/6894.73326 
  CALL SURFACE(Control%Dim, TYPE1 = 3, F_PRESSURE = Fluid%OutVars(1,:), P_FORCE = SProp%pLoad, P_DISPLACEMENT = SProp%Disp, C_FORCE = SCase%pLoad, C_DISPLACEMENT = SCase%Disp)
  CALL SURFACE(Control%Dim, TYPE1 = 103, FILENUM = iter)

  !press = (Time+dTmin)*80000.0*6894.73326
  !write (*,*) 'pressure:', press/6894.73326 
  !press = (Time+dTmin)*2000000.0
  !call Press2Load(Control%Dim, press, SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, &
  !                                    SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad)  
									  
!  fn = "./output/solid/load/load0000.plt"
!    write ( fn(25:28), '(I4.4)' ) iter
!    open (Unit=20, File=fn, STATUS='replace', ACTION='write')
!    Write (20,'(A)') 'TITLE="Check load"'
!    if (Control%Dim == 2) then
!        Write (20,*) 'VARIABLES="x", "y", "fx", "fy"'
!        Write (20,'(A,I5,A,I5, A)') 'ZONE T = "Propellant", N =', SProp%nBNode , ', E =', SProp%nBFace, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
!    elseif (Control%Dim == 3) then
!        Write (20,*) 'VARIABLES="x", "y", "z", "fx", "fy", "fz"'
!        Write (20,'(A,I5,A,I5, A)') 'ZONE T = "Propellant", N =', SProp%nBNode , ', E =', SProp%nBFace, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
!    endif
!    do i = 1, SProp%nBNode
!        write (20,*) SProp%xyz(:,i), SProp%pLoad(:,i)
!    enddo
!    do i = 1, SProp%nBFace
!        write (20,*) SProp%f2n(:,i)
!    enddo
!	if (Control%Dim == 2) then
!		Write (20,'(A,I5,A,I5, A)') 'ZONE T = "Case", N =', SCase%nBNode , ', E =', SCase%nBFace, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
!	elseif (Control%Dim == 3) then
!        Write (20,'(A,I5,A,I5, A)') 'ZONE T = "Case", N =', SCase%nBNode , ', E =', SCase%nBFace, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
!    endif
!    do i = 1, SCase%nBNode
!        write (20,*) SCase%xyz(:,i), SCase%pLoad(:,i)
!    enddo
!    do i = 1, SCase%nBFace
!        write (20,*) SCase%f2n(:,i)
!    enddo
!    close(20)
  ! Solid Lagrangian Description Stage
  CALL FEM_GXL(Control%Dim, 2, iter_solid, dTmin, n_remesh, remesh_num, &
               SProp%Patch, SProp%Remesh_flag, SProp%corner_index, SProp%edge_length, &
               SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, SProp%pLoad, SProp%Disp, &
               SProp%ridge_group_num, SProp%ridge_num, SProp%ridge, &
               SCase%nBNode, SCase%nBFace, SCase%xyz, SCase%f2n, SCase%pLoad, SCase%Disp)
			   
  !call surface_disp_interpolation(Control%Dim, maxval(SProp%Patch), iter, Fluid%nBNode, Fluid%xyz)

  !b_rate = 0.00004
  !write (*,*) 'SRMS >> move surface(burning_rate: ', b_rate, ')'
  !if (Control%Dim == 2) then
  !  CALL surface_move(SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, b_rate)
  !elseif (Control%Dim == 3) then
    !CALL surface_move_3D_temp(SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, b_rate)
  !  CALL surface_move_3D(iter_solid, SProp%nBNode, SProp%nBFace, SProp%xyz, SProp%f2n, b_rate)
  !endif
  
  ! Surface Moving with Deformation
  CALL SURFACE(Control%Dim, TYPE1 = 1, TYPE2 = 0, &
               F_POINT_NUM = Fluid%nBNode, F_POINT = Fluid%xyz, F_EDGE_NUM = Fluid%nBFace, F_EDGE = Fluid%f2n, F_LOC = Fluid%Patch,&
               P_POINT_NUM = SProp%nBNode, P_POINT = SProp%xyz, P_EDGE_NUM = SProp%nBFace, P_EDGE = SProp%f2n, P_DISPLACEMENT = SProp%Disp, P_LOC = SProp%Patch, &
               C_POINT_NUM = SCase%nBNode, C_POINT = SCase%xyz, C_EDGE_NUM = SCase%nBFace, C_EDGE = SCase%f2n, C_DISPLACEMENT = SCase%Disp)
  CALL SURFACE(Control%Dim, TYPE1 = 101, FILENUM = iter )
  CALL SURFACE(Control%Dim, TYPE1 = 104, FILENUM = iter )

  ! Solve Fluid Region
  CALL FLUS(Control%Dim,2,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%Patch,Fluid%BcFlag,Fluid%InVars,Fluid%OutVars,Fluid%Ref_Pressure)

  ! Fluid Surface Remeshing
  CALL SURFACE(Control%Dim, TYPE1 = 0, TYPE2 = 1, F_POINT_NUM = Fluid%nBNode, F_POINT = Fluid%xyz, F_EDGE_NUM = Fluid%nBFace, F_EDGE = Fluid%f2n, F_LOC = Fluid%Patch, FLAG = Fluid%n_remesh)
  IF( Fluid%n_remesh == .TRUE. ) THEN
    ! Fluid Remesh
    CALL FLUS(Control%Dim,3,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%Patch,Fluid%BcFlag,Fluid%InVars,Fluid%OutVars,Fluid%Ref_Pressure)
    ! Burn Remesh
    CALL SPBS(Control%Dim,3,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%BcFlag,Fluid%OutVars,Fluid%InVars,Fluid%Prop_Dens,Fluid%Ref_Pressure)
  END IF

  ! Solver Burn Region
  CALL SPBS(Control%Dim,2,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%BcFlag,Fluid%OutVars,Fluid%InVars,Fluid%Prop_Dens,Fluid%Ref_Pressure)

  ! Surface Moving with Burning Rate
  !Fluid%InVars(1,:) = Fluid%InVars(1,:) * regression_rate_x
  Fluid%InVars(1,:) = 0.005
  CALL SURFACE(Control%Dim, TYPE1 = 0, TYPE2 = 0, TIMESTEP = dTmin, &
               F_POINT_NUM = Fluid%nBNode, F_POINT = Fluid%xyz, F_EDGE_NUM = Fluid%nBFace, F_EDGE = Fluid%f2n, F_LOC = Fluid%Patch, F_B_RATE = Fluid%InVars(1,:), F_BCFLAG = Fluid%BcFlag, &
               P_POINT_NUM = SProp%nBNode, P_POINT = SProp%xyz, P_LOC = SProp%Patch, C_POINT = SCase%xyz)
  CALL SURFACE(Control%Dim, TYPE1 = 102, FILENUM = iter)

  ! Calculate Mass Flow Rate
  Fluid%InVars(1,:) = Fluid%InVars(1,:) * Fluid%Prop_Dens

  ! Time Update
  Time = Time + dTmin

  ! Time Check
  !IF( Time >= Control%TargetTime-1d-15 ) EXIT
  IF (iter >= Control%IterMax) EXIT

  ! Post Fluid / Burn Region
  IF( MOD(iter,Control%OutputIteration) == 0 ) THEN
    CALL POST
  END IF

END DO


! Post Final
CALL POST

! Finalize Fluid Region
CALL FLUS(Control%Dim,5,Time,Time,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%Patch,Fluid%BcFlag,Fluid%InVars,Fluid%OutVars,Fluid%Ref_Pressure)
! Finalize Burn Region
CALL SPBS(Control%Dim,5,Time,dTmin,Fluid%nBNode,Fluid%nBFace,Fluid%xyz,Fluid%f2n,Fluid%BcFlag,Fluid%OutVars,Fluid%InVars,Fluid%Prop_Dens,Fluid%Ref_Pressure)
! Finalize Surface
CALL SURFACE(Control%Dim, TYPE1 = -2)

! Print Finish Status
IF( iter > Control%IterMax ) iter = Control%IterMax
WRITE(*,*)
WRITE(*,*) 'SRMS >> ======================================================'
WRITE(*,*) 'SRMS >> Iteration          =', iter
WRITE(*,*) 'SRMS >> System Time        =', Time
WRITE(*,*) 'SRMS >> ======================================================'
WRITE(*,*)

END PROGRAM
