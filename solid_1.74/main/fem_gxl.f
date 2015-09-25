! Copyright(2006~2013) LASSCOM: Large Structural System Computing Lab, DGU, Seoul

subroutine FEM_GXL(Dim, solid_flag, current_step, delta_t, n_remesh, remesh, &
                   PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                   PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                   PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                   CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)

use file_info_house
use system_info_house
use system_domain

implicit none

! In or Out Variables
integer, intent(in   ) :: Dim, solid_flag
real,    intent(in   ) :: delta_t
integer, intent(inout) :: current_step
integer, intent(inout) :: n_remesh, remesh

INTEGER, INTENT(inout) :: PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_RIDGE_GROUP_NUM
INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_PATCH(:)
INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_CORNER(:)
INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_RIDGE_NUM(:)
INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_RIDGE(:,:)
REAL(8), INTENT(inout), ALLOCATABLE :: PROPEL_EDGE_LENGTH(:)
INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_EDGE(:,:)
INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_REMESH_FLAG(:)
REAL(8), INTENT(inout), ALLOCATABLE :: PROPEL_POINT(:,:)
REAL(8), INTENT(inout), ALLOCATABLE :: PROPEL_DISP(:,:)
REAL(8), INTENT(inout), ALLOCATABLE :: PROPEL_PLOAD(:,:)
        
INTEGER, INTENT(inout) :: CASE_POINT_NUM, CASE_EDGE_NUM
INTEGER, INTENT(inout), ALLOCATABLE :: CASE_EDGE(:,:)
REAL,    INTENT(inout), ALLOCATABLE :: CASE_POINT(:,:)
REAL,    INTENT(inout), ALLOCATABLE :: CASE_DISP(:,:)
REAL,    INTENT(inout), ALLOCATABLE :: CASE_PLOAD(:,:)
    
integer, parameter :: min_step_order = -3, max_step_order = 1
integer :: converge_status
real :: delta_t0, cur_time, temp_delt, pre_time, error
logical :: last_step
real, allocatable :: prop_inc_pload(:,:), prop_pre_pload(:,:)
real, allocatable :: case_inc_pload(:,:), case_pre_pload(:,:)

INTERFACE
  SUBROUTINE fsi_interface(Dim, action, current_step, delta_t, n_remesh, remesh, &
                           PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                           PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                           PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                           CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
    INTEGER, INTENT(in   ) :: Dim, action
    REAL(8), INTENT(in   ) :: delta_t
    INTEGER, INTENT(inout) :: current_step, n_remesh, remesh

    INTEGER, INTENT(inout) :: PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_RIDGE_GROUP_NUM
    INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_PATCH(:)
    INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_CORNER(:)
    INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_RIDGE_NUM(:)
    INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_RIDGE(:,:)
    REAL(8), INTENT(inout), ALLOCATABLE :: PROPEL_EDGE_LENGTH(:)
    INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_EDGE(:,:)
    INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_REMESH_FLAG(:)
    REAL(8), INTENT(inout), ALLOCATABLE :: PROPEL_POINT(:,:)
    REAL(8), INTENT(inout), ALLOCATABLE :: PROPEL_DISP(:,:)
    REAL(8), INTENT(inout), ALLOCATABLE :: PROPEL_PLOAD(:,:)
        
    INTEGER, INTENT(inout) :: CASE_POINT_NUM, CASE_EDGE_NUM
    INTEGER, INTENT(inout), ALLOCATABLE :: CASE_EDGE(:,:)
    REAL(8), INTENT(inout), ALLOCATABLE :: CASE_POINT(:,:)
    REAL(8), INTENT(inout), ALLOCATABLE :: CASE_DISP(:,:)
    REAL(8), INTENT(inout), ALLOCATABLE :: CASE_PLOAD(:,:)
  END SUBROUTINE
END INTERFACE

! Solid_flag = 0 : Open major files
!            = 1 : Initialize & build
!                  Initial stage (problem type: static, quasi-static, dynamic, nonlinearity)
!			       Preparing stage (solver type, sparse, symmetricity options)
!                  Convey data(Solid -> Surface)
!            = 11: Convey underomed coordinate(Solid -> Surface)
!            = 2 : Running solver stage (step(s))
!            = 3 : Movement stage
!            = 4 : Eulerian description stage
!            = 5 : equal 1 with remesh data
!            = 6 : close files

select case(solid_flag)
case(0)
    n_remesh = 0
    call write_input_filename(Dim)
    call fem_gxl_initialize(Dim, 0, solid_flag, remesh)
case(1)
    call fem_gxl_initialize(Dim, 1, solid_flag, remesh)
    call fem_gxl_foro(Dim, 0, n_remesh, remesh, current_step, delta_t, 1)
    call fsi_interface(Dim, -1, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(2) 
    ! Lagrangian description stage
    !if (Dim == 2) then 
        ! initialze 'timestep_order' & 'delt' of time_domain module
        !call fsi_interface(Dim, 0, current_step, delta_t, n_remesh, remesh, &
        !                    PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
        !                    PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
        !                    CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
        allocate (prop_inc_pload(Dim,PROPEL_POINT_NUM), prop_pre_pload(Dim,PROPEL_POINT_NUM))
        if (CASE_POINT_NUM /= 0) allocate (case_inc_pload(Dim,CASE_POINT_NUM), case_pre_pload(Dim,CASE_POINT_NUM))
        
        if (n_remesh == 0) then
            call difference_pload(Dim, 1, PROPEL_POINT_NUM, PROPEL_PLOAD, prop_inc_pload)
            prop_pre_pload = PROPEL_PLOAD 
            if (CASE_POINT_NUM /= 0) then
                call difference_pload(Dim, 2, CASE_POINT_NUM, CASE_PLOAD, case_inc_pload)
                case_pre_pload = CASE_PLOAD
            endif
        else
            prop_inc_pload = 0.0
            prop_pre_pload = PROPEL_PLOAD 
            if (CASE_POINT_NUM /= 0) then
                case_inc_pload = 0.0
                case_pre_pload = CASE_PLOAD
            endif
        endif
        current_step = current_step + 1
        if (delt == 0.0) then
            delta_t0 = delta_t
        else
            delta_t0 = delt
        endif
        cur_time = 0.0
        converge_status = 0
        call timestep_control(0, converge_status, min_step_order, max_step_order, delta_t0)
        write (*,*) 'init_step:', current_step
        write(*,'(/A,I6,5X,A)') ' >> FEM_GXL: computation start for time step #: ', current_step, solver_kind
        if (auto_time_stepping) write (*,'(A)') '  cur_time,  converge_status,  pre_time,   delt,      error_norm,   system_time'
        time_step_loop: do   ! Loop for time steps
            if (auto_time_stepping) then
                last_step = .FALSE.
                pre_time = cur_time
                if (delt > delta_t) call timestep_control(2, 1, min_step_order, max_step_order, delta_t)
                
                if ( delta_t*0.999 <= cur_time + delt) then
                    temp_delt = delt
                    call timestep_control(2, 1, min_step_order, max_step_order, delta_t-cur_time)
                    cur_time = delta_t
                    last_step = .TRUE.
                else
                    cur_time = cur_time + delt
                endif
                PROPEL_PLOAD = (prop_pre_pload-prop_inc_pload) + prop_inc_pload*cur_time/delta_t
                if (CASE_POINT_NUM /= 0) CASE_PLOAD = (case_pre_pload-case_inc_pload) + case_inc_pload*cur_time/delta_t
            else
                PROPEL_PLOAD = prop_pre_pload
                if (CASE_POINT_NUM /= 0) CASE_PLOAD = case_pre_pload
            endif
                
            call fsi_interface(Dim, 0, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)

            call fem_gxl_prima(Dim, current_step, converge_status, error, last_step)
            
            if (auto_time_stepping) then
                write (*,'(ES12.5,5X,I3,5X,4(ES12.5,2X))') cur_time, converge_status, pre_time, delt, error, delta_t
                call timestep_control(1, converge_status, min_step_order, max_step_order, delta_t0)
                if (converge_status > 0) then
                    if (last_step) then
                        call timestep_control(2, 1, min_step_order, max_step_order, temp_delt)
                        exit time_step_loop
                    endif
                else
                    cur_time = pre_time
                endif
            else
                if (converge_status > 0) then
                    exit time_step_loop
                else
                    STOP 'fem_gxl_prima: Fail to converge!'   !-------------------------------------------->>>
                endif
            endif
        enddo time_step_loop
        write (*,*)
        PROPEL_PLOAD = prop_pre_pload
        deallocate (prop_inc_pload, prop_pre_pload)
        if (CASE_POINT_NUM /= 0) then
            CASE_PLOAD = case_pre_pload
            deallocate (case_inc_pload, case_pre_pload)
        endif
        call fsi_interface(Dim, -2, current_step, delta_t, n_remesh, remesh, &
                            PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                            PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                            PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                            CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
    !elseif (Dim == 3) then
    !    call fem_gxl_prima(Dim, current_step, converge_status)
    !    
    !    call fsi_interface(Dim, -2, current_step, delta_t, n_remesh, remesh, &
    !                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
    !                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
    !                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
    !endif
case(3)
    ! movement stage
    call fsi_interface(Dim, 1, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(4)
    ! Eulerian description stage
    call fsi_interface(Dim, 3, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(5)
    call fem_gxl_initialize(Dim, 1, solid_flag, remesh)
    call fem_gxl_foro(Dim, 0, n_remesh, remesh, current_step, delta_t, 2) 
    call fsi_interface(Dim, -1, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(11)
    call fsi_interface(Dim, -3, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(12)
    call fsi_interface(Dim, 12, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(13)
    call fsi_interface(Dim, 13, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(14)
    call fsi_interface(Dim, 14, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(100)
    call fsi_interface(Dim, 100, current_step, delta_t, n_remesh, remesh, &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
case(6) 
    call close_files
case default
    STOP 'fem_gxl : invalid solid_flag number!'
end select    

contains

subroutine write_input_filename(Dim)

implicit none

integer :: Dim

if (Dim == 2) then
    filename = './input/solid.inp'
elseif (Dim == 3) then
    filename = './input/solid_3D.inp'
endif

end subroutine write_input_filename

subroutine close_files

implicit none 

close(unit_number)
close(out1_unit_number)
close(out2_unit_number)
close(unit_num_node) 
close(unit_num_elem) 
close(unit_num_cont)
close(unit_num_inte)
if ( remesh /= 0 ) close(unit_num_remesh)

end subroutine close_files

end subroutine fem_gxl




