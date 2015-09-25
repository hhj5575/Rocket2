subroutine fsi_interface(Dim, action, current_step, delta_t, n_remesh, remesh,  &
                        PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)

! action = -1: execute 'get_fsi_info' only
!        =  0: getting data(ni, ni_flag, inte_coord) only
!        =  1: all but 'action = -1'

use file_info_house
use system_info_house
use fsi_domain

implicit none

! In or Out Variables
integer, intent(in   ) :: Dim, action
real,    intent(in   ) :: delta_t
integer, intent(inout) :: current_step
integer, intent(inout) :: n_remesh, remesh

INTEGER, INTENT(in):: PROPEL_RIDGE_GROUP_NUM
INTEGER, INTENT(inout):: PROPEL_POINT_NUM, PROPEL_EDGE_NUM
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

integer :: num_rn, ng, cont_node(PROPEL_POINT_NUM), pre_n_remesh

!INTERFACE
!  SUBROUTINE alloc_fsi_info(PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_EDGE, PROPEL_POINT, PROPEL_DISP, PROPEL_PLOAD, &
!                            CASE_POINT_NUM, CASE_EDGE_NUM, CASE_EDGE, CASE_POINT, CASE_DISP, CASE_PLOAD)
!
!    INTEGER, INTENT(inout) :: PROPEL_POINT_NUM, PROPEL_EDGE_NUM
!    INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_EDGE(:,:)
!    REAL,    INTENT(inout), ALLOCATABLE :: PROPEL_POINT(:,:)
!    REAL,    INTENT(inout), ALLOCATABLE :: PROPEL_DISP(:,:)
!    REAL,    INTENT(inout), ALLOCATABLE :: PROPEL_PLOAD(:)
!        
!    INTEGER, INTENT(inout) :: CASE_POINT_NUM, CASE_EDGE_NUM
!    INTEGER, INTENT(inout), ALLOCATABLE :: CASE_EDGE(:,:)
!    REAL,    INTENT(inout), ALLOCATABLE :: CASE_POINT(:,:)
!    REAL,    INTENT(inout), ALLOCATABLE :: CASE_DISP(:,:)
!    REAL,    INTENT(inout), ALLOCATABLE :: CASE_PLOAD(:)
!  END SUBROUTINE
!END INTERFACE

if (Dim == 2) then
    ng = region(1)%num_groups ! physical domain module
    num_rn = num_of_regions ! global variable from physical_domain module
elseif (Dim == 3) then
    ng = sub_mesh_region
    num_rn = divided_mesh_region
endif

if (action <= -2) then
    call convey_fsi_info(Dim, action, n_remesh, PROPEL_RIDGE_GROUP_NUM, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_CORNER, &
                         CASE_POINT_NUM, CASE_POINT, CASE_DISP)
    if (action == -2) then
        write (*,*) '>> Update pre pload'  
        call fsi_data_update(Dim, 1, PROPEL_POINT_NUM, PROPEL_PLOAD)
        if (num_rn == 2) call fsi_data_update(Dim, 2, CASE_POINT_NUM, CASE_PLOAD)
    endif
elseif (action == -1) then
    call alloc_fsi_info(Dim, PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                        PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                        CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)
    if (n_remesh == 0) then
        call convey_fsi_info(Dim, action, n_remesh, PROPEL_RIDGE_GROUP_NUM, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_CORNER, &
                             CASE_POINT_NUM, CASE_POINT, CASE_DISP)
    else
        call convey_fsi_info(Dim, action, n_remesh, PROPEL_RIDGE_GROUP_NUM, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_CORNER, &
                             CASE_POINT_NUM, CASE_POINT, CASE_DISP, &
                             PROPEL_RIDGE_NUM = PROPEL_RIDGE_NUM(1:PROPEL_RIDGE_GROUP_NUM), PROPEL_RIDGE = PROPEL_RIDGE(1:PROPEL_RIDGE_GROUP_NUM,1:1000))
    endif
    !write (*,*) ' >> Adjust deformed coordinate'
elseif (action == 0) then
    if (solver_flag) then
        call write_fsi_info(Dim, current_step, PROPEL_POINT_NUM, PROPEL_PLOAD, CASE_POINT_NUM, CASE_PLOAD)
    endif
    !call write_fsi_info_bak(Dim, current_step, PROPEL_POINT_NUM, PROPEL_PLOAD, CASE_POINT_NUM, CASE_PLOAD)
elseif (action == 1) then
    call surface_control(Dim, num_rn, ng, current_step, n_remesh, &
                        PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_PATCH, PROPEL_REMESH_FLAG, CASE_POINT_NUM, CASE_POINT, CASE_DISP)
elseif (action == 3) then
    if (n_remesh >= 3) then
        write (*,'(A)') " >> Start rearrange_edge"
        call rearrange_edge(Dim, current_step, ng, PROPEL_REMESH_FLAG, PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_EDGE, &
                            PROPEL_PATCH, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                            PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM(1:PROPEL_RIDGE_GROUP_NUM), PROPEL_RIDGE(1:PROPEL_RIDGE_GROUP_NUM,1:1000))
		n_remesh = 2
    endif
    call ale_control(Dim, num_rn, delta_t, current_step, remesh, n_remesh, ng, &
                     PROPEL_REMESH_FLAG, PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_EDGE, PROPEL_PATCH)
    if (n_remesh == 0) then
        call convey_fsi_info(Dim, -2, n_remesh, PROPEL_RIDGE_GROUP_NUM, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_CORNER, &
                             CASE_POINT_NUM, CASE_POINT, CASE_DISP)
    endif
elseif (action == 12) then
    call convey_pair_fsi_info(Dim, 1, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_PLOAD, PROPEL_DISP)
elseif (action == 13) then
    call convey_pair_fsi_info(Dim, 2, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_PLOAD, PROPEL_DISP)
elseif (action == 14) then
    call convey_pair_fsi_info(Dim, 3, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_PLOAD, PROPEL_DISP)
elseif (action == 100) then
    !if (solver_flag) then
        if (current_step /= 0 .and. mod(current_step, abs(output_interval2)) == 0) then
            call make_backup_file(Dim, current_step, remesh)
        endif
    !endif
endif

end subroutine fsi_interface
