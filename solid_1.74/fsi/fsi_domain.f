module fsi_domain 

use physical_domain
use time_domain

implicit none
save

type fsi_inte_type
    integer :: fsi_ni = 0, old_fsi_ni = 0, fsi_load_type, fsi_inte_point(4)
    real :: delta_t, tol, max_b_rate, fsi_inte_coord(2,2)
    integer, allocatable :: fsi_inte_node(:,:), old_fsi_inte(:,:), B_flag(:), ni_flag(:)
    real, allocatable :: B_rate(:), fsi_coord(:,:)
    logical, allocatable :: new_fsi_node(:)
    
    integer :: prop_point_num, case_point_num
    integer, allocatable :: prop_adjn(:), prop_pair(:,:)
    real, allocatable :: prop_point(:,:), case_point(:,:)
    real, allocatable :: prop_disp(:,:), case_disp(:,:)
    real, allocatable :: prop_pre_pload(:,:), case_pre_pload(:,:)
end type

type(fsi_inte_type) :: fsi_inte_set
logical :: fsi_node_change = .FALSE.

! variables about divided mesh for 3D
integer :: divided_mesh_region = 1,  sub_mesh_region = 1, remesh_type_number = 1

contains

!===========================================================================
subroutine write_fsi_info_bak(Dim, current_step, PROPEL_POINT_NUM, PROPEL_PLOAD, CASE_POINT_NUM, CASE_PLOAD)

implicit none

integer, intent(in) :: Dim, current_step, PROPEL_POINT_NUM, CASE_POINT_NUM
real, intent(inout) :: PROPEL_PLOAD(PROPEL_POINT_NUM), CASE_PLOAD(CASE_POINT_NUM)

integer :: i, j, k, adjn, num, seq, cp, rn, flag, nBNode, tol_rn
integer :: case_adjn, ngroup, max_adjn, inte_rn(2)
real :: dis, temp_node(2,2), temp_load(2), vec(2), temp_vec(2), temp_val
integer, allocatable :: case_adj(:), prop_sp(:), prop_adjn(:), prop_adj(:,:)
integer, allocatable :: inte_load_num(:,:), load_point(:,:)
real, allocatable :: inte_load(:), load_value(:), pload(:)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375
logical :: check
!-------------------------------------------------------------------------------
tol_rn = num_of_regions
nBNode = fsi_inte_set%fsi_ni

allocate (inte_load_num(3,nBNode*2), inte_load(nBNode*2), pload(nBNode))
inte_load_num = 0
inte_load = 0.0

ngroup = region(1)%num_groups
allocate (prop_adjn(ngroup), prop_sp(ngroup))
do i = 1, ngroup
    call set_skin(1, i, prop_sp(i), prop_adjn(i))
enddo
max_adjn = maxval(prop_adjn(:))
allocate (prop_adj(ngroup, max_adjn))
do i = 1, ngroup
    prop_adj(i,1:prop_adjn(i)) = region(1)%skin_nodes(prop_sp(i):prop_sp(i)+prop_adjn(i)-1)
enddo

if (tol_rn >= 2) then
    call set_skin(2, 1, num, case_adjn)
    allocate (case_adj(case_adjn))
    case_adj = region(2)%skin_nodes(num:num+case_adjn-1)
endif

pload = 0.0
do i = 1, nBNode
    rn = fsi_inte_set%fsi_inte_node(1,i)
    cp = fsi_inte_set%fsi_inte_node(2,i)
    !write (*,*) 'RN, CP: ', rn, cp
    if (tol_rn < rn) then
        !write (*,*) 'Error(write_fsi_info), Total number of region < Interface region'
        !write (*,*) 'i:', i, '. rn, cp:', fsi_inte_set%fsi_inte_node(1:2,i)
        !stop
        fsi_inte_set%fsi_inte_node(1,i) = 2
        rn = 2
    endif
    if (rn == 1) then
        num = 0
        seq = 0
        do j = 1, ngroup
            do k = 1, prop_adjn(j)
                if ( prop_adj(j,k) == cp ) then
                    seq = k
                    exit
                endif
            enddo
            if (seq /= 0) then
                pload(i) = PROPEL_PLOAD(prop_adjn(j)-seq+1+num)
                exit
            endif
            num = num + prop_adjn(j)
        enddo
        if (seq == 0) then
            write (*,*) 'RN, CP: ', rn, cp
            STOP 'Error : find_adj_seq: cannot find sequence of node(cp) in adj!!!'
        endif
    elseif (rn == 2) then
        call find_adj_seq(cp, case_adjn, case_adj, seq)
        pload(i) = CASE_PLOAD(CASE_POINT_NUM-seq+1)
    endif
enddo

do i = 1, nBNode
    temp_load = 0.0
    flag = 0
    do j = 1, 2
        check = .TRUE.
        if (j == 1 .AND. i /= 1) then
            inte_rn(2) = fsi_inte_set%fsi_inte_node(1,i-1)
            cp = fsi_inte_set%fsi_inte_node(2,i-1)
            temp_node(:,1) = region(inte_rn(2))%node_coord(:,cp) + history(inte_rn(2))%data(region(inte_rn(2))%node(4:5,cp),1)
            inte_rn(1) = fsi_inte_set%fsi_inte_node(1,i)
            cp = fsi_inte_set%fsi_inte_node(2,i)
            temp_node(:,2) = region(inte_rn(1))%node_coord(:,cp) + history(inte_rn(1))%data(region(inte_rn(1))%node(4:5,cp),1)
            temp_val = (pload(i)+(pload(i)+pload(i-1))*0.5)*0.5
        elseif (j == 2 .AND. i /= nBNode) then
            inte_rn(1) = fsi_inte_set%fsi_inte_node(1,i)
            cp = fsi_inte_set%fsi_inte_node(2,i)
            temp_node(:,1) = region(inte_rn(1))%node_coord(:,cp) + history(inte_rn(1))%data(region(inte_rn(1))%node(4:5,cp),1)
            inte_rn(2) = fsi_inte_set%fsi_inte_node(1,i+1)
            cp = fsi_inte_set%fsi_inte_node(2,i+1)
            temp_node(:,2) = region(inte_rn(2))%node_coord(:,cp) + history(inte_rn(2))%data(region(inte_rn(2))%node(4:5,cp),1)
            temp_val = (pload(i)+(pload(i)+pload(i+1))*0.5)*0.5
        else
            check = .FALSE.
        endif
        temp_val = PROPEL_PLOAD(1)
        !if (inte_rn(1) == 1 .and. inte_rn(2) == 2) then
        !    !write (*,*) 'rn:', inte_rn, 'cp:', fsi_inte_set%fsi_inte_node(2,i)
        !    check = .FALSE.
        !    flag = 1
        !endif
        
        if (check) then
            call calc_len(temp_node(:,1), temp_node(:,2), dis)
            vec(:) = ((temp_node(:,2)-temp_node(:,1))/dis)
            temp_vec(1) = -vec(2)
            temp_vec(2) = vec(1)
            temp_load = temp_load + temp_vec*dis*0.5*temp_val
        endif
    enddo
    if (flag == 1) then
        temp_load = temp_load * 2.0
    endif
    !write (*,*) fsi_inte_set%fsi_inte_node(1:2,i), temp_load
    do  j = 1, 2
        inte_load_num(1,(i-1)*2+j) = fsi_inte_set%fsi_inte_node(1,i)  ! rn number
        inte_load_num(2,(i-1)*2+j) = fsi_inte_set%fsi_inte_node(2,i)  ! node number
        inte_load_num(3,(i-1)*2+j) = j   ! dof number
        inte_load((i-1)*2+j) = temp_load(j)
    enddo
enddo

do rn = 1, tol_rn
    if (rn == 1) then
        adjn = sum(prop_adjn(:))
    else
        adjn = case_adjn
    endif
    
    allocate (load_point(2,adjn*2), load_value(adjn*2))
    num = 0
    do i = 1, nBNode*2
        if (inte_load_num(1,i) == rn) then
            num = num + 1
            load_point(1:2,num) = inte_load_num(2:3,i)
            load_value(num) = inte_load(i)
            ! test
            !load_value(num) = load_value(num) * 1.0/sqrt(sqrt(1.0-region(rn)%node_coord(2,load_point(1,num))))
        endif
    enddo
    call update_fsi_load(rn, Dim, num, load_point(:,1:num), load_value(1:num))
    deallocate (load_point, load_value)
enddo

! B_flag
! case = 0, propellent = 1
!do i = 1, nBNode
!    fsi_inte_set%B_flag(i) = b_flag(i)
!enddo 


deallocate (inte_load_num, inte_load)
deallocate (prop_adjn, prop_sp, prop_adj)
if ( tol_rn >= 2 ) deallocate(case_adj)

end subroutine write_fsi_info_bak

!===========================================================================
subroutine write_fsi_info(Dim, current_step, PROPEL_POINT_NUM, PROPEL_PLOAD, CASE_POINT_NUM, CASE_PLOAD)

implicit none

integer, intent(in) :: Dim, current_step, PROPEL_POINT_NUM, CASE_POINT_NUM
real, intent(in) :: PROPEL_PLOAD(Dim,PROPEL_POINT_NUM), CASE_PLOAD(Dim,CASE_POINT_NUM)

integer :: tol_rn, rn, i, j, nn, ne, status, num, adjn, bn, bfn
integer, allocatable :: adj(:), load_point(:,:), elem(:,:), bound_node(:), bound_face(:,:)
real, allocatable :: load_value(:), temp_pload(:,:), node(:,:), load(:,:)
character(len=40) :: fn

!-------------------------------------------------------------------------------
!fn = "./output/solid/load/tec_load0000000.plt"
!write (fn(29:35), '(I7.7)' ) current_step
!open (Unit=37, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)
!Write (37,'(A)') 'TITLE="Check Load"'
!if (Dim == 2) then
!    Write (37,'(A,A)') 'VARIABLES= "x", "y", "fx", "fy"'
!elseif (Dim == 3) then
!    Write (37,'(A,A)') 'VARIABLES= "x", "y", "z", "fx", "fy", "fz"'
!    Write (37,'(A,I5,A,I5, A)') 'ZONE T = "Propellant", N =', SProp%nBNode , ', E =', SProp%nBFace, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
!endif
!write (*,*) ' >> write_fsi_info'
if (Dim == 2) then
    tol_rn = num_of_regions
    do rn = 1, tol_rn
        nn = region(rn)%num_nodes
        ne = region(rn)%num_elements
        allocate(node(Dim,nn), load(Dim,nn), elem((Dim-1)*4,ne))
        node = region(rn)%node_coord
        elem = region(rn)%element(4:3+(Dim-1)*4,:)
        load = 0.0
        !Write (37,'(A,I1,A,I5,A,I5, A)') 'ZONE T = "Region', rn, '", N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
    
        adjn = region(rn)%num_skin
        allocate(adj(adjn), load_point(Dim,adjn*2), load_value(adjn*Dim))
        allocate(temp_pload(Dim,adjn))
        adj = region(rn)%skin_nodes
        temp_pload = 0.0
        num = 0
        do i = adjn, 1, -1
            num = num + 1
            if (rn == 1) then
                temp_pload(1,num) = PROPEL_PLOAD(2,i)
                temp_pload(2,num) = PROPEL_PLOAD(1,i)
            elseif (rn == 2) then
                temp_pload(1,num) = CASE_PLOAD(2,i)
                temp_pload(2,num) = CASE_PLOAD(1,i)
            endif
        enddo

        num = 0
        do i = 1, adjn
            do j = 1, 2
                num = num + 1
                load_point(1:2,num) = (/ adj(i), j /)
                load_value(num) = temp_pload(j,i)
                load(j,adj(i)) = temp_pload(j,i)
            enddo
            if (rn == 1) then
                fsi_inte_set%prop_point(:,i) = region(rn)%node_coord(:,adj(i)) + history(rn)%data(region(rn)%node(4:5,adj(i)),1)
            elseif (rn == 2) then
                fsi_inte_set%case_point(:,i) = region(rn)%node_coord(:,adj(i)) + history(rn)%data(region(rn)%node(4:5,adj(i)),1)
            endif
        enddo
    
        call update_fsi_load(rn, Dim, num, load_point(:,1:num), load_value(1:num))
        deallocate(adj, load_point, load_value)
        deallocate(temp_pload)
    
        !do i = 1, nn
        !    write (37,*) node(:,i), load(:,i)
        !enddo
        !do i = 1, ne
        !    write (37,*) elem(:,i)
        !enddo
        
        deallocate (node, load, elem)
    enddo
elseif (Dim == 3) then
    !tol_rn = num_of_regions
    nn = region(1)%num_nodes
    allocate(load_point(2,nn*Dim), load_value(nn*Dim), temp_pload(Dim,nn))
    num = 0
    load_value = 0
    load_point = 0
    do i = 1, nn
        do j = 1, Dim
            load_point(1:2,(i-1)*Dim+j) = (/ i, j /)
        enddo
    enddo
    do rn = 1, divided_mesh_region
        temp_pload = 0.0
        call get_skin_surface_info1(rn, bn, bfn)
        allocate (bound_node(bn), bound_face(4,bfn))
        call get_skin_surface_info2(rn, bn, bound_node, bfn, bound_face)
        
        if (rn == 1) then
            temp_pload(:,1:PROPEL_POINT_NUM) = PROPEL_PLOAD(:,1:PROPEL_POINT_NUM)
        elseif (rn == 2) then
            temp_pload(:,1:CASE_POINT_NUM) = CASE_PLOAD(:,1:CASE_POINT_NUM)
        endif

        do i = 1, bn
            do j = 1, Dim
                num = (bound_node(i)-1)*Dim+j
                !load_point(1:2,num) = (/ bound_node(i), j /)
                load_value(num) = load_value(num) + temp_pload(j,i)
                !write (*,*) load_point(1:2,num), load_value(num)
            enddo
        enddo
        deallocate (bound_node, bound_face)
    enddo
    
    fn = "./output/solid/load/pload000000.plt"
    write ( fn(26:31), '(I6.6)' ) current_step
    open (Unit=20, File=fn, STATUS='replace', ACTION='write')
    Write (20,'(A)') 'TITLE="Check load"'
    Write (20,*) 'VARIABLES="x", "y", "z", "fx", "fy", "fz"'
    Write (20,'(A,I5,A,I5, A)') 'ZONE T = "Propellant", N =', nn , ', E =', region(1)%num_elements, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
    do i = 1, nn
        write (20,*) region(1)%node_coord(:,i), load_value((i-1)*3+1:(i-1)*3+3)
    enddo
    do i = 1, region(1)%num_elements
        write (20,*) region(1)%element(4:11,i)
    enddo
    close(20)
    ! case one domain
    call update_fsi_load(1, Dim, nn*Dim, load_point, load_value)
        
    !deallocate (bound_node, bound_face, load_point, load_value, temp_pload)
    deallocate (load_point, load_value, temp_pload)
endif
!close(37)

end subroutine write_fsi_info
!===========================================================================

subroutine fsi_data_update(Dim, rn, nn, pload)

implicit none

integer, intent(in) :: Dim, rn, nn
real, intent(in) :: pload(Dim,nn)

if (rn == 1) then
    fsi_inte_set%prop_pre_pload = pload
else
    fsi_inte_set%case_pre_pload = pload
endif
    
end subroutine fsi_data_update
!===========================================================================
subroutine find_number_group(cp, ngroup, max_adjn, adjn, adj, ng, seq)
                
implicit none

integer, intent(in) :: cp, ngroup, max_adjn
integer, intent(in) :: adjn(ngroup), adj(ngroup,max_adjn)
integer, intent(inout) :: ng, seq

integer :: i, j
logical :: check

check = .FALSE.
do i = 1, ngroup
    do j = 1, adjn(i)
        if ( adj(i,j) == cp ) then
            ng = i
            seq = j
            check = .TRUE.
            exit
        endif
    enddo
    if ( check ) exit
enddo    
    
if ( check == .FALSE. ) stop 'Error : find_mesh_set_num(fsi_lib)!!' 

end subroutine find_number_group


subroutine find_adj_seq(cp, adjn, adj, seq)

implicit none

integer, intent(in) :: cp, adjn, adj(adjn)
integer, intent(inout) :: seq

integer :: i, j

seq = 0
do i = 1, adjn
    if ( adj(i) == cp ) then
        seq = i
        exit
    endif
enddo

j = 0
!if (seq == 0) STOP 'Error : find_adj_seq: cannot find sequence of node in adj!!!'
if (seq == 0) then
    write (*,*) 'cp:', cp
    write (*,*) adj(j)
endif

end subroutine find_adj_seq


subroutine convey_fsi_info(Dim, action, n_remesh, PRGN, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_CORNER, &
                           CASE_POINT_NUM, CASE_POINT, CASE_DISP, PROPEL_RIDGE_NUM, PROPEL_RIDGE)

use remesh_domain
implicit none

INTEGER, INTENT(in)    :: Dim, action, n_remesh, PRGN
INTEGER, INTENT(inout) :: PROPEL_POINT_NUM
INTEGER, INTENT(INOUT) :: PROPEL_CORNER(PROPEL_POINT_NUM)
REAL,    INTENT(inout) :: PROPEL_POINT(Dim,PROPEL_POINT_NUM)
REAL,    INTENT(inout) :: PROPEL_DISP(Dim,PROPEL_POINT_NUM)
INTEGER, INTENT(inout) :: CASE_POINT_NUM
REAL,    INTENT(inout) :: CASE_POINT(Dim,CASE_POINT_NUM)
REAL,    INTENT(inout) :: CASE_DISP(Dim,CASE_POINT_NUM)

INTEGER, OPTIONAL :: PROPEL_RIDGE_NUM(PRGN)
INTEGER, OPTIONAL :: PROPEL_RIDGE(PRGN,1000)

integer :: i, cp, num, cont_node(PROPEL_POINT_NUM), CASE_CORNER(CASE_POINT_NUM)
real :: spatial_coord(2)
real :: prop_deform(2,PROPEL_POINT_NUM), case_deform(2,CASE_POINT_NUM)
real :: pdisp(2,PROPEL_POINT_NUM), cdisp(2,CASE_POINT_NUM)

if (Dim == 2) then
    !if (action == -2) then
    !    if (num_of_regions == 2) then
    !        do i = 1, PROPEL_POINT_NUM
    !            cp = region(1)%skin_nodes(i)
    !            prop_deform(:,i) = region(1)%node_coord(:,cp)
    !            pdisp(:,i) = history(1)%data(region(1)%node(4:5,cp),1)
    !        enddo
    !        do i = 1, CASE_POINT_NUM
    !            cp = region(2)%skin_nodes(i)
    !            case_deform(:,i) = region(2)%node_coord(:,cp)
    !            cdisp(:,i) = history(2)%data(region(2)%node(4:5,cp),1)
    !        enddo
    !    else
    !        do i = 1, PROPEL_POINT_NUM
    !            cp = region(1)%skin_nodes(i)
    !            prop_deform(:,i) = region(1)%node_coord(:,cp)+history(1)%data(region(1)%node(4:5,cp),1)
    !        enddo
    !        cont_node = 0
    !    endif
    !endif
    num = 0
    do i = PROPEL_POINT_NUM, 1, -1
        cp = region(1)%skin_nodes(i)
        num = num + 1
        if (action == -1) then
            PROPEL_POINT(2,num) = region(1)%node_coord(1,cp)
            PROPEL_POINT(1,num) = region(1)%node_coord(2,cp)
            PROPEL_DISP(2,num) = history(1)%data(region(1)%node(4,cp),1)
            PROPEL_DISP(1,num) = history(1)%data(region(1)%node(5,cp),1)
        elseif (action == -2) then
            PROPEL_DISP(2,num) = history(1)%data(region(1)%node(4,cp),1)
            PROPEL_DISP(1,num) = history(1)%data(region(1)%node(5,cp),1)
        elseif (action == -3) then
            PROPEL_POINT(2,num) = region(1)%node_coord(1,cp)
            PROPEL_POINT(1,num) = region(1)%node_coord(2,cp)
        endif
    enddo

    if (num_of_regions == 2) then
        num = 0
        do i = CASE_POINT_NUM, 1, -1
            cp = region(2)%skin_nodes(i)
            num = num + 1
            if (action == -1) then
                CASE_POINT(2,num) = region(2)%node_coord(1,cp)
                CASE_POINT(1,num) = region(2)%node_coord(2,cp)
                CASE_DISP(2,num) = history(2)%data(region(2)%node(4,cp),1)
                CASE_DISP(1,num) = history(2)%data(region(2)%node(5,cp),1)
            elseif (action == -2) then
                CASE_DISP(2,num) = history(2)%data(region(2)%node(4,cp),1)
                CASE_DISP(1,num) = history(2)%data(region(2)%node(5,cp),1)
            elseif (action == -3) then
                CASE_POINT(2,num) = region(2)%node_coord(1,cp)
                CASE_POINT(1,num) = region(2)%node_coord(2,cp)
            endif
        enddo
    endif
elseif (Dim == 3) then
    if (action <= -1) then
        call get_skin_surface_info3(1, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_CORNER, action)
        if (divided_mesh_region == 2) then
            CASE_CORNER = 0
            call get_skin_surface_info3(2, CASE_POINT_NUM, CASE_POINT, CASE_DISP, CASE_CORNER, action)
        endif
    endif
    if (action == -1) then
        if (n_remesh /= 0) then
            !call revise_add_data_for_remesh(nn, bn_coord, mesh_set(rn)%bound_lv)
            call revise_b_edge_info(sub_mesh_region, PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_CORNER, &
                                    PRGN, PROPEL_RIDGE_NUM, PROPEL_RIDGE)
            call free_geo_domain(2)
        endif
        call create_remesh_domain(sub_mesh_region, 2)
    endif
endif

end subroutine convey_fsi_info



subroutine check_new_fsi_node

implicit none

integer :: i, j, rn, point

do i = 1, fsi_inte_set%fsi_ni
    rn = fsi_inte_set%fsi_inte_node(1,i)
    point = fsi_inte_set%fsi_inte_node(2,i)
    do j = 1, fsi_inte_set%old_fsi_ni
        if (rn == fsi_inte_set%old_fsi_inte(1,j) .AND. point == fsi_inte_set%old_fsi_inte(1,j)) then
            fsi_inte_set%new_fsi_node(i) = .TRUE.
            exit
        endif
    enddo
enddo

end subroutine check_new_fsi_node



!=============================================================================

subroutine find_inte_coord(nn, coord)

implicit none

integer, intent(in) :: nn
real, intent(in) :: coord(2,nn)

integer :: i

do i = 1, 2
    if (fsi_inte_set%fsi_inte_point(i*2-1) == 1) then
        fsi_inte_set%fsi_inte_coord(:,i) = coord(:,fsi_inte_set%fsi_inte_point(2*i))
    endif
enddo

end subroutine find_inte_coord



subroutine alloc_fsi_info(Dim, PROPEL_PATCH, PROPEL_REMESH_FLAG, PROPEL_CORNER, PROPEL_EDGE_LENGTH, &
                          PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_EDGE, PROPEL_PLOAD, PROPEL_DISP, &
                          CASE_POINT_NUM, CASE_EDGE_NUM, CASE_POINT, CASE_EDGE, CASE_PLOAD, CASE_DISP)

implicit none

INTEGER, INTENT(in   ) :: Dim
INTEGER, INTENT(inout) :: PROPEL_POINT_NUM, PROPEL_EDGE_NUM
INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_PATCH(:)
INTEGER, INTENT(inout), ALLOCATABLE :: PROPEL_CORNER(:)
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

integer :: i, j, num, sp, lp, ng, sn, nn, prop_num, f2n_num, bn, bfn, num_rn
logical :: tempflag
integer, allocatable :: conn(:), bound_node(:)
real, allocatable :: prop_coord(:,:)

tempflag = .false.
if (allocated(PROPEL_POINT)) then
    prop_num = PROPEL_POINT_NUM
    allocate (prop_coord(Dim,prop_num))
    prop_coord = PROPEL_POINT
    do i = 1, region(1)%num_groups
        if (PROPEL_REMESH_FLAG(i) == 2) THEN
            tempflag = .TRUE.
            exit
        endif
    enddo
    DEALLOCATE(PROPEL_EDGE, PROPEL_POINT, PROPEL_DISP, PROPEL_PLOAD)
    DEALLOCATE(PROPEL_PATCH, PROPEL_REMESH_FLAG)
    DEALLOCATE(PROPEL_CORNER, PROPEL_EDGE_LENGTH)
    deallocate(fsi_inte_set%prop_disp, fsi_inte_set%prop_adjn)
    deallocate(fsi_inte_set%prop_point, fsi_inte_set%prop_pair)
endif
if (allocated(CASE_POINT)) then
    DEALLOCATE(CASE_EDGE, CASE_POINT, CASE_DISP, CASE_PLOAD)
    deallocate (fsi_inte_set%case_disp)
    deallocate (fsi_inte_set%case_point)
endif

ng = region(1)%num_groups
f2n_num = (Dim-1)*2
if (Dim == 2) then
    num_rn = num_of_regions
    num = region(1)%num_skin
    PROPEL_POINT_NUM = num
    PROPEL_EDGE_NUM = num
    fsi_inte_set%prop_point_num = num
elseif (Dim == 3) then
    num_rn = divided_mesh_region
    ng = sub_mesh_region
    call get_skin_surface_info1(1, bn, bfn)
    PROPEL_POINT_NUM = bn
    PROPEL_EDGE_NUM = bfn
    fsi_inte_set%prop_point_num = bn
endif
allocate (PROPEL_EDGE(f2n_num,PROPEL_EDGE_NUM), PROPEL_POINT(Dim,PROPEL_POINT_NUM))
allocate (PROPEL_DISP(Dim,PROPEL_POINT_NUM), PROPEL_PLOAD(Dim,PROPEL_POINT_NUM))
allocate (PROPEL_PATCH(PROPEL_POINT_NUM), PROPEL_REMESH_FLAG(ng))
allocate (PROPEL_CORNER(PROPEL_POINT_NUM), PROPEL_EDGE_LENGTH(PROPEL_EDGE_NUM))
allocate (fsi_inte_set%prop_disp(Dim,PROPEL_POINT_NUM), fsi_inte_set%prop_adjn(ng))
allocate (fsi_inte_set%prop_point(Dim,PROPEL_POINT_NUM), fsi_inte_set%prop_pair(Dim,PROPEL_POINT_NUM))
allocate (fsi_inte_set%prop_pre_pload(Dim,PROPEL_POINT_NUM))
PROPEL_EDGE = 0
PROPEL_POINT = 0.0
PROPEL_DISP = 0.0
PROPEL_PLOAD = 0.0
PROPEL_PATCH = 0
PROPEL_REMESH_FLAG = 0
PROPEL_CORNER = 0
PROPEL_EDGE_LENGTH = 0.0
fsi_inte_set%prop_point = 0.0
fsi_inte_set%prop_disp = 0.0
fsi_inte_set%prop_adjn = 0
fsi_inte_set%prop_pair = 0
fsi_inte_set%prop_pre_pload = 0.0
    
fsi_inte_set%prop_adjn(1) = 1
if (Dim == 2) then
    num = 1
    do i = ng, 2, -1
        if (i == ng) then
            sn = region(1)%num_skin - region(1)%skin_group_pointer(i) + 1
        else
            sn = region(1)%skin_group_pointer(i+1) - region(1)%skin_group_pointer(i)
        endif
        fsi_inte_set%prop_adjn(num+1) = fsi_inte_set%prop_adjn(num) + sn
        num = num + 1
    enddo
    
    do i = 1, ng
        sp = fsi_inte_set%prop_adjn(i)
        if (i == ng) then
            lp = PROPEL_POINT_NUM
        else
            lp = fsi_inte_set%prop_adjn(i+1) - 1
        endif
        do j = sp, lp
            if (j == lp) then
                PROPEL_EDGE(:,j) = (/ j, sp /)
            else
                PROPEL_EDGE(:,j) = (/ j, j+1 /)
            endif
            PROPEL_PATCH(j) = i
        enddo
    enddo
        if (tempflag) then
        call fsi_set_prop_pair(Dim, prop_num, prop_coord)
    else
        do i = 1, PROPEL_POINT_NUM
            fsi_inte_set%prop_pair(:,i) = (/ i, i/)
        enddo
    endif
elseif (Dim == 3) then
    write (*,*) 'PROPEL_POINT_NUM:', PROPEL_POINT_NUM
    allocate (bound_node(bn))
    call get_skin_surface_info2(1, bn, bound_node, bfn, PROPEL_EDGE)
    deallocate (bound_node)
    call get_skin_surface_info4(1, ng, bn, PROPEL_PATCH)
endif

if (num_rn == 2) then
    if (Dim == 2) then
        num = region(2)%num_skin
        CASE_POINT_NUM = num
        CASE_EDGE_NUM = num
        fsi_inte_set%case_point_num = num
    elseif (Dim == 3) then
        call get_skin_surface_info1(2, bn, bfn)
        CASE_POINT_NUM = bn
        CASE_EDGE_NUM = bfn
        fsi_inte_set%case_point_num = bn
    endif
    allocate (CASE_EDGE(f2n_num,CASE_EDGE_NUM), CASE_POINT(Dim,CASE_POINT_NUM))
    allocate (CASE_DISP(Dim,CASE_POINT_NUM), CASE_PLOAD(Dim,CASE_POINT_NUM))
    allocate (fsi_inte_set%case_disp(Dim,CASE_POINT_NUM))
    allocate (fsi_inte_set%case_point(Dim,CASE_POINT_NUM))
    allocate (fsi_inte_set%case_pre_pload(Dim,CASE_POINT_NUM))
    CASE_EDGE = 0
    CASE_POINT = 0.0
    CASE_DISP = 0.0
    CASE_PLOAD = 0.0
    fsi_inte_set%case_point = 0.0
    fsi_inte_set%case_disp = 0.0
    fsi_inte_set%case_pre_pload = 0.0
    
    if (Dim == 2) then
        do i = 1, CASE_POINT_NUM
            if (i == CASE_POINT_NUM) then
                CASE_EDGE(:,i) = (/ i, 1 /)
            else
                CASE_EDGE(:,i) = (/ i, i+1 /)
            endif
        enddo
    elseif (Dim == 3) then
        write (*,*) 'CASE_POINT_NUM:', CASE_POINT_NUM
        allocate (bound_node(bn))
        call get_skin_surface_info2(2, bn, bound_node, bfn, CASE_EDGE)
        deallocate (bound_node)
    endif
else
    CASE_POINT_NUM = 1
    CASE_EDGE_NUM = 1
    allocate (CASE_EDGE(f2n_num,CASE_EDGE_NUM), CASE_POINT(Dim,CASE_POINT_NUM))
    allocate (CASE_DISP(Dim,CASE_POINT_NUM), CASE_PLOAD(Dim,CASE_POINT_NUM))
    allocate (fsi_inte_set%case_disp(Dim,CASE_POINT_NUM))
    allocate (fsi_inte_set%case_point(Dim,CASE_POINT_NUM))
    allocate (fsi_inte_set%case_pre_pload(Dim,CASE_POINT_NUM)) 
    CASE_EDGE = 0
    CASE_POINT = 0.0
    CASE_DISP = 0.0
    CASE_PLOAD = 0.0
    fsi_inte_set%case_disp = 0.0 
    fsi_inte_set%case_point = 0.0
    fsi_inte_set%case_pre_pload = 0.0
endif
    
end subroutine alloc_fsi_info
                  
                          
                          
subroutine fsi_set_prop_pair(Dim, prop_num, prop_coord)

implicit none

integer, intent(in) :: DIm, prop_num
real, intent(in) :: prop_coord(2,prop_num)

integer :: i, j, temp_num, skin_num
real :: dis, min_dis
real, allocatable :: skin_coord(:,:)
character(len=50) :: filename

skin_num = region(1)%num_skin
allocate (skin_coord(2,skin_num))
if (skin_num /= prop_num) then
    write (*,*) 'prop_num:', skin_num, prop_num
    Stop "Number of skin node(Solid) /= NUmber of propel point(Surface)"
endif
do i = 1, skin_num
    skin_coord(1,skin_num-i+1) = region(1)%node_coord(2,region(1)%skin_nodes(i))
    skin_coord(2,skin_num-i+1) = region(1)%node_coord(1,region(1)%skin_nodes(i))
enddo

fsi_inte_set%prop_pair = 0
do i = 1, skin_num
    min_dis = 10e8
    do j = 1, prop_num
        call calc_len(skin_coord(:,i), prop_coord(:,j), dis)
        if (dis < 10e-8) then
            fsi_inte_set%prop_pair(:,i) = (/ i, j /)
            exit
        elseif (min_dis > dis) then
            min_dis = dis
            temp_num = j
        endif
    enddo
    if (fsi_inte_set%prop_pair(1,i) == 0) then
        !write (*,'(A)') 'Warning: Cannot find matching point'
        !write (*,'(A,I5,A,I5,A,ES15.7,A)') 'Temporory:', region(1)%skin_nodes(i), '<->', temp_num, '(Dis:', min_dis, ')'
        fsi_inte_set%prop_pair(:,i) = (/ i, temp_num /)
    endif
enddo

!filename = "./output/solid/ale/check_pair000.plt"
!!write ( filename(30:32), '(I3.3)' ) n_remesh
!open (Unit=30, File=filename, STATUS='replace', ACTION='write')
!Write (30,'(A,A,A)') 'TITLE="Check fsi pair"'
!Write (30,*) 'VARIABLES="x", "y", "Solid", "Surface"'
!Write (30,'(A,I5,A,I5, A)') 'ZONE N =', skin_num , ', E =', skin_num, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
!do i = 1, skin_num
!    write (30,'(2ES15.7,2I5)') skin_coord(:,i), fsi_inte_set%prop_pair(:,i)
!enddo
!do i = 1, skin_num
!    if (i == skin_num) then
!        write (30,*) i, 1
!    else
!        write (30,*) i, i+1
!    endif
!enddo
!close(30)
deallocate (skin_coord)

end subroutine fsi_set_prop_pair

subroutine convey_pair_fsi_info(Dim, flag, prop_num, prop_coord, prop_load, prop_disp)

implicit none

integer, intent(in) :: Dim, flag, prop_num
real, intent(inout) :: prop_coord(2,prop_num), prop_load(2,prop_num), prop_disp(2,prop_num)

integer :: i, pair(2)
real :: temp_val(2,prop_num)

if (flag == 1) then
    temp_val = prop_coord
    do i = 1, prop_num
        pair = fsi_inte_set%prop_pair(:,i)
        prop_coord(:,pair(1)) = temp_val(:,pair(2))
    enddo
elseif (flag == 2) then
    temp_val = prop_load
    do i = 1, prop_num
        pair = fsi_inte_set%prop_pair(:,i)
        prop_load(:,pair(1)) = temp_val(:,pair(2))
    enddo
elseif (flag == 3) then
    temp_val = prop_disp
    do i = 1, prop_num
        pair = fsi_inte_set%prop_pair(:,i)
        prop_disp(:,pair(2)) = temp_val(:,pair(1))
    enddo
endif
          
end subroutine convey_pair_fsi_info
!=============================================================================

subroutine free_fsi_domain

if (allocated(fsi_inte_set%fsi_inte_node)) then
    deallocate(fsi_inte_set%fsi_inte_node, fsi_inte_set%B_flag, fsi_inte_set%ni_flag)
    deallocate(fsi_inte_set%B_rate, fsi_inte_set%fsi_coord, fsi_inte_set%new_fsi_node)
endif
if (allocated(fsi_inte_set%old_fsi_inte)) deallocate(fsi_inte_set%old_fsi_inte)
deallocate(fsi_inte_set%prop_pre_pload)
if (allocated(fsi_inte_set%case_pre_pload)) deallocate(fsi_inte_set%case_pre_pload)

call free_history
call free_pdomain

end subroutine free_fsi_domain
!=============================================================================

end module fsi_domain