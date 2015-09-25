subroutine update_surface(Dim, case_check, PROPEL_POINT_NUM, PROPEL_POINT, CASE_POINT_NUM, CASE_POINT, &
                          nn, coord, ng, PROPEL_PATCH, PROPEL_REMESH_FLAG, max_move_dis) 
  
use surface_domain

implicit none

integer, intent(in) :: Dim, nn
real, intent(inout) :: coord(Dim,nn), max_move_dis
integer, intent(in) :: PROPEL_POINT_NUM, CASE_POINT_NUM, ng, PROPEL_REMESH_FLAG(ng)
integer, intent(in) :: PROPEL_PATCH(PROPEL_POINT_NUM)
real, intent(in) :: PROPEL_POINT(Dim,PROPEL_POINT_NUM), CASE_POINT(Dim,CASE_POINT_NUM)
logical, intent(in) :: case_check

integer :: i, j, k, temp_cp, cp, check_arr_num, num, num_group, sp, adjn, num_patch, num_skin
integer :: check_arr(PROPEL_POINT_NUM), temp_patch(PROPEL_POINT_NUM), temp_remesh_flag(ng)
real :: vec(Dim), dis, min_dis, temp_coord(2,nn)
logical :: check, separate_check

move_flag = 0
if (Dim == 2) then
    num_patch = maxval(PROPEL_PATCH)
    !separate_check = .false.
    !do i = 1, num_patch
    !    if (PROPEL_REMESH_FLAG(i) == 4) then
    !        separate_check = .true.
    !        exit
    !    endif
    !enddo
    
    abla_flag = 0
    if (case_check) then
        mesh_set(2)%nn = region(2)%num_nodes
        mesh_set(2)%ne = region(2)%num_elements
        if (allocated(mesh_set(2)%coord)) deallocate(mesh_set(2)%coord, mesh_set(2)%ele)
        allocate(mesh_set(2)%coord(Dim,mesh_set(2)%nn), mesh_set(2)%ele(4,mesh_set(2)%ne))
        mesh_set(2)%coord = region(2)%node_coord(1:2,:)
        mesh_set(2)%ele = region(2)%element(4:7, :)
        
        max_move_dis = 0.0
        do i = 1, region(2)%num_skin
            cp = region(2)%skin_nodes(i)
            do j = 1, inte_set%cn
                if (inte_set%c_nodes(j) == cp) then
                    vec = mesh_set(2)%coord(:,cp) - CASE_POINT(:,i)
                    dis = sqrt(vec(1)**2 + vec(2)**2)
                    if (max_move_dis < dis) max_move_dis = dis
                    mesh_set(2)%coord(:,cp) = CASE_POINT(:,i)
                    exit
                endif
            enddo
        enddo
        
        if (max_move_dis > mesh_set(2)%min_d(1)*10e-8) abla_flag = 1
    endif
    
    do i = 1, PROPEL_POINT_NUM
        !if (separate_check == .false. .and. PROPEL_PATCH(i) == ng+1) then
        !    PROPEL_PATCH(i) = ng
        !endif
        temp_patch(PROPEL_POINT_NUM-i+1) = num_patch - PROPEL_PATCH(i) + 1
    enddo
    do i = 1, num_patch
        temp_remesh_flag(num_patch-i+1) = PROPEL_REMESH_FLAG(i)
    enddo
    
    max_move_dis = 0.0
    if (sum(PROPEL_REMESH_FLAG(1:num_patch)) == 0) then
        do i = 1, region(1)%num_skin
            cp = region(1)%skin_nodes(i)
            vec = coord(:,cp) - PROPEL_POINT(:,i)
            dis = sqrt(vec(1)**2 + vec(2)**2)
            if (max_move_dis < dis) max_move_dis = dis
            coord(:,cp) = PROPEL_POINT(:,i)
        enddo
    else
        !write (*,*) 'PROPEL_REMESH_FLAG: ', temp_Remesh_flag
        !write (*,*) 'PROPEL_PATCH'
        !write (*,*) temp_patch(1:PROPEL_POINT_NUM)
        !num = 0
        !do  ! update coordinate
            !num = num + 1
            !if (num > PROPEL_POINT_NUM) exit  ! update coordinate
            !if (temp_remesh_flag(temp_patch(num)) == 0) then
            !    num_group = temp_patch(num)
            !    num = num - 1
            !    call set_skin(1, num_group, sp, adjn)
            !    !write (*,*) 'number_group, sp, adjn: ', num_group, sp, adjn
            !    do j = sp, sp+adjn-1
            !        num = num + 1
            !        cp = region(1)%skin_nodes(j)
            !        coord(:,cp) = PROPEL_POINT(:,num)
            !    enddo
            !endif
        !enddo  ! update coordinate
        temp_coord = coord
        num_skin = region(1)%num_skin
        do i = 1, PROPEL_POINT_NUM
            if (temp_remesh_flag(temp_patch(i)) == 0) then
                min_dis = 10e8
                num = 0
                do j = 1, num_skin
                    cp = region(1)%skin_nodes(j)
                    call calc_len(temp_coord(:,cp), PROPEL_POINT(:,i), dis)
                    if (min_dis > dis) then
                        min_dis = dis
                        num = cp
                    endif
                enddo
                if (num == 0) then
                    stop 'Update_surface(Error): Cannot find the near PROPEL_POINT at skin_node'
                else
                    coord(:,num) = PROPEL_POINT(:,i)
                endif
            endif
        enddo
    endif
    if (max_move_dis > mesh_set(1)%min_d(1)*10e-8) move_flag = 1
elseif (Dim == 3) then
    max_move_dis = 0.0
    do i = 1, mesh_set(1)%bn
        cp = mesh_set(1)%bound_node(i)
        vec = coord(:,cp) - PROPEL_POINT(:,i)
        dis = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
        if (max_move_dis < dis) max_move_dis = dis
        coord(:,cp) = PROPEL_POINT(:,i)
    enddo
    if (max_move_dis > 10e-6) move_flag = 1
endif


end subroutine update_surface
    