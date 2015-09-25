subroutine surface_control(Dim, num_rn, ng, current_step, n_remesh, &
                           PROPEL_POINT_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_PATCH, PROPEL_REMESH_FLAG, &
                           CASE_POINT_NUM, CASE_POINT, CASE_DISP)

use surface_domain
use file_info_house
use system_info_house

implicit none

integer, intent(in) :: Dim, num_rn, ng, current_step
integer, intent(inout) :: n_remesh, PROPEL_REMESH_FLAG(ng)
integer, intent(in) :: PROPEL_POINT_NUM, CASE_POINT_NUM, PROPEL_PATCH(PROPEL_POINT_NUM)
real, intent(in) :: PROPEL_POINT(Dim,PROPEL_POINT_NUM), PROPEL_DISP(Dim,PROPEL_POINT_NUM)
real, intent(in) :: CASE_POINT(Dim,CASE_POINT_NUM), CASE_DISP(Dim,CASE_POINT_NUM)

integer :: nn, ne, ele_node, adjn, del_pn, count, i, rn, sp(2), lp(2), sub_nn, sub_ne
real :: tol_skew, tol_kurt, max_move_dis
real :: temp_prop_coord(Dim,PROPEL_POINT_NUM), temp_case_coord(Dim,CASE_POINT_NUM)
real, allocatable :: coord(:,:), sub_coord(:,:)
integer, allocatable :: el_data(:,:), sub_elem(:,:)
integer, allocatable :: temp_adj(:), del_p(:) 
real :: env_d, skew, kurt, area, tol_env_d, det_jacobi, tol_jacobi
logical :: remesh_check, thin_check(ng), non_matching_check(ng), delete_subdomain_check

! remesh_flag : 0 - don`t remesh
!               1 - seperation remesh
!               2 - the other remesh(evaluation, length between boundary nodes)

move_flag = 0
if (Dim == 2) then
    remesh_flag = 0
    rn = 1
    mesh_set(rn)%flag = 0
    nn = region(rn)%num_nodes
    ne = region(rn)%num_elements
    ele_node = 4
    allocate(coord(2,nn), el_data(4,ne))
    coord = region(1)%node_coord(1:2, :)
    el_data = region(rn)%element(4:7, :)
    
    !call move_surface(num_rn, rn, nn, coord, remesh_flag(rn), mp, move_flag)
    do i = 1, PROPEL_POINT_NUM
        temp_prop_coord(1,PROPEL_POINT_NUM-i+1) = PROPEL_POINT(2,i) ! - PROPEL_DISP(2,i)
        temp_prop_coord(2,PROPEL_POINT_NUM-i+1) = PROPEL_POINT(1,i) ! - PROPEL_DISP(1,i)
    enddo
    
    !write (*,*) 'current_step:', current_step
    !write (*,*) 'SUM(PROPEL_REMESH_FLAG): ', sum(PROPEL_REMESH_FLAG(1:ng))
    if (ablation_flag) then
        do i = 1, CASE_POINT_NUM
            temp_case_coord(1,CASE_POINT_NUM-i+1) = CASE_POINT(2,i) !- CASE_DISP(2,i)
            temp_case_coord(2,CASE_POINT_NUM-i+1) = CASE_POINT(1,i) !- CASE_DISP(1,i)
        enddo
    endif
    
    !call write_tecplot(PROPEL_POINT_NUM, temp_prop_coord, ne, el_data, rn, current_step, 7)
    thin_check = .false.
    do i = 1, ng
        if (PROPEL_REMESH_FLAG(i) == 3) then
            write (*,'(A,I4)') ' >>    Delete subdomain: ', i
            exit
        endif
        if (PROPEL_REMESH_FLAG(i) == 11) then
        !if (thin_chk(i) == .false.) then
            write (*,'(A)') ' >>    Thin check: True'
			thin_check(i) = .true.
            PROPEL_REMESH_FLAG(i) = 0
            exit
        endif
    enddo
    call update_surface(Dim, ablation_flag, PROPEL_POINT_NUM, temp_prop_coord, CASE_POINT_NUM, temp_case_coord, &
                        nn, coord, ng, PROPEL_PATCH, PROPEL_REMESH_FLAG, max_move_dis) 
    
    if (n_remesh == 2 .or. current_step <= 1) then
        accumulate_move_dis = max_move_dis
    else
        accumulate_move_dis = accumulate_move_dis + max_move_dis
    endif
    
    !write (*,*) 'SUM(PROPEL_REMESH_FLAG): ', sum(PROPEL_REMESH_FLAG(1:ng))
    !if (is_coarse_problem) then
    !    call check_non_matching_disp(ng, non_matching_check, current_step)
    !    do i = 1, ng
    !        if (non_matching_check(i)) then
    !            remesh_flag(rn) = 2
    !            exit
    !        endif
    !    enddo
    !endif
    if (move_flag == 1 .and. sum(PROPEL_REMESH_FLAG(1:ng)) == 0 .and. remesh_flag(rn) == 0)  then
        !if (num_rn >= 2) then
        !    print *, " >> Adjust contact points"
        !    call adjust_contact(nn, coord)
        !endif
        call coord_interporation(Dim, rn, ng, nn, coord, thin_check)
    endif
    call update_mesh_info(Dim, 1, nn, ne, coord, el_data, 1)
    
    remesh_check = .false.
    if (move_flag == 1) then
        remesh_check = .true.
        if (sum(PROPEL_REMESH_FLAG(1:ng)) /= 0) then
            remesh_check = .false.
        elseif (0.3*mesh_set(1)%cri_d > accumulate_move_dis) then
            remesh_check = .false.
        else
            do i = 1, ng
                if (thin_check(i) .or. non_matching_check(i)) then
                    remesh_check = .false.
                    exit
                endif
            enddo
        endif
    endif
    
    !call write_tecplot( nn, coord, ne, el_data, 1, current_step, 4 )
    !call write_tecplot( region(2)%num_nodes, region(2)%node_coord(1:2,:), region(2)%num_elements, region(2)%element(4:7,:), 2, current_step, 4 )
    !call write_tecplot( nn, coord, ne, el_data, 1, current_step, 4 )
    
    write (*,*) 
    write (*,'(A,I7,A)') ' >> Current step: ', current_step, ', surface control info.'
    write (*,'(A,F12.7)') ' >>    Maximum moving distance:', max_move_dis
    write (*,'(A,F12.7)') ' >> Accumulate moving distance:', accumulate_move_dis
    write (*,'(A,F12.7)') ' >>   Criteria moving distance:', 0.3*mesh_set(1)%cri_d
    write (*,*) 
    if (move_flag == 1) then
        write (*,'(A,I2,A)') ' >>    Move flag: ', move_flag, ' (Run eulerian stage)'
    else
        write (*,'(A,I2,A)') ' >>    Move flag: ', move_flag, ' (Do not run eulerian stage)'
    endif
    if (remesh_check) then
        write (*,'(A)') " >> Remesh check: Yes"
    else
        write (*,'(A)') " >> Remesh check: No"
        do i = 1, ng
            if (thin_check(ng)) write (*,'(A,I2,A)') "   > Thin check(", ng, "): Yes "
        enddo
        do i = 1, ng
            if (non_matching_check(ng)) write (*,'(A,I2,A)') "   > Non-matching disp check(", ng, "): Yes "
        enddo
    endif
    write (*,'(A,L)') " >> ablation_flag:", ablation_flag

    !remesh check
    !if (current_step == 2) then
    !    move_flag = 1
    !    remesh_check = .true.
    !endif
    n_remesh = 0
    if (move_flag == 1) then
        if ( remesh_flag(rn) == 1 ) then
            write (*,*) '>> Separate the object!!'
            allocate (surface_set(1)%coord(2,nn))
            surface_set(1)%nn = nn
            surface_set(1)%coord = coord
        else
            if (remesh_check) then
                do i = 1, ng
                    allocate (surface_set(i)%coord(2,nn))
                    surface_set(i)%nn = nn
                    surface_set(i)%coord = coord
        
                    call set_skin(rn, i, sp(1), adjn)
                    allocate (temp_adj(adjn), del_p(adjn))
                    temp_adj = region(rn)%skin_nodes(sp(1):sp(1)+adjn-1)
                    del_pn = 0;  del_p = 0
        
                    ! check total length
                    !call calc_tol_length(region(rn)%num_nodes, region(rn)%node_coord, adjn, temp_adj, tol_dis)
                    !write (*,*) i, ' total length :', tol_dis, mesh_set(2)%ave_d(1)*5.
                    !if ( num_rn >= 2 ) then
                    !    limit_dis = mesh_set(2)%ave_d(1)*5.
                    !else 
                    !    limit_dis = mesh_set(1)%ave_d(i)*12.
                    !endif
                    !if (tol_dis < limit_dis) then
                    !    !write (*,*) 'Delete subdomain : ', mn
                    !    remesh_flag(rn) = 2
                    !    mesh_set(rn)%flag(i) = 4
                    !endif
        
                    ! check sharp of corner point
                    !if (mesh_set(rn)%flag(i) /= 4) then
                    !    call sharp_check(rn, i, nn, coord, adjn, temp_adj)   
                    !endif
                    !if (mesh_set(rn)%flag(i) == 0 ) then
                    !    ! check the angle between boudary nodes
                    !    !call check_boundary_angle(rn, nn, coord, adjn, temp_adj, del_pn, del_p)
                    !    ! check the length between boudary nodes
                    !    call check_boundary_length(rn, i, nn, coord, adjn, temp_adj, del_pn, del_p)
                    !endif
            
                    if(mesh_set(rn)%flag(i) == 0) then
                        call evaluation_mesh(rn, i, nn, coord, env_d, skew, kurt, area)
                        !if (mod(current_step,80) == 0) skew = 65.0
                        if (num_rn == 1) then
                            tol_env_d = 1.25*mesh_set(rn)%ave_d(i)
                            tol_skew = 65.0
                        else
                            tol_env_d = 1.5*mesh_set(rn)%ave_d(i)
                            if (ng >= 2) then
                                tol_skew = 70.0
                            else
                                tol_skew = 65.0
                            endif
                        endif
                        tol_kurt = 0.1
                        write (*,'(A,I2,A)') '             <<Evaluate ', i, ' mesh>>'
                        write (*,'(A)') '-------------------------------------------------'
                        write (*,'(A)') '      |  env_d  |    Max   |    Min   |    Min   |'
                        write (*,'(A)') '      |         | Skewness | Kurtosis |   area   |'
                        write (*,'(A)') '-------------------------------------------------'
                        write (*,'(A,F7.4,A,F7.2,A)') ' Cri  | ', tol_env_d, ' |  ',   tol_skew, ' |        0.10 |             |'
                        write (*,'(A,F7.4,A,F5.2,A,ES11.4,A,ES11.4,A)') ' Eval | ',env_d,' |    ',skew,' |    ',kurt,' |',area,' |'
                        write (*,'(A)') '-------------------------------------------------'
                
                        !if (current_step == 2) skew = 100.0
                        if ((skew >= tol_skew .OR. env_d >= tol_env_d) .AND. rn == 1) then
                    
                            if (env_d >= tol_env_d) then
                                write (*,*) mesh_set(rn)%ave_d(i)
                                write (*,'(5x,A)') 'Env_d : Remesh srart!!!'
                            elseif ( skew >= tol_skew ) then
                                write (*,'(5x,A)') 'Skewness : Remesh srart!!!'
                            else
                                write (*,'(5x,A)') 'Kurtosis : Remesh srart!!!'
                            endif
                            remesh_flag(rn) = 2
                            mesh_set(rn)%flag(i) = 1
                            call adjust_corner(rn, i, nn, coord, adjn, temp_adj, ng, 1)
                        endif
                        write (*,*) " "
                    endif
                    deallocate (temp_adj, del_p)
                enddo
            endif
        endif
    endif
    call find_inte_coord(nn, coord)

    deallocate(coord, el_data)
elseif (Dim == 3) then
    remesh_flag = 0
    rn = 1
    mesh_set(rn)%flag = 0
    nn = mesh_set(rn)%nn
    ne = mesh_set(rn)%ne
    ele_node = 8
    allocate(coord(Dim,nn), el_data(ele_node,ne))
    coord = mesh_set(rn)%node_coord
    el_data = mesh_set(rn)%element
    
    temp_prop_coord = PROPEL_POINT ! - PROPEL_DISP
    call update_surface(Dim, ablation_flag, PROPEL_POINT_NUM, temp_prop_coord, CASE_POINT_NUM, temp_case_coord, &
                        nn, coord, ng, PROPEL_PATCH, PROPEL_REMESH_FLAG, max_move_dis) 
    
    if (n_remesh == 2 .or. current_step <= 1) then
        accumulate_move_dis = max_move_dis
    else
        accumulate_move_dis = accumulate_move_dis + max_move_dis
    endif
    
    delete_subdomain_check = .FALSE.
    do i = 1, ng
        if (PROPEL_REMESH_FLAG(i) == 3) then
            delete_subdomain_check = .TRUE.
            exit
        endif
    enddo    
    write (*,*) 'move_flag, remesh_flag:', move_flag, PROPEL_REMESH_FLAG(1:ng)
    !if (move_flag == 1 .and. sum(PROPEL_REMESH_FLAG(1:ng)) == 0) then
    if (move_flag == 1) then
        call coord_interporation(Dim, rn, ng, nn, coord, thin_check)
        if (divided_mesh_region >= 2 .and. delete_subdomain_check == .FALSE.) then
            if (remesh_type_number == 1) then
                call ale_region_smoothing(nn, coord)
            elseif (remesh_type_number == 2) then
                call ale_region_smoothing2(nn, coord, 2)
            endif
        endif
        call write_tecplot_3D(nn, coord, ne, el_data, 1, current_step, 4)
    endif
    call update_mesh_info(Dim, 1, nn, ne, coord, el_data, 1)
    
    remesh_check = .false.
    if (move_flag == 1) then
        remesh_check = .true.
        do i = 1, ng
            if (PROPEL_REMESH_FLAG(i) /= 0) then
                remesh_check = .false.
            endif
        enddo
        if (remesh_check) then
            if (0.3*mesh_set(1)%cri_d > accumulate_move_dis) then
                remesh_check = .false.
            else
                if (delete_subdomain_check) then 
                    remesh_check = .false.
                endif
            endif
        endif
    endif
    
    n_remesh = 0
    write (*,*) 
    write (*,'(A,I7,A)') ' >> Current step: ', current_step, ', surface control info.'
    if (remesh_check) then
        write (*,'(A)') " >> Remesh check: Yes"
    else
        write (*,'(A)') " >> Remesh check: No"
        do i = 1, ng
            if (thin_check(ng)) write (*,'(A,I2,A)') "   > Thin check(", ng, "): Yes "
        enddo
    endif
    write (*,'(A,L)') " >> ablation_flag:", ablation_flag
    if (remesh_check) then
        tol_jacobi = 0.2
        do i = 1, ng
            sp = mesh_set(rn)%sub_num(i,1:2)
            lp = mesh_set(rn)%sub_num(i+1,1:2) - 1
            sub_nn = lp(1) - sp(1) + 1
            sub_ne = lp(2) - sp(2) + 1
            allocate(sub_coord(3,sub_nn), sub_elem(8,sub_ne))
            sub_coord = coord(:,sp(1):lp(1))
            sub_elem = el_data(:,sp(2):lp(2))
            if (ng /= 1) sub_elem = sub_elem - (mesh_set(rn)%sub_num(i,1)-1)
            call calc_elem_jacobian(sub_nn, sub_ne, sub_coord, sub_elem, det_jacobi)
            write (*,'(A,I2,A,I2,A)') '<<Evaluate ', rn, ' region, ', i, ' sub domain 3D mesh>>'
            write (*,'(A)') ' ------------------------------------------'
            write (*,'(A)') '|                     Cri          Eval    |'
            write (*,'(A,F7.4,6X,F7.4,A)') '| det_jacobi  =    ', tol_jacobi, det_jacobi,  '    |'
            write (*,'(A)') ' ------------------------------------------'
            !if (current_step == 5 .or. current_step == 10 .and. i == 1) det_jacobi = 0.1
            if (det_jacobi <= tol_jacobi) then
                remesh_flag(rn) = 2
                mesh_set(rn)%flag(i) = 1    
                !call quad2tri(current_step, rn, i)
            endif
            deallocate(sub_coord, sub_elem)
        enddo
        if (move_flag == 1) then
            write (*,'(A,I2,A)') ' >>    Move flag: ', move_flag, ' (Run eulerian stage)'
        else
            write (*,'(A,I2,A)') ' >>    Move flag: ', move_flag, ' (Do not run eulerian stage)'
        endif
        deallocate (coord, el_data)
    endif
endif
!if (current_step == 5) stop
if (remesh_flag(rn) > 0) then
    if (Dim == 2) then
        n_remesh = 2
    elseif (Dim == 3) then
        n_remesh = 3
    endif
elseif (sum(PROPEL_REMESH_FLAG(1:ng)) /= 0) then
    n_remesh = 3
else
    n_remesh = 0
endif
write (*,'(A,I)') " >> remesh number:", n_remesh
write (*,'(A)') " >> Complete surface_control"
contains

!=====================================================================================
subroutine evaluation_mesh(rn, ng, nn, node, env_d, skew, kurt, tol_area)

implicit none

integer, intent(in) :: rn, nn, ng
real, intent(in) :: node(2,nn)
real, intent(inout) :: env_d, skew, kurt, tol_area

integer :: i, j, ne, adjn, ne1, ne2, sp, ngroup
real :: angle, s, area_ne, min_area, length, ave_d
real :: temp_p(2,5), dis(3), area(4)
integer, allocatable :: ele(:,:), adj(:)
real, allocatable :: skew_ne(:), kurt_ne(:)
!-----------------------------------------------------------------------------------------
ngroup = region(rn)%num_groups
ne1 = region(rn)%element_group(ng)
if ( ng == ngroup ) then
    ne2 = region(rn)%num_elements
else
    ne2 = region(rn)%element_group(ng+1) - 1
endif
ne = ne2 - ne1 + 1

allocate ( ele(4,ne), skew_ne(ne), kurt_ne(ne) )
ele = 0;  skew_ne = 0.;   kurt_ne = 0.;  tol_area = 0.
ele = region(rn)%element(4:7,ne1:ne2)
ave_d = mesh_set(rn)%ave_d(ng)
env_d = 0.0

call set_skin(rn, ng, sp, adjn)
allocate( adj(adjn) )
adj = region(rn)%skin_nodes(sp:sp+adjn-1) 

do i = 1, adjn
    if (i == adjn) then
        call calc_len(node(:,adj(i)), node(:,adj(1)), length)
    else
        call calc_len(node(:,adj(i)), node(:,adj(i+1)), length)
    endif

    if ( env_d < abs(length - ave_d) ) env_d = abs(length - ave_d)
enddo

do i = 1, ne
    temp_p = 0.
    do j = 1, 4
        if ( j /= 4 ) then
            temp_p(1,j) = (node(1,ele(j,i))+node(1,ele(j+1,i)))*0.5
            temp_p(2,j) = (node(2,ele(j,i))+node(2,ele(j+1,i)))*0.5
        else
            temp_p(1,j) = (node(1,ele(j,i))+node(1,ele(1,i)))*0.5
            temp_p(2,j) = (node(2,ele(j,i))+node(2,ele(1,i)))*0.5
        endif
    enddo
    call calc_cross_point(temp_p(:,1), temp_p(:,3), temp_p(:,2), temp_p(:,4), temp_p(:,5))
    call calc_angle(temp_p(:,1), temp_p(:,5), temp_p(:,4), angle)
    skew_ne(i) = abs(90. - angle)
    
    !call calc_cross_point(node(:,ele(1,i)), node(:,ele(3,i)), node(:,ele(2,i)), node(:,ele(4,i)), temp_p(:,5))
    !do j = 1, 4
    !    dis = 0.
    !    if ( j /= 4 ) then
    !        call calc_len(node(:,ele(j,i)), node(:,ele(j+1,i)), dis(1))
    !        call calc_len(node(:,ele(j,i)), temp_p(:,5), dis(2))
    !        call calc_len(node(:,ele(j+1,i)), temp_p(:,5), dis(3))
    !    else
    !        call calc_len(node(:,ele(j,i)), node(:,ele(1,i)), dis(1))
    !        call calc_len(node(:,ele(j,i)), temp_p(:,5), dis(2))
    !        call calc_len(node(:,ele(1,i)), temp_p(:,5), dis(3))
    !    endif
    !
    !    s = (dis(1)+dis(2)+dis(3))*0.5
    !    area(j) = sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
    !enddo
    !area_ne = sum(area(1:4))
    do j = 1, 4
        temp_p(:,j) = node(:,ele(j,i))
    enddo
    call calc_area(4, temp_p(:,1:4), area_ne)
    min_area = minval(area(1:4))
    kurt_ne(i) = 4. * (min_area/area_ne)
    tol_area = tol_area + area_ne
enddo

if ( ne == 0 ) then
    skew = 0.
    kurt = 0.
    tol_area = 0.
else
    skew = maxval(skew_ne(:))
    kurt = minval(kurt_ne(:))
    !write (0822,'(I4,3X,4F8.3)') current_step, skew, sum(skew_ne(:))/float(ne), kurt, sum(kurt_ne(:))/float(ne)
endif

deallocate (ele, skew_ne, kurt_ne)
deallocate(adj)

end subroutine evaluation_mesh



end subroutine surface_control

!************************************
! surface -> solid function start
subroutine surface_coord_update(ng, PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_REMESH_FLAG)

use surface_domain
implicit none 

INTEGER, INTENT(in) :: ng
INTEGER, INTENT(in) :: PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_REMESH_FLAG(ng)
REAL,    INTENT(in) :: PROPEL_POINT(2,PROPEL_POINT_NUM), PROPEL_DISP(2,PROPEL_POINT_NUM)

integer :: i, j, sp, lp, adjn, num, accumulate_adjn
integer, allocatable :: adj(:)

remesh_flag(1) = 2
accumulate_adjn = 0
do i = 1, ng
    if (allocated(surface_set(i)%coord)) deallocate(surface_set(i)%coord)
    allocate(surface_set(i)%coord(2,PROPEL_POINT_NUM))
    surface_set(i)%nn = PROPEL_POINT_NUM
    num = 0
    do j= PROPEL_POINT_NUM, 1, -1
        num = num + 1
        !surface_set(i)%coord(1,num) = PROPEL_POINT(2,j) - fsi_inte_set%prop_disp(2,i)
        !surface_set(i)%coord(2,num) = PROPEL_POINT(1,j) - fsi_inte_set%prop_disp(1,i)
        surface_set(i)%coord(1,num) = PROPEL_POINT(2,j) - PROPEL_DISP(2,i)
        surface_set(i)%coord(2,num) = PROPEL_POINT(1,j) - PROPEL_DISP(1,i)
    enddo
    
    sp = accumulate_adjn + 1
    if (i == ng) then
        lp = PROPEL_POINT_NUM
    else
        lp = PROPEL_POINT_NUM - fsi_inte_set%prop_adjn(i+1) + 1
    endif
    
    adjn = lp - sp + 1
    accumulate_adjn = accumulate_adjn + adjn
    allocate(adj(adjn))
    num = 0
    do j = sp, lp
        num = num + 1
        adj(num) = j
    enddo
    
    remesh_flag(1) = 1
    mesh_set(1)%flag(i) = PROPEL_REMESH_FLAG(i)
    deallocate(adj)
enddo
close(30)
end subroutine surface_coord_update



subroutine rearrange_edge(Dim, current_step, ng, PROPEL_REMESH_FLAG, PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_EDGE, &
          PROPEL_PATCH, PROPEL_CORNER, PROPEL_EDGE_LENGTH, PROPEL_RIDGE_GROUP_NUM, PROPEL_RIDGE_NUM, PROPEL_RIDGE)

use surface_domain
use remesh_domain
implicit none

integer, intent(in) :: Dim, current_step, ng, PROPEL_REMESH_FLAG(ng)
integer, intent(in) :: PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_RIDGE_GROUP_NUM
integer, intent(in) :: PROPEL_CORNER(PROPEL_POINT_NUM), PROPEL_PATCH(PROPEL_EDGE_NUM)
integer, intent(in) :: PROPEL_RIDGE_NUM(PROPEL_RIDGE_GROUP_NUM), PROPEL_RIDGE(PROPEL_RIDGE_GROUP_NUM,1000)
integer, intent(inout) :: PROPEL_EDGE((Dim-1)*2,PROPEL_EDGE_NUM)
real, intent(in) :: PROPEL_EDGE_LENGTH(PROPEL_EDGE_NUM)
real, intent(inout) :: PROPEL_POINT(Dim,PROPEL_POINT_NUM), PROPEL_DISP(Dim,PROPEL_POINT_NUM)

integer :: i, j, k, q, n, rn, adjn, sp, ep, lp, np, nn, ne, num, cn1, cn2, num_skin, tol_cn, cnum, count, output_flag, remesh_type
integer :: num_patch, prop_bn, cur_ng, temp_num(10)
integer :: target_adjn, target_cn, tol_corn_group(ng)
integer :: temp_corner(PROPEL_POINT_NUM)
real :: temp_coord(2,PROPEL_POINT_NUM), temp_edge_length(PROPEL_EDGE_NUM), dis, cp(2), p(2)
character(len=60) :: fn
logical :: cont_check, odd_check

integer, allocatable :: prop_edge(:,:), prop_patch(:), prop_corner(:)
real, allocatable :: prop_coord(:,:), prop_edge_length(:)
integer, allocatable :: tol_corn(:), adj(:), cont_node(:), target_adj(:), target_cp(:)
real,    allocatable :: tol_edge(:,:), tol_coord(:,:), adj_coord(:,:)

output_flag = 1
remesh_type = 2
num_patch = maxval(PROPEL_PATCH)

rn = 1
remesh_flag(rn) = 2
temp_coord = PROPEL_POINT !- PROPEL_DISP

if (Dim == 2) then
    num_skin =  region(rn)%num_skin
    allocate (tol_corn(num_skin), tol_edge(2,num_skin), tol_coord(2,num_skin))
    tol_cn = 0
    tol_corn_group = 0
    tol_corn = 0
    tol_edge = 0.0
    tol_coord = 0.0
    odd_check = .false.
    write (*,*) 'rearrange_edge: check1(num_patch):', num_patch
    do i = 1, num_patch
        ! corner data
        cn1 = mesh_set(rn)%corn_group(i)
        if ( i == num_patch ) then
            cn2 = mesh_set(rn)%cn
        else
            cn2 = mesh_set(rn)%corn_group(i+1) - 1
        endif
        write (*,*) 'rearrange_edge: check2(ng):', i
        if (PROPEL_REMESH_FLAG(num_patch-i+1) == 0) then
            cnum = cn2 - cn1 + 1
            tol_corn(tol_cn+1:tol_cn+cnum) = mesh_set(rn)%corn(cn1:cn2)
            tol_edge(:,tol_cn+1:tol_cn+cnum) = mesh_set(rn)%edge(:,cn1:cn2)
            tol_coord(:,tol_cn+1:tol_cn+cnum) = mesh_set(rn)%corn_coord(:,cn1:cn2)
            tol_corn_group(i) = tol_cn + 1
            tol_cn = tol_cn + cnum
        else
            ! find first corner node
            adjn = 0
            do j = 1, PROPEL_EDGE_NUM
                if (num_patch-PROPEL_PATCH(j)+1 == i) then
                    adjn = adjn + 1
                    if (adjn == 1) then
                        sp = PROPEL_EDGE(1,j)
                        ne = j
                    endif
                endif
            enddo
            allocate (adj_coord(2,adjn))
            write (*,*) 'rearrange_edge: check3'
            if (output_flag == 1) then
                fn = './output/solid/remesh/check_rearrange0000000_0.plt'
                write (fn(38:44), '(I7.7)' ) current_step
                open(unit=30, file=fn, status='REPLACE', action='WRITE')
                Write (30,'(A)') 'TITLE="Check random edge"'
                Write (30,*) 'VARIABLES="x", "y", "corner"'
                Write (30,'(A,I5,A,I5, A)') 'ZONE N =', PROPEL_POINT_NUM , ', E =', PROPEL_EDGE_NUM, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
                do j = 1, PROPEL_POINT_NUM
                    write (30,*) PROPEL_POINT(:,j), PROPEL_CORNER(j)
                enddo
                do j = 1, PROPEL_EDGE_NUM
                    write (30,*) PROPEL_EDGE(:,j)
                enddo
                close(30)
            endif
        
            count = 1
            np = sp
            !write (*,*) 'sp:', sp
            rearrange_adj: do
                adj_coord(:,count) = temp_coord(:,np)
                temp_corner(count) = PROPEL_CORNER(np)
                temp_edge_length(count) = PROPEL_EDGE_LENGTH(ne)
                np = PROPEL_EDGE(2,ne)
                if (np == sp) then
                    if (count /= adjn) then
                        write (*,*) 'surface_control.f90 - Error(rearrange_edge): count /= adjn'
                        stop
                    endif
                    exit rearrange_adj
                endif
                do j = 1, PROPEL_EDGE_NUM
                    if (PROPEL_EDGE(1,j) == np) then
                        count = count + 1
                        ne = j
                        exit
                    endif
                enddo
            enddo rearrange_adj
        
            if (output_flag == 1) then
                fn = './output/solid/remesh/check_rearrange0000000_1.plt'
                write (fn(38:44), '(I7.7)' ) current_step
                open(unit=30, file=fn, status='REPLACE', action='WRITE')
                Write (30,'(A)') 'TITLE="Check random edge"'
                Write (30,*) 'VARIABLES="x", "y", "corner", "edge_length"'
                Write (30,'(A,I5,A,I5, A)') 'ZONE N =', adjn , ', E =', adjn, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
                do j = 1, adjn
                    write (30,*) adj_coord(:,j), temp_corner(j), temp_edge_length(j)
                enddo
                do j = 1, adjn
                    if (j == adjn) then
                        write (30,*) j, 1
                    else
                        write (30,*) j, j+1
                    endif
                enddo
                close(30)
            endif

            ! translate data by x=y graph line
            allocate (prop_edge(2,adjn), prop_patch(adjn), prop_corner(adjn))
            allocate (prop_coord(2,adjn+1), prop_edge_length(adjn))
            num = 0
            do j = adjn, 1, -1
                num = num + 1
                prop_coord(1,num) = adj_coord(2,j)
                prop_coord(2,num) = adj_coord(1,j)
                prop_corner(num) = temp_corner(j)
            enddo
            num = 0
            do j = adjn-1, 1, -1
                num = num + 1
                prop_edge_length(num) = temp_edge_length(j)
            enddo
            prop_edge_length(adjn) = temp_edge_length(adjn)

            if (remesh_type == 2) then
                allocate (adj(adjn))
                do j = 1, adjn
                    adj(j) = j
                enddo
                ! find the contact points
                call set_skin(2, 1, sp, target_adjn)
                allocate (target_adj(target_adjn), target_cp(target_adjn))
                allocate (cont_node(adjn))
                target_adj = region(2)%skin_nodes(sp:sp+target_adjn-1)
                sp = inte_set%cont_node(2);  ep = inte_set%cont_node(4)
                call find_contact_point(sp, ep, target_adjn, target_adj, 2, target_cn, target_cp)
    
                cont_node = 0
                do j = 1, adjn
                    cp = prop_coord(:,adj(j))
                    do k = 1, target_cn-1
                        call calc_len(region(2)%node_coord(:,target_cp(k)), region(2)%node_coord(:,target_cp(k+1)), dis)
                        call near_point(region(2)%node_coord(:,target_cp(k)), region(2)%node_coord(:,target_cp(k+1)), cp, p, dis*0.01, cont_check)
                        if (cont_check) then
                            cont_node(j) = 1
                            exit
                        endif
                    enddo
                enddo
                deallocate(target_adj, target_cp)
            
                call odd_to_even_node(adjn, prop_coord, cont_node, odd_check)
                deallocate (cont_node, adj)
            
                allocate (adj(adjn))
                do j = 1, adjn
                    adj(j) = j
                enddo
                nn = adjn
                call adjust_corner(1, i, nn, prop_coord(:,1:nn), adjn, adj, num_patch, 2)
                deallocate(adj)
            endif
            if (allocated(surface_set(i)%coord)) deallocate(surface_set(i)%coord)
            allocate(surface_set(i)%coord(2,adjn))
            surface_set(i)%nn = adjn
            surface_set(i)%coord(:,1:adjn) = prop_coord(:,1:adjn)
            if (output_flag == 1) then
                fn = './output/solid/remesh/check_rearrange0000000_2.plt'
                write (fn(38:44), '(I7.7)' ) current_step
                open(unit=30, file=fn, status='REPLACE', action='WRITE')
                Write (30,'(A)') 'TITLE="Check random edge"'
                Write (30,*) 'VARIABLES="x", "y"'
                Write (30,'(A,I5,A,I5, A)') 'ZONE N =', adjn , ', E =', adjn, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
                do j = 1, surface_set(i)%nn
                    write (30,*) surface_set(i)%coord(:,j)
                enddo
                do j = 1, adjn
                    if (j == adjn) then
                        write (30,*) j, 1
                    else
                        write (30,*) j, j+1
                    endif
                enddo
                close(30)
            endif
            mesh_set(rn)%ave_d(i) = 0.0
            do j = 1, adjn
                if (i == adjn) then
                    call calc_len(surface_set(i)%coord(:,j), surface_set(i)%coord(:,1), dis)
                else
                    call calc_len(surface_set(i)%coord(:,j), surface_set(i)%coord(:,j+1), dis)
                endif
                mesh_set(rn)%ave_d(i) = mesh_set(rn)%ave_d(i) + dis
            enddo
            mesh_set(rn)%ave_d(i) = mesh_set(rn)%ave_d(i)/float(adjn)
            !
            if (remesh_type == 1) then
            ! find corner info
                tol_corn_group(i) = tol_cn + 1
                do j = 1, adjn
                    if (prop_corner(j) == 2) then
                        tol_cn = tol_cn + 1
                        tol_corn(tol_cn) = j
                        tol_coord(:,tol_cn) = prop_coord(:,j)
                    endif
                enddo
        
                ! find corner_edge info
                do j = tol_corn_group(i), tol_cn
                    tol_edge(1,j) = prop_edge_length(tol_corn(j))
                    if (j == tol_corn_group(i)) then
                        if (tol_corn(j) == 1) then
                            tol_edge(2,tol_cn) = prop_edge_length(adjn)
                        else
                            tol_edge(2,tol_cn) = prop_edge_length(tol_corn(j)-1)
                        endif
                    else
                        tol_edge(2,j-1) = prop_edge_length(tol_corn(j)-1)
                    endif
                enddo
            endif
        
            deallocate (adj_coord)
            deallocate (prop_edge, prop_patch, prop_corner)
            deallocate (prop_coord, prop_edge_length)
            !mesh_set(rn)%flag(i) = PROPEL_REMESH_FLAG(ng-i+1)
            if (PROPEL_REMESH_FLAG(num_patch-i+1) == 3) then
                mesh_set(rn)%flag(i) = 4
            else
                if (odd_check) then
                    mesh_set(rn)%flag(i) = 3
                else
                    mesh_set(rn)%flag(i) = 2
                endif
            endif
        endif
    enddo

    if (remesh_type == 1) then
        deallocate (mesh_set(rn)%corn, mesh_set(rn)%edge, mesh_set(rn)%corn_coord)
        mesh_set(rn)%cn = tol_cn
        allocate (mesh_set(rn)%corn(mesh_set(rn)%cn))
        allocate (mesh_set(rn)%edge(2,mesh_set(rn)%cn))
        allocate (mesh_set(rn)%corn_coord(2,mesh_set(rn)%cn))
        mesh_set(rn)%corn = tol_corn(1:tol_cn)
        mesh_set(rn)%edge = tol_edge(:,1:tol_cn)
        mesh_set(rn)%corn_coord = tol_coord(:,1:tol_cn)
        mesh_set(rn)%corn_group = tol_corn_group
    endif
elseif (Dim == 3) then
    num_patch = sub_mesh_region   ! temporory
    !if (allocated(surface_set(num_patch)%p2)) deallocate(surface_set(num_patch)%p2)
    cur_ng = 0 
    do i = 1, num_patch
        ! update corner data
        sp = mesh_set(1)%sub_num(i,3)
        lp = mesh_set(1)%sub_num(i+1,3)
        prop_bn = lp-sp
        
        allocate (prop_corner(prop_bn))
        prop_corner = 0
        num = 0
        do j = sp, lp-1
            if (PROPEL_CORNER(j) == 3 .or. PROPEL_CORNER(j) == 6) then
                num = num + 1
                !prop_corner(num) = j - (sp-1)
                prop_corner(num) = j
            endif
        enddo
            
        ! temporory
        open (unit=37, file='tsm_edge.dat', status='OLD', action='READ')
        read (37,*) num
        read (37,*) prop_corner(1:num)
        ! temporory
        
        write (*,*) 'PROPEL_REMESH_FLAG:', PROPEL_REMESH_FLAG(i)
        write (*,*) 'prop_corner_num:', num
        write (*,*) prop_corner(1:num)
        
        !mesh_set(rn)%flag(i) = PROPEL_REMESH_FLAG(i)
        cur_ng = cur_ng + 1
        if (PROPEL_REMESH_FLAG(i) /= 0) then
            if (PROPEL_REMESH_FLAG(i) == 3) then
                mesh_set(rn)%flag(i) = 3
                cur_ng = cur_ng - 1
            else
                if (allocated(surface_set(i)%tri_p2)) deallocate(surface_set(i)%tri_p2)
                allocate(surface_set(i)%tri_p2(num))
                surface_set(i)%tri_p2 = 0
                surface_set(i)%tri_np2 = num
                do j = 1, num
                    surface_set(i)%tri_p2(j) = prop_corner(j) - (sp-1)
                enddo
        
                ! update ridge data
                if (allocated(surface_set(i)%ridge_num)) then
                    deallocate(surface_set(i)%ridge_num, surface_set(i)%ridge, surface_set(i)%ridge_conn)
                endif
                write (*,*) 'PROPEL_RIDGE_GROUP_NUM:', PROPEL_RIDGE_GROUP_NUM
                allocate(surface_set(i)%ridge_num(PROPEL_RIDGE_GROUP_NUM))
                allocate(surface_set(i)%ridge(500,PROPEL_RIDGE_GROUP_NUM))
                allocate(surface_set(i)%ridge_conn(PROPEL_RIDGE_GROUP_NUM))
                
                !surface_set(i)%ridge_group_num = 0
                !num = 0
                !do j = 1, propel_ridge_group_num
                !    write (*,*) j, 'ridge_num:', propel_ridge_num(j)
                !    write (*,*) j, 'ridge_first_node:', propel_ridge(j,1)
                !    if (sp <= propel_ridge(j,1) .and. propel_ridge(j,1) < lp) then
                !        num = num + 1
                !        surface_set(i)%ridge_num(num) = propel_ridge_num(j)
                !        surface_set(i)%ridge_conn(num) = j
                !        do k = 1, propel_ridge_num(j)
                !            surface_set(i)%ridge(k,num) = propel_ridge(j,k) - (sp-1)
                !        enddo
                !        !write (*,*) 'ridge_num:', propel_ridge_num(j)
                !        write (*,*) surface_set(i)%ridge(1:propel_ridge_num(j),num)
                !    endif
                !enddo
                !surface_set(i)%ridge_group_num = num
                !write (*,*) 'surface_set(i)%ridge_group_num:', num
            
                ! temporory
                read (37,*) surface_set(i)%ridge_group_num
                num = surface_set(i)%ridge_group_num
                read (37,*) surface_set(i)%ridge_num(1:num)
                do j = 1, num
                    n = MOD(surface_set(i)%ridge_num(j), 10)
		            count = 1
		            do k = 1, surface_set(i)%ridge_num(j)/10                                 ! each line: 10 input data
			            read(37,*) (surface_set(i)%ridge(q,j), q=count, count+9)   ! limitation: field(50) -> 5 input lines
			            count = count + 10
		            end do
		            if (n > 0) then
		              read(37,*) (surface_set(i)%ridge(k,j), k=count, count+n-1)  ! read last line if less than 10 data exist there
		            endif
                enddo
                close(37)
                do j = 1, num
                    surface_set(i)%ridge_conn(j) = j
                enddo
                write (*,'(100I8,1X)') surface_set(i)%ridge_num(1:num)
                do j = 1, num
                    write (*,'(/,(10I9))') surface_set(i)%ridge(1:surface_set(i)%ridge_num(i),i)
                enddo
              !  ! temporory          
              !  
                call quad2tri(current_step, 1, i)
                mesh_set(rn)%flag(i) = 2
            endif
        else
            Hexa(cur_ng)%np2 = num
            if (allocated(Hexa(cur_ng)%bp2_coord)) deallocate(Hexa(cur_ng)%bp2_coord)
            allocate (Hexa(cur_ng)%bp2_coord(3,num))
            do j = 1, num
                Hexa(cur_ng)%bp2_coord(:,j) = PROPEL_POINT(:,prop_corner(j))
            enddo
            
            if (allocated(Hexa(cur_ng)%bedn)) then
                deallocate(Hexa(cur_ng)%bedn, Hexa(cur_ng)%control_edge_coord, Hexa(cur_ng)%b_edge_conn)
            endif
            write (*,*) 'PROPEL_RIDGE_GROUP_NUM:', PROPEL_RIDGE_GROUP_NUM
            allocate(Hexa(cur_ng)%bedn(PROPEL_RIDGE_GROUP_NUM))
            allocate(Hexa(cur_ng)%control_edge_coord(3,200,PROPEL_RIDGE_GROUP_NUM))
            allocate(Hexa(cur_ng)%b_edge_conn(PROPEL_RIDGE_GROUP_NUM))
            
            Hexa(cur_ng)%control_edge_coord = 0.0
            num = 0
            do j = 1, PROPEL_RIDGE_GROUP_NUM
                write (*,*) j, 'ridge_num:', PROPEL_RIDGE_NUM(j)
                write (*,*) j, 'ridge_first_node:', PROPEL_RIDGE(j,1)
                if (sp <= PROPEL_RIDGE(j,1) .and. PROPEL_RIDGE(j,1) < lp) then
                    num = num + 1
                    Hexa(cur_ng)%bedn(num) = PROPEL_RIDGE_NUM(j)
                    Hexa(cur_ng)%b_edge_conn(num) = j
                    do k = 1, PROPEL_RIDGE_NUM(j)
                        cnum = PROPEL_RIDGE(j,k)
                        Hexa(cur_ng)%control_edge_coord(:,k,num) = PROPEL_POINT(:,cnum)
                    enddo
                endif
            enddo
            Hexa(cur_ng)%ben = num
              
            ! temporory
          !  read (37,*) Hexa(cur_ng)%ben
          !  num = Hexa(cur_ng)%ben
          !  read (37,*) Hexa(cur_ng)%bedn(1:num)
          !  do j = 1, num
          !      n = MOD(Hexa(cur_ng)%bedn(j), 10)
		        !count = 1
		        !do k = 1, Hexa(cur_ng)%bedn(j)/10                                 ! each line: 10 input data
			       ! read(37,*) temp_num   ! limitation: field(50) -> 5 input lines
          !          do q = 1, 10
          !              Hexa(cur_ng)%control_edge_coord(:,count-1+q,j) = PROPEL_POINT(:,temp_num(q))
          !          enddo
			       ! count = count + 10
          !      enddo
		        !if (n > 0) then
		        !    read(37,*) temp_num(1:n)  ! read last line if less than 10 data exist there
          !          do q = 1, n
          !          Hexa(cur_ng)%control_edge_coord(:,count-1+q,j) = PROPEL_POINT(:,temp_num(q))
          !          enddo
		        !endif
          !  enddo
          !  close(37)
          !  do j = 1, num
          !      Hexa(cur_ng)%b_edge_conn(j) = j
          !  enddo
          !  write (*,'(100I8,1X)') Hexa(cur_ng)%bedn(1:num)
            ! temporory
              
            write (*,*) ' Hexa(i)%ben:', num
            mesh_set(rn)%flag(i) = 0
        endif
        deallocate (prop_corner)
    enddo
endif
!if (output_flag == 1) then
!    fn = './output/solid/remesh/check_rearrange0000000_3.plt'
!    write (fn(38:44), '(I7.7)' ) current_step
!    open (Unit=30, File=fn, STATUS='replace', ACTION='write')
!    do i = 1, ng
!        if (mesh_set(rn)%flag(i) == 1) then
!            cn1 = mesh_set(rn)%corn_group(i)
!            if ( i == ng ) then
!                cn2 = mesh_set(rn)%cn
!            else
!                cn2 = mesh_set(rn)%corn_group(i+1) - 1
!            endif
!            do j = cn1, cn2
!                write (30,*) surface_set(i)%coord(:,mesh_set(rn)%corn(j))
!            enddo
!            exit
!        endif
!    enddo
!    close(30)
!endif

end subroutine rearrange_edge

          
subroutine odd_to_even_node(adjn, prop_coord, prop_cont_node, odd_check)

implicit none 

integer, intent(inout) :: adjn, prop_cont_node(adjn)
real,intent(inout) :: prop_coord(2,adjn+1)
logical, intent(inout) :: odd_check

integer :: i, j, num, sp, lp, flag, memory_j
integer :: clu_group_number, temp_clu, temp_clu_num, clu_num
integer :: this(3), clu_group(adjn), cont_node(adjn)
real :: dis, ang, cri_ang, max_dis, tol_ave_dis, clu_ave_dis, remain_dis, temp_dis
real :: clu_dis(adjn), coord(2,adjn+1), temp_coord(2,adjn), vec(2), temp_cp(2)

cri_ang = 1.0
! find first_node
cont_node = 0
do i = 1, adjn
    this = (/ i-1, i, i+1 /)
    if (i == 1) this(1) = adjn
    if (i == adjn) this(3) = 1
    
    call calc_angle(prop_coord(:,this(1)), prop_coord(:,this(2)), prop_coord(:,this(3)), ang)
    call calc_len(prop_coord(:,this(2)), prop_coord(:,this(3)), dis)
    if (abs(ang-180.0) > cri_ang) then
        clu_dis(1) = dis
        num = 0
        do j = i, adjn
            num = num + 1
            coord(:,num) = prop_coord(:,j)
            cont_node(num) = prop_cont_node(j)
        enddo
        do j = 1, i-1
            num = num + 1
            coord(:,num) = prop_coord(:,j)
            cont_node(num) = prop_cont_node(j)
        enddo
        exit
    endif
enddo

if (2*(adjn/2) /= adjn ) then
    odd_check = .true.
    write (*,*) 'number of propellant node is odd number: ', adjn
    ! set cluster
    clu_group_number = 1
    clu_group(clu_group_number) = 1
    do i = 2, adjn
        this = (/ i-1, i, i+1 /)
        if (i == adjn) this(3) = 1
    
        call calc_angle(coord(:,this(1)), coord(:,this(2)), coord(:,this(3)), ang)
        call calc_len(coord(:,this(2)), coord(:,this(3)), dis)
    
        if (abs(ang-180.0) > cri_ang) then
            clu_group_number = clu_group_number + 1
            clu_dis(clu_group_number) = dis
            clu_group(clu_group_number) = i
        else
            clu_dis(clu_group_number) = clu_dis(clu_group_number) + dis
        endif
    enddo
    clu_group(clu_group_number+1) = adjn+1
    max_dis = 0.0
    do i = 1, clu_group_number
        sp = clu_group(i)
        lp = clu_group(i+1)-1
        num = 0
        do j = sp, lp
            num = num + cont_node(j)
        enddo
        !write (*,*) 'clu_num:', i, num
        if (num < 2 .and. max_dis < clu_dis(i)) then
            temp_clu = i
            max_dis = clu_dis(i)
        endif
    enddo
    !write (*,*) 
    !write (*,*) 'temp_clu:', temp_clu
    tol_ave_dis = sum(clu_dis(1:clu_group_number))/float(adjn)
    sp = clu_group(temp_clu)
    lp = clu_group(temp_clu+1)-1
    clu_num = lp - sp + 1
    clu_ave_dis = clu_dis(temp_clu)/float(clu_num)
    if (tol_ave_dis > clu_ave_dis) then
        flag = -1
        do i = lp, adjn-1
            coord(:,i) = coord(:,i+1)
        enddo
    else
        flag = 1
        do i = adjn, lp+1, -1
            coord(:,i+1) = coord(:,i)
        enddo
    endif
    !do i = 1, clu_group_number
    !    write (*,*) 'clu_number:', i
    !    write (*,*) 'sp, lp, clu_num:', clu_group(i), clu_group(i+1)-1, clu_group(i+1)-clu_group(i)
    !    write (*,*) 'clu_dis:', clu_dis(i)
    !    write (*,*) 
    !enddo
    !write (*,*) 'total_clu_number:', clu_group_number
    !write (*,*) 'temp_clu:', temp_clu
    !write (*,*) 'sp, lp, clu_num:', sp, lp, clu_num
    !write (*,*) 'clu_flag:', flag
    
    do i = 1, clu_num
        temp_coord(:,i) = coord(:,sp-1+i)
    enddo
    if (temp_clu == clu_group_number) then
        temp_coord(:,clu_num+1) = coord(:,1)
    else
        temp_coord(:,clu_num+1) = coord(:,lp+1)
    endif
    temp_clu_num = clu_num + 1
    
    clu_num = clu_num + flag
    clu_ave_dis = clu_dis(temp_clu)/float(clu_num)
    temp_cp = temp_coord(:,1)
    memory_j = 2
    do i = sp+1, lp+flag
        remain_dis = clu_ave_dis
        do j = memory_j, temp_clu_num
            call calc_len(temp_cp, temp_coord(:,j), temp_dis)
            if (remain_dis < temp_dis) then
                vec = (temp_coord(:,j)-temp_cp(:))/temp_dis
                coord(:,i) = temp_cp(:)+vec(:)*remain_dis
                temp_cp = coord(:,i)
                memory_j = j
                exit
            else
                temp_cp = temp_coord(:,j)
                remain_dis = remain_dis - temp_dis
            endif
        enddo  
    enddo
    
    adjn = adjn + flag
    prop_coord(:,1:adjn) = coord(:,1:adjn)
else
    prop_coord(:,1:adjn) = coord
endif

end subroutine odd_to_even_node

! surface -> solid function end
!************************************
