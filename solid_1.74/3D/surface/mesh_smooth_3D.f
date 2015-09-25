subroutine ale_region_smoothing(nn, coord)

use surface_domain
implicit none

integer, intent(in) :: nn
real, intent(in) :: coord(3,nn)

integer :: i, j, k, q, iter, ng, num_group, case_bn, prop_bn, case_nn, case_ne, num, count, prop_bfn
integer :: p2_csn, csn(3), cp(3), i_seq(3)
real :: tol_dis, dis(2), vec(3), pl(4), temp_pl(4), p(3,3), op(3), yz_coord(2,3)
real :: ang, mean_dis
logical :: check
integer, allocatable :: case_elem(:,:)
real, allocatable :: case_coord(:,:), cs_nodes(:,:,:)
real, allocatable :: seep_vec(:,:,:), seep_dis(:,:,:), seep_vec2(:,:), seep_dis2(:,:)

!num_group = region(1)%num_groups
!case_nn = region(2)%num_nodes
!case_ne = region(2)%num_elements
!allocate(case_elem(8,case_ne), case_coord(3,case_nn))
!case_coord = region(2)%node_coord
!case_elem = region(2)%element
num_group = sub_mesh_region
case_nn = mesh_set(2)%nn
case_ne = mesh_set(2)%ne
allocate(case_elem(8,case_ne), case_coord(3,case_nn))
case_coord = mesh_set(2)%node_coord
case_elem = mesh_set(2)%element

do ng = 1, num_group
    csn(1) = move_set(ng)%csn(1,1)
    p2_csn = move_set(ng)%p2_csn - 2
    allocate (seep_vec2(3,move_set(ng)%p2cn), seep_dis2(p2_csn,move_set(ng)%p2cn))

    ! smoothing ale region nodes at matching nodes
    
    ! calculation seep vec
    yz_coord(:,2) = (/ 0.0, 0.0 /)
    yz_coord(:,3) = (/ 1.0, 0.0 /)
    do i = 1, move_set(ng)%p2cn
        yz_coord(:,1) = case_coord(2:3,move_set(ng)%prop2case(2,i))
        call calc_angle(yz_coord(:,1), yz_coord(:,2), yz_coord(:,3), ang)
        seep_vec2(:,i) = (/ 0.0, cos(ang*pi/180.0), sin(ang*pi/180.0) /)
    enddo

    do j = 1, move_set(ng)%p2cn
        do i = 2, move_set(ng)%p2c_nn(j)-1
            cp(1:2) = (/ move_set(ng)%prop2case(i,j), move_set(ng)%prop2case(i+1,j) /)
            call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
            seep_dis2(i-1,j) = dis(1)
        enddo
    enddo

    ! move ale region nodes matching propellant nodes
    do i = 1, move_set(ng)%p2cn
        cp(1:2) = move_set(ng)%prop2case(1:2,i)
        case_coord(:,cp(2)) = coord(:,cp(1))
    enddo

    do j = 1, move_set(ng)%p2cn
        do i = 2, move_set(ng)%p2c_nn(j)-1
            cp(1:2) = (/ move_set(ng)%prop2case(i,j), move_set(ng)%prop2case(i+1,j) /)
            vec = seep_vec2(:,j)
            dis(1) = seep_dis2(i-1,j)
            case_coord(:,cp(2)) = case_coord(:,cp(1)) + vec*dis(1)
        enddo
    enddo
    deallocate (seep_dis2, seep_vec2)

    do iter = 1, 2
        csn = move_set(ng)%csn(:,iter)
        ! calc seep vector
        allocate (seep_vec(3,csn(2),csn(3)-1), seep_dis(csn(1)-1,csn(2),csn(3)-1))
        allocate (cs_nodes(csn(1),csn(2),csn(3)))
        if (iter == 1) then
            cs_nodes = move_set(ng)%cs_nodes1
        else
            cs_nodes = move_set(ng)%cs_nodes2
        endif
    
        do i = 1, csn(2)
            do j = 1, csn(3)-1
                cp(1:2) = (/ cs_nodes(1,i,j), cs_nodes(2,i,j) /)
                call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
                seep_vec(:,i,j) = (case_coord(:,cp(2))-case_coord(:,cp(1)))/dis(1)
            enddo
        enddo

        ! calc seep distance
        do i = 1, csn(1)-1
            do j = 1, csn(2)
                do k = 1, csn(3)-1
                    cp(1:2) = (/ cs_nodes(i,j,k), cs_nodes(i+1,j,k) /)
                    call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
                    seep_dis(i,j,k) = dis(1)
                enddo
            enddo
        enddo
    
        ! smoothing ale region nodes excluding matching nodes
        do i = 1, csn(2)
            num = 0
            tol_dis = 0.0
            do j = 1, csn(3)-1
                cp(1:2) = (/ cs_nodes(1,i,j), cs_nodes(1,i,j+1) /)
                call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
                num = num + 1
                tol_dis = tol_dis + dis(1)
            enddo
            mean_dis = tol_dis/float(num)
            do j = 1, num-1
                cp = (/ cs_nodes(1,i,j), cs_nodes(1,i,j+1), cs_nodes(1,i,j+2) /)
                call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
                if (dis(1) <= mean_dis) then
                    call calc_length_3D(case_coord(:,cp(2)), case_coord(:,cp(3)), dis(2))
                    vec = case_coord(:,cp(3)) - case_coord(:,cp(2))
                    vec = vec/dis(2)
                else
                    vec = case_coord(:,cp(1)) - case_coord(:,cp(2))
                    vec = vec/dis(1)
                endif
                vec = vec * abs(mean_dis-dis(1))
                case_coord(:,cp(2)) = case_coord(:,cp(2)) + vec
            enddo
        enddo
    
        ! move inner_nodes and surface nodes on the other side of matching surface
        do i = 1, csn(1)-1
            do j = 1, csn(2)
                do k = 1, csn(3)-1
                    cp(1:2) = (/ cs_nodes(i,j,k), cs_nodes(i+1,j,k) /)
                    vec = seep_vec(:,j,k)
                    dis(1) = seep_dis(i,j,k)
                    case_coord(:,cp(2)) = case_coord(:,cp(1)) + vec*dis(1)
                enddo
            enddo
        enddo
        deallocate (seep_vec, seep_dis, cs_nodes)
    enddo
enddo

call update_mesh_info(3, 2, case_nn, case_ne, case_coord, case_elem, 1)
call write_tecplot_3D(case_nn, case_coord, case_ne, case_elem, 2, 1, 4)

deallocate (case_elem, case_coord)

end subroutine ale_region_smoothing


subroutine ale_region_smoothing2(nn, coord, flag)

use surface_domain
implicit none

integer, intent(in) :: nn, flag
real, intent(in) :: coord(3,nn)

integer :: i, j, k, q, iter, ng, num_group, rn_nn, case_bn, prop_bn, case_nn, case_ne, num, count, prop_bfn, p2_csn
integer :: csn(3), cp(3), i_seq(3)
real :: tol_dis, dis(2), mean_dis, vec(3), pl(4), temp_pl(4), p(3,5), op(3), ang, t
logical :: check, parallel
integer, allocatable :: case_elem(:,:), prop_b_face(:,:), prop_p2c_node(:), pre_face(:,:)
real, allocatable :: case_coord(:,:), cs_nodes(:,:,:), rn_coord(:,:)
real, allocatable :: seep_vec(:,:,:), seep_dis(:,:,:), seep_vec2(:,:), seep_dis2(:,:)

!num_group = region(1)%num_groups
!case_nn = region(2)%num_nodes
!case_ne = region(2)%num_elements
!allocate(case_elem(8,case_ne), case_coord(3,case_nn))
!case_coord = region(2)%node_coord
!case_elem = region(2)%element
num_group = sub_mesh_region

case_nn = mesh_set(2)%nn
case_ne = mesh_set(2)%ne
allocate(case_elem(8,case_ne), case_coord(3,case_nn))
case_coord = mesh_set(2)%node_coord
case_elem = mesh_set(2)%element
prop_bfn = mesh_set(1)%bfn
allocate (prop_b_face(4,prop_bfn), prop_p2c_node(mesh_set(1)%nn))
prop_b_face = mesh_set(1)%bound_face
op = (/ 0.0, 0.0, 0.0 /)

prop_p2c_node = 0
do ng = 1, num_group
    do i = 1, move_set(ng)%p2cn
        prop_p2c_node(move_set(ng)%prop2case(1,i)) = 1
    enddo
enddo
p(:,4) = (/ 10.0, 0.0, 0.0 /)
do ng = 1, num_group
    csn(1) = move_set(ng)%csn(1,1)
    p2_csn = move_set(ng)%p2_csn - 2
    allocate (seep_vec2(3,move_set(ng)%p2cn), seep_dis2(p2_csn,move_set(ng)%p2cn))
    seep_dis2 = 0.0
    
    ! smoothing ale region nodes at matching nodes
    do i = 1, move_set(ng)%p2cn
        num = move_set(ng)%prop2case(1,i)
        do j = 1, mesh_set(1)%bn
            if (mesh_set(1)%bound_node(j) == num) then
                num = j
                exit
            endif
        enddo
        p(:,2) = mesh_set(1)%node_coord(:,mesh_set(1)%bound_node(num))
        !p(:,2) = mesh_set(1)%node_coord(:,num)
        pl = 0.0
        do j = 1, prop_bfn
            do k = 1, 4
                if (prop_b_face(k,j) == num) then
                    count = 0
                    do q = 1, 4
                        count = count + prop_p2c_node(mesh_set(1)%bound_node(prop_b_face(q,j)))
                    enddo
                    if (count == 4) then
                        if (k == 1) then
                            p(:,1) = mesh_set(1)%node_coord(:,mesh_set(1)%bound_node(prop_b_face(4,j)))
                            p(:,3) = mesh_set(1)%node_coord(:,mesh_set(1)%bound_node(prop_b_face(2,j)))
                        elseif (k == 4) then
                            p(:,1) = mesh_set(1)%node_coord(:,mesh_set(1)%bound_node(prop_b_face(3,j)))
                            p(:,3) = mesh_set(1)%node_coord(:,mesh_set(1)%bound_node(prop_b_face(1,j)))
                        else
                            p(:,1) = mesh_set(1)%node_coord(:,mesh_set(1)%bound_node(prop_b_face(k-1,j)))
                            p(:,3) = mesh_set(1)%node_coord(:,mesh_set(1)%bound_node(prop_b_face(k+1,j)))
                        endif
                        call calc_plane(p(:,1), p(:,2), p(:,3), temp_pl)
                        pl = pl - temp_pl
                        exit
                    endif
                endif
            enddo
        enddo
        
        !cp(1:2) = (/ move_set(ng)%prop2case(2,i), move_set(ng)%prop2case(3,i) /)
        !call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
        !seep_vec2(:,i) = (case_coord(:,cp(2))-case_coord(:,cp(1)))/dis(1)
        p(:,1) = p(:,2)
        p(1,1) = 0.0
        call calc_plane(p(:,1), op, p(:,4), temp_pl)
        p(:,1) = p(:,2) + pl(1:3)
        p(:,3) = p(:,1)+10.0*temp_pl(1:3)
        p(:,1) = p(:,1)-10.0*temp_pl(1:3)
        call cross_point_plane(temp_pl, p(:,1), p(:,3), p(:,5), t, parallel)
        p(:,2) = p(:,5) - p(:,2)
        dis(1) = sqrt(p(1,2)**2.0+p(2,2)**2.0+p(3,2)**2.0)
        seep_vec2(:,i) = p(:,2)/dis(1)
    enddo

    !call cross_point_plane(pl(:,1), p(:,1), p(:,2), cp, t, parallel)
    !do i = 2, csn(1)
    do j = 1, move_set(ng)%p2cn
        do i = 2, move_set(ng)%p2c_nn(j) - 1
            cp(1:2) = (/ move_set(ng)%prop2case(2,j), move_set(ng)%prop2case(i+1,j) /)
            call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
            vec = case_coord(:,cp(2)) - case_coord(:,cp(1))
            call calc_angle_cos(vec, op, seep_vec2(:,j), ang)
            seep_dis2(i-1,j) = cos(ang*pi/180.)*dis(1)
        enddo
    enddo

    ! move ale region nodes matching propellant nodes
    do i = 1, move_set(ng)%p2cn
        cp(1:2) = move_set(ng)%prop2case(1:2,i)
        case_coord(:,cp(2)) = coord(:,cp(1))
    enddo

    do j = 1, move_set(ng)%p2cn
        do i = 2, move_set(ng)%p2c_nn(j)-1
            cp(1:2) = (/ move_set(ng)%prop2case(2,j), move_set(ng)%prop2case(i+1,j) /)
            vec = seep_vec2(:,j)
            dis(1) = seep_dis2(i-1,j)
            case_coord(:,cp(2)) = case_coord(:,cp(1)) + vec*dis(1)
        enddo
    enddo
    if (flag == 1) then
        allocate (move_set(ng)%seep_dis(p2_csn,move_set(ng)%p2cn))
    endif
    move_set(ng)%seep_dis = seep_dis2
    
    !call check_match_face_info(ng, p2_csn, move_set(ng)%p2cn, move_set(ng)%prop2case(2,:), seep_dis2, move_set(ng)%p2c_nn)
    deallocate (seep_dis2, seep_vec2)

    do iter = 1, 2
        csn = move_set(ng)%csn(:,iter)
        ! calc seep vector
        allocate (seep_vec(3,csn(2),csn(3)-1), seep_dis(csn(1)-1,csn(2),csn(3)-1))
        allocate (cs_nodes(csn(1),csn(2),csn(3)))
        if (iter == 1) then
            cs_nodes = move_set(ng)%cs_nodes1
        else
            cs_nodes = move_set(ng)%cs_nodes2
        endif
    
        do i = 1, csn(2)
            do j = 1, csn(3)-1
                cp(1:2) = (/ cs_nodes(1,i,j), cs_nodes(2,i,j) /)
                call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
                seep_vec(:,i,j) = (case_coord(:,cp(2))-case_coord(:,cp(1)))/dis(1)
            enddo
        enddo

        ! calc seep distance
        do i = 1, csn(1)-1
            do j = 1, csn(2)
                do k = 1, csn(3)-1
                    cp(1:2) = (/ cs_nodes(i,j,k), cs_nodes(i+1,j,k) /)
                    call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
                    seep_dis(i,j,k) = dis(1)
                enddo
            enddo
        enddo
    
        ! smoothing ale region nodes excluding matching nodes
        do i = 1, csn(2)
            num = 0
            tol_dis = 0.0
            do j = 1, csn(3)-1
                cp(1:2) = (/ cs_nodes(1,i,j), cs_nodes(1,i,j+1) /)
                call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
                num = num + 1
                tol_dis = tol_dis + dis(1)
            enddo
            mean_dis = tol_dis/float(num)
            do j = 1, num-1
                cp = (/ cs_nodes(1,i,j), cs_nodes(1,i,j+1), cs_nodes(1,i,j+2) /)
                call calc_length_3D(case_coord(:,cp(1)), case_coord(:,cp(2)), dis(1))
                if (dis(1) <= mean_dis) then
                    call calc_length_3D(case_coord(:,cp(2)), case_coord(:,cp(3)), dis(2))
                    vec = case_coord(:,cp(3)) - case_coord(:,cp(2))
                    vec = vec/dis(2)
                else
                    vec = case_coord(:,cp(1)) - case_coord(:,cp(2))
                    vec = vec/dis(1)
                endif
                vec = vec * abs(mean_dis-dis(1))
                case_coord(:,cp(2)) = case_coord(:,cp(2)) + vec
            enddo
        enddo
    
        ! move inner_nodes and surface nodes on the other side of matching surface
        do i = 1, csn(1)-1
            do j = 1, csn(2)
                do k = 1, csn(3)-1
                    cp(1:2) = (/ cs_nodes(i,j,k), cs_nodes(i+1,j,k) /)
                    vec = seep_vec(:,j,k)
                    dis(1) = seep_dis(i,j,k)
                    case_coord(:,cp(2)) = case_coord(:,cp(1)) + vec*dis(1)
                enddo
            enddo
        enddo
        deallocate (seep_vec, seep_dis, cs_nodes)
    enddo
enddo

if (flag == 1) then
    rn_nn = region(1)%num_nodes
    allocate (rn_coord(3,rn_nn))
    rn_coord = region(1)%node_coord
    do i = 1, case_nn
        rn_coord(:,mesh_set(2)%surf2ori(i)) = case_coord(:,i)
    enddo
    call update_coord(1, 3, rn_nn, rn_coord )
    deallocate (rn_coord)
endif
call update_mesh_info(3, 2, case_nn, case_ne, case_coord, case_elem, 1)
call write_tecplot_3D(case_nn, case_coord, case_ne, case_elem, 2, 1, 4)

deallocate (prop_b_face, prop_p2c_node)
deallocate (case_elem, case_coord)

end subroutine ale_region_smoothing2


subroutine check_match_face_info(ng, p2_csn, p2cn, prop2case, seep_dis, p2c_nn)

use surface_domain
implicit none

integer, intent(in) :: ng, p2_csn, p2cn, p2c_nn(p2cn), prop2case(p2cn)
real, intent(in) :: seep_dis(p2_csn,p2cn)

integer :: i, j, k, conn, count, num, case_nn, temp_node(p2cn)
logical :: check
integer, allocatable :: temp_b_face(:,:)

allocate(temp_b_face(4,move_set(ng)%match_bfn))
temp_node = 0
do i = 1, move_set(ng)%match_bfn
    do j = 1, 4
        num = move_set(ng)%match_b_face(j,i)
        conn = 0
        do k = 1, p2cn 
            if (prop2case(k) == num) then
                conn = k
                exit
            endif
        enddo
        if (conn /= 0) then
            temp_node(conn) = num
            temp_b_face(j,i) = conn
        else
            write (*,*) 'Cannot find the node connected with match_b_face', num
        endif
    enddo
enddo

open (Unit=113, File='check_ale_match_face.plt', STATUS='replace', ACTION='write')
Write (113,'(A)') 'TITLE="Boundary face"'
Write (113,'(A,A)') 'VARIABLES="x", "y", "z", "node", "p2c_nn", "dis1", "dis2", "dis3", "dis4", "dis5",'
Write (113,'(A, I5,A,I5, A)') 'ZONE N =', p2cn, ', E =', move_set(ng)%match_bfn,  ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
do i = 1, p2cn
    write (113,*) mesh_set(2)%node_coord(:,temp_node(i)), temp_node(i), p2c_nn(i), seep_dis(1:5,i)
enddo
do i = 1, move_set(ng)%match_bfn
    write (113,*) temp_b_face(:,i)
enddo
close(113)
stop

end subroutine check_match_face_info