subroutine remesh_ale_region1(num_group, mesh_flag, temp_arr, new_nn, new_ne, new_coord, new_elem, mat_array_num, mat_array)

use remesh_domain
use surface_domain
implicit none

integer, intent(in) :: num_group, temp_arr, mat_array_num, mesh_flag(num_group)
integer, intent(inout) :: new_nn, new_ne, new_elem(8,temp_arr), mat_array(mat_array_num)
real,    intent(inout) :: new_coord(3,temp_arr)

integer :: iter, i, j, k, q, temp_k, nn, ng, ne, num, temp_nn, temp_ne, count, layer_elem, mat_num
integer :: pre_nn, pre_ne, p2c_ne, p2cn, prop_bn, prop_bfn, match_tri_ne, mat_layer_num, cur_ng
integer :: face(4), del_num(2), elem(8,2), mat_layer(3,3), div_num(2), csn(3), face_nn(4,6,2), oppo_face(4,6), sub_num(3)
integer, allocatable :: temp_elem(:,:), prop_b_node(:), prop_b_face(:,:), match_tri_elem(:,:)
integer, allocatable :: node_check(:), conn(:), prop_conn(:), del_node(:,:), side_nodes(:,:,:), mat_num_mesh(:)
integer, allocatable :: cs_nodes(:,:,:), temp_cs_nodes(:), elem_check(:)
real :: theta, p(3,7), dis(2), pl(4), t, area, yz_coord(2,3), vec(3)
real, allocatable :: temp_coord(:,:), match_tri_coord(:,:), seep_vec(:,:), seep_dis(:), cs_nodes_112(:,:,:)
logical :: check(2)
!real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

!num_group = sub_mesh_region
allocate(cs_nodes_112(3,2,num_group))

mat_layer_num = 3
mat_layer(:,1) = (/ 1, 3, 7 /)
mat_layer(:,2) = (/ 2, 6, 7 /)
mat_layer(:,3) = (/ 1, 2, 3 /)
mat_layer(:,3) = mat_layer(:,3) + num_group

oppo_face(:,1) = (/ 5, 6, 7, 8 /)
oppo_face(:,2) = (/ 1, 4, 3, 2 /)
oppo_face(:,3) = (/ 4, 8, 7, 3 /)
oppo_face(:,4) = (/ 1, 5, 8, 4 /)
oppo_face(:,5) = (/ 2, 6, 5, 1 /)
oppo_face(:,6) = (/ 2, 3, 7, 6 /)

nn = mesh_set(2)%nn
ne = mesh_set(2)%ne
allocate (node_check(nn), conn(nn))
allocate (temp_coord(3,nn*2), temp_elem(8,ne*2), mat_num_mesh(ne*2))
allocate (side_nodes(move_set(1)%csn(1,1)*move_set(1)%csn(2,1),2,num_group))
mat_num_mesh = 0

! Find nodes & elements need to keep
node_check = 1
conn = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 2, move_set(ng)%csn(1,1)+1
            do j = 1, move_set(ng)%p2cn
                node_check(move_set(ng)%prop2case(i,j)) = 0
            enddo
        enddo
    endif
enddo

!num = 0
!do i = 1, inte_set%csn(1,1)
!    do j = 1, inte_set%csn(2,1)
!        node_check(inte_set%cs_nodes1(i,j,inte_set%csn(3,1))) = 1
!        num = num + 1
!        side_nodes(num,1) = inte_set%cs_nodes1(i,j,inte_set%csn(3,1))
!    enddo
!enddo
!num = 0
!do i = 1, inte_set%csn(1,2)
!    do j = 1, inte_set%csn(2,2)
!        node_check(inte_set%cs_nodes2(i,j,inte_set%csn(3,2))) = 1
!        num = num + 1
!        side_nodes(num,2) = inte_set%cs_nodes2(i,j,inte_set%csn(3,2))
!    enddo
!enddo

!temp_nn = 0
!do i = 1, nn
!    if (node_check(i) == 1) then
!        temp_nn = temp_nn + 1
!        temp_coord(:,temp_nn) = mesh_set(2)%node_coord(:,i)
!        conn(i) = temp_nn
!    endif
!enddo
!do i = 1, 2
!    side_nodes(:,i) = conn(side_nodes(:,i))
!enddo

do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 2, move_set(ng)%csn(3,1)
            do j = 1, move_set(ng)%csn(1,1)
                do k = 1, move_set(ng)%csn(2,1)
                    node_check(move_set(ng)%cs_nodes1(j,k,i)) = 0
                enddo
            enddo
        enddo
        do i = 2, move_set(ng)%csn(3,2)
            do j = 1, move_set(ng)%csn(1,2)
                do k = 1, move_set(ng)%csn(2,2)
                    node_check(move_set(ng)%cs_nodes2(j,k,i)) = 0
                enddo
            enddo
        enddo
    endif
enddo

temp_nn = 0
do i = 1, nn
    if (node_check(i) == 1) then
        temp_nn = temp_nn + 1
        temp_coord(:,temp_nn) = mesh_set(2)%node_coord(:,i)
        conn(i) = temp_nn
    endif
enddo

temp_ne = 0
do i = 1, ne
    num = sum(node_check(mesh_set(2)%element(:,i)))
    if (num == 8) then
        temp_ne = temp_ne + 1
        temp_elem(:,temp_ne) = conn(mesh_set(2)%element(:,i))
        elem(:,1) = mesh_set(2)%surf2ori(mesh_set(2)%element(:,i))
        do j = 1, region(1)%num_elements
            elem(:,2) = region(1)%element(4:11,j)
            call find_count(8, elem(:,1), 8, elem(:,2), count)
            if (count == 8) then
                mat_num_mesh(temp_ne) = region(1)%element(3,j)
                exit
            endif
        enddo
    endif
enddo
deallocate (node_check, conn)

! make mesh at region(3) - smoothing region
temp_k = temp_ne
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do iter = 1, 2
            q = 0
            csn = move_set(ng)%csn(:,iter)
            allocate (cs_nodes(csn(1),csn(2),csn(3)))
            if (iter == 1) then
                cs_nodes = move_set(ng)%cs_nodes1
            else
                cs_nodes = move_set(ng)%cs_nodes2
            endif
            pre_nn = temp_nn
            pre_ne = temp_ne
            p(:,1) = mesh_set(2)%node_coord(:,cs_nodes(1,1,1))
            p(:,2) = mesh_set(2)%node_coord(:,cs_nodes(1,1,csn(3)))
            call calc_length_3D(p(:,1), p(:,2), dis(1))
            div_num(1) = int(dis(1)/remesh_dis)+2
            layer_elem = 0
            do i = 1, csn(1)
                if (i /= csn(1)) then
                    do j = 1, mat_layer_num
                        if (mat_layer(j,1) <= i .and. mat_layer(j,2) >= i) then
                            mat_num = mat_layer(j,3)
                            exit
                        endif
                    enddo
                endif
    
                count = temp_nn
                allocate (temp_cs_nodes(csn(2)))
                do j = 1, csn(2)
                    p(:,1) = mesh_set(2)%node_coord(:,cs_nodes(i,j,1))
                    p(:,2) = mesh_set(2)%node_coord(:,cs_nodes(i,j,csn(3)))
                    call find_match_node_num(p(:,1), pre_nn, temp_coord(:,1:pre_nn), temp_cs_nodes(j), remesh_dis*0.01, -1)
                    vec = (p(:,2)-p(:,1))/float(div_num(1))
                    do k = 1, div_num(1)
                        temp_nn = temp_nn + 1
                        temp_coord(:,temp_nn) = p(:,1) + k*vec
                        if (i == 1 .and. j == 1 .and. k == 1) cs_nodes_112(:,iter,ng) = temp_coord(:,temp_nn)
                        if (k == div_num(1)) then
                            q = q + 1
                            side_nodes(q,iter,ng) = temp_nn
                        endif
                    enddo
                enddo
        
                if (i /= csn(1)) then
                    do j = 2, csn(2)
                        do k = 1, div_num(1)
                            temp_ne = temp_ne + 1
                            mat_num_mesh(temp_ne) = mat_num
                            if (i == 1) layer_elem = layer_elem + 1
                            if (k == 1) then                    
                                temp_elem(1:4,temp_ne) = (/ temp_cs_nodes(j-1), count+k+(j-2)*div_num(1), &
                                                            count+k+(j-1)*div_num(1), temp_cs_nodes(j) /)
                            else
                                temp_elem(1:4,temp_ne) = (/ count+k-1+(j-2)*div_num(1), count+k+(j-2)*div_num(1), &
                                                            count+k+(j-1)*div_num(1), count+k-1+(j-1)*div_num(1) /)
                            endif
                        enddo
                    enddo
                endif
                if (i /= 1) then
                    if (i == csn(1)) then
                        num = 0
                        do j = 2, csn(2)
                            do k = 1, div_num(1)
                                num = num + 1
                                if (k == 1) then                    
                                    temp_elem(5:8,pre_ne+(i-2)*layer_elem+num) = (/ temp_cs_nodes(j-1), count+k+(j-2)*div_num(1), &
                                                                                    count+k+(j-1)*div_num(1), temp_cs_nodes(j) /)
                                else
                                    temp_elem(5:8,pre_ne+(i-2)*layer_elem+num) = (/ count+k-1+(j-2)*div_num(1), count+k+(j-2)*div_num(1), &
                                                                                    count+k+(j-1)*div_num(1), count+k-1+(j-1)*div_num(1) /)
                                endif
                            enddo
                        enddo
                    else
                        num = 0
                        do j = 2, csn(2)
                            do k = 1, div_num(1)
                                num = num + 1
                                temp_elem(5:8,pre_ne+(i-2)*layer_elem+num) = temp_elem(1:4,pre_ne+(i-1)*layer_elem+num)
                            enddo
                        enddo
                    endif
                endif
                deallocate(temp_cs_nodes)
            enddo
            deallocate(cs_nodes)
        enddo
    endif
enddo

! reaarange connectivity number
!pre_ne = temp_k
!allocate(elem_check(temp_ne))
!elem_check = 0
!elem_check(1:pre_ne) = 1
!do
!    do i = pre_ne+1, temp_ne
!        if (elem_check(i) == 0) then
!            elem(:,1) = temp_elem(:,i)
!            call make_face_nn(8, 4, elem(:,1), face_nn(:,:,1))
!            do j = temp_ne, 1, -1
!                if (elem_check(j) == 1) then
!                    elem(:,2) = temp_elem(:,j)
!                    call find_count(8, temp_elem(:,i), 8, temp_elem(:,j), count)
!                    if (count == 4) then
!                        call make_face_nn(8, 4, temp_elem(:,j), face_nn(:,:,2))
!                        !write (*,'(8I6)') temp_elem(:,i)
!                        !do k = 1, 6
!                        !    write (*,*) face_nn(:,k,1)
!                        !enddo
!                        !write (*,*) 
!                        !write (*,'(8I6)') temp_elem(:,j)
!                        !do k = 1, 6
!                        !    write (*,*) face_nn(:,k,2)
!                        !enddo
!                        do k = 1, 2
!                            do q = 1, 6
!                                if (k == 1) then
!                                    call find_count(4, face_nn(:,q,k), 8, elem(:,2), count)
!                                else
!                                    call find_count(4, face_nn(:,q,k), 8, elem(:,1), count)
!                                endif
!                                if (count == 4) then
!                                    face(k) = q
!                                    exit
!                                endif
!                            enddo
!                        enddo
!                        do k = 1, 4
!                            if (face_nn(1,face(1),1) == face_nn(k,face(2),2)) then
!                                check(1) = .false.
!                                if (k == 4) then
!                                    if (face_nn(2,face(1),1) == face_nn(1,face(2),2)) then
!                                        check(2) = .true.
!                                    endif
!                                else
!                                    if (face_nn(2,face(1),1) == face_nn(k+1,face(2),2)) then
!                                        check(1) = .true.
!                                    endif
!                                endif
!                                exit
!                            endif
!                        enddo
!                        !write (*,*) check(1)
!                        !write (*,*) face_nn(:,face(1),1)
!                        !write (*,*) face_nn(:,face(2),2)
!                        if (check(1)) then
!                            temp_elem(1,i) = face_nn(1,face(1),1)
!                            temp_elem(5,i) = elem(oppo_face(1,face(1)),1)
!                            do k = 4, 2, -1
!                                temp_elem(6-k,i) = face_nn(k,face(1),1)
!                                temp_elem(10-k,i) = elem(oppo_face(k,face(1)),1)
!                            enddo
!                            
!                        else
!                            temp_elem(1:4,i) = face_nn(:,face(1),1)
!                            do k = 1, 4
!                                temp_elem(k+4,i) = elem(oppo_face(k,face(1)),1)
!                            enddo
!                        endif
!                        elem_check(i) = 1
!                    endif
!                endif
!            enddo
!        endif        
!    enddo
!    if (sum(elem_check) == temp_ne) exit
!enddo

prop_bn = 0
prop_bfn = 0
num = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        num = num + 1
        prop_bn = prop_bn + Hexa(num)%bn
        prop_bfn = prop_bfn + Hexa(num)%bfn
    endif
enddo
allocate (prop_b_node(prop_bn), prop_b_face(4,prop_bfn), prop_conn(new_nn))
sub_num = 0
num = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        num = num + 1
        write (*,*) 'group number:', num
        prop_b_node(sub_num(1)+1:sub_num(1)+Hexa(num)%bn) = Hexa(num)%bound_node + sub_num(3)
        do j = 1, Hexa(num)%bfn
            prop_b_face(:,sub_num(2)+j) = Hexa(num)%bound_face(:,j) + sub_num(3)
        enddo
        sub_num(1) = sub_num(1) + Hexa(num)%bn
        sub_num(2) = sub_num(2) + Hexa(num)%bfn
        sub_num(3) = sub_num(3) + Hexa(num)%nn
    endif
enddo
do i = 1, prop_bn
    prop_conn(prop_b_node(i)) = i
enddo

match_tri_ne = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        match_tri_ne = match_tri_ne + move_set(ng)%match_bfn*2
    endif
enddo
allocate (match_tri_elem(3,match_tri_ne))
open (Unit=20, File='./output/solid/remesh/check_match_surf.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="check_match_surf"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6, A)') 'ZONE N=', mesh_set(2)%nn , ', E=', match_tri_ne, ', DATAPACKING = POINT, ZONETYPE = FETRIANGLE'
do i = 1, mesh_set(2)%nn
    write (20,*) mesh_set(2)%node_coord(:,i)
enddo
num = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 1, move_set(ng)%match_bfn
            face = move_set(ng)%match_b_face(:,i)
            p(:,1:4) = mesh_set(2)%node_coord(:,face)
            call calc_length_3D(p(:,1), p(:,3), dis(1))
            call calc_length_3D(p(:,2), p(:,4), dis(2))
            if (dis(1) < dis(2)) then
                match_tri_elem(:,num+1) = (/ face(1), face(2), face(3) /)
                match_tri_elem(:,num+2) = (/ face(3), face(4), face(1) /)
            else
                match_tri_elem(:,num+1) = (/ face(1), face(2), face(4) /)
                match_tri_elem(:,num+2) = (/ face(2), face(3), face(4) /)
            endif
            write (20,'(3(I7,2X))') match_tri_elem(:,num+1)
            write (20,'(3(I7,2X))') match_tri_elem(:,num+2)
            num = num + 2
        enddo
    endif
enddo
close(20)
    
allocate (node_check(prop_bn), conn(prop_bn))
node_check = 0
conn = 0
pre_nn = temp_nn
do i = 1, prop_bn
    p(:,1) = new_coord(:,prop_b_node(i))
    do j = 1, match_tri_ne
        p(:,2) = mesh_set(2)%node_coord(:,match_tri_elem(1,j))
        p(:,3) = mesh_set(2)%node_coord(:,match_tri_elem(2,j))
        p(:,4) = mesh_set(2)%node_coord(:,match_tri_elem(3,j))
        call calc_plane(p(:,2), p(:,3), p(:,4), pl)
        dis(1) = abs(pl(1)*p(1,1) + pl(2)*p(2,1) + pl(3)*p(3,1) + pl(4))
        check = .false.
        if (dis(1) < remesh_dis*0.05) then
            p(:,5) = p(:,1) + pl(1:3)*remesh_dis
            p(:,6) = p(:,1) - pl(1:3)*remesh_dis
            call cross_point_plane(pl, p(:,5), p(:,6), p(:,7), t, check(1))
            if (check(1) == .false.) then  ! 평면과 평행하지 않으면...
                call check_inner_area(p(:,2:4), p(:,7), check(2), area)
                if (check(2)) then
                    node_check(i) = 1
                    temp_nn = temp_nn + 1
                    temp_coord(:,temp_nn) = p(:,1)
                    conn(i) = temp_nn
                    exit
                endif
            endif
        endif
        if (check(2)) exit
    enddo
enddo

p2cn = sum(node_check)
p2c_ne = 0
pre_ne = temp_ne
do i = 1, prop_bfn
    face = prop_conn(prop_b_face(:,i))
    num = sum(node_check(face))
    if (num == 4) then
        temp_ne = temp_ne + 1
        p2c_ne = p2c_ne + 1
        temp_elem(1:4,temp_ne) = conn(face)
    endif
enddo

! calculation seep dis
allocate (seep_dis(move_set(1)%csn(1,1)-1))
p(:,1) = mesh_set(2)%node_coord(:,move_set(1)%cs_nodes1(1,1,1))
do i = 1, move_set(1)%csn(1,1)-1
    p(:,2) = mesh_set(2)%node_coord(:,move_set(1)%cs_nodes1(i+1,1,1))
    call calc_length_3D(p(:,1), p(:,2), seep_dis(i))
enddo

! calculation seep vec
allocate (seep_vec(3,p2cn))
yz_coord(:,2) = (/ 0.0, 0.0 /)
yz_coord(:,3) = (/ 1.0, 0.0 /)
do i = pre_nn+1, temp_nn
    yz_coord(:,1) = temp_coord(2:3,i)
    call calc_angle(yz_coord(:,1), yz_coord(:,2), yz_coord(:,3), theta)
    seep_vec(:,i-pre_nn) = (/ 0.0, cos(theta*pi/180.0), sin(theta*pi/180.0) /)
enddo

do i = 1, move_set(1)%csn(1,1)-1
    do j = 1, p2cn
        p(:,1) = temp_coord(:,pre_nn+j)
        vec = seep_vec(:,j)
        dis(1) = seep_dis(i)
        temp_nn = temp_nn + 1
        temp_coord(:,temp_nn) = p(:,1) + vec*dis(1)
    enddo
    do j = 1, mat_layer_num
        if (mat_layer(j,1) <= i .and. mat_layer(j,2) >= i) then
            num = mat_layer(j,3)
            exit
        endif
    enddo
    !write (*,*) i, 'mat_num:', num
    do j = 1, p2c_ne
        temp_elem(5:8,pre_ne+j) = temp_elem(1:4,pre_ne+j)+p2cn
        mat_num_mesh(pre_ne+j) = num
        if (i /= move_set(1)%csn(1,1)-1) then
            temp_ne = temp_ne + 1
            temp_elem(1:4,temp_ne) = temp_elem(5:8,pre_ne+j)
        endif
    enddo
    pre_ne = temp_ne - p2c_ne
enddo
deallocate (seep_dis, seep_vec)
open (Unit=20, File='./output/solid/remesh/check_prop_surf.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check prop surface mesh"'
Write (20,*) 'VARIABLES="x", "y", "z", "num"'
Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Prop_surf", N = ', prop_bn, ', E = ', prop_bfn, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
do i = 1, prop_bn
    write (20,*) new_coord(:,prop_b_node(i)), node_check(i)
enddo
do i = 1, prop_bfn
    write (20,*) prop_conn(prop_b_face(:,i))
enddo
close(20)

!open (Unit=20, File='./output/solid/remesh/check_case_mesh1.plt', STATUS='replace', ACTION='write')
!Write (20,'(A)') 'TITLE="Check case mesh"'
!Write (20,*) 'VARIABLES="x", "y", "z"'
!Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Case", N = ', temp_nn, ', E = ', temp_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
!do i = 1, temp_nn
!    write (20,*) temp_coord(:,i)
!enddo
!do i = 1, temp_ne
!    write (20,*) temp_elem(:,i)
!enddo
!close(20)

! summation with case nodes & case nodes
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 1, 2
            do j = 1, move_set(1)%csn(1,1)*move_set(1)%csn(2,1)
                p(:,1) = temp_coord(:,side_nodes(j,i,ng))
                dis(2) = 10e8
                del_num = 0
                do k = pre_nn+1, temp_nn
                    p(:,2) = temp_coord(:,k)
                    call calc_length_3D(p(:,1), p(:,2), dis(1))
                    if (dis(1) < remesh_dis*0.01) then
                        del_num = (/ k, side_nodes(j,i,ng) /)
                        exit
                    elseif (dis(2) > dis(1)) then
                        temp_k = k
                        dis(2) = dis(1)
                    endif
                enddo
                if (del_num(1) == 0) then
                    del_num = (/ temp_k, side_nodes(j,i,ng) /)
                endif
                do k = del_num(1), temp_nn-1
                    temp_coord(:,k) = temp_coord(:,k+1)
                enddo
                temp_nn = temp_nn - 1
                do k = 1, temp_ne
                    do q = 1, 8
                        if (temp_elem(q,k) > del_num(1)) then
                            temp_elem(q,k) = temp_elem(q,k) - 1
                        elseif (temp_elem(q,k) == del_num(1)) then
                            temp_elem(q,k) = del_num(2)
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif
enddo
deallocate (side_nodes)

open (Unit=20, File='./output/solid/remesh/check_case_mesh.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Case", N = ', temp_nn, ', E = ', temp_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do i = 1, temp_nn
    write (20,*) temp_coord(:,i)
enddo
do i = 1, temp_ne
    write (20,*) temp_elem(:,i)
enddo
close(20)

pre_nn = new_nn
pre_ne = new_ne
new_coord(:,new_nn+1:new_nn+temp_nn) = temp_coord(:,1:temp_nn)
new_elem(:,new_ne+1:new_ne+temp_ne) = temp_elem(:,1:temp_ne)+new_nn
new_nn = new_nn + temp_nn
new_ne = new_ne + temp_ne
deallocate (temp_coord, temp_elem)

! summation with propellant nodes & case nodes
do i = pre_nn, 1, -1
    check(1) = .FALSE.
    do j = 1, prop_bn
        if (prop_b_node(j) == i) then
            if (node_check(j) == 1) then
                check(1) = .TRUE.
            endif
            exit
        endif
    enddo
    
    if (check(1)) then
        p(:,1) = new_coord(:,i)
        dis(2) = 10e8
        del_num = 0
        do j = pre_nn+1, new_nn
            p(:,2) = new_coord(:,j)
            call calc_length_3D(p(:,1), p(:,2), dis(1))
            if (dis(1) < remesh_dis*0.01) then
                !del_num = (/ i, j-1 /)
                del_num = (/ i, j /)
                exit
            elseif (dis(2) > dis(1)) then
                temp_k = j
                dis(2) = dis(1)
            endif
        enddo
        if (del_num(1) == 0) then
            !del_num = (/ i, temp_k-1 /)
            del_num = (/ i, temp_k /)
        endif
        do j = del_num(2), new_nn-1
            new_coord(:,j) = new_coord(:,j+1)
        enddo
        new_nn = new_nn - 1
        do j = 1, new_ne
            do k = 1, 8
                if (new_elem(k,j) == del_num(2)) then
                    new_elem(k,j) = del_num(1)
                elseif (new_elem(k,j) > del_num(2)) then
                    new_elem(k,j) = new_elem(k,j) - 1
                endif
            enddo
        enddo
    endif
enddo
deallocate (node_check)

do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 1, 2
            do j = 1, move_set(ng)%csn(2,1)+1
                if (i == 1) then
                    if (j == move_set(ng)%csn(2,1)+1) then
                        !p(:,1) = mesh_set(2)%node_coord(:,inte_set%cs_nodes1(1,1,2))
                        p(:,1) = cs_nodes_112(:,i,ng)
                    else
                        p(:,1) = mesh_set(2)%node_coord(:,move_set(ng)%cs_nodes1(1,j,1))
                    endif
                else
                    if (j == move_set(ng)%csn(2,1)+1) then
                        !p(:,1) = mesh_set(2)%node_coord(:,inte_set%cs_nodes2(1,1,2))
                        p(:,1) = cs_nodes_112(:,i,ng)
                    else
                        p(:,1) = mesh_set(2)%node_coord(:,move_set(ng)%cs_nodes2(1,j,1))
                    endif
                endif
                dis(2) = 10e8
                temp_k = 0
                do k = 1, new_nn
                    p(:,2) = new_coord(:,k)
                    call calc_length_3D(p(:,1), p(:,2), dis(1))
                    if (dis(1) < remesh_dis*0.01) then
                        temp_k = k
                        exit
                    elseif (dis(2) > dis(1)) then
                        temp_k = k
                        dis(2) = dis(1)
                    endif
                enddo
                if (i == 1) then
                    if (j == move_set(ng)%csn(2,1)+1) then
                        move_set(ng)%cs_nodes1(1,1,2) = temp_k
                    else
                        move_set(ng)%cs_nodes1(1,j,1) = temp_k
                    endif
                else
                    if (j == move_set(ng)%csn(2,1)+1) then
                        move_set(ng)%cs_nodes2(1,1,2) = temp_k
                    else
                        move_set(ng)%cs_nodes2(1,j,1) = temp_k
                    endif
                endif
            enddo
        enddo
    endif
enddo
deallocate(cs_nodes_112)

! element rearrange according mat_num
allocate (temp_elem(8,new_ne))
temp_elem = new_elem
pre_ne = 0
cur_ng = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        cur_ng = cur_ng + 1
        mat_array(cur_ng) = Hexa(ng)%ne + pre_ne
        pre_ne = pre_ne + Hexa(ng)%ne
    endif
enddo
do i = 1, mat_layer_num
    num = 0
    do j = 1, temp_ne
        if (mat_num_mesh(j) == i+num_group) then
            num = num + 1
            new_elem(:,pre_ne+num) = temp_elem(:,j+mat_array(cur_ng))
        endif
    enddo
    pre_ne = pre_ne + num
    mat_array(i+cur_ng) = pre_ne
enddo
deallocate (temp_elem)
!write (*,*) inte_set%cs_nodes1(1,:,1), inte_set%cs_nodes1(1,1,2)
!write (*,*) inte_set%cs_nodes2(1,:,1), inte_set%cs_nodes2(1,1,2)

open (Unit=20, File='./output/solid/remesh/final_mesh.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
do i = 1, mat_array_num
    if (i == 1) then
        k = 1
    else
        k = mat_array(i-1)+1
    endif
    q = mat_array(i)
    num = q - k + 1
    Write (20,'(A,I2,A,I6,A,I6,A)') 'ZONE T = "mat_num', i, '", N = ', new_nn, ', E = ', num, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
    do j = 1, new_nn
        write (20,*) new_coord(:,j)
    enddo
    do j = k, q
        write (20,*) new_elem(:,j)
    enddo
enddo
close(20)

write (*,*) '>> Hexa_mesh_rearrange'
call Hexa_mesh_rearrange(new_nn, new_coord, new_ne, new_elem)

open (Unit=20, File='./output/solid/remesh/final_mesh2.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
do i = 1, mat_array_num
    if (i == 1) then
        k = 1
    else
        k = mat_array(i-1)+1
    endif
    q = mat_array(i)
    num = q - k + 1
    Write (20,'(A,I2,A,I6,A,I6,A)') 'ZONE T = "mat_num', i, '", N = ', new_nn, ', E = ', num, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
    do j = 1, new_nn
        write (20,*) new_coord(:,j)
    enddo
    do j = k, q
        write (20,*) new_elem(:,j)
    enddo
enddo
close(20)

end subroutine remesh_ale_region1
    
    
    
subroutine remesh_ale_region2(num_group, mesh_flag, temp_arr, new_nn, new_ne, new_coord, new_elem, mat_array_num, mat_array)

use remesh_domain
use surface_domain
implicit none

integer, intent(in) :: num_group, temp_arr, mat_array_num, mesh_flag(num_group)
integer, intent(inout) :: new_nn, new_ne, new_elem(8,temp_arr), mat_array(mat_array_num)
real,    intent(inout) :: new_coord(3,temp_arr)

integer :: iter, i, j, k, q, temp_k, nn, ng, ne, num, temp_nn, temp_ne, count, layer_elem, mat_num
integer :: pre_nn, pre_ne, pre_ne2, p2c_ne, p2cn, prop_bn, prop_bfn, match_tri_ne, mat_layer_num, cur_ng
integer :: face(4), del_num(2), elem(8,2), mat_layer(3,3), div_num(2), csn(3), face_nn(4,6,2), oppo_face(4,6), sub_num(3)
integer, allocatable :: temp_elem(:,:), prop_b_node(:), prop_b_face(:,:), match_tri_elem(:,:)
integer, allocatable :: node_check(:), conn(:), prop_conn(:), del_node(:,:), side_nodes(:,:,:), mat_num_mesh(:)
integer, allocatable :: cs_nodes(:,:,:), temp_cs_nodes(:), elem_check(:)
real :: theta, min_dis, p(3,7), dis(3), pl(4), temp_pl(4), t, area, yz_coord(2,3), vec(3), op(3)
real, allocatable :: temp_coord(:,:), match_tri_coord(:,:), seep_vec(:,:), seep_dis(:,:), cs_nodes_112(:,:,:)
logical :: check(2), parallel, near_check

integer :: ori_p2cn, min_p2c_nn, case_nn, p2_csn, temp_j, end_num, sum_p2c_node_num
integer, allocatable :: match_conn(:), p2c_nn(:), end_pos(:), elem_pos(:,:), sum_p2c_node(:,:)
real, allocatable :: temp_val(:,:)

!num_group = sub_mesh_region
allocate(cs_nodes_112(3,2,num_group))

op = (/ 0.0, 0.0, 0.0 /)
mat_layer_num = 3
mat_layer(:,1) = (/ 1, 3, 7 /)
mat_layer(:,2) = (/ 2, 6, 7 /)
mat_layer(:,3) = (/ 1, 2, 3 /)
mat_layer(:,3) = mat_layer(:,3) + num_group

oppo_face(:,1) = (/ 5, 6, 7, 8 /)
oppo_face(:,2) = (/ 1, 4, 3, 2 /)
oppo_face(:,3) = (/ 4, 8, 7, 3 /)
oppo_face(:,4) = (/ 1, 5, 8, 4 /)
oppo_face(:,5) = (/ 2, 6, 5, 1 /)
oppo_face(:,6) = (/ 2, 3, 7, 6 /)

nn = mesh_set(2)%nn
ne = mesh_set(2)%ne
allocate (node_check(nn), conn(nn))
allocate (temp_coord(3,nn*2), temp_elem(8,ne*2), mat_num_mesh(ne*2))
allocate (side_nodes(move_set(1)%csn(1,1)*move_set(1)%csn(2,1),2,num_group))
mat_num_mesh = 0

! Find nodes & elements need to keep
node_check = 2
conn = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do j = 1, move_set(ng)%p2cn
            do i = 2, move_set(ng)%p2c_nn(j)
                node_check(move_set(ng)%prop2case(i,j)) = 1
            enddo
        enddo
    endif
enddo

!num = 0
!do i = 1, inte_set%csn(1,1)
!    do j = 1, inte_set%csn(2,1)
!        node_check(inte_set%cs_nodes1(i,j,inte_set%csn(3,1))) = 1
!        num = num + 1
!        side_nodes(num,1) = inte_set%cs_nodes1(i,j,inte_set%csn(3,1))
!    enddo
!enddo
!num = 0
!do i = 1, inte_set%csn(1,2)
!    do j = 1, inte_set%csn(2,2)
!        node_check(inte_set%cs_nodes2(i,j,inte_set%csn(3,2))) = 1
!        num = num + 1
!        side_nodes(num,2) = inte_set%cs_nodes2(i,j,inte_set%csn(3,2))
!    enddo
!enddo

!temp_nn = 0
!do i = 1, nn
!    if (node_check(i) == 1) then
!        temp_nn = temp_nn + 1
!        temp_coord(:,temp_nn) = mesh_set(2)%node_coord(:,i)
!        conn(i) = temp_nn
!    endif
!enddo
!do i = 1, 2
!    side_nodes(:,i) = conn(side_nodes(:,i))
!enddo

do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 2, move_set(ng)%csn(3,1)
            do j = 1, move_set(ng)%csn(1,1)
                do k = 1, move_set(ng)%csn(2,1)
                    node_check(move_set(ng)%cs_nodes1(j,k,i)) = 0
                enddo
            enddo
        enddo
        do i = 2, move_set(ng)%csn(3,2)
            do j = 1, move_set(ng)%csn(1,2)
                do k = 1, move_set(ng)%csn(2,2)
                    node_check(move_set(ng)%cs_nodes2(j,k,i)) = 0
                enddo
            enddo
        enddo
    endif
enddo

open (Unit=20, File='./output/solid/remesh/check_case_mesh0.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z", "node_check"'
Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Case", N = ', mesh_set(2)%nn, ', E = ', mesh_set(2)%ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do i = 1, mesh_set(2)%nn
    write (20,*) mesh_set(2)%node_coord(:,i), node_check(i)
enddo
do i = 1, mesh_set(2)%ne
    write (20,*) mesh_set(2)%element(:,i)
enddo
close(20)

temp_nn = 0
do i = 1, nn
    if (node_check(i) == 2) then
        temp_nn = temp_nn + 1
        temp_coord(:,temp_nn) = mesh_set(2)%node_coord(:,i)
        conn(i) = temp_nn
    endif
enddo

temp_ne = 0
allocate (sum_p2c_node(2,mesh_set(2)%bn))
sum_p2c_node_num = 0
sum_p2c_node = 0
do i = 1, ne
    num = sum(node_check(mesh_set(2)%element(:,i)))
    if (num == 12) then
        do j = 1, 8
            count = mesh_set(2)%element(j,i)
            if (node_check(count) == 1) then
                check(1) = .TRUE.
                do k = 1, sum_p2c_node_num
                    if (sum_p2c_node(1,k) == count) then
                        check(1) = .FALSE.
                        temp_k = k
                        exit
                    endif
                enddo
                if (check(1)) then
                    temp_nn = temp_nn + 1
                    temp_coord(:,temp_nn) = mesh_set(2)%node_coord(:,count)
                    conn(count) = temp_nn
                    sum_p2c_node_num = sum_p2c_node_num + 1
                    sum_p2c_node(:,sum_p2c_node_num) = (/ count, temp_nn /)
                endif
            endif
        enddo
    endif
    
    if (num == 12 .or. num == 16) then
        temp_ne = temp_ne + 1
        temp_elem(:,temp_ne) = conn(mesh_set(2)%element(:,i))
        elem(:,1) = mesh_set(2)%surf2ori(mesh_set(2)%element(:,i))
        do j = 1, region(1)%num_elements
            elem(:,2) = region(1)%element(4:11,j)
            call find_count(8, elem(:,1), 8, elem(:,2), count)
            if (count == 8) then
                mat_num_mesh(temp_ne) = region(1)%element(3,j)
                exit
            endif
        enddo
    endif
enddo
deallocate (node_check, conn)

open (Unit=20, File='./output/solid/remesh/check_case_mesh1.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Case", N = ', temp_nn, ', E = ', temp_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do i = 1, temp_nn
    write (20,*) temp_coord(:,i)
enddo
do i = 1, temp_ne
    write (20,*) temp_elem(:,i)
enddo
close(20)

! make mesh at region(3) - smoothing region
temp_k = temp_ne
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do iter = 1, 2
            q = 0
            csn = move_set(ng)%csn(:,iter)
            allocate (cs_nodes(csn(1),csn(2),csn(3)))
            if (iter == 1) then
                cs_nodes = move_set(ng)%cs_nodes1
            else
                cs_nodes = move_set(ng)%cs_nodes2
            endif
            pre_nn = temp_nn
            pre_ne = temp_ne
            p(:,1) = mesh_set(2)%node_coord(:,cs_nodes(1,1,1))
            p(:,2) = mesh_set(2)%node_coord(:,cs_nodes(1,1,csn(3)))
            call calc_length_3D(p(:,1), p(:,2), dis(1))
            div_num(1) = int(dis(1)/remesh_dis)+2
            layer_elem = 0
            do i = 1, csn(1)
                if (i /= csn(1)) then
                    do j = 1, mat_layer_num
                        if (mat_layer(j,1) <= i .and. mat_layer(j,2) >= i) then
                            mat_num = mat_layer(j,3)
                            exit
                        endif
                    enddo
                endif
    
                count = temp_nn
                allocate (temp_cs_nodes(csn(2)))
                do j = 1, csn(2)
                    p(:,1) = mesh_set(2)%node_coord(:,cs_nodes(i,j,1))
                    p(:,2) = mesh_set(2)%node_coord(:,cs_nodes(i,j,csn(3)))
                    call find_match_node_num(p(:,1), pre_nn, temp_coord(:,1:pre_nn), temp_cs_nodes(j), remesh_dis*0.01, -1)
                    vec = (p(:,2)-p(:,1))/float(div_num(1))
                    do k = 1, div_num(1)
                        temp_nn = temp_nn + 1
                        temp_coord(:,temp_nn) = p(:,1) + k*vec
                        if (i == 1 .and. j == 1 .and. k == 1) cs_nodes_112(:,iter,ng) = temp_coord(:,temp_nn)
                        if (k == div_num(1)) then
                            q = q + 1
                            side_nodes(q,iter,ng) = temp_nn
                        endif
                    enddo
                enddo
        
                if (i /= csn(1)) then
                    do j = 2, csn(2)
                        do k = 1, div_num(1)
                            temp_ne = temp_ne + 1
                            mat_num_mesh(temp_ne) = mat_num
                            if (i == 1) layer_elem = layer_elem + 1
                            if (k == 1) then                    
                                temp_elem(1:4,temp_ne) = (/ temp_cs_nodes(j-1), count+k+(j-2)*div_num(1), &
                                                            count+k+(j-1)*div_num(1), temp_cs_nodes(j) /)
                            else
                                temp_elem(1:4,temp_ne) = (/ count+k-1+(j-2)*div_num(1), count+k+(j-2)*div_num(1), &
                                                            count+k+(j-1)*div_num(1), count+k-1+(j-1)*div_num(1) /)
                            endif
                        enddo
                    enddo
                endif
                if (i /= 1) then
                    if (i == csn(1)) then
                        num = 0
                        do j = 2, csn(2)
                            do k = 1, div_num(1)
                                num = num + 1
                                if (k == 1) then                    
                                    temp_elem(5:8,pre_ne+(i-2)*layer_elem+num) = (/ temp_cs_nodes(j-1), count+k+(j-2)*div_num(1), &
                                                                                    count+k+(j-1)*div_num(1), temp_cs_nodes(j) /)
                                else
                                    temp_elem(5:8,pre_ne+(i-2)*layer_elem+num) = (/ count+k-1+(j-2)*div_num(1), count+k+(j-2)*div_num(1), &
                                                                                    count+k+(j-1)*div_num(1), count+k-1+(j-1)*div_num(1) /)
                                endif
                            enddo
                        enddo
                    else
                        num = 0
                        do j = 2, csn(2)
                            do k = 1, div_num(1)
                                num = num + 1
                                temp_elem(5:8,pre_ne+(i-2)*layer_elem+num) = temp_elem(1:4,pre_ne+(i-1)*layer_elem+num)
                            enddo
                        enddo
                    endif
                endif
                deallocate(temp_cs_nodes)
            enddo
            deallocate(cs_nodes)
        enddo
    endif
enddo

open (Unit=20, File='./output/solid/remesh/check_case_mesh2.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Case", N = ', temp_nn, ', E = ', temp_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do i = 1, temp_nn
    write (20,*) temp_coord(:,i)
enddo
do i = 1, temp_ne
    write (20,*) temp_elem(:,i)
enddo
close(20)

! reaarange connectivity number
!pre_ne = temp_k
!allocate(elem_check(temp_ne))
!elem_check = 0
!elem_check(1:pre_ne) = 1
!do
!    do i = pre_ne+1, temp_ne
!        if (elem_check(i) == 0) then
!            elem(:,1) = temp_elem(:,i)
!            call make_face_nn(8, 4, elem(:,1), face_nn(:,:,1))
!            do j = temp_ne, 1, -1
!                if (elem_check(j) == 1) then
!                    elem(:,2) = temp_elem(:,j)
!                    call find_count(8, temp_elem(:,i), 8, temp_elem(:,j), count)
!                    if (count == 4) then
!                        call make_face_nn(8, 4, temp_elem(:,j), face_nn(:,:,2))
!                        !write (*,'(8I6)') temp_elem(:,i)
!                        !do k = 1, 6
!                        !    write (*,*) face_nn(:,k,1)
!                        !enddo
!                        !write (*,*) 
!                        !write (*,'(8I6)') temp_elem(:,j)
!                        !do k = 1, 6
!                        !    write (*,*) face_nn(:,k,2)
!                        !enddo
!                        do k = 1, 2
!                            do q = 1, 6
!                                if (k == 1) then
!                                    call find_count(4, face_nn(:,q,k), 8, elem(:,2), count)
!                                else
!                                    call find_count(4, face_nn(:,q,k), 8, elem(:,1), count)
!                                endif
!                                if (count == 4) then
!                                    face(k) = q
!                                    exit
!                                endif
!                            enddo
!                        enddo
!                        do k = 1, 4
!                            if (face_nn(1,face(1),1) == face_nn(k,face(2),2)) then
!                                check(1) = .false.
!                                if (k == 4) then
!                                    if (face_nn(2,face(1),1) == face_nn(1,face(2),2)) then
!                                        check(2) = .true.
!                                    endif
!                                else
!                                    if (face_nn(2,face(1),1) == face_nn(k+1,face(2),2)) then
!                                        check(1) = .true.
!                                    endif
!                                endif
!                                exit
!                            endif
!                        enddo
!                        !write (*,*) check(1)
!                        !write (*,*) face_nn(:,face(1),1)
!                        !write (*,*) face_nn(:,face(2),2)
!                        if (check(1)) then
!                            temp_elem(1,i) = face_nn(1,face(1),1)
!                            temp_elem(5,i) = elem(oppo_face(1,face(1)),1)
!                            do k = 4, 2, -1
!                                temp_elem(6-k,i) = face_nn(k,face(1),1)
!                                temp_elem(10-k,i) = elem(oppo_face(k,face(1)),1)
!                            enddo
!                            
!                        else
!                            temp_elem(1:4,i) = face_nn(:,face(1),1)
!                            do k = 1, 4
!                                temp_elem(k+4,i) = elem(oppo_face(k,face(1)),1)
!                            enddo
!                        endif
!                        elem_check(i) = 1
!                    endif
!                endif
!            enddo
!        endif        
!    enddo
!    if (sum(elem_check) == temp_ne) exit
!enddo
!deallocate(elem_check)

prop_bn = 0
prop_bfn = 0
num = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        num = num + 1
        prop_bn = prop_bn + Hexa(num)%bn
        prop_bfn = prop_bfn + Hexa(num)%bfn
    endif
enddo
allocate (prop_b_node(prop_bn), prop_b_face(4,prop_bfn), prop_conn(new_nn))
sub_num = 0
num = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        num = num + 1
        write (*,*) 'group number:', num
        prop_b_node(sub_num(1)+1:sub_num(1)+Hexa(num)%bn) = Hexa(num)%bound_node + sub_num(3)
        do j = 1, Hexa(num)%bfn
            prop_b_face(:,sub_num(2)+j) = Hexa(num)%bound_face(:,j) + sub_num(3)
        enddo
        sub_num(1) = sub_num(1) + Hexa(num)%bn
        sub_num(2) = sub_num(2) + Hexa(num)%bfn
        sub_num(3) = sub_num(3) + Hexa(num)%nn
    endif
enddo
do i = 1, prop_bn
    prop_conn(prop_b_node(i)) = i
enddo

match_tri_ne = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        match_tri_ne = match_tri_ne + move_set(ng)%match_bfn*2
    endif
enddo
allocate (match_tri_elem(3,match_tri_ne))
open (Unit=20, File='./output/solid/remesh/check_match_surf0.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="check_match_surf"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6, A)') 'ZONE N=', mesh_set(2)%nn , ', E=', match_tri_ne, ', DATAPACKING = POINT, ZONETYPE = FETRIANGLE'
do i = 1, mesh_set(2)%nn
    write (20,*) mesh_set(2)%node_coord(:,i)
enddo
num = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 1, move_set(ng)%match_bfn
            face = move_set(ng)%match_b_face(:,i)
            p(:,1:4) = mesh_set(2)%node_coord(:,face)
            call calc_length_3D(p(:,1), p(:,3), dis(1))
            call calc_length_3D(p(:,2), p(:,4), dis(2))
            if (dis(1) < dis(2)) then
                match_tri_elem(:,num+1) = (/ face(1), face(2), face(3) /)
                match_tri_elem(:,num+2) = (/ face(3), face(4), face(1) /)
            else
                match_tri_elem(:,num+1) = (/ face(1), face(2), face(4) /)
                match_tri_elem(:,num+2) = (/ face(2), face(3), face(4) /)
            endif
            write (20,'(3(I7,2X))') match_tri_elem(:,num+1)
            write (20,'(3(I7,2X))') match_tri_elem(:,num+2)
            num = num + 2
        enddo
    endif
enddo
close(20)
    
allocate (node_check(prop_bn), conn(prop_bn))
node_check = 0
conn = 0
pre_nn = temp_nn
do i = 1, prop_bn
    p(:,1) = new_coord(:,prop_b_node(i))
    do j = 1, match_tri_ne
        p(:,2) = mesh_set(2)%node_coord(:,match_tri_elem(1,j))
        p(:,3) = mesh_set(2)%node_coord(:,match_tri_elem(2,j))
        p(:,4) = mesh_set(2)%node_coord(:,match_tri_elem(3,j))
        call calc_length_3D(p(:,2), p(:,3), dis(1))
        call calc_length_3D(p(:,3), p(:,4), dis(2))
        call calc_length_3D(p(:,2), p(:,4), dis(3))
        dis(2) = (dis(1)+dis(2)+dis(3))/3.0
        call calc_plane(p(:,2), p(:,3), p(:,4), pl)
        dis(1) = abs(pl(1)*p(1,1) + pl(2)*p(2,1) + pl(3)*p(3,1) + pl(4))
        check = .false.
        !if (dis(1) < remesh_dis*0.05) then
        if (dis(1) < dis(2)*0.05) then
            p(:,5) = p(:,1) + pl(1:3)*remesh_dis
            p(:,6) = p(:,1) - pl(1:3)*remesh_dis
            call cross_point_plane(pl, p(:,5), p(:,6), p(:,7), t, check(1))
            if (check(1) == .false.) then  
                call check_inner_area(p(:,2:4), p(:,7), check(2), area)
                if (check(2)) then
                    node_check(i) = 1
                    temp_nn = temp_nn + 1
                    temp_coord(:,temp_nn) = p(:,1)
                    conn(i) = temp_nn
                    exit
                endif
            endif
        endif
        if (check(2)) exit
    enddo
enddo

allocate(elem_pos(4,prop_bfn))
p2cn = sum(node_check)
p2c_ne = 0
pre_ne = temp_ne
do i = 1, prop_bfn
    face = prop_conn(prop_b_face(:,i))
    num = sum(node_check(face))
    if (num == 4) then
        temp_ne = temp_ne + 1
        p2c_ne = p2c_ne + 1
        temp_elem(1:4,temp_ne) = conn(face)
        elem_pos(:,p2c_ne) = temp_elem(1:4,temp_ne) - pre_nn
    endif
enddo

open (Unit=20, File='./output/solid/remesh/check_match_surf1.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="check_match_surf"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6, A)') 'ZONE N=', p2cn , ', E=', p2c_ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
do i = 1, p2cn
    write (20,*) temp_coord(:,i+pre_nn)
enddo
do i = 1, p2c_ne
    write (20,*) elem_pos(:,i)
enddo
close(20)

! calculation seep dis
p2_csn = move_set(1)%p2_csn-2
ori_p2cn = maxval(move_set(1:num_group)%p2cn)
case_nn = mesh_set(2)%nn
allocate (seep_dis(p2_csn,p2cn), p2c_nn(p2cn), end_pos(p2cn))
allocate (temp_val(p2_csn,4), match_conn(case_nn))
seep_dis = 0.0
match_conn = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 1, move_set(ng)%p2cn
            match_conn(move_set(ng)%prop2case(2,i)) = i
        enddo
    endif
enddo

do i = 1, p2cn
    p(:,1) = temp_coord(:,i+pre_nn)
    temp_j = 0
    do ng = 1, num_group
        if (mesh_flag(ng) /= 3) then
            do j = 1, move_set(ng)%p2cn
                p(:,2) = mesh_set(2)%node_coord(:,move_set(ng)%prop2case(2,j))
                call calc_length_3D(p(:,1), p(:,2), dis(1))
                if (dis(1) <= remesh_dis*0.02) then
                    temp_j = j
                    exit
                endif
            enddo
            if (temp_j == 0) then
                do j = 1, move_set(ng)%match_bfn
                    face = move_set(ng)%match_b_face(:,j)
                    p(:,2:5) = mesh_set(2)%node_coord(:,face)
                    !call calc_length_3D(p(:,2), p(:,4), dis(2))
                    !call calc_length_3D(p(:,3), p(:,5), dis(3))
                    !if (dis(2) < dis(3)) dis(2) = dis(3)
                    call near_plane(p(:,2:5), p(:,1), remesh_dis*0.2, near_check) 
                    if (near_check) then
                        temp_j = j
                        face = match_conn(face)
                        min_p2c_nn = minval(move_set(ng)%p2c_nn(face)) - 2
                        do k = 1, 4
                            temp_val(1:min_p2c_nn,k) = move_set(ng)%seep_dis(1:min_p2c_nn,face(k))
                        enddo
                        
                        ! interpolation
                        p2c_nn(i) = min_p2c_nn
                        call rotation_n_interpolation(i, p(:,1:5), p2_csn, min_p2c_nn, temp_val, seep_dis(1:min_p2c_nn,i), remesh_dis*0.01)
                        !if (i == 1137) then
                        !    write (*,*) 'face  :', face
                        !    write (*,*) 'b_face:', move_set(ng)%match_b_face(:,temp_j)
                        !    write (*,*) seep_dis(1:3,i)
                        !endif
                        exit
                    endif
                enddo
            else
                p2c_nn(i) = move_set(ng)%p2c_nn(temp_j)-2
                seep_dis(1:p2c_nn(i),i) = move_set(ng)%seep_dis(1:p2c_nn(i),temp_j)
            endif
            if (temp_j /= 0) exit
        endif
    enddo
    if (temp_j == 0) then
        write (*,*) 'Cannot find the point near matching surface at ', i, i+pre_nn, 'point'
        write (*,'(3(F15.10,2X))') p(:,1)
        stop
    endif
enddo
deallocate (match_conn)

do ng = 1, num_group
    do i = 1, remesh_add(ng)%fix_point_num
        num = prop_conn(Hexa(ng)%fixed_point_group(1,i))
        num = conn(num) - pre_nn
        p(:,1) = temp_coord(:,num+pre_nn)
        min_dis = 10e8
        temp_j = 0
        do j = 1, move_set(ng)%p2cn
            p(:,2) = mesh_set(2)%node_coord(:,move_set(ng)%prop2case(2,j))
            call calc_length_3D(p(:,1), p(:,2), dis(1))
            if (dis(1) <= remesh_dis*0.02) then
                temp_j = j
                exit
            elseif (dis(1) < min_dis) then
                temp_j = j
                min_dis = dis(1)
            endif
        enddo
        p2c_nn(num) = move_set(ng)%p2c_nn(temp_j)-2
        seep_dis(1:p2c_nn(num),num) = move_set(ng)%seep_dis(1:p2c_nn(num),temp_j)
        do j = 2, Hexa(ng)%fixed_point_group_num
            num = prop_conn(Hexa(ng)%fixed_point_group(j,i))
            num = conn(num) - pre_nn
            p2c_nn(num) = move_set(ng)%p2c_nn(temp_j)-2
            seep_dis(1:p2c_nn(num),num) = move_set(ng)%seep_dis(1:p2c_nn(num),temp_j)
        enddo
    enddo
enddo

! calculation seep vec
allocate (seep_vec(3,p2cn))
p(:,4) = (/ 10.0, 0.0, 0.0 /)
do i = pre_nn+1, temp_nn
    p(:,2) = temp_coord(:,i)
    pl = 0.0
    do j = 1, p2c_ne
        do k = 1, 4
            if (elem_pos(k,j) == i-pre_nn) then
                if (k == 1) then
                    p(:,1) = temp_coord(:,elem_pos(4,j)+pre_nn)
                    p(:,3) = temp_coord(:,elem_pos(2,j)+pre_nn)
                elseif (k == 4) then
                    p(:,1) = temp_coord(:,elem_pos(3,j)+pre_nn)
                    p(:,3) = temp_coord(:,elem_pos(1,j)+pre_nn)
                else
                    p(:,1) = temp_coord(:,elem_pos(k-1,j)+pre_nn)
                    p(:,3) = temp_coord(:,elem_pos(k+1,j)+pre_nn)
                endif
                call calc_plane(p(:,1), p(:,2), p(:,3), temp_pl)
                pl = pl - temp_pl
                exit
            endif
        enddo
    enddo
    
    p(:,1) = p(:,2)
    p(1,1) = 0.0        
    call calc_plane(p(:,1), op, p(:,4), temp_pl)
    p(:,1) = p(:,2) + pl(1:3)
    p(:,3) = p(:,1)+10.0*temp_pl(1:3)
    p(:,1) = p(:,1)-10.0*temp_pl(1:3)
    call cross_point_plane(temp_pl, p(:,1), p(:,3), p(:,5), t, parallel)
    p(:,2) = p(:,5) - p(:,2)
    dis(1) = sqrt(p(1,2)**2.0+p(2,2)**2.0+p(3,2)**2.0)
    seep_vec(:,i-pre_nn) = p(:,2)/dis(1)
enddo
    
open (Unit=20, File='./output/solid/remesh/check_match_surf2.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="check_match_surf"'
Write (20,*) 'VARIABLES="x", "y", "z", "p2c_nn", "dis1", "dis2", "vec1", "vec2", "vec3",'
Write (20,'(A,I6,A,I6, A)') 'ZONE N=', p2cn , ', E=', p2c_ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
do i = 1, p2cn
    write (20,*) temp_coord(:,i+pre_nn), p2c_nn(i), seep_dis(1:2,i), seep_vec(:,i)
enddo
do i = 1, p2c_ne
    write (20,*) elem_pos(:,i)
enddo
close(20)


deallocate(conn)
allocate(conn(p2cn))
do i = 1, p2_csn
    end_pos = 0
    conn = 0
    do j = 1, p2cn
        if (p2c_nn(j) >= i) then
            p(:,1) = temp_coord(:,pre_nn+j)
            vec = seep_vec(:,j)
            dis(1) = seep_dis(i,j)
            temp_nn = temp_nn + 1
            temp_coord(:,temp_nn) = p(:,1) + vec*dis(1)
            if (p2c_nn(j) /= i) end_pos(j) = 1
            conn(j) = temp_nn
        endif
    enddo
    do j = 1, mat_layer_num
        if (mat_layer(j,1) <= i .and. mat_layer(j,2) >= i) then
            num = mat_layer(j,3)
            exit
        endif
    enddo
    !write (*,*) i, 'mat_num:', num
    pre_ne2 = temp_ne
    count = 0
    do j = 1, p2c_ne
        check(1) = .TRUE.
        do k = 1, 4
            if (conn(elem_pos(k,j)) == 0) then
                check(1) = .FALSE.
                exit
            endif
        enddo
        if (check(1)) then
            count = count + 1
            end_num = 0
            do k = 1, 4
                temp_elem(4+k,pre_ne+count) = conn(elem_pos(k,j))
                end_num = end_num + end_pos(elem_pos(k,j))
            enddo
            mat_num_mesh(pre_ne+count) = num
            if (end_num == 4) then
                temp_ne = temp_ne + 1
                temp_elem(1:4,temp_ne) = temp_elem(5:8,pre_ne+count)
            endif
        endif
    enddo
    pre_ne = pre_ne2
enddo

open (Unit=20, File='./output/solid/remesh/check_case_mesh3.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Case", N = ', temp_nn, ', E = ', temp_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do i = 1, temp_nn
    write (20,*) temp_coord(:,i)
enddo
do i = 1, temp_ne
    write (20,*) temp_elem(:,i)
enddo
close(20)

deallocate (seep_dis, seep_vec)
!open (Unit=20, File='./output/solid/remesh/check_prop_surf.plt', STATUS='replace', ACTION='write')
!Write (20,'(A)') 'TITLE="Check prop surface mesh"'
!Write (20,*) 'VARIABLES="x", "y", "z", "num"'
!Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Prop_surf", N = ', prop_bn, ', E = ', prop_bfn, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
!do i = 1, prop_bn
!    write (20,*) new_coord(:,prop_b_node(i)), node_check(i)
!enddo
!do i = 1, prop_bfn
!    write (20,*) prop_conn(prop_b_face(:,i))
!enddo
!close(20)

!open (Unit=20, File='./output/solid/remesh/check_case_mesh1.plt', STATUS='replace', ACTION='write')
!Write (20,'(A)') 'TITLE="Check case mesh"'
!Write (20,*) 'VARIABLES="x", "y", "z"'
!Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Case", N = ', temp_nn, ', E = ', temp_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
!do i = 1, temp_nn
!    write (20,*) temp_coord(:,i)
!enddo
!do i = 1, temp_ne
!    write (20,*) temp_elem(:,i)
!enddo
!close(20)

! summation with case nodes & case nodes
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 1, 2
            do j = 1, move_set(1)%csn(1,1)*move_set(1)%csn(2,1)
                p(:,1) = temp_coord(:,side_nodes(j,i,ng))
                dis(2) = 10e8
                del_num = 0
                do k = pre_nn+1, temp_nn
                    p(:,2) = temp_coord(:,k)
                    call calc_length_3D(p(:,1), p(:,2), dis(1))
                    if (dis(1) < remesh_dis*0.01) then
                        del_num = (/ k, side_nodes(j,i,ng) /)
                        exit
                    elseif (dis(2) > dis(1)) then
                        temp_k = k
                        dis(2) = dis(1)
                    endif
                enddo
                if (del_num(1) == 0) then
                    del_num = (/ temp_k, side_nodes(j,i,ng) /)
                endif
                do k = del_num(1), temp_nn-1
                    temp_coord(:,k) = temp_coord(:,k+1)
                enddo
                temp_nn = temp_nn - 1
                do k = 1, temp_ne
                    do q = 1, 8
                        if (temp_elem(q,k) > del_num(1)) then
                            temp_elem(q,k) = temp_elem(q,k) - 1
                        elseif (temp_elem(q,k) == del_num(1)) then
                            temp_elem(q,k) = del_num(2)
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif
enddo
deallocate (side_nodes)

if (sum_p2c_node_num /= 0) then
    do i = 1, sum_p2c_node_num
        p(:,1) = temp_coord(:,sum_p2c_node(2,i))
        dis(2) = 10e8
        del_num = 0
        do k = pre_nn+1, temp_nn
            p(:,2) = temp_coord(:,k)
            call calc_length_3D(p(:,1), p(:,2), dis(1))
            if (dis(1) < remesh_dis*0.01) then
                del_num = (/ k, sum_p2c_node(2,i) /)
                exit
            elseif (dis(2) > dis(1)) then
                temp_k = k
                dis(2) = dis(1)
            endif
        enddo
        if (del_num(1) == 0) then
            del_num = (/ temp_k, sum_p2c_node(2,i) /)
        endif
        do k = del_num(1), temp_nn-1
            temp_coord(:,k) = temp_coord(:,k+1)
        enddo
        temp_nn = temp_nn - 1
        do k = 1, temp_ne
            do q = 1, 8
                if (temp_elem(q,k) > del_num(1)) then
                    temp_elem(q,k) = temp_elem(q,k) - 1
                elseif (temp_elem(q,k) == del_num(1)) then
                    temp_elem(q,k) = del_num(2)
                endif
            enddo
        enddo
    enddo
    deallocate (sum_p2c_node)
endif

open (Unit=20, File='./output/solid/remesh/check_case_mesh.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Case", N = ', temp_nn, ', E = ', temp_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do i = 1, temp_nn
    write (20,*) temp_coord(:,i)
enddo
do i = 1, temp_ne
    write (20,*) temp_elem(:,i)
enddo
close(20)

pre_nn = new_nn
pre_ne = new_ne
new_coord(:,new_nn+1:new_nn+temp_nn) = temp_coord(:,1:temp_nn)
new_elem(:,new_ne+1:new_ne+temp_ne) = temp_elem(:,1:temp_ne)+new_nn
new_nn = new_nn + temp_nn
new_ne = new_ne + temp_ne
deallocate (temp_coord, temp_elem)

! summation with propellant nodes & case nodes
do i = pre_nn, 1, -1
    check(1) = .FALSE.
    do j = 1, prop_bn
        if (prop_b_node(j) == i) then
            if (node_check(j) == 1) then
                check(1) = .TRUE.
            endif
            exit
        endif
    enddo
    
    if (check(1)) then
        p(:,1) = new_coord(:,i)
        dis(2) = 10e8
        del_num = 0
        do j = pre_nn+1, new_nn
            p(:,2) = new_coord(:,j)
            call calc_length_3D(p(:,1), p(:,2), dis(1))
            if (dis(1) < remesh_dis*0.01) then
                !del_num = (/ i, j-1 /)
                del_num = (/ i, j /)
                exit
            elseif (dis(2) > dis(1)) then
                temp_k = j
                dis(2) = dis(1)
            endif
        enddo
        if (del_num(1) == 0) then
            !del_num = (/ i, temp_k-1 /)
            del_num = (/ i, temp_k /)
        endif
        do j = del_num(2), new_nn-1
            new_coord(:,j) = new_coord(:,j+1)
        enddo
        new_nn = new_nn - 1
        do j = 1, new_ne
            do k = 1, 8
                if (new_elem(k,j) == del_num(2)) then
                    new_elem(k,j) = del_num(1)
                elseif (new_elem(k,j) > del_num(2)) then
                    new_elem(k,j) = new_elem(k,j) - 1
                endif
            enddo
        enddo
    endif
enddo
deallocate (node_check)

do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        do i = 1, 2
            do j = 1, move_set(ng)%csn(2,1)+1
                if (i == 1) then
                    if (j == move_set(ng)%csn(2,1)+1) then
                        !p(:,1) = mesh_set(2)%node_coord(:,inte_set%cs_nodes1(1,1,2))
                        p(:,1) = cs_nodes_112(:,i,ng)
                    else
                        p(:,1) = mesh_set(2)%node_coord(:,move_set(ng)%cs_nodes1(1,j,1))
                    endif
                else
                    if (j == move_set(ng)%csn(2,1)+1) then
                        !p(:,1) = mesh_set(2)%node_coord(:,inte_set%cs_nodes2(1,1,2))
                        p(:,1) = cs_nodes_112(:,i,ng)
                    else
                        p(:,1) = mesh_set(2)%node_coord(:,move_set(ng)%cs_nodes2(1,j,1))
                    endif
                endif
                dis(2) = 10e8
                temp_k = 0
                do k = 1, new_nn
                    p(:,2) = new_coord(:,k)
                    call calc_length_3D(p(:,1), p(:,2), dis(1))
                    if (dis(1) < remesh_dis*0.01) then
                        temp_k = k
                        exit
                    elseif (dis(2) > dis(1)) then
                        temp_k = k
                        dis(2) = dis(1)
                    endif
                enddo
                if (i == 1) then
                    if (j == move_set(ng)%csn(2,1)+1) then
                        move_set(ng)%cs_nodes1(1,1,2) = temp_k
                    else
                        move_set(ng)%cs_nodes1(1,j,1) = temp_k
                    endif
                else
                    if (j == move_set(ng)%csn(2,1)+1) then
                        move_set(ng)%cs_nodes2(1,1,2) = temp_k
                    else
                        move_set(ng)%cs_nodes2(1,j,1) = temp_k
                    endif
                endif
            enddo
        enddo
    endif
enddo
deallocate(cs_nodes_112)

! element rearrange according mat_num
allocate (temp_elem(8,new_ne))
temp_elem = new_elem
pre_ne = 0
cur_ng = 0
do ng = 1, num_group
    if (mesh_flag(ng) /= 3) then
        cur_ng = cur_ng + 1
        mat_array(cur_ng) = Hexa(ng)%ne + pre_ne
        pre_ne = pre_ne + Hexa(ng)%ne
    endif
enddo
do i = 1, mat_layer_num
    num = 0
    do j = 1, temp_ne
        if (mat_num_mesh(j) == i+num_group) then
            num = num + 1
            new_elem(:,pre_ne+num) = temp_elem(:,j+mat_array(cur_ng))
        endif
    enddo
    pre_ne = pre_ne + num
    mat_array(i+cur_ng) = pre_ne
enddo
deallocate (temp_elem)
!write (*,*) inte_set%cs_nodes1(1,:,1), inte_set%cs_nodes1(1,1,2)
!write (*,*) inte_set%cs_nodes2(1,:,1), inte_set%cs_nodes2(1,1,2)

open (Unit=20, File='./output/solid/remesh/final_mesh.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
do i = 1, mat_array_num
    if (i == 1) then
        k = 1
    else
        k = mat_array(i-1)+1
    endif
    q = mat_array(i)
    num = q - k + 1
    Write (20,'(A,I2,A,I6,A,I6,A)') 'ZONE T = "mat_num', i, '", N = ', new_nn, ', E = ', num, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
    do j = 1, new_nn
        write (20,*) new_coord(:,j)
    enddo
    do j = k, q
        write (20,*) new_elem(:,j)
    enddo
enddo
close(20)

write (*,*) '>> Hexa_mesh_rearrange'
call Hexa_mesh_rearrange(new_nn, new_coord, new_ne, new_elem)

open (Unit=20, File='./output/solid/remesh/final_mesh2.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check case mesh"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6,A)') 'ZONE N = ', new_nn, ', E = ', new_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do j = 1, new_nn
    write (20,*) new_coord(:,j)
enddo
do j = 1, new_ne
    write (20,*) new_elem(:,j)
enddo
close(20)

end subroutine remesh_ale_region2


subroutine rotation_n_interpolation(ci, p, p2_csn, min_p2c_nn, temp_val, seep_dis, tol)

implicit none

integer, intent(in) :: ci, p2_csn, min_p2c_nn
real, intent(in) :: tol, p(3,5), temp_val(p2_csn,4)
real, intent(inout) :: seep_dis(min_p2c_nn)

integer :: i, j, err
real :: f1, f2, det_j, u , v
real :: xi, eta , norm, temp, ang_tol, move_angle
real :: move_vec(3), temp_p(3,5), temp_node(2,4), tn(3), jacobi(4), target_node(2)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

ang_tol = 0.05
move_vec = p(:,2)
do i = 1, 5
    temp_p(:,i) = p(:,i) - move_vec
enddo

! 1. y axis rotation
temp_node(:,1) = (/ temp_p(1,3), temp_p(3,3) /)
temp_node(:,2) = (/ 0.0, 0.0 /)
temp_node(:,3) = (/ 1.0, 0.0 /)
call calc_angle(temp_node(:,1), temp_node(:,2), temp_node(:,3), move_angle)

if ( (move_angle > 0.+tol) .AND. (move_angle < 360.-tol) ) then
	!write (*,'(A,F6.2)') "Rotation at y axis : ", move_angle
	do i = 1, 5
		tn(:) = temp_p(:,i)
		temp_p(1,i) = tn(1)*cos(move_angle*pi/180.) + tn(3)*sin(move_angle*pi/180.)
		temp_p(3,i) = -tn(1)*sin(move_angle*pi/180.) + tn(3)*cos(move_angle*pi/180.)
	enddo
endif

! 2. z axis rotation
temp_node(:,1) = (/ temp_p(1,3), temp_p(2,3) /)
temp_node(:,2) = (/ 0.0, 0.0 /)
temp_node(:,3) = (/ 1.0, 0.0 /)
call calc_angle(temp_node(:,1), temp_node(:,2), temp_node(:,3), move_angle)

if ( (move_angle > 0.+tol) .AND. (move_angle < 360.-tol) ) then
	!write (*,'(A,F6.2)') "Rotation at z axis : ", move_angle
    move_angle = -move_angle
	do i = 1, 5
		tn(:) = temp_p(:,i)
		temp_p(1,i) = tn(1)*cos(move_angle*pi/180.) - tn(2)*sin(move_angle*pi/180.)
		temp_p(2,i) = tn(1)*sin(move_angle*pi/180.) + tn(2)*cos(move_angle*pi/180.)
	enddo
endif

! 3. x axis rotation
temp_node(:,1) = (/ temp_p(2,5), temp_p(3,5) /)
temp_node(:,2) = (/ 0.0, 0.0 /)
temp_node(:,3) = (/ 1.0, 0.0 /)
call calc_angle(temp_node(:,1), temp_node(:,2), temp_node(:,3), move_angle)

if ( (move_angle > 0.+tol) .AND. (move_angle < 360.-tol) ) then
	!write (*,'(A,F6.2)') "Rotation at x axis : ", move_angle
    move_angle = -move_angle
	do i = 1, 5
		tn(:) = temp_p(:,i)
		temp_p(2,i) = tn(2)*cos(move_angle*pi/180.) - tn(3)*sin(move_angle*pi/180.)
		temp_p(3,i) = +tn(2)*sin(move_angle*pi/180.) + tn(3)*cos(move_angle*pi/180.)
	enddo
endif

!if (ci == 1137) then
!    do i = 1, 5
!        write (*,'(3(F15.10,2X))') temp_p(:,i)
!    enddo
!    write (*,*) 
!endif

! 4. calculation u, v
target_node = temp_p(1:2,1)
do i = 1, 4
    temp_node(:,i) = temp_p(1:2,i+1)
enddo
err = 0 
xi = 0.0 ; eta = 0.0
j = 0 
do 
    f1 = ( (1-xi)*(1-eta)*temp_node(1,1) + (1+xi)*(1-eta)*temp_node(1,2)+ &
    (1+xi)*(1+eta)*temp_node(1,3)+(1-xi)*(1+eta)*temp_node(1,4) )*0.25-target_node(1)

    f2 = ((1-xi)*(1-eta)*temp_node(2,1) + (1+xi)*(1-eta)*temp_node(2,2)+ &
    (1+xi)*(1+eta)*temp_node(2,3) + (1-xi)*(1+eta)*temp_node(2,4))*0.25-target_node(2)
    
    jacobi(1) = (-1+eta)*temp_node(1,1) +(1-eta)*temp_node(1,2) &
            + (1+eta)*temp_node(1,3) +(-1-eta)*temp_node(1,4)
    jacobi(3) = (-1+eta)*temp_node(2,1) +(1-eta)*temp_node(2,2) &
            + (1+eta)*temp_node(2,3) +(-1-eta)*temp_node(2,4)
    jacobi(2) = (-1+xi)*temp_node(1,1) +(-1-xi)*temp_node(1,2) &
            + (1+xi)*temp_node(1,3) +(1-xi)*temp_node(1,4)
    jacobi(4) = (-1+xi)*temp_node(2,1) +(-1-xi)*temp_node(2,2) &
            + (1+xi)*temp_node(2,3) +(1-xi)*temp_node(2,4)

    jacobi = jacobi * 0.25
    det_j = jacobi(1)*jacobi(4)-jacobi(2)*jacobi(3)
    temp = jacobi(1)
    jacobi(1) = jacobi(4)
    jacobi(2) = -jacobi(2)
    jacobi(3) = -jacobi(3)
    jacobi(4) = temp
    jacobi = jacobi / det_j

    u = jacobi(1) * f1 + jacobi(2) * f2
    v = jacobi(3) * f1 + jacobi(4) * f2
    xi = xi - u
    eta = eta - v

    norm = sqrt( u**2 + v**2 )
    if ( norm < 1.0e-8 ) exit

    j = j+ 1 
    if (j > 100 ) then
        err = 1
        exit
    endif
enddo

seep_dis = 0.0
if (err /= 1) then
    do i = 1, min_p2c_nn
        seep_dis(i) = ((1-xi)*(1-eta)*temp_val(i,1) + (1+xi)*(1-eta)*temp_val(i,2)+ &
                      (1+xi)*(1+eta)*temp_val(i,3)+(1-xi)*(1+eta)*temp_val(i,4))*0.25
    enddo
endif

!if (ci == 1137) then
!    write (*,*) 'xi, eta:', xi, eta
!    write (*,'(4(F15.10,2X))') temp_val(1,1), temp_val(1,2), temp_val(1,3), temp_val(1,4)
!    write (*,*) 'seep_dis(1):', seep_dis(1)
!endif

end subroutine rotation_n_interpolation