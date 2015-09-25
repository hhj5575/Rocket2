subroutine make_hexa_mesh(ng, cur_ng, add_nn, add_coord, fsi_flag)
    
use remesh_domain
implicit none

integer, intent(in) :: ng, cur_ng, fsi_flag, add_nn
real, intent(in) :: add_coord(3,20)

integer :: nn, ne, sub_nn, sub_ne
integer :: elem(8,50000)
real :: node(3,50000)

integer :: i, j, k, sd, count, ce, num, status, memory, temp_k, bn, bpn, max_num
integer :: hexa_line(2,12), div_num(3), line(2), line2(2)
real :: vec(3,3), vec_dis(3), p(3,6), cp(3), op(3), val(3), temp_base(3,3), temp_val(3), elem_dis(4), pl(4)
real :: tol_dis, dis, ratio, remainder_dis, min_dis, area, min_area, t
integer, allocatable :: del_node(:,:), line_node(:,:), proj_pn(:), line_node_num(:), bound_lv(:), hexa_ben(:,:)
real, allocatable :: upper_node(:,:), hexa_ben_coord(:,:,:)
logical :: inner_check, pl_check, check
character(len=50) :: fn

nn = 0
ne = 0
hexa_line(1,:) = (/ 1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 3, 4 /)
hexa_line(2,:) = (/ 2, 4, 6, 8, 3, 1, 7, 5, 5, 6, 7, 8 /)

fn = './output/solid/remesh/check_sub_hexa_00.plt'
write (fn(38:39),'(I2.2)') cur_ng
open (Unit=30, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)
Write (30,'(A)') 'TITLE="check_sub_hexa"'
do sd = 1, remesh_add(ng)%sub_hexa_num
    do i = 1, 2
        if (remesh_add(ng)%sub_hexa(hexa_line(i,1),sd) > Tri%nn) then
            num = remesh_add(ng)%sub_hexa(hexa_line(i,1),sd) - Tri%nn
            p(:,i) = add_coord(:,num)
        else
            p(:,i) = Tri%node(:,remesh_add(ng)%sub_hexa(hexa_line(i,1),sd))
        endif
    enddo
    vec(:,1) = p(:,2) - p(:,1)
    call calc_length_3D(p(:,1), p(:,2), vec_dis(1))
    if (remesh_add(ng)%sub_hexa(hexa_line(2,2),sd) > Tri%nn) then
        num = remesh_add(ng)%sub_hexa(hexa_line(2,2),sd) - Tri%nn
        p(:,3) = add_coord(:,num)
    else
        p(:,3) = Tri%node(:,remesh_add(ng)%sub_hexa(hexa_line(2,2),sd))
    endif
    if (remesh_add(ng)%sub_hexa(hexa_line(1,2),sd) > Tri%nn) then
        num = remesh_add(ng)%sub_hexa(hexa_line(1,2),sd) - Tri%nn
        p(:,4) = add_coord(:,num)
    else
        p(:,4) = Tri%node(:,remesh_add(ng)%sub_hexa(hexa_line(1,2),sd))
    endif
    vec(:,2) = p(:,4) - p(:,3)
    call calc_length_3D(p(:,3), p(:,4), vec_dis(2))
    div_num = 0
    do i = 1, 3
        do j = 1, Tri%begn
            do k = 1, Tri%b_edge_group_num(j)
                ce = Tri%b_edge_group(k,j)
                line = (/ Tri%b_edge(1,ce), Tri%b_edge(Tri%bedn(ce),ce) /)
                line2 = remesh_add(ng)%sub_hexa(hexa_line(1:2,i*4-3),sd)
                call find_count(2, line2, 2, line, count)
                if (count == 2) then
                    div_num(i) = Tri%b_edge_divnum(j)
                    !write (*,*) line, div_num(i)
                    exit
                endif
            enddo
            if (div_num(i) /= 0) exit
        enddo
    enddo
    
    sub_nn = 0
    sub_ne = 0
    do i = 1, div_num(1)+1
        p(:,5) = p(:,1) + vec(:,1)*(float(i-1)/float(div_num(1)))
        p(:,6) = p(:,3) + vec(:,2)*(float(i-1)/float(div_num(1)))
        vec(:,3) = p(:,6)-p(:,5)
        do j = 1, div_num(2)+1
            sub_nn = sub_nn + 1
            node(:,sub_nn+nn) = p(:,5) + vec(:,3)*(float(j-1)/float(div_num(2)))
        enddo
    enddo
    
    allocate (upper_node(3,sub_nn))
    if (remesh_add(ng)%sub_hexa(hexa_line(1,3),sd) > Tri%nn) then
        num = remesh_add(ng)%sub_hexa(hexa_line(1,3),sd) - Tri%nn
        p(:,1) = add_coord(:,num)
    else
        p(:,1) = Tri%node(:,remesh_add(ng)%sub_hexa(hexa_line(1,3),sd))
    endif
    if (remesh_add(ng)%sub_hexa(hexa_line(2,3),sd) > Tri%nn) then
        num = remesh_add(ng)%sub_hexa(hexa_line(2,3),sd) - Tri%nn
        p(:,2) = add_coord(:,num)
    else
        p(:,2) = Tri%node(:,remesh_add(ng)%sub_hexa(hexa_line(2,3),sd))
    endif
    vec(:,1) = p(:,2) - p(:,1)
    call calc_length_3D(p(:,1), p(:,2), vec_dis(1))
    
    if (remesh_add(ng)%sub_hexa(hexa_line(2,4),sd) > Tri%nn) then
        num = remesh_add(ng)%sub_hexa(hexa_line(2,4),sd) - Tri%nn
        p(:,3) = add_coord(:,num)
    else
        p(:,3) = Tri%node(:,remesh_add(ng)%sub_hexa(hexa_line(2,4),sd))
    endif
    if (remesh_add(ng)%sub_hexa(hexa_line(1,4),sd) > Tri%nn) then
        num = remesh_add(ng)%sub_hexa(hexa_line(1,4),sd) - Tri%nn
        p(:,4) = add_coord(:,num)
    else
        p(:,4) = Tri%node(:,remesh_add(ng)%sub_hexa(hexa_line(1,4),sd))
    endif
    vec(:,2) = p(:,4) - p(:,3)
    call calc_length_3D(p(:,3), p(:,4), vec_dis(2))
    num = 0
    do i = 1, div_num(1)+1
        p(:,5) = p(:,1) + vec(:,1)*(float(i-1)/float(div_num(1)))
        p(:,6) = p(:,3) + vec(:,2)*(float(i-1)/float(div_num(1)))
        vec(:,3) = p(:,6)-p(:,5)
        do j = 1, div_num(2)+1
            num = num + 1
            upper_node(:,num) = p(:,5) + vec(:,3)*(float(j-1)/float(div_num(2)))
        enddo
    enddo
    
    do i = 1, div_num(3)
        do j = 1, num
            vec(:,3) = upper_node(:,j) - node(:,j+nn)
            sub_nn = sub_nn + 1
            node(:,sub_nn+nn) = node(:,j+nn) + vec(:,3)*(float(i)/float(div_num(3)))
        enddo
    enddo
    do k = 1, div_num(3)
        do j = 1, div_num(1)
            do i = 1, div_num(2)
                sub_ne = sub_ne + 1
                elem(1,sub_ne+ne) = (j-1)*(div_num(2)+1)+(k-1)*num+i + nn
                elem(2,sub_ne+ne) = (j-1)*(div_num(2)+1)+(k-1)*num+i+1 + nn
                elem(3,sub_ne+ne) = j*(div_num(2)+1)+(k-1)*num+i+1 + nn
                elem(4,sub_ne+ne) = j*(div_num(2)+1)+(k-1)*num+i + nn
                elem(5:8,sub_ne+ne) = elem(1:4,sub_ne+ne)+num
            enddo
        enddo
    enddo
    !Write (30,'(A,I3,A,I5,A,I5, A)') 'ZONE T = "sub', sd, '", N=', sub_nn , ', E=', sub_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
    !do i = 1, sub_nn
	   ! write(30,*) node(:,i+nn)
    !enddo
    !do i = 1, sub_ne
	   ! write (30,'(8I6)') elem(:,i+ne)-nn
    !enddo
    nn = nn + sub_nn
    ne = ne + sub_ne
    deallocate (upper_node)
enddo

allocate(del_node(2,nn))
del_node = 0
do i = 1, nn
    if (del_node(1,i) == 0) then
        do j = i+1, nn
            call calc_length_3D(node(:,i), node(:,j), dis)
            if (dis < 10e-7) then
                del_node(1,j) = 1
                del_node(2,j) = i
            endif
        enddo
    endif
enddo
do i = nn, 1, -1
    if (del_node(1,i) == 1) then
        do j = i, nn-1
            node(:,j) = node(:,j+1)
        enddo
        do j = 1, ne
            do k = 1, 8
                if (elem(k,j) == i) then
                    elem(k,j) = del_node(2,i)
                elseif (elem(k,j) >= i) then
                    elem(k,j) = elem(k,j) - 1
                endif
            enddo
        enddo
    endif
enddo
nn = nn - sum(del_node(1,:))
deallocate(del_node)

call Hexa_mesh_rearrange(nn, node(:,1:nn), ne, elem(:,1:ne))

Hexa(cur_ng)%nn = nn
Hexa(cur_ng)%ne = ne
allocate (Hexa(cur_ng)%node(3,nn), Hexa(cur_ng)%elem(8,ne), Hexa(cur_ng)%node_lv(2,nn))
Hexa(cur_ng)%node = node(:,1:nn)
Hexa(cur_ng)%elem = elem(:,1:ne)
Hexa(cur_ng)%node_lv = 0
call remesh_get_boundary_info_3D(cur_ng, 0)
bn = Hexa(cur_ng)%bn

num = maxval(Tri%b_edge_divnum(1:Tri%begn))+1
allocate (line_node(num,Tri%ben), line_node_num(Tri%ben))
line_node_num = 0
line_node = 0
do i = 1, Tri%ben
    !if (Tri%bedn(i) > 2) then
        call tri_or_add_node(Tri%nn, Tri%node, add_nn, add_coord(:,1:add_nn), Tri%b_edge(1,i), p(:,1))
        call tri_or_add_node(Tri%nn, Tri%node, add_nn, add_coord(:,1:add_nn), Tri%b_edge(Tri%bedn(i),i), p(:,2))
        !p(:,1) = Tri%node(:,Tri%b_edge(1,i))
        !p(:,2) = Tri%node(:,Tri%b_edge(Tri%bedn(i),i))
        vec(:,1) = p(:,2)-p(:,1)
    
        div_num(1) = 0
        do j = 1, Tri%begn
            do k = 1, Tri%b_edge_group_num(j)
                if (Tri%b_edge_group(k,j) == i) then
                    div_num(1) = Tri%b_edge_divnum(j)
                    !write (*,*) i, 'edge div_num:', div_num(1)
                    exit
                endif
            enddo
            if (div_num(1) /= 0) exit
        enddo
        line_node_num(i) = div_num(1)+1
        min_dis = 10e8
        temp_k = 0
        do j = 0, div_num(1)
            p(:,3) = p(:,1) + vec(:,1)*(float(j)/float(div_num(1)))
            do k = 1, nn
                call calc_length_3D(node(:,k), p(:,3), dis)
                if (dis < 10e-8) then
                    line_node(j+1,i) = k
                    exit
                elseif (min_dis > dis) then
                    min_dis = dis
                    temp_k = k
                endif
            enddo
            if (line_node(j+1,i) == 0) then
                if (temp_k == 0) stop 'Cannot find line node!! (make_hexa_mesh)'
                line_node(j+1,i) = temp_k
            endif
        enddo
    if (Tri%bedn(i) > 2) then
        tol_dis = 0.0
        do j = 1, Tri%bedn(i)-1
            call tri_or_add_node(Tri%nn, Tri%node, add_nn, add_coord(:,1:add_nn), Tri%b_edge(j,i), p(:,1))
            call tri_or_add_node(Tri%nn, Tri%node, add_nn, add_coord(:,1:add_nn), Tri%b_edge(j+1,i), p(:,2))
            !p(:,1) = Tri%node(:,Tri%b_edge(j,i))
            !p(:,2) = Tri%node(:,Tri%b_edge(j+1,i))
            call calc_length_3D(p(:,1), p(:,2), dis)
            tol_dis = tol_dis + dis
        enddo
        ratio = tol_dis/float(div_num(1))
        memory = 1
        !cp = Tri%node(:,Tri%b_edge(1,i))
        call tri_or_add_node(Tri%nn, Tri%node, add_nn, add_coord(:,1:add_nn), Tri%b_edge(1,i), cp)
        do j = 1, div_num(1)-1
            remainder_dis = ratio
            p(:,1) = cp
            call tri_or_add_node(Tri%nn, Tri%node, add_nn, add_coord(:,1:add_nn), Tri%b_edge(memory+1,i), p(:,2))
            !p(:,2) = Tri%node(:,Tri%b_edge(memory+1,i))
            vec(:,1) = p(:,2) - p(:,1)
            call calc_length_3D(p(:,1), p(:,2), dis)
            if (remainder_dis < dis) then
                cp = p(:,1) + vec(:,1)/dis*remainder_dis
                node(:,line_node(j+1,i)) = cp
            else
                remainder_dis = remainder_dis - dis
                do k = memory+1, Tri%bedn(i)-1
                    call tri_or_add_node(Tri%nn, Tri%node, add_nn, add_coord(:,1:add_nn), Tri%b_edge(k,i), p(:,1))
                    call tri_or_add_node(Tri%nn, Tri%node, add_nn, add_coord(:,1:add_nn), Tri%b_edge(k+1,i), p(:,2))
                    !p(:,1) = Tri%node(:,Tri%b_edge(k,i))
                    !p(:,2) = Tri%node(:,Tri%b_edge(k+1,i))
                    call calc_length_3D(p(:,1), p(:,2), dis)
                    if (dis < remainder_dis) then
                        remainder_dis = remainder_dis - dis
                    else
                        node(:,line_node(j+1,i)) = p(:,1) + (p(:,2)-p(:,1))*remainder_dis/dis
                        cp = node(:,line_node(j+1,i))
                        memory = k
                        exit
                    endif
                enddo
            endif
        enddo
    endif
enddo
!write (*,*) 'Tri%ben:', Tri%ben
!do i = 1, Tri%ben
!    write (*,*) 'line_node_num', i, ':', line_node_num(i)
!    write (*,*) line_node(1:line_node_num(i),i)
!    write (*,*)
!enddo

Write (30,'(A,I5,A,I5, A)') 'ZONE N=', nn , ', E=', ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do i = 1, nn
	 write(30,*) Hexa(cur_ng)%node(:,i)
enddo
do i = 1, ne
	 write (30,'(8I6)') Hexa(cur_ng)%elem(:,i)
enddo
close(30)

allocate (proj_pn(bn), bound_lv(nn))
proj_pn = 0
bound_lv = 0
call set_proj_pn(cur_ng, bn, proj_pn, nn, bound_lv, Tri%ben, num, line_node_num, line_node)
Hexa(cur_ng)%node = node(:,1:nn)

fn = './output/solid/remesh/check_sub_proj_00.plt'
write (fn(38:39),'(I2.2)') cur_ng
open (Unit=30, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)
Write (30,'(A)') 'TITLE="check_sub_proj"'
Write (30,*) 'VARIABLES="x", "y", "z", "proj_pn","lv"'
Write (30,'(A,I5,A,I5, A)') 'ZONE N=', nn , ', E=', ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
do i = 1, nn
    num = Hexa(cur_ng)%bound_conn(i)
    if (num /= 0) num = proj_pn(num)        
	 write(30,*) node(:,i), num, bound_lv(i)
enddo
do i = 1, ne
	 write (30,'(8I6)') elem(:,i)
enddo
close(30)

call projection_point_to_plane(cur_ng, bn, proj_pn, nn, node(:,1:nn))

Hexa(cur_ng)%node = node
Hexa(cur_ng)%node_lv(1,1:nn) = bound_lv(1:nn)
do i = 1, Tri%np2
    min_dis = 10e8
    temp_k = 0
    do j = 1, nn
        call calc_length_3D(Tri%node(:,Tri%bp2(i)), Hexa(cur_ng)%node(:,j), dis)
        if (dis < 10e-8) then
            !Hexa(ng)%node_lv(1,j) = 3
            temp_k = j
            exit
        elseif (min_dis > dis) then
            temp_k = j
            min_dis = dis
        endif
    enddo
    !if (temp_k /= 0) then
        Hexa(cur_ng)%node_lv(1,temp_k) = 3
        if (fsi_flag == 2) Hexa(cur_ng)%bp2_coord(:,i) = node(:,temp_k)
    !endif
enddo

if (remesh_add(ng)%fix_point_num /= 0) call adjust_fix_coord(cur_ng)
call surface_smoothing(cur_ng, 0)
call projection_point_to_plane(cur_ng, bn, proj_pn, nn, Hexa(ng)%node(:,1:nn))
call inner_node_smoothing(cur_ng)
call write_tecplot_3D(Hexa(cur_ng)%nn, Hexa(cur_ng)%node(:,1:Hexa(cur_ng)%nn), Hexa(cur_ng)%ne, Hexa(cur_ng)%elem(:,1:Hexa(cur_ng)%ne), 1, cur_ng, 2) 

if (fsi_flag == 2) then
    max_num = (maxval(Tri%b_edge_divnum)+1)*3
    !allocate (Hexa(ng)%b_edge(max_num,Hexa(ng)%ben))
    !Hexa(ng)%b_edge = 0
    allocate (Hexa_ben(max_num,Hexa(cur_ng)%ben), hexa_ben_coord(3,max_num,Hexa(cur_ng)%ben))
    Hexa(cur_ng)%bedn = 0
    Hexa_ben = 0
    hexa_ben_coord = 0.0
    do i = 1, Hexa(cur_ng)%ben
        do j = 1, Hexa(cur_ng)%cen(i)
            num = Hexa(cur_ng)%control_edge_number(j,i)
            if (j == 1) then
                Hexa(cur_ng)%bedn(i) = Hexa(cur_ng)%bedn(i) + 1
                Hexa_ben(Hexa(cur_ng)%bedn(i),i) = line_node(1,num)
                hexa_ben_coord(:,Hexa(cur_ng)%bedn(i),i) = Hexa(cur_ng)%node(:,line_node(1,num))
            endif
            do k = 2, line_node_num(num)
                Hexa(cur_ng)%bedn(i) = Hexa(cur_ng)%bedn(i) + 1
                Hexa_ben(Hexa(cur_ng)%bedn(i),i) = line_node(k,num)
                hexa_ben_coord(:,Hexa(cur_ng)%bedn(i),i) = Hexa(cur_ng)%node(:,line_node(k,num))
            enddo
        enddo
    enddo
    max_num = maxval(Hexa(cur_ng)%bedn)
    allocate (Hexa(cur_ng)%control_edge_coord(3,max_num,Hexa(cur_ng)%ben))
    Hexa(cur_ng)%control_edge_coord = 0.0
    write (*,*) 'number of group:', cur_ng
    write (*,*) 'Hexa(cur_ng)%ben:', Hexa(cur_ng)%ben
    do i = 1, Hexa(cur_ng)%ben
        !write (*,*) 'Hexa(cur_ng)%bedn(', i, '):', Hexa(cur_ng)%bedn(i)
        do j = 1, Hexa(cur_ng)%bedn(i)
            Hexa(cur_ng)%control_edge_coord(:,j,i) = hexa_ben_coord(:,j,i)
        enddo
        !write (*,*) Hexa_ben(1:Hexa(cur_ng)%bedn(i),i)
        !write (*,*)
    enddo
    deallocate(Hexa_ben, Hexa_ben_coord)
endif

deallocate(line_node, line_node_num)
deallocate (proj_pn, bound_lv)

end subroutine make_hexa_mesh