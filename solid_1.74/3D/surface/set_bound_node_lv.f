subroutine set_bound_node_lv(nn, bn, bfn, node, bound_node, bound_face, bound_lv)

implicit none

integer, intent(in) :: nn, bn, bfn
integer, intent(in) :: bound_node(bn), bound_face(4,bfn)
real, intent(in) :: node(3,nn)
integer, intent(inout) :: bound_lv(bn)

integer :: rn, i, j, k, q, cn, temp_k, count, ln, ele_nodes, this(2), temp_bound_face(4,bfn)
real :: cri_ang
integer, allocatable :: temp_line(:,:), line(:,:), node_lv(:)
real, allocatable :: line_ang(:)

rn = 1
allocate(temp_line(4,bfn*3))

do i = 1, bfn
    temp_bound_face(:,i) = bound_node(bound_face(:,i))
enddo
ln = 0
temp_line = 0
count = 0
do i = 1, bfn
    do j = 1, 4
        if (j == 4) then
            this = (/ temp_bound_face(4,i), temp_bound_face(1,i) /)
        else
            this = (/ temp_bound_face(j,i), temp_bound_face(j+1,i) /)
        endif
        
        do k = 1, ln
            count = 0
            do q = 1, 2
                if (temp_line(q,k) == this(1)) count = count + 1
                if (temp_line(q,k) == this(2)) count = count + 1
            enddo
            if (count == 2) then
                temp_k = k
                exit
            endif
        enddo
        if (count == 2) then
            temp_line(4,temp_k) = i
        else
            ln = ln + 1
            temp_line(1:2,ln) = this
            temp_line(3,ln) = i
        endif
    enddo
enddo

allocate (line(4,ln), line_ang(ln), node_lv(nn))
line = temp_line(:,1:ln)
line_ang= 0.0
ele_nodes = 4
call calc_line_angle(ele_nodes, nn, bfn, ln, node, temp_bound_face, line, line_ang)
write (*,*) '>> complete calculation line angle'

cri_ang = 24.
call find_edge_node_lv(nn, node_lv, ln, line, line_ang, cri_ang)
do i = 1, bn
    cn = bound_node(i)
    bound_lv(i) = node_lv(cn)
enddo

deallocate (line, line_ang, node_lv)

end subroutine set_bound_node_lv


subroutine calc_line_angle(ele_nodes, nn, ne, ln, node, elem, line, line_ang)
    
implicit none

integer, intent(in) :: ele_nodes, nn, ne, ln, elem(ele_nodes,ne), line(4,ln)
real, intent(in) :: node(3,nn)
real, intent(inout) :: line_ang(ln)

integer :: i, j
real :: temp_node(3,3), elem_pl(4,ne), op(3), ang

do i = 1, ne
    do j = 1, 3
        temp_node(:,j) = node(:,elem(j,i))
    enddo
    call calc_plane(temp_node(:,1), temp_node(:,2), temp_node(:,3), elem_pl(:,i))
enddo

op = (/ 0.0, 0.0, 0.0 /)
do i = 1, ln
    if (line(4,i) /= 0) then
        do j = 1, 2
            temp_node(:,j) = elem_pl(1:3,line(j+2,i))
        enddo
        call calc_angle_cos(temp_node(:,1), op, temp_node(:,2), ang)
        line_ang(i) = ang
    endif
enddo

end subroutine calc_line_angle



subroutine find_edge_node_lv(nn, node_lv, ln, line, line_ang, cri_ang)
    
implicit none

integer, intent(in) :: nn, ln, line(4,ln)
real, intent(in) :: line_ang(ln), cri_ang
integer, intent(inout) :: node_lv(nn)

integer :: i, num, count, temp_ned, ned
integer :: temp_edge(ln), this(4), edge(ln)
integer, allocatable :: edge_node(:)

temp_ned = 0
do i = 1, ln
    if (line_ang(i) >= cri_ang) then
        call find_count(4, this, 2, line(1:2,i), count)
        if (count /= 2) then
            temp_ned = temp_ned + 1
            temp_edge(temp_ned) = i
        endif
    endif
enddo

call delete_alone_line(ln, line, temp_ned, temp_edge, ned, edge)
allocate(edge_node(ned*2))
do i = 1, ned
    !write (*,*) 'edge sp, lp:', line(1:2,edge(i))
    edge_node(i*2-1:i*2) = line(1:2,edge(i))
enddo
call descending_sort_int(ned*2, edge_node)

node_lv = 1
num = 0
count = 1
do i = 1, ned*2
    if (num /= edge_node(i)) then
        if (count >= 3) then
            node_lv(num) = 3
        elseif (num /= 0) then
            node_lv(num) = 2
        endif
        num = edge_node(i)
        count = 1
    else
        count = count + 1
        if (i == ned*2 .AND. count >= 3) then
            node_lv(num) = 3
        elseif (i == ned*2 .AND. count < 3) then
            node_lv(num) = 2
        endif
    endif
enddo

end subroutine find_edge_node_lv


subroutine delete_alone_line(ln, line, inp_ned, inp_edge, ned, edge)

implicit none

integer, intent(in) :: ln, inp_ned, line(4,ln), inp_edge(ln)
integer, intent(inout) :: ned, edge(ln)

integer :: i, j, k, temp_ned, del_num, count, num
integer :: temp_edge(ln), target_edge(2), ball_edge(2)

temp_ned = inp_ned
temp_edge = 0
temp_edge(1:temp_ned) = inp_edge(1:temp_ned)
delete_line: do 
    ned = 0
    edge = 0
    del_num = 0
    
    do i = 1, temp_ned
        num = 0
        target_edge = line(1:2,temp_edge(i))
        do j = 1, 2
            do k = 1, temp_ned
                if (i /= k) then
                    ball_edge = line(1:2,temp_edge(k))
                    call find_count(1, target_edge(j), 2, ball_edge, count)
                    if (count == 1) then
                        num = num + 1
                        exit
                    endif
                endif
            enddo
        enddo

        if (num == 2) then
            ned = ned + 1
            edge(ned) = temp_edge(i)
        else
            del_num = del_num + 1
        endif
    enddo
    if (del_num == 0) then
        exit delete_line
    else
        temp_ned = ned
        temp_edge = 0
        temp_edge(1:ned) = edge(1:ned)
    endif    
enddo delete_line
    
end subroutine delete_alone_line
