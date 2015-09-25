subroutine find_other_two_node(line, ele, o_node)

implicit none

integer, intent(in) :: line(2), ele(3,2)
integer, intent(inout) :: o_node(2)

integer :: i, j

do i = 1, 2
    do j = 1, 3
        if (ele(j,i) /= line(1) .AND. ele(j,i) /= line(2)) then
            o_node(i) = ele(j,i)
            exit
        endif
    enddo
enddo

end subroutine find_other_two_node

!****************************************************************************

subroutine find_inner_line(ele_nodes, ne, elem, ln, line, inner_line)

implicit none
integer, intent(in) :: ele_nodes, ne, elem(ele_nodes,ne), ln
integer, intent(inout) :: inner_line(ln), line(2,ln)

integer :: i, j, k, q, o, nn, mesh_type, count, cp, sn, num
integer :: face_num, pre_ele, lp, j_r, iter, sfn, temp_k
integer :: this(3), temp(2), temp_line_face(4), elem_line(2,6), face_nn(3,4)
logical :: check
integer, allocatable :: line_hn(:), line_tfn(:), line_face(:,:,:)
real, allocatable :: node(:,:)

do i = 1, ln
    if (line(1,i) > line(2,i)) then
        temp(1) = line(1,i)
        line(1,i) = line(2,i)
        line(2,i) = temp(1)
    endif
enddo
mesh_type = 6
face_num = 3
! find line_tn, line_tfn
allocate (line_hn(ln), line_tfn(ln), line_face(30,face_num,ln))
line_hn = 0
line_tfn = 0
line_face = 0
do i = 1, ne
    call extract_line_from_elem(ele_nodes, mesh_type, elem(:,i), elem_line)
    call make_face_nn(ele_nodes, face_num, elem(:,i), face_nn)
    do j = 1, mesh_type
        if (elem_line(1,j) > elem_line(2,j)) then
            temp(1) = elem_line(1,j)
            elem_line(1,j) = elem_line(2,j)
            elem_line(2,j) = temp(1)
        endif
        do k = 1, ln
            if (elem_line(1,j) == line(1,k)) then
                call find_count(2, elem_line(:,j), 2, line(:,k), count)
                if (count == 2) then
                    line_tfn(k) = line_tfn(k) + 1
                    exit  ! (k)
                endif
            endif
        enddo  ! (k)
    enddo
    
    iter = ele_nodes/2+2
    do j = 1, iter
        do k = 1, face_num
            if (k == face_num) then
                this(1:2) = (/ face_nn(k,j), face_nn(1,j) /)
            else
                this(1:2) = (/ face_nn(k,j), face_nn(k+1,j) /)
            endif
            if (this(1) > this(2)) then
                temp(1) = this(1)
                this(1) = this(2)
                this(2) = temp(1)
            endif
            do q = 1, ln
                if (this(1) == line(1,q)) then
                    call find_count(2, this, 2, line(:,q), count)
                    if (count == 2) then
                        check = .true.
                        do o = 1, line_hn(q)
                            temp_line_face(1:face_num) = line_face(o,:,q)
                            call find_count(face_num, temp_line_face, face_num, face_nn(:,j),count)
                            if (count == face_num) then
                                check = .false.
                                exit
                            endif
                        enddo
                        if (check) then
                            line_hn(q) = line_hn(q) + 1
                            line_face(line_hn(q),:,q) = face_nn(:,j)
                            exit
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
enddo
deallocate (line_face) 
    
! find inner line
inner_line = 1
do i = 1, ln
    if ( line_hn(i) == line_tfn(i) ) then
        inner_line(i) = 0
    endif
enddo
deallocate (line_hn, line_tfn)

end subroutine find_inner_line

subroutine tecplot_edge(nn, node, ln, line, ned, edge)

implicit none

integer :: nn, ne, ln, line(4,ln), ned, edge(ln)
real :: node(3,nn)

integer :: i, j, status

open (Unit=20, File='./output/solid/remesh/check_edge.plt', STATUS='replace', ACTION='write', IOSTAT=status)
Write (20,'(A)') 'TITLE="3D surface check!!!"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6, A)') 'ZONE N=', nn , ', E=', ned, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'        

do i = 1, nn
    write(20,*) node(:,i)
enddo
do i = 1, ned
    write (20,'(3I7)') line(1:2,edge(i))
enddo

end subroutine tecplot_edge

subroutine tecplot_point(nn, node, ne, elem, bnn, bound_node)

implicit none

integer :: nn, ne, elem(3,ne), bnn, bound_node(bnn)
real :: node(3,nn)

integer :: i, j, status, num, memory_i

open (Unit=20, File='./output/solid/remesh/check_point.plt', STATUS='replace', ACTION='write', IOSTAT=status)
Write (20,'(A)') 'TITLE="3D surface check!!!"'
Write (20,*) 'VARIABLES="x", "y", "z", "b_node"'
Write (20,'(A,I6,A,I6, A)') 'ZONE N=', nn , ', E=', ne, ', DATAPACKING = POINT, ZONETYPE = FETRIANGLE'        
!memory_i = 1
do i = 1, nn
    num = 0
    do j = 1, bnn
        if (bound_node(j) == i) then
            num = 1
            exit
        endif
    enddo
    write(20,*) node(:,i), num
enddo
do i = 1, ne
    write (20,'(3I7)') elem(:,i)
enddo

end subroutine tecplot_point


subroutine tecplot_face(nn, node, ne, elem, bfn, temp_bfn, bound_fn, bound_face, flag)

implicit none

integer :: nn, ne, elem(3,ne), bfn, temp_bfn, bound_fn(bfn), bound_face(temp_bfn,bfn), flag
real :: node(3,nn)

integer :: i, j, k, status, num, memory_i

if (flag == 1) then
    open (Unit=20, File='./output/solid/remesh/check_face1.plt', STATUS='replace', ACTION='write', IOSTAT=status)
else
    open (Unit=20, File='./output/solid/remesh/check_face2.plt', STATUS='replace', ACTION='write', IOSTAT=status)
endif
    
Write (20,'(A)') 'TITLE="3D surface check!!!"'
Write (20,*) 'VARIABLES="x", "y", "z"'
do k = 1, bfn
    Write (20,'(A,I4,A,I6,A,I6, A)') 'ZONE T = "face_num', k, '", N=', nn , ', E=', bound_fn(k), ', DATAPACKING = POINT, ZONETYPE = FETRIANGLE'        
    do i = 1, nn
        write(20,*) node(:,i)
    enddo
    do j = 1, bound_fn(k)
        write (20,'(3I7)') elem(:,bound_face(j,k))
    enddo
enddo

end subroutine tecplot_face

subroutine tecplot_surface(filename, nn, ne, node, elem) 

implicit none

integer :: nn, ne, elem(3,ne)
real :: node(3,nn)
character(len=40) :: filename 

integer :: i, status

open (Unit=20, File=filename, STATUS='replace', ACTION='write', IOSTAT=status)
Write (20,'(A)') 'TITLE="3D surface check!!!"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I6,A,I6, A)') 'ZONE N=', nn , ', E=', ne, ', DATAPACKING = POINT, ZONETYPE = FETRIANGLE'        
do i = 1, nn
    write(20,*) node(:,i)
enddo
do i = 1, ne
    write (20,'(3I7)') elem(:,i)
enddo
close(20)

end subroutine tecplot_surface


subroutine adjust_boundary_info(ln, line, np2, p2, tcen, cen, control_edge, ned, edge, control_edge2, max_conn, p2cn, p2_conn)

implicit none

integer, intent(in) :: ln, line(4,ln), np2, tcen, cen(ln), control_edge(ln,200), max_conn
integer, intent(inout) :: p2(2,ln), ned, edge(ln), control_edge2(ln,200), p2cn(np2), p2_conn(max_conn,np2)

integer :: i, j, k, this(2), this2(2), count

ned = 0
p2cn = 0
do i = 1, tcen
    do j = 1, cen(i)-1
        this = control_edge(j:j+1,i)
        do k = 1, ln
            this2 = line(1:2,k)
            call find_count(2, this, 2, line(1:2,k), count)
            if (count == 2) then
                ned = ned + 1
                edge(ned) = k
                control_edge2(j,i) = k
                exit
            endif
        enddo
    enddo
    do j = 1, 2
        if (j == 1) then
            this = (/ control_edge(1,i), control_edge(2,i) /)
        else
            this = (/ control_edge(cen(i),i), control_edge(cen(i)-1,i) /)
        endif
        do k = 1, np2
            if (p2(1,k) == this(1)) then
                p2cn(k) = p2cn(k) + 1
                p2_conn(p2cn(k),k) = this(2)
                exit
            endif
        enddo
    enddo
enddo
do i = 1, np2
    p2(2,i) = p2cn(i)
enddo

end subroutine adjust_boundary_info