subroutine find_edge(ln, line, line_ang, ned, np2, edge, p2, cri_ang)
    
implicit none

integer, intent(in) :: ln, line(4,ln)
real, intent(in) :: line_ang(ln), cri_ang
integer, intent(inout) :: ned, np2, edge(ln), p2(2,ln)

integer :: i, num, count, temp_ned
integer :: temp_edge(ln), this(4)
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

num = 0
count = 1
do i = 1, ned*2
    if (num /= edge_node(i)) then
        if (count >= 3) then
            np2 = np2 + 1
            p2(1,np2) = num
            p2(2,np2) = count
            !write (*,*) 'p2:', p2(:,np2)
        endif
        num = edge_node(i)
        count = 1
    else
        count = count + 1
        if (i == ned*2 .AND. count >= 3) then
            np2 = np2 + 1
            p2(1,np2) = num
            p2(2,np2) = count
            !write (*,*) 'p2:', p2(:,np2)
        endif
    endif
enddo

end subroutine find_edge


subroutine remesh_find_edge(ln, line, line_ang, nn, node, ne, elem, ned, np2, p2, edge, cri_ang, check_edge)

implicit none

integer, intent(in) :: ln, line(4,ln), np2, nn, ne, elem(3,ne)
real, intent(in) :: line_ang(ln), cri_ang, node(3,nn)
integer, intent(inout) :: ned, edge(ln), p2(2,np2)
logical, intent(out) :: check_edge

integer :: i, j, k, temp_k, temp_ned, count, num, seq, pre_num, conn_edge_num, status
integer :: temp_edge(ln), temp_line(2), slp(2), temp_p2(np2), this(ln)
integer :: np2_near(2), excute_p2_num(np2), excute_p2(5,np2), conn_edge(2,6)
real :: p(3,3), ang, max_ang
logical :: check, check_np2
integer, allocatable :: edge_node(:)

temp_p2 = p2(1,:)
temp_ned = 0
excute_p2 = 0
excute_p2_num = 0
ned = 0
edge = 0
p2(2,:) = 0
do i = 1, ln
    if (line_ang(i) >= cri_ang) then
        temp_ned = temp_ned + 1
        temp_edge(temp_ned) = i
    endif
enddo

do i = 1, temp_ned
    temp_line = line(1:2,temp_edge(i))
    slp = 0
    check_np2 = .false.
    do j = 1, np2
        if (temp_line(1) == temp_p2(j)) then
            slp = temp_line
            check_np2 = .true.
        elseif (temp_line(2) == temp_p2(j)) then
            slp = (/ temp_line(2), temp_line(1) /)
            check_np2 = .true.
        endif
        
        check = .true.
        do k = 1, excute_p2_num(j)
            if (slp(2) == excute_p2(k,j)) then
                check = .false.
                exit
            endif
        enddo

        if (slp(1) /= 0 .and. check) then
            num = 1
            this(1) = temp_edge(i)
            np2_near = (/ j, slp(2) /)
            pre_num = i
            find_another_np2: do 
                call check_lp_np2(slp(2), np2, temp_p2, seq)
                if (seq /= 0) then
                    excute_p2_num(np2_near(1)) = excute_p2_num(np2_near(1)) + 1
                    excute_p2(excute_p2_num(np2_near(1)),np2_near(1)) = np2_near(2)
                    excute_p2_num(seq) = excute_p2_num(seq) + 1
                    excute_p2(excute_p2_num(seq),seq) = slp(1)
                    
                    edge(ned+1:ned+num) = this(1:num)
                    ned = ned + num
                    exit find_another_np2
                endif
                
                conn_edge = 0
                conn_edge_num = 0
                do k = 1, temp_ned
                    if (k /= pre_num) then
                        temp_line = line(1:2,temp_edge(k))
                        if (temp_line(1) == slp(2)) then
                            conn_edge_num = conn_edge_num+ 1
                            conn_edge(:,conn_edge_num) = (/ temp_line(2), k /)
                        elseif (temp_line(2) == slp(2)) then
                            conn_edge_num = conn_edge_num+ 1
                            conn_edge(:,conn_edge_num) = (/ temp_line(1), k /)
                        endif
                    endif
                enddo
                if (conn_edge_num == 0) then
                    exit find_another_np2
                else
                    max_ang = 0.0
                    p(:,1) = node(:,slp(1))
                    p(:,2) = node(:,slp(2))
                    do k = 1, conn_edge_num
                        p(:,3) = node(:,conn_edge(1,k))
                        call calc_angle_cos(p(:,1), p(:,2), p(:,3), ang)
                        if (max_ang < ang) then
                            max_ang = ang
                            temp_k = k
                        endif
                    enddo
                    slp = (/ slp(2), conn_edge(1,temp_k) /)
                    num = num + 1
                    this(num) = temp_edge(conn_edge(2,temp_k))
                    pre_num = conn_edge(2,temp_k)
                endif
            enddo find_another_np2
        endif
        if (check_np2) exit
    enddo
enddo

check_edge = .true.
do i = 1, np2
    p2(2,i) = excute_p2_num(i)
    if (excute_p2_num(i) < 3) then
        check_edge = .false.
        exit
    endif
enddo

if (check_edge) then
    allocate(edge_node(ned*2))
    do i = 1, ned
        edge_node(i*2-1:i*2) = line(1:2,edge(i))
    enddo
    call descending_sort_int(ned*2, edge_node)
    open (Unit=20, File='check_remesh_edge.plt', STATUS='replace', ACTION='write', IOSTAT=status)
    Write (20,'(A)') 'TITLE="3D mesh generation!!!"'
    Write (20,*) 'VARIABLES="x", "y", "z", "edge"'
    Write (20,'(A,I6,A,I6, A)') 'ZONE N=', nn , ', E=', ne, ', DATAPACKING = POINT, ZONETYPE = FETRIANGLE'        
    do i = 1, nn
        num = 0
        do j = 1, ned*2
            if (edge_node(j) == i) then
                num = 1
                exit
            endif
        enddo
        Write (20,'(3f13.8,2x,I2)') node(:,i), num
    enddo
    do i = 1, ne
        write (20,*) elem(:,i)
    enddo
endif

end subroutine remesh_find_edge

subroutine remesh_find_edge_bak(ln, line, line_ang, ned, np2, edge, p2, tri_np2, tri_p2, cri_ang, check)
    
implicit none

integer, intent(in) :: ln, line(4,ln), tri_np2, tri_p2(tri_np2)
real, intent(in) :: line_ang(ln), cri_ang
integer, intent(inout) :: ned, np2, edge(ln), p2(2,ln)
logical, intent(inout) :: check

integer :: i, num, count, temp_ned
integer :: temp_edge(ln), this(4)
integer, allocatable :: edge_node(:), temp_p2(:)

this(1:2) = (/ 2045, 2047 /)
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

num = 0
count = 1
do i = 1, ned*2
    if (num /= edge_node(i)) then
        if (count >= 3) then
            np2 = np2 + 1
            p2(1,np2) = num
            p2(2,np2) = count
            !write (*,*) 'p2:', p2(:,np2)
        endif
        num = edge_node(i)
        count = 1
    else
        count = count + 1
        if (i == ned*2 .AND. count >= 3) then
            np2 = np2 + 1
            p2(1,np2) = num
            p2(2,np2) = count
            !write (*,*) 'p2:', p2(:,np2)
        endif
    endif
enddo
allocate (temp_p2(np2))
temp_p2 = p2(1,1:np2)
write (*,*) 'np2:', np2
write (*,*) temp_p2

!check = .FALSE.
!if (tri_np2 > 1) then
!    if (tri_np2 == np2) then
!        call find_count(tri_np2, tri_p2, np2, temp_p2, count)
!        if (count == tri_np2) then
!            check = .TRUE.
!        endif
!    endif
!else
    check = .TRUE.
!endif

end subroutine remesh_find_edge_bak


subroutine check_lp_np2(lp, np2, p2, seq)
        
implicit none
        
integer, intent(in) :: lp, np2, p2(np2)
integer, intent(out) :: seq

integer :: i

seq = 0
do i= 1, np2
    if (p2(i) == lp) then
        seq = i
        exit
    endif
enddo

end subroutine check_lp_np2