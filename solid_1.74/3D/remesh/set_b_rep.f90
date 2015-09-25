subroutine set_b_rep(ng, cur_ng, nn, node, ne, elem, add_nn, add_coord, ln, line, line_ang, np2, p2, ned, edge, tcen, cen, control_edge, fsi_flag)
                    
use remesh_domain
    
implicit none

integer, intent(in) :: ng, cur_ng, nn, ne, ln, fsi_flag, add_nn
integer, intent(in) :: elem(3,ne), line(4,ln)
real, intent(in) :: node(3,nn), line_ang(ln), add_coord(3,20)
integer, intent(inout) :: np2, p2(2,ln), ned, edge(ln), tcen, cen(ln), control_edge(ln,200)

integer :: i, j, k, q, o1, o2, temp_k, sp, lp, pp, cl, ten, max_conn, num, a, b, status, memory_i, flag, pre_np2
integer :: temp_ele, current_ele, temp_e, current_edge, seq, temp(5), except, check_num, temp_node, pre_tcen, count
integer :: e_num(2), p2cn(np2), temp_edge(ln), temp_edge2(ln), edge_check(ln), this(4), temp_elem(3,2)
integer :: control_edge2(ln,200), other_node(2), except_ele(ne)
integer :: temp_tcen, temp_cen(ln), temp_control_edge(ln,200)
real :: temp_pl(4), current_pl(4), p(3,3), val(2), temp_val, op(3)
integer,allocatable :: p2_conn(:,:), p2_conn_use(:,:), temp_p2(:)
logical :: check, corner_check
character(len=35) :: fn

integer :: bnn, bln, bfn, euler_bfn, temp_bfn
integer, allocatable :: bound_fn(:), bound_face(:,:), ce_check(:), bound_pl(:,:), bound_pln(:)
integer, allocatable :: bound_contain_num(:,:), bound_contain_node(:,:,:), temp_pln(:), temp_pl_node(:,:)

! manually add variable
integer :: add_np2
!----------------------------------------------------------------------------
op = (/ 0.0, 0.0, 0.0 /)
max_conn = 6
allocate(p2_conn(max_conn,np2), p2_conn_use(max_conn,np2))
p2_conn = 0;  p2cn = 0;   p2_conn_use = 0
control_edge2 = 0
edge_check = 0
write (*,*) 'fsi_flag:', fsi_flag
if (fsi_flag == 1) then
    if (np2 == 0) then  ! control point가 0인 경우
        
    else
        do i = 1, np2
            do j = 1, p2(2,i)
                sp = 0;  lp = 0
                ten = 0 ; temp_edge = 0;  temp_edge2 = 0
                do k = 1, ned
                    e_num = line(1:2,edge(k))
                    if (e_num(1) == p2(1,i) .OR. e_num(2) == p2(1,i)) then
                        cl = k
                        edge_check(k) = 1
                        sp = p2(1,i)
                        if (e_num(1) == p2(1,i)) then
                            lp = e_num(2)
                        else
                            lp = e_num(1)
                        endif
                    
                        check = .TRUE.
                        do q = 1, p2cn(i)
                            if (p2_conn(q,i) == lp) then
                                check = .FALSE.
                                exit
                            endif
                        enddo
                    
                        if (check) then
                            ten = 1
                            temp_edge(ten) = sp
                            temp_edge2(ten) = edge(k)
                            p2cn(i) = p2cn(i) + 1
                            p2_conn(p2cn(i),i) = lp
                            exit
                        else
                            lp = 0
                        endif
                    endif
                enddo

                if (lp /= 0) then
                    !write (*,*) 'sp, lp:', sp, lp
                    a = 0
                    do 
                        a = a +1
                        if (a == 1000) stop 'loop over(a)'
                        check = .TRUE.
                        do k = 1, np2
                            if (p2(1,k) == lp) then
                                ten = ten + 1
                                temp_edge(ten) = lp
                                p2cn(k) = p2cn(k) + 1
                                p2_conn(p2cn(k),k) = temp_edge(ten-1)
                                check = .FALSE.
                                exit
                            endif
                        enddo
                    
                        if (check) then
                            do k = 1, ned
                                if (cl /= k) then
                                    e_num = line(1:2,edge(k))
                                    if (e_num(1) == lp .OR. e_num(2) == lp) then
                                        !write (*,*) 'e_num:', e_num(1:2)
                                        edge_check(k) = 1
                                        cl = k
                                        ten = ten + 1
                                        temp_edge(ten) = lp
                                        temp_edge2(ten) = edge(k)
                                        if (e_num(1) == lp) then
                                            lp = e_num(2)
                                        else
                                            lp = e_num(1)
                                        endif
                                        exit
                                    endif
                                endif
                            enddo
                        else
                            tcen = tcen + 1
                            cen(tcen) = ten
                            control_edge(1:ten,tcen) = temp_edge(1:ten)
                            control_edge2(1:ten-1,tcen) = temp_edge2(1:ten-1)
                            !write (*,*) tcen, 'edge point number:', ten
                            !write (*,*) 'edge point:', temp_edge(1:ten)
                            !write (*,*) 'edge line:', temp_edge2(1:ten-1)
                            !write (*,*) 
                            exit
                        endif
                    enddo
                endif
            enddo
        enddo
    endif
    pre_tcen = tcen
    open (unit=37, file='edge.dat', status='REPLACE', action='WRITE')
    write (37,*) np2
    write (37,'(100I8,1X)') p2(1,1:np2)
    write (37,*) tcen
    write (37,'(100I8,1X)') cen(1:tcen)
    do i = 1, tcen
        write (37,'(/,(10I9))') control_edge(1:cen(i),i)
    enddo
    close(37)
    call tecplot_point(nn, node, ne, elem, np2, p2(1,1:np2))
    stop
elseif (fsi_flag == 2) then
    write (*,*) '>> adjust_boundary_info'
    call adjust_boundary_info(ln, line, np2, p2, tcen, cen, control_edge, ned, edge, control_edge2, max_conn, p2cn, p2_conn)
    Hexa(cur_ng)%np2 = np2
    Hexa(cur_ng)%ben = tcen
    allocate (Hexa(cur_ng)%bp2(np2), Hexa(cur_ng)%bp2_coord(3,np2), Hexa(cur_ng)%bedn(tcen))
    allocate (Hexa(cur_ng)%cen(tcen), Hexa(cur_ng)%control_edge_number(tcen,tcen))
    Hexa(cur_ng)%cen = 1
    do i = 1, tcen
        Hexa(cur_ng)%control_edge_number(1,i) = i
    enddo
    pre_tcen = tcen
endif
do i = 1, tcen
    !write (*,'(/,(10I9))') control_edge(1:cen(i),i)
    write (*,'(2I9)') control_edge(1,i), control_edge(cen(i),i)
enddo
do i = 1, np2
    write (*,*) 'p2cn:', p2cn(i)
    write (*,*) p2(1,i), ': ', p2_conn(1:p2cn(i),i)
enddo

call tecplot_edge(nn, node, ln, line, ned, edge)
!=========================================================================================

bfn = 0
euler_bfn = 2 + tcen - np2
temp_bfn = ne/2
allocate (bound_fn(euler_bfn), bound_face(temp_bfn,euler_bfn), ce_check(tcen))
allocate (bound_pln(euler_bfn), bound_pl(np2*4,euler_bfn))
bound_fn = 0;  bound_face = 0
bound_pln = 0;  bound_pl = 0

ce_check = 0
memory_i = 1
b = 0
! p2_conn_use = 0: don`t use
!             = 1: use
do
    if (bfn == euler_bfn) exit
    i = memory_i
    ! p2에서 first & next point 설정
    sp = p2(1,i)
    current_edge = 0
    if (bfn == 0) then
        check_num = 0
    else
        check_num = 2
    endif
    p_conn: do j = 1, p2cn(i)
        if (p2_conn_use(j,i) == check_num) then
            lp = p2_conn(j,i)
            do k = 1, tcen
                if (ce_check(k) /= 2) then
                    if (control_edge(1,k) == sp .AND. control_edge(2,k) == lp) then
                        current_edge = k
                        temp_edge(1:cen(k)) = control_edge(1:cen(k),k)
                        temp_edge2(1:cen(k)-1) = control_edge2(1:cen(k)-1,k)
                        p2_conn_use(j,i) = p2_conn_use(j,i) + 1
                        exit p_conn !p_conn
                    elseif (control_edge(cen(k),k) == sp .AND. control_edge(cen(k)-1,k) == lp) then
                        current_edge = k
                        num = 0
                        do q = cen(k), 1, -1
                            num = num + 1
                            temp_edge(num) = control_edge(q,k)
                            if (q /= 1) temp_edge2(num) = control_edge2(q-1,k)
                        enddo
                        p2_conn_use(j,i) = p2_conn_use(j,i) + 1
                        exit p_conn !p_conn
                    endif
                endif
            enddo
        endif
    enddo p_conn
    
    ! simple closed curve 검색
    if (current_edge /= 0) then
        bfn = bfn + 1
        current_ele = line(3,temp_edge2(1))
        bfn_check: do j = 1, bfn-1
            do k = 1, bound_fn(j)
                if (bound_face(k,j) == current_ele) then
                    current_ele = line(4,temp_edge2(1))
                    exit bfn_check
                endif
            enddo
        enddo bfn_check
        bound_fn(bfn) = bound_fn(bfn) + 1
        bound_face(bound_fn(bfn),bfn) = current_ele
        bound_pln(bfn) = bound_pln(bfn) + 1
        bound_pl(bound_pln(bfn),bfn) = sp
        !write (*,*) 'bfn:', bfn, ', sp, lp, ce:', sp, lp, current_ele
        !write (*,*) 'line:',line(:,temp_edge2(1))
    
        p(:,1) = node(:,elem(1,current_ele))
        p(:,2) = node(:,elem(2,current_ele))
        p(:,3) = node(:,elem(3,current_ele))
        call calc_plane(p(:,1), p(:,2), p(:,3), current_pl)
        a = 0
        scc: do  
            a = a + 1
            ce_check(current_edge) = ce_check(current_edge) + 1
            !write (*,*) cen(current_edge)
            do j = 2, cen(current_edge)
                check = .true.
                p_conn2: do k = 1, np2
                    if (p2(1,k) == lp) then
                        temp_k = k
                        check = .false.
                        do q = 1, p2cn(k)
                            if (p2_conn(q,k) == pp) then
                                p2_conn_use(q,k) = p2_conn_use(q,k) + 2
                                except = q
                                exit p_conn2
                            endif
                        enddo
                    endif
                enddo p_conn2
                
                if (lp == sp) then
                    if (bound_face(1,bfn) == bound_face(bound_fn(bfn),bfn)) bound_fn(bfn) = bound_fn(bfn) - 1
                    !write (*,*) bfn, ' boundary face'
                    !write (*,*) bound_face(1:bound_fn(bfn),bfn)
                    exit scc
                endif
                if (check) then   ! lp가 p2가 아닐일 때
                    num = temp_edge2(j)
                    !write (*,*) line(1:2,num)
                    temp_elem = elem(:,line(3:4,num))
                    call find_other_two_node(line(1:2,num), temp_elem, other_node)
                    val = 0.
                    !do k = 1, 2
                    !    do q = 1, 3
                    !        val(k) = val(k) + current_pl(q)*node(q,other_node(k))
                    !    enddo
                    !    val(k) = val(k) + current_pl(4)
                    !    write (*,*) other_node(k), val(k)
                    !enddo
                    do k = 1, 2
                        p(:,1) = node(:,elem(1,line(2+k,num)))
                        p(:,2) = node(:,elem(2,line(2+k,num)))
                        p(:,3) = node(:,elem(3,line(2+k,num)))
                        call calc_plane(p(:,1), p(:,2), p(:,3), temp_pl)
                        call calc_angle_cos(temp_pl(1:3), op, current_pl(1:3), val(k))
                    enddo
                    if (abs(val(1)) > abs(val(2))) then
                        temp_ele = line(4,num)
                        temp_node = other_node(2)
                    else
                        temp_ele = line(3,num)
                        temp_node = other_node(1)
                    endif       
                    !write (*,*) lp, temp_node
                    bound_fn(bfn) = bound_fn(bfn) + 1
                    bound_face(bound_fn(bfn),bfn) = temp_ele
                    p(:,1) = node(:,elem(1,temp_ele))
                    p(:,2) = node(:,elem(2,temp_ele))
                    p(:,3) = node(:,elem(3,temp_ele))
                    call calc_plane(p(:,1), p(:,2), p(:,3), current_pl)
                    pp = lp
                    lp = temp_edge(j+1)
                else ! lp가 p2일 때, 선분이 꺽임.
                    pp = lp
                    bound_pln(bfn) = bound_pln(bfn) + 1
                    bound_pl(bound_pln(bfn),bfn) = lp
                    if (j == 2) then
                        num = temp_edge2(1)
                        temp_elem = elem(:,line(3:4,num))
                        call find_other_two_node(line(1:2,num), temp_elem, other_node)
                        val = 0.
                        do k = 1, 2
                            p(:,1) = node(:,elem(1,line(2+k,num)))
                            p(:,2) = node(:,elem(2,line(2+k,num)))
                            p(:,3) = node(:,elem(3,line(2+k,num)))
                            call calc_plane(p(:,1), p(:,2), p(:,3), temp_pl)
                            call calc_angle_cos(temp_pl(1:3), op, current_pl(1:3), val(k))
                        enddo
                        if (abs(val(1)) > abs(val(2))) then
                            temp_ele = line(4,num)
                            temp_node = other_node(2)
                        else
                            temp_ele = line(3,num)
                            temp_node = other_node(1)
                        endif       
                    endif
                    ! temp - 1: temp_ele, 2: temp_lp, 3: temp_edge, 4: seq, 5: p2_conn_use
                    call find_next_point(p2cn(temp_k), p2_conn(1:p2cn(temp_k),temp_k), current_pl, nn, node, ne, elem, temp_ele, pp, temp_node, temp)
                    do q = 1, tcen
                        if (ce_check(q) /= 2) then
                            if (control_edge(1,q) == pp .AND. control_edge(2,q) == temp(2)) then
                                temp(3:4) = (/ q, 1 /)  ! 정방향
                                exit
                            elseif (control_edge(cen(q),q) == pp .AND. control_edge(cen(q)-1,q) == temp(2)) then
                                temp(3:4) = (/ q, 2 /)  ! 역방향
                                exit
                            endif
                        endif
                    enddo
                    if (temp(4) == 1) then
                        num = control_edge2(1,temp(3))
                    else
                        num = control_edge2(cen(temp(3))-1,temp(3))
                    endif
                    p(:,1) = node(:,elem(1,temp(1)))
                    p(:,2) = node(:,elem(2,temp(1)))
                    p(:,3) = node(:,elem(3,temp(1)))                    
                    call calc_plane(p(:,1), p(:,2), p(:,3), current_pl)
                    lp = temp(2)
                    current_edge = temp(3)
                    if (bound_face(bound_fn(bfn),bfn) /= temp(1)) then
                        bound_fn(bfn) = bound_fn(bfn) + 1
                        bound_face(bound_fn(bfn),bfn) = temp(1)
                    endif
                    if (temp(4) == 1) then
                        temp_edge(1:cen(temp(3))) = control_edge(1:cen(temp(3)),temp(3))
                        temp_edge2(1:cen(temp(3))-1) = control_edge2(1:cen(temp(3))-1,temp(3))
                    else 
                        num = 0
                        do k = cen(temp(3)), 1, -1
                            num = num + 1
                            temp_edge(num) = control_edge(k,temp(3))
                            if (k /= 1) temp_edge2(num) = control_edge2(k-1,temp(3))
                        enddo
                    endif
                    if (pp /= sp) p2_conn_use(temp(5),temp_k) = p2_conn_use(temp(5),temp_k) + 1
                    exit
                endif
            enddo
            !write (*,*) 'pp, lp:', pp, lp
            !write (*,*) 'temp:', temp
            !write (*,*) 
            if (a == 10) stop 'loop over(a)'
        enddo scc
    else
        memory_i = memory_i + 1
    endif
    b = b + 1
    if (b == 100) stop 'loop over(b)'
enddo
max_conn = maxval(bound_pln)
allocate (bound_contain_node(np2*2,max_conn,euler_bfn))
allocate (bound_contain_num(max_conn,euler_bfn))
allocate (temp_pln(euler_bfn), temp_pl_node(max_conn,euler_bfn))
bound_contain_node = 0
bound_contain_num = 0
temp_pln = bound_pln
temp_pl_node = 0
do i = 1, bfn
    temp_pl_node(1:bound_pln(i),i) = bound_pl(1:bound_pln(i),i)
enddo

pre_np2 = np2
add_np2 = remesh_add(ng)%np2
if (add_np2 /= 0) then
    p2(:,np2+1:np2+add_np2) = remesh_add(ng)%p2
endif
np2 = np2 + add_np2
temp_tcen = tcen
temp_cen = cen
temp_control_edge = control_edge
control_edge2 = 0
tcen = 0
cen = 0
control_edge = 0
Hexa(cur_ng)%cen = 0
Hexa(cur_ng)%control_edge_number = 0
do i = 1, remesh_add(ng)%ln
    !if (i == 69) write (*,*) 'remesh_line:', remesh_add(ng)%line(:,i) 
    sp = 0;  lp = 0;  flag = 0
    corner_check = .FALSE.
    do j = 1, temp_tcen
        call find_count(2, remesh_add(ng)%line(:,i), temp_cen(j), temp_control_edge(1:temp_cen(j),j), count)
        if (count == 2) then
            corner_check = .TRUE.
            flag = 3
            cl = j
            do k = 1, temp_cen(j)
                if (remesh_add(ng)%line(1,i) == temp_control_edge(k,j) .or. remesh_add(ng)%line(2,i) == temp_control_edge(k,j)) then
                    if (sp == 0) then
                        sp = k
                    else
                        lp = k
                        exit
                    endif
                endif
            enddo
            exit
        endif
    enddo
    if (corner_check == .FALSE.) then
        do j = 1, temp_tcen
            do k = 1, temp_cen(j)
                if (remesh_add(ng)%line(1,i) == temp_control_edge(k,j)) then
                    flag = 1
                    exit
                elseif (remesh_add(ng)%line(2,i) == temp_control_edge(k,j)) then
                    flag = 2
                    exit
                endif
            enddo
            if (flag /= 0) then
                cl = j
                exit
            endif
        enddo
    endif
    !write (*,'(A,5(I4,1X))') 'i, flag, cl, sp, lp:', i, flag, cl, sp, lp
    if (flag <= 2) then
        if (flag == 1 .or. flag == 2) then
            this(1:2) = (/ temp_control_edge(1,cl), temp_control_edge(temp_cen(cl),cl) /)
            call find_count(1, remesh_add(ng)%line(flag,i), 2, this(1:2), count)
            if (count == 1) flag = 0
        endif
        tcen = tcen + 1
        cen(tcen) = 2
        control_edge(1:2,tcen) = remesh_add(ng)%line(:,i)
    elseif (flag == 3) then
        this(1:2) = (/ temp_control_edge(1,cl), temp_control_edge(temp_cen(cl),cl) /)
        do k = 1, 2
            call find_count(1, remesh_add(ng)%line(k,i), 2, this(1:2), count)
            if (count == 1) flag = flag - k
        enddo

        tcen = tcen + 1
        cen(tcen) = lp - sp + 1
        control_edge(1:cen(tcen),tcen) = temp_control_edge(sp:lp,cl)
            
        Hexa(cur_ng)%cen(cl) = Hexa(cur_ng)%cen(cl) + 1
        Hexa(cur_ng)%control_edge_number(Hexa(cur_ng)%cen(cl),cl) = tcen
        control_edge2(Hexa(cur_ng)%cen(cl),cl) = 0      ! 0 = don`t include sp, lp 
                                                        ! 1 = include sp
                                                        ! 2 = include lp
                                                        ! 3 = include sp, lp
        if (flag == 0) then
            control_edge2(Hexa(cur_ng)%cen(cl),cl) = 3
        elseif (flag == 1 .or. flag == 2) then
            call find_count(2, remesh_add(ng)%line(:,i), 1, this(1), count)
            if (count == 1) then
                control_edge2(Hexa(cur_ng)%cen(cl),cl) = 1
            else
                control_edge2(Hexa(cur_ng)%cen(cl),cl) = 2
            endif
        endif
    endif
    write (*,'(I4,A,2I5,A,2(I7,2X))') tcen, ' tcen num, flag:', cen(tcen), flag, ', line_num:', control_edge(1,tcen), control_edge(cen(tcen),tcen)
    if (flag /= 0) then
        do q = 1, bfn
            do k = 1, bound_pln(q)
                if (k == bound_pln(q)) then
                    this(3:4) = (/ bound_pl(k,q), bound_pl(1,q) /)
                else
                    this(3:4) = (/ bound_pl(k,q), bound_pl(k+1,q) /)
                endif
                call find_count(2, this(1:2), 2, this(3:4), num)
                if (num == 2) then
                    !do o1 = bound_pln(q), k+1, - 1
                    !    bound_pl(o1+1,q) = bound_pl(o1,q)
                    !enddo
                    !bound_pl(k+1,q) = remesh_add(ng)%line(flag,i)
                    !bound_pln(q) = bound_pln(q) + 1
                    if (flag == 3) then
                        bound_contain_node(bound_contain_num(k,q)+1:bound_contain_num(k,q)+2,k,q) = remesh_add(ng)%line(:,i)
                        bound_contain_num(k,q) = bound_contain_num(k,q) + 2
                    else
                        bound_contain_node(bound_contain_num(k,q)+1,k,q) = remesh_add(ng)%line(flag,i)
                        bound_contain_num(k,q) = bound_contain_num(k,q) + 1
                    endif
                    exit
                endif
            enddo
        enddo
    endif
enddo

! rearrange bound surface point
bound_pln = 0
bound_pl = 0
do i = 1, bfn
    do j = 1, temp_pln(i)
        bound_pln(i) = bound_pln(i) + 1
        bound_pl(bound_pln(i),i) = temp_pl_node(j,i)
        if (bound_contain_num(j,i) /= 0) then
            memory_i = temp_pl_node(j,i)
            this(1) = temp_pl_node(j,i)
            if (j == temp_pln(i)) then
                this(2) = temp_pl_node(1,i)
            else
                this(2) = temp_pl_node(j+1,i)
            endif
            do k = 1, temp_tcen
                this(3:4) = (/ temp_control_edge(1,k), temp_control_edge(temp_cen(k),k) /)
                call find_count(2, this(1:2), 2, this(3:4), count)
                if (count == 2) then
                    num = k
                    if (this(1) == this(3)) then
                        sp = 2;  lp = temp_cen(k)-1;  a = 1
                    else
                        sp = temp_cen(k)-1;  lp = 2;  a = -1
                    endif
                    exit
                endif
            enddo
            do k = sp, lp, a
                do q = 1, bound_contain_num(j,i)
                    if (temp_control_edge(k,num) == bound_contain_node(q,j,i)) then
                        if (memory_i /= temp_control_edge(k,num)) then
                            bound_pln(i) = bound_pln(i) + 1
                            bound_pl(bound_pln(i),i) = temp_control_edge(k,num)
                            memory_i = temp_control_edge(k,num)
                        endif
                    endif
                enddo
            enddo
        endif
    enddo
    !write (*,*) i, ' bound_pln:', bound_pln(i)
    !write (*,*) bound_pl(1:bound_pln(i),i)
enddo
deallocate (bound_contain_node, bound_contain_num)
deallocate (temp_pln, temp_pl_node)

num = 0
do i = 1, tcen
    num = num + cen(i)-1
enddo
open (Unit=20, File='./output/solid/remesh/check_control_edge.plt', STATUS='replace', ACTION='write', IOSTAT=status)
Write (20,'(A)') 'TITLE="3D surface control_edge!!!"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I4,A,I6,A,I6, A)') 'ZONE N=', nn+add_nn , ', E=', num, ', DATAPACKING = POINT, ZONETYPE = FELINESEG' 
do i = 1, nn
    write(20,*) node(:,i)
enddo
do i = 1, add_nn
    write(20,*) add_coord(:,i)
enddo
do i = 1, tcen
    do j = 1, cen(i)-1
        write (20,'(2I7)') control_edge(j:j+1,i) 
    enddo
enddo
close(20)

call rearrange_hexa_b_edge(cur_ng, ln, tcen, cen, control_edge, temp_tcen, temp_cen, temp_control_edge, control_edge2)
!write (*,*) 'Hexa(ng)%ben:', Hexa(ng)%ben
!do i = 1, Hexa(ng)%ben
!    write (*,*) 'Hexa(ng)%cen', i, ':', Hexa(ng)%cen(i)
!    write (*,*) Hexa(ng)%control_edge_number(1:Hexa(ng)%cen(i),i)
!    write (*,*)
!enddo

call tecplot_point(nn, node, ne, elem, np2, p2(1,1:np2))
call tecplot_face(nn, node, ne, elem, bfn, temp_bfn, bound_fn, bound_face, 1)

except_ele = 0
do i = 1, bfn
    except_ele(bound_face(1:bound_fn(i),i)) = 1
enddo
do i = 1, bfn
    call find_face_ele(ne, elem, except_ele, temp_bfn, bound_fn(i), bound_face(:,i))
enddo
call tecplot_face(nn, node, ne, elem, bfn, temp_bfn, bound_fn, bound_face, 2)

!write (*,*) 'tcen:', tcen
!num = 0
!do i = 1, tcen
!    write (*,*) i, 'line num:', cen(i)
!    write (*,*) control_edge(1:cen(i),i)
!    !write (*,*) control_edge(1,i), control_edge(cen(i),i)
!    num = num + cen(i) - 1
!enddo

allocate (temp_p2(np2))
temp_p2 = p2(1,1:pre_np2)
call update_boundary_info(ng, pre_np2, temp_p2, add_nn, add_coord, ln, tcen, cen, control_edge, bfn, temp_bfn, bound_fn, bound_face, bound_pln, bound_pl)

deallocate(p2_conn, p2_conn_use, temp_p2)
deallocate (bound_fn, bound_face, ce_check)
deallocate (bound_pln, bound_pl)

end subroutine set_b_rep

!====================================================================================                

subroutine set_proj_pn(ng, bn, proj_pn, nn, bound_lv, ben, col, line_node_num, line_node)

use remesh_domain

implicit none

integer, intent(in) :: ng, bn, nn, ben, col, line_node_num(ben), line_node(col,ben)
integer, intent(inout) :: proj_pn(bn), bound_lv(nn)

integer :: i, j, k, q, temp_k, ccn, count, coe, fp, num_nci, exc_num, sp, lp, pre_sp, cnn, first_elem, inner_point
integer :: chk_line_num, end_line_num, end_elem_num, chk_num, chk_elem_num, inner_node_num, near_node_num, near_elem_num
integer :: iter
integer :: this(4), num(2), nci(8), except(4), line(2), elem(4), near_node(13), near_elem(6), temp(3)
integer :: chk_face(4,2), chk_point(2), chk_line(2,Hexa(ng)%nn), end_line(2,Hexa(ng)%nn*2), inner_node(Hexa(ng)%nn), conn(Hexa(ng)%nn)
integer, allocatable :: closed_curve(:), line_info(:,:), end_elem(:), chk_elem(:)
logical :: check
real :: p(3,3), pl(4,2), op(3), alpha(2), ang

!write (*,*) line_node_num
allocate (end_elem(Hexa(ng)%bfn), chk_elem(Hexa(ng)%bfn))
conn = Hexa(ng)%bound_conn
do i = 1, Tri%bfn
    write (*,*) '---- ', i, ' plane ------'
    write (*,*) Tri%bound_pl(1:Tri%bound_bpln(i),i)
    ccn = 0
    allocate(line_info(2,Tri%bound_bpln(i)))
    line_info = 0
    do j = 1, Tri%bound_bpln(i)
        if (j == Tri%bound_bpln(i)) then
            this(1:2) = (/ Tri%bound_pl(j,i), Tri%bound_pl(1,i) /)
        else
            this(1:2) = (/ Tri%bound_pl(j,i), Tri%bound_pl(j+1,i) /)
        endif
        do k = 1, Tri%ben
            if (line_node_num(k) /= 0) then
                this(3:4) = (/ Tri%b_edge(1,k), Tri%b_edge(Tri%bedn(k),k) /)
                call find_count(2, this(1:2), 2, this(3:4), count)
                if (count == 2) then
                    if (this(1) == this(3)) then
                        line_info(1,j) = 1   ! line: 3 -> 4
                    else
                        line_info(1,j) = 2   ! line: 4 -> 3
                    endif
                    line_info(2,j) = k
                    exit
                endif
            endif
        enddo
        if (line_info(2,j) == 0) then
            write (*,*) 'bound_pl:', this(1:2)
            stop
        endif            
        ccn = ccn + line_node_num(line_info(2,j)) - 1
    enddo
    
    allocate (closed_curve(ccn))
    count = 0
    do j = 1, Tri%bound_bpln(i)
        num(1) = line_node_num(line_info(2,j))
        if (line_info(1,j) == 1) then
            coe = 1
            sp = 1
            lp = num(1)-1
        else
            coe = -1
            sp = num(1)
            lp = 2
        endif
        do k = sp, lp, coe
            count = count + 1
            closed_curve(count) = line_node(k,line_info(2,j))
        enddo
    enddo
    !write (*,*) 'ccn, count:', ccn, count
    !write (*,*) 'closed curve:', closed_curve(1:count)
    bound_lv(closed_curve(1:count)) = 2
    end_line_num = ccn
    end_line = 0
    do j = 1, ccn
        if (j == ccn) then
            end_line(:,j) = (/ closed_curve(j), closed_curve(1) /)
        else
            end_line(:,j) = (/ closed_curve(j), closed_curve(j+1) /)
        endif
    enddo
    
    !fp = closed_curve(1)
    !sp = closed_curve(2)
    !lp = closed_curve(ccn)
    !this(1:2) = (/ line_info(2,1), line_info(2,Tri%bound_bpln(i)) /)
    !call find_nci(Hexa(ng)%bfn, Hexa(ng)%bound_face(:,1:Hexa(ng)%bfn), fp, num_nci, nci)
    !exc_num = 0
    !do j = 1, num_nci
    !    do k = 1, Tri%ben
    !        if (line_node_num(k) /= 0) then
    !            if (k /= this(1) .and. k /= this(2)) then
    !                if (nci(j) == line_node(2,k) .or. nci(j) == line_node(line_node_num(k)-1,k)) then
    !                    exc_num = exc_num + 1
    !                    except(exc_num) = nci(j)
    !                    exit
    !                endif
    !            endif
    !        endif
    !    enddo
    !enddo
    !
    !this(1) = fp
    !this(2) = closed_curve(2)
    !near_elem_num =  0
    !do j = 1, Hexa(ng)%bfn
    !    call find_count(1, this(1), 4, Hexa(ng)%bound_face(:,j), count)
    !    if (count == 1) then
    !        near_elem_num = near_elem_num + 1
    !        near_elem(near_elem_num) = j
    !    endif
    !enddo
    !
    !near_node_num = 3
    !do j = 1, near_elem_num
    !    elem = Hexa(ng)%bound_face(:,near_elem(j))
    !    call find_count(2, this(1:2), 4, elem, count)
    !    if (count == 2) then
    !        do k = 1, 4
    !            if (elem(k) == this(1)) then
    !                temp = (/ k+1, k+2, k+3 /)
    !                if (temp(1) > 4) temp(1) = temp(1) - 4
    !                if (temp(2) > 4) temp(2) = temp(2) - 4
    !                if (temp(3) > 4) temp(3) = temp(3) - 4
    !
    !                if (elem(temp(1)) == this(2)) then
    !                    near_node(1:3) = elem(temp)
    !                else
    !                    near_node(1:3) = (/ elem(temp(3)), elem(temp(2)), elem(temp(1)) /)
    !                endif
    !                exit
    !            endif
    !        enddo
    !        exit
    !    endif
    !enddo
    !
    !iter = 0
    !do 
    !    iter = iter + 1
    !    if (iter > near_node_num-1) then
    !        write (*,*) 'iteration is too much: ', iter
    !        stop 
    !    endif
    !    
    !    this(2) = near_node(near_node_num)
    !    this(3) = near_node(near_node_num-2)
    !    do j = 1, near_elem_num
    !        elem = Hexa(ng)%bound_face(:,near_elem(j))
    !        call find_count(2, this(1:2), 4, elem, count)
    !        if (count == 2) then
    !            call find_count(1, this(3), 4, elem, count)
    !            if (count == 0) then
    !                do k = 1, 4
    !                    if (elem(k) == this(1)) then
    !                        temp = (/ k+1, k+2, k+3 /)
    !                        if (temp(1) > 4) temp(1) = temp(1) - 4
    !                        if (temp(2) > 4) temp(2) = temp(2) - 4
    !                        if (temp(3) > 4) temp(3) = temp(3) - 4
    !                
    !                        if (elem(temp(1)) == this(2)) then
    !                            near_node(near_node_num+1:near_node_num+2) = elem(temp(2:3))
    !                        else
    !                            near_node(near_node_num+1:near_node_num+2) = (/ elem(temp(2)), elem(temp(1)) /)
    !                        endif
    !                        near_node_num = near_node_num + 2
    !                        exit
    !                    endif
    !                enddo
    !                exit
    !            endif
    !        endif
    !    enddo
    !    if (near_node(near_node_num) == near_node(1)) exit
    !enddo
    !write (*,*) 'exc_node:', except(1:exc_num)
    !write (*,*) 'near_node:', near_node(1:near_node_num-1)
    !do j = 2, near_node_num-1
    !    call find_count(1, near_node(j), exc_num, except(1:exc_num), count)
    !    if (count == 1) then
    !        this(1:3) = (/ fp, closed_curve(2), near_node(near_node_num-2) /)
    !        exit
    !    elseif (near_node(j) == closed_curve(ccn)) then
    !        this(1:3) = (/ fp, closed_curve(2), near_node(2) /)
    !        exit
    !    endif
    !enddo
    !
    !do k = 1, near_elem_num
    !    temp_k = near_elem(k)
    !    call find_count(3, this(1:3), 4, Hexa(ng)%bound_face(:,temp_k), count)
    !    if (count == 3) then
    !        first_elem = temp_k
    !        exit
    !    endif
    !enddo
    fp = closed_curve(1)
    sp = closed_curve(2)
    lp = closed_curve(ccn)
    p(:,1) = Hexa(ng)%node(:,fp)
    p(:,2) = Hexa(ng)%node(:,sp)
    p(:,3) = Hexa(ng)%node(:,lp)
    call calc_plane(p(:,1), p(:,2), p(:,3), pl(:,1))
    this(1:2) = (/ fp, sp /)
    num(1) = 0
    do j = 1, Hexa(ng)%bfn
        call find_count(2, this(1:2), 4, Hexa(ng)%bound_face(:,j), count)
        if (count == 2) then
            num(1) = num(1) + 1
            elem(num(1)) = j
            do k = 1, 4
                if (Hexa(ng)%bound_face(k,j) == this(1)) then
                    temp = (/ k+1, k+2, k+3 /)
                    if (temp(1) > 4) temp(1) = temp(1) - 4
                    if (temp(2) > 4) temp(2) = temp(2) - 4
                    if (temp(3) > 4) temp(3) = temp(3) - 4
                    p(:,1) = Hexa(ng)%node(:,Hexa(ng)%bound_face(k,j))
                    p(:,3) = Hexa(ng)%node(:,Hexa(ng)%bound_face(temp(2),j))
                    if (Hexa(ng)%bound_face(temp(1),j) == this(2)) then
                        p(:,2) = Hexa(ng)%node(:,Hexa(ng)%bound_face(temp(1),j))
                    else
                        p(:,2) = Hexa(ng)%node(:,Hexa(ng)%bound_face(temp(3),j))
                    endif
                    call calc_plane(p(:,1), p(:,2), p(:,3), pl(:,2))
                    call calc_angle_cos(pl(1:3,1), op, pl(1:3,2), alpha(num(1)))
                    exit
                endif
            enddo
        endif
        if (num(1) == 2) exit
    enddo
    !write (*,*) elem(1), 'element angle:', alpha(1)
    !write (*,*) elem(2), 'element angle:', alpha(2)
    if (alpha(1) <= 10e-4) then
        first_elem = elem(1)
    elseif (alpha(2) <= 10e-4) then
        first_elem = elem(2)
    else
        if (alpha(1) < alpha(2)) then
            first_elem = elem(1)
        else
            first_elem = elem(2)
        endif
    endif
        
    end_elem_num = 1
    end_elem(1) = first_elem
    !write (*,*) remesh%conn_node(Hexa(ng)%bound_face(:,first_elem))
    chk_line_num = 0
    do j = 1, 4
        if (j == 4) then
            line = (/ Hexa(ng)%bound_face(j,first_elem), Hexa(ng)%bound_face(1,first_elem) /)
        else
            line = (/ Hexa(ng)%bound_face(j,first_elem), Hexa(ng)%bound_face(j+1,first_elem) /)
        endif
        
        check = .true.
        do k = 1, end_line_num
            call find_count(2, line, 2, end_line(:,k), count)
            if (count == 2) then
                check = .false.
                exit
            endif
        enddo
        if (check) then
            chk_line_num = chk_line_num + 1
            chk_line(:,chk_line_num) = line
        endif
    enddo
    
    num(1) = 0
    find_remainder_elem: do 
        num(1) = num(1) + 1
        if (num(1) == bn) stop 'interation is too much'
        
        chk_elem_num = 0
        chk_elem = 0
        do j = 1, chk_line_num
            line = chk_line(:,j)
            do k = 1, Hexa(ng)%bfn
                call find_count(2, line, 4, Hexa(ng)%bound_face(:,k), count)
                if (count == 2) then
                    call array_include_check(k, end_elem_num, end_elem(1:end_elem_num), check)
                    if (check == .false.) then
                        !write (*,*) k, chk_elem_num, chk_elem(1:chk_elem_num)
                        call array_include_check(k, chk_elem_num, chk_elem(1:chk_elem_num), check)
                        if (check == .false.) then
                            chk_elem_num = chk_elem_num + 1
                            chk_elem(chk_elem_num) = k
                        endif
                        exit
                    endif
                endif
            enddo
            end_line_num = end_line_num + 1
            end_line(:,end_line_num) = line
        enddo
        
        chk_line_num = 0
        chk_line = 0
        if (chk_elem_num == 0) then
            exit find_remainder_elem
        else
            do j = 1, chk_elem_num
                elem = Hexa(ng)%bound_face(:,chk_elem(j))
                do k = 1, 4
                    if (k == 4) then
                        line = (/ elem(k), elem(1) /)
                    else
                        line = (/ elem(k), elem(k+1) /)
                    endif
                    check = .true.
                    do q = 1, end_line_num
                        call find_count(2, line, 2, end_line(:,q), count)
                        if (count == 2) then
                            check = .false.
                            exit
                        endif
                    enddo
                    if (check) then
                        check = .true.
                        do q = 1, chk_line_num
                            call find_count(2, line, 2, chk_line(:,q), count)
                            if (count == 2) then
                                check = .false.
                                exit
                            endif
                        enddo
                        if (check) then
                            chk_line_num = chk_line_num + 1
                            chk_line(:,chk_line_num) = line
                        endif
                    endif
                enddo
                end_elem_num = end_elem_num + 1
                end_elem(end_elem_num) = chk_elem(j)
            enddo
        endif
    enddo find_remainder_elem
    !write (*,*) 'end_elem_num:', end_elem_num
    !write (*,*) end_elem(1:end_elem_num)
    !write (*,*) remesh%conn_node(Hexa(ng)%bound_face(:,end_elem(1)))
    
    inner_node_num = 0
    inner_node = 0
    do j = 1, end_elem_num
        do k = 1, 4
            fp = Hexa(ng)%bound_face(k,end_elem(j))
            call array_include_check(fp, ccn, closed_curve, check)
            if (check == .false.) then
                call array_include_check(fp, inner_node_num, inner_node(1:inner_node_num), check)
                if (check == .false.) then
                    inner_node_num = inner_node_num + 1
                    inner_node(inner_node_num) = fp
                endif
            endif
        enddo
    enddo
    !inner_node(1:inner_node_num) = remesh%conn_node(inner_node(1:inner_node_num))
    proj_pn(conn(inner_node(1:inner_node_num))) = i
    bound_lv(inner_node(1:inner_node_num)) = 1
    !write (*,*) 'inner_node_num:', inner_node_num
    !write (*,*) inner_node(1:inner_node_num)
    !write (*,*) 
    deallocate(line_info, closed_curve)
enddo      

deallocate (end_elem, chk_elem)


end subroutine set_proj_pn

subroutine check_repetition(num, arr, a, check)

implicit none

integer :: num, arr(num), a
logical :: check

integer :: i

check = .false.
do i = 1, num
    if (arr(i) == a) then
        check = .true.
        exit
    endif
enddo

end subroutine check_repetition



subroutine find_next_point(p2cn, p2_conn, current_pl, nn, node, ne, elem, ce, cp, cp_ele_node, info)

implicit none

integer, intent(in) :: p2cn, p2_conn(p2cn), nn, ne, elem(3,ne), ce, cp
integer, intent(inout) :: cp_ele_node, info(5)
real, intent(in) :: current_pl(4), node(3,nn)

integer :: i, j, k, count(3), limit, num, remainder, temp_ce, temp(3)
real :: dis, min_dis
logical :: check

!write (*,*) 'ce, cp, cp_ele_node:', ce, cp, cp_ele_node
!write (*,*) 'p2_conn:', p2_conn

check = .true.
do j = 1, p2cn
    if (p2_conn(j) == cp_ele_node) then
        info(1:2) = (/ ce, cp_ele_node /)
        info(5) = j
        check = .false.
        exit 
    endif
enddo

temp_ce = ce
if (check) then
    limit = 12
    num = 0
    info = 0
    find_ele: do 
        if (num >= limit) stop 'error(find_next_point): cannot find next node'
        num = num + 1
        min_dis = 10e8
        temp = 0
        do i = 1, ne
            if (i /= temp_ce) then
                count = 0
                do j = 1, 3
                    if (elem(j,i) == cp .or. elem(j,i) == cp_ele_node) count(j) = 1
                enddo
                if (sum(count) == 2) then
                    do j = 1, 3
                        if (count(j) == 0) then
                            remainder = elem(j,i)
                            exit
                        endif
                    enddo
                    do j = 1, p2cn
                        if (p2_conn(j) == remainder) then
                            dis = 0.0
                            do k = 1, 3
                                dis = dis + current_pl(k)*node(k,remainder)
                            enddo
                            dis = dis + current_pl(4)
                            !write (*,'(A,I7,A,F12.7,A,F12.7)') 'element:',i, ', min_dis:', min_dis, ', dis:', dis
                            if (min_dis > dis) then
                                min_dis = dis
                                temp = (/ i, remainder, j /)
                            endif
                            exit
                        endif
                    enddo
                    if (temp(1) == 0) then
                        cp_ele_node = remainder
                        temp_ce = i
                        exit
                    endif
                endif
            endif
        enddo
        if (temp(1) /= 0) then
            info(1:2) = temp(1:2)
            info(5) = temp(3)
            exit find_ele
        endif
    enddo find_ele !find_ele
endif
!write (*,*) 'info:', info

                            
end subroutine find_next_point


subroutine find_face_ele(ne, elem, except_ele, temp_bfn, bound_fn, bound_face)

implicit none

integer, intent(in) :: ne, elem(3,ne), temp_bfn
integer, intent(inout) :: bound_fn, bound_face(temp_bfn), except_ele(ne)

integer :: i, j, k, q, count, line(2), check_ele(ne), pre_num, pre_fn, ele_count

pre_fn = 0
do
    count = 0
    pre_num = bound_fn
    do i = pre_fn+1, bound_fn
        do j = 1, 3
            line(1) = elem(j,bound_face(i))
            if (j == 3) then
                line(2) = elem(1,bound_face(i))
            else
                line(2) = elem(j+1,bound_face(i))
            endif
            do k = 1, ne
                if (except_ele(k) == 0) then
                    ele_count = 0
                    do q = 1, 3
                        if (elem(q,k) == line(1) .or. elem(q,k) == line(2)) ele_count = ele_count + 1
                    enddo
                    if (ele_count == 2) then
                        bound_fn = bound_fn + 1
                        bound_face(bound_fn) = k
                        count = count + 1
                        except_ele(k) = 1
                        exit
                    endif
                endif
            enddo
        enddo
    enddo
    !write (*,*) 'pre_num, bound_fn, count:', pre_num, bound_fn, count
    if (count == 0) then
        exit
    else
        pre_fn = pre_num
    endif
enddo

end subroutine find_face_ele