module remesh_domain

implicit none
save
    
Type mesh_type
    integer :: remesh_flag  ! = 0: not remesh, /= 0: remesh
    integer :: nn, ne, sn, sfn
    integer, allocatable :: elem(:,:)
    real, allocatable :: node(:,:)
    
    ! boundary representation info
    integer :: bfn, ben, bn, inn, np2
    integer, allocatable :: bp2(:), bound_fn(:), bedn(:), b_edge(:,:), bound_conn(:), bound_node(:)
    integer, allocatable :: bound_face(:,:), node_lv(:,:), bound_bpln(:), bound_pl(:,:)
    integer, allocatable :: inner_node(:), inn_nci_num(:), inn_nci(:,:)
    
    ! manually variable 
    integer :: begn, fixed_point_group_num
    integer, allocatable :: b_edge_group_num(:), b_edge_group(:,:), b_edge_divnum(:), cen(:), control_edge_number(:,:)
    integer, allocatable :: b_edge_conn(:), fixed_point_group(:,:)
    real, allocatable :: control_edge_coord(:,:,:), bp2_coord(:,:)
    
    real :: tol = 0.01
End type

Type add_data
    integer :: ori_np2 = 0,  np2 = 0, ln = 0, sub_hexa_num = 0, fix_line_num, fix_point_num
    integer, allocatable :: p2(:,:), line(:,:), sub_hexa(:,:), ori_p2(:), p2_lv(:), fix_line(:,:), p2_flag(:), p2_node(:,:), fix_p_axis(:)
    real, allocatable :: ori_p2_coord(:,:), p2_coord(:,:), fix_coord(:,:)
End Type

type (mesh_type) :: Tri
type (mesh_type), allocatable :: Hexa(:)
type (add_data), allocatable :: remesh_add(:)
!type (add_data) :: remesh_add
integer :: arr_number
real :: remesh_dis = 0.004

contains

subroutine create_remesh_domain(ng, flag)

implicit none

integer, intent(in) :: ng, flag

if (flag == 2) then
    allocate (Hexa(ng))
elseif (flag == 3) then
    allocate (remesh_add(ng))
endif


end subroutine create_remesh_domain


subroutine update_boundary_info(ng, np2, bp2, add_nn, add_coord, ln, tcen, cen, edge, bfn, temp_bfn, bound_fn, bound_face, bound_pln, bound_pl)

implicit none
integer, intent(in) :: ng, np2, bp2(np2), ln, tcen, cen(ln), edge(ln,ln), add_nn
integer, intent(in) :: bfn, temp_bfn, bound_fn(bfn), bound_face(temp_bfn,bfn), bound_pln(bfn), bound_pl(np2*4,bfn)
real, intent(in) :: add_coord(3,20)
integer :: i, num

Tri%np2 = np2
allocate (Tri%bp2(np2))
Tri%bp2 = bp2
Tri%ben = tcen
num = maxval(cen(1:tcen))
allocate (Tri%bedn(tcen), Tri%b_edge(num,tcen))
Tri%bedn = 0
Tri%b_edge = 0
do i = 1, tcen
    Tri%bedn(i) = cen(i)
    Tri%b_edge(1:cen(i),i) = edge(1:cen(i),i)
enddo

Tri%bfn = bfn
num = maxval(bound_fn(1:bfn))
allocate (Tri%bound_fn(bfn), Tri%bound_face(num,bfn))
num = maxval(bound_pln(1:bfn))
allocate (Tri%bound_bpln(bfn), Tri%bound_pl(num,bfn))
Tri%bound_pl = 0
Tri%bound_face = 0
do i = 1, bfn
    Tri%bound_fn(i) = bound_fn(i)
    Tri%bound_face(1:bound_fn(i),i) = bound_face(1:bound_fn(i),i)
    Tri%bound_bpln(i) = bound_pln(i)
    Tri%bound_pl(1:bound_pln(i),i) = bound_pl(1:bound_pln(i),i)
enddo

call set_edge_level(ng, add_nn, add_coord)

end subroutine update_boundary_info

subroutine set_edge_level(ng,  add_nn, add_coord)

implicit none

integer, intent(in) :: ng, add_nn
real, intent(in) :: add_coord(3,20)

integer :: i, j, k, q, lv, count, tln, tlgn, ce, temp_k, pre_tln
integer :: line(2), line2(2), this(4)
integer :: line_group_num(Tri%ben/3), line_group(Tri%ben/2,Tri%ben/3)
integer :: line_flag(Tri%ben), sub_hexa_flag(remesh_add(ng)%sub_hexa_num)
integer :: hexa_line(2,12), temp_line(Tri%ben)
real :: tol_dis, dis, b_edge_dis(Tri%ben), p(3,2)
logical :: check
integer :: fixed_sub_hexa_num, fixed_line(2)
integer, allocatable :: fixed_num(:), fixed_line_num(:)
!
hexa_line(1,:) = (/ 1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 3, 4 /)
hexa_line(2,:) = (/ 2, 4, 6, 8, 3, 1, 7, 5, 5, 6, 7, 8 /)
!fixed_sub_hexa_num = remesh_add(ng)%fix_line(1)
!fixed_line = remesh_add(ng)%fix_line(2:3)
!fixed_line = remesh_add(ng)%sub_hexa(fixed_line,fixed_sub_hexa_num)
!fixed_num = remesh_add(ng)%fix_line(4)
!fixed_line_num = 0
allocate (fixed_num(remesh_add(ng)%fix_line_num), fixed_line_num(remesh_add(ng)%fix_line_num))
fixed_num = 0
fixed_line_num = 0
do i = 1, remesh_add(ng)%fix_line_num
    fixed_sub_hexa_num = remesh_add(ng)%fix_line(1,i)
    fixed_line = remesh_add(ng)%fix_line(2:3,i)
    fixed_line = remesh_add(ng)%sub_hexa(fixed_line,fixed_sub_hexa_num)
    remesh_add(ng)%fix_line(2:3,i) = fixed_line
enddo
!
!allocate (Tri%b_edge_lv(Tri%ben))
!Tri%b_edge_lv = 0
!do i = 1, Tri%ben
!    line = (/ Tri%b_edge(1,i), Tri%b_edge(Tri%bedn(i),i) /)
!    lv = 0
!    do j = 1, Tri%sub_hexa_num
!        call find_count(2, line, 8, Tri%sub_hexa(:,j), count)
!        if (count == 2) lv = lv + 1
!    enddo
!    !write (*,*) i, 'edge level:', lv
!    Tri%b_edge_lv(i) = lv
!enddo
line_flag = 0
tlgn = 0
do i = 1, Tri%ben
    if (line_flag(i) == 0) then
        line_flag(i) = 1
        tln = 1
        temp_line(tln) = i
        sub_hexa_flag = 0
        do 
            pre_tln = tln
            do j = 1, remesh_add(ng)%sub_hexa_num
                if (sub_hexa_flag(j) == 0) then
                    temp_k = 0
                    do k = 1, 12
                        this(1:2) = remesh_add(ng)%sub_hexa(hexa_line(:,k),j)
                        do q = 1, tln
                            ce = temp_line(q)
                            this(3:4) = (/ Tri%b_edge(1,ce), Tri%b_edge(Tri%bedn(ce),ce) /)
                            call find_count(2, this(1:2), 2, this(3:4), count)
                            if (count == 2) then
                                temp_k = k
                                exit
                            endif
                        enddo
                        if (temp_k /= 0) then
                            exit
                        endif
                    enddo
                    
                    if (temp_k /= 0) then
                        if (temp_k <= 4) then
                            this = (/ 1, 2, 3, 4 /)
                        elseif (temp_k >= 5 .and. temp_k <= 8) then
                            this = (/ 5, 6, 7, 8 /)
                        else
                            this = (/ 9, 10, 11, 12 /)
                        endif
                        do k = 1, 4
                            line = remesh_add(ng)%sub_hexa(hexa_line(:,this(k)),j)
                            check = .TRUE.
                            do q = 1, tln
                                ce = temp_line(q)
                                line2 = (/ Tri%b_edge(1,ce), Tri%b_edge(Tri%bedn(ce),ce) /)
                                call find_count(2, line, 2, line2, count)
                                if (count == 2) then
                                    check = .FALSE.
                                    exit
                                endif
                            enddo
                            if (check) then
                                tln = tln + 1
                                do q = 1, Tri%ben
                                    if (line_flag(q) == 0) then
                                        line2 = (/ Tri%b_edge(1,q), Tri%b_edge(Tri%bedn(q),q) /)
                                        call find_count(2, line, 2, line2, count)
                                        if (count == 2) then
                                            line_flag(q) = 1
                                            temp_line(tln) = q
                                            !write (*,*) tln, line2
                                            exit
                                        endif
                                    endif
                                enddo
                            endif
                        enddo
                        sub_hexa_flag(j) = 1
                    endif
                endif
            enddo
            if (pre_tln == tln) exit
        enddo
        !write (*,*) 'tln:', tln
        !write (*,*) temp_line(1:tln)
        tlgn = tlgn + 1
        line_group_num(tlgn) = tln
        line_group(1:tln,tlgn) = temp_line(1:tln)
    endif
enddo
tln = maxval(line_group_num(1:tlgn))
allocate (Tri%b_edge_group_num(tlgn), Tri%b_edge_group(tln,tlgn))
Tri%begn = tlgn
do i = 1, tlgn
    Tri%b_edge_group_num(i) = line_group_num(i)
    Tri%b_edge_group(1:line_group_num(i),i) = line_group(1:line_group_num(i),i)
enddo

b_edge_dis = 0.0
do i = 1, Tri%ben
    tol_dis = 0.0
    do j = 1, Tri%bedn(i)-1
        line = Tri%b_edge(j:j+1,i)
        do k = 1, 2
            if (line(k) > Tri%nn) then
                p(:,k) = add_coord(:,line(k)-Tri%nn)
            else
                p(:,k) = Tri%node(:,line(k))
            endif
        enddo
        call calc_length_3D(p(:,1), p(:,2), dis)
        tol_dis = tol_dis + dis
    enddo
    do j = 1, remesh_add(ng)%fix_line_num
        if (fixed_line_num(j) == 0) then
            line = (/ Tri%b_edge(1,i), Tri%b_edge(Tri%bedn(i),i) /)
            fixed_line = remesh_add(ng)%fix_line(2:3,j)
            call find_count(2, line, 2, fixed_line, count)
            if (count == 2) then
                fixed_num(j) = remesh_add(ng)%fix_line(4,j)
                fixed_line_num(j) = i
            endif
        endif
    enddo
    b_edge_dis(i) = tol_dis
    !write (*,*) i, 'b_edge_dis:', b_edge_dis(i)
enddo

allocate (Tri%b_edge_divnum(tlgn))
Tri%b_edge_divnum = 0.0
do i = 1, tlgn
    check = .FALSE.
    tol_dis = 0.0
    do j = 1, Tri%b_edge_group_num(i)
        tol_dis = tol_dis + b_edge_dis(Tri%b_edge_group(j,i))
        do k = 1, remesh_add(ng)%fix_line_num
            if (Tri%b_edge_group(j,i) == fixed_line_num(k)) then
                check = .TRUE.
                temp_k = k
                exit
            endif
        enddo
        if (check) exit
    enddo
    if (check) then
        Tri%b_edge_divnum(i) = fixed_num(temp_k)
    else
        tol_dis = tol_dis/float(Tri%b_edge_group_num(i))
        Tri%b_edge_divnum(i) = nint(tol_dis/remesh_dis)
        if (Tri%b_edge_divnum(i) <= 2) Tri%b_edge_divnum(i) = 2 ! Minimum element number at each layer is 2
    endif
    write (*,*) 'Tri%b_edge_divnum:', Tri%b_edge_divnum(i)
enddo
   

end subroutine set_edge_level



subroutine add_data_for_remesh(ng, np2, p2, add_nn, add_coord)

implicit none

integer, intent(in) :: ng, np2, p2(2,np2)
integer, intent(inout) :: add_nn
real, intent(inout) :: add_coord(3,20)

integer :: i, j, temp_j, check_p2(np2), ori_num, line(2), sub_hexa(8), node_num(3)
real :: dis, min_dis, p(3,2)
logical :: check

!if (remesh_add(ng)%ori_np2 == 0) then
!    remesh_add(ng)%ori_np2 = 8
!    allocate (remesh_add(ng)%ori_p2(remesh_add(ng)%ori_np2))
!    !remesh_add%ori_p2 = (/ 97, 120, 121, 125, 2105, 2139, 2141, 2149 /)   ! ACM 
!    if (ng == 1) then
!        remesh_add(ng)%ori_p2 = (/ 1, 12, 49, 60, 714, 736, 780, 802 /)   ! BCM 1 subdomain
!    elseif (ng == 2) then
!        remesh_add(ng)%ori_p2 = (/ 1, 4, 301, 304, 690, 696, 1142, 1148 /)   ! BCM 2 subdomain
!    endif
!endif
ori_num = remesh_add(ng)%ori_np2
check_p2 = 0
do i = 1, remesh_add(ng)%ori_np2
    check = .FALSE.
    do j = 1, np2
        if (check_p2(j) == 0) then
            if (remesh_add(ng)%ori_p2(i) == p2(1,j)) then
                check = .TRUE.
                check_p2(j) = 1
                exit
            endif
        endif
    enddo
    if (check == .FALSE.) then
        p(:,1) = Tri%node(:,remesh_add(ng)%ori_p2(i))
        min_dis = 10e8
        do j = 1, np2
            if (check_p2(j) == 0) then
                p(:,2) = Tri%node(:,p2(1,j))
                call calc_length_3D(p(:,1), p(:,2), dis)
                if (dis < min_dis) then
                    min_dis = dis 
                    temp_j = j
                endif
            endif
        enddo
        remesh_add(ng)%ori_p2(i) = p2(1,temp_j)
        check_p2(temp_j) = 1
    endif
enddo

add_nn = 0
if (remesh_add(ng)%np2 /= 0) then
    do i = 1, remesh_add(ng)%np2
        if (remesh_add(ng)%p2_flag(i) == 2) then
            node_num = remesh_add(ng)%p2_node(:,i)
            do j = 1, 3
                if (node_num(j) <= remesh_add(ng)%ori_np2) then
                    p(j,1) = Tri%node(j,remesh_add(ng)%ori_p2(node_num(j)))
                else
                    node_num(j) = node_num(j) - remesh_add(ng)%ori_np2
                    p(j,1) = Tri%node(j,remesh_add(ng)%p2(1,node_num(j)))
                endif
                write (*,*) 'node_num:', node_num
                write (*,*) p(:,1)
            enddo
            add_nn = add_nn + 1
            add_coord(:,add_nn) = p(:,1)
            remesh_add(ng)%p2(1,i) = Tri%nn+add_nn
        endif
    enddo
endif

do i = 1, remesh_add(ng)%ln
    line = remesh_add(ng)%line(:,i)
    do j = 1, 2
        if (line(j) <= ori_num) then
            remesh_add(ng)%line(j,i) = remesh_add(ng)%ori_p2(line(j))
        else
            remesh_add(ng)%line(j,i) = remesh_add(ng)%p2(1,line(j)-ori_num)
        endif
    enddo
enddo
do i = 1, remesh_add(ng)%sub_hexa_num
    sub_hexa = remesh_add(ng)%sub_hexa(:,i)
    do j = 1, 8
        if (sub_hexa(j) <= ori_num) then
            remesh_add(ng)%sub_hexa(j,i) = remesh_add(ng)%ori_p2(sub_hexa(j))
        else
            remesh_add(ng)%sub_hexa(j,i) = remesh_add(ng)%p2(1,sub_hexa(j)-ori_num)
        endif
    enddo
enddo

write (*,*) 'remesh_add%ori_p2:', remesh_add(ng)%ori_p2(:) 
write (*,*)
if (remesh_add(ng)%ln /= 0) then
    write (*,*) 'remesh_add%line(1,:):', remesh_add(ng)%line(1,:)
    write (*,*) 'remesh_add%line(2,:):', remesh_add(ng)%line(2,:)
    write (*,*)
endif
do i = 1, remesh_add(ng)%sub_hexa_num
    write (*,'(I2,A,8(I6,1X))') i, ' remesh_add%sub_hexa:', remesh_add(ng)%sub_hexa(:,i)
enddo

!if (remesh_add%np2 == 0) then
!    remesh_add%np2 = 8
!    remesh_add%ln = 14
!    remesh_add%sub_hexa_num = 4
!    allocate (remesh_add%p2(2,remesh_add%np2), remesh_add%p2_lv(remesh_add%np2))
!    allocate (remesh_add%line(2,remesh_add%ln), remesh_add%sub_hexa(8,remesh_add%sub_hexa_num))
!    
!    remesh_add%p2(1,:) = (/ 2252, 2605, 2191, 2432, 210, 505, 151, 359 /)
!    remesh_add%p2(2,:) = 3
!    remesh_add%p2_lv = (/ 
!    
!    remesh_add%line(1,:) = (/ 2252, 2252, 2252, 2191, 2191, 210, 210, 210, 151, 151, 2252, 2605, 2191, 2432 /)
!    remesh_add%line(2,:) = (/ 2605, 2191, remesh_add%ori_p2(6), remesh_add%ori_p2(7), 2432, 505, 151, &
!                              remesh_add%ori_p2(2), remesh_add%ori_p2(3), 359, 210, 505, 151, 359 /)
!    
!    remesh_add%sub_hexa_num = 4
!    remesh_add%sub_hexa(:,1) = (/ remesh_add%ori_p2(5), remesh_add%ori_p2(6), 2252, 2605, remesh_add%ori_p2(1), remesh_add%ori_p2(2), 210, 505 /)
!    remesh_add%sub_hexa(:,2) = (/ 2252, 2191, 2432, 2605, 210, 151, 359, 505 /)
!    remesh_add%sub_hexa(:,3) = (/ 2252, remesh_add%ori_p2(6), remesh_add%ori_p2(7), 2191, 210, remesh_add%ori_p2(2), remesh_add%ori_p2(3), 151 /)
!    remesh_add%sub_hexa(:,4) = (/ 2432, 2191, remesh_add%ori_p2(7), remesh_add%ori_p2(8), 359, 151, remesh_add%ori_p2(3), remesh_add%ori_p2(4) /)
!    !remesh_add%sub_hexa(:,1) = (/ remesh_add%ori_p2(5), 598, 502, remesh_add%ori_p2(6), remesh_add%ori_p2(1), 123, 41, remesh_add%ori_p2(2) /)
!    !remesh_add%sub_hexa(:,2) = (/ 598, 743, 510, 502, 123, 226, 49, 41 /)
!    !remesh_add%sub_hexa(:,3) = (/ 502, 510, remesh_add%ori_p2(7), remesh_add%ori_p2(6), 41, 49, remesh_add%ori_p2(3), remesh_add%ori_p2(2) /)
!    !remesh_add%sub_hexa(:,4) = (/ 743, remesh_add%ori_p2(8), remesh_add%ori_p2(7), 510, 226, remesh_add%ori_p2(4), remesh_add%ori_p2(3), 49 /)
!endif

end subroutine add_data_for_remesh


!subroutine revise_add_data_for_remesh(new_nn, new_coord, new_bound_lv)
!
!implicit none
!
!integer, intent(in) :: new_nn, new_bound_lv(new_nn)
!real, intent(in) :: new_coord(3,new_nn)
!
!integer :: i, j, k, temp_num, pn, cbn, piar_nodes(2,remesh_add%ln*2)
!real :: min_dis, dis, p(3,2)
!logical :: check
!
!! match the Hexa%bound_node & Add data points
!do i = 1, remesh_add%ori_np2
!    p(:,1) = Tri%node(:,remesh_add%ori_p2(i))
!    min_dis = 10e8
!    temp_num = 0
!    do j = 1, new_nn
!        p(:,2) = new_coord(:,j)
!        call calc_length_3D(p(:,1), p(:,2), dis)
!        if (dis < remesh_dis*0.001) then
!            temp_num = j
!            exit
!        elseif (min_dis > dis) then
!            temp_num = j
!            min_dis = dis
!        endif
!    enddo
!    remesh_add%ori_p2(i) = temp_num
!enddo
!
!pn = 0
!piar_nodes = 0
!do i = 1, remesh_add%np2
!    p(:,1) = Tri%node(:,remesh_add%p2(1,i))
!    min_dis = 10e8
!    temp_num = 0
!    do j = 1, new_nn
!        if (new_bound_lv(j) >= 2) then
!            p(:,2) = new_coord(:,j)
!            call calc_length_3D(p(:,1), p(:,2), dis)
!            if (dis < remesh_dis*0.001) then
!                temp_num = j
!                exit
!            elseif (min_dis > dis) then
!                temp_num = j
!                min_dis = dis
!            endif
!        endif
!    enddo
!    pn = pn + 1
!    piar_nodes(:,pn) = (/ remesh_add%p2(1,i), temp_num /)
!    remesh_add%p2(1,i) = temp_num
!enddo
!
!do i = 1, remesh_add%ln
!    do j = 1, 2
!        cbn = remesh_add%line(j,i)
!        check = .FALSE.
!        do k = 1, pn
!            if (piar_nodes(1,k) == cbn) then
!                temp_num = piar_nodes(2,k)
!                check = .TRUE.
!                exit
!            endif
!        enddo
!        if (check) then
!            remesh_add%line(j,i) = temp_num
!        else
!            p(:,1) = Tri%node(:,cbn)
!            min_dis = 10e8
!            temp_num = 0
!            do k = 1, new_nn
!                p(:,2) = new_coord(:,k)
!                call calc_length_3D(p(:,1), p(:,2), dis)
!                if (dis < remesh_dis*0.001) then
!                    temp_num = k
!                    exit
!                elseif (min_dis > dis) then
!                    temp_num = k
!                    min_dis = dis
!                endif
!            enddo
!            pn = pn + 1
!            piar_nodes(:,pn) = (/ cbn, temp_num /)
!            remesh_add%line(j,i) = temp_num
!        endif
!    enddo
!enddo
!
!do i = 1, remesh_add%sub_hexa_num
!    do j = 1, 8
!        cbn = remesh_add%sub_hexa(j,i)
!        check = .FALSE.
!        do k = 1, pn
!            if (piar_nodes(1,k) == cbn) then
!                temp_num = piar_nodes(2,k)
!                check = .TRUE.
!                exit
!            endif
!        enddo
!        if (check) then
!            remesh_add%sub_hexa(j,i) = temp_num
!        else
!            p(:,1) = Tri%node(:,cbn)
!            min_dis = 10e8
!            temp_num = 0
!            do k = 1, new_nn
!                p(:,2) = new_coord(:,k)
!                call calc_length_3D(p(:,1), p(:,2), dis)
!                if (dis < remesh_dis*0.001) then
!                    temp_num = k
!                    exit
!                elseif (min_dis > dis) then
!                    temp_num = k
!                    min_dis = dis
!                endif
!            enddo
!            pn = pn + 1
!            piar_nodes(:,pn) = (/ cbn, temp_num /)
!            remesh_add%sub_hexa(j,i) = temp_num
!        endif
!    enddo
!enddo
!
!write (*,*) 'remesh_add%p2:', remesh_add(ng)%ori_p2(:) 
!write (*,*) 
!write (*,*) 'remesh_add%p2:', remesh_add(ng)%p2(1,:) 
!write (*,*)
!write (*,*) 'remesh_add%line(1,:):', remesh_add(ng)%line(1,:)
!write (*,*) 'remesh_add%line(2,:):', remesh_add(ng)%line(2,:)
!write (*,*)
!write (*,*) 'remesh_add%sub_hexa(1):', remesh_add%sub_hexa(:,1)
!write (*,*) 'remesh_add%sub_hexa(2):', remesh_add%sub_hexa(:,2)
!write (*,*) 'remesh_add%sub_hexa(3):', remesh_add%sub_hexa(:,3)
!write (*,*) 'remesh_add%sub_hexa(4):', remesh_add%sub_hexa(:,4)
!
!end subroutine revise_add_data_for_remesh

subroutine revise_b_edge_info(num_group, bnn, b_coord, b_corner, brgn, b_ridge_num, b_ridge)

implicit none

integer, intent(in) :: num_group, bnn, brgn
integer, intent(inout) :: b_corner(bnn), b_ridge_num(brgn), b_ridge(brgn,1000)
real, intent(in) :: b_coord(3,bnn)

integer :: i, j, k, temp_num, max_num, ng, cl, temp_ridge(1000), num
real :: dis, min_dis, cp(3), tol
integer, allocatable :: line_nodes(:), p2(:)
character(len=80) :: fn

write (*,*) 'num_group:', num_group
tol = remesh_dis*0.001
!fn = './output/solid/remesh/check_ridge_00.plt'
fn = './output/solid/remesh/check_ridge.plt'
!write (fn(35:36),'(I2.2)') ng
open (Unit=30, File=fn, STATUS='replace', ACTION='write')
num = 0
do ng = 1, num_group
    do i = 1, Hexa(ng)%ben
        num = num + Hexa(ng)%bedn(i) - 1
    enddo
enddo
Write (30,'(A,A,A)') 'TITLE="Check b_ridge"'
Write (30,*) 'VARIABLES="x", "y", "z"'
Write (30,'(A,I6,A,I6,A)') 'ZONE T = "ridge", N = ', bnn, ',E = ', num, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
do i = 1, bnn
    write (30,'(3ES15.7)') b_coord(:,i)
enddo

b_corner = 0
do ng = 1, num_group
    write (*,*) ng, ' Hexa(ng)%ben:', Hexa(ng)%ben
    !if (Hexa(ng)%remesh_flag /= 0) then
        max_num = maxval(Hexa(ng)%bedn(1:Hexa(ng)%ben))
        write (*,*) 'max_num, np2:', max_num, Hexa(ng)%np2
        allocate (p2(Hexa(ng)%np2), line_nodes(max_num))
        do i = 1, Hexa(ng)%np2
            min_dis = 10e8
            cp = Hexa(ng)%bp2_coord(:,i)
            temp_num = 0
            do j = 1, bnn
                call calc_length_3D(cp, b_coord(:,j), dis)
                if (dis < tol) then
                    temp_num = j
                    exit
                elseif (dis < min_dis) then
                    min_dis = dis
                    temp_num = j
                endif
            enddo
            p2(i) = temp_num
            b_corner(temp_num) = 3
        enddo
        write (*,*) 
        write (*,*) 'p2:', p2
        
        do i = 1, Hexa(ng)%ben
            cl = Hexa(ng)%b_edge_conn(i)
            b_ridge_num(cl) = Hexa(ng)%bedn(i)
            do j = 1, Hexa(ng)%bedn(i)
                min_dis = 10e8
                cp = Hexa(ng)%control_edge_coord(:,j,i)
                temp_num = 0
                do k = 1, bnn
                    call calc_length_3D(cp, b_coord(:,k), dis)
                    if (dis < tol) then
                        temp_num = k
                        exit
                    elseif (dis < min_dis) then
                        min_dis = dis
                        temp_num = k
                    endif
                enddo
                b_ridge(cl,j) = temp_num
            enddo
            !write (*,*) 
            write (*,*) cl, ' rideg number:', b_ridge_num(cl)
            temp_ridge = b_ridge(cl,:)
            write (*,*) temp_ridge(1:b_ridge_num(cl))
            do j = 1, Hexa(ng)%bedn(i)-1
                write (30,*) temp_ridge(j:j+1)
            enddo
        enddo
        deallocate (p2, line_nodes)
    !endif
enddo
close(30)

end subroutine revise_b_edge_info


subroutine rearrange_hexa_b_edge(ng, ln, tcen, cen, control_edge, temp_tcen, temp_cen, temp_control_edge, control_edge2)

implicit none
integer, intent(in) :: ng, ln, tcen, cen(ln), control_edge(ln,200), control_edge2(ln,200)
integer, intent(in) :: temp_tcen, temp_cen(ln), temp_control_edge(ln,200)

integer :: i, j, k, num, en(ln), check_edge(ln), pre_num

do i = 1, Hexa(ng)%ben
    !write (*,*) i, ' Pre Hexa(ng)%cen:', Hexa(ng)%cen(i)
    !write (*,*) Hexa(ng)%control_edge_number(1:Hexa(ng)%cen(i),i)
    if (Hexa(ng)%cen(i) == 0) then
        write (*,*) i, 'Hexa(ng)%cen:', Hexa(ng)%cen(i)
        write (*,*) temp_control_edge(1:temp_cen(i),i)
        stop 'Hexa control edge has not Tri control edge(rearrange_hexa_b_edge-remesh_domain)'
    elseif (Hexa(ng)%cen(i) /= 1) then
        en(1:Hexa(ng)%cen(i)) = Hexa(ng)%control_edge_number(1:Hexa(ng)%cen(i),i)
        Hexa(ng)%control_edge_number(:,i) = 0
        check_edge = 0
        num = 0 
        do j = 1, Hexa(ng)%cen(i)
            if (control_edge2(j,i) == 1) then
                num = num + 1
                Hexa(ng)%control_edge_number(num,i) = en(j)  ! frist control_edge
                check_edge(j) = 1
                pre_num = control_edge(cen(en(j)),en(j))
                exit
            endif
        enddo
        find_remainder_edge: do
            do j = 1, Hexa(ng)%cen(i)
                if (check_edge(j) == 0) then
                    if (pre_num == control_edge(1,en(j))) then
                        num = num + 1
                        Hexa(ng)%control_edge_number(num,i) = en(j)
                        if (control_edge2(j,i) == 2) then
                            exit find_remainder_edge
                        else
                            pre_num = control_edge(cen(en(j)),en(j))
                            exit
                        endif
                    endif
                endif
            enddo
            if (num == Hexa(ng)%cen(i)) exit find_remainder_edge
        enddo find_remainder_edge
    endif
    !write (*,*) i, ' Post Hexa(ng)%cen:', Hexa(ng)%cen(i)
    !write (*,*) Hexa(ng)%control_edge_number(1:Hexa(ng)%cen(i),i)
    !write (*,*)
enddo

end subroutine rearrange_hexa_b_edge


subroutine projection_point_to_plane(cur_ng, bn, proj_pn, nn, node)

implicit none
integer, intent(in) :: cur_ng, bn, nn, proj_pn(bn)
real, intent(inout) :: node(3,nn)

integer :: i, j, k, bpn, temp_k
real :: min_area, area, min_dis, dis, t
real :: op(3), temp_val(3), temp_base(3,3), pl(4), p(3,2), cp(3), val(3), elem_dis(4)
logical :: pl_check, check, inner_check

do i = 1, Tri%bfn
    do j = 1, bn
        bpn = Hexa(cur_ng)%bound_node(j)
        if (proj_pn(j) == i) then
            if (bpn == 23581) write (*,*) node(:,bpn)
            op = node(:,bpn)
            min_area = 10e8
            min_dis = 10e8
            temp_k = 0
            pl_check = .false.
            do k = 1, Tri%bound_fn(i)
                temp_base(:,1) = Tri%node(:,Tri%elem(1,Tri%bound_face(k,i)))
                temp_base(:,2) = Tri%node(:,Tri%elem(2,Tri%bound_face(k,i)))
                temp_base(:,3) = Tri%node(:,Tri%elem(3,Tri%bound_face(k,i)))
                call calc_plane(temp_base(:,1), temp_base(:,2), temp_base(:,3), pl)
                dis = abs(pl(1)*op(1) + pl(2)*op(2) + pl(3)*op(3) + pl(4))
                p(:,1) = op(:) + pl(1:3)*6.0*sqrt(3.0)
                p(:,2) = op(:) - pl(1:3)*6.0*sqrt(3.0)
                call cross_point_plane(pl, p(:,1), p(:,2), cp, t, check)
                if (check == .false.) then
                    call check_inner_area(temp_base, cp, inner_check, area)
                    if (inner_check) then
                        if (min_dis > dis) then
                            min_dis = dis
                            temp_val = cp
                            pl_check = .true.
                        endif
                    else
                        val = cp
                        call move_point_near_ele(Tri%nn, Tri%node, val, cp, Tri%elem(:,Tri%bound_face(k,i)))
                        call calc_length_3D(val, cp, dis)
                        call calc_length_3D(temp_base(:,1), temp_base(:,2), elem_dis(1))
                        call calc_length_3D(temp_base(:,2), temp_base(:,3), elem_dis(2))
                        call calc_length_3D(temp_base(:,3), temp_base(:,1), elem_dis(3))
                        elem_dis(4) = maxval(elem_dis(1:3))
                        if (dis <= elem_dis(4)*sqrt(3.0)) then
                            if (min_area > area .and. pl_check == .false.) then
                                min_area = area
                                temp_val = cp
                                temp_k = k 
                            endif
                        endif
                    endif
                endif
            enddo
            if (pl_check) then
                node(:,bpn) = temp_val
            else
                if (temp_k == 0) then
                    write (*,*) 'Cannot find the near boundary plane at ', bpn, 'node'
                    write (*,*) 'area:', min_area, ',dis:', min_dis
                    stop
                else
                    call move_point_near_ele(Tri%nn, Tri%node, temp_val, cp, Tri%elem(:,Tri%bound_face(temp_k,i)))
                    node(:,bpn) = cp
                endif
            endif
            if (bpn == 23581) write (*,*) node(:,bpn)
        endif
    enddo
enddo

end subroutine projection_point_to_plane



subroutine free_geo_domain(flag)

implicit none
integer, intent(in) :: flag 

! mesh
if (flag == 1) then
    Tri%nn = 0 ;  Tri%ne = 0;  Tri%sn = 0;  Tri%sfn = 0
    deallocate (Tri%node, Tri%elem)
    Tri%np2 = 0;  Tri%ben = 0;  Tri%bfn = 0
    deallocate (Tri%bp2)
    deallocate (Tri%bedn, Tri%b_edge)
    deallocate (Tri%bound_fn, Tri%bound_face)
    deallocate (Tri%bound_bpln, Tri%bound_pl)
    deallocate (Tri%b_edge_group_num, Tri%b_edge_group, Tri%b_edge_divnum)
elseif (flag == 2) then
    if (allocated(Hexa)) deallocate (Hexa)
elseif (flag == 3) then
    if (allocated(remesh_add)) deallocate (remesh_add)
endif

end subroutine free_geo_domain


end module
