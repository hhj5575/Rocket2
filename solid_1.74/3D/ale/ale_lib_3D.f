subroutine calc_overlap_volume(ce, cq, target_ce, target_cubic, ball_cubic, overlap_volume)
    
implicit none

integer, intent(in) :: ce, cq, target_ce
real, intent(in) :: target_cubic(3,8), ball_cubic(3,8)
real, intent(inout) :: overlap_volume

integer :: i, j, k, q, count, overlap_num, pre_num, next_num, nci_num, tri_num, num, iter
integer :: line(2,12), face_nn(4,6), node_check(2,8), tri(3,2)
integer :: face_node_num(12), face_nodes(8,12), temp_num, temp_nodes(48), this(2), temp_face(8)
real :: diff_volume, t, tol, dis, volume, volume_tol, dis_tol
real :: overlap_coord(3,24), cp(3), cross_p(3), line_coord(3,2), face_coord(3,4), pl(4), tet_coord(3,4), tri_coord(3,3)
logical :: check(4), face_check
character(len=80) :: fn

integer, allocatable :: nci(:,:), tri_line(:,:), tri_check(:)

line(1,:) = (/ 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8 /)
line(2,:) = (/ 2, 3, 4, 1, 5, 6, 7, 8, 6, 7, 8, 5 /)
face_nn(1,:) = (/ 1, 5, 1, 2, 3, 1 /)
face_nn(2,:) = (/ 2, 8, 5, 6, 7, 4 /)
face_nn(3,:) = (/ 3, 7, 6, 7, 8, 8 /)
face_nn(4,:) = (/ 4, 6, 2, 3, 4, 5 /)
tol = 0.0
do i = 1, 12
    line_coord = ball_cubic(:,line(:,i))
    call calc_length_3D(line_coord(:,1), line_coord(:,2), dis)
    tol = tol + dis
    line_coord = target_cubic(:,line(:,i))
    call calc_length_3D(line_coord(:,1), line_coord(:,2), dis)
    tol = tol + dis
enddo
tol = tol/24.0
volume_tol = tol*tol*tol
dis_tol = 0.001

face_node_num = 0
face_nodes = 0

overlap_num = 0
overlap_coord = 0.0
node_check = 0
! check each points lie other cubic
do i = 1, 8
    cp = ball_cubic(:,i)
    call inout_check_cubic2(ce, cp, target_cubic, check(1), volume_tol)
    if (check(1)) then
        overlap_num = overlap_num + 1
        overlap_coord(:,overlap_num) = cp
        node_check(1,i) = overlap_num
        do j = 1, 6
            do k = 1, 4
                if (face_nn(k,j) == i) then
                    face_node_num(j) = face_node_num(j) + 1
                    face_nodes(face_node_num(j),j) = overlap_num
                    exit
                endif
            enddo
        enddo
    endif
enddo
do i = 1, 8
    cp = target_cubic(:,i)
    call inout_check_cubic2(ce, cp, ball_cubic, check(1), volume_tol)
    if (check(1)) then
        check(2) = .TRUE.
        do j = 1, overlap_num
            call calc_length_3D(overlap_coord(:,j), cp, dis)
            if (dis < tol*dis_tol) then
                check(2) = .FALSE.
                num = j
                exit
            endif
        enddo
        if (check(2)) then
            overlap_num = overlap_num + 1
            overlap_coord(:,overlap_num) = cp
            num = overlap_num
        endif
        node_check(2,i) = num
        do j = 1, 6
            do k = 1, 4
                if (face_nn(k,j) == i) then
                    face_node_num(j+6) = face_node_num(j+6) + 1
                    face_nodes(face_node_num(j+6),j+6) = num
                    exit
                endif
            enddo
        enddo
    endif
enddo
!if (ce == 155 .and. cq == 1) then
!    if (target_ce == 4 .or. target_ce == 9) then
        !write (*,*) 'ball_coord'
        !do i = 1, 8
        !    write (*,*) ball_cubic(:,i)
        !enddo
        !write (*,*) 
        !write (*,*) 'target_coord'
        !do i = 1, 8
        !    write (*,*) target_cubic(:,i)
        !enddo
        !write (*,*) 'point overlap_num:', overlap_num
!    endif
!endif

! check that line cross the other cubic
tri(:,1) = (/ 1, 2, 3 /)
tri(:,2) = (/ 3, 4, 1 /)
do i = 1, 12
    line_coord = ball_cubic(:,line(:,i))
    do j = 1, 6
        face_coord = target_cubic(:,face_nn(:,j))
        check = .FALSE.
        do q = 1, 2
            call calc_plane(face_coord(:,tri(1,q)), face_coord(:,tri(2,q)), face_coord(:,tri(3,q)), pl)
            call cross_point_plane(pl, line_coord(:,1), line_coord(:,2), cross_p, t, check(1))
            if (check(1) == .FALSE.) then
                tri_coord = face_coord(:,tri(:,q))
                call check_inner_face2(tri_coord, cross_p, check(2))
                if (check(2)) then
                    check(3) = .TRUE.
                    do k = 1, overlap_num
                        call calc_length_3D(cross_p, overlap_coord(:,k), dis)
                        if (dis < tol*dis_tol) then
                            check(3) = .FALSE.
                            num = k
                            exit
                        endif
                    enddo
                    if (check(3)) then
                        overlap_num = overlap_num + 1
                        overlap_coord(:,overlap_num) = cross_p
                        num = overlap_num
                    endif
                
                    do k = 1, 6
                        call find_count(2, line(:,i), 4, face_nn(:,k), count)
                        if (count == 2) then
                            call check_face_node(num, face_node_num(k), face_nodes(1:face_node_num(k),k), check(4))
                            if (check(4)) then
                                face_node_num(k) = face_node_num(k) + 1
                                face_nodes(face_node_num(k),k) = num
                            endif
                        endif
                    enddo
                    call check_face_node(num, face_node_num(j+6), face_nodes(1:face_node_num(j+6),j+6), check(4))
                    if (check(4)) then
                        face_node_num(j+6) = face_node_num(j+6) + 1
                        face_nodes(face_node_num(j+6),j+6) = num
                    endif
                endif
            endif
            if (check(3)) exit
        enddo
    enddo
enddo
do i = 1, 12
    line_coord = target_cubic(:,line(:,i))
    do j = 1, 6
        face_coord = ball_cubic(:,face_nn(:,j))
        check = .FALSE.
        do q = 1, 2
            call calc_plane(face_coord(:,tri(1,q)), face_coord(:,tri(2,q)), face_coord(:,tri(3,q)), pl)
            call cross_point_plane(pl, line_coord(:,1), line_coord(:,2), cross_p, t, check(1))
            if (check(1) == .FALSE.) then
                tri_coord = face_coord(:,tri(:,q))
                call check_inner_face2(tri_coord, cross_p, check(2))
                if (check(2)) then
                    check(3) = .TRUE.
                    do k = 1, overlap_num
                        call calc_length_3D(cross_p, overlap_coord(:,k), dis)
                        if (dis < tol*dis_tol) then
                            check(3) = .FALSE.
                            num = k
                            exit
                        endif
                    enddo
                    if (check(3)) then
                        overlap_num = overlap_num + 1
                        overlap_coord(:,overlap_num) = cross_p
                        num = overlap_num
                    endif
                
                    do k = 1, 6
                        call find_count(2, line(:,i), 4, face_nn(:,k), count)
                        if (count == 2) then
                            call check_face_node(num, face_node_num(k+6), face_nodes(1:face_node_num(k+6),k+6), check(4))
                            if (check(4)) then
                                face_node_num(k+6) = face_node_num(k+6) + 1
                                face_nodes(face_node_num(k+6),k+6) = num
                            endif
                        endif
                    enddo
                    call check_face_node(num, face_node_num(j), face_nodes(1:face_node_num(j),j), check(4))
                    if (check(4)) then
                        face_node_num(j) = face_node_num(j) + 1
                        face_nodes(face_node_num(j),j) = num
                    endif
                endif
            endif
            if (check(3)) exit
        enddo
    enddo
enddo

!if (ce == 155 .and. cq == 1) then
!    if (target_ce == 4 .or. target_ce == 9) then
        !write (*,*) 'overlap_num:', overlap_num
        !do i = 1, overlap_num
        !    write (*,'(I6,A,3F12.7)') i, ': ', overlap_coord(:,i)
        !enddo
        !write (*,*) 
        !do i = 1, 12
        !    write (*,*) i, face_nodes(1:face_node_num(i),i)
        !enddo
        !write (*,*) 
!    endif
!endif

overlap_volume = 0.0
if (overlap_num == 4) then
    tet_coord = overlap_coord(:,1:4)
    call calc_tet_volume(tet_coord, overlap_volume)
elseif (overlap_num > 4) then
    allocate (nci(8,overlap_num))
    nci = 0
    tri_num = 0
    do i = 1, overlap_num
        temp_num = 0
        do j = 1, 12
            do k = 1, face_node_num(j)
                if (i == face_nodes(k,j)) then
                    temp_nodes(temp_num+1:temp_num+face_node_num(j)) = face_nodes(1:face_node_num(j),j)
                    temp_num = temp_num + face_node_num(j)
                    exit
                endif
            enddo
        enddo
        call sortc(temp_nodes(1:temp_num), temp_num)
        !write (*,'(100(I4,1X))') i, temp_nodes(1:temp_num)
    
        nci_num = 0
        pre_num = 0
        count = 0
        do j = 1, temp_num
            if (pre_num == temp_nodes(j)) then
                count = count + 1
                if (j == temp_num .and. count >= 2 .and. pre_num /= i) then
                    nci_num = nci_num + 1
                    nci(nci_num,i) = pre_num
                    tri_num = tri_num + 1
                endif
            else
                if (count >= 2 .and. pre_num /= i) then
                    nci_num = nci_num + 1
                    nci(nci_num,i) = pre_num
                    tri_num = tri_num + 1
                endif
                pre_num = temp_nodes(j)
                count = 1
            endif
        enddo        
        !write (*,*) i, nci(1:nci_num,i)
    enddo

    allocate (tri_line(2,tri_num))
    tri_num = 0
    tri_line = 0
    do i = 1, overlap_num
        this(1) = i
        do j = 1, 8
            if (nci(j,i) == 0) then
                exit
            else
                this(2) = nci(j,i)
                check(1) = .TRUE.
                do k = 1, tri_num
                    call find_count(2, this, 2, tri_line(:,k), count)
                    if (count == 2) then
                        check(1) = .FALSE.
                        exit
                    endif
                enddo
                if (check(1)) then
                    tri_num = tri_num + 1
                    tri_line(:,tri_num) = this
                endif
            endif
        enddo
    enddo

    face_check = .TRUE.
    allocate (tri_check(tri_num))
    do i = 1, 12
        if (face_node_num(i) > 2) then
            num = face_node_num(i)
            temp_face(1:num) = face_nodes(1:num,i)
            temp_num = 0
            pre_num = 0
            next_num = 0
            iter = 0
            tri_check = 0
            do
                iter = iter + 1
                do j = 1, tri_num
                    if (j /= pre_num .and. tri_check(j) == 0) then
                        call find_count(2, tri_line(:,j), num, temp_face(1:num), count)
                        if (count == 2) then
                            if (pre_num == 0) then
                                face_nodes(temp_num+1:temp_num+2,i) = tri_line(:,j)
                                temp_num = temp_num + 2
                                pre_num = j
                                next_num = tri_line(2,j)
                                exit
                            else
                                if (tri_line(1,j) == next_num .or. tri_line(2,j) == next_num) then
                                    if (next_num == tri_line(1,j)) then
                                        next_num = tri_line(2,j)
                                    else
                                        next_num = tri_line(1,j)
                                    endif
                                    face_nodes(temp_num+1,i) = next_num
                                    temp_num = temp_num + 1
                                    pre_num = j
                                    exit
                                endif
                            endif
                        endif
                    endif
                enddo
                if (num == temp_num) exit
                if (iter == 20) then
                    !write (*,*)'Error(calc_overlap_volume): iteration is over than 20 at ', i, ce
                    !write (*,*) 'face_nodes:',  face_nodes(1:temp_num,i)
                    !write (*,*) 
                    !fn = './output/solid/remesh/overlap/check_overlap_line_000000_00.plt'
                    !write (fn(50:55),'(I6.6)') ce
                    !write (fn(57:58),'(I2.2)') cq
                    !open (Unit=30, File=fn, STATUS='replace', ACTION='write')
                    !Write (30,'(A,A,A)') 'TITLE="Check line(3D)"'
                    !Write (30,*) 'VARIABLES="x", "y", "z"'
                    !Write (30,'(A,I5,A,I5, A)') 'ZONE N =', overlap_num, ', E =', tri_num, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
                    !do j = 1, overlap_num
                    !    write (30,'(3ES15.7)') overlap_coord(:,j)
                    !enddo
                    !do j = 1, tri_num
                    !    write (30,*) tri_line(:,j)
                    !enddo
                    !close(30)
                    face_node_num(i) = 0
                    !face_check = .FALSE.
                    exit
                endif
            enddo
        endif
        !if (face_check == .FALSE.) exit
    enddo
    deallocate (tri_check)
    
    overlap_volume = 0.0
    if (overlap_num /= 0) then
        tet_coord(:,4) = 0.0
        do i = 1, overlap_num
            tet_coord(:,4) = tet_coord(:,4) + overlap_coord(:,i)
        enddo
        tet_coord(:,4) = tet_coord(:,4)/float(overlap_num)

        do i = 1, 12
            if (face_node_num(i) > 2) then
                tet_coord(:,1) = overlap_coord(:,face_nodes(1,i))
                do j = 1, face_node_num(i)-2
                    tet_coord(:,2) = overlap_coord(:,face_nodes(j+1,i))
                    tet_coord(:,3) = overlap_coord(:,face_nodes(j+2,i))
                    call calc_tet_volume(tet_coord, volume)
                    overlap_volume = overlap_volume + volume
                enddo
            endif    
        enddo
    endif
    deallocate(nci, tri_line)
endif

end subroutine calc_overlap_volume

subroutine check_face_node(cp, n, face_nodes, check)

implicit none

integer, intent(in) :: cp, n, face_nodes(n)
logical, intent(out) :: check

integer :: i

check = .TRUE.
do i = 1, n
    if (cp == face_nodes(i)) then
        check = .FALSE.
        exit
    endif
enddo

end subroutine check_face_node


subroutine overlap_check(ball_coord, target_coord, check)

implicit none

real, intent(in) :: ball_coord(3,8), target_coord(3,8)
logical, intent(out) :: check

integer :: i
real :: cp(3,2), dis(2), max_dis(2)

cp = 0.0
do i = 1, 8
    cp(:,1) = cp(:,1) + ball_coord(:,i)
    cp(:,2) = cp(:,2) + target_coord(:,i)
enddo
cp = cp * 0.125

max_dis = 0.0
do i = 1, 8
    call calc_length_3D(cp(:,1), ball_coord(:,i), dis(1))
    call calc_length_3D(cp(:,2), target_coord(:,i), dis(2))
    if (max_dis(1) < dis(1)) max_dis(1) = dis(1)
    if (max_dis(2) < dis(2)) max_dis(2) = dis(2)
enddo
call calc_length_3D(cp(:,1), cp(:,2), dis(1))
dis(2) = max_dis(1)+max_dis(2)

if (dis(1) <= dis(2)) then
    check = .TRUE.
else
    check = .FALSE.
endif

end subroutine overlap_check

subroutine find_same_number(na, a, nb, b, count)

implicit none

integer, intent(in) :: na, nb, b(nb)
integer, intent(inout) :: a(na), count

integer :: i, j, c(na)

c = a
count = 0
a = 0
do i = 1, na
    do j = 1, nb
        if (c(i) == b(j)) then
            count = count + 1
            a(count) = c(i)
            exit
        endif
    enddo
enddo

end subroutine find_same_number

subroutine find_nci_elem_hexa(cp, elem, nci)

implicit none

integer, intent(in) :: cp, elem(8)
integer, intent(inout) :: nci(3)

integer :: i, seq

do i = 1, 8
    if (cp == elem(i)) then
        seq = i
        exit
    endif
enddo
if (seq == 1) then
    nci = (/ 2, 4, 5 /)
elseif (seq == 2) then
    nci = (/ 1, 3, 6 /)
elseif (seq == 3) then
    nci = (/ 2, 4, 7 /)
elseif (seq == 4) then
    nci = (/ 1, 3, 8 /)
elseif (seq == 5) then
    nci = (/ 1, 6, 8 /)
elseif (seq == 6) then
    nci = (/ 2, 5, 7 /)
elseif (seq == 7) then
    nci = (/ 3, 6, 8 /)
elseif (seq == 8) then
    nci = (/ 4, 5, 7 /)
endif
nci = elem(nci)

end subroutine find_nci_elem_hexa


subroutine calc_quad_cp_dis(quad_nn, quad_node, quad_ne, quad_ele, quad_cp, quad_dis)

implicit none

integer, intent(in) :: quad_nn, quad_ne, quad_ele(8,quad_ne)
real, intent(in) :: quad_node(3,quad_nn)
real, intent(inout) :: quad_dis(quad_ne), quad_cp(3,quad_ne)

integer :: i, j
real :: elem_coord(3,8), dis

do i = 1, quad_ne
    elem_coord = quad_node(:,quad_ele(:,i))

    quad_cp(:,i) = 0.0
    do j = 1, 8
        quad_cp(:,i) = quad_cp(:,i) + elem_coord(:,j)
    enddo
    quad_cp(:,i) = quad_cp(:,i) * 0.125

    quad_dis(i) = 0.0
    do j = 1, 8
        call calc_length_3D(quad_cp(:,i), elem_coord(:,j), dis)
        if (quad_dis(i) < dis) quad_dis(i) = dis
    enddo
enddo

end subroutine calc_quad_cp_dis


subroutine calc_quad_volume(quad_nn, quad_node, quad_ne, quad_ele, quad_vol)

implicit none

integer, intent(in) :: quad_nn, quad_ne, quad_ele(8,quad_ne)
real, intent(in) :: quad_node(3,quad_nn)
real, intent(inout) :: quad_vol(quad_ne)

integer :: i
real :: target_coord(3,8)

do i = 1, quad_ne
    target_coord = quad_node(:,quad_ele(:,i))
    call calc_cubic_volume(target_coord, quad_vol(i))
enddo

end subroutine calc_quad_volume



subroutine find_other_revol_nodes(group_num, ne, elem, revol_info, revol_nodes)

implicit none

integer, intent(in) :: group_num, ne, elem(8,ne), revol_info(3)
integer, intent(inout) :: revol_nodes(revol_info(1), revol_info(3))

integer :: i, j, k, q, pre_ne, pre_layer_ne, count, cp, ce, rv_num, nci(3)
integer :: pre_elem(revol_info(1)), revol_surf_nodes(revol_info(1)), pre_layer_elem(revol_info(1)*2)
logical :: check

rv_num = revol_info(1)
revol_surf_nodes = revol_nodes(:,1)
pre_layer_ne = 0
pre_layer_elem = 0
do i = 2, revol_info(3)
    pre_ne = 0
    pre_elem = 0
    do j = 1, rv_num
        cp = revol_nodes(j,i-1)
        ce = 0
        
        ! find the element including revol_nodes(j,i-1)
        do k = 1, pre_ne
            call find_count(1, revol_nodes(j,i-1), 8, elem(:,pre_elem(k)), count)
            if (count == 1) then
                ce = pre_elem(k)
                call find_nci_elem_hexa(cp, elem(:,ce), nci)
                exit
            endif
        enddo
        
        if (ce == 0) then
            do k = 1, ne
                call find_count(1, revol_nodes(j,i-1), 8, elem(:,k), count)
                if (count == 1) then
                    check = .TRUE.
                    if (i /= 2) then
                        do q = 1, pre_layer_ne
                            if (pre_layer_elem(q) == k) then
                                check = .FALSE.
                                exit
                            endif
                        enddo
                    endif
                    if (check) then
                        call find_count(rv_num, revol_surf_nodes, 8, elem(:,k), count)
                        if (count == 4) then
                            ce = k
                            pre_ne = pre_ne + 1
                            pre_elem(pre_ne) = k
                            call find_nci_elem_hexa(cp, elem(:,k), nci)
                            exit
                        endif
                    endif
                endif
            enddo
        endif
        
        if (ce == 0) then
            write (*,*) 'check point:', cp
            stop 'Cannot find the element including check point'
        endif
        !write (*,*) 'ce, cp:' , ce, cp
        !write (*,'(8(I6,2X))') elem(:,ce)
        do k = 1, 3
            call find_count(1, nci(k), rv_num, revol_surf_nodes, count)
            if (count == 0) then
                revol_nodes(j,i) = nci(k)
                exit
            endif
        enddo
    enddo
    
    if (i /= revol_info(3)) then
        !pre_layer_ne = pre_ne
        !pre_layer_elem(1:pre_ne) = pre_elem(1:pre_ne)
        pre_layer_ne = 0
        do j = 1, ne
            call find_count(rv_num, revol_surf_nodes, 8, elem(:,j), count)
            if (count == 4) then
                pre_layer_ne = pre_layer_ne + 1
                pre_layer_elem(pre_layer_ne) = j
            endif
        enddo
        revol_surf_nodes = revol_nodes(:,i)
    endif
    !write (*,*) i, 'revol_nodes:', revol_nodes(:,i)
    !write (*,*)
enddo


end subroutine find_other_revol_nodes