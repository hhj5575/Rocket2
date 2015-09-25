Subroutine remesh_get_boundary_info_3D(ng, flag)
    
use remesh_domain
use qsort_c_module

implicit none

integer, intent(in) :: ng, flag

integer :: nn, ne, rn, bn, bfn, hn, tfn, temp_k, count, num_line, memory_j, adjn, chk_num
integer :: i, j, k, q, o, sp, lp, inn, cn
integer :: nci(4), face_nn(4,6), temp_face(4,20), over(3), this(2), elem_line(2,12), temp(2), temp_line_face(4)
integer, allocatable :: elem(:,:), line(:,:), temp_line(:,:), bound_node(:), bound_face(:,:), bound_lv(:), bound_conn(:)
integer, allocatable :: del_line(:), line_sum(:), line_hn(:), line_tfn(:), adj(:), ln(:), line_face(:,:,:)
real, allocatable :: coord(:,:)
logical :: check(2), line_check
integer :: time_values_0(8), time_values_1(8)
real :: t0, t1, t2, start_time, finish_time, event_time
!----------------------------------------------------------------------------
if (flag == 1) deallocate(Hexa(ng)%bound_node, Hexa(ng)%bound_face, Hexa(ng)%bound_conn)
nn = Hexa(ng)%nn
ne = Hexa(ng)%ne
allocate (elem(8,ne), coord(3,nn))
coord = Hexa(ng)%node
elem = Hexa(ng)%elem
    
! find line
allocate (temp_line(2,ne*12), del_line(ne*12), line_sum(ne*12))
del_line = 1
temp_line = 0
do i = 1, ne
    call extract_line_from_elem(8, 12, elem(:,i), elem_line)
    do j = 1, 12
        if (elem_line(1,j) > elem_line(2,j)) then
            temp_line(1,(i-1)*12+j) = elem_line(2,j)
            temp_line(2,(i-1)*12+j) = elem_line(1,j)
        else
            temp_line(:,(i-1)*12+j) = elem_line(:,j)
        endif
    enddo
enddo
    
! descending sort
call QsortC(temp_line(:,1:ne*12), 2)
!do i = 1, ne*12-1
!    do j = i+1, ne*12
!        if ( temp_line(1,i) > temp_line(1,j) ) then
!            temp = temp_line(:,i)
!            temp_line(:,i) = temp_line(:,j)
!            temp_line(:,j) = temp
!        endif
!    enddo
!enddo
    
do i = 1, ne*12-1
    if (del_line(i) == 1) then
        do j = i+1, ne*12
            if (del_line(j) == 1 .and. temp_line(1,i) == temp_line(1,j)) then
                call find_count(2, temp_line(:,i), 2, temp_line(:,j), count)
                if (count == 2) del_line(j) = 0
            elseif (temp_line(1,i) /= temp_line(1,j)) then
                exit
            endif
        enddo
    endif
enddo
    
deallocate (line_sum)
num_line = sum(del_line)
allocate (line(2,num_line))
num_line = 0
do i = 1, ne*12
    if (del_line(i) == 1) then
        num_line = num_line + 1
        if (temp_line(1,i) > temp_line(2,i)) then
            line(1,num_line) = temp_line(2,i)
            line(2,num_line) = temp_line(1,i)
        else
            line(:,num_line) = temp_line(:,i)
        endif
        !line_sum(num_line) = temp_line(1,i)+temp_line(2,i)
    endif
enddo
deallocate (temp_line, del_line)
    
allocate (adj(nn+1))
adj(1) = 1
chk_num = 1
do i = 2, num_line
    if (chk_num /= line(1,i)) then
        if (chk_num+1 /= line(1,i)) then
            do j = chk_num+1, line(1,i)-1
                adj(j) = i
            enddo
        endif
        chk_num = line(1,i)
        adj(chk_num) = i
    endif
enddo
do i = chk_num+1, nn+1
    adj(i) = num_line + 1
enddo
    
! find line_tn, line_tfn
allocate (line_hn(num_line), line_tfn(num_line), line_face(20,4,num_line))
line_hn = 0
line_tfn = 0
line_face = 0
do i = 1, ne
    call extract_line_from_elem(8, 12, elem(:,i), elem_line)
    call make_face_nn(8, 4, elem(:,i), face_nn)
    do j = 1, 12
        if (elem_line(1,j) > elem_line(2,j)) then
            temp(1) = elem_line(1,j)
            elem_line(1,j) = elem_line(2,j)
            elem_line(2,j) = temp(1)
        endif
        sp = adj(elem_line(1,j))
        lp = adj(elem_line(1,j)+1)-1
        do k = sp, lp
            call find_count(2, elem_line(:,j), 2, line(:,k), count)
            if (count == 2) then
                line_tfn(k) = line_tfn(k) + 1
                exit
            endif
        enddo
    enddo
        
    do j = 1, 6
        do k = 1, 4
            if (k == 4) then
                this = (/ face_nn(k,j), face_nn(1,j) /)
            else
                this = (/ face_nn(k,j), face_nn(k+1,j) /)
            endif
            if (this(1) > this(2)) then
                temp(1) = this(1)
                this(1) = this(2)
                this(2) = temp(1)
            endif
            sp = adj(this(1))
            lp = adj(this(1)+1)-1
            do q = sp, lp
                call find_count(2, this, 2, line(:,q), count)
                if (count == 2) then
                    check(1) = .true.
                    do o = 1, line_hn(q)
                        temp_line_face = line_face(o,:,q)
                        call find_count(4,temp_line_face,4,face_nn(:,j),count)
                        if (count == 4) then
                            check(1) = .false.
                            exit
                        endif
                    enddo
                    if (check(1)) then
                        line_hn(q) = line_hn(q) + 1
                        line_face(line_hn(q),:,q) = face_nn(:,j)
                        exit
                    endif
                endif
            enddo
        enddo
    enddo
enddo
deallocate (line_face) 
    
! find bound_node
allocate (bound_node(nn), bound_conn(nn))
bound_conn = 0
bn = 0
do i = 1, num_line
    if ( line_hn(i) /= line_tfn(i) ) then
        do j = 1, 2
            check(1) = .TRUE.
            do k = 1, bn 
                if (bound_node(k) == line(j,i)) then
                    check(1) = .FALSE.
                    exit
                endif
            enddo
            if ( check(1) ) then
                bn = bn + 1
                bound_node(bn) = line(j,i)
            endif
        enddo
    endif
enddo
!write (*,*) rn, 'domain, num_line, num_bound_node: ', num_line, bn
call ascending_sort_int(bn, bound_node(1:bn))
do i = 1, bn
    bound_conn(bound_node(i)) = i
enddo
allocate (Hexa(ng)%bound_node(bn), Hexa(ng)%bound_conn(nn))
Hexa(ng)%bn = bn
Hexa(ng)%bound_node(:) = bound_node(1:bn)
Hexa(ng)%bound_conn = bound_conn
deallocate (line, adj, line_hn, line_tfn, bound_conn)
    
! find boundary faces
write (*,*) 'find boundary faces'
bfn = 0
allocate (bound_face(4,ne*4))
do i = 1, ne
    call make_face_nn(8, 4, elem(:,i), face_nn)
    do j = 1, 6
        call find_count(bn, bound_node(1:bn), 4, face_nn(:,j), count)
        if (count == 4) then
            check(1) = .true.
            do k = 1, bfn
                call find_count(4, bound_face(:,k), 4, face_nn(:,j), count)
                if (count >= 3) then
                    check(1) = .false.
                    temp_k = k
                    exit
                endif
            enddo
            
            if (check(1)) then
                bfn = bfn + 1
                bound_face(:,bfn) = face_nn(:,j)
            else
                do k = temp_k, bfn-1
                    bound_face(:,k) = bound_face(:,k+1)
                enddo
                bfn = bfn - 1
            endif
        endif
    enddo
enddo

write (*,*) 'check_bound_face'
call check_bound_face(nn, coord, bfn, bound_face(:,1:bfn))
allocate (Hexa(ng)%bound_face(4,bfn))
Hexa(ng)%bfn = bfn
Hexa(ng)%bound_face = bound_face(:,1:bfn)

deallocate (elem, coord, bound_node, bound_face)

End subroutine remesh_get_boundary_info_3D

!============================================================================

subroutine array_include_check(chk_val, arr_num, arr, check)

implicit none

integer, intent(in) :: chk_val, arr_num, arr(arr_num)
logical, intent(out) :: check

integer :: i

check = .false.
do i = 1, arr_num
    if (arr(i) == chk_val) then
        check = .true.
        exit
    endif
enddo

end subroutine array_include_check

subroutine move_point_near_ele(nn, coord, temp, cp, ele)

implicit none

integer :: nn, ele(3)
real :: coord(3,nn), cp(3), temp(3)

integer :: i, this(2), temp_i
real :: p(3,3), np(3), min_dis, dis, t, temp_cp(3)
logical :: check

p(:,1) = coord(:,ele(1))
p(:,2) = coord(:,ele(2))
p(:,3) = coord(:,ele(3))

min_dis = 10e+8
temp_i = 0
do i = 1, 3
    if (i == 3) then
        this = (/ 3, 1 /)
    else
        this = (/ i, i+1 /)
    endif
    call calc_perpendicular_point(temp, p(:,this(1)), p(:,this(2)), np, t, check)
    if (t <= 1.01 .AND. t >= -0.01) then
        call calc_length_3D(temp, np, dis)
        if (min_dis > dis) then
            min_dis = dis
            temp_cp = np
            temp_i = i
        endif
    else
        call calc_length_3D(temp, p(:,this(1)), dis)
        if (min_dis > dis) then
            min_dis = dis
            temp_cp = p(:,this(1))
            temp_i = i
        endif
        call calc_length_3D(temp, p(:,this(2)), dis)
        if (min_dis > dis) then
            min_dis = dis
            temp_cp = p(:,this(2))
            temp_i = i
        endif
    endif
enddo
cp = temp_cp

end subroutine move_point_near_ele

!****************************************************************************
!
!  PROGRAM : calc_perpendicular_point
!
!  PURPOSE : 임의의 점에서 선분으로 수선의 위치 계산
!
!****************************************************************************
subroutine calc_perpendicular_point(a, b, c, d, t, judgment)

implicit none
!a: 임의의 점, b, c: 선분을 이루는 두 점
real, intent(in) :: a(3), b(3), c(3)
real, intent(inout) :: d(3), t
logical, intent(inout) :: judgment

real :: l(3), vec(3), tol

!l(1) = (c(1)-b(1))**2+(c(2)-b(2))**2+(c(3)-b(3))**2
!l(2) = (a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2
!l(3) = (a(1)-c(1))**2+(a(2)-c(2))**2+(a(3)-c(3))**2
!t = ((l(2)-l(3))/l(1) + 1)*0.5
tol = 0.01

vec = c - b
t = vec(1)*(a(1)-b(1))+vec(2)*(a(2)-b(2))+vec(3)*(a(3)-b(3))
t = t / (vec(1)**2+vec(2)**2+vec(3)**2)

if (t < 1. + tol .AND. t > 0. - tol) then
    judgment = .TRUE.
else
    judgment = .FALSE.
endif
d = b+vec*t

end subroutine calc_perpendicular_point


Subroutine check_inner_area(p, cp, check, area)

implicit none

real, intent(in) :: p(3,3), cp(3)
real, intent(inout) :: area
logical, intent(inout) :: check

integer :: i, i_r
integer :: pl_check(3)
real :: sum_area, ori_area
real :: dis(3), s, tol
!----------------------------------------------------------------------------
call calc_length_3D(p(:,1), p(:,2), dis(1))
call calc_length_3D(p(:,2), p(:,3), dis(2))
call calc_length_3D(p(:,3), p(:,1), dis(3))
s = (dis(1)+dis(2)+dis(3))/2.
if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
    ori_area = sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
else
    ori_area = 0.0
endif

sum_area = 0.
do i = 1, 3
    i_r = i+1
    if ( i == 3 ) i_r = 1
    call calc_length_3D(p(:,i), p(:,i_r), dis(1))
    call calc_length_3D(p(:,i), cp, dis(2))
    call calc_length_3D(p(:,i_r), cp, dis(3))
    s = (dis(1)+dis(2)+dis(3))/2.
    if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
        sum_area = sum_area + sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
    endif
enddo

tol = 0.05
area = abs(sum_area-ori_area)
if ( abs(sum_area-ori_area) < ori_area*tol ) then
    check = .TRUE.
else
    check = .FALSE.
endif

end subroutine check_inner_area

!****************************************************************************

subroutine find_match_node2(num, cp, target_nn, target_coord, tol)

implicit none

integer, intent(in) :: target_nn
real, intent(in) :: cp(3), target_coord(3,target_nn), tol
integer, intent(inout) :: num

integer :: i, temp_num
real :: min_dis, dis

min_dis = 10e8
temp_num = 0
num = 0
do i = 1, target_nn
    call calc_length_3D(cp, target_coord(:,i), dis)
    if (dis < tol) then
        num = i
        exit
    elseif (min_dis > dis) then
        temp_num = i
        min_dis = dis
    endif
enddo
if (num == 0) then
    num = temp_num
endif

end subroutine find_match_node2


subroutine find_cons_number(unit_num, unit_num_inte, cons_num, fix_cons_num, ori_nn, ori_coord, tol_nn, tol_coord, prop_cons, tol)

implicit none

integer, intent(in) :: unit_num, unit_num_inte, tol_nn, ori_nn
real, intent(in) :: tol_coord(3,tol_nn), ori_coord(3,ori_nn), tol
integer, intent(inout) :: cons_num, fix_cons_num, prop_cons(4,tol_nn)

integer :: i, j, k, q, cons_set_num
real :: pl(4), cp(3), dis, p(3,3)
integer, allocatable :: cons(:,:)
real, allocatable :: cons_set_coord(:,:)
logical :: check
character(len=5 ) :: keyword

rewind(unit_num_inte)
do 
    read(unit_num_inte,*) keyword
    if (keyword == '*reco') then
        backspace(unit_num_inte)
        read(unit_num_inte,*) keyword, fix_cons_num, cons_set_num
        allocate (cons_set_coord(3,cons_set_num*3), cons(3,cons_set_num))
        do i = 1, cons_set_num
            read(unit_num_inte,*) cons(:,i)
            do j = 1, 3
                read(unit_num_inte,*) cons_set_coord(:,(i-1)*3+j)
            enddo
        enddo
        exit
    endif
enddo

read(unit_num,*) keyword, cons_num
cons_num = fix_cons_num
do i = 1, cons_num
    read(unit_num,*) prop_cons(:,i)
    call find_match_node2(prop_cons(1,i), ori_coord(:,prop_cons(1,i)), tol_nn, tol_coord, tol)
enddo

do i = 1, cons_set_num
    do j = 1, 3
        p(:,j) = cons_set_coord(:,(i-1)*3+j)
    enddo
    call calc_plane(p(:,1), p(:,2), p(:,3), pl)
    do j = 1, tol_nn
        cp = tol_coord(:,j)
        dis = abs(pl(1)*cp(1) + pl(2)*cp(2) + pl(3)*cp(3) + pl(4))
        if (dis < tol) then
            check = .TRUE.
            do k = 1, cons_num
                if (prop_cons(1,k) == j) then
                    do q = 1, 3
                        if (cons(q,i) == 1) then
                            prop_cons(q+1,k) = 1
                        endif
                    enddo
                    check = .FALSE.
                    exit
                endif
            enddo
            if (check) then
                cons_num = cons_num + 1
                prop_cons(1,cons_num) = j
                prop_cons(2:4,cons_num) = cons(:,i)
            endif
        endif
    enddo
enddo
deallocate (cons_set_coord, cons)

end subroutine find_cons_number


subroutine write_inte_file_3D(unit_num_inte, inte_unit_num, p2_csn, csn, num_rn, num_group, temp_val, ori_ng, mesh_flag)

use surface_domain
use remesh_domain
implicit none

integer, intent(in) :: unit_num_inte, inte_unit_num, num_rn, num_group, ori_ng
integer, intent(in) :: p2_csn, csn(2), temp_val(csn(2)+1,2*num_group), mesh_flag(ori_ng)

character(len=5) :: keyword
integer :: i, j, rn, num, ng, pre_bn, cp, status, flag, del_ng
integer :: cons_num(2), cons(3), line_num(8)
real :: p(3)

rn = 1
write (inte_unit_num,*) '*divd, ', num_rn
write (inte_unit_num,*) '*subd, ', num_group
write (inte_unit_num,*) '*rety, ', remesh_type_number
write (inte_unit_num,*) '*redi, ', remesh_dis
write (inte_unit_num,*) ' '

write (inte_unit_num,*) '*p2cn, ', p2_csn
write (inte_unit_num,*) '*cont, ', csn
do i = 1, num_group
    write (inte_unit_num,'(100(I8))') temp_val(:,(i-1)*2+1)
    write (inte_unit_num,'(100(I8))') temp_val(:,(i-1)*2+2)
enddo
write (inte_unit_num,*) ' '

rewind(unit_num_inte)
do 
    read(unit_num_inte,*) keyword
    if (keyword == '*reco') then
        backspace(unit_num_inte)
        read(unit_num_inte,*) keyword, cons_num
        write (inte_unit_num,*) keyword, cons_num
        do i = 1, cons_num(2)
            read(unit_num_inte,*) cons
            write (inte_unit_num,*) cons
            do j = 1, 3
                read(unit_num_inte,*) p
                write (inte_unit_num,*) p
            enddo
        enddo
        write (inte_unit_num,*)
        exit
    endif
enddo

del_ng = 0
rewind(unit_num_inte)
do 
    read(unit_num_inte,*,iostat=status) keyword
    if (status == -1) then
        exit
    elseif (keyword == '*end') then
        write (inte_unit_num,*) keyword
        exit
    elseif (keyword == '*remp') then
        backspace(unit_num_inte)
        read(unit_num_inte,*) keyword, ng
        if (mesh_flag(ng) /= 3) then
            write (inte_unit_num,*) keyword, ng-del_ng
            read(unit_num_inte,*) keyword, num
            write (inte_unit_num,*) keyword, num
            pre_bn = mesh_set(1)%sub_num(ng,3) - 1
            write (*,*) '*op2n'
            do i = 1, num
                read(unit_num_inte,*) p
                cp = mesh_set(rn)%bound_node(remesh_add(ng)%ori_p2(i)+pre_bn)
                write (inte_unit_num,'(3(F15.8,2X))') mesh_set(rn)%node_coord(:,cp)
                write (*,'(3(F15.8,2X))') mesh_set(rn)%node_coord(:,cp)
            enddo
            read(unit_num_inte,*) keyword, num
            write (inte_unit_num,*) keyword, num
            write (*,*) '*ap2n'
            do i = 1, num
                read(unit_num_inte,*) flag
                backspace(unit_num_inte)
                if (flag == 1) then
                    read(unit_num_inte,*) flag, p, line_num(1:2)
                    cp = mesh_set(rn)%bound_node(remesh_add(ng)%p2(1,i)+pre_bn)
                    write (inte_unit_num,'(I3,3(F15.8,2X),2(I3,1X))') flag, mesh_set(rn)%node_coord(:,cp), line_num(1:2)
                    write (*,'(3(F15.8,2X))') mesh_set(rn)%node_coord(:,cp)
                elseif (flag == 2) then
                    read(unit_num_inte,*) flag, line_num(1:5)
                    write (inte_unit_num,'(I3,3(I7,1X),I3,1X,I3)') flag, line_num(1:5)
                endif
            enddo
            read(unit_num_inte,*) keyword, num
            write (inte_unit_num,*) keyword, num
            do i = 1, num
                read(unit_num_inte,*) line_num
                write (inte_unit_num,'(8(I4,1X))') line_num
            enddo
            read(unit_num_inte,*) keyword, num
            write (inte_unit_num,*) keyword, num
            do i = 1, num
                read(unit_num_inte,*) line_num(1:4)
                write (inte_unit_num,'(4(I4,1X))') line_num(1:4)
            enddo
            read(unit_num_inte,*) keyword, num
            write (inte_unit_num,*) keyword, num
            do i = 1, num
                read(unit_num_inte,*) line_num(1:3)
                write (inte_unit_num,'(I3,1X,I7,1X,I7)') line_num(1:3)
            enddo
        else
            del_ng = del_ng + 1
        endif
    endif
enddo
        
end subroutine write_inte_file_3D


subroutine find_match_node_num(cp, nn, temp_coord, num, tol, flag)

implicit none

integer, intent(   in) :: nn, flag
integer, intent(inout) :: num
real,    intent(   in) :: cp(3), temp_coord(3,nn), tol

integer :: sp, lp, i, temp_i
real :: dis, min_dis 

if (flag == 1) then
    sp = 1
    lp = nn
elseif (flag == -1) then
    sp = nn
    lp = 1
endif

min_dis = 10e8
num = 0
temp_i = 0
do i = sp, lp, flag
    call calc_length_3D(cp, temp_coord(:,i), dis)
    if (dis < tol) then
        num = i
        exit
    elseif (min_dis > dis) then
        temp_i = i
        min_dis = dis
    endif
enddo
if (num == 0) then
    num = temp_i
endif

end subroutine find_match_node_num

subroutine Hexa_mesh_rearrange(nn, node, ne, elem)

implicit none

integer, intent(in) :: nn, ne
real, intent(in) :: node(3,nn)
integer, intent(inout) :: elem(8,ne)

integer :: i, temp_elem(8)
real :: vec(3,4), dot

do i = 1, ne
    vec(:,1) = node(:,elem(2,i)) - node(:,elem(1,i))
    vec(:,2) = node(:,elem(4,i)) - node(:,elem(1,i))
    vec(:,3) = node(:,elem(5,i)) - node(:,elem(1,i))
    
    vec(1,4) = vec(2,1)*vec(3,2) - vec(2,2)*vec(3,1)
    vec(2,4) = -(vec(1,1)*vec(3,2) - vec(1,2)*vec(3,1))
    vec(3,4) = vec(1,1)*vec(2,2) - vec(1,2)*vec(2,1)
    
    dot = vec(1,3)*vec(1,4)+vec(2,3)*vec(2,4)+vec(3,3)*vec(3,4)
    if (dot <= 0.0) then
        !write (*,*) i, ': inner product of element is minus'
        temp_elem = elem(:,i)
        elem(2,i) = temp_elem(4)
        elem(4,i) = temp_elem(2)
        elem(6,i) = temp_elem(8)
        elem(8,i) = temp_elem(6)
    endif
enddo

end subroutine Hexa_mesh_rearrange


subroutine tri_or_add_node(Tri_nn, Tri_node, add_nn, add_node, num, p)

implicit none

integer, intent(in) :: tri_nn, add_nn, num
real, intent(in) :: tri_node(3,tri_nn), add_node(3,add_nn)
real, intent(inout) :: p(3)

if (num > tri_nn) then
    p = add_node(:,num-tri_nn)
else
    p = tri_node(:,num)
endif

end subroutine tri_or_add_node