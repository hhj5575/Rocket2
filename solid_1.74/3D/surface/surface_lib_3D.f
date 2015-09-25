Subroutine get_boundary_info_3D(flag, rn, num_group)
    
use surface_domain
use qsort_c_module

implicit none

integer, intent(in) :: flag, rn, num_group

integer :: nn, ne, bn, bfn, hn, ng, tfn, temp_k, count, num_line, memory_j, adjn, chk_num, sp(2), lp(2)
integer :: i, j, k, q, o, inn, cn
integer :: nci(4), face_nn(4,6), temp_face(4,20), over(3), this(2), elem_line(2,12), temp(2), temp_line_face(4)
integer, allocatable :: line(:,:), temp_line(:,:), bound_node(:), bound_face(:,:), bound_lv(:)
integer, allocatable :: del_line(:), line_sum(:), line_hn(:), line_tfn(:), adj(:), ln(:), line_face(:,:,:)
integer, allocatable :: elem(:,:), inn_nci_num(:), inn_nci(:,:), inner_node(:)
real, allocatable :: coord(:,:)
integer :: rn_bn, rn_bfn, rn_inn
integer, allocatable :: rn_bound_node(:), rn_bound_face(:,:), rn_bound_lv(:), inv_bound_node(:)
integer, allocatable :: rn_inner_node(:), rn_inn_nci_num(:), rn_inn_nci(:,:)
logical :: check(2), line_check
integer :: time_values_0(8), time_values_1(8)
real :: t0, t1, t2, start_time, finish_time, event_time
!----------------------------------------------------------------------------
!t0 = 0.0
!t1 = 0.0
!call cpu_time(t0)
!call date_and_time(values=time_values_0)
!start_time = time_values_0(5)*3600 + time_values_0(6)*60 + time_values_0(7)*1 + time_values_0(8)*0.001 

if (flag == 1) deallocate(mesh_set(rn)%bound_node, mesh_set(rn)%bound_face)
allocate(rn_bound_node(mesh_set(rn)%nn))
allocate(rn_bound_face(4,mesh_set(rn)%ne*6))
allocate(rn_bound_lv(mesh_set(rn)%nn))
allocate(inv_bound_node(mesh_set(rn)%nn))
allocate(rn_inner_node(mesh_set(rn)%nn))
allocate(rn_inn_nci_num(mesh_set(rn)%nn))
allocate(rn_inn_nci(20,mesh_set(rn)%nn))
rn_bn = 0;  rn_bfn = 0;  rn_inn = 0
rn_bound_lv = 0
do ng = 1, num_group
    if (rn == 1) then
        sp = mesh_set(rn)%sub_num(ng,1:2)
        lp = mesh_set(rn)%sub_num(ng+1,1:2) - 1
        nn = lp(1) - sp(1) + 1
        ne = lp(2) - sp(2) + 1
        write (*,*) 'rn, nn, ne:', rn, nn, ne
        allocate (elem(8,ne), coord(3,nn))
        coord = mesh_set(rn)%node_coord(:,sp(1):lp(1))
        elem = mesh_set(rn)%element(:,sp(2):lp(2))
        if (ng /= 1) elem = elem - (mesh_set(rn)%sub_num(ng,1)-1)
        !call write_tecplot_3D(nn, coord, ne, elem, rn, ng, 7)
    else
        nn = mesh_set(rn)%nn
        ne = mesh_set(rn)%ne
        write (*,*) 'rn, nn, ne:', rn, nn, ne
        allocate (elem(8,ne), coord(3,nn))
        coord = mesh_set(rn)%node_coord
        elem = mesh_set(rn)%element
    endif
    
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
            sp(1) = adj(elem_line(1,j))
            lp(1) = adj(elem_line(1,j)+1)-1
            do k = sp(1), lp(1)
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
                sp(1) = adj(this(1))
                lp(1) = adj(this(1)+1)-1
                do q = sp(1), lp(1)
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
    allocate (bound_node(nn))
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
    
    inv_bound_node = 0
    do i = 1, bn
        inv_bound_node(bound_node(i)) = i+rn_bn
    enddo
    if (rn == 1) then
        rn_bound_node(rn_bn+1:rn_bn+bn) = bound_node(1:bn)+mesh_set(rn)%sub_num(ng,1)-1
    else
        rn_bound_node(rn_bn+1:rn_bn+bn) = bound_node(1:bn)
    endif
    rn_bn = rn_bn + bn
    if (rn == 1) mesh_set(rn)%sub_num(ng+1,3) = rn_bn+1
    deallocate (line, adj, line_hn, line_tfn)
    
    ! find boundary faces
    write (*,*) 'find boundary faces'
    bfn = 0
    allocate (bound_face(4,ne*6))
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
    
    !if (rn == 1) then
    !    write (*,*) 'rearrange_elem'
    !    call rearrange_elem(bfn, bound_face(:,1:bfn))
    !    call check_bound_face(nn, coord, bfn, bound_face(:,1:bfn))
    !endif
    do i = 1, bfn
        if (rn == 1) then
            !rn_bound_face(:,rn_bfn+i) = inv_bound_node(bound_face(:,i))+mesh_set(rn)%sub_num(ng,3)-1
            rn_bound_face(:,rn_bfn+i) = inv_bound_node(bound_face(:,i))
        else
            rn_bound_face(:,rn_bfn+i) = inv_bound_node(bound_face(:,i))
        endif
    enddo
    rn_bfn = rn_bfn + bfn
    if (rn == 1) mesh_set(rn)%sub_num(ng+1,4) = rn_bfn+1
    
    deallocate (elem, coord, bound_node, bound_face)
enddo

allocate (mesh_set(rn)%bound_node(rn_bn))
mesh_set(rn)%bn = rn_bn
mesh_set(rn)%bound_node = rn_bound_node(1:rn_bn)

allocate (mesh_set(rn)%bound_face(4,rn_bfn))
mesh_set(rn)%bfn = rn_bfn
mesh_set(rn)%bound_face = rn_bound_face(:,1:rn_bfn)

if (rn == 1) then
    allocate (bound_lv(rn_bn))
    call set_bound_node_lv(mesh_set(rn)%nn, rn_bn, rn_bfn, mesh_set(rn)%node_coord, rn_bound_node(1:rn_bn), rn_bound_face(:,1:rn_bfn), bound_lv)
    !rn_bound_lv(mesh_set(rn)%sub_num(ng,3):mesh_set(rn)%sub_num(ng+1,3)-1) = bound_lv
    !deallocate (bound_lv)
        
    inn = mesh_set(rn)%nn - rn_bn
    !allocate(line_sum(mesh_set(rn)%nn))
    allocate(inner_node(inn))
    allocate(inn_nci_num(inn))
    allocate(inn_nci(20,inn))
    inn_nci_num = 0
    inn_nci = 0
    inn = 0
    !line_sum = 1
    do i = 1, mesh_set(rn)%nn
        check(1) = .true.
        do j = 1, rn_bn
            if (i == rn_bound_node(j)) then
                check(1) = .false.
                !line_sum(i) = 0
                exit
            endif
        enddo
        if (check(1)) then
            inn = inn + 1
            inner_node(inn) = i
        endif
    enddo
    
    do i = inn, 1, -1
        cn = inner_node(i)
        do j = 1, mesh_set(rn)%ne
            do k = 1, 8
                if (mesh_set(rn)%element(k,j) == cn) then
                    inn_nci_num(i) = inn_nci_num(i) + 1
                    inn_nci(inn_nci_num(i),i) = j
                    exit
                endif
            enddo
        enddo
    enddo
    !rn_inner_node(rn_inn+1:rn_inn+inn) = inner_node(1:inn)+mesh_set(rn)%sub_num(ng,1)-1
    !rn_inn_nci_num(rn_inn+1:rn_inn+inn) = inn_nci_num(1:inn)
    !rn_inn_nci(:,rn_inn+1:rn_inn+inn) = inn_nci(:,1:inn)+mesh_set(rn)%sub_num(ng,1)-1
    !rn_inn = rn_inn + inn
    !mesh_set(rn)%sub_num(ng+1,5) = rn_inn+1
    !deallocate(inner_node, inn_nci_num, inn_nci)
    
    allocate(mesh_set(rn)%bound_lv(rn_bn))
    !mesh_set(rn)%bound_lv = rn_bound_lv(1:rn_bn)
    mesh_set(rn)%bound_lv = bound_lv(1:rn_bn)
    
    allocate(mesh_set(rn)%inner_node(inn))
    allocate(mesh_set(rn)%inn_nci_num(inn))
    allocate(mesh_set(rn)%inn_nci(20,inn))
    !mesh_set(rn)%inn = rn_inn
    !mesh_set(rn)%inner_node = rn_inner_node(1:rn_inn)
    !mesh_set(rn)%inn_nci_num = rn_inn_nci_num(1:rn_inn)
    !mesh_set(rn)%inn_nci = rn_inn_nci(:,1:rn_inn)
    mesh_set(rn)%inn = inn
    mesh_set(rn)%inner_node = inner_node(1:inn)
    mesh_set(rn)%inn_nci_num = inn_nci_num(1:inn)
    mesh_set(rn)%inn_nci = inn_nci(:,1:rn_inn)
    write (*,*) 'sub_nn:', mesh_set(rn)%sub_num(1:ng,1)
    write (*,*) 'sub_ne:', mesh_set(rn)%sub_num(1:ng,2)
    write (*,*) 'sub_bn:', mesh_set(rn)%sub_num(1:ng,3)
    write (*,*) 'sub_bfn:', mesh_set(rn)%sub_num(1:ng,4)
    !write (*,*) 'sub_inn:', mesh_set(rn)%sub_num(1:ng,5)
    open (Unit=37, File='check_bound_lv.plt', STATUS='replace', ACTION='write')
    Write (37,'(A)') 'TITLE="Check face"'
    Write (37,'(A,A)') 'VARIABLES= "x", "y", "z", "lv"'
    Write (37,'(A,I5,A,I5, A)') 'ZONE N =', mesh_set(rn)%bn , ', E =', mesh_set(rn)%bfn, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
    do i = 1, rn_bn
        write (37,*) mesh_set(rn)%node_coord(:,rn_bound_node(i)), bound_lv(i)
    enddo
    do i = 1, rn_bfn
        write (37,*) rn_bound_face(:,i)
    enddo
    close(37)
    deallocate(bound_lv, inner_node, inn_nci_num, inn_nci)
endif
deallocate (rn_bound_node, rn_bound_face, inv_bound_node)
deallocate (rn_bound_lv, rn_inner_node, rn_inn_nci_num, rn_inn_nci)

!if (rn == 1) then
!    open (Unit=37, File='check_face_rn01.plt', STATUS='replace', ACTION='write')
!else
!    open (Unit=37, File='check_face_rn02.plt', STATUS='replace', ACTION='write')
!endif
!Write (37,'(A)') 'TITLE="Check face"'
!Write (37,'(A,A)') 'VARIABLES= "x", "y", "z"'
!Write (37,'(A,I5,A,I5, A)') 'ZONE N =', mesh_set(rn)%bn , ', E =', mesh_set(rn)%bfn, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
!do i = 1, mesh_set(rn)%bn
!    write (37,*) mesh_set(rn)%node_coord(:,mesh_set(rn)%bound_node(i))
!enddo
!do i = 1, mesh_set(rn)%bfn
!    write (37,*) mesh_set(rn)%bound_face(:,i)
!enddo
    
End subroutine get_boundary_info_3D
!============================================================================

subroutine make_face_nn(elem_nodes, fn, ele, face_nn)

implicit none

integer, intent(in) :: elem_nodes, fn, ele(elem_nodes)
integer, intent(inout) :: face_nn(fn,elem_nodes/2+2)
!------------------------------------------------------------------------------

if (elem_nodes == 8) then
    face_nn(:,1) = (/ ele(1), ele(2), ele(3), ele(4) /)
    face_nn(:,2) = (/ ele(5), ele(8), ele(7), ele(6) /)
    face_nn(:,3) = (/ ele(1), ele(5), ele(6), ele(2) /)
    face_nn(:,4) = (/ ele(2), ele(6), ele(7), ele(3) /)
    face_nn(:,5) = (/ ele(3), ele(7), ele(8), ele(4) /)
    face_nn(:,6) = (/ ele(1), ele(4), ele(8), ele(5) /)
elseif (elem_nodes == 4) then
    face_nn(:,1) = (/ ele(1), ele(2), ele(3) /)
    face_nn(:,2) = (/ ele(1), ele(2), ele(4) /)
    face_nn(:,3) = (/ ele(2), ele(3), ele(4) /)
    face_nn(:,4) = (/ ele(3), ele(1), ele(4) /)
endif

end subroutine make_face_nn
!============================================================================

subroutine find_count(na, a, nb, b, count)

implicit none

integer, intent(in) :: na, nb
integer, intent(in) :: a(na), b(nb)
integer, intent(inout) :: count

integer :: i, j, num
!----------------------------------------------------------------------------
count = 0
do i = 1, na
    do j = 1, nb
        if ( a(i) == b(j) ) then
            count = count + 1 
            exit
        endif
    enddo
enddo

end subroutine find_count
!============================================================================

subroutine ascending_sort_int(na, a)

implicit none

integer :: na, a(na)

integer :: i, j, temp
!----------------------------------------------------------------------------
do i = 1, na
    do j = i+1, na
        if ( a(i) < a(j) ) then
            temp = a(i)
            a(i) = a(j)
            a(j) = temp
        endif
    enddo
enddo

end subroutine ascending_sort_int
!============================================================================
subroutine extract_line_from_elem(elem_nodes, mesh_type, elem, elem_line)

implicit none

integer, intent(in) :: elem_nodes, mesh_type, elem(elem_nodes)
integer, intent(inout) :: elem_line(2,mesh_type)

integer :: i, num

if (elem_nodes == 4) then
    num = 3
else
    num = 4
endif
do i = 1, num
    elem_line(:,i) = (/ i, i+1 /)
    if (i == num) elem_line(2,i) = 1
enddo
if (elem_nodes == 4) then
    do i = 1, num
        elem_line(:,i+num) = (/ i, 4 /)
    enddo
else
    do i = 1, num
        elem_line(:,i+num) = (/ i+num, i+1+num /)
        elem_line(:,i+2*num) = (/ i, i+num /)
        if (i == num) elem_line(2,i+num) = 5
    enddo
endif

do i = 1, mesh_type
    elem_line(:,i) = elem(elem_line(:,i))
enddo


end subroutine extract_line_from_elem

!****************************************************************************
subroutine near_plane(b_p, np, tol, check)

implicit none

real, intent(in) :: b_p(3,4), np(3), tol
logical, intent(inout) :: check

real :: dis, pl(4), p(3,2), cp(3), t
logical :: parallel

check = .false.
call calc_plane(b_p(:,1), b_p(:,2), b_p(:,3), pl)
dis = abs(pl(1)*np(1)+pl(2)*np(2)+pl(3)*np(3)+pl(4))
dis = dis/sqrt(pl(1)**2.0+pl(2)**2.0+pl(3)**2.0)

if (dis <= tol) then
    p(:,1) = np + pl(1:3)*tol*10.0
    p(:,2) = np - pl(1:3)*tol*10.0
    call cross_point_plane(pl, p(:,1), p(:,2), cp, t, parallel)
    if ( parallel == .FALSE. ) then
        call check_inner_face(b_p, cp, check)
    endif
endif

end subroutine near_plane
!****************************************************************************
Subroutine calc_plane(a, b, c, pl)

implicit none

real, intent(in) :: a(3), b(3), c(3)
real, intent(inout) :: pl(4)

real :: normal
integer :: i, j

!----------------------------------------------------------------------------

pl(1) = (b(2)-a(2))*(c(3)-a(3)) - (b(3)-a(3))*(c(2)-a(2))
pl(2) = -( (b(1)-a(1))*(c(3)-a(3)) - (b(3)-a(3))*(c(1)-a(1)) )
pl(3) = (b(1)-a(1))*(c(2)-a(2)) - (b(2)-a(2))*(c(1)-a(1))
pl(4) = -(pl(1)*a(1) + pl(2)*a(2) + pl(3)*a(3))

if ( (pl(1) /= 0.0) .OR. (pl(2) /= 0.0) .OR. (pl(3) /= 0.0) ) then
	normal = sqrt(pl(1)**2 + pl(2)**2 + pl(3)**2)
	pl = pl / normal
endif

end subroutine calc_plane
!****************************************************************************
!
!  PROGRAM : Cross_point
!
!  PURPOSE : 평면과 벡터(선분)의 교차점
!
!****************************************************************************
Subroutine cross_point_plane(pl, cri_coord, ref_coord, cross_p, t, parallel)

implicit none

real, intent(in) :: pl(4), cri_coord(3), ref_coord(3)
real, intent(inout) :: cross_p(3)
logical, intent(inout) :: parallel

real :: t, under_t, tol, vec(3), dis
!----------------------------------------------------------------------------
call calc_length_3D(cri_coord, ref_coord, dis)
tol = dis * 0.01
vec = ref_coord - cri_coord

t = -(pl(1)*cri_coord(1) + pl(2)*cri_coord(2) + pl(3)*cri_coord(3) + pl(4))
under_t = (pl(1)*vec(1) + pl(2)*vec(2) + pl(3)*vec(3))

parallel = .FALSE.
if ( under_t <= 0. + tol .AND. under_t >= 0. - tol ) then
    parallel = .TRUE.
	cross_p = (/ 0.0, 0.0, 0.0 /)
else
    t = t / under_t
    if ( t > 1 + tol .OR. t < 0 - tol ) parallel = .TRUE.
    cross_p = vec*t + cri_coord
endif

end subroutine cross_point_plane
!****************************************************************************
!
!  PROGRAM : inner_check_area
!
!  PURPOSE : 평면의 내부 node인지 검사(넓이로 검사)
!
!****************************************************************************
Subroutine check_inner_face(b_p, cp, check)

implicit none

real, intent(in) :: b_p(3,4), cp(3)
logical, intent(inout) :: check

integer :: i, i_r, pl_check(3), this(3,2)
real :: sum_area, ori_area
real :: dis(3), s
!----------------------------------------------------------------------------
this(:,1) = (/ 1, 2, 3 /)
this(:,2) = (/ 3, 4, 1 /)
ori_area = 0.0
do i = 1, 2
    call calc_length_3D(b_p(:,this(1,i)), b_p(:,this(2,i)), dis(1))
    call calc_length_3D(b_p(:,this(2,i)), b_p(:,this(3,i)), dis(2))
    call calc_length_3D(b_p(:,this(3,i)), b_p(:,this(1,i)), dis(3))
    s = (dis(1)+dis(2)+dis(3))/2.
    if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
        ori_area = ori_area + sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
    endif
enddo

sum_area = 0.
do i = 1, 4
    i_r = i+1
    if ( i == 4 ) i_r = 1
    call calc_length_3D(b_p(:,i), b_p(:,i_r), dis(1))
    call calc_length_3D(b_p(:,i), cp, dis(2))
    call calc_length_3D(b_p(:,i_r), cp, dis(3))
    s = (dis(1)+dis(2)+dis(3))/2.
    if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
        sum_area = sum_area + sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
    endif
enddo

if ( abs(sum_area-ori_area) < ori_area*0.01 ) then
    check = .TRUE.
else
    check = .FALSE.
endif

end subroutine check_inner_face
!****************************************************************************
Subroutine check_inner_face2(b_p, cp, check)

implicit none

real, intent(in) :: b_p(3,3), cp(3)
logical, intent(inout) :: check

integer :: i, j, j_r, pl_check(3), this(3,2)
real :: sum_area, ori_area
real :: dis(3), s
!----------------------------------------------------------------------------
ori_area = 0.0
call calc_length_3D(b_p(:,1), b_p(:,2), dis(1))
call calc_length_3D(b_p(:,2), b_p(:,3), dis(2))
call calc_length_3D(b_p(:,3), b_p(:,1), dis(3))
s = (dis(1)+dis(2)+dis(3))/2.
if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
    ori_area = ori_area + sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
endif

sum_area = 0.
do j = 1, 3
    j_r = j+1
    if ( j == 3 ) j_r = 1
    call calc_length_3D(b_p(:,j), b_p(:,j_r), dis(1))
    call calc_length_3D(b_p(:,j), cp, dis(2))
    call calc_length_3D(b_p(:,j_r), cp, dis(3))
    s = (dis(1)+dis(2)+dis(3))/2.
    if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
        sum_area = sum_area + sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
    endif
enddo

if ( abs(sum_area-ori_area) < ori_area*0.001 ) then
    check = .TRUE.
else
    check = .FALSE.
endif

end subroutine check_inner_face2
!==============================================================================
subroutine make_eci(ele, eci)

implicit none

integer,intent(in) :: ele(8)
integer,intent(inout) :: eci(3,8)
!---------------------------------------------------------------------------------
eci(:,1) = (/ ele(4), ele(2), ele(5) /)
eci(:,2) = (/ ele(1), ele(3), ele(6) /)
eci(:,3) = (/ ele(2), ele(4), ele(7) /)
eci(:,4) = (/ ele(3), ele(1), ele(8) /)
eci(:,5) = (/ ele(8), ele(6), ele(1) /)
eci(:,6) = (/ ele(5), ele(7), ele(2) /)
eci(:,7) = (/ ele(6), ele(8), ele(3) /)
eci(:,8) = (/ ele(7), ele(5), ele(4) /)

end subroutine make_eci
!****************************************************************************
!
!  PROGRAM: angle_cos
!
!  PURPOSE: Calculate the angle between three points by using the cos method
!
!****************************************************************************
Subroutine calc_angle_cos(a, b, c, angle)

implicit none

Real, intent(in) :: a(3), b(3), c(3)
Real, intent(out) :: angle

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

real :: dis(3)

call calc_length_3D(a, b, dis(1))
call calc_length_3D(b, c, dis(2))
call calc_length_3D(a, c, dis(3))

angle = 0.
if ( dis(1) == 0. .OR. dis(2) == 0. ) then
	angle = 0.
else
	angle = (dis(1)**2 + dis(2)**2 - dis(3)**2) / (2*dis(1)*dis(2))
    if ( angle < -1. ) then
        angle = -1.
    elseif ( angle > 1. ) then
        angle = 1.
    endif
    
	angle = acos(angle)
	angle = angle*(180/pi)
endif

end subroutine
! 세점의 각도 계산 종료
!****************************************************************************
!
!  PROGRAM: descending_sort_int
!
!  PURPOSE: 정수의 내림차순 정렬
!
!****************************************************************************
subroutine descending_sort_int(na, a)

implicit none

integer :: na, a(na)

integer :: i, j, temp
!----------------------------------------------------------------------------
do i = 1, na
    do j = i+1, na
        if ( a(i) > a(j) ) then
            temp = a(i)
            a(i) = a(j)
            a(j) = temp
        endif
    enddo
enddo

end subroutine descending_sort_int


subroutine find_nci(bfn, bound_face, sp, nci_num, nci)

implicit none

integer :: bfn, bound_face(4,bfn), sp, nci_num, nci(8)
integer :: i, j, k, q, this(2)
logical :: check

nci_num = 0
do i = 1, bfn
    do j = 1, 4
        if (bound_face(j,i) == sp) then
            if (j == 1) then
                this = (/ bound_face(4,i), bound_face(2,i) /)
            elseif (j == 4) then
                this = (/ bound_face(3,i), bound_face(1,i) /)
            else
                this = (/ bound_face(j-1,i), bound_face(j+1,i) /)
            endif
            do k = 1, 2
                check = .true.
                do q = 1, nci_num
                    if (this(k) == nci(q)) then
                        check = .false.
                        exit
                    endif
                enddo
                if (check) then
                    nci_num = nci_num + 1
                    nci(nci_num) = this(k)
                endif
            enddo
        endif
    enddo
    if (nci_num == 8) exit
enddo

end subroutine find_nci


!****************************************************************************
!
!  PROGRAM : rearrange_ele
!
!  PURPOSE : 표면의 삼각형 element의 외적 방향 재설정
!
!****************************************************************************
Subroutine rearrange_elem(prop_bfn, prop_face)

implicit none

integer, intent(in) :: prop_bfn
integer, intent(inout) :: prop_face(4,prop_bfn)

integer :: i, j, k, q, line(2), count, check_ele(prop_bfn), qele(4), this(4)
integer :: qne, pre_qne, temp_qne, temp_ele(prop_bfn), ce
real :: p(3,3), pl(4)

check_ele = 1
!line = prop_face(1:2,1)
line = (/ prop_face(2,1), prop_face(1,1) /)
find_first_ele: do i = 1, prop_bfn
    call find_count(2, line, 4, prop_face(:,i), count)
    if (count == 2) then
        check_ele(i) = 0
        qne = 1
        temp_ele(1) = i
        qele = prop_face(:,i)
        prop_face(1:2,i) = line
        do j = 1, 4
            if (line(2) == qele(j)) then
                if (j == 1) then
                    this = (/ 1, 2, 3, 4 /)
                elseif (j == 2) then
                    this = (/ 2, 3, 4, 1 /)
                elseif (j == 3) then
                    this = (/ 3, 4, 1, 2 /)
                elseif (j == 4) then
                    this = (/ 4, 1, 2, 3 /)
                endif
                if (qele(this(2)) == line(1)) then
                    prop_face(3,i) = qele(this(4))
                else
                    prop_face(3,i) = qele(this(2))
                endif
                prop_face(4,i) = qele(this(3))
                exit find_first_ele
            endif
        enddo
    endif
enddo find_first_ele
if (qne == 0) stop 'Error(rearrange_ele): Cannot find element matching the first plane'

temp_qne = qne
pre_qne = 1
whole_loop: do
    do i = pre_qne, temp_qne
        ce = temp_ele(i)
        do j = 1, 4
            if (j == 4) then
                line = (/ prop_face(1,ce), prop_face(4,ce) /)
            else
                line = (/ prop_face(j+1,ce), prop_face(j,ce) /)
            endif
            do k = 1, prop_bfn
                if (check_ele(k) == 1) then
                    call find_count(2, line, 4, prop_face(:,k), count)
                    if (count == 2) then
                        check_ele(k) = 0
                        qne = qne + 1
                        temp_ele(qne) = k
                        qele = prop_face(:,k)
                        prop_face(1:2,k) = line
                        do q = 1, 4
                            if (line(2) == qele(q)) then
                                if (q == 1) then
                                    this = (/ 1, 2, 3, 4 /)
                                elseif (q == 2) then
                                    this = (/ 2, 3, 4, 1 /)
                                elseif (q == 3) then
                                    this = (/ 3, 4, 1, 2 /)
                                elseif (q == 4) then
                                    this = (/ 4, 1, 2, 3 /)
                                endif
                                if (qele(this(2)) == line(1)) then
                                    prop_face(3,k) = qele(this(4))
                                else
                                    prop_face(3,k) = qele(this(2))
                                endif
                                prop_face(4,k) = qele(this(3))
                                exit
                            endif
                        enddo
                        exit
                    endif
                endif
            enddo
        enddo
    enddo
    if (qne == pre_qne-1) exit
    pre_qne = temp_qne + 1
    temp_qne = qne
    !write (*,*) pre_tne, temp_tne, tne
enddo whole_loop

end subroutine rearrange_elem

!****************************************************************************
!
!  PROGRAM : rearrange_bound_face
!
!  PURPOSE : 표면의 사각형 face의 외적 방향 재설정
!
!****************************************************************************
Subroutine rearrange_bound_face(prop_bfn, prop_face)

implicit none

integer, intent(in) :: prop_bfn
integer, intent(inout) :: prop_face(4,prop_bfn)

integer :: i, j, k, q, tnfn, count, cf, num, temp_tnfn, ro_dir
integer :: tnf(prop_bfn), temp_tnf(prop_bfn), line(2), temp(4), lv(4)
real :: p(3,3), pl(4)
logical :: check, check_face(prop_bfn)

check_face = .true.
!do i = 1, prop_bfn
!    lv = (1,bound%bound_face(:,i))
!    check = .true.
!    do j = 1, 4
!        if (lv(j) /= 1) then
!            check = .false.
!            exit
!        endif
!    enddo
!    if (check) then
!        tnfn = 1
!        tnf(tnfn) = i
!        check_face(i) = .false.
!        exit
!    endif
!enddo

tnfn = 1
tnf(tnfn) = 1
check_face(1) = .false.
whole_loop: do
    temp_tnfn = 0
    do i = 1, tnfn
        cf = tnf(i)
        do j = 1, 4
            if (j == 4) then
                line = (/ prop_face(1,cf), prop_face(j,cf) /)
            else
                line = (/ prop_face(j+1,cf), prop_face(j,cf) /)
            endif
            do k = 1, prop_bfn
                if (check_face(k)) then
                    temp = prop_face(:,k)
                    call find_count(2, line, 4, temp, count)
                    if (count == 2) then
                        temp_tnfn = temp_tnfn + 1
                        temp_tnf(temp_tnfn) = k
                        ro_dir = 1
                        do q = 1, 4
                            if (prop_face(q,k) == line(1)) then
                                if (q == 4) then
                                    if (prop_face(1,k) /= line(2)) then
                                        ro_dir = 2
                                    endif
                                else
                                    if (prop_face(q+1,k) /= line(2)) then
                                        ro_dir = 2
                                    endif
                                endif
                                exit
                            endif
                        enddo
                        if (ro_dir == 2) then
                            prop_face(:,k) = (/ temp(1), temp(4), temp(3), temp(2) /)
                            !call calc_plane(p(:,1), p(:,2), p(:,3), pl)
                            !write (*,*) temp
                            !write (*,*) bound%bound_face(:,k)
                            !write (*,'(4f13.8)') pl
                            !write (*,*) 
                        endif
                        check_face(k) = .false.
                    endif
                endif
            enddo
        enddo
    enddo
    if (temp_tnfn == 0) exit
    tnfn = temp_tnfn
    tnf(1:tnfn) = temp_tnf(1:tnfn)
enddo whole_loop


end subroutine rearrange_bound_face

subroutine check_bound_face(nn, coord, bfn, bound_face)

implicit none

integer, intent(in) :: nn, bfn, bound_face(4,bfn)
real, intent(in) :: coord(3,nn)

integer :: i, j
real :: pl(3,nn), plane(4), p(3,3)

pl = 0.0
do i = 1, bfn
    p = coord(:,bound_face(1:3,i))
    call calc_plane(p(:,1), p(:,2), p(:,3), plane)
    do j = 1, 4
        pl(:,bound_face(j,i)) = pl(:,bound_face(j,i)) + plane(1:3)
    enddo
enddo

open (Unit=20, File='bound_face_check.plt', STATUS='replace', ACTION='write')
Write (20,'(A)') 'TITLE="Check bound face"'
Write (20,*) 'VARIABLES="x", "y", "z", "pl_x", "pl_y", "pl_z"'
Write (20,'(A,I5,A,I5, A)') 'ZONE N =', nn , ', E =', bfn, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
do i = 1, nn
    write (20,*) coord(:,i), pl(:,i)
enddo
do i = 1, bfn
    write (20,*) bound_face(:,i)
enddo
close(20)

end subroutine check_bound_face



subroutine inout_check_cubic(cp, cubic_coord, check, diff_volume, tol)

implicit none

real, intent(in) :: cp(3), cubic_coord(3,8), tol
real, intent(inout) :: diff_volume
logical, intent(inout) :: check

integer :: i, j
integer :: this(4,5), face_nn(3,12)
real :: ea, cubic_vol(2), coord(3,4), volume, center_p(3), dis, max_dis

center_p = 0.0
do i = 1, 8
    center_p = center_p + cubic_coord(:,i)*0.125
enddo
max_dis = 0.0
do i = 1, 8
    call calc_length_3D(cubic_coord(:,i), center_p, dis)
    if (dis > max_dis) max_dis = dis
enddo
call calc_length_3D(cp, center_p, dis)
!write (*,*) dis, max_dis
if (dis <= max_dis*1.01) then
    this(:,1) = (/ 1, 2, 3, 6 /)
    this(:,2) = (/ 1, 3, 8, 6 /)
    this(:,3) = (/ 1, 3, 4, 8 /)
    this(:,4) = (/ 1, 6, 8, 5 /)
    this(:,5) = (/ 3, 8, 6, 7 /)
    cubic_vol = 0.0

    do i = 1, 5
        coord = cubic_coord(:,this(:,i))
        !vec(:,1) = coord(:,1)-coord(:,4)
        !vec(:,2) = coord(:,2)-coord(:,4)
        !vec(:,3) = coord(:,3)-coord(:,4)
        !
        !vec(1,4) = vec(2,2)*vec(3,3)-vec(3,2)*vec(2,3)
        !vec(2,4) = -(vec(1,2)*vec(3,3)-vec(3,2)*vec(1,3))
        !vec(3,4) = vec(1,2)*vec(2,3)-vec(2,2)*vec(1,3)
        !
        !cubic_vol(1) = cubic_vol(1) + abs(vec(1,1)*vec(1,4)+vec(2,1)*vec(2,4)+vec(3,1)*vec(3,4))/6.0
        call calc_tet_volume(coord, volume)
        cubic_vol(1) = cubic_vol(1) + volume
    enddo

    face_nn(:,1) = (/ 1, 2, 3 /)
    face_nn(:,2) = (/ 3, 4, 1 /)
    face_nn(:,3) = (/ 2, 1, 5 /)
    face_nn(:,4) = (/ 5, 6, 2 /)
    face_nn(:,5) = (/ 3, 2, 6 /)
    face_nn(:,6) = (/ 6, 7, 3 /)
    face_nn(:,7) = (/ 3, 7, 8 /)
    face_nn(:,8) = (/ 8, 4, 3 /)
    face_nn(:,9) = (/ 4, 8, 5 /)
    face_nn(:,10) = (/ 5, 1, 4 /)
    face_nn(:,11) = (/ 5, 8, 7 /)
    face_nn(:,12) = (/ 7, 6, 5 /)

    coord(:,4) = cp
    do i = 1, 12
        do j = 1, 3
            coord(:,j) = cubic_coord(:,face_nn(j,i))
        enddo
        !vec(:,1) = coord(:,1)-coord(:,4)
        !vec(:,2) = coord(:,2)-coord(:,4)
        !vec(:,3) = coord(:,3)-coord(:,4)
        !
        !vec(1,4) = vec(2,2)*vec(3,3)-vec(3,2)*vec(2,3)
        !vec(2,4) = -(vec(1,2)*vec(3,3)-vec(3,2)*vec(1,3))
        !vec(3,4) = vec(1,2)*vec(2,3)-vec(2,2)*vec(1,3)
        !
        !cubic_vol(2) = cubic_vol(2) + abs(vec(1,1)*vec(1,4)+vec(2,1)*vec(2,4)+vec(3,1)*vec(3,4))/6.0
        call calc_tet_volume(coord, volume)
        cubic_vol(2) = cubic_vol(2) + volume
    enddo

    diff_volume = abs(cubic_vol(1)-cubic_vol(2))
    ea = abs((cubic_vol(1)-cubic_vol(2))/cubic_vol(1))
    check = .true.
    if (ea > tol) check = .false.
else
    diff_volume = 10e8
    check = .false.
endif
    
!write (*,*) 'cubic_volume:', cubic_vol, check

end subroutine inout_check_cubic

subroutine inout_check_cubic2(ce, cp, cubic_coord, check, tol)

implicit none

integer, intent(in) :: ce
real, intent(in) :: cp(3), cubic_coord(3,8), tol
logical, intent(inout) :: check

integer :: i, j, k
integer :: this(4,5)
real :: det_mat(4,4), temp_mat(4,4), det(5), val
logical :: degen_elem

this(:,1) = (/ 1, 2, 3, 6 /)
this(:,2) = (/ 1, 6, 3, 8 /)
this(:,3) = (/ 1, 3, 4, 8 /)
this(:,4) = (/ 1, 6, 8, 5 /)
this(:,5) = (/ 3, 8, 6, 7 /)

check = .FALSE.
do i = 1, 5
    det_mat = 1.0
    do j = 1, 4
        det_mat(2:4,j) = cubic_coord(:,this(j,i))
    enddo
    call calc_determinat_4(det_mat, det(1))
    !det(1) = DMGT(eps,4,det_mat)
    val = det(1)
    if (abs(det(1)) <= tol*0.0001) then
        !write (*,*) ce, ' element has the degenerated tetrahedron '
        !stop 
        degen_elem = .FALSE.
    else
        degen_elem = .TRUE.
    endif
    if (degen_elem) then
        do j = 1, 4
            temp_mat = det_mat
            temp_mat(2:4,j) = cp
            call calc_determinat_4(temp_mat, det(j+1))
            val = val - det(j+1)
        enddo
        if (val/det(1) < 0.001) then
            check = .TRUE.
            do j = 2, 5
                if (det(j)/det(1) < -0.001) then
                    check = .FALSE.
                    exit
                endif
            enddo
        endif
    endif
    if (check) then
        exit
    endif
enddo

end subroutine inout_check_cubic2


subroutine calc_tet_volume(coord, volume)

implicit none

real, intent(in) :: coord(3,4)
real, intent(out) :: volume

real :: vec(3,4)

vec(:,1) = coord(:,1)-coord(:,4)
vec(:,2) = coord(:,2)-coord(:,4)
vec(:,3) = coord(:,3)-coord(:,4)
    
vec(1,4) = vec(2,2)*vec(3,3)-vec(3,2)*vec(2,3)
vec(2,4) = -(vec(1,2)*vec(3,3)-vec(3,2)*vec(1,3))
vec(3,4) = vec(1,2)*vec(2,3)-vec(2,2)*vec(1,3)
    
volume = abs(vec(1,1)*vec(1,4)+vec(2,1)*vec(2,4)+vec(3,1)*vec(3,4))/6.0

end subroutine calc_tet_volume

subroutine get_skin_surface_info1(rn, bn, bfn)

use surface_domain
implicit none

integer, intent(in) :: rn
integer, intent(inout) :: bn, bfn

bn = mesh_set(rn)%bn
bfn = mesh_set(rn)%bfn
    
end subroutine get_skin_surface_info1

subroutine get_skin_surface_info2(rn, bn, bound_node, bfn, bound_face)

use surface_domain
implicit none

integer, intent(in) :: rn, bn, bfn
integer, intent(inout) :: bound_node(bn), bound_face(4,bfn)

!integer :: i, nn, conn(mesh_set(rn)%nn)

!bound_node = mesh_set(rn)%bound_node
bound_node = mesh_set(rn)%surf2ori(mesh_set(rn)%bound_node)
!conn = 0
!do i = 1, bn
!    conn(mesh_set(rn)%bound_node(i)) = i
!enddo
!do i = 1, bfn
!    bound_face(:,i) = conn(mesh_set(rn)%bound_face(:,i))
!enddo
bound_face = mesh_set(rn)%bound_face

end subroutine get_skin_surface_info2

subroutine get_skin_surface_info3(rn, bn, bound_point, bound_disp, bound_corner, flag)

use surface_domain
use system_domain

integer, intent(in) :: rn, bn, flag
integer, intent(inout) :: bound_corner(bn)
real, intent(inout) :: bound_point(3,bn), bound_disp(3,bn)

integer :: i, cp

if (flag == -1 .or. flag == -3) then
    do i = 1, bn
        !cp = mesh_set(rn)%bound_node(i)
        !bound_point(:,i) = region(rn)%node_coord(:,cp)
        cp = mesh_set(rn)%surf2ori(mesh_set(rn)%bound_node(i))
        bound_point(:,i) = region(1)%node_coord(:,cp)
        if (rn == 1) bound_corner(i) = mesh_set(rn)%bound_lv(i)
    enddo
endif
if (flag == -1 .or. flag == -2) then
    write (*,*) 'allocated(sdomain(1)%u_cartesian):', allocated(sdomain(1)%u_cartesian)
    if (allocated(sdomain(1)%u_cartesian)) then
        do i = 1, bn
            !cp = mesh_set(rn)%bound_node(i)
            cp = mesh_set(rn)%surf2ori(mesh_set(rn)%bound_node(i))
            !bound_disp(:,i) = history(rn)%data(region(rn)%node(4:6,cp),1)
            !bound_disp(:,i) = sdomain(1)%u_cartesian(sdomain(1)%dof_index(region(1)%node(4:6,cp)))
            bound_disp(:,i) = sdomain(1)%u_cartesian(region(1)%node(4:6,cp))
        enddo
    else
        do i = 1, bn
            bound_disp(:,i) = 0.0
        enddo
    endif
endif

end subroutine get_skin_surface_info3

subroutine get_skin_surface_info4(rn, ng, bn, prop_patch)

use surface_domain
implicit none

integer, intent(in) :: rn, bn, ng
integer, intent(inout) :: prop_patch(bn)

integer :: i, sp, lp

do i = 1, ng
    sp = mesh_set(rn)%sub_num(i,3)
    lp = mesh_set(rn)%sub_num(i+1,3) - 1
    prop_patch(sp:lp) = i   
enddo

end subroutine get_skin_surface_info4

subroutine calc_cubic_volume(cubic_coord, volume)

implicit none

real, intent(in) :: cubic_coord(3,8)
real, intent(inout) :: volume

integer :: i, this(4,5)
real :: tet_volume, coord(3,4)

this(:,1) = (/ 1, 2, 3, 6 /)
this(:,2) = (/ 1, 3, 8, 6 /)
this(:,3) = (/ 1, 3, 4, 8 /)
this(:,4) = (/ 1, 6, 8, 5 /)
this(:,5) = (/ 3, 8, 6, 7 /)

volume = 0.0
do i = 1, 5
    coord = cubic_coord(:,this(:,i))
    call calc_tet_volume(coord, tet_volume)
    volume = volume + tet_volume
enddo

end subroutine calc_cubic_volume


subroutine quad2tri(current_step, rn, ng)

use surface_domain
implicit none

integer, intent(in) :: current_step, rn, ng

integer :: i, num, status, nn, ne, sp(2), lp(2)
real    :: p(3,4), dis(2)
character(len=45) :: filename
integer, allocatable :: elem(:,:)
real, allocatable :: coord(:,:)

!nn = mesh_set(ng)%bn
!ne = mesh_set(ng)%bfn
sp = mesh_set(rn)%sub_num(ng,3:4)
lp = mesh_set(rn)%sub_num(ng+1,3:4) - 1
nn = lp(1) - sp(1) + 1
ne = lp(2) - sp(2) + 1
allocate (coord(3,nn), elem(4,ne))
num = 0
do i = sp(1), lp(1)
    num = num + 1
    coord(:,num) = mesh_set(rn)%node_coord(:,mesh_set(rn)%bound_node(i))
enddo
elem = mesh_set(rn)%bound_face(:,sp(2):lp(2))
if (ng /= 1) elem = elem - (mesh_set(rn)%sub_num(ng,3)-1)

surface_set(ng)%tri_nn = nn
surface_set(ng)%tri_ne = ne*2
allocate (surface_set(ng)%tri_coord(3,nn))
allocate (surface_set(ng)%tri_elem(3,ne*2))

filename = './output/solid/remesh/remesh3D_00_0001.plt'
write (filename(32:33), '(I2.2)' ) ng
write (filename(35:38), '(I4.4)' ) current_step
open (Unit=2001, File=filename, ACTION='write', IOSTAT=status)
Write (2001,'(A)') 'TITLE="3D surface check!!!"'
Write (2001,*) 'VARIABLES="x", "y", "z"'
Write (2001,'(A,I6,A,I6, A)') 'ZONE N=', nn , ', E=', ne*2, ', DATAPACKING = POINT, ZONETYPE = FETRIANGLE'        
do i = 1, nn
    write (2001,'(3ES15.7)') coord(:,i)
enddo

surface_set(ng)%tri_coord = coord
do i = 1, ne
    num = (i-1)*2
    p = coord(:,elem(:,i))
    call calc_length_3D(p(:,1), p(:,3), dis(1))
    call calc_length_3D(p(:,2), p(:,4), dis(2))
    if (dis(1) < dis(2)) then
        surface_set(ng)%tri_elem(:,num+1) = (/ elem(1,i), elem(2,i), elem(3,i) /)
        surface_set(ng)%tri_elem(:,num+2) = (/ elem(3,i), elem(4,i), elem(1,i) /)
    else
        surface_set(ng)%tri_elem(:,num+1) = (/ elem(1,i), elem(2,i), elem(4,i) /)
        surface_set(ng)%tri_elem(:,num+2) = (/ elem(2,i), elem(3,i), elem(4,i) /)
    endif
    write (2001,'(3(I7,2X))') surface_set(ng)%tri_elem(:,num+1)
    write (2001,'(3(I7,2X))') surface_set(ng)%tri_elem(:,num+2)
enddo
close(2001)

end subroutine quad2tri


!****************************************************************************
!
!  PROGRAM : Cross
!
!  PURPOSE : 벡터의 외적
!
!****************************************************************************
Subroutine calc_vec_cross(a, b, c, ijk)

implicit none

Real, intent(in) :: a(3), b(3), c(3)
Real, intent(inout) :: ijk(3)

Real :: vec(3,2)
integer :: i, j

!----------------------------------------------------------------------------

vec(:,1) = b - a
vec(:,2) = c - b
ijk(1) = vec(2,1)*vec(3,2) - vec(3,1)*vec(2,2)
ijk(2) = -(vec(1,1)*vec(3,2) - vec(3,1)*vec(1,2))
ijk(3) = vec(1,1)*vec(2,2) - vec(2,1)*vec(1,2)

end subroutine calc_vec_cross


!****************************************************************************
!
!  PROGRAM : line_inner_check
!
!  PURPOSE : 임의의 점이 선분안에 포함되어있는지 검사(길이검사)
!
!****************************************************************************
! 
subroutine line_inner_check(p, cp, tol, check)

implicit none

real :: p(3,2), cp(3), tol
logical :: check

real :: dis(3)

call calc_length_3D(p(:,1), p(:,2), dis(1))
call calc_length_3D(p(:,1), cp, dis(2))
call calc_length_3D(p(:,2), cp, dis(3))

if (abs(dis(1)-dis(2)-dis(3)) < tol) then
    check = .TRUE.
else
    check = .FALSE.
endif

end subroutine line_inner_check


subroutine calc_determinat_4(a, det)

implicit none

real, intent(in) :: a(4,4)
real, intent(inout) :: det

det = 0.0
det = det + a(1,1)*a(2,2)*a(3,3)*a(4,4) - a(1,1)*a(2,2)*a(4,3)*a(3,4)
det = det - a(1,1)*a(2,3)*a(3,2)*a(4,4) + a(1,1)*a(2,3)*a(4,2)*a(3,4)
det = det + a(1,1)*a(2,4)*a(3,2)*a(4,3) - a(1,1)*a(2,4)*a(4,2)*a(3,3)

det = det - a(1,2)*a(2,1)*a(3,3)*a(4,4) + a(1,2)*a(2,1)*a(4,3)*a(3,4)
det = det + a(1,2)*a(2,3)*a(3,1)*a(4,4) - a(1,2)*a(2,3)*a(4,1)*a(3,4)
det = det - a(1,2)*a(2,4)*a(3,1)*a(4,3) + a(1,2)*a(2,4)*a(4,1)*a(3,3)

det = det + a(1,3)*a(2,1)*a(3,2)*a(4,4) - a(1,3)*a(2,1)*a(4,2)*a(3,4)
det = det - a(1,3)*a(2,2)*a(3,1)*a(4,4) + a(1,3)*a(2,2)*a(4,1)*a(3,4)
det = det + a(1,3)*a(2,4)*a(3,1)*a(4,2) - a(1,3)*a(2,4)*a(4,1)*a(3,2)

det = det - a(1,4)*a(2,1)*a(3,2)*a(4,3) + a(1,4)*a(2,1)*a(4,2)*a(3,3)
det = det + a(1,4)*a(2,2)*a(3,1)*a(4,3) - a(1,4)*a(2,2)*a(4,1)*a(3,3)
det = det - a(1,4)*a(2,3)*a(3,1)*a(4,2) + a(1,4)*a(2,3)*a(4,1)*a(3,2)

end subroutine calc_determinat_4