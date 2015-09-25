subroutine set_ale_region(unit_num_inte, Dim)

use surface_domain
use remesh_domain
implicit none

integer, intent(in) :: unit_num_inte, Dim 

integer :: i, j, k, q, o, temp_j, sp, ng, num, lp, iter, n, count, pre_num, user_case_nn, num_group, status
integer :: csn(Dim,2), this(3), seq(2), face(4), nci1(3), nci2(3), line(2), sub_hexa(8), node_num(3)
integer :: case_bn, case_nn, case_ne, prop_bn, prop_nn, p2cn, p2_csn, case_bfn, match_bfn, temp_ln, flag
real :: tol_dis, min_dis, dis, p(Dim,4), np(Dim), tol(5)
logical :: check
character(len=5) :: keyword
integer, allocatable :: case_b_nodes(:), case_elem(:,:), case_b_face(:,:), user2case(:)
integer, allocatable :: prop_b_nodes(:), prop2case(:,:), match_b_face(:,:), prop_b_lv(:), p2c_nn(:)
integer, allocatable :: cs_nodes1(:,:,:), cs_nodes2(:,:,:)
integer, allocatable :: temp_val(:), near_elem_num(:), near_elem(:,:), check_elem1(:), check_elem2(:)
integer, allocatable :: temp_line(:,:)
real, allocatable :: case_coord(:,:), prop_coord(:,:), conn_off(:)

!ng = region(1)%num_groups
num_group = sub_mesh_region
allocate (move_set(num_group))

! add data for remeshing
rewind(unit_num_inte)
do 
    read(unit_num_inte,*,iostat=status) keyword
    if (keyword == '*end' .or. status == -1) then
        exit
    elseif (keyword == '*remp') then
        backspace(unit_num_inte)
        read(unit_num_inte,*) keyword, ng
        
        prop_nn = mesh_set(1)%nn
        sp = mesh_set(1)%sub_num(ng,3)
        lp = mesh_set(1)%sub_num(ng+1,3)
        prop_bn = lp-sp
        allocate (prop_coord(Dim,prop_nn), prop_b_nodes(prop_bn), prop_b_lv(prop_bn))
        prop_coord = mesh_set(1)%node_coord
        prop_b_nodes = mesh_set(1)%bound_node(sp:lp-1)
        prop_b_lv = mesh_set(1)%bound_lv(sp:lp-1)
        
        read(unit_num_inte,*) keyword, num
        remesh_add(ng)%ori_np2 = num
        allocate (remesh_add(ng)%ori_p2_coord(3,num), remesh_add(ng)%ori_p2(num))
        remesh_add(ng)%ori_p2 = 0
        do i = 1, num
            read(unit_num_inte,*) remesh_add(ng)%ori_p2_coord(:,i)
            min_dis = 10e8
            do j = 1, prop_bn
                np = prop_coord(:,prop_b_nodes(j))
                call calc_length_3D(np, remesh_add(ng)%ori_p2_coord(:,i), dis)
                if (dis < tol(ng)*10e-4) then
                    remesh_add(ng)%ori_p2(i) = j
                    exit
                elseif (min_dis > dis) then
                    temp_j = j
                    min_dis = dis
                endif
            enddo
            if (remesh_add(ng)%ori_p2(i) == 0) then
                remesh_add(ng)%ori_p2(i) = temp_j
            endif
        enddo
        !write (*,*) 'remesh_add(ng)%ori_p2:', remesh_add(ng)%ori_p2(1:num)
        read(unit_num_inte,*) keyword, num
        remesh_add(ng)%np2 = num
        if (num /= 0) then
            allocate (remesh_add(ng)%p2_coord(3,num), remesh_add(ng)%p2(2,num), remesh_add(ng)%p2_lv(num))
            allocate (remesh_add(ng)%p2_flag(num), remesh_add(ng)%p2_node(3,num))
            remesh_add(ng)%p2 = 0
            remesh_add(ng)%p2_node = 0
            do i = 1, num
                read(unit_num_inte,*) flag
                backspace(unit_num_inte)
                if (flag == 1) then
                    read(unit_num_inte,*) remesh_add(ng)%p2_flag(i), remesh_add(ng)%p2_coord(:,i), remesh_add(ng)%p2(2,i), remesh_add(ng)%p2_lv(i)
                elseif (flag == 2) then
                    read(unit_num_inte,*) remesh_add(ng)%p2_flag(i), node_num, remesh_add(ng)%p2(2,i), remesh_add(ng)%p2_lv(i)
                    remesh_add(ng)%p2_node(:,i) = node_num
                    do j = 1, 3
                        if (node_num(j) <= remesh_add(ng)%ori_np2) then
                            remesh_add(ng)%p2_coord(j,i) = remesh_add(ng)%ori_p2_coord(j,node_num(j))
                        else
                            node_num(j) = node_num(j) - remesh_add(ng)%ori_np2
                            remesh_add(ng)%p2_coord(j,i) = remesh_add(ng)%p2_coord(j,node_num(j))
                        endif
                    enddo
                endif
                if (remesh_add(ng)%p2_flag(i) == 1) then
                    min_dis = 10e8
                    do j = 1, prop_bn
                        if (prop_b_lv(j) == remesh_add(ng)%p2_lv(i)) then
                            np = prop_coord(:,prop_b_nodes(j))
                            call calc_length_3D(np, remesh_add(ng)%p2_coord(:,i), dis)
                            if (dis < tol(ng)*10e-4) then
                                remesh_add(ng)%p2(1,i) = j
                                exit
                            elseif (min_dis > dis) then
                                temp_j = j
                                min_dis = dis
                            endif
                        endif
                    enddo
                    if (remesh_add(ng)%p2(1,i) == 0) then
                        remesh_add(ng)%p2(1,i) = temp_j
                    endif
                endif
            enddo
            write (*,*) 'remesh_add(ng)%p2:', remesh_add(ng)%p2(1,1:num)
        endif
        
        !write (*,*) 
        !read(unit_num_inte,*) keyword, num
        !remesh_add(ng)%ln = num
        !if (num /= 0) then
        !    allocate (remesh_add(ng)%line(2,num))
        !    do i = 1, num
        !        read(unit_num_inte,*) remesh_add(ng)%line(:,i)
        !    enddo
        !endif
        read(unit_num_inte,*) keyword, num
        remesh_add(ng)%sub_hexa_num = num
        allocate (remesh_add(ng)%sub_hexa(8,num))
        do i = 1, num
            read(unit_num_inte,*) remesh_add(ng)%sub_hexa(:,i)
            write (*,'(9(I6,1X))') i, remesh_add(ng)%sub_hexa(:,i)
        enddo
        allocate (temp_line(2,num*12))
        temp_ln = 0
        call find_remesh_add_line(remesh_add(ng)%ori_np2, remesh_add(ng)%np2, &
                                  remesh_add(ng)%sub_hexa_num, remesh_add(ng)%sub_hexa, temp_ln, temp_line)
        remesh_add(ng)%ln = temp_ln
        allocate (remesh_add(ng)%line(2,temp_ln))
        remesh_add(ng)%line = temp_line(:,1:temp_ln)
        deallocate (temp_line)
        
        remesh_add(ng)%fix_line_num = 0
        read (unit_num_inte,*) keyword, num
        write (*,*) keyword, num
        if (num /= 0) then
            allocate (remesh_add(ng)%fix_line(4,num))
            remesh_add(ng)%fix_line_num = num
            do i = 1, num
                read(unit_num_inte,*) remesh_add(ng)%fix_line(:,i)
                remesh_add(ng)%fix_line(4,i) = remesh_add(ng)%fix_line(4,i) - 1
                write (*,*) 'fix_line:', remesh_add(ng)%fix_line(:,i)
            enddo
        endif
        
        read (unit_num_inte,*) keyword, num
        write (*,*) keyword, num
        if (num /= 0) then
            allocate (remesh_add(ng)%fix_coord(6,num), remesh_add(ng)%fix_p_axis(num))
            remesh_add(ng)%fix_point_num = num
            do i = 1, num
                read(unit_num_inte,*) remesh_add(ng)%fix_p_axis(i), this(1:2)
                remesh_add(ng)%fix_coord(1:3,i) = region(1)%node_coord(:,this(1))
                remesh_add(ng)%fix_coord(4:6,i) = region(1)%node_coord(:,this(2))
            enddo
        endif
        
        deallocate (prop_coord, prop_b_nodes, prop_b_lv)
    endif
enddo

if (divided_mesh_region == 2) then
    case_nn = mesh_set(2)%nn
    case_ne = mesh_set(2)%ne
    case_bn = mesh_set(2)%bn
    case_bfn = mesh_set(2)%bfn
    allocate(case_elem((Dim-1)*4,case_ne), case_coord(Dim,case_nn))
    allocate(case_b_nodes(case_bn), case_b_face((Dim-1)*2,case_bfn))
    allocate(conn_off(case_bfn))
    case_coord = mesh_set(2)%node_coord
    case_elem = mesh_set(2)%element
    case_b_nodes = mesh_set(2)%bound_node
    do i = 1, case_bfn
        case_b_face(:,i) = case_b_nodes(mesh_set(2)%bound_face(:,i))
    enddo
    user_case_nn = region(1)%num_nodes
    allocate(user2case(user_case_nn))
    user2case = 0
    do i = 1, case_nn
        user2case(mesh_set(2)%surf2ori(i)) = i
    enddo

    rewind(unit_num_inte)
    do 
        read(unit_num_inte,*) keyword
        if (keyword == '*p2cn') then
            backspace(unit_num_inte)
            exit
        endif
    enddo

    do i = 1, case_bfn
        tol_dis = 0.0
        do j = 1, 4
            if (j == 4) then
                call calc_length_3D(case_coord(:,case_b_face(j,i)), case_coord(:,case_b_face(1,i)), dis)
            else
                call calc_length_3D(case_coord(:,case_b_face(j,i)), case_coord(:,case_b_face(j+1,i)), dis)
            endif
            tol_dis = tol_dis + dis
        enddo
        conn_off(i) = tol_dis*0.25
    enddo

    ! set case ale region
    read(unit_num_inte,*) keyword, p2_csn
    read(unit_num_inte,*) keyword, csn(1:2,1)
    csn(1:2,2) = csn(1:2,1)
    write (*,*) 'csn:', csn(1:2,1)

    do ng = 1, num_group
        prop_nn = mesh_set(1)%nn
        sp = mesh_set(1)%sub_num(ng,3)
        lp = mesh_set(1)%sub_num(ng+1,3)
        prop_bn = lp-sp
        allocate (prop_coord(Dim,prop_nn), prop_b_nodes(prop_bn))
        prop_coord = mesh_set(1)%node_coord
        prop_b_nodes = mesh_set(1)%bound_node(sp:lp-1)

        allocate (cs_nodes1(csn(1,1),csn(2,1),case_bn/4), cs_nodes2(csn(1,1),csn(2,1),case_bn/4))
        allocate (temp_val(csn(2,1)))
        cs_nodes1 = 0
        cs_nodes2 = 0
        read(unit_num_inte,*) temp_val, cs_nodes1(1,1,2)
        !write (*,*) temp_val
        cs_nodes1(1,:,1) = temp_val
        read(unit_num_inte,*) temp_val, cs_nodes2(1,1,2)
        cs_nodes2(1,:,1) = temp_val
        do i = 1, csn(2,1)
            cs_nodes1(1,i,1) = user2case(cs_nodes1(1,i,1))
            cs_nodes2(1,i,1) = user2case(cs_nodes2(1,i,1))
        enddo
        cs_nodes1(1,1,2) = user2case(cs_nodes1(1,1,2))
        cs_nodes2(1,1,2) = user2case(cs_nodes2(1,1,2))
        deallocate (temp_val)
        !write (*,*) cs_nodes1(1,:,1)

        allocate (prop2case(p2_csn,prop_bn), p2c_nn(prop_bn))
        p2cn = 0
        p2c_nn = 0
        prop2case = 0
        ! find propellant surface nodes in matching surface
        do i = 1, prop_bn
            np = prop_coord(:,prop_b_nodes(i))
            do j = 1, case_bfn
                p = case_coord(:,case_b_face(:,j))
                call near_plane(p, np, conn_off(j)*0.01, check)
                if (check) then
                    p2cn = p2cn + 1
                    prop2case(1,p2cn) = prop_b_nodes(i)
                    p2c_nn(p2cn) = p2c_nn(p2cn) + 1
                    exit
                endif
            enddo
        enddo

        ! match propellant surface nodes to case surface node
        tol(ng) = minval(conn_off)
        do i = 1, p2cn
            do j = 1, case_bn
                call calc_length_3D(prop_coord(:,prop2case(1,i)), case_coord(:,case_b_nodes(j)), dis)
                if (dis < tol(ng)*10e-4) then
                    prop2case(2,i) = case_b_nodes(j)
                    p2c_nn(i) = p2c_nn(i) + 1
                    exit
                endif
            enddo
        enddo

        ! set matching surface(prop & case)
        allocate (match_b_face(4,case_bfn), temp_val(p2cn))
        temp_val = prop2case(2,1:p2cn)
        match_bfn = 0
        match_b_face = 0
        do i = 1, case_bfn
            face = case_b_face(:,i)
            call find_count(4, face, p2cn, temp_val, count)
            if (count == 4) then
                match_bfn = match_bfn + 1
                match_b_face(:,match_bfn) = face
            endif
        enddo
    
        open (Unit=113, File='check_ale_match_face.plt', STATUS='replace', ACTION='write')
        Write (113,'(A)') 'TITLE="Boundary face"'
        Write (113,'(A,A)') 'VARIABLES= "x", "y", "z"'
        Write (113,'(A, I5,A,I5, A)') 'ZONE N =', case_nn, ', E =', match_bfn,  ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
        do i = 1, case_nn
            write (113,*) case_coord(:,i)
        enddo
        do i = 1, match_bfn
            write (113,*) match_b_face(:,i)
        enddo
        close(113)
    
        ! find case inner node
        allocate (near_elem_num(case_nn), near_elem(20,case_nn))
        allocate (check_elem1(case_ne), check_elem2(case_ne))
        check_elem1 = 0
        check_elem2 = 0
        near_elem_num = 0
        near_elem = 0
        do i = 1, case_ne
            do j = 1, 8
                near_elem_num(case_elem(j,i)) = near_elem_num(case_elem(j,i)) + 1
                near_elem(near_elem_num(case_elem(j,i)), case_elem(j,i)) = i
            enddo
        enddo
        do i = 3, p2_csn
            check_elem2 = 0
            do j = 1, p2cn
                if (p2c_nn(j) == i-1) then
                    !if (temp_val(j) == 0) then
                    !    write (*,*) i, p2c_nn(j)
                    !    write (*,*) prop2case(1:i-1,j)
                    !endif
                    nci1 = 0
                    num = 0
                    do k = 1, near_elem_num(temp_val(j))
                        if (check_elem1(near_elem(k,temp_val(j))) == 0) then
                            if (nci1(1) == 0) then
                                call find_nci_elem_hexa(temp_val(j), case_elem(:,near_elem(k,temp_val(j))), nci1)
                                num = 3
                            else
                                call find_nci_elem_hexa(temp_val(j), case_elem(:,near_elem(k,temp_val(j))), nci2)
                                call find_same_number(num, nci1, 3, nci2, count)
                                num = count
                            endif
                            check_elem2(near_elem(k,temp_val(j))) = 1
                        endif
                    enddo
                    if (num == 1) then
                        call find_count(1, nci1(1), p2cn, temp_val, count)
                        if (count == 0) then
                            prop2case(i,j) = nci1(1)
                            p2c_nn(j) = p2c_nn(j) + 1
                        endif
                    elseif (num >= 2) then
                        do k = 1, num
                            call find_count(1, nci1(k), p2cn, temp_val, count)
                            if (count == 0) then
                                prop2case(i,j) = nci1(k)
                                p2c_nn(j) = p2c_nn(j) + 1
                                exit
                            endif
                        enddo
                        !if (count == 1) stop 'Cannot find connected node!!'
                    elseif (num == 0) then
                        !stop 'Cannot find connected node!!'
                    endif
                endif
            enddo
            if (i /= p2_csn) then
                temp_val = prop2case(i,1:p2cn)
                do j = 1, case_ne
                    if (check_elem2(j) == 1) check_elem1(j) = 1
                enddo
            endif
        enddo
        temp_val = prop2case(2,1:p2cn)

        num = case_bn/4
        call find_case_ale_region(csn(:,1), num, cs_nodes1, case_bfn, case_b_face, case_ne, case_elem, p2cn, temp_val)
        call find_case_ale_region(csn(:,2), num, cs_nodes2, case_bfn, case_b_face, case_ne, case_elem, p2cn, temp_val)

        deallocate (temp_val)
        deallocate (near_elem_num, near_elem)
        deallocate (check_elem1, check_elem2)
    
        ! save surface domain
        allocate (move_set(ng)%csn(Dim,2))
        allocate (move_set(ng)%cs_nodes1(csn(1,1),csn(2,1),csn(3,1)), move_set(ng)%cs_nodes2(csn(1,2),csn(2,2),csn(3,2)))
        allocate (move_set(ng)%prop2case(p2_csn,p2cn), move_set(ng)%p2c_nn(p2cn))
        allocate (move_set(ng)%match_b_face(4,match_bfn))
        move_set(ng)%csn = csn
        move_set(ng)%p2cn = p2cn
        move_set(ng)%p2_csn = p2_csn
        do i = 1, csn(1,1)
            do j = 1, csn(2,1)
                move_set(ng)%cs_nodes1(i,j,:) = cs_nodes1(i,j,1:csn(3,1))
                !if (i == 1) write (*,*) move_set(ng)%cs_nodes1(i,j,:)
            enddo
        enddo
        do i = 1, csn(1,2)
            do j = 1, csn(2,2)
                move_set(ng)%cs_nodes2(i,j,:) = cs_nodes2(i,j,1:csn(3,2))
            enddo
        enddo
        move_set(ng)%prop2case = prop2case(:,1:p2cn)
        move_set(ng)%p2c_nn = p2c_nn(1:p2cn)

        move_set(ng)%match_bfn = match_bfn
        move_set(ng)%match_b_face = match_b_face(:,1:match_bfn)
        deallocate (cs_nodes1, cs_nodes2, prop2case, p2c_nn)

        deallocate(prop_coord, prop_b_nodes)
        deallocate(match_b_face)
    enddo
    deallocate(case_elem, case_coord)
    deallocate(case_b_nodes, case_b_face)
    deallocate(conn_off)
    
    rewind(unit_num_inte)
    do 
        read(unit_num_inte,*,iostat=status) keyword
        if (keyword == '*end' .or. status == -1) then
            exit
        elseif (keyword == '*revo') then
            case_ne = region(1)%num_elements
            allocate (case_elem(8,case_ne))
            case_elem = region(1)%element(4:11,:)
            
            backspace(unit_num_inte)
            read(unit_num_inte,*) keyword, revol_region_num
            allocate(revol_set(revol_region_num))
            do i = 1, revol_region_num
                read(unit_num_inte,*) keyword, num, revol_set(i)%revol_info(1:3)
                allocate (revol_set(i)%revol_nodes(revol_set(i)%revol_info(1),revol_set(i)%revol_info(3)))
                revol_set(i)%revol_nodes = 0
                n = MOD(revol_set(i)%revol_info(1), 10)
                count = 1
		        do j = 1, revol_set(i)%revol_info(1)/10                                 ! each line: 10 input data
			        read(unit_num_inte,*) (revol_set(i)%revol_nodes(k,1), k=count, count+9)   ! limitation: field(50) -> 5 input lines
			        count = count + 10
                end do
                if (n > 0) then
		            read(unit_num_inte,*) (revol_set(i)%revol_nodes(k,1), k=count, count+n-1)  ! read last line if less than 10 data exist there
                endif
                call find_other_revol_nodes(i, case_ne, case_elem, revol_set(i)%revol_info, revol_set(i)%revol_nodes)
                do j = 1, revol_set(i)%revol_info(3)
                    revol_set(i)%revol_nodes(:,j) = user2case(revol_set(i)%revol_nodes(:,j))
                    !if (j == 1) write (*,*) revol_set(i)%revol_nodes(:,j)
                enddo
            enddo
            deallocate (case_elem)
        endif
    enddo
    deallocate (user2case)
endif

end subroutine set_ale_region
!==========================================================================================

subroutine find_opposite_node(face, seq, elem, num)

implicit none

integer, intent(in) :: face(4), seq, elem(8)
integer, intent(inout) :: num

integer :: i, j, k, eci(3,8)
logical :: check

call make_eci(elem, eci)
do i= 1, 8
    if (elem(i) == face(seq)) then
        do j = 1, 3
            check = .TRUE.
            do k = 1, 4
                if (face(k) == eci(j,i)) then
                    check = .FALSE.
                    exit
                endif
            enddo
            if (check) then
                num = eci(j,i)
                exit
            endif
        enddo
        exit
    endif
enddo

end subroutine find_opposite_node


subroutine find_case_ale_region(csn, temp_csn, cs_nodes, case_bfn, case_b_face, case_ne, case_elem, p2cn, prop2case)

implicit none

integer, intent(in) :: temp_csn, case_bfn, case_b_face(4,case_bfn), case_ne, case_elem(8,case_ne), p2cn, prop2case(p2cn)
integer, intent(inout) :: csn(3), cs_nodes(csn(1),csn(2),temp_csn)

integer :: pre_num, i, j, k, q, o, count, num, seq(2), this(3), face_seq(4,2), face(4)
logical :: check

csn(3) = 1
!write (*,*) cs_nodes(1,:,csn(3))
!write (*,*) 
pre_num = 0
find_case_ale_region1: do 
    csn(3) = csn(3) + 1
    ! find first node(cs_nodes(1,1,csn(3))
    if (csn(3) /= 2) then
        this = (/ cs_nodes(1,1,csn(3)-1), cs_nodes(1,2,csn(3)-1), pre_num /)
        do j = 1, case_bfn
            call find_count(2, this(1:2), 4, case_b_face(:,j), count)
            if (count == 2) then
                call find_count(3, this, 4, case_b_face(:,j), count)
                if (count == 2) then
                    do k = 1, 4
                        if (case_b_face(k,j) == this(1)) then
                            if (k == 1) then
                                seq = (/ 2, 4 /)
                            elseif (k == 4) then
                                seq = (/ 3, 1 /)
                            else
                                seq = (/ k-1, k+1 /)
                            endif
                            if (case_b_face(seq(1),j) == this(2)) then
                                cs_nodes(1,1,csn(3)) = case_b_face(seq(2),j)
                            else
                                cs_nodes(1,1,csn(3)) = case_b_face(seq(1),j)
                            endif
                            exit ! k
                        endif
                    enddo
                    exit ! j
                endif
            endif
        enddo
    endif

    ! find other nodes(cs_nodes(1,2:csn(2),csn(3))
    do i = 2, csn(2) 
        this = (/ cs_nodes(1,i-1,csn(3)-1), cs_nodes(1,i,csn(3)-1), cs_nodes(1,i-1,csn(3)) /)
        do j = 1, case_bfn
            call find_count(3, this, 4, case_b_face(:,j), count)
            if (count == 3) then
                do k = 1, 4
                    num = case_b_face(k,j)
                    do q = 1, 3
                        if (case_b_face(k,j) == this(q)) then
                            num = 0
                            exit  ! q
                        endif
                    enddo
                    if (num /= 0) then
                        cs_nodes(1,i,csn(3)) = num
                        exit ! k
                    endif
                enddo
                exit ! j
            endif
        enddo
    enddo
    pre_num = cs_nodes(1,1,csn(3)-1)
    !write (*,*) 'csn(3):', csn(3)
    !write (*,*) cs_nodes(1,:,csn(3))
    !write (*,*) 'pre_num:', pre_num
    call find_count(1, cs_nodes(1,csn(2),csn(3)), p2cn, prop2case, count)
    !if (cs_nodes(1,csn(2),csn(3)) == end_nodes(csn(2))) exit find_case_ale_region1
    if (count == 1) exit find_case_ale_region1
enddo find_case_ale_region1

do i = 2, csn(1)
    do j = 1, csn(2)-1
        face_seq(:,1) = (/ j, j+1, j, j+1 /)
        do k = 1, csn(3)-1
            face_seq(:,2) = (/ k, k, k+1, k+1 /)
            face = (/ cs_nodes(i-1,j,k), cs_nodes(i-1,j+1,k), cs_nodes(i-1,j,k+1), cs_nodes(i-1,j+1,k+1) /)
            !if (i == 4 .and. j == 1 .and. k == 1) write (*,*) 'face:', face
            do q = 1, case_ne
                call find_count(4, face, 8, case_elem(:,q), count)
                if (count == 4) then
                    check = .true.
                    if (i /= 2) then
                        do o = 1, 8
                            if (cs_nodes(i-2,j,k) == case_elem(o,q)) then
                                check = .false.
                                exit
                            endif
                        enddo
                    endif
                    !if (i == 4 .and. j == 1 .and. k == 1) then
                    !    write (*,*) 'Check:', check
                    !    write (*,'(A,8(I6,1X))') 'elem:', case_elem(:,q)
                    !endif
                    if (check) then
                        do o = 1, 4
                            if (cs_nodes(i,face_seq(o,1),face_seq(o,2)) == 0) then
                                call find_opposite_node(face, o, case_elem(:,q), num)
                                cs_nodes(i,face_seq(o,1),face_seq(o,2)) = num
                            endif
                        enddo
                        !write (*,*) 'cs_nodes(i):', cs_nodes(i,j,k), cs_nodes(i,j+1,k), cs_nodes(i,j,k+1), cs_nodes(i,j+1,k+1)
                        exit ! q
                    endif
                endif
            enddo ! q
        enddo ! k
    enddo ! j
enddo  ! i
!write (*,*) 'csn(3):', csn(3)
!do j = 1, csn(3)
!    do i = 1, 1 !csn(2)
!        write (*,*) cs_nodes(:,i,j)
!    enddo
!enddo
!write (*,*)

end subroutine find_case_ale_region

subroutine find_remesh_add_line(ori_np2, add_np2, sub_hexa_num, sub_hexa, temp_ln, temp_line)

implicit none

integer, intent(in) :: ori_np2, add_np2, sub_hexa_num, sub_hexa(8,sub_hexa_num)
integer, intent(inout) :: temp_ln, temp_line(2,sub_hexa_num*12)

integer :: i, j, k, np2, count, line(2), hexa_line(2,12), num, node_lv(2)
integer, allocatable :: p2_lv(:)
logical :: check

hexa_line(1,:) = (/ 1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 3, 4 /)
hexa_line(2,:) = (/ 2, 4, 6, 8, 3, 1, 7, 5, 5, 6, 7, 8 /)

temp_ln = 0
open(unit=35, file='hexa_line.dat')
do i = 1, sub_hexa_num
    do j = 1, 12
        line = sub_hexa(hexa_line(:,j),i)
        check = .TRUE.
        do k = 1, temp_ln
            call find_count(2, line, 2, temp_line(:,k), count)
            if (count == 2) then
                check = .FALSE.
                exit
            endif
        enddo
        if (check) then
            temp_ln = temp_ln + 1
            temp_line(:,temp_ln) = line
            write (35,*) temp_line(:,temp_ln)
        endif
    enddo
enddo
close(35)

end subroutine find_remesh_add_line