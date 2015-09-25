subroutine surface_smoothing(ng, flag)

use remesh_domain
    
implicit none

integer, intent(in) :: ng, flag
integer :: i, j, k, cbn, ten, root, num_nci, num
integer :: nele(8,Hexa(ng)%nn), nele_num(Hexa(ng)%nn), surf_ele(4), nci(8), this(4), temp_ele(4,8), near_node(2)
real :: t, u, v, quad_val, area, tol_area, dot_a, dot_b, sd, suv, ang
real :: p(3,4), su(3), sv(3), sc(3), q1(4), q2(4), dis(4), ori_coord(3,Hexa(ng)%nn), prop_line(3,3)
real :: alpha_ii(4) ! (1) = alpha21, (2) = alpha43, (3) = alpha31, (4) = alpha42
real :: alpha_i(4) ! (1) = -alpha21-alpha31, (2) = -alpha21-alpha42
                   ! (3) = -alpha31-alpha43, (4) = -alpha43-alpha42
real :: pl(3), temp_pl(4), mp(3), vec(3), move_vec(3), val, cp(3)
real :: DNRM2
logical :: check

!if (flag == 1) then
!    call get_boundary_node(1)
!    call get_boundary_face(1)
!endif

quad_val = sqrt(3.0)/6.0
q1 = (/ -1, 1, -1, 1 /)
q2 = (/ -1, -1, 1, 1 /)

ori_coord = Hexa(ng)%node(:,1:Hexa(ng)%nn)
do i = 1, Hexa(ng)%bn
    cbn = Hexa(ng)%bound_node(i)
    if (Hexa(ng)%node_lv(1,cbn) == 1 .or. Hexa(ng)%node_lv(1,cbn) == 2) then
        ten = 0
        do j = 1, Hexa(ng)%bfn
            do k = 1, 4
                if (Hexa(ng)%bound_face(k,j) == cbn) then
                    ten = ten + 1
                    nele(ten,i) = j
                    exit
                endif
            enddo
        enddo
        nele_num(i) = ten
    endif
enddo

do root = 1, 4
    do i = 1, Hexa(ng)%bn
        cbn = Hexa(ng)%bound_node(i)
        if (Hexa(ng)%node_lv(1,cbn) == 2) then
            ten = nele_num(i) 
            do j = 1, nele_num(i)
                temp_ele(:,j) = Hexa(ng)%bound_face(:,nele(j,i))
            enddo
            call find_nci(ten, temp_ele(:,1:ten), cbn, num_nci, nci)
            
            num = 0
            do j = 1, num_nci
                if (Hexa(ng)%node_lv(1,nci(j)) == 2 .or. Hexa(ng)%node_lv(1,nci(j)) == 3) then
                    num = num + 1
                    near_node(num) = nci(j)
                endif
                if (num == 2) exit
            enddo
            
            p(:,1) = Hexa(ng)%node(:,cbn)
            prop_line(:,1) = ori_coord(:,cbn)
            do j = 1, 2
                call calc_length_3D(p(:,1), Hexa(ng)%node(:,near_node(j)), dis(j))
                prop_line(:,j+1) = ori_coord(:,near_node(j))
            enddo
            call calc_angle_cos(prop_line(:,2), prop_line(:,1), prop_line(:,3), ang)
            if (ang > 160.0) then
                dis(3) = (dis(1)+dis(2))*0.5
                if (dis(1) > dis(2)) then
                    num = 1
                else
                    num = 2
                endif
                dis(3) = dis(num) - dis(3)
                vec = Hexa(ng)%node(:,near_node(num)) - p(:,1)
                !vec = ori_coord(:,near_node(num)) - p(:,1)
                val = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
                vec = vec/val
                Hexa(ng)%node(:,cbn) = Hexa(ng)%node(:,cbn) + vec*dis(3)
                p(:,1) = Hexa(ng)%node(:,cbn)
                do j = 1, 2
                    call calc_perpendicular_point(p(:,1), prop_line(:,1), prop_line(:,j+1), p(:,2), t, check)
                    if (check) then
                        Hexa(ng)%node(:,cbn) = p(:,2)
                        exit
                    endif
                enddo
            endif
        endif
    enddo
enddo

do root = 1, 4
    do i = 1, Hexa(ng)%bn
        cbn = Hexa(ng)%bound_node(i)
        !if (Hexa(ng)%node_lv(1,cbn) == 1 .or. Hexa(ng)%node_lv(1,cbn) == 2) then
        if (Hexa(ng)%node_lv(1,cbn) == 1) then
            cp = Hexa(ng)%node(:,cbn)
            ten = nele_num(i) 
            do j = 1, nele_num(i)
                temp_ele(:,j) = Hexa(ng)%bound_face(:,nele(j,i))
            enddo
            pl = 0.0
            tol_area = 0.0
            mp = 0.0
            do j = 1, ten
                surf_ele = temp_ele(:,j)
                p = Hexa(ng)%node(:,surf_ele)
                call calc_plane(p(:,1), p(:,2), p(:,3), temp_pl)
                pl = pl + temp_pl(1:3)
                area = 0.0
                do k = 1, 4
                    u = 0.5 + q1(k)*quad_val
                    v = 0.5 + q2(k)*quad_val
                    su = (1-v)*(p(:,2)-p(:,1))+v*(p(:,4)-p(:,3))
                    sv = (1-u)*(p(:,3)-p(:,1))+u*(p(:,4)-p(:,2))
                    sc(1) = su(2)*sv(3)-su(3)*sv(2)
                    sc(2) = su(3)*sv(1)-su(1)*sv(3)
                    sc(3) = su(1)*sv(2)-su(2)*sv(1)
                    sd = DOT_PRODUCT(su, sv)
                    dot_a = DOT_PRODUCT(su, su)
                    dot_b = DOT_PRODUCT(sv, sv)
                    suv = dot_a * dot_b - sd**2
                    area = area + suv
                enddo
                area = area * 0.25
                tol_area = tol_area + area
                !mp = mp + (p(:,1)+p(:,2)+p(:,3)+p(:,4))*0.25*area
                do k = 1, 4
                    if (surf_ele(k) == cbn) then
                        if (k == 1) then
                            mp = mp + (p(:,4)+p(:,2))*0.5*area
                        elseif (k == 4) then
                            mp = mp + (p(:,3)+p(:,1))*0.5*area
                        else
                            mp = mp + (p(:,k-1)+p(:,k+1))*0.5*area
                        endif
                        exit
                    endif
                enddo
            enddo
            pl = pl / ten
            mp = mp / tol_area
            !mp = mp / ten
            move_vec = mp-cp
            val = pl(1)*move_vec(1) + pl(2)*move_vec(2) + pl(3)*move_vec(3)
            move_vec = move_vec - val*pl
            Hexa(ng)%node(:,cbn) = cp + move_vec
        endif
    enddo
enddo

end subroutine surface_smoothing


subroutine inner_node_smoothing(ng)

use remesh_domain
implicit none
integer, intent(in) :: ng

integer :: i, j, k, q, inn, num, cn, count, root
integer :: nci(16), face_nn(4,6), this(3)
integer, allocatable :: inner_node(:), inn_nci(:,:), inn_nci_num(:)
real :: cp(3)
logical :: check

inn = Hexa(ng)%nn - Hexa(ng)%bn
allocate (inner_node(inn), inn_nci_num(inn), inn_nci(16,inn))
inn_nci_num = 0
inn_nci = 0
inn = 0
do i = 1, Hexa(ng)%nn
    check = .true.
    do j = 1, Hexa(ng)%bn
        if (i == Hexa(ng)%bound_node(j)) then
            check = .false.
            exit
        endif
    enddo
    if (check) then
        inn = inn + 1
        inner_node(inn) = i
    endif
enddo

do i = inn, 1, -1
    cn = inner_node(i)
    do j = 1, Hexa(ng)%ne
        do k = 1, 8
            if (Hexa(ng)%elem(k,j) == cn) then
                inn_nci_num(i) = inn_nci_num(i) + 1
                inn_nci(inn_nci_num(i),i) = j
                exit
            endif
        enddo
    enddo
enddo

do root = 1, 5
    do i = 1, inn
        cn = inner_node(i)
        nci = inn_nci(:,i)
        num = inn_nci_num(i)
        cp = 0.0
        count = 0
        do j = 1, num
            if (nci(j) > Hexa(ng)%ne) then
                write (*,*) nci(j), Hexa(ng)%ne
                stop
            endif
            call make_face_nn(8, 4, Hexa(ng)%elem(:,nci(j)), face_nn)
            do k = 1, 6
                check = .false.
                do q = 1, 4
                    if (face_nn(q,k) == cn) then
                        if (q == 1) then
                            this = (/ face_nn(2,k), face_nn(4,k), face_nn(3,k) /)
                        elseif (q == 2) then
                            this = (/ face_nn(3,k), face_nn(1,k), face_nn(4,k) /)
                        elseif (q == 3) then
                            this = (/ face_nn(4,k), face_nn(2,k), face_nn(1,k) /)
                        else
                            this = (/ face_nn(1,k), face_nn(3,k), face_nn(2,k) /)
                        endif
                        check = .true.
                        exit
                    endif
                enddo
                if (check) then
                    count = count + 2
                    cp = cp + Hexa(ng)%node(:,this(1)) + Hexa(ng)%node(:,this(2)) ! + Hexa(ng)%node(:,this(3))
                endif
            enddo
        enddo
        cp = cp / real(count)
        !if (cn == 20912) then
        !    write (*,'(I6,A,3F12.7)') cn, ' node:', Hexa(ng)%node(:,cn)
        !    write (*,'(I6,3F12.7)') count, cp
        !endif
        Hexa(ng)%node(:,cn) = cp
    enddo
enddo
deallocate (inner_node, inn_nci_num, inn_nci)

end subroutine inner_node_smoothing



subroutine adjust_fix_coord(ng)

use remesh_domain
!use surface_domain

implicit none
integer, intent(in) :: ng

integer :: i, j, k, cp, nci(3), axis, count, pre_num(1), temp_num, num
integer, allocatable :: fixed_point_group(:,:)
real :: min_dis, max_dis, dis, ang
real :: p(3,4), pre_p(3), yz_p(2,2)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

allocate (fixed_point_group(Hexa(ng)%nn,remesh_add(ng)%fix_point_num))
fixed_point_group = 0
do i = 1, remesh_add(ng)%fix_point_num
    axis = remesh_add(ng)%fix_p_axis(i)
    ! find first point
    min_dis = 10e8
    cp = 0
    p(:,1) = remesh_add(ng)%fix_coord(1:3,i)
    do j = 1, Hexa(ng)%nn
        if (Hexa(ng)%node_lv(1,j) == 2) then
            p(:,2) = Hexa(ng)%node(:,j)
            call calc_length_3D(p(:,1), p(:,2), dis)
            if (dis < min_dis) then
                min_dis = dis
                cp = j
            endif
        endif
    enddo
    
    p(:,3) = (/ 0.0, Hexa(ng)%node(2,cp), Hexa(ng)%node(3,cp) /)
    p(:,4) = (/ 0.0, 0.0, 0.0 /)
    pre_p = Hexa(ng)%node(:,cp)
    Hexa(ng)%node(:,cp) = p(:,1)
    Hexa(ng)%node_lv(1,cp) = 3
    fixed_point_group(1,i) = cp
    
    pre_num(1) = 0
    num = 1
    !do j = 2, csn
    do
        num = num + 1
        find_elem: do j = 1, Hexa(ng)%ne
            do k = 1, 8
                if (Hexa(ng)%elem(k,j) == cp) then
                    call find_nci_elem_hexa(cp, Hexa(ng)%elem(:,j), nci)
                    call find_count(1, pre_num(1), 3, nci, count)
                    if (count == 0) exit find_elem
                endif
            enddo
        enddo find_elem
        pre_num(1) = cp
        
        yz_p(:,1) = pre_p(2:3)
        max_dis = 0.0
        temp_num = 0
        do k = 1, 3
            yz_p(:,2) = Hexa(ng)%node(2:3,nci(k))
            !call calc_len(yz_p(:,1), yz_p(:,2), dis)  
            dis = abs(yz_p(2,1) - yz_p(2,2))
            if (dis > max_dis) then
                max_dis = dis
                cp = nci(k)
            endif
        enddo
        p(1,2) = 0.0
        p(2:3,2) = Hexa(ng)%node(2:3,cp)
        call calc_angle_cos(p(:,3), p(:,4), p(:,2), ang)

        Hexa(ng)%node(1,cp) = p(1,1)
        Hexa(ng)%node(2,cp) = p(2,3)*cos(ang*pi/180.) - p(3,3)*sin(ang*pi/180.)
        Hexa(ng)%node(3,cp) = p(2,3)*sin(ang*pi/180.) + p(3,3)*cos(ang*pi/180.)
        fixed_point_group(num,i) = cp
        
        if (Hexa(ng)%node_lv(1,cp) == 2) then
            Hexa(ng)%node(:,cp) = remesh_add(ng)%fix_coord(4:6,i)
            Hexa(ng)%node_lv(1,cp) = 3
            exit
        endif
        Hexa(ng)%node_lv(1,cp) = 3
    enddo
enddo
Hexa(ng)%fixed_point_group_num = num
allocate (Hexa(ng)%fixed_point_group(num,remesh_add(ng)%fix_point_num))
do i = 1, remesh_add(ng)%fix_point_num
    Hexa(ng)%fixed_point_group(:,i) = fixed_point_group(1:num,i)
enddo
deallocate(fixed_point_group)


end subroutine adjust_fix_coord