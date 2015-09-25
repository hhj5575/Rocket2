!=====================================================================================
subroutine inner_smooth(rn, sn, s_nodes, nn, node)

use physical_domain

implicit none

integer, intent(in) :: rn, nn, sn, s_nodes(sn)
real, intent(inout) :: node(2,nn)

integer :: i, j, k, soft_n, roof_i, ne, ng, sp, ep, inc, adjn
integer :: nni(16,nn), nni_n(nn), soft(2,16), bound(nn)
integer :: nei_n(nn), nei1(7,nn), nei2(7,nn), nei3(7,nn)
real :: soft_sum(2), ang(3), vec(2), nndis(16), coe
logical :: s_check(nn)
integer, allocatable :: ele(:,:)

!print *, " >> inner smoothing!!!"
nni = 0;  nni_n = 0; nei_n=0; nei1 = 0; nei2 = 0; nei3 = 0
ne = region(rn)%num_elements
ng = region(rn)%num_groups
bound = 0
s_check = .FALSE.
do i = 1, ng
    call set_skin(rn, i, sp, adjn)
    do j = sp, sp+adjn-1
        bound(region(rn)%skin_nodes(j)) = 1
    enddo
enddo
allocate ( ele(4,ne) )
ele(:,1:ne) = region(rn)%element(4:7,:)

do i = 1, nn
    if ( bound(i) == 0 ) then
        do j = 1, sn
            if (i == s_nodes(j)) then
                s_check(i) = .TRUE.
                exit
            endif
        enddo
        
        if (s_check(i)) then
            soft_n = 0
            soft_sum(1) = 0.
            soft_sum(2) = 0.

            do j = 1, ne
                do k = 1, 4
                    if (ele(k,j) == i) then
                        soft_n = soft_n + 1
                        nei_n(i) = nei_n(i) + 1
                        if (k == 1) then
                            soft(1, soft_n) = ele(2,j)
                            soft(2, soft_n) = ele(4,j)
                            nei1(nei_n(i), i) = ele(2,j)
                            nei2(nei_n(i), i) = ele(3,j)
                            nei3(nei_n(i), i) = ele(4,j)
                        elseif (K == 2) then
                            soft(1, soft_n) = ele(k+1,j)
                            soft(2, soft_n) = ele(k-1,j)
                            nei1(nei_n(i), i) = ele(3,j)
                            nei2(nei_n(i), i) = ele(4,j)
                            nei3(nei_n(i), i) = ele(1,j)
                        elseif (k == 3) then
                            soft(1, soft_n) = ele(k+1,j)
                            soft(2, soft_n) = ele(k-1,j)
                            nei1(nei_n(i), i) = ele(4,j)
                            nei2(nei_n(i), i) = ele(1,j)
                            nei3(nei_n(i), i) = ele(2,j)
                        elseif (k == 4) then
                            soft(1, soft_n) = ele(1,j)
                            soft(2, soft_n) = ele(3,j)
                            nei1(nei_n(i), i) = ele(1,j)
                            nei2(nei_n(i), i) = ele(2,j)
                            nei3(nei_n(i), i) = ele(3,j)                    
                        endif
                        soft_n = soft_n + 1
                        if ( k == 3 ) then
                            soft(1, soft_n) = ele(1,j)
                        elseif ( k == 4 ) then
                            soft(1, soft_n) = ele(2,j)
                        else
                            soft(1, soft_n) = ele(k+2,j)
                        endif
                    endif
                enddo
            enddo
            nni_n(i) = soft_n
            nni(1:nni_n(i),i) = soft(1,1:nni_n(i))
        endif
    endif
enddo

do roof_i = 1, 6
    do i = 1, nn
        if ( bound(i) == 0 .AND. s_check(i)) then
        !if ( bound(i) == 0) then
            soft_sum = 0.0
            do j = 1, nei_n(i)
                soft_sum = soft_sum + node(:,nei1(j,i))
            enddo
            node(:,i) = soft_sum(:)/(nei_n(i)*1.0)
        endif
    enddo
    !do i = 1, nn
    !    if ( bound(i) == 0) then
    !        soft_sum = 0.0
    !        do j = 1, nei_n(i)
    !            soft_sum = soft_sum + (node(:,nei1(j,i))+node(:,nei3(j,i))-node(:,nei2(j,i)))
    !        enddo
    !        node(:,i) = soft_sum(:)/(nei_n(i)*1.0)
    !    endif
    !enddo
            !soft_sum = 0.
            !do j = 1, nei_n(i)
            !    call calc_angle(node(:,nei1(j,i)), node(:,nei2(j,i)), node(:,i), ang(1))
            !    call calc_angle(node(:,i), node(:,nei2(j,i)), node(:,nei3(j,i)), ang(2))
            !    ang(3) = 0.5*(ang(1)-ang(2))
            !    vec = node(:,i) - node(:,nei2(j,i))
            !    
            !    soft_sum(1) = soft_sum(1) + (node(1,nei2(j,i)) + vec(1)*cos(ang(3)*pi/180.) - vec(2)*sin(ang(3)*pi/180.))
            !    soft_sum(2) = soft_sum(2) + (node(2,nei2(j,i)) + vec(1)*sin(ang(3)*pi/180.) + vec(2)*cos(ang(3)*pi/180.))
            !    
            !    call calc_angle(node(:,i), node(:,nei1(j,i)), node(:,nei2(j,i)), ang(1))
            !    ang(3) = 90.0 - ang(1)
            !    vec = node(:,i) - node(:,nei1(j,i))
            !    soft_sum(1) = soft_sum(1) + (node(1,nei1(j,i)) + vec(1)*cos(ang(3)*pi/180.) - vec(2)*sin(ang(3)*pi/180.))
            !    soft_sum(2) = soft_sum(2) + (node(2,nei1(j,i)) + vec(1)*sin(ang(3)*pi/180.) + vec(2)*cos(ang(3)*pi/180.))
            !    
            !    call calc_angle(node(:,nei2(j,i)), node(:,nei3(j,i)), node(:,i), ang(1))
            !    ang(3) = ang(1) - 90.0
            !    vec = node(:,i) - node(:,nei3(j,i))
            !    soft_sum(1) = soft_sum(1) + (node(1,nei3(j,i)) + vec(1)*cos(ang(3)*pi/180.) - vec(2)*sin(ang(3)*pi/180.))
            !    soft_sum(2) = soft_sum(2) + (node(2,nei3(j,i)) + vec(1)*sin(ang(3)*pi/180.) + vec(2)*cos(ang(3)*pi/180.))
            !enddo
            !nei_n(i) = 0
            !nni_n(i) = 0
            !node(:,i) = soft_sum(:)/(nni_n(i)+nei_n(i))
            !node(:,i) = soft_sum(:)/nei_n(i)
            !node(:,i) = soft_sum(:)/nni_n(i)
    !    endif
    !enddo
enddo

end subroutine inner_smooth

subroutine surface_smooth_3D(rn, nn, node)

use surface_domain
implicit none

integer, intent(in) :: rn, nn
real, intent(inout) :: node(3,nn)

integer :: i, j, k, q, ten, cbn, loop, num, num_nci
integer :: temp_ele(4,8), near_node(2), nci(8), surf_ele(4), this(3,2)
integer :: nele(8,nn), nele_num(nn)
real :: t, u, v, quad_val, area, tol_area, dot_a, dot_b, sd, suv, val, temp_area, s
real :: p(3,4), su(3), sv(3), sc(3), q1(4), q2(4), dis(3)
real :: pl(3), temp_pl(4), mp(3), vec(3), move_vec(3), cp(3)

quad_val = sqrt(3.0)/6.0
q1 = (/ -1.0, 1.0, -1.0, 1.0 /)
q2 = (/ -1.0, -1.0, 1.0, 1.0 /)

this(:,1) = (/ 1, 2, 3 /)
this(:,2) = (/ 3, 4, 1 /)

do i = 1, mesh_set(rn)%bn
    cbn = mesh_set(rn)%bound_node(i)
    if (mesh_set(rn)%bound_lv(i) == 1 .or. mesh_set(rn)%bound_lv(i) == 2) then
        ten = 0
        do j = 1, mesh_set(rn)%bfn
            do k = 1, 4
                if (mesh_set(rn)%bound_face(k,j) == i) then
                    ten = ten + 1
                    nele(ten,i) = j
                    exit
                endif
            enddo
        enddo
        nele_num(i) = ten
    endif
enddo

do loop = 1, 4
    do i = 1, mesh_set(rn)%bn
        cbn = mesh_set(rn)%bound_node(i)
        if (mesh_set(rn)%bound_lv(i) == 2) then
            ten = nele_num(i) 
            do j = 1, nele_num(i)
                temp_ele(:,j) = mesh_set(rn)%bound_face(:,nele(j,i))
            enddo
            call find_nci(ten, temp_ele(:,1:ten), i, num_nci, nci)
            
            num = 0
            do j = 1, num_nci
                if (mesh_set(rn)%bound_lv(nci(j)) == 2 .or. mesh_set(rn)%bound_lv(nci(j)) == 3) then
                    num = num + 1
                    near_node(num) = mesh_set(rn)%bound_node(nci(j))
                endif
                if (num == 2) exit
            enddo
            
            p(:,1) = node(:,cbn)
            do j = 1, 2
                call calc_length_3D(p(:,1), node(:,near_node(j)), dis(j))
            enddo
            dis(3) = (dis(1)+dis(2))*0.5
            if (dis(1) > dis(2)) then
                num = 1
            else
                num = 2
            endif
            dis(3) = dis(num) - dis(3)
            vec = node(:,near_node(num)) - p(:,1)
            val = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
            if (val > 0.000001) then
                vec = vec/val
            endif
            node(:,cbn) = node(:,cbn) + vec*dis(3)
        endif
    enddo
enddo

do loop = 1, 6
    do i = 1, mesh_set(rn)%bn
        cbn = mesh_set(rn)%bound_node(i)
        if (mesh_set(rn)%bound_lv(i) == 1) then
            cp = node(:,cbn)
            ten = nele_num(i) 
            do j = 1, nele_num(i)
                temp_ele(:,j) = mesh_set(rn)%bound_face(:,nele(j,i))
            enddo
            pl = 0.0
            tol_area = 0.0
            mp = 0.0
            do j = 1, ten
                surf_ele = mesh_set(rn)%bound_node(temp_ele(:,j))
                p = node(:,surf_ele)
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
                    if (suv <= 0.0) then
                        write (*,*) 'area:', suv
                        write (*,*) p(:,1)
                        write (*,*) p(:,2)
                        write (*,*) p(:,3)
                        write (*,*) p(:,4)
                        temp_area = 0.0
                        do q = 1, 2
                            call calc_length_3D(p(:,this(1,q)), p(:,this(2,q)), dis(1))
                            call calc_length_3D(p(:,this(2,q)), p(:,this(3,q)), dis(2))
                            call calc_length_3D(p(:,this(3,q)), p(:,this(1,q)), dis(3))
                            s = (dis(1)+dis(2)+dis(3))/2.
                            if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
                                temp_area = temp_area + sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
                            endif
                        enddo
                        write (*,*) 'temp_area:', temp_area
                        stop
                    endif
                    area = area + suv
                enddo
                area = area * 0.25
                tol_area = tol_area + area
                mp = mp + (p(:,1)+p(:,2)+p(:,3)+p(:,4))*0.25*area
            enddo
            pl = pl / float(ten)
            mp = mp / tol_area
            move_vec = mp-cp
            val = pl(1)*move_vec(1) + pl(2)*move_vec(2) + pl(3)*move_vec(3)
            move_vec = move_vec - val*pl
            node(:,cbn) = cp + move_vec
        endif
    enddo
enddo

end subroutine surface_smooth_3D

subroutine inner_smooth_3D(rn, nn, node)

use surface_domain

implicit none

integer, intent(in) :: rn, nn
real, intent(inout) :: node(3,nn)

integer :: i, j, k, q, inn, num, cn, count, root, bn, ne, ng
integer :: nci(20), face_nn(4,6), this(3)
real :: cp(3)
logical :: check

ne = mesh_set(rn)%ne
bn = mesh_set(rn)%bn
inn = mesh_set(rn)%inn

do root = 1, 6
    do i = 1, inn
        cn = mesh_set(rn)%inner_node(i)
        nci = mesh_set(rn)%inn_nci(:,i)
        num = mesh_set(rn)%inn_nci_num(i)
        cp = 0.0
        count = 0
        do j = 1, num
            if (nci(j) > ne) then
                write (*,*) nci(j), ne
                stop
            endif
            !call make_face_nn(8, 4, region(rn)%element(4:11,nci(j)), face_nn)
            call make_face_nn(8, 4, mesh_set(rn)%element(:,nci(j)), face_nn)
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
                    cp = cp + node(:,this(1)) + node(:,this(2)) ! + node(:,this(3))
                endif
            enddo
        enddo
        cp = cp / real(count)
        node(:,cn) = cp
    enddo
enddo

end subroutine inner_smooth_3D


Subroutine bound_smooth(bn, adj, nndis, sn, s_nodes, nn, node)

implicit none

integer, intent(in) :: bn, nn, sn
integer, intent(in) :: adj(bn), s_nodes(sn)
real, intent(in) :: nndis(bn)
real, intent(inout) :: node(2,nn)

integer :: i, j, k, this(3), cp, sp, ep, inc, iter, count
real :: alpha, beta, ratio(2), cri_ang(2), dis(3), vec(2), tol_dis(2), cur_dis(2)
logical :: check, s_check, end_check
integer :: re_adj(bn)
real :: re_nndis(bn)

cri_ang = (/ 178.0, 182.0 /)
do i = 1, bn
    if (i == bn) then
        this(:) = (/ i-1, i, 1 /)
    elseif (i == 1) then
        this(:) = (/ bn, i, i+1 /)
    else
        this(:) = (/ i-1, i, i+1 /)
    endif
    call calc_angle(node(:,adj(this(1))), node(:,adj(this(2))), node(:,adj(this(3))), alpha)
    
    if ((alpha <= cri_ang(1)) .or. (alpha >= cri_ang(2))) then
        count = 0
        do j = i, bn
            count = count + 1
            re_adj(count) = adj(j)
            re_nndis(count) = nndis(j)
        enddo
        do j = 1, i-1
            count = count + 1
            re_adj(count) = adj(j)
            re_nndis(count) = nndis(j)
        enddo
        exit
    endif
enddo

cp = 1
sp = 1
tol_dis = 0.0
end_check = .false.
b_smooth: do ! check skin node
    cp = cp + 1
    if (cp > bn) exit ! b_smooth
    
    if (cp == bn) then
        this(:) = (/ cp-1, cp, 1 /)
    elseif (cp == 1) then
        this(:) = (/ bn, cp, cp+1 /)
    else
        this(:) = (/ cp-1, cp, cp+1 /)
    endif
    call calc_angle(node(:,re_adj(this(1))), node(:,re_adj(this(2))), node(:,re_adj(this(3))), alpha)
    call calc_len(node(:,re_adj(this(2))), node(:,re_adj(this(1))), dis(1))
    tol_dis(1) = tol_dis(1) + re_nndis(this(1))
    tol_dis(2) = tol_dis(2) + dis(1)
    if (this(3) == 1) then
        call calc_len(node(:,re_adj(this(3))), node(:,re_adj(this(2))), dis(1))
        tol_dis(1) = tol_dis(1) + re_nndis(this(2))
        tol_dis(2) = tol_dis(2) + dis(1)
    endif
    
    if ( ((alpha <= cri_ang(1)) .or. (alpha >= cri_ang(2))) .or. this(3) == 1) then
        ep = cp
        if (this(3) == 1) then
            ep = cp + 1
            end_check = .true.
        endif

        cur_dis = 0.0
        do i = sp+1, ep-1
            if (i == bn) then
                this(:) = (/ i-1, i, 1 /)
            elseif (i == 1) then
                this(:) = (/ bn, i, i+1 /)
            else
                this(:) = (/ i-1, i, i+1 /)
            endif
            cur_dis(1) = cur_dis(1) + re_nndis(this(1))
            dis(1) = tol_dis(2)*cur_dis(1)/tol_dis(1)
            dis(1) = dis(1) - cur_dis(2)
            call calc_len(node(:,re_adj(this(2))), node(:,re_adj(this(1))), dis(2))
            call calc_len(node(:,re_adj(this(3))), node(:,re_adj(this(2))), dis(3))
            !write (*,'(A,I,1X,3(F12.7,2X))') 'dis:', i, dis
            if (dis(1) <= dis(2)) then
                vec(:) = (node(:,re_adj(this(2)))-node(:,re_adj(this(1))))/dis(2)
                if (dis(1) <= dis(2)*0.15) then
                    node(:,re_adj(this(2))) = node(:,re_adj(this(1))) + vec(:)*dis(2)*0.15
                    dis(1) = dis(2)*0.15
                else
                    node(:,re_adj(this(2))) = node(:,re_adj(this(1))) + vec(:)*dis(1)
                endif
            else
                vec(:) = (node(:,re_adj(this(3)))-node(:,re_adj(this(2))))/dis(3)
                if ((dis(1)-dis(2)) >= dis(3)*0.85) then
                    node(:,re_adj(this(2))) = node(:,re_adj(this(2))) + vec(:)*dis(3)*0.85
                    dis(1) = dis(3)*0.85 + dis(2)
                else
                    node(:,re_adj(this(2))) = node(:,re_adj(this(2))) + vec(:)*(dis(1)-dis(2))
                endif
            endif
            cur_dis(2) = cur_dis(2) + dis(1)
        enddo
        sp = cp
        tol_dis = 0.0
        if (end_check) exit ! b_smooth
    endif
enddo b_smooth

        
!do iter = 1, 20
!    if (mod(iter,2) == 1) then
!        sp = 1;  ep = bn;  inc = 1
!    else
!        sp = bn;  ep = 1;  inc = -1
!    endif
!    do i = sp, ep, inc
!        !s_check = .FALSE.
!        !do j = 1, sn
!        !    if (s_nodes(j) == adj(i)) then
!        !        s_check = .TRUE.
!        !        exit
!        !    endif
!        !enddo
!    
!        !if (s_check) then
!            if (i == bn) then
!                this(:) = (/ i-1, i, 1 /)
!            elseif (i == 1) then
!                this(:) = (/ bn, i, i+1 /)
!            else
!                this(:) = (/ i-1, i, i+1 /)
!            endif
!            call calc_angle(node(:,adj(this(1))), node(:,adj(this(2))), node(:,adj(this(3))), alpha)
!        
!            if ( (alpha >= cri_ang(1)) .AND. (alpha <= cri_ang(2)) ) then
!                ratio(1) = nndis(this(1))/(nndis(this(1))+nndis(this(2)))
!                call calc_len(node(:,adj(this(2))), node(:,adj(this(1))), dis(1))
!                call calc_len(node(:,adj(this(2))), node(:,adj(this(3))), dis(2))
!                ratio(2) = dis(1)/(dis(1)+dis(2))
!                if (ratio(1) < ratio(2)) then
!                    vec(:) = (node(:,adj(this(1)))-node(:,adj(this(2))))/dis(1)
!                    dis(3) = dis(1) - ratio(1)*(dis(1)+dis(2))
!                else
!                    vec(:) = (node(:,adj(this(3)))-node(:,adj(this(2))))/dis(2)
!                    dis(3) = ratio(1)*(dis(1)+dis(2)) - dis(1)
!                endif
!                node(:,adj(this(2))) = node(:,adj(this(2))) + vec(:)*dis(3)
!            endif
!        !endif
!    enddo
!enddo

end subroutine bound_smooth
