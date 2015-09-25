! 점과 직선 사이의 거리
Subroutine dist_between_point_line(a, b, c, off, dis)
            
implicit none

real, intent(in) :: a(2), b(2), c(2), off
real, intent(inout) :: dis

real :: val(4)

if ( abs(b(1) - a(1)) <= off * 0.01) then
    dis = abs( c(1) - a(1) ) 
else
    val(1) = -(b(2)-a(2))/(b(1)-a(1))
    val(2) = 1
    val(3) = -(b(2)+val(1)*b(1))

    dis = abs( (val(1)*c(1)+val(2)*c(2)+val(3))/sqrt(val(1)**2+val(2)**2) )
endif

end subroutine dist_between_point_line
    

Subroutine near_point(a, b, c, cp, off, check)
            
implicit none

real, intent(in) :: a(2), b(2), c(2), off
real, intent(inout) :: cp(2)
logical, intent(inout) :: check

real :: val(4), dis, t

check = .FALSE.
if ( abs(b(1) - a(1)) <= off * 0.01) then
    dis = abs( c(1) - a(1) ) 
else
    val(1) = -(b(2)-a(2))/(b(1)-a(1))
    val(2) = 1
    val(3) = -(b(2)+val(1)*b(1))

    dis = abs( (val(1)*c(1)+val(2)*c(2)+val(3))/sqrt(val(1)**2+val(2)**2) )
endif

if ( dis <= off ) then
    val(1) = b(1)-a(1)
    val(2) = a(1)-c(1)
    val(3) = b(2)-a(2)
    val(4) = a(2)-c(2)

    t = -(val(1)*val(2) + val(3)*val(4))
    t = t/(val(1)**2+val(3)**2)
    if ( t >= -0.001 .AND. t <= 1.001 ) then
        check = .TRUE.
        cp(:) = a(:) + t*(b(:)-a(:))
    endif   
endif

end Subroutine near_point
!===============================================================================
! 점의 수가 짝수이면 홀수로 바꿈
subroutine make_even_number(nn, node)

implicit none

integer, intent(inout) :: nn
real, intent(inout) :: node(2,nn)

integer :: i, i_r, temp_num
real :: max_len, seg, dis, vector(2)

max_len = 0.
do i = 1, nn-1
    i_r = i+1
    if ( i == nn-1 ) i_r = 1
    call calc_len(node(:,i), node(:,i_r), dis)
    
    if ( max_len < dis ) then
        max_len = dis
        temp_num = i 
    endif
enddo
	
do i = nn-1, temp_num + 1, -1
	node(:,i+1) = node(:,i)
enddo
seg = max_len * 0.5

i_r = temp_num + 1
if ( temp_num == nn-1 ) i_r = 1
vector(:) = ( node(:,i_r) - node(:,temp_num) ) / max_len
node(:,temp_num+1) = node(:,temp_num) + (seg * vector(:))

end subroutine

!===============================================================================
subroutine two_point_ave(node)

implicit none

real :: node(2,3)

real :: dis(3), vec(2), ang
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

call calc_len(node(:,1), node(:,2), dis(1))
call calc_len(node(:,2), node(:,3), dis(2))
call calc_angle(node(:,1), node(:,2), node(:,3), ang)
dis(3) = 0.5*(dis(1)+dis(2))

if (dis(1) .LE. dis(3)) then
	vec(:) = ((node(:,3)-node(:,2))/dis(2)) * (dis(3)-dis(1))! * cos((beta - 180) / 2. * pi / 180.)
else
	vec(:) = ((node(:,2)-node(:,1))/dis(1)) * (dis(3)-dis(1)) * cos((ang - 180.) / 2. * pi / 180.)
endif
node(:,2) = node(:,2) + vec(:)

end subroutine
!================================================================================
! 경계의 총 길이
subroutine calc_tol_len(adjn, adj, nn, node, tol_dis)

implicit none

integer :: adjn, nn, adj(adjn)
real :: tol_dis, node(2,nn)

integer :: i, i_r
real :: dis

tol_dis = 0
do i = 1, adjn
    i_r = i + 1
    if ( i == adjn ) i_r = 1
    
    call calc_len(node(:,adj(i)), node(:,adj(i_r)), dis)
    tol_dis = tol_dis + dis
enddo

end subroutine
!================================================================================
! 형상의 넓이
subroutine calc_tol_area(nn, node, ne, ele, tol_area)

implicit none

integer :: ne, nn, ele(4,ne)
real :: tol_area, node(2,nn)

integer :: i
real :: area, radial(2)

tol_area = 0.0
do i = 1, ne
    radial = 0.0
    radial(1) = node(2,ele(1,i))*node(1,ele(2,i)) + node(2,ele(2,i))*node(1,ele(3,i))
    radial(1) = radial(1) + node(2,ele(3,i))*node(1,ele(4,i)) + node(2,ele(4,i))*node(1,ele(1,i))
    radial(2) = node(1,ele(1,i))*node(2,ele(2,i)) + node(1,ele(2,i))*node(2,ele(3,i))
    radial(2) = radial(2) + node(1,ele(3,i))*node(2,ele(4,i)) + node(1,ele(4,i))*node(2,ele(1,i))
    area = abs(radial(1)-radial(2)) * 0.5
    tol_area = tol_area + area
enddo

end subroutine
!================================================================================

subroutine load_interporation( load_coord, load_value, nn, coord, sp, bn, inte_load_num, inte_load_value, tol, load_check )

implicit none

integer, intent(in) :: nn, sp, bn
real, intent(in) :: load_coord(2), load_value, coord(2,nn), tol
integer, intent(inout) :: inte_load_num(2)
real, intent(inout) :: inte_load_value(2)
logical, intent(inout) :: load_check

integer :: i, i_r
real :: cp(2), dis(3)
logical :: check
!--------------------------------------------------------------------------------
load_check = .FALSE.
do i = sp, sp+bn-1
    i_r = i + 1
    if ( i == bn ) i_r = sp
    
    call near_point(coord(:,i), coord(:,i_r), load_coord, cp, tol, check)
    if ( check ) then
        inte_load_num(:) = (/ i, i_r /)
        call calc_len( coord(:,i), coord(:,i_r), dis(1) )
        call calc_len( coord(:,i), cp(:), dis(2) )
        call calc_len( coord(:,i_r), cp(:), dis(3) )
        inte_load_value(1) = dis(3)/dis(1) * load_value
        inte_load_value(2) = dis(2)/dis(1) * load_value
        load_check = .TRUE.
        exit
    endif
enddo

end subroutine

!================================================================================

subroutine sep_load_interporation( load_coord, load_value, coord, nn, sn, bn, inte_load_num, inte_load_value, tol, check )

implicit none

integer, intent(in) :: nn, sn, bn
real, intent(in) :: load_coord(2), load_value, coord(2,nn), tol
integer, intent(inout) :: inte_load_num(2)
real, intent(inout) :: inte_load_value(2)
logical, intent(inout) :: check

integer :: i, i_r
real :: cp(2), dis(3)
logical :: find_check, near_check
!--------------------------------------------------------------------------------
find_check = .FALSE.
do i = sn, bn
    i_r = i + 1
    if ( i == bn ) i_r = sn
    
    call near_point(coord(:,i), coord(:,i_r), load_coord, cp, tol, near_check)
    if ( near_check ) then
        inte_load_num(:) = (/ i, i_r /)
        call calc_len( coord(:,i), coord(:,i_r), dis(1) )
        call calc_len( coord(:,i), cp(:), dis(2) )
        call calc_len( coord(:,i_r), cp(:), dis(3) )
        inte_load_value(1) = dis(3)/dis(1) * load_value
        inte_load_value(2) = dis(2)/dis(1) * load_value
        find_check = .TRUE.
        exit
    endif
enddo

if ( find_check == .FALSE. ) then
    check = .FALSE.
else
    check = .TRUE.
endif

end subroutine

!================================================================================

subroutine get_num_contact ( unit_num, num_contact )

implicit none

integer :: unit_num, num_contact

integer :: a(4)
character :: keyword

rewind( unit=unit_num ) 
read (unit_num,*) keyword 
read (unit_num,*) keyword, a(:), num_contact
rewind( unit=unit_num ) 

end subroutine

!================================================================================

! 직선 사이에 점을 포함여부 확인
Subroutine contain_point(a, b, c, off, check)

implicit none

real, intent(in) :: a(2), b(2), c(2), off
logical, intent(inout) :: check

real :: val(4)

check = .FALSE.
call calc_len(a, c, val(1))
call calc_len(b, c, val(2))
call calc_len(a, b, val(3))
if ( abs(val(3)-val(1)-val(2)) <= off ) then
    check = .TRUE.
else
    check = .FALSE.
endif

end Subroutine

!=======================================================================================


subroutine calc_tol_length(nn, node, adjn, adj, tol_dis)

implicit none

integer, intent(in) :: nn, adjn, adj(adjn)
real, intent(in) :: node(2,nn)
real :: tol_dis

integer :: i, i_r
real :: dis

tol_dis = 0
do i = 1, adjn
    i_r = i + 1
    if ( i == adjn ) i_r = 1
    
    call calc_len(node(:,adj(i)), node(:,adj(i_r)), dis)
    tol_dis = tol_dis + dis
enddo

end subroutine


Subroutine calc_cross_point(l_a, l_b, l_c, l_d, cross)

real, intent(in) :: l_a(2), l_b(2), l_c(2), l_d(2)
real, intent(out) :: cross(2)

real :: x, dis, tol

tol = 10e-4
call calc_len(l_a, l_b, dis)
if (abs(((l_d(2)-l_c(2)) * (l_b(1)-l_a(1))) - ((l_d(1)-l_c(1)) * (l_b(2)-l_a(2)))) < dis*tol) then
    x = 0.5
else
    x = (((l_d(1) - l_c(1)) * (l_a(2) - l_c(2))) - ((l_d(2) - l_c(2)) * (l_a(1) - l_c(1)))) / (((l_d(2) - l_c(2)) * (l_b(1) - l_a(1))) - ((l_d(1) - l_c(1)) * (l_b(2) - l_a(2))))
endif
cross = l_a + x * (l_b - l_a)

end subroutine calc_cross_point



subroutine ut_jacobi(coord, ut, vec, cur_ut, error)

implicit none

real, intent(in) :: coord(2,4), ut(2,4), vec(2)
real, intent(out) :: cur_ut(2)
integer, intent(out) :: error

real :: f1 ,f2 , det_j , u , v
real :: xi, eta ,norm, tol, temp
real :: jacobi(4)
integer :: j

error = 0 
xi = 0.0 ; eta = 0.0
tol = 1.0e-6
j = 0 
do 
    f1 = ((1-xi)*(1-eta)*coord(1,1)+(1+xi)*(1-eta)*coord(1,2)+ &
   (1+xi)*(1+eta)*coord(1,3)+(1-xi)*(1+eta)*coord(1,4))*0.25-vec(1)
    f2 = ((1-xi)*(1-eta)*coord(2,1)+(1+xi)*(1-eta)*coord(2,2)+ &
   (1+xi)*(1+eta)*coord(2,3)+(1-xi)*(1+eta)*coord(2,4))*0.25-vec(2)
   
    jacobi(1) = (-1+eta)*coord(1,1) +(1-eta)*coord(1,2) &
           + (1+eta)*coord(1,3) +(-1-eta)*coord(1,4)
    jacobi(3) = (-1+eta)*coord(2,1) +(1-eta)*coord(2,2) &
           + (1+eta)*coord(2,3) +(-1-eta)*coord(2,4)
    jacobi(2) = (-1+xi)*coord(1,1) +(-1-xi)*coord(1,2) &
           + (1+xi)*coord(1,3) +(1-xi)*coord(1,4)
    jacobi(4) = (-1+xi)*coord(2,1) +(-1-xi)*coord(2,2) &
           + (1+xi)*coord(2,3) +(1-xi)*coord(2,4)

    jacobi = jacobi * 0.25
    det_j = jacobi(1)*jacobi(4)-jacobi(2)*jacobi(3)
    temp = jacobi(1)
    jacobi(1) = jacobi(4)
    jacobi(2) = -jacobi(2)
    jacobi(3) = -jacobi(3)
    jacobi(4) = temp
    jacobi = jacobi / det_j

    u = jacobi(1) * f1 + jacobi(2) * f2
    v = jacobi(3) * f1 + jacobi(4) * f2
    xi = xi - u
    eta = eta - v

    norm = sqrt( u**2 + v**2 )
    if (norm < tol) exit

    j = j + 1 
    if (j > 100) then
        error = 1
        exit
    endif
enddo

cur_ut(1) = ( (1-xi)*(1-eta)*ut(1,1) + (1+xi)*(1-eta)*ut(1,2)+ &
(1+xi)*(1+eta)*ut(1,3)+(1-xi)*(1+eta)*ut(1,4) )*0.25
cur_ut(2) = ((1-xi)*(1-eta)*ut(2,1) + (1+xi)*(1-eta)*ut(2,2)+ &
(1+xi)*(1+eta)*ut(2,3) + (1-xi)*(1+eta)*ut(2,4))*0.25

end subroutine ut_jacobi



subroutine matrixinv(a,b,n)
! subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
integer :: i,j,k,l,n,irow
real:: big,a(n,n),b(n,n),dum

!build the identity matrix
do i = 1,n
    do j = 1,n
        b(i,j) = 0.0
    end do
    b(i,i) = 1.0
end do

do i = 1,n ! this is the big loop over all the columns of a(n,n)
    ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot 
    ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
    big = a(i,i)
    do j = i,n
        if (a(j,i).gt.big) then
            big = a(j,i)
            irow = j
        end if
    end do
        ! interchange lines i with irow for both a() and b() matrices
    if (big.gt.a(i,i)) then
        do k = 1,n
            dum = a(i,k) ! matrix a()
            a(i,k) = a(irow,k)
            a(irow,k) = dum
            dum = b(i,k) ! matrix b()
            b(i,k) = b(irow,k)
            b(irow,k) = dum
        end do
    end if
    ! divide all entries in line i from a(i,j) by the value a(i,i); 
    ! same operation for the identity matrix
    dum = a(i,i)
    do j = 1,n
        a(i,j) = a(i,j)/dum
        b(i,j) = b(i,j)/dum
    end do
        ! make zero all entries in the column a(j,i); same operation for indent()
    do j = i+1,n
        dum = a(j,i)
        do k = 1,n
            a(j,k) = a(j,k) - dum*a(i,k)
            b(j,k) = b(j,k) - dum*b(i,k)
        end do
    end do
end do

! substract appropiate multiple of row j from row j-1 
do i = 1,n-1
    do j = i+1,n
        dum = a(i,j)
        do l = 1,n
            a(i,l) = a(i,l)-dum*a(j,l)
            b(i,l) = b(i,l)-dum*b(j,l)
        end do
    end do
end do
 
end subroutine matrixinv



subroutine write_adj(nn, nodes, ne, ele, sp, adjn, adj)

implicit none

integer, intent(in) :: nn, nodes(nn), ne, ele(4,ne), sp
integer, intent(inout) :: adjn, adj(nn)

integer :: i, j, k, cn, soft_n, bn, sum_val
integer :: bound(nn), soft(2,8), bound_conn(2,nn)

bound_conn = 0
bound = 0
bn = 0
do i = 1, nn
    cn = nodes(i)
	soft_n = 0
    soft = 0
	do j = 1, ne
		do k = 1, 4
			if (cn == ele(k,j)) then 
				soft_n = soft_n + 1
				if (k .EQ. 4) then
					soft(1, soft_n) = ele(1,j)
					soft(2, soft_n) = ele(3,j)
				else if (k .EQ. 1) then
					soft(1, soft_n) = ele(2,j)
					soft(2, soft_n) = ele(4,j)
				else
					soft(1, soft_n) = ele(k+1,j)
					soft(2, soft_n) = ele(K-1,j)
				endif
			endif
		enddo
	enddo
	if (sum(soft(1,1:soft_n)) /= sum(soft(2,1:soft_n))) then
        sum_val = 0
        do j = 1, soft_n 
            do k = 1, soft_n
                if ( soft(1,j) == soft(2,k) ) then
                    sum_val = sum_val + soft(1,j)
                    exit
                endif
            enddo
        enddo
        bound_conn(1,i) = sum(soft(1,1:soft_n)) - sum_val
        bound_conn(2,i) = sum(soft(2,1:soft_n)) - sum_val
        bn = bn + 1
        bound(bn) = cn
    endif
enddo

adjn = 1
adj(1) = sp
do
    cn = adj(adjn)
    do i = 1, nn
        if ( nodes(i) == cn ) exit
    enddo
    if ( sp == bound_conn(1,i) ) exit
    
    adjn = adjn + 1
    adj(adjn) = bound_conn(1,i)
enddo

end subroutine write_adj


subroutine thin_check(check, nn, coord, ng, adjn, adj, cont_slp, min_d)

implicit none

integer, intent(in) :: nn, ng, adjn, adj(adjn), cont_slp(4,ng)
real, intent(in) :: coord(2,nn), min_d
logical, intent(inout) :: check

integer :: i, j, k, o, count, cont_num, check_num 
integer :: this(2), seq(2), check_nodes(6,2), cont_nodes(adjn), ang_node(3)
real :: p(2,3), cp(2), tol, ang, cri_ang
logical :: end_check, near_check

cri_ang = 165.0
tol = min_d * 2.5
do i = 1, ng
    this = cont_slp(1:2,i)
    call find_count(adjn, adj, 2, this, count)
    if (count == 2) then
        ! find non-matching nodes
        cont_num = 0
        do j = 1, adjn
            if (adj(j) == this(1)) then
                cont_num = cont_num + 1
                cont_nodes(cont_num) = adj(j)
                seq(1) = j
                exit
            endif
        enddo
        end_check = .false.
        do j = seq(1)+1, adjn
            cont_num = cont_num + 1
            cont_nodes(cont_num) = adj(j)
            if (adj(j) == this(2)) then
                end_check = .true.
                seq(2) = j
                exit
            endif
        enddo
        if (end_check == .false.) then
            do j = 1, adjn
                cont_num = cont_num + 1
                cont_nodes(cont_num) = adj(j)
                if (adj(j) == this(2)) then
                    seq(2) = j
                    exit
                endif
            enddo
        endif
        
        if (seq(1) == 1) then
            check_nodes(:,1) = (/ adj(adjn-5), adj(adjn-4), adj(adjn-3), adj(adjn-2), adj(adjn-1), adj(adjn) /)
        elseif (seq(1) == 2) then
            check_nodes(:,1) = (/ adj(adjn-4), adj(adjn-3), adj(adjn-2), adj(adjn-1), adj(adjn), adj(1) /)
        elseif (seq(1) == 3) then
            check_nodes(:,1) = (/ adj(adjn-3), adj(adjn-2), adj(adjn-1), adj(adjn), adj(1), adj(2) /)
        elseif (seq(1) == 4) then
            check_nodes(:,1) = (/ adj(adjn-2), adj(adjn-1), adj(adjn), adj(1), adj(2), adj(3) /)
        elseif (seq(1) == 5) then
            check_nodes(:,1) = (/ adj(adjn-1), adj(adjn), adj(1), adj(2), adj(3), adj(4) /)
        elseif (seq(1) == 6) then
            check_nodes(:,1) = (/ adj(adjn), adj(1), adj(2), adj(3), adj(4), adj(5) /)
        else
            check_nodes(:,1) = (/ adj(seq(1)-6), adj(seq(1)-5), adj(seq(1)-4), adj(seq(1)-3), adj(seq(1)-2), adj(seq(1)-1) /)
        endif
        
        if (seq(2) == adjn) then
            check_nodes(:,2) = (/ adj(1), adj(2), adj(3), adj(4), adj(5), adj(6) /)
        elseif (seq(2) == adjn-1) then
            check_nodes(:,2) = (/ adj(adjn), adj(1), adj(2), adj(3), adj(4), adj(3) /)
        elseif (seq(2) == adjn-2) then
            check_nodes(:,2) = (/ adj(adjn-1), adj(adjn), adj(1), adj(2), adj(3), adj(2) /)
        elseif (seq(2) == adjn-3) then
            check_nodes(:,2) = (/ adj(adjn-2), adj(adjn-1), adj(adjn), adj(1), adj(2), adj(1) /)
        elseif (seq(2) == adjn-4) then
            check_nodes(:,2) = (/ adj(adjn-3), adj(adjn-2), adj(adjn-1), adj(adjn), adj(1), adj(2) /)
        elseif (seq(2) == adjn-5) then
            check_nodes(:,2) = (/ adj(adjn-4), adj(adjn-3), adj(adjn-2), adj(adjn-1), adj(adjn), adj(1) /)
        else
            check_nodes(:,2) = (/ adj(seq(2)+1), adj(seq(2)+2), adj(seq(2)+3), adj(seq(2)+4), adj(seq(2)+5), adj(seq(2)+6) /)
        endif
        write (*,*) 'cont_nodes:', cont_nodes(1:cont_num)
        do j = 1, 2
            write (*,*) 'check_node:', check_nodes(:,j)
            check_num = 0
            do k = 2, 5
                p(:,1) = coord(:,check_nodes(k-1,j))
                p(:,2) = coord(:,check_nodes(k+1,j))
                p(:,3) = coord(:,check_nodes(k,j))
                call calc_angle(p(:,1), p(:,3), p(:,2), ang)
                !if (ang <= cri_ang) then
                    do o = 1, cont_num-1
                        p(:,1) = coord(:,cont_nodes(o))
                        p(:,2) = coord(:,cont_nodes(o+1))
                        call near_point(p(:,1), p(:,2), p(:,3), cp, tol, near_check)
                        if (near_check) then
                            check_num = check_num + 1
                            exit
                        endif
                    enddo
                !endif
                if (check_num /= k-1) exit
            enddo
            if (check_num == 4) then
                check = .false.
                exit
            endif
        enddo
        
        exit
    endif
enddo

write (*,*) 'min_dis:', min_d
write (*,*) 'check thin surface:', check

end subroutine thin_check



subroutine find_contact_point(sp, lp, target_adjn, target_adj, flag, cn, cp)

implicit none

integer, intent(in) :: sp, lp, target_adjn, target_adj(target_adjn), flag
integer, intent(inout) :: cn, cp(target_adjn)

integer :: i, j
logical :: check

cn = 0    
check = .FALSE.
do i = 1, target_adjn
    if ( sp == target_adj(i) ) then
        if (flag == 2) then
            cn = cn + 1
            if (i == 1) then
                cp(cn) = target_adj(target_adjn)
            else
                cp(cn) = target_adj(i-1)
            endif
        endif
        cn = cn + 1
        cp(cn) = target_adj(i)
        check = .TRUE.
    elseif ( lp == target_adj(i) .AND. check == .TRUE. ) then
        cn = cn + 1
        cp(cn) = target_adj(i)
        if (flag == 2) then
            cn = cn + 1
            if (i == target_adjn) then
                cp(cn) = target_adj(1)
            else
                cp(cn) = target_adj(i+1)
            endif
        endif
        check = .FALSE.
        exit
    elseif ( check ) then
        cn = cn + 1
        cp(cn) = target_adj(i)
    endif
enddo
if ( check ) then
    do i = 1, target_adjn
        if ( lp == target_adj(i) ) then
            cn = cn + 1
            cp(cn) = target_adj(i)
            if (flag == 2) then
                cn = cn + 1
                if (i == target_adjn) then
                    cp(cn) = target_adj(1)
                else
                    cp(cn) = target_adj(i+1)
                endif
            endif
            exit
        else
            cn = cn + 1
            cp(cn) = target_adj(i)
        endif
    enddo
endif

end subroutine find_contact_point

subroutine read_mat_element(unit_num_elem, num_materials, sub_mesh_region, nn, ne, coord, elem, conn, sub_nn, sub_ne, rn, ng)

implicit none

integer, intent(in) :: unit_num_elem, num_materials, sub_mesh_region, nn, ne, rn, ng
integer, intent(inout) :: sub_nn, sub_ne, elem(8,ne), conn(nn)
real, intent(inout) :: coord(3,nn)

integer :: i, j, status, ele_num, mat_num, num, mat_rn
integer :: temp_elem(8,ne), temp_elem_num(ne), temp_conn(nn)
integer, allocatable :: temp_node(:), del_node(:), rn_node(:)
real :: temp_coord(3,nn)
character(len=5) :: keyword
logical :: check

temp_elem = elem
temp_coord = coord
elem = 0
sub_ne = 0
temp_elem_num = 0
rewind(unit=unit_num_elem)
find_mat_elem: do
    read(unit_num_elem,*, iostat=status) keyword
    if (keyword == '*elem' .or. keyword == '*spri') then
        backspace(unit_num_elem)
        read(unit_num_elem,*, iostat=status) keyword, mat_rn, ele_num, mat_num
        check = .false.
        if (rn == 1 .and. mat_num == ng) then
            check = .true.
        elseif (rn == 2 .and. mat_num > sub_mesh_region) then
            check = .true.
        endif
        if (check) then
            if (ele_num > 200 .and. ele_num < 400) then
                do 
                    read(unit_num_elem,*, iostat=status) keyword
                    backspace(unit_num_elem)
                    if (status == -1 .or. keyword == '*elem' .or. keyword == '*spri' .or. keyword == '*end') then
                        backspace(unit_num_elem)
                        exit
                    else
                        sub_ne = sub_ne + 1
                        read(unit_num_elem,*, iostat=status) temp_elem_num(sub_ne), (elem(j,sub_ne), j=1,8)
                    endif
                enddo
            endif
            if (rn == 1 .and. mat_num == ng) then
                exit find_mat_elem
            elseif (rn == 2 .and. mat_num == num_materials) then
                exit find_mat_elem
            endif
        endif
    endif
enddo find_mat_elem

allocate (del_node(sub_ne*8), temp_node(sub_ne*8), rn_node(nn))
temp_node = 0
del_node = 0
sub_nn = 0
rn_node = 0
do i = 1, sub_ne
    temp_node((i-1)*8+1:(i-1)*8+8) = elem(:,i)
enddo
do i = 1, sub_ne*8
    if (del_node(i) == 0) then
        sub_nn = sub_nn + 1
        rn_node(sub_nn) = temp_node(i)
        
        do j = i+1, sub_ne*8
            if (temp_node(i) == temp_node(j)) then
                del_node(j) = 1
            endif
        enddo
    endif
enddo
conn = 0
temp_conn = 0
do i = 1, sub_nn
    coord(:,i) = temp_coord(:,rn_node(i))
    conn(i) = rn_node(i)
    temp_conn(rn_node(i)) = i
enddo
do i = 1, sub_ne
    elem(:,i) = temp_conn(elem(:,i))
enddo
write (*,*) 'rn, ng, sub_nn, sub_ne:', rn, ng, sub_nn, sub_ne
call write_tecplot_3D(sub_nn, coord(:,1:sub_nn), sub_ne, elem(:,1:sub_ne), rn, ng, 7) 

end subroutine read_mat_element