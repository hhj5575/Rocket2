subroutine set_local_info( nn, ndim, global, nc, connect, size, local )

implicit none

integer :: nn, ndim , size, nc
integer,intent(in) :: connect(nc) 
real :: global(ndim, nn)
real :: local(ndim, size)
integer :: i

local = 0.0
if (size .EQ. 7 ) then 

	do i=1,3
		local(:, i) = global(:, connect(i))
		local(:, size) = local(:, size) + global(:, connect(i))
	enddo
	do i=4,6
		local(:, i) = global(:, connect(i))
	enddo
	local(:, 7) = local(:, 7) / 3

elseif ( size .EQ. 9 ) then

	do i=1,4
		local(:, i) = global(:, connect(i))
	enddo
	do i=1,3
		local(:, i+4) = ( local(:, i) + local(:, i+1) ) *0.5
	enddo
	local(:, 8) = ( local(:, 1) + local(:, 4) ) *0.5
	local(:, 9) = ( local(:, 5) + local(:, 7) ) *0.5

else 

	do i=1,size
		local(:, i) = global(:, connect(i))
	enddo
			
endif

end subroutine



subroutine calc_avg_vel( ndim, nsize , node, velocity,  surface_size,avg )

implicit none

integer :: ndim, nsize, surface_size
integer :: info(2,surface_size)
real :: node( ndim, nsize) , velocity(ndim, nsize) 
real :: avg(surface_size)

real :: point_a(ndim), point_b(ndim), normal(ndim) 
real :: length
real :: avg_vel
integer :: i

select case (nsize)

	case (3)
		info(1,:) = (/ 1,2,3 /)
		info(2,:) = (/ 2,3,1 /)
	case (4)
		info(1,:) = (/ 1,2,3,4 /)
		info(2,:) = (/ 2,3,4,1 /)
	case (7) 
		info(1,:) = (/ 1,4,2,5,3,6,4,6,7 /)
		info(2,:) = (/ 4,2,5,3,6,1,7,7,5 /)
	case (9) 
		info(1,:) = (/ 1,5,2,6,3,7,4,8,5,9,8,9 /)
		info(2,:) = (/ 5,2,6,3,7,4,8,1,9,7,9,6 /)

end select
avg = 0.0
do i= 1, surface_size
	point_a = node( :, info( 1, i) )
	point_b = node( :, info( 2, i) )
	call get_normal(point_a, point_b, normal, length)
	avg(i) = avg_vel( velocity( :, info( 1, i) ), velocity( :, info( 2, i) ) , normal  ) * length
enddo
end subroutine



logical function sharing_node( node , element , nelem )

implicit none

integer :: nelem
integer :: node(2) , element(nelem) 
integer :: i, j

sharing_node = .FALSE.

do i=1,nelem
	if ( node(1) .EQ. element(i) ) then 
		do j=1,nelem
			if ( node(2) .EQ. element(j) ) then
				sharing_node = .TRUE. 
				exit
			endif
		enddo
	endif 
	if ( sharing_node ) exit
enddo 

end function

subroutine calc_sub_area( coord, m  , area , n  )

implicit none 
integer :: m , n 
real :: coord(2,m)
real :: area(n)

area = 0.0
if ( n .EQ. 4 ) then

	area(1)=( coord(1,1)*coord(2,5) + coord(1,5)*coord(2,9)+ &
            coord(1,9)*coord(2,8) + coord(1,8)*coord(2,1)- &
            coord(1,5)*coord(2,1) - coord(1,9)*coord(2,5)- &
            coord(1,8)*coord(2,9) - coord(1,1)*coord(2,8) ) *0.5
         
	area(2)=( coord(1,5)*coord(2,2) + coord(1,2)*coord(2,6 )+ &
            coord(1,6)*coord(2,9) + coord(1,9)*coord(2,5)- &
            coord(1,2)*coord(2,5) - coord(1,6)*coord(2,2)- &
            coord(1,9)*coord(2,6) - coord(1,5)*coord(2,9) ) *0.5

	area(3)=( coord(1,9)*coord(2,6) + coord(1,6)*coord(2,3)+ &
            coord(1,3)*coord(2,7) + coord(1,7)*coord(2,9)- &
            coord(1,6)*coord(2,9) - coord(1,3)*coord(2,6)- &
            coord(1,7)*coord(2,3) - coord(1,9)*coord(2,7) ) *0.5
         
	area(4)=( coord(1,8)*coord(2,9) + coord(1,9)*coord(2,7)+ &
            coord(1,7)*coord(2,4) + coord(1,4)*coord(2,8)- &
            coord(1,9)*coord(2,8) - coord(1,7)*coord(2,9)- &
            coord(1,4)*coord(2,7) - coord(1,8)*coord(2,4) ) *0.5

elseif ( n .EQ. 3 ) then 

	area(1)=( coord(1,1)*coord(2,4) + coord(1,4)*coord(2,7)+ &
            coord(1,7)*coord(2,6) + coord(1,6)*coord(2,1)- &
            coord(1,4)*coord(2,1) - coord(1,7)*coord(2,4)- &
            coord(1,6)*coord(2,7) - coord(1,1)*coord(2,6) ) *0.5
         
	area(2)=( coord(1,4)*coord(2,2) + coord(1,2)*coord(2,5)+ &
            coord(1,5)*coord(2,7) + coord(1,7)*coord(2,4)- &
            coord(1,2)*coord(2,4) - coord(1,5)*coord(2,2)- &
            coord(1,7)*coord(2,5) - coord(1,4)*coord(2,7) ) *0.5

	area(3)=( coord(1,5)*coord(2,3) + coord(1,3)*coord(2,6)+ &
            coord(1,6)*coord(2,7) + coord(1,7)*coord(2,5)- &
            coord(1,3)*coord(2,5) - coord(1,6)*coord(2,3)- &
            coord(1,7)*coord(2,6) - coord(1,5)*coord(2,7) ) *0.5

elseif ( m .EQ. 4 ) then

	area(1)=( coord(1,1)*coord(2,2) + coord(1,2)*coord(2,3)+ &
            coord(1,3)*coord(2,4) + coord(1,4)*coord(2,1)- &
            coord(1,2)*coord(2,1) - coord(1,3)*coord(2,2)- &
            coord(1,4)*coord(2,3) - coord(1,1)*coord(2,4) ) *0.5

elseif ( m .EQ. 3 ) then

	area(1)=( coord(1,1)*coord(2,2) + coord(1,2)*coord(2,3)+ &
            coord(1,3)*coord(2,1) - coord(1,2)*coord(2,1)- &
            coord(1,3)*coord(2,2) - coord(1,1)*coord(2,3) ) *0.5

else
	STOP "Size of dimension is wrong!!"
endif


end subroutine 		





real function avg_vel(a, b, normal )

real, intent(in) :: a(2) , b(2) , normal(2)
real :: temp(2)

temp = ( a + b  ) * 0.5 
avg_vel = temp(1)*normal(1) + temp(2)*normal(2)

end function



subroutine get_normal(a, b, normal, length)

implicit none

real , intent(in) :: a(2), b(2)
real :: normal(2), length
real :: temp(2)

temp = b-a
length = sqrt( temp(1)**2 + temp(2)**2 ) 
normal(1) = temp(2) /length
normal(2) = -temp(1) / length
end subroutine

subroutine get_position( vec , len1, val , len2  )

implicit none 

integer, intent(in) :: len1, len2 
integer, intent(in) :: vec(len1) 

integer,intent(inout) :: val(len2)

integer :: i, j, value



do i=1,len2
	value = val(i)
	do j=1,len1
		if ( vec(j) == value ) then
			val(i) = j
			exit
		endif
	enddo
enddo



end subroutine



subroutine set_quad_num( position, num_quad )

implicit none
integer, intent(in) :: num_quad
integer, intent(inout) :: position(2)

integer :: quad_point(10)
integer :: num1 , num2 
integer :: temp, i, j
integer :: surface_num

logical :: positive


num1 = position(1)
num2 = position(2)

temp =  num2 - num1

if (temp > 0 ) then
	positive = .TRUE.
else 
	positive = .FALSE.
endif

temp = abs(temp) 	

if ( temp == 1 ) then
	surface_num = num1
else
	surface_num = num2 
endif

do i=1,num_quad
	quad_point(i) = i
enddo
quad_point(num_quad+1) = 1

if (positive) then
	position(1) = quad_point(surface_num)
	position(2) = quad_point(surface_num+1)
else
	position(2) = quad_point(surface_num) 
	position(1) = quad_point(surface_num+1)
endif

end subroutine 


Subroutine calc_len(pa, pb, length)

implicit none

real, intent(in) :: pa(1:2,1), pb(1:2,1)
real, intent(out) :: length

length = sqrt( ( (pa(1,1) - pb(1,1)) ** 2) + ((pa(2,1) - pb(2,1)) ** 2 ) )

end subroutine



subroutine inelement( coord , vec, in_element, diff_area, tol )

implicit none

real, intent(in) ::  coord(2,4), vec(2)
real, intent(inout) :: diff_area
logical, intent(inout) :: in_element

real :: x, y, direction
integer :: i, j
real :: nrm, tol, area, calc_length, heron, local_area
real :: for_tri(2,5), for_local(2, 4) , len(3)

in_element = .FALSE.

area = abs(coord(1,1)*coord(2,2) + coord(1,2)*coord(2,3)+ &
            coord(1,3)*coord(2,4) + coord(1,4)*coord(2,1)- &
            coord(1,2)*coord(2,1) - coord(1,3)*coord(2,2)- &
            coord(1,4)*coord(2,3) - coord(1,1)*coord(2,4))*0.5

for_tri(:,1:4) = coord 
for_tri(:,5) = coord(:,1)
local_area = 0.0 


do j = 1, 4 
        for_local(:,1 ) = vec
        for_local(:,2:3 ) = for_tri(:, j:j+1) 
        for_local(:,4 ) = vec
        
        !do i= 1 , 3 
        !    len(i) = calc_length( for_local(:,i) , for_local(:,i+1))
        !enddo 
        !
        !heron = sum( len )*0.5 
        !if (heron-len(1)>tol .AND. heron-len(2)>tol .AND. heron-len(3)>tol) then
        !    local_area  = local_area + sqrt(heron*(heron-len(1))*(heron-len(2))*(heron-len(3)) ) 
        !endif
        local_area = local_area + 0.5*abs(for_local(1,1)*for_local(2,2)+for_local(1,2)*for_local(2,3)+for_local(1,3)*for_local(2,1)- &
                            for_local(1,2)*for_local(2,1)-for_local(1,3)*for_local(2,2)-for_local(1,1)*for_local(2,3))
enddo 

!print *, abs(area-local_area) / area
diff_area = abs(area-local_area)/area
if ((abs(area-local_area)/area) < tol) then
    in_element = .TRUE.
endif

end subroutine



subroutine outelement(coord, in_coord, vec, ut, vt, at, error)

implicit none

integer, intent(inout) :: error
real,intent(in) ::  coord(2,4) , ut(2,4) , vt(2,4) , at(2,4)
real,intent(inout) :: in_coord(2), vec(2,3)

integer :: i, j, min_i, i_r, err
real :: min_len, len, f1 ,f2 ,det_j, u, v, xi, eta, norm, temp, dis
real :: temp_coord(2,4), check_node(2), jacobi(4), cp(2), temp_cp(2)
logical :: in_element, check

min_len = 10e8
min_i = 0
temp_cp = 0.0
do i = 1, 4
    i_r = i + 1
    if ( i == 4 ) i_r = 1
    call calc_len(coord(:,i), coord(:,i_r), dis)
    call near_point(coord(:,i), coord(:,i_r), in_coord, cp, dis*0.01, check )
    if ( check ) then
        call calc_len(in_coord, cp, dis)
        if (dis < min_len) then
            temp_cp = cp
            min_i = i
        endif
    endif
enddo

if (min_i /= 0) then
    in_coord = temp_cp
else
    len = sqrt((coord(1,1)-in_coord(1))**2 + (coord(2,1)-in_coord(2))**2)
    min_len = len
    min_i = 1
    do i = 2, 4
	    len = sqrt((coord(1,i)-in_coord(1))**2 + (coord(2,i)-in_coord(2))**2)
	    if ( min_len > len ) then
		    min_len = len
		    min_i = i
	    endif
    enddo
    in_coord = coord(:,min_i)
endif

call goto_NC( coord, in_coord, ut, vt, at, vec, error )

end subroutine outelement


subroutine outelement_3D(coord, in_coord, vec, ut, vt, at, error)

implicit none

integer, intent(inout) :: error
real,intent(in) ::  coord(3,8) , ut(3,8) , vt(3,8) , at(3,8)
real,intent(inout) :: in_coord(3), vec(3,3)

integer :: i, face_nn(4,6), line(2,12)
real :: t, min_dis, dis
real :: face_coord(3,4), pl(4,2), p(3,2), cp(3), temp_cp(3), pre_p(3)
logical :: parallel, check

face_nn(1,:) = (/ 1, 5, 1, 2, 3, 1 /)
face_nn(2,:) = (/ 2, 8, 5, 6, 7, 4 /)
face_nn(3,:) = (/ 3, 7, 6, 7, 8, 8 /)
face_nn(4,:) = (/ 4, 6, 2, 3, 4, 5 /)

line(1,:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4 /)
line(2,:) = (/ 2, 3, 4, 1, 6, 7, 8, 5, 5, 6, 7, 8 /)

pre_p = in_coord

min_dis = 10e8
do i = 1, 6
    face_coord = coord(:,face_nn(:,i))
    call calc_plane(face_coord(:,1), face_coord(:,2), face_coord(:,3), pl(:,1))
    call calc_plane(face_coord(:,3), face_coord(:,4), face_coord(:,1), pl(:,2))
    pl(:,1) = (pl(:,1)+pl(:,2))*0.5
    dis = abs(pl(1,1)*in_coord(1)+pl(2,1)*in_coord(2)+pl(3,1)*in_coord(3)+pl(4,1))
    dis = dis/sqrt(pl(1,1)**2.0+pl(2,1)**2.0+pl(3,1)**2.0)
    
    if (dis < min_dis) then
        p(:,1) = in_coord + pl(1:3,1)*dis*10.0
        p(:,2) = in_coord - pl(1:3,1)*dis*10.0
        call cross_point_plane(pl(:,1), p(:,1), p(:,2), cp, t, parallel)
        if ( parallel == .FALSE. ) then
            call check_inner_face2(face_coord, cp, check)
            if (check) then
                min_dis = dis
                temp_cp = cp
            endif
        endif
    endif
enddo
do i = 1, 12
    p = coord(:,line(:,i))
    call calc_perpendicular_point(in_coord, p(:,1), p(:,2), cp, t, check)
    if (check) then
        call calc_length_3D(in_coord, cp, dis)
        if (dis < min_dis) then
            min_dis = dis
            temp_cp = cp
        endif
    endif
enddo
do i = 1, 8
    call calc_length_3D(in_coord, coord(:,i), dis)
    if (dis < min_dis) then
        min_dis = dis
        temp_cp = coord(:,i)
    endif
enddo
in_coord = temp_cp

!write (*,'(3(F11.6,1X),A,3(F11.6,1X))') pre_p, '->', in_coord
call goto_NC_3D(coord, in_coord, ut, vt, at, vec, error )

end subroutine outelement_3D


logical function in_element( coord , vec, tol )

real, intent(in) :: coord(2,4), vec(2)
real, intent(in) :: tol
real :: x, y, direction
integer :: i, j
real :: nrm, tol_area, area, local_area, calc_length, heron
real :: for_tri(2,5), for_local(2, 4) , len(3)

in_element = .FALSE.

area =( coord(1,1)*coord(2,2) + coord(1,2)*coord(2,3)+ &
            coord(1,3)*coord(2,4) + coord(1,4)*coord(2,1)- &
            coord(1,2)*coord(2,1) - coord(1,3)*coord(2,2)- &
            coord(1,4)*coord(2,3) - coord(1,1)*coord(2,4) ) *0.5

tol_area = tol*area

for_tri(:,1:4) = coord 
for_tri(:,5) = coord(:,1)
local_area = 0.0 


do j = 1, 4 
        for_local(:,1 ) = vec
        for_local(:,2:3 ) = for_tri(:, j:j+1) 
        for_local(:,4 ) = vec
        
        do i= 1 , 3 
                len(i) = calc_length( for_local(:,i) , for_local(:,i+1))
        enddo 
        
        heron = sum( len )*0.5 


        local_area  = local_area + sqrt(heron*(heron-len(1))*(heron-len(2))*(heron-len(3)) ) 

enddo 
!print *, abs(area-local_area) / area 
if ( ( abs(area-local_area) ) < tol_area ) then
        in_element = .TRUE.
endif
end function 


subroutine check_inner_node(nn, coord, ne, ele, ni, inner_node)

implicit none

integer, intent(in) :: nn, ne, ni, ele(4,ne)
real, intent(in) :: coord(2,nn)
integer, intent(inout) :: inner_node(ni)

integer :: i, j, k, q, r, num, nei_n, count
integer :: nei(2,7), nen(7)

inner_node = 0
count = 0
do i = 1, nn
    nei_n = 0
    do j = 1, ne
        do k = 1, 4
            if (ele(k,j) == i) then
                nei_n = nei_n + 1
                if (k == 4) then
                    nei(1, nei_n) = ele(1,j)
                    nei(2, nei_n) = ele(3,j)
                else if (k == 1) then
                    nei(1, nei_n) = ele(2,j)
                    nei(2, nei_n) = ele(4,j)
                else
                    nei(1, nei_n) = ele(k+1,j)
                    nei(2, nei_n) = ele(K-1,j)
                endif
            endif
        enddo
    enddo
        
    if (sum(nei(1,1:nei_n)) == sum(nei(2,1:nei_n))) then
        count = count + 1
        inner_node(count) = i
    endif
enddo
    
end subroutine

!==============================================================================

subroutine goto_NC( coord, temp_coord, ut, vt, at, vec, error )

implicit none

integer, intent(inout) :: error
real,intent(in) ::  coord(2,4), temp_coord(2), ut(2,4) , vt(2,4) , at(2,4)
real,intent(inout) :: vec(2,3)

real :: f1 ,f2 , det_j , u , v
real :: xi, eta , norm, tol , temp
real :: jacobi(4)
integer :: i , j, err

err = 0 
xi = 0.0 ; eta = 0.0
tol = 1.0e-7
j = 0 
do 
    f1 = ( (1-xi)*(1-eta)*coord(1,1) + (1+xi)*(1-eta)*coord(1,2)+ &
    (1+xi)*(1+eta)*coord(1,3)+(1-xi)*(1+eta)*coord(1,4) )*0.25-temp_coord(1)

    f2 = ((1-xi)*(1-eta)*coord(2,1) + (1+xi)*(1-eta)*coord(2,2)+ &
    (1+xi)*(1+eta)*coord(2,3) + (1-xi)*(1+eta)*coord(2,4))*0.25-temp_coord(2)
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
    if ( norm < tol ) exit

    j = j+ 1 
    if (j > 100 ) then
        err = 1
        exit
    endif
enddo

vec = 0.0
if (err /= 1) then
    vec(1,1) = ( (1-xi)*(1-eta)*ut(1,1) + (1+xi)*(1-eta)*ut(1,2)+ &
    (1+xi)*(1+eta)*ut(1,3)+(1-xi)*(1+eta)*ut(1,4) )*0.25
    vec(2,1) = ((1-xi)*(1-eta)*ut(2,1) + (1+xi)*(1-eta)*ut(2,2)+ &
    (1+xi)*(1+eta)*ut(2,3) + (1-xi)*(1+eta)*ut(2,4))*0.25

    vec(1,2) = ( (1-xi)*(1-eta)*vt(1,1) + (1+xi)*(1-eta)*vt(1,2)+ &
    (1+xi)*(1+eta)*vt(1,3)+(1-xi)*(1+eta)*vt(1,4) )*0.25
    vec(2,2) = ((1-xi)*(1-eta)*vt(2,1) + (1+xi)*(1-eta)*vt(2,2)+ &
    (1+xi)*(1+eta)*vt(2,3) + (1-xi)*(1+eta)*vt(2,4))*0.25

    vec(1,3) = ( (1-xi)*(1-eta)*at(1,1) + (1+xi)*(1-eta)*at(1,2)+ &
    (1+xi)*(1+eta)*at(1,3)+(1-xi)*(1+eta)*at(1,4) )*0.25
    vec(2,3) = ((1-xi)*(1-eta)*at(2,1) + (1+xi)*(1-eta)*at(2,2)+ &
    (1+xi)*(1+eta)*at(2,3) + (1-xi)*(1+eta)*at(2,4))*0.25
endif

end subroutine


subroutine goto_NC_3D(coord, temp_coord, ut, vt, at, vec, err )

implicit none

integer, intent(inout) :: err
real,intent(in) :: coord(3,8), temp_coord(3), ut(3,8) , vt(3,8) , at(3,8)
real,intent(inout) :: vec(3,3)

real :: f1 ,f2, f3, det_j , u , v, w
real :: xi, eta, zeta, norm, tol, temp, coe
real :: jacobi(3,3), jacobi_inv(3,3), val(8,3)
integer :: i, j, k, errorflag, iter

val(:,1) = (/ -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0 /)
val(:,2) = (/ -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0 /)
val(:,3) = (/ -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 /)

err = 0 
xi = 0.0 ; eta = 0.0; zeta = 0.0
tol = 1.0e-4
iter = 0 
do
    f1 = ((1.0-xi)*(1.0-eta)*(1.0-zeta)*coord(1,1) + (1.0+xi)*(1.0-eta)*(1.0-zeta)*coord(1,2) + &
          (1.0+xi)*(1.0+eta)*(1.0-zeta)*coord(1,3) + (1.0-xi)*(1.0+eta)*(1.0-zeta)*coord(1,4) + &
          (1.0-xi)*(1.0-eta)*(1.0+zeta)*coord(1,5) + (1.0+xi)*(1.0-eta)*(1.0+zeta)*coord(1,6) + &
          (1.0+xi)*(1.0+eta)*(1.0+zeta)*coord(1,7) + (1.0-xi)*(1.0+eta)*(1.0+zeta)*coord(1,8))  &
          *0.125-temp_coord(1)

    f2 = ( (1.0-xi)*(1.0-eta)*(1.0-zeta)*coord(2,1) + (1+xi)*(1.0-eta)*(1.0-zeta)*coord(2,2) + &
           (1.0+xi)*(1.0+eta)*(1.0-zeta)*coord(2,3) + (1.0-xi)*(1.0+eta)*(1.0-zeta)*coord(2,4) + &
           (1.0-xi)*(1.0-eta)*(1.0+zeta)*coord(2,5) + (1.0+xi)*(1.0-eta)*(1.0+zeta)*coord(2,6) + &
           (1.0+xi)*(1.0+eta)*(1.0+zeta)*coord(2,7) + (1.0-xi)*(1.0+eta)*(1.0+zeta)*coord(2,8))  &
           *0.125-temp_coord(2)
        
    f3 = ( (1.0-xi)*(1.0-eta)*(1.0-zeta)*coord(3,1) + (1.0+xi)*(1.0-eta)*(1.0-zeta)*coord(3,2) + &
           (1.0+xi)*(1.0+eta)*(1.0-zeta)*coord(3,3) + (1.0-xi)*(1.0+eta)*(1.0-zeta)*coord(3,4) + &
           (1.0-xi)*(1.0-eta)*(1.0+zeta)*coord(3,5) + (1.0+xi)*(1.0-eta)*(1.0+zeta)*coord(3,6) + &
           (1.0+xi)*(1.0+eta)*(1.0+zeta)*coord(3,7) + (1.0-xi)*(1.0+eta)*(1.0+zeta)*coord(3,8))  &
           *0.125-temp_coord(3)

    jacobi = 0.0 
    do i = 1, 3
        do j = 1, 3
            do k = 1, 8
                coe = val(k,i)
                if (i == 1) then
                    jacobi(j,i) = jacobi(j,i)+coe*(1.0+val(k,2)*eta)*(1.0+val(k,3)*zeta)*coord(j,k)
                elseif (i == 2) then
                    jacobi(j,i) = jacobi(j,i)+coe*(1.0+val(k,1)*xi)*(1.0+val(k,3)*zeta)*coord(j,k)
                else
                    jacobi(j,i) = jacobi(j,i)+coe*(1.0+val(k,1)*xi)*(1.0+val(k,2)*eta)*coord(j,k)
                endif
            enddo
        enddo
    enddo
    jacobi = jacobi * 0.125
    jacobi_inv = 0.0
    !call inverse_mat(jacobi, jacobi_inv, 3)
    !jacobi_inv = jacobi
    !call invert(jacobi_inv,1,3)
    call FINDInv(jacobi, jacobi_inv, 3, errorflag)
    
    if (errorflag < 0 .or. iter > 20) then
        err = 1
        exit
    endif

    u = jacobi_inv(1,1) * f1 + jacobi_inv(1,2) * f2 + jacobi_inv(1,3) * f3
    v = jacobi_inv(2,1) * f1 + jacobi_inv(2,2) * f2 + jacobi_inv(2,3) * f3
    w = jacobi_inv(3,1) * f1 + jacobi_inv(3,2) * f2 + jacobi_inv(3,3) * f3
    xi = xi - u
    eta = eta - v
    zeta = zeta - w

    norm = sqrt( u**2.0 + v**2.0 + w**2.0)
    if ( norm < tol ) exit

    iter = iter + 1 
enddo

vec = 0.0
if (err /= 1) then
    do i = 1, 8
        do j = 1, 3
            vec(j,1) = vec(j,1)+((1+val(i,1)*xi)*(1+val(i,2)*eta)*(1+val(i,3)*zeta))*ut(j,i)
            vec(j,2) = vec(j,2)+((1+val(i,1)*xi)*(1+val(i,2)*eta)*(1+val(i,3)*zeta))*vt(j,i)
            vec(j,3) = vec(j,3)+((1+val(i,1)*xi)*(1+val(i,2)*eta)*(1+val(i,3)*zeta))*at(j,i)
        enddo
    enddo
    vec = vec * 0.125
endif
       
end subroutine goto_NC_3D


subroutine get_quad_coord( coord )

implicit none

real, intent(inout) ::  coord(2,4)
real :: quad(2,4)
real :: point 
real :: position(2,4)
real :: xi , eta 
integer :: i 

point = 1.0 / sqrt(3.0) 
position(1, :) = (/ -point, point, point ,-point /)
position(2, :) = (/ -point, -point, point ,point /)

do i=1,4
        xi = position(1, i)
        eta = position(2, i)

        quad(1,i) = ( (1-xi)*(1-eta)*coord(1,1) + (1+xi)*(1-eta)*coord(1,2)+ &
       (1+xi)*(1+eta)*coord(1,3)+(1-xi)*(1+eta)*coord(1,4) )*0.25
        quad(2,i) = ((1-xi)*(1-eta)*coord(2,1) + (1+xi)*(1-eta)*coord(2,2)+ &
       (1+xi)*(1+eta)*coord(2,3) + (1-xi)*(1+eta)*coord(2,4))*0.25
enddo

coord = quad

end subroutine 



subroutine set_local_rotation(ori_nn, ori_coord, new_nn, new_coord, ori_element, new_element, this, temp_element) 

implicit none

integer, intent(in) :: ori_nn, new_nn, ori_element(4), new_element(4)
real, intent(in) :: ori_coord(2,ori_nn), new_coord(2,new_nn)
integer, intent(out) :: this(4), temp_element(4)

integer :: i, temp_i
real :: dis, min_dis

temp_element = 0
min_dis = 10e8
do i = 1, 4
    call calc_len(ori_coord(:,ori_element(1)), new_coord(:,new_element(i)), dis)
    if (min_dis > dis) then
        temp_i = i
        min_dis = dis
    endif
enddo

if (temp_i == 1) then
    this = (/ 1, 2, 3, 4 /)
elseif (temp_i == 2) then
    this = (/ 2, 3, 4, 1 /)
elseif (temp_i == 3) then
    this = (/ 3, 4, 1, 2  /)
elseif (temp_i == 4) then
    this = (/ 4, 1, 2, 3 /)
endif

do i = 1, 4
    temp_element(i) = new_element(this(i))
enddo
end subroutine set_local_rotation


subroutine cholesky_ale(A, b, ndim, xdim)


! Solve Ax = b  using LAPACK
!    
! ndim: dimension of A
! LAPACK: DPOTRF, DPOTRS

integer, intent(in) :: ndim , xdim
real :: A(ndim,ndim), b(ndim, xdim)
real :: x(ndim, xdim )
integer :: info


call DPOTRF('U', ndim, A, ndim, info)
call DPOTRS('U', ndim, xdim, A, ndim, b, ndim, info)
if (info < 0) stop 'cholesky: LAPACK DPOTRS error- illegal argument value'
        
end subroutine cholesky_ale


subroutine theta_interpolation(nn, ne, coord, elem, theta, node_check, ut, flag)

implicit none

integer, intent(in) :: nn, ne, elem(4,ne)
real, intent(in) :: coord(2,nn), theta(ne)
integer, intent(inout) :: node_check(nn), flag
real, intent(inout) :: ut(2,nn)

integer :: i, j, k, q, uxi, cp, ce, iter, dof, e, this_num, num, conn_num, tn, elem_num
integer :: nodes(nn), elem_chk(ne), conn_ele(6,nn), conn_ele_num(nn), this(6), conn_elem(6), iter_node_check(nn)
real :: g, xsj, ur, ratio, vol_sum, vec_dis, sum_rate, sum_theta(2), temp_val, tol
real :: sg(3,4), elem_coord(2,4), lr(4), lz(4), elem_ut(2,4), ori_ut(2), interp_ut(2), cri_dis(nn), coe(6)
real :: shp(3,4,4), detf_a(2), detf_b(2,2), poly_coe(3,5), dvol(4), val(2,6), dis(2), pre_theta(ne), aft_theta(ne)
real :: elem_vol(6), f(2,2), detf(4), vec(2), pre_ut(2,nn), conn_ele_theta(6), conn_ele_rate(6), p(2), a(2)
real :: mat_a(2,2), vec_b(2), det_a, zero_tolerance
logical :: check

zero_tolerance = 10e-8

call calc_ori_theta(nn, ne, coord, elem, ut, pre_theta)
!open (Unit=20, File='interp_pre.plt', STATUS='replace', ACTION='write')
!Write (20,'(A,A,A)') 'TITLE="interpolation previous"'
!Write (20,*) 'VARIABLES="x", "y", "ux", "uy"'
!Write (20,'(A,I5,A,I5, A)') 'ZONE N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
!do i = 1, nn
!    write (20,*) coord(:,i), ut(:,i)
!enddo
!do i = 1, ne
!    write (20,*) elem(:,i)
!enddo
!close(20)

lr = (/ -1.0,1.0,1.0,-1.0 /)
lz = (/-1.0,-1.0,1.0,1.0 /)
g = 1.0/sqrt(3.0)		! sqrt(1/3)
do i = 1, 4
    sg(:,i) = (/ g*lr(i), g*lz(i), 1.d0 /)
end do

conn_ele_num = 0
conn_ele = 0
do i = 1, ne
    do j = 1, 4
        cp = elem(j,i)
        conn_ele_num(cp) = conn_ele_num(cp) + 1
        conn_ele(conn_ele_num(cp),cp) = i
    enddo
enddo

pre_ut = ut
do i = 1, nn
    if (node_check(i) == 1) then
        cri_dis(i) = sqrt(ut(1,i)**2.0 + ut(2,i)**2.0) * 2.0
    endif
enddo

do iter = 1, 1
    !pre_ut = ut
    !sub_iter: do 
    elem_chk = 0
    do e = 1, ne
        iter_node_check = 1
        do i = 1, ne
            if (elem_chk(i) == 0) then
                num = sum(node_check(elem(:,i)))
                if (num == 0) then
                    elem_chk(i) = 1
                elseif (num == 1) then
                    conn_num = 1
                    conn_elem(1) = i
                    do j = 1, 4
                        if (node_check(elem(j,i)) == 1) then
                            tn = elem(j,i)
                            exit
                        endif
                    enddo
                    
                    elem_num = 1
                    do j = 1, ne
                        if (j /= i) then
                            num = sum(node_check(elem(:,j)))
                            do k = 1, 4
                                if (elem(k,j) == tn) then
                                    elem_num = elem_num + 1
                                    if (num == 1) then
                                        conn_num = conn_num + 1
                                        conn_elem(conn_num) = j
                                    endif
                                    exit  !(k)
                                endif
                            enddo  !(k)
                            !if (conn_num == 2) exit  !(j)
                        endif
                    enddo  !(j)
                    val = 0.0
                    do j = 1, conn_num
                        do dof = 1, 2
                            ce = conn_elem(j)
                            elem_coord = coord(:,elem(:,ce))
                            elem_ut = ut(:,elem(:,ce))
                            !write (*,'(A,2F15.10)') 'elem_coord1:', elem_coord(:,1)
                            !write (*,'(A,2F15.10)') 'elem_coord2:', elem_coord(:,2)
                            !write (*,'(A,2F15.10)') 'elem_coord3:', elem_coord(:,3)
                            !write (*,'(A,2F15.10)') 'elem_coord4:', elem_coord(:,4)
                            !write (*,*) 
                            !write (*,'(A,2F15.10)') 'elem_ut1:', elem_ut(:,1)
                            !write (*,'(A,2F15.10)') 'elem_ut2:', elem_ut(:,2)
                            !write (*,'(A,2F15.10)') 'elem_ut3:', elem_ut(:,3)
                            !write (*,'(A,2F15.10)') 'elem_ut4:', elem_ut(:,4)
                            do k = 1, 4
                                call shapef(sg(1,k), sg(2,k), elem_coord ,shp(:,:,k), xsj, 2, .false.)
                                dvol(k) = xsj*sg(3,k)
                                if (elem(k,ce) == tn) uxi = k
                                
                                !write (*,'(A,F15.10)') 'dvol:', dvol(k)
                                !write (*,'(A,4(F12.7,1X))') 'shp1:', shp(1,1,k), shp(1,2,k), shp(1,3,k), shp(1,4,k)
                                !write (*,'(A,4(F12.7,1X))') 'shp2:', shp(2,1,k), shp(2,2,k), shp(2,3,k), shp(2,4,k)
                                !write (*,'(A,4(F12.7,1X))') 'shp3:', shp(3,1,k), shp(3,2,k), shp(3,3,k), shp(3,4,k)
                                !write (*,*) 
                            enddo
                            !write (*,*) 'uxi:', uxi
                    
                            ! calculate 1st-order polynomial coefficient
                            do k = 1, 4  
                                detf_a = 0.0
                                detf_b = 0.0
                                f = 0.0
                                do q = 1, 4 
                                    if (q == uxi) then
                                        if (dof == 1) then
                                            detf_a(1) = shp(1,q,k)
                                            detf_a(2) = shp(2,q,k)
                                            detf_b(2,1) = detf_b(2,1) + elem_ut(2,q)*shp(1,q,k)
                                            detf_b(2,2) = detf_b(2,2) + elem_ut(2,q)*shp(2,q,k)
                                        elseif (dof == 2) then
                                            detf_a(1) = shp(1,q,k)
                                            detf_a(2) = shp(2,q,k)
                                            detf_b(1,1) = detf_b(1,1) + elem_ut(1,q)*shp(1,q,k)
                                            detf_b(1,2) = detf_b(1,2) + elem_ut(1,q)*shp(2,q,k)
                                        endif
                                    else
                                        detf_b(1,1) = detf_b(1,1) + elem_ut(1,q)*shp(1,q,k)
                                        detf_b(1,2) = detf_b(1,2) + elem_ut(1,q)*shp(2,q,k)
                                        detf_b(2,1) = detf_b(2,1) + elem_ut(2,q)*shp(1,q,k)
                                        detf_b(2,2) = detf_b(2,2) + elem_ut(2,q)*shp(2,q,k)
                                    endif
                                    f(1,1) = f(1,1) + elem_ut(1,q)*shp(1,q,k)
                                    f(1,2) = f(1,2) + elem_ut(1,q)*shp(2,q,k)
                                    f(2,1) = f(2,1) + elem_ut(2,q)*shp(1,q,k)
                                    f(2,2) = f(2,2) + elem_ut(2,q)*shp(2,q,k)
                                enddo
                                detf_b(1,1) = detf_b(1,1) + 1.0
                                detf_b(2,2) = detf_b(2,2) + 1.0
                                if (dof == 1) then
                                    poly_coe(1,k) = detf_a(1)*detf_b(2,2) - detf_a(2)*detf_b(2,1)
                                elseif (dof == 2) then
                                    poly_coe(1,k) = detf_a(2)*detf_b(1,1) - detf_a(1)*detf_b(1,2)
                                endif
                                poly_coe(2,k) = detf_b(1,1)*detf_b(2,2) - detf_b(1,2)*detf_b(2,1)
                                !poly_coe(3,k) = detf_b(1,1)*detf_b(2,2) - detf_b(1,2)*detf_b(2,1)
                                f(1,1) = f(1,1) + 1.0
                                f(2,2) = f(2,2) + 1.0
                                detf(k) = f(1,1)*f(2,2) - f(1,2)*f(2,1)
                                !write (*,'(A,2F15.10)') 'f(1,:):', f(1,1), f(1,2)
                                !write (*,'(A,2F15.10)') 'f(2,:):', f(2,1), f(2,2)
                                !write (*,'(A,F15.10)') 'detf:', detf(k)
                                !write (*,*) 
                                !write (*,'(A,2F15.10)') '     detf_a:', detf_a(1:2)
                                !write (*,'(A,2F15.10)') 'detf_b(1,:):', detf_b(1,1), detf_b(1,2)
                                !write (*,'(A,2F15.10)') 'detf_b(2,:):', detf_b(2,1), detf_b(2,2)
                                !write (*,'(A,2F15.10)') '   poly_coe:', poly_coe(1:2,k)
                                !write (*,*) 
                            enddo
                            
                            ! calculate theta
                            poly_coe(:,5) = 0.0
                            do k = 1, 4
                                ratio = dvol(k)/sum(dvol)
                                poly_coe(1:2,5) = poly_coe(1:2,5) + ratio*poly_coe(1:2,k)
                                !elem_vol(j) = elem_vol(j) + abs(dvol(k)*detf(k))
                            enddo
                            poly_coe(2,5) = poly_coe(2,5) - theta(ce)  ! temporary
                            if (abs(poly_coe(1,5)) <= zero_tolerance) then
                                val(dof,j) = 0.0
                            else
                                val(dof,j) = -(poly_coe(2,5)/poly_coe(1,5))
                            endif
                        enddo
                    enddo
                    
                    if (conn_num == 1) then
                        vec = val(:,1) - ut(:,tn)
                        tol = sqrt(ut(1,tn)**2.0+ut(2,tn)**2.0)
                        tol = max(tol*10e-6,1e-10) 
                        if (abs(vec(1)) <= tol .or. abs(vec(2)) <= tol) then
                            vec = 0.0
                        else
                            dis(1) = sqrt(vec(1)**2.0+vec(2)**2.0)
                            dis(2) = abs(vec(1)*vec(2))/dis(1)
                            dis(2) = sqrt(vec(2)**2.0-dis(2)**2.0)
                            ratio = dis(2)/dis(1)
                            vec(1) = vec(1)*ratio
                            vec(2) = vec(2) * (1-ratio)
                        endif
                    elseif (conn_num == 2) then
                        temp_val = abs(val(1,1))+abs(val(2,1))+abs(val(1,2))+abs(val(2,2))
                        if (temp_val <= zero_tolerance) then
                            vec = 0.0
                        elseif (abs(val(1,1)) <= zero_tolerance .or. abs(val(1,2)) <= zero_tolerance) then
                            vec(1) = 0.0
                            vec(2) = (val(2,1)+val(2,2))*0.5
                        elseif (abs(val(2,1)) <= zero_tolerance .or. abs(val(2,2)) <= zero_tolerance) then
                            vec(1) = (val(1,1)+val(1,2))*0.5
                            vec(2) = 0.0
                        else
                            !a(1) = -(val(2,1)/val(1,1))
                            !a(2) = -(val(2,2)/val(1,2))
                            !b(1) = val(2,1)
                            !b(2) = val(2,2)
                            !x = -(b(1)-b(2))/(a(1)-a(2))
                            val(:,1) = (val(:,1) - ut(:,tn))
                            val(:,2) = (val(:,2) - ut(:,tn))
                            vec(1) = (val(2,1)-val(2,2))/(val(2,1)/val(1,1)-val(2,2)/val(1,2))
                            vec(2) = -(val(2,1)/val(1,1))*vec(1) + val(2,1)
                        endif
                        if (elem_num >= 5) then
                            vec = vec * 0.5
                        endif
                    else
                        !write (*,*) conn_num, ' :', conn_elem(1:conn_num)
                        !write (*,*) 'tn:', tn
                        coe = 0.0
                        temp_val = 0.0
                        do j = 1, conn_num
                            val(:,j) = (val(:,j) - ut(:,tn))
                            call calc_coefficient(val(:,j), coe)
                            temp_val = temp_val + abs(val(1,j)) + abs(val(2,j))
                        enddo
                        if (temp_val <= zero_tolerance) then
                            vec = 0.0
                        else
                            mat_a(1,1) = coe(1)*2.0
                            mat_a(1,2) = coe(3)
                            mat_a(2,1) = coe(3)
                            mat_a(2,2) = coe(2)*2.0
                            vec_b(1) = -coe(4)
                            vec_b(2) = -coe(5)
                            det_a = mat_a(1,1)*mat_a(2,2)-mat_a(1,2)*mat_a(2,1)
                            temp_val = mat_a(1,1)
                            WRITE (*,*) 'DET_A:', det_a
                            if (abs(det_a) <= zero_tolerance) then
                                vec = 0.0
                            else
                                mat_a(1,1) = mat_a(2,2)/det_a
                                mat_a(2,2) = temp_val/det_a
                                mat_a(1,2) = -mat_a(1,2)/det_a
                                mat_a(2,1) = -mat_a(2,1)/det_a
                                vec(1) = mat_a(1,1)*vec_b(1)+mat_a(1,2)*vec_b(2)
                                vec(2) = mat_a(2,1)*vec_b(1)+mat_a(2,2)*vec_b(2)
                                !vec = vec * 0.5
                            endif
                        endif
                    endif
                    ori_ut = ut(:,tn)
                    ut(:,tn) = ut(:,tn) + vec * 0.5
                    do j = 1, conn_num
                        elem_chk(conn_elem(j)) = 1
                    enddo
                    iter_node_check(tn) = 0
                endif
            endif
        enddo
        do i = 1, nn
            if (iter_node_check(i) == 0) node_check(i) = 0
        enddo
        if (sum(elem_chk) == ne) exit  ! (e)
    enddo  ! (e)
enddo

!do iter = 1, 1
!    !sub_iter: do 
!    elem_chk = 0
!    do i = 1, nn
!        if (node_check(i) == 1) then
!            check = .true.
!            do j = 1, conn_ele_num(i)
!                if (elem_chk(conn_ele(j,i)) == 1) then
!                    check = .false.
!                    exit
!                endif
!            enddo
!            
!            if (check) then
!                allocate (mat_A(2,conn_ele_num(i)), t_mat_A(conn_ele_num(i),2), vec_b(conn_ele_num(i)))
!                mat_A = 0.0
!                vec_b = 0.0
!                vec = 0.0
!                elem_vol = 0.0
!                do j = 1, conn_ele_num(i)
!                    ce = conn_ele(j,i)
!                    elem_coord = coord(:,elem(:,ce))
!                    elem_ut = ut(:,elem(:,ce))
!                    !if (i == 217 .and. ce == 76) then
!                    !    write (*,'(A,2F12.7)') 'elem_coord1:', elem_coord(:,1)
!                    !    write (*,'(A,2F12.7)') 'elem_coord2:', elem_coord(:,2)
!                    !    write (*,'(A,2F12.7)') 'elem_coord3:', elem_coord(:,3)
!                    !    write (*,'(A,2F12.7)') 'elem_coord4:', elem_coord(:,4)
!                    !    write (*,*) 
!                    !    write (*,'(A,2F12.7)') 'elem_ut1:', elem_ut(:,1)
!                    !    write (*,'(A,2F12.7)') 'elem_ut2:', elem_ut(:,2)
!                    !    write (*,'(A,2F12.7)') 'elem_ut3:', elem_ut(:,3)
!                    !    write (*,'(A,2F12.7)') 'elem_ut4:', elem_ut(:,4)
!                    !endif
!                    do k = 1, 4
!                        call shapef(sg(1,k), sg(2,k), elem_coord ,shp(:,:,k), xsj, 2, .false.)
!                        dvol(k) = xsj*sg(3,k)
!                        if (elem(k,ce) == i) uxi = k
!                        
!                        !if (i == 217 .and. ce == 76) then
!                        !    write (*,'(A,F15.10)') 'dvol:', dvol(j)
!                        !    write (*,'(A,4(F12.7,1X))') 'shp1:', shp(1,1,k), shp(1,2,k), shp(1,3,k), shp(1,4,k)
!                        !    write (*,'(A,4(F12.7,1X))') 'shp2:', shp(2,1,k), shp(2,2,k), shp(2,3,k), shp(2,4,k)
!                        !    write (*,'(A,4(F12.7,1X))') 'shp3:', shp(3,1,k), shp(3,2,k), shp(3,3,k), shp(3,4,k)
!                        !    write (*,*) 
!                        !endif
!                    enddo
!                    
!                    !ur = elem_ut(2,uxi)/elem_ut(1,uxi)
!                    !write (*,*) 'uxi:', uxi
!                    !write (*,*) 'ur:', ur
!           
!                    ! calculate 1st-order polynomial coefficient
!                    do k = 1, 4  
!                        detf_a = 0.0
!                        detf_b = 0.0
!                        f = 0.0
!                        do q = 1, 4 
!                            if (q == uxi) then
!                                detf_a(1) = shp(1,q,k)
!                                detf_a(2) = shp(2,q,k)
!                            else
!                                detf_b(1,1) = detf_b(1,1) + elem_ut(1,q)*shp(1,q,k)
!                                detf_b(1,2) = detf_b(1,2) + elem_ut(1,q)*shp(2,q,k)
!                                detf_b(2,1) = detf_b(2,1) + elem_ut(2,q)*shp(1,q,k)
!                                detf_b(2,2) = detf_b(2,2) + elem_ut(2,q)*shp(2,q,k)
!                            endif
!                            f(1,1) = f(1,1) + elem_ut(1,q)*shp(1,q,k)
!                            f(1,2) = f(1,2) + elem_ut(1,q)*shp(2,q,k)
!                            f(2,1) = f(2,1) + elem_ut(2,q)*shp(1,q,k)
!                            f(2,2) = f(2,2) + elem_ut(2,q)*shp(2,q,k)
!                        enddo
!                        detf_b(1,1) = detf_b(1,1) + 1.0
!                        detf_b(2,2) = detf_b(2,2) + 1.0
!                        poly_coe(1,k) = detf_a(1)*detf_b(2,2)-detf_a(2)*detf_b(2,1)
!                        poly_coe(2,k) = detf_a(2)*detf_b(1,1)-detf_a(1)*detf_b(1,2)
!                        poly_coe(3,k) = detf_b(1,1)*detf_b(2,2) - detf_b(1,2)*detf_b(2,1)
!                        f(1,1) = f(1,1) + 1.0
!                        f(2,2) = f(2,2) + 1.0
!                        detf(k) = f(1,1)*f(2,2) - f(1,2)*f(2,1)
!                    enddo
!            
!                    ! calculate theta
!                    do k = 1, 4
!                        ratio = dvol(k)/sum(dvol)
!                        mat_A(1:2,j) = mat_A(1:2,j) + ratio*poly_coe(1:2,k)
!                        vec_b(j) = vec_b(j) + ratio*poly_coe(3,k)
!                        elem_vol(j) = elem_vol(j) + abs(dvol(k)*detf(k))
!                    enddo
!                    vec_b(j) = theta(ce) - vec_b(j)
!                enddo
!                t_mat_A = transpose(mat_A)
!                val = MATMUL( t_mat_A, mat_A )
!                vec_b = MATMUL( t_mat_A, vec_b )
!                call general_lu(val, vec_b, vec, 2, 1)
!                ori_ut = ut(:,i)
!                ut(:,i) = ut(:,i) + (vec(:) - ut(:,i))/2.0
!                !ut(:,i) = vec
!                !if (i >= 300 .and. i <= 500) write (*,'(I6,3(2F11.7,1X))') i, ori_ut, ut(:,i), ori_ut-ut(:,i)
!                deallocate (mat_A, t_mat_A, vec_b)
!            endif
!        endif
!    enddo
!enddo

!open (Unit=20, File='interp_aft.plt', STATUS='replace', ACTION='write')
!Write (20,'(A,A,A)') 'TITLE="interpolation after"'
!Write (20,*) 'VARIABLES="x", "y", "ux", "uy"'
!Write (20,'(A,I5,A,I5, A)') 'ZONE N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
!do i = 1, nn
!    write (20,*) coord(:,i), ut(:,i)
!enddo
!do i = 1, ne
!    write (20,*) elem(:,i)
!enddo
!close(20)
call calc_ori_theta(nn, ne, coord, elem, ut, aft_theta)
sum_theta = 0.0 
do i = 1, ne
    sum_theta(1) = sum_theta(1) + abs(theta(i) - pre_theta(i))
    sum_theta(2) = sum_theta(2) + abs(theta(i) - aft_theta(i))
enddo

if (sum_theta(1) < sum_theta(2)) then
    flag = 2
else
    flag = 1
endif    
!write (*,'(A,3F12.7)') 'thet(pre, after):', sum_theta, flag

end subroutine theta_interpolation


subroutine calc_ori_theta(nn, ne, coord, elem, ut, theta)

implicit none 

integer, intent(in) :: nn, ne, elem(4,ne)
real, intent(in) :: coord(2,nn), ut(2,nn)
real, intent(inout) :: theta(ne)

integer :: i, j, k
real :: g, xsj, ratio
real :: sg(3,4), elem_coord(2,4), lr(4), lz(4), elem_ut(2,4)
real :: shp(3,4,4), dvol(4), f(2,2), detf(4)

lr = (/ -1.0,1.0,1.0,-1.0 /)
lz = (/-1.0,-1.0,1.0,1.0 /)
g = 1.0/sqrt(3.0)		! sqrt(1/3)
do i = 1, 4
    sg(:,i) = (/ g*lr(i), g*lz(i), 1.d0 /)
end do

theta = 0.0
do i = 1, ne
    elem_coord = coord(:,elem(:,i))
    elem_ut = ut(:,elem(:,i))
    do j = 1, 4
        call shapef(sg(1,j), sg(2,j), elem_coord ,shp(:,:,j), xsj, 2, .false.)
        dvol(j) = xsj*sg(3,j)
    enddo

    do j = 1, 4  
        f = 0.0
        do k = 1, 4 
            f(1,1) = f(1,1) + elem_ut(1,k)*shp(1,k,j)
            f(1,2) = f(1,2) + elem_ut(1,k)*shp(2,k,j)
            f(2,1) = f(2,1) + elem_ut(2,k)*shp(1,k,j)
            f(2,2) = f(2,2) + elem_ut(2,k)*shp(2,k,j)
        enddo
        f(1,1) = f(1,1) + 1.0
        f(2,2) = f(2,2) + 1.0
        detf(j) = f(1,1)*f(2,2) - f(1,2)*f(2,1)
    enddo
            
    ! calculate theta
    do j = 1, 4
        ratio = dvol(j)/sum(dvol)
        theta(i) = theta(i) + ratio*detf(j)
    enddo
enddo

end subroutine calc_ori_theta


subroutine ori_to_new_theta(ori_nn, ori_ne, ori_coord, ori_ele, ori_theta, new_nn, new_ne, new_coord, new_ele, new_theta, tol)

implicit none

integer, intent(in) :: ori_nn, ori_ne, new_nn, new_ne, ori_ele(4,ori_ne), new_ele(4,new_ne)
real, intent(in) :: tol, ori_coord(2,ori_nn), new_coord(2,new_nn), ori_theta(ori_ne)
real, intent(inout) :: new_theta(new_ne)

integer :: e, i, j, k, temp_ne, cycle_num, temp_k
integer :: temp_ele(ori_ne), temp_pn(ori_ne)
real :: ele_area, local_area, min_area, ori_ele_area
real :: ele_coord(2,4), ori_ele_coord(2,4), vec(2,3), node(2,3), ratio(2), ori_area(ori_ne)
logical :: check(2)

do e = 1, ori_ne
    ele_coord = ori_coord(:,ori_ele(:,e))
    call calc_area(4, ele_coord, ori_area(e))
enddo
    
new_theta = 0.0
cycle_num = 40
do e = 1, new_ne
    ele_coord = new_coord(:,new_ele(:,e))
    call calc_area(4, ele_coord, ele_area)
    
    temp_ne = 0;  temp_ele = 0;  temp_pn = 0
    vec(:,1) = ele_coord(:,2) - ele_coord(:,1)
    vec(:,2) = ele_coord(:,3) - ele_coord(:,4)
    do i = 1, cycle_num
        ratio(1) = (1/float(cycle_num))*0.5 + (i-1)*(1/float(cycle_num))

        node(:,1) = ele_coord(:,1) + vec(:,1)*ratio(1)
        node(:,2) = ele_coord(:,4) + vec(:,2)*ratio(1)
        vec(:,3) = node(:,2) - node(:,1)
        do j = 1, cycle_num
            ratio(2) = (1/float(cycle_num))*0.5 + (j-1)*(1/float(cycle_num))
            node(:,3) = node(:,1) + vec(:,3)*ratio(2)
        
            ! node(:,3)가 어떤 quad_element에 속해 있는지 검사
            check(1) = .FALSE.
            do k = 1, temp_ne
                ori_ele_coord = ori_coord(:,ori_ele(:,temp_ele(k)))
                call inout_check(4, ori_ele_coord, node(:,3), check(1))
                if (check(1)) then
                    temp_pn(k) = temp_pn(k) + 1
                    exit
                endif
            enddo
            
            if (check(1) == .FALSE.) then
                min_area = 10e8
                find_ori_ele: do k = 1, ori_ne
                    ori_ele_coord = ori_coord(:,ori_ele(:,k))
                    call inelement(ori_ele_coord, node(:,3), check(2), local_area, tol*1.0e-4)
                    if (check(2)) then
                        temp_ne = temp_ne + 1
                        temp_ele(temp_ne) = k
                        temp_pn(temp_ne) = 1
                        exit find_ori_ele
                    else
                        if (min_area > local_area) then
                            temp_k = k
                            min_area = local_area
                        endif
                    endif
                enddo find_ori_ele
                if (check(2) == .false.) then
                    temp_ne = temp_ne + 1
                    temp_ele(temp_ne) = temp_k
                    temp_pn(temp_ne) = 1
                endif
            endif
        enddo
    enddo
    
    do i = 1, temp_ne
        ratio(1) = float(temp_pn(i))/float(cycle_num**2)
        new_theta(e) = new_theta(e) + ori_theta(temp_ele(i))*ratio(1)
        !write (*,*) 'ele_num, ratio: ', temp_ele(i), ratio(1)
    enddo
    !do i = 1, temp_ne
    !    ratio(1) = (float(temp_pn(i))/float(cycle_num**2))*(ele_area/ori_area(temp_ele(i)))
    !    new_theta(e) = new_theta(e) + ori_theta(temp_ele(i))*ratio(1)
    !    !write (*,*) 'ele_num, ratio: ', temp_ele(i), ratio(1)
    !enddo
    !if (e < 200) write (*,*) 'ele_num:', e, new_theta(e)
enddo


end subroutine ori_to_new_theta


subroutine calc_coefficient(val, coe)

implicit none

real, intent(   in) :: val(2)
real, intent(inout) :: coe(6)

real :: a, b, c, d

if (val(1) < 10e-8) then
    coe(2) = coe(2) + 1.0
    coe(5) = coe(5) - 2.0*val(2)
    coe(6) = coe(6) + val(2)**2.0
else
    a = -(val(2)/val(1))
    b = -1.0
    c = val(2)
    d = a**2.0 + b**2.0

    coe(1) = coe(1) + a**2.0/d    ! x^2
    coe(2) = coe(2) + b**2.0/d    ! y^2
    coe(3) = coe(3) + 2.0*a*b/d   ! xy
    coe(4) = coe(4) + 2.0*a*c/d   ! x
    coe(5) = coe(5) + 2.0*b*c/d   ! y
    coe(6) = coe(6) + c**2.0/d    ! 1
endif

end subroutine calc_coefficient


  subroutine inverse_mat(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer :: n
real :: a(n,n), c(n,n)
real :: L(n,n), U(n,n), b(n), d(n), x(n)
real :: coeff
integer :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do


! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do

end subroutine inverse_mat

  
subroutine elem_to_quadelem(dim, ele_node, nn, ne, coord, elem, quad_nn, quad_ne, quad_node, quad_ele)

implicit none

integer, intent(in) :: dim, ele_node, nn, ne, elem(ele_node,ne), quad_nn, quad_ne
real, intent(in) :: coord(dim,nn)
integer, intent(inout) :: quad_ele(ele_node,quad_ne)
real, intent(inout) :: quad_node(dim,quad_nn)

integer :: i, j, k, size, face_nn(4,6)

face_nn(1,:) = (/ 1, 5, 1, 2, 3, 1 /)
face_nn(2,:) = (/ 2, 8, 5, 6, 7, 4 /)
face_nn(3,:) = (/ 3, 7, 6, 7, 8, 8 /)
face_nn(4,:) = (/ 4, 6, 2, 3, 4, 5 /)

if (dim == 2) then
    size = 9
elseif (dim == 3) then
    size = 27
endif

do i = 1, ne 
    quad_node(:,(i-1)*size+1:(i-1)*size+ele_node) = coord(:,elem(:,i))
        
    if (dim == 2) then
        do j = 1, 3
            quad_node(:,(i-1)*9+4+j) = (coord(:,elem(j,i))+coord(:,elem(j+1,i)))*0.5
        enddo
        quad_node(:,(i-1)*9+8) = (coord(:,elem(4,i))+coord(:,elem(1,i)))*0.5
        quad_node(:,(i-1)*9+9) = (quad_node(:,(i-1)*9+5)+quad_node(:,(i-1)*9+7))*0.5
    
        quad_ele(:,(i-1)*4+1) = (/ (i-1)*9+1, (i-1)*9+5, (i-1)*9+9, (i-1)*9+8 /)
        quad_ele(:,(i-1)*4+2) = (/ (i-1)*9+5, (i-1)*9+2, (i-1)*9+6, (i-1)*9+9 /)
        quad_ele(:,(i-1)*4+3) = (/ (i-1)*9+6, (i-1)*9+3, (i-1)*9+7, (i-1)*9+9 /)
        quad_ele(:,(i-1)*4+4) = (/ (i-1)*9+7, (i-1)*9+4, (i-1)*9+8, (i-1)*9+9 /)
    elseif (dim == 3) then
        do j = 1, 3
            quad_node(:,(i-1)*size+8+j) = (coord(:,elem(j,i))+coord(:,elem(j+1,i)))*0.5
            quad_node(:,(i-1)*size+12+j) = (coord(:,elem(j+4,i))+coord(:,elem(j+5,i)))*0.5
            quad_node(:,(i-1)*size+16+j) = (coord(:,elem(j,i))+coord(:,elem(j+4,i)))*0.5
        enddo
        quad_node(:,(i-1)*size+12) = (coord(:,elem(4,i))+coord(:,elem(1,i)))*0.5
        quad_node(:,(i-1)*size+16) = (coord(:,elem(8,i))+coord(:,elem(5,i)))*0.5
        quad_node(:,(i-1)*size+20) = (coord(:,elem(4,i))+coord(:,elem(8,i)))*0.5
        quad_node(:,(i-1)*size+21) = 0.0
        do j = 1, 8
            quad_node(:,(i-1)*size+21) = quad_node(:,(i-1)*size+21) + coord(:,elem(j,i))*0.125
        enddo
        
        do j = 1, 6
            quad_node(:,(i-1)*size+21+j) = 0.0
            do k = 1, 4
                quad_node(:,(i-1)*size+21+j) = quad_node(:,(i-1)*size+21+j) + coord(:,elem(face_nn(k,j),i))*0.25
            enddo
        enddo
        
        quad_ele(:,(i-1)*8+1) = (/ (i-1)*27+12, (i-1)*27+1 , (i-1)*27+9 , (i-1)*27+22, &
                                   (i-1)*27+27, (i-1)*27+17, (i-1)*27+24, (i-1)*27+21 /)
        quad_ele(:,(i-1)*8+2) = (/ (i-1)*27+9 , (i-1)*27+2 , (i-1)*27+10, (i-1)*27+22, &
                                   (i-1)*27+24, (i-1)*27+18, (i-1)*27+25, (i-1)*27+21 /)
        quad_ele(:,(i-1)*8+3) = (/ (i-1)*27+10, (i-1)*27+3 , (i-1)*27+11, (i-1)*27+22, &
                                   (i-1)*27+25, (i-1)*27+19, (i-1)*27+26, (i-1)*27+21 /)
        quad_ele(:,(i-1)*8+4) = (/ (i-1)*27+11, (i-1)*27+4 , (i-1)*27+12, (i-1)*27+22, &
                                   (i-1)*27+26, (i-1)*27+20, (i-1)*27+27, (i-1)*27+21 /)
        do j = 1, 4
            quad_ele(1:4,(i-1)*8+j+4) = quad_ele(5:8,(i-1)*8+j)
        enddo
        quad_ele(5:8,(i-1)*8+5) = (/ (i-1)*27+16, (i-1)*27+5 , (i-1)*27+13, (i-1)*27+23 /)
        quad_ele(5:8,(i-1)*8+6) = (/ (i-1)*27+13, (i-1)*27+6 , (i-1)*27+14, (i-1)*27+23 /)
        quad_ele(5:8,(i-1)*8+7) = (/ (i-1)*27+14, (i-1)*27+7 , (i-1)*27+15, (i-1)*27+23 /)
        quad_ele(5:8,(i-1)*8+8) = (/ (i-1)*27+15, (i-1)*27+8 , (i-1)*27+16, (i-1)*27+23 /)
    endif
enddo

end subroutine elem_to_quadelem


SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	REAL :: m
	REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END SUBROUTINE FINDinv