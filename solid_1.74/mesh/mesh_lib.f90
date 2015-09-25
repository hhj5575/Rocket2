Subroutine calc_angle(a, b, c, angle)

implicit none

Real, intent(in) :: a(2), b(2), c(2)
Real, intent(out) :: angle

Real :: a_x(4)
Real :: a_y(4)
Real :: ang(2)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

a_x(1) = a(1) - b(1)
a_y(1) = a(2) - b(2)  

if ((a_x(1) .GE. 0) .AND. (a_y(1) .EQ. 0)) then
	
	ang(1) = 0

else if ((a_x(1) .LT. 0) .AND. (a_y(1) .EQ. 0)) then

	ang(1) = 180

else if ((a_x(1) .EQ. 0) .AND. (a_y(1) .GT. 0)) then

	ang(1) = 90

else if ((a_x(1) .EQ. 0) .AND. (a_y(1) .LT. 0)) then

	ang(1) = 270

else if ((a_x(1) .GT. 0) .AND. (a_y(1) .GT. 0)) then

	ang(1) = atan(a_y(1)/a_x(1))
	ang(1) = ang(1) * (180 / pi)
	
else if ((a_x(1) .LT. 0) .AND. (a_y(1) .GT. 0)) then

	ang(1) = atan(a_y(1)/a_x(1))
	ang(1) = ang(1) * (180 / pi)
	ang(1) = 180 + ang(1)

else if ((a_x(1) .LT. 0) .AND. (a_y(1) .LT. 0)) then

	ang(1) = atan(a_y(1)/a_x(1))
	ang(1) = ang(1) * (180 / pi)
	ang(1) = 180 + ang(1)

else if ((a_x(1) .GT. 0) .AND. (a_y(1) .LT. 0)) then

	ang(1) = atan(a_y(1)/a_x(1))
	ang(1) = ang(1) * (180 / pi)
	ang(1) = 360 + ang(1)

endif

a_x(2) = c(1) - b(1)
a_y(2) = c(2) - b(2)

if ((a_x(2) .GE. 0) .AND. (a_y(2) .EQ. 0)) then
	
	ang(2) = 0

else if ((a_x(2) .LT. 0) .AND. (a_y(2) .EQ. 0)) then

	ang(2) = 180

else if ((a_x(2) .EQ. 0) .AND. (a_y(2) .GT. 0)) then

	ang(2) = 90

else if ((a_x(2) .EQ. 0) .AND. (a_y(2) .LT. 0)) then

	ang(2) = 270

else if ((a_x(2) .GT. 0) .AND. (a_y(2) .GT. 0)) then

	ang(2) = atan(a_y(2)/a_x(2))
	ang(2) = ang(2) * (180 / pi)
	
else if ((a_x(2) .LT. 0) .AND. (a_y(2) .GT. 0)) then

	ang(2) = atan(a_y(2)/a_x(2))
	ang(2) = ang(2) * (180 / pi)
	ang(2) = 180 + ang(2)

else if ((a_x(2) .LT. 0) .AND. (a_y(2) .LT. 0)) then

	ang(2) = atan(a_y(2)/a_x(2))
	ang(2) = ang(2) * (180 / pi)
	ang(2) = 180 + ang(2)

else if ((a_x(2) .GT. 0) .AND. (a_y(2) .LT. 0)) then

	ang(2) = atan(a_y(2)/a_x(2))
	ang(2) = ang(2) * (180 / pi)
	ang(2) = 360 + ang(2)

endif

angle = ang(1) - ang(2)

if (angle .LT. 0) then

	angle = 360 + angle

endif

end subroutine
!==========================================================================
Subroutine calc_length_bak(p_a, p_b, length)

implicit none

real, intent(in) :: p_a(2), p_b(2)
real, intent(out) :: length

length = sqrt( ( (p_a(1) - p_b(1)) ** 2) + ((p_a(2) - p_b(2)) ** 2 ) )

end subroutine
!==========================================================================
Subroutine calc_length_3D(p_a, p_b, length)

implicit none

real, intent(in) :: p_a(3), p_b(3)
real, intent(out) :: length

length = sqrt(((p_a(1)-p_b(1))**2.0) + ((p_a(2)-p_b(2))**2.0) + ((p_a(3)-p_b(3))**2.0))

end subroutine
!==========================================================================
subroutine calc_area(pn, coord, area)

implicit none

integer,intent(in) :: pn
real, intent(in) :: coord(2,pn)
real, intent(out) :: area

integer :: i

area = 0.0
do i = 1, pn
    if (i == pn) then
        area = area + coord(1,i)*coord(2,1) - coord(1,1)*coord(2,i)
    else
        area = area + coord(1,i)*coord(2,i+1) - coord(1,i+1)*coord(2,i)
    endif
enddo
area = abs(area) * 0.5

end subroutine calc_area
!==========================================================================
Subroutine inout_check(pn, p, cp, in_out)
	
implicit none

integer, intent(in) :: pn	
real, intent(in) :: p(2,pn), cp(2)
logical, intent(out) :: in_out

integer :: i, i_r
real :: total_ang, ang

total_ang = 0
do i = 1, pn
	i_r = i + 1
	if ( i == pn ) i_r = 1

	call calc_angle( p(:,i_r), cp, p(:,i), ang )

	if ( ang < 180. ) then
		total_ang = total_ang + ang
	else
		total_ang = total_ang + ( ang - 360. )
	endif
enddo

if ( (abs(total_ang) <= 1.) ) then
	in_out = .FALSE.
else
	in_out = .TRUE.
endif

end subroutine

subroutine projection_s_cont(point_a, point_b, test, s)

implicit none

real, intent(in) :: point_a(2), point_b(2), test(2)
real, intent(out) :: s
real :: vector1(2), vector2(2), norm

vector1 = point_b - point_a
vector2 = test - point_a
norm = DOT_PRODUCT(vector1, vector1)
if (norm == 0.0) STOP 'projection_s: zero norm denominator!'

s = DOT_PRODUCT(vector1, vector2)/norm

end subroutine projection_s_cont


subroutine calc_divided_number(length, a0, am, m, scale_factor)

implicit none

integer, intent(inout) :: m
real, intent(in) :: length, a0, am
real, intent(inout) :: scale_factor

real :: r, sn, an

if (a0 * 0.6 >= length) then
    m = 1
    scale_factor = 1.0
else
    r = ((am-a0)/length) + 1.0
    m = 1
    sn = 0.0
    find_m: do  ! find m
        an = a0*(r**(m-1))
        if (an < a0 * 0.01) then
            write (*,'(A,I4,3(A,ES15.7))') 'm: ', m, ', a0: ', a0, ', an: ', an, ', sn:', sn
            stop 'divide_point: an is small !!'
        endif
        sn = sn + an
        if (sn*1.02 >= length) then
            scale_factor = length/sn
            exit find_m
        else
            m = m + 1
        endif
    enddo find_m  ! find m
endif

!write (*,'(A,I5,A,3ES15.7)') 'm: ', m, ', Scale factor:', scale_factor, a0, am

end subroutine calc_divided_number



subroutine adjust_scale_factor(length, a0, am, m, scale_factor)

implicit none

integer, intent(in) :: m
real, intent(in) :: length, a0, am
real, intent(inout) :: scale_factor

integer :: i
real :: sn, r

sn = 0.0
r = ((am-a0)/length) + 1.0
do i = 0, m - 1
    sn = sn + a0*(r**i)
enddo
scale_factor = length/sn


end subroutine adjust_scale_factor