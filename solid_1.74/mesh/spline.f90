!****************************************************************************
!
!  PROGRAM: spline
!
!****************************************************************************
subroutine spline(pn, bs, b, count, p, check, d)

implicit none

integer :: pn, bs, count, check
real :: d
real :: b(2,pn*500), p(2,count)

integer :: i, j, m
real :: t, tt, t_t, tension, alpha, length, coe, sum_angle
real :: h(6), o(2), vp(2,4)

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

!============================================================================
!°î·üÀ» Á¤ÇÏ´Â º¯¼ö ÁöÁ¤
tension = 0.6
!do i = 1, count
!	write (*,*) p(i,:)
!enddo
!write (*,*) " "
if (count < 3) then
	write (*,*) "Not enough curve points"
	stop
elseif ( check == 2 ) then
	do i = 1, count-1
		call calc_len(p(:,i), p(:,i+1), length)
		if ( length > 1.3 * d ) then
			m = 2 * nint(length / (2*d))
		else
			m = 0
		endif

		do j = 0, m-1
			bs = bs + 1
			t = m
			t = (1/t) * j
		
			tt = t*t
			t_t = 1-t
			h(1) = (1+2*t) * t_t * t_t
			h(2) = t * t_t * t_t
			h(3) = tt * (3-2*t)
			h(4) = tt * (t-1)
						
			if (i == 1) then
				vp(:,1) = tension * (p(:,i+1) - p(:,i))
				vp(:,2) = tension * (p(:,i+2) - p(:,i))
				b(:,bs) = h(1)*p(:,i) + h(2)*vp(:,1) + h(3)*p(:,i+1) + h(4)*vp(:,2)
			elseif (i < count - 1) then
				vp(:,2) = tension * (p(:,i+1) - p(:,i-1))
				vp(:,3) = tension * (p(:,i+2) - p(:,i))
				b(:,bs) = h(1)*p(:,i) + h(2)*vp(:,2) + h(3)*p(:,i+1) + h(4)*vp(:,3)
			elseif (i == count - 1) then
				vp(:,3) = tension * (p(:,count) - p(:,count-2))
				vp(:,4) = tension * (p(:,count) - p(:,count-1))
				b(:,bs) = h(1)*p(:,count-1) + h(2)*vp(:,3) + h(3)*p(:,count) + h(4)*vp(:,4)
			endif
		enddo
	enddo
elseif ( check == 3 ) then
	if ( count /= 3 ) then
		write (*,*) "Too much curve points"
		stop
	else
		h(1) = p(1,2)-p(1,1)
		h(2) = p(1,3)-p(1,1)
		h(3) = p(1,2)**2-p(1,1)**2
		h(4) = p(1,3)**2-p(1,1)**2
		h(5) = p(2,2)**2-p(2,1)**2
		h(6) = p(2,3)**2-p(2,1)**2

		o(2) = h(1)*(h(4)+h(6)) - h(2)*(h(3)+h(5))
		o(2) = o(2) / (2*((p(2,1)-p(2,2))*h(2) - (p(2,1)-p(2,3))*h(1)))

		o(1) = (h(3)+h(5)+2*o(2)*(p(2,1)-p(2,2))) / (2*h(1))
		
		!r = sqrt( (p(1,1)-o(1))**2 + (p(1,2)-o(2))**2 )

		!write (*,*) o(1), o(2) ! r
		sum_angle = 0.
		do i = 1, 2
			call calc_angle(p(:,i), o, p(:,i+1), alpha)
			sum_angle = sum_angle + alpha
		enddo
		coe = 1.
		if ( sum_angle > 360. ) coe = -1.
			
		do i = 1, count - 1 
			call calc_len(p(:,i), p(:,i+1), length)
			length = 1.1 * length
			if ( length > 1.3 * d ) then
				m = 2 * nint(length / (2*d))
			else
				m = 0
			endif
			
			call calc_angle(p(:,i), o, p(:,i+1), alpha)
			if (coe == -1.) alpha = 360 - alpha
			alpha = alpha / m
			vp(:,1) = p(:,i) - o(:)
			
			do j = 0, m-1
				bs = bs + 1
				b(1,bs) = o(1)+(cos(alpha*j*pi/180.)*vp(1,1) + coe*sin(alpha*j*pi/180.)*vp(2,1))
				b(2,bs) = o(2)+((-1.)*coe*sin(alpha*j*pi/180.)*vp(1,1) + cos(alpha*j*pi/180.)*vp(2,1))
			enddo
		enddo
	endif
endif

end subroutine spline