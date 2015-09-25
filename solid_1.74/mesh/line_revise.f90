subroutine line_revise(on, op, opn, d)

implicit none

integer :: on, opn(on)
real :: op(on*500,5), d

integer :: i, j, k, q, num, op_num, count, max_opn, m
integer :: temp_opn(on)
real :: length, t, tt, t_t, tension, coe, alpha, sum_angle
real :: temp_op(on*500,5), h(6), vp(4,2), o(2)
integer, allocatable :: this(:)
real, allocatable :: p(:,:)

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375
!--------------------------------------------------------------------------
temp_opn = opn;   temp_op = op
opn = 0;   op = 0;   num = 0;   op_num = 0;   tension = 0.6

max_opn = temp_opn(1)
do i = 2, on
	if ( max_opn < temp_opn(i) ) max_opn = temp_opn(i)
enddo	
allocate( this(max_opn), p(max_opn,2) )
this = 0;   count = 0

do i = 1, on
	do j = 1, temp_opn(i)
		num = num + 1
		if ( num > sum(temp_opn(1:i)) ) then
			num = num - 1
			exit
		endif

		if ( int(temp_op(num,5)) == 1 ) then
			op_num = op_num + 1
			op(op_num,1) = i
			op(op_num,2) = op_num - sum(opn(1:i-1))
			op(op_num,3:4) = temp_op(num,3:4)
		elseif ( int(temp_op(num,5)) == 2 ) then
			count = 1
			this(count) = num
			p(count,1:2) = temp_op(num,3:4)
			do 
				if ( int(temp_op(num+count,5)) /= 2 ) then
					this(count+1) = num+count
					count = count + 1
					p(count,1:2) = temp_op(this(count),3:4)
					exit
				else
					this(count+1) = num+count
					if ( num+count == temp_opn(i) ) then
						this(count+1) = sum(temp_opn(1:i-1)) + 1
						count = count + 1
						p(count,1:2) = temp_op(this(count),3:4)
						exit
					endif
					count = count + 1
					p(count,1:2) = temp_op(this(count),3:4)
				endif
			enddo

			do k = 1, count - 1
				call calc_len(p(k,:), p(k+1,:), length)
				if ( length > 1.3 * d ) then
					m = 2 * nint(length / (2*d))
				else
					m = 0
				endif	
			
				do q = 0, m-1
					op_num = op_num + 1
					t = m
					t = (1/t) * q
				
					tt = t*t
					t_t = 1-t
					h(1) = (1+2*t) * t_t * t_t
					h(2) = t * t_t * t_t
					h(3) = tt * (3-2*t)
					h(4) = tt * (t-1)

					if (k == 1) then
						vp(1,:) = tension * (p(k+1,:) - p(k,:))
						vp(2,:) = tension * (p(k+2,:) - p(k,:))
						op(op_num,1) = i
						op(op_num,2) = op_num - sum(opn(1:i-1)) 
						op(op_num,3:4) = h(1)*p(k,:) + h(2)*vp(1,:) + h(3)*p(k+1,:) + h(4)*vp(2,:)
					elseif (k < count - 1) then
						vp(2,:) = tension * (p(k+1,:) - p(k-1,:))
						vp(3,:) = tension * (p(k+2,:) - p(k,:))
						op(op_num,1) = i
						op(op_num,2) = op_num - sum(opn(1:i-1)) 
						op(op_num,3:4) = h(1)*p(k,:) + h(2)*vp(2,:) + h(3)*p(k+1,:) + h(4)*vp(3,:)
					elseif (k == count - 1) then
						vp(3,:) = tension * (p(count,:) - p(count-2,:))
						vp(4,:) = tension * (p(count,:) - p(count-1,:))
						op(op_num,1) = i
						op(op_num,2) = op_num - sum(opn(1:i-1)) 
						op(op_num,3:4) = h(1)*p(count-1,:) + h(2)*vp(3,:) + h(3)*p(count,:) + h(4)*vp(4,:)
					endif
				enddo
			enddo
			num = num + count - 2
			p = 0
		elseif ( int(temp_op(num,5)) == 3 ) then
			count = 3
			this(1:3) = (/num, num+1, num+2 /)
			if (this(2) == temp_opn(i)) this(3) = sum(temp_opn(1:i-1)) + 1
			p(1:3,:) = temp_op(this(1:3),3:4)
						
			h(1) = p(2,1)-p(1,1)
			h(2) = p(3,1)-p(1,1)
			h(3) = p(2,1)**2-p(1,1)**2
			h(4) = p(3,1)**2-p(1,1)**2
			h(5) = p(2,2)**2-p(1,2)**2
			h(6) = p(3,2)**2-p(1,2)**2

			o(2) = h(1)*(h(4)+h(6)) - h(2)*(h(3)+h(5))
			o(2) = o(2) / (2*((p(1,2)-p(2,2))*h(2) - (p(1,2)-p(3,2))*h(1)))

			o(1) = (h(3)+h(5)+2*o(2)*(p(1,2)-p(2,2))) / (2*h(1))

			sum_angle = 0.
			do k = 1, 2
				call calc_angle(p(k,:), o, p(k+1,:), alpha)
				sum_angle = sum_angle + alpha
			enddo
			coe = 1.
			if ( sum_angle > 360. ) coe = -1.

			do k = 1, count - 1 
				call calc_len(p(k,:), p(k+1,:), length)
				length = 1.1 * length
				if ( length > 1.3 * d ) then
					m = 2 * nint(length / (2*d))
				else
					m = 0
				endif
				
				call calc_angle(p(k,:), o, p(k+1,:), alpha)
				if (coe == -1.) alpha = 360 - alpha
				alpha = alpha / m
				vp(1,:) = p(k,:) - o(:)
				
				do q = 0, m-1
					op_num = op_num + 1
					op(op_num,1) = i
					op(op_num,2) = op_num - sum(opn(1:i-1)) 
					op(op_num,3) = o(1)+(cos(alpha*q*pi/180.)*vp(1,1) + coe*sin(alpha*q*pi/180.)*vp(1,2))
					op(op_num,4) = o(2)+((-1.)*coe*sin(alpha*q*pi/180.)*vp(1,1) + cos(alpha*q*pi/180.)*vp(1,2))
				enddo
			enddo

			num = num + 1
		endif
	enddo
	opn(i) = op_num - sum(opn(1:i-1))
enddo

num = 0 
!do i = 1, on
!	do j = 1, opn(i)
!		num = num + 1
!		write (*,*) op(num,1:4)
!	enddo
!enddo
!write (*,*)

end subroutine
