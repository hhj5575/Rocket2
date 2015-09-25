!****************************************************************************
!
!  PROGRAM : len_smooth
!
!  PURPOSE : 경계 절점의 길이 유연화
!
!****************************************************************************
Subroutine len_smooth(pn, p_num, point, nn, np)

implicit none

integer, intent(in) :: pn, nn, p_num(pn)
real, intent(inout) :: point(2,pn), np(2,nn)

integer :: i, j, k, q, o, i_l, i_r, stop_p
real :: len_sum, len_num, length, alpha, beta, len_f, len_a, cri_ang
real :: seg, seg_sum, seg_num, seg_ext
real :: vector(2)

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

cri_ang = 145.
k = 0 
do i = 1, pn
	k = k + 1
	if (k .GT. pn) then
		exit
	endif
	
	len_sum = 0.
	len_num = 0.
	q = k 
	do j = q, pn
		
		k = k + 1
		if (k == pn) then
			i_l = pn - 1	
			i_r = 1
		else
			i_l = k - 1
			i_r = k + 1 
		endif
				
		call calc_angle(point(:,i_l), point(:,k), point(:,i_r), alpha)
		call calc_len(point(:,k), point(:,i_l), length)
								
		if (k == pn) then
			len_sum = len_sum + length
			call calc_len(point(:,k), point(:,1), length)
			len_sum = len_sum + length
			len_num = len_num + 2.
			exit
		elseif ( (alpha >= cri_ang) .AND. (alpha <= 185.) ) then
			len_sum = len_sum + length
			len_num = len_num + 1.
		elseif (k - 1 == q) then
			len_sum = len_sum + length
			len_num = len_num + 1.
		else
			len_sum = len_sum + length
			len_num = len_num + 1.
			k = k - 1
			exit
		endif
		
	enddo
	seg = len_sum / len_num
	
	if ( len_num >= 5 ) then
		seg_sum = 0;  seg_num = 0; seg_ext = 0
		do j = q, k - 1
			seg_num = seg_num + 1
			call calc_len(point(:,j), point(:,j+1), length)
			seg_ext = seg_ext + (length - seg)
			seg_sum = seg_sum + length
			if ( ( seg_num /= 1 ) .AND. ( seg_ext > seg ) ) then
				k = j + 1
				seg = seg_sum / seg_num
				exit
			endif
		enddo
	endif
	
	do j = q + 1, k
		
		if (j == pn) then
			i_r = 1
			i_l = j - 1
		else
			i_r = j + 1
			i_l = j - 1
		endif
		
		call calc_len(point(:,i_l), point(:,j), len_f)
		call calc_len(point(:,j), point(:,i_r), len_a)
		call calc_angle(point(:,i_l), point(:,j), point(:,i_r), beta)
        
		if ( (seg - len_f) .LT. len_a ) then

			!if ( j == k ) then
			!	if (len_a .LE. seg) then
			!		vector(1) = ((point(i_l,2) - point(j,2)) / len_a) * (seg - len_a)
			!		vector(2) = ((point(i_l,3) - point(j,3)) / len_a) * (seg - len_a)
			!	else
			!		vector(1) = ((point(j,2) - point(i_r,2)) / len_a) * (seg - len_a)
			!		vector(2) = ((point(j,3) - point(i_r,3)) / len_a) * (seg - len_a)
			!	endif
			!	point(j,2) = point(j,2) + vector(1)
			!	point(j,3) = point(j,3) + vector(2)
			!	np(nint(point(j,1)),2) = point(j,2)
			!	np(nint(point(j,1)),3) = point(j,3)
			!else
				if (len_f .LE. seg) then
					
					vector(1) = ((point(1,i_r) - point(1,j)) / len_a) * (seg - len_f)! * cos((beta - 180) / 2. * pi / 180.)
					vector(2) = ((point(2,i_r) - point(2,j)) / len_a) * (seg - len_f)! * cos((beta - 180) / 2. * pi / 180.)

				else
					
					vector(1) = ((point(1,j) - point(1,i_l)) / len_f) * (seg - len_f)! * cos((beta - 180) / 2. * pi / 180.)
					vector(2) = ((point(2,j) - point(2,i_l)) / len_f) * (seg - len_f)! * cos((beta - 180) / 2. * pi / 180.)

				endif
				
				point(:,j) = point(:,j) + vector(:)
				np(:,p_num(j)) = point(:,j)
			!endif
			
		else 
				
			stop_p = j
			len_sum = 0.
			len_num = 0.

			do o = q+1, stop_p + 1
				call calc_len(point(:,o), point(:,o-1), length)
				len_sum = len_sum + length
				len_num = len_num + 1.
			enddo

			seg = len_sum / len_num

			do o = q + 1, stop_p

				if (o == pn) then
					i_r = 1
					i_l = o - 1
				else
					i_r = o + 1
					i_l = o - 1
				endif
				
				call calc_len(point(:,i_l), point(:,o), len_f)
				call calc_len(point(:,o), point(:,i_r), len_a)
				call calc_angle(point(:,i_l), point(:,o), point(:,i_r), beta)

				if (len_f .LE. seg) then
					
					vector(1) = ((point(1,i_r) - point(1,o)) / len_a) * (seg - len_f)! * cos((beta - 180) / 2. * pi / 180.)
					vector(2) = ((point(2,i_r) - point(2,o)) / len_a) * (seg - len_f)! * cos((beta - 180) / 2. * pi / 180.)

				else
					
					vector(1) = ((point(1,o) - point(1,i_l)) / len_f) * (seg - len_f)! * cos((beta - 180) / 2. * pi / 180.)
					vector(2) = ((point(2,o) - point(2,i_l)) / len_f) * (seg - len_f)! * cos((beta - 180) / 2. * pi / 180.)

				endif
				
				point(:,o) = point(:,o) + vector(:)
				np(:,p_num(o)) = point(:,o)
			enddo
			
			k = stop_p - 1
			exit
					
		endif
	
	enddo

enddo

end subroutine len_smooth