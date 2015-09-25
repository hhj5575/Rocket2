subroutine point_summation(fan, thin, face_cn, face_cp, face_n, face_p, face_nn, face_tn)

implicit none

integer, intent(inout) :: fan, face_n(fan), face_nn(fan,2000), face_cn(fan), face_tn, thin(fan)
real, intent(inout) :: face_cp(2,500), face_p(2,2000)

integer :: fan_i, i, j, k, count, sn, num, sp,lp, tn, face_count
integer :: max_i, max_val, temp_i(2), cf, max_line, max_sp, max_lp, pct(2)
integer :: pre_bs, post_bs, p_case, temp_k(3), this(fan), temp_thin(fan), temp(2)
real :: cri, max_dis, line_dis, seg
real :: temp_p(2,2000), temp_n(fan), p(2,2), dis(2), vec(2)
logical :: p_pre, p_post, post_last

cri = 0.000001
sn = 0
do i = 1, fan
	!write (*,*) face_n(i)
	!do j = 1, face_n(i)
	!	write (*,*) face_p(sn+j,:)
	!enddo
	!write (*,*) face_cn(i)
	!do j = 1, face_cn(i)
	!	write (*,*) face_cp(sn+j,:)
	!enddo
	sn = sn + face_cn(i)
enddo

temp_thin = thin
do i = 1, fan
	max_val = temp_thin(1)
	max_i = 1
	do j = 1, fan
		if ( max_val < temp_thin(j) ) then
			max_val = temp_thin(j)
			max_i = j
		endif
	enddo
	this(i) = max_i
	thin(i) = temp_thin(max_i)
	temp_thin(max_i) = -1
enddo

write (*,*) thin(:)
write (*,*) this(:)

count = face_n(this(1))
temp_p = face_p
temp_n = face_n
face_p = 0.0
face_n = 0
face_n(1) = count
tn = sum(temp_n(1:this(1)-1))
do i = 1, count
	face_nn(1,i) = i
	face_p(:,i) = temp_p(:,tn+i)
enddo

do fan_i = 2, fan
	cf = this(fan_i)
	face_count = 0
	do i = 1, face_cn(cf)
		sn = sum(face_cn(1:cf-1))
		p(:,1) = face_cp(:,sn+i)
		if ( i == face_cn(cf) ) then
			p(:,2) = face_cp(:,sn+1)
		else
			p(:,2) = face_cp(:,sn+i+1)
		endif

		p_pre = .FALSE.;  p_post = .FALSE.;  temp_k = 0
		do j = 1, fan_i-1
			temp_i = 0
			do k = 1, face_n(j)
				num = face_nn(j,k)
				call calc_len(p(:,1), face_p(:,num), dis(1))
				call calc_len(p(:,2), face_p(:,num), dis(2))
				!if ( (p(1,1) == face_p(num,1)) .AND. (p(1,2) == face_p(num,2)) ) then
				if ( dis(1) <= cri ) then
					p_pre = .TRUE.
					pre_bs = num
					temp_i(1) = j
					temp_k(1) = k
					if ( p_post  == .TRUE. .AND. temp_i(1) == temp_i(2) ) then
						pct(1) = num
						pct(2) = post_bs
						temp_k(3) = abs(temp_k(1)-temp_k(2))
					endif
					!pre_first = .FALSE.
					!if ( k == 1 ) pre_first = .TRUE.
				!elseif ( (p(2,1) == face_p(num,1)) .AND. (p(2,2) == face_p(num,2)) ) then
				elseif ( dis(2) <= cri ) then
					p_post = .TRUE.
					temp_i(2) = j
					temp_k(2) = k
					if ( p_pre == .TRUE. .AND. temp_i(1) == temp_i(2) ) then 
						pct(1) = pre_bs
						if ( k == face_n(j) ) then
							post_last = .TRUE.
							pct(2) = num + 1
						else 
							post_last = .FALSE.
							pct(2) = face_nn(j,k+1)
						endif
						if ( temp_k(1) == 1 ) then
							temp_k(3) = face_n(j)+1 - temp_k(2)
						else
							temp_k(3) = abs(temp_k(1)-temp_k(2))
						endif
					endif
					if ( k == face_n(j) ) then
						post_last = .TRUE.
						post_bs = num + 1
					else 
						post_last = .FALSE.
						post_bs = face_nn(j,k+1)
					endif
				endif
			enddo
		enddo
	
		!if ( i == 1 ) then
		!	p_first = .FALSE.
		!	if (p_pre) p_first = .TRUE.
		!endif
			
		if ( p_pre == .FALSE. ) then 
			p_case = 1
		elseif ( (p_pre) .AND. (p_post == .FALSE.) ) then
			p_case = 2
		elseif ( (p_pre) .AND. (p_post) ) then
			p_case = 3
		endif
		
		!if ( i == face_cn(fan_i) ) then
		!	if ( (p_case == 3) .AND. (p_first == .FALSE. ) ) p_case = 2
		!endif

		tn = sum(temp_n(1:cf-1))
		num = 0
		sp = 0;  lp = 0
		
		do j = 1, temp_n(cf)
			call calc_len(p(:,1), temp_p(:,tn+j), dis(1))
			call calc_len(p(:,2), temp_p(:,tn+j), dis(2))
			if ( dis(1) <= cri ) then
				sp = tn+j
				num = num + 1
			elseif ( dis(2) <= cri ) then
				lp = tn+j
				if ( i == face_cn(cf) ) lp = tn+temp_n(cf)+1
				num = num + 1
			endif
			if ( num == 2 ) exit
		enddo
		if ( sp == 0 .OR. lp == 0 ) then
			write (*,*) 'Match Error!!'
			stop
		endif
		!write (*,*) 
		!write (*,*) face_n(1:fan_i-1)
		!write (*,*) p(1,:)
		!write (*,*) p(2,:)
		!write (*,*) p_case, sp, lp
		if ( p_case == 1 ) then
			do j = sp, lp-1
				count = count + 1
				face_count = face_count + 1
				face_nn(fan_i,face_count) = count
				face_p(:,count) = temp_p(:,j)
				if ( j == sp ) then
					temp(1) = face_count
				elseif ( j == lp-1 ) then
					temp(2) = face_count
				endif
			enddo
		elseif ( p_case == 2 ) then ! 시작점 겹치고 끝점 겹치지 않음
			face_count = face_count + 1
			face_nn(fan_i,face_count) = pre_bs
			temp(1) = face_count
			do j = sp+1,lp-1
				count = count + 1
				face_count = face_count + 1
				face_nn(fan_i,face_count) = count
				face_p(:,count) = temp_p(:,j)
				if ( j == lp-1 ) temp(2) = face_count
			enddo
		elseif ( p_case == 3 ) then ! 시작점, 끝점 둘 다 겹침
			!if ( temp_k(3) < abs(lp-sp) ) then
				!do j = count, post_bs+num-1, -1
				!	face_p(j+num,:) = face_p(j,:)
				!enddo
				!do j = 1, num
				!	face_p(post_bs+j-1,:) = temp_p(lp-j,:)
				!enddo
				!do j = face_n(temp_i), temp_k(2)+num,-1
				!	face_nn(temp_i,j+num) = face_nn(temp_i,j)
				!enddo
				!do j = 1, face_n(temp_i)
				!	if ( face_nn(temp_i,j) >= post_bs+num-1 ) then
				!		face_nn(temp_i,j) = face_nn(temp_i,j) + 1
				!	endif
				!enddo
				!do j = 1, num
				!	face_nn(temp_i,temp_k(2)+j) = post_bs-1+j
				!enddo
				!face_n(temp_i) = face_n(temp_i) + num
				
				!count = count + num
				!if ( post_bs <= pre_bs ) pre_bs = pre_bs + num
			!endif
			!if ( pre_first ) then
			if ( temp_k(3) == 0 ) then
				face_count = face_count + 1
				face_nn(fan_i,face_count) = pre_bs
				num = abs(lp-sp)
				temp(1) = face_count
				do j = sp+1,lp-1
					count = count + 1
					face_count = face_count + 1
					face_nn(fan_i,face_count) = count
					face_p(:,count) = temp_p(:,j)
					if ( j == lp-1 ) temp(2) = face_count
				enddo
			else
				face_count = face_count + 1
				face_nn(fan_i,face_count) = pct(1)
				num = temp_k(3)
				do j = pct(2)+num-2, pct(2),-1
					face_count = face_count + 1
					face_nn(fan_i,face_count) = j
				enddo
			endif
		endif
		call calc_len(p(:,1), p(:,2), line_dis)

		if ( i == 1 ) then
			max_line = 1
			max_dis = line_dis
			max_sp = temp(1);  max_lp = temp(2)
		else
			if ( max_dis < line_dis .AND. p_case /= 3 ) then
				max_line = i
				max_dis = line_dis
				max_sp = temp(1);  max_lp = temp(2)
			endif
		endif
	enddo
	
	if ( 2 * (face_count/2) /= face_count ) then 
		num = max_lp-max_sp+1
		seg = max_dis / (num+1)
		p(:,1) = face_p(:,face_nn(fan_i,max_sp))
		p(:,2) = face_p(:,face_nn(fan_i,max_lp))
		call calc_len( p(:,1), p(:,2), line_dis )
		vec(:) = (p(:,2)-p(:,1))/line_dis
		do j = count, face_nn(fan_i,max_lp)+1, -1
			face_p(:,j+1) = face_p(:,j)
		enddo
		do j = 1, num+1
			face_p(:,face_nn(fan_i,max_sp)+j) = p(:,1) + (seg*vec(:)*j)
		enddo
		count = count + 1
		do j = face_count, max_lp+1,-1
			face_nn(fan_i,j+1) = face_nn(fan_i,j)
		enddo
		face_count = face_count + 1
		do j = 1, face_count
			if ( face_nn(fan_i,j) >= face_nn(fan_i,max_lp)+1 ) face_nn(fan_i,j) = face_nn(fan_i,j)+1				
		enddo
		face_nn(fan_i,max_lp+1) = face_nn(fan_i,max_lp)+1
	endif
	face_n(fan_i) = face_count
enddo
face_tn = count

end subroutine
