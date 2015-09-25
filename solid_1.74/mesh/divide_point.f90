subroutine divide_point(pn, p, d, line, bs, b, edge)

implicit none

integer, intent(in) :: pn, line(pn)
integer, intent(inout) :: bs
real, intent(in) :: d, p(2,pn), edge(2,pn)
real, intent(inout) :: b(2,pn*500)

integer :: i, j, i_r, m, status, count, num
integer :: temp_start, temp_end, temp_m, temp_i, ex_num
integer :: this(pn)
real :: vector(2), length, temp_len, a0, sn, r, scale_factor, an, am

ex_num = pn
bs = 0
temp_len = 0.0
temp_i = 0
i = 0

do num = 1, ex_num
	i = i + 1
	if (i > ex_num) exit

	if ( line(i) == 1 ) then
		i_r = i + 1
		if (i == ex_num) i_r = 1

		length = sqrt(((p(1,i_r) - p(1,i)) ** 2) + ((p(2,i_r) - p(2,i)) ** 2))
		vector(:) = p(:,i_r) - p(:,i)
		
        a0 = edge(1,i)
        am = edge(2,i)
        if (a0 * 0.6 >= length) then
            m = 0
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
        !write (*,'(A,I5,A,I5,A,3ES15.7)') 'i: ', i, ', m: ', m, ', Scale factor:', scale_factor, a0, am
        
		if ( length > temp_len ) then
			temp_m = m
			temp_len = length
			temp_i = i
		endif
			
		if (m > 0) then
			vector(:) = vector(:) / length
            sn = 0.0
            bs = bs + 1
            b(:,bs) = p(:,i)
            if ( temp_i == i ) temp_start = bs
			do j = 0, m - 2
				bs = bs + 1
                sn = sn + a0*(r**j)
				b(:,bs) = p(:,i) + (sn * vector(:))*scale_factor
				if ( temp_i == i ) then
					if ( j == m - 2 ) temp_end = bs
				endif
			enddo
		else
			bs = bs + 1
			b(:,bs) = p(:,i)
			
			if ( temp_i == i ) then
				temp_start = bs;  temp_end = bs
			endif
		endif
	elseif ( line(i) == 2 ) then
		count = 1
		this(count) = i
		do 
			if ( line(i+count) /= 2 ) then
				this(count+1) = i+count
				count = count + 1
				exit
			else
				this(count+1) = i+count
				if ( i+count == ex_num ) then
					this(count+1) = 1
					count = count + 1
					exit
				endif
				count = count + 1
			endif
		enddo

		call spline( pn, bs, b, count, p(:,this(1:count)), line(i), d )
		i = i + count - 2
	elseif ( line(i) == 3 ) then
		count = 3
		this(1:3) = (/i, i+1, i+2 /)
		if (this(2) == ex_num) this(3) = 1
		
		call spline( pn, bs, b, count, p(:,this(1:count)), line(i), d )
		i = i + 1
	endif	
enddo

if ( 2 * (bs / 2) /= bs ) then
    write (*,*) "bs : ", bs, temp_start, temp_i, temp_m
    
    do i = bs, temp_end + 1, -1
	    b(:,i+1) = b(:,i)
    enddo
    bs = bs + 1
    
    i = temp_i
    i_r = i + 1
    if (i == ex_num) i_r = 1
    
    a0 = edge(1,i)
    am = edge(2,i)
    m = temp_m + 1
    length = sqrt(((p(1,i_r) - p(1,i)) ** 2) + ((p(2,i_r) - p(2,i)) ** 2))
    r = ((am-a0)/length) + 1.0
    vector(:) = (p(:,i_r) - p(:,i))/length
    
    sn = 0.0
    do j = 0, m - 1
        sn = sn + a0*(r**j)
    enddo
    scale_factor = length/sn
    
    sn = 0.0
    b(:,temp_start) = p(:,i)
	do j = 0, m - 2
        sn = sn + a0*(r**j)
		b(:,temp_start+j+1) = p(:,i) + (sn * vector(:))*scale_factor
	enddo
endif

! check surface node
open (Unit=30, File='./output/solid/remesh/check_div_point.plt', STATUS='replace', ACTION='write', IOSTAT=status)
do i = 1, bs
    Write (30,*) b(:,i)
enddo
close(30)

end subroutine divide_point

    
    
subroutine divide_point_bak(pn, p, ins, ins_num, d, enn, line, edge, face_i, fan, face_n, face_p)

implicit none

integer, intent(in) :: pn, ins_num, face_i, fan, line(pn)
integer, intent(inout) :: enn, face_n(fan)
real, intent(in) :: d, p(2,pn), edge(2,pn)
real, intent(inout) :: face_p(2,2000)
logical, intent(in) :: ins

integer :: i, j, i_r, m, bs, count, num
integer :: temp_start, temp_end, temp_m, temp_i, ex_num
integer :: this(pn)
real :: seg, vector(2), length, temp_len
real, allocatable :: b(:,:)

allocate ( b(2,pn*500) )

ex_num = pn
if ( ins ) ex_num = pn - ins_num

bs = 0
temp_len = 0
i = 0

do num = 1, ex_num
	i = i + 1
	if (i > ex_num) exit

	if ( line(i) == 1 ) then
		i_r = i + 1
		if (i == ex_num) i_r = 1

		length = sqrt(((p(1,i_r) - p(1,i)) ** 2) + ((p(2,i_r) - p(2,i)) ** 2))
		vector(:) = p(:,i_r) - p(:,i)
		if ( length > 1.3 * d ) then
			m = 2 * nint(length / (2*d))
		else
			m = 0
		endif

		if ( length > temp_len ) then
			temp_m = m
			temp_len = length
			temp_i = i
		endif
			
		if (m > 0) then
			seg = length / m
			vector(:) = vector(:) / length

			do j = 0, m - 1
				bs = bs + 1
				b(:,bs) = p(:,i) + (j * seg * vector(:))
				
				if ( temp_i == i ) then
					if ( j == 0 ) temp_start = bs
					if ( j == m - 1 ) temp_end = bs
				endif
			enddo
		else
			bs = bs + 1
			b(:,bs) = p(:,i)
			
			if ( temp_i == i ) then
				temp_start = bs;  temp_end = bs
			endif
		endif
	elseif ( line(i) == 2 ) then
		count = 1
		this(count) = i
		do 
			if ( line(i+count) /= 2 ) then
				this(count+1) = i+count
				count = count + 1
				exit
			else
				this(count+1) = i+count
				if ( i+count == ex_num ) then
					this(count+1) = 1
					count = count + 1
					exit
				endif
				count = count + 1
			endif
		enddo

		call spline( pn, bs, b, count, p(:,this(1:count)), line(i), d )
		i = i + count - 2
	elseif ( line(i) == 3 ) then
		count = 3
		this(1:3) = (/i, i+1, i+2 /)
		if (this(2) == ex_num) this(3) = 1
		
		call spline( pn, bs, b, count, p(:,this(1:count)), line(i), d )
		i = i + 1
	endif	
enddo

if ( 2 * (bs / 2) /= bs ) then
    write (*,*) "bs : ", bs, temp_start
	if ( temp_len / temp_m >= d ) then
		do i = bs, temp_end + 1, -1
			b(:,i+1) = b(:,i)
		enddo
		bs = bs + 1
		seg = temp_len / (temp_m + 1) 

		i_r = temp_i + 1
		if ( temp_i == ex_num ) i_r = 1
		vector(:) = ( p(:,i_r) - p(:,temp_i) ) / temp_len
		do i = 0, temp_m
			b(:,temp_start + i) = p(:,temp_i) + (i * seg * vector(:))
		enddo
	elseif ( temp_len / temp_m < d ) then
		do i = temp_end + 1, bs
			b(:,i-1) = b(:,i)
		enddo 
		bs = bs - 1
		seg = temp_len / (temp_m - 1) 
			
		i_r = temp_i + 1
		if ( temp_i == ex_num ) i_r = 1
		vector(:) = ( p(:,i_r) - p(:,temp_i) ) / temp_len
		do i = 0, temp_m - 2
			b(:,temp_start + i) = p(:,temp_i) + (i * seg * vector(:))
		enddo
	endif
endif

temp_len = 0
i = ex_num
if ( ins ) then  ! <== 내부 점이 존재 한다면...
	
	do num = ex_num + 1, pn
		i = i + 1
		if (i > pn) exit

		if ( line(i) == 1 ) then
			i_r = i + 1
			if ( i == pn ) i_r = ex_num + 1

			length = sqrt(((p(1,i_r) - p(1,i)) ** 2) + ((p(2,i_r) - p(2,i)) ** 2))
			vector(:) = p(:,i_r) - p(:,i)
			if ( length > 1.3 * d ) then
				m = 2 * nint(length / (2*d))
			else
				m = 0
			endif
				
			if ( length > temp_len ) then
				temp_m = m
				temp_len = length
				temp_i = i
			endif
				
			if (m > 0) then

				seg = length / m
				
				vector(:) = vector(:) / length
				
				do j = 0, m - 1
					
					bs = bs + 1
					b(:,bs) = p(:,i) + (j * seg * vector(:))
					
					if ( temp_i == i ) then
						if ( j == 0 ) temp_start = bs
						if ( j == m - 1 ) temp_end = bs
					endif
					
				enddo

			else

				bs = bs + 1
				b(:,bs) = p(:,i)
				
				if ( temp_i == i ) then
					temp_start = bs;  temp_end = bs
				endif
											
			endif
		elseif ( line(i) == 2 ) then
			count = 1
			this(count) = i
			do 
				if ( line(i+count) /= 2 ) then
					this(count+1) = i+count
					count = count + 1
					exit
				else
					this(count+1) = i+count
					if ( i+count == pn ) then
						this(count+2) = ex_num+1
						count = count + 2
						exit
					endif
					count = count + 1
				endif
			enddo
						
			call spline( pn, bs, b, count, p(:,this(1:count)), line(i), d )
			i = i + count - 2
		elseif ( line(i) == 3 ) then
			count = 3
			this(1:3) = (/i, i+1, i+2 /)
			if (this(2) == pn) this(3) = ex_num + 1
			
			call spline( pn, bs, b, count, p(:,this(1:count)), line(i), d )
			i = i + 1
		endif

	enddo

	if ( 2 * (bs / 2) /= bs ) then
		if ( temp_len / temp_m >= d ) then
			do i = bs, temp_end + 1, -1
				b(:,i+1) = b(:,i)
			enddo
			bs = bs + 1
			seg = temp_len / (temp_m + 1) 

			i_r = temp_i + 1
			if ( temp_i == pn ) i_r = 1
			vector(:) = ( p(:,i_r) - p(:,temp_i) ) / temp_len
			do i = 0, temp_m
				b(:,temp_start + i) = p(:,temp_i) + (i * seg * vector(1))
			enddo
		elseif ( temp_len / temp_m < d ) then
			do i = temp_end + 1, bs
				b(:,i-1) = b(:,i)
			enddo 
			bs = bs - 1
			seg = temp_len / (temp_m - 1) 
				
			i_r = temp_i + 1
			if ( temp_i == pn ) i_r = 1
			vector(:) = ( p(:,i_r) - p(:,temp_i) ) / temp_len
			do i = 0, temp_m - 2
				b(:,temp_start + i) = p(:,temp_i) + (i * seg * vector(:))
			enddo
		endif
	endif

endif
enn = bs

face_n(face_i) = bs
face_p(:,sum(face_n(1:face_i-1))+1:sum(face_n(1:face_i))) = b(:,1:bs)

end subroutine divide_point_bak
