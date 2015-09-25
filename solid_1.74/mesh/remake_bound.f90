!****************************************************************************
!
!  PROGRAM : Remake boundary
!
!  PURPOSE : 경계선 재설정
!
!****************************************************************************
subroutine remake_bound( anes, pn, corner_node, nn, node, d, val, ins, ins_num, enn )

implicit none

integer, intent(in) :: anes, pn, ins_num
integer, intent(inout) :: nn 
integer, intent(out) :: enn
real, intent(in) :: d, val, corner_node(pn,2)
real, intent(inout) :: node(anes,3)
logical, intent(in) :: ins

integer :: i, j, i_r, j_r, k, m, in, temp_i, temp_start, temp_end, temp_m, cri, ex_num, temp_nn
integer :: corner_num(pn)
real :: vector(2), temp_p(2,2)
real :: length, seg, temp_len, dis(2)
logical :: cross

ex_num = pn
if ( ins ) ex_num = pn - ins_num

cri = nint(d / val)
temp_len = 0

do i = 1, ex_num

	if (i == ex_num) then ! <- 마지막 선이라면...
	
		call calc_len(corner_node(1,1:2), corner_node(i,1:2), length)
		vector(1) = corner_node(1,1) - corner_node(i,1)
		vector(2) = corner_node(1,2) - corner_node(i,2)
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
			vector(1) = vector(1) / length
			vector(2) = vector(2) / length
			
			do j = 0, m - 1
				
				nn = nn + 1
				node(nn,1) = nn
				node(nn,2) = corner_node(i,1) + (j * seg * vector(1))
				node(nn,3) = corner_node(i,2) + (j * seg * vector(2))
				if (j == 0 ) corner_num(i) = nn

				if ( temp_i == i ) then
					if ( j == 0 ) temp_start = nn
					if ( j == m - 1 ) temp_end = nn
				endif
				
			enddo

		else

			nn = nn + 1
			node(nn,1) = nn
			node(nn,2:3) = corner_node(i,1:2)
			corner_num(i) = nn

			if ( temp_i == i ) then
				temp_start = nn;  temp_end = nn
			endif
					
		endif

	else
		
		call calc_len(corner_node(i,1:2), corner_node(i+1,1:2), length)
		vector(1) = corner_node(i+1,1) - corner_node(i,1)
		vector(2) = corner_node(i+1,2) - corner_node(i,2)
		if ( length > 1.3 * d ) then
			m = 2 * nint(length / (2*d))
		else
			m = 0
		endif
		if ( length > temp_len ) then
			temp_len = length
			temp_m = m
			temp_i = i
		endif
				
		if (m > 0) then

			seg = length / m
			vector(1) = vector(1) / length
			vector(2) = vector(2) / length
						
			do j = 0, m - 1
				
				nn = nn + 1
				node(nn,1) = nn
				node(nn,2) = corner_node(i,1) + (j * seg * vector(1))
				node(nn,3) = corner_node(i,2) + (j * seg * vector(2))
				if (j == 0 ) corner_num(i) = nn

				if ( temp_i == i ) then
					if ( j == 0 ) temp_start = nn
					if ( j == m - 1 ) temp_end = nn
				endif
								
			enddo

		else

			nn = nn + 1
			node(nn,1) = nn
			node(nn,2:3) = corner_node(i,1:2)
			corner_num(i) = nn

			if ( temp_i == i ) then
				temp_start = nn;  temp_end = nn
			endif
											
		endif

	endif
	
enddo

in = 0; temp_nn = nn * 2
do i = 1, temp_nn
	in = in + 1
	i_r = in + 1
	if ( in == nn ) i_r = 1

	temp_p(1,1:2) = (node(i_r,2:3) + node(in,2:3)) * 0.5 
	call calc_len( node(i_r,2:3), node(in,2:3), dis(1) )
	vector(1) = (node(i_r,2) - node(in,2)) * 1.2
	vector(2) = (node(i_r,3) - node(in,3)) * 1.2
	
	temp_p(2,1) = temp_p(1,1) - vector(2)
	temp_p(2,2) = temp_p(1,2) + vector(1)

	do j = 1, nn
		if ( (j < in - 6) .OR. (j > in + 5) ) then
			j_r = j + 1
			if ( j == nn ) j_r = 1
			call calc_cross( temp_p(1,1:2), temp_p(2,1:2), node(j,2:3), node(j_r,2:3), cross )
			if ( cross .EQV. .FALSE. ) then
				do k = nn, in+1, -1
					node(k+1,2:3) = node(k,2:3)
				enddo
				nn = nn + 1
				node(nn,1) = nn
				
				node(in+1,2:3) = temp_p(1,1:2)
				if ( temp_start > in ) temp_start = temp_start + 1
				if ( temp_end > in ) temp_end = temp_end + 1

				exit

			endif
		endif
	enddo

	if ( in == nn ) exit
enddo


temp_m = temp_end - temp_start + 1
if ( 2 * (nn / 2) /= nn ) then  ! <= 경계 절점이 짝수가 아니라면...
	if ( temp_len / temp_m >= d ) then
		do i = nn, temp_end + 1, -1
			node(i+1,1) = node(i,1) + 1
			node(i+1,2:3) = node(i,2:3)
		enddo
		nn = nn + 1
		seg = temp_len / (temp_m + 1) 

		i_r = temp_i + 1
		if ( temp_i == ex_num ) i_r = 1
		vector(1) = ( corner_node(i_r,1) - corner_node(temp_i,1) ) / temp_len
		vector(2) = ( corner_node(i_r,2) - corner_node(temp_i,2) ) / temp_len
		do i = 0, temp_m
			node(temp_start + i,1) = temp_start + i
			node(temp_start + i,2) = corner_node(temp_i,1) + (i * seg * vector(1))
			node(temp_start + i,3) = corner_node(temp_i,2) + (i * seg * vector(2))
		enddo
	elseif ( temp_len / temp_m < d ) then
		do i = temp_end + 1, nn
			node(i-1,2:3) = node(i,2:3)
		enddo 
		nn = nn - 1
		seg = temp_len / (temp_m - 1) 
			
		i_r = temp_i + 1
		if ( temp_i == ex_num ) i_r = 1
		vector(1) = ( corner_node(i_r,1) - corner_node(temp_i,1) ) / temp_len
		vector(2) = ( corner_node(i_r,2) - corner_node(temp_i,2) ) / temp_len
		do i = 0, temp_m - 2
			node(temp_start + i,1) = temp_start + i
			node(temp_start + i,2) = corner_node(temp_i,1) + (i * seg * vector(1))
			node(temp_start + i,3) = corner_node(temp_i,2) + (i * seg * vector(2))
		enddo
	endif
endif

temp_len = 0
enn = nn
if ( ins ) then  ! <== 내부 점이 존재 한다면...
	
	do i = ex_num + 1, pn
		
		i_r = i + 1
		if ( i == pn ) i_r = ex_num + 1

		call calc_len( corner_node(i_r,1:2), corner_node(i,1:2), length )
		vector(1) = corner_node(i_r,1) - corner_node(i,1)
		vector(2) = corner_node(i_r,2) - corner_node(i,2)
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
			
			vector(1) = vector(1) / length
			vector(2) = vector(2) / length
			
			do j = 0, m - 1
				
				nn = nn + 1
				node(nn,1) = nn
				node(nn,2) = corner_node(i,1) + (j * seg * vector(1))
				node(nn,3) = corner_node(i,2) + (j * seg * vector(2))
				
				if ( temp_i == i ) then
					if ( j == 0 ) temp_start = nn
					if ( j == m - 1 ) temp_end = nn
				endif
				
			enddo

		else

			nn = nn + 1
			node(nn,1) = nn
			node(nn,2) = corner_node(i,1)
			node(nn,3) = corner_node(i,2)
			
			if ( temp_i == i ) then
				temp_start = nn;  temp_end = nn
			endif
										
		endif

	enddo

	if ( 2 * (nn / 2) /= nn ) then
		if ( temp_len / temp_m >= d ) then
			do i = nn, temp_end + 1, -1
				node(i+1,1) = node(i,1) + 1
				node(i+1,2:3) = node(i,2:3)
			enddo
			nn = nn + 1
			seg = temp_len / (temp_m + 1) 

			i_r = temp_i + 1
			if ( temp_i == pn ) i_r = 1
			vector(1) = ( corner_node(i_r,1) - corner_node(temp_i,1) ) / temp_len
			vector(2) = ( corner_node(i_r,2) - corner_node(temp_i,2) ) / temp_len
			do i = 0, temp_m
				node(temp_start + i,1) = temp_start + i
				node(temp_start + i,2) = corner_node(temp_i,1) + (i * seg * vector(1))
				node(temp_start + i,3) = corner_node(temp_i,2) + (i * seg * vector(2))
			enddo
		elseif ( temp_len / temp_m < d ) then
			do i = temp_end + 1, nn
				node(i-1,2:3) = node(i,2:3)
			enddo 
			nn = nn - 1
			seg = temp_len / (temp_m - 1) 
				
			i_r = temp_i + 1
			if ( temp_i == pn ) i_r = 1
			vector(1) = ( corner_node(i_r,1) - corner_node(temp_i,1) ) / temp_len
			vector(2) = ( corner_node(i_r,2) - corner_node(temp_i,2) ) / temp_len
			do i = 0, temp_m - 2
				node(temp_start + i,1) = temp_start + i
				node(temp_start + i,2) = corner_node(temp_i,1) + (i * seg * vector(1))
				node(temp_start + i,3) = corner_node(temp_i,2) + (i * seg * vector(2))
			enddo
		endif
	endif

endif

end subroutine
!============================================================================

!============================================================================
!****************************************************************************
!
!  PROGRAM : Remake divided boundary
!
!  PURPOSE : 갈라진 경계선 재설정
!
!****************************************************************************
subroutine remake_divided_bound( anes, pn, corner_node, nn, node, d, val, on, opn, op_nn)

implicit none

integer, intent(in) :: anes, on, pn
integer, intent(inout) :: nn ,opn(on), op_nn(on,pn*1000)
real, intent(in) :: d, val, corner_node(pn,2)
real, intent(inout) :: node(anes,3)

integer :: i, j, i_r, m, temp_i, temp_start, temp_end, temp_m, cri, on_i, count
integer :: temp_opn(on), temp_op_nn(on,pn)
real :: vector(2), length, seg, temp_len

cri = nint(d / val)
temp_len = 0
temp_op_nn(1:on,1:pn) = op_nn(1:on,1:pn)
temp_opn = opn
op_nn = 0
nn = 0

do on_i = 1, on
	count = 0
	do i = 1, temp_opn(on_i)
		i_r = i + 1
		if (i == temp_opn(on_i)) i_r = 1
		
		call calc_len(corner_node(temp_op_nn(on_i,i),1:2), corner_node(temp_op_nn(on_i,i_r),1:2), length)
		vector(1:2) = corner_node(temp_op_nn(on_i,i_r),1:2) - corner_node(temp_op_nn(on_i,i),1:2)
		
		if ( length > 1.3 * d ) then
			m = 2 * nint(length / (2*d))
		else
			m = 0
		endif
			
		if ( length > temp_len ) then
			temp_m = m
			temp_len = length
			temp_i = temp_op_nn(on_i,i)
		endif

		if (m > 0) then

			seg = length / m
			vector = vector / length
							
			do j = 0, m - 1
				
				count = count + 1
				nn = nn + 1
				node(nn,1) = nn
				node(nn,2) = corner_node(temp_op_nn(on_i,i),1) + (j * seg * vector(1))
				node(nn,3) = corner_node(temp_op_nn(on_i,i),2) + (j * seg * vector(2))
				op_nn(on_i,count) = nn
				
				if ( temp_i == i ) then
					if ( j == 0 ) temp_start = nn
					if ( j == m - 1 ) temp_end = nn
				endif
				
			enddo

		else
			
			count = count + 1
			nn = nn + 1
			node(nn,1) = nn
			node(nn,2:3) = corner_node(temp_op_nn(on_i,i),1:2)
			op_nn(on_i,count) = nn
			
			if ( temp_i == i ) then
				temp_start = nn;  temp_end = nn
			endif
					
		endif

		opn(on_i) = count
			
	enddo

!	in = 0; temp_nn = opn(on_i) * 2
!	do i = 1, temp_nn
!		in = in + 1
!		i_r = in + 1
!		if ( in == opn(on_i) ) i_r = 1
!
!		temp_p(1,1:2) = (node(op_nn(on_i,in),2:3) + node(op_nn(on_i,i_r),2:3)) * 0.5 
!		call calc_len( node(op_nn(on_i,i_r),2:3), node(op_nn(on_i,in),2:3), dis(1) )
!		vector(1) = (node(op_nn(on_i,i_r),2) - node(op_nn(on_i,in),2)) * 1.2
!		vector(2) = (node(op_nn(on_i,i_r),3) - node(op_nn(on_i,in),3)) * 1.2
!		
!		temp_p(2,1) = temp_p(1,1) - vector(2)
!		temp_p(2,2) = temp_p(1,2) + vector(1)
!
!		do j = 1, opn(on_i)
!			if ( (j < in - 6) .OR. (j > in + 5) ) then
!				j_r = j + 1
!				if ( j == nn ) j_r = 1  <= 여기서 부터 수정
!				call calc_cross( temp_p(1,1:2), temp_p(2,1:2), node(op_nn(on_i,j),2:3), node(op_nn(on_i,j_r),2:3), cross )
!				if ( cross .EQV. .FALSE. ) then
!					do k = opn(on_i), in+1, -1
!						node(k+1,2:3) = node(k,2:3)
!					enddo
!					nn = nn + 1
!					node(nn,1) = nn
!					
!					node(in+1,2:3) = temp_p(1,1:2)
!					if ( temp_start > in ) temp_start = temp_start + 1
!					if ( temp_end > in ) temp_end = temp_end + 1
!
!					exit
!
!				endif
!			endif
!		enddo
!
!		if ( in == nn ) exit
!	enddo

	temp_m = temp_end - temp_start + 1
	if ( 2 * (opn(on_i) / 2) /= opn(on_i) ) then  ! <= 경계 절점이 짝수가 아니라면...
		if ( temp_len / temp_m >= d ) then
			do i = opn(on_i), temp_end + 1, -1
				node(op_nn(on_i,i+1),1) = node(op_nn(on_i,i),1) + 1
				node(op_nn(on_i,i+1),2:3) = node(op_nn(on_i,i),2:3)
				op_nn(on_i,i) = op_nn(on_i,i) + 1
			enddo
			nn = nn + 1
			opn(on_i) = opn(on_i) + 1
			seg = temp_len / (temp_m + 1) 

			i_r = temp_i + 1
			if ( temp_i == opn(on_i) ) i_r = temp_op_nn(on_i,1)
			vector(1) = ( corner_node(i_r,1) - corner_node(temp_i,1) ) / temp_len
			vector(2) = ( corner_node(i_r,2) - corner_node(temp_i,2) ) / temp_len
			do i = 0, temp_m
				node(temp_start + i,1) = temp_start + i
				node(temp_start + i,2) = corner_node(temp_i,1) + (i * seg * vector(1))
				node(temp_start + i,3) = corner_node(temp_i,2) + (i * seg * vector(2))
			enddo
		elseif ( temp_len / temp_m < d ) then
			do i = temp_end + 1, opn(on_i)
				node(i-1,2:3) = node(i,2:3)
			enddo 
			nn = nn - 1
			opn(on_i) = opn(on_i) -1
			seg = temp_len / (temp_m - 1) 
				
			i_r = temp_i + 1
			if ( temp_i == opn(on_i) ) i_r = temp_op_nn(on_i,1)
			vector(1) = ( corner_node(i_r,1) - corner_node(temp_i,1) ) / temp_len
			vector(2) = ( corner_node(i_r,2) - corner_node(temp_i,2) ) / temp_len
			do i = 0, temp_m - 2
				node(temp_start + i,1) = temp_start + i
				node(temp_start + i,2) = corner_node(temp_i,1) + (i * seg * vector(1))
				node(temp_start + i,3) = corner_node(temp_i,2) + (i * seg * vector(2))
			enddo
		endif
	endif
enddo

end subroutine

