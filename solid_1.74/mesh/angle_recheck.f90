!****************************************************************************
!
!  PROGRAM : Angle check
!
!  PURPOSE : 각도 검사 후 봉합 
!		     - 각도가 0 ~ 60, 290 ~ 360 도 일때 봉합
!			 - 봉합 할 수 있는 점이 없을때까지 반복
!
!****************************************************************************
Subroutine angle_recheck(anes, temp_bs, b_num, bs, nn, ne, bs_num, b, temp_b_num, temp_b, node, ele, dpn, div, ec)

implicit none

integer, intent(in) :: anes, temp_bs, ne, dpn
integer, intent(inout) :: bs, b_num, nn, ele(5,anes), div(dpn,2), bs_num(temp_bs), temp_b_num(temp_bs)
real, intent(inout) :: b(2,temp_bs), node(2,anes), temp_b(2,temp_bs)
logical, intent(out) :: ec

integer :: i, j, k, q, o, i_l, i_r
integer :: ang_check, ang_check_num, change_i, del_i, mid_i, cross_judg
real :: beta, ang(2), temp_b_nodes(2,bs)

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

q = 0
k = 1
ec = .FALSE.

do ang_check = 1, bs
	
	ang_check_num = 0
	do i = 1, bs
		if (i+ang_check_num > bs) then
			exit
		elseif (i == 1) then ! <- 처음일 때...

			call calc_angle(b(:,bs), b(:,i), b(:,i+1), beta)
			
			if (((beta >= 0) .AND. (beta <= 60)) .OR. ((beta >= 290) .AND. (beta <= 360))) then
				change_i = 0
                del_i = 0
                mid_i = 0
				do j = 1, nn	
					if (j == bs_num(i+1)) then
						change_i = j
					elseif (j == bs_num(bs)) then
						del_i = j
					elseif (j == bs_num(i)) then
						mid_i = j
					endif
                enddo
                
                if (change_i == 0 .or. del_i == 0 .or. mid_i == 0) exit
                do j = 1, bs
				    temp_b_nodes(:,j) = b(:,j)
                enddo
				call calc_ele_cross(anes, nn, ne, ele, change_i, del_i, i, bs, node, temp_b_nodes, cross_judg)
                do j = 1, bs
				    b(:,j) = temp_b_nodes(:,j)
                enddo
                
				if ( cross_judg == 0 ) then			
					! 유동절점 봉합
					b(i,1) = change_i
					b(:,i) = (b(:,bs) + b(:,i+1)) / 2.0
					
					do j = i + 1, bs - 2
						b(:,j) = b(:,j+1)
					enddo
					bs = bs - 2

					! node번호 수정
					
					node(:,change_i) = b(:,i)
									
					do j = del_i, nn - 1
						node(:,j) = node(:,j+1)
					enddo
					nn = nn - 1
					
					! element번호 수정

					do j = 1, ne 
						do o = 2, 5
							if (ele(o,j) .EQ. del_i) then
								ele(o,j) = change_i
							elseif (ele(o,j) .GT. del_i) then
								ele(o,j) = ele(o,j) - 1
							endif
						enddo
					enddo

					! div의 node번호 수정
					
					do j = 1, dpn
						if ( div(j,2) == del_i ) then
							div(j,2) = change_i
						elseif ( div(j,2) > del_i ) then
							div(j,2) = div(j,2) - 1
						endif
					enddo
					
					ang_check_num = ang_check_num + 1
				
				endif

			endif

		elseif (i == bs) then ! <- 마지막일 때...

			call calc_angle(b(:,i-1), b(:,i), b(:,1), beta)
			
			if (((beta >= 0) .AND. (beta <= 60)) .OR. ((beta >= 290) .AND. (beta <= 360))) then
				!write (*,*) change_i, del_i, mid_i
                change_i = 0
                del_i = 0
                mid_i = 0
				do j = 1, nn	
					if (j == bs_num(1)) then
						change_i = j
					elseif (j == bs_num(i-1)) then
						del_i = j
					elseif (j == bs_num(i)) then
						mid_i = j
					endif
				enddo
				!write (*,*) change_i, del_i, mid_i
                
                if (change_i == 0 .or. del_i == 0 .or. mid_i == 0) exit
                do j = 1, bs
				    temp_b_nodes(:,j) = b(:,j)
                enddo
				call calc_ele_cross(anes, nn, ne, ele, change_i, del_i, i, bs, node, temp_b_nodes, cross_judg)
				do j = 1, bs
				    b(:,j) = temp_b_nodes(:,j)
                enddo

				if ( cross_judg == 0 ) then
					! 유동절점 봉합
					bs_num(1) = change_i
					b(:,1) = (b(:,i-1) + b(:,1)) / 2.0

					bs = bs - 2

					! node번호 수정

					node(:,change_i) = b(:,1)

					do j = del_i, nn - 1
						node(:,j) = node(:,j+1)
					enddo
					nn = nn - 1
					
					! element번호 수정

					do j = 1, ne 
						do o = 2, 5
							if (ele(o,j) == del_i) then
								ele(o,j) = change_i
							elseif (ele(o,j) .GT. del_i) then
								ele(o,j) = ele(o,j) - 1			
							endif
						enddo
					enddo

					! div의 node번호 수정
					
					do j = 1, dpn
						if ( div(j,2) == del_i ) then
							div(j,2) = change_i
						elseif ( div(j,2) > del_i ) then
							div(j,2) = div(j,2) - 1
						endif
					enddo

					ang_check_num = ang_check_num + 1

				endif

			endif

		else ! <- 처음과 마지막이 아닐 때...

			call calc_angle(b(:,i-1), b(:,i), b(:,i+1), beta)
			
			if (((beta >= 0) .AND. (beta <= 60)) .OR. ((beta >= 290) .AND. (beta <= 360))) then
				change_i = 0
                del_i = 0
                mid_i = 0
				do j = 1, nn	
					if (j == bs_num(i-1)) then
						change_i = j
					elseif (j == bs_num(i+1)) then
						del_i = j
					elseif (j == bs_num(i)) then
						mid_i = j
					endif
                enddo
				
                if (change_i == 0 .or. del_i == 0 .or. mid_i == 0) exit
                do j = 1, bs
				    temp_b_nodes(:,j) = b(:,j)
                enddo
				call calc_ele_cross(anes, nn, ne, ele, change_i, del_i, i, bs, node, temp_b_nodes, cross_judg)
				do j = 1, bs
				    b(:,j) = temp_b_nodes(:,j)
                enddo
				!call calc_ele_cross(anes, nn, ne, ele, change_i, del_i, i, bs, node, b(:,1:bs), cross_judg)

				if ( cross_judg == 0 ) then		
					! 유동절점 봉합
					bs_num(i-1) = change_i
					b(:,i-1) = (b(:,i-1) + b(:,i+1)) / 2.0
					
					do j = i, bs - 2
						if ( bs_num(j+2) > del_i ) then
							bs_num(j) = bs_num(j+2) - 1
						else 
							bs_num(j) = bs_num(j+2)
						endif
						b(:,j) = b(:,j+2)
					enddo
					bs = bs - 2
					
					! node번호 수정

					node(:,change_i) = b(:,i-1)

					do j = del_i, nn
						node(:,j) = node(:,j+1)
					enddo
					nn = nn - 1
								
					! element번호 수정

					do j = 1, ne 
						do o = 2, 5
							if (ele(o,j) == del_i) then
								ele(o,j) = change_i
							elseif (ele(o,j) > del_i) then
								ele(o,j) = ele(o,j) - 1			
							endif
						enddo
					enddo
					
					! div의 node번호 수정
					
					do j = 1, dpn
						if ( div(j,2) == del_i ) then
							div(j,2) = change_i
						elseif ( div(j,2) > del_i ) then
							div(j,2) = div(j,2) - 1
						endif
					enddo

					ang_check_num = ang_check_num + 1
				
				endif

			endif
			 
		endif
		
		if (bs <= 2) then
			exit
		endif
	
	enddo
	
	if (bs <= 6) then

		ec = .TRUE.
	
		if ( bs == 6 ) then
			do j = 1, bs
				if ( j == 1 ) then
					i_l = bs
					i_r = j + 1
				elseif ( j == bs ) then
					i_l = j - 1
					i_r = 1
				else
					i_l = j - 1
					i_r = j + 1
				endif
				
				call calc_angle(b(:,i_l), b(:,j), b(:,i_r), beta)
				if  ( ( (beta >= 0.) .AND. (beta <= 45.) ) .OR. ( (beta .GT. 290.) .AND. (beta <= 360.) ) ) then
					ec = .FALSE.
					exit
				endif
			enddo

		elseif (bs == 4) then
			do j = 1, 2
				if ( j == 1 ) then
					call calc_angle( b(:,4), b(:,1), b(:,2), ang(1) )
					call calc_angle( b(:,2), b(:,3), b(:,4), ang(2) )
					if  ( ( (ang(1) >= 0.) .AND. (ang(1) <= 45.) ) .OR. ( (ang(1) .GT. 290.) .AND. (ang(1) <= 360.) ) ) then
					if  ( ( (ang(2) >= 0.) .AND. (ang(2) <= 45.) ) .OR. ( (ang(2) .GT. 290.) .AND. (ang(2) <= 360.) ) ) then
						ec = .FALSE.
						exit
					endif
					endif
				else
					call calc_angle( b(:,1), b(:,2), b(:,3), ang(1) )
					call calc_angle( b(:,3), b(:,4), b(:,1), ang(2) )
					if  ( ( (ang(1) >= 0.) .AND. (ang(1) <= 45.) ) .OR. ( (ang(1) .GT. 290.) .AND. (ang(1) <= 360.) ) ) then
					if  ( ( (ang(2) >= 0.) .AND. (ang(2) <= 45.) ) .OR. ( (ang(2) .GT. 290.) .AND. (ang(2) <= 360.) ) ) then
						ec = .FALSE.
					endif
					endif						
				endif
			enddo
		endif

		if ( ec ) then
			!write (*,*) "---=== 남은 경계절점과 유동절점이 6이하가 되었습니다.<angle(2)> ===---"
			do j = 1, bs
                temp_b_num(j) = bs_num(j)
				temp_b(:,j) = b(:,j)
				!write (*,*) temp_b(j,1:3)
			enddo
			b_num = bs
			!write (*,'(A,I5)') "b_num :", b_num
			exit
		endif

	endif

	if (ang_check_num == 0) then
		exit
	endif
enddo

end subroutine