!****************************************************************************
!
!  PROGRAM : Clock check
!
!  PURPOSE : 6개이하의 절점 검사
!			 1. 남은 절점들이 시계방향성을 가지는 검사
!			 2. 남은 절점들 간의 교차 검사(보류)		
!
!****************************************************************************
Subroutine clock_check(anes, b_num, node, ne, temp_b_num, temp_b, ele, iter_check)

implicit none

integer, intent(in) :: anes, b_num, ne, ele(5,anes), temp_b_num(b_num)
real, intent(in) :: temp_b(2,b_num)
real, intent(inout) :: node(2,anes)
logical, intent(inout) :: iter_check

integer :: i, i_l, i_r, j, j_r, k, ele_num, r_p(6,2), count
integer :: r_ele(6,2), temp_ele(10)
real :: sum_angle, alpha, p(2)
logical :: cross, over, check

count = 0
cross = .TRUE.
!write (*,*) 'b_num:', b_num
!do i = 1, b_num
!    write (*,*) temp_b_num(i), temp_b(:,i)
!enddo
do i = 1, b_num - 1
	do j = i + 1, b_num
		if ( j == b_num ) then
			j_r = 1
		else
			j_r = j + 1
		endif

		call calc_cross(temp_b(:,i), temp_b(:,i+1), temp_b(:,j), temp_b(:,j_r), cross)
		if ( cross .EQV. .FALSE. ) exit
	enddo
enddo

over = .FALSE.
if ( cross ) then
	sum_angle = 0
	do i = 1, b_num
		if ( i == 1 ) then
			i_l = b_num;   i_r = 2
		elseif (i == b_num ) then
			i_l = i-1;     i_r = 1
		else
			i_l = i-1;     i_r = i+1
		endif
		
		call calc_angle(temp_b(:,i_l), temp_b(:,i), temp_b(:,i_r), alpha)

		sum_angle = sum_angle + alpha
	enddo
	
	if ( b_num == 4 ) then
		if ( sum_angle > 370. ) over = .TRUE.
	elseif ( b_num == 6 ) then
		if ( sum_angle > 740. ) over = .TRUE.
	endif

	if ( over ) then
		do i = 1, b_num	
			ele_num = 0
			do j = 1, ne
				do k = 2, 5
					if ( ele(k,j) == temp_b_num(i) ) then
						ele_num = ele_num + 1
                        if (ele_num > 10) then
                            iter_check = .FALSE.
                            stop
                            exit
                        endif
						temp_ele(ele_num) = j
					endif
                enddo
                if (iter_check == .false.) exit
            enddo
			if (iter_check == .false.) exit
            
			if ( ele_num  == 2 ) then 
				count = count + 1
				r_p(count,1) = temp_b_num(i)
				r_ele(count,1:2) = temp_ele(1:2)
			elseif ( ele_num > 8 ) then
			    iter_check = .FALSE.
			    exit
			endif

		enddo 
	endif
			
endif

if ( iter_check ) then
    do i = 1, count
        check = .FALSE.
	    do j = 2, 5
		    do k = 2, 5
			    if ( (ele(j,r_ele(i,1)) == ele(k,r_ele(i,2))) .AND. (ele(j,r_ele(i,1)) /= r_p(i,1)) ) then
				    r_p(i,2) = ele(j,r_ele(i,1))
				    check = .TRUE.
				    exit
			    endif
		    enddo
		    if ( check ) exit
	    enddo
	    
    	if ( check == .FALSE. ) then
            write (*,*) 'Fail - clock_ckeck subroutine'
    	    iter_check = .FALSE.
    	    exit
    	endif
    	
	    cross = .TRUE.
	    do j = 1, b_num
		    j_r = j + 1
		    if ( j == b_num ) j_r = 1
    			
		    call calc_cross(temp_b(:,j), temp_b(:,j_r), node(:,r_p(i,1)), node(:,r_p(i,2)), cross)
    		
		    if ( cross .EQV. .FALSE. ) then
			    call calc_cross_p( temp_b(:,j), temp_b(:,j_r), node(:,r_p(i,1)), node(:,r_p(i,2)), p )
			    node(:,r_p(i,1)) = ( p(:) + node(:,r_p(i,2)) ) * 0.5
			    exit
		    endif
	    enddo
	    if ( cross ) then
		    node(:,r_p(i,1)) = ( node(:,r_p(i,1)) + node(:,r_p(i,2)) ) * 0.5
	    endif
    enddo
endif

end subroutine clock_check