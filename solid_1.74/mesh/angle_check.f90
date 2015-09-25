!****************************************************************************
!
!  PROGRAM : Angle check
!
!  PURPOSE : 각도 검사 후 봉합 
!		     - 각도가 0 ~ 60, 290 ~ 360 도 일때 봉합
!			 - 봉합 할 수 있는 점이 없을때까지 반복
!
!****************************************************************************
Subroutine angle_check(paving_i, anes, temp_bs, pbn, nn, ne, n_b, b_num, me, node, ele, temp_b_num, temp_b, ec, iter_check)

implicit none

integer, intent(in) :: paving_i, anes, temp_bs, n_b, ne
integer, intent(inout) :: pbn, nn, b_num, ele(5,anes), temp_b_num(temp_bs)
Real, intent(inout) :: me(2,temp_bs), node(2,anes), temp_b(2,temp_bs)
logical, intent(out) :: ec
logical, intent(inout) :: iter_check

integer :: i, j, k, q, o, i_l, i_r, iter, this(3), flag, tol_iter, temp_j
integer :: ang_check, ang_check_num, change_i, del_i, mid_i, cross_judg
real :: beta, change_i_num, sum_ang, dis, ang(2), dist(3), val(2)
logical :: judg

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

iter = 0
tol_iter = pbn
ang_check_loop: do 
    iter = iter + 1
    flag = 0
    do i = 1, pbn
        this = (/ i-1, i, i+1 /)
        if (i == 1) this(1) = pbn
        if (i == pbn) this(3) = 1
        
        call calc_angle(me(:,this(1)), me(:,this(2)), me(:,this(3)), beta)
        call calc_len(me(:,this(1)), me(:,this(2)), dist(1))
        call calc_len(me(:,this(3)), me(:,this(2)), dist(2))
        if (abs(dist(2)) < 10e-16) then
            val(1) = 1.0
        else
            val(1) = dist(1)/dist(2)
        endif
        if (abs(dist(1)) < 10e-16) then
            val(2) = 1.0
        else
            val(2) = dist(2)/dist(1)
        endif
        if (dist(1) < 10e-8 .or. val(1) < 0.1) then
            do j = 1, nn
                call calc_len(node(:,j), me(:,this(1)), dist(3))
                if (dist(3) < 10e-8) then
                    temp_j = j
                    exit
                endif
            enddo
            if (this(1) == 1) then
                me(:,this(1)) = (me(:,this(2))+me(:,pbn))*0.5
            else
                me(:,this(1)) = (me(:,this(2))+me(:,this(1)-1))*0.5
            endif
            node(:,temp_j) = me(:,this(1))
        elseif (dist(2) < 10e-8 .or. val(2) < 0.1) then
            do j = 1, nn
                call calc_len(node(:,j), me(:,this(3)), dist(3))
                if (dist(3) < 10e-8) then
                    temp_j = j
                    exit
                endif
            enddo
            if (this(3) == pbn) then
                me(:,this(3)) = (me(:,this(2))+me(:,1))*0.5
            else
                me(:,this(3)) = (me(:,this(2))+me(:,this(3)+1))*0.5
            endif
            node(:,temp_j) = me(:,this(3))
        endif
        if (((beta >= 0.0) .AND. (beta <= 60.0)) .OR. ((beta >= 290.0) .AND. (beta <= 360.0))) then
            !write (*,*) 'beta:', beta
            change_i_num = 0
            change_i = 0;  del_i = 0;  mid_i = 0
			do j = 1, nn
                do k = 1, 3
                    call calc_len(node(:,j), me(:,this(k)), dist(k))
                enddo                
                
                if (i == 1 .or. i == pbn) then
                    if (dist(3) < 10e-8) then
					    change_i_num = change_i_num + 1
					    if (change_i_num == 1) then
						    change_i = j
					    elseif (change_i_num == 2) then
						    del_i = j
					    endif
				    elseif (dist(1) < 10e-8) then
					    del_i = j
                    elseif (dist(2) < 10e-8) then
					    mid_i = j
				    endif
                else
                    if (dist(1) < 10e-8) then
					    change_i_num = change_i_num + 1
					    if (change_i_num == 1) then
						    change_i = j
					    elseif (change_i_num == 2) then
						    del_i = j
					    endif
				    elseif (dist(3) < 10e-8) then
					    del_i = j
                    elseif (dist(2) < 10e-8) then
					    mid_i = j
                    endif
                endif
                if (change_i /= 0 .and. mid_i /= 0 .and. del_i /= 0) exit
            enddo
            
			if ( change_i > nn .OR. change_i < 1 .OR. del_i > nn .OR. del_i < 1 .OR. mid_i > nn .OR. mid_i < 1 ) then
                write (*,*) 'Fail - angle_check subroutine'
				iter_check = .FALSE.
				exit
            else
                call calc_one_ele(ne, ele(:,1:ne), change_i, del_i, mid_i, judg)
				call calc_ele_cross(anes, nn, ne, ele, change_i, del_i, i, pbn, node, me(:,1:pbn), cross_judg)
            endif
				
		    if ( ( judg ) .AND. ( cross_judg == 0 ) ) then
                !write (*,*) 'i, ci, mi, di:', i, change_i, mid_i, del_i
			    ! 유동절점 봉합
                if (i == 1) then
                    me(:,this(2)) = (me(:,this(1)) + me(:,this(3)))*0.5
					do j = i+1, pbn-2
						me(:,j) = me(:,j+1)
                    enddo
                    node(:,change_i) = me(:,i)
                elseif (i == pbn) then
                    me(:,this(3)) = (me(:,this(1)) + me(:,this(3)))*0.5
                    node(:,change_i) = me(:,1)
                else
			        me(:,this(1)) = (me(:,this(1)) + me(:,this(3)))*0.5
			        do j = i, pbn - 2
				        me(:,j) = me(:,j+2)
                    enddo
                    node(:,change_i) = me(:,this(1))
                endif
			    pbn = pbn - 2
					
			    ! node번호 수정
			    do j = del_i, nn - 1
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
                flag = 1
                exit
            endif
        endif
    enddo
    !write (*,*) 'iter, flag:', iter, flag
    if (iter > tol_iter) iter_check = .false.
    
    if ( iter_check == .FALSE. .or. flag == 0) then
	    exit ang_check_loop
	elseif (pbn <= 6) then
		ec = .TRUE.
		if ( pbn == 6 ) then
		    sum_ang = 0.
			do j = 1, pbn
				if ( j == 1 ) then
					i_l = pbn
					i_r = j + 1
				elseif ( j == pbn ) then
					i_l = j - 1
					i_r = 1
				else
					i_l = j - 1
					i_r = j + 1
				endif
				
				call calc_angle(me(:,i_l), me(:,j), me(:,i_r), beta)
				sum_ang = sum_ang + beta
				
			    if  ( ( (beta >= 0.) .AND. (beta <= 45.) ) .OR. ( (beta .GT. 290.) .AND. (beta <= 360.) ) ) then
				    ec = .FALSE.
			    endif
			    
			    if ( sum_ang > 725. ) then
			        ec = .TRUE.
			        exit
			    endif
			enddo
		elseif (pbn == 4) then
			do j = 1, 2
				if ( j == 1 ) then
					call calc_angle( me(:,4), me(:,1), me(:,2), ang(1) )
					call calc_angle( me(:,2), me(:,3), me(:,4), ang(2) )
					if  ( ( (ang(1) >= 0.) .AND. (ang(1) <= 45.) ) .OR. ( (ang(1) >= 290.) .AND. (ang(1) <= 360.) ) ) then
					if  ( ( (ang(2) >= 0.) .AND. (ang(2) <= 45.) ) .OR. ( (ang(2) >= 290.) .AND. (ang(2) <= 360.) ) ) then
						ec = .FALSE.
						exit
					endif
					endif
				else
					call calc_angle( me(:,1), me(:,2), me(:,3), ang(1) )
					call calc_angle( me(:,3), me(:,4), me(:,1), ang(2) )
					if  ( ( (ang(1) >= 0.) .AND. (ang(1) <= 45.) ) .OR. ( (ang(1) >= 290.) .AND. (ang(1) <= 360.) ) ) then
					if  ( ( (ang(2) >= 0.) .AND. (ang(2) <= 45.) ) .OR. ( (ang(2) >= 290.) .AND. (ang(2) <= 360.) ) ) then
						ec = .FALSE.
					endif
					endif						
				endif
			enddo
		endif				
		
		if ( ec ) then
			!write (*,*) "---=== 남은 경계절점과 유동절점이 6이하가 되었습니다.<angle(1)> ===---"
			do j = 1, pbn
				do k = 1, nn
				    call calc_len(node(:,k), me(:,j), dis) 
					if ( dis < 10e-6 ) then
						temp_b_num(j) = k
						temp_b(:,j) = me(:,j)
						exit
					endif
				enddo
			enddo
			b_num = pbn
			!write (*,'(A,I5)') "b_num :", b_num
			exit ang_check_loop
		endif
	endif
enddo ang_check_loop

end subroutine angle_check
    
    
!****************************************************************************
!
!  PROGRAM : Angle check
!
!  PURPOSE : 각도 검사 후 봉합 
!		     - 각도가 0 ~ 60, 290 ~ 360 도 일때 봉합
!			 - 봉합 할 수 있는 점이 없을때까지 반복
!
!****************************************************************************
!Subroutine angle_check_bak(paving_i, anes, temp_bs, pbn, nn, ne, n_b, b_num, me, node, ele, temp_b_num, temp_b, ec, iter_check)
!
!implicit none
!
!integer, intent(in) :: paving_i, anes, n_b, ne
!integer, intent(inout) :: pbn, nn, b_num, ele(5,anes), temp_b_num(temp_bs)
!Real, intent(inout) :: me(2,pbn), node(2,anes), temp_b(2,temp_bs)
!logical, intent(out) :: ec
!logical, intent(inout) :: iter_check
!
!integer :: i, j, k, q, o, i_l, i_r
!integer :: ang_check, ang_check_num, change_i, del_i, mid_i, cross_judg
!real :: beta, change_i_num, sum_ang, dis, ang(2), dist(3)
!logical :: judg
!
!real, parameter :: pi = 3.141592653589793238462643383279502884197169399375
!q = 0
!k = 1
!ec = .FALSE.
!do ang_check = 1, pbn
!    
!	ang_check_num = 0
!	do i = 1, pbn
!		
!		if (i + ang_check_num > pbn) then
!			
!			exit
!
!		elseif (i == 1) then ! <- 처음일 때...
!
!			call calc_angle(me(:,pbn), me(:,i), me(:,i+1), beta)
!								
!			if (((beta >= 0) .AND. (beta <= 60)) .OR. ((beta >= 290) .AND. (beta <= 360))) then
!				
!				change_i_num = 0
!				do j = n_b, nn
!					if ((node(1,j) == me(1,i+1)) .AND. (node(2,j) == me(2,i+1))) then
!						change_i_num = change_i_num + 1
!						if (change_i_num == 1) then
!							change_i = j
!						elseif (change_i_num == 2) then
!							del_i = j
!						endif
!					elseif ((node(1,j) == me(1,pbn)) .AND. (node(2,j) == me(2,pbn))) then
!						del_i = j
!					elseif ((node(1,j) == me(1,i)) .AND. (node(2,j) == me(2,i))) then
!						mid_i = j
!					endif
!				enddo
!				
!				if ( change_i > nn .OR. change_i < 1 .OR. del_i > nn .OR. del_i < 1 .OR. mid_i > nn .OR. mid_i < 1 ) then
!                    write (*,*) 'Fail - angle_check subroutine'
!				    iter_check = .FALSE.
!				    exit
!				else
!				    call calc_ele_cross(anes, nn, ne, ele, change_i, del_i, i, pbn, node, me(:,1:pbn), cross_judg)
!				endif
!				
!				if ( cross_judg == 0 ) then
!					
!					! 유동절점 봉합
!					me(:,i) = (me(:,pbn) + me(:,i+1)) / 2.0
!									
!					do j = i + 1, pbn - 2
!						me(:,j) = me(:,j+1)
!					enddo
!					pbn = pbn - 2
!
!					! node번호 수정
!					
!					node(:,change_i) = me(:,i)
!					
!					do j = del_i, nn - 1
!						node(:,j) = node(:,j+1)
!					enddo
!					nn = nn - 1
!					
!					! element번호 수정
!
!					do j = 1, ne
!						do o = 2, 5
!							if (ele(o,j) == del_i) then
!								ele(o,j) = change_i
!							elseif (ele(o,j) > del_i) then
!								ele(o,j) = ele(o,j) - 1
!							endif
!						enddo
!					enddo
!					
!					ang_check_num = ang_check_num + 1
!
!				endif
!			
!			endif
!
!		elseif (i == pbn) then ! <- 마지막일 때...
!			
!			call calc_angle( me(:,i-1), me(:,i), me(:,1), beta)
!			
!			if ( ((beta >= 0.) .AND. (beta <= 60.)) .OR. ((beta >= 290.) .AND. (beta <= 360.)) ) then
!				
!				change_i_num = 0
!				do j = n_b, nn	
!					if ((node(1,j) == me(1,1)) .AND. (node(2,j) == me(2,1))) then
!						change_i_num = change_i_num + 1
!						if (change_i_num == 1) then
!							change_i = j
!						elseif (change_i_num == 2) then
!							del_i = j
!						endif
!					elseif ( (node(1,j) == me(1,i-1)) .AND. (node(2,j) == me(2,i-1)) ) then
!						del_i = j
!					elseif ( (node(1,j) == me(1,i)) .AND. (node(2,j) == me(2,i)) ) then
!						mid_i = j
!					endif
!				enddo
!				
!				if ( (me(1,i-1) == me(1,1)) .AND. (me(2,i-1) == me(2,1)) ) then
!					change_i = n_b+1
!				endif
!                
!                if ( change_i > nn .OR. change_i < 1 .OR. del_i > nn .OR. del_i < 1 .OR. mid_i > nn .OR. mid_i < 1 ) then
!                    write (*,*) 'Fail - angle_check subroutine'
!				    iter_check = .FALSE.
!				    exit
!				else
!				    call calc_ele_cross(anes, nn, ne, ele, change_i, del_i, i, pbn, node, me(:,1:pbn), cross_judg)
!				endif
!				
!				if ( cross_judg == 0 ) then
!					
!					! 유동절점 봉합
!					me(:,1) = (me(:,i-1) + me(:,1)) / 2.0
!
!					pbn = pbn - 2
!					
!					! node번호 수정
!
!					node(:,change_i) = me(:,1)
!					do j = del_i, nn - 1
!						node(:,j) = node(:,j+1)
!					enddo
!					
!					nn = nn - 1
!
!					! element번호 수정
!
!					do j = 1, ne
!						do o = 2, 5
!							if (ele(o,j) == del_i) then
!								ele(o,j) = change_i
!							elseif (ele(o,j) > del_i) then
!								ele(o,j) = ele(o,j) - 1			
!							endif
!						enddo
!					enddo
!
!					ang_check_num = ang_check_num + 1
!
!				endif
!
!			endif
!
!		else ! <- 처음과 마지막이 아닐 때...
!			
!			call calc_angle(me(:,i-1), me(:,i), me(:,i+1), beta)
!			
!			if (((beta >= 0) .AND. (beta <= 60)) .OR. ((beta >= 290) .AND. (beta <= 360))) then
!				
!				change_i_num = 0
!				do j = 1, nn
!                    call calc_len(node(:,j), me(:,i-1), dist(1))
!                    call calc_len(node(:,j), me(:,i), dist(2))
!                    call calc_len(node(:,j), me(:,i+1), dist(3))
!					!if ((node(1,j) == me(1,i-1)) .AND. (node(2,j) == me(2,i-1))) then
!                    if (dist(1) < 10e-8) then
!						change_i_num = change_i_num + 1
!						if (change_i_num == 1) then
!							change_i = j
!						elseif (change_i_num == 2) then
!							del_i = j
!						endif
!					!elseif ((node(1,j) == me(1,i+1)) .AND. (node(2,j) == me(2,i+1))) then
!                    elseif (dist(3) < 10e-8) then
!						del_i = j
!                    !elseif ((node(1,j) == me(1,i)) .AND. (node(2,j) == me(2,i))) then
!                    elseif (dist(2) < 10e-8) then
!						mid_i = j
!					endif
!				enddo
!				
!				if ( change_i > nn .OR. change_i < 1 .OR. del_i > nn .OR. del_i < 1 .OR. mid_i > nn .OR. mid_i < 1 ) then
!                    write (*,*) 'Fail - angle_check subroutine'
!				    iter_check = .FALSE.
!				    exit
!				else
!				    call calc_ele_cross(anes, nn, ne, ele, change_i, del_i, i, pbn, node, me(:,1:pbn), cross_judg)
!				
!				    call calc_one_ele(ne, ele(:,1:ne), change_i, del_i, mid_i, judg)
!				endif
!				
!				if ( ( judg ) .AND. ( cross_judg == 0 ) ) then
!					
!					! 유동절점 봉합
!					me(:,i-1) = (me(:,i-1) + me(:,i+1)) / 2.0
!					
!					do j = i, pbn - 2
!						me(:,j) = me(:,j+2)
!					enddo
!					pbn = pbn - 2
!					
!					! node번호 수정
!
!					node(:,change_i) = me(:,i-1)
!					do j = del_i, nn - 1
!						node(:,j) = node(:,j+1)
!					enddo
!					nn = nn - 1
!					
!					! element번호 수정
!
!					do j = 1, ne
!						do o = 2, 5
!							if (ele(o,j) == del_i) then
!								ele(o,j) = change_i
!							elseif (ele(o,j) > del_i) then
!								ele(o,j) = ele(o,j) - 1			
!							endif
!						enddo
!					enddo
!
!					ang_check_num = ang_check_num + 1
!
!				endif
!			
!			endif
!			 
!		endif
!		
!		if ( pbn <= 2 ) exit
!		
!	enddo
!	
!	if ( iter_check == .FALSE. ) then
!	
!	    exit
!
!	elseif (ang_check_num == 0) then
!		
!		exit
!
!	elseif (pbn <= 6) then
!		
!		ec = .TRUE.
!
!		if ( pbn == 6 ) then
!		    sum_ang = 0.
!			do j = 1, pbn
!				if ( j == 1 ) then
!					i_l = pbn
!					i_r = j + 1
!				elseif ( j == pbn ) then
!					i_l = j - 1
!					i_r = 1
!				else
!					i_l = j - 1
!					i_r = j + 1
!				endif
!				
!				call calc_angle(me(:,i_l), me(:,j), me(:,i_r), beta)
!				sum_ang = sum_ang + beta
!				
!			    if  ( ( (beta >= 0.) .AND. (beta <= 45.) ) .OR. ( (beta .GT. 290.) .AND. (beta <= 360.) ) ) then
!				    ec = .FALSE.
!			    endif
!			    
!			    if ( sum_ang > 725. ) then
!			        ec = .TRUE.
!			        exit
!			    endif
!			enddo
!		elseif (pbn == 4) then
!			do j = 1, 2
!				if ( j == 1 ) then
!					call calc_angle( me(:,4), me(:,1), me(:,2), ang(1) )
!					call calc_angle( me(:,2), me(:,3), me(:,4), ang(2) )
!					if  ( ( (ang(1) >= 0.) .AND. (ang(1) <= 45.) ) .OR. ( (ang(1) >= 290.) .AND. (ang(1) <= 360.) ) ) then
!					if  ( ( (ang(2) >= 0.) .AND. (ang(2) <= 45.) ) .OR. ( (ang(2) >= 290.) .AND. (ang(2) <= 360.) ) ) then
!						ec = .FALSE.
!						exit
!					endif
!					endif
!				else
!					call calc_angle( me(:,1), me(:,2), me(:,3), ang(1) )
!					call calc_angle( me(:,3), me(:,4), me(:,1), ang(2) )
!					if  ( ( (ang(1) >= 0.) .AND. (ang(1) <= 45.) ) .OR. ( (ang(1) >= 290.) .AND. (ang(1) <= 360.) ) ) then
!					if  ( ( (ang(2) >= 0.) .AND. (ang(2) <= 45.) ) .OR. ( (ang(2) >= 290.) .AND. (ang(2) <= 360.) ) ) then
!						ec = .FALSE.
!					endif
!					endif						
!				endif
!			enddo
!		endif				
!		
!		if ( ec ) then
!			!write (*,*) "---=== 남은 경계절점과 유동절점이 6이하가 되었습니다.<angle(1)> ===---"
!			do j = 1, pbn
!				do k = 1, nn
!				    call calc_len(node(1:2,k), me(1:2,j), dis) 
!					if ( dis < 10e-6 ) then
!						temp_b(1,j) = k
!						temp_b(2:3,j) = me(:,j)
!						exit
!					endif
!				enddo
!			enddo
!			b_num = pbn
!			!write (*,'(A,I5)') "b_num :", b_num
!			exit
!		endif
!
!	endif
!	
!enddo
!!write (*,*) "각도 검사 후의 절점"
!!write (*,'(A,I3)') "pbn : ", pbn
!
!end subroutine