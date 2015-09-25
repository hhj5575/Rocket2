!****************************************************************************
!
!  PROGRAM : Cross_Near
!
!  PURPOSE : 투영된 절점의 교차와 근접 확인 후 봉합
!
!****************************************************************************
Subroutine cross_near(temp_bs, n_b, anes, nn, ne, pbn, node, ele, me, cross_num, near_num, cross, near, wi)

implicit none

integer, intent(in) :: anes, temp_bs, n_b, wi
integer, intent(inout) :: nn, ne, pbn, cross_num, near_num, ele(5,anes), cross(4,temp_bs), near(4,temp_bs)
real, intent(inout) :: node(2,anes), me(2,temp_bs)

integer :: i, j, k, o, change_i, del_i, change_i_num, extra_num
integer :: temp_i, temp_j, near_check_recheck, near_check, num, temp_pbn
integer :: this(4)
real :: length, cri, cur_d, temp_dis, tol, cri_ang
real :: alpha(2), beta(2), me_dis(pbn), val(2)
logical :: check

cross = 0.;     near = 0
cross_num = 0;  near_num = 0
cri_ang = 45.
!----------------------------------------------------------------------------------------
me_dis = 0.0
temp_pbn = pbn
do i = 1, pbn
    this(1:3) = (/ i-1, i, i+1 /)
    if (i == 1) this(1) = pbn
    if (i == pbn) this(3) = 1
    call calc_len(me(:,this(1)), me(:,this(2)), val(1))
    call calc_len(me(:,this(2)), me(:,this(3)), val(2))
    me_dis(i) = min(val(1),val(2))
enddo
!do i = 1, pbn-1
!    do j = i+1, pbn
!        if (me_dis(i) > me_dis(j)) then
!            temp_dis = me_dis(i)
!            me_dis(i) = me_dis(j)
!            me_dis(j) = temp_dis
!        endif
!    enddo
!enddo
num = int(pbn*0.2)+1
tol = sum(me_dis(1:num))/real(num)
!----------------------------------------------------------------------------------------
! 근접 절점 봉합 - 거리가 1.5 * d 이하이면 점과 점을 연결하여 봉합
cri = 1.5
near_check = 0
if ( wi == 2) cri = 1.5 
if (cross_num == 0) then !<- 임시로
do i = 1, pbn
	
	near_check_recheck = near_check
	near_check = 0
	check = .FALSE.
	
	if (near_num == 0) then
	    check = .TRUE.
	else
	    if (near(1,near_num) /= 0 .AND. near(3,near_num) /= 0) then
	        check = .TRUE.
	    endif
    endif
    cur_d = me_dis(i)
    
    if ( check ) then
        this(1:2) = (/ i-1, i+1 /)
		if (i == 1) then
			k = 2
            this(1) = pbn
		elseif (i == 2) then
			k = 1
		else 
			k = 0
            if (i == pbn) this(2) = 1
        endif
        
		do j = i + 5, temp_pbn - k
            this(3:4) = (/ j-1, j+1 /)
            if (j == 1) this(3) = pbn
            if (j == pbn) this(4) = 1
            !cur_d = 0.25*sum(me_dis(this))
            
			call calc_len(me(:,i), me(:,j), length)
			temp_i = i
			temp_j = j
        
			if (length .LE. cri * cur_d) then
				!write (*,'(2(I,A,2F12.7))') i, ': ', me(:,i), j, ': ', me(:,j)
				change_i_num = 0
				do o = n_b + 1, nn	
                    call calc_len(node(:,o), me(:,i), val(1))
                    call calc_len(node(:,o), me(:,j), val(2))
                    if (val(1) < tol*10e-8) then
					!if ((node(1,o) == me(1,i)) .AND. (node(2,o) == me(2,i))) then
						change_i_num = change_i_num + 1
						if (change_i_num == 1) then
							change_i = o
						elseif (change_i_num == 2) then
							del_i = o
						endif
					!elseif ((node(1,o) == me(1,j)) .AND. (node(2,o) == me(2,j))) then
                    elseif (val(2) < tol*10e-8) then
						del_i = o
					endif
				enddo
				
				!extra_num = 0
				!do o = i, j
				!	extra_num = extra_num + 1
				!enddo
                call calc_angle( me(:,this(1)), me(:,i), me(:,j), alpha(1) )
				call calc_angle( me(:,j), me(:,i), me(:,this(2)), alpha(2) )
				call calc_angle( me(:,this(3)), me(:,j), me(:,i), beta(1) )
				call calc_angle( me(:,i), me(:,j), me(:,this(4)), beta(2) )
				extra_num = j-i+1
				if ( ((extra_num/2) * 2 == extra_num) .AND. (extra_num >= 4) ) then !<- 짝수라면...
                    if ( (alpha(1) > cri) .AND. (alpha(1) < 180.0+cri) .AND. (beta(1) > cri) .AND. (beta(1) < 180.0+cri) .and. (alpha(2) > cri) .AND. (alpha(2) < 180.0+cri) .AND. (beta(2) > cri) .AND. (beta(2) < 180.0+cri) ) then
					    if ( near_num > 1 ) then
						    !extra_num = 0
						    !do o = i, j
						    !	extra_num = extra_num + 1
						    !enddo	
                            extra_num = j-i+1
                        endif
					    ! 근접 봉합 node번호 check
					    near_num = near_num + 1
					    near(1,near_num) = change_i
					    near(2,near_num) = del_i
					    !write (*,*) 'near info: ', near_num, near(1:2,near_num)
					    exit
                    endif
				endif
			endif
		enddo
	elseif ((near(1,near_num) /= 0) .AND. (near(3,near_num) == 0)) then
		
		if (extra_num == 4) then
			near(3,near_num) = change_i
			near(4,near_num) = del_i
            temp_pbn = temp_j
		else
			temp_j = temp_j - 1
			call calc_len(me(:,i), me(:,temp_j), length)
			!this = (/ i-1, i, temp_j-1, temp_j /)
   !         if (i == pbn) this(1) = 1
   !         if (temp_j == pbn) this(3) = 1
   !         cur_d = 0.25*sum(me_dis(this))
            
			change_i_num = 0
			do o = n_b + 1, nn	
                call calc_len(node(:,o), me(:,i), val(1))
                call calc_len(node(:,o), me(:,temp_j), val(2))
                if (val(1) < tol*10e-8) then
				!if ((node(1,o) == me(1,i)) .AND. (node(2,o) == me(2,i))) then
					change_i_num = change_i_num + 1
					if (change_i_num == 1) then
						change_i = o
					elseif (change_i_num == 2) then
						del_i = o
					endif
				!elseif ((node(1,o) == me(1,temp_j)) .AND. (node(2,o) == me(2,temp_j))) then
                elseif (val(2) < tol*10e-8) then
					del_i = o
				endif
			enddo

			if (length .LE. cri * cur_d) then
				if (i - 1 == temp_i) then
					ne = ne + 1
					ele(:,ne) = (/ ne, near(1,near_num), change_i, del_i, near(2,near_num) /)
				else
					ne = ne + 1
					ele(:,ne) = (/ ne, ele(3,ne-1), change_i, del_i, ele(4,ne-1) /)
				endif
				extra_num = extra_num - 2
			else 
				if (i - 1 == temp_i) then
					near(3:4,near_num) = near(1:2,near_num)
				else 
					near(3:4,near_num) = ele(3:4,ne)
                endif
                temp_pbn = temp_j
			endif
		endif			
	endif
enddo

endif !<- 임시로

!if (near_num /= 0) then
!    write (*,*) "근접 하여 봉합된 절점 수"
!    write (*,'(A,I3)') "near_num = ", near_num
!    do i = 1, near_num
!	    write (*,*) near(1:4,i)
!    enddo
!    write(*,*) "----------------------------------------------------------------"
!endif

!========================================================================================
end subroutine