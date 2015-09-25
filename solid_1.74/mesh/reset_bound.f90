!****************************************************************************
!
!  PROGRAM : Reset bound
!
!  PURPOSE : 교차와 근접 검사로 나누어진 영역 분할
!
!****************************************************************************
Subroutine reset_bound(anes, temp_bs, pbn, n_b, nn, b_num, cross_num, near_num, node, me, cross, near, temp_b_num, temp_b, div, div_num, dpn)

implicit none

integer, intent(in) :: anes, pbn, n_b, nn, cross_num, near_num, temp_bs
integer, intent(in) :: cross(4,cross_num), near(4,near_num)
integer, intent(inout) :: b_num, div_num, dpn, div(10*temp_bs,2), temp_b_num(temp_bs)
real, intent(in) :: node(2,anes), me(2,pbn)
real, intent(inout) :: temp_b(2,temp_bs)

integer :: i, j, reset_b_size, count, reset_b_num(temp_bs)
real :: reset_b(2,temp_bs)

! 경계선 재설정 : 임시 경계로 옮김 -> 크로스 확인 -> 근접 확인 -> 마무리 경계 절점 설정
! 1. 임시 경계로 옮김
reset_b_num = 0
reset_b = 0.0
do i = 1, pbn 
	do j = n_b + 1, nn 
		if ((node(1,j) == me(1,i)) .AND. (node(2,j) == me(2,i))) then
			reset_b_num(i) = j
			exit
		endif
	enddo

	reset_b(:,i) = me(:,i)
enddo
reset_b_size = pbn

! 2. 크로스 확인
if (cross_num .NE. 0) then
	do i = 1, (cross_num/2) + 1

		b_num = 0

		if (i == 1) then

			do j = 1, reset_b_size
				if ((reset_b_num(j) <= cross(1,i)) .OR. (reset_b_num(j) >= cross(2,i))) then
					count = count + 1
                    temp_b_num(count) = reset_b_num(j)
					temp_b(:,count) = reset_b(:,j)
				endif
			enddo
			b_num = count
			
		elseif (i == cross_num/2 + 1) then

			div_num = div_num + 1
			do j = 1, reset_b_size
				if ((reset_b_num(j) >= cross(3,i-1)) .AND. (reset_b_num(j) <= cross(4,i-1))) then
					count = count + 1
					dpn = dpn + 1
					div(dpn,1) = div_num
					div(dpn,2) = reset_b_num(j) !, reset_b(j,2), reset_b(j,3) /)
				endif
			enddo
						
		else

			div_num = div_num + 1
			do j = 1, reset_b_size
				if ((reset_b_num(j) >= cross(3,i-1) .AND. reset_b_num(j) <= cross(1,i)) .OR. (reset_b_num(j) >= cross(2,i) .AND. reset_b_num(j) <= cross(4,i-1))) then
					count = count + 1
					dpn = dpn + 1
					div(dpn,1) = div_num
					div(dpn,2) = reset_b_num(j)!, reset_b(j,2), reset_b(j,3) /)
				endif	
			enddo
			
		endif

	enddo
endif

!write (*,*) 'reset_b_size:', reset_b_size
!write (*,*) 'reset_b:', nint(reset_b(1,1:reset_b_size))
!write (*,*) near(3:4,near_num)
! 3. 근접 확인
if (near_num .NE. 0) then
	do i = 1, near_num + 1

		count = 0

		if (i == 1) then

			do j = 1, reset_b_size
				if ((reset_b_num(j) <= near(1,i)) .OR. (reset_b_num(j) >= near(2,i))) then
					count = count + 1
					temp_b_num(count) = reset_b_num(j)
					temp_b(:,count) = reset_b(:,j)
				endif
			enddo
			b_num = count

		elseif (i == near_num + 1) then

			div_num = div_num + 1
			do j = 1, reset_b_size
				if ((reset_b_num(j) >= near(3,i-1)) .AND. (reset_b_num(j) <= near(4,i-1))) then
					count = count + 1
					dpn = dpn + 1
					div(dpn,1) = div_num
					div(dpn,2) = reset_b_num(j)!, reset_b(j,2), reset_b(j,3) /)
				endif
			enddo
			
		else

			div_num = div_num + 1
			do j = 1, reset_b_size
				if ((reset_b_num(j) >= near(3,i-1) .AND. reset_b_num(j) <= near(1,i)) .OR. (reset_b_num(j) >= near(2,i) .AND. reset_b_num(j) <= near(4,i-1))) then
					count = count + 1
					dpn = dpn + 1
					div(dpn,1) = div_num
					div(dpn,2) = reset_b_num(j)!, reset_b(j,2), reset_b(j,3) /)
				endif	
			enddo
			
		endif

	enddo

endif

! 4. 마무리 경계 절점 설정
if ((cross_num == 0) .AND. (near_num == 0)) then
	do i = 1, reset_b_size
        temp_b_num(i) = reset_b_num(i)
        temp_b(:,i) = reset_b(:,i)
	enddo
	b_num = reset_b_size 
endif

! 경계선 재설정 종료
!========================================================================================
end subroutine