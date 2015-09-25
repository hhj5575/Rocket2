!****************************************************************************
!
!  PROGRAM : End_check
!
!  PURPOSE : 남은 절점의 수가 6개 이하이면 종료
!
!****************************************************************************
Subroutine end_check(anes, b_num, temp_b_num, temp_b, ne, ele, nn, node)

implicit none

integer, intent(in) :: b_num, anes, temp_b_num(b_num)
integer, intent(inout) :: ne, nn, ele(5,anes)
real, intent(in) :: temp_b(2,b_num)
real, intent(inout) :: node(2,anes)

if (b_num <= 6) then
	
	if (b_num == 4) then !<= 4개 보다 작으면...

		ne = ne + 1
		ele(:,ne) = (/ ne, temp_b_num(1), temp_b_num(2), temp_b_num(3), temp_b_num(4) /)
				
	elseif (b_num == 6) then !<= 6개 보다 작으면...
		
		!do i = 1, 6 
		!	write (*,*) temp_b(i,:)
		!enddo
		call make_end_ele( anes, temp_b_num(1:6), temp_b(:,1:6), ne, ele, nn, node )
	endif

endif

end subroutine end_check
!======================================================================================
!======================================================================================
! Last Element생성
Subroutine make_end_ele(anes, temp_b_num, temp_b, ne, ele, nn, node)

implicit none

integer, intent(in) :: anes, temp_b_num(6)
real, intent(in) :: temp_b(2,6)
integer, intent(inout) :: ne, nn, ele(5,anes)
real, intent(inout) :: node(2,anes)

integer :: i, j, k, q, num, sum_case, elem(4), this(2)
real :: min_dis
real :: dis(3)
logical :: check

call calc_len(temp_b(:,1), temp_b(:,4), dis(1))
call calc_len(temp_b(:,2), temp_b(:,5), dis(2))
call calc_len(temp_b(:,3), temp_b(:,6), dis(3))

min_dis = dis(1)
do i = 2, 3
	if (dis(i) .LT. min_dis) then
		min_dis = dis(i)
	endif
enddo

if (min_dis == dis(1)) then

	ne = ne + 1
	ele(:,ne) = (/ ne, temp_b_num(1), temp_b_num(2), temp_b_num(3), temp_b_num(4) /)
	
	ne = ne + 1
    ele(:,ne) = (/ ne, temp_b_num(4), temp_b_num(5), temp_b_num(6), temp_b_num(1) /)

elseif (min_dis == dis(2)) then

	ne = ne + 1
    ele(:,ne) = (/ ne, temp_b_num(2), temp_b_num(3), temp_b_num(4), temp_b_num(5) /)
	
	ne = ne + 1
    ele(:,ne) = (/ ne, temp_b_num(5), temp_b_num(6), temp_b_num(1), temp_b_num(2) /)
    
elseif (min_dis == dis(3)) then
	
	ne = ne + 1
    ele(:,ne) = (/ ne, temp_b_num(3), temp_b_num(4), temp_b_num(5), temp_b_num(6) /)
	
	ne = ne + 1
    ele(:,ne) = (/ ne, temp_b_num(6), temp_b_num(1), temp_b_num(2), temp_b_num(3) /)
    
endif

check = .false.
do q = 1, 2
    num = 1
    if (q == 2) num = 0 
    
    do i = ne-num, ne
        elem = ele(2:5,i)
        call check_diamond_element(elem, nn, node(:,1:nn), sum_case, check)
        if (check) then
            ! sum_case = 1-sum of node 2, 4, 2-sum of node 1, 3
            if (sum_case == 1) then
                this = (/ elem(2), elem(4) /)
                if (elem(2) > elem(4)) this = (/ elem(4), elem(2) /)
            elseif (sum_case == 2) then
                this = (/ elem(1), elem(3) /)
                if (elem(1) > elem(3)) this = (/ elem(3), elem(1) /)
            else
                stop 'error(end_check): check is true but sum_case is zero'
            endif
            
            node(:,this(1)) = (node(:,this(1))+node(:,this(2)))*0.5
            do j = this(2), nn-1
                node(:,j) = node(:,j+1)
            enddo
            nn = nn - 1
            do j = 1, ne
                do k = 2, 5
                    if (ele(k,j) == this(2)) then
                        ele(k,j) = this(1)
                    elseif (ele(k,j) > this(2)) then
                        ele(k,j) = ele(k,j) - 1
                    endif
                enddo
            enddo
            if (q == 1 .and. i == ne-1) ele(2:5,ne-1) = ele(2:5,ne)
            ne = ne - 1
            
            exit
        endif
    enddo ! end i
enddo ! end q

end subroutine

! Element생성 종료
!======================================================================================

subroutine check_diamond_element(elem, nn, node, sum_case, check)

implicit none

integer, intent(in) :: nn, elem(4)
real, intent(in) :: node(2,nn)
integer, intent(out) :: sum_case
logical, intent(out) :: check

integer :: i, this(3)
real :: ang(4), temp_node(2,3), cri_ang(2)

cri_ang = (/ 40., 140. /)
check = .false.
sum_case = 0
do i = 1, 4
    this = (/ i-1, i, i+1 /)
    if (i == 1) this(1) = 4
    if (i == 4) this(3) = 1
    
    temp_node = node(:,elem(this))
    call calc_angle(temp_node(:,1), temp_node(:,2), temp_node(:,3), ang(i))
enddo

if (ang(1) <= cri_ang(1) .and. ang(2) >= cri_ang(2) .and. ang(3) <= cri_ang(1) .and. ang(4) >= cri_ang(2)) then
    check = .true.
    sum_case = 1
elseif (ang(1) >= cri_ang(2) .and. ang(2) <= cri_ang(1) .and. ang(3) >= cri_ang(2) .and. ang(4) <= cri_ang(1)) then
    check = .true.
    sum_case = 2
endif
  

end subroutine check_diamond_element