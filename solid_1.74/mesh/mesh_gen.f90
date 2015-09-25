!****************************************************************************
!
!  PROGRAM : Mesh_gen
!
!  PURPOSE : Mesh 생성 
!
!****************************************************************************
Subroutine mesh_gen( temp_bs, in_bs_num, in_b, d, anes, nn, node, ne, ele, wi, iter_check )
                     
implicit none

integer, intent(in) :: temp_bs, anes, wi, in_bs_num(temp_bs)
integer, intent(inout) :: nn, ne, ele(5,anes)
real, intent(in) :: d, in_b(2,temp_bs)
real, intent(inout) :: node(2,anes)
logical, intent(inout) :: iter_check

integer :: i, j, ln, bs, en, pbn, paving_i, ic, n_b, b_num, cross_num, near_num
integer :: div_num, run_div_num, count_div, count_div_num, dpn
integer :: ep(temp_bs), cross(4,temp_bs), near(4,temp_bs), div(10*temp_bs,2), bs_num(temp_bs)
integer :: temp_b_num(temp_bs), end_b_num(6)
real :: max_x, max_y
real :: me(2,temp_bs), temp_b(2,temp_bs), end_b(2,6), b(2,temp_bs)
logical :: ij, ec
integer, allocatable :: temp_div(:,:)

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

bs = temp_bs
bs_num = in_bs_num
b = in_b

max_x = b(1,1)
max_y = b(2,1)
do i = 2, temp_bs
	if (b(1,i) .GT. max_x) then
		max_x = b(1,i)
	elseif (b(2,i) .GT. max_y) then
		max_y = b(2,i)
	endif
enddo

!if ( max_x .GE. max_y ) then
!	ln = nint(max_x / d) 
!else
!	ln = nint(max_y / d)
!endif
ln = 100

!write (*,*) "----------------------------------------------"
!write (*,'(A,I3)') "Size of boundary : ", bs
!do i = 1, bs
!	write (*,*) b(i,:)
!enddo

!========================= 초기 설정 =========================
div_num = 0;  count_div_num = 0;  run_div_num = 0;  
div = 0.;      dpn = 0;  count_div = 0;
ij = .FALSE.
!=============================================================

call first_continue_end_node(anes, bs, bs_num(1:bs), b(:,1:bs), nn, ne, node, ele)
do paving_i = 1, ln
	!write (*,'(I3,A)') paving_i, "번째 paving."
	! Move boundary nodes into inside
    temp_b_num = 0
    temp_b = 0.0
	call paving(temp_bs, d, bs, bs_num(1:bs), b(:,1:bs), pbn, me, en, ep, iter_check)
	if (iter_check == .false.) exit
    
	! Make elements
	n_b = nn
	call make_ele(paving_i, temp_bs, nn, node, n_b, anes, ne, ele, bs, bs_num(1:bs), b(:,1:bs), pbn, me, en, ep, ic, b_num, end_b_num, end_b, iter_check)
    
	if ( ( ic == 1 ) .OR. ( b_num > 6 ) ) then
		if ( b_num <= 6 ) then
			do i = 1, pbn
				nn = nn + 1
				node(:,nn) = me(:,i)
				end_b_num(i) = nn
				end_b(:,i) = node(:,nn)
            enddo
			!write (*,*) "---=== 남은 경계절점과 유동절점이 6이하가 되었습니다.<Element> ===---"
			!write (*,'(A,I3)') "bs : ", b_num
			ij = .TRUE.
			goto 200
		else
			do i = 1, pbn
				nn = nn + 1
				node(:,nn) = me(:,i)
			enddo
		endif
	else
		ij = .TRUE.
		goto 200
	endif
    if (iter_check == .false.) exit
    
	! Angle check
    b_num = 0
	call angle_check(paving_i, anes, temp_bs, pbn, nn, ne, n_b, b_num, me, node, ele, temp_b_num, temp_b, ec, iter_check)

    if (iter_check == .false.) exit
    
    if (pbn > 6) then
	    call con_node_cross(anes, nn, ne, n_b, node, ele, pbn, me(:,1:pbn))
	    if ( ec ) goto 200
    endif
    
	! Cross & near check
	call cross_near(temp_bs, n_b, anes, nn, ne, pbn, node, ele, me, cross_num, near_num, cross, near, wi)
    
	! Reset boundary nodes
	call reset_bound(anes, temp_bs, pbn, n_b, nn, b_num, cross_num, near_num, node, me(:,1:pbn), cross(:,1:cross_num), near(:,1:near_num), temp_b_num, temp_b, div, div_num, dpn)
    
	! Print boundary nodes
	!if (div_num > run_div_num) then
	!	write(*,*) "----------------------------------------"
	!	write (*,*) "분활된 경계 절점"
	!	do i = run_div_num + 1, div_num
	!		write (*,'(I3,A)') i+1, "번째 분활된 경계 절점"
	!		do j = 1, dpn
	!			if ( div(j,1) == i ) then
	!				write (*,'(I5,$)') div(j,2)
	!			endif
	!		enddo
	!		write(*,'(/,A)') "----------------------------------------"
	!	enddo
	!endif
	
	! Check the end of paving
	200 continue
	if ( ij ) then
		do i = 1, b_num
			temp_b(:,i) = end_b(:,i)
		enddo
	endif
	ij = .FALSE.

	!write (*,*) b_num
	if ( b_num <= 6 ) then
		!write (*,'(I3,A)') run_div_num+1, "번째 경계면 Paving 종료"
		!write(*,*) "----------------------------------------------------------------"
        call clock_check(anes, b_num, node, ne, temp_b_num(1:b_num), temp_b(:,1:b_num), ele, iter_check)
	    if (iter_check == .false.) exit
        
		call end_check( anes, b_num, temp_b_num(1:b_num), temp_b(:,1:b_num), ne, ele, nn, node )
		
		! Soften inner nodes
		if ( wi == 2 ) call inner_soft_end(anes, n_b, nn, ne, node, ele, iter_check)
		if (iter_check == .false.) exit
        
		if (run_div_num == div_num) then
			exit
		else
			run_div_num = run_div_num + 1
			bs = 0 
			do i = 1, dpn
				if (div(i,1) == run_div_num) then
					bs = bs + 1
					do j = 1, nn
						if ( div(i,2) == j ) then
							bs_num(bs) = j
							b(:,bs) = node(:,j)
						!	write (*,*) bs, b(bs,1:3)
						endif
					enddo
				endif
			enddo
			!write(*,*) "----------------------------------------------------------------"		
			!write (*,'(I3,A)') run_div_num+1, "번째 분할된 경계면 Paving 시작"
			!do i = 1, bs
			!	write (*,*) b(i,1:3)
			!enddo
			!write(*,*) "----------------------------------------------------------------"
		endif
	else
		bs = b_num
		do i = 1, bs 
            bs_num(i) = temp_b_num(i)
			b(:,i) = temp_b(:,i)
		enddo
	endif
	
	! Soften inner nodes
	call inner_soft(anes, n_b, nn, ne, bs, bs_num(1:bs), node, ele, iter_check)
	if (iter_check == .false.) exit
    
	! Soften boundary nodes
	! 1. 각도 재검사 (교차와 근접 절점의 봉합 후 생기는 예각 수정)
    allocate ( temp_div(dpn,2) )
    temp_div = div(1:dpn,:)
	call angle_recheck(anes, temp_bs, b_num, bs, nn, ne, bs_num, b, temp_b_num, temp_b, node, ele, dpn, temp_div, ec)
    div(1:dpn,:) = temp_div(1:dpn,:)
    deallocate ( temp_div )
	if ( ec ) goto 200
	
	! 2. 튀어나온 절점 안으로 넣기 
	call bound_soft(anes, bs, bs_num(1:bs), b(:,1:bs), d, ne, nn, ele, node)
	
	! 3. 경계 절점의 유연화
	call len_smooth(bs, bs_num(1:bs), b(:,1:bs), nn, node(:,1:nn))

	! 4. 연속 끌점점의 처리 
	call continue_end_smooth(anes, temp_bs, bs, bs_num(1:bs), b(:,1:bs), nn, ne, node, ele)
	
	!write(*,*) "boundary point"
	!write(*,'(A,I3)') "b_size = ", bs
	!do i = 1, bs 
	!	write (*,*) b(i,:)
	!enddo
	!write(*,*) "----------------------------------------------------------------"
	
enddo ! <- 전체 루프문 종료

end subroutine
