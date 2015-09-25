!****************************************************************************
!
!  PROGRAM : first_continue_end_node
!
!****************************************************************************
Subroutine first_continue_end_node(anes, bs, bs_num, b, nn, ne, node, ele)

implicit none

integer, intent(in) :: anes
integer, intent(inout) :: bs, bs_num(bs), nn, ne, ele(5,anes)
real, intent(inout) :: b(2,bs), node(2,anes)

integer :: i, this(4), num
real :: alpha, beta, length(2)

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

do 
    num = 0
    do i = 1, bs
	    if ( i == 1 ) then
		    this = (/ bs-1, bs, 1, 2 /)
        elseif ( i == 2 ) then
            this = (/ bs, 1, 2, 3 /)
        elseif ( i == bs ) then
            this = (/ i-2, i-1, i, 1 /)
        else
            this = (/ i-2, i-1, i, i+1 /)
	    endif
	
	    call calc_angle(b(:,this(1)), b(:,this(2)), b(:,this(3)), alpha)
	    call calc_angle(b(:,this(2)), b(:,this(3)), b(:,this(4)), beta)
	    if ( ( alpha < 125. ) .AND. ( beta < 125. ) .AND. ( alpha + beta > 180. ) .AND. ( alpha + beta < 220. ) ) then
		    call calc_len(b(:,this(1)), b(:,this(2)), length(1))
		    call calc_len(b(:,this(3)), b(:,this(4)), length(2))
            
		    if ( length(1) <= length(2) ) then
			    nn = nn + 1
			    node(:,nn) = ( b(:,this(1)) + b(:,this(3)) ) * 0.5
			
			    ne = ne + 1
			    ele(:,ne) = (/ ne, bs_num(this(1)), bs_num(this(2)), bs_num(this(3)), nn /)

			    bs_num(this(2)) = nn
			    b(:,this(2)) = node(:,nn)
		    elseif ( length(1) > length(2) ) then
			    nn = nn + 1
			    node(:,nn) = ( b(:,this(2)) + b(:,this(4)) ) * 0.5
			
			    ne = ne + 1
			    ele(:,ne) = (/ ne, bs_num(this(2)), bs_num(this(3)), bs_num(this(4)), nn /)

    		    bs_num(this(3)) = nn
			    b(:,this(3)) = node(:,nn)
            endif
            num = 1
            exit
	    endif
    enddo		
    
    if (num == 0) exit
enddo

end subroutine first_continue_end_node


!****************************************************************************
!
!  PROGRAM : Continue_end_smooth
!
!  PURPOSE : 연속되는 끝절점의 처리
!
!****************************************************************************
Subroutine continue_end_smooth(anes, temp_bs, bs, bs_num, b, nn, ne, node, ele)

implicit none

integer, intent(in) :: anes, temp_bs
integer, intent(inout) :: bs, bs_num(temp_bs), nn, ne, ele(5,anes)
real, intent(inout) :: b(2,temp_bs), node(2,anes)

integer :: i, i_l, i_r,i_l_l
real :: alpha, beta, length(2), vector(2,2)

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

do i = 1, bs
	if ( i == 1 ) then
		i_l = bs;  i_l_l = bs - 1;  i_r = 2
	elseif ( i == 2 ) then
		i_l = 1;  i_l_l = bs;  i_r = 3
	elseif ( i == bs ) then
		i_l = i - 1;  i_l_l = i - 2;  i_r = 1
	else
		i_l = i - 1;  i_l_l = i - 2;  i_r = i + 1
	endif
	
	call calc_angle(b(:,i_l_l), b(:,i_l), b(:,i), alpha)
	call calc_angle(b(:,i_l), b(:,i), b(:,i_r), beta)
	
	if ( ( alpha < 125. ) .AND. ( beta < 125. ) .AND. ( alpha + beta > 190. ) .AND. ( alpha + beta < 220. ) ) then
		call calc_len(b(:,i_l), b(:,i_l_l), length(1))
		call calc_len(b(:,i), b(:,i_r), length(2))

		if ( length(1) <= length(2) ) then

			nn = nn + 1
			node(:,nn) = ( b(:,i_l_l) + b(:,i) ) * 0.5
			
			ne = ne + 1
			ele(:,ne) = (/ ne, bs_num(i_l_l), bs_num(i_l), bs_num(i), nn /)

			vector(:,1) = b(:,i_l) - b(:,i_l_l)
			vector(:,2) = b(:,i_l) - b(:,i)

			b(:,i_l_l) = b(:,i_l_l) + (vector(:,1) * 0.15)
			b(:,i) = b(:,i) + (vector(:,2) * 0.15)

			node(:,bs_num(i_l_l)) = b(:,i_l_l)
			node(:,bs_num(i)) = b(:,i)
			
			bs_num(i_l) = nn
			b(:,i_l) = node(:,nn)
							
		elseif ( length(1) > length(2) ) then

			nn = nn + 1
			node(:,nn) = ( b(:,i_l) + b(:,i_r) ) * 0.5
			
			ne = ne + 1
			ele(:,ne) = (/ ne, bs_num(i_l), bs_num(i), bs_num(i_r), nn /)

			vector(:,1) = b(:,i) - b(:,i_l)
			vector(:,2) = b(:,i) - b(:,i_r)

			b(:,i_l) = b(:,i_l) + (vector(:,1) * 0.15)
			b(:,i_r) = b(:,i_r) + (vector(:,2) * 0.15)

			node(:,bs_num(i_l)) = b(:,i_l)
			node(:,bs_num(i_r)) = b(:,i_r)

			bs_num(i) = nn
			b(:,i) = node(:,nn)
			
		endif

	endif
enddo		

end subroutine continue_end_smooth