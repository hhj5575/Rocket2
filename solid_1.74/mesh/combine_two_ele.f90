!****************************************************************************
!
!  PROGRAM : Combines two elements
!
!  PURPOSE : 6개 이상 연결된 절점의 분할
!
!****************************************************************************
Subroutine combines_two_ele(anes, pre_nn, nn, ne, node, ele, iter_check)

implicit none

integer, intent(in) :: anes, pre_nn
integer, intent(inout) :: ne, nn, ele(5,anes)
real, intent(inout) :: node(2,anes)
logical, intent(inout) :: iter_check

integer :: i, j, k, ele_sum, temp_ele(20), l_ele, u_ele, iter
integer :: temp_kl, temp_ku, temp_kl_r, temp_ku_r
logical :: check

iter = 0
do 
    iter = iter + 1
    if (iter > nn - pre_nn) then
        write (*,*) 'Error: subroutine combines_two_ele'
        iter_Check = .FALSE.
        exit
    endif    
    check = .TRUE.
    do i = pre_nn+1, nn
	    ele_sum = 0
	    do j = 1, ne
		    do k = 2, 5
			    if ( i == ele(k,j) ) then
				    ele_sum = ele_sum + 1
				    temp_ele(ele_sum) = j
			    endif
		    enddo
	    enddo
	
	    if ( ele_sum == 2 ) then
		
		    l_ele = temp_ele(1); u_ele = temp_ele(2)
		    if ( temp_ele(2) < temp_ele(1) ) then
			    l_ele = temp_ele(2)
			    u_ele = temp_ele(1)
		    endif
		    !write (*,*) i, l_ele, u_ele
		    !write (*,*) ele(1:5,l_ele)
		    !write (*,*) ele(1:5,u_ele)

		    do j = 2, 5
			    if ( ele(j,l_ele) == i ) then
				    temp_kl = j
			    endif
			    if ( ele(j,u_ele) == i ) then
				    temp_ku = j
			    endif
		    enddo

		    temp_kl = temp_kl + 2;  temp_ku = temp_ku + 2
		    if ( temp_kl >= 6 ) temp_kl = temp_kl - 4
		    if ( temp_ku >= 6 ) temp_ku = temp_ku - 4
		    temp_kl_r = temp_kl + 1;  temp_ku_r = temp_ku + 1
		    if ( temp_kl_r == 6 ) temp_kl_r = 2
		    if ( temp_ku_r == 6 ) temp_ku_r = 2
		
		    ele(:,l_ele) = (/ l_ele, ele(temp_kl,l_ele), ele(temp_kl_r,l_ele), ele(temp_ku,u_ele), ele(temp_ku_r,u_ele) /)
						
		    do j = i, nn - 1
			    node(:,j) = node(:,j+1)
		    enddo
		    nn = nn - 1
		
		    do j = 1, ne 
			    do k = 2, 5
				    if ( ele(k,j) >= i ) then
					    ele(k,j) = ele(k,j) - 1
				    endif
			    enddo
		    enddo
		    do j = u_ele, ne - 1
			    ele(2:5,j) = ele(2:5,j+1)
		    enddo
		    ne = ne - 1
            
            check = .FALSE.
	        exit
	    !elseif ( ele_sum >= 9 ) then
	    !    iter_check = .FALSE.
	    !    exit
	    endif
    enddo
    if (check) exit
enddo

end subroutine combines_two_ele