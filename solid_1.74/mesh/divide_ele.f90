!****************************************************************************
!
!  PROGRAM : Divide element
!
!  PURPOSE : 6개 이상 연결된 절점의 분할
!
!****************************************************************************
! 1개의 절점에서 6개 이상의 Element가 생성될 때...
Subroutine divide_over_ele(anes, nn, ne, node, ele, pre_nn, pre_ne, iter_check)

implicit none

integer, intent(in) :: anes, pre_nn, pre_ne
integer, intent(inout) :: ne, nn, ele(5,anes)
real, intent(inout) :: node(2,anes)
logical, intent(inout) :: iter_check

integer :: i, j, k, o, j_r, num, count, ele_num, temp_i, b_num, ele_con_num, ele_con, ele_i, temp_j
integer :: div_ele(10,5), ele_right_node(10,3), temp_ele_right_node(10,3), b_u(11), b_l(11), ele_u(5), ele_l(5)
integer :: ele_nodes(16), this(2), elem(4), pre_num, div_ele_num, iter, fn, temp_nodes(16), ce
real :: min_dis, sum_x(2), sum_y(2), min_ang, ang, vec(2)
real :: dis(8) 

do 
    div_ele_num = 0
    do i = pre_nn+1, nn
	    ele_num = 0
	    do j = pre_ne, ne
		    do k = 2, 5
			    if (ele(k,j) == i) then
				    ele_num = ele_num + 1
				    div_ele(ele_num,1:5) = ele(:,j)
				    if (k == 5) then
					    ele_right_node(ele_num,1) = ele(2,j)
					    ele_right_node(ele_num,2) = 2
				    else
					    ele_right_node(ele_num,1) = ele(k+1,j)
					    ele_right_node(ele_num,2) = k + 1
				    endif
			    endif
		    enddo
	    enddo
	
	    do j = 1, ele_num
		    ele_con_num = 0
		    do k = 1, ne
			    do o = 2, 5
				    if (ele(o,k) == ele_right_node(j,1)) then
					    ele_con_num = ele_con_num + 1
				    endif
			    enddo
		    enddo
		    ele_right_node(j,3) = ele_con_num
        enddo
    
        b_u = 0
        b_l = 0
        !if (ele_num >= 6) write (*,*) i, ' node ele_num:', ele_num
	    if (ele_num == 6) then	
		    ele_con = 0
		    do j = 1, 3
			    if ( (ele_right_node(j,3) < 5) .AND. (ele_right_node(j+3,3) < 5) ) then
                    ele_con = j
                    exit
                endif
		    enddo
		    if (ele_con == 0) then
			    ele_num = 0
			    goto 300
		    endif

		    do j = 1, ele_num
			    call calc_len(node(:,i), node(:,ele_right_node(j,1)), dis(j))
		    enddo

		    temp_ele_right_node = ele_right_node
		
		    min_dis = 0
		    do j = 1, ele_num
			    j_r = j + 3
			    if (j_r >= 7) j_r = j_r - 6

			    if ( (temp_ele_right_node(j_r,3) .LE. 4) .AND. (temp_ele_right_node(j,3) .LE. 4) ) then
				
				    if ( min_dis == 0 ) then
					    min_dis = dis(j)
					    temp_i = j
				    else
					    if ( dis(j) < min_dis ) then
						    min_dis = dis(j)
						    temp_i = j
					    endif
				    endif
			    endif
            enddo
            temp_i = ele_con
        
		    b_num = 1
		    b_u(b_num) = temp_ele_right_node(temp_i,1)
		    ele_u(1) = temp_i

		    do j = 1, 3
			
			    if (j == 1) then

				    if (ele_right_node(temp_i,2) == 4) then
					    b_u(b_num+1) = div_ele(temp_i, 5)
					    b_u(b_num+2) = div_ele(temp_i, 2)
					    b_num = b_num + 2					
				    elseif (ele_right_node(temp_i,2) == 5) then
					    b_u(b_num+1) = div_ele(temp_i, 2)
					    b_u(b_num+2) = div_ele(temp_i, 3)
					    b_num = b_num + 2
				    else
					    b_u(b_num+1) = div_ele(temp_i, ele_right_node(temp_i,2)+1)
					    b_u(b_num+2) = div_ele(temp_i, ele_right_node(temp_i,2)+2)
					    b_num = b_num + 2
				    endif
				
			    elseif (j == 2) then
				
				    do k = 1, 6 
					    do o = 2, 5
						    if ((div_ele(k,o) == b_u(3)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k	
							    ele_u(2) = k
							    if (o == 4) then
								    b_u(b_num+1) = div_ele(k,5)
								    b_u(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_u(b_num+1) = div_ele(k,2)
								    b_u(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_u(b_num+1) = div_ele(k,o+1)
								    b_u(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 5) then
						    exit
					    endif
				    enddo

			    else
				    do k = 1, 6 
					    do o = 2, 5
						    if ((div_ele(k,o) == b_u(5)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k
							    ele_u(3) = k
							    if (o == 4) then
								    b_u(b_num+1) = div_ele(k,5)
								    b_u(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_u(b_num+1) = div_ele(k,2)
								    b_u(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_u(b_num+1) = div_ele(k,o+1)
								    b_u(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 7) then
						    exit
					    endif
				    enddo
			
			    endif

		    enddo
		
		    b_num = 1
		    b_l(1) = b_u(7)
		    do j = 1, 3
			
			    if (j == 1) then

				    do k = 1, 6 
					    do o = 2, 5
						    if ((div_ele(k,o) == b_l(1)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k	
							    ele_l(1) = k
							    if (o == 4) then
								    b_l(b_num+1) = div_ele(k,5)
								    b_l(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_l(b_num+1) = div_ele(k,2)
								    b_l(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_l(b_num+1) = div_ele(k,o+1)
								    b_l(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 3) then
						    exit
					    endif
				    enddo
				
			    elseif (j == 2) then

				    do k = 1, 6 
					    do o = 2, 5
						    if ((div_ele(k,o) == b_l(3)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k	
							    ele_l(2) = k
							    if (o == 4) then
								    b_l(b_num+1) = div_ele(k,5)
								    b_l(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_l(b_num+1) = div_ele(k,2)
								    b_l(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_l(b_num+1) = div_ele(k,o+1)
								    b_l(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 5) then
						    exit
					    endif
				    enddo

			    else
				
				    do k = 1, 6 
					    do o = 2, 5
						    if ((div_ele(k,o) == b_l(5)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k
							    ele_l(3) = k
							    if (o == 4) then
								    b_l(b_num+1) = div_ele(k,5)
								    b_l(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_l(b_num+1) = div_ele(k,2)
								    b_l(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_l(b_num+1) = div_ele(k,o+1)
								    b_l(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 7) then
						    exit
					    endif
				    enddo

			    endif

		    enddo
		    sum_x = 0.
		    sum_y = 0.
		    do j = 1, b_num
                if (b_u(j) <= 0 .OR. b_u(j) > nn .OR. b_l(j) <= 0 .OR. b_l(j) > nn) then
                    iter_check = .FALSE.
                    exit
                endif
			    sum_x(1) = sum_x(1) + node(1,b_u(j))
			    sum_y(1) = sum_y(1) + node(2,b_u(j))
			    sum_x(2) = sum_x(2) + node(1,b_l(j))
			    sum_y(2) = sum_y(2) + node(2,b_l(j))
		    enddo
            if (iter_check == .FALSE.) exit
		    node(1,i) = sum_x(1) / 7.
		    node(2,i) = sum_y(1) / 7.
		    nn = nn + 1
		    node(1,nn) = (sum_x(2) + node(1,i)) / 8.
		    node(2,nn) = (sum_y(2) + node(2,i)) / 8.

		    do j = 1, 3
                if (ele_l(j) == 0 .or. ele_l(j) > 10) then
                    iter_check = .false.
                    exit
                endif
			    do k = 2, 5
				    if (div_ele(ele_l(j),k) == i) then
					    ele(k,div_ele(ele_l(j),1)) = nn
				    endif
			    enddo
		    enddo
            if (iter_check == .FALSE.) exit
        
		    ne = ne + 1
		    ele(1,ne) = ne
		    ele(2,ne) = b_u(1)
		    ele(3,ne) = i
		    ele(4,ne) = b_u(7)
		    ele(5,ne) = nn
            div_ele_num = 1
            exit
		    300 continue
		    !write (*,*) "-------------------------------------------------------------"
		    !if (ele_num == 0) then
		    !	write (*,'(I5,A)') i, " node에서 6개의 Element로 연결되었지만 분할 할 수 없습니다." 
		    !else
		    !	write (*,'(I5,A)') i, " node에서 6개의 Element로 연결되어 분할하였습니다." 
		    !	write (*,'(A,I5)') " 생성된 node : ", nn
		    !	write (*,'(A,I5)') " 생성된 Element : ", ne
		    !	write (*,*) ele(2:5,ne)  
		    !endif
		    !write (*,*) "-------------------------------------------------------------"
    !===========================================================================================

    !===========================================================================================
	    elseif (ele_num == 7) then
		    ele_con = 0
		    ele_con_num = 0 
		    temp_ele_right_node = ele_right_node
		    do ele_i = 1, 7
			
			    do j = 1, ele_num
				    if (temp_ele_right_node(j,3) .LE. 4) then
					    min_dis = dis(j)
					    temp_i = j
					    exit
				    endif
			    enddo
			    do j = temp_i+1, ele_num
				    if ((dis(j) .LT. min_dis) .AND. (temp_ele_right_node(j,3) .LE. 4)) then
					    min_dis = dis(j)
					    temp_i = j
				    endif
			    enddo
		    
		        if ( temp_i > ele_num .OR. temp_i < 1 ) then
		            iter_check = .FALSE.
                    exit
	            endif
	        
			    b_num = 1
			    b_u(b_num) = temp_ele_right_node(temp_i,1)
			    ele_u(1) = temp_i
			
			    do j = 1, 4
				
				    if (j == 1) then

					    if (ele_right_node(temp_i,2) == 4) then
						    b_u(b_num+1) = div_ele(temp_i, 5)
						    b_u(b_num+2) = div_ele(temp_i, 2)
						    b_num = b_num + 2					
					    elseif (ele_right_node(temp_i,2) == 5) then
						    b_u(b_num+1) = div_ele(temp_i, 2)
						    b_u(b_num+2) = div_ele(temp_i, 3)
						    b_num = b_num + 2
					    else
						    b_u(b_num+1) = div_ele(temp_i, ele_right_node(temp_i,2)+1)
						    b_u(b_num+2) = div_ele(temp_i, ele_right_node(temp_i,2)+2)
						    b_num = b_num + 2
					    endif

				    elseif (j == 2) then

					    do k = 1, 7 
						    do o = 2, 5
							    if ((div_ele(k,o) == b_u(3)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
								    temp_i = k	
								    ele_u(2) = k
								    if (o == 4) then
									    b_u(b_num+1) = div_ele(k,5)
									    b_u(b_num+2) = div_ele(k,2)
									    b_num = b_num + 2					
								    elseif (o == 5) then
									    b_u(b_num+1) = div_ele(k,2)
									    b_u(b_num+2) = div_ele(k,3)
									    b_num = b_num + 2
								    else
									    b_u(b_num+1) = div_ele(k,o+1)
									    b_u(b_num+2) = div_ele(k,o+2)
									    b_num = b_num + 2
								    endif
								    exit
							    endif
						    enddo
						    if (b_num == 5) then
							    exit
						    endif
					    enddo

				    elseif (j == 3) then
					
					    do k = 1, 7 
						    do o = 2, 5
							    if ((div_ele(k,o) == b_u(5)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
								    temp_i = k
								    ele_u(3) = k
								    if (o == 4) then
									    b_u(b_num+1) = div_ele(k,5)
									    b_u(b_num+2) = div_ele(k,2)
									    b_num = b_num + 2					
								    elseif (o == 5) then
									    b_u(b_num+1) = div_ele(k,2)
									    b_u(b_num+2) = div_ele(k,3)
									    b_num = b_num + 2
								    else
									    b_u(b_num+1) = div_ele(k,o+1)
									    b_u(b_num+2) = div_ele(k,o+2)
									    b_num = b_num + 2
								    endif
								    exit
							    endif
						    enddo
						    if (b_num == 7) then
							    exit
						    endif
					    enddo
				
					    if (ele_right_node(temp_i,3) .LE. 4) then
						    ele_con = 3
						    temp_j = temp_i
					    else 
						    ele_con = 0
					    endif
				
				    else
					
					    do k = 1, 7 
						    do o = 2, 5
							    if ((div_ele(k,o) == b_u(7)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
								    temp_i = k
								    ele_u(4) = k
								    if (o == 4) then
									    b_u(b_num+1) = div_ele(k,5)
									    b_u(b_num+2) = div_ele(k,2)
									    b_num = b_num + 2					
								    elseif (o == 5) then
									    b_u(b_num+1) = div_ele(k,2)
									    b_u(b_num+2) = div_ele(k,3)
									    b_num = b_num + 2
								    else
									    b_u(b_num+1) = div_ele(k,o+1)
									    b_u(b_num+2) = div_ele(k,o+2)
									    b_num = b_num + 2
								    endif
								    exit
							    endif
						    enddo
						    if (b_num == 9) then
							    exit
						    endif
					    enddo

					    if (ele_right_node(temp_i,3) .LE. 4) then
						    ele_con = 4
						    b_num = 1
						    b_l(b_num) = b_u(9)
					    elseif (ele_con == 3) then					
						    b_num = 1
						    b_l(b_num) = b_u(7)
						    temp_i = temp_j
					    else
						    do k = 1, ele_num
							    if (temp_ele_right_node(k,1) == b_u(1)) then
								    do o = k, ele_num - 1
									    temp_ele_right_node(k,1:3) = temp_ele_right_node(k+1,1:3)
								    enddo
								    exit
							    endif
						    enddo
						    ele_num = ele_num - 1
					    endif

				    endif

			    enddo

			    if ((ele_con == 3) .OR. (ele_con == 4)) then
				    exit
			    endif

		    enddo
		
		    if ( iter_check == .FALSE. ) exit
		
		    if (ele_num == 0) then

			    goto 400

		    endif

		    do j = ele_con + 1, 7
			
			    if (j == ele_con + 1) then

				    do k = 1, 7 
					    do o = 2, 5
						    if ((div_ele(k,o) == b_l(1)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k	
							    ele_l(1) = k
							    if (o == 4) then
								    b_l(b_num+1) = div_ele(k,5)
								    b_l(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_l(b_num+1) = div_ele(k,2)
								    b_l(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_l(b_num+1) = div_ele(k,o+1)
								    b_l(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 3) then
						    exit
					    endif
				    enddo

			    elseif (j == ele_con + 2) then

				    do k = 1, 7 
					    do o = 2, 5
						    if ((div_ele(k,o) == b_l(3)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k	
							    ele_l(2) = k
							    if (o == 4) then
								    b_l(b_num+1) = div_ele(k,5)
								    b_l(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_l(b_num+1) = div_ele(k,2)
								    b_l(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_l(b_num+1) = div_ele(k,o+1)
								    b_l(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 5) then
						    exit
					    endif
				    enddo

			    elseif (j == ele_con + 3) then
				    do k = 1, 7
					    do o = 2, 5
						    if ((div_ele(k,o) == b_l(5)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k
							    ele_l(3) = k
							    if (o == 4) then
								    b_l(b_num+1) = div_ele(k,5)
								    b_l(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_l(b_num+1) = div_ele(k,2)
								    b_l(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_l(b_num+1) = div_ele(k,o+1)
								    b_l(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 7) then
						    exit
					    endif
				    enddo
				
			    elseif (j == ele_con + 4) then
				
				    do k = 1, 7 
					    do o = 2, 5
						    if ((div_ele(k,o) == b_l(7)) .AND. (div_ele(k,1) .NE. div_ele(temp_i,1))) then
							    temp_i = k
							    ele_l(4) = k
							    if (o == 4) then
								    b_l(b_num+1) = div_ele(k,5)
								    b_l(b_num+2) = div_ele(k,2)
								    b_num = b_num + 2					
							    elseif (o == 5) then
								    b_l(b_num+1) = div_ele(k,2)
								    b_l(b_num+2) = div_ele(k,3)
								    b_num = b_num + 2
							    else
								    b_l(b_num+1) = div_ele(k,o+1)
								    b_l(b_num+2) = div_ele(k,o+2)
								    b_num = b_num + 2
							    endif
							    exit
						    endif
					    enddo
					    if (b_num == 9) then
						    exit
					    endif
				    enddo

			    endif

		    enddo
		
		    sum_x = 0.
		    sum_y = 0.
		    do j = 1, 2 * ele_con + 1
                if (b_u(j) <= 0 .or. b_u(j) > nn) then
                    iter_check = .FALSE.
                    exit
                endif
			    sum_x(1) = sum_x(1) + node(1,b_u(j))
			    sum_y(1) = sum_y(1) + node(2,b_u(j))
		    enddo
		    do j = 1, 2 * (7 - ele_con) + 1
                if (b_l(j) <= 0 .or. b_l(j) > nn) then
                    iter_check = .FALSE.
                    exit
                endif
			    sum_x(2) = sum_x(2) + node(1,b_l(j))
			    sum_y(2) = sum_y(2) + node(2,b_l(j))
            enddo
            if (iter_check == .FALSE.) exit

		    node(1,i) = sum_x(1) / (2 * ele_con + 1)
		    node(2,i) = sum_y(1) / (2 * ele_con + 1)
		
		    nn = nn + 1
		    node(1,nn) = (sum_x(2) + node(1,i)) / (2 * (7 - ele_con) + 2)
		    node(2,nn) = (sum_y(2) + node(2,i)) / (2 * (7 - ele_con) + 2)
		
		    do j = 1, (7 - ele_con)
			    do k = 2, 5
				    if (div_ele(ele_l(j),k) == i) then
					    ele(k,div_ele(ele_l(j),1)) = nn
				    endif
			    enddo
		    enddo

		    ne = ne + 1
		    ele(1,ne) = ne
		    ele(2,ne) = b_u(1)
		    ele(3,ne) = i
		    ele(4,ne) = b_u(2*ele_con+1)
		    ele(5,ne) = nn
	        div_ele_num = 1
            exit
		    400 continue
		    !write (*,*) "-------------------------------------------------------------"
		    !if (ele_num == 0) then
		    !	write (*,*) "분할 할 수 없는 곳입니다."
		    !else
		    !	write (*,'(I5,A)') i, " node에서 7개의 Element로 연결되어 분할하였습니다." 
		    !	write (*,'(A,I5)') " 생성된 node : ", nn
		    !	write (*,'(A,I5)') " 생성된 Element : ", ne
		    !	write (*,*) ele(2:5,ne)  
		    !endif			
		    !write (*,*) "-------------------------------------------------------------"
        elseif (ele_num == 8) then
            !iter_check = .FALSE.
            !exit
            this(1) = i
            num = 0
            pre_num = div_ele(1,1)
            do j = 1, 4
                if (div_ele(1,j+1) == i) then
                    if (j == 1) then
                        this(2) = div_ele(1,5)
                    else
                        this(2) = div_ele(1,j)
                    endif
                    exit
                endif
            enddo
            
            iter = 0
            do 
                iter = iter + 1
                if (iter == 10) stop
                do j = 1, 8
                    if (pre_num /= div_ele(j,1)) then
                        elem = div_ele(j,2:5)
                        call find_count(2, this, 4, elem, count)
                        if (count == 2) then
                            pre_num = div_ele(j,1)
                            do k = 1, 4
                                if (this(2) == elem(k)) then
                                    if (k == 3) then
                                        ele_nodes(num+1:num+2) = (/ elem(4), elem(1) /)
                                    elseif (k == 4) then
                                        ele_nodes(num+1:num+2) = (/ elem(1), elem(2) /)
                                    else
                                        ele_nodes(num+1:num+2) = elem(k+1:k+2)
                                    endif
                                    this(2) = ele_nodes(num+2)
                                    num = num + 2
                                    exit
                                endif
                            enddo
                            exit
                        endif
                    endif
                enddo
                if (num == 16) exit
            enddo
            
            min_ang = 360.0
            do j = 1, 16
                if (j == 1) then
                    this = (/ ele_nodes(16), ele_nodes(2) /)
                elseif (j == 16) then
                    this = (/ ele_nodes(15), ele_nodes(1) /)
                else
                    this = (/ ele_nodes(j-1), ele_nodes(j+1) /)
                endif
                call calc_angle(node(:,this(1)), node(:,ele_nodes(j)), node(:,this(2)), ang)
                if (min_ang > ang) then
                    min_ang = ang
                    fn = j
                endif
            enddo
            temp_nodes = ele_nodes
            num = 0
            do j = fn, 16
                num = num + 1
                ele_nodes(num) = temp_nodes(j)
            enddo
            do j = 1, fn-1
                num = num + 1
                ele_nodes(num) = temp_nodes(j)
            enddo
            vec = node(:,ele_nodes(1+8)) - node(:,ele_nodes(1))
            node(:,nn+1) = node(:,ele_nodes(1)) + vec*0.25
            node(:,nn+1+4) = node(:,ele_nodes(1)) + vec*0.75
            vec = node(:,ele_nodes(3+8)) - node(:,ele_nodes(3))
            node(:,nn+2) = node(:,ele_nodes(3)) + vec*0.25
            node(:,nn+2+4) = node(:,ele_nodes(3)) + vec*0.75
            vec = node(:,ele_nodes(5+8)) - node(:,ele_nodes(5))
            node(:,nn+3) = node(:,ele_nodes(5)) + vec*0.25
            node(:,nn+3+4) = node(:,ele_nodes(5)) + vec*0.75
            vec = node(:,ele_nodes(7+8)) - node(:,ele_nodes(7))
            node(:,nn+4) = node(:,ele_nodes(7)) + vec*0.25
            node(:,nn+4+4) = node(:,ele_nodes(7)) + vec*0.75
            node(:,i) = 0.0
            do j = 1, 16
                node(:,i) = node(:,i) + node(:,ele_nodes(j))
            enddo
            node(:,i) = node(:,i)/16.0
            do j = 1, 4
                ce = div_ele(j*2-1,1)
                ele(1,ce) = ce
                ele(2:5,ce) = (/ ele_nodes((j-1)*4+2), ele_nodes((j-1)*4+3), nn+(j*2), nn+(j*2)-1 /)
            enddo
            do j = 1, 4
                ce = div_ele(j*2,1)
                ele(1,ce) = ce
                if (j == 4) then
                    ele(2:5,ce) = (/ ele_nodes((j-1)*4+3), ele_nodes((j-1)*4+4), nn+1, nn+(j*2) /)
                else
                    ele(2:5,ce) = (/ ele_nodes((j-1)*4+3), ele_nodes((j-1)*4+4), nn+(j*2)+1, nn+(j*2) /)
                endif
            enddo
            do j = 1, 8
                ele(1,ne+j) = ne+j
            enddo
            ele(2:5,ne+1) = (/ ele_nodes(16), ele_nodes(1), ele_nodes(2), nn+1 /)
            ele(2:5,ne+2) = (/ ele_nodes(4), ele_nodes(5), ele_nodes(6), nn+3 /)
            ele(2:5,ne+3) = (/ ele_nodes(8), ele_nodes(9), ele_nodes(10), nn+5 /)
            ele(2:5,ne+4) = (/ ele_nodes(12), ele_nodes(13), ele_nodes(14), nn+7 /)
            ele(2:5,ne+5) = (/ nn+8, nn+1, nn+2, i /)
            ele(2:5,ne+6) = (/ nn+2, nn+3, nn+4, i /)
            ele(2:5,ne+7) = (/ nn+4, nn+5, nn+6, i /)
            ele(2:5,ne+8) = (/ nn+6, nn+7, nn+8, i /)
            ne = ne + 8
            nn = nn + 8
            div_ele_num = 1
            exit
        elseif (ele_num >= 9) then
            iter_check = .FALSE.
            exit
	    endif
    enddo
    
    if (div_ele_num == 0) then
        exit
    elseif ( iter_check == .FALSE. ) then
        exit
    endif
enddo

end subroutine divide_over_ele
!======================================================================================

subroutine divide_bound_node(anes, pre_nn, nn, ne, node, ele, iter_check)

implicit none

integer, intent(in) :: anes, pre_nn
integer, intent(inout) :: ne, nn, ele(5,anes)
real, intent(inout) :: node(2,anes)
logical, intent(inout) :: iter_check

integer :: i, j, k, q, count, bound_type, iter
integer :: this(3), line(2), ele_nodes(6), elem(2)
real :: cri_ang(2), ang
logical :: check

cri_ang = (/ 170.0, 185.0 /)
iter = 0
do 
    iter = iter + 1
    bound_type = 0
    if (iter > pre_nn) then
        write (*,*) 'Error: subroutine divide_bound_node'
        iter_check = .FALSE.
        exit
    endif
    
    check = .TRUE.
    do i = 1, pre_nn
        this = (/ i-1, i, i+1 /)
        if (i == 1) this(1) = pre_nn
        if (i == pre_nn) this(3) = 1
        
        do j = 1, ne
            call find_count(3, this, 4, ele(2:5,j), count)
            if (count == 3) then
                elem(1) = j
                call calc_angle(node(:,this(1)), node(:,this(2)), node(:,this(3)), ang)
                if (ang >= cri_ang(1) .and. ang <= cri_ang(2)) then
                    line(1) = this(1)
                    do k = 2, 5
                        if (ele(k,j) /= this(1) .and. ele(k,j) /= this(2) .and. ele(k,j) /= this(3)) then
                            line(2) = ele(k,j)
                            exit
                        endif
                    enddo
                    
                    do k = 1, ne
                        if (k /= j) then
                            call find_count(2, line, 4, ele(2:5,k), count)
                            if (count == 2) then
                                elem(2) = k
                                if (this(1) == 1) then
                                    call calc_angle(node(:,pre_nn), node(:,this(1)), node(:,this(2)), ang)
                                else
                                    call calc_angle(node(:,this(1)-1), node(:,this(1)), node(:,this(2)), ang)
                                endif
                                if (ang >= cri_ang(1) .and. ang <= cri_ang(2)) then
                                    bound_type = 1
                                else
                                    bound_type = 2
                                endif
                                ele_nodes(1:3) = this
                                ele_nodes(4) = line(2)
                                do q = 2, 5
                                    if (ele(q,k) == line(2)) then
                                        if (q == 4) then
                                            ele_nodes(5:6) = (/ ele(q+1,k), ele(2,k) /)
                                        elseif (q == 5) then
                                            ele_nodes(5:6) = (/ ele(2,k), ele(3,k) /)
                                        else
                                            ele_nodes(5:6) = (/ ele(q+1,k), ele(q+2,k) /)
                                        endif
                                        exit
                                    endif
                                enddo
                                exit
                            endif
                        endif
                    enddo
                    if (bound_type == 1) then
                        node(:,nn+1) = (node(:,ele_nodes(1))+node(:,ele_nodes(5)))*0.5
                        node(:,nn+2) = (node(:,ele_nodes(2))+node(:,ele_nodes(4)))*0.5
                        ele(2:5,elem(1)) = (/ ele_nodes(1), nn+1, ele_nodes(5), ele_nodes(6) /)
                        ele(2:5,elem(2)) = (/ ele_nodes(1), ele_nodes(2), nn+2, nn+1 /)
                        ele(:,ne+1) = (/ ne+1, nn+1, nn+2, ele_nodes(4), ele_nodes(5) /)
                        ele(:,ne+2) = (/ ne+2, ele_nodes(2), ele_nodes(3), ele_nodes(4), nn+2 /)
                        nn = nn + 2
                        ne = ne + 2
                    else
                        node(:,nn+1) = (node(:,ele_nodes(2))+node(:,ele_nodes(4))+node(:,ele_nodes(6)))/3.0
                        ele(2:5,elem(1)) = (/ ele_nodes(1), ele_nodes(2), nn+1, ele_nodes(6) /)
                        ele(2:5,elem(2)) = (/ ele_nodes(2), ele_nodes(3), ele_nodes(4), nn+1 /)
                        ele(:,ne+1) = (/ ne+1, ele_nodes(4), ele_nodes(5), ele_nodes(6), nn+1 /)
                        nn = nn + 1
                        ne = ne + 1
                    endif
                    check = .FALSE.
                endif
            endif
        enddo
        if (check == .FALSE.) exit      
    enddo
    if (check) exit
enddo


end subroutine divide_bound_node