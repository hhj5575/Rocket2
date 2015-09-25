!*******    *********************************************************************
!
!  PROGRAM : Make elements
!
!  PURPOSE : 내부로 투영 된 절점과 경계절점으로 Element 생성
!
!****************************************************************************
subroutine make_ele(ci, temp_bs, nn, node, n_b, anes, ne, ele, bs, bs_num, b, pbn, me, en, ep, ic, b_num, end_b_num, end_b, iter_check)

implicit none

integer, intent(in) :: ci, bs, anes, en, n_b, temp_bs, bs_num(bs)
integer, intent(inout) :: nn, pbn, ne, ic, ele(5,anes), ep(bs), b_num, end_b_num(6)
real, intent(in) :: b(2,bs)
real, intent(inout) :: me(2,temp_bs), node(2,anes), end_b(2,6)
logical, intent(inout) :: iter_check

integer :: i, j, k, m, q, o, c, i_r, i_r_r, q_r, ele_r, temp_i
integer :: e_f, d_s, pbc, end_check, inner_num
real :: alpha, beta, dgr, criterion_d, cri_d, dis(3), tol_dis
logical :: cross(2)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

q = 0
k = 0
j = 0
m = 0
e_f = 0
d_s = 0
ic = 0
pbc = 0
criterion_d = 0.6
if (ci == 1) criterion_d = 0.30
cri_d = criterion_d * 1.4
do i = 1, bs
	!if ((pi == 2) .AND. (i == 3)) exit
	end_check = 0	
	if (i+j+k+m > bs) then
		exit
	elseif (i == 1) then !<-- 처음 절점일때는
		
		call calc_angle(b(:,bs), b(:,1), b(:,2), alpha)
		call calc_angle(b(:,1), b(:,2), b(:,3), beta)
		call calc_angle(b(:,2), b(:,3), b(:,4), dgr)
        call calc_len(b(:,bs), b(:,1), dis(1))
        call calc_len(b(:,1), b(:,2), dis(2))
        tol_dis = (dis(1)+dis(2))*0.5
        if (alpha < 180.0) then
            cri_d = (sin(alpha*pi/180.)*0.4+1.0)*criterion_d
        else
            cri_d = criterion_d * 1.4
        endif
        		
		if ((alpha >= 0.) .AND. (alpha <= 110.)) then !<- 끝점(end point) 일 때... 

			ne = ne + 1

			ele(:,ne) = (/ ne, bs_num(1), bs_num(2), n_b+1, bs_num(bs) /)
					
			pbc = pbc + 1;  k = 1;  q = 0
								
		elseif ((alpha > 110.) .AND. (alpha <= 240.)) then ! <- 측면절점(corner point)일 때...
		
			if ((beta >= 0.) .AND. (beta <= 110.)) then ! <- 처음:측면 + 다음:끝 일 때...

				ne = ne + 1
				
				ele(:,ne) = (/ ne, bs_num(1), bs_num(2), bs_num(3), n_b+1 /)
				
				pbc = pbc + 1
				q = 1; j = j + 1

			elseif ((beta > 110.) .AND. (beta <= 240.)) then ! <- 처음:측면 + 다음:측면 일 때...

				call calc_len(me(:,pbn), me(:,1), dis(1))
				call calc_len(me(:,1), me(:,2), dis(2))
				dis(3) = dis(1) + dis(2)
				
				do o = 1, en
					if (ep(o) == 1) then
						end_check = end_check + 1
					else
						end_check = end_check + 0
					endif
				enddo			

				if ( (dis(1) > cri_d*tol_dis) .OR. (dis(2) > cri_d*tol_dis) .OR. (end_check == 1) ) then !<- 거리가 길거나, 끝점 일 때...

					call calc_len(me(:,1), me(:,2), dis(1))
					call calc_len(me(:,2), me(:,3), dis(2))
					dis(3) = dis(1) + dis(2)
					
					end_check = 0 
					do o = 1, en
						if (ep(o) == 2) then
							end_check = end_check + 1
						else
							end_check = end_check + 0
						endif
					enddo

					if ( (dis(1) > cri_d*tol_dis) .OR. (dis(2) > cri_d*tol_dis) .OR. (end_check == 1) ) then !<- 거리가 길거나, 끝점 일 때...

						ne = ne + 1

						ele(:,ne) = (/ ne, bs_num(1), bs_num(2), n_b+2, n_b+1 /)
						
						pbc = pbc + 2; q = 0

					else                                               !<- 거리도 짧고, 끝점이 아닐 때...
						
						me(:,1) = (me(:,1) + me(:,2) + me(:,3)) / 3.

						do o = 2, pbn - 2
							me(:,o) = me(:,o+2)
						enddo
						pbn = pbn - 2	
						
						do o = 1, en
							if (ep(o) > 2) then
								ep(o) = ep(o) - 2
							endif
						enddo

						ne = ne + 1							
	
						ele(:,ne) = (/ ne, bs_num(1), bs_num(2), bs_num(3), n_b+1 /)
						
						pbc = pbc + 1;  q = 1;  j = j + 1
					
						if ((dgr > 240.) .AND. (dgr <= 360.)) then
							d_s = 1
						else 
							d_s = 0
						endif

					endif

				else !<- 거리도 짧고, 끝점이 아닐 때...
					
					me(:,1) = (me(:,pbn) + me(:,1) + me(:,2)) / 3.

					do o = 2, pbn - 2
						me(:,o) = me(:,o+1)
					enddo
					pbn = pbn - 2	
					do o = 1, en
						if (ep(o) > 1) then
							ep(o) = ep(o) - 1
						endif
					enddo
										
					ne = ne + 1

					ele(:,ne) = (/ ne, bs_num(1), bs_num(2), n_b+1, bs_num(bs) /)
					
					pbc = pbc + 1;  q = 0;  m = 1
					
					if ((dgr > 240.) .AND. (dgr <= 360.)) then
						d_s = 1
					else 
						d_s = 0
					endif
					
				endif

			elseif ((beta > 240.) .AND. (beta < 315.)) then ! <- 처음:측면 + 다음:코너 일 때...

				call calc_len(me(:,pbn), me(:,1), dis(1))
				call calc_len(me(:,1), me(:,2), dis(2))
				dis(3) = dis(1) + dis(2)
				
				do o = 1, en
					if (ep(o) == 1) then
						end_check = end_check + 1
					else
						end_check = end_check + 0
					endif
				enddo
			
				if ( ( (dis(1) > cri_d*tol_dis) .AND. (dis(2) > cri_d*tol_dis) ) .OR. (end_check == 1) ) then !<- 거리가 길거나, 끝점 일 때...
				
					ne = ne + 1

					ele(:,ne) = (/ ne, bs_num(1), bs_num(2), n_b+2, n_b+1 /)
					
					pbc = pbc + 2

				else !<- 거리도 짧고, 끝점이 아닐 때...
					
					me(:,1) = (me(:,pbn) + me(:,1) + me(:,2)) / 3.

					do o = 2, pbn - 2
						me(:,o) = me(:,o+1)
					enddo
					pbn = pbn - 2	
					do o = 1, en
						if (ep(o) > 1) then
							ep(o) = ep(o) - 2
						endif
					enddo
										
					ne = ne + 1							

					ele(:,ne) = (/ ne, bs_num(1), bs_num(2), n_b+1, bs_num(bs) /)
					
					pbc = pbc + 1;  k = 1
					
				endif

				ne = ne + 1

				ele(:,ne) = (/ ne, bs_num(2), n_b+4, n_b+3, ele(4+q,ne-1) /)
				
				pbc = pbc + 2;  q = -1

			elseif ((beta >= 315.) .AND. (beta < 360.)) then ! <- 처음:측면 + 다음:전환 일 때...

				call calc_len(me(:,pbn), me(:,1), dis(1))
				call calc_len(me(:,1), me(:,2), dis(2))
				dis(3) = dis(1) + dis(2)
				
				do o = 1, en
					if (ep(o) == 1) then
						end_check = end_check + 1
					else
						end_check = end_check + 0
					endif
				enddo
			
				if ( (dis(1) > cri_d*tol_dis) .OR. (dis(2) > cri_d*tol_dis) .OR. (end_check == 1) ) then !<- 거리가 길거나, 끝점 일 때...
				
					ne = ne + 1

					ele(:,ne) = (/ ne, bs_num(1), bs_num(2), n_b+2, n_b+1 /)
					
					pbc = pbc + 2

				else !<- 거리도 짧고, 끝점이 아닐 때...
					
					me(:,1) = (me(:,pbn) + me(:,1) + me(:,2)) / 3.

					do o = 2, pbn - 2
						me(:,o) = me(:,o+1)
					enddo
					pbn = pbn - 2	
					do o = 1, en
						if (ep(o) > 1) then
							ep(o) = ep(o) - 2
						endif
					enddo
										
					ne = ne + 1							

					ele(:,ne) = (/ ne, bs_num(1), bs_num(2), n_b+1, bs_num(bs) /)
					
					pbc = pbc + 1;  k = 1
					
				endif

				ne = ne + 1

				ele(:,ne) = (/ ne, bs_num(2), n_b+4, n_b+3, ele(4+q,ne-1) /)
				
				ne = ne + 1

				ele(:,ne) = (/ ne, bs_num(2), n_b+6, n_b+5, n_b+4 /)
				
				pbc = pbc + 4;  q = -1

			endif

		elseif ((alpha > 240.) .AND. (alpha < 360.)) then ! <- 코너절점(corner point)일 때...
			
			ne = ne + 1
			
			ele(:,ne) = (/ ne, bs_num(1), n_b+3, n_b+2, n_b+1 /)
			
			pbc = pbc + 3

			call calc_len(me(:,3), me(:,4), dis(1))
			call calc_len(me(:,4), me(:,5), dis(2))
			dis(3) = dis(1) + dis(2)
			
			do o = 1, en
				if (ep(o) == 4) then
					end_check = end_check + 1
				else
					end_check = end_check + 0
				endif
			enddo
		
			if ( (dis(1) > cri_d*tol_dis) .OR. (dis(2) > cri_d*tol_dis) .OR. (end_check == 1) ) then !<- 거리가 길거나, 끝점 일 때...
			
				ne = ne + 1

				ele(:,ne) = (/ ne, bs_num(1), bs_num(2), n_b+4, n_b+3 /)

				pbc = pbc + 1

			else !<- 거리도 짧고, 끝점이 아닐 때...
				
				me(:,3) = (me(:,3) + me(:,4) + me(:,5)) / 3.

				do o = 4, pbn - 2
					me(:,o) = me(:,o+1)
				enddo
				pbn = pbn - 2	
				do o = 1, en
					if (ep(o) > 4) then
						ep(o) = ep(o) - 2
					endif
				enddo
									
				ne = ne + 1							

				ele(:,ne) = (/ ne, bs_num(1), bs_num(2), bs_num(3), n_b+3 /)
				
				q = 1
				j = j + 1
				
			endif
						
		!elseif ((alpha >= 315.) .AND. (alpha < 360.)) then ! <- 전환절점(reversal point)일 때...
  !
		!	ne = ne + 1
  !
		!	ele(:,ne) = (/ ne, bs_num(1)), n_b+3, n_b+2, n_b+1 /)
		!	
		!	ne = ne + 1
  !
		!	ele(:,ne) = (/ ne, bs_num(1)), n_b+5, n_b+4, n_b+3 /)
		!	
		!	pbc = pbc + 4
		!	
		!	call calc_len(me(:,5), me(:,6), dis(1))
		!	call calc_len(me(:,6), me(:,7), dis(2))
		!	dis(3) = dis(1) + dis(2)
		!	
		!	do o = 1, en
		!		if (ep(o) == 6) then
		!			end_check = end_check + 1
		!		else
		!			end_check = end_check + 0
		!		endif
		!	enddo
		!
		!	if ((dis(3) > 1.2*tol_dis) .OR. (end_check == 1)) then !<- 거리가 길거나, 끝점 일 때...
		!	
		!		ne = ne + 1
  !
		!		ele(:,ne) = (/ ne, bs_num(1)), bs_num(2)), n_b+6, n_b+5 /)
		!		
		!		pbc = pbc + 1
  !
		!	else !<- 거리도 짧고, 끝점이 아닐 때...
		!		
		!		me(:,5) = (me(:,5) + me(:,6) + me(:,7)) / 3.
  !
		!		do o = 4, pbn - 2
		!			me(:,o) = me(:,o+1)
		!		enddo
		!		pbn = pbn - 2	
		!		do o = 1, en
		!			if (ep(o) > 6) then
		!				ep(o) = ep(o) - 2
		!			endif
		!		enddo
		!								
		!		ne = ne + 1							
  !
		!		ele(:,ne) = (/ ne, bs_num(1)), bs_num(2)), bs_num(3)), n_b+5 /)
		!		
		!		q = 1
		!		j = 1
		!		
		!	endif

		endif
			
	elseif ((i+j+m == bs) .AND. (k .NE. 1)) then !<- 마지막 절점 일 때...

		if ( m == 1 ) then

			call calc_cross(b(:,bs-1), b(:,bs), me(:,1), me(:,pbn), cross(1))
			call calc_cross(b(:,bs), me(:,1), me(:,pbn), b(:,bs-1), cross(2))

			if ( ( cross(1) .EQV. .FALSE. ) .OR. ( cross(2) .EQV. .FALSE. ) ) then
			
				me(:,1) = ( 2 * me(:,1) + b(:,bs) ) / 3

			endif
			ne = ne + 1
			ele(:,ne) = (/ ne, bs_num(bs-1), bs_num(bs), n_b+1, ele(4+q,ne-1) /)

			j = j + 2;  ic = 1

		else

			call calc_cross(b(:,i+j), b(:,1), me(:,1), me(:,pbn), cross(1))
			call calc_cross(b(:,1), me(:,1), me(:,pbn), b(:,i+j), cross(2))

			if ( ( cross(1) .EQV. .FALSE. ) .OR. ( cross(2) .EQV. .FALSE. ) ) then
			
				me(:,1) = ( 2 * me(:,1) + b(:,1) ) / 3

			endif
			ne = ne + 1
			ele(:,ne) = (/ ne, bs_num(i+j), bs_num(1), n_b+1, ele(4+q,ne-1) /)

			j = j + 2;  ic = 1

		endif
		
	else                                           !<- 처음과 마지막 절점 아닐 때...
		if (pbn == 1) then

			ne = ne + 1

			ele(:,ne) = (/ ne, bs_num(i+j), bs_num(i+j+1), bs_num(i+j+2), ele(4+q,ne-1) /)
			
			j = j + 1;  q = 1;  d_s = 0;  e_f = 0

		else
		
			if (pbc == pbn) then
				i_r = 1
				ele_r = n_b + 1
				if (pbn <= 3) then
					ic = 0	
				else
					ic = 1
				endif
			else
				i_r = pbc + 1
				ele_r = ele(4+q,ne) + 1
            endif
			
            call calc_len(b(:,i+j), b(:,i+j+1), dis(1))
			if (i+j+1 == bs-1) then
				call calc_angle(b(:,i+j), b(:,i+j+1), b(:,i+j+2), alpha)
				call calc_angle(b(:,i+j+1), b(:,i+j+2), b(:,1), beta)
                call calc_len(b(:,i+j+1), b(:,i+j+2), dis(2))
			elseif (i+j+1 == bs) then
				call calc_angle(b(:,i+j), b(:,i+j+1), b(:,1), alpha)
				call calc_angle(b(:,i+j+1), b(:,1), b(:,2), beta)
                call calc_len(b(:,i+j+1), b(:,1), dis(2))
			else
				call calc_angle(b(:,i+j), b(:,i+j+1), b(:,i+j+2), alpha)
				call calc_angle(b(:,i+j+1), b(:,i+j+2), b(:,i+j+3), beta)
                call calc_len(b(:,i+j+1), b(:,i+j+2), dis(2))
            endif
            tol_dis = (dis(1)+dis(2))*0.5
            if (alpha < 180.0) then
                cri_d = (sin(alpha*pi/180.)*0.4+1.0)*criterion_d
            else
                cri_d = 1.4*criterion_d
            endif
			
			if ((alpha >= 0.) .AND. (alpha <= 110.)) then ! <- 끝절점(point point)일 때...
				
				ne = ne + 1 
				ele(1:3,ne) = (/ ne, bs_num(i+j), bs_num(i+j+1) /)
				if ( (i+j+1 == bs) .AND. (k == 1) )then
					ele(4,ne) = n_b + 1
				elseif (i+j+1 == bs) then
					ele(4,ne) = bs_num(1)
				else	
					ele(4,ne) = bs_num(i+j+2)
				endif
				ele(5,ne) = ele(4+q,ne-1)
		
				j = j + 1;  q = 1;  d_s = 0

				if ((beta >= 240.) .AND. (beta <= 360.)) then
					e_f = 1
				else
					e_f = 0
				endif

			elseif ((alpha > 110.) .AND. (alpha < 240.)) then ! <- 측면절점(corner point)일 때...
				
				if (pbc == pbn - 1) then
					i_r_r = 1
				else 
					i_r_r = i_r + 1
				endif

				call calc_len(me(:,pbc), me(:,i_r), dis(1))
				call calc_len(me(:,i_r), me(:,i_r_r), dis(2))
				dis(3) = dis(1) + dis(2)
				do o = 1, en
					if (ep(o) == i_r) then
						end_check = end_check + 1
					else
						end_check = end_check + 0
					endif
				enddo

				if (pbn == 2) then
					
					me(:,1) = (me(:,1) + me(:,2)) / 2.

					pbn = 1

					ne = ne + 1

					ele(:,ne) = (/ ne, bs_num(i+j), bs_num(i+j+1), bs_num(i+j+2), ele(4+q,ne-1) /)
										
					j = j + 1;  q = 1;  d_s = 0;  e_f = 0	
					
				elseif ( (dis(1) > cri_d*tol_dis) .OR. (dis(2) > cri_d*tol_dis) .OR. (end_check == 1) ) then  !<- 거리가 길거나, 끝점 일 때...
				
					ne = ne + 1

					ele(:,ne) = (/ ne, bs_num(i+j), bs_num(i+j+1), ele_r, ele(4+q,ne-1) /)
										
					if (pbc .NE. pbn) then
						pbc = pbc + 1
					endif
                    
					if ( ( ele_r == n_b + 1) .AND. ( pbc == pbn ) ) then
						call calc_cross(b(:,i+j), b(:,i+j+1), me(:,1), me(:,pbn), cross(1))
						call calc_cross(b(:,i+j+1), me(:,1), me(:,pbn), b(:,i+j), cross(2))

						if ( ( cross(1) .EQV. .FALSE. ) .OR. ( cross(2) .EQV. .FALSE. ) ) then
		
							me(:,1) = ( 2 * me(:,1) + b(:,1) ) / 3

						endif
					endif
						
					q = 0;  d_s = 0;  e_f = 0
										
				else !<- 거리도 짧고, 끝점이 아닐 때...
					
					if (pbc == pbn - 1) then
						
						temp_i = ele(4+q,ne)
						
						me(:,1) = (me(:,i_r-1) + me(:,i_r) + me(:,i_r_r)) / 3.
						
						pbn = pbn - 2
					
						do o = 1, ne
							do c = 2, 5
								if ( ele(c,o) == temp_i ) then
									ele(c,o) = n_b + 1
								endif
							enddo
                        enddo
						
                        if (i+j+1 > bs) then
                            iter_check = .false.
                            exit
                        endif
                        
						ne = ne + 1							
						ele(1:3,ne) = (/ ne, bs_num(i+j), bs_num(i+j+1) /)
						if ( i+j+1 == bs ) then
							ele(4,ne) = bs_num(1)
						else
							ele(4,ne) = bs_num(i+j+2)
						endif
						ele(5,ne) = n_b + 1

						j = j + 1;  q = 1;  e_f = 0;  ic = 1

					else

						me(:,i_r-1) = (me(:,i_r-1) + me(:,i_r) + me(:,i_r_r)) / 3.

						do o = i_r, pbn - 2
							me(:,o) = me(:,o+2)
						enddo
						pbn = pbn - 2	
						do o = 1, en
							if (ep(o) > i_r) then
								ep(o) = ep(o) - 2
							endif
                        enddo
                        
                        if (i+j+1 > bs .or. i+j+2 > bs) then
                            iter_check = .false.
                            exit
                        endif
                        
						ne = ne + 1							

						ele(:,ne) = (/ ne, bs_num(i+j), bs_num(i+j+1), bs_num(i+j+2), ele(4+q,ne-1) /)
						
						q = 1;  e_f = 0
						if ((beta > 240.) .AND. (beta <= 360.)) then
							d_s = 1
						else
							d_s = 0
							j = j + 1
						endif

					endif
					
				endif
				
			elseif ((alpha >= 240.) .AND. (alpha < 360.)) then ! <- 코너절점(corner point)일 때...
				
				if ((d_s .NE. 1) .AND. (e_f .NE. 1)) then
				
					ne = ne + 1

					ele(:,ne) = (/ ne, bs_num(i+j), bs_num(i+j+1), ele_r, ele(4+q,ne-1) /)
					
					pbc = pbc + 1;  q = 0

				endif

				ne = ne + 1;  q_r = ele(4+q,ne-1)

				ele(:,ne) = (/ ne, bs_num(i+j+1), q_r+2, q_r+1, q_r /)
				
				pbc = pbc + 2;  q = -1;  d_s = 0;  e_f = 0

			!elseif ((alpha >= 315.) .AND. (alpha < 360.)) then ! <- 전환절점(reversal point)일 때...
   !
			!	if ((d_s .NE. 1) .AND. (e_f .NE. 1)) then
			!	
			!		ne = ne + 1
   !
			!		ele(:,ne) = (/ ne, bs_num(i+j)), bs_num(i+j+1)), ele_r, ele(4+q,ne-1) /)
			!		
			!		pbc = pbc + 1;  q = 0
   !
			!	endif
   !
			!	ne = ne + 1;  q_r = ele(4+q,ne-1)
   !
			!	ele(:,ne) = (/ ne, bs_num(i+j+1)), q_r+2, q_r+1, q_r /)
			!	
			!	pbc = pbc + 2;  q = -1
   !
			!	ne = ne + 1;  q_r = ele(4+q,ne-1)
   !
			!	ele(:,ne) = (/ ne, bs_num(i+j+1)), q_r+2, q_r+1, q_r /)
			!	
			!	pbc = pbc + 2;  q = -1;  d_s = 0;  e_f = 0

			endif
		endif
	endif

	if (ic .NE. 1) then
		inner_num = bs + 1 - k - i - j - m + pbc
		if (inner_num <= 6) then
			!write (*,*) "---=== 남은 경계절점과 유동절점이 6이하가 되었습니다.<Element> ===---"
			b_num = 0
			!write (*,*) i+j, bs , pbc
			do o = i + j + 1, bs
				b_num = b_num + 1
                end_b_num(b_num) = bs_num(o)
				end_b(:,b_num) = b(:,o)
			enddo
			if ( (k .NE. 1) .AND. (m .NE. 1) )then
				b_num = b_num + 1
				end_b_num(b_num) = bs_num(1)
				end_b(:,b_num) = b(:,1)
			endif
			do o = 1, pbc 
				b_num = b_num + 1
				end_b_num(b_num) = nn + o
				end_b(:,b_num) = me(:,o)
				node(:,nn+o) = me(:,o)
			enddo
			nn = nn + pbc
			
			!write (*,'(A,I3)') "bs : ", b_num
			exit
		endif
	endif
	b_num = pbn
enddo

end subroutine make_ele
