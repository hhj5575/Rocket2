!****************************************************************************
!
!  PROGRAM : parition
!
!  PURPOSE : 180도 이상되는 절점을 찾아 전체 구간을 분할
!
!****************************************************************************
subroutine partition(nn, anes, node, bn, db, dbs, dn, first_nn)

implicit none

integer, intent(inout) :: nn, dn, first_nn
integer, intent(in) :: bn, anes
integer, intent(inout) :: db(first_nn,bn), dbs(first_nn)
real, intent(inout) :: node(2,anes)

integer :: i, j, k, m, bs, i_l, i_r, j_l, j_r, k_r, ins_dn, cp, rmn, fp, lp, knn, en
integer :: temp_dbs, temp_dn, temp_kp, temp_case
integer :: kp(nn), temp_db(bn), ep(nn), part_case(2)
real :: length, seg, angle, len_a, len_f, len_sum, div_d, temp_d, dis, cri, cri_angle
real :: scale_factor, a0, am, r, sn
real :: vector(2), alpha(2), beta(2), temp_ang(2), val(2), temp_val(2)
logical :: cross, check_ep

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

dn = 1;   kp = 0;   bs = 0;   ins_dn = 0
knn = 0;   en = 0;   cri = 58.
cri_angle = 220.

do i = 1, dbs(1)
		
	if ( i == 1 ) then
		i_l = dbs(1)
		i_r = i + 1
	elseif ( i == dbs(1) ) then
		i_l = i - 1
		i_r = 1
	else
		i_l = i - 1
		i_r = i + 1
	endif
		
	call calc_angle(node(:,db(1,i_l)), node(:,db(1,i)), node(:,db(1,i_r)), angle )

	if ( angle > cri_angle ) then
		
		if ( knn == 0 ) then
			knn = knn + 1
			kp(knn) = db(1,i)
		elseif ( kp(knn)+1 /= i ) then
			knn = knn + 1
			kp(knn) = db(1,i)
		endif
	
	elseif ( angle <= 110. ) then
		
		if ( i == 1 ) then
			ep(en+1:en+2) = (/ db(1,dbs(1)), db(1,i_r) /)
			en = en + 2
		elseif (i == dbs(1) ) then
			ep(en+1:en+2) = (/ db(1,i_l), db(1,1) /)
			en = en + 2
		else
			ep(en+1:en+2) = (/ db(1,i_l), db(1,i_r) /)
			en = en + 2
		endif
	endif

enddo

!do i = 1, knn
	!write (*,'(I4,$)') kp(i)
!enddo

!write (*,'(/,A)') "--------------"

do i = 1, knn
	cp = 0
	len_f = 0
    part_case = 0
    temp_kp = 0
	do k = 1, dn
		do j = 1, dbs(k)
			if ( db(k,j) == kp(i) ) then
				temp_kp = j
				temp_dn = k
                exit
			endif
        enddo
        if (temp_kp /= 0) exit
	enddo
	
	if ( ( dn == 1 ) .OR. ( db(temp_dn,2) /= kp(i) ) )then	
	
		if ( temp_kp == 1) then
			i_l = db(temp_dn, dbs(temp_dn))
			i_r = db(temp_dn,temp_kp + 1)
		elseif( temp_kp == dbs(temp_dn) ) then
			i_l = db(temp_dn,temp_kp - 1)
			i_r = db(temp_dn,1)
		else
			i_l = db(temp_dn,temp_kp - 1)
			i_r = db(temp_dn,temp_kp + 1)
		endif
		
		do j = 1, dbs(temp_dn)
			if ( j == 1 ) then
				j_l = db(temp_dn, dbs(temp_dn))
				j_r = db(temp_dn,j + 1)
			elseif( j == dbs(temp_dn) ) then
				j_l = db(temp_dn,j - 1)
				j_r = db(temp_dn,1)
			else
				j_l = db(temp_dn,j - 1)
				j_r = db(temp_dn,j + 1)
			endif

			check_ep = .TRUE.
			do k = 1, en
				if ( ep(k) == db(temp_dn,j) ) then
					check_ep = .FALSE.
					exit
				endif
			enddo
            if (temp_kp == j) check_ep = .FALSE.
            
			if ( check_ep ) then

				call calc_angle( node(:,j_l), node(:,db(temp_dn,j)), node(:,db(temp_dn,temp_kp)), alpha(1) )
				call calc_angle( node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), node(:,j_r), beta(1) )
				call calc_angle( node(:,i_l), node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), alpha(2) )
				call calc_angle( node(:,db(temp_dn,j)), node(:,db(temp_dn,temp_kp)), node(:,i_r), beta(2) )
				call calc_len( node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), dis )
		
				if ( ((alpha(1) > cri) .AND. (alpha(1) < 180.0+cri)) .AND. ((beta(1) > cri) .AND. (beta(1) < 180.0+cri)) ) then
				if ( ((alpha(2) > cri) .AND. (alpha(2) < 180.0+cri)) .AND. ((beta(2) > cri) .AND. (beta(2) < 180.0+cri)) ) then
				!if ( dis >= d ) then 
					
					rmn = 0
					len_sum = 0
					if ( temp_kp > j ) then
						do k = 1, j 
							rmn = rmn + 1
							if ( k /= j ) then
								call calc_len(node(:,db(temp_dn,k)), node(:,db(temp_dn,k+1)), length)
								len_sum = len_sum + length
							endif
						enddo
						do k = temp_kp, dbs(temp_dn)
							rmn = rmn + 1
							if ( k == dbs(temp_dn) ) then
								call calc_len(node(:,db(temp_dn,k)), node(:,db(temp_dn,1)), length)
							else
								call calc_len(node(:,db(temp_dn,k)), node(:,db(temp_dn,k+1)), length)
							endif
							len_sum = len_sum + length
						enddo
					else
						do k = temp_kp, j 
							rmn = rmn + 1
							if ( k /= j ) then
								call calc_len(node(:,db(temp_dn,k)), node(:,db(temp_dn,k+1)), length)
								len_sum = len_sum + length
							endif
						enddo
					endif
					temp_d = len_sum / (rmn - 1)

					do k = 1, dbs(temp_dn)
						if( k == dbs(temp_dn) ) then
							k_r = db(temp_dn,1)
						else
							k_r = db(temp_dn,k + 1)
						endif
						
						call calc_cross( node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), node(:,db(temp_dn,k)), node(:,k_r), cross )
						
						if ( cross .EQV. .FALSE. ) then
							exit
						endif				
                    enddo
                    
                    temp_case = 0
					if ( ( rmn /= (rmn / 2) * 2 ) .AND. (rmn >= 3) .AND. ( cross ) ) then       ! odd number
						temp_case = 1
                    elseif ( ( rmn == (rmn / 2) * 2 ) .AND. (rmn >= 4) .AND. ( cross ) ) then   ! even number
                        temp_case = 2
                    endif
                    if (temp_case /= 0) then
						!len_a = sqrt( ((node(db(temp_dn,temp_kp),1) - node(db(temp_dn,j),1)) ** 2) + ((node(db(temp_dn,temp_kp),2) - node(db(temp_dn,j),2)) ** 2) )
						call calc_len(node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), len_a)
						call calc_len(node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,i_l)), val(1))
                        call calc_len(node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,i_r)), val(2))
                        temp_val(1) = min(val(1), val(2))
                        call calc_len(node(:,db(temp_dn,j)), node(:,db(temp_dn,j_l)), val(1))
                        call calc_len(node(:,db(temp_dn,j)), node(:,db(temp_dn,j_r)), val(2))
                        temp_val(2) = min(val(1), val(2))
						if ( len_f == 0. ) then
							len_f = len_a 
							cp = j
							temp_ang(1) = alpha(2)
							temp_ang(2) = beta(2)
							div_d = temp_d
                            a0 = temp_val(1)
                            am = temp_val(2)
                            part_case(1) = temp_case
						else
							if ( len_a < len_f ) then
								len_f = len_a
								cp = j
								temp_ang(1) = alpha(2)
								temp_ang(2) = beta(2)
								div_d = temp_d
                                a0 = temp_val(1)
                                am = temp_val(2)
                                part_case(1) = temp_case
							endif
						endif
					endif
				!endif	
				endif
				endif
			endif
		enddo
						
		if ( cp .NE. 0 ) then
            write (*,*)
            write (*,*) "========================================================"
		    write (*,'(I3,A)') i, "th partition"
		    write (*,'(A,I7,A,I7)') 'Partition  Start Info - Seq: ', temp_kp, ', Node num:', db(temp_dn,temp_kp)
            write (*,'(A,I7,A,I7)') 'Partition    End Info - Seq: ', cp, ', Node num:', db(temp_dn,cp)
            write (*,'(A,I7,A,I7)') 'Partition Domain Info - Num: ', temp_dn, ',    Count:', dbs(temp_dn)
            
            call calc_divided_number(len_f, a0, am, m, scale_factor)
            if (part_case(1) == 1) then
                write (*,'(A,I3,A)') 'Partition case: ', part_case(1), ' - Number of remainder node is odd number'
            else
                write (*,'(A,I3,A)') 'Partition case: ', part_case(1), ' - Number of remainder node is even number'
            endif
            ! part_case = 1: number of remainder node is odd number
            !           = 2: number of remainder node is even number
            !write (*,'(A,I5)') 'Previous m = ', m
            if (part_case(1) == 1) then  ! odd number
                !m = 2 * nint( len_f / (2 * div_d) )
                if (m /= 2*(m/2)) then ! if m(edge number) is odd number
                    if (scale_factor >= 1.0 .or. m == 1) then
                        m = m + 1
                    else
                        m = m - 1
                    endif
                    call adjust_scale_factor(len_f, a0, am, m, scale_factor)
                endif
            else                         ! even number
                !m = 2 * nint( len_f / (2 * div_d) ) + 1
                if (m == 2*(m/2)) then ! if m is even number
                    if (scale_factor >= 1.0) then
                        m = m + 1
                    elseif (m /= 0) then
                        m = m - 1
                    endif
                    call adjust_scale_factor(len_f, a0, am, m, scale_factor)
                endif
            endif
            !write (*,'(A,I5)') '  Adjust m = ', m
			
            vector(:) = node(:,db(temp_dn,cp)) - node(:,db(temp_dn,temp_kp))
			if (m == 0 .or. m == 1) then
        !        if (m == 0) then
				    !nn = nn + 1
				    !!node(:,nn) = ( node(:,db(temp_dn,cp)) + node(:,db(temp_dn,temp_kp)) ) * 0.5
        !            node(:,nn) = node(:,db(temp_dn,temp_kp)) + vector*a0/len_f
				    !fp = nn;  lp = nn
        !        else
        !            fp = 0;  lp = 0
        !        endif
                fp = 0;  lp = 0
			else
				!seg = len_f / m
				!vector(:) = vector(:) / len_f
				sn = 0.0
                r = ((am-a0)/len_f) + 1.0
                do j = 0, m - 2
                    sn = sn + a0*(r**j)
                    nn = nn + 1
                    node(:,nn) = node(:,db(temp_dn,temp_kp)) + (sn * vector(:))*scale_factor
                    if ( j == 0 ) then
						fp = nn
						if ( j == m-2) then
							lp = nn
						endif
					elseif( j == m-2 ) then
						lp = nn
					endif
                enddo
    		endif
			
			temp_dbs = dbs(temp_dn)
			do j = 1, temp_dbs
				temp_db(j) = db(temp_dn, j)
			enddo
			dbs(temp_dn) = 0
			if ( temp_kp > cp ) then
				
				do j = 1, cp
					dbs(temp_dn) = dbs(temp_dn) + 1
					db(temp_dn,dbs(temp_dn)) = temp_db(j)
                enddo
                if (fp /= 0) then
				    do j = lp, fp, -1
					    dbs(temp_dn) = dbs(temp_dn) + 1
					    db(temp_dn,dbs(temp_dn)) = j 
                    enddo
                endif
				do j = temp_kp, temp_dbs
					dbs(temp_dn) = dbs(temp_dn) + 1
					db(temp_dn,dbs(temp_dn)) = temp_db(j)
				enddo

				dn = dn + 1
				do j = cp, temp_kp
					dbs(dn) = dbs(dn) + 1
					db(dn,dbs(dn)) = temp_db(j)
                enddo
                if (fp /= 0) then
				    do j = fp, lp
					    dbs(dn) = dbs(dn) + 1
					    db(dn,dbs(dn)) = j
                    enddo
                endif
				
				if ( temp_ang(1) > 190. ) then
					temp_dn = dn
				elseif ( temp_ang(2) > 190. ) then
					temp_dn = temp_dn
				else 
					temp_dn = 0
				endif
															
			else 
			
				do j = 1, temp_kp
					dbs(temp_dn) = dbs(temp_dn) + 1
					db(temp_dn,dbs(temp_dn)) = temp_db(j)
                enddo
                if (fp /= 0) then
				    do j = fp, lp
					    dbs(temp_dn) = dbs(temp_dn) + 1
					    db(temp_dn,dbs(temp_dn)) = j
                    enddo
                endif
				do j = cp, temp_dbs
					dbs(temp_dn) = dbs(temp_dn) + 1
					db(temp_dn,dbs(temp_dn)) = temp_db(j)
				enddo

				dn = dn + 1
				do j = temp_kp, cp
					dbs(dn) = dbs(dn) + 1
					db(dn,dbs(dn)) = temp_db(j)
                enddo
                if (fp /= 0) then
				    do j = lp, fp, -1
					    dbs(dn) = dbs(dn) + 1
					    db(dn,dbs(dn)) = j
                    enddo
                endif

				if ( temp_ang(1) > 190.  ) then
					temp_dn = temp_dn
				elseif ( temp_ang(2) > 190.  ) then
					temp_dn = dn
				else 
					temp_dn = 0
				endif

			endif
			
		else

			temp_dn = 0

		endif
		
		if ( temp_dn .NE. 0 ) then   !<- if the angle is more than 180 at keypoint
			
			cp = 0
			len_f = 0
			do j = 1, dbs(temp_dn)
				if ( db(temp_dn,j) == kp(i) ) then
					temp_kp = j
				endif
			enddo

			if ( temp_kp == 1) then
				i_l = db(temp_dn, dbs(temp_dn))
				i_r = db(temp_dn,temp_kp + 1)
			elseif( temp_kp == dbs(temp_dn) ) then
				i_l = db(temp_dn,temp_kp - 1)
				i_r = db(temp_dn,1)
			else
				i_l = db(temp_dn,temp_kp - 1)
				i_r = db(temp_dn,temp_kp + 1)
			endif

			do j = 1, dbs(temp_dn)
				if ( j == 1 ) then
					j_l = db(temp_dn, dbs(temp_dn))
					j_r = db(temp_dn,j + 1)
				elseif( j == dbs(temp_dn) ) then
					j_l = db(temp_dn,j - 1)
					j_r = db(temp_dn,1)
				else
					j_l = db(temp_dn,j - 1)
					j_r = db(temp_dn,j + 1)
				endif

				check_ep = .TRUE.
				do k = 1, en
					if ( ep(k) == db(temp_dn,j) ) then
						check_ep = .FALSE.
						exit
					endif
				enddo

				if ( check_ep ) then
				
					call calc_angle( node(:,j_l), node(:,db(temp_dn,j)), node(:,db(temp_dn,temp_kp)), alpha(1) )
					call calc_angle( node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), node(:,j_r), beta(1) )
					call calc_angle( node(:,i_l), node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), alpha(2) )
					call calc_angle( node(:,db(temp_dn,j)), node(:,db(temp_dn,temp_kp)), node(:,i_r), beta(2) )
					call calc_len( node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), dis )

					if ( ((alpha(1) > cri) .AND. (alpha(1) < 180.0+cri)) .AND. ((beta(1) > cri) .AND. (beta(1) < 180.0+cri)) ) then
					if ( ((alpha(2) > cri) .AND. (alpha(2) < 180.0+cri)) .AND. ((beta(2) > cri) .AND. (beta(2) < 180.0+cri)) ) then
					!if ( dis >= d ) then
						
						rmn = 0
                        len_sum = 0
						if ( temp_kp > j ) then
							do k = 1, j 
								rmn = rmn + 1
								if ( k /= j ) then
									call calc_len(node(:,db(temp_dn,k)), node(:,db(temp_dn,k+1)), length)
									len_sum = len_sum + length	
								endif
							enddo
							do k = temp_kp, dbs(temp_dn)
								rmn = rmn + 1
								if ( k == dbs(temp_dn) ) then
									call calc_len(node(:,db(temp_dn,k)), node(:,db(temp_dn,1)), length)
								else
									call calc_len(node(:,db(temp_dn,k)), node(:,db(temp_dn,k+1)), length)
								endif
								len_sum = len_sum + length
							enddo
						else
							do k = temp_kp, j 
								rmn = rmn + 1
								if ( k /= j ) then
									call calc_len(node(:,db(temp_dn,k)), node(:,db(temp_dn,k+1)), length)
									len_sum = len_sum + length
								endif
							enddo
						endif
						temp_d = len_sum / (rmn - 1)

						do k = 1, dbs(temp_dn)
							if( k == dbs(temp_dn) ) then
								k_r = db(temp_dn,1)
							else
								k_r = db(temp_dn,k + 1)
							endif
							
							call calc_cross( node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,j)), node(:,db(temp_dn,k)), node(:,k_r), cross )
							
							if ( cross .EQV. .FALSE. ) then
								exit
							endif				
                        enddo
                        
                        temp_case = 0
						if ( ( rmn .NE. (rmn / 2) * 2 ) .AND. (rmn >= 3) .AND. ( cross ) ) then
							temp_case = 1
                        elseif ( ( rmn == (rmn / 2) * 2 ) .AND. (rmn >= 4) .AND. ( cross ) ) then
                            temp_case = 2
                        endif
                        if (temp_case /= 0) then
							!len_a = sqrt( ((node(db(temp_dn,j),1) - node(db(temp_dn,temp_kp),2)) ** 1) + ((node(db(temp_dn,j),2) - node(db(temp_dn,temp_kp),2)) ** 2) )
							call calc_len(node(:,db(temp_dn,j)), node(:,db(temp_dn,temp_kp)), len_a)
                            call calc_len(node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,i_l)), val(1))
                            call calc_len(node(:,db(temp_dn,temp_kp)), node(:,db(temp_dn,i_r)), val(2))
                            temp_val(1) = min(val(1), val(2))
                            call calc_len(node(:,db(temp_dn,j)), node(:,db(temp_dn,j_l)), val(1))
                            call calc_len(node(:,db(temp_dn,j)), node(:,db(temp_dn,j_r)), val(2))
                            temp_val(2) = min(val(1), val(2))
							if ( len_f == 0 ) then
								len_f = len_a 
								cp = j
                                div_d = temp_d
                                a0 = temp_val(1)
                                am = temp_val(2)
                                part_case(2) = temp_case
							else
								if ( len_a < len_f ) then
									len_f = len_a 
									cp = j
                                    div_d = temp_d
                                    a0 = temp_val(1)
                                    am = temp_val(2)
                                    part_case(2) = temp_case
								endif
							endif
						endif
					!endif
					endif
					endif

				endif
				
			enddo
			
			if ( cp .NE. 0 ) then
                
                write (*,*)
                write (*,*) "========================================================"
		        write (*,'(I3,A)') i, "th, created more than 180"
		        write (*,'(A,I7,A,I7)') 'Partition  Start Info - Seq: ', temp_kp, ', Node num:', db(temp_dn,temp_kp)
                write (*,'(A,I7,A,I7)') 'Partition    End Info - Seq: ', cp, ', Node num:', db(temp_dn,cp)
                write (*,'(A,I7,A,I7)') 'Partition Domain Info - Num: ', temp_dn, ',    Count:', dbs(temp_dn)
            
                call calc_divided_number(len_f, a0, am, m, scale_factor)
                if (part_case(2) == 1) then
                    write (*,'(A,I3,A)') 'Partition case: ', part_case(2), ' - Number of remainder node is odd number'
                else
                    write (*,'(A,I3,A)') 'Partition case: ', part_case(2), ' - Number of remainder node is even number'
                endif

                if (part_case(2) == 1) then  ! odd number
                    !m = 2 * nint( len_f / (2 * div_d) )
                    if (m /= 2*(m/2)) then ! if m(edge number) is odd number
                        if (scale_factor >= 1.0 .or. m == 1) then
                            m = m + 1
                        else
                            m = m - 1
                        endif
                        call adjust_scale_factor(len_f, a0, am, m, scale_factor)
                    endif
                else                         ! even number
                    !m = 2 * nint( len_f / (2 * div_d) ) + 1
                    if (m == 2*(m/2)) then ! if m is even number
                        if (scale_factor >= 1.0) then
                            m = m + 1
                        elseif (m /= 0) then
                            m = m - 1
                        endif
                        call adjust_scale_factor(len_f, a0, am, m, scale_factor)
                    endif
                endif
                !write (*,'(A,I5)') '  Adjust m = ', m
			
                vector(:) = node(:,db(temp_dn,cp)) - node(:,db(temp_dn,temp_kp))
			    if (m == 0 .or. m == 1) then
            !        if (m == 0) then
				        !nn = nn + 1
				        !!node(:,nn) = ( node(:,db(temp_dn,cp)) + node(:,db(temp_dn,temp_kp)) ) * 0.5
            !            node(:,nn) = node(:,db(temp_dn,temp_kp)) + vector*a0/len_f
				        !fp = nn;  lp = nn
            !        else
            !            fp = 0;  lp = 0
            !        endif
                    fp = 0;  lp = 0
			    else
				    !seg = len_f / m
				    !vector(:) = vector(:) / len_f
				    sn = 0.0
                    r = ((am-a0)/len_f) + 1.0
                    do j = 0, m - 2
                        sn = sn + a0*(r**j)
                        nn = nn + 1
                        node(:,nn) = node(:,db(temp_dn,temp_kp)) + (sn * vector(:))*scale_factor
                        if ( j == 0 ) then
						    fp = nn
						    if ( j == m-2) then
							    lp = nn
						    endif
					    elseif( j == m-2 ) then
						    lp = nn
					    endif
                    enddo
                endif
                
				!write (*,*) fp, lp
				temp_dbs = dbs(temp_dn)
				do j = 1, temp_dbs
					temp_db(j) = db(temp_dn, j)
				enddo

				dbs(temp_dn) = 0
				if ( temp_kp > cp ) then
								
					do j = 1, cp
						dbs(temp_dn) = dbs(temp_dn) + 1
						db(temp_dn,dbs(temp_dn)) = temp_db(j)
                    enddo
                    if (fp /= 0) then
					    do j = lp, fp, -1
						    dbs(temp_dn) = dbs(temp_dn) + 1
						    db(temp_dn,dbs(temp_dn)) = j
                        enddo
                    endif
					do j = temp_kp, temp_dbs
						dbs(temp_dn) = dbs(temp_dn) + 1
						db(temp_dn,dbs(temp_dn)) = temp_db(j)
					enddo

					dn = dn + 1
					do j = cp, temp_kp
						dbs(dn) = dbs(dn) + 1
						db(dn,dbs(dn)) = temp_db(j)
                    enddo
                    if (fp /= 0) then
					    do j = fp, lp
						    dbs(dn) = dbs(dn) + 1
						    db(dn,dbs(dn)) = j
					    enddo
					endif
					!if ( temp_ang(1) > cri_angle ) then
					!	temp_dn = dn
					!elseif ( temp_ang(2) > cri_angle ) then
					!	temp_dn = temp_dn
					!else 
						temp_dn = 0
					!endif
																
				else 
				
					do j = 1, temp_kp
						dbs(temp_dn) = dbs(temp_dn) + 1
						db(temp_dn,dbs(temp_dn)) = temp_db(j)
                    enddo
                    if (fp /= 0) then
					    do j = fp, lp
						    dbs(temp_dn) = dbs(temp_dn) + 1
						    db(temp_dn,dbs(temp_dn)) = j
                        enddo
                    endif
					do j = cp, temp_dbs
						dbs(temp_dn) = dbs(temp_dn) + 1
						db(temp_dn,dbs(temp_dn)) = temp_db(j)
					enddo

					dn = dn + 1
					do j = temp_kp, cp
						dbs(dn) = dbs(dn) + 1
						db(dn,dbs(dn)) = temp_db(j)
                    enddo
                    if (fp /= 0) then
					    do j = lp, fp, -1
						    dbs(dn) = dbs(dn) + 1
						    db(dn,dbs(dn)) = j
                        enddo
                    endif

					!if ( temp_ang(1) > cri_angle ) then
					!	temp_dn = temp_dn
					!elseif ( temp_ang(2) > cri_angle ) then
					!	temp_dn = dn
					!else 
						temp_dn = 0
					!endif

				endif

			endif

		endif

	endif
enddo
open(unit=3,file='./output/solid/remesh/partition_check.plt')
do i = 1,nn
	write (3,*) node(:,i)
enddo
close(3)

end subroutine partition
