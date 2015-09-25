!****************************************************************************
!
!  PROGRAM : Paving
!
!  PURPOSE : 경계 절점의 내부로 투영	
!
!****************************************************************************
subroutine paving(temp_bs, d, bs, b_num, cur_b, pbn, me, en, ep, iter_check)

implicit none

integer, intent(in) :: bs, temp_bs
integer, intent(inout) :: b_num(bs), pbn, en, ep(bs)
real, intent(in) :: d
real, intent(inout) :: me(2,temp_bs), cur_b(2,bs)
logical, intent(inout) :: iter_check

integer :: i, j, k, i_l_l, i_l, i_r, i_r_r, conn_num
integer :: e_first, e_second, e_f, this(2), conn_b_me(2,bs), conn_b(bs)
real :: alpha, length, angle_a, angle_b, line_d, first_line_d
real :: vector(4), unit_vector(2), dis, seg(3), b(2,bs)
real :: line_dis(2)
logical :: judg 

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

j = 0
k = 0
pbn = 0
e_f = 0
e_first = 0
e_second = 0
en = 0
me = 0
conn_b_me = 0
conn_num = 0
line_dis = 0.0

! find starting point
call find_b_start_point(bs, b_num, cur_b, b, conn_b)

line_d = d
do i = 1, bs ! <- paving 경계 전체 루프문 시작
	if (i+j > bs) then
		exit	
	elseif (i+j == 1) then   ! <- 처음 절점
		i_l = bs
		i_r = i + j + 1
	elseif (i+j == bs) then

		if (k == 1) then
	
			call calc_angle(b(:,i+j-1), b(:,i+j), b(:,1), alpha) 
			length = sqrt((b(1,1) - b(1,i+j)) ** 2 + (b(2,1) - b(2,i+j)) ** 2)
			vector(1:2) = b(:,1) - b(:,i+j)
			unit_vector(:) = vector(1:2) / length
					
			if ((alpha > 240.) .AND. (alpha < 360.)) then
							
				!pbn = pbn + 1
				!vector(1) = unit_vector(1) * d
				!vector(2) = unit_vector(2) * d
				!vector(3) = unit_vector(1) * d * sqrt(2.)
				!vector(4) = unit_vector(2) * d * sqrt(2.)
				!me(1,pbn) = b(1,i+j) + (cos(alpha*pi/270.) * vector(1) - sin(alpha*pi/270.) * vector(2))
				!me(2,pbn) = b(2,i+j) + (sin(alpha*pi/270.) * vector(1) + cos(alpha*pi/270.) * vector(2))
				!me(1,pbn+1) = b(1,i+j) + (cos(alpha*pi/360.) * vector(3) - sin(alpha*pi/360.) * vector(4))
				!me(2,pbn+1) = b(2,i+j) + (sin(alpha*pi/360.) * vector(3) + cos(alpha*pi/360.) * vector(4))
				!pbn = pbn + 1
                
                this = (/ i+j, 1 /)
                pbn = pbn + 1
                vector(1:2) = b(:,i+j) - b(:,1)
                me(:,pbn) = b(:,i+j) + vector(1:2)
                vector(1:2) = b(:,i+j) - b(:,i_l)
                me(:,pbn+1) = me(:,pbn) + vector(1:2)
                pbn = pbn + 1
		        
                line_d = sqrt((vector(1))**2 + (vector(2))**2)
                line_d = line_d * cos(abs(270-alpha)*pi/180.)
                alpha = alpha - 180
                
                if (conn_num /= 0) then
                    call adjust_line_d(1, first_line_d, alpha, temp_bs, bs, this, conn_num, conn_b_me, line_dis, b, me, iter_check)
                    if (iter_check == .false.) exit
                endif
                
			!elseif ((alpha >= 315.) .AND. (alpha < 360.)) then
			!	
			!	pbn = pbn + 1
			!	vector(1) = unit_vector(1) * d
			!	vector(2) = unit_vector(2) * d
			!	vector(3) = unit_vector(1) * d * sqrt(2.)
			!	vector(4) = unit_vector(2) * d * sqrt(2.)
			!	me(1,pbn) = b(1,i+j) + (cos(alpha*pi/240.) * vector(1) - sin(alpha*pi/240.) * vector(2))
			!	me(2,pbn) = b(2,i+j) + (sin(alpha*pi/240.) * vector(1) + cos(alpha*pi/240.) * vector(2))
			!	me(1,pbn+1) = b(1,i+j) + (cos(alpha*pi/288.) * vector(1) - sin(alpha*pi/288.) * vector(2))
			!	me(2,pbn+1) = b(2,i+j) + (sin(alpha*pi/288.) * vector(1) + cos(alpha*pi/288.) * vector(2))
			!	me(1,pbn+2) = b(1,i+j) + (cos(alpha*pi/360.) * vector(3) - sin(alpha*pi/360.) * vector(4))
			!	me(2,pbn+2) = b(2,i+j) + (sin(alpha*pi/360.) * vector(3) + cos(alpha*pi/360.) * vector(4))
			!	me(1,pbn+3) = b(1,i+j) + (cos(alpha*pi/480.) * vector(1) - sin(alpha*pi/480.) * vector(2))
			!	me(2,pbn+3) = b(2,i+j) + (sin(alpha*pi/480.) * vector(1) + cos(alpha*pi/480.) * vector(2))
			!	pbn = pbn + 3

			endif

			exit

		else 
			i_l = i + j - 1
			i_r = 1
		endif
	else 
		i_l = i + j - 1
		i_r = i + j + 1
	endif

	call calc_angle(b(:,i_l), b(:,i+j), b(:,i_r), alpha)
	length = sqrt((b(1,i_r) - b(1,i+j)) ** 2 + (b(2,i_r) - b(2,i+j)) ** 2)
	call calc_len( b(:,i+j), b(:,i_l), seg(1) )
	call calc_len( b(:,i+j), b(:,i_r), seg(2) )
	seg(3) = (seg(1) + seg(2)) * 0.5
	!if ( seg(3) > d ) seg(3) = d
	vector(1) = b(1,i_r) - b(1,i+j)
	vector(2) = b(2,i_r) - b(2,i+j)
	unit_vector(1) = vector(1) / length
	unit_vector(2) = vector(2) / length
	
	if ((alpha >= 0.) .AND. (alpha <= 110.)) then 
		this = (/ i+j, i_r /)
		if (i+j == 1) then
			
			call calc_angle(b(:,1), b(:,2), b(:,3), angle_a)
			call calc_angle(b(:,bs-1), b(:,bs), b(:,1), angle_b)
			
			if ( angle_b < angle_a) angle_a = angle_b
			
			if ( ( angle_a > 90 ) .AND. ( angle_a <= 160 ) ) then
				dis = sin((angle_a-65)*pi/180) * sin(alpha*pi/180) * sqrt(2.) * seg(3)
			else 
				dis = sin(alpha*pi/180) * sqrt(2.) * seg(3)
			endif

			pbn = pbn + 1
			!vector(1) = unit_vector(1) * dis
			!vector(2) = unit_vector(2) * dis
            !me(1,pbn) = b(2,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2))
			!me(2,pbn) = b(3,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2))
            
            vector(1:2) = b(:,i_r)+b(:,i_l)-(2.0*b(:,i+j))
            length = sqrt(vector(1)**2.0 + vector(2)**2.0)
            vector(1:2) = vector(1:2)/length*dis
            me(1,pbn) = b(1,i+j) + vector(1)
            me(2,pbn) = b(2,i+j) + vector(2)
            
            vector(1:2) = b(:,i_l) - b(:,i+j)
            line_d = sqrt((vector(1))**2.0 + (vector(2))**2.0)
            line_d = line_d * cos(abs(90-alpha)*pi/180.)
            
			j = j + 1;  k = 1;  e_first = 1;  e_f = 1

			en = en + 1
			ep(en) = pbn

		elseif (i+j == bs) then

			call calc_angle(b(:,i+j), b(:,1), b(:,2), angle_a)
			call calc_angle(b(:,i+j-2), b(:,i+j-1), b(:,i+j), angle_b)

			if ( angle_b < angle_a) angle_a = angle_b
			
			if ( ( angle_a > 90 ) .AND. ( angle_a <= 160 ) ) then
				dis = sin((angle_a-65)*pi/180) * sin(alpha*pi/180) * sqrt(2.) * seg(3)
			else 
				dis = sin(alpha*pi/180) * sqrt(2.) * seg(3)
			endif
		
			vector(1) = unit_vector(1) * dis
			vector(2) = unit_vector(2) * dis
			
			if (e_f == 1) then
				
				!me(1,1) = (b(2,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2)) + me(1,pbn)) / 2.
				!me(2,1) = (b(3,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2)) + me(2,pbn)) / 2.
                vector(1:2) = b(:,i_r)+b(:,i_l)-(2.0*b(:,i+j))
                length = sqrt(vector(1)**2.0 + vector(2)**2.0)
                vector(1:2) = vector(1:2)/length*dis
                me(1,pbn) = (b(1,i+j) + vector(1) + me(1,pbn)) / 2.
                me(2,pbn) = (b(2,i+j) + vector(2) + me(2,pbn)) / 2.
                
                vector(1:2) = b(:,i_l) - b(:,i+j)
                line_d = sqrt((vector(1))**2 + (vector(2))**2)
                line_d = line_d * cos(abs(90-alpha)*pi/180.)

			elseif (e_second == 1) then
			
				!me(1,1) = (b(2,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2)) + me(1,1)) / 2.
				!me(2,1) = (b(3,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2)) + me(2,1)) / 2.
                vector(1:2) = b(:,i_r)+b(:,i_l)-(2.0*b(:,i+j))
                length = sqrt(vector(1)**2.0 + vector(2)**2.0)
                vector(1:2) = vector(1:2)/length*dis
                me(1,1) = (b(1,i+j) + vector(1) + me(1,1)) / 2.
                me(2,1) = (b(2,i+j) + vector(2) + me(2,1)) / 2.
                
                vector(1:2) = b(:,i_l) - b(:,i+j)
                line_d = sqrt((vector(1))**2 + (vector(2))**2)
                line_d = line_d * cos(abs(90-alpha)*pi/180.)
			
			else
			
				!me(1,1) = b(2,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2))
				!me(2,1) = b(3,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2))
                vector(1:2) = b(:,i_r)+b(:,i_l)-(2.0*b(:,i+j))
                length = sqrt(vector(1)**2.0 + vector(2)**2.0)
                vector(1:2) = vector(1:2)/length*dis
                me(1,1) = b(1,i+j) + vector(1)
                me(2,1) = b(2,i+j) + vector(2)
                
                vector(1:2) = b(:,i_l) - b(:,i+j)
                line_d = sqrt((vector(1))**2 + (vector(2))**2)
                line_d = line_d * cos(abs(90-alpha)*pi/180.)
                
				en = en + 1
				ep(en) = 1

			endif

			pbn = pbn - 1

			call calc_cross(b(:,i+j-1), me(:,1), b(:,i+j-2), me(:,pbn), judg)

			if (judg .EQV. .FALSE.) then
			
				me(1,1) = (1.5 * me(1,pbn) + b(1,i+j)) / 2.5
				me(2,1) = (1.5 * me(2,pbn) + b(2,i+j)) / 2.5

			endif 

		else 
			
			if (i+j == bs - 1) then
				i_r_r = 1
				i_l_l = i + j - 2
			elseif (i+j == 2) then
				i_r_r = i + j + 2
				i_l_l = bs
			else 
				i_r_r = i + j + 2
				i_l_l = i + j - 2
			endif
			
			call calc_angle(b(:,i+j), b(:,i_r), b(:,i_r_r), angle_a)
			call calc_angle(b(:,i_l_l), b(:,i_l), b(:,i+j), angle_b)
			if ( angle_b < angle_a) then
				angle_a = angle_b
			endif
			
			if ( ( angle_a >= 90 ) .AND. ( angle_a <= 160 ) ) then
				dis = sin((angle_a-65)*pi/180) * sin(alpha*pi/180) * sqrt(2.) * seg(3)
			else 
				dis = sin(alpha*pi/180) * sqrt(2.) * seg(3)
			endif
						
			vector(1) = unit_vector(1) * dis
			vector(2) = unit_vector(2) * dis

			if (e_f == 1) then
				
				!me(1,pbn) = (b(2,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2)) + me(1,pbn)) / 2.
				!me(2,pbn) = (b(3,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2)) + me(2,pbn)) / 2.
                vector(1:2) = b(:,i_r)+b(:,i_l)-(2.0*b(:,i+j))
                length = sqrt(vector(1)**2.0 + vector(2)**2.0)
                vector(1:2) = vector(1:2)/length*dis
                me(1,pbn) = (b(1,i+j) + vector(1) + me(1,pbn)) / 2.
                me(2,pbn) = (b(2,i+j) + vector(2) + me(2,pbn)) / 2.
                
                vector(1:2) = b(:,i_l) - b(:,i+j)
                !vector(1:2) = me(:,pbn) - b(:,i_r)
                line_d = sqrt((vector(1))**2 + (vector(2))**2)
                line_d = line_d * cos(abs(90-alpha)*pi/180.)
					
			elseif (i == 2) then

				!me(1,pbn) = b(2,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2))
				!me(2,pbn) = b(3,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2))
                vector(1:2) = b(:,i_r)+b(:,i_l)-(2.0*b(:,i+j))
                length = sqrt(vector(1)**2.0 + vector(2)**2.0)
                vector(1:2) = vector(1:2)/length*dis
                me(1,1) = b(1,i+j) + vector(1)
                me(2,2) = b(2,i+j) + vector(2)
                
                vector(1:2) = b(:,i_l) - b(:,i+j)
                line_d = sqrt((vector(1))**2 + (vector(2))**2)
                line_d = line_d * cos(abs(90-alpha)*pi/180.)

				e_second = 1

				en = en + 1
				ep(en) = pbn
			
			elseif ((i+j == bs - 1) .AND. (e_first == 1)) then
			
				!me(1,1) = (b(2,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2)) + me(1,1)) / 2.
				!me(2,1) = (b(3,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2)) + me(2,1)) / 2.
                vector(1:2) = b(:,i_r)+b(:,i_l)-(2.0*b(:,i+j))
                length = sqrt(vector(1)**2.0 + vector(2)**2.0)
                vector(1:2) = vector(1:2)/length*dis
                me(1,1) = (b(1,i+j) + vector(1) + me(1,1)) / 2.
                me(2,1) = (b(2,i+j) + vector(2) + me(2,1)) / 2.
                
                vector(1:2) = b(:,i_l) - b(:,i+j)
                !vector(1:2) = me(:,1) - b(:,i_r)
                line_d = sqrt((vector(1))**2 + (vector(2))**2)
                line_d = line_d * cos(abs(90-alpha)*pi/180.)
                
				pbn = pbn - 1

			else

				!me(1,pbn) = b(2,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2))
				!me(2,pbn) = b(3,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2))
                vector(1:2) = b(:,i_r)+b(:,i_l)-(2.0*b(:,i+j))
                length = sqrt(vector(1)**2.0 + vector(2)**2.0)
                vector(1:2) = vector(1:2)/length*dis
                me(1,pbn) = b(1,i+j) + vector(1)
                me(2,pbn) = b(2,i+j) + vector(2)
                
                vector(1:2) = b(:,i_l) - b(:,i+j)
                line_d = sqrt((vector(1))**2 + (vector(2))**2)
                line_d = line_d * cos(abs(90-alpha)*pi/180.)                
                
				en = en + 1
				ep(en) = pbn

			endif

			if ((i+j == 2) .OR. ((i+j == bs - 1) .AND. (e_first == 1)) .OR. (e_f == 1)) then

			else

				call calc_cross(b(:,i+j-1), me(:,pbn), b(:,i+j-2), me(:,pbn-1), judg)

				if (judg .EQV. .FALSE.) then
				
					me(1,pbn) = (1.5 * me(1,pbn-1) + b(1,i+j)) / 2.5
					me(2,pbn) = (1.5 * me(2,pbn-1) + b(2,i+j)) / 2.5
	
				endif 

			endif
						
			j = j + 1
			e_f = 1	

        endif
        
        if (conn_num == 0) then
            line_dis(1) = line_d
            if (this(1) == 1) then
                first_line_d = sqrt((b(1,i_r)-b(1,this(1)))**2 + (b(2,i_r)-b(2,this(1)))**2)
                first_line_d = first_line_d * cos(abs(90-alpha)*pi/180.) 
            endif
        else
            call adjust_line_d(0, line_d, alpha, temp_bs, bs, this, conn_num, conn_b_me, line_dis, b, me, iter_check)
            if (iter_check == .false.) exit
        endif
                
    elseif ((alpha > 110.) .AND. (alpha <= 240.)) then

		pbn = pbn + 1
        if ( pbn > temp_bs ) then
            !write (*,*) 'pbn, temp_bs:', pbn, temp_bs
            write (*,*) 'Fail - paving subroutine'
		    iter_check = .FALSE.
		    exit
		endif
        
        conn_num = conn_num + 1
        conn_b_me(:,conn_num) = (/ i+j, pbn /)
		!vector(1) = unit_vector(1) * seg(3) / sqrt( sin(alpha*pi/360.) )
		!vector(2) = unit_vector(2) * seg(3) / sqrt( sin(alpha*pi/360.) )
        vector(1) = unit_vector(1) * line_d / sin(alpha*pi/360.)
		vector(2) = unit_vector(2) * line_d / sin(alpha*pi/360.)
		me(1,pbn) = b(1,i+j) + (cos(alpha*pi/360.) * vector(1) - sin(alpha*pi/360.) * vector(2))
		me(2,pbn) = b(2,i+j) + (sin(alpha*pi/360.) * vector(1) + cos(alpha*pi/360.) * vector(2))
        this = (/ i+j, i_r /)

		if (e_f == 1) then
			
			call calc_cross(b(:,i+j), me(:,pbn), b(:,i+j-1), me(:,pbn-1), judg)
			if (judg .EQV. .FALSE.) then
				me(:,pbn-1) = (1.5 * me(:,pbn) + b(:,i+j-2)) / 2.5
			endif
		
		elseif ( (i+j == bs - 1) .AND. (e_first == 1) ) then

			call calc_cross(b(:,i+j), me(:,pbn), b(:,i+j+1), me(:,1), judg)
			if (judg .EQV. .FALSE.) then
				me(:,1) = (1.5 * me(:,pbn) + b(:,1)) / 2.5
            endif 
            if (conn_num /= 0) then
                call adjust_line_d(1, first_line_d, alpha, temp_bs, bs, this, conn_num, conn_b_me, line_dis, b, me, iter_check)
                if (iter_check == .false.) exit
            endif

		elseif ( (i+j == bs) .AND. (e_second == 1) ) then

			call calc_cross(b(:,i+j), me(:,pbn), b(:,1), me(:,1), judg)

			if (judg .EQV. .FALSE.) then
				
				me(:,1) = (1.5 * me(:,pbn) + b(:,2)) / 2.5

			endif
			
		else 
			
			if ( i+j .NE. 1) then
				call calc_cross(b(:,i+j-1), me(:,pbn-1), b(:,i+j), me(:,pbn), judg)	
                if (judg .EQV. .FALSE.) then
					me(:,pbn) = (1.5 * me(:,pbn) + b(:,i+j)) / 2.5
				endif
			endif

		endif

		e_f = 0 

    elseif ((alpha > 240.) .AND. (alpha < 360.)) then
        !vector(1) = unit_vector(1) * seg(3)
		!vector(2) = unit_vector(2) * seg(3)
		!vector(3) = unit_vector(1) * seg(3) * sqrt(2.)
		!vector(4) = unit_vector(2) * seg(3) * sqrt(2.)
		!me(1,pbn) = b(1,i+j) + (cos(alpha*pi/270.) * vector(1) - sin(alpha*pi/270.) * vector(2))
		!me(2,pbn) = b(2,i+j) + (sin(alpha*pi/270.) * vector(1) + cos(alpha*pi/270.) * vector(2))
		!me(1,pbn+1) = b(1,i+j) + (cos(alpha*pi/360.) * vector(3) - sin(alpha*pi/360.) * vector(4))
		!me(2,pbn+1) = b(2,i+j) + (sin(alpha*pi/360.) * vector(3) + cos(alpha*pi/360.) * vector(4))
		!me(1,pbn+2) = b(1,i+j) + (cos(alpha*pi/540.) * vector(1) - sin(alpha*pi/540.) * vector(2))
		!me(2,pbn+2) = b(2,i+j) + (sin(alpha*pi/540.) * vector(1) + cos(alpha*pi/540.) * vector(2))
		!pbn = pbn + 2
		!e_f = 0 
        
        this = (/ i+j, i_r /)
        pbn = pbn + 1
        vector(1:2) = b(:,i+j) - b(:,i_r)
        me(:,pbn) = b(:,i+j) + vector(1:2)
        vector(1:2) = b(:,i+j) - b(:,i_l)
        me(:,pbn+1) = me(:,pbn) + vector(1:2)
        me(:,pbn+2) = b(:,i+j) + vector(1:2)
        pbn = pbn + 2
		e_f = 0 
        line_d = sqrt((vector(1))**2 + (vector(2))**2)
        line_d = line_d * cos(abs(270-alpha)*pi/180.)
        alpha = alpha - 180
        
        if (conn_num == 0) then
            line_dis(1) = line_d
            if (this(1) == 1) then
                first_line_d = sqrt((b(1,i_r)-b(1,this(1)))**2 + (b(2,i_r)-b(2,this(1)))**2)
                first_line_d = first_line_d * cos(abs(90-alpha)*pi/180.) 
            endif
        else
            call adjust_line_d(0, line_d, alpha, temp_bs, bs, this, conn_num, conn_b_me, line_dis, b, me, iter_check)
            if (iter_check == .false.) exit
        endif

	!elseif ((alpha >= 315.) .AND. (alpha < 360.)) then

		!pbn = pbn + 1
		!vector(1) = unit_vector(1) * seg(3)
		!vector(2) = unit_vector(2) * seg(3)
		!vector(3) = unit_vector(1) * seg(3) * sqrt(2.)
		!vector(4) = unit_vector(2) * seg(3) * sqrt(2.)
		!me(1,pbn) = b(1,i+j) + (cos(alpha*pi/240.) * vector(1) - sin(alpha*pi/240.) * vector(2))
		!me(2,pbn) = b(2,i+j) + (sin(alpha*pi/240.) * vector(1) + cos(alpha*pi/240.) * vector(2))
		!me(1,pbn+1) = b(1,i+j) + (cos(alpha*pi/288.) * vector(1) - sin(alpha*pi/288.) * vector(2))
		!me(2,pbn+1) = b(2,i+j) + (sin(alpha*pi/288.) * vector(1) + cos(alpha*pi/288.) * vector(2))
		!me(1,pbn+2) = b(1,i+j) + (cos(alpha*pi/360.) * vector(3) - sin(alpha*pi/360.) * vector(4))
		!me(2,pbn+2) = b(2,i+j) + (sin(alpha*pi/360.) * vector(3) + cos(alpha*pi/360.) * vector(4))
		!me(1,pbn+3) = b(1,i+j) + (cos(alpha*pi/480.) * vector(1) - sin(alpha*pi/480.) * vector(2))
		!me(2,pbn+3) = b(2,i+j) + (sin(alpha*pi/480.) * vector(1) + cos(alpha*pi/480.) * vector(2))
		!me(1,pbn+4) = b(1,i+j) + (cos(alpha*pi/720.) * vector(1) - sin(alpha*pi/720.) * vector(2))
		!me(2,pbn+4) = b(2,i+j) + (sin(alpha*pi/720.) * vector(1) + cos(alpha*pi/720.) * vector(2))
		!pbn = pbn + 4
		!e_f = 0
        
	endif
enddo
do i = 1, bs
    b_num(i) = conn_b(i)
    cur_b(:,i) = b(:,i)
enddo

!write (*,*) "----------------------------------------------"
!write (*,'(A,I3)') "끝절점(en) : ", en
!write (*,'(A,I3)') "투영된 절점(pbn) : ", pbn
!do i = 1, pbn
!	write (*,*) me(i,:)
!enddo
!write (*,*) "----------------------------------------------"

end subroutine

subroutine find_b_start_point(bs, b_num, cur_b, b, conn_b)

implicit none

integer, intent(in) :: bs, b_num(bs)
real, intent(in) :: cur_b(2,bs)
real, intent(inout) :: b(2,bs)
integer, intent(inout) :: conn_b(bs)

integer :: i, st, num, this(3)
real :: ang
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

st = 1
do i = 1, bs
    this = (/ i-1, i, i+1 /)
    if (i == 1) this(1) = bs
    if (i == bs) this(3) = 1
    
    call calc_angle(cur_b(:,this(1)), cur_b(:,this(2)), cur_b(:,this(3)), ang)
    
    if ((ang >= 0.) .AND. (ang <= 110.)) then
        st = i
        exit
    endif
enddo

num = 0
b = 0.0
conn_b = 0
do i = st, bs
    num = num + 1
    conn_b(num) = b_num(i)
    b(:,num) = cur_b(:,i)
enddo
do i = 1, st-1
    num = num + 1
    conn_b(num) = b_num(i)
    b(:,num) = cur_b(:,i)
enddo
if (num /= bs) stop 'error-find_b_start_point(paving): new bs is not equal original bs'

end subroutine find_b_start_point


subroutine adjust_line_d(flag, line_d, alpha, temp_bs, bs, this, conn_num, conn_b_me, line_dis, b, me, iter_check)

implicit none

!flag = 1: last(line_d = first_line_d)
!     = 0: not last
integer, intent(in) :: flag, temp_bs, bs, this(2)
integer, intent(inout) :: conn_num, conn_b_me(2,bs)
real, intent(in) :: line_d, alpha
real, intent(inout) :: me(2,temp_bs), b(2,bs), line_dis(2)
logical, intent(inout) :: iter_check

integer :: i, num
real :: ratio, length, vector(2)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

if (flag == 1) then
    line_dis(2) = line_d
    num = conn_num+1
else
    line_dis(2) = sqrt((b(1,this(2))-b(1,this(1)))**2 + (b(2,this(2))-b(2,this(1)))**2)
    line_dis(2) = line_dis(2) * cos(abs(90-alpha)*pi/180.) 
    num = conn_num
endif
ratio = (line_dis(2) - line_dis(1))/(num)

if (abs(line_dis(1)) <= 10e-8) then
    iter_check = .false.
else
    do i = 1, num-1
        vector = me(:,conn_b_me(2,i)) - b(:,conn_b_me(1,i))
        length = sqrt((vector(1))**2 + (vector(2))**2)
        !vector = (vector/length)*(ratio*i+line_dis(1))
        vector = vector*(ratio*i/line_dis(1)+1)
        me(:,conn_b_me(2,i)) = b(:,conn_b_me(1,i)) + vector
    enddo
    conn_num = 0
    conn_b_me = 0
    if (flag == 0) line_dis(1) = line_d
endif

end subroutine adjust_line_d