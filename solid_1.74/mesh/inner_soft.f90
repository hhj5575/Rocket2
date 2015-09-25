!****************************************************************************
!
!  PROGRAM : Inner soft
!
!  PURPOSE : 내부 절점의 유연화
!
!****************************************************************************
subroutine inner_soft(anes, n_b, nn, ne, bs, b_num, node, ele, iter_check)

implicit none

integer, intent(in) :: anes, n_b, nn, ne, bs, ele(5,anes), b_num(bs)
real, intent(inout) :: node(2,anes)
logical, intent(inout) :: iter_check

integer :: i, j, k, soft_n, soft(2,12)
real :: soft_sum(2)
logical :: check

inner_soft_do: do i = n_b + 1, nn
	check = .TRUE.
	do j = 1, bs
		if ( i == b_num(j) ) then
			check = .FALSE.
			exit
		endif
	enddo

	if ( check ) then
		soft_n = 0
		soft_sum(1) = 0.
		soft_sum(2) = 0.

		do j = 1, ne
			do k = 2, 5
				if (i == ele(k,j)) then 
					soft_n = soft_n + 1
                    if (soft_n > 12) then
                        iter_check = .false.
                        exit inner_soft_do
                    endif
					if (k .EQ. 5) then
						soft(1, soft_n) = ele(2,j)
						soft(2, soft_n) = ele(4,j)
					else if (k .EQ. 2) then
						soft(1, soft_n) = ele(3,j)
						soft(2, soft_n) = ele(5,j)
					else
						soft(1, soft_n) = ele(k+1,j)
						soft(2, soft_n) = ele(K-1,j)
					endif
				endif
            enddo
		enddo

		if (sum(soft(1,1:soft_n)) == sum(soft(2,1:soft_n))) then
			
			do j = 1, soft_n
				soft_sum(1) = soft_sum(1) + node(1,soft(1,j))
				soft_sum(2) = soft_sum(2) + node(2,soft(1,j))
			enddo

			node(:,i) = soft_sum(:) / soft_n
						
		endif
		
		if ( soft_n >= 8 ) then
            write (*,*) 'Fail - smoothing subroutine'
		    iter_check = .FALSE.
		    exit
		endif
		
	endif

enddo inner_soft_do

! 유연화 종료
end subroutine inner_soft

subroutine inner_soft_end(anes, n_b, nn, ne, node, ele, iter_check)

implicit none

integer, intent(in) :: anes, nn, ne, n_b, ele(5,anes)
real, intent(inout) :: node(2,anes)
logical, intent(inout) :: iter_check

integer :: i, j, k, soft_n, soft(2,12)
real :: soft_sum(2)

do i = n_b + 1, nn
	soft_n = 0
	soft_sum(1) = 0.
	soft_sum(2) = 0.

	do j = 1, ne
		do k = 2, 5
			if (ele(k,j) == i) then 
				soft_n = soft_n + 1
				if (k == 5) then
					soft(1, soft_n) = ele(2,j)
					soft(2, soft_n) = ele(4,j)
				else if (k == 2) then
					soft(1, soft_n) = ele(3,j)
					soft(2, soft_n) = ele(5,j)
				else
					soft(1, soft_n) = ele(k+1,j)
					soft(2, soft_n) = ele(K-1,j)
				endif
			endif
		enddo
	enddo
	
	if (sum(soft(1,1:soft_n)) == sum(soft(2,1:soft_n))) then
		
		do j = 1, soft_n
			soft_sum(:) = soft_sum(:) + node(:,soft(1,j))
		enddo

		node(:,i) = soft_sum(:) / soft_n
					
	endif
	
	if ( soft_n >= 8 ) then
        write (*,*) 'Fail - smoothing subroutine'
	    iter_check = .FALSE.
	    exit
	endif
	
enddo

! 유연화 종료
end subroutine inner_soft_end

subroutine inner_soft_object(anes, pnn, pne, nn, ne, node, ele, iter_check)

implicit none

integer, intent(in) :: anes, nn, ne, pnn, pne, ele(5,anes)
real, intent(inout) :: node(2,anes)
logical, intent(inout) :: iter_check

integer :: i, j, k, soft_n, roof_i
integer :: nni(16,nn), nni_n(nn), soft(2,16)
integer :: nei_n(nn), nei1(8,nn), nei2(8,nn), nei3(8,nn)
real :: soft_sum(2)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375
logical :: check(nn)

nni = 0;  nni_n = 0; nei_n=0; nei1 = 0; nei2 = 0; nei3 = 0
check = .FALSE.

do i = pnn + 1, nn
	soft_n = 0
	soft_sum(1) = 0.
	soft_sum(2) = 0.
	
	do j = pne, ne
	    do k = 2, 5
		    if (ele(k,j) == i) then
                soft_n = soft_n + 1
                nei_n(i) = nei_n(i) + 1
                if (k == 2) then
                    soft(1, soft_n) = ele(3,j)
                    soft(2, soft_n) = ele(5,j)
                    nei1(nei_n(i), i) = ele(3,j)
                    nei2(nei_n(i), i) = ele(4,j)
                    nei3(nei_n(i), i) = ele(5,j)
                elseif (k == 3) then
                    soft(1, soft_n) = ele(k+1,j)
                    soft(2, soft_n) = ele(k-1,j)
                    nei1(nei_n(i), i) = ele(4,j)
                    nei2(nei_n(i), i) = ele(5,j)
                    nei3(nei_n(i), i) = ele(2,j)
                elseif (k == 4) then
                    soft(1, soft_n) = ele(k+1,j)
                    soft(2, soft_n) = ele(k-1,j)
                    nei1(nei_n(i), i) = ele(5,j)
                    nei2(nei_n(i), i) = ele(2,j)
                    nei3(nei_n(i), i) = ele(3,j)
                elseif (k == 5) then
                    soft(1, soft_n) = ele(2,j)
                    soft(2, soft_n) = ele(4,j)
                    nei1(nei_n(i), i) = ele(2,j)
                    nei2(nei_n(i), i) = ele(3,j)
                    nei3(nei_n(i), i) = ele(4,j)                    
                endif
                soft_n = soft_n + 1
                if ( k == 4 ) then
                    soft(1:2, soft_n) = ele(2,j)
                elseif ( k == 5 ) then
                    soft(1:2, soft_n) = ele(3,j)
                else
                    soft(1:2, soft_n) = ele(k+2,j)
                endif
            endif
        enddo
    enddo
	
    if (sum(soft(1,1:soft_n)) == sum(soft(2,1:soft_n))) then
		check(i) = .TRUE.
		nni_n(i) = soft_n
		nni(1:nni_n(i),i) = soft(1,1:nni_n(i))
    else
        write (*,*) 'Fail - smoothing subroutine'
        iter_check = .FALSE.
        exit
	endif
enddo

if ( iter_check ) then
    do roof_i = 1, 6
        do i = pnn + 1, nn
            soft_sum = 0.0
            do j = 1, nei_n(i)
                soft_sum = soft_sum + node(:,nei1(j,i))
            enddo
            node(:,i) = soft_sum(:)/(nei_n(i)*1.0)
        enddo
        !do i = pnn + 1, nn
        !    soft_sum = 0.0
        !    do j = 1, nei_n(i)
        !        soft_sum = soft_sum + (node(:,nei1(j,i))+node(:,nei3(j,i))-node(:,nei2(j,i)))
        !    enddo
        !    node(:,i) = soft_sum(:)/(nei_n(i)*1.0)
        !enddo
    enddo
    !do roof_i = 1, 6
	   ! do i = 
		  !  if ( check(i) ) then
			 !   soft_sum = 0.
			 !   do j = 1, nni_n(i)
				!    soft_sum(:) = soft_sum(:) + node(:,nni(j,i))
    !            enddo	
    !            node(:,i) = soft_sum(:)/nni_n(i)
		  !  endif
	   ! enddo
    !enddo
endif

! 유연화 종료
end subroutine inner_soft_object

subroutine plane_soft(nn, ne, node, ele)

implicit none

integer, intent(in) :: nn, ne, ele(5,ne)
real, intent(inout) :: node(2,nn)

integer :: i, j, k, soft_n, soft(2,8)
real :: soft_sum(2)

do i = 1, nn
	soft_n = 0
	soft_sum(1) = 0.
	soft_sum(2) = 0.

	do j = 1, ne
		do k = 2, 5
			if (ele(k,j) == i) then 
				soft_n = soft_n + 1
				if (k == 5) then
					soft(1, soft_n) = ele(2,j)
					soft(2, soft_n) = ele(4,j)
				else if (k == 2) then
					soft(1, soft_n) = ele(3,j)
					soft(2, soft_n) = ele(5,j)
				else
					soft(1, soft_n) = ele(k+1,j)
					soft(2, soft_n) = ele(K-1,j)
				endif
			endif
		enddo
	enddo
	
	if (sum(soft(1,1:soft_n)) == sum(soft(2,1:soft_n))) then
		
		do j = 1, soft_n
			soft_sum(:) = soft_sum(:) + node(:,soft(1,j))
		enddo

		node(:,i) = soft_sum(:) / soft_n
					
	endif
	
enddo

! 유연화 종료
end subroutine