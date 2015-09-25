!****************************************************************************
!
!  PROGRAM : One element
!
!  PURPOSE : 1. 3개의 점이 하나의 Element로 구성되어있는지 확인
!			 2. 1개의 Element를 구성하는 선분의 교차 검사
!
!****************************************************************************
! <1>
Subroutine calc_one_ele(ne, ele, a, b, c, judg)

implicit none

integer, intent(in) :: ne, a, b, c
integer, intent(in) :: ele(5,ne)
logical, intent(inout) :: judg

integer :: i, j, sum

do i = 1, ne
	sum = 0
	do j = 2, 5	
		if ( (ele(j,i) == a) .OR. (ele(j,i) == b) .OR. (ele(j,i) == c) ) then
			sum = sum + 1
		endif
	enddo
	
	if (sum == 3) then
		judg = .FALSE.  !<- 3 points consist of an element
		exit
	else
		judg = .TRUE. !<- 3 points don`t consist of an element
	endif

enddo

end Subroutine

! <2>
Subroutine calc_ele_cross(anes, nn, ne, ele, ci, di, i, ps, node, point, cross_judg)

implicit none

integer, intent(in) :: nn, ne, ci, di, i, ps, anes
integer, intent(in) :: ele(5,anes)
real, intent(inout) :: node(2,anes), point(2,ps)
integer, intent(inout) :: cross_judg

integer :: j, o, ci_r, di_r, count, i_r, i_r_r, i_l, i_l_l
integer :: a(4), ele_sum(2)
logical :: ja, jb

real :: tn(2,nn)

tn(:,1:nn) = node(:,1:nn)

tn(1,ci) = ( tn(1,ci) + tn(1,di) ) * 0.5
tn(2,ci) = ( tn(2,ci) + tn(2,di) ) * 0.5
tn(:,di) = tn(:,ci)

if ( i == 1 ) then
	ci_r = ci + 1;  di_r = di - 1;  i_r = i + 1;  i_r_r = i + 2;  i_l = ps  ;i_l_l = ps - 1
elseif ( i == 2 ) then
	ci_r = ci + 1;  di_r = di - 1;  i_r =i - 1;  i_r_r =ps;  i_l = i + 1;  i_l_l = i + 2
elseif ( i == ps ) then
	ci_r = ci + 1;  di_r = di - 1;  i_r = 1;  i_r_r = 2;  i_l = i - 1;  i_l_l = i - 2
elseif ( i == ps - 1 ) then
	ci_r = ci + 1;  di_r = di - 1;  i_r = i - 1;  i_r_r = i - 2;  i_l = i + 1;  i_l_l = 1
else
	!ci_r = ci - 1;  di_r = di + 1;  i_r = i - 1;  i_r_r = i - 2;  i_l = i + 1;  i_l_l = i + 2
    ci_r = ci + 1;  di_r = di - 1;  i_r = i - 1;  i_r_r = i - 2;  i_l = i + 1;  i_l_l = i + 2
endif

cross_judg = 0
count = 0
do j = 1, ne
	ele_sum = 0
	do o = 2, 5
		if ( ci == ele(o,j) ) then
			ele_sum(1) = ele_sum(1) + 1
		elseif ( ci_r == ele(o,j) ) then
			ele_sum(1) = ele_sum(1) + 1
		elseif ( di == ele(o,j) ) then
			ele_sum(2) = ele_sum(2) + 1
		elseif ( di_r == ele(o,j) ) then
			ele_sum(2) = ele_sum(2) + 1
		endif
	enddo
		
	if (ele_sum(1) == 2) then
		a = (/ ele(2,j), ele(3,j), ele(4,j), ele(5,j) /)
		call calc_cross(tn(:,a(1)), tn(:,a(2)), tn(:,a(3)), tn(:,a(4)), ja)
		if ( ja .EQV. .FALSE. ) cross_judg = cross_judg + 1 
		call calc_cross(tn(:,a(2)), tn(:,a(3)), tn(:,a(4)), tn(:,a(1)), jb)
		if ( jb .EQV. .FALSE. ) cross_judg = cross_judg + 1 
		count = count + 1
		
		if ( ( ja .EQV. .FALSE. ) .OR. ( jb .EQV. .FALSE. ) ) then
			point(:,i_r) = ( point(:,i_r_r) + point(:,i) ) * 0.5
			tn(:,ci) = point(:,i_r)
			node(:,ci) = point(:,i_r)
		endif
	elseif (ele_sum(2) == 2) then
		a = (/ ele(2,j), ele(3,j), ele(4,j), ele(5,j) /)
		call calc_cross(tn(:,a(1)), tn(:,a(2)), tn(:,a(3)), tn(:,a(4)), ja)
		if ( ja .EQV. .FALSE. ) cross_judg = cross_judg + 1 
		call calc_cross(tn(:,a(2)), tn(:,a(3)), tn(:,a(4)), tn(:,a(1)), jb)
		if ( jb .EQV. .FALSE. ) cross_judg = cross_judg + 1 
		count = count + 1
		
		if ( ( ja .EQV. .FALSE. ) .OR. ( jb .EQV. .FALSE. ) ) then
			point(:,i_l) = ( point(:,i_l_l) + point(:,i) ) * 0.5
			tn(:,di) = point(:,i_l)
			node(:,di) = point(:,i_l)
		endif
	endif

	if ( count == 2 ) then
		exit
	endif

enddo

end Subroutine
