!****************************************************************************
!
!  PROGRAM : In-Out
!
!  PURPOSE : 한 점의 위치가 다각형 내외부의 어느쪽인지 검사
!
!****************************************************************************
Subroutine inout(pn, p, cp, in_out)
	
implicit none

integer, intent(in) :: pn	
real, intent(in) :: p(pn,2), cp(2)
logical, intent(out) :: in_out

integer :: i, i_r
real :: total_ang, ang

write (*,*) cp
total_ang = 0
do i = 1, pn
	i_r = i + 1
	if ( i == pn ) i_r = 1

	call calc_angle( p(i_r,1:2), cp(1:2), p(i,1:2), ang )
	write (*,*) i, ang 

	if ( ang < 180 ) then
		total_ang = total_ang + ang
	else
		total_ang = total_ang + ( ang - 360 )
	endif
enddo

write (*,*) total_ang

if ( (total_ang >= 359) .AND. (total_ang <= 361.) ) then
	in_out = .TRUE.
else
	in_out = .FALSE.
endif

end subroutine



