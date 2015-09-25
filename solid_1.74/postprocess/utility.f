subroutine find_inter_skin( nn, nb, st, ed, num_node, skin_info, bound_node ) 

implicit none

integer, intent(in) :: nn, nb, st, ed, skin_info(nn)
integer, intent(out) :: num_node
integer, intent(inout) :: bound_node(nn)

integer :: i, position 
logical :: find_next

find_next = .TRUE.

do i = 1, nb
	if ( skin_info(i) == st ) then
		position = i 
		exit
	endif
enddo

bound_node(1) = skin_info(position)
num_node = 1 
do i = position+1, nb
	num_node = num_node + 1 
	bound_node(num_node) = skin_info(i)
	if ( skin_info(i) == ed ) then
		find_next = .FALSE.
		exit
	endif
enddo

if ( find_next ) then 
	do i= 1, position-1
		num_node = num_node + 1 
		bound_node(num_node) = skin_info(i)
		if ( skin_info(i) == ed ) then
			find_next = .FALSE.
			exit
		endif
	enddo 
endif 

end subroutine

!=========================================================================

subroutine find_inter_skin_normal( nn, nb, st, ed, num_node, skin_info, skin_normal, bound_normal ) 

implicit none

integer, intent(in) :: nn, nb, st, ed, skin_info(nn)
integer, intent(out) :: num_node
real, intent(in) :: skin_normal(2,nn)
real, intent(inout) :: bound_normal(2,nn)

integer :: i, position 
logical :: find_next

find_next = .TRUE.

do i = 1, nb
	if ( skin_info(i) == st ) then
		position = i 
		exit
	endif
enddo

bound_normal(:,1) = skin_normal(:,position)
num_node = 1 
do i = position+1, nb
	num_node = num_node + 1 
	bound_normal(:,num_node) = skin_normal(:,i)
	if ( skin_info(i) == ed ) then
		find_next = .FALSE.
		exit
	endif
enddo

if ( find_next ) then 
	do i= 1, position-1
		num_node = num_node + 1 
		bound_normal(:,num_node) = skin_normal(:,i)
		if ( skin_info(i) == ed ) then
			find_next = .FALSE.
			exit
		endif
	enddo 
endif 

bound_normal = -(bound_normal)

end subroutine

!=========================================================================

subroutine normalize_old( vec , val , len )

implicit none

real,intent(in) :: vec( 2, 3 )
real,intent(out) :: val(2)
real,intent(out) :: len

real :: x1 , y1, x2 , y2, nrm1 , nrm2


y1 = vec(1,2) -vec( 1, 1 ) ;  y2 = vec( 1, 3 ) - vec( 1, 2 )
x1 = vec( 2, 1 ) - vec( 2, 2 );  x2 =  vec( 2, 2 ) - vec( 2, 3 )

nrm1 = sqrt( x1**2 + y1**2 ) ;   nrm2 = sqrt( x2**2 + y2**2 )

x1 = x1 / nrm1 ; y1 = y1 / nrm1  ;  x2 = x2 / nrm2;  y2 = y2 / nrm2

len =  sqrt( (x1+x2)**2 + (y1+y2)**2 )

val(1) = ( x1 + x2 ) / len
val(2) = ( y1 + y2 ) / len

len = ( nrm1 + nrm2 ) * 0.5


end subroutine

subroutine sortc (array, n)
!
! Purpose : To sort a character array into ascending order using a
!           selection sort.
!
implicit none
    
integer, intent(in) :: n
integer, intent(inout) :: array(n)
    
integer :: i, j, iptr, temp
    
! Sort the array
outer: do i = 1, n-1
    iptr = i 
    inner: do j = i+1, n 
        minval: if ( array(j) < array(iptr) ) then
            iptr = j 
        endif minval
    enddo inner
        
    swap: if ( i /= iptr ) then
        temp        = array(i)
        array(i)    = array(iptr)
        array(iptr) = temp
    endif swap
enddo outer
    
end subroutine sortc