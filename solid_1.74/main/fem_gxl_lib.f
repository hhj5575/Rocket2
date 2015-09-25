subroutine find_star(unit_num, state_flag)

implicit none

integer, intent(in) :: unit_num
integer, intent(out) :: state_flag
integer :: stat
character(len=1) :: key


state_flag = 0

do
	read(unit_num, *, iostat=stat) key

	if (stat == -1) then
		state_flag = -1
		exit   ! exit from do-while: meet EOF
	elseif (key == '*') then
		backspace(unit_num)
		state_flag = 1
		EXIT   ! exit from do-while
	endif
enddo

end subroutine find_star






subroutine find_starword(unit_num, keyword, state_flag)

implicit none

integer, intent(in) :: unit_num
character(len=4), intent(in) :: keyword
integer, intent(out) :: state_flag
integer :: stat
character(len=5) :: key


state_flag = 0

do
	read(unit_num, *, iostat=stat) key

	if (stat == -1) then
		state_flag = -1
		exit   ! exit from do-while: meet EOF
	elseif (key(1:1) == '*') then
    if (key(2:5) == keyword) then
      backspace(unit_num)
      state_flag = 1
      EXIT   ! exit from do-while
    endif
	endif
enddo

end subroutine find_starword






subroutine find_loc(info, dim, val)

implicit none

integer,intent(in) :: dim
integer,intent(in) :: info(dim)
integer,intent(inout) :: val
integer :: i

do i = 1, dim
	if ( info(i) .EQ. val ) then
		val = i 
		exit
	endif 
enddo

end subroutine find_loc






subroutine write_matrix(a, dim1, dim2)

integer :: dim1,dim2, i
real :: a(dim1,dim2)

open(unit=201, file='matrix.txt', status='REPLACE', action='WRITE')

do i=1,dim1
write(201,*) a(i,:)
enddo

end subroutine write_matrix








subroutine add_char(name, new_name, flag)


character(35),intent(in) :: name
integer,intent(in) :: flag
character(35), intent(out) :: new_name

character(5) :: flag1, flag2, flag3, flag4
integer :: a,b


new_name = name
flag1 = '.node'
flag2 = '.elem'
flag3 = '.cont'
flag4 = '.inte'

a = len_trim(name)
b = a + 5
a = a + 1

if (flag == 1) then
  new_name(a:b) = flag1
elseif (flag ==2) then
  new_name(a:b) = flag2
elseif (flag ==3) then
  new_name(a:b) = flag3
elseif (flag ==4) then
  new_name(a:b) = flag4
else
  STOP 'add_char: flag error!'
endif

end subroutine add_char







subroutine logic_to_num(tf_vec, m,  result, n, count)

integer, intent(in) :: m, n
logical, intent(in) :: tf_vec(m)
integer, intent(out) :: result(n) , count 
integer :: i 

count = 0 
result = 0
do i = 1, m 
	if ( tf_vec(i) ) then
		count = count + 1 
		result(count) = i 
	endif 
enddo

end subroutine logic_to_num







subroutine matrix_output_debug(rn, ele_num, iter, m, n, matrix)

use file_info_house !------------------------------------------>>>>
use element_specification

implicit none

integer, intent(in) :: rn, ele_num, iter, m, n
real, intent(in) :: matrix(m,n)
integer :: i, j, k
integer, parameter :: ndim = 2, num_elem_nodes = 4, ndof = 2
real :: xl(ndim, num_elem_nodes), ul(ndof, 5*num_elem_nodes)
integer, parameter :: max_n = 100 ! max number of columns of the matrix

if (n > max_n) STOP 'matrix_output_debug: max_n = 100'

k = iter + 1    ! to avoid 'zero' index in k_debug

write(out3_unit_number,'(A,I2,A)') 'k_debug(:,:,', k, ') = ['
do i = 1, m
  write(out3_unit_number,100) (matrix(i,j), j=1,n)														
end do
write(out3_unit_number,'(A/)') '];'


call fetch_xl(rn, ele_num, num_elem_nodes, xl, ndim, .TRUE.)

write(out3_unit_number,'(A,I2,A)') 'xl(:,:,', k, ') = ['
do i = 1, ndim
  write(out3_unit_number,100) (xl(i,j), j=1,num_elem_nodes)														
end do
write(out3_unit_number,'(A/)') '];'

call fetch_ul(rn, ele_num, num_elem_nodes, ul, ndof, .TRUE.)

write(out3_unit_number,'(A,I2,A)') 'ul(:,:,', k, ') = ['
do i = 1, ndof
  write(out3_unit_number,100) (ul(i,j), j=1,num_elem_nodes)														
end do
write(out3_unit_number,'(A/)') '];'


100 format (100ES15.7)

end subroutine matrix_output_debug
