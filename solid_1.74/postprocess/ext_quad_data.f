subroutine ext_quad_data(Dim, unit_num, num_var, tol_ne, temp_quad_data)

implicit none

integer :: Dim, tol_ne, num_var, unit_num
real :: temp_quad_data(num_var,tol_ne*(Dim-1)*4)

integer :: i, roof, count, status
integer :: mat, ele, quad
character(len=5) :: keyword
!======================================================================
rewind(unit_num)
count = 0
temp_quad_data = 0.0
do 
	read (unit_num,*, iostat=status ) keyword
    if ( status == -1 ) then
        exit
    elseif (keyword == '*numb') then
        read (unit_num,*, iostat=status) keyword
        do 
            read (unit_num,*,iostat=status) mat
            backspace(unit_num)
            if (status == 59 .or. status == -1 .or. keyword == '*numb') then
                exit
            else
                count = count + 1
                read (unit_num,*,iostat=status) mat, ele, quad, temp_quad_data(:,count)
            endif
        enddo
    endif
enddo

end subroutine ext_quad_data