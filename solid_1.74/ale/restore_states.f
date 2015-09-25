subroutine restore_states(Dim, flag, num_rn, remesh_num, step_num ) 

use physical_domain
use time_domain
use materials
use system_domain
use file_info_house

implicit none

integer, intent(in) :: Dim, flag, num_rn, remesh_num, step_num
integer :: unit_num1, unit_num2 , num_var , col, i, j, rn, nn, status, mat_num, a, b
real :: vec(Dim,3)
real, allocatable  :: value(:,:)
integer :: position(Dim), mat_model_num
character(len=40) :: mat_file,disp_file

if (flag == 1) then
    mat_file = main_input_file
    disp_file = main_input_file
    a = len_trim(main_input_file)
    b = a + 5
    a = a + 1
    mat_file(a:b) = '.mate'
    disp_file(a:b) = '.disp'    
elseif (flag == 2) then
    mat_file = './remesh/mate000.dat'
    disp_file = './remesh/disp000.dat'
    write( mat_file(14:16) , '(I3.3)') remesh_num
    write( disp_file(14:16) , '(I3.3)') remesh_num
endif

unit_num1 = 6000
unit_num2 = 7000
open(unit=unit_num2, file=disp_file , STATUS='old',IOSTAT=status,form='unformatted')

do rn = 1, num_rn
    nn = region(rn)%num_nodes
    do i= 1, nn
        do j = 1, Dim
	        read(unit_num2) vec(j,1)
        enddo
        do j = 1, Dim
	        read(unit_num2) vec(j,2)
        enddo
        do j = 1, Dim
	        read(unit_num2) vec(j,3)
        enddo
	    do j = 1, Dim
	        position(j) = region(rn)%node(3+j,i)
        enddo
        if (Dim == 2) then
	        call update_motion(rn, Dim, position, vec)
        elseif (Dim == 3) then
            call update_motion_3D(rn, Dim, position, vec)
        endif
    enddo
enddo
close( unit_num2 )

open(unit=unit_num1, file=mat_file , STATUS='old',IOSTAT=status, form='unformatted')
do mat_num = 1, number_of_materials
    num_var=material(mat_num)%num_row_state
    col=material(mat_num)%num_col_state
    mat_model_num = material(mat_num)%model_num
    if (mat_model_num /= 31) then
        allocate( value( num_var, col) )
        do i = 1, col
            do j = 1, num_var
                read (unit_num1) value(j,i)
            enddo
        enddo
        call update_states(mat_num, num_var, col, value )
        deallocate ( value )
    endif
enddo
close( unit_num1 ) 

end subroutine restore_states