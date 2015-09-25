subroutine remesh_routine_3D(unit_num, num_rn, remesh_num, current_step, ng, &
                             prop_remesh_flag, prop_num, prop_face_num, prop_point, prop_disp, prop_face, prop_patch)
    
use surface_domain
use remesh_domain
use materials
use file_info_house
implicit none

integer, intent(in) :: unit_num, num_rn, current_step, ng, prop_num, prop_face_num
integer, intent(in) :: prop_face(4,prop_face_num), prop_patch(prop_num)
real,    intent(in) :: prop_point(3,prop_num), prop_disp(3, prop_num)
integer, intent(inout) :: remesh_num, prop_remesh_flag(ng)

integer :: i, j, num, rn, num_patch, temp_arr, temp_num, cons_num, fix_cons_num, number_of_regions, mat_array_num
integer :: cur_ng, cur_num, new_mat_num, spring_num
integer :: input_val(10), cbou_info(2), cons(3)
real :: ref_temp, area, density, fac, theta
real :: damping(2), field(10), cbou_coord(6), cp(3)
character(len=10) :: str
character(len=120) :: line_str
integer, allocatable :: temp_val(:,:), mat_array(:)

integer :: tri_nn, tri_ne, tol_nn, tol_ne, prop_nn, prop_ne
integer, allocatable :: tol_elem(:,:), prop_cons(:,:), saving_point(:)
real,    allocatable :: tol_coord(:,:)

! file variables
integer :: main_unit_num, node_unit_num, ele_unit_num , cont_unit_num, inte_unit_num, status
character(len=30 ) :: re_main_file, re_node_file, re_ele_file, re_inte_file, re_mat_file, re_disp_file
integer :: unit_num1, unit_num2
character(len=5 ) :: keyword 
logical :: check

remesh_num = remesh_num + 1 
!num_patch = maxval(prop_patch)

re_main_file = './remesh/remesh000.main'
re_node_file = './remesh/remesh000.node'
re_ele_file= './remesh/remesh000.elem'
re_inte_file = './remesh/remesh000.inte'
re_mat_file = './remesh/mate000.dat'
re_disp_file = './remesh/disp000.dat' 

write(re_main_file(16:18) , '(I3.3)') remesh_num
write(re_node_file(16:18) , '(I3.3)') remesh_num
write(re_ele_file(16:18) , '(I3.3)') remesh_num
write(re_inte_file(16:18) , '(I3.3)') remesh_num
write(re_mat_file(14:16) , '(I3.3)') remesh_num
write(re_disp_file(14:16) , '(I3.3)') remesh_num

main_unit_num = 5999 
node_unit_num = 5000
ele_unit_num = 5001
!cont_unit_num = 5002
inte_unit_num = 5003
unit_num1 = 5004
unit_num2 = 5005

open(unit=main_unit_num, file=re_main_file , STATUS='replace', IOSTAT=status)
open(unit=node_unit_num, file=re_node_file , STATUS='replace', IOSTAT=status)
open(unit=ele_unit_num, file=re_ele_file , STATUS='replace', IOSTAT=status)
open(unit=inte_unit_num, file=re_inte_file , STATUS='replace', IOSTAT=status)
open(unit=unit_num1, file=re_mat_file , STATUS='replace', IOSTAT=status, form='unformatted')
open(unit=unit_num2, file=re_disp_file , STATUS='replace', IOSTAT=status, form='unformatted')

write(*,*) "Start remeshing!!!"
temp_arr = region(1)%num_nodes*5
allocate (tol_coord(3,temp_arr), tol_elem(8,temp_arr))
tol_nn = 0
tol_ne = 0
tol_coord = 0.0
tol_elem = 0

cur_ng = ng
mat_array_num = number_of_materials
if (material(number_of_materials)%model_num == 31) mat_array_num = mat_array_num - 1
        
do rn = 1, divided_mesh_region
    if (rn == 1) then
        write(*,*) " >> Region 1 remeshing"
        if (remesh_flag(rn) >= 1) then
            cur_num = 0
            do i = 1, ng
                prop_remesh_flag(i) = mesh_set(rn)%flag(i)
                write (*,*) i, 'sub domain mesh_flag:', mesh_set(rn)%flag(i)
                Hexa(i)%remesh_flag = mesh_set(rn)%flag(i)
                if (mesh_set(rn)%flag(i) /= 0) then
                    if (mesh_set(rn)%flag(i) == 3) then
                        mat_array_num = mat_array_num - 1
                        cur_ng = cur_ng - 1
                    else
                        cur_num = cur_num + 1
                        num = surface_set(i)%ridge_group_num
                        call remesh_3D_main(i, cur_num, temp_arr, surface_set(i)%tri_nn, surface_set(i)%tri_ne, &
                                            surface_set(i)%tri_coord, surface_set(i)%tri_elem, surface_set(i)%tri_np2, surface_set(i)%tri_p2, &
                                            surface_set(i)%ridge_group_num, surface_set(i)%ridge_num(1:num), surface_set(i)%ridge(:,1:num), surface_set(i)%ridge_conn(1:num), &
                                            tol_nn, tol_ne, tol_coord, tol_elem)
                    endif
                elseif (mesh_set(rn)%flag(i) == 0) then
                    cur_num = cur_num + 1
                    call meshless_3D(rn, i, cur_num, temp_arr, tol_nn, tol_ne, tol_coord, tol_elem)
                endif
            enddo
        endif
    elseif (rn == 2) then
        write(*,*) " >> Region 2 remeshing"
        allocate (mat_array(mat_array_num))
        mat_array = 0
        write (*,*) 'prop_remesh_flag:', prop_remesh_flag
        if (remesh_type_number == 1) then
            call remesh_ale_region1(ng, prop_remesh_flag, temp_arr, tol_nn, tol_ne, tol_coord, tol_elem, mat_array_num, mat_array)
        elseif (remesh_type_number == 2) then
            call remesh_ale_region2(ng, prop_remesh_flag, temp_arr, tol_nn, tol_ne, tol_coord, tol_elem, mat_array_num, mat_array)
        endif
    endif
enddo

rn = 1
! main input file
number_of_regions = 1
rewind(unit=unit_num) 
read(unit_num,*) keyword
if (keyword(1:1) == '!') then
    read(unit_num,*) str
else
    backspace(unit_num)
    read(unit_num,*) str
endif

temp_num = mat_array_num
if (material(number_of_materials)%model_num == 31) temp_num = temp_num + 1
write(main_unit_num, '(A,A,I3,A,I3,A,I3,A,I3)') str, ', ', number_of_regions, ', ', temp_num, ', ', num_load_sets, ', 0, ', inte_set%cont_num
input_val(1) = region(rn)%num_node_dof
input_val(2) = region(rn)%num_dim
input_val(3) = region(rn)%nodes_element
input_val(4) = region(rn)%num_quad_pts
if (rn == 1) then
    !input_val(5) = num_patch
    input_val(5) = region(rn)%num_groups
elseif (rn == 2) then
    input_val(5) = region(rn)%num_groups
endif

write (main_unit_num, *) '*doma, ', rn
write (main_unit_num, '(5x,I2,2x,I2,2x,I6,2x,I6,2x,I2,2x,I2,2x,I2)') input_val(1:2), tol_nn, tol_ne, input_val(3:5)

!write (*,*) size(load_set), num_load_sets
spring_num = 0
do
    read(unit_num, *, iostat=status) keyword
    if (keyword == '*doma') then
        backspace(unit_num)
        read(unit_num,*) keyword, num
        if (num == rn) then
            do
                read(unit_num, *, iostat=status) keyword
                if (keyword == '*spri') then
                    backspace(unit_num)
                    read(unit_num,*) str, spring_num
                    write (main_unit_num,*) str, spring_num
	            elseif (keyword == '*cbou') then
                    backspace(unit_num)
                    read(unit_num,*) keyword, num
                    write (main_unit_num,*) keyword, num
		            do i= 1, num
			            read (unit_num, *) cbou_coord, cbou_info
                        write (main_unit_num,*) cbou_coord, cbou_info
                    enddo
                elseif (keyword == '*cons') then
        		    backspace(unit_num)
                    allocate (prop_cons(4,tol_nn))
                    prop_cons = 0
                    call find_cons_number(unit_num, unit_num_inte, cons_num, fix_cons_num, region(rn)%num_nodes, region(rn)%node_coord, &
                                          tol_nn, tol_coord, prop_cons, remesh_dis*0.005)
                    write (main_unit_num,'(A,I6)') '*cons, ', cons_num
                    do i = 1, cons_num
                        write (main_unit_num,'(I7,3(A,I2))') prop_cons(1,i), ', ', prop_cons(2,i), ', ', prop_cons(3,i), ', ', prop_cons(4,i)
                    enddo
                    open (Unit=20, File='./output/solid/remesh/check_cons.plt', STATUS='replace', ACTION='write')
                    Write (20,'(A)') 'TITLE="Check case mesh"'
                    Write (20,*) 'VARIABLES="x", "y", "z", "cons_x", "cons_y", "cons_z"'
                    Write (20,'(A,I6,A,I6,A)') 'ZONE T = "Mesh", N = ', tol_nn, ', E = ', tol_ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
                    do i = 1, tol_nn
                        cons = 0
                        do j = 1, cons_num
                            if (prop_cons(1,j) == i) then
                                cons = prop_cons(2:4,j)
                                exit
                            endif
                        enddo
                        write (20,*) tol_coord(:,i), cons
                    enddo
                    do i = 1, tol_ne
                        write (20,*) tol_elem(:,i)
                    enddo
                    close(20)
                    deallocate(prop_cons)
                elseif (keyword == '*mpcn') then
                    backspace(unit_num)
                    read(unit_num,*) keyword, num
                    write (main_unit_num,*) keyword, num
                    do i = 1, num
                        read(unit_num,*) input_val(1:2), theta, input_val(3)
                        write (main_unit_num,*) input_val(1:2), theta, input_val(3)
                        do j = 1, 3
                            read(unit_num,*) field(1:3)
                            write (main_unit_num,*) field(1:3)
                        enddo
                    enddo                    
                elseif (keyword == '*doma' .OR. keyword == '*mate' .OR. keyword == '*tloa') then
                    backspace(unit_num)
                    exit
                elseif (status == -1) then
                    exit
                endif
            enddo
            exit
        endif
    elseif (status == -1) then
        exit
    endif
enddo           
write (main_unit_num, '(A,I3,A)') '*load, 1, ', load_set(1)%load_type, ', 1'
write (main_unit_num, *) ' 1, 1, 0.0'
write (main_unit_num, *)
    
! node & element file
! write information of nodes
write (*,*) 'Make: solid_3D.inp.node'
write(node_unit_num,*) '*node, 1'
do i = 1, tol_nn
    do j = 1, 3
        if (abs(tol_coord(j,i)) <= 10e-10) tol_coord(j,i) = 0.0
    enddo
    write(node_unit_num,'(I8,2X,3(ES20.10,2X))') i, tol_coord(:,i)
enddo

! write information of elements
write (*,*) 'Make: solid_3D.inp.elem'
rewind(unit_num_elem)
do i = 1, mat_array_num
    if (i == 1) then
        input_val(1) = 1
    else
        input_val(1) = mat_array(i-1)+1
    endif
    input_val(2) = mat_array(i)
    do 
        read(unit_num_elem,*) keyword
        if (keyword == '*elem') then
            backspace(unit_num_elem)
            read(unit_num_elem,*) str, input_val(3:5), area, density, fac, damping
            if (input_val(5) <= ng) then
                if (prop_remesh_flag(input_val(5)) /= 3) then
                    write(ele_unit_num,'(A,3(I3,1X),5(F12.7,1X))') str, input_val(3:5), area, density, fac, damping
                    do j = input_val(1), input_val(2)
                        write (ele_unit_num,'(I7,8(A,I7))') j,', ',tol_elem(1,j),', ',tol_elem(2,j),', ',tol_elem(3,j),', ',tol_elem(4,j), &
                                                            ', ',tol_elem(5,j),', ',tol_elem(6,j),', ',tol_elem(7,j),', ',tol_elem(8,j)
                    enddo
                    exit
                endif
            else
                write(ele_unit_num,'(A,3(I3,1X),5(F12.7,1X))') str, input_val(3:4), i, area, density, fac, damping
                do j = input_val(1), input_val(2)
                    write (ele_unit_num,'(I7,8(A,I7))') j,', ',tol_elem(1,j),', ',tol_elem(2,j),', ',tol_elem(3,j),', ',tol_elem(4,j), &
                                                        ', ',tol_elem(5,j),', ',tol_elem(6,j),', ',tol_elem(7,j),', ',tol_elem(8,j)
                enddo
                exit
            endif
        endif
    enddo
enddo
if (material(number_of_materials)%model_num == 31) then
    allocate(saving_point(region(1)%num_nodes))
    saving_point = 0
    rewind(unit_num_elem)
    do 
        read(unit_num_elem,*) keyword
        if (keyword == '*spri') then
            backspace(unit_num_elem)
            read(unit_num_elem,*) str, input_val(3:5), area
            write(ele_unit_num,'(A,3(I3,1X),F12.7)') str, input_val(3:4), mat_array_num+1, area
            do i = 1, spring_num
                read(unit_num_elem,*) num, input_val(1:2)
                do j = 1, 2
                    if (saving_point(input_val(j)) == 0) then
                        cp = region(1)%node_coord(:,input_val(j))
                        call find_match_node_num(cp, tol_nn, tol_coord(:,1:tol_nn), input_val(4+j), remesh_dis*0.001, 1)
                    else
                        input_val(4+j) = saving_point(input_val(j))
                    endif
                enddo
                write (ele_unit_num,'(I7,2(A,I7))') num,', ',input_val(5),', ',input_val(6)
            enddo
            exit
        endif
    enddo
    deallocate(saving_point)
endif    

! main file (common input variables)
do 
    read(unit_num,*,iostat=status) keyword
    if (status == -1) exit
    if ( keyword == '*tloa' ) then
        backspace(unit=unit_num)
        read(unit_num,*) keyword, temp_num, ref_temp
        write (main_unit_num,*) keyword, temp_num, ref_temp
        do i = 1, temp_num
            read(unit_num,*) field(1:2)
            write (main_unit_num,*) field(1:2)
        enddo
        write (main_unit_num,*)
    elseif ( keyword == '*mate' ) then
        backspace(unit=unit_num)
        write (main_unit_num,*)
        new_mat_num = 0
        do i = 1, number_of_materials
            read(unit_num,*) keyword, input_val(1:3)
            check = .true.
            if (input_val(1) <= ng) then
                if (prop_remesh_flag(input_val(1)) == 3) then
                    check = .false.
                endif
            endif
            if (check) then
                new_mat_num = new_mat_num + 1
                write (main_unit_num,*) keyword, new_mat_num, input_val(2:3)
            endif
            num = MOD(input_val(3), 10)
            do j = 1, input_val(3)/10
                read(unit_num,*) field(:)
                if (check) write (main_unit_num,'(10(ES16.9,1X))') field(:)
            enddo
            if ( num /= 0 ) then
                read(unit_num,*) field(1:num)
                if (check) write (main_unit_num,'(10(ES16.9,1X))') field(1:num)
            endif
        enddo
    elseif (keyword == '*end') then
        write (main_unit_num,*) '*end'
        write (main_unit_num,*) ' '
    elseif (keyword == '*reme') then
        write (main_unit_num,'(A,I3)') '*reme, ', remesh_num
    elseif (keyword == '*init') then
        write (main_unit_num,'(A,I7)') '*initstep, ', current_step
    elseif (keyword == '*stop') then
        write (main_unit_num,*) '*stop'
        exit
    else
        backspace(unit=unit_num)
        read(unit_num,'(A)') line_str
        write (main_unit_num,'(A)') line_str
    endif
enddo
rewind(unit=unit_num)
    
! write contact file
!if (num_rn >= 2) call check_contact_node(num_patch, prop_nn, prop_coord, unit_num_cont, cont_unit_num, 0)

! write interface file
write (*,*) 'Make: solid_3D.inp.inte'
allocate (temp_val(move_set(1)%csn(2,1)+1,2*cur_ng))
num = 0
do i = 1, ng
    if (prop_remesh_flag(i) /= 3) then
        num = num + 1
        temp_val(1:move_set(1)%csn(2,1),(num-1)*2+1) = move_set(i)%cs_nodes1(1,:,1)
        temp_val(1+move_set(1)%csn(2,1),(num-1)*2+1) = move_set(i)%cs_nodes1(1,1,2)
        temp_val(1:move_set(1)%csn(2,1),(num-1)*2+2) = move_set(i)%cs_nodes2(1,:,1)
        temp_val(1+move_set(1)%csn(2,1),(num-1)*2+2) = move_set(i)%cs_nodes2(1,1,2)
    endif
enddo
call write_inte_file_3D(unit_num_inte, inte_unit_num, move_set(1)%p2_csn, move_set(1)%csn(1:2,1), num_rn, cur_ng, temp_val, &
                        ng, prop_remesh_flag)

allocate (mesh_set(rn)%mat_array(mat_array_num))
mesh_set(rn)%mat_array = mat_array
call update_mesh_info(3, rn, tol_nn, tol_ne, tol_coord(:,1:tol_nn), tol_elem(:,1:tol_ne), 2)

deallocate (temp_val)
write (*,*) 'Free: remesh domain(Remesh_add)'
call free_geo_domain(3)

write (*,*) 'Close: remesh files'
close(node_unit_num)
close(ele_unit_num)
close(main_unit_num)
close(inte_unit_num)
deallocate (tol_coord, tol_elem)
write (*,*) 'End: subroutine remesh_routine_3D'

end subroutine remesh_routine_3D
                             
                             
                             
                             
subroutine meshless_3D(rn, ng, cur_ng, temp_arr, new_nn, new_ne, new_coord, new_elem)

use surface_domain
use remesh_domain
use materials
use file_info_house
implicit none

integer, intent(in) :: rn, ng, cur_ng, temp_arr
integer, intent(inout) :: new_nn, new_ne, new_elem(8,temp_arr)
real,    intent(inout) :: new_coord(3,temp_arr)

integer :: i, sp(4), lp(4), sub_nn, sub_ne, sub_bn, sub_bfn
integer, allocatable :: elem(:,:), bound_face(:,:)

sp = mesh_set(rn)%sub_num(ng,1:4)
lp = mesh_set(rn)%sub_num(ng+1,1:4) - 1
sub_nn = lp(1) - sp(1) + 1
sub_ne = lp(2) - sp(2) + 1
sub_bn = lp(3) - sp(3) + 1
sub_bfn = lp(4) - sp(4) + 1
new_coord(:,new_nn+1:new_nn+sub_nn) = mesh_set(rn)%node_coord(:,sp(1):lp(1))
allocate (elem(8,sub_ne))
elem = mesh_set(rn)%element(:,sp(2):lp(2))
if (ng /= 1) elem = elem - (mesh_set(rn)%sub_num(ng,1)-1)
new_elem(:,new_ne+1:new_ne+sub_ne) = elem + new_nn

new_nn = new_nn + sub_nn
new_ne = new_ne + sub_ne

! update mesh data to Hexa module
write (*,*) 'Don`t remesh: update mesh data to Hexa module', ng, cur_ng
write (*,*) 'nn, ne, bn, bfn:', sub_nn, sub_ne, sub_bn, sub_bfn
Hexa(cur_ng)%nn = sub_nn
Hexa(cur_ng)%ne = sub_ne
Hexa(cur_ng)%bn = sub_bn
Hexa(cur_ng)%bfn = sub_bfn
allocate (Hexa(cur_ng)%node(3,sub_nn), Hexa(cur_ng)%elem(8,sub_ne))
allocate (Hexa(cur_ng)%bound_node(sub_bn), Hexa(cur_ng)%bound_face(4,sub_bfn))
allocate (bound_face(4,sub_bfn))
Hexa(cur_ng)%node = mesh_set(rn)%node_coord(:,sp(1):lp(1))
Hexa(cur_ng)%elem = elem
Hexa(cur_ng)%bound_node = mesh_set(rn)%bound_node(sp(3):lp(3)) - (mesh_set(rn)%sub_num(ng,1)-1)
bound_face = mesh_set(rn)%bound_face(:,sp(4):lp(4))
do i = 1, sub_bfn
    bound_face(:,i) = mesh_set(rn)%bound_node(bound_face(:,i)) - (mesh_set(rn)%sub_num(ng,1)-1)
enddo
Hexa(cur_ng)%bound_face = bound_face
write (*,*) Hexa(cur_ng)%bound_node(1:6)

deallocate (elem, bound_face)

end subroutine meshless_3D
                             