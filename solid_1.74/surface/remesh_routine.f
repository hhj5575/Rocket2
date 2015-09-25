subroutine remesh_routine(unit_num, num_rn, unit_num_cont, remesh_num, current_step, ng, prop_remesh_flag, prop_num, prop_edge_num, prop_point, prop_disp, prop_edge, prop_patch)

use surface_domain
use materials

implicit none

integer, intent(in) :: unit_num, num_rn, unit_num_cont, current_step, ng, prop_num, prop_edge_num
integer, intent(in) :: prop_edge(2,prop_edge_num), prop_patch(prop_num)
real, intent(in) :: prop_point(2,prop_num), prop_disp(2, prop_num)
integer, intent(inout) :: remesh_num, prop_remesh_flag(ng)

integer :: i, j, rn, nn, ne, st, num, nb, num_patch, number, num_groups
integer :: tol_nn, tol_ne, ori_nn, ori_ne
integer :: max_node_dof, max_dim, nodes_elem, num_quad_pts, temp_num, cons_num
integer :: prop_nn, prop_ne, case_nn, mat_ele_num
integer :: cons(3), temp_var(10), cbou_info(2), input_val(10)
real :: ref_temp
real :: damping(2), field(10), cbou_coord(4)

real, allocatable :: coord(:,:), prop_coord(:,:)
real, allocatable :: case_coord(:,:)
integer , allocatable :: element(:,:), prop_ele(:,:)
integer :: type_num, mat_num, count, del_flag
real :: area, density, fac

integer :: main_unit_num, node_unit_num, ele_unit_num , cont_unit_num, inte_unit_num, status
integer :: unit_num1, unit_num2 
character(len=30 ) :: main_file, node_file, ele_file , cont_file, inte_file, mat_file,disp_file
character(len=5 ) :: keyword 
character(len=120) :: line_str

remesh_num = remesh_num + 1 
num_patch = maxval(prop_patch)

main_file = './remesh/remesh000.main'
node_file = './remesh/remesh000.node'
ele_file= './remesh/remesh000.elem'
cont_file = './remesh/remesh000.cont'
inte_file = './remesh/remesh000.inte'
mat_file = './remesh/mate000.dat'
disp_file = './remesh/disp000.dat' 

write( main_file(16:18) , '(I3.3)') remesh_num
write( node_file(16:18) , '(I3.3)') remesh_num
write( ele_file(16:18) , '(I3.3)') remesh_num
write( cont_file(16:18) , '(I3.3)') remesh_num
write( inte_file(16:18) , '(I3.3)') remesh_num
write( mat_file(14:16) , '(I3.3)') remesh_num
write( disp_file(14:16) , '(I3.3)') remesh_num

main_unit_num = 5999 
node_unit_num = 5000
ele_unit_num = 5001
cont_unit_num = 5002
inte_unit_num = 5003
unit_num1 = 5004
unit_num2 = 5005

open(unit=main_unit_num, file=main_file , STATUS='replace', IOSTAT=status)
open(unit=node_unit_num, file=node_file , STATUS='replace', IOSTAT=status)
open(unit=ele_unit_num, file=ele_file , STATUS='replace', IOSTAT=status)
open(unit=cont_unit_num, file=cont_file , STATUS='replace', IOSTAT=status)
open(unit=inte_unit_num, file=inte_file , STATUS='replace', IOSTAT=status)
open(unit=unit_num1, file=mat_file , STATUS='replace', IOSTAT=status, form='unformatted')
open(unit=unit_num2, file=disp_file , STATUS='replace', IOSTAT=status, form='unformatted')

nn = 0;  ne = 0; nb = 0;  del_flag = 0
write(*,*) "Start remeshing!!!"

! main file
write(main_unit_num, '(A,I3,A,I3,A,I3,A,I3)') 'Remesh, ', num_rn, ', ', number_of_materials, ', ', num_load_sets, ', 0, ', inte_set%cont_num

do rn = 1, num_rn
    !node, ele input file
    ori_nn = region(rn)%num_nodes
    ori_ne = region(rn)%num_elements
    allocate( coord( 2,ori_nn*3), element(4,ori_ne*3) )
    coord(:,1:ori_nn) = region(rn)%node_coord
    element(:,1:ori_ne) = region(rn)%element(4:7,:)
    
    if ( rn == 1 .and. remesh_flag(rn) >= 1 ) then
        tol_nn = 0; tol_ne = 0
        coord = 0.0;  element = 0
        write (*,*) 'num_patch:', num_patch
        do i = 1, num_patch
            write (*,'(A,I2,A,I2,A,I3)') 'mesh_set(', rn, ')%flag(', i, ') :', mesh_set(rn)%flag(i)
            prop_remesh_flag(num_patch-i+1) = mesh_set(rn)%flag(i)
            !if ( mesh_set(rn)%flag(i) == 1 .OR. mesh_set(rn)%flag(i) == 2 ) then
            if (mesh_set(rn)%flag(i) /= 0) then
                call mesh_gen_main(current_step, rn, num_patch, i, tol_nn, tol_ne, ori_nn, coord, ori_ne, element, prop_remesh_flag(num_patch-i+1))
                !call adjust_PROPEL_POINT(i, tol_nn, st, surface_set(i)%adjn, coord, prop_num, prop_edge_num, prop_point, prop_disp, prop_edge, prop_patch)
            !elseif ( mesh_set(rn)%flag(i) == 3 ) then
            !    call mesh_gen_part(current_step, rn, i, tol_nn, tol_ne, ori_nn, coord, ori_ne, element)
            elseif ( mesh_set(rn)%flag(i) == 0 ) then
                call meshless(rn, i, tol_nn, tol_ne, ori_nn, coord, ori_ne, element)
            endif
        enddo
        nn = tol_nn;  ne = tol_ne
	else
		!call write_input(ori_nn, ori_ne, coord, element, node_unit_num, ele_unit_num, 1)
		nn = ori_nn;  ne = ori_ne
    endif
    write (*,*) rn, 'region of remesing end'

    ! main input file
    max_node_dof = region(rn)%num_node_dof
    max_dim = region(rn)%num_dim
    nodes_elem = region(rn)%nodes_element
    num_quad_pts = region(rn)%num_quad_pts
    if (rn == 1) then
        num_groups = num_patch
    elseif (rn == 2) then
        num_groups = region(rn)%num_groups
    endif

    write (main_unit_num, *) '*doma, ', rn
    write (main_unit_num, '(5x,I2,2x,I2,2x,I6,2x,I6,2x,I2,2x,I2,2x,I2)') max_node_dof, max_dim, nn, ne, nodes_elem, num_quad_pts, num_groups
    if ( num_groups > 1 .AND. rn == 1 ) then
        do i = 1, num_groups
            !if ( mesh_set(rn)%flag(i) /= 4 ) then
                write (main_unit_num, *) i, surface_set(i)%num_nodes, surface_set(i)%num_elements
            !endif
        enddo
    endif

    !write (*,*) size(load_set), num_load_sets
    if ( rn == 1 ) then
        rewind(unit=unit_num)
        do
            read(unit_num, *, iostat=status) keyword
            if (keyword == '*doma') then
                backspace(unit_num)
                read(unit_num,*) keyword, number
                if (number == rn) then
                    do
                        read(unit_num, *, iostat=status) keyword
	                    if (keyword == '*cbou') then
                            backspace(unit_num)
                            read(unit_num,*) keyword, number
                            write (main_unit_num,*) keyword, number
		                    do i= 1, number
			                    read (unit_num, *) cbou_coord, cbou_info
                                write (main_unit_num,*) cbou_coord, cbou_info
                            enddo
                            exit
                        elseif (keyword == '*doma' .OR. keyword == '*mate' .OR. status == -1) then
                            exit
                        endif
                    enddo
                    exit
                endif
            elseif (status == -1) then
                exit
            endif
        enddo            
        write (main_unit_num, *) '*load, 1, ', load_set(1)%load_type, ', 1'
        write (main_unit_num, *) ' 1, 1, 0.0'
        write (main_unit_num, *)
    elseif ( rn == 2 ) then
        rewind(unit=unit_num)
        do 
            read(unit_num, *, iostat=status) keyword
            if (keyword == '*doma') then
                backspace(unit_num)
                read(unit_num,*) keyword, number
                if (number == rn) then
                    do
                        read(unit_num, *, iostat=status) keyword
	                    if (keyword == '*cons') then
        		            backspace(unit_num)
        		            read(unit_num,*) keyword, cons_num
                            write (main_unit_num,*) keyword, cons_num
                            do i = 1, cons_num
                                read(unit_num,*) cons(:) 
                                write (main_unit_num,*) cons(1), ', ', cons(2), ', ', cons(3)
                            enddo
                            exit        		    
        		        elseif (keyword == '*doma' .OR. keyword == '*mate' .OR. keyword == '*tloa' .OR. status == -1) then
        		            exit
        		        endif
                    enddo
                    exit
                endif
            elseif (keyword == '*mate' .OR. keyword == '*tloa' .OR. status == -1) then
                exit
            endif
        enddo
        write (main_unit_num, *) '*load, 1, ', load_set(2)%load_type, ', 2'
        write (main_unit_num, *) ' 1, 1, 0.0'
    endif        
    
    if ( rn == 1 ) then
        ! write information of nodes
        write (node_unit_num,'(A,I2)') '*node, ', rn
        prop_nn = nn;  prop_ne = ne
        allocate ( prop_coord(2,nn), prop_ele(4,ne) )
        prop_coord = coord(:,1:nn)
        prop_ele = element(:,1:ne)
        
        ! remesh contact
        if (num_rn >= 2) call remesh_contact(num_groups, prop_nn, prop_coord)
        do i = 1, prop_nn
            write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i,', ',prop_coord(1,i),', ',prop_coord(2,i)
        enddo
        
        ! write information of elements
        keyword = '*elem'
        type_num = region(rn)%element(2,1)
        mat_num = region(rn)%element(3,1)
        area = region(rn)%element_area(1)
        density = region(rn)%element_density(1)
        fac = region(rn)%consis_mass_ratio(1)
        damping(:) = region(rn)%damping_factor(1,:)
        write(ele_unit_num,'(A,2x,I2,2x,I3,2x,I2,2x,F8.4,2x,F14.6,2x,E10.3,2x,F5.2,2x,F5.2)') keyword, rn, type_num, mat_num, area, density, fac, damping(:)
        do i = 1, ne
            write (ele_unit_num,'(I5,A,I5,A,I5,A,I5,A,I5)') i,', ',element(1,i),', ',element(2,i),', ',element(3,i),', ',element(4,i)
        enddo
    elseif ( rn == 2 ) then
        case_nn = nn;
        allocate ( case_coord(2,nn) )
        case_coord = coord(:,1:nn)
        write (node_unit_num,'(A,I2)') '*node, ', rn
        do i = 1, case_nn
            write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i,', ',case_coord(1,i),', ',case_coord(2,i)
        enddo
        
        ! write information of elements
        keyword = '*elem'
        count = 0
        do mat_num = 2, number_of_materials
            mat_ele_num = material(mat_num)%num_col_state / num_quad_pts
            type_num = region(rn)%element(2,count+1)
            area = region(rn)%element_area(count+1)
            density = region(rn)%element_density(count+1)
            fac = region(rn)%consis_mass_ratio(count+1)
            damping(:) = region(rn)%damping_factor(count+1,:)
            write(ele_unit_num,'(A,2x,I2,2x,I3,2x,I2,2x,F8.4,2x,F14.6,2x,E10.3,2x,F5.2,2x,F5.2)') keyword, rn, type_num, mat_num, area, density, fac, damping(:)
            do i = 1, mat_ele_num
                count = count + 1
                write (ele_unit_num,'(I5,A,I5,A,I5,A,I5,A,I5)') count,', ',element(1,count),', ',element(2,count),', ',element(3,count),', ',element(4,count)
            enddo
        enddo
        if ( count /= ne ) STOP "Remesh_routine : ne is not summation of mat_ele_number"
    endif
    
    call update_mesh_info(2, rn, nn, ne, coord(:,1:nn), element(:,1:ne), 2)
    deallocate( coord, element )
enddo

write (main_unit_num,*) ' '
rewind(unit=unit_num)
do 
    read(unit_num,*) keyword
    if ( keyword == '*tloa' .or. keyword == '*mate') then
        backspace(unit=unit_num)
        exit
    endif
enddo

do 
    read(unit_num,*) keyword
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
        do i = 1, number_of_materials
            read(unit_num,*) keyword, input_val(1:3)
            write (main_unit_num,*) keyword, input_val(1:3)
            num = MOD(input_val(3), 10)
            do j = 1, input_val(3)/10
                read(unit_num,*) field(:)
                write (main_unit_num,'(10(ES16.9,1X))') field(:)
            enddo
            if ( num /= 0 ) then
                read(unit_num,*) field(1:num)
                write (main_unit_num,'(10(ES16.9,1X))') field(1:num)
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

! write cont_file
if (num_rn >= 2) call check_contact_node(num_patch, prop_nn, prop_coord, unit_num_cont, cont_unit_num, 0)

!write inte_file
if (num_rn >= 2) then
    call write_inte_file(del_flag, num_rn, prop_nn, prop_coord, case_nn, case_coord, inte_unit_num)
else
    call write_inte_file(del_flag, num_rn, prop_nn, prop_coord, prop_nn, prop_coord, inte_unit_num)
endif
    
close(node_unit_num)
close(ele_unit_num)
close(main_unit_num)
close(cont_unit_num)
close(inte_unit_num)

contains

subroutine adjust_PROPEL_POINT(ng, tol_nn, st, adjn, coord, prop_num, prop_edge_num, prop_point, prop_disp, prop_edge, prop_patch)

implicit none

integer, intent(in) :: ng, tol_nn, st, adjn, prop_num, prop_edge_num
integer, intent(in) :: prop_patch(prop_num), prop_edge(2,prop_edge_num)
real, intent(in) :: prop_point(2,prop_num), prop_disp(2,prop_num)
real, intent(inout) :: coord(2,tol_nn)

integer :: i, j, this(2), temp_patch(prop_num)
real :: cp(2), np(2), dis, temp_coord(2,prop_num)
logical :: cont_check

num = 0
do i = prop_num, 1, -1
    num = num + 1
    temp_coord(1,num) = prop_point(2,i) ! - prop_disp(2,i)
    temp_coord(2,num) = prop_point(1,i) ! - prop_disp(1,i)
    temp_patch(num) = prop_patch(i)
enddo

!open (Unit=20, File='check_PROPEL_surface.plt', STATUS='replace', ACTION='write')
!Write (20,'(A)') 'TITLE="Check random edge"'
!Write (20,*) 'VARIABLES="x", "y"'
!Write (20,'(A,I5,A,I5, A)') 'ZONE N =', prop_num , ', E =', prop_edge_num, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
!do i = 1, prop_num
!    write (20,*) temp_coord(:,i)
!enddo
!do i = 1, prop_edge_num
!    write (20,*) prop_num-prop_edge(1,i)+1, prop_num-prop_edge(2,i)+1
!enddo
!close(20)

!call set_skin(1, ng, sp, lp)
!lp = sp + lp - 1
!memory_j = sp
!do i = st+1, st+adjn
!    cp = coord(:,i)
!    do j = memory_j, lp
!        this = (/ j, j+1 /)
!        if (j == lp) this(2) = sp
!        call calc_len(temp_coord(:,this(1)), temp_coord(:,this(2)), dis)
!        call near_point(temp_coord(:,this(1)), temp_coord(:,this(2)), cp, np, dis*0.01, cont_check)
!        if (cont_check) then
!            coord(:,i) = np
!            memory_j = j
!            exit
!        endif
!    enddo 
!enddo

do i = st+1, st+adjn
    cp = coord(:,i)
    do j = 1, prop_edge_num
        if (temp_patch(j) == ng) then
            this = (/ prop_num-prop_edge(1,j)+1, prop_num-prop_edge(2,j)+1 /)
            call calc_len(temp_coord(:,this(1)), temp_coord(:,this(2)), dis)
            call near_point(temp_coord(:,this(1)), temp_coord(:,this(2)), cp, np, dis*0.01, cont_check)
            if (cont_check) then
                coord(:,i) = np
                exit
            endif
        endif
    enddo
enddo

end subroutine adjust_PROPEL_POINT

end subroutine remesh_routine
