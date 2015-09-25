subroutine separation_routine(unit_num, num_rn, unit_num_cont, remesh_num, current_step)

use surface_domain
use materials

implicit none

integer, intent(in) :: unit_num, num_rn, unit_num_cont, current_step
integer, intent(inout) :: remesh_num

integer :: i, j, k, num, rn, rg, sp, lp, count, adjn, nn, ne, sn, current_domain, number, temp_bn
integer :: tol_nn, tol_ne, temp_bs, temp_nn, temp_ne, temp_num, node_num, cons_num, del_flag
integer :: type_num, mat_num, max_node_dof, max_dim, nodes_elem, num_quad_pts
integer :: load_type_number, ndim, coarse_count, the_dof, direct, sep_node, nb, ori_nn
integer :: roof, remain, case_adjn, new_adjn, ord_adjn, new_nn, new_ne, num_groups
integer :: prop_nn, prop_ne, mat_ele_num
integer :: this(2), inte_load_node(2), cons(3), val(5), bn(2), temp_var(10), cbou_info(2)
real :: tol, load_value, dis, temp_d, fac, area, density, tol_area
real :: cont_coord(2,(inte_set%cont_num+1)*2), load_coord(2)
real :: target_node(2,2), inte_load_value(2), cp(2), vec(2), temp_sep_coord(2,2), damping(2), field(10), cbou_coord(4)
logical :: check, frist_check, non_check, near_check, temp_check(2), sep_check(2)

integer, allocatable ::temp_adj(:), tadj(:), temp_ele(:,:), element(:,:), adj(:), new_adj(:), ord_adj(:)
integer, allocatable :: temp_inte_node(:,:), case_adj(:), case_bound(:), new_ele(:,:)
integer, allocatable :: prop_ele(:,:)
real, allocatable :: temp_b(:,:), temp_coord(:,:), coord(:,:), temp_inte_value(:), new_coord(:,:)
real, allocatable :: prop_coord(:,:)

integer :: main_unit_num, node_unit_num, ele_unit_num , cont_unit_num, inte_unit_num, status
integer :: unit_num1 , unit_num2
character(len=30 ) :: main_file, node_file, ele_file , cont_file, inte_file, mat_file, disp_file, fn
character(len=5 ) :: keyword
!-------------------------------------------------------------------------------------
remesh_num = remesh_num + 1 

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
open(unit=unit_num1, file=mat_file , STATUS='replace', IOSTAT=status,form='unformatted')
open(unit=unit_num2, file=disp_file , STATUS='replace', IOSTAT=status,form='unformatted')

tol = mesh_set(num_rn)%ave_d(1) * 0.4
!call set_skin(2, 1, sp, case_adjn)
!allocate (case_adj(case_adjn+1))
!case_adj(1:case_adjn) = region(2)%skin_nodes(sp:sp+case_adjn-1)
!case_adj(case_adjn+1) = region(2)%skin_nodes(sp)
sep_check = .true.
del_flag = 0
write(*,*) "Start separation!!!"

! main file
write(main_unit_num, '(A,I3,A,I3,A,I3,A,I3)') 'Remesh, ', num_rn, ', ', number_of_materials, ', ', num_load_sets, ', 0, ', inte_set%cont_num
do rn = 1, num_rn
    !node, ele input file
    nn = region(rn)%num_nodes
    ne = region(rn)%num_elements
    allocate( coord( 2, nn), element(4,ne) )
    allocate( new_coord(2,nn*4), new_ele(4,ne*4) )
    if ( rn == 1 ) then
        coord = surface_set(1)%coord
    else
        coord = region(rn)%node_coord
    endif
    element = region(rn)%element(4:7,:)
    
    if ( rn == 1 ) then
        deallocate (surface_set)
        allocate (surface_set(2))
        temp_d = mesh_set(rn)%ave_d(1)
	    coord(:,sep_num(3)) = sep_point(:)
        temp_check = .FALSE.

        call set_skin(1, 1, sp, adjn)
        allocate ( tadj(adjn), temp_adj(adjn) )
        temp_adj = region(1)%skin_nodes(sp:sp+adjn-1)
        if ( sep_seq(3) > sep_seq(1) ) then !움직이는 선들의 번호가 접촉면의 번호보다 클 때~
            
            ! 첫 번째 형상
            tadj = 0;  num = 0
            do i = 1, sep_seq(1)
                num = num + 1
                tadj(num) = temp_adj(i)
            enddo
            call calc_len(coord(:,sep_num(1)), coord(:,sep_num(3)), dis)
            if ( dis < 0.35*temp_d ) then
                num = num - 1
                !temp_check(1) = .TRUE.
            endif
            vec = coord(:,tadj(num)) - coord(:,sep_num(3))
            temp_sep_coord(:,1) = coord(:,sep_num(3)) + vec(:)/10.0
            do i = sep_seq(3), adjn
                num = num + 1
                tadj(num) = temp_adj(i)
            enddo
            allocate ( ord_adj(num) ) 
            ord_adjn = num
            ord_adj(1:num) = tadj(1:num)
            
            ! 두 번째 형상
            tadj = 0;  num = 0
            call calc_len(coord(:,sep_num(2)), coord(:,sep_num(3)), dis)
            num = num + 1
            tadj(num) = sep_num(3)
            if ( dis >= 0.35*temp_d ) then
                num = num + 1
                tadj(num) = sep_num(2)
                !temp_check(2) = .TRUE.
            endif
            do i = sep_seq(2)+1, sep_seq(3)-1
                num = num + 1
                tadj(num) = temp_adj(i)
            enddo
            vec = coord(:,tadj(2)) - coord(:,sep_num(3))
            temp_sep_coord(:,2) = coord(:,sep_num(3)) + vec(:)/10
            allocate ( new_adj(num) )
            new_adjn = num
            new_adj = tadj(1:num)

        else !움직이는 선들의 번호가 접촉면의 번호보다 작을 때~
            ! 첫 번째 형상
            tadj = 0;  num = 0
            do i = 1, sep_seq(3)
                num = num + 1
                tadj(num) = temp_adj(i)
            enddo
            call calc_len(coord(:,sep_num(2)), coord(:,sep_num(3)), dis)
            if ( dis < 0.35*temp_d ) then
                num = num - 1
                temp_check(1) = .TRUE.
            endif
            vec = coord(:,sep_num(2)) - coord(:,sep_num(3))
            temp_sep_coord(:,1) = coord(:,sep_num(3)) + vec(:)/10
            do i = sep_seq(2), adjn
                num = num + 1
                tadj(num) = temp_adj(i)
            enddo
            
            allocate ( ord_adj(num) ) 
            ord_adjn = num
            ord_adj(1:num) = tadj(1:num)
            
            ! 두 번째 형상
            tadj = 0;  num = 0
            call calc_len(coord(:,sep_num(1)), coord(:,sep_num(3)), dis)
            num = num + 1
            tadj(num) = sep_num(3)
            if ( dis >= 0.35*temp_d ) then
                temp_check(2) = .TRUE.
            endif
            do i = sep_seq(3)+1, sep_seq(1)-1
                num = num + 1
                tadj(num) = i
            enddo
            if (temp_check(2)) then
                num = num + 1
                tadj(num) = sep_num(1)
            endif
            vec = coord(:,sep_num(1)) - coord(:,sep_num(3))
            temp_sep_coord(:,2) = coord(:,sep_num(3)) + vec(:)/10

            allocate ( new_adj(num) )
            new_adjn = num
            new_adj = tadj(1:num)
        endif

        ! 요소망 생성
        tol_nn = 0; tol_ne = 0
        new_coord = 0.0;  new_ele = 0
        do i = 1, 2
            if ( i == 1 ) then
                temp_bs = ord_adjn
            else
                temp_bs = new_adjn
            endif
            ! temp_bs의 홀수, 짝수 검사
            if ( 2 * (temp_bs / 2) /= temp_bs ) then
                temp_bs = temp_bs + 1
                allocate ( temp_b(2,temp_bs) )
                temp_b = 0.
                do j = 1, temp_bs-1
                    if ( i == 1 ) then
                        temp_b(:,j) = coord(:,ord_adj(j))
                        !call calc_len(temp_b(:,j), coord(:,sep_num(3)), dis)
                        !if ( dis < 10e-6 ) temp_b(:,j) = temp_sep_coord(:,1)
                        if (ord_adj(j) == sep_num(3)) then
                            temp_b(:,j) = temp_sep_coord(:,1)
                        endif
                    else
                        temp_b(:,j) = coord(:,new_adj(j))
                        !call calc_len(temp_b(:,j), coord(:,sep_num(3)), dis)
                        !if ( dis < 10e-6 ) temp_b(:,j) = temp_sep_coord(:,2)
                        if (new_adj(j) == sep_num(3)) then
                            temp_b(:,j) = temp_sep_coord(:,2)
                        endif
                    endif
                enddo
                call make_even_number(temp_bs, temp_b)
            else
                allocate ( temp_b(2,temp_bs) )
                temp_b = 0.
                do j = 1, temp_bs
                    if ( i == 1 ) then
                        temp_b(:,j) = coord(:,ord_adj(j))
                        call calc_len(temp_b(:,j), coord(:,sep_num(3)), dis)
                        if ( dis < 10e-6 ) temp_b(:,j) = temp_sep_coord(:,1)
                    else
                        temp_b(:,j) = coord(:,new_adj(j))
                        call calc_len(temp_b(:,j), coord(:,sep_num(3)), dis)
                        if ( dis < 10e-6 ) temp_b(:,j) = temp_sep_coord(:,2)
                    endif
                enddo
            endif
            
            ! 분리된 형상의 넓이 검사
            ! mesh_set(2)%ave_d(1)^2*5.0 이하이면 mesh를 생성하지 않음(skip)
            tol_area = mesh_set(num_rn)%ave_d(1)**2.0 * sqrt(5.0)
            vec = 0.
            do j = 1, temp_bs
                if (j == temp_bs) then
                    vec(1) = vec(1) + temp_b(1,j)*temp_b(2,1)
                    vec(2) = vec(2) + temp_b(1,1)*temp_b(2,j)
                else
                    vec(1) = vec(1) + temp_b(1,j)*temp_b(2,j+1)
                    vec(2) = vec(2) + temp_b(1,j+1)*temp_b(2,j)
                endif
            enddo
            area = 0.5*abs(vec(1)-vec(2))

            if (area < tol_area) then
                write (*,*) i, ' object is smaller than mesh_set(2)%ave_d(1)^2*sqrt(5.0)'
                write (*,*) 'object area, tolerance area:', area, tol_area
                sep_check(i) = .false.
                del_flag = 1
            else
                surface_set(i)%num_nodes = tol_nn + 1
                surface_set(i)%num_elements = tol_ne + 1
                call sep_mesh_gen(i, temp_bs, temp_b, temp_d, tol_nn, tol_ne, nn, ne, new_coord, new_ele, new_adjn)
            endif
            deallocate (temp_b)
        enddo
        
        ! coord, element 재설정
        deallocate ( coord, element )
        allocate ( coord(2,tol_nn), element(4,tol_ne) ) 
        coord = new_coord(:,1:tol_nn)
        element = new_ele(:,1:tol_ne)
        nn = tol_nn;  ne = tol_ne
        deallocate ( new_coord, new_ele)
        deallocate ( tadj, temp_adj, new_adj )
    else
		!call write_input(nn, ne, coord, element, node_unit_num, ele_unit_num, 2)
	endif
    write (*,*) rn, ': end remeshing'
    
    ! main input file
    max_node_dof = region(rn)%num_node_dof
    max_dim = region(rn)%num_dim
    nodes_elem = region(rn)%nodes_element
    num_quad_pts = region(rn)%num_quad_pts
    num_groups = region(rn)%num_groups
    if ( rn == 1 ) then
        if (sep_check(1) .AND. sep_check(2)) num_groups = num_groups + 1
    endif
    write (main_unit_num, *) '*doma, ', rn
    write (main_unit_num, '(5x,I2,2x,I2,2x,I6,2x,I6,2x,I2,2x,I2,2x,I2)') max_node_dof, max_dim, nn, ne, nodes_elem, num_quad_pts, num_groups
    if ( num_groups > 1 ) then
        do i = 1, num_groups
            write (main_unit_num, *) i, surface_set(i)%num_nodes, surface_set(i)%num_elements
        enddo
    endif
    
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
        write (main_unit_num, *) '  1, 1, 0.0' 
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
        		        elseif (keyword == '*doma' .OR. keyword == '*mate' .OR. status == -1) then
        		            exit
        		        endif
                    enddo
                    exit
                endif
            elseif (keyword == '*mate' .OR. status == -1) then
                exit
            endif
        enddo
        write (main_unit_num, *) '*load, 1, ', load_set(2)%load_type, ', 2'
        write (main_unit_num, *) ' 1, 1, 0.0'
        write (main_unit_num, *)
    endif
   
    if ( rn == 1 ) then
        allocate ( prop_coord(2,nn), prop_ele(4,ne) ) 
        prop_nn = nn;  prop_ne = ne
        prop_coord = coord
        prop_ele = element
        
        write (node_unit_num,'(A,I2)') '*node, ', rn
        do i = 1, nn
            write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i,', ',coord(1,i),', ',coord(2,i)
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
        write (node_unit_num,'(A,I2)') '*node, ', rn
        do i = 1, nn
            write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i,', ',coord(1,i),', ',coord(2,i)
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
do 
    read(unit_num,*) keyword
    if ( keyword == '*mate' ) then
        backspace(unit=unit_num)
        do i = 1, number_of_materials
            read(unit_num,*) keyword, temp_var(1:3)
            write (main_unit_num,*) keyword, temp_var(1:3)
            num = MOD(temp_var(3), 10)
            do j = 1, temp_var(3)/10
                read(unit_num,*) field(:)
                write (main_unit_num,*) field(:)
            enddo
            if ( num /= 0 ) then
                read(unit_num,*) field(1:num)
                write (main_unit_num,*) field(1:num)
            endif
        enddo
        exit
    endif
enddo
rewind(unit=unit_num)
write (main_unit_num,*) '*end'
write (main_unit_num,*) '*step, ', current_step
write (main_unit_num,*) '*reme, ', remesh_num

! write cont_file
if (num_rn >= 2) then
    if (sep_check(1) == .false. .OR. sep_check(2) == .false.) then
        if (sep_check(1) == .false.) then
            allocate (surface_set(1)%adj(surface_set(2)%adjn))
            surface_set(1)%adjn = surface_set(2)%adjn
            surface_set(1)%adj = surface_set(2)%adj
        endif
        !call check_contact_node(prop_nn, prop_coord, unit_num_cont, cont_unit_num, 0)
    else
        !call check_contact_node(prop_nn, prop_coord, unit_num_cont, cont_unit_num, 1)
    endif
endif

!write inte_file
call write_inte_file(del_flag, num_rn, prop_nn, prop_coord, prop_nn, prop_coord, inte_unit_num)

close( node_unit_num )
close( ele_unit_num )
close( main_unit_num )
close( cont_unit_num )

end subroutine