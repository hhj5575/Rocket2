subroutine ale_control(Dim, num_rn, delta_t, current_step, remesh, n_remesh, ng, PROPEL_REMESH_FLAG, PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_EDGE, PROPEL_PATCH)                                                                          
                                                                                                
use ale_domain
use file_info_house
use system_info_house

implicit none

integer, intent(in) :: Dim, num_rn, current_step, ng, PROPEL_POINT_NUM, PROPEL_EDGE_NUM
integer, intent(in) :: PROPEL_EDGE(2,PROPEL_EDGE_NUM), PROPEL_PATCH(PROPEL_POINT_NUM)
real, intent(in) :: delta_t, PROPEL_POINT(2,PROPEL_POINT_NUM), PROPEL_DISP(2,PROPEL_POINT_NUM)
integer, intent(inout) :: remesh, n_remesh, PROPEL_REMESH_FLAG(ng)

integer :: rn, i, id
! remesh_flag : 0 - don`t remesh
!               1 - seperation remesh
!               2 - the other remesh(evaluation, length between boundary nodes)
 
! remesh or separation
if (Dim == 2) then
    if (n_remesh == 2) then
        write (*,*) " >> Remesh routine"
        !remesh_energy_file = './output/solid/remesh/energy_chk000.dat'
        !write (remesh_energy_file(33:35), '(I3.3)') remesh+1
        !open (Unit=unit_num_remesh_energy, File=remesh_energy_file, STATUS='replace', ACTION='write')
        !write (unit_num_remesh_energy,'(A,I)') 'Remesh:', remesh+1
        !write (unit_num_remesh_energy,'(A,ES15.7)') 'Before:', sta_energy
        !close(unit_num_remesh_energy)
        remesh_energy_flag = .TRUE.
        !call set_fsi_inte(num_rn, 2, 1)
        if (remesh_flag(1) == 1) then
            call separation_routine(unit_number, num_rn, unit_num_cont, remesh, current_step)
        else
            call remesh_routine(unit_number, num_rn, unit_num_cont, remesh, current_step, ng, &
                                PROPEL_REMESH_FLAG, PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_EDGE, PROPEL_PATCH)
        endif
    
        do rn = 1, num_rn
            call remesh_interpolation2(Dim, rn, delta_t, solver_flag)
        enddo
        call free_fsi_domain
        deallocate (loadf, dispf)    ! 'system_info_house'
    elseif (n_remesh == 0) then
        if (move_flag >= 1) then
            write (*,'(A)') " >> FVM motion: Burning"
            call fvm_routine_bak(Dim, 1, solver_flag)
            !call fvm_routine(Dim, 1, delta_t)
    
            call update_corner_info(1)
            if (solver_flag) then
                if (num_rn >= 2) then
                    !call post_build_interface(unit_num_cont)
                    call post_build_interface(current_step)
                    id = 1  ! temporary setting
                    if (is_coarse_problem) then
                        if (coarse_schur) then
                            call update_coarse_disp
                        elseif (.NOT. coarse_schur) then
                        !if (.NOT. coarse_schur) then
                            !call set_global_sparse(id)
                        endif
                    endif
                endif
            endif
        endif
        if (abla_flag == 1) then
            write (*,'(A)') " >> FVM motion: Ablation"
            call fvm_routine_case(Dim, delta_t)
        endif
        write (*,*) 
        call free_mesh_info(1)
        do i = 1, ng
            if (allocated(surface_set(i)%coord)) deallocate(surface_set(i)%coord)
        enddo
        call set_fsi_inte(num_rn, 1, 1, current_step)
        if (fsi_inte_set%fsi_ni /= fsi_inte_set%old_fsi_ni) then
            n_remesh = 1
        endif
    endif
    !write (*,'(A,2ES25.15)') '8(ale_control)  :', region(1)%node_coord(2,8), region(1)%node_coord(1,8)
elseif (Dim == 3) then
    if (n_remesh == 2) then
        write (*,*) " >> Remesh routine"
        call remesh_routine_3D(unit_number, num_rn, remesh, current_step, ng, &
                               PROPEL_REMESH_FLAG, PROPEL_POINT_NUM, PROPEL_EDGE_NUM, PROPEL_POINT, PROPEL_DISP, PROPEL_EDGE, PROPEL_PATCH)
        
        !do rn = 1, num_rn
            !call remesh_interpolation2(Dim, rn, delta_t, solver_flag)
        !enddo
        write (*,*) " >> Remesh interpolation"
        call remesh_interpolation2(Dim, 1, delta_t, solver_flag)
        write (*,*) " >> Free domain varialbes"
        call free_fsi_domain
        deallocate (loadf, dispf)    ! 'system_info_house'
    elseif (n_remesh <= 1) then
        if (move_flag >= 1) then
            write (*,'(A)') " >> FVM motion: Burning"
            call fvm_routine_bak(Dim, 1, solver_flag)
        endif
    endif
endif

end subroutine ale_control

subroutine make_backup_file(Dim, current_step, remesh)

use ale_domain
use remesh_domain
use file_info_house
use system_info_house
implicit none

integer, intent(in) :: Dim, current_step, remesh

integer :: rn, i, j, ng, nn, ne, bn, sp, lp, num, mat_num, row, col, status, seq(2)
integer :: mat_ele_num, ele_nodes, flag
real :: area, density, fac, ratio(2), damping(2), zero, p(3)
logical :: existence

integer, allocatable :: b_nodes(:)
real, allocatable :: coord(:,:)

! backup file data
integer :: node_unit_num,  cont_unit_num, unit_num1, unit_num2, unit_num3
integer :: main_unit_num, elem_unit_num, inte_unit_num, type_num
character(len=40) :: bak_main_file, bak_elem_file, bak_inte_file
character(len=40) :: bak_node_file, bak_cont_file, bak_mate_file, bak_disp_file
character(len=200) :: char
character(len=5) :: keyword

bak_main_file = './restart/solid0000000.inp'
bak_elem_file = './restart/solid0000000.inp.elem'
bak_inte_file = './restart/solid0000000.inp.inte'
bak_node_file = './restart/solid0000000.inp.node'
bak_mate_file = './restart/solid0000000.inp.mate'
bak_disp_file = './restart/solid0000000.inp.disp'
write (bak_main_file(16:22), '(I7.7)') current_step
write (bak_node_file(16:22), '(I7.7)') current_step
write (bak_elem_file(16:22), '(I7.7)') current_step
write (bak_inte_file(16:22), '(I7.7)') current_step
write (bak_mate_file(16:22), '(I7.7)') current_step
write (bak_disp_file(16:22), '(I7.7)') current_step

main_unit_num = 5999 
node_unit_num = 5001
elem_unit_num = 5002
inte_unit_num = 5003
unit_num1 = 5005
unit_num2 = 5006

open(unit=main_unit_num, file=bak_main_file , STATUS='replace', IOSTAT=status)
open(unit=node_unit_num, file=bak_node_file , STATUS='replace', IOSTAT=status)
open(unit=elem_unit_num, file=bak_elem_file , STATUS='replace', IOSTAT=status)
open(unit=inte_unit_num, file=bak_inte_file , STATUS='replace', IOSTAT=status)
open(unit=unit_num1, file=bak_mate_file , STATUS='replace', IOSTAT=status, form='unformatted')
open(unit=unit_num2, file=bak_disp_file , STATUS='replace', IOSTAT=status,form='unformatted')

if (Dim == 2) then
    bak_cont_file = './restart/solid0000000.inp.cont'
    write (bak_cont_file(16:22), '(I7.7)') current_step
    cont_unit_num = 5004
    open(unit=cont_unit_num, file=bak_cont_file , STATUS='replace', IOSTAT=status)
endif
ele_nodes = (Dim-1)*4

do rn = 1, num_of_regions
    nn = region(rn)%num_nodes
    ne = region(rn)%num_elements
    
    ! write node file
    write (node_unit_num,'(A,I2)') '*node, ', rn
    do i = 1, nn
        if (Dim == 2) then
            write (node_unit_num,'(I6,2(A,F16.10))') i,', ',region(rn)%node_coord(1,i),', ', region(rn)%node_coord(2,i)
        elseif (Dim == 3) then
            write (node_unit_num,'(I6,3(A,F16.10))') i,', ',region(rn)%node_coord(1,i),', ', region(rn)%node_coord(2,i),', ', region(rn)%node_coord(3,i)
        endif
    enddo

    ! write element file
    keyword = '*elem'
    num = 0
    mat_ele_num = 0
    do i = 1, ne
        mat_num = region(rn)%element(3,i)
        if (mat_ele_num /= mat_num) then
            type_num = region(rn)%element(2,i)
            area = region(rn)%element_area(i)
            density = region(rn)%element_density(i)
            fac = region(rn)%consis_mass_ratio(i)
            damping(:) = region(rn)%damping_factor(i,:)
            write(elem_unit_num,'(A,2x,I2,2x,I3,2x,I2,2x,F8.4,2x,F14.6,2x,E10.3,2x,F5.2,2x,F5.2)') keyword, rn, type_num, mat_num, area, density, fac, damping(:)
            mat_ele_num = mat_num
        endif
        write (elem_unit_num,'(I6,2X,8(I7,2X))') i, region(rn)%element(4:3+ele_nodes,i)
    enddo
    if (region(rn)%num_springs /= 0) then
        type_num = region(rn)%spring(2,1)
        mat_num = region(rn)%spring(3,1)
        area = region(rn)%spring_area(1)
        write(elem_unit_num,'(A,2x,I2,2x,I3,2x,I2,2x,F8.4)') '*spring ' , rn, type_num, mat_num, area
        do i = 1, region(rn)%num_springs
            write (elem_unit_num,'(I6,2X,2(I7,2X))') i, region(rn)%spring(4:5,i)
        enddo
    endif
    !else
    !    count = 0
    !    do mat_num = 2, number_of_materials
    !        mat_ele_num = material(mat_num)%num_col_state / region(rn)%num_quad_pts
    !        type_num = region(rn)%element(2,count+1)
    !        area = region(rn)%element_area(count+1)
    !        density = region(rn)%element_density(count+1)
    !        fac = region(rn)%consis_mass_ratio(count+1)
    !        damping(:) = region(rn)%damping_factor(count+1,:)
    !        write(ele_unit_num,'(A,2x,I2,2x,I3,2x,I2,2x,F8.4,2x,F14.6,2x,E10.3,2x,F5.2,2x,F5.2)') keyword, rn, type_num, mat_num, area, density, fac, damping(:)
    !        do i = 1, mat_ele_num
    !            count = count + 1
    !            write (ele_unit_num,'(I5,2X,4(I6,2X))') count, region(rn)%element(4:7,count)
    !        enddo
    !    enddo
    !    if ( count /= ne ) STOP "make_backup_file : ne is not summation of mat_ele_number"
    !endif
    
    ! write disp file
    zero = 0.0
    if (solver_flag) then
        if (Dim == 2) then
            do i=1,nn
                write(unit_num2) history(rn)%data(region(rn)%node(4,i), 1)
	            write(unit_num2) history(rn)%data(region(rn)%node(5,i), 1)
    	
	            write(unit_num2) history(rn)%data(region(rn)%node(4,i), 2)
	            write(unit_num2) history(rn)%data(region(rn)%node(5,i), 2)
    	
	            write(unit_num2) history(rn)%data(region(rn)%node(4,i), 3)
	            write(unit_num2) history(rn)%data(region(rn)%node(5,i), 3)
            enddo
        elseif (Dim == 3) then
            do i = 1, nn
                write(unit_num2) sdomain(rn)%u_cartesian(region(rn)%node(4,i))
                write(unit_num2) sdomain(rn)%u_cartesian(region(rn)%node(5,i))
                write(unit_num2) sdomain(rn)%u_cartesian(region(rn)%node(6,i))
            
	            write(unit_num2) zero
	            write(unit_num2) zero
                write(unit_num2) zero
    	
	            write(unit_num2) zero
	            write(unit_num2) zero
                write(unit_num2) zero
            enddo
        endif
    endif
enddo

! write material file
if (solver_flag) then
    do mat_num = 1, number_of_materials
        row = material(mat_num)%num_row_state
        col = material(mat_num)%num_col_state
        do i = 1, col
            do j = 1, row
                write (unit_num1) material(mat_num)%state1(j,i)
            enddo
        enddo
    enddo
endif
     
! write main file
seq = 0
do i = 1, rest_set%line_num
    keyword = rest_set%line_str(i)(1:5)
    if (keyword == '*cons') then
        write (main_unit_num,'(A)') rest_set%line_str(i)
        seq(1) = seq(1) + 1
        do j = 1, rest_set%cons_num(seq(1))
            if (seq(1) == 1) then
                write (main_unit_num,*) rest_set%cons1(:,j)
            elseif (seq(1) == 2) then
                write (main_unit_num,*) rest_set%cons2(:,j)
            elseif (seq(1) == 3) then
                write (main_unit_num,*) rest_set%cons3(:,j)
            endif
        enddo
    elseif (keyword == '*load') then
        write (main_unit_num,'(A)') rest_set%line_str(i)
        seq(2) = seq(2) + 1
        do j = 1, rest_set%load_num(1,seq(2))
            if (seq(2) == 1) then
                write (main_unit_num,*) rest_set%load_info1(:,j), rest_set%load1(j)
            elseif (seq(2) == 2) then
                write (main_unit_num,*) rest_set%load_info2(:,j), rest_set%load2(j)
            elseif (seq(2) == 3) then
                write (main_unit_num,*) rest_set%load_info3(:,j), rest_set%load3(j)
            endif
        enddo
    elseif (keyword == '*reme') then
        write (main_unit_num,'(A,I3)') '*reme, ', remesh
    elseif (keyword == '*init') then
        write (main_unit_num,'(A,I7)') '*initstep, ', current_step
    else
        write (main_unit_num,'(A)') rest_set%line_str(i)
    endif
enddo

! write inte file
if (Dim == 3) then
    rewind (unit_num_inte)
    do 
        read(unit_num_inte,*,iostat=status) keyword
        backspace(unit_num_inte)
        if (keyword == '*remp') then
            read(unit_num_inte,*) keyword, ng
            write (inte_unit_num,'(A,1X,I3)') keyword, ng
        
            nn = mesh_set(1)%nn
            sp = mesh_set(1)%sub_num(ng,3)
            lp = mesh_set(1)%sub_num(ng+1,3)
            bn = lp-sp
            allocate (coord(Dim,nn), b_nodes(bn))
            coord = mesh_set(1)%node_coord
            b_nodes = mesh_set(1)%bound_node(sp:lp-1)
        
            read(unit_num_inte,*,iostat=status) keyword, num
            write (inte_unit_num,'(A,1X,I5)') keyword, num
            do i = 1, num
                read(unit_num_inte,*) p
                write (inte_unit_num,'(3(F15.8,2X))') coord(:,b_nodes(remesh_add(ng)%ori_p2(i)))
            enddo
            read(unit_num_inte,*,iostat=status) keyword, num
            write (inte_unit_num,'(A,1X,I5)') keyword, num
            do i = 1, num
                read(unit_num_inte,*) flag
                backspace(unit_num_inte)
                if (flag == 1) then
                    read(unit_num_inte,*) flag, p, seq(1:2)
                    write (inte_unit_num,'(I2,2X,3(F15.8,2X),2(I2,2X))') flag, coord(:,b_nodes(remesh_add(ng)%p2(1,i))), seq(1:2)
                elseif (flag == 2) then
                    read(unit_num_inte,'(A)', iostat=status) char
                    write (inte_unit_num,'(A)') char
                endif
            enddo
            read(unit_num_inte,*,iostat=status) keyword, num
            write (inte_unit_num,'(A,1X,I5)') keyword, num
            do i = 1, num
                read(unit_num_inte,'(A)', iostat=status) char
                write (inte_unit_num,'(A)') char
            enddo
            read(unit_num_inte,*,iostat=status) keyword, num
            write (inte_unit_num,'(A,1X,I5)') keyword, num
            do i = 1, num
                read(unit_num_inte,'(A)', iostat=status) char
                write (inte_unit_num,'(A)') char
            enddo
            deallocate (coord, b_nodes)
        else
            read(unit_num_inte,'(A)', iostat=status) char
            if (status /= 0) exit
            write (inte_unit_num,'(A)') char
        endif
    enddo
elseif (Dim == 2) then
    rewind (unit_num_inte)
    do
        read(unit_num_inte,'(A)', iostat=status) char
        if (status /= 0) exit
        write (inte_unit_num,'(A)') char
    enddo
endif

! write contact file
if (move_flag > 1 .and. Dim == 2) then
    inquire (file=bak_cont_file, EXIST=existence)
    if (existence == .false.) then
        open(unit=cont_unit_num, file=bak_cont_file , STATUS='replace', IOSTAT=status)
        rewind(unit=unit_num_cont)
        do 
            read(unit_num_cont, '(A)', iostat=status) char
            if (status /= 0) exit
            write(cont_unit_num, '(A)') char
        enddo
        close(cont_unit_num)
    endif
endif

close(main_unit_num)
close(node_unit_num)
close(elem_unit_num)
close(inte_unit_num)
close(unit_num1)
close(unit_num2)

end subroutine make_backup_file