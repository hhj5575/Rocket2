!=======================================================================

subroutine fvm_routine_bak(Dim, rn, solver_flag)

use ale_domain

implicit none 

integer, intent(in) :: DIm, rn
logical, intent(in) :: solver_flag

integer :: num_dim, ele_node, num_quad , st, ed
integer :: element((Dim-1)*4)

integer :: num_face, size, val_size, ce, num
integer :: i, j, k , num_mat, col , row , error, temp_ne, num_group
integer :: position(Dim), this((Dim-1)*4)

real :: vec(Dim,3), in_coord(Dim), temp_dis, temp_cp(Dim)
real :: temp_at(Dim,(Dim-1)*4), temp_vt(Dim,(Dim-1)*4), temp_ut(Dim,(Dim-1)*4)
real :: dist, min_val, local_val, tol
logical :: in_element

integer :: quad_nn, quad_ne
! mesh data
integer :: nn, ne, temp_quad_ne, mest_set_nn, mat_ne, mat_num, mat_model_num
integer, allocatable :: ele(:,:), temp_quad_ele(:), conn(:), mat_elem_num(:)
real, allocatable :: coord(:,:), ori_coord(:,:), at(:,:), vt(:,:), ut(:,:), temp_coord(:,:)
!quad data
integer, allocatable :: quad_ele(:,:), new_quad_ele(:,:)
real, allocatable :: quad_node(:,:), quad_dis(:), quad_cp(:,:), quad_vol(:)
real, allocatable :: new_quad_node(:,:), new_quad_dis(:), new_quad_cp(:,:)
real,allocatable :: new_value(:,:), temp_value(:)
real :: local_node(Dim,9)

!nn = mesh_set(rn)%nn
!ne = mesh_set(rn)%ne
nn = region(rn)%num_nodes
ne = region(rn)%num_elements
ele_node = (Dim-1)*4
allocate(coord(Dim,nn), ele(ele_node,ne))
coord = 0.0
ele = region(rn)%element(4:3+ele_node,:)
if (Dim == 2) then
    coord = mesh_set(rn)%coord
elseif (Dim == 3) then
    do i = 1, divided_mesh_region
        do j = 1, mesh_set(i)%nn
            coord(:,mesh_set(i)%surf2ori(j)) = mesh_set(i)%node_coord(:,j)
        enddo
    enddo
endif
num_dim = Dim

if (solver_flag) then
    allocate(ori_coord(num_dim,nn), at(num_dim,nn), vt(num_dim,nn), ut(num_dim,nn))
    ele_node = region(rn)%nodes_element
    num_quad = region(rn)%num_quad_pts
    mat_num = 1
    num_mat = pdomain_mat_db%num_materials

    row = material(1)%num_row_state
    col = material(1)%num_col_state
    if (num_dim == 2) then
        size = 9
        do i=1,nn
            do j = 1, num_dim
	            ut(j,i) =  history(rn)%data(region(rn)%node(3+j,i), 1)
	            vt(j,i) =  history(rn)%data(region(rn)%node(3+j,i), 2)
                at(j,i) =  history(rn)%data(region(rn)%node(3+j,i), 3)
            enddo
        enddo
    elseif (num_dim == 3) then
        size = 27
        do i = 1,nn
            ut(:,i) = sdomain(rn)%u_cartesian(region(rn)%node(4:6,i))
            vt(:,i) = 0.0
            at(:,i) = 0.0
        enddo
    endif
    
    !ori_coord = region(rn)%node_coord(1:2, :) + ut
    ori_coord = region(rn)%node_coord(1:num_dim, :)
    
    num_group = region(rn)%num_groups
    if (num_dim == 2) then
        tol = sum(mesh_set(rn)%ave_d(1:num_group))/float(num_group)
    elseif (num_dim == 3) then
        tol = 10e-2                 ! temporory
    endif
    
    allocate (temp_coord(Dim,ele_node))
    st = 1;  ed = ne
    do  i = 1, nn
        if (num_dim == 2) then
	        call calc_len(ori_coord(:,i), coord(:,i), dist)
        elseif (num_dim == 3) then
            call calc_length_3D(ori_coord(:,i), coord(:,i), dist)
        endif
	    if (dist > tol*10e-4) then
		    vec = 0.0
		    in_coord = coord(:,i)
		    min_val = 10e8
		    do j= st, ed
			    do k = 1, ele_node
				    temp_coord(:,k) = ori_coord(:,ele(k,j))
                enddo
                if (num_dim == 2) then
			        call inelement(temp_coord, in_coord, in_element, local_val, 1.0e-03)
                elseif (num_dim == 3) then
                    call inout_check_cubic(in_coord, temp_coord, in_element, local_val, 1.0e-03)
                endif
			    if (in_element) then
				    do k = 1, ele_node
					    temp_ut(:,k) = ut(:,ele(k,j))
					    temp_vt(:,k) = vt(:,ele(k,j))
					    temp_at(:,k) = at(:,ele(k,j))
                    enddo
                    
                    if (num_dim == 2) then
				        call goto_NC(temp_coord, in_coord, temp_ut, temp_vt, temp_at, vec, error )
                    elseif (num_dim == 3) then
                        call goto_NC_3D(temp_coord, in_coord, temp_ut, temp_vt, temp_at, vec, error )
                    endif
                    
                    if (error == 1) then
                        !open (Unit=31, File='goto_NC_err.dat', STATUS='replace', ACTION='write')
                        !do k = 1, ele_node
                        !    write (31,*) temp_coord(:,k)
                        !enddo
                        !write (31,*)
                        !write (31,*) in_coord
                        !close(31)
				    elseif ( error /= 1 ) then
                        do k = 1, num_dim
				            position(k)= region(rn)%node(3+k,i)
                        enddo
                        if (num_dim == 2) then
				            call update_motion(rn, num_dim, position, vec)
                        elseif (num_dim == 3) then
                            call update_motion_3D(rn, num_dim, position, vec)
                        endif
					    !coord(:,i) = coord(:,i) - vec(:,1)
					    exit
				    endif
			    else
				    if ( min_val > local_val ) then
					    min_val = local_val
					    temp_ne = j
				    endif
			    endif
			    if ( j == ed ) then
                    !write (*,'(I6,1X,4(F11.6,1X),I6,1X,A)') i, coord(:,i), min_val, temp_ne, ':: Cannot Found!'
				    do k = 1, ele_node
					    temp_coord(:,k) = ori_coord(:,ele(k,temp_ne))
					    temp_ut(:,k) = ut(:,ele(k,temp_ne))
					    temp_vt(:,k) = vt(:,ele(k,temp_ne))
					    temp_at(:,k) = at(:,ele(k,temp_ne))
                    enddo
                    if (num_dim == 2) then
				        call outelement(temp_coord, in_coord, vec, temp_ut, temp_vt, temp_at, error)
                    elseif (num_dim == 3) then
                        !write (*,*) i, ' node: outelement_3D. ', temp_ne
                        call outelement_3D(temp_coord, in_coord, vec, temp_ut, temp_vt, temp_at, error)
                    endif
                    
                    if (error == 1) stop 'Error: Not convergence at subroutine goto_NC'
                    do k = 1, num_dim
				        position(k)= region(rn)%node(3+k,i)
                    enddo
                    if (num_dim == 2) then
				        call update_motion(rn, num_dim, position, vec)
                    elseif (num_dim == 3) then
                        call update_motion_3D(rn, num_dim, position, vec)
                    endif
				    !coord(:,i) = coord(:,i) - vec(:,1)
				    exit
                endif
            enddo
	    else
		    !coord(:,i) = coord(:,i) - ut(:,i)
	    endif
    enddo
    call update_coord(rn, num_dim, nn, coord )

    ! state variables update
    if (num_dim == 2) then
        if ( row .NE. 0 ) then 
            allocate(temp_value(row)) 
            allocate( new_value(row, col)) 
            temp_value = 0.0
            new_value = 0.0 
            quad_nn = ne*size;  quad_ne = ne*ele_node
            allocate (quad_node(Dim,quad_nn), quad_ele(ele_node,quad_ne))
            call elem_to_quadelem(num_dim, ele_node, nn, ne, ori_coord, ele, quad_nn, quad_ne, quad_node, quad_ele)
            !if (num_dim == 2) then
            !    call write_tecplot(quad_nn, quad_node, quad_ne, quad_ele, 1, 1, 1) 
            !elseif (num_dim == 3) then
            !    call write_tecplot_3D(quad_nn, quad_node, quad_ne, quad_ele, 1, 1, 1) 
            !endif
       
            do i = 1, ne
                element = ele(:,i)
                call set_local_info( nn, num_dim, coord, ele_node, element, size, local_node )
                temp_quad_ne = 1 + ale_set(rn)%ele_adj(i+1) - ale_set(rn)%ele_adj(i)
                allocate (temp_quad_ele(temp_quad_ne)) 
                temp_quad_ele(1) = i
                temp_quad_ele(2:temp_quad_ne) = ale_set(rn)%ele_adjn(ale_set(rn)%ele_adj(i):ale_set(rn)%ele_adj(i+1)-1)
                do j = 1, 4
                    this = (/ j, j+4, 9, j+3 /)
                    if (j == 1) this(4) = 8
                    temp_coord = local_node(:,this)
                    !write (*,*) 'quad number : ', (i-1)*4+j
                    call area_interpolation(rn, mat_num, temp_quad_ne, temp_quad_ele, quad_nn, quad_node, quad_ne, quad_ele, temp_coord, row, temp_value, 10)
                    new_value(:,(i-1)*4+j) = temp_value
                    !write (*,*) 
                enddo
                deallocate (temp_quad_ele) 
            enddo
        endif
    elseif (num_dim == 3) then
        quad_nn = ne*size;  quad_ne = ne*ele_node
        allocate (quad_node(Dim,quad_nn), quad_ele(ele_node,quad_ne))
        allocate (new_quad_node(Dim,quad_nn), new_quad_ele(ele_node,quad_ne))
        allocate (quad_dis(quad_ne), new_quad_dis(quad_ne))
        allocate (quad_cp(3,quad_ne), new_quad_cp(3,quad_ne))
        allocate (quad_vol(quad_ne))
        call elem_to_quadelem(num_dim, ele_node, nn, ne, ori_coord, ele, quad_nn, quad_ne, quad_node, quad_ele)
        call elem_to_quadelem(num_dim, ele_node, nn, ne, coord, ele, quad_nn, quad_ne, new_quad_node, new_quad_ele)
        !call calc_quad_cp_dis(quad_nn, quad_node, quad_ne, quad_ele, quad_cp, quad_dis)
        !call calc_quad_cp_dis(quad_nn, new_quad_node, quad_ne, new_quad_ele, new_quad_cp, new_quad_dis)
        call calc_quad_volume(quad_nn, quad_node, quad_ne, quad_ele, quad_vol)
        do mat_num = 1, number_of_materials
            mat_ne = 0
            do i = 1, ne
                if (region(rn)%element(3,i) == mat_num) mat_ne = mat_ne +1
            enddo
            allocate (mat_elem_num(mat_ne))
            mat_ne = 0
            do i = 1, ne
                if (region(rn)%element(3,i) == mat_num) then
                    mat_ne = mat_ne + 1
                    mat_elem_num(mat_ne) = i
                endif
            enddo
            row = material(mat_num)%num_row_state
            col = material(mat_num)%num_col_state
            mat_model_num = material(mat_num)%model_num
            allocate(temp_value(row)) 
            allocate(new_value(row,col)) 
            temp_value = 0.0
            new_value = 0.0
            if ( row .NE. 0 .and. mat_model_num /= 31) then 
                do i = 1, mat_ne
                    ce = mat_elem_num(i)
                    !write (*,*) '=========', ce, 'element ================'
                    temp_quad_ne = ale_set(rn)%ele_adjn(ce+1) - ale_set(rn)%ele_adjn(ce)
                    allocate (temp_quad_ele(temp_quad_ne))
                    !temp_quad_ele = ale_set(rn)%ele_adj(ale_set(rn)%ele_adjn(ce):ale_set(rn)%ele_adjn(ce+1)-1)
                    num = 0
                    do j = ale_set(rn)%ele_adjn(ce), ale_set(rn)%ele_adjn(ce+1)-1
                        if (mat_num == region(rn)%element(3,ale_set(rn)%ele_adj(j))) then
                            num = num + 1 
                            temp_quad_ele(num) = ale_set(rn)%ele_adj(j)
                        endif
                    enddo
                    temp_quad_ne = num
                    do j = 1, 8
                        temp_coord = new_quad_node(:,quad_ele(:,(ce-1)*8+j))
                        temp_cp = new_quad_cp(:,(ce-1)*8+j)
                        temp_dis = new_quad_dis((ce-1)*8+j)
                        !call vol_interpolation_3D(rn, ce, j, mat_num, 1, temp_quad_ne, temp_quad_ele(1:temp_quad_ne), temp_cp, temp_dis, &
                        !                          quad_nn, quad_node, quad_ne, quad_ele, quad_cp, quad_dis, quad_vol, temp_coord, row, temp_value)
                        call vol_interpolation_3D_temp(rn, ce, j, mat_num, quad_ne, quad_vol, temp_coord, row, temp_value)
                        new_value(:,(i-1)*8+j) = temp_value
                    enddo
                    deallocate (temp_quad_ele)
                enddo
            endif
            call update_states(mat_num, row, col, new_value)

            deallocate(new_value, temp_value)
            deallocate(mat_elem_num)
        enddo
        deallocate(quad_node, quad_ele)
        deallocate (quad_dis, quad_cp, quad_vol)
        deallocate (new_quad_cp, new_quad_dis)
    endif 
    deallocate(temp_coord)
    deallocate(ori_coord, at, vt, ut)
else
    call update_coord(rn, num_dim, nn, coord)
endif

deallocate(coord, ele)

end subroutine fvm_routine_bak
 

subroutine area_interpolation(rn, mat_num, temp_quad_ne, temp_quad_ele, quad_nn, quad_node, quad_ne, quad_ele, temp_node, row, temp_value, cycle_num)

use materials
use physical_domain

implicit none

integer, intent(in) :: rn, mat_num, row, quad_nn, quad_ne, quad_ele(4,quad_ne), cycle_num
integer, intent(in) :: temp_quad_ne, temp_quad_ele(temp_quad_ne)
real, intent(in) :: quad_node(2,quad_nn), temp_node(2,4)
real, intent(inout) :: temp_value(row)

integer :: i, j, k, q, temp_ne, ori_ne, pos, val_point, num
real :: area, local_area, quad_area
real :: quad_coord(2,4), vec(2,3), node(2,3), ratio(2), ori_val(row)
integer :: temp_ele(30), temp_pn(30)
logical :: check

temp_ne = 0;  temp_ele = 0;  temp_pn = 0
vec(:,1) = temp_node(:,2) - temp_node(:,1)
vec(:,2) = temp_node(:,3) - temp_node(:,4)
call calc_area(4, temp_node, area)
do i = 1, cycle_num
    ratio(1) = (1/float(cycle_num))*0.5 + (i-1)*(1/float(cycle_num))

    node(:,1) = temp_node(:,1) + vec(:,1)*ratio(1)
    node(:,2) = temp_node(:,4) + vec(:,2)*ratio(1)
    vec(:,3) = node(:,2) - node(:,1)
    do j = 1, cycle_num
        ratio(2) = (1/float(cycle_num))*0.5 + (j-1)*(1/float(cycle_num))
        node(:,3) = node(:,1) + vec(:,3)*ratio(2)
        
        ! node(:,3)가 어떤 quad_element에 속해 있는지 검사
        check = .FALSE.
        do k = 1, temp_ne
            quad_coord = quad_node(:,quad_ele(:,temp_ele(k)))
            call inout_check(4, quad_coord, node(:,3), check)
            if (check) then
                temp_pn(k) = temp_pn(k) + 1
                exit
            endif
        enddo
        if (check == .FALSE.) then
            !find_quad_ele: do k = 1, temp_quad_ne
            !    num = temp_quad_ele(k)
            !    do q = (num-1)*4+1, num*4
                do q = 1, quad_ne
                    quad_coord = quad_node(:,quad_ele(:,q))
                    call inout_check(4, quad_coord, node(:,3), check)
                    if (check) then
                        temp_ne = temp_ne + 1
                        temp_ele(temp_ne) = q
                        temp_pn(temp_ne) = 1
                        !exit find_quad_ele
                        exit
                    endif
                enddo
            !enddo find_quad_ele
        endif
    enddo
enddo

temp_value = 0.0
do i = 1, temp_ne
    !write (*,'(A,I6,3X,A,I6,3X,A,F12.7)') 'element :', temp_ele(i), 'point number : ', temp_pn(i), 'area :', area/float(cycle_num**2) * temp_pn(i)
    pos = mod(temp_ele(i),4)
    if (pos == 0) then
        ori_ne = temp_ele(i)/4
        pos = 4
    else
        ori_ne = temp_ele(i)/4 + 1
    endif
    local_area = area/float(cycle_num**2) * float(temp_pn(i))
    !local_area = float(temp_pn(i))/float(cycle_num**2) 
    quad_coord = quad_node(:,quad_ele(:,temp_ele(i)))
    call calc_area(4, quad_coord, quad_area)
    
    call find_point_number(rn, ori_ne, pos, val_point)
    call get_variable( mat_num , val_point, row, ori_val )
    temp_value = temp_value + (ori_val*(local_area/quad_area))
    !if (row <= 10) then
    !    temp_value = temp_value + ori_val*local_area/area
    !else
    !    temp_value(1:10) = temp_value(1:10) + ori_val*local_area/area
    !    temp_value(11:row) = temp_value(11:row) + (ori_val*(local_area/quad_area))
    !endif
    
    !write (*,*) 'local_area, quad_area : ', local_area, quad_area
enddo

end subroutine area_interpolation


subroutine surface_disp_interpolation(Dim, num_rn, current_step, target_nn, target_coord)

use surface_domain

implicit none

integer, intent(in) :: Dim, num_rn, current_step, target_nn
real, intent(inout) :: target_coord(Dim,target_nn)

integer :: i, j, k, nn, ne, error, temp_ne, cn, cadjn, ele_node, num(2)
real :: local_val, min_val, in_coord(3), vec(3,3)
real :: temp_coord(Dim,(Dim-1)*4), temp_ut(Dim,(Dim-1)*4),vt(Dim,(Dim-1)*4), at(Dim,(Dim-1)*4)
logical :: in_element
character(len=40) :: fn

integer, allocatable :: elem(:,:)
real, allocatable :: coord(:,:), ut(:,:)

!fn = './output/solid/ale/fluid_surf_00_000000.plt'
!write (fn(31:32), '(I2.2)' ) rn
!write (fn(34:39), '(I6.6)' ) current_step
!open (Unit=30, File=fn, STATUS='replace', ACTION='write')

ele_node = (Dim-1)*4
nn = 0;  ne = 0
do i = 1, num_rn
    nn = nn + region(i)%num_nodes
    ne = ne + region(i)%num_elements
enddo
allocate(coord(Dim,nn), elem(ele_node,ne))
nn = 0;  ne = 0
do i = 1, num_rn
    num(1) = region(i)%num_nodes
    num(2) = region(i)%num_elements
    coord(:,nn+1:nn+num(1)) = region(i)%node_coord
    elem(:,ne+1:ne+num(2)) = region(i)%element(4:3+ele_node,:)
    nn = nn + num(1)
    ne = ne + num(2)
enddo

allocate (ut(Dim,nn))
nn = 0
do i = 1, num_rn
    num(1) = region(i)%num_nodes
    do j = 1, num(1)
        do k = 1, Dim
	        ut(k,nn+j) =  history(i)%data(region(i)%node(3+k,j), 1)
        enddo
    enddo
    nn = nn + num(1)
enddo
vt = 0.0
at = 0.0

do i = 1, target_nn
    in_coord = target_coord(:,i)
    min_val = 10e8
    do j = 1, ne
        do k = 1, ele_node
			temp_coord(:,k) = coord(:,elem(k,j))
        enddo
        call inout_check_cubic(in_coord, temp_coord, in_element, local_val, 1.0e-03)
        if (in_element) then
			do k = 1, ele_node
				temp_ut(:,k) = ut(:,elem(k,j))
            enddo
            if (Dim == 2) then
				call goto_NC(temp_coord, in_coord, temp_ut, vt, at, vec, error )
            elseif (Dim == 3) then
                call goto_NC_3D(temp_coord, in_coord, temp_ut, vt, at, vec, error )
            endif
			if ( error /= 1 ) then
                !target_disp(:,i) = vec(:,1)
                target_coord(:,i) = target_coord(:,i) + vec(:,1)
				exit
            endif
        else
			if ( min_val > local_val ) then
				min_val = local_val
				temp_ne = j
			endif
        endif
        if ( j == ne ) then
			do k = 1, ele_node
				temp_coord(:,k) = coord(:,elem(k,temp_ne))
				temp_ut(:,k) = ut(:,elem(k,temp_ne))
            enddo
            if (Dim == 2) then
				call outelement(temp_coord, in_coord, vec, temp_ut, vt, at, error)
            elseif (Dim == 3) then
                !write (*,*) i, ' node: outelement_3D. ', temp_ne
                call outelement_3D(temp_coord, in_coord, vec, temp_ut, vt, at, error)
            endif
                    
            if (error == 1) stop 'Error: Not convergence at subroutine goto_NC'
            target_coord(:,i) = target_coord(:,i) + vec(:,1)
			exit
        endif
    enddo
enddo

end subroutine surface_disp_interpolation


subroutine fvm_routine_case(Dim, dt)

use ale_domain
implicit none

integer, intent(in) :: Dim
real, intent(in) :: dt

integer :: i, j, k, rn, error
integer :: position(2)
 
real :: vec(2,3), temp_coord(2,4), temp_at(2,4), temp_vt(2,4), temp_ut(2,4), in_coord(2)
real :: dist, min_area, local_area, tol
logical :: in_element, check

! mesh data
integer :: nn, ne, temp_quad_ne, mat_ne, mat_nn, cur_ne, temp_ne, cn, ce, elem_num, node_num
integer, allocatable :: ele(:,:), c_nodes(:), c_elem(:)
real, allocatable :: coord(:,:), ori_coord(:,:), at(:,:), vt(:,:), ut(:,:)

if (Dim == 2) then
    rn = 2
    nn = mesh_set(rn)%nn
    ne = mesh_set(rn)%ne
    allocate(coord(Dim,nn), ele((Dim-1)*4,ne))
    coord = mesh_set(rn)%coord
    ele = mesh_set(rn)%ele
    !
    !allocate(ori_coord(2,nn), at(2,nn), vt(2,nn), ut(2,nn))
    !do i=1,nn
	   ! ut(1,i) =  history(rn)%data(region(rn)%node(4,i), 1)
	   ! ut(2,i) =  history(rn)%data(region(rn)%node(5,i), 1)
    !
	   ! vt(1,i) =  history(rn)%data(region(rn)%node(4,i), 2)
	   ! vt(2,i) =  history(rn)%data(region(rn)%node(5,i), 2)
	   !
	   ! at(1,i) =  history(rn)%data(region(rn)%node(4,i), 3)
	   ! at(2,i) =  history(rn)%data(region(rn)%node(5,i), 3)
    !enddo
    !ori_coord = region(rn)%node_coord(1:2,:)
    !
    !cn = inte_set%cn
    !ce = inte_set%ce
    !allocate (c_nodes(cn), c_elem(ce))
    !c_nodes = inte_set%c_nodes
    !c_elem = inte_set%c_elem
    !call inner_smooth(2, cn, c_nodes, nn, coord)
    !tol = mesh_set(rn)%ave_d(1)
    !do i = 1, cn
    !    node_num = c_nodes(i)
	   ! call calc_len(ori_coord(:,node_num), coord(:,node_num), dist)
	   ! if (dist > tol*10e-5) then
		  !  vec= 0.0
		  !  in_coord = coord(:,node_num)
		  !  min_area = 10e8
		  !  do j= 1, ce
    !            elem_num = c_elem(j)
			 !   do k = 1, 4
				!    temp_coord(:,k) = ori_coord(:,ele(k,elem_num))
			 !   enddo
			 !   call inelement(temp_coord, in_coord, in_element, local_area, 1.0e-03)
			 !   if (in_element) then
				!    do k = 1, 4
				!	    temp_ut(:,k) = ut(:,ele(k,j))
				!	    temp_vt(:,k) = vt(:,ele(k,j))
				!	    temp_at(:,k) = at(:,ele(k,j))
				!    enddo
    !    
				!    call goto_NC(temp_coord, in_coord, temp_ut, temp_vt, temp_at, vec, error )
    !                if ( error /= 1 ) then
				!	    position(1)= region(rn)%node(4,node_num)
				!	    position(2)= region(rn)%node(5,node_num)
    !    
				!	    call  update_motion( rn, 2, position , vec )
				!	    exit
				!    endif
			 !   else
				!    if ( min_area > local_area ) then
				!	    min_area = local_area
				!	    temp_ne = elem_num
				!    endif
			 !   endif
			 !   if ( j == ce ) then
				!    !print *, i , coord(:,i), temp_ne, ':: Cannot Found!! '
				!    do k = 1, 4
				!	    temp_coord(:,k) = ori_coord(:,ele(k,temp_ne))
				!	    temp_ut(:,k) = ut(:,ele(k,temp_ne))
				!	    temp_vt(:,k) = vt(:,ele(k,temp_ne))
				!	    temp_at(:,k) = at(:,ele(k,temp_ne))
				!    enddo
				!    call outelement(temp_coord, in_coord, vec, temp_ut, temp_vt, temp_at, error)
    !    
				!    position(1)= region(rn)%node(4,node_num)
				!    position(2)= region(rn)%node(5,node_num)
				!    call update_motion( rn, 2, position , vec )
				!    !coord(:,i) = coord(:,i) - vec(:,1)
				!    exit
				!    !STOP 'ERROR!!!!!!!!!!'
			 !   endif
    !        enddo
	   ! endif
    !enddo
    call update_coord(rn, Dim, nn, coord )
    deallocate(coord, ele)
endif       
            
end subroutine fvm_routine_case



subroutine vol_interpolation_3D(rn, ce, cq, mat_num, flag, temp_quad_ne, temp_quad_ele, temp_cp, temp_dis, &
                                quad_nn, quad_node, quad_ne, quad_ele, quad_cp, quad_dis, quad_vol, ball_coord, row, temp_value)

use materials
use physical_domain

implicit none

integer, intent(in) :: rn, ce, cq, mat_num, flag, row, quad_nn, quad_ne, quad_ele(8,quad_ne)
integer, intent(in) :: temp_quad_ne, temp_quad_ele(temp_quad_ne)
real, intent(in) :: quad_node(3,quad_nn), ball_coord(3,8), temp_cp(3), temp_dis, quad_cp(3,quad_ne), quad_dis(quad_ne), quad_vol(quad_ne)
real, intent(inout) :: temp_value(row)

integer :: i, j, k, num, val_point, ne, temp_num
real :: overlap_volume, accumulate_volume, ball_volume, target_volume, dis(2)
real :: target_coord(3,8), target_value(row), temp_elem(quad_ne)
logical :: check(2)
character(len=80) :: fn
    
call calc_cubic_volume(ball_coord, ball_volume)

accumulate_volume = 0.0
temp_value = 0.0
if (flag == 1) then
    ne = temp_quad_ne
elseif (flag == 2) then
    ne = quad_ne/8
endif
temp_num = 0
temp_elem = 0
do i = 1, ne
    check = .TRUE.
    if (flag == 2) then
        if (mat_num /= region(rn)%element(3,i)) then
            check(1) = .FALSE.
        endif
    endif
    
    if (check(1)) then
        do j = 1, 8
            if (flag == 1) then
                num = (temp_quad_ele(i)-1)*8+j
            elseif (flag == 2) then
                num = (i-1)*8+j
            endif
            target_coord = quad_node(:,quad_ele(:,num))
            
            call calc_length_3D(temp_cp, quad_cp(:,num), dis(1))
            dis(2) = temp_dis+quad_dis(num)
            if (dis(1) > dis(2)) then
                check(2) = .FALSE.
            endif
            
            if (check(2)) then
                temp_num = temp_num + 1
                temp_elem(temp_num) = num
                
                call calc_overlap_volume(ce, cq, temp_num, ball_coord, target_coord, overlap_volume)
                accumulate_volume = accumulate_volume + overlap_volume
                !if (ce == 160 .and. cq == 2) then
                !    write (*,*) 'overlap_volume:', temp_num, overlap_volume
                !    write (*,*) target_value*overlap_volume/target_volume
                !endif
                if (abs(overlap_volume) > 10e-8*target_volume) then
                    if (flag == 1) then
                        call find_point_number(rn, temp_quad_ele(i), j, val_point)
                    elseif (flag == 2) then
                        call find_point_number(rn, i, j, val_point)
                    endif
                    call get_variable(mat_num, val_point, row, target_value)
                    temp_value = temp_value + target_value*overlap_volume/quad_vol(num)
                endif
            endif
        enddo
    endif
enddo

!if (flag == 2) then
!    fn = './output/solid/remesh/overlap/check_overlap_elem_000000_000.plt'
!    write (fn(50:55),'(I6.6)') ce
!    write (fn(57:59),'(I3.3)') cq
!    open (Unit=30, File=fn, STATUS='replace', ACTION='write')
!    Write (30,'(A,A,A)') 'TITLE="Check overlap element(3D)"'
!    Write (30,*) 'VARIABLES="x", "y", "z"'
!    Write (30,'(A,I6,A,I6,A)') 'ZONE T = "ball_elem", N = ', temp_num*8, ',E = ', temp_num, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
!    do i = 1, temp_num
!        do j = 1, 8
!            write (30,'(3ES15.7)') quad_node(:,quad_ele(j,temp_elem(i)))
!        enddo
!    enddo
!    do i = 1, temp_num
!        write (30,*) (j, j=(i-1)*8+1,i*8)
!    enddo
!    Write (30,'(A)') 'ZONE T = "ball_elem", N = 8, E = 1, DATAPACKING = POINT, ZONETYPE = FEBRICK'
!    do k = 1, 8
!        write (30,'(3ES15.7)') ball_coord(:,k)
!    enddo
!    write (30,*) '1, 2, 3, 4, 5, 6, 7, 8'
!    close(30)
    !write (*,*)
    !write (*,*) '      ball_volume:', ball_volume
    !write (*,*) 'accumulate_volume:', accumulate_volume
!    stop
!endif

end subroutine vol_interpolation_3D

                                
subroutine vol_interpolation_3D_temp(rn, ce, cq, mat_num, quad_ne, quad_vol, ball_coord, row, temp_value)
                                     
use materials
use physical_domain
implicit none

integer, intent(in) :: rn, ce, cq, mat_num, quad_ne, row
real, intent(in) :: quad_vol(quad_ne), ball_coord(3,8)
real, intent(inout) :: temp_value(row)

integer :: num, val_point
real :: ball_volume, target_value(row)

num = (ce-1)*8+cq
call calc_cubic_volume(ball_coord, ball_volume)

call find_point_number(rn, ce, cq, val_point)
                    
call get_variable(mat_num, val_point, row, target_value)
temp_value = target_value*ball_volume/quad_vol(num)

end subroutine vol_interpolation_3D_temp