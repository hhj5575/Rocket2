subroutine remesh_interpolation1(rn, dt, solver_flag)
! 2nd order polynominal interpolation
use ale_domain

implicit none

integer, intent(in) :: rn
real, intent(in) :: dt
logical, intent(in) :: solver_flag

integer :: i, j, k , ele_node, num_quad, mat_num, num_var, col, count
integer :: st, ed , num_in, ar_size, error, node_num , num ,target_num 
integer :: ori_nn, ori_ne, new_nn, new_ne, temp_ne, flag
real :: calc_length, new_length , local_area, min_area, length
logical :: in_element

integer :: num_bound, unit_num1, unit_num2, node_unit_num

real, allocatable :: at(:,:), vt(:,:), ut( :,:), new_at(:,:),new_vt(:,:), new_ut( :,:)  , quad_coord(:,:)
real, allocatable :: ori_at(:,:), ori_vt(:,:), ori_ut(:,:), ori_theta(:), new_theta(:)
real,allocatable :: new_var(:,:), new_value(:,: ),  poly(:, :),quad(:,:),target(:,:) , variable(:,:) , temp(:,:) , pre(:,:), temp_var(:)
real :: tol, vec(2,3), node(2),  new_node(2), local_coord(2,4),  poly_mat(9,9),ele_coord(2,4)
real :: temp_ut(2,4), temp_vt(2,4), temp_at(2,4), temp_coord(2)

integer :: el_data(4), ele_nodes(4)
integer, allocatable :: target_xadj(:) , target_adj(:), target_quad(:),target_node(:)
integer,allocatable :: inner_node(:), ori_ele(:,:), new_ele(:,:), node_check(:)
real, allocatable :: ori_coord(:,:), new_coord(:,:)

unit_num1 = 5004
unit_num2 = 5005
node_unit_num = 5000

ori_nn = region(rn)%num_nodes
ori_ne = region(rn)%num_elements
ele_node = region(rn)%nodes_element
num_quad = region(rn)%num_quad_pts
num_bound = region(rn)%num_skin
tol = mesh_set(rn)%ave_d(1)
!open(unit=0830, file='check_state.dat')
!open(unit=0831, file='check_tdata.dat')

if ( rn == 1 ) then
    mat_num = 1
    num_var = material(mat_num)%num_row_state
    col = material(mat_num)%num_col_state
    new_nn = mesh_set(rn)%nn
    new_ne = mesh_set(rn)%ne
    allocate( new_at(2,new_nn), new_vt(2,new_nn), new_ut(2,new_nn) )
    allocate( new_value(new_ne*4, num_var )  )
    new_ut = 0.0
    new_vt = 0.0
    new_at = 0.0
    new_value= 0.0
    
    if (solver_flag) then
        allocate(at(2,ori_nn), vt(2,ori_nn), ut(2,ori_nn))
        allocate(ori_at(2,ori_nn), ori_vt(2,ori_nn), ori_ut(2,ori_nn))
        allocate(ori_coord(2,ori_nn), ori_ele(4,ori_ne))
        allocate(new_coord(2,new_nn), new_ele(4,new_ne))
        new_coord = mesh_set(rn)%coord
        new_ele = mesh_set(rn)%ele

        do i=1,ori_nn
            ut(1,i) =  history(rn)%data(  region(rn)%node(4,i) , 1 )
            ut(2,i) =  history(rn)%data(  region(rn)%node(5,i) , 1 )
    	
            vt(1,i) =  history(rn)%data(  region(rn)%node(4,i) , 2 )
            vt(2,i) =  history(rn)%data(  region(rn)%node(5,i) , 2 )
    	
            at(1,i) =  history(rn)%data(  region(rn)%node(4,i) , 3 )
            at(2,i) =  history(rn)%data(  region(rn)%node(5,i) , 3 )
        enddo
    
        !ori_coord = region(rn)%node_coord(1:2, :) + ut    
        ori_coord = region(rn)%node_coord(1:2, :)
        ori_ele = region(rn)%element(4:7,:)
    
        allocate(quad_coord(2,new_ne*4), target_node(new_ne*4) )

        num_in = ori_nn - num_bound
        allocate( inner_node(num_in) )    
        call check_inner_node(ori_nn, ori_coord, ori_ne, ori_ele, num_in, inner_node)

        do i = 1 ,new_ne
            ele_nodes = new_ele(:,i)
            ele_coord = new_coord(:, ele_nodes )
            call get_quad_coord( ele_coord )
            st = 4*(i-1)+1
            quad_coord( :, st:(st+3) ) = ele_coord 
        enddo

        do i = 1 , new_ne*4
            new_node = quad_coord(:,i)
            length = 1.0e9
            do j= 1, num_in
	            new_length = calc_length( new_node , ori_coord(:,inner_node(j)))

        !print *, i , inner_node(j), length, new_length 

	            if ( new_length < length ) then
		            length = new_length
		            target_node(i) = inner_node(j)
	            endif
            enddo
        enddo 

        allocate( target_xadj( num_in+1) )
        allocate( target_adj( new_nn*4) )

        count = 0 
        do i = 1, num_in
            node_num = inner_node(i)
            target_xadj(i) = count + 1 
            do j =1, new_ne*4
	            if ( node_num  == target_node(j) ) then
		            count = count + 1
		            target_adj(count) = j
	            endif
            enddo
            target_xadj(i+1) = count + 1 
        enddo

        do i = 1, num_in
            node_num = inner_node(i)
            node = ori_coord(:,node_num)
            count = 1  
            st = ale_set(rn)%ele_adj(node_num)
            ed = ale_set(rn)%ele_adj(node_num+1) -1
            num = ale_set(rn)%ele_adj(node_num+1) - ale_set(rn)%ele_adj(node_num)
            ar_size = num_quad*num
            allocate(quad(2, ar_size), variable( ar_size, num_var ), temp_var(num_var), stat = error)
            variable = 0.0
            do j =st,ed
	            num = ale_set(rn)%ele_adjn(j)
	            el_data = region(rn)%element(4:7, num)
	            local_coord = ori_coord(:,el_data)
	            call get_quad_coord( local_coord )
	            quad(:, count:(count+3) ) = local_coord
	            el_data = region(rn)%element(8:11, num)
                    
	            do k=1,4
                    temp_var = variable(count+k-1,:)
                    !temp_var = 0.0
                    call get_variable(mat_num, el_data(k),num_var,temp_var)
                    variable(count+k-1,:) = temp_var
	            enddo
                !do k=1,4
			    !    call get_variable( mat_num , el_data(k),num_var,variable( count+k-1, :)  )
		        !enddo
	            count = count + 4
            enddo
                    
            quad(1, :) = node(1) - quad(1, :)
            quad(2, :) = node(2) - quad(2, :)
            
            allocate( pre( 9, ar_size ), temp(ar_size, 9), new_var(9 , num_var) ,stat= error)
            do j= 1, ar_size
	            pre(1,j) = 1.0
	            pre(2,j) = quad(1, j)
	            pre(3,j) = quad(2, j)
	            pre(4,j) = quad(1, j)**2
	            pre(5,j) = quad(1, j)*quad(2, j)
	            pre(6,j) = quad(2, j)**2
	            pre(7,j) = (quad(1,j)**2) * quad(2,j)
	            pre(8,j) = quad(1, j)*( quad(2, j)**2)
	            pre(9,j) = ( quad(1, j)*quad(2, j) )**2
            enddo   

            do j = 1,9
	            temp(:, j) = pre( j,:)
            enddo
            
            new_var = MATMUL( pre, variable )
            
            deallocate( variable, temp_var )

            poly_mat = MATMUL( pre, temp )           
            deallocate( pre,temp  )

            call cholesky_ale(poly_mat, new_var, 9, num_var)
            
            st = target_xadj(i)
            ed = target_xadj(i+1) -1
            num = ed - st + 1

            allocate( target(2, num), target_quad(num) )
            target = 0.0
            target_num = 0
            
            target_quad = target_adj(st:ed)
            target = quad_coord(:, target_quad ) 

            allocate( poly( num , 9 ), stat=error )
            allocate( variable(num , num_var ),stat=error )
        
            target(1, :) = node(1) - target(1, :)
            target(2, :) = node(2) - target(2, :)

            do j= 1, num
	            poly(j,1) = 1.0
	            poly(j,2) = target(1,j)
	            poly(j,3) = target(2,j)
	            poly(j,4) = target(1,j)**2
	            poly(j,5) = target(1,j)*target(2,j)
	            poly(j,6) = target(2, j)**2
	            poly(j,7) = (target(1,j)**2)*target(2,j)
	            poly(j,8) = target(1,j)*(target(2,j)**2)
	            poly(j,9) = (target(1,j)*target(2,j))**2
            enddo   
            
            variable = MATMUL( poly, new_var )      

            new_value(target_quad,:) = variable
            
            deallocate( quad )
            deallocate( poly , variable, new_var, target, target_quad  )
        enddo
    
        !call write_tecplot(ori_nn, ori_coord, ori_ne, ori_ele, 1, 1, 5 ) 
        do i = 1, new_nn
            vec = 0.0
            temp_coord = new_coord(:,i)
            min_area = 10e8
            do j = 1, ori_ne
	            el_data = region(rn)%element(4:7, j)
                local_coord = ori_coord(:,el_data)

                call inelement(local_coord, temp_coord, in_element, local_area, tol*1.0e-4)
                if ( in_element )then
                    temp_ut = ut(:,el_data)
                    temp_vt = vt(:,el_data)
                    temp_at = at(:,el_data)

                    call goto_NC(local_coord, temp_coord, temp_ut, temp_vt, temp_at, vec, error )
                    new_at(:,i) = new_at(:,i) + vec(:,3)
                    new_vt(:,i) = new_vt(:,i) + vec(:,2)
                    new_ut(:,i) = new_ut(:,i) + vec(:,1)
                    !if (i == 9) new_ut(:,i) = (/ 0.040278901645480, 0.005967272745175 /)
                    exit
                else
                    if ( min_area > local_area ) then
                        min_area = local_area
                        temp_ne = j
                    endif
                endif

	            if ( j == ori_ne ) then
		            !print *, i , new_coord(:,i), temp_ne, ':: Cannot Found!! '
                    !STOP 'ERROR!!!!!!!!!!'
                    el_data = region(rn)%element(4:7, temp_ne)
                    local_coord = ori_coord(:,el_data)
                
                    temp_ut = ut(:,el_data)
                    temp_vt = vt(:,el_data)
                    temp_at = at(:,el_data)
                    call outelement(local_coord, temp_coord, vec, temp_ut, temp_vt, temp_at, error)
                    new_at(:,i) = new_at(:,i) + vec(:,3)
                    new_vt(:,i) = new_vt(:,i) + vec(:,2)
                    new_ut(:,i) = new_ut(:,i) + vec(:,1)
                    new_coord(:,i) = temp_coord
                    exit
	            endif
            enddo
        enddo
    
        allocate(ori_theta(ori_ne), new_theta(new_ne))
        allocate(node_check(new_nn))
        node_check = mesh_set(rn)%boundary
        call calc_ori_theta(ori_nn, ori_ne, ori_coord, ori_ele, ut, ori_theta)
        call ori_to_new_theta(ori_nn, ori_ne, ori_coord, ori_ele, ori_theta, new_nn, new_ne, new_coord, new_ele, new_theta, tol)
        call theta_interpolation(new_nn, new_ne, new_coord, new_ele, new_theta, node_check, new_ut, flag)
        deallocate (ori_theta, new_theta, node_check)
    
        deallocate(ut , vt , at )
        deallocate(ori_coord, ori_ele, new_coord, new_ele)
    endif

    do i = 1 , new_ne*4
        do j = 1 , num_var 
	        write (unit_num1) new_value(i,j)
        enddo
    enddo
    
    !new_coord = new_coord - new_ut
    do i= 1, new_nn
        write(unit_num2) new_ut(1,i)
        write(unit_num2) new_ut(2,i)
        write(unit_num2) new_vt(1,i)
        write(unit_num2) new_vt(2,i)
        write(unit_num2) new_at(1,i)
        write(unit_num2) new_at(2,i)
    enddo 
    
    !write (node_unit_num,'(A,I2)') '*node, ', rn
    !do i = 1, new_nn
    !    write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i,', ',new_coord(1,i),', ',new_coord(2,i)
    !enddo
    
    deallocate(new_value)
    deallocate(new_ut , new_vt , new_at )
else
    do i=1,ori_nn
	    write(unit_num2) history(rn)%data(  region(rn)%node(4,i) , 1 )
	    write(unit_num2) history(rn)%data(  region(rn)%node(5,i) , 1 )
    	
	    write(unit_num2) history(rn)%data(  region(rn)%node(4,i) , 2 )
	    write(unit_num2) history(rn)%data(  region(rn)%node(5,i) , 2 )
    	
	    write(unit_num2) history(rn)%data(  region(rn)%node(4,i) , 3 )
	    write(unit_num2) history(rn)%data(  region(rn)%node(5,i) , 3 )
    enddo
    
    !write (node_unit_num,'(A,I2)') '*node, ', rn
    !do i = 1, ori_nn
    !    write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i,', ',region(rn)%node_coord(1,i),', ', region(rn)%node_coord(2,i)
    !enddo
    
    do mat_num = 2, number_of_materials
        num_var=material(mat_num)%num_row_state
        col=material(mat_num)%num_col_state
        do i = 1, col
            do j = 1, num_var
                write (unit_num1) material(mat_num)%state1(j,i)
            enddo
        enddo
    enddo
endif

if (rn == num_of_regions) then
    !close(node_unit_num)
    call close_all_domain_file
endif
!close(0830)
!close(0831)

end subroutine remesh_interpolation1
    

subroutine remesh_interpolation2(Dim, rn, delta_t, solver_flag)

use ale_domain

implicit none

integer, intent(in) :: Dim, rn
real, intent(in) :: delta_t
logical, intent(in) :: solver_flag

integer :: i, j, k, ele_node, num_quad, mat_num, num_var,col, mat_model_num
integer :: ori_nn, ori_ne, new_nn, new_ne, temp_ne, flag, error
real :: length, local_val, min_val
logical :: in_element, check

integer :: num_bound, unit_num1, unit_num2 

real, allocatable :: at(:,:), vt(:,:), ut(:,:), new_at(:,:),new_vt(:,:), new_ut(:,:), pre_ut(:,:)
real :: tol, vec(Dim,3), local_coord(Dim,(Dim-1)*4), temp_cp(3), temp_dis
real :: temp_ut(Dim,(Dim-1)*4), temp_vt(Dim,(Dim-1)*4), temp_at(Dim,(Dim-1)*4), temp_coord(Dim)

integer :: size, num_face, val_size, num_dim, row, num_group, temp_quad_ne, ce, new_mat_num
integer :: el_data((Dim-1)*4), this(4)
integer,allocatable :: ori_ele(:,:), new_ele(:,:), temp_quad_ele(:), node_check(:)
real, allocatable :: ori_coord(:,:), new_coord(:,:), ori_theta(:), new_theta(:)

!quad data
integer :: quad_nn, quad_ne, element((Dim-1)*4), new_quad_nn, new_quad_ne
integer, allocatable :: quad_ele(:,:), new_quad_ele(:,:)
real, allocatable :: quad_node(:,:), quad_dis(:), quad_cp(:,:), quad_vol(:)
real, allocatable :: new_quad_node(:,:), new_quad_dis(:), new_quad_cp(:,:)
real,allocatable :: new_value(:,:), temp_value(:)
real :: temp_node(Dim,(Dim-1)*4), local_node(2,9)

unit_num1 = 5004
unit_num2 = 5005

ori_nn = region(rn)%num_nodes
ori_ne = region(rn)%num_elements
ele_node = region(rn)%nodes_element
num_quad = region(rn)%num_quad_pts
num_bound = region(rn)%num_skin
num_group = region(rn)%num_groups
tol = sum(mesh_set(rn)%ave_d(1:num_group))/float(num_group)
num_dim = Dim

!if ( num_quad .EQ. 1 ) then
!	size = ele_node
!	num_face = size 
!	val_size = num_face+1
!elseif ( ( num_quad .EQ. 3) .AND. ( ele_node .EQ. 3 )  ) then
!	size = 7
!	num_face = 9 
!	val_size = num_face
!elseif ( ( num_quad .EQ. 4) .AND. ( ele_node .EQ. 4 )  ) then
!	size = 9
!	num_face = 12
!	val_size = num_face
!endif
if (num_dim == 2) then
    size = 9
elseif (num_dim == 3) then
    size = 27
endif

if ( rn == 1 ) then
    new_nn = mesh_set(rn)%nn
    new_ne = mesh_set(rn)%ne
    mat_num = 1
    row = material(mat_num)%num_row_state
    col = material(mat_num)%num_col_state
    allocate(new_at(num_dim,new_nn), new_vt(num_dim,new_nn), new_ut(num_dim, new_nn), pre_ut(num_dim,new_nn))
    new_at = 0.0
    new_vt = 0.0
    new_ut = 0.0
    
    if (solver_flag) then
        allocate(at(num_dim,ori_nn), vt(num_dim,ori_nn), ut(num_dim,ori_nn))
        allocate(ori_coord(num_dim,ori_nn), ori_ele(ele_node,ori_ne))
        allocate(new_coord(num_dim,new_nn), new_ele(ele_node,new_ne))
        if (num_dim == 2) then
            new_coord = mesh_set(rn)%coord
            new_ele = mesh_set(rn)%ele
        elseif (num_dim == 3) then
            new_coord = mesh_set(rn)%node_coord
            new_ele = mesh_set(rn)%element
        endif
        ori_coord = region(rn)%node_coord(1:num_dim, :)
        ori_ele = region(rn)%element(4:3+ele_node,:)

        if (num_dim == 2) then
            do i=1,ori_nn
                do j = 1, num_dim
	                ut(j,i) =  history(rn)%data(region(rn)%node(3+j,i), 1)
	                vt(j,i) =  history(rn)%data(region(rn)%node(3+j,i), 2)
                    at(j,i) =  history(rn)%data(region(rn)%node(3+j,i), 3)
                enddo
            enddo
        elseif (num_dim == 3) then
            do i = 1,ori_nn
                ut(:,i) = sdomain(rn)%u_cartesian(region(rn)%node(4:6,i))
                vt(:,i) = 0.0
                at(:,i) = 0.0
            enddo
        endif
        write (*,*) 'start displacement interpolation'
        do i =1, new_nn
            vec = 0.0
            temp_coord = new_coord(:,i)
            min_val = 10e8
            do j = 1, ori_ne
	            el_data = region(rn)%element(4:3+ele_node, j)
                local_coord = ori_coord(:,el_data)
                
                if (num_dim == 2) then
                    call inelement(local_coord, temp_coord, in_element, local_val, tol*1.0e-4)
                elseif (num_dim == 3) then
                    call inout_check_cubic(temp_coord, local_coord, in_element, local_val, 1.0e-03)
                endif
                if ( in_element )then
                    temp_ut = ut(:,el_data)
                    temp_vt = vt(:,el_data)
                    temp_at = at(:,el_data)
                
                    if (num_dim == 2) then
				        call goto_NC(local_coord, temp_coord, temp_ut, temp_vt, temp_at, vec, error )
                    elseif (num_dim == 3) then
                        call goto_NC_3D(local_coord, temp_coord, temp_ut, temp_vt, temp_at, vec, error )
                    endif
                    
                    if (error /= 1) then
                        new_at(:,i) = vec(:,3)
                        new_vt(:,i) = vec(:,2)
                        new_ut(:,i) = vec(:,1)
                        exit
                    endif
                else
                    if ( min_val > local_val ) then
                        min_val = local_val
                        temp_ne = j
                    endif
                endif

	            if ( j == ori_ne ) then
		            !print *, i , new_coord(:,i), temp_ne, ':: Cannot Found!! '
                    !STOP 'ERROR!!!!!!!!!!'
                    el_data = region(rn)%element(4:3+ele_node, temp_ne)
                    local_coord = ori_coord(:,el_data)
                
                    temp_ut = ut(:,el_data)
                    temp_vt = vt(:,el_data)
                    temp_at = at(:,el_data)
                    
                    if (num_dim == 2) then
				        call outelement(local_coord, temp_coord, vec, temp_ut, temp_vt, temp_at, error)
                    elseif (num_dim == 3) then
                        !write (*,*) i, ' node: outelement_3D. ', temp_ne
                        call outelement_3D(local_coord, temp_coord, vec, temp_ut, temp_vt, temp_at, error)
                    endif
                    if (error == 1) stop 'Error: Not convergence at subroutine goto_NC'
                    new_at(:,i) = vec(:,3)
                    new_vt(:,i) = vec(:,2)
                    new_ut(:,i) = vec(:,1)
                    exit
	            endif
            enddo
        enddo
        
        if (num_dim == 2) then
            pre_ut = new_ut
            allocate(ori_theta(ori_ne), new_theta(new_ne))
            allocate(node_check(new_nn))
            node_check = mesh_set(rn)%boundary
            call calc_ori_theta(ori_nn, ori_ne, ori_coord, ori_ele, ut, ori_theta)
            call ori_to_new_theta(ori_nn, ori_ne, ori_coord, ori_ele, ori_theta, new_nn, new_ne, new_coord, new_ele, new_theta, tol)
            call theta_interpolation(new_nn, new_ne, new_coord, new_ele, new_theta, node_check, new_ut, flag)
            deallocate (ori_theta, new_theta, node_check)
        endif
        !if (flag == 2) new_ut = pre_ut
        write (*,*) 'End displacement interpolation'
        deallocate(ut, vt, at)

        ! state variables update by using area interpolation
        ele_node = region(rn)%nodes_element
        num_quad = region(rn)%num_quad_pts
        row = material(1)%num_row_state
        col = material(1)%num_col_state
 
        quad_nn = ori_ne*size;  quad_ne = ori_ne*ele_node
        allocate (quad_node(Dim,quad_nn), quad_ele(ele_node,quad_ne))
        call elem_to_quadelem(num_dim, ele_node, ori_nn, ori_ne, ori_coord, ori_ele, quad_nn, quad_ne, quad_node, quad_ele)
        if (num_dim == 2) then
            if ( row .NE. 0 ) then 
                allocate(temp_value(row))
                allocate(new_value(row, new_ne*ele_node)) 
                new_value = 0.0 
            
                temp_quad_ne = 1
                allocate (temp_quad_ele(temp_quad_ne)) 
                temp_quad_ele(1) = i
                do i = 1, new_ne
                    element = new_ele(:,i)
                    call set_local_info( new_nn, num_dim, new_coord, ele_node, element, size, local_node )
                    !temp_quad_ne = 1 + ale_set(rn)%ele_adj(i+1) - ale_set(rn)%ele_adj(i)
                    !temp_quad_ele(2:temp_quad_ne) = ale_set(rn)%ele_adjn(ale_set(rn)%ele_adj(i):ale_set(rn)%ele_adj(i+1)-1)
                    do j = 1, 4
                        this = (/ j, j+4, 9, j+3 /)
                        if (j == 1) this(4) = 8
                        temp_node = local_node(:,this)
                        call area_interpolation(rn, mat_num, temp_quad_ne, temp_quad_ele, quad_nn, quad_node, quad_ne, quad_ele, temp_node, row, temp_value, 25)
                        new_value(:,(i-1)*4+j) = temp_value
                    enddo
                enddo
                do i = 1, new_ne*ele_node
                    do j = 1, row 
	                    write (unit_num1) new_value(j,i)
                    enddo
                enddo 
                deallocate(new_value)
                deallocate(temp_quad_ele)
                deallocate (temp_value)
            endif
        elseif (num_dim == 3) then
            num_group = sub_mesh_region
            new_quad_nn = new_ne*size;  new_quad_ne = new_ne*ele_node
            allocate (new_quad_node(Dim,new_quad_nn), new_quad_ele(ele_node,new_quad_ne))
            allocate (quad_dis(quad_ne), new_quad_dis(new_quad_ne))
            allocate (quad_cp(3,quad_ne), new_quad_cp(3,new_quad_ne))
            allocate (quad_vol(quad_ne))
            call elem_to_quadelem(num_dim, ele_node, new_nn, new_ne, new_coord, new_ele, new_quad_nn, new_quad_ne, new_quad_node, new_quad_ele)
            call calc_quad_cp_dis(quad_nn, quad_node, quad_ne, quad_ele, quad_cp, quad_dis)
            call calc_quad_cp_dis(new_quad_nn, new_quad_node, new_quad_ne, new_quad_ele, new_quad_cp, new_quad_dis)
            call calc_quad_volume(quad_nn, quad_node, quad_ne, quad_ele, quad_vol)
            temp_quad_ne = 1
            allocate (temp_quad_ele(temp_quad_ne)) 
            temp_quad_ele(1) = i
            new_mat_num = 0
            do mat_num = 1, number_of_materials
                check = .TRUE.
                if (mat_num <= sub_mesh_region) then
                    if (mesh_set(1)%flag(mat_num) == 3) then
                        check = .FALSE.
                    endif
                endif
                if (check) then
                    new_mat_num = new_mat_num + 1
                    write (*,*) '=========', mat_num, 'material  ================'
                    row = material(mat_num)%num_row_state
                    mat_model_num = material(mat_num)%model_num
                    if ( row .NE. 0 .and. mat_model_num /= 31) then 
                        if (mat_num == 1) then
                            this(1:2) = (/ 1, mesh_set(rn)%mat_array(new_mat_num) /)
                        else
                            this(1:2) = (/ mesh_set(rn)%mat_array(new_mat_num-1)+1, mesh_set(rn)%mat_array(new_mat_num) /)
                        endif
                        this(3) = this(2)-this(1)+1
                        write (*,*) 'sp, lp, ne:', this(1:3)
                        allocate(temp_value(row))
                        allocate(new_value(row, this(3)*ele_node)) 
                        new_value = 0.0
                        temp_value = 0.0
                        if (mat_num <= number_of_materials-1) then
                            do i = 1, this(3)
                                !write (*,*) i+this(1)-1, ' element interpolation'
                                ce = this(1)+i-1
                                do j = 1, 8
                                    temp_node = new_quad_node(:,new_quad_ele(:,(ce-1)*8+j))
                                    temp_cp = new_quad_cp(:,(ce-1)*8+j)
                                    temp_dis = new_quad_dis((ce-1)*8+j)
                                    call vol_interpolation_3D(rn, ce, j, mat_num, 2, temp_quad_ne, temp_quad_ele, temp_cp, temp_dis, &
                                                              quad_nn, quad_node, quad_ne, quad_ele, quad_cp, quad_dis, quad_vol, temp_node, row, temp_value)
                                    new_value(:,(i-1)*8+j) = temp_value
                                enddo
                            enddo
                        endif
                        do j = 1, this(3)*ele_node
                            do k = 1, row 
	                            write (unit_num1) new_value(k,j)
                            enddo
                        enddo
                        deallocate (new_value, temp_value)
                    endif
                endif
            enddo
            
            deallocate (new_quad_node, new_quad_ele)
            deallocate (quad_dis, quad_cp, quad_vol)
            deallocate (new_quad_cp, new_quad_dis)
            deallocate (temp_quad_ele)
        endif
        deallocate(quad_node, quad_ele)
        deallocate(ori_coord, ori_ele, new_coord, new_ele)
    else
        new_mat_num = 0
        do mat_num = 1, number_of_materials
            check = .TRUE.
            if (mat_num <= sub_mesh_region) then
                if (mesh_set(1)%flag(mat_num) == 3) then
                    check = .FALSE.
                endif
            endif
            if (check) then
                new_mat_num = new_mat_num + 1
                row = material(mat_num)%num_row_state
                mat_model_num = material(mat_num)%model_num
                if ( row .NE. 0 .and. mat_model_num /= 31) then 
                    if (mat_num == 1) then
                        this(1:2) = (/ 1, mesh_set(rn)%mat_array(new_mat_num) /)
                    else
                        this(1:2) = (/ mesh_set(rn)%mat_array(new_mat_num-1)+1, mesh_set(rn)%mat_array(new_mat_num) /)
                    endif
                    this(3) = this(2)-this(1)+1
                    allocate(temp_value(row))
                    allocate(new_value(row, this(3)*ele_node)) 
                    new_value = 0.0
                    temp_value = 0.0
                    do j = 1, this(3)*ele_node
                        do k = 1, row 
	                        write (unit_num1) new_value(k,j)
                        enddo
                    enddo
                    deallocate (new_value, temp_value)
                endif
            endif
        enddo
    endif
    open (Unit=30, File='remesh_disp.plt', STATUS='replace', ACTION='write')
    do i= 1, new_nn
        write (30,*) new_ut(:,i)
    enddo
    close(30)
    do i= 1, new_nn
        do j = 1, num_dim
            write(unit_num2) new_ut(j,i)
        enddo
        do j = 1, num_dim
            write(unit_num2) new_vt(j,i)
        enddo
        do j = 1, num_dim
            write(unit_num2) new_at(j,i)
        enddo
    enddo
    deallocate(new_ut , new_vt , new_at )
else
    do i=1,ori_nn
	    write(unit_num2) history(rn)%data(  region(rn)%node(4,i) , 1 )
	    write(unit_num2) history(rn)%data(  region(rn)%node(5,i) , 1 )
    	
	    write(unit_num2) history(rn)%data(  region(rn)%node(4,i) , 2 )
	    write(unit_num2) history(rn)%data(  region(rn)%node(5,i) , 2 )
    	
	    write(unit_num2) history(rn)%data(  region(rn)%node(4,i) , 3 )
	    write(unit_num2) history(rn)%data(  region(rn)%node(5,i) , 3 )
    enddo
    do mat_num = 2, number_of_materials
        num_var=material(mat_num)%num_row_state
        col=material(mat_num)%num_col_state
        do i = 1, col
            do j = 1, num_var
                write (unit_num1) material(mat_num)%state1(j,i)
            enddo
        enddo
    enddo
endif    
if (rn == num_of_regions) call close_all_domain_file

end subroutine remesh_interpolation2