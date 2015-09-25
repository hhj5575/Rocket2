module ale_domain 

use materials
use coarse_domain
use surface_domain
use output_domain

implicit none
save

type ale_type
  integer, allocatable :: xadj(:), adjncy(:), ele_adj(:), ele_adjn(:)
end type

type restart_type
    integer :: line_num, cons_num(3), load_num(3,3)
    integer, allocatable :: cons1(:,:), cons2(:,:), cons3(:,:)
    integer, allocatable :: load_info1(:,:), load_info2(:,:), load_info3(:,:)
    real, allocatable :: load1(:), load2(:), load3(:)
    character(len=180) :: line_str(200)
end type

type(ale_type), allocatable :: ale_set(:)
type(restart_type) :: rest_set

integer :: pre_remesh_step = 0

contains

!=======================================================================
subroutine check_ale_set(dim, rn) 

implicit none

integer, intent(in) :: dim, rn

integer :: nn, ne, ele_node, ele_adjn, count, pre_num
integer, allocatable :: ele(:,:), ele_num(:), near_ele(:,:), temp_ele(:), ele_adj(:)

integer :: i, j, k, numflag, type

ele_node = (dim-1)*4
nn = region(rn)%num_nodes
ne = region(rn)%num_elements
allocate ( ele(ele_node,ne) ) 
ele = region(rn)%element(4:3+ele_node,:)

allocate( ale_set(rn)%xadj(ne+1),  ale_set(rn)%adjncy(ele_node*ne) )
if (dim == 2) then
    allocate( ale_set(rn)%ele_adj(nn+1),  ale_set(rn)%ele_adjn(ele_node*nn) )
    ale_set(rn)%ele_adj= 0 ;  ale_set(rn)%ele_adjn = 0
    
    numflag = 1
    type = 4
    call METIS_MeshTODual(ne, nn, ele, type, numflag, ale_set(rn)%xadj, ale_set(rn)%adjncy )

    count = 0
    do i = 1, nn
        ale_set(rn)%ele_adj(i) = count + 1
        do j = 1, ne
            do k = 1, 4
                if ( i == ele(k,j) ) then
                    count = count + 1
                    ale_set(rn)%ele_adjn(count) = j
                    exit
                endif
            enddo
        enddo
    enddo
    ale_set(rn)%ele_adj(nn+1)= count + 1
    write (*,*) ale_set(rn)%ele_adj
else
    allocate (ele_num(nn), near_ele(20,nn), temp_ele(160), ele_adj(60))
    allocate( ale_set(rn)%ele_adjn(ne+1),  ale_set(rn)%ele_adj(30*ne) )
    ale_set(rn)%ele_adjn= 0 ;  ale_set(rn)%ele_adj = 0
    ele_num = 0
    near_ele = 0
    do i = 1, ne
        do j = 1, 8
            ele_num(ele(j,i)) = ele_num(ele(j,i)) + 1
            near_ele(ele_num(ele(j,i)),ele(j,i)) = i
        enddo
    enddo
    ale_set(rn)%ele_adjn(1) = 1
    do i = 1, ne
        count = 0
        do j = 1, 8
            temp_ele(count+1:count+ele_num(ele(j,i))) = near_ele(1:ele_num(ele(j,i)),ele(j,i))
            count = count + ele_num(ele(j,i))
        enddo
        call sortc(temp_ele(1:count), count)
        pre_num = 0
        ele_adjn = 0
        do j = 1, count
            if (pre_num /= temp_ele(j)) then
                ele_adjn = ele_adjn + 1
                ele_adj(ele_adjn) = temp_ele(j)
                pre_num = temp_ele(j)
            endif
        enddo
        ale_set(rn)%ele_adj(ale_set(rn)%ele_adjn(i):ale_set(rn)%ele_adjn(i)+ele_adjn-1) = ele_adj(1:ele_adjn)
        ale_set(rn)%ele_adjn(i+1) = ale_set(rn)%ele_adjn(i)+ele_adjn
    enddo
    deallocate (ele_num, near_ele, temp_ele, ele_adj)
endif
deallocate (ele) 

end subroutine check_ale_set

!=======================================================================

subroutine fvm_routine(Dim, rn, dt)

implicit none 

integer, intent(in) :: Dim, rn
real,intent(in) :: dt

integer :: nn, ne
integer :: num_dim, ele_node, num_quad , st, ed  
integer :: element( 4 )

integer ::  num_face, size,val_size, bound
integer :: i, j, k , num_mat, col , row , error, temp_ne, num_group
integer :: position(2)
real :: vec(2,3), tol

real :: temp_coord(2,4), temp_at(2,4), temp_vt(2,4), temp_ut(2,4), in_coord(2)
real,allocatable :: area(:), coord(:,:), ori_coord(:,:), at(:,:), vt(:,:), ut(:,:), new_ut(:,:)
real,allocatable ::  avg(:) , velocity(:,:), local_velocity(:,:), local_node(:,:) , value(:,:)
real :: dist, min_area, local_area
real,allocatable :: new_value(:,: ), temp_value(:,:), ele(:,:), temp_val(:)
logical :: in_element
character(len=50) :: fn

ele_node = region(rn)%nodes_element
num_quad = region(rn)%num_quad_pts
num_dim = Dim
nn = mesh_set(rn)%nn
ne = mesh_set(rn)%ne
allocate(coord(2,nn), ele(4,ne))
allocate(ori_coord(2,nn), at(2,nn), vt(2,nn), ut(2,nn), new_ut(2,nn))
coord = mesh_set(rn)%coord
ele = mesh_set(rn)%ele

num_mat = pdomain_mat_db%num_materials

if ( num_mat .NE. 1 ) then
!    STOP "Too many Materials  "
endif

row = material(1)%num_row_state
col = material(1)%num_col_state

allocate( area( num_quad),  velocity(2,nn)  ) 
allocate( temp_value( row, num_quad), temp_val(row) ) 
temp_value = 0.0 

if ( num_quad .EQ. 1 ) then
	size = ele_node
	num_face = size 
	val_size = num_face+1
elseif ( ( num_quad .EQ. 3) .AND. ( ele_node .EQ. 3 )  ) then
	size = 7
	num_face = 9 
	val_size = num_face
elseif ( ( num_quad .EQ. 4) .AND. ( ele_node .EQ. 4 )  ) then
	size = 9
	num_face = 12
	val_size = num_face
endif

do i=1,nn
	ut(1,i) =  history(rn)%data(  region(rn)%node(4,i) , 1 )
	ut(2,i) =  history(rn)%data(  region(rn)%node(5,i) , 1 )
    
	vt(1,i) =  history(rn)%data(  region(rn)%node(4,i) , 2 )
	vt(2,i) =  history(rn)%data(  region(rn)%node(5,i) , 2 )
	
	at(1,i) =  history(rn)%data(  region(rn)%node(4,i) , 3 )
	at(2,i) =  history(rn)%data(  region(rn)%node(5,i) , 3 )
enddo
!ori_coord = region(rn)%node_coord(1:2, :) + ut
ori_coord = region(rn)%node_coord(1:2, :)

num_group = region(rn)%num_groups
tol = sum(mesh_set(rn)%ave_d(1:num_group))/float(num_group)
do  i = 1, nn
	call calc_len(ori_coord(:,i), coord(:,i), dist)
	if (dist > tol*10e-4) then
		st = 1;  ed = ne
		vec= 0.0
		in_coord = coord(:,i)
		min_area = 10e8
		do j= st , ed
			do k = 1, 4
				temp_coord(:,k) = ori_coord(:,ele(k,j))
			enddo
			call inelement(temp_coord, in_coord, in_element, local_area, 1.0e-03)
			if (in_element) then
				do k = 1, 4
					temp_ut(:,k) = ut(:,ele(k,j))
					temp_vt(:,k) = vt(:,ele(k,j))
					temp_at(:,k) = at(:,ele(k,j))
				enddo

				call goto_NC(temp_coord, in_coord, temp_ut, temp_vt, temp_at, vec, error )
                
				if ( error /= 1 ) then
					position(1)= region(rn)%node(4,i)
					position(2)= region(rn)%node(5,i)

					call  update_motion( rn, 2, position , vec )
                    new_ut(:,i) = vec(:,1)
					!coord(:,i) = coord(:,i) - vec(:,1)
					exit
				endif
			else
				if ( min_area > local_area ) then
					min_area = local_area
					temp_ne = j
				endif
			endif
			if ( j == ed ) then
				!print *, i , coord(:,i), temp_ne, ':: Cannot Found!! '
				do k = 1, 4
					temp_coord(:,k) = ori_coord(:,ele(k,temp_ne))
					temp_ut(:,k) = ut(:,ele(k,temp_ne))
					temp_vt(:,k) = vt(:,ele(k,temp_ne))
					temp_at(:,k) = at(:,ele(k,temp_ne))
				enddo
				call outelement(temp_coord, in_coord, vec, temp_ut, temp_vt, temp_at, error)

				position(1)= region(rn)%node(4,i)
				position(2)= region(rn)%node(5,i)
				call update_motion( rn, 2, position , vec )
                new_ut(:,i) = vec(:,1)
				!coord(:,i) = coord(:,i) - vec(:,1)
				exit
				!STOP 'ERROR!!!!!!!!!!'
			endif
        enddo
	else
		!coord(:,i) = coord(:,i) - ut(:,i)
	endif
enddo
call update_coord(rn, Dim, nn, coord )

if ( row .NE. 0 ) then 
  velocity = ( coord - ori_coord ) /dt 

  allocate( avg(num_face) , local_velocity(num_dim,size), local_node(num_dim,size) , value(row,val_size) ) 
  allocate( new_value(row, col )  ) 
  new_value = 0.0 

  do i=1,ne
! Set local node & veloocity info ( 9 point ) 
        temp_value = 0.0 
        call get_ele_connect( rn , i , ele_node , element )
        call set_local_info( nn, num_dim, coord, ele_node, element,size, local_node ) 
        call set_local_info( nn, num_dim, velocity, ele_node, element,size, local_velocity )
        call calc_avg_vel( num_dim, size, local_node, local_velocity,num_face,  avg )
        call calc_sub_area( local_node, size, area, num_quad )
        call find_near_value( rn, ele_node, num_quad, i, val_size,row, bound, value )
        call get_ele_quad( rn , i , num_quad , element )
        call solve_flux( num_quad, num_face, val_size, value, avg, row, temp_value, area, dt )        
        new_value(:,element) = temp_value
  enddo          !element end

  call update_states(rn, row, col, new_value )

  deallocate(avg,local_velocity,local_node,value) 
  deallocate(new_value) 
endif 

size = int(time/dt)
fn = "./output/solid/ale/2D_ale0000000.plt"
write ( fn(26:32), '(I7.7)' ) size
open (Unit=20, File=fn, STATUS='replace', ACTION='write')
Write (20,'(A,A,A)') 'TITLE="ale_check"'
Write (20,*) 'VARIABLES="x", "y", dx, dy'
Write (20,'(A,I5,A,I5, A)') 'ZONE N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
do i = 1, nn
    write (20,*) coord(:,i), new_ut(:,i)
enddo
do i = 1, ne
    write (20,*) ele(:,i)
enddo
close(20)

deallocate(area,velocity) 
deallocate(temp_value) 
deallocate (coord, ele)
deallocate(ori_coord, at, vt, ut, new_ut)

end subroutine fvm_routine

!============================================================================

subroutine find_near_value( rn, ele_node, num_quad,  my_num, vsize, num_stat, bound, values )

implicit none

integer :: rn, my_num,ele_node, num_face, num_quad, vsize, num_stat, bound


real :: values(num_stat, vsize) 
real :: temp_val(num_stat)

integer :: ele_num, num_near_ele 
integer :: i, j, val_point
integer :: my_element(2,ele_node) , near_element(10)
integer, allocatable :: near_elem(:)
integer :: position(2) , quad_table(10) 
integer :: element(4) 
logical :: sharing_node , check
integer :: mat_num

! near element Variable num 
!   face  :    1  2  3
!info(1,:) = (/ 4, 1, 2, 3 /)
!info(2,:) = (/ 5, 4, 1, 2 /)

bound = 0
call get_ele_connect( rn , my_num , ele_node , element )
my_element(1,:) = element(:)
my_element(2,1:(ele_node-1)) = my_element(1,2:ele_node)
my_element(2,ele_node) = my_element(1,1)
values = 0.0
quad_table = 0 
mat_num=1
num_face=ele_node ; ele_num=0 ; near_element = 0 
num_near_ele = ale_set(rn)%xadj(my_num+1) - ale_set(rn)%xadj(my_num)
allocate( near_elem(num_near_ele) )

do i=1,num_quad
        quad_table(i) = i
enddo
quad_table(num_quad+1) = 1

i =  ale_set(rn)%xadj(my_num)
j =  ale_set(rn)%xadj(my_num+1) - 1 

near_elem = ale_set(rn)%adjncy(i:j)
 
do i=1,num_near_ele
	ele_num = near_elem(i)
	call get_ele_connect( rn , ele_num , ele_node , element )
        
	do j=1,ele_node
		check = sharing_node( my_element(:,j) , element, ele_node )
		if ( check ) then
			near_element(j) = ele_num
			exit
		endif
	enddo
enddo 

if ( num_quad .EQ. 4 ) then 

	do i=1 , num_quad
		call find_point_number(rn,my_num,i,val_point)
		call get_variable( mat_num , val_point, num_stat, temp_val )
        values(:,i) = temp_val
    enddo
    
	do i=1,ele_node
		ele_num = near_element(i)
		if ( ele_num .EQ. 0 ) then
			call find_point_number(rn,my_num, quad_table(i) ,val_point)
			call get_variable( mat_num , val_point, num_stat, temp_val )
            values(:, 2*(i-1)+5 ) = temp_val
			call find_point_number(rn,my_num, quad_table(i+1) ,val_point)
			call get_variable( mat_num , val_point, num_stat, temp_val )
            values(:, 2*(i-1)+6 ) = temp_val
            bound = 1
		else
			position = my_element(:,i)
			call get_ele_connect(rn, ele_num, ele_node , element)
			call get_position( element , ele_node , position , 2)
			call set_quad_num( position, num_quad ) 
			call find_point_number(rn,ele_num, position(1) ,val_point)
			call get_variable( mat_num , val_point, num_stat, temp_val )
            values(:, 2*(i-1)+5 ) = temp_val
			call find_point_number(rn,ele_num, position(2) ,val_point)
			call get_variable( mat_num , val_point, num_stat, temp_val )
            values(:, 2*(i-1)+6 ) = temp_val 
		endif
	enddo

endif

deallocate(near_elem)

end subroutine find_near_value


subroutine solve_flux( num_quad, num_face, vsize, value, velocity, num_stat, new_value, area, dt) 

implicit none

integer,intent(in) :: num_quad, num_face, vsize,num_stat
real,intent(in) :: dt
real,intent(in) :: value(num_stat, vsize), velocity(num_face), area(num_quad)
real,intent(inout) :: new_value(num_stat, num_quad)
integer,allocatable :: ix_velocity(:,:), ix_near(:,:), direction(:,:)
integer :: i, j , local_face
!real, allocatable :: direction(:,:)
!real :: factor

if ( num_face < 4 ) then
	local_face =num_face 
else
	local_face = 4
endif

allocate( direction(num_quad,local_face) , ix_velocity(num_quad,local_face) ) 
allocate( ix_near(num_quad,local_face) )

if ( num_face .EQ. 9 ) then 

        direction(1,:)=(/ 1, 1, -1, 1 /)
        direction(2,:)=(/ 1, 1, -1, -1 /)
        direction(3,:)=(/ 1, 1, 1, 1 /)

        ix_velocity(1,:) = (/ 1, 7, 8, 6 /)
        ix_velocity(2,:) = (/ 2, 3, 9, 7 /)
        ix_velocity(3,:) = (/ 4, 5, 8, 9 /)

        ix_near(1,:)= (/ 4, 2, 3, 9 /)
        ix_near(2,:)= (/ 5, 6, 3, 1 /)
        ix_near(3,:)= (/ 7, 8, 1, 2 /)

elseif ( num_face .EQ. 12 ) then 
    
        !if (bound == 0) then
            !factor = 0.25
            !direction(1,:)=(/ 1-factor, 1.0+factor, -1.0-factor, 1-factor /)
            !direction(2,:)=(/ 1-factor, 1-factor, -1.0-factor, -1.0-factor /)
            !direction(3,:)=(/ 1.0+factor, 1-factor, 1-factor, -1.0-factor /)
            !direction(4,:)=(/ 1.0+factor, 1.0+factor, 1-factor, 1-factor /)
            direction(1,:)=(/ 1, 1, -1, 1 /)
            direction(2,:)=(/ 1, 1, -1, -1 /)
            direction(3,:)=(/ 1, 1, 1, -1 /)
            direction(4,:)=(/ 1, 1, 1, 1 /)
        !elseif (bound == 1) then
        !    factor = 0.75
        !    direction(1,:)=(/ 1-factor, 1, -1, 1-factor /)
        !    direction(2,:)=(/ 1-factor, 1-factor, -1, -1 /)
        !    direction(3,:)=(/ 1, 1-factor, 1-factor, -1 /)
        !    direction(4,:)=(/ 1, 1, 1-factor, 1-factor /)
        !    !write (*,*) ' boundary element : ', my_num
        !else
        !    STOP 'Error : ale_domain(solve_flux) - set wrong boundary!!!'
        !endif
                            
        ix_velocity(1,:) = (/ 1, 9, 12, 8 /)
        ix_velocity(2,:) = (/ 2, 3, 10, 9 /)
        ix_velocity(3,:) = (/ 10, 4, 5, 11 /)
        ix_velocity(4,:) = (/ 12, 11, 6, 7 /)
                    
        ix_near(1,:)= (/ 5, 2, 4, 12 /)
        ix_near(2,:)= (/ 6, 7, 3, 1 /)
        ix_near(3,:)= (/ 2, 8, 9, 4 /)
        ix_near(4,:)= (/ 1, 3, 10, 11 /)

elseif ( num_face .EQ. 4 ) then 

        direction(1,:)=(/ 1, 1, 1, 1 /)
        ix_velocity(1,:) = (/ 1, 2, 3, 4 /)     
        ix_near(1,:)= (/  2, 3, 4, 5 /) 
        
elseif ( num_face .EQ. 3 ) then 

        direction(1,:)=(/ 1, 1, 1 /)
        ix_velocity(1,:) = (/ 1, 2, 3 /)        
        ix_near(1,:)= (/ 2, 3, 4 /)     
        
endif 

do i = 1,num_quad
	do j=1,local_face 
		new_value(:,i)=new_value(:,i)+0.5*(value(:,ix_near(i,j))-value(:,i))*velocity(ix_velocity(i,j))*direction(i,j)
    enddo
	new_value(:,i) = dt * new_value(:,i) / area(i) + value(:,i) 
enddo

 deallocate(ix_velocity, ix_near, direction)

end subroutine solve_flux

subroutine back_up_main_file(Dim, unit_num)

implicit none
integer, intent(in) :: Dim, unit_num

integer :: i, num(3), status
integer, allocatable :: int_val(:,:)
real, allocatable :: real_val(:)
character(len=5) :: keyword

rewind(unit_num)
num = 0
rest_set%cons_num = 0
rest_set%load_num = 0
do 
    num(1) = num(1) + 1
    read(unit_num,'(A)', iostat=status) rest_set%line_str(num(1))
    if (status /= 0) exit
    keyword = rest_set%line_str(num(1))(1:5)
    if (keyword == '*stop') then
        exit
    elseif (keyword == '*cons') then
        num(2) = num(2) + 1
        backspace(unit_num)
        read(unit_num,*) keyword, rest_set%cons_num(num(2))
        allocate (int_val(Dim+1,rest_set%cons_num(num(2))))
        do i = 1, rest_set%cons_num(num(2))
            read(unit_num,*) int_val(:,i)
        enddo            
        if (num(2) == 1) then
            allocate(rest_set%cons1(Dim+1,rest_set%cons_num(num(2))))
            rest_set%cons1 = int_val
        elseif (num(2) == 2) then
            allocate(rest_set%cons2(Dim+1,rest_set%cons_num(num(2))))
            rest_set%cons2 = int_val
        elseif (num(2) == 3) then
            allocate(rest_set%cons3(Dim+1,rest_set%cons_num(num(2))))
            rest_set%cons3 = int_val
        endif
        deallocate (int_val)
    elseif (keyword == '*load') then
        num(3) = num(3) + 1
        backspace(unit_num)
        read(unit_num,*) keyword, rest_set%load_num(:,num(3))
        allocate (int_val(2,rest_set%load_num(1,num(3))))
        allocate (real_val(rest_set%load_num(1,num(3))))
        do i = 1, rest_set%load_num(1,num(3))
            read(unit_num,*) int_val(:,i), real_val(i)
        enddo
        if (num(3) == 1) then
            allocate(rest_set%load_info1(2,rest_set%load_num(1,num(3))))
            allocate(rest_set%load1(rest_set%load_num(1,num(3))))
            rest_set%load_info1 = int_val
            rest_set%load1 = real_val
        elseif (num(3) == 2) then
            allocate(rest_set%load_info2(2,rest_set%load_num(1,num(3))))
            allocate(rest_set%load2(rest_set%load_num(1,num(3))))
            rest_set%load_info2 = int_val
            rest_set%load2 = real_val
        elseif (num(3) == 3) then
            allocate(rest_set%load_info3(2,rest_set%load_num(1,num(3))))
            allocate(rest_set%load3(rest_set%load_num(1,num(3))))
            rest_set%load_info3 = int_val
            rest_set%load3 = real_val
        endif
        deallocate (int_val, real_val)
    endif
enddo
rest_set%line_num = num(1)

end subroutine back_up_main_file

!====================================================================================================

subroutine close_all_domain_file

implicit none

integer :: i, unit_num1 = 5004, unit_num2 = 5005

close(unit_num1) 
close(unit_num2) 

call free_material
if (number_of_contacts > 0) call free_coarse_domain
call free_sdomain
call free_ale_domain
call free_output_domain
call free_surface_domain

do i = 1, 3
    if (rest_set%cons_num(i) /= 0) then
        if (i == 1) then
            deallocate (rest_set%cons1)
        elseif (i == 2) then
            deallocate (rest_set%cons2)
        elseif (i == 3) then
            deallocate (rest_set%cons3)
        endif
    endif
    if (rest_set%load_num(1,i) /= 0) then
        if (i == 1) then
            deallocate (rest_set%load_info1, rest_set%load1)
        elseif (i == 2) then
            deallocate (rest_set%load_info2, rest_set%load2)
        elseif (i == 3) then
            deallocate (rest_set%load_info3, rest_set%load3)
        endif
    endif
enddo

end subroutine 
!=========================================================================================================

subroutine free_ale_domain

deallocate( ale_set )

end subroutine free_ale_domain

end module
