module physical_domain

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (October 28, 2006)
! Modification History: December 05, 2006 (by Jeeho Lee)
!                       October 31, 2009
!                       August 20, 2010 (pressure & temperature loads)
!                       July 07, 2015 (add 3D version features)

implicit none
save

type pdomain_type
  integer :: num_node_dof, num_dim, num_nodes, num_elements, nodes_element, num_quad_pts
  integer :: num_springs = 0
	integer :: num_dofs = 0, num_constr = 0, num_mpc_group = 0, max_bandwidth = 0
	integer :: num_skin = 0, num_skin_contact_nodes = 0
  integer :: num_groups = 1
  integer :: num_skin_surfaces = 0, num_inner_faces = 0
  integer :: size_str_output
	integer, allocatable :: node(:,:), node_number_table(:), node_group(:) ! starting node number in each group
	integer, allocatable :: element(:,:), element_group(:) ! starting element number in each group
  integer, allocatable :: spring(:,:)
  integer, allocatable :: element_face(:,:,:)  ! element_face: 3D only
  
	real :: body_force(3) = 0.0, grfactor(3) = 0.0
	real, allocatable :: node_coord(:,:)
	real, allocatable :: element_area(:), element_density(:), consis_mass_ratio(:), damping_factor(:,:)
	real, allocatable :: element_result(:,:,:) 
  real, allocatable :: spring_area(:)

	integer, allocatable :: active_dof(:), dof_to_array(:)
	integer, allocatable :: disp_bc(:)  ! inactive (constrained) nominal DOF number
  logical(1), allocatable :: dof_activity_flag(:)  ! TRUE: active DOF, FALSE: constrained DOF
	logical(1), allocatable :: rigid_constraints(:) ! TRUE: unchangeable rigid-body constrained nominal DOF

  real, allocatable :: mpc_theta(:)
  integer, allocatable :: mpc_element_addr(:), mpc_element_nodes(:), mpc_plane(:) ! MPC address in each element and nodal data info
  integer, allocatable ::  mpc_element_group(:)
  logical :: mpc_existence = .FALSE.

  integer :: skin_pressure_ends(2)  ! 2D skin pressure start(1) & end numbers(2) (skin point number)
	integer, allocatable :: skin_nodes(:), skin_node_group(:), skin_surfaces(:,:)  ! skin_surfaces: 3D only
  integer, allocatable :: skin_contact_nodes(:), skin_group_pointer(:), skin_node_map(:)  ! starting skin node number in each group
  real, allocatable :: skin_pressure(:)
	real, allocatable :: skin_normal(:,:), surface_normal(:,:)  ! surface_normal: 3D skin surface normal vector
  real :: skin_length
  logical(1), allocatable :: skin_load_set(:), surface_pressure_participation(:)

  logical :: dof_changed = .FALSE.
  logical(1), allocatable :: node_occupancy(:), element_participation(:), spring_participation(:)
  logical(1), allocatable :: element_face_participation(:,:), contact_face_participation(:,:)  ! element_face_participation: 3D only
end type

type pdomain_mat_db_type
  integer :: num_materials
	integer, allocatable :: mat_num_model_num(:)
end type

type load_set_type
	integer :: region_info = 0, load_type = 0, num_load = 0
	integer, allocatable :: load_bc(:), load_node(:)
	real, allocatable :: load(:)
  logical :: surface_pressure_existence = .FALSE.
end type

type disp_set_type
	integer :: region_info = 0
	real, allocatable :: constr(:)
end type

type temp_load_type 
	integer :: num_temp 
  real :: reference_temperature ! stress-free temperature
	real,allocatable :: step_time(:)
	real,allocatable :: temp_value(:)
end type

! global variables -----------------------------
integer, protected :: num_of_regions = 1, num_load_sets = 0, num_disp_sets = 0, total_load_sets = 0,total_disp_sets = 0
logical :: is_coarse_problem = .FALSE.

integer, parameter :: max_num_elements_one_node = 8 ! max number of elements sharing one node
real, parameter :: Numerical_Zero = 1.0e-10
real, parameter :: PI = 3.141592653589793238462643383279502884197169399375
! ----------------------------------------------
type (pdomain_type), allocatable, protected :: region(:)
type (pdomain_mat_db_type), protected :: pdomain_mat_db

type (load_set_type), allocatable, protected :: load_set(:)
type (disp_set_type), allocatable, protected :: disp_set(:), disp_set_backup(:)
real, allocatable, protected :: load_factor_db(:,:)

type (temp_load_type), protected :: temperature 


contains
!==============================================================================

subroutine create_domain(number)

implicit none

integer, intent(in) :: number

allocate(region(number))   
num_of_regions = number

end subroutine create_domain





subroutine set_domain(n, max_node_dof, max_dim, max_nodes, max_elements, nodes_elem, num_quad_pts, ngroup)

implicit none

integer, intent(in) :: n, max_node_dof, max_dim, max_nodes, max_elements, nodes_elem, num_quad_pts, ngroup
integer :: i


if (max_node_dof < 1) STOP 'set_domain: range error in max_node_dof!'
if (max_dim > 3 .or. max_dim < 2) STOP 'set_domain: range error in max_dim!'

i = max_node_dof + 3  ! node#, # elements sharing the node, flag if it is actually used by any element

allocate(region(n)%node(i,max_nodes))
allocate(region(n)%node_coord(max_dim,max_nodes))
allocate(region(n)%node_occupancy(max_nodes))
allocate(region(n)%skin_node_map(max_nodes))
region(n)%node = 0
region(n)%node_occupancy = .FALSE.
region(n)%skin_node_map = 0

if (max_elements > 0) then
  i = 3 + nodes_elem + num_quad_pts   ! num_quad_pts: number of material history address data
  allocate(region(n)%element(i,max_elements))
  allocate(region(n)%element_area(max_elements))
  allocate(region(n)%element_density(max_elements))
  allocate(region(n)%consis_mass_ratio(max_elements))
  allocate(region(n)%damping_factor(max_elements,2))
  allocate(region(n)%element_participation(max_elements))
  region(n)%element_participation = .TRUE.
  region(n)%consis_mass_ratio = 0.0
  if (max_dim == 3) then
    region(n)%size_str_output = 12
  else
    region(n)%size_str_output = 10
  endif
  allocate(region(n)%element_result(num_quad_pts,region(n)%size_str_output,max_elements))
elseif (max_elements < 0) then
  STOP  'set_domain: max_elements cannot be negative!'
endif

allocate(region(n)%node_group(ngroup))
allocate(region(n)%element_group(ngroup))
allocate(region(n)%skin_group_pointer(ngroup))


region(n)%num_node_dof = max_node_dof
region(n)%num_dim = max_dim
region(n)%num_nodes = max_nodes
region(n)%num_elements = max_elements
region(n)%nodes_element = nodes_elem
region(n)%num_quad_pts = num_quad_pts
region(n)%num_groups = ngroup

region(n)%max_bandwidth = 0
region(n)%node_group = 1
region(n)%element_group = 1
region(n)%skin_group_pointer = 1

region(n)%num_springs = 0

if (max_dim > 2) then ! 3D space region - Hexa elements with 4-node face
  allocate(region(n)%element_face(6,6,max_elements))
  allocate(region(n)%element_face_participation(6,max_elements))
  region(n)%element_face_participation = .TRUE. ! default: external surface
endif


end subroutine set_domain





subroutine set_spring_elements(n, max_springs)

! 2-node spring (link) element

implicit none

integer, intent(in) :: n, max_springs
integer :: i

if (max_springs > 0) then
  i = 3 + 2 + 1    ! nodes_elem = 2, # of material history address data = 1
  allocate(region(n)%spring(i,max_springs))
  allocate(region(n)%spring_area(max_springs))

  allocate(region(n)%spring_participation(max_springs))
  region(n)%spring_participation = .TRUE.

  region(n)%num_springs = max_springs
endif

end subroutine set_spring_elements





subroutine create_boundary_set(nl, nd, num_contacts)

implicit none

integer, intent(in) :: nl, nd, num_contacts
integer :: i

num_load_sets = nl
num_disp_sets = nd
total_load_sets = num_load_sets + num_contacts*2
total_disp_sets = num_disp_sets + num_contacts*2

allocate(load_set(total_load_sets))
allocate(disp_set(total_disp_sets))
allocate(disp_set_backup(total_disp_sets))

do i = 1, num_of_regions
  allocate(region(i)%skin_load_set(total_load_sets))
  region(i)%skin_load_set = .FALSE.
end do

end subroutine create_boundary_set






subroutine write_str_output(rn, elem_num, nquad, str_out, num_output)

implicit none

integer, intent(in) :: rn, elem_num, nquad, num_output
real, intent(in) :: str_out(nquad,num_output)


! for stress/strain output
region(rn)%element_result(1:nquad,1:num_output,elem_num) = str_out(:,:)

end subroutine write_str_output






subroutine set_material_db(num_materials,mat_num,model_num)

implicit none

integer, intent(in) :: num_materials,mat_num,model_num


if (.not.(allocated(pdomain_mat_db%mat_num_model_num))) then
  allocate(pdomain_mat_db%mat_num_model_num(num_materials))
  pdomain_mat_db%num_materials = num_materials
endif

pdomain_mat_db%mat_num_model_num(mat_num) = model_num

end subroutine set_material_db







subroutine build_nodes(unit_num, rn)

implicit none

integer,intent(in) :: unit_num, rn
integer :: max_dim, number, counter, i, j, status, max_nodes, current_domain, max_user_number
real :: field(50)
character(5) :: keyword

do
  read(unit_num,*, iostat=status) keyword
  if (status == -1 .or. keyword == '*end ') then
    write(*,*) '*** EOF in node input file ! -----------------------------'
    EXIT   ! exit from do-while loop
  elseif (keyword == '*node') then
    backspace(unit_num)
    read(unit_num, *, iostat=status) keyword, current_domain
    if (rn == current_domain) then
      max_dim = region(rn)%num_dim
      max_nodes = region(rn)%num_nodes
      counter = 0

      do j = 1, max_nodes
        read(unit_num, *, iostat=status) number, (field(i), i=1,max_dim)
        if (status == -1) EXIT

        counter = counter + 1

        region(rn)%node(1,counter) = number ! number: user node number / counter: internal node number
        region(rn)%node(2,counter) = 0  ! node(2,:): dof bandwidth of the node
        region(rn)%node(3,counter) = 0  ! node(3,:): number of 2D/3D solid elements using the node
        region(rn)%node_coord(1:max_dim,counter) = field(1:max_dim) !nodal coordinates

      end do
      write(*,*) 'Domain Region #:', current_domain, '     Number of nodes:', counter
    endif
  endif
end do


! internal_user_node_number_table (previously subroutine)
max_user_number = MAXVAL(region(rn)%node(1,:))
allocate(region(rn)%node_number_table(max_user_number))
region(rn)%node_number_table = -123

do i = 1, region(rn)%num_nodes
  region(rn)%node_number_table(region(rn)%node(1,i)) = i
end do


end subroutine build_nodes
	






subroutine build_elements(unit_num, rn, elem_type, mat_num, area, density, cfac, damping, mat_counter, counter, elem_flag)

! Last modification: July, 2015

! elem_flag  = .TRUE.  : 2D/3D elements
!            = .FALSE. : spring elements

implicit none

integer, intent(in) :: unit_num, rn, elem_type, mat_num
real, intent(in) :: area, density, cfac, damping(:)
logical, intent(in) :: elem_flag
integer, intent(inout) :: mat_counter, counter
character(len=5) :: keyword
integer, allocatable :: connectivity(:)
integer :: i, j, n, number, num_elem, status


if (elem_flag) then	
  num_elem = region(rn)%num_elements
  n = region(rn)%nodes_element
else
  num_elem = region(rn)%num_springs
  n = 2
endif

allocate(connectivity(n))

do
  read(unit_num,*, iostat=status) keyword
  backspace(unit_num)
  if (status == -1 .OR. keyword == '*elem' .OR. keyword == '*spri' .OR. keyword == '*end ') then
    EXIT   ! exit from do-while loop
  else

    read(unit_num,*, iostat=status) number, (connectivity(i), i=1,n)

    if (elem_flag) then
      counter = counter + 1
      mat_counter = mat_counter + 1
      call write_elem(rn,counter,number,elem_type,mat_num,mat_counter,area,density,cfac,damping,connectivity,n)
    else
      counter = counter + 1
      mat_counter = mat_counter + 1
            
      call write_spring(rn,counter,number,elem_type,mat_num,mat_counter,area,connectivity)
    endif
    call assign_dof_element(rn,connectivity,n,elem_flag)  ! element-by-element
  endif
end do			
	
deallocate(connectivity)
write(*,*) '  Number of DOF:', region(rn)%num_dofs

end subroutine build_elements 





subroutine write_elem(rn,p,number,type_num,mat_num,mat_counter,area,density,cfac,damping,connectivity,ndim)

! Last modification: July, 2015

! rn: region number
! p: internal element number
! number: user element number
! type_num: element type number
! mat_num: material number
! area: sectional area for 1D, thickness for 2D elements, or a physical property
! density
! cfac: consistence mass factor
! connectivity(1:ndim): element connectivity in user(in)/internal(out) node number

implicit none

integer, intent(in) :: rn, p, number, type_num, mat_num, ndim
real, intent(in) :: area, density, cfac, damping(:)
integer, intent(inout) :: connectivity(ndim), mat_counter
integer :: n_dim, i, j, k, ne, nnd, nq, elem(ndim)
integer :: ix_face(6,4), ix1(4), ix2(7), start, node1, node2, p_node1, p_node2
integer, allocatable :: v(:)


n_dim = region(rn)%num_dim
ne = region(rn)%nodes_element
! nnd = region(rn)%num_nodes
nq = region(rn)%num_quad_pts

region(rn)%element(1,p) = number
region(rn)%element(2,p) = type_num
region(rn)%element(3,p) = mat_num
region(rn)%element_area(p) = area
region(rn)%element_density(p) = density
region(rn)%consis_mass_ratio(p) = cfac
region(rn)%damping_factor(p,:) = damping(:)

elem = connectivity

! write element connectivity: using internal node numbers
! allocate(v(nnd))
! v = region(rn)%node(1,:)
! call find_match(v, nnd, elem, connectivity, ne)
connectivity = region(rn)%node_number_table(elem)

region(rn)%element(4:ne+3,p) = connectivity

do i = 1, nq
   region(rn)%element(ne+3+i,p) = mat_counter - 1 + i
end do

mat_counter = mat_counter - 1 + nq


if (n_dim > 2) then ! 3D case - Hexa only
  if (ne /= 8) STOP 'write_elem: for 3D hexa element only!'
  ix_face(1,:) = [1, 4, 3, 2]
  ix_face(2,:) = [5, 6, 7, 8]
  ix_face(3,:) = [1, 2, 6, 5]
  ix_face(4,:) = [2, 3, 7, 6]
  ix_face(5,:) = [3, 4, 8, 7]
  ix_face(6,:) = [4, 1, 5, 8]

  do i = 1, 6 ! over the current face
    ix2(1:4) = connectivity(ix_face(i,:))
    ix2(5:7) = ix2(1:3) ! make one cycle
    start = MINLOC(ix2(1:4),1)
    ix1 = ix2(start:(start+3))
    region(rn)%element_face(1:4,i,p) = ix1
    region(rn)%element_face(5:6,i,p) = 0   ! partner element & face info (0: no partner yet)
    node1 = ix1(1)
    node2 = ix1(3)
    match: do j = 1, p-1 ! over element
      do k = 1, 6   ! over face
        if (region(rn)%element_face_participation(k,j)) then
          p_node1 = region(rn)%element_face(1,k,j)
          p_node2 = region(rn)%element_face(3,k,j)
          if ((node1 == p_node1)) then
            if ((node2 == p_node2)) then  ! bingo - found the matched face!
              region(rn)%element_face(5,k,j) = p  ! record the current info into the partner
              region(rn)%element_face(6,k,j) = i

              region(rn)%element_face(5,i,p) = j  ! record the partner info into the current
              region(rn)%element_face(6,i,p) = k

              region(rn)%element_face_participation(k,j) = .FALSE.  ! internal face - partner
              region(rn)%element_face_participation(i,p) = .FALSE.  ! internal face - current

              region(rn)%num_inner_faces = region(rn)%num_inner_faces + 2
              EXIT match
            endif
          endif
        endif
      end do  ! k
    end do match ! j

  end do  ! i
endif ! 3D case

! deallocate(v)

end subroutine write_elem






subroutine write_spring(rn,p,number,type_num,mat_num,mat_counter,area,connectivity)

! rn: region number
! p: internal spring element number
! number: user spring element number
! type_num: element type number
! mat_num: material number
! area: sectional area for 1D, thickness for 2D elements, or a physical property
! connectivity(1:ndim): element connectivity in user(in)/internal(out) node number

implicit none

integer, intent(in) :: rn, p, number, type_num, mat_num
real, intent(in) :: area
integer, intent(inout) :: connectivity(2), mat_counter
integer :: elem(2)


region(rn)%spring(1,p) = number
region(rn)%spring(2,p) = type_num
region(rn)%spring(3,p) = mat_num
region(rn)%spring_area(p) = area

elem = connectivity

connectivity = region(rn)%node_number_table(elem)

region(rn)%spring(4:5,p) = connectivity
region(rn)%spring(6,p) = mat_counter


end subroutine write_spring






subroutine read_boundary(unit_num, rn)

implicit none

integer,intent(in) :: unit_num, rn
integer :: i, j, k, n, nn, nsize, sbound_type, sboundary(3)
integer :: number, status, ndof, ndim, nelem, n_nodes, group_number, n_group, plane_num
integer,allocatable :: boundary(:,:), constr_node(:), found(:)
integer, allocatable :: positions(:,:)
real :: degree
real,allocatable :: point(:,:), points(:,:,:)
character(5) :: keyword
logical, allocatable :: constr(:), bound_info(:,:), check(:)


nn = region(rn)%num_nodes
ndof = region(rn)%num_node_dof
ndim = region(rn)%num_dim
nelem = region(rn)%num_elements

allocate(constr(nn))
allocate(bound_info(ndof, nn)) 

constr = .FALSE.
bound_info = .FALSE.


do   ! do-while in a domain
	read(unit_num, *, iostat=status) keyword

	select case (keyword)
	case ('*cons')  ! 2D & 3D
		backspace(unit_num)
		read(unit_num,*) keyword, number
		allocate(boundary(ndof,number), constr_node(number), found(number))
		boundary = 0

		do i = 1, number
			read (unit_num, *) constr_node(i), (boundary(j,i), j= 1, ndof)
		end do
		
    found = region(rn)%node_number_table(constr_node)
    
		do i = 1, number
			constr(found(i)) = .TRUE.
			do j = 1, ndof
				k = boundary(j,i)
				if (k /= 0) bound_info(j, found(i)) = .TRUE.
			end do  ! j
!      write(*,*) '*cons:', found(i)

		end do  ! i

		deallocate(boundary, constr_node, found)


	case ('*cbou')  ! 2D boundary along the line defined by two-end points
		backspace(unit_num)
		read(unit_num,*) keyword, number
		allocate(boundary(ndof,number))
		allocate(point(4,number))
		boundary = 0

		do i = 1, number
			read (unit_num,*) (point(j,i), j=1,4), (boundary(j,i), j=1,ndof) ! 2D only
		end do

		call set_cbound(rn, number, point, boundary, nn, ndof, constr, bound_info)

		deallocate(boundary, point)


	case ('*pbou')  ! 2D surface boundary based on 2 skin nodal coordiantes: point(1:2,:) & point(3:4,:)
		backspace(unit_num)
		read(unit_num,*) keyword, group_number, number
		allocate(boundary(ndof,number))
		allocate(point(4,number))
		boundary = 0

		do i = 1, number
			read (unit_num,*) (point(j,i), j=1,4), (boundary(j,i), j=1,ndof)
		end do

		call set_pbound(rn, group_number, number, point, boundary, nn, ndof, constr, bound_info)

		deallocate(boundary, point)


	case ('*nbou') ! 2D surface boundary based on 2 skin nodal numbers
		backspace(unit_num)
		read(unit_num,*) keyword, group_number, number
		allocate(boundary(ndof,number))
		allocate(positions(2,number))
		boundary = 0

		do i = 1, number
			read (unit_num,*) (positions(j,i), j=1,2), (boundary(j,i), j= 1,ndof)
		end do

		call set_nbound(rn, group_number, number, positions, boundary, nn, ndof, constr, bound_info)

		deallocate(boundary, positions)


	case ('*sbou')  ! 3D surface boundary based on 4( or 3) points rectangular palne
		backspace(unit_num)
    sboundary = 0

		read(unit_num,*) keyword, sbound_type, number, (sboundary(j), j= 1,ndof)

    if (ndim /= 3) STOP 'read_boundary: sbou(nd) is only for 3D!'

    if ((sbound_type == 0) .AND. (number < 3)) then
      STOP 'read_boundary: sbound_type == 0 needs 3 points!'
    elseif ((sbound_type /= 0) .AND. (number < 4)) then
      STOP 'read_boundary: sbound_type == 1 needs 4 points (CCW)!'
    endif

		allocate(point(ndim,number))

    do i = 1, number
      read (unit_num,*) (point(j,i), j=1,ndim)
    end do

		call set_sbound(rn, sbound_type, number, point, sboundary, nn, ndof, constr, bound_info)

		deallocate(point)


	case ('*mpcn')  ! multi-point constraints: 2D & 3D
		backspace(unit_num)
		read(unit_num,*) keyword, n_group
    write(*,'(/A,I5)') ' read MPC data: number of MPC groups =', n_group
    region(rn)%num_mpc_group = n_group  ! number of MPC groups

    allocate(region(rn)%mpc_plane(n_group)) ! MPC principal plane number
    allocate(region(rn)%mpc_theta(n_group))
    allocate(check(n_group))
    check = .FALSE.
    allocate(points(ndim,3,n_group))

    do n = 1, n_group
      read(unit_num,*) group_number, plane_num, degree, number
      if (group_number > n_group) STOP 'read_boundary: MPC group number is larger than number of MPC groups!'
      if ((plane_num > 3) .OR. (plane_num < 1)) STOP 'read_boundary: MPC plane number is out of range [1,3]!'
      if (number < 2) STOP 'read_boundary: MPC needs at least 2 points to define inclined line or plane!'
      if ((ndim > 2) .AND. (number /= 3)) then
        write(*,*) '   Number of MPC line or plane points =', number, '  MPC group number =', group_number
        STOP 'set_mpc: MPC 3D plane needs 3 points!'
      endif

      do i = 1, number
        read (unit_num,*) (points(j,i,group_number), j=1,ndim)
      end do

      region(rn)%mpc_plane(group_number) = plane_num
      region(rn)%mpc_theta(group_number) = (PI/180.0)*degree  ! transform from degree (user input) to radian
      check(group_number) = .TRUE.
    end do ! n: read over MPC group

    if (.NOT. ALL(check)) STOP 'read_boundary: MPC must be prescribed for all MPC groups!'

    call set_mpc(rn, ndim, ndim, points)  ! must be called only one time in the same region

    deallocate(points, check)

	case default
		backspace(unit_num)
		EXIT  !------------------------------------------------------------>>>
	end select

end do  ! do-while in a domain

! write(*,'(A/,(3L5))') '   >> contraints', ((bound_info(i,j), i=1,ndof), j=1,nn)
! assign DOF numbers for constrained nodes ---------------------------------------------------------

do j = 1, nn    ! loop over all nodes (maybe inefficient!)
	if (constr(j)) then ! for constrained node only
		do k = 1, ndof
			if (bound_info(k,j)) then ! for constrained DOF at this constrained node
				region(rn)%node(3+k, j) = -region(rn)%node(3+k, j)   ! temporally negative DOF numbers for inactive DOF
			endif
		end do
	endif
end do

deallocate(constr, bound_info)


end subroutine read_boundary






subroutine set_pbound(rn, group, num, axis, cnst, nn, ndof, constr, bound)

! apply constraints along a boundary surface using the coordinates of two given nodal points
! 2D
! Input:
!       num: number of nodal sets
!       axis(1:2,:): 1st nodal coordinates in the plane
!       axis(3:4,:): 2nd nodal coordinates in the plane
!       cnst(:,:): constraints (=1: fixed)

implicit none

integer,intent(in) :: rn, group, ndof, num, nn 
integer,intent(in) :: cnst(ndof,num)
real,intent(in) :: axis(4,num)
logical, intent(inout) :: constr(nn), bound(ndof,nn)

integer :: k, ndim, num_bound, st, ed, counter, found(2)
integer, allocatable :: bound_node(:), skin_node(:)
real :: cons_coord(2,2), v1(2), distance(2), TOL, DNRM2
real, allocatable :: coord(:,:)


num_bound = region(rn)%num_skin	
allocate(bound_node(num_bound)) 
allocate(skin_node(num_bound)) 
allocate(coord(2, num_bound)) 
coord = region(rn)%node_coord(:,region(rn)%skin_nodes(:))

counter = 0

do k = 1, num 
	cons_coord(:,1) = axis(1:2, k)
	cons_coord(:,2) = axis(3:4, k)

  v1 = cons_coord(:,2) - cons_coord(:,1)
  TOL = 1.0e-3*DNRM2(2,v1,1)

	call find_match_coord(2, cons_coord, 2, coord, num_bound, found, distance)
	
	if ((distance(1) < TOL ).AND. (distance(2) < TOL)) then
    st = region(rn)%skin_nodes(found(1))
    ed = region(rn)%skin_nodes(found(2))

    call find_skin(rn, group, st, ed, counter, bound_node, skin_node, 0) 
	
    if (cnst(1,k) == 1) bound(1, bound_node(1:counter)) = .TRUE.
    if (cnst(2,k) == 1) bound(2, bound_node(1:counter)) = .TRUE.
    constr(bound_node(1:counter)) = .TRUE.
	endif

end do  ! k

deallocate(bound_node, skin_node, coord) 
 
end subroutine set_pbound







subroutine set_nbound(rn, group, num, point, cnst, nn, ndof, constr, bound)

! apply constraints along a boundary surface using numbering data of two given nodal points
! 2D
! Input:
!       num: number of nodal sets
!       point(1:2,:): starting & end nodal numbers
!       cnst(:,:): constraints (=1: fixed)

implicit none

integer,intent(in) :: rn, group, nn, ndof, num, point(2,num), cnst(ndof,num)
logical,intent(inout) :: constr(nn), bound(ndof,nn)

integer :: k, num_bound, counter, node_numbers(2)
integer, allocatable :: bound_node(:), skin_node(:)


num_bound = region(rn)%num_skin
allocate(bound_node(num_bound))
allocate(skin_node(num_bound)) 

counter = 0

do k = 1, num

  node_numbers = region(rn)%node_number_table(point(1:2,k))

	call find_skin(rn, group, node_numbers(1), node_numbers(2), counter, bound_node, skin_node, 0)

	constr(bound_node(1:counter)) = .TRUE.
	if (cnst(1,k) == 1) bound(1, bound_node(1:counter)) = .TRUE.
	if (cnst(2,k) == 1) bound(2, bound_node(1:counter)) = .TRUE.

end do  ! k

deallocate(bound_node, skin_node) 

end subroutine set_nbound






subroutine set_cbound(rn, num, points, cnst, nn, ndof, constr, bound)

! new version: Jan. 2015
! 2D
! search part is an embedded version of subroutine 'search_nodes_along_line'
! can handle multiple constraint-lines
!
! Input:
!       num: number of constraint-lines
!       nn: number of all nodes in region 'rn'

implicit none

integer, intent(in) :: rn , nn, ndof, num
integer, intent(in) :: cnst(ndof,num) 
real, intent(in) :: points(4,num)
logical, intent(inout) :: constr(nn), bound(ndof,nn)

integer ::  i, j, k
real :: test_point(2), vec1(2), vec2(2), v1(2), v2(2)
real :: theta, norm1, norm2, DNRM2
logical :: on_points = .FALSE.
real, parameter :: TOL = 1.0e-5


do i = 1, num
  vec1 = points(1:2, i)
  vec2 = points(3:4, i)

  do j = 1, nn  ! loop over all nodes --------------------------

    test_point = region(rn)%node_coord(1:2,j)
    v1 = vec1 - test_point
    v2 = vec2 - test_point
    norm1 = DNRM2(2,v1,1)
    norm2 = DNRM2(2,v2,1)

    on_points = .FALSE.
    if ((norm1 <= TOL*norm2) .OR. (norm2 <= TOL*norm1)) then
      on_points = .TRUE.
    else
      call vector_angle(v1, v2, 2, theta)
      if (ABS(theta - PI) < TOL) on_points = .TRUE.
    endif

    if (on_points) then
      constr(j) = .TRUE.
      do k = 1, ndof
        if (cnst(k,i) == 1) then
          bound(k,j) = .TRUE.
        endif
      end do ! k
    endif

  end do ! j ----------------------------------------------------
  
end do ! i

end subroutine set_cbound






subroutine set_sbound(rn, sbound_type, num, points, cnst, nn, ndof, constr, bound)

! apply constraints over boundary surface faces using 3~4 points
! 3D
! Input:
!       num: number of nodal sets
!       point(1:2,:): starting & end nodal numbers
!       cnst(:,:): constraints (=1: fixed)

implicit none

integer,intent(in) :: rn, sbound_type, num, nn, ndof, cnst(ndof)
real, intent(in) :: points(3,4)
logical,intent(inout) :: constr(nn), bound(ndof,nn)

integer :: i, j, k, num_surfaces, face_number, elem_number
integer, allocatable :: face_ix(:,:)
real :: corners(3,4)
logical, allocatable :: bound_face(:)
real, parameter :: angle = 90


num_surfaces = region(rn)%num_skin_surfaces

if (sbound_type == 0) then
  corners(:,1) = points(:,1)
  corners(:,2) = points(:,3)
  corners(:,3) = points(:,2) + points(:,3) - points(:,1)
  corners(:,4) = points(:,2)

else
  corners = points(1:3,1:4)
endif

allocate(face_ix(4,num_surfaces))
allocate(bound_face(num_surfaces))

do i = 1, num_surfaces
  face_number = region(rn)%skin_surfaces(1,i)
  elem_number = region(rn)%skin_surfaces(2,i)
  face_ix(:,i) = region(rn)%element_face(1:4,face_number,elem_number)
end do

call search_contact_faces(sbound_type, num_surfaces, face_ix, region(rn)%num_nodes, region(rn)%node_coord, corners, angle, bound_face)

write(*,'(/A)') '*sbound:'

do i = 1, num_surfaces
  if (bound_face(i)) then
    constr(face_ix(1:4,i)) = .TRUE.
    if (cnst(1) == 1) bound(1, face_ix(1:4,i)) = .TRUE.
    if (cnst(2) == 1) bound(2, face_ix(1:4,i)) = .TRUE.
    if (cnst(3) == 1) bound(3, face_ix(1:4,i)) = .TRUE.
    do j = 1, 4
      write(*,'(A,2I9,3L3)') '   constrained node(internal & user), dof: ', face_ix(j,i), region(rn)%node(1,face_ix(j,i)), (bound(k, face_ix(j,i)), k=1,3)
    end do
  endif
end do  ! i

write(*,'(A)') '*sbound set done!'
deallocate(bound_face, face_ix)

end subroutine set_sbound






subroutine set_mpc(rn, ndim, npoint, points)

implicit none

integer, intent(in) :: rn, ndim, npoint
real, intent(in) :: points(ndim,npoint,*)

integer :: n, nn, ne, nelem, i, j, nsize, addr, old_addr
integer :: elem_num, face_num, nface, ngroup, group_num
integer, allocatable :: mpc_nodes(:), face_ix(:,:), mpc_group(:), mpc_group_packed(:)
real :: corners(3,4)
real, allocatable :: node_coord(:,:)
logical, allocatable :: found_flag(:), checkboard(:,:), in_or_out(:)
real, parameter :: angle = 90


ne = region(rn)%nodes_element
nelem = region(rn)%num_elements
ngroup = region(rn)%num_mpc_group

if (ndim == 2) then
  nn = region(rn)%num_skin
  allocate(found_flag(nn))
  allocate(mpc_group(nn))
  allocate(node_coord(ndim,nn))
  node_coord(:,:) = region(rn)%node_coord(:,region(rn)%skin_nodes)  ! MPC is checked only over skin nodes

  call search_nodes_along_line(ngroup, ndim, points, nn, node_coord, found_flag, mpc_group)

  n = COUNT(found_flag)

  if (n > 0) then

    region(rn)%mpc_existence = .TRUE.

    allocate(mpc_nodes(n))
    allocate(mpc_group_packed(n))
    allocate(region(rn)%mpc_element_addr(nelem+1))
    allocate(region(rn)%mpc_element_group(nelem))
    region(rn)%mpc_element_group = 0

    mpc_nodes = PACK(region(rn)%skin_nodes, found_flag)
    mpc_group_packed = PACK(mpc_group, found_flag)

    write(*,'(/A,I5)') ' > number of MPC groups =', ngroup
    write(*,*) '   number of nodes along MPC line =', n
    write(*,'(A/,(8I10/))') '    nodes along MPC line : ', mpc_nodes


    allocate(checkboard(n,nelem))
    checkboard = .FALSE.

    do i = 1, n ! over MPC nodes
      do j = 1, nelem ! over all elements
        checkboard(i,j) = ANY(region(rn)%element(4:3+ne,j) == mpc_nodes(i))
      end do ! j
    end do ! i

    nsize = COUNT(checkboard)
    allocate(region(rn)%mpc_element_nodes(nsize))

    region(rn)%mpc_element_addr(1) = 1
    addr = 1
    old_addr = 1
    do i = 1, nelem
      do j = 1, n
        if (checkboard(j,i)) then
          region(rn)%mpc_element_nodes(addr) = mpc_nodes(j)
          region(rn)%mpc_element_group(i) = mpc_group_packed(j) ! can be overwritten by later coming node info
          addr = addr + 1
        endif
      end do ! j
      region(rn)%mpc_element_addr(i+1) = addr
      old_addr = addr
    end do ! i

    deallocate(checkboard, mpc_nodes, mpc_group, mpc_group_packed)
    deallocate(node_coord, found_flag)

  else  ! found no node: fail to search nodes along the MPC line
    write(*,*) '   *** Found NO node: fail to search nodes along the MPC line!'
  endif

elseif (ndim == 3) then

  if (ngroup > 1) STOP 'set_mpc: presently only 1 MPC group for 3D!'
  group_num = 1  ! presently only one group for 3D

  nface = region(rn)%num_skin_surfaces
  allocate(face_ix(4,nface))
  allocate(found_flag(nface))

  corners(:,1) = points(:,1,group_num)
  corners(:,2) = points(:,3,group_num)
  corners(:,3) = points(:,2,group_num) + points(:,3,group_num) - points(:,1,group_num)
  corners(:,4) = points(:,2,group_num)
  write(*,'(A/4(3E11.3/))') '  MPC plane corners', corners

  do i = 1, nface
    elem_num = region(rn)%skin_surfaces(2,i)
    face_num = region(rn)%skin_surfaces(1,i)
    face_ix(1:4,i) = region(rn)%element_face(1:4,face_num,elem_num)
  end do

  call search_contact_faces(0, nface, face_ix, region(rn)%num_nodes, region(rn)%node_coord, corners, angle, found_flag)

  n = COUNT(found_flag)

  if (n > 0) then

    region(rn)%mpc_existence = .TRUE.

    allocate(region(rn)%mpc_element_addr(nelem+1))
    allocate(region(rn)%mpc_element_group(nelem))
    region(rn)%mpc_element_group = 0

    nsize = 4*n   ! 4 nodes per a face and only one face is under MPC per one element
    allocate(region(rn)%mpc_element_nodes(nsize))

    write(*,'(/A,I5)') ' > number of MPC groups =', ngroup
    write(*,'(A,I10/)') '   number of faces on MPC plnae =', n


    region(rn)%mpc_element_addr(1) = 1
    addr = 1
    old_addr = 1
    do i = 1, nelem

      do j = 1, nface
        if (found_flag(j)) then
          elem_num = region(rn)%skin_surfaces(2,j)
          if (i == elem_num) then
            face_num = region(rn)%skin_surfaces(1,j)
            region(rn)%mpc_element_nodes(addr:(addr+3)) = region(rn)%element_face(1:4,face_num,elem_num)
            addr = addr + 4
          endif
        endif
      end do ! j
      region(rn)%mpc_element_addr(i+1) = addr
      if (addr > old_addr) then
        region(rn)%mpc_element_group(i) = group_num  ! presently only one group for 3D
        write(*,'(A,3I10)') '   element, face, MPC group on MPC plane: ', elem_num, face_num, group_num
      endif
      old_addr = addr
    end do ! i
    
    write(*,'(A,I10/)') '   MPC address length =', addr

  else  ! found no node: fail to search faced on the MPC plane
    write(*,*) '   *** Found NO face: fail to search faces on the MPC Plane!'
  endif

  deallocate(face_ix, found_flag)

else
  STOP 'set_mpc: space dimension must be 2 or 3!'
endif


end subroutine set_mpc







subroutine set_disp_bc(sn, rn)

! sn: disp load set number
! rn: region number

implicit none

integer, intent(in) :: sn, rn
integer :: m1, m2, size


m1 = region(rn)%num_dofs      ! # of active dof
m2 = region(rn)%num_constr    ! # of constrained dof
size = m1 + m2

if (allocated(disp_set(sn)%constr)) then
  deallocate(disp_set(sn)%constr)
  deallocate(disp_set_backup(sn)%constr)
endif
allocate(disp_set(sn)%constr(size))
allocate(disp_set_backup(sn)%constr(size))
disp_set(sn)%constr = 0.0
disp_set(sn)%region_info = rn

end subroutine set_disp_bc






subroutine write_disp_bc(sn, rn, node_dof, disp, ndim)

! sn: disp load set number
! rn: region number
! node_dof(2,ndim): user node # (1,:) & DOF # (2,:)
! disp(ndim): prescribed displacement
! user node number is transformed to internal node number

implicit none

integer, intent(in) :: sn, rn, node_dof(2,ndim), ndim
real, intent(in) :: disp(ndim)
integer :: i, m1, node, dof, inactive_dof(ndim)


m1 = region(rn)%num_dofs

do i = 1, ndim
  node = region(rn)%node_number_table(node_dof(1,i))  ! transform to internal node number
  dof = node_dof(2,i)
  inactive_dof(i) = region(rn)%node(3+dof,node)
end do

disp_set(sn)%constr(inactive_dof) = disp

write(*,*) 'write disp dof:', inactive_dof

end subroutine write_disp_bc






subroutine update_disp_bc(sn, rn, ix, disp, ndim, flag)

! similar to 'write_disp_bc' but dof number input
! 2 dof per node
! sn: disp load set number
! rn: region number
! disp(ndim): prescribed displacement
! flag = 1: ix is nominal/user dof numbering
!      = 2: ix is internal system dof numbering

implicit none

integer, intent(in) :: sn, rn, ix(ndim), ndim, flag
real, intent(in) :: disp(ndim)
integer :: i, m1, node, dof, inactive_dof(ndim)

if (flag == 1) then
  disp_set(sn)%constr(ix) = disp
else
  inactive_dof = ix - region(rn)%num_dofs
  inactive_dof = region(rn)%disp_bc(inactive_dof)
  disp_set(sn)%constr(inactive_dof) = disp
endif  
  
end subroutine update_disp_bc






subroutine backup_disp_bc(flag)

! similar to 'write_disp_bc' but dof number input
! sn: disp load set number
! rn: region number
! flag =  1: back-up
!      = -1: recover

implicit none

integer, intent(in) :: flag
integer :: i

if (flag == 1) then
  do i = 1, total_disp_sets
    if (allocated(disp_set(i)%constr)) disp_set_backup(i)%constr = disp_set(i)%constr
  end do
elseif (flag == -1) then
  do i = 1, total_disp_sets
    if (allocated(disp_set(i)%constr)) disp_set(i)%constr = disp_set_backup(i)%constr
  end do
endif  
  
end subroutine backup_disp_bc






subroutine assign_dof_element(rn, connectivity, ndim, elem_flag)

! element-by-element dof numbers assigning
! assign DOF numbers except constrained nodes, which are already dof-asigned
! ndim: number of nodes for this element
! elem_flag  = .TRUE.  : quad elements
!            = .FALSE. : spring elements

implicit none

integer, intent(in) :: rn, ndim, connectivity(ndim)
logical, intent(in) :: elem_flag
integer :: i, j, n_node_dof, node_num, dof_number
logical :: virgin

dof_number = region(rn)%num_dofs     ! initially zero

n_node_dof = region(rn)%num_node_dof

do i = 1, ndim
	node_num = connectivity(i)
	virgin = .NOT. region(rn)%node_occupancy(node_num)
	if (virgin) then
		do  j = 4, n_node_dof+3
      dof_number = dof_number + 1
      region(rn)%node(j,node_num) = dof_number ! Nominal DOF #
		end do
    region(rn)%node_occupancy(node_num) = .TRUE.  ! DOF-assigned node
	endif		
  if (elem_flag) region(rn)%node(3,node_num) = region(rn)%node(3,node_num) + 1 ! node(3,:): number of 2D/3D solid elements using the node

  region(rn)%node(2,node_num) = region(rn)%node(2,node_num) + ndim*n_node_dof ! dof bandwidth of the node
end do

region(rn)%num_dofs = dof_number 


end subroutine assign_dof_element






subroutine assign_global_dof(rn)

! called from subroutine 'build_model'

implicit none

integer, intent(in) :: rn
integer :: i, j, k, counter, b_counter
integer :: nn, ndof, m1, m2, total_num_dof


nn = region(rn)%num_nodes
ndof = region(rn)%num_node_dof
allocate(region(rn)%dof_activity_flag(nn*ndof))

region(rn)%dof_activity_flag =.TRUE.
b_counter = 0

do i = 1, nn
	do j = 1, ndof
		k = region(rn)%node(3+j, i)   ! assigned negative DOF numbers in 'assign_dof_element'
		if (k < 0) then
			region(rn)%node(3+j, i) = -k 
      region(rn)%dof_activity_flag(-k) = .FALSE.
      b_counter = b_counter + 1
    elseif (k == 0) then
      write(*,*) '*** assign_global_dof: No DOF assigned for node:', i,' dof:', j
      STOP 
		endif
	end do
end do

region(rn)%num_constr = b_counter ! number of constrained DOFs 
m2 = region(rn)%num_constr
region(rn)%num_dofs = region(rn)%num_dofs - m2
m1 = region(rn)%num_dofs
total_num_dof = m1 + m2
write(*,*) '  number of active dofs = ', m1
write(*,*) '  number of constrained dofs = ', m2 

if ((nn*ndof) /= total_num_dof) STOP '*** assign_global_dof: Error in total_num_dof !'

! set active DOF index array (active_dof) & inactive DOF index array (disp_bc)

allocate(region(rn)%active_dof(m1))
allocate(region(rn)%disp_bc(m2))

allocate(region(rn)%dof_to_array(total_num_dof))
allocate(region(rn)%rigid_constraints(total_num_dof))

counter = 1
b_counter = 1
do i = 1, total_num_dof
  if (region(rn)%dof_activity_flag(i)) then  ! active dof
    region(rn)%active_dof(counter) = i
    region(rn)%dof_to_array(i) = counter
    counter = counter + 1
  else
    region(rn)%disp_bc(b_counter) = i
    region(rn)%dof_to_array(i) = b_counter + m1
    b_counter = b_counter + 1
  endif  
end do 

region(rn)%rigid_constraints = .NOT.(region(rn)%dof_activity_flag)

end subroutine assign_global_dof







subroutine dof_transformer(rn, flag, dof_numbers, ndim)

! total number of dofs is preserved!
! flag =  1: move to active dofs from constrained dofs
!      = -1: discharge from active dofs, & move to constrained dofs
! dof_numbers: nominal (physical) dof numbers to move

implicit none

integer, intent(in) :: rn, flag, dof_numbers(ndim), ndim
integer :: i, j, m1, m2, total_num_dof, counter, b_counter
integer, allocatable :: squad(:)


! allocate(squad(ndim))
! squad = region(rn)%dof_to_array(dof_numbers)

if (flag == 1) then   ! dofs move to active dof part
  counter = 0
  do i = 1, ndim
    j = dof_numbers(i)
    if ((.NOT.region(rn)%rigid_constraints(j)).AND.(region(rn)%dof_activity_flag(j) == .FALSE.)) then
      region(rn)%dof_activity_flag(j) = .TRUE.
      counter = counter + 1
    endif
  end do

  m1 = region(rn)%num_dofs + counter
  m2 = region(rn)%num_constr - counter
  total_num_dof = m1 + m2

  if (counter > 0) then    
    deallocate(region(rn)%active_dof, region(rn)%disp_bc)
    allocate(region(rn)%active_dof(m1))
    allocate(region(rn)%disp_bc(m2))
    region(rn)%dof_to_array = 0
  
! renumbering DOFs
    counter = 1
    b_counter = 1
    do i = 1, total_num_dof
      if (region(rn)%dof_activity_flag(i)) then  ! active dof
        region(rn)%active_dof(counter) = i
        region(rn)%dof_to_array(i) = counter
        counter = counter + 1
      else
        region(rn)%disp_bc(b_counter) = i
        region(rn)%dof_to_array(i) = b_counter + m1
        b_counter = b_counter + 1
      endif  
    end do 
    
    region(rn)%num_dofs = m1
    region(rn)%num_constr = m2
    call write_dof_changed(rn, .TRUE.)
  endif
  
elseif (flag == -1) then   ! dofs move to constrained dof part
 
  counter = 0
  do i = 1, ndim
    if (region(rn)%dof_activity_flag(dof_numbers(i)) == .FALSE.) then
      counter = counter + 1   ! already constrained dof !!!
    else
      region(rn)%dof_activity_flag(dof_numbers(i)) = .FALSE.  
    endif
  end do

  m1 = region(rn)%num_dofs - ndim + counter
  m2 = region(rn)%num_constr + ndim - counter
  total_num_dof = m1 + m2

  if (counter < ndim) then    
    deallocate(region(rn)%active_dof, region(rn)%disp_bc)
    allocate(region(rn)%active_dof(m1))
    allocate(region(rn)%disp_bc(m2))
    region(rn)%dof_to_array = 0
  
! renumbering DOFs
    counter = 1
    b_counter = 1
    do i = 1, total_num_dof
      if (region(rn)%dof_activity_flag(i)) then  ! active dof
        region(rn)%active_dof(counter) = i
        region(rn)%dof_to_array(i) = counter
        counter = counter + 1
      else
        region(rn)%disp_bc(b_counter) = i
        region(rn)%dof_to_array(i) = b_counter + m1   
        b_counter = b_counter + 1
      endif  
    end do   
  
    region(rn)%num_dofs = m1
    region(rn)%num_constr = m2
    call write_dof_changed(rn, .TRUE.)
  endif
  
else
  STOP 'dof_tranformer: invalid flag number !'
endif  

write(*,*) '> dof_transformer: region number =', rn
write(*,*) '     new numbers of active dof & constrained dof:', m1, m2
! write(*,*) 'dof_to_array:', region(rn)%dof_to_array
write(*,*)

end subroutine dof_transformer







subroutine write_dof_changed(rn, toggle)

! write the protected variable, dof_changed

implicit none

integer, intent(in) :: rn
logical, intent(in) :: toggle

region(rn)%dof_changed = toggle

end subroutine write_dof_changed







subroutine read_dof_number(rn, user_nodes, dofs, n_dofs, ndim)

! rn: region number
! user_nodes: array having user node numbers to look up the dof numbers
! dofs(nodal dof number,node number): result array
! ndim: number of nodes to read

implicit none

integer, intent(in) :: rn, user_nodes(ndim), n_dofs, ndim
integer, intent(out) :: dofs(n_dofs,ndim)
integer :: internal_nodes(ndim)


internal_nodes = region(rn)%node_number_table(user_nodes)   ! transform to internal node numbers

dofs(:,1:ndim) = region(rn)%node(4:,internal_nodes(1:ndim))

end subroutine read_dof_number










subroutine read_loads(unit_num, rn, flag_bc, flag_load, input_flag, gr_unit_numbers, gr_motion_flags)

implicit none

integer,intent(in) :: unit_num, rn, gr_unit_numbers(3)
logical,intent(in) :: flag_bc, flag_load
logical,intent(out) :: gr_motion_flags(3)
integer, intent(inout) :: input_flag

integer :: i, j, k, n, load_type_number
integer, allocatable :: load_point(:,:)

integer :: number, status, set_number, num_dim
real :: value, angle, ref_temperature
real, allocatable :: points(:,:), field(:), load_value(:), temp_time(:), temp_load(:)

character(5) :: keyword
character(len=20) :: gr_file_name


num_dim = region(rn)%num_dim


do   ! do-while in a domain region -------------------------------------------------------------
		
  read(unit_num, *, iostat=status) keyword

  if (status == -1 .OR. keyword == '*end ') then
    EXIT   ! exit from do-while loop
  elseif  (keyword == '*doma' .OR. keyword == '*mate') then
    backspace(unit_num)
    EXIT
  endif

  select case (keyword)
		
!   =======================================================

  case ('*disp') ! this input keyword cannot be repeated for the same disp load set number
    if (flag_bc) STOP 'build_model: input order error!'
    backspace(unit_num)
    read(unit_num,*) keyword, number, set_number
    if (set_number > num_disp_sets) STOP 'build_model: disp load set number error!'
      
    i = disp_set(set_number)%region_info
    if (i > 0)	STOP 'build_model:  this disp load set number is already defined!'
      
    call set_disp_bc(set_number, rn)
      
    allocate(load_point(2,number), load_value(number))
      
    n = region(rn)%num_node_dof
    do i = 1, number
      read (unit_num, *)  j, k, value
      if ((1 > k).OR.(n < k)) STOP 'build_model: invalid dof # in *disp input!'
      load_point(1,i) = j   ! node number (user node number)
      load_point(2,i) = k   ! DOF number
      load_value(i) = value   ! prescribed disp value
      value = 0.0
    end do
      
    call  write_disp_bc(set_number, rn, load_point, load_value, number)
      
    deallocate(load_point, load_value)
      
  case ('*load') ! this input keyword cannot be repeated for the same disp load set number
    if (input_flag < 2 .OR. flag_load) STOP 'build_model: input order error!'
    backspace(unit_num)
    read(unit_num,*) keyword, number, load_type_number, set_number
    if (set_number > num_load_sets) STOP 'build_model: load set number error!'

! load_type_number: physical_domain variable
! load_type_number == 2: axisymmetric (nonlinear)
!                  == 3: axisymmetric loads (linear)
    
    i = load_set(set_number)%region_info
    if (i > 0)	STOP 'build_model:  this load set number is already defined!'

    allocate(load_point(2,number), load_value(number))
    
    n = region(rn)%num_node_dof
    do i = 1, number
      read (unit_num,*) j, k, value
      if ((1 > k).or.(n < k)) STOP 'build_model: invalid dof # in *load input!'
      load_point(1,i) = j   ! node number (user node number)
      load_point(2,i) = k   ! DOF number
      load_value(i) = value   ! load value
      value = 0.0
    end do
    
    input_flag = 4

    call write_load_bc(set_number, rn, load_point, load_value, number, load_type_number)
      
    deallocate(load_point, load_value)
    
			
  case ('*pres') ! this input keyword cannot be repeated for the same disp load set number
    if (input_flag < 2 .OR. flag_load) STOP 'build_model:input order error!'
    backspace(unit_num)
    read(unit_num,*) keyword, number, load_type_number, set_number, angle, value
    if (set_number > num_load_sets) STOP 'build_model: load set number error!'
    if ((num_dim > 2) .AND. (number /= 4)) STOP 'buildmodel: 3D pressure plane needs 4 points(CCW)!'
    if ((num_dim > 2) .AND. ((load_type_number == 2).OR.(load_type_number == 3))) STOP 'buildmodel: load_type_number(2/3) is not for 3D pressure!'

! load_type_number: physical_domain variable      
! load_type_number == 2: axisymmetric (nonlinear)
!                  == 3: axisymmetric loads (linear)

!    if (load_set(set_number)%region_info > 0)	write(*,*) '***** warning: pressure load set number is already defined!'

    allocate(load_point(3,number), load_value(number), points(num_dim,number))
                    
    do i = 1, number
      read (unit_num,*) (points(j,i), j=1,num_dim)    ! coordinates of pressured domain
    end do

    input_flag = 4

    call write_surface_load(rn, num_dim, number, points, angle, value, load_type_number, set_number)

    deallocate(load_point, load_value, points)

  case ('*tloa') ! temperature load: this input keyword cannot be repeated in the whole input
    backspace(unit_num)
    read(unit_num,*) keyword, number, ref_temperature
    allocate(temp_time(number), temp_load(number))
    do i= 1, number
          read (unit_num,*) temp_time(i), temp_load(i)
    end do

    call write_temp_load(ref_temperature, number, temp_time, temp_load)

    deallocate(temp_time, temp_load)
  case ('*grou')

    backspace(unit_num)
    read(unit_num, *, iostat=status) keyword, i, value, gr_file_name
    
    select case(i)
    case (1)
      open(unit=gr_unit_numbers(1), file=gr_file_name, status='OLD', action='READ') 		
      call write_gr_motion(rn, 1, value)
      write(*,*) '>> FEMULA: open ground motion file: x1-direction:  ', gr_file_name
      gr_motion_flags(1) = .TRUE.
    case (2)
      if (num_dim < 2) STOP 'build_model: ground motion dimension(2) exceeds the max problem dimension!'        
      open(unit=gr_unit_numbers(2), file=gr_file_name, status='OLD', action='READ') 		
      call write_gr_motion(rn, 2, value)
      write(*,*) '>> FEMULA: open ground motion file: x2-direction:  ', gr_file_name
      gr_motion_flags(2) = .TRUE.
    case (3)
      if (num_dim < 3) STOP 'build_model: ground motion dimension(3) exceeds the max problem dimension!'
      open(unit=gr_unit_numbers(3), file=gr_file_name, status='OLD', action='READ') 		
      call write_gr_motion(rn, 3, value)	
      write(*,*) '>> FEMULA: open ground motion file: x3-direction:  ', gr_file_name
      gr_motion_flags(3) = .TRUE.
    end select  ! i

  end select
end do   ! do-while in a domain region ---------------------------------------------------------


end subroutine read_loads








subroutine write_load_bc(sn, rn, load_point, load_value, nline, load_type_number)

! sn: load set number
! rn: region number
! load_point(1:2,nline): load posiions (1: node #, 2: dof #)
! load_value(nline): load values
! load_type_number == 2: axisymmetric (nonlinear)
!                  == 3: axisymmetric loads (linear)

! user node numbers are transformed to internal node numbers

implicit none

integer, intent(in) :: rn, nline, load_type_number, sn
integer, intent(in) :: load_point(2,nline)
real, intent(in) :: load_value(nline)
integer :: i, k, counter
integer :: the_dof, load_nodes(nline), found(nline)
integer :: bc(nline), bc_node(nline)
real :: cn(nline)


load_set(sn)%region_info = rn
load_set(sn)%load_type = load_type_number

load_nodes = load_point(1,:)
found = region(rn)%node_number_table(load_nodes)    ! transform to internal node number

counter = 0
do i = 1, nline
  if (found(i) /= 0 ) then 
    k = load_point(2,i)   ! DOF number at one point
    the_dof = region(rn)%node(3+k,found(i))
    if (the_dof > 0) then
      counter = counter + 1
      bc(counter) = the_dof
      cn(counter) = load_value(i)
      bc_node(counter) = found(i)
    endif
  endif
end do

! if (counter < 1) STOP 'write_constraints: boundary condition counter is less than 1'

load_set(sn)%num_load = counter
load_set(sn)%region_info = rn ! region number for this load set

allocate(load_set(sn)%load_bc(counter))
load_set(sn)%load_bc(1:counter) = bc(1:counter)  ! DOF numbers for this load set

allocate(load_set(sn)%load(counter))
load_set(sn)%load(1:counter) = cn(1:counter)   ! load values along with DOF numbers

allocate(load_set(sn)%load_node(counter))
load_set(sn)%load_node(1:counter) = bc_node(1:counter)   ! node numbers along with DOF numbers


end subroutine write_load_bc







subroutine write_surface_load(rn, num_dim, n_points, points, angle, pressure, load_type_number, sn)

! Last modification: Aug, 2015 by Jeeho Lee
! 2D & 3D
! Input:
!       rn: region number
!       num_dim: dimension/dofs at one skin node (2 or 3)
!       nline: number of pressure loading lines
!       points(i,j): pressure point coordinates(i) for starting(j=1) and end(j=2) points
!       pressure: pressure value
!       load_type_number == 2: axisymmetric (nonlinear)
!                  == 3: axisymmetric loads (linear)
!       sn: load set number

implicit none

integer, intent(in) :: rn, num_dim, n_points, load_type_number, sn
real, intent(in) :: points(num_dim,n_points), angle, pressure

integer :: i, n_skin, dofs_skin, counter, pressure_point_st, pressure_point_end, group
integer :: num_surfaces, face_number, elem_number
integer, allocatable :: snode(:), bound_node(:), skin_node(:), face_nodes(:,:)
real :: tol
logical, allocatable :: bound_face(:)


n_skin = region(rn)%num_skin
dofs_skin = num_dim*n_skin

if (num_dim == 2) then

  allocate(bound_node(n_skin))
  allocate(skin_node(n_skin))

  tol = 0.1*region(rn)%skin_length ! tolerance based on minimum distance between two skin nodes
  call find_surface_node(rn, n_skin, points, tol, group, pressure_point_st, pressure_point_end)

  call find_skin(rn, group, pressure_point_st, pressure_point_end, counter, bound_node, skin_node, 0)

  allocate(snode(counter))
  snode = skin_node(1:counter)    ! skin point numbers
  region(rn)%skin_pressure(snode) = pressure  ! overwriting the pressure value !!
  region(rn)%skin_pressure_ends(1) = snode(1)
  region(rn)%skin_pressure_ends(2) = snode(counter)

  deallocate(bound_node, skin_node, snode)

elseif (num_dim == 3) then

  num_surfaces = region(rn)%num_skin_surfaces

  allocate(face_nodes(4,num_surfaces))
  allocate(bound_face(num_surfaces))

  do i = 1, num_surfaces
    face_number = region(rn)%skin_surfaces(1,i)
    elem_number = region(rn)%skin_surfaces(2,i)
    face_nodes(:,i) = region(rn)%element_face(1:4,face_number,elem_number)
  end do

  call search_contact_faces(1, num_surfaces, face_nodes, region(rn)%num_nodes, region(rn)%node_coord, points, angle, bound_face)

  forall (i = 1:num_surfaces, bound_face(i)) region(rn)%skin_pressure(i) = pressure

  region(rn)%surface_pressure_participation = bound_face

  do i = 1, num_surfaces
    face_number = region(rn)%skin_surfaces(1,i)
    elem_number = region(rn)%skin_surfaces(2,i)
    if (bound_face(i)) write(*,'(A,2I8,I3,E13.5)') '   * surface, element, face numbers & pressure:', i, elem_number,face_number, region(rn)%skin_pressure(i)
  end do
  write(*,'(A,I9/)') '  Number of pressure faces =', COUNT(bound_face)

  deallocate(bound_face, face_nodes)

endif


if (.NOT. load_set(sn)%surface_pressure_existence) then
  if (allocated(load_set(sn)%load_bc)) then  ! deallocate if not a pressure load set previously
    deallocate(load_set(sn)%load_bc)
    deallocate(load_set(sn)%load)
    deallocate(load_set(sn)%load_node)
  endif

  load_set(sn)%num_load = dofs_skin
  allocate(load_set(sn)%load_bc(dofs_skin))
  allocate(load_set(sn)%load(dofs_skin))
  allocate(load_set(sn)%load_node(dofs_skin))

  load_set(sn)%surface_pressure_existence = .TRUE.
  load_set(sn)%region_info = rn
  load_set(sn)%load_type = load_type_number

  region(rn)%skin_load_set(sn) = .TRUE.

else
  if (rn /= load_set(sn)%region_info) STOP 'write_surface_load: region number is inconsistent with load set number!'
endif

end subroutine write_surface_load








subroutine surface_load(flag, sn, rn)

! Last modification: Aug, 2015 by Jeeho Lee
! Transform pressure to nodal loads based on averaged gradient at a nodal point
! 2D & 3D version
! Input:
!       flag = 0: based on the original config
!            = 1: based on the deformed config (x0 + u)
!       sn: load set number
!       rn: region number

implicit none

integer, intent(in) :: flag, sn, rn

integer :: load_rn, gn, i, j, k, jj, n_skin, node, num_dim, num_surfaces, face_number, elem_number
integer, allocatable :: bound_face(:), face_nodes(:,:)
real :: load(3), pressure, normal_vector(3), s(4), nodal_areas(4), corners(3,4), trans(3,3), point(3)
logical :: in_or_out


load_rn = load_set(sn)%region_info

if (rn == load_rn) then

  load_set(sn)%load = 0.0
  num_dim = region(rn)%num_dim
  n_skin = region(rn)%num_skin

  if (num_dim == 2) then

    do i = 1, n_skin
      node = region(rn)%skin_nodes(i)
      pressure = region(rn)%skin_pressure(i)  ! point-wise pressure

      if (i == region(rn)%skin_pressure_ends(1)) then ! start point (left only)
        load(1:2) = pressure *region(rn)%skin_normal(1:2,i)
      elseif (i == region(rn)%skin_pressure_ends(2)) then ! end point (right only)
        load(1:2) = pressure *region(rn)%skin_normal(3:4,i)
      else
        load(1:2) = pressure *(region(rn)%skin_normal(1:2,i) + region(rn)%skin_normal(3:4,i))
      endif

      jj = 2*i-1
      load_set(sn)%load_node(jj:(jj+1)) = node                        ! node number
      load_set(sn)%load_bc(jj:(jj+1))   = region(rn)%node(4:5,node)   ! global DOF numbers
      load_set(sn)%load(jj:(jj+1))      = load(1:2)                   ! values

    end do  ! i

  elseif (num_dim == 3) then

    num_surfaces = region(rn)%num_skin_surfaces

    allocate(face_nodes(4,num_surfaces))
    allocate(bound_face(num_surfaces))

    do i = 1, num_surfaces  ! over skin surfaces

      face_number = region(rn)%skin_surfaces(1,i)
      elem_number = region(rn)%skin_surfaces(2,i)
      face_nodes(:,i) = region(rn)%element_face(1:4,face_number,elem_number)

      if (region(rn)%surface_pressure_participation(i)) then

        corners = region(rn)%node_coord(1:3,face_nodes(:,i))
        call surface_trans_matrix(corners, trans)

        point = 0.25*(corners(:,1) + corners(:,2) + corners(:,3) + corners(:,4))
        call point_position(rn, corners, trans, point, s, nodal_areas, in_or_out)
        if (.NOT. in_or_out) then
          write(*,'(A)') 'surface_load: the center point is not located in the element! corners:'
          write(*,'(4(3E12.3/))') (corners(:,j), j=1,4)
          write(*,'(A,2I7,I3,3E12.3/)') '  surface, elem, face numbers & point:', i, elem_number, face_number, point
          STOP
        endif
        pressure = region(rn)%skin_pressure(i)

        if (flag == 0) then
          normal_vector = trans(3,:)
        else
          normal_vector = region(rn)%surface_normal(1:3,i)
        endif

        do j = 1, 4 ! over nodes of an element
          node = face_nodes(j,i)
          load(1:3) = pressure * s(j) * normal_vector

          jj = 3*region(rn)%skin_node_map(node) - 2
          load_set(sn)%load_node(jj:(jj+2)) = node                                ! node number
          load_set(sn)%load_bc(jj:(jj+2)) = region(rn)%node(4:6,node)             ! global DOF numbers
          load_set(sn)%load(jj:(jj+2)) = load_set(sn)%load(jj:(jj+2)) + load(1:3) ! values
        end do  ! j

      else  ! assign DOF and node numbers for no-pressure surface load DB to avoid error

        do j = 1, 4 ! over nodes of an element
          node = face_nodes(j,i)

          jj = 3*region(rn)%skin_node_map(node) - 2
          load_set(sn)%load_node(jj:(jj+2)) = node
          load_set(sn)%load_bc(jj:(jj+2)) = region(rn)%node(4:6,node)
        end do  ! j

      endif
    end do  ! i

    deallocate(bound_face, face_nodes)

  endif ! num_dim
  

  if (flag == 0) then
    write(*,'(/A,3I7)') ' > load set & region, number of skin nodes = ', sn, rn, n_skin
    do i = 1, n_skin*num_dim
      if (ABS(load_set(sn)%load(i)) > 0) then
        node = load_set(sn)%load_node(i)
        write(*,'(A,3I9,E15.5)') '   pressure (internal & user) node, dof, load: ', node, region(rn)%node(1,node), load_set(sn)%load_bc(i), load_set(sn)%load(i)
      endif
    end do
  endif

endif ! load_rn


end subroutine surface_load







subroutine update_load_bc(sn, rn, node_dof, load, ndim)

! sn: load set number
! rn: region number
! node_dof(2,ndim): node # (1,:) & DOF # (2,:)
! load(ndim): load vector

implicit none

integer, intent(in) :: sn, rn, node_dof(2,ndim), ndim
real, intent(in) :: load(ndim)
integer :: i, m1, node, dof, load_dof(ndim)


if (allocated(load_set(sn)%load_bc)) then
  deallocate(load_set(sn)%load_bc)
  deallocate(load_set(sn)%load)
endif

do i = 1, ndim
  node = node_dof(1,i)
  dof = node_dof(2,i)
  load_dof(i) = region(rn)%node(3+dof,node)
end do

! followings are 'load_dof' ndim-size vector operation
! load_dof = region(rn)%dof_to_array(load_dof)


load_set(sn)%num_load = ndim
load_set(sn)%region_info = rn ! region number for this load set

allocate(load_set(sn)%load_bc(ndim))
load_set(sn)%load_bc = load_dof  ! DOF numbers for this load set

allocate(load_set(sn)%load(ndim))
load_set(sn)%load = load   ! load values along with DOF numbers

end subroutine update_load_bc






subroutine add_load_bc(flag, sn, rn, number, node, load, ndim, size)

! sn: load set number
! rn: region number
! node_dof(2,ndim): node # (1,:) & DOF # (2,:)
! load(ndim): load vector

implicit none

integer, intent(in) :: flag, sn, rn, number, node, ndim, size
real, intent(in) :: load(ndim)



if (flag == 0) then   ! setup
  if (allocated(load_set(sn)%load_bc)) then
    deallocate(load_set(sn)%load_bc)
    deallocate(load_set(sn)%load)
  endif
  allocate(load_set(sn)%load_bc(size))    ! 2D case
  allocate(load_set(sn)%load(size))
  load_set(sn)%num_load = size
  load_set(sn)%region_info = rn

else
  load_set(sn)%load_bc(number:number+ndim-1) = region(rn)%node(4:3+ndim,node)
  if (flag == -1) then
    load_set(sn)%load(number:number+ndim-1) = 0.0  
  elseif (flag == 1) then
    load_set(sn)%load(number:number+ndim-1) = load 
  endif
endif

end subroutine add_load_bc




	


subroutine loadfactor(unit_num_factor,load_factor_file)

implicit none 

integer, intent(in) :: unit_num_factor
character(len=20) :: load_factor_file

integer :: i, number, status
real :: time, factor

open(unit=unit_num_factor, file=load_factor_file, status='OLD', action='READ')
read(unit_num_factor, *, iostat=status) number
allocate(load_factor_db(number,2))

do i = 1, number
  read(unit_num_factor, *, iostat=status) time, factor
	load_factor_db(i,1) = time
  load_factor_db(i,2) = factor
end do

end subroutine loadfactor







subroutine write_temp_load(ref_temperature, dim, time, time_load)

implicit none 

real, intent(in) :: ref_temperature
integer,intent(in) :: dim
integer :: load_size ,  i 
real, intent(in) :: time(dim), time_load(dim)
real :: tol 


tol = 1.0e4*Numerical_Zero

if (time(1) < tol) then 
  load_size = dim + 1
else 
  load_size = dim
endif

temperature%reference_temperature = ref_temperature

temperature%num_temp= load_size
allocate(temperature%step_time(load_size), temperature%temp_value(load_size)) 

if (dim /= load_size) then 
  temperature%step_time(1) = 0.0
  temperature%temp_value(1) = 0.0
  temperature%step_time(2:load_size) = time
  temperature%temp_value(2:load_size) = time_load
else 
  temperature%step_time = time
  temperature%temp_value = time_load
endif

end subroutine write_temp_load






subroutine get_temperature(now, result, ref_temperature)

implicit none

real, intent(in) :: now
real, intent(out) :: result, ref_temperature
integer :: load_size , i 
real :: slope , delta_t, tol
logical :: found 
 
 
tol = 1.0e4*Numerical_Zero
found = .FALSE.

load_size = temperature%num_temp
do i = 2 , load_size
  if (abs(temperature%step_time(i) -  now) < tol) then
    result = temperature%temp_value(i) 
    found = .TRUE.
    EXIT !------------------------------------------->>
  elseif (temperature%step_time(i) > now) then
    slope =  (temperature%temp_value(i)-temperature%temp_value(i-1))/(temperature%step_time(i)-temperature%step_time(i-1))
    delta_t = now - temperature%step_time(i-1)
    result = slope * delta_t + temperature%temp_value(i-1)
    found = .TRUE.
    EXIT !------------------------------------------->> 
  endif
end do   

! if (found) write(*,*) '+++ time & temperature: ', now, result 
ref_temperature = temperature%reference_temperature

end subroutine get_temperature





subroutine write_gr_motion(rn, dof, factor)

implicit none

integer, intent(in) :: rn, dof
real, intent(in) :: factor


region(rn)%grfactor(dof) = factor


end subroutine write_gr_motion







subroutine find_surface_node(rn, max_finding, coord, tol, group, start_point, end_point)

! 2D only

implicit none

integer, intent(in) :: rn, max_finding
real, intent(in) :: coord(2,2), tol
integer, intent(out) :: group, start_point, end_point

integer :: i, j, n, n_skin, node_found1(max_finding), node_found2(max_finding), n_nodes1, n_nodes2
integer, allocatable :: group_cand1(:), group_cand2(:)
real, allocatable :: skin_coord(:,:)
logical :: found

n_skin = region(rn)%num_skin
allocate(skin_coord(2, n_skin)) 
skin_coord = region(rn)%node_coord(1:2,region(rn)%skin_nodes(1:n_skin))

! start point --------------------------------
call find_match_node(tol, coord(:,1), skin_coord, n_skin, node_found1, max_finding, n_nodes1)
if (n_nodes1 < 1) then
  write(*,'(A,I5/)') 'find_surface_node: No node is found. Check the 1st input point! region = ', rn
  STOP
endif
allocate(group_cand1(n_nodes1))
call skin_group_finder(rn, node_found1, n_nodes1, group_cand1)

! end point ----------------------------------
call find_match_node(tol, coord(:,2), skin_coord, n_skin, node_found2, max_finding, n_nodes2)
if (n_nodes2 < 1)  then
  write(*,'(A,I5/)') 'find_surface_node: No node is found. Check the 2nd input point! region =', rn
  STOP
endif
allocate(group_cand2(n_nodes2))
call skin_group_finder(rn, node_found2, n_nodes2, group_cand2)

!---------------------------------------------
if ((n_nodes1 == 1).AND.(n_nodes2 == 1)) then
  if (group_cand1(1) == group_cand2(1)) then
    start_point = node_found1(1)
    end_point = node_found2(1)
    group = group_cand1(1)
  else
    STOP 'find_surface_node: no group match(1)!'
  endif
else
  found = .FALSE.
  gfinder: do i = 1, n_nodes1
    do j = 1, n_nodes2
      if (group_cand1(i) == group_cand2(j)) then
        start_point = node_found1(i)
        end_point = node_found2(j)
        group = group_cand1(i)
        found = .TRUE.
        EXIT gfinder
      endif
    end do
  end do gfinder
  if (.NOT. found) STOP 'find_surface_node: no group match(2)!'
endif

start_point = region(rn)%skin_nodes(start_point)  ! transform from surface number to node number
end_point = region(rn)%skin_nodes(end_point)

write(*,*) '########### group, start & end skin nodes:', group, start_point, end_point
deallocate(skin_coord, group_cand1, group_cand2)

end subroutine find_surface_node








subroutine find_point_number(rn,elem_num,quad_num,point_address)

implicit none

integer, intent(in) :: rn, elem_num, quad_num
integer, intent(out) :: point_address
integer :: n

n = 3 + region(rn)%nodes_element + quad_num
point_address = region(rn)%element(n,elem_num)

end subroutine find_point_number





subroutine find_point_number_spring(rn,spring_num,point_address)

! 'find_point_number' for spring (connection) elements

implicit none

integer, intent(in) :: rn, spring_num
integer, intent(out) :: point_address
integer :: n

n = 3 + (2 + 1)   ! always # of spring nodes = 2, quad_num = 1
point_address = region(rn)%spring(n,spring_num)

end subroutine find_point_number_spring







subroutine set_contact_group(rn, num_contact_group)

implicit none

integer, intent(in) :: rn, num_contact_group
integer :: num_surfaces


num_surfaces = region(rn)%num_skin_surfaces

if (allocated(region(rn)%contact_face_participation)) deallocate(region(rn)%contact_face_participation)
allocate(region(rn)%contact_face_participation(num_surfaces, num_contact_group))
region(rn)%contact_face_participation = .FALSE.


end subroutine set_contact_group






subroutine extract_nodes_from_faces(rn, contact_group, nn, nodes, n_nodes)

! July 2015, coded by Jeeho Lee
!
! to extract nodes from faces

implicit none

integer,intent(in) :: rn, contact_group, nn
integer, intent(out) :: n_nodes, nodes(nn)

integer :: i, num_nodes, nsize
integer :: num_surfaces, face_number, elem_number
logical, allocatable :: check(:)


num_nodes = region(rn)%num_nodes
allocate(check(num_nodes))
check = .FALSE.

nsize = COUNT(region(rn)%contact_face_participation(:,contact_group))

num_surfaces = region(rn)%num_skin_surfaces

do i = 1, num_surfaces
  if (region(rn)%contact_face_participation(i,contact_group)) then
    face_number = region(rn)%skin_surfaces(1,i)
    elem_number = region(rn)%skin_surfaces(2,i)
    check(region(rn)%element_face(1:4,face_number,elem_number)) = .TRUE.
  endif
end do

n_nodes = 0
do i = 1, num_nodes
  if (check(i)) then
    n_nodes = n_nodes + 1
    nodes(n_nodes) = i
  endif
end do

deallocate(check)

end subroutine extract_nodes_from_faces







subroutine find_contact_surface(rn, contact_group, corners)

! July 2015, coded by Jeeho Lee
!
! to write the surface of 3D solid object

implicit none

integer,intent(in) :: rn, contact_group
real, intent(in) :: corners(3,4)

integer :: i
integer :: num_surfaces, face_number, elem_number
integer, allocatable :: face_ix(:,:)
real, parameter :: angle = 90


num_surfaces = region(rn)%num_skin_surfaces

allocate(face_ix(4,num_surfaces))

do i = 1, num_surfaces
  face_number = region(rn)%skin_surfaces(1,i)
  elem_number = region(rn)%skin_surfaces(2,i)
  face_ix(:,i) = region(rn)%element_face(1:4,face_number,elem_number)
end do

call search_contact_faces(1, num_surfaces, face_ix, region(rn)%num_nodes, region(rn)%node_coord, corners, angle, region(rn)%contact_face_participation(:,contact_group))

deallocate(face_ix)


end subroutine find_contact_surface








subroutine write_surface(rn)

! June 2015, coded by Jeeho Lee
! 3D
! write the surface of 3D solid object

implicit none

integer,intent(in) :: rn

integer :: ndim, nn, ne, nodes_elem, num_surfaces
integer :: i, j, counter
logical, allocatable :: surface_nodes(:)


nodes_elem = region(rn)%nodes_element
if (nodes_elem /= 8) STOP 'write_surface: currently 8-node brick element only!'

ndim = region(rn)%num_dim 
if (ndim /= 3) STOP 'write_skin: Work only for 3D!'

nn = region(rn)%num_nodes
ne = region(rn)%num_elements

allocate(surface_nodes(nn))
surface_nodes = .FALSE.

num_surfaces = 6*ne - region(rn)%num_inner_faces
if (num_surfaces < 6) STOP 'write_surface: number of external surfaces cannot be smaller than 6!'

allocate(region(rn)%skin_surfaces(2,num_surfaces))
allocate(region(rn)%surface_normal(3,num_surfaces))
allocate(region(rn)%surface_pressure_participation(num_surfaces))
allocate(region(rn)%skin_pressure(num_surfaces))
region(rn)%skin_pressure = 0.0
region(rn)%surface_pressure_participation = .FALSE.

counter = 0
do i = 1, ne
  do j = 1, 6
    if (region(rn)%element_face_participation(j,i)) then
      counter = counter + 1
      region(rn)%skin_surfaces(1,counter) = j ! face number
      region(rn)%skin_surfaces(2,counter) = i ! eleemnt number
      surface_nodes(region(rn)%element_face(1:4,j,i)) = .TRUE.
    endif
  end do  ! j
end do  ! i

if (counter /= num_surfaces) then
  write(*,*) '  counter, num_surfaces =', counter, num_surfaces
  STOP 'write_surface: number of surfaces is not consistent with the counter!'
endif

region(rn)%num_skin_surfaces = num_surfaces

counter = 0
do i = 1, nn
  if (surface_nodes(i)) then
    counter = counter + 1
    region(rn)%skin_node_map(i) = counter ! map from the nodal number to surface nodal number
  endif
end do
region(rn)%num_skin = counter

write(*,*) '  > number of inner faces:', region(rn)%num_inner_faces
write(*,*) '  > number of surfaces:', num_surfaces


deallocate(surface_nodes)

end subroutine write_surface







subroutine write_skin(unit_num, rn)

! Modified: Jan 25, 2013 by Jeeho Lee
!
! The subroutine works only for 2D & Q4 elements!

implicit none

integer,intent(in) :: unit_num, rn

integer :: i, j, k, ndim, nn, ne, nodes_elem, num_bound, status, node_number, ndof
integer :: node_num, sort_flag, number, ngroup, nn_group, ne_group, nn1, nn2, ne1, ne2
integer :: st_node, ed_node
integer, allocatable :: bound_nodes_all(:), bound_nodes(:), nodes(:), elements(:,:)
real, allocatable :: coord(:,:), u(:,:)
real :: length


nodes_elem = region(rn)%nodes_element
if (nodes_elem /= 4) STOP 'write_skin: currently 4-node element only!'

ndim = region(rn)%num_dim 
if (ndim /= 2) STOP 'write_skin: Work only for 2D!'

nn = region(rn)%num_nodes
ne = region(rn)%num_elements
ndof =region(rn)%num_node_dof
ngroup = region(rn)%num_groups

allocate(bound_nodes_all(nn)) 

num_bound = 0
do i = 1, ngroup
  region(rn)%skin_group_pointer(i) = num_bound + 1

  if (i < ngroup) then
    nn1 = region(rn)%node_group(i)
    nn2 = region(rn)%node_group(i+1) - 1
    ne1 = region(rn)%element_group(i)
    ne2 = region(rn)%element_group(i+1) - 1
  else
    nn1 = region(rn)%node_group(i)
    nn2 = nn
    ne1 = region(rn)%element_group(i)
    ne2 = ne    
  endif
  nn_group = nn2 - nn1 + 1
  ne_group = ne2 - ne1 + 1
  allocate(nodes(nn_group))
  allocate(elements(4,ne_group))
  allocate(coord(ndim,nn_group))
  allocate(bound_nodes(nn_group))

  nodes = region(rn)%node(3,nn1:nn2)
  coord = region(rn)%node_coord(:,nn1:nn2)
  elements = region(rn)%element(4:7,ne1:ne2) - nn1 + 1
  sort_flag = -1  ! default sorting start point

  call mesh_surface(ndim, nn_group, nodes, coord, ne_group, nodes_elem, elements, sort_flag, number, bound_nodes) ! bound_nodes: internal node numbers

  bound_nodes_all(num_bound+1:num_bound+number) = bound_nodes(1:number) + nn1 - 1
  num_bound = num_bound + number

  deallocate(nodes, elements, coord, bound_nodes)
end do

region(rn)%num_skin = num_bound     ! number of boundary nodes
allocate(region(rn)%skin_nodes(num_bound))
allocate(region(rn)%skin_pressure(num_bound))
allocate(region(rn)%skin_normal(2*ndim,num_bound))
allocate(region(rn)%skin_node_group(num_bound))


region(rn)%skin_nodes(:) = bound_nodes_all(1:num_bound)

region(rn)%skin_node_group = 1  ! default group number: 1
region(rn)%skin_pressure = 0.0
region(rn)%skin_normal = 0.0 

write(*,*) ' > number of boundary nodes:', region(rn)%num_skin
write(*,*) ' > boundary nodes:', region(rn)%skin_nodes


do i = 1, ngroup
  if (i < ngroup) then
    st_node = region(rn)%skin_group_pointer(i)
    ed_node = region(rn)%skin_group_pointer(i+1) - 1
  else
    st_node = region(rn)%skin_group_pointer(i)
    ed_node = region(rn)%num_skin
  endif
  num_bound = ed_node - st_node + 1

  region(rn)%skin_node_group(st_node:ed_node) = i ! assign the group number to each skin node

  allocate(u(ndim, num_bound))
  u = 0.0   ! dummy u
  call update_skin_normal(0, rn, i, ndim, num_bound, u, length)
  deallocate(u)
  if (i == 1) then
    region(rn)%skin_length = length
  elseif (length < region(rn)%skin_length) then
    region(rn)%skin_length = length
  endif

end do

deallocate(bound_nodes_all)

end subroutine write_skin





subroutine update_skin_normal(flag, rn, group, ndim, num_skin, u, min_len)

! 2D only
! flag = 0: based on the original config
!      = 1: based on the deformed config (x0 + u)

implicit none

integer, intent(in) :: rn, group, ndim, num_skin, flag
real, intent(in) :: u(ndim,num_skin)
real, intent(out) :: min_len
integer :: i, n1, n2, n
real :: coord(ndim,3), vec(ndim), norm_vec(ndim), len
real :: skin_coord(ndim, num_skin) 
logical :: masu

if (ndim /= 2) STOP 'update_skin_normal: Work only for 2D!'

masu = .TRUE.
min_len = 1.0e2*Numerical_Zero


n1 = region(rn)%skin_group_pointer(group) ! first node
n2 = n1 + num_skin - 1            ! last node

if (flag == 0 ) then
  skin_coord = region(rn)%node_coord(:,region(rn)%skin_nodes(n1:n2))
else
  skin_coord = region(rn)%node_coord(:,region(rn)%skin_nodes(n1:n2)) + u
endif 

do i = 1, num_skin
	coord(:,2) = skin_coord(:,i)    ! current skin point
  
	if (i == num_skin) then    
		coord(:,3) = skin_coord(:,1)  ! starting point
	else
		coord(:,3) = skin_coord(:,i+1)
	endif
  
  coord(:,2) = coord(:,3) - coord(:,2)  ! tangntial vector along the edge (i, i+1)
  call normalize(coord(:,2), 2, vec, len) ! normalized tangntial vector with its length
  if (len > Numerical_Zero) then
    if (masu) then
      min_len = len
      masu = .FALSE.
    elseif (min_len > len) then
      min_len = len
    endif
  else
    write(*,*) 'update_skin_normal: WARNING! Skin node is too close to its next!', i
  endif
    
  norm_vec(1) = vec(2)
  norm_vec(2) = - vec(1)
  norm_vec = 0.5*len*norm_vec ! scaled by half of the edge length for pressure computation
  
  n = n1 + i - 1
	region(rn)%skin_normal(1:2,n) = norm_vec    ! normal vector along left side for i-th point

  if (i < num_skin) then
    region(rn)%skin_normal(3:4,n+1) = norm_vec  ! normal vector along right side for (i+1)-th point
  else
    region(rn)%skin_normal(3:4,n1) = norm_vec   ! normal vector at the first node
  endif
end do 

end subroutine update_skin_normal






subroutine update_surface_normal(rn, n, normal)

! Last modification: Aug, 2015 by Jeeho Lee
! 3D version

integer, intent(in) :: rn, n
real, intent(in) :: normal(3,n)


region(rn)%surface_normal(1:3,1:n) = normal


end subroutine update_surface_normal






subroutine find_skin(rn, group, st, ed, num_nodes, bound_node, skin_node, address) 

! bound_node: boundary point node numbers (internal nodal number)
! skin_node: boundary point order numbers [1:number of boundary nodes]

implicit none

integer, intent(in) :: rn, group, st, ed, address
integer, intent(out) :: num_nodes, bound_node(:), skin_node(:)

integer :: i, j, k, position, n1, n2, num_nodes_1
logical :: find_next


write(*,'(A,4I8)') 'rn, group, st, ed:', rn, group, st, ed

n1 = region(rn)%skin_group_pointer(group) ! first node

if (group < region(rn)%num_groups) then  
  n2 = region(rn)%skin_group_pointer(group+1) - 1 ! last node
else
  n2 = region(rn)%num_skin    ! last node in the last group
endif
find_next = .TRUE.

!-----------------------------------------
position = 0
do i = n1, n2
	if (region(rn)%skin_nodes(i) == st) then
		position = i 
		EXIT
	endif
end do
if (position == 0) STOP 'find_skin: fail to find starting position!'

skin_node(1+address) = position
bound_node(1+address) = region(rn)%skin_nodes(position)

num_nodes = 1 
do i = position+1, n2
	num_nodes = num_nodes + 1 
  skin_node(num_nodes+address) = i
	bound_node(num_nodes+address) = region(rn)%skin_nodes(i)
	if (region(rn)%skin_nodes(i) == ed) then
		find_next = .FALSE.
		EXIT
	endif
end do
write(*,*) 'skin_node1:', bound_node(1+address:num_nodes+address)

! if the end point exists beyond the last boundary point: search from the skin loop beginnng

if (find_next) then 
  num_nodes_1 = num_nodes + 1
	do i = n1, position-1
		num_nodes = num_nodes + 1
    skin_node(num_nodes+address) = i
		bound_node(num_nodes+address) = region(rn)%skin_nodes(i)
		if (region(rn)%skin_nodes(i) == ed) then
			EXIT
		endif
	end do
  write(*,*) 'skin_node2:', bound_node(num_nodes_1+address:num_nodes+address)
endif 

write(*,*) 'number of skin nodes=',  num_nodes

end subroutine find_skin






subroutine skin_group_finder(rn, skin_node, n_nodes, group)

! group finder for a given surface point(node)

implicit none

integer, intent(in) :: rn, skin_node(n_nodes), n_nodes
integer, intent(out) :: group(n_nodes)
integer :: n, i, j, position 


n = region(rn)%num_groups
do j = 1, n_nodes
  if (n == 1) then
    group(j) = 1
  else
    group(j) = n
    do i = 2, n
      if (skin_node(j) < region(rn)%skin_group_pointer(i)) then
        group(j) = i - 1
        EXIT
      endif
    end do
  endif
end do

end subroutine skin_group_finder






subroutine skin_neighbor_finder(rn, group, node, node_f, node_b)

! node finder for forward and backward surface points of a given surface point(node)

implicit none

integer, intent(in) :: rn, group, node
integer, intent(out) :: node_f, node_b
integer :: i, position, n1, n2


! in case of failure to find nodes, return zeros
node_f = 0
node_b = 0

n1 = region(rn)%skin_group_pointer(group) ! first node

if (group < region(rn)%num_groups) then  
    n2 = region(rn)%skin_group_pointer(group+1) - 1 ! last node
else
    n2 = region(rn)%num_skin    ! last node in the last group
endif

do i = n1, n2
	if (region(rn)%skin_nodes(i) == node) then
		position = i 
		EXIT
	endif
end do

if (position == n1) then ! head case
  node_f = region(rn)%skin_nodes(position+1)
  node_b = region(rn)%skin_nodes(n2)
elseif (position == n2) then  ! tail case
  node_f = region(rn)%skin_nodes(n1)
  node_b = region(rn)%skin_nodes(position-1)
else
  node_f = region(rn)%skin_nodes(position+1)
  node_b = region(rn)%skin_nodes(position-1)
endif

end subroutine skin_neighbor_finder





subroutine surface_node_length(rn, group, node, flag, length)


! January 2015
! flag =  1: use only forward node
!      = -1: use only backward node
!      = otherwise: use both backward and forward nodes (default)

implicit none

integer, intent(in) :: rn, group, node, flag
real, intent(out) :: length
integer :: node_f, node_b, n
real :: l1, l2

n = region(rn)%num_dim

call skin_neighbor_finder(rn, group, node, node_f, node_b)

if (node_f*node_f == 0) STOP 'surface_length: the node is not on the surface!'

if (flag == 1) then
  call nodal_distance(rn, node, node_f, n, length)

elseif (flag == -1) then
  call nodal_distance(rn, node, node_b, n, length)

else
  call nodal_distance(rn, node, node_f, n, l1)
  call nodal_distance(rn, node, node_b, n, l2)
  length = 0.5*(l1 + l2)
endif

end subroutine surface_node_length






subroutine nodal_distance(rn, node1, node2, n, length)

! January 2015
! Calculate length between two nodes

integer, intent(in) :: rn, node1, node2, n
real, intent(out) :: length
real :: coord1(3) , coord2(3), calc_length_n

coord1(1:n) = region(rn)%node_coord(1:n, node1)     
coord2(1:n) = region(rn)%node_coord(1:n, node2) 
   
length = calc_length_n(coord1 , coord2, n) 

end subroutine nodal_distance








subroutine write_skin_contact(rn, bound_nodes, ndim)

implicit none

integer, intent(in) :: rn, ndim, bound_nodes(ndim)


if (allocated(region(rn)%skin_contact_nodes)) deallocate(region(rn)%skin_contact_nodes)

allocate(region(rn)%skin_contact_nodes(ndim))
region(rn)%skin_contact_nodes = bound_nodes
region(rn)%num_skin_contact_nodes = ndim

end subroutine write_skin_contact






subroutine write_subregion_group(flag, rn, group, node_pointer, elem_pointer)

! flag: for future use

implicit none

integer, intent(in) :: flag, rn, group, node_pointer, elem_pointer

region(rn)%node_group(group) = node_pointer
region(rn)%element_group(group) = elem_pointer

end subroutine write_subregion_group







subroutine element_output(unit_num, output_control, rn, output_rn, current_i, output_interval)

! output_interval: if negative, only current time step data are written (overwrite!)

implicit none

integer, intent(in) :: unit_num, rn, output_rn, current_i, output_control, output_interval
integer :: i, j, k, ne, n_quad, mat_num, num_prst, num_field, max_str_output
character(len=8) :: name_2d(8)=(/'stress11','stress22','stress33','stress12','state001','state002','state003','state004'/)
character(len=8) :: name_3d(12)=(/'stress11','stress22','stress33','stress12','stress23','stress31', &
																  'state001','state002','state003','state004','state005','state006'/)
character(len=8) :: prst_name(3)=(/'prncpl_1','prncpl_2','prncpl_3'/)
real :: prst(3), vector(6)

if ((output_rn < 1).OR.(rn == output_rn)) then  !--------------------------------------------

  max_str_output = region(rn)%size_str_output
! 'output_control' is smaller than 20, nothing will be done in this subroutine

if ((output_control >= 20).and.(output_control < 30)) then   ! 2D cases
!	write(*,'(A,I3)') ' * 2D output for element stress/strain are written!: cnt # = ', output_control
	num_field = 8
	if (max_str_output < num_field) STOP 'element_output: num_field exceeds maximum!'
	if (output_control == 21) then
		num_prst = 2
	else
		num_prst = 0
	endif


	if (mod(current_i,abs(output_interval)) == 0) then
  
  	write(unit_num,100)  '*number of stress and state output fields =', (num_field+num_prst),'  @ current time step =', current_i
    write(unit_num,101)  'mat# elem# quad#', (name_2d(k), k=1,num_field), (prst_name(k), k=1,num_prst)
    
		ne = region(rn)%num_elements
		n_quad = region(rn)%num_quad_pts
		do i = 1, ne
			mat_num = region(rn)%element(3,i)
			do j = 1, n_quad
				if (output_control == 21) then
					vector(1:4) = region(rn)%element_result(j,1:4,i)
					call spect_2d(vector(1:4), prst(1:2))
				endif
				write(unit_num,102) mat_num, i, j, (region(rn)%element_result(j,k,i), k=1,num_field), &
														(prst(k), k=1,num_prst)
			end do
		end do
	endif
	
elseif (output_control >= 30) then   ! 3D cases
!	write(*,'(A,I3)') ' * 3D output for element stress/strain are written!: cnt # = ', output_control
	num_field = 12
	if (max_str_output < num_field) STOP 'element_output: num_field exceeds maximum!'
	if (output_control == 31) then
		num_prst = 3
	else
		num_prst = 0
	endif
	if (mod(current_i,abs(output_interval)) == 0) then
  
  	write(unit_num,100)  '*number of stress and state output fields =', (num_field+num_prst),'  @ current time step =', current_i
    write(unit_num,101)  'mat# elem# quad#', (name_3d(k), k=1,num_field), (prst_name(k), k=1,num_prst)
		ne = region(rn)%num_elements
		n_quad = region(rn)%num_quad_pts
		do i = 1, ne
			mat_num = region(rn)%element(3,i)
			do j = 1, n_quad
				if (output_control == 31) then
					vector(1:6) = region(rn)%element_result(j,1:6,i)
					call spect_3d(vector, prst)
				endif
				write(unit_num,102) mat_num, i, j, (region(rn)%element_result(j,k,i), k=1,num_field), &
														(prst(k), k=1,num_prst)
			end do
		end do
	endif
endif
if (output_interval < 0) then
	rewind(unit=unit_num)
else
	write(unit_num,*)
endif
! write(*,*) ' --------------------------------------------------------------'

endif  !-------------------------------------------------------------------------------------

100 format (A,I4,A,I9)
101 format (A,100(7X,A))
102 format (I4,I6,I6,2X,100ES15.7)

end subroutine element_output


!==============================================================================

subroutine update_coord(rn, num_dim, nn, coord )

implicit none

integer, intent(in) :: rn , num_dim, nn
real, intent(in) :: coord(num_dim, nn)

region(rn)%node_coord(1:num_dim,:) = coord

end subroutine update_coord


subroutine update_constr(sn, num, disp)

implicit none

integer, intent(in) :: sn, num
real, intent(in) :: disp

disp_set(sn)%constr(num) = disp

end subroutine update_constr



subroutine get_ele_connect( rn , el_num , nel , conn )

implicit none

integer , intent(in) :: rn, el_num, nel
integer , intent(out) :: conn( nel )

integer :: i


do i = 1, nel
    conn( i ) =  region(rn)%element(3+i , el_num )
enddo

end subroutine get_ele_connect



subroutine get_ele_quad( rn , el_num , nqd , conn )

implicit none

integer , intent(in) :: rn, el_num, nqd
integer , intent(out) :: conn( nqd )

integer :: i ,n

n = 3 + region(rn)%nodes_element

do i = 1, nqd
    conn( i ) =  region(rn)%element(n+i , el_num )
enddo

end subroutine get_ele_quad



subroutine fsi_load_set(sn, ni, inte_info)

integer, intent(in) :: sn, ni, inte_info(ni)

integer :: i, ln

ln = 0
do i = 1, ni
    if ( inte_info(i) == 1 ) ln = ln + 1
enddo        

if (allocated(load_set(sn)%load_bc)) then
    deallocate(load_set(sn)%load_bc, load_set(sn)%load, load_set(sn)%load_node)
endif

load_set(sn)%num_load = ln*2
allocate(load_set(sn)%load_bc(load_set(sn)%num_load))
allocate(load_set(sn)%load(load_set(sn)%num_load))
allocate(load_set(sn)%load_node(load_set(sn)%num_load))
load_set(sn)%load = 0

end subroutine fsi_load_set



subroutine update_fsi_load(rn, Dim, ndim, load_point, load_value)

implicit none

integer, intent(in) :: rn, Dim, ndim, load_point(2,ndim)
real, intent(in) :: load_value(ndim)

integer :: i, j, k, count, n_nodes, coarse_count, ori_num_load, num_load
integer :: the_dof, load_nodes(ndim), found(ndim)
integer, allocatable :: nodes(:), ori_load_bc(:), ori_load_node(:)
real, allocatable :: ori_load(:)
integer :: bc(ndim), bc_node(ndim)
real :: cn(ndim)
logical :: check
!-------------------------------------------------------------------------------
n_nodes = region(rn)%num_nodes
allocate(nodes(n_nodes), ori_load_bc(n_nodes*Dim), ori_load(n_nodes*Dim), ori_load_node(n_nodes*Dim))
ori_num_load = load_set(rn)%num_load
ori_load_bc(1:ori_num_load) = load_set(rn)%load_bc(1:ori_num_load)
ori_load(1:ori_num_load) = load_set(rn)%load(1:ori_num_load)
ori_load_node(1:ori_num_load) = load_set(rn)%load_node(1:ori_num_load)

!===========================================================
! After fsi load subtracted from original load, add new fsi load 
! But when currently run fsi, original load is zero
num_load = 0
ori_load_bc = 0
ori_load = 0
ori_load_node = 0
load_set(rn)%load_bc = 0
load_set(rn)%load = 0
load_set(rn)%load_node = 0 
!==============================================================

nodes = region(rn)%node(1,:)
load_nodes = load_point(1,:)
found = region(rn)%node_number_table(load_nodes) 
!call find_my_match(nodes, n_nodes, load_nodes, found, ndim)

do i = 1, ndim
    if (found(i) /= 0 ) then 
        k = load_point(2,i)   ! DOF number at one point
        the_dof = region(rn)%node(3+k,found(i))
        coarse_count = region(rn)%node(2,found(i))
        if (the_dof > 0) then
            check = .TRUE.
            do j = 1, ori_num_load
                if ( ori_load_bc(j) == the_dof .AND. ori_load_node(j) == found(i)) then
                    ori_load(j) = ori_load(j) + load_value(i)/coarse_count
                    check = .FALSE.
                    exit
                endif
            enddo
            if ( check ) then
                if (abs(load_value(i)) > 10e-8) then
                    num_load = num_load + 1
                    ori_load_bc(num_load) = the_dof
                    ori_load_node(num_load) = found(i)
                    ori_load(num_load) = load_value(i)/coarse_count
                    !write (*,*) num_load, ori_load_bc(num_load), ori_load_node(num_load), ori_load(num_load)
                endif
            endif
        endif
    endif
enddo
if (num_load == 0) then
    num_load = 1
    ori_load_bc(1) = 1
    ori_load_node(1) = 1
    ori_load(1) = 0.0
endif

if (load_set(rn)%num_load /= num_load) then
    load_set(rn)%num_load = num_load
    deallocate (load_set(rn)%load_bc, load_set(rn)%load, load_set(rn)%load_node)
    allocate (load_set(rn)%load_bc(num_load))
    allocate (load_set(rn)%load(num_load))
    allocate (load_set(rn)%load_node(num_load))
endif
load_set(rn)%load_bc(1:num_load) = ori_load_bc(1:num_load)  ! DOF numbers for this load set
load_set(rn)%load(1:num_load) = ori_load(1:num_load)   ! load values along with DOF numbers
load_set(rn)%load_node(1:num_load) = ori_load_node(1:num_load)   ! node numbers along with DOF numbers

deallocate ( nodes, ori_load_bc, ori_load, ori_load_node )

end subroutine update_fsi_load



subroutine free_pdomain

implicit none 

if (allocated(load_set)) deallocate(load_set)
if (allocated(disp_set)) deallocate (disp_set)
if (allocated(disp_set_backup)) deallocate (disp_set_backup)
deallocate(region) 

pdomain_mat_db%num_materials = 0 
deallocate( pdomain_mat_db%mat_num_model_num ) 

if ( temperature%num_temp /= 0 )  deallocate( temperature%step_time , temperature%temp_value )
if ( size(load_factor_db) /= 0 ) deallocate( load_factor_db )

end subroutine free_pdomain


subroutine set_skin(rn, group, sp, adjn)

integer, intent(in) :: rn, group
integer, intent(out) :: sp, adjn

integer :: ngroup

ngroup = region(rn)%num_groups
sp = region(rn)%skin_group_pointer(group)
if (group < ngroup) then  
    adjn = region(rn)%skin_group_pointer(group+1) - region(rn)%skin_group_pointer(group)
else
    adjn = region(rn)%num_skin + 1 - region(rn)%skin_group_pointer(group)  
endif

end subroutine set_skin


subroutine free_skin_contact_nodes( rn )

implicit none

integer, intent(in) :: rn

deallocate ( region(rn)%skin_contact_nodes )

end subroutine free_skin_contact_nodes



subroutine free_constr(sn)

implicit none

integer :: sn

if (allocated(disp_set(sn)%constr)) then
    deallocate (disp_set(sn)%constr)
    deallocate (disp_set_backup(sn)%constr)
endif

end subroutine free_constr
   
!==============================================================================

end module physical_domain
