
!----------------- matching & finding algorithms -------------------------------------------------------------


subroutine find_match(v, n, u, x, ndim)

! v(n): target vector
! u(ndim): screen vector
! x(ndim): result vector containing array index number of v where its value matches u
! Note: it may cause this warning error: fort: (1):" In call to FIND_MATCH, an array temporary was created for argument #1"

implicit none

integer, intent(in) :: v(n), n, u(ndim), ndim
integer, intent(out) :: x(ndim)
integer :: i, j, m
logical :: toggle

do i = 1, ndim
	toggle = .FALSE.
	m = u(i)

  if (m < n) then 
    do j = m, n
      if (m == v(j)) then
        x(i) = j
        toggle = .TRUE.
        EXIT
      endif
    end do
	
    if (.NOT. toggle) then
      do j = m-1, 1, -1
        if (m == v(j)) then
          x(i) = j
          toggle = .TRUE.
          EXIT
        endif
      end do
    endif
  else
    do j = 1, n
      if (m == v(j)) then
        x(i) = j
        toggle = .TRUE.
        EXIT
      endif
    end do		
  endif

	if (.NOT. toggle) STOP 'find_match: fail to find a match!'
end do

end subroutine find_match






subroutine find_my_match(v, n, u, x, ndim)

! v(n): target vector
! u(ndim): screen vector
! x(ndim): result vector containing array index number of v where its value matches u

implicit none

integer, intent(in) :: v(n), n, u(ndim), ndim
integer, intent(out) :: x(ndim)
integer :: i, j, m, count
logical :: toggle

x = 0
do i = 1, ndim
	toggle = .FALSE.
	m = u(i)
  do j = 1, n
    if (m == v(j)) then
			x(i) = j
      toggle = .TRUE.
      exit
    endif
  end do
end do

end subroutine find_my_match








subroutine find_match_coord(ndim, target_coord, num_target, coord, nn, node_num, distance)

! Aug. 2015
! finding minimum distance point in a nodal point system: multi-target algorithm
! 2D & 3D
! Input:
!       ndim: spatial dimension
!       num_target: number of target point(s)
!       target_coord: target coordinate to be compared with
!       coord: coordinates of nodal points of a system
! Output:
!       node_num: found point number(s) in a nodal system
!       distance: distance between the found point in a nodal system and target point(s)

implicit none

integer, intent(in) :: ndim, num_target, nn
real, intent(in) :: target_coord(ndim,num_target), coord(ndim,nn)
integer, intent(out) :: node_num(num_target)
real, intent(out) :: distance(num_target)

integer :: i, j
real :: DNRM2               ! BLAS Library
real :: vector(ndim,nn), len(nn)


do i = 1, num_target
  do j = 1, nn
    vector(:,j) = target_coord(:,i) - coord(:,j)
    len(j) = DNRM2(ndim, vector(1:ndim,j), 1)
  end do
  node_num(i) = MINLOC(len,1)
  distance(i) = len(node_num(i))
enddo

end subroutine find_match_coord








subroutine find_match_node(tol, target_coord, coord, nn, node_num, max_n, n_nodes)

! Aug. 2015
! finding matched points based on the given tolerance
! 2D only
! only for one target point
! cf) subroutine 'search_nodes_along_line'
! Input:
!       target_coord: target coordinate to be compared with
!       coord: object coordinate to be referenced
!       max_n: maximum number of finding 'node_num'
! Output:
!       node_num: found point number(s) in a nodal system
!       n_nodes: number of found points

implicit none

integer, intent(in) :: nn, max_n 
real, intent(in) :: tol, target_coord(2), coord(2,nn)
integer, intent(out) :: node_num(max_n), n_nodes

integer :: j
real :: DNRM2               ! BLAS Library
real :: vector(2,nn), len(nn)
logical :: found(nn)


found  = .FALSE.

do j = 1, nn
  vector(:,j) = target_coord(:) - coord(:,j)
  len(j) = DNRM2(2, vector(:,j), 1)
  node_num(j) = j
  if (len(j) < tol) found(j) = .TRUE.
end do

n_nodes = COUNT(found)
node_num(1:n_nodes) = PACK(node_num, found)


end subroutine find_match_node









subroutine ranking_by_value(flag, n, values, order)

! Based on Selection Sort Algorithm (for small data sets less than 1000)
! flag = 1: give ranking by value (lower ranking order number for smaller value)

implicit none

integer, intent(in) :: flag, n
real,intent(in) :: values(n)
integer,intent(out) :: order(n)

integer :: i, j, addr
logical :: mask(n)

mask = .TRUE. ! unoccupied now

do i = 1, n
  addr = MINLOC(values, 1, mask)
  order(i) = addr
  mask(addr) = .FALSE.
end do

end subroutine ranking_by_value






!----------------- vector computation ------------------------------------------------------------------------


subroutine normalize(u, dim, v, len)

! Normalize the given vector u(dim)

implicit none

integer, intent(in) :: dim
real,intent(in) :: u(dim)
real,intent(out) :: v(dim), len
real :: DNRM2

len = DNRM2(dim,u,1)
if (len <= 0.0) STOP 'normalize: length is zero!'
v = u/len

end subroutine normalize






subroutine cross_product(a, b, c)

! c = a X b

implicit none

real, intent(in) :: a(3), b(3)
real, intent(out) :: c(3)

c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = a(3)*b(1) - a(1)*b(3)
c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine cross_product








subroutine set_line(point1, point2, line)

implicit none

real, intent(in) :: point1(2), point2(2)
real, intent(out) :: line(3)

line(1) = point2(2) - point1(2)  ! y2 - y1 
line(2) = point2(1) - point1(1)  ! x2 - x1 
line(3) = line(2)*point1(2) - line(1)*point1(1)
line(2) = -line(2)

end subroutine set_line







real function calc_dist(line , point)

implicit none

real, intent(in) :: line(3) , point(2)

calc_dist = abs(line(1)*point(1)+line(2)*point(2)+line(3)) / SQRT(ABS(line(1)**2+line(2)**2))

end function calc_dist








subroutine projection_s(point_a, point_b, test, s)

! distance from point_a to 'test'

implicit none

real, intent(in) :: point_a(2), point_b(2), test(2)
real, intent(out) :: s
real :: vector1(2), vector2(2), norm

vector1 = point_b - point_a
vector2 = test - point_a
norm = DOT_PRODUCT(vector1, vector1)
if (norm == 0.0) STOP 'projection_s: zero norm denominator!'

s = DOT_PRODUCT(vector1, vector2)/norm

end subroutine projection_s








subroutine search_nodes_along_line(ngroup, ndim, points, nn, node_coord, found_flag, group_number)

! multi line-group search over ngroup
! cf) subroutine 'set_cbound' in 'physical_domain'

implicit none

integer, intent(in) :: ngroup, ndim, nn
real, intent(in) :: points(ndim,2,ngroup), node_coord(ndim,nn)
logical, intent(out) :: found_flag(nn)
integer, intent(out) :: group_number(nn)

integer ::  i, n
real :: test_point(ndim), vec1(ndim), vec2(ndim), v1(ndim), v2(ndim)
real :: theta, norm1, norm2, DNRM2
real, parameter :: PI = 3.141592653589793238462643383279502884197169399375
real, parameter :: TOL = 1.0e-5


group_number = 0
found_flag = .FALSE.

do n = 1, ngroup
  vec1 = points(:,1,n)
  vec2 = points(:,2,n)

  do i = 1, nn  ! loop over given nodes

    test_point = node_coord(:,i)
    v1 = vec1 - test_point
    v2 = vec2 - test_point
    norm1 = DNRM2(ndim,v1,1)
    norm2 = DNRM2(ndim,v2,1)

    if ((norm1 <= TOL*norm2) .OR. (norm2 <= TOL*norm1)) then
      found_flag(i) = .TRUE.
      group_number(i) = n
    else
      call vector_angle(v1, v2, ndim, theta)

      if (ABS(theta - PI) < TOL) then
        found_flag(i) = .TRUE.
        group_number(i) = n
      endif
    endif
  end do ! i
end do ! n

end subroutine search_nodes_along_line







subroutine face_normal_vector(node_coord, normal)

! for 4-node surface element only
! 3D version
! using node 2-1 vector to node 4-1 vector (CCW +)
! Input:
!         node_coord(3,node_order): coordinates of each node of one 4-node surface plane
! Output:
!         normal(3): normal (unit) vector


implicit none

real, intent(in) :: node_coord(3,4)
real, intent(out) :: normal(3)
real :: vec1(3), vec2(3), length
real :: DNRM2

vec1 = node_coord(:,2) - node_coord(:,1)
vec2 = node_coord(:,4) - node_coord(:,1)

call cross_product(vec1, vec2, normal)

length = DNRM2(3,normal,1)
if (length <= 0.0) STOP 'face_normal_vector: vector length is zero!'
normal = normal/length

end subroutine face_normal_vector







subroutine vector_angle(v1, v2, ndim, theta)

! theta: angle from v1 to v2 (CCW +)
! 2D and 3D

implicit none

integer, intent(in) :: ndim
real, intent(in) :: v1(ndim), v2(ndim)
real, intent(out) :: theta
real :: u1(ndim), u2(ndim), u_normal(ndim), cos_normal, vector(ndim)
real, parameter :: PI = 3.141592653589793238462643383279502884197169399375
real :: norm, DNRM2


norm = DNRM2(ndim, v1, 1)
! if (norm <= 0) STOP 'vector_angle: invalid v1 norm denominator!'
if (norm <= 0) norm = 1.0
u1 = v1/norm

norm = DNRM2(ndim, v2, 1)
! if (norm <= 0) STOP 'vector_angle: invalid v2 norm denominator!'
if (norm <= 0) norm = 1.0
u2 = v2/norm


if (ndim == 2) then

  u_normal(1) = -u1(2)    ! u_normal: normal vector to u1 in plane
  u_normal(2) = u1(1)
  norm = u1(1)*u2(1) + u1(2)*u2(2)

  if (norm >= 1.0) then ! to prevent numerical-zero error
    theta = 0.0
  elseif (norm <= -1.0) then
    theta = PI
  else
    theta = ACOS(norm)
  endif

  cos_normal= u_normal(1)*u2(1) + u_normal(2)*u2(2)
  if (cos_normal < 0.0) then
    theta = 2.0*PI - theta
  endif

elseif (ndim == 3) then

  call cross_product(u1, u2, vector)
  vector = vector/DNRM2(ndim, vector, 1)
  call cross_product(vector, u1, u_normal)    ! u_normal: normal vector to u1 in plane (u1 x u2)
  norm = u1(1)*u2(1) + u1(2)*u2(2) + u1(3)*u2(3)

  if (norm >= 1.0) then ! to prevent numerical-zero error
    theta = 0.0
  elseif (norm <= -1.0) then
    theta = PI
  else
    theta = ACOS(norm)
  endif

  cos_normal= u_normal(1)*u2(1) + u_normal(2)*u2(2) + u_normal(3)*u2(3)
  if (cos_normal < 0.0) then
    theta = 2.0*PI - theta
  endif

else
  STOP 'vector_angle: vector dimension must be 2 or 3!'
endif

end subroutine vector_angle






!----------------- surface geometry computation --------------------------------------------------------------


subroutine surface_trans_matrix(corners, trans)

! Coded by Jeeho Lee (July, 2015)
!
! 3D
! construct transformation matrix based on 4 corner points
!
! Input:  
!         corners:
! Output:
!         trans:

implicit none

real, intent(in) :: corners(3,4)
real,intent(out) :: trans(3,3)

real :: v1(3), v2(3), normal(3)
real :: DNRM2


v1 = corners(:,2) - corners(:,1)
v1 = v1/DNRM2(3,v1,1)
v2 = corners(:,4) - corners(:,1)
v2 = v2/DNRM2(3,v2,1)

call cross_product(v1, v2, normal)  ! v1 and v2 are not necessarily perpendicular each other
normal = normal/DNRM2(3, normal, 1)
call cross_product(normal, v1, v2)  ! trans(2,:) is automatically unit vector

trans(1,:) = v1
trans(2,:) = v2
trans(3,:) = normal

end subroutine surface_trans_matrix











subroutine search_contact_faces(flag, nface, face_ix, nn, coord, corners, angle, in_or_out)

! Coded by Jeeho Lee (July, 2015)
!
! 3D
! check if the face is inside or outside of RECTANGULAR represented by 4 corner points
!
! Input:  
!         flag = 0: check if the face is on the test plane (corners) without direction check
!         nface: number of faces
!         face_ix: node number of faces
!         coord(i,j): coordinates(i=1:3) of j_th node
!         nn: number of nodes
!         corners: coordiantes of 4 corner points (CCW +)
!         angle: angle in degree for direction check
! Output:
!         in_or_out = TRUE(in) or FALSE(out) result for each face

implicit none

integer, intent(in) :: flag, nface, nn, face_ix(4,nface)
real, intent(in) :: coord(3,nn), corners(3,4), angle
logical, intent(out) :: in_or_out(nface)

integer :: i,j
real :: origin(3), v1(3), v2(3), v3(3), normal(3), vector(3), new_corners(3,4)
real :: box_21(2), box_34(2), point(2), vec2(2), theta, theta2, trans(3,3)
real :: DNRM2, norm, norm1, norm2, depth_tolerance, angle_cos
real, parameter :: PI = 3.141592653589793238462643383279502884197169399375
logical :: on_Plane_problem = .FALSE.


if (flag == 0) then
  on_Plane_problem = .TRUE. ! find points excatly on the plane
else
  on_Plane_problem = .FALSE.
  angle_cos = COS(angle*PI/180.0)
endif

call surface_trans_matrix(corners, trans)

! write(*,'(A/,3(3F10.3/))') 'trans=', trans

new_corners = MATMUL(trans,corners) ! new_corners: corners in the new plane (its own plane)

if (.NOT. on_Plane_problem) then ! do not need to check when flag=0
  ! check all corner points are in the same plane
  vector = new_corners(:,3) - new_corners(:,1)
  norm = DNRM2(3,vector,1)
  if (ABS(vector(3)) > 0.01*norm) then

    write(*,'(A)') 'search_contact_faces: four corner points are NOT in the same plane!'
    write(*,'(A/, 4(3E11.3/))') '   old corners = ', corners
    write(*,'(A/, 4(3E11.3/))') '   new_corners = ', new_corners
    STOP
  endif
endif

box_21 = new_corners(1:2,2) - new_corners(1:2,1)
norm1 = DNRM2(2,box_21,1)
box_21 = box_21/norm1

box_34 = new_corners(1:2,4) - new_corners(1:2,3)
norm2 = DNRM2(2,box_34,1)
box_34 = box_34/norm2

depth_tolerance = 1.0e-5*(norm1+norm2)

in_or_out = .FALSE.

if (on_Plane_problem) then    ! assume rectangular corner points

 do i = 1, nface
    v1 = 0.5*(coord(:,face_ix(3,i)) + coord(:,face_ix(1,i)))  ! mid-point
    v1 = v1 - corners(:,1)
    v1 = MATMUL(trans,v1)

    if (ABS(v1(3)) < depth_tolerance) then
      point = v1(1:2) ! on-the-plane point for in/out test
      vec2 = point - new_corners(1:2,1)
      vec2 = vec2/DNRM2(2,vec2,1)
      call vector_angle(box_21, vec2, 2, theta)

      if (theta < 0.5*PI) then
        vec2 = point - new_corners(1:2,3)
        vec2 = vec2/DNRM2(2,vec2,1)
        call vector_angle(box_34, vec2, 2, theta)
        if (theta < 0.5*PI) in_or_out(i) = .TRUE.
      endif
    endif
  end do

else
  normal = trans(3,:)

  write(*,'(A,4F12.3)') '* normal:', angle_cos, normal

  do i = 1, nface
    origin = coord(:,face_ix(1,i))
    v1 = coord(:,face_ix(2,i)) - origin
    v1 = v1/DNRM2(3,v1,1)
    v2 = coord(:,face_ix(4,i)) - origin
    v2 = v2/DNRM2(3,v2,1)

    call cross_product(v1, v2, v3)
    v3 = v3/DNRM2(3,v3,1)
    norm = v3(1)*normal(1) + v3(2)*normal(2) + v3(3)*normal(3)  ! inner product

    if (norm > angle_cos) then  ! direction check
      v1 = 0.5*(coord(:,face_ix(3,i)) + coord(:,face_ix(1,i)))  ! mid-point
!      v1 = v1 - origin
      v1 = MATMUL(trans,v1)

      point = v1(1:2) ! projected point for in/out test
      vec2 = point - new_corners(1:2,1)
      vec2 = vec2/DNRM2(2,vec2,1)
      call vector_angle(box_21, vec2, 2, theta)

      theta2 = 0.0
      if (theta < 0.5*PI) then
        vec2 = point - new_corners(1:2,3)
        vec2 = vec2/DNRM2(2,vec2,1)
        call vector_angle(box_34, vec2, 2, theta2)
        if (theta2 < 0.5*PI) in_or_out(i) = .TRUE.
      endif
    endif

!    write(*,'(A,I3,7F9.3,L3)') 'angle check:', i, norm, v3, (coord(j,face_ix(3,i)), j=1,3), in_or_out(i)

  end do  ! i

endif


end subroutine search_contact_faces









subroutine point_position(rn, corners, trans, point, s, nodal_areas, in_or_out)

! July 2015
! 3D version
! find a position factor s(1:4) on a 4-node face

implicit none

integer, intent(in) :: rn
real, intent(in) :: corners(3,4), trans(3,3), point(3)
real, intent(out) :: s(4), nodal_areas(4)
logical, intent(out) :: in_or_out

integer :: i, j, number
real :: origin(3), new_corners(3,4), v1(3), normal(3)
real :: point_2d(2), theta(4), vector(2,5), norm(5), s0(5), area(0:6), sum_area
real :: DNRM2, norm1, tol, depth_tolerance
real, parameter :: PI = 3.141592653589793238462643383279502884197169399375

s = 0.0
in_or_out = .FALSE.

new_corners = MATMUL(trans,corners)
! write(*,'(A/,4(3E12.3/))') 'new_corners:', (new_corners(:,j), j=1,4)

v1 = point - corners(:,1)
v1 = MATMUL(trans,v1)
norm1 = DNRM2(12,new_corners,1)
depth_tolerance = 0.25*(norm1)
tol = 1.0e-5*norm1

if (ABS(v1(3)) < depth_tolerance) then
  in_or_out = .TRUE.

  point_2d = v1(1:2) + new_corners(1:2,1)  ! on-the-plane point for in/out test

  number = 0
  do i = 1, 4
    vector(:,i) = new_corners(1:2,i) - point_2d
    norm(i) = DNRM2(2, vector(:,i), 1)
    if (norm(i) < tol) then
      number = i
      EXIT  !--------------------------------------->>>
    endif
  end do


  if (number > 0) then  ! the tested point is at the i-th corner
    s(number) = 1.0
  else
    vector(:,5) = vector(:,1)
    norm(5) = norm(1)

    do i = 1, 4
      call vector_angle(vector(:,i), vector(:,(i+1)), 2, theta(i))

      if (theta(i) > PI) then
        in_or_out = .FALSE.
        EXIT  !--------------------------------------->>>
      endif
    end do

    if (in_or_out) then
      do i = 1, 4
        area(i) = 0.5*norm(i)*norm(i+1)*SIN(theta(i))
      end do
      area(5) = area(1)
      area(6) = area(2)
      area(0) = area(4)

      number = 0
      do i = 1, 4
        if (area(i) < tol) then
          number = i
          EXIT  !--------------------------------------->>>
        endif
      end do

      if (number > 0) then  ! the tested point lies on the edge line
        i = number
        sum_area = area(i+1) + area(i-1)
        s0 = 0.0
        s0(i)   = area(i+1)/sum_area
        s0(i+1) = area(i-1)/sum_area
        s = s0(1:4)
        s(1) = s0(1) + s0(5)
      else                  ! the tested point lies at the inside
        sum_area = area(1)*area(2) + area(2)*area(3) + area(3)*area(4) + area(4)*area(1)
        sum_area = 1.0/sum_area
        do i = 1, 4
          s(i) = sum_area*(area(i+1)*area(i+2))
        end do
      endif

    endif ! in_or_out
  endif ! numbat at the i-th corner

else
  STOP 'point_position: depth exceeds depth_tolerance!'
endif ! depth_tolerance

do i = 1, 4
  nodal_areas(i) = 0.5*(area(i) + area(i-1))
end do


end subroutine point_position








subroutine mesh_surface(ndim, nn, nodes, coord, ne, nodes_elem, elements, start, num_bound, bound_nodes) 

! Coded by Jeeho Lee (Dec. 18, 2012)
! Element surface tracing algorithm (by Jeeho Lee)
! CCW element connectivity is assumed
! METIS 4.x is required
!   nodes: array for number of elements using each node
!   start: MUST be internal node number, if less than 1, default start point

implicit none

integer, intent(in) :: ndim, nn, ne, nodes_elem, start
integer, intent(in) :: nodes(nn), elements(nodes_elem,ne)
real, intent(in) :: coord(ndim,nn)
integer, intent(out) :: num_bound
integer, intent(out) :: bound_nodes(nn)

integer :: i, j, k, counter, node_num, current_node, next_node, previous_node
integer :: mesh_type, numflag, egg, neigh_size, num_test_elem
integer :: adj(nn+1), adjn(4*nn)
integer, allocatable :: boundary(:), neigh_nodes(:)
real :: v0(ndim), v1(ndim), theta0, theta
logical :: boundary_flag(nn), egg_flag
logical, allocatable ::mask_neighbor(:)

boundary_flag = .FALSE.
numflag = 1   ! Fortan format
mesh_type = 4

call METIS_MeshTONodal(ne, nn, elements, mesh_type, numflag, adj, adjn)

num_bound = 0

do i = 1, nn  ! over all nodes
  k = adj(i+1) - adj(i)
  if (k > nodes(i)) then
    boundary_flag(i) = .TRUE.
    num_bound = num_bound + 1
  endif
end do

allocate(boundary(num_bound))

counter = 0
do i = 1, nn  ! over all nodes
  if (boundary_flag(i)) then
    counter = counter + 1
    boundary(counter) = i
  endif
  if (counter == num_bound) EXIT
end do

!-----------------------------------------------------------------------
! starting node setting
current_node = boundary(1)  ! default starting point
if (start > 0) then
  do i = 1, num_bound
    if (start == boundary(i)) then
      current_node = start
      EXIT
    endif
  end do
endif

bound_nodes(1) = current_node
write(*,*) '*********** current_node:', current_node
!-----------------------------------------------------------------------

neigh_size = adj(current_node+1) - adj(current_node) ! number of neighbor (adjacent) nodes
allocate(neigh_nodes(neigh_size))
neigh_nodes = adjn(adj(current_node):(adj(current_node+1)-1))
allocate(mask_neighbor(neigh_size))
mask_neighbor = .TRUE.

num_test_elem = nodes(current_node)

counter = 0
finder: do j = 1, ne  ! over all elements
  do k = 1, nodes_elem
    node_num = elements(k,j)    ! CCW connectivity is assumed
    if (node_num == current_node) then
      counter = counter + 1
      if (k == 1) then  ! clockwise search to exclude that node
        previous_node = elements(nodes_elem,j)
      else
        previous_node = elements(k-1,j)
      endif
      do i = 1, neigh_size
        if (previous_node == neigh_nodes(i)) mask_neighbor(i) = .FALSE.
      end do      
    endif
    if (counter == num_test_elem) EXIT finder
  end do  ! k
end do finder ! j

if (COUNT(mask_neighbor) /= 1) stop 'mesh_surface: found more/less than one next_node!'
  
do i = 1, neigh_size
  if (mask_neighbor(i)) then
    next_node = neigh_nodes(i)
    EXIT
  endif
end do

write(*,*) '*********** next_node:', next_node

deallocate(neigh_nodes)
deallocate(mask_neighbor)

bound_nodes(2) = next_node 

!-----------------------------------------------------------------------

theta0 = 0.0
theta = 0.0
    
counter = 2
do
  previous_node = current_node
  current_node = next_node  
  j = adj(current_node+1) - adj(current_node)
  egg = 0
  egg_flag = .FALSE.
  do k = 1, j
    node_num = adjn(adj(current_node)+k-1)
    if (boundary_flag(node_num).AND.(node_num /= previous_node)) then
      if (egg == 0) then
        egg = node_num
        egg_flag = .TRUE.

      elseif (egg_flag) then  ! when there are two surface points
        v0 = coord(:,previous_node) - coord(:,current_node)
        v1 = coord(:,egg) - coord(:,current_node)
        call vector_angle(v0, v1, ndim, theta0)
        v1 = coord(:,node_num) - coord(:,current_node)
        call vector_angle(v0, v1, ndim, theta)
        if (theta < theta0) then
          theta0 = theta
          egg = node_num
        endif
        egg_flag = .FALSE.
      else  ! when there are more than two surface points
        v1 = coord(:,node_num) - coord(:,current_node)
        call vector_angle(v0, v1, ndim, theta)
        if (theta < theta0) then
          theta0 = theta
          egg = node_num
        endif
      endif
    endif
  end do
  next_node = egg
  counter = counter + 1
  bound_nodes(counter) = next_node

  if (counter == num_bound) EXIT
end do    

end subroutine mesh_surface








subroutine point_to_surface(point, p0, p_f, p_b, ndim, indicator, convexity, zone, projection)

! 2D only
! check 'point' with p0 (center point), p_f(forward pt), p_b(backward pt): CCW(+)
! indicator = .true.  : inside
!           = .false. : outside
! convexity: TRUE if convex
! projection: projection height & length
!             (1): height (negative if inside) (2): length (forward: +, backward: -)

implicit none

integer, intent(in) :: ndim
real, intent(in) :: point(ndim), p0(ndim), p_f(ndim), p_b(ndim)
logical, intent(out) :: indicator, convexity
integer, intent(out) :: zone
real, intent(out) :: projection(2)
real :: v0(ndim), v1(ndim), v2(ndim), theta1, theta2, length, PI_over_two, tester
real, parameter :: PI = 3.141592653589793238462643383279502884197169399375
real :: DNRM2

if (ndim /= 2) stop 'vector_angle: currently 2D only!'

PI_over_two = 0.5*PI

v0 = point - p0
v1 = p_b - p0
v2 = p_f - p0

call vector_angle(v1, v2, ndim, theta1)

call vector_angle(v1, v0, ndim, theta2)

zone = 0
length = DNRM2(ndim, v0, 1)

if (theta2 > theta1) then
  indicator = .TRUE.    ! inside
  
  if (theta1 > PI) then
    convexity = .TRUE.    ! convex
    if (theta2 > theta1 + PI_over_two) then
      zone = -1
      projection(1) = length*SIN(theta2)
      projection(2) = -length*COS(theta2)
    elseif (theta2 < 3.0*PI_over_two) then
      zone = 1
      theta2 = theta1 - theta2
      projection(1) = length*SIN(theta2)
      projection(2) = length*COS(theta2)   
    else
      zone = 0
      projection(1) = length*SIN(theta2)
      tester = length*SIN(theta1 - theta2)
      if (projection(1) > tester) then  ! near backward zone(-1)
        projection(2) = -length*COS(theta2)
        zone = -1
      else                              ! near forward zone(+1)
        projection(1) = tester
        projection(2) = length*COS(theta1 - theta2)  
        zone = 1
      endif
    endif ! zones
  else
    convexity = .FALSE.   ! concave  
    if (theta2 > 3.0*PI_over_two) then
      zone = -1
      projection(1) = length*SIN(theta2)
      projection(2) = -length*COS(theta2)
    elseif (theta2 < theta1 + PI_over_two) then
      zone = 1
      theta2 = theta1 - theta2
      projection(1) = length*SIN(theta2)
      projection(2) = length*COS(theta2)   
    else
      zone = 0
      projection(1) = length
      projection(2) = 0.0
      if (theta2 > (0.5*theta1 + PI)) then  ! near backward zone(-1)
        zone = -1
      else
        zone = 1
      endif
    endif ! zones
  endif ! convex-concave
    
else
  indicator = .FALSE.   ! outside

  if (theta1 > PI) then
    convexity = .TRUE.    ! convex
    if (theta2 < PI_over_two) then
      zone = -1
      projection(1) = length*SIN(theta2)
      projection(2) = -length*COS(theta2)
    elseif (theta2 > theta1 - PI_over_two) then
      zone = 1
      theta2 = theta1 - theta2
      projection(1) = length*SIN(theta2)
      projection(2) = length*COS(theta2)   
    else
      zone = 0
      projection(1) = length
      projection(2) = 0.0
      if (theta2 < 0.5*theta1) then ! near backward zone(-1)
        zone = -1
      else
        zone = 1
      endif
    endif ! zones
  else
    convexity = .FALSE.   ! concave
    if (theta2 < theta1 - PI_over_two) then
      zone = -1
      projection(1) = length*SIN(theta2)
      projection(2) = -length*COS(theta2)
    elseif (theta2 > PI_over_two) then
      zone = 1
      theta2 = theta1 - theta2
      projection(1) = length*SIN(theta2)
      projection(2) = length*COS(theta2)   
    else
      zone = 0
      projection(1) = length*SIN(theta2)
      tester = length*SIN(theta1 - theta2)
      if (projection(1) < tester) then  ! near backward zone(-1)
        projection(2) = -length*COS(theta2)
        zone = -1
      else                              ! near forward zone(+1)
        projection(1) = tester
        projection(2) = length*COS(theta1 - theta2)  
        zone = 1
      endif
    endif ! zones
  endif ! convex-concave
endif ! inside-outside

end subroutine point_to_surface
