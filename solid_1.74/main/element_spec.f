module element_specification

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! Coded by Jeeho Lee (Nov 22, 2009)
! Modification History:	May 2015


use elements
use time_domain

implicit none
save

real :: thickness=1.0, body(3), density, consistent_ratio=0.0, damping_factor(2)=0.0
real :: dt, ctan(3), gr_accl(3)

contains
!==============================================================================


subroutine element_pre(action, elem_type, rn, ele_num, ul, xl, tl, tl0)

! called by: element_interface

implicit none

integer, intent(in) :: action, elem_type, rn, ele_num
real, intent(out) :: ul(:,:), xl(:,:), tl(:), tl0

integer :: num_elem_nodes, ndof, ndim
real :: temperature = 0.0


call element_area(rn, ele_num, thickness, .TRUE.)
call element_density(rn, ele_num, density)
call element_consis_mass(rn, ele_num, consistent_ratio)
call element_damping(rn, ele_num, damping_factor)
call element_body_force(rn, body)

ndim = region(rn)%num_dim
num_elem_nodes = region(rn)%nodes_element
ndof = region(rn)%num_node_dof

call fetch_xl(rn, ele_num, num_elem_nodes, xl, ndim, .TRUE.)
call fetch_ul(rn, ele_num, num_elem_nodes, ul, ndof, .TRUE.)

!----------------------------------------------------------

call get_temperature(time, temperature, tl0)          ! time: from 'time_domain' module
! call fetch_tl(rn, ele_num, tl, num_elem_nodes, ndof)
tl = temperature   ! Note: tl(nst) is a vector array; uniform temperature assumption

! Time related data
dt = delt    ! delt: from 'time_domain' module
ctan = intgrn_para

call element_gr_accel(rn, ground_accel_data, gr_accl)  ! ground_accel_data(1:3): 'time_domain' module
 
end subroutine element_pre






subroutine spring_pre(action, elem_type, rn, ele_num, ul, xl, tl, tl0)

! called by: element_interface
! for spring elements

implicit none

integer, intent(in) :: action, elem_type, rn, ele_num
real, intent(out) :: ul(:,:), xl(:,:), tl(:), tl0
integer :: num_elem_nodes, ndof, ndim
real :: temperature = 0.0

call element_area(rn, ele_num, thickness, .FALSE.)
density = 0.0
consistent_ratio = 0.0
body = 0.0

ndim = region(rn)%num_dim
num_elem_nodes = 2
ndof = region(rn)%num_node_dof
call fetch_xl(rn, ele_num, num_elem_nodes, xl, ndim, .FALSE.)
call fetch_ul(rn, ele_num, num_elem_nodes, ul, ndof, .FALSE.)
!----------------------------------------------------------

call get_temperature(time, temperature, tl0)          ! time: from 'time_domain' module
! call fetch_tl(rn, ele_num, tl, num_elem_nodes, ndof)
tl = temperature   ! Note: tl(nst) is a vector array; uniform temperature assumption

! Time related data
dt = delt    ! delt: from 'time_domain' module
ctan = intgrn_para
call element_gr_accel(rn, ground_accel_data, gr_accl)  ! ground_accel_data(1:3): 'time_domain' module
end subroutine spring_pre






subroutine fetch_xl(rn, ele_num, elem_nodes, xl, dim, elem_flag)

! dim: Dimension of a node

implicit none

integer, intent(in) ::  rn, ele_num, dim, elem_nodes
logical, intent(in) :: elem_flag
real, intent(out) :: xl(dim, elem_nodes)
real :: coord

call element_coord(rn, ele_num, xl, dim, elem_flag)

end subroutine fetch_xl





subroutine fetch_ul(rn, ele_num, num_elem_nodes, ul, ndof, elem_flag)

! num_elem_nodes = 2: spring elements
! ndof: Number of DOF per a node

implicit none

integer, intent(in) ::  rn, ele_num, ndof, num_elem_nodes
logical, intent(in) :: elem_flag
real, intent(out) :: ul(ndof, 5*num_elem_nodes)
integer :: nodes(num_elem_nodes)
integer :: ix(ndof*num_elem_nodes), n2, n3, n4, n5, shape(2)

n2 = 2*num_elem_nodes
n3 = 3*num_elem_nodes
n4 = 4*num_elem_nodes
n5 = 5*num_elem_nodes

call element_dof(rn, ele_num, ix, ndof, num_elem_nodes, elem_flag)

shape(1) = ndof
shape(2) = num_elem_nodes
ul(:,1:num_elem_nodes) = reshape(history(rn)%data(ix,1),shape)    ! disp vector at time n+1
ul(:,num_elem_nodes+1:n2) = ul(:,1:num_elem_nodes) - reshape(history(rn)%data(ix,4),shape) ! disp increment vector between t_n+1 & t_n

ul(:,n2+1:n3) = 0.0                      
ul(:,n3+1:n4) = reshape(history(rn)%data(ix,2),shape)   ! velocity vector
ul(:,n4+1:n5) = reshape(history(rn)%data(ix,3),shape)   ! accel vector

end subroutine fetch_ul





subroutine trans_ul(u, nix, trans)

implicit none

integer, intent(in) ::  nix
real, intent(in) :: trans(nix,nix)
real, intent(inout) :: u(nix, 5)

u = MATMUL(trans,u)   ! all disp, velocity, acceleration data


end subroutine trans_ul


!==============================================================================
end module element_specification