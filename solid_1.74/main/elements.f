module elements

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (October 28, 2006)
! Modification History:		December, 2006 (by Jeeho Lee)
!                         November, 2014

use physical_domain

! constant parameters used over element models -----------------------------
! integer, parameter :: num_str_output = 12  --> use 'max_str_output' in module 'physical_domain'


contains   !===================================================================


subroutine element_dof(rn, elem_n, ix, nodal_dof, num_nodes, elem_flag)

! last modification: Nov. 2014
! num_nodes = 2: spring elements
!           > 2: quad or triangular elements

implicit none

integer, intent(in) :: rn, elem_n, nodal_dof, num_nodes
logical, intent(in) :: elem_flag
integer, intent(out) :: ix(nodal_dof*num_nodes)
integer :: i, node, inode(nodal_dof,num_nodes)

if (elem_flag) then
  do i = 1, num_nodes
    node = region(rn)%element(3+i,elem_n)
    inode(1:nodal_dof,i) = region(rn)%node(4:nodal_dof+3,node)
  end do
else ! spring element
  do i = 1, 2
    node = region(rn)%spring(3+i,elem_n)
    inode(1:nodal_dof,i) = region(rn)%node(4:nodal_dof+3,node)
  end do
endif

ix = pack(inode, .true.)

end subroutine element_dof




subroutine element_coord(rn, elem_n, coord, dim, elem_flag)

! rn: region number
! elem_n: element number
! coord(coord#,node#): coordinate of each node of element

implicit none

integer, intent(in) :: rn, elem_n, dim
logical, intent(in) :: elem_flag
real, intent(out) :: coord(:,:)
integer :: i, j, n
integer, allocatable :: nodes(:)


if (elem_flag) then
  n = region(rn)%nodes_element
  allocate(nodes(n))
  nodes = region(rn)%element(4:3+n,elem_n)

  do i = 1, n
    coord(:,i) = region(rn)%node_coord(1:dim,nodes(i))
  end do

  deallocate(nodes)
else
  n = 2
  allocate(nodes(n))
  nodes = region(rn)%spring(4:3+n,elem_n)

  do i = 1, n
    coord(:,i) = region(rn)%node_coord(1:dim,nodes(i))
  end do

  deallocate(nodes)
endif

end subroutine element_coord




subroutine element_material(rn, ele_num, mat_num, mat_model_num, elem_flag)

implicit none

integer, intent(in) :: rn, ele_num
logical, intent(in) :: elem_flag
integer, intent(out) :: mat_num, mat_model_num


if (elem_flag) then
  mat_num = region(rn)%element(3,ele_num)
else
  mat_num = region(rn)%spring(3,ele_num)
endif
mat_model_num = pdomain_mat_db%mat_num_model_num(mat_num)

end subroutine element_material




subroutine element_area(rn, ele_num, area, elem_flag)

implicit none

integer, intent(in) :: rn, ele_num
logical, intent(in) :: elem_flag
real, intent(out) :: area


if (elem_flag) then
  area = region(rn)%element_area(ele_num)
else
  area = region(rn)%spring_area(ele_num)
endif

end subroutine element_area




subroutine element_density(rn, ele_num, density)

implicit none

integer, intent(in) :: rn, ele_num
real, intent(out) :: density


density = region(rn)%element_density(ele_num)

end subroutine element_density





subroutine element_consis_mass(rn, ele_num, cfac)

implicit none

integer, intent(in) :: rn, ele_num
real, intent(out) :: cfac


cfac = region(rn)%consis_mass_ratio(ele_num)

end subroutine element_consis_mass





subroutine element_damping(rn, ele_num, damping_factor)

implicit none

integer, intent(in) :: rn, ele_num
real, intent(out) :: damping_factor(2)


damping_factor(:) = region(rn)%damping_factor(ele_num,:)

end subroutine element_damping





subroutine element_body_force(rn, body)

implicit none

integer, intent(in) :: rn
real, intent(out) :: body(3)

body = region(rn)%body_force


end subroutine element_body_force





subroutine element_gr_accel(rn, gr_accl_data, gr_accl)

implicit none

integer, intent(in) :: rn
real, intent(in) :: gr_accl_data(3)
real, intent(out) :: gr_accl(3)
integer :: i

do i = 1, 3
  gr_accl(i) = region(rn)%grfactor(i)*gr_accl_data(i)
enddo

end subroutine element_gr_accel



!==============================================================================
end module elements
