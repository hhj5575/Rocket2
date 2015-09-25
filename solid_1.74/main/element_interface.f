subroutine element_interface(action, iter, reg_num, ele_num, elem_type, pe, ke, ndf, nix, nel, nquad, ndim, n_out, elem_flag, trans, uel, trans_flag)

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! Coded by Jeeho Lee (October 28, 2006)
! Modification History:	December 05, 2006
!                       May 2015

! External subroutine
! elem_type numbering guide: 1~99: linear elements
!                            100~499: nonlinear elements


! elem_type 

! action = 0: Create state variables (nonlinear cases)
!        = 1: Assemble k & p
! ndof: # of DOF of one node
! nix: # of all DOF of one element
! nel: # of nodes of one element
! ndim: space dimension
! n_out: size of stress/strain output fields
! elem_flag = .FALSE. : spring (connector) element

use element_specification

implicit none

external :: cohesive1_2d

integer, intent(inout) :: action
integer, intent(in) :: reg_num, ele_num, elem_type, ndf, nix, nel, nquad, ndim, iter, n_out
logical, intent(in) :: elem_flag, trans_flag
real, intent(in) :: trans(nix,nix)
real, intent(out) :: pe(nix), ke(nix,nix), uel(ndf,nel)

real :: xl(ndim,nel), ul(ndf,5*nel), tl(nel), tl0, str_out(nquad,n_out)



if (elem_flag) then
  call element_pre(action, elem_type, reg_num, ele_num, ul, xl, tl, tl0)
else
  call spring_pre(action, elem_type, reg_num, ele_num, ul, xl, tl, tl0)
endif

if (trans_flag) then
  call trans_ul(ul, nix, trans)
endif

uel = ul(:,1:nel)    ! NOT input data. This is for disp output in the cartesian coordinate!

select case (elem_type)

case (112:113)   ! 2D/3D spring element
  call spring_link(action, iter, 0, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, ndf, ndim, nix, nel, nquad)
  str_out = 0.0

case (210)   ! Plane stress (pro_type = 4) for test
  call nl_elmt2d(action, iter, 4, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)
case (211)   ! Plane strain: Small strain/Standard formulation
  call nl_elmt2d(action, iter, 1, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)
case (212)   ! Axisymmetric: Small strain/Mixed formulation
  call nl_elmt2d(action, iter, 2, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)

case (213)   ! 3D: Small strain/Standard formulation (pro_type = 3)
  call nl_elmt3d(action, iter, 3, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)

case (221)   ! Plane strain: Small strain/Mixed formulation
  call nl_2d_mixed(action, iter, 1, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)
case (222)   ! Axisymmetric: Small strain/Mixed formulation
  call nl_2d_mixed(action, iter, 2, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)
  
case (311)   ! Plane strain: Large Deformation/Standard formulation
  call nl_2d_fn(action, iter, 1, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)
case (312)   ! Axisymmetric: Large Deformation/Standard formulation
  call nl_2d_fn(action, iter, 2, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)
  
case (313)   ! 3D: Large Deformation/Mixed formulation (pro_type = 3)
  call nl_3d_fn(action, iter, 3, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)

case (321)   ! Plane strain: Large Deformation/Mixed formulation
  call nl_2d_fn_mx(action, iter, 1, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)
case (322)   ! Axisymmetric: Large Deformation/Mixed formulation
  call nl_2d_fn_mx(action, iter, 2, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)

case (323)   ! 3D: Large Deformation/Mixed formulation (pro_type = 3)
  call nl_3d_fn_mx(action, iter, 3, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, str_out, ndf, ndim, nix, nel, nquad, n_out)

case (412)   ! Axisymmetric-cohesive truss element
  call cohesive1_2d(action, iter, 1, reg_num, ele_num, ul, xl ,tl, tl0, ke, pe, ndf, ndim, nix, nel, nquad)
  str_out = 0.0
case default
	STOP 'element_interface: element type number is invalid!'
end select


! for stress/strain output
if ((action > 0).AND.(elem_flag)) then
  call write_str_output(reg_num, ele_num, nquad, str_out, n_out)
endif


end subroutine element_interface
