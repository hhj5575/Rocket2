subroutine material_initialize(action, mat_num, mat_model, ndim)


! Coded by Jeeho Lee (October 31, 2009)
! Modification History:
!        July 2015 (2D and 3D)

! External subroutine
! material_type numbering guide: 1~99: 

! material_type = 1: 
!			= 2: Plastic-damage Model(Lee & Fenves) - Small strain formulation
!			= 3: Viscoelastic Model - Small strain formulation


! action = 0: set up material array (d-array) & Create state variables (nonlinear cases)
!        = 1: assemble k & p (tangent stiffness matrix & residual vector)

! ndim = 4: 2D
!      = 6: 3D

include 'material_models_list'

implicit none

integer, intent(in) :: action, mat_num, mat_model, ndim

integer :: nothing = 0, nothing2 = 0  ! dummy values
real :: eps(ndim), vps(ndim), str(ndim), res(ndim), tang(ndim,ndim), secd(ndim,ndim), noreal
real :: delt, chleng, ta, tl0, dmg_var(2)
real :: f(3,3),finv(3,3),detf,bb(6),cc(6)

if (action == 0) then
  write(*,*) 'mat model #:  ',mat_model

  select case (mat_model)
  case(1)
    call linearelastic(action,nothing,mat_num,ta,tl0,eps,str,tang,ndim)
  case(2:3)
    call viscoelastic(action,nothing,nothing,mat_num,ta,tl0,delt,eps,str,tang,secd,ndim)	
  case(12)
    call mooney_rivlin_2(action,nothing,nothing,nothing2,nothing2,mat_num,ta,tl0,delt,f,finv,detf,bb,cc,str,tang,ndim)
  case(14:15)
    call mooney_rivlin_2(action,nothing,nothing,nothing2,nothing2,mat_num,ta,tl0,delt,f,finv,detf,bb,cc,str,tang,ndim)
  case(31)
    call nonlinear_uniaxial(action,nothing,nothing,nothing2,mat_num,ta,tl0,delt,noreal,eps,str,tang,secd,ndim)
  case(51)
    call ozupek(action,nothing,nothing,nothing2,nothing2,mat_num,ta,tl0,delt,f,finv,detf,bb,cc,str,tang,ndim)
!  case(21)
!    call cohesive1(action,nothing,nothing,mat_num,ta,tl0,delt,noreal,noreal,noreal,1)
  case default
    STOP 'material_initialize: material type number is invalid!'
  end select
endif

end subroutine material_initialize
