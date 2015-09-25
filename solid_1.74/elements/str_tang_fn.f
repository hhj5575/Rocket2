subroutine str_tang_fn(action,iter,reg_num,ele_num,there,mat_num,material_model,ndm,ndf,nel &
                      ,ta,tl0,delt,xl,ul,f,finv,detf,sig,tang,ndim,str_data,n_out)

!   Compute the stress and tangent stiffness at the quadrature point
!   - Standard & Mixed method
!   - Plane strain & Axisymmetric & 3D cases

!     Inputs:
!       f(3,3)       - Deformation gradient
!       finv(3,3)    - Inverse deformation gradient
!       df(3,3)      - Incremental deformation gradient (Not used here)
!       detf         - Determinant of deformation gradient
!     Outputs
!       tang(6,6)    - Expanded material moduli for mixed computation
!       sig(ndim)      - Cauchy stress values at t_n+1


include '../material_models_list'

implicit none

integer, intent(in) :: action, reg_num, ele_num, there, mat_num, material_model,ndm, ndf, nel, ndim, iter, n_out
real, intent(in) ::  xl(ndm,*), ul(ndf,*), ta, tl0, delt, detf(*), f(3,3), finv(3,3)
real, intent(out) :: sig(ndim), tang(ndim,ndim), str_data(n_out)

integer :: ii, i, j, ndim2
real :: bb(6), cc(6)


ndim2 = 2*ndim

tang = 0.0
sig = 0.0

!     Compute Left Cauchy-Green deformation tensor

bb(1) = f(1,1)*f(1,1) + f(1,2)*f(1,2) + f(1,3)*f(1,3)
bb(2) = f(2,1)*f(2,1) + f(2,2)*f(2,2) + f(2,3)*f(2,3)
bb(3) = f(3,1)*f(3,1) + f(3,2)*f(3,2) + f(3,3)*f(3,3)
bb(4) = f(1,1)*f(2,1) + f(1,2)*f(2,2) + f(1,3)*f(2,3)
bb(5) = f(2,1)*f(3,1) + f(2,2)*f(3,2) + f(2,3)*f(3,3)
bb(6) = f(1,1)*f(3,1) + f(1,2)*f(3,2) + f(1,3)*f(3,3)

! write(*,*) 'Left Cauchy: ', bb

!     Compute Right Cauchy-Green deformation tensor

cc(1) = f(1,1)*f(1,1) + f(2,1)*f(2,1) + f(3,1)*f(3,1)
cc(2) = f(1,2)*f(1,2) + f(2,2)*f(2,2) + f(3,2)*f(3,2)
cc(3) = f(1,3)*f(1,3) + f(2,3)*f(2,3) + f(3,3)*f(3,3)
cc(4) = f(1,1)*f(1,2) + f(2,1)*f(2,2) + f(3,1)*f(3,2)
cc(5) = f(1,2)*f(1,3) + f(2,2)*f(2,3) + f(3,2)*f(3,3)
cc(6) = f(1,1)*f(1,3) + f(2,1)*f(2,3) + f(3,1)*f(3,3)

! write(*,*) 'Right Cauchy: ', cc

select case (material_model)
case(12)  ! mooney-rivlin material model without viscoelasticity(optional)

  call mooney_rivlin_2(action,1,reg_num,ele_num,there,mat_num,ta,tl0,delt,f,finv,detf(1),bb,cc,sig,tang,ndim)

case(14)  ! mooney-rivlin material model with deviatoric viscoelasticity(optional)

  call mooney_rivlin_2(action,3,reg_num,ele_num,there,mat_num,ta,tl0,delt,f,finv,detf(1),bb,cc,sig,tang,ndim)
		
case(15)  ! mooney-rivlin material model with volumetric & deviatoric viscoelasticity(optional)

  call mooney_rivlin_2(action,4,reg_num,ele_num,there,mat_num,ta,tl0,delt,f,finv,detf(1),bb,cc,sig,tang,ndim)

case(51)  ! Ozupek material model

  call ozupek(action,iter,reg_num,ele_num,there,mat_num,ta,tl0,delt,f,finv,detf(1),bb,cc,sig,tang,ndim)

case default
  write(*,*) 'Material model #: ', material_model
  stop 'str_tang_fn: material model number is invalid!'
end select


if (n_out >= ndim2) then
  str_data(1:ndim) = sig(1:ndim)
  str_data((ndim+1):ndim2) = cc(1:ndim)  ! Green strain tensor computed in each material model
endif

end subroutine str_tang_fn
