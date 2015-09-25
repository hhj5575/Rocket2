subroutine str_tang_infn(action,iter,pro_type,there,mat_num,material_model,ndm,ndf,nel &
                      ,ta,tl0,delt,xl,ul,eps,sig,tang,secd,ndim,str_data,n_out)

!   compute the stress and strain at the point


include '../material_models_list'

implicit none

integer, intent(in) :: action, pro_type, there, mat_num, material_model,ndm, ndf, nel, iter, ndim, n_out
real, intent(in) ::  ta, tl0, delt, xl(ndm,*), ul(ndf,*), eps(ndim)
real, intent(out) :: sig(ndim), tang(ndim,ndim), secd(ndim,ndim), str_data(n_out)

integer :: j, ndim2
real :: dd(ndim,ndim), strain(ndim), dmg_var(2)


ndim2 = 2*ndim

str_data = 0.0
tang = 0.0
secd = 0.0
dmg_var = 0.0
strain = eps

select case (material_model)
case(1)
  call linearelastic(action,pro_type,mat_num,ta,tl0,strain,sig,dd,ndim)
case(2)
  call viscoelastic(action,1,there,mat_num,ta,tl0,delt,strain,sig,dd,secd,ndim)
case(3)
  call viscoelastic(action,2,there,mat_num,ta,tl0,delt,strain,sig,dd,secd,ndim)

case default
  write(*,*) 'Material model #: ', material_model
  STOP 'str_tang_infn: material model number is invalid!'
end select

tang(1:ndim,1:ndim) = dd

if (n_out >= ndim2) then
  str_data(1:ndim) = sig(1:ndim)
  str_data((ndim+1):ndim2) = strain(1:ndim)
  if ((ndim < 5) .AND. (n_out >= 10)) str_data(9:10) = dmg_var(1:2)
endif

end subroutine str_tang_infn
