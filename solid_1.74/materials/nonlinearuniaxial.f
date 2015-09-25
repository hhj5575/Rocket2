module nonlinear_uniaxial_material

use materials

implicit none

integer, parameter :: max_size_d = 30, max_size_state = 5

contains
!==============================================================================

subroutine nonlinear_uniaxial(action,iter,pro_type,there,mat_num,ta,tl0,delt,length,eps,str,dd,secd,ndim)

implicit none

integer, intent(in) :: mat_num, action, ndim, iter, pro_type, there
real, intent(in) :: ta, tl0, delt, length
real, intent(out) :: str(ndim), dd(ndim,ndim), secd(ndim,ndim)
real, intent(inout) :: eps(ndim)
	  
integer :: size_d, size_state
real :: d_array(max_size_d), state(max_size_state), state_now(max_size_state)


select case (action)
case(0)
	  		
	if (.NOT. material(mat_num)%flag) then
	
		call make_d_nl_uniaxial(material(mat_num)%para,d_array)
    size_d = max_size_d
		call set_d_array(1,mat_num,size_d,d_array)
	endif 

  size_state = 1

  call number_material_state_var(mat_num,size_state)
  
case default
  size_d = material(mat_num)%size_d
	call set_d_array(2,mat_num,size_d,d_array)

! The subrountine in 'materials' to download history data
  size_state = material(mat_num)%num_row_state
!	call material_history(1,mat_num,there,state(1:size_state),size_state)
	call material_history(0,mat_num,there,state_now(1:size_state),size_state)

	call contact_link(iter,pro_type,d_array,length,eps,state_now(1),str,dd)
  secd = 0.0

! The subrountine in 'materials' to upload history data
	call material_history(2,mat_num,there,state_now(1:size_state),size_state)

end select


end subroutine nonlinear_uniaxial



	  	  
subroutine make_d_nl_uniaxial(props,d)

	  
implicit none

integer, parameter :: nprops = 4

real, intent(in) :: props(nprops)
real, intent(out) :: d(:)


d = 0.0	

d(1) = props(1) ! e
d(2) = props(2) ! order
d(3) = props(3) ! link_contact type: transform to Integer when used
d(4) = props(4) ! clearance


end subroutine make_d_nl_uniaxial





subroutine contact_link(iter,pro_type,d,length,eps,state,str,dd)

implicit none

integer, intent(in) :: iter, pro_type
real, intent(in) :: d(:), length, eps(:)
real, intent(out) :: str(:), dd(:,:)
real, intent(inout) :: state(1)

real, parameter :: transient = 1.0e-8
integer :: link_contact_type
real :: e, order, clearance, dl, new_dl, strain, prev_strain


prev_strain = state(1)
strain = eps(1)
e = d(1)
order = d(2)
link_contact_type = d(3)
clearance = d(4)

if ((clearance > length) .OR. (link_contact_type == 0)) clearance = length

dl = strain*length

if (link_contact_type /= 0) then
  new_dl = dl + length - clearance
  strain = new_dl/length
endif

if (iter == 0) prev_strain = strain

if (prev_strain < 0.0) then
  e = e**order
  dd(1,1) = e
  str(1) = dd(1,1)*strain
!  write(*,*) '+++++++++++++++++++++++++++ contact +++++++++++++++++++++++++'
!  write(*,*) 'strain, eps =', strain, eps(1)
!  write(*,*) 'dd =', dd(1,1)
!  write(*,*) 'str =', str(1)
else
  order = 1.0
  dd(1,1) = transient*e*order*EXP(-order*strain)
  str(1) = transient*e*(1.0 - EXP(-order*strain))
endif

state(1) = prev_strain


end subroutine contact_link



! module contains end ---------------------------------------------------------	  
  	  
end module nonlinear_uniaxial_material
