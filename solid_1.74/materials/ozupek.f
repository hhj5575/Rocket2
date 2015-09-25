module ozupek_material

!     Ozupek Model
!
!     Jeeho Lee @ DGU
!     Last Modification: April, 2015
	
use materials

implicit none

integer, parameter :: max_size_d = 100, max_size_state = 100

contains
!==============================================================================
            
            
subroutine ozupek(action,nstep,reg_num,ele_num,there,mat_num,ta,tl0, delt,f,finv,detf,bb,cc,str,tang,ndim)

implicit none
	  
integer, intent(in) :: reg_num, ele_num, there, mat_num, action, nstep, ndim
real, intent(in) :: delt, ta, tl0, f(3,3), finv(3,3), bb(6), detf
real, intent(out) :: str(ndim), tang(ndim,ndim)
real, intent(inout) :: cc(6)   ! cc: Right Cauchy-Green in, Green strain out
	  
integer :: npara, size_d, size_state
real :: d_array(max_size_d), state(max_size_state)


select case (action) 
case(0)
  
	if (.NOT. material(mat_num)%flag) then
  
    npara = material(mat_num)%size_para
		call make_d_ozupek(npara,material(mat_num)%para,d_array,size_d)
		call set_d_array(1,mat_num,size_d,d_array)

    size_state = max_size_state ! max_size_state: same as prescribed by '*DEPVAR' in Abaqus

    call number_material_state_var(mat_num,size_state) 
    write(*,*) 'maximum size_state & size_state: ', max_size_state, size_state
    if (size_state > max_size_state) stop 'ozupek: size_state exceeds maximum setting!'
  else
    write(*,*) 'Warning!!! ozupek: material initialize flag is not FALSE!'
	endif 
		
case default
  size_d = material(mat_num)%size_d
	call set_d_array(2,mat_num,size_d,d_array)
		
! The subrountine in 'materials' to download history data
  size_state = material(mat_num)%num_row_state
	call material_history(1,mat_num,there,state(1:size_state),size_state)	  
	  
! Call material model
  call ozupek_wrapper(nstep,size_d,d_array,delt,ta,f,size_state,state,ndim,str,tang)

! The subrountine in 'materials' to upload history data
	call material_history(2,mat_num,there,state(1:size_state),size_state)
	  
end select
	  
end subroutine ozupek



	  	  
subroutine make_d_ozupek(nprops,props,d,size_d)

implicit none

integer, intent(in) :: nprops
real, intent(in) :: props(nprops)
real, intent(out) :: d(:)
integer, intent(out) :: size_d

size_d = nprops
d(1:size_d) = props(1:nprops)


end subroutine make_d_ozupek

            
            
            
            
            
subroutine ozupek_wrapper(nstep,nprops,props,dtime,temp,dfgrd1,nstatev,statev,ntens,stress,ddsdde)

!     Ozupek model wrapper to call Abaqus UMAT

implicit none

integer, intent(in) :: nstep, ntens, nprops, nstatev
real, intent(in) :: props(nprops), dtime, dfgrd1(3,3)
real, intent(out) :: stress(ntens), ddsdde(ntens,ntens)
real, intent(inout) :: statev(nstatev)
integer, parameter :: ndi = 3
integer :: nshr, kinc
integer :: noel, npt, layer, kspt, kstep
real :: sse, spd, scd, rpl ,ddsddt(ntens), drplde(ntens), drpldt, stran(ntens), dstran(ntens), time(2)
real :: dtemp, predef(1), dpred(1), coords(3),drot(3,3), pnewdt, celent, dfgrd0(3,3), temp
character (len=8) :: cmname


nshr = ntens - ndi
kinc = nstep + 1

call UMAT_01(stress,statev,ddsdde,sse,spd,scd, &
  rpl,ddsddt,drplde,drpldt, &
  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
  ndi,nshr,ntens,nstatev,props,nprops,coords,drot,pnewdt, &
  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)


end subroutine ozupek_wrapper


! module contains end -----------------------------------------------------	  
  	  
end module ozupek_material
