module linear_elastic_material

use materials

implicit none

integer, parameter :: max_size_d = 30

contains
!==============================================================================

subroutine linearelastic(action, pro_type, mat_num, ta, tl0, eps, str, dd, ndim)

implicit none


integer, intent(in) :: mat_num, action, ndim, pro_type
real, intent(in) :: ta, tl0
real, intent(out) :: str(ndim), dd(ndim,ndim)
real, intent(inout) :: eps(ndim)
	  
integer :: size_d
real :: d_array(max_size_d)


select case (action)
case(0)
	  		
	if (.NOT. material(mat_num)%flag) then
	
		call make_d(material(mat_num)%para,d_array)
    size_d = max_size_d
		call set_d_array(1,mat_num,size_d,d_array)
	endif 
  
  call number_material_state_var(mat_num,1)
  
case default
  size_d = material(mat_num)%size_d
	call set_d_array(2,mat_num,size_d,d_array)
		
	call le_2d(pro_type,d_array,ta,tl0,eps,str,dd)
	
end select


end subroutine linearelastic



	  	  
subroutine make_d(props,d)

!
!    Linear Elastic Model
!
!     Last Modification: Feb, 2010
!              By  Jeeho Lee
!
	  
implicit none

integer, parameter :: nprops = 3

real, intent(in) :: props(nprops)
real, intent(out) :: d(:)

real :: e, nu


d = 0.0	

d(1)  = props(1)    ! e
d(2) = props(2)    ! nu
d(5) = props(3)    ! thermal expansion coefficient

e = d(1)
nu = d(2)
d(6) = e/(3.d0*(1.d0 - 2.d0*nu))       ! K (bulk modulus)
       
!      d(21) = 1.d0/e
!      d(22) = -nu*d(21)
!      d(23) = 2.d0*(1.d0+nu)*d(21)
d(24) = e/(1.d0-nu*nu)
d(25) = nu*d(24)
d(26) = e/(2.d0+2.d0*nu)   ! G (mu)
! d(27) = e/(2.d0-2.d0*nu)
! d(28) = (1.d0-2.d0*nu)/(1.d0-nu)	
      
d(29) = e*nu/(1.d0-nu-2.d0*nu*nu)   ! lamda
d(30) = d(29) + 2.d0*d(26)

end subroutine make_d




subroutine le_2d(pro_type,d,ta,tl0,eps,str,dd)

implicit none

integer, intent(in) :: pro_type
real, intent(in) :: d(:), ta, tl0
real, intent(out) :: str(:), dd(:,:)
real, intent(inout) :: eps(:)

integer :: i
real :: strt, epst, K, alpha

alpha = d(5)     ! thermal expansion coefficient (line)

if (pro_type == 0) then
  dd(1,1) = d(1)
  epst = alpha*(ta - tl0)
  str(1) = dd(1,1)*(eps(1)-epst)
else

  K = d(6)

  if (pro_type == 3) then  ! 3D problem
    call eltang(1,d,dd)
  elseif (pro_type == 4) then ! Plane stress problem
    call eltang(2,d,dd)
  else                   ! Plane strain problem (& Axi-symmetric problem)
    call eltang(3,d,dd)
  endif

  str = matmul(dd,eps)

  epst = alpha*(ta - tl0)
  strt = 3.0*K*epst
  do i = 1,3
    str(i) = str(i) - strt
    eps(i) = eps(i) - epst   ! mechanical strain for post-processing
  end do ! i

  !write(*,*) 'eps : ', eps
  !write(*,*) 'str : ', str

endif


end subroutine le_2d




subroutine eltang(id,d,dd)

!   id = 1: 3D
!      = 2: Plane Stress
!      = 3: Plane Strain (& Axi-Symmetry)

implicit none

integer, intent(in) :: id
real, intent(in) :: d(:)
real, intent(out) :: dd(:,:)
integer :: i,j,dim
real :: diag1,diag2,offdia


if (id == 1) then
  dim = 6
else
  dim = 4
endif

dd = 0.0

if (id == 2) then

  diag1 = d(24)
  offdia = d(25)
  diag2 = d(26)

  dd(1,1) = diag1
  dd(2,2) = diag1
  dd(1,2) = offdia
  dd(2,1) = offdia
  dd(4,4) = diag2    ! tangent matrix format: use 4th row & column
  
else   ! 3D or Plane strain  (& Axi-Symmetry)
  diag1 = d(30)
  offdia = d(29)
  diag2 = d(26)
   
  do i = 1, 3
    dd(i,i) = diag1
  enddo

  do j = 4, dim
    dd(j,j) = diag2
  enddo
    
  dd(1,2) = offdia
  dd(1,3) = offdia
  dd(2,3) = offdia
  dd(2,1) = offdia
  dd(3,1) = offdia
  dd(3,2) = offdia

endif

end subroutine eltang


! module contains end ---------------------------------------------------------	  
  	  
end module linear_elastic_material
