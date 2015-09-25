subroutine spring_link(action,iter,pro_type,reg_num,spring_num,ul,xl,tl,tl0,s,p,ndf,ndm,nst,nel,nquad)

! uniaxial spring link element in 2D or 3D space

!   action = 0 : initialization stage
!          = 1 : computational stage (residuals, tangent stiffness)
!          = 2 : dynamic problems sloving stage

!    pro_type: spring (link) in 2D/3D = 0

use element_specification

! Global variables from module 'elements' & 'element_specification':
!  - Parameters: Numerical_Zero
!  - Geometry: thickness, body(3), gr_accl(3), density, consistent_ratio
!  - Time related: dt, ctan(3)

implicit none

integer, intent(in) :: action, pro_type, reg_num, spring_num, ndf, ndm, nst, nel, nquad, iter
real, intent(in) :: xl(ndm,nel), ul(ndf,5*nel), tl(nel), tl0
real, intent(out) :: s(nst,nst), p(nst)

integer, parameter :: ndim = 1, max_lint = 1, max_nel = 2

integer :: mat_num, material_model, point_address
integer :: nea, nea2      
integer :: i, j, k, l, i1, j1, j2, lint

real :: sig(ndim), eps(ndim), dd(ndim,ndim), q(2), ke(2,2), secd(ndim,ndim)
      	  
real :: lfac, cfac, mass(nel,nel), vl(ndf,nel), accl(ndf,nel), ct1, ct2, ct3
real :: ug, infv(ndf), damfac_k, damfac_m, damp(nst,nst), m(nst,nst), factor, dvol
real :: area, rr, ta, length, elong_length, coord(ndm,nel), vect(ndm), cos_theta(ndm), trans(2,2*ndm)
real :: DNRM2, temp(2*ndm,2)

logical :: mat_damping = .FALSE., flgrm = .FALSE., dyna = .FALSE.
logical :: ther = .FALSE.

!--------------------------------------------------------------------------------------------


if (pro_type /= 0) STOP 'spring_link: pro_type is not for spring link!'
if (ndm /= ndf) STOP 'spring_link: ndm is not consistent with ndf!'

if (action == 2) then
	dyna = .TRUE.
else
	dyna = .FALSE.
endif

dyna = .FALSE.  ! presently dynamic option is off
ther = .FALSE.  ! presently thermal option is off

nea = 4*nel
nea2 = 3*nel

! vl = ul(1:ndf,nea2+1:nea)
vl = 0.0
s = 0.0

!     read data from module 'element_specification'
area = thickness
rr = density

accl = 0.0
infv = 0.0	  
	  
if (dyna) then
!     read data from module 'element_specification'
	cfac = consistent_ratio
	
	damfac_m = damping_factor(1)   ! mass matrix proportional damping factor	
	damfac_k = damping_factor(2)   ! stiffness matrix proportional damping factor
	
	accl = ul(1:ndf,nea+1:nea+nel)

  lfac = 1.d0 ! lumped mass only

!  Zero mass & damping matrix
  mass = 0.0
  damp = 0.0

	ct1 = area*ctan(1)
	ct2 = area*ctan(2)
	ct3 = area*ctan(3)

! Ground motion part (default: off)
	flgrm = .TRUE.
  if (flgrm) then
		infv(1:ndm) = gr_accl(1:ndm)  ! ndm = ndf
	else
		infv = 0.0
	endif

else
	ct1 = area
	ct2 = 0.0
	ct3 = 0.0
	mat_damping = .FALSE.

endif

call element_material(reg_num,spring_num, mat_num, material_model, .FALSE.)

coord = xl(:,1:2)   ! two-end nodes
vect = coord(:,2) - coord(:,1)
length = DNRM2(ndm, vect, 1)

! compute the strain
coord = xl(:,1:2) + ul(:,1:2)
vect = coord(:,1) - coord(:,2)
elong_length = DNRM2(ndm, vect, 1)

if (length > Numerical_Zero) then
  eps = elong_length/length - 1.0
else
  STOP 'spring_link: zero length!'
endif

if (ther) then
  ta = 0.5*(tl(1) + tl(2))
else
  ta = tl0
endif


! Get the address in history db using the subrountine in 'physical_domain'

call find_point_number_spring(reg_num,spring_num,point_address)

call str_tang_spring(action,iter,pro_type,point_address,mat_num,material_model,ndm,ndf,nel &
                    ,ta,tl0,delt,xl,ul,length,eps,sig,dd,secd,ndim)

factor = dd(ndim,ndim)*area/length ! component factor for tangent stiffness matrix
q(1) = -sig(ndim)*area        ! axial force for end-node 1
q(2) =  sig(ndim)*area        ! axial force for end-node 2

!  Compute accelerations of relative motion and ground motion
if (dyna) then
  dvol = 0.5*area*length
  j1 = 1
  j2 = ndf
  do j = 1, nel
    p(j1:j2) = dvol*(body(:) -rr*(infv(:) + accl(:,j)))
    j1 = j1 + ndf
    j2 = j2 + ndf
  end do
else
  p = 0.0
endif


cos_theta = vect/length

trans = 0.0
trans(1,1:ndm) = cos_theta
trans(2,ndm+1:2*ndm) = cos_theta

ke(1,:) = [factor, -factor]
ke(2,:) = [-factor, factor]

temp = TRANSPOSE(trans)
p = p + MATMUL(temp, q)
temp = MATMUL(temp, ke)

s = MATMUL(temp, trans)


end subroutine spring_link






subroutine str_tang_spring(action,iter,pro_type,there,mat_num,material_model,ndm,ndf,nel &
                      ,ta,tl0,delt,xl,ul,length,eps,sig,tang,secd,ndim)

!   compute the stress and strain at the point


include '../material_models_list'

implicit none

integer, intent(in) :: action, pro_type, there, mat_num, material_model,ndm, ndf, nel, iter, ndim
real, intent(in) ::  ta, tl0, delt, xl(ndm,*), ul(ndf,*), length, eps(ndim)
real, intent(out) :: sig(ndim), tang(ndim,ndim), secd(ndim,ndim)

real :: dd(ndim,ndim), strain(ndim)


tang = 0.0
secd = 0.0

strain = eps

select case (material_model)
case(1)
  call linearelastic(action,pro_type,mat_num,ta,tl0,strain,sig,dd,ndim)
!case(2)
!  call viscoelastic(action,1,there,mat_num,ta,tl0,delt,strain,sig,dd,secd,ndim)
!case(3)
!  call viscoelastic(action,2,there,mat_num,ta,tl0,delt,strain,sig,dd,secd,ndim)
case(31)
  call nonlinear_uniaxial(action,iter,pro_type,there,mat_num,ta,tl0,delt,length,strain,sig,dd,secd,ndim)
	 
case default
  write(*,*) 'Material model #: ', material_model
  STOP 'str_tang_spring: material model number is invalid!'
end select

tang(1:ndim,1:ndim) = dd


end subroutine str_tang_spring
