subroutine cohesive1_2d(action,iter,pro_type,reg_num,spring_num,ul,xl,tl,tl0,s,p,ndf,ndm,nst,nel,nquad)


!   action = 0 : initialization stage
!          = 1 : computational stage (residuals, tangent stiffness)
!          = 2 : dynamic problems sloving stage

!    pro_type: plane strain = 1, axisymmetric = 2ß

use element_specification

! Global variables (previously COMMON block data) from module 'element_specification':
!
!  - Geometry: thickness, body(3), gr_accl(3), density, consistent_ratio
!  - Time related: dt, ctan(3)
!  - Ouput: num_str_output

implicit none

integer, intent(in) :: action, pro_type, reg_num, spring_num, ndf, ndm, nst, nel, nquad, iter
real, intent(in) :: xl(ndm,nel), ul(ndf,5*nel), tl(nel), tl0
real, intent(out) :: s(nst,nst), p(nst)

integer :: lint, mat_num, material_model, there, nodes(2)
real :: area, length, coord(ndm,nel), ta, sin_value, cos_value, ke(4,5), str(2), tang(2), tran(4,4), temp(4,4)
real :: DNRM2, vect(4),vnormal(2), forc(2),bl2(3,2),br2(3,2)
real :: db(3,2),bl(3,2),br(3,2),p1(ndf,ndm),p2(ndf,ndm), db2(3,2)
integer :: i, i1, ii, j1, jj, j

!--------------------------------------------------------------------------------------------

area = 1.0 ! thickness
vnormal(1) = 0   ! need to get it from the input file later. shear bond(0,1),
vnormal(2) = 1   ! normal (1,0)


call element_material(reg_num,spring_num, mat_num, material_model, .FALSE.)

! coord = xl + ul(1:ndm,1:nel)   
coord = xl(:,1:2)
vect(1) = vnormal(1)    !!coord(:,1) - coord(:,2)
vect(2) = vnormal(2)
length = DNRM2(2, vect, 1)
vect = vect/length
vect(3) = -vect(2)
vect(4) =  vect(1)

!write(*,*),vnormal(1),vnormal(2)
!write(*,*),vect(1),vect(2)

ta = 0.5*(tl(1) + tl(2))

!   Get the address in history db using the subrountine in 'physical_domain'
lint = 1
call find_point_number_spring(reg_num,spring_num,there)



call str_tang_chsv1_2d(action,iter,pro_type,there,mat_num,material_model,ndm,ndf,nel,ta,tl0,delt,xl,ul &
                         ,length,str,tang,vect)



! multiply tangent modulus by length and area
tang(1) = tang(1)*1.*1.
tang(2) = tang(2)*1.*1.


! compute strain-displacement matrix
do i = 1,ndm

 bl(i,1) =  -vect(i)  !(xl(i,1) - xl(i,2))*xlen2
 br(i,1) =  bl(i,1)
 bl2(i,1) = -vect(i+ndm)
 br2(i,1) = bl2(i,1)

!           Set remaining terms and multiply by elastic modulus

bl(i,2) = -bl(i,1)
br(i,2) = -br(i,1)
bl2(i,2) = -bl2(i,1)
br2(i,2) = -br2(i,1)

db(i,1) =  bl(i,1)*tang(1)
db(i,2) = -db(i,1)
db2(i,1) = bl2(i,1)*tang(2)
db2(i,2) = -db2(i,1)

end do ! i

! form internal force
forc(1) = str(1)*1.*1.
forc(2) = str(2)*1.*1.

do i = 1,ndm

p1(i,1) = - bl(i,1)*forc(1)
p1(i,2) = - bl(i,2)*forc(1)
p2(i,1) = -bl2(i,1)*forc(2)
p2(i,2) = -bl2(i,2)*forc(2)

end do ! i

p(1) = p1(1,1) + p2(1,1)
p(2) = p1(2,1) + p2(2,1)
p(3) = p1(1,2) + p2(1,2)
p(4) = p1(2,2) + p2(2,2)



i1 = 0
do ii = 1,2
 j1 = 0
 do jj = 1,2
  do i = 1,ndm
   do j = 1,ndm
     s(i+i1,j+j1) = db(i,ii)*br(j,jj) + db2(i,ii)*br2(j,jj)
   end do ! j
  end do ! i
  j1 = j1 + ndf
 end do ! jj
 i1 = i1 + ndf
end do ! ii

!---------------------------------------------------------------
nodes = region(reg_num)%spring(4:5,spring_num)

call surface_node_length(reg_num, 1, nodes(1), 0, length)

s = length*s
p = length*p
!---------------------------------------------------------------

end subroutine cohesive1_2d






subroutine str_tang_chsv1_2d(action,iter,pro_type,there,mat_num,material_model,ndm,ndf,nel,ta,tl0,delt,xl,ul &
                         ,length,sig,tang,vect)

!   compute the stress and strain at the point


use cohesive_material

implicit none

integer, parameter :: ndim = 1

integer, intent(in) :: action, pro_type, there, mat_num, material_model, ndm, ndf, nel, iter
real, intent(in) ::  xl(ndm,*), ul(ndf,*), ta, tl0, delt, length
real, intent(out) :: sig(ndm), tang(2)
integer :: i 

integer :: npara, size_d, size_state
real :: d_array(max_size_d), state(max_size_state)
real :: coord(2,2), elong_length, eps(2), dd(2), str(2), dx(2), du(2)
real :: DNRM2, vect(4),vn(4)

tang(1) = 0.0
tang(2) = 0.0
size_d = material(mat_num)%size_d
call set_d_array(2,mat_num,size_d,d_array)
! coord = xl(1:2,1:2) + ul(1:2,1:2)

vn(1:4) = vect(1:4)
eps(1) = 0.0d0
eps(2) = 0.0

do i = 1,ndm
du(i) = ul(i,2) - ul(i,1)
eps(1) = eps(1)  +  du(i)*vn(i)
eps(2) = eps(2)  +  du(i)*vn(i+ndm)
end do ! i
!eps(1) = abs(eps(1))
!eps(2) = abs(eps(2))


select case (material_model)
case(21)
  call cohesive1(action,pro_type,there,mat_num,ta,tl0,delt,eps,str,dd,ndim)

case default
  write(*,*) 'Material model #: ', material_model
  stop 'call str_tang_chsv1_2d: material model number is invalid!'
end select

tang(1) = dd(1)
tang(2) = dd(2)
sig(1) = str(1)
sig(2) = str(2)


end subroutine str_tang_chsv1_2d
			






subroutine material_initialize_chsv1_2d(action, mat_num, mat_model)


! action = 0: set up material array (d-array) & Create state variables (nonlinear cases)
!        = 1: assemble k & p (tangent stiffness matrix & residual vector)

use cohesive_material

implicit none

integer, parameter :: ndim = 1
integer, intent(in) :: action, mat_num, mat_model

integer :: nothing  ! dummy value
real :: noreal(2)  ! dummy value
real :: delt, ta, tl0

if (action == 0) then
  write(*,*) 'mat model #:  ',mat_model

  select case (mat_model)
  case(21)
    call cohesive1(action,nothing,nothing,mat_num,ta,tl0,delt,noreal,noreal,noreal,1)
  case default
    STOP 'material_initialize: material type number is invalid!'
  end select
endif

end subroutine material_initialize_chsv1_2d

