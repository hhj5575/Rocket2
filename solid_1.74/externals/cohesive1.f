module cohesive_material

!     Cohesive Material Model
!
	
use materials

implicit none

integer, parameter :: max_size_d = 10, max_size_state = 6  ! max_size_d >= size_d
integer, parameter :: pre_props = 2, num_internal_props = 0

contains
!==============================================================================
            
subroutine cohesive1(action,pro_type,there,mat_num,ta,tl0,delt,eps,str,dd,ndim)

implicit none
	  
integer, intent(in) :: pro_type, there, mat_num, action, ndim
real, intent(in) :: delt, ta, tl0, eps(2)
real, intent(out) :: str(2), dd(2)
	  
integer :: npara, size_d, size_state
real :: d_array(max_size_d), state(max_size_state)



if (ndim /= 1) STOP 'cohesive1: invalid ndim!'

select case (action) 
case(0)
  
	if (material(mat_num)%flag .EQV. .FALSE.) then
  
    npara = material(mat_num)%size_para

		call make_d_cohesive(npara,material(mat_num)%para,d_array,size_d)

		call set_d_array(1,mat_num,size_d,d_array)

    size_state = 6
    call number_material_state_var(mat_num,size_state)

	else
    write(*,*) 'Warning!!! cohesive: material initialize flag is not false!'
	endif 
		
case default

  size_d = material(mat_num)%size_d
	call set_d_array(2,mat_num,size_d,d_array)

		
! The subrountine in 'materials' to download history data
  size_state = 6
	call material_history(1,mat_num,there,state(1:size_state),size_state)

	  
	call cohesive1_mat(pro_type,size_d,d_array,ta,tl0,state,eps,str,dd)


! The subrountine in 'materials' to upload history data
	call material_history(2,mat_num,there,state(1:size_state),size_state)

	  
end select
	  
end subroutine cohesive1





	  	  
subroutine make_d_cohesive(nprops,props,d,size_d)

implicit none

integer, intent(in) :: nprops
real, intent(in) :: props(nprops)
real, intent(out) :: d(:)
integer, intent(out) :: size_d
		
	  
d(1) = props(1)   ! Young's modulus
d(2) = props(2)    ! Poisson's ratio
d(3) = props(3)    ! yield limit
d(4) = props(4)    ! isotropic hardening modulus( < 0 for softening)
d(5) = props(5)   ! thermal expansion coefficient
d(6) = props(6)   ! Young's modulus (2nd material)
d(7) = props(7)    ! Poisson's ratio
d(8) = props(8)    ! yield limit
d(9) = props(9)    ! isotropic hardening modulus( < 0 for softening)
d(10) = props(10)   ! thermal expansion coefficient

size_d = 10 ! pre_props + num_internal_props ! num_internal_props: none


end subroutine make_d_cohesive

          
            
              
                  
            
            
            
subroutine cohesive1_mat(pro_type,size_d,d,ta,ta0,state,eps,str,dd)

implicit none

integer, intent(in) :: pro_type, size_d ! not used here
real, intent(in) :: d(:), ta, ta0, eps(2)
real, intent(out) :: str(2), dd(2)
real, intent(inout) :: state(6)
logical :: noconv

real :: e, alpha, epst, lambda, alp, epn, epp, cc, kep

real :: gamma, nu, sy, h, xlen, ff, yld, sigtr, count,sig, dsig, dlam
real :: epstr, dyld, nn, rs, rf, det, etol


etol = 1.0d-15
e = d(1)
nu = d(2)
sy = d(3)
h  = d(4)   ! isotropic hardening
alpha = d(5)


dd(1) = e
epn = state(1)
epp = state(2)


if (state(3).le. 0 .or. state(3) .gt. 1.0 ) then
state(3) = 1.0
endif


!    compute trial stress

   sigtr = d(1)*(eps(1) - epn)


!     Check yield

ff    =  sigtr
yld   = d(3)*exp(d(4)*epp)
kep = e*h*yld/(e+h*yld)
! write(*,*)ff, yld

!     Plastic state

if (ff.gt.yld) then

lambda =  0.0d0
sig    =  sigtr
cc     =  1./d(1)
epstr  =  cc*sigtr



!         Compute Newton parameters

noconv = .true.
count =  0

do while (noconv)

count = count +1
ff    = abs(sig)
yld   = d(3)*exp(d(4)*epp)
dyld  = h*yld
nn    = (sig)/ff
Rs    = epstr - cc*sig - lambda*nn
rf    = yld - ff
det = 1.0/(dyld*cc + 1.0d0)

! increments
dsig  =  det*(dyld*Rs + nn*rf)
dlam  =  det*(  nn*Rs - cc*rf)



! update variables 
lambda = lambda + dlam
sig    =  sig + dsig

epn    = state(1)  + lambda*nn
epp    = state(2)  + lambda


!         Check convergence
if(abs(dlam).lt.etol*abs(lambda)) then   ! Newton converge
 if(abs(dsig).lt.etol*abs(sig))  then
 noconv = .false.
 endif
elseif(count.gt.25) then                     ! Count limit
 write(*,*) "no convergence!", lambda,dlam
 noconv = .false.
endif

enddo

!       Elasto-plastic modulus

dd(1) = kep

!       Update history variables

state(1) = epn
state(2) = epp
state(3) = 1.d0 - epn/eps(1)   ! damage


!     Elastic state

else

!       Damaged Elastic modulus
if (state(3).gt.0 .and. eps(1).ge.0) then
 dd(1) = d(1)*state(3)
 sig = dd(1)*eps(1)
!       elastic modulus
elseif (eps(1).ge.0) then
 dd(1) = d(1)
 sig = dd(1)*eps(1)
!       contact status
else
 dd(1) = d(1)*2.0d0   ! for compression (contact)
 sig = dd(1)*eps(1)
endif

!
endif
!

str(1) = sig


!------------------------------------------------
e = d(6)
nu = d(7)
sy = d(8)
h  = d(9)   ! isotropic hardening
alpha = d(10)


dd(2) = e
epn = state(4)
epp = state(5)


if (state(6).le. 0 .or. state(6) .gt. 1.0 ) then
state(6) = 1.0
endif


!    compute trial stress

sigtr = e*(eps(2) - epn)


!     Check yield

ff    =  sigtr
yld   = sy*exp(h*epp)
kep = e*h*yld/(e+h*yld)
!write(*,*)ff, yld

!     Plastic state

if (ff.gt.yld) then

lambda =  0.0d0
sig    =  sigtr
cc     =  1./e
epstr  =  cc*sigtr



!         Compute Newton parameters

noconv = .true.
count =  0

do while (noconv)

count = count +1
ff    = abs(sig)
yld   = sy*exp(h*epp)
dyld  = h*yld
nn    = (sig)/ff
Rs    = epstr - cc*sig - lambda*nn
rf    = yld - ff
det = 1.0/(dyld*cc + 1.0d0)

! increments
dsig  =  det*(dyld*Rs + nn*rf)
dlam  =  det*(  nn*Rs - cc*rf)



! update variables
lambda = lambda + dlam
sig    =  sig + dsig

epn    = state(4)  + lambda*nn
epp    = state(5)  + lambda


!         Check convergence
if(abs(dlam).lt.etol*abs(lambda)) then   ! Newton converge
if(abs(dsig).lt.etol*abs(sig))  then
noconv = .false.
endif
elseif(count.gt.25) then                     ! Count limit
write(*,*) "no convergence!", lambda,dlam
noconv = .false.
endif

enddo

!       Elasto-plastic modulus

dd(2) = kep

!       Update history variables

state(4) = epn
state(5) = epp
state(6) = 1.d0 - epn/eps(2)   ! damage


!     Elastic state

else

!       Damaged Elastic modulus
if (state(6).gt.0 .and. eps(2).ge.0) then
dd(2) = e*state(6)
sig = dd(2)*eps(2)
!       elastic modulus
elseif (eps(2).ge.0) then
dd(2) =  e
sig = dd(2)*eps(2)
!       contact status
else
dd(2) = e*2.0d0   ! for compression (contact)
sig = dd(2)*eps(2)
endif

!
endif
!

str(2) = sig


end subroutine cohesive1_mat




! module contains end -----------------------------------------------------

end module cohesive_material
