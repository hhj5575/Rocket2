module mnrv_2_material

!     MR2 Model
!
!     Last Modification: Feb, 2010
!                        May, 2013
!                        March, 2014
!                        June, 2015
	
use materials

implicit none

integer, parameter :: max_size_d = 50, max_size_state = 112  ! for maximum ndim = 6 & max 15-set visco data
integer, parameter :: max_mu_data_set = 15, pre_props = 10, num_internal_props = 2
real, parameter :: Numerical_zero = 1.0e-10, Nu_incomp = 0.49999999999

contains
!==============================================================================
            
            
subroutine mooney_rivlin_2(action,mv_flag,reg_num,ele_num,there,mat_num,ta,tl0, delt,f,finv,detf,bb,cc,str,tang,ndim)

implicit none
	  
integer, intent(in) :: reg_num, ele_num, there, mat_num, action, mv_flag, ndim
real, intent(in) :: delt, ta, tl0, f(3,3), finv(3,3), bb(6), detf
real, intent(out) :: str(ndim), tang(ndim,ndim)
real, intent(inout) :: cc(6)   ! cc: Right Cauchy-Green in, Green strain out
	  
integer :: npara, size_d, size_state
real :: d_array(max_size_d), state(max_size_state)


select case (action) 
case(0)
  
	if (.NOT. material(mat_num)%flag) then
  
    npara = material(mat_num)%size_para
		call make_d_mnrv2(npara,material(mat_num)%para,d_array,size_d)
		call set_d_array(1,mat_num,size_d,d_array)

    if (size_d > (pre_props + num_internal_props)) then
      size_state = (ndim+1) + (size_d - (pre_props + num_internal_props))*(ndim+1)/2
    else
      size_state = 1
    endif
    
    call number_material_state_var(mat_num,size_state) 
    write(*,*) 'maximum size_state & size_state: ', max_size_state, size_state
    if (size_state > max_size_state) STOP 'mnrv: size_state exceeds maximum setting!'
  else
    write(*,*) 'Warning!!! mooney_rivlin: material initialize flag is not false!'
	endif 
		
case default
  size_d = material(mat_num)%size_d
	call set_d_array(2,mat_num,size_d,d_array)
		
! The subrountine in 'materials' to download history data
  size_state = material(mat_num)%num_row_state
	call material_history(1,mat_num,there,state(1:size_state),size_state)	  
	  
! Call Mooney-Rivlin material model
  call general_mnrv(mv_flag,reg_num,ele_num,size_d,d_array,ta,tl0,delt,f,finv,detf,bb,cc,state,ndim,str,tang)
! The subrountine in 'materials' to upload history data
	call material_history(2,mat_num,there,state(1:size_state),size_state)
	  
end select
	  
end subroutine mooney_rivlin_2



	  	  
subroutine make_d_mnrv2(nprops,props,d,size_d)

implicit none

integer, intent(in) :: nprops
real, intent(in) :: props(nprops)
real, intent(out) :: d(:)
integer, intent(out) :: size_d

integer :: i, j, num_visco_para  
real, parameter :: prop_limit = 1.0		

		
if (pre_props > nprops) stop 'make_d_mnrv2: too few data (minimum 10)!' 
    
d = 0.0
	  
d(1) = props(1)    ! Young's modulus
d(2) = props(2)    ! Poisson's ratio
if (d(2) > Nu_incomp) d(2) = Nu_incomp  ! Incompressible case
d(3) = d(1)/(2.d0*(1.d0 + d(2)))        ! G modulus
d(4) = d(1)/(3.d0*(1.d0 - 2.d0*d(2)))   ! K modulus

d(5) = props(3)    ! thermal expansion coefficient
d(6) = props(4)    ! c factor (if c = 0: Neo-Hookean) for C10 & C01
d(7) = props(5)    ! C20
d(8) = props(6)    ! C11
d(9) = props(7)    ! C02

d(10) = props(8)   ! T_ref for WLF model
d(11) = props(9)   ! c1 in WLF
d(12) = props(10)  ! c2 in WLF

size_d = pre_props + num_internal_props ! num_internal_props: G & K

if (nprops > pre_props) then    ! MR with viscosity case
  j = nprops - pre_props  
  if ((j/2 > max_mu_data_set) .OR. (mod(j,2) == 1)) STOP 'make_d_mnrv2: visco pair error!'
  
  num_visco_para = 0
  j = size_d - pre_props
  do i = (pre_props+1), nprops, 2
    d(i+j) = props(i)
    d(i+j+1) = props(i+1)
    if (props(i) > prop_limit) d(i+j) = prop_limit
    num_visco_para = num_visco_para + 2
  end do
  size_d = size_d + num_visco_para
endif

end subroutine make_d_mnrv2

            
            
            
            
            
subroutine general_mnrv(mv_flag,reg_num,ele_num,size_d,d,ta,tl0,dt,f,finv,detf0,b1,cc,state,ntm,sig,aa)

!		mv_flag = 1: hyperelasticity without viscoelasticity
! 	mv_flag = 2: NOT used presently
!		mv_flag = 3: deviatoric viscoelasticity
!		mv_flag = 4: deviatoric & volumetric viscoelasticity

implicit  none

integer, intent(in) :: mv_flag, reg_num, ele_num, size_d, ntm
real, intent(in) :: d(:), ta, tl0, dt, f(3,3), finv(3,3), b1(6) , detf0
real, intent(out) :: sig(ntm), aa(ntm,ntm)
real, intent(inout) :: cc(6), state(*)  ! cc: Right Cauchy-Green in, Green strain out

logical :: incompressible_flag, visco_flag, vol_visco_flag
integer :: i,j,k,l, ii,jj, fn_type, nv, n_pointer
integer :: n, n_order, np, nq
real :: detf, detfi, j23, trbb3,trbb2, bdi, mub1,mub2,mub3
real :: nub1,nub2,nub3,nub4, i1b,i2b
real :: u, up, upp, press, c_factor
real :: pk2(6), i1_c(6), i2_c(6), bb(6),bd(6),b2(6),ad(6,6)
real :: K_bulk, G, nu, gg, dt_shift
real :: t_alpha, t_beta, J_th, J_el, K_incomp
real :: c20, c11, c02, c_two1, c_two2, del_i1b(6), del_i2b(6), del_j(6), dd_i1b(6,6), dd_i2b(6,6)
real :: tensor1(3,3), tensor2(3,3), tensor3(3,3)
real :: c_pq(0:2,0:2),D_p(1:3),fac1, fac2, fac3, fac4
real, parameter :: one3 = 1.0/3.0, two3 = 2.0/3.0, one6 = 1.0/6.0


if (detf0 < Numerical_zero) then
  write(*,'(/A,2I8,E15.5/)') 'general_mnrv: detf is negative! region, element, detf:', reg_num, ele_num, detf0
  STOP
endif

n_order = 2

nu = d(2)
G = d(3)
K_bulk = d(4)
t_alpha = d(5)   ! t_alpha: thermal expansion coeff.
c_factor = d(6)

c_pq = 0.0

nub1 = c_factor*G         ! 2*C01
mub1 = (1.0-c_factor)*G   ! 2*C10

c_pq(0,1) = 0.5*nub1    ! C01
c_pq(1,0) = 0.5*mub1    ! C10
c_pq(2,0) = d(7)
c_pq(1,1) = d(8)
c_pq(0,2) = d(9)

D_p(1) = 2.0/K_bulk
D_p(2) = D_p(1)
D_p(3) = D_p(1)

if ((size_d > (pre_props + num_internal_props)).AND.(dt > 0.0).AND.(mv_flag >= 3)) then
  visco_flag = .TRUE.
else
  visco_flag = .FALSE.
endif

if ((mv_flag == 4).AND.visco_flag) then
  vol_visco_flag = .TRUE.
else
  vol_visco_flag = .FALSE.
endif

incompressible_flag = .FALSE.
if ((nu + Numerical_zero) > Nu_incomp) then   ! incompressible case
  incompressible_flag = .TRUE.
endif


detf = detf0

detfi = 1.d0/detf
j23   = detfi**two3
do i = 1,6
  bb(i) = b1(i) * j23
end do


!     Compute bb:bb

b2(1) = bb(1)*bb(1) + bb(4)*bb(4) + bb(6)*bb(6)
b2(2) = bb(4)*bb(4) + bb(2)*bb(2) + bb(5)*bb(5)
b2(3) = bb(6)*bb(6) + bb(5)*bb(5) + bb(3)*bb(3)
b2(4) = bb(1)*bb(4) + bb(4)*bb(2) + bb(6)*bb(5)
b2(5) = bb(4)*bb(6) + bb(2)*bb(5) + bb(5)*bb(3)
b2(6) = bb(6)*bb(1) + bb(5)*bb(4) + bb(3)*bb(6)
trbb2 =  bb(1)**2 + bb(2)**2 + bb(3)**2 + (bb(4)**2 + bb(5)**2 + bb(6)**2)*2.0


i1b = bb(1) + bb(2) + bb(3)
i2b = 0.5*(i1b*i1b - trbb2)

trbb3 = i1b * one3
do i = 1,3
  del_i1b(i) = bb(i) - trbb3
  del_i1b(i+3) = bb(i+3)
end do


do i = 1, ntm
  del_i2b(i) = i1b*bb(i) - b2(i)
end do

del_j = 0.0

do i = 1, 3
  del_i2b(i) = del_i2b(i) - two3*i2b
  del_j(i) = 0.5*detf
end do



!     Compute stress tensor

sig = 0.0

do n = 1, n_order
  do np = 0, n
    nq = n - np
    if (np > 0) then
      fac1 = 2.0*c_pq(np,nq) * np * (i1b-3.0)**(np-1) * (i2b-3.0)**nq
      sig = sig + fac1*del_i1b(1:ntm)
    endif
    if (nq > 0) then
      fac1 = 2.0*c_pq(np,nq) * nq * (i1b-3.0)**np * (i2b-3.0)**(nq-1)
      sig = sig + fac1*del_i2b(1:ntm)
    endif
  end do
end do



!     Compute tangent stiffness tensor


! build the push-forwarded dd_i1b

dd_i1b = 0.0

do i = 1, 3
  dd_i1b(i,i) = one3*i1b
  do j = 1, 3
    dd_i1b(i,j) = dd_i1b(i,j) - one3*(del_i1b(j) + bb(i))
  end do
  do j = 4, ntm
    dd_i1b(i,j) = dd_i1b(i,j) - one3*del_i1b(j)
  end do
end do

do i = 4, ntm
  dd_i1b(i,i) = one6*i1b
  do j = 1, 3
    dd_i1b(i,j) = dd_i1b(i,j) - one3*bb(i)
  end do
end do



! build the push-forwarded dd_i2b


dd_i2b = 0.0

do i = 1, 3
  dd_i2b(i,i) = two3*i2b
  do j = 1, 3
    dd_i2b(i,j) = dd_i2b(i,j) - two3*del_i2b(j) - one3*i1b*bb(i) + two3*b2(i) + bb(i)*del_i1b(j)
  end do
  do j = 4, ntm
    dd_i2b(i,j) = dd_i2b(i,j) - two3*del_i2b(j) + bb(i)*del_i1b(j)
  end do
end do

do i = 4, ntm
  dd_i2b(i,i) = one3*i2b
  do j = 1, 3
    dd_i2b(i,j) = dd_i2b(i,j) - one3*i1b*bb(i) + two3*b2(i) + bb(i)*del_i1b(j)
  end do
  do j = 4, ntm
    dd_i2b(i,j) = dd_i2b(i,j) + bb(i)*del_i1b(j)
  end do
end do

call cross_over(bb, -1.0, dd_i2b)


! build tangent stiffness tensor

aa = 0.0

do n = 1, n_order
  do np = 0, n
    nq = n - np
    if (np > 1) then
      fac1 = 4.0*c_pq(np,nq) * np * (np-1)* (i1b-3.0)**(np-2) * (i2b-3.0)**nq
      forall (i=1:ntm, j=1:ntm) aa(i,j) = aa(i,j) + fac1*del_i1b(i)*del_i1b(j)
    endif
    if (np > 0) then
      fac1 = 4.0*c_pq(np,nq) * np * (i1b-3.0)**(np-1) * (i2b-3.0)**nq
      forall (i=1:ntm, j=1:ntm) aa(i,j) = aa(i,j) + fac1*dd_i1b(i,j)
    endif
    if ((np > 0).AND.(nq > 0)) then
      fac1 = 4.0*c_pq(np,nq) * np * nq * (i1b-3.0)**(np-1) * (i2b-3.0)**(nq-1)
      forall (i=1:ntm, j=1:ntm) aa(i,j) = aa(i,j) + fac1*(del_i1b(i)*del_i2b(j)+del_i1b(j)*del_i2b(i))
    endif
    if (nq > 1) then
      fac1 = 4.0*c_pq(np,nq) * nq * (nq-1)* (i1b-3.0)**np * (i2b-3.0)**(nq-2)
      forall (i=1:ntm, j=1:ntm) aa(i,j) = aa(i,j) + fac1*del_i2b(i)*del_i2b(j)
    endif
    if (nq > 0) then
      fac1 = 4.0*c_pq(np,nq) * nq * (i1b-3.0)**np * (i2b-3.0)**(nq-1)
      forall (i=1:ntm, j=1:ntm) aa(i,j) = aa(i,j) + fac1*dd_i2b(i,j)
    endif
  end do
end do


!     Compute viscoelasticity (Deviatoric only)

if (visco_flag .AND. (.NOT.(vol_visco_flag))) then
  nv = (size_d - (pre_props + num_internal_props))/2
  call viscoe_fn(vol_visco_flag,nv,d,ta,dt,dt_shift,f,finv,detf,state(1),state(ntm+1),ntm,sig,aa)
endif



!     Compute pressure term

! Thermal expansion

t_beta = 1.0
J_th = (1.0 + t_beta*t_alpha*(ta - tl0))**3
J_el = detf/J_th


do n = 1, 1
  fac2 = 2.0 * n * J_el/D_p(n)
  fac1 = fac2 * (J_el-1.0)**(2*n-1)   ! pressure (tension-positive term)
  forall (i=1:3) sig(i) = sig(i) + fac1
  fac3 = fac1 + fac2 * (2*n-1) * J_el * (J_el-1.0)**(2*n-2)
  fac4 = 2.0 * fac1
  forall (i=1:3, j=1:3) aa(i,j) = aa(i,j) + fac3
  forall (i=1:3) aa(i,i) = aa(i,i) - fac4
  forall (i=4:ntm) aa(i,i) = aa(i,i) - fac1
end do


  
!     Compute viscoelasticity (Volumetric & Deviatoric)
  
if (vol_visco_flag) then
  nv = (size_d - (pre_props + num_internal_props))/2
  call viscoe_fn(vol_visco_flag,nv,d,ta,dt,dt_shift,f,finv,detf,state(1),state(ntm+1),ntm,sig,aa)
endif

!     Compute spatial Cauchy stress and tangent modulus

sig = detfi*sig
aa = detfi*aa


! Green strain tensor for post-processing: mechanical part only
cc = 0.5*cc
do i = 1,3
	cc(i) = cc(i) - 0.5 - t_alpha*(ta - tl0)
end do ! i

end subroutine general_mnrv








subroutine viscoe_fn(vol_visco_flag,nv,d,ta,dt,dt_shift,fi,finv,xj,sn,hn,ntm,tau,aa)

!  Purpose:  Viscoelasticity model for 3D/2D finite deformation
!           - 2D cases: ntm = 4
!     History data:
!       sn(ntm)
!       hn(ntm,nv)

!		vol_visco_flag = .TRUE.: deviatoric & volumetric viscoelasticity

implicit  none

logical, intent(in) :: vol_visco_flag
integer, intent(in) :: nv, ntm
real, intent(in) :: ta, dt, d(:), fi(3,3),finv(3,3), xj
real, intent(out) :: dt_shift, aa(ntm,ntm)
real, intent(inout) :: tau(:), sn(ntm), hn(ntm,*)

integer :: i, j, n, nd
real :: gg, xj_23, xji23, trsmh, gammai, gexpi
real :: tau0(ntm), PI2(ntm), smh(ntm), bigh(ntm), snold(ntm), ht(ntm)
real :: t_ref_wlf, c1_wlf, c2_wlf
real, parameter ::  one3 = 1.0/3.0, two3 = 2.0/3.0
real, parameter :: one(6) = [1.0,1.0,1.0,0.0,0.0,0.0]



!     Save elastic stress for updates
do i = 1, ntm
  tau0(i) = tau(i)
end do

!     Compute viscoelastic parameters

if (nv > 0) then

  xj_23 = xj**(two3)   ! Jn**(1-1/3)
  xji23 = 1.0/xj_23

  call push_pull(finv, tau0, PI2, 1, ntm)

!       Update history parameters for 2nd Piola-Kirchhoff stress

  do i = 1, ntm
    snold(i) = sn(i)
    sn(i)    = PI2(i)*xj_23
    bigh(i)  = 0.0
  end do

!       Sum over series of terms in viscoelasticity

  dt_shift = dt   ! copy for recovering when return
  
  t_ref_wlf = d(10)    ! T_ref for WLF model

	if (t_ref_wlf > -273.15) then
    c1_wlf = d(11)
    c2_wlf = d(12)
		call wlf(c1_wlf,c2_wlf,t_ref_wlf,ta,dt_shift)
	endif
	if (dt_shift <= 0.0) stop 'viscfd: dt and/or dt_shift is zero!'

  gg = 1.0
  nd = pre_props + num_internal_props  ! module global variables
  do n = 1, nv
    gammai = d(2*n+nd-1)
    gexpi  = exp(-0.5*dt_shift/d(2*n+nd))

    do i = 1, ntm
      ht(i)   = gexpi*(gexpi*hn(i,n) - snold(i))
      hn(i,n) = ht(i) + gexpi*sn(i)
      bigh(i) = bigh(i) + gammai*ht(i)
      
    end do ! i
    gg = gg - gammai*(1.0 - gexpi)
  end do ! n


!       Push forward to compute stress deviator

  call push_pull(fi, bigh, smh, 1, ntm)

  do i = 1, ntm
    smh(i) = smh(i)*xji23
  end do
  
  if (ntm == 1) then
    trsmh = 0.0
  else

    if (vol_visco_flag) then
      trsmh = 0.0
    else
      trsmh = one3*(smh(1) + smh(2) + smh(3))
      do i = 1, 3
        smh(i) = smh(i) - trsmh
      end do
    endif
  endif

!       Update Kirchhoff stress
  do i = 1, ntm
    tau(i) = gg*tau0(i) + smh(i)
  end do

!       Update tangent
  if (vol_visco_flag) then
    do i = 1, ntm
      do j = 1, ntm
        aa(i,j) = gg*aa(i,j) - two3*(smh(i)*one(j) + one(i)*smh(j))
      end do
    end do

  else

    do i = 1, ntm
      do j = 1, ntm
        aa(i,j) = gg*aa(i,j) - two3*((smh(i)+trsmh*one(i))*one(j) + one(i)*smh(j))
      end do
      aa(i,i) = aa(i,i) + trsmh*(1.0 + one(i))
    end do
  endif

endif

end subroutine viscoe_fn








subroutine push_pull(f, tensor1, tensor2, rank, ndim)

implicit none

integer, intent(in) :: rank, ndim
real, intent(in) :: f(3,3), tensor1(ndim,ndim)
real, intent(out) :: tensor2(ndim,ndim)

integer, parameter :: ix(2,9) = RESHAPE((/1,1,2,2,3,3,1,2,2,3,3,1,2,1,3,2,1,3/), (/2,9/))
integer :: i, j
real :: ff(ndim,9), ff_trans(9,ndim), matrix(ndim,9), tensor0(9,9)


if (ndim < 4) STOP 'push_pull: ndim is smaller than 4!'

ff = 0.0
tensor0 = 0.0

do i = 1, ndim
  do j = 1, 9
    ff(i,j) = f(ix(1,i),ix(1,j))*f(ix(2,i),ix(2,j))
  end do
end do


if (rank == 1) then
! expansion of tensor in vector form from symmetric size 6 to full size 9
  tensor0(1:ndim,1) = tensor1(1:ndim,1)
  do i = 4, ndim
    tensor0(i+3,1) = tensor1(i,1)
  end do
  tensor2(1:ndim,1) = MATMUL(ff(1:ndim,1:9),tensor0(1:9,1))

elseif (rank == 2) then
! expansion of tensor in matrix form from symmetric size 6X6 to full size 9X9
  tensor0(1:ndim,1:ndim) = tensor1(1:ndim,1:ndim)
  do i = 4, ndim
    tensor0(i+3,1:ndim) = tensor1(i,1:ndim)
  end do
  do i = 4, ndim
    tensor0(1:ndim,i+3) = tensor1(1:ndim,i)
  end do

  matrix = MATMUL(ff,tensor0)

  ff_trans = TRANSPOSE(ff)
  tensor2 = MATMUL(matrix,ff_trans)
else
  STOP 'push_pull: rank is wrong!'
endif

end subroutine push_pull





subroutine tensorizer(vector, length, tensor, ndim)

! matrix form recovery for symmetric tensor in vector form

implicit none

integer, intent(in) :: length, ndim
real, intent(in) :: vector(length)
real, intent(out) :: tensor(ndim,ndim)

integer :: i


if ((ndim < 2).OR.(length < 4)) STOP 'tensorizer: minimum vector size is 4 & minimum tensor matrix is 2X2!'

do i = 1, ndim
  tensor(i,i) = vector(i)
end do

tensor(1,2) = vector(4)
tensor(2,1) = tensor(1,2)

if (length == 6) then
  tensor(1,3) = vector(6)
  tensor(2,3) = vector(5)
  tensor(3,1) = tensor(1,3)
  tensor(3,2) = tensor(2,3)
endif

end subroutine tensorizer






subroutine cross_over(v, factor, tensor)

! tensor(i,j) = tensor(i,j) + 0.5*(b(I,k)*b(J,L) + b(I,L)*b(J,K))
!   i,j = 1..6
!   I,J,K,L = 1..9
!   b(1:9,1:9) <- v(1:6) : symmetrically tensorized

implicit none

real, intent(in) :: v(6), factor
real, intent(inout) :: tensor(6,6)

integer :: i, j, k, l, ii, jj
integer, parameter :: ix(3,3) = RESHAPE((/1,4,6,4,2,5,6,5,3/), (/3,3/))


do l = 1,3
  do k = l,3
    jj = ix(k,l)
    do j = 1,3
      do i = j,3
        ii = ix(i,j)
        tensor(ii,jj) = tensor(ii,jj) + 0.5*factor*(v(ix(i,k))*v(ix(j,l)) &
                        + v(ix(i,l))*v(ix(j,k)))
      end do ! i
    end do ! j
  end do ! k
end do ! l

end subroutine cross_over








subroutine wlf(c1,c2,t_ref,ta,dt)

! WLF equation shift for temperature consideration in viscoelastic model
!   - output: modified dt

real, intent(in) :: c1 , c2, t_ref, ta
real, intent(inout) :: dt
real :: log_at
real, parameter :: Numerical_zero = 1.0e-10


if (ABS(c2 + ta - t_ref) > Numerical_Zero) then
  log_at = -c1*(ta - t_ref)/(c2 + ta - t_ref)
  dt = dt/(10.0**log_at)
endif

end subroutine wlf


! module contains end -----------------------------------------------------	  
  	  
end module mnrv_2_material