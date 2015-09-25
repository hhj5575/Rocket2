module viscoelastic_material

!     Viscoelastic Model
!
!     Last Modification: Feb, 2010
!                        May, 2013
!                        March, 2014
	
use materials

implicit none

integer, parameter :: max_size_d = 50, max_size_state = 90  ! for maximum ndim = 6 & max 15-set visco data
integer, parameter :: max_mu_data_set = 15, pre_props = 10, num_internal_props = 2

contains
!==============================================================================
            
subroutine viscoelastic(action,vol_flag,there,mat_num,ta,tl0,delt,eps,str,tang,secd,ndim)

implicit none
	  
integer, intent(in) :: there, mat_num, action, ndim, vol_flag
real, intent(in) :: delt, ta, tl0
real, intent(out) :: str(ndim), tang(ndim,ndim), secd(ndim,ndim)
real, intent(inout) :: eps(ndim)
	  
integer :: npara, size_d, size_state, nv, n_pointer
real :: d_array(max_size_d), state(max_size_state)


select case (action) 
case(0)
  
	if (.NOT. material(mat_num)%flag) then
  
    npara = material(mat_num)%size_para
		call make_d_viscoe(npara,material(mat_num)%para,d_array,size_d)
		call set_d_array(1,mat_num,size_d,d_array)
    
    if (size_d > (pre_props + num_internal_props)) then
      nv = (size_d - (pre_props + num_internal_props))/2
      n_pointer = ndim + 1 + ndim*nv
      size_state = n_pointer + 1
    else
      size_state = 1
    endif
    
    call number_material_state_var(mat_num,size_state) 
    write(*,*) 'maximum size_state & size_state: ', max_size_state, size_state
    if (size_state > max_size_state) STOP 'viscoelastic: size_state exceeds maximum setting!'
	else
    write(*,*) 'Warning!!! viscoelstic: material initialize flag is not FALSE!'
	endif 
		
case default
  size_d = material(mat_num)%size_d
	call set_d_array(2,mat_num,size_d,d_array)
		
! The subrountine in 'materials' to download history data
  size_state = material(mat_num)%num_row_state
	call material_history(1,mat_num,there,state(1:size_state),size_state)	  
	  
! Call viscoelastic material model
  nv = (size_d - (pre_props + num_internal_props))/2

	call viscoe_linear(vol_flag,nv,size_d,d_array,ta,tl0,delt,eps,state(1),state(ndim+1),state(size_state-1),state(size_state),ndim,str,tang,secd)

! The subrountine in 'materials' to upload history data
	call material_history(2,mat_num,there,state(1:size_state),size_state)
	  
end select
	  
end subroutine viscoelastic




	  	  
subroutine make_d_viscoe(nprops,props,d,size_d)

implicit none

integer, intent(in) :: nprops
real, intent(in) :: props(nprops)
real, intent(out) :: d(:)
integer, intent(out) :: size_d

integer :: i, j, num_visco_para  
		

if ((pre_props + 2) > nprops) stop 'Viscoelastic: too few data (minimum 12)!'

d = 0.0
	  
d(1) = props(1)    ! Young's modulus
d(2) = props(2)    ! Poisson's ratio

d(3) = d(1)/(2.d0*(1.d0 + d(2)))        ! G modulus
d(4) = d(1)/(3.d0*(1.d0 - 2.d0*d(2)))   ! K modulus

d(5) = props(3)    ! thermal expansion coefficient

! d(6) = props(4)    ! any physical value is not assigned
! d(7) = props(5)    ! any physical value is not assigned
! d(8) = props(6)    ! any physical value is not assigned
! d(9) = props(7)    ! any physical value is not assigned

d(10) = props(8)   ! T_ref for WLF model
d(11) = props(9)   ! c1 in WLF
d(12) = props(10)  ! c2 in WLF

size_d = pre_props + num_internal_props ! num_internal_props: G & K


if (nprops > pre_props) then    ! viscocity case
  j = nprops - pre_props
  if ((j/2 > max_mu_data_set) .OR. (mod(j,2) == 1)) STOP 'Viscoelastic: visco pair error!'
  
  num_visco_para = 0
  j = size_d - pre_props
  do i = (pre_props+1), nprops
    d(i+j) = props(i)
    num_visco_para = num_visco_para + 1
  end do
  size_d = size_d + num_visco_para
endif

end subroutine make_d_viscoe

            
            
            
            
subroutine viscoe_linear(vol_flag,nv,size_d,d,ta,tl0,dt,eps,en,qi,thetan,qtheta,ntm,sig,dd,dr)

!      Linear viscoelastic material model

!      Inputs:
!         d(*)    - Material parameters
!         ta      - Temperature
!         tl0     - Stress-free (reference) temperature
!         eps(ntm)  - Strains at t_n+1
!         en(ntm)   - Strains at t_n
!         qi(ntm,nv)   - Viscoelastic strain
!         ntm     - Number of terms

!      Outputs:
!         sig(ntm)  - Stresses
!         eps(ntm)  - Mechanical Strains at t_n+1 for post-processing
!         dd(6,6) - Viscoelastic moduli
!         dr(6,6) - Rayleigh moduli (Elastic stiffness)

implicit  none

integer, intent(in) :: nv, ntm, vol_flag, size_d
real, intent(in) :: d(:), ta, tl0, dt
real, intent(inout) :: eps(ntm), en(ntm), qi(ntm,nv), thetan, qtheta(nv)
real, intent(out) :: sig(ntm), dd(ntm,ntm), dr(ntm,ntm)

integer :: i, j, n, nd
logical :: vol_visco_flag
real :: gfac, exp_n, mu_0, mu_n, dq_n, dtau, theta
real :: alpha, e, nu, G, Gg, K, Kg, G2, Gg2
real :: hvisc, hh, press, ee(ntm), epst
real :: dt_shift, t_ref_wlf, c1_wlf, c2_wlf
real, parameter :: one2 = 1.0/2.0, two3 = 2.0/3.0

! Set elastic parameters for G (mu) and lambda

e = d(1)
nu = d(2)
G = d(3)
K = d(4)
alpha = d(5)   ! thermal expansion coefficient (line)

dd = 0.0
dr = 0.0

if (vol_flag == 2) then
  vol_visco_flag = .TRUE.
else
  vol_visco_flag = .FALSE.
endif

! write(*,'(/A,6E11.3)') 'eps = ', eps
! write(*,'(A,6E11.3)')  'en  = ', en
! write(*,'(A)')  'qi (in)= '
! call print_internal_var(qi,ntm, nv)

! Compute volumetric strain and deviatoric components

theta = (eps(1) + eps(2) + eps(3))/3.0
do i = 1, 3
  ee(i) = eps(i) - theta
end do ! i
do i = 4, ntm
  ee(i) = one2*eps(i)
end do ! i

! Set properties for integrating the q_i terms

hh = 0.0

sig = 0.0
mu_0 = 0.0
gfac = 0.0

dt_shift = dt   ! copy for recovering when return
  
t_ref_wlf = d(10)    ! T_ref for WLF model

if (t_ref_wlf > -273.15) then
  c1_wlf = d(11)
  c2_wlf = d(12)
  call wlf(c1_wlf,c2_wlf,t_ref_wlf,ta,dt_shift)
endif
if (dt_shift <= 0.0) STOP 'viscfd: dt and/or dt_shift is zero!'

nd = pre_props + num_internal_props   ! module global variables
do n = 1, nv
  mu_n = d(2*n+nd-1)
  dtau = dt_shift/d(2*n+nd)
  exp_n = exp(-dtau)
  call vis_func(dtau,exp_n,hvisc)
  dq_n = hvisc         ! delta_h
  gfac = gfac + mu_n*dq_n
  mu_0 = mu_0 + mu_n

! Update history and compute viscoelastic deviatoric stress
  do i = 1, ntm
    qi(i,n) = exp_n*qi(i,n) + dq_n*(ee(i) - en(i))
    sig(i)  = sig(i) + mu_n*qi(i,n)
  end do
  qtheta(n) = exp_n*qtheta(n) + dq_n*(theta - thetan)
  hh = hh + mu_n*qtheta(n)
end do ! n


! Finish updates and save the strains
mu_0 = 1.d0 - mu_0
gfac = gfac + mu_0

G2 = 2.0*G
do i = 1, ntm
  sig(i) = G2*(mu_0*ee(i) + sig(i))
  en(i)  = ee(i)
end do ! i

if (.NOT. vol_visco_flag) then
  mu_0 = 1.0
  hh = 0.0
  Kg = K
else
  Kg = K*gfac
endif


epst = alpha*(ta - tl0)

press = 3.0*K*(mu_0*(theta - epst) + hh)
thetan = theta
      
! Add elastic bulk term

do i = 1, 3
  sig(i) = sig(i) + press
	eps(i) = eps(i) - epst   ! mechanical strain for post-processing
end do ! i

! Set tangent parameters
Gg = G*gfac
Kg = Kg - two3*Gg
K  = K - two3*G
Gg2 = G2*gfac

do j = 1, 3
  do i = 1, 3
    dd(i,j) = Kg
    dr(i,j) = K
  end do ! i
  dd(j,j) = dd(j,j) + Gg2
  dr(j,j) = dr(j,j) + G2
end do ! j

do i = 4, ntm
  dd(i,i) = Gg
  dr(i,i) = G
end do ! i

! write(*,'(A)')  'qi (out)= '
! call print_internal_var(qi,ntm, nv)

end subroutine viscoe_linear





subroutine vis_func(x,expx,hvisc)

!     Purpose: Compute integration factor for viscoelastic material
!              H_visc(x) = [ 1 - exp(-x)]/x

!     Inputs:
!       x        - Current relaxation time parameter
!       expx     - Exponential of x

!     Outputs:
!       hvisc    - Viscoelastic function

implicit none

real, intent(in) :: x, expx
real, intent(out) :: hvisc


if (x < 1.0d-04) then
    hvisc = 1.d0 - 0.5d0*x*(1.d0 - x/3.d0*(1.d0 - 0.25d0*x*(1.d0 - 0.2d0*x)))
else
    hvisc = (1.d0 - expx)/x
endif

end subroutine vis_func






subroutine wlf(c1,c2,t_ref,ta,dt)

! WLF equation shift for temperature consideration in viscoelastic model
!   - output: modified dt

implicit none

real, intent(in) :: c1 , c2, t_ref, ta
real, intent(inout) :: dt
real :: log_at
real, parameter :: Numerical_zero = 1.0e-10


if (ABS(c2 + ta - t_ref) > Numerical_Zero) then
  log_at = -c1*(ta - t_ref)/(c2 + ta - t_ref)
  dt = dt/(10.0**log_at)
endif

end subroutine wlf




subroutine print_internal_var(qi,ntm, nv)

implicit none

integer, intent(in) :: ntm, nv
real, intent(in) :: qi(ntm,nv)
integer :: i, n


do i = 1, ntm
  write(*,'(10E11.3)')  (qi(i,n), n=1,nv)
end do
write(*,'(/)')

end subroutine print_internal_var

! module contains end -----------------------------------------------------	  
  	  
end module viscoelastic_material