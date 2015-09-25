subroutine nl_2d_fn_mx(action,iter,pro_type,reg_num,ele_num,ul,xl,tl,tl0,s,p,str_out,ndf,ndm,nst,nel,nquad,n_out)

! Coded by Jeeho Lee (April 2010)
!		- Dynamic part (add Mass & Ground acceleration): Aug 2011
!		- Bad deformation exit signal: March 2014

!     Plane/axisymmetric large deformation element routine - Mixed formulation
!     4-node element formulation
!     Present version: Statics problems only

!     pro_type: plane strain = 1, axisymmetric = 2	  

use element_specification

! Global variables from module 'elements' & 'element_specification':
!
!  - Geometry: thickness, body(3), gr_accl(3), density, consistent_ratio
!  - Time related: dt, ctan(3)


implicit none

integer, intent(inout) :: action
integer, intent(in) :: pro_type, reg_num, ele_num, ndf, ndm, nst, nel, nquad, iter, n_out
real, intent(in) :: xl(ndm,nel), ul(ndf,5*nel), tl(nel), tl0
real, intent(out) :: s(nst,nst), p(nst), str_out(nquad,n_out)

integer, parameter :: ndim = 4, max_lint = 4
integer ::  mat_num, material_model, point_address
      
integer :: i, i1, j, jj, j1, k, k1, l, lint, nea
integer :: npm
real :: xsj
real :: dsigtr, mpress, dpress, dtheta, fac, bdb, bd3
real :: sg(3,max_lint), r(ndf,nel), xx(2,max_lint), shp(3,16,max_lint), dvol0(max_lint), dvol(max_lint)
real :: bbd(3,7), bbar(2,16,max_lint), tang(ndim,ndim), ad(6,6,max_lint), dd(7,7), sigm(4), sigl(6,max_lint)
real :: detf(2,max_lint), df(9,max_lint), fi(9,2,max_lint), finv(9,max_lint)
real :: phi(6,max_lint), theta(2,max_lint), hh(6,6), vl(3)
real :: press(max_lint), pbar(max_lint)
real :: irad(max_lint), ta(max_lint)
real :: xr(2,max_lint), ur(2,max_lint), xu(2,16), el(4,7), ru(3,16), shpr(16,max_lint), shpbar(2,16,max_lint)

real :: b1, b2, thick
real :: str_data(n_out)

real :: mass(nel,nel), al(2), accl(ndf,nel), infv(2), lfac, cfac, rr, m(nst,nst)
real :: aj0, ct1, ct3

logical :: flgrm = .false., dyna = .false.
logical :: ther, quad, geo_flag



if ((nel > 4) .or. (nel < 4)) stop 'nl_2d_fn_mx: Too many nodes for MIXED mode - valid only for 4 node element!'

!if (nel == 4 .or. nel == 8 .or. nel == 9) then            
!  quad = .true.
!else
!  quad = .false.
!endif
quad = .true.

if (action == 2) then
	dyna = .true.
else
	dyna = .false.
endif

ther = .false.
geo_flag = .true.

!       Set element quadrature order
npm = 1
l   = 2
nea = 4*nel
         
s = 0.0
p = 0.0
r = 0.0

!   read data from module 'element_specification'
thick = thickness
b1 = body(1)
b2 = body(2)
rr = density

accl = 0.0
infv = 0.0	  

if (dyna) then  ! dynamic case	

	cfac = consistent_ratio
	accl = ul(1:ndf,nea+1:nea+nel)

  lfac = 1.d0 - cfac

	mass = 0.0

	ct1 = thick*ctan(1)
	ct3 = thick*ctan(3)
	
	
	! Ground motion part (default: off)
	flgrm = .true.
  if (flgrm) then
		infv(1) = gr_accl(1)
		infv(2) = gr_accl(2)
	else
		infv = 0.0
	endif

else
	ct1 = thick
	ct3 = 0.0
  
endif



!     Compute current geometry (useless now)

do j = 1, nel
  do i = 1, 2
    xu(i,j) = xl(i,j) + ul(i,j)
  end do ! i
end do ! j


call element_material(reg_num, ele_num, mat_num, material_model, .TRUE.)


if (quad) then
  l = 2
  call int2d(l,lint,sg,nel)
else    ! triangular elements
  if (nel == 3) then
    if (pro_type == 1) then
        l =  1
    else
        l = -3
    endif
    call tint2d(l,lint,el)
  else if (nel == 6 .OR. nel == 7 ) then
    l = 7
    call tint2d(l,lint,el)
  endif
endif

if (lint > max_lint) then
  write(*,*) 'nl_2d_fn_mx: too many quadrature points!: lint =', lint
	STOP
endif

phi = 0.0

do l = 1, lint

!         Shape functions and derivatives

  call shp2d(sg(1,l),xl,shp(1,1,l),xsj,ndm,nel,.FALSE.)
  dvol0(l) = xsj*sg(3,l)

!         Mixed volume effect and temperature projection

  ta(l) = 0.0

!         Compute coordinates at gauss points

  xr(:,l) = 0.0
  ur(:,l) = 0.0
  do i = 1, nel
    xr(1,l) = xr(1,l) + xl(1,i)*shp(3,i,l)
    xr(2,l) = xr(2,l) + xl(2,i)*shp(3,i,l)
    ur(1,l) = ur(1,l) + ul(1,i)*shp(3,i,l)
    ur(2,l) = ur(2,l) + ul(2,i)*shp(3,i,l)
            
    ta(l)   = ta(l) + shp(3,i,l)*tl(i)
  end do


!         Compute volumetric strain from displacements

  if (pro_type == 2) then   ! Axisymmetric case
    dvol0(l) = dvol0(l)*xr(1,l)
    do i = 1,nel
      shpr(i,l)  = shp(3,i,l)/(xr(1,l) + ur(1,l))
    end do
  else                     ! Plane case
    do i = 1,nel
      shpr(i,l)  = 0.0
    end do
  endif

 !         Set the pressure functions
  phi(1,l) = 1.0
  
end do ! l

vl(3) = 0.0


!       Mixed model for volumetric and temperature response

!         Compute f, finv, df, detf and shp at conf t-n+1

call kine2d(pro_type,shp,xl,ul,fi,finv,df,detf,ndm,ndf,nel,lint)

!		- Bad deformation exit signal
if (MINVAL(detf(1,1:lint)) < Numerical_Zero) then
  action = -1
  RETURN  !-------------------------------------------------------------------->>>>>
endif

!         Compute volume at current state

do l = 1, lint
  dvol(l) = dvol0(l)*detf(1,l)
end do

!         Mixed model for volumetric response

call bbar2m(phi,shp,shpr,dvol,detf,lint,nel,npm,hh,theta,shpbar)



!  write(*,*) 'before: ', l,fi(1,1,l),finv(1,l),theta(1,l),sigl(1,l),ad(1,1,l)

!         Compute mixed model deformation gradient

call fbar(fi,detf,theta,lint)

!         Compute Cauchy stresses and spatial tangent tensor at t-n+1

do l = 1, lint
!           Compute coordinates in reference and current configuration  

!     Get the address in history db using the subrountine in 'physical_domain'
  call find_point_number(reg_num,ele_num,l,point_address)

  call str_tang_fn(action,iter,reg_num,ele_num,point_address,mat_num,material_model,ndm,ndf,nel &
                     ,ta(l),tl0,dt,xl,ul,fi(1,1,l),finv(1,l),theta(1,l),sigl(1,l),tang,ndim,str_data,n_out)
  ad(:,:,l) = 0.0
  ad(1:ndim,1:ndim,l)  = tang
  str_out(l,:) = str_data
end do ! l

press(1) = 0.0
do l = 1,lint
!               Modify volume element and integrate pressure over reference volume
  press(1) = press(1) + dvol0(l)*(sigl(1,l) + sigl(2,l) + sigl(3,l))/3.0
  dvol(l)  = dvol0(l)*theta(1,l)
end do ! l

! Divide pressure by reference volume


press(1) = press(1) * hh(1,1)
do l = 2, lint
  press(l) = press(1)
end do ! l


! Compute mixed pressure


do l = 1, lint

! Compute mixed stress and multiply by volume element

  dsigtr  = press(l)*detf(1,l)/theta(1,l) - (sigl(1,l)+sigl(2,l)+sigl(3,l))/3.0
  sigm(1) =  sigl(1,l) + dsigtr
  sigm(2) =  sigl(2,l) + dsigtr
  sigm(3) =  sigl(3,l) + dsigtr
  sigm(4) =  sigl(4,l)

  do i = 1, ndim
    sigm(i) = sigm(i)*dvol(l)
  end do
              
	if (dyna) then			
		al = 0.0
			        
		do j = 1,nel
			do k = 1,2
				al(k) = al(k) + shp(3,j,l)*(accl(k,j) + infv(k))
			enddo		
		enddo
		al = cfac*rr*al		
		aj0 = lfac*rr				
	else
		aj0 = 0.0
		al = 0.0
	endif
							
! Compute residual

  do i = 1, nel
    ru(1,i) = shp(1,i,l)*sigm(1) + shp(2,i,l)*sigm(4)
    ru(2,i) = shp(1,i,l)*sigm(4) + shp(2,i,l)*sigm(2)

    r(1,i)  = r(1,i) - ru(1,i) - shpr(i,l)*sigm(3) + (b1 -al(1) - aj0*(infv(1)+accl(1,i)))*dvol0(l)*shp(3,i,l)
    r(2,i)  = r(2,i) - ru(2,i) + (b2 -al(2) - aj0*(infv(2)+accl(2,i)))*dvol0(l)*shp(3,i,l)
  end do

! Multiply tangent moduli by volume element
  
! Part 1: Geometric tangent matrix
  if (geo_flag) then
  
    i1 = 0
    do i = 1, nel
      bd3 = shpr(i,l)*sigm(3)
      j1 = 0
      do j = 1,nel
        bdb = ru(1,j)*shp(1,i,l) + ru(2,j)*shp(2,i,l) ! i,j are symmetric
        s(i1+1,j1+1) = s(i1+1,j1+1) + bdb + bd3*shpr(j,l)
        s(i1+2,j1+2) = s(i1+2,j1+2) + bdb

        j1 = j1 + ndf
      end do ! j

      i1 = i1 + ndf
    end do ! i
  endif
  
! Part 2: Material tangent matrix

!  Modify tangent moduli for stress factors

  mpress = press(l)*detf(1,l)/theta(1,l)
  dpress = (sigl(1,l) + sigl(2,l) + sigl(3,l))/3.0
  

  call dmatmx(ad(1,1,l),dd)

  call dmatdx(dd,sigl(1,l),dpress,mpress)

!  Multiply tangent moduli by volume element

  dd = dvol(l)*dd


 !  Compute row terms

  i1    = 0
  do i = 1, nel
!                 Compute bmat-t * dd * dvol

    do jj = 1, 7

      bbd(1,jj) =   shp(1,i,l)*dd(1,jj)  &
                +   shpr( i,l)*dd(3,jj)  &
                +   shp(2,i,l)*dd(4,jj)  &
                + shpbar(1,i,l)*dd(7,jj)

      bbd(2,jj) =    shp(2,i,l)*dd(2,jj)  &
                +    shp(1,i,l)*dd(4,jj)  &
                + shpbar(2,i,l)*dd(7,jj)
    end do ! jj

    j1 = 0
    do j = 1, nel  ! can use 'i' instead of 'nel', if s(:,:) is definitely symmetric

      do jj = 1, 2
        s(i1+jj,j1+1) = s(i1+jj,j1+1)  &
                      + bbd(jj,1)*shp(1,j,l)  &
                      + bbd(jj,4)*shp(2,j,l)  &
                      + bbd(jj,3)*shpr(j,l)  &
                      + bbd(jj,7)*shpbar(1,j,l)

        s(i1+jj,j1+2) = s(i1+jj,j1+2)  &
                      + bbd(jj,2)*shp(2,j,l)  &
                      + bbd(jj,4)*shp(1,j,l)  &
                      + bbd(jj,7)*shpbar(2,j,l)
      end do ! jj
      j1 = j1 + ndf
    end do ! j
    i1 = i1 + ndf
  end do ! i

	if (dyna) then

		do i = 1, nel
			aj0 = shp(3,i,l)*dvol0(l)*rr
		
!		Compute lumped mass matrix
			mass(i,i) = mass(i,i) + lfac*aj0
			                          
!		Compute consistent mass matrix
			do j = 1, nel
				mass(i,j) = mass(i,j) + cfac*aj0*shp(3,j,l)
			end do ! j
		end do ! i		
    
	endif  ! dyna
	
end do ! l


if (dyna) then
!  Expand 'mass' matrix (scalar mass: no directional mass) to the general m(nst,nst) matrix

	m = 0.0
	j1 = 1
	do j = 1, ndf*nel, ndf
      k1 = 1
      do k = 1,ndf*nel,ndf
        m(j,k) = mass(j1,k1)
        m(j+1,k+1) = mass(j1,k1)
        k1 = k1 + 1
      end do
      j1 = j1 + 1
	end do
	
!  Add mass matrix(mass) & damping matrix(damp) to stiffness matrix(s)
      
	s = ct1*s + ct3*m

else
	s = ct1*s   ! thickness
endif  ! dyna


j = 1
do i = 1, nel
  p(j) = r(1,i)*thick
  p(j+1) = r(2,i)*thick
  j = j + 2
enddo


! write(*,*) 'p: ', p
! write(*,*) 's: ', s

end subroutine nl_2d_fn_mx
