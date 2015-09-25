subroutine nl_2d_fn(action,iter,pro_type,reg_num,ele_num,ul,xl,tl,tl0,s,p,str_out,ndf,ndm,nst,nel,nquad,n_out)

! Coded by Jeeho Lee (May 2010)
!		- Dynamic part (add Mass & Ground acceleration): Aug 2011
!		- Bad deformation exit signal: March 2014

!     Plane/axisymmetric large deformation element routine - Displacement formulation
!     4-node element formulation
!     Present version: Statics problems only

!     pro_type: plane strain = 1, axisymmetric = 2	  

use element_specification

! Global variables from module 'elements' & 'element_specification':
!  - Parameters: Numerical_Zero
!  - Geometry: thickness, body(3), gr_accl(3), density, consistent_ratio
!  - Time related: dt, ctan(3)


implicit  none

integer, intent(inout) :: action
integer, intent(in) :: pro_type, reg_num, ele_num, ndf, ndm, nst, nel, nquad, iter, n_out
real, intent(in) :: xl(ndm,nel), ul(ndf,5*nel), tl(nel), tl0
real, intent(out) :: s(nst,nst), p(nst), str_out(nquad,n_out)

integer, parameter :: ndim = 4, max_lint = 9, max_nel = 9
integer :: mat_num, material_model, point_address
      
integer :: i, i1, j, jj, j1, k, k1, l, lint, nea
real :: xsj
real :: bdb,bd3
real :: xx1,xx2,xx3
real :: r(ndf,nel), r1(3,16), xcur(3)
real :: sg(3,max_lint), bbd(6,3), bb(6), shp(3,16,max_lint), shpr(16), dvol0(max_lint), dvol(max_lint)
real :: detf(2,max_lint), df(3,3,max_lint), fi(9,2,max_lint), finv(3,3,max_lint)

real :: dd(ndim,ndim), sigv(ndim)
real :: xr(2,max_lint), ur(2,max_lint), el(4,7), ta(max_lint)

real :: b1, b2, thick
real :: str_data(n_out)

real :: mass(nel,nel), al(2), accl(ndf,nel), infv(2), lfac, cfac, rr, m(nst,nst)
real :: aj0, ct1, ct3, stress_mean(4)

logical :: flgrm = .false., dyna = .false.
logical :: ther, quad, geo_flag



if (nel > max_nel) stop 'nl_2d_fn: # of nodes per element exceeds maximum!'

if (nel == 4 .or. nel == 8 .or. nel == 9) then            
  quad = .true.
else
  quad = .false.
endif

if (action == 2) then
	dyna = .true.
else
	dyna = .false.
endif

ther = .false.
geo_flag = .true.


! Set element quadrature order
l   = 2
nea = 4*nel

s = 0.0
p = 0.0
r = 0.0
shpr = 0.0

! Read data from module 'element_specification'
thick = thickness
b1 = body(1)
b2 = body(2)
rr = density

accl = 0.0
infv = 0.0	  

if (dyna) then

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


call element_material(reg_num, ele_num, mat_num, material_model, .TRUE.)


!  Compute tangent stiffness and residual force vector

!     Set quadrature points and weights

if (quad) then
	if (nquad == 4) then
		l = 2
	elseif (nquad == 1) then
		l = 1
	else
		l = 3
	endif	
 
  call int2d(l,lint,sg,nel)
else    ! triangular elements
  if (nel == 3) then
    if (pro_type == 1) then
        l =  1
    else
        l = -3
    endif
    call tint2d(l,lint,el)
  else if (nel == 6 .or. nel == 7 ) then
    l = 7
    call tint2d(l,lint,el)
  endif
endif

if (lint > max_lint) then
  write(*,*) 'nl_2d_fn: Too many quadrature points!: lint =', lint
	stop
endif

!     COMPUTE TANGENT STIFFNESS AND RESIDUAL FORCE VECTOR

!     Compute shape functions and derivatives in reference configuration

do l = 1, lint

  if (quad) then
    call shp2d(sg(1,l),xl,shp(1,1,l),xsj,ndm,nel,.false.)
    dvol0(l) = xsj*sg(3,l)  ! sg(3,l): quadrature weights
  else
    call trishp(el(1,l),xl,ndm,nel-4,xsj,shp(1,1,l))
    dvol0(l) = xsj*el(4,l)  ! el(4,l): quadrature weights
  endif

!       Compute coordinates at gauss points

  ta(l) = 0.0

  xr(:,l) = 0.0
  ur(:,l) = 0.0
  do i = 1, nel
    xr(1,l) = xr(1,l) + xl(1,i)*shp(3,i,l)
    xr(2,l) = xr(2,l) + xl(2,i)*shp(3,i,l)
    ur(1,l) = ur(1,l) + ul(1,i)*shp(3,i,l)
    ur(2,l) = ur(2,l) + ul(2,i)*shp(3,i,l)
    
    ta(l)   = ta(l) + shp(3,i,l)*tl(i)
  end do ! i
  

!       Axisymmetric volume

  if (pro_type == 2) then     ! Axisymmetric case
    dvol0(l) = dvol0(l)*xr(1,l)
  endif
end do ! l
      
xcur(3) = 0.0d0

!     Compute deformation gradient and determinant; transform shape
!     functions to current configuration.

call kine2d(pro_type,shp,xl,ul,fi,finv,df,detf,ndm,ndf,nel,lint)

!		- Bad deformation exit signal
if (MINVAL(detf(1,1:lint)) < Numerical_Zero) then
  action = -1
  RETURN  !-------------------------------------------------------------------->>>>>
endif


!     Set the loop limits and consistent/lumped mass factor



!     Compute proper rank mass effects for 3-node triangle

!      if(cfac.gt.0.0d0 .and. nel.eq.3) then
!        call masst3(stype,cfac,d(4),xl,ul(1,1,5),r,s)
!        cfac = 0.0d0


!     LOOP OVER GAUSS POINTS

do l = 1, lint

!       Set reference coordinates
  xcur(1) = xr(1,l) + ur(1,l)
  xcur(2) = xr(2,l) + ur(2,l)

!       Check for axisymmetry
  if (pro_type == 2) then
    do i = 1,nel
      shpr(i) = shp(3,i,l)/xcur(1)
    end do ! i
  else
    do i = 1,nel
      shpr(i) = 0.0
    end do ! i
  end if

!       Compute Cauchy stresses and spatial tangent tensor

!     Get the address in history db using the subrountine in 'physical_domain'
  call find_point_number(reg_num,ele_num,l,point_address)

  call str_tang_fn(action,iter,reg_num,ele_num,point_address,mat_num,material_model,ndm,ndf,nel &
                  ,ta(l),tl0,dt,xl,ul,fi(1,1,l),finv(1,1,l),detf(1,l),sigv,dd,ndim,str_data,n_out)
  str_out(l,:) = str_data

!       Compute volume and mass factor

  dvol(l) = dvol0(l)*detf(1,l)


!         Multiply tangent moduli and stresses by volume element

  do i = 1, ndim
    sigv(i) = sigv(i)*dvol(l)
    do j = 1, ndim
      dd(i,j) = dd(i,j)*dvol(l)
    end do ! j
  end do ! i


	if (dyna) then			
		al = 0.0
		
		do j = 1, nel
			do k = 1, 2
				al(k) = al(k) + shp(3,j,l)*(accl(k,j) + infv(k))
			enddo		
		enddo
		al = cfac*rr*al		
		aj0 = lfac*rr			        	
	else
		aj0 = 0.0
		al = 0.0
	endif


!         COMPUTE STRESS DIVERGENCE AND INERTIA TERMS

  do i = 1,nel

!           Stress divergence term (used in geometric stiffness)

    r1(1,i) = shp(1,i,l)*sigv(1) + shp(2,i,l)*sigv(4)
    r1(2,i) = shp(1,i,l)*sigv(4) + shp(2,i,l)*sigv(2)

    r(1,i) = r(1,i) - r1(1,i) - shpr(i)*sigv(3) + (b1 -al(1) - aj0*(infv(1)+accl(1,i)))*dvol0(l)*shp(3,i,l)
    r(2,i) = r(2,i) - r1(2,i) + (b2 -al(2) - aj0*(infv(2)+accl(2,i)))*dvol0(l)*shp(3,i,l)
  end do ! i
	


!         COMPUTE K (s(nst,nst) = K)

!           PART 1. - Geometric part.

  if(geo_flag) then

    i1  = 0
    do i = 1, nel
      bd3 = shpr(i)*sigv(3)
      j1  = 0
      do j = 1, nel
        bdb  = r1(1,i)*shp(1,j,l) + r1(2,i)*shp(2,j,l)
        s(i1+1,j1+1) = s(i1+1,j1+1) + bdb + bd3*shpr(j)
        s(i1+2,j1+2) = s(i1+2,j1+2) + bdb
        
        j1 = j1 + ndf
      end do ! j
      
      i1 = i1 + ndf
    end do ! i
  endif

!           PART 2. - Tangent modulus part (based upon dd-array)

  i1 = 0
  do i  = 1, nel
!             Compute bmat-t * dd * dvol
    do jj = 1, ndim
      bbd(jj,1) = shp(1,i,l)*dd(1,jj) + shpr(i)*dd(3,jj) + shp(2,i,l)*dd(4,jj)
      bbd(jj,2) = shp(1,i,l)*dd(4,jj) + shp(2,i,l)*dd(2,jj)
    end do ! jj

!             Compute tangent stiffness
    j1 = 0
    do j  = 1, nel

      do k = 1, 2
        s(i1+k,j1+1) = s(i1+k,j1+1) + bbd(1,k)*shp(1,j,l) + bbd(3,k)*shpr(j)  &
                      + bbd(4,k)*shp(2,j,l)
        s(i1+k,j1+2) = s(i1+k,j1+2) + bbd(4,k)*shp(1,j,l) + bbd(2,k)*shp(2,j,l)
      end do ! k

      j1 = j1 + ndf
    end do ! j
    i1 = i1 + ndf
  end  do ! i

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
      do k = 1, ndf*nel, ndf
        m(j,k) = mass(j1,k1)
        m(j+1,k+1) = mass(j1,k1)
        k1 = k1 + 1
      end do
      j1 = j1 + 1
	end do
	
!  Add mass matrix(mass) to stiffness matrix(s)
	s = ct1*s + ct3*m
	
else
	s = ct1*s   ! thickness
endif  ! dyna

j = 1
do i = 1,nel
  p(j) = r(1,i)*thick
  p(j+1) = r(2,i)*thick
  j = j + 2
enddo

!write(*,*) 'p: ', p
!write(*,*) 's: ', s


end subroutine nl_2d_fn
