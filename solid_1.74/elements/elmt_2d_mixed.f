subroutine nl_2d_mixed(action,iter,pro_type,reg_num,ele_num,ul,xl,tl,tl0,s,p,str_out,ndf,ndm,nst,nel,nquad,n_out)

! Coded by Jeeho Lee (April 2010)
!		- Dynamic part (add Mass & Ground acceleration): Aug 2011


!     Plane/axisymmetric linear element routine - B-bar formulation
!     4-node element formulation
!     pro_type: plane strain = 1, axisymmetric = 2	  

use element_specification

! Global variables (previously COMMON block data) from module 'element_specification':
!
!  - Geometry: thickness, body(3), gr_accl(3), density, consistent_ratio
!  - Time related: dt, ctan(3)

implicit none

integer, intent(in) :: action, pro_type, reg_num, ele_num, ndf, ndm, nst, nel, nquad, iter, n_out
real, intent(in) :: xl(ndm,nel), ul(ndf,5*nel), tl(nel), tl0
real, intent(out) :: s(nst,nst), p(nst), str_out(nquad,n_out)

integer, parameter :: ndim = 4, max_lint = 4, max_nel = 9
integer ::  mat_num, material_model, point_address
      
integer :: i,i1,j,jj,j1,k,k1,l,lint, nea
integer :: ncp,ncs,npm
real :: epp,xsj0,xsj
real :: dsigtr, mpress, dtheta, fac
real :: sg(3,max_lint), r(ndf,nel), xx(2,max_lint), shp(3,16,max_lint)
real :: bbd(3,7), aa(6,6,max_lint), dd(7,7), ad(ndim,ndim), secd(ndim,ndim), dvol(max_lint)
real :: sigm(6), sigl(6,max_lint), bpra(3), bbar(2,16,max_lint)
real :: phi(6,max_lint), theta(3,max_lint), hh(6,6)
real :: press(max_lint), pbar(max_lint), eps(6,max_lint)
real :: irad(max_lint), el(4,7), ta(max_lint), epsv(max_lint)
real :: xr(2)

real :: b1, b2, thick
real :: str_data(n_out)

real :: mass(nel,nel), al(2), accl(ndf,nel), infv(2), lfac, cfac, rr, m(nst,nst)
real :: aj0, ct1, ct3

logical :: flgrm = .false., dyna = .false.
logical :: ther, quad



if ((nel > 4) .or. (nel < 4)) stop 'nl_2d_mixed: Too many nodes for MIXED mode - valid only for 4 node element!'

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

!		Set element quadrature order
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


call element_material(reg_num, ele_num, mat_num, material_model, .TRUE.)

if (quad) then
  l = 2
  call int2d(l,lint,sg,nel)
else    ! triangular elements
  if (nel == 3) then
    if(pro_type == 1) then
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
  write(*,*) 'nl_2d_fn: Too many quadrature points!', lint
endif

! Numerical integration at Gauss points

do l = 1,lint

!         Shape functions and derivatives
  if(quad) then
    call shp2d(sg(1,l),xl,shp(1,1,l),xsj0,ndm,nel,.false.)
    dvol(l) = xsj0*sg(3,l)  ! sg(3,l): quadrature weights
  else
    call trishp(el(1,l),xl,ndm,nel-4,xsj0,shp(1,1,l))
    dvol(l) = xsj0*el(4,l)
  endif
  
!         Mixed volume effect and temperature projection

  ta(l) = 0.0

!         Compute coordinates

  xx(1,l) = 0.0d0
  xx(2,l) = 0.0d0
  do i = 1,nel
    xx(1,l) = xx(1,l) + shp(3,i,l)*xl(1,i)
    xx(2,l) = xx(2,l) + shp(3,i,l)*xl(2,i)
  end do ! i

!         Compute volumetric strain from displacements

  if (pro_type == 2) then
    dvol(l) = dvol(l)*xx(1,l)
    irad(l) = 1.d0/xx(1,l)
    ncp     = 2
    ncs     = 4
  else
    irad(l) = 0.0d0
    ncp     = 2
    ncs     = 4
  endif
  
  theta(:,l) = 0.0

  do i = 1,nel
    fac        = shp(1,i,l) + shp(3,i,l)*irad(l)
    theta(1,l) = theta(1,l) + fac*ul(1,i) + shp(2,i,l)*ul(2,i)
    theta(2,l) = 0.0
    theta(3,l) = 0.0

    ta(l)      = ta(l) + shp(3,i,l)*tl(i)
  end do ! i
!         Set the pressure functions
  phi(1,l) = 1.0
  
end do ! l



!       Mixed model for volumetric and temperature response

call bbar2s(phi,shp,dvol,lint,nel,npm,hh,irad,theta,bbar)


!       Compute strains and stresses at quadrature points

do l = 1,lint
  call strn2m(shp(1,1,l),xl,ul,theta(1,l),irad(l),ndm,ndf,nel,eps(1,l))

  epsv(l) = theta(1,l)
  
!     Get the address in history db using the subrountine in 'physical_domain'
  call find_point_number(reg_num,ele_num,l,point_address)

  call str_tang_infn(action,iter,pro_type,point_address,mat_num,material_model,ndm,ndf,nel &
                     ,ta(l),tl0,dt,xl,ul,eps(1,l),sigl(1,l),ad,secd,ndim,str_data,n_out)

  aa(:,:,l) = 0.0
  aa(1:ndim,1:ndim,l) = ad
  str_out(l,:) = str_data
     
!         Volumetric stress

  pbar(l) = (sigl(1,l) + sigl(2,l) + sigl(3,l))/3.0

end do ! l

!       Integrate constant pressure over volume

mpress    = 0.0
do l = 1,lint
  mpress  = mpress  + pbar(l)*dvol(l)
end do ! l

!         Divide pressure by volume

press(1) = mpress * hh(1,1)
do l = 2,lint
  press(l) = press(1)
end do ! l


!       Compute mixed stress and multiply by volume element
do l = 1,lint
  dsigtr    =  press(l)  - pbar(l)
  sigl(1,l) =  sigl(1,l) + dsigtr
  sigl(2,l) =  sigl(2,l) + dsigtr
  sigl(3,l) =  sigl(3,l) + dsigtr
end do ! l

! Residual and Tangent computations

do l = 1,lint

  do j =1,ndim
    sigm(j) = sigl(j,l)*dvol(l)
  enddo

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
  do j = 1,nel
    r(1,j) = r(1,j) + (b1 -al(1) - aj0*(infv(1)+accl(1,j)))*shp(3,j,l)*dvol(l)  &
						- shp(1,j,l)*sigm(1) - shp(2,j,l)*sigm(4) - shp(3,j,l)*sigm(3)*irad(l)
    r(2,j) = r(2,j) + (b2 -al(2) - aj0*(infv(2)+accl(2,j)))*shp(3,j,l)*dvol(l)  &
						- shp(2,j,l)*sigm(2) - shp(1,j,l)*sigm(4)
  end do

!  Compute mixed tangent stiffness matrix
!      Multiply tangent moduli by volume element

  call dmatmx(aa(1,1,l), dd)
  
  do i = 1,7
    do j = 1,7
      dd(i,j) = dd(i,j)*dvol(l)
    end do ! j
  end do ! i

!               Compute row terms
  i1 = 0
  do i = 1,nel
!                 Compute bmat-t * dd * dvol
    do jj = 1,7

      bbd(1,jj) = shp(1,i,l)*dd(1,jj) + shp(2,i,l)*dd(4,jj)  &
                  + shp(3,i,l)*dd(3,jj)*irad(l) + bbar(1,i,l)*dd(7,jj)

      bbd(2,jj) = shp(2,i,l)*dd(2,jj) + shp(1,i,l)*dd(4,jj)  &
                  + bbar(2,i,l)*dd(7,jj)

      bbd(3,jj) =  shp(2,i,l)*dd(5,jj) + (shp(1,i,l)-shp(3,i,l)*irad(l))*dd(6,jj)
    end do ! jj
			                          
!		Compute: B_bar_trans * D * alpha * j * w
    j1 = 0
    do j = 1,nel
			do jj = 1,ncp ! ncp = 2
           s(i1+jj,j1+1) = s(i1+jj,j1+1)  &
                          + bbd(jj,1)*shp(1,j,l)  &
                          + bbd(jj,4)*shp(2,j,l)  &
                          + bbd(jj,3)*shp(3,j,l)*irad(l)  &
                          + bbd(jj,7)*bbar(1,j,l)

           s(i1+jj,j1+2) = s(i1+jj,j1+2)  &
                          + bbd(jj,2)*shp(2,j,l)  &
                          + bbd(jj,4)*shp(1,j,l)  &
                          + bbd(jj,7)*bbar(2,j,l)
			end do ! jj				
			j1 = j1 + ndf
    end do ! j
    i1 = i1 + ndf
  end do ! i


	if (dyna) then
		do i = 1,nel
		
!		Compute lumped mass matrix
			aj0 = shp(3,i,l)*dvol(l)*rr
			mass(i,i) = mass(i,i) + lfac*aj0
			                          
!		Compute consistent mass matrix
			do j = 1,nel
				mass(i,j)    = mass(i,j) + cfac*aj0*shp(3,j,l)				
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
	
!  Compute the effective stiffness matrix
		  
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


end subroutine nl_2d_mixed
