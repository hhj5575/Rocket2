subroutine nl_elmt2d(action,iter,pro_type,reg_num,ele_num,ul,xl,tl,tl0,s,p,str_out,ndf,ndm,nst,nel,nquad,n_out)

! Coded by Jeeho Lee (Nov 2009)
!		- Axi-symmetry part (Dynamic) included: Mar 2011
!		- Include Ground acceleration: June 2011
				
!  Compute the effective residual vector(p) and tangent stiffness matrix(s)
!    - Dynamic/Static procedures	
!    - 2X2 Integration
!    - Plane stress/strain element
!    - Axi-symmetric element

!     pro_type: plane strain = 1, axisymmetric = 2	  

use element_specification
! Global variables from module 'elements' & 'element_specification':
!
!  - Geometry: thickness, body(3), gr_accl(3), density, consistent_ratio, damping_factor(2)
!  - Time related: dt, ctan(3)

implicit none
      
integer, intent(in) :: action, pro_type, reg_num, ele_num, ndf, ndm, nst, nel, nquad, iter, n_out
real, intent(in) :: xl(ndm,nel), ul(ndf,5*nel), tl(nel), tl0
real, intent(out) :: s(nst,nst), p(nst), str_out(nquad,n_out)

integer, parameter :: ndim = 4, max_lint = 9, max_nel = 9	  
integer :: mat_num, material_model, point_address

integer :: nea, nea2      
integer :: i,j,k,l,i1,j1,k1,lint
real :: dv, rn, zn, un
real :: a11, a12, a21, a22, a41, a42
real :: bd11, bd21, bd12, bd22, bd13, bd23, bd14, bd24

real :: dot, shp(3,max_nel), sg(3,max_lint), el(4,7), tg(16), dvol, xsj0, shpr(max_nel)
real :: sig(4), eps(4), dd(4,4), secd(ndim,ndim)
      	  
real :: lfac, cfac, mass(nel,nel), vl(ndf,nel), accl(ndf,nel), al(2), aj0, aj1, aj2, aj3
real :: ug, infv(2), damfac_k, damfac_m, damp(nst,nst), m(nst,nst)
real :: z11,z12,z21,z22,z41,z42,ct1,ct2,ct3
real :: thick, b1, b2, rr, radius, ta
real :: str_data(n_out)

logical :: mat_damping = .false., flgrm = .false., dyna = .false.
logical :: ther, quad
	


if (nel > max_nel) stop 'nl_elmt2d: # of nodes per element exceeds maximum!'

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

l   = 2
nea = 4*nel
nea2 = 3*nel

vl = ul(1:ndf,nea2+1:nea)
      
s = 0.0
p = 0.0

!     read data from module 'element_specification'
thick = thickness
b1 = body(1)
b2 = body(2)
rr = density

accl = 0.0
infv = 0.0	  
	  
if (dyna) then
!     read data from module 'element_specification'
	cfac = consistent_ratio
	
	damfac_m = damping_factor(1)   ! mass matrix proportional damping factor	
	damfac_k = damping_factor(2)   ! stiffness matrix proportional damping factor
	
	accl = ul(1:ndf,nea+1:nea+nel)

  lfac = 1.d0 - cfac

!  Zero mass & damping matrix
  mass = 0.0
  damp = 0.0

	ct1 = thick*ctan(1)
	ct2 = thick*ctan(2)
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
	ct2 = 0.0
	ct3 = 0.0
	mat_damping = .false.

endif


call element_material(reg_num, ele_num, mat_num, material_model, .TRUE.)

mat_damping = .false.
if(material_model == 5) then
	mat_damping = .true.
endif

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
  write(*,*) 'nl_elmt_2d: Too many quadrature points!', lint
endif


! Numerical integration at Gauss points

do l = 1,lint

!		Shape functions and derivatives
  if(quad) then
    call shp2d(sg(1,l),xl,shp,xsj0,ndm,nel,.false.)
    xsj0 = xsj0*sg(3,l)  ! sg(3,l): quadrature weights
  else
    call trishp(el(1,l),xl,ndm,nel-4,xsj0,shp)
    xsj0 = xsj0*el(4,l)  ! el(4,l): quadrature weights
  endif


!   Get the address in history db using the subrountine in 'physical_domain'
  call find_point_number(reg_num,ele_num,l,point_address)

  eps = 0.0
  ta  = 0.0
  radius = 0.0   ! for axi-symmetry

  do j  = 1, nel
    radius = radius + shp(3,j)*xl(1,j)
    eps(1) = eps(1) + shp(1,j)*ul(1,j)
    eps(2) = eps(2) + shp(2,j)*ul(2,j)
    eps(3) = eps(3) + shp(3,j)*ul(1,j)   ! for axi-symmetry
    eps(4) = eps(4) + shp(1,j)*ul(2,j) + shp(2,j)*ul(1,j)
    ta     = ta + shp(3,j)*tl(j)    ! temperature at this quadrature point
  enddo

  if (pro_type == 2) then   ! for axi-symmetry
    eps(3) = eps(3)/radius
  else                      ! for plane strain
    radius = 1.0
    eps(3) = 0.0 
  endif

  call str_tang_infn(action,iter,pro_type,point_address,mat_num,material_model,ndm,ndf,nel &
                  ,ta,tl0,dt,xl,ul,eps,sig,dd,secd,ndim,str_data,n_out)

  str_out(l,:) = str_data
  
  dvol = xsj0*radius    ! radius = 1.0 for plane strain/stress cases
  
  if (pro_type == 2) then
    do j = 1, nel
      shpr(j) = shp(3,j)/radius
    end do ! j
  else
    xsj0 = 0.0
  endif

!  Compute accelerations of relative motion and ground motion
  if (dyna) then      
    al = 0.0          
    do j = 1,nel
      do k = 1,2
        al(k) = al(k) + shp(3,j)*(accl(k,j)+infv(k))
      enddo
    enddo		      
    al = cfac*rr*al		
    aj0 = lfac*rr		
  else
    aj0 = 0.0
    al = 0.0
  endif		
    
!   Residual computation
    
  j1 = 1
  do j = 1,nel
    rn = shp(1,j)*dvol
    zn = shp(2,j)*dvol
    un = shp(3,j)*dvol

    p(j1) = p(j1) + thick*(-rn*sig(1) -sig(4)*zn -shp(3,j)*xsj0*sig(3)  &
           + (b1 -al(1) -aj0*(infv(1)+accl(1,j)))*un) 
    p(j1+1) = p(j1+1) + thick*(-zn*sig(2) -sig(4)*rn  &
           + (b2 -al(2) -aj0*(infv(2)+accl(2,j)))*un) 
    j1 = j1 + ndf
  end do  ! j


!   Stiffness computation


  j1 = 1
  do j = 1, nel
    aj1 = shp(1,j)*dvol
    aj2 = shp(2,j)*dvol
    aj3 = shp(3,j)*xsj0

!			Compute B_trans * D * j * w
    bd11 = aj1*dd(1,1) + aj3*dd(3,1) + aj2*dd(4,1)
    bd21 = aj2*dd(2,1) + aj1*dd(4,1)
    bd12 = aj1*dd(1,2) + aj3*dd(3,2) + aj2*dd(4,1)
    bd22 = aj2*dd(2,2) + aj1*dd(4,2)
    bd13 = aj1*dd(1,3) + aj3*dd(3,3) + aj2*dd(4,3)
    bd23 = aj2*dd(2,3) + aj1*dd(4,3)
    bd14 = aj1*dd(1,4) + aj3*dd(3,4) + aj2*dd(4,4)
    bd24 = aj2*dd(2,4) + aj1*dd(4,4)

    k1 = 1
    do k = 1, nel
      s(j1,k1)     = s(j1,k1)     + bd11*shp(1,k) + bd14*shp(2,k) + bd13*shpr(k)
      s(j1,k1+1)   = s(j1,k1+1)   + bd12*shp(2,k) + bd14*shp(1,k)
      s(j1+1,k1)   = s(j1+1,k1)   + bd21*shp(1,k) + bd24*shp(2,k) + bd23*shpr(k)
      s(j1+1,k1+1) = s(j1+1,k1+1) + bd22*shp(2,k) + bd24*shp(1,k)

      k1 = k1 + ndf
    end do ! k
    j1 = j1 + ndf
  end do ! j



  if (dyna) then
  
    dv = dvol*rr

    do j = 1, nel
    
!			Compute lumped mass matrix
      aj0 = shp(3,j)*dv
      mass(j,j) = mass(j,j) + lfac*aj0
!				write(*,*) 'j, lfac, aj0, dv, rr =', j, lfac, aj0, dv, rr
    
!			Compute consistent mass matrix
      do k = 1, nel
        mass(j,k)    = mass(j,k) + cfac*aj0*shp(3,k)
      end do ! k
    end do ! j

!   Compute material damping matrix                
    if (mat_damping) then
      j1 = 1
      do j = 1,nel
        rn  = shp(1,j)*dvol
        zn  = shp(2,j)*dvol

        z11 = secd(1,1)*rn + secd(1,4)*zn
        z21 = secd(2,1)*rn + secd(2,4)*zn
        z41 = secd(4,1)*rn + secd(4,4)*zn
        z12 = secd(1,2)*zn + secd(1,4)*rn
        z22 = secd(2,2)*zn + secd(2,4)*rn
        z42 = secd(4,2)*zn + secd(4,4)*rn

        i1 = 1
        do i = 1,nel
          rn = shp(1,i)
          zn = shp(2,i)			

          damp(i1  ,j1  ) = damp(i1  ,j1  ) + rn*z11 + zn*z41
          damp(i1  ,j1+1) = damp(i1  ,j1+1) + rn*z12 + zn*z42
          damp(i1+1,j1+1) = damp(i1+1,j1+1) + zn*z22 + rn*z42
          damp(i1+1,j1  ) = damp(i1+1,j1  ) + zn*z21 + rn*z41

          i1 = i1 + ndf
        end do	! i		  
        j1 = j1 + ndf
      end do	! j
    endif  ! mat_damping
    
  endif  ! dyna
		
end do  ! l


if (dyna) then
!  Expand 'mass' matrix (scalar mass: no directional mass) to the general m(nst,nst) matrix

	m = 0.0
	j1 = 1
	do j = 1,ndf*nel,ndf		  
      k1 = 1
      do k = 1,ndf*nel,ndf
        m(j,k) = mass(j1,k1)
        m(j+1,k+1) = mass(j1,k1)
        k1 = k1 + 1
      end do
      j1 = j1 + 1
	end do


! Add viscous damping (Stiffness proportional) with material damping - matrix operation (nst,nst)

	damp = damp + damfac_k*s + damfac_m*m


!  Add residual from damping to p
	j1 = 1 
	do j = 1,nel
		i1 = 1
		do i = 1,nel
			p(j1) = p(j1) - thick*(damp(j1,i1)*vl(1,i) + damp(j1,i1+1)*vl(2,i))
			p(j1+1) = p(j1+1) - thick*(damp(j1+1,i1)*vl(1,i) + damp(j1+1,i1+1)*vl(2,i))
				 
			i1 = i1 + 2
		enddo
		j1 = j1 + 2
	enddo		

!  Add mass matrix(mass) & damping matrix(damp) to stiffness matrix(s)
		  
	s = ct1*s + ct2*damp + ct3*m

else
	s = ct1*s   ! thickness
endif  ! dyna

!write(*,*) 'ct =', ct1, ct2, ct3
!write(*,*) 'mass m=', m
					
end subroutine nl_elmt2d
