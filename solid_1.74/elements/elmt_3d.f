subroutine nl_elmt3d(action,iter,pro_type,reg_num,ele_num,ul,xl,tl,tl0,s,p,str_out,ndf,ndm,nst,nel,nquad,n_out)

! Coded by Jeeho Lee (June 2015)

!  Compute the effective residual vector(p) and tangent stiffness matrix(s)
!    - Dynamic/Static procedures	
!    - 1X1X1 or 2X2X2 Integration
!    - 3D element

!     pro_type: 3D problems = 3

use element_specification

! Global variables from module 'elements' & 'element_specification':
!
!  - Geometry: thickness, body(3), gr_accl(3), density, consistent_ratio, damping_factor(2)
!  - Time related: dt, ctan(3)

implicit none
      
integer, intent(in) :: action, pro_type, reg_num, ele_num, ndf, ndm, nst, nel, nquad, iter, n_out
real, intent(in) :: xl(ndm,nel), ul(ndf,5*nel), tl(nel), tl0
real, intent(out) :: s(nst,nst), p(nst), str_out(nquad,n_out)

integer, parameter :: ndim = 6, max_lint = 8, max_nel = 8

integer :: mat_num, material_model, point_address

integer :: nea, nea2      
integer :: i, j, k, l, i1, j1, k1, lint
real :: xn, yn, zn
real :: a11, a12, a13, a21, a22, a23, a31, a32, a33, a41, a42, a43, a51, a52, a53, a61, a62, a63

real :: dot, shp(4,max_nel,max_lint), sg(4,max_lint), el(4,7), tg(16), dvol, xsj(max_lint), dv(max_lint)
real :: sig(6), eps(6), dd(6,6), secd(ndim,ndim)
      	  
real :: lfac, cfac, mass(nel,nel), vl(ndf,nel), accl(ndf,nel), al(3), aj0, aj1, aj2, aj3, aj4
real :: ug, infv(3), damfac_k, damfac_m, damp(nst,nst), m(nst,nst)
real :: z11,z12,z21,z22,z41,z42,ct1,ct2,ct3
real :: thick, b1, b2, b3, rr, ta
real :: str_data(n_out)

logical :: mat_damping = .FALSE., flgrm = .FALSE., dyna = .FALSE.
logical :: ther, hexa
	

if (pro_type /= 3) STOP 'nl_elmt3d: pro_type is not for 3D!'
if (nel > max_nel) STOP 'nl_elmt3d: # of nodes per element exceeds maximum!'
if (ndf /= 3) STOP 'nl_elmt3d: nodal dof is not 3!'

if (nel == 8) then
  hexa = .TRUE.
else
  STOP 'nl_elmt3d: hexahedron only presently!'
  hexa = .FALSE.
endif

if (action == 2) then
	dyna = .TRUE.
else
	dyna = .FALSE.
endif


ther = .FALSE.

nea = 4*nel
nea2 = 3*nel

vl = ul(1:ndf,nea2+1:nea)
      
s = 0.0
p = 0.0

!     read data from module 'element_specification'
thick = 1.0   ! 3D
b1 = body(1)
b2 = body(2)
b3 = body(3)
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
	flgrm = .TRUE.
  if (flgrm) then
		infv(1) = gr_accl(1)
		infv(2) = gr_accl(2)
    infv(3) = gr_accl(3)
	else
		infv = 0.0
	endif

else
	ct1 = thick
	ct2 = 0.0
	ct3 = 0.0
	mat_damping = .FALSE.

endif


call element_material(reg_num, ele_num, mat_num, material_model, .TRUE.)

mat_damping = .FALSE.
! if(material_model == 5) then
! 	mat_damping = .TRUE.
! endif

if (hexa) then
	if (nquad == 8) then
		l = 2
	elseif (nquad == 1) then
		l = 1
	else
		STOP 'nl_elmt_3d: invalid nquad!'
	endif	
		
  call int3d(l,lint,sg,nel)
endif


if (lint > max_lint) then
  write(*,*) 'nl_elmt_3d: too many quadrature points!', lint
endif


! Numerical integration at Gauss points

do l = 1, lint
!		Shape functions and derivatives
  call shp3d(sg(1,l),xl,shp(1,1,l),xsj(l),ndm,nel)
  dv(l) = xsj(l)*sg(4,l)
end do  ! l


do l = 1, lint
!   Get the address in history db using the subrountine in 'physical_domain'
  call find_point_number(reg_num,ele_num,l,point_address)

!   compute the strain
  eps = 0.0
  ta  = 0.0

  do j  = 1, nel
    eps(1) = eps(1) + shp(1,j,l)*ul(1,j)
    eps(2) = eps(2) + shp(2,j,l)*ul(2,j)
    eps(3) = eps(3) + shp(3,j,l)*ul(3,j)
    eps(4) = eps(4) + shp(1,j,l)*ul(2,j) + shp(2,j,l)*ul(1,j)
    eps(5) = eps(5) + shp(2,j,l)*ul(3,j) + shp(3,j,l)*ul(2,j)
    eps(6) = eps(6) + shp(3,j,l)*ul(1,j) + shp(1,j,l)*ul(3,j)
    ta     = ta + shp(4,j,l)*tl(j)    ! temperature at this quadrature point
  enddo

  call str_tang_infn(action,iter,pro_type,point_address,mat_num,material_model,ndm,ndf,nel &
                     ,ta,tl0,dt,xl,ul,eps,sig,dd,secd,ndim,str_data,n_out)
                     
  str_out(l,:) = str_data

  dvol = dv(l)

!  Compute accelerations of relative motion and ground motion
  if (dyna) then
    al = 0.0
    do j = 1, nel
      do k = 1, 3
        al(k) = al(k) + shp(4,j,l)*(accl(k,j)+infv(k))
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
  do j = 1, nel
    aj1 = shp(1,j,l)*dvol
    aj2 = shp(2,j,l)*dvol
    aj3 = shp(3,j,l)*dvol
    aj4 = shp(4,j,l)*dvol

    p(j1)   = p(j1)   + thick*(-aj1*sig(1)  -aj2*sig(4)  -aj3*sig(6)  &
						 + (b1 -al(1) -aj0*(infv(1)+accl(1,j)))*aj4)
    p(j1+1) = p(j1+1) + thick*(-aj1*sig(4)  -aj2*sig(2)  -aj3*sig(5)  &
						 + (b2 -al(2) -aj0*(infv(2)+accl(2,j)))*aj4)
    p(j1+2) = p(j1+2) + thick*(-aj1*sig(6)  -aj2*sig(5)  -aj3*sig(3)  &
						 + (b3 -al(3) -aj0*(infv(3)+accl(3,j)))*aj4)
    j1 = j1 + ndf

  end do  ! j


!   Stiffness computation


  j1 = 1
  do j = 1, nel

    xn  = shp(1,j,l)*dvol
    yn  = shp(2,j,l)*dvol
    zn  = shp(3,j,l)*dvol

    a11 = dd(1,1)*xn + dd(1,4)*yn + dd(1,6)*zn
    a21 = dd(2,1)*xn + dd(2,4)*yn + dd(2,6)*zn
    a31 = dd(3,1)*xn + dd(3,4)*yn + dd(3,6)*zn
    a41 = dd(4,1)*xn + dd(4,4)*yn + dd(4,6)*zn
    a51 = dd(5,1)*xn + dd(5,4)*yn + dd(5,6)*zn
    a61 = dd(6,1)*xn + dd(6,4)*yn + dd(6,6)*zn
    a12 = dd(1,2)*yn + dd(1,4)*xn + dd(1,5)*zn
    a22 = dd(2,2)*yn + dd(2,4)*xn + dd(2,5)*zn
    a32 = dd(3,2)*yn + dd(3,4)*xn + dd(3,5)*zn
    a42 = dd(4,2)*yn + dd(4,4)*xn + dd(4,5)*zn
    a52 = dd(5,2)*yn + dd(5,4)*xn + dd(5,5)*zn
    a62 = dd(6,2)*yn + dd(6,4)*xn + dd(6,5)*zn
    a13 = dd(1,3)*zn + dd(1,5)*yn + dd(1,6)*xn
    a23 = dd(2,3)*zn + dd(2,5)*yn + dd(2,6)*xn
    a33 = dd(3,3)*zn + dd(3,5)*yn + dd(3,6)*xn
    a43 = dd(4,3)*zn + dd(4,5)*yn + dd(4,6)*xn
    a53 = dd(5,3)*zn + dd(5,5)*yn + dd(5,6)*xn
    a63 = dd(6,3)*zn + dd(6,5)*yn + dd(6,6)*xn

    i1 = 1
    do i = 1,nel

      xn   = shp(1,i,l)
      yn   = shp(2,i,l)
      zn   = shp(3,i,l)

      s(i1  ,j1  ) = s(i1  ,j1  ) + xn*a11 + yn*a41 + zn*a61
      s(i1  ,j1+1) = s(i1  ,j1+1) + xn*a12 + yn*a42 + zn*a62
      s(i1  ,j1+2) = s(i1  ,j1+2) + xn*a13 + yn*a43 + zn*a63
      s(i1+1,j1  ) = s(i1+1,j1  ) + yn*a21 + xn*a41 + zn*a51
      s(i1+1,j1+1) = s(i1+1,j1+1) + yn*a22 + xn*a42 + zn*a52
      s(i1+1,j1+2) = s(i1+1,j1+2) + yn*a23 + xn*a43 + zn*a53
      s(i1+2,j1  ) = s(i1+2,j1  ) + zn*a31 + yn*a51 + xn*a61
      s(i1+2,j1+1) = s(i1+2,j1+1) + zn*a32 + yn*a52 + xn*a62
      s(i1+2,j1+2) = s(i1+2,j1+2) + zn*a33 + yn*a53 + xn*a63

      i1 = i1 + ndf
    end do ! i

    j1 = j1 + ndf
  end do ! j



  if (dyna) then
		
			dvol = dvol*rr

			do j = 1, nel
			
!			Compute lumped mass matrix
				aj0 = shp(4,j,l)*dvol
				mass(j,j) = mass(j,j) + lfac*aj0
!				write(*,*) 'j, lfac, aj0, dvol, rr =', j, lfac, aj0, dvol, rr
			
!			Compute consistent mass matrix
				do k = 1, nel
					mass(j,k)    = mass(j,k) + cfac*aj0*shp(4,k,l)
				end do ! k
			end do ! j
  endif  ! dyna
		
end do  ! l


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
    end do  ! k
    j1 = j1 + 1
	end do  ! j


!  Compute the effective stiffness matrix
		  
	s = ct1*s + ct3*m

else
	s = ct1*s   ! thickness
endif  ! dyna

!write(*,*) 'ct =', ct1, ct2, ct3
!write(*,*) 'mass m=', m
!write(*,'(24E10.1)') s


end subroutine nl_elmt3d
