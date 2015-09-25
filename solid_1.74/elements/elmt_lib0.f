subroutine int2d(l,lint,sg,nel)

!-----------------------------------------------------------------
!      Purpose: Form Gauss points and weights for two dimensions

!      Inputs:
!         l       - Number of points/direction

!      Outputs:
!         lint    - Total number of points
!         sg(3,*) - Array of points and weights
!-----------------------------------------------------------------

implicit none

integer :: i,j,k,l,lint, lr(9),lz(9),lw(9),nel
real :: g,h,sg(3,*),ss(5),ww(5)

data lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
data lw/4*25,4*40,64/


!     Set number of total points

lint = l*l

!     5 pt. integration

if(l == 0) then
        lint = 5
        g    = sqrt(0.6)		! sqrt(0.6)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 5.0/9.0
        end do ! i

        sg(1,5) = 0.0d0
        sg(2,5) = 0.0d0
        sg(3,5) = 16.0/9.0		! 16/9

!     1x1 integration

elseif(l == 1) then
        sg(1,1) = 0.d0
        sg(2,1) = 0.d0
        if(nel.eq.3) sg(2,1) = -1.0/3.0
        sg(3,1) = 4.d0

!     2x2 integration

elseif(l == 2) then
        g = 1.0/sqrt(3.0)		! sqrt(1/3)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 1.d0
        end do ! i

!     3x3 integration

elseif(l == 3) then
        g = sqrt(0.6)		! sqrt(0.6)
        h = 1.d0
        h = h/81.d0
        do i = 1,9
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = h*lw(i)
        end do ! i

else
	stop 'int2d: Quadrature points exceed the maximum(3X3)!'
        
endif

end subroutine int2d


	  
	  
	  

subroutine shp2d(ss,xl,shp,xsj,ndm,nel,flg)

!-----------------------------------------------------------------
!      Purpose: Computes shape function and derivatives for
!               quadrilateral elements

!      Inputs:
!         ss(2)     - Natural coordinates for point
!          xl(ndm,*) - Nodal coordinates for element
!         ndm       - Spatial dimension of mesh
!         nel       - Number of nodes on element
!         flg       - Flag, compute global x/y derivatives if false,
!                           else derivatives are w/r natural coords.

!      Outputs:
!         shp(3,*)  - Shape functions and derivatives at point
!                     shp(1,i) = dN_i/dx or dN_i/dxi_1
!                     shp(2,i) = dN_i/dy or dN_i/dxi_2
!                     shp(3,i) = N_i
!         xsj       - Jacobian determinant at point
!-----------------------------------------------------------------

implicit none

logical :: flg
integer :: ndm,nel, i,j,k
real :: xsj, temp
real :: shp(3,9),xl(ndm,*), s(4),t(4),xs(3,2),sx(2,2),ss(2)

data s/-0.5d0,0.5d0,0.5d0,-0.5d0/,t/-0.5d0,-0.5d0,0.5d0,0.5d0/  
! (data s: set values of half natural coords at nodes)



if (nel > 9) stop 'shp2d: exceed maximum element nodes(9)!'

if (nel == 4 .and. .not.flg) then
!     Form 4-node quadrilateral shape functions
        call shapef(ss(1),ss(2),xl,shp,xsj,ndm,flg)
        
else
	do i = 1,4
		shp(3,i) = (0.5d0+s(i)*ss(1))*(0.5d0+t(i)*ss(2))
		shp(1,i) = s(i)*(0.5d0+t(i)*ss(2))
		shp(2,i) = t(i)*(0.5d0+s(i)*ss(1))
	end do ! i
	
!       Form triangle by adding third and fourth together

	if (nel == 3) then
		do i = 1,3
			shp(i,3) = shp(i,3)+shp(i,4)
		end do ! i
	end if

!       Add quadratic terms if necessary

	if (nel > 4) then
		call shap2(ss(1),ss(2),shp,nel)
	endif

!       Construct jacobian and its inverse

	xs = 0.0
	do i = 1,ndm
		do j = 1,2
			do k = 1,nel
				xs(i,j) = xs(i,j) + xl(i,k)*shp(j,k)
			end do ! k
		end do ! j
	end do ! i
        
	if(ndm == 2) then
		xsj = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
	elseif(ndm == 3) then
		xsj = sqrt((xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1))**2  &
					+ (xs(3,1)*xs(1,2)-xs(3,2)*xs(1,1))**2  &
					+ (xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1))**2)
	endif
        
	if (.not.flg) then
		if(xsj == 0.0d0) then
			temp = 1.0d0
		else
			temp = 1.d0/xsj
		endif
          
		sx(1,1) = xs(2,2)*temp
		sx(2,2) = xs(1,1)*temp
		sx(1,2) =-xs(1,2)*temp
		sx(2,1) =-xs(2,1)*temp

!         Form global derivatives

		do i = 1,nel
			temp     = shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
			shp(2,i) = shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
			shp(1,i) = temp
		end do ! i

!         Return center node in hierarchical form for 8-nodes

		if (nel == 8) then
			temp     = shp(1,9)*sx(1,1)+shp(2,9)*sx(2,1)
			shp(2,9) = shp(1,9)*sx(1,2)+shp(2,9)*sx(2,2)
			shp(1,9) = temp
		endif
          
	endif
endif


end subroutine shp2d





subroutine shapef(s,t,xl,shp,xsj,ndm,flg)

!-----------------------------------------------------------------
!      Purpose: Shape function routine for 4-node isoparametric
!               quadrilaterals

!      Inputs:
!         s,t       - Natural coordinates of point
!         xl(ndm,*) - Nodal coordinates for element
!         ndm       - Spatial dimension of mesh
!         flg       - Flag, Compute global derivatives if true,
!                           else compute derivatives w/r natural coords.

!      Outputs:
!         shp(3,*)  - Shape functions and derivatives at point
!                     shp(1,i) = dN_i/dx  or dN_i/dxi_1
!                     shp(2,i) = dN_i/dy  or dN_i/dxi_2
!                     shp(3,i) = N_i
!         xsj       - Jacobian determinant at point
!-----------------------------------------------------------------

implicit none

logical :: flg
integer :: ndm
real :: s,t, xsj,xsj1, sh,th,sp,tp,sm,tm
real :: xo,xs,xt, yo,ys,yt
real :: xsm,xsp,xtm,xtp, ysm,ysp,ytm,ytp
real :: xl(ndm,4),shp(3,4)

!     Set up interpolations

      sh = 0.5d0*s
      th = 0.5d0*t
      sp = 0.5d0 + sh
      tp = 0.5d0 + th
      sm = 0.5d0 - sh
      tm = 0.5d0 - th
      shp(3,1) =   sm*tm
      shp(3,2) =   sp*tm
      shp(3,3) =   sp*tp
      shp(3,4) =   sm*tp

!     Set up natural coordinate functions (times 4)

      xo =  xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4)
      xs = -xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4) + xo*t
      xt = -xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4) + xo*s
      yo =  xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4)
      ys = -xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4) + yo*t
      yt = -xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4) + yo*s

!     Compute jacobian (times 16)

      xsj1 = xs*yt - xt*ys

!     Divide jacobian by 16 (multiply by .0625)

      xsj = 0.0625d0*xsj1
      if(.not.flg) then
        if(xsj1.eq.0.0d0) then
          xsj1 = 1.0d0
        else
          xsj1 = 1.0d0/xsj1
        endif

!       Divide functions by jacobian

        xs  = (xs+xs)*xsj1
        xt  = (xt+xt)*xsj1
        ys  = (ys+ys)*xsj1
        yt  = (yt+yt)*xsj1

!       Multiply by interpolations

        ytm =  yt*tm
        ysm =  ys*sm
        ytp =  yt*tp
        ysp =  ys*sp
        xtm =  xt*tm
        xsm =  xs*sm
        xtp =  xt*tp
        xsp =  xs*sp

!       Compute shape functions

        shp(1,1) = - ytm+ysm
        shp(1,2) =   ytm+ysp
        shp(1,3) =   ytp-ysp
        shp(1,4) = - ytp-ysm
        shp(2,1) =   xtm-xsm
        shp(2,2) = - xtm-xsp
        shp(2,3) = - xtp+xsp
        shp(2,4) =   xtp+xsm
      endif

end subroutine shapef
	  




subroutine shap2(s,t,shp,nel)

!-----------------------------------------------------------------
!      Purpose: Adds quadratic functions to quadrilaterals for any
!               non-zero mid-side or central node (Q8 & Q9 only)

!      Inputs:
!         s,t      - Natural coordinates
!         nel      - Maximum number of local node on element <= 9

!      Outputs:
!         shp(3,9) - Shape functions and derivatives w/r natural coords
!                    shp(1,i) = dN_i/dxi_1
!                    shp(2,i) = dN_i/dxi_2
!                    shp(3,i) = N_i
!-----------------------------------------------------------------

implicit none

integer :: i, j, k, l, nel
real :: s, t, s2, t2
real :: shp(3,9)


if (nel < 8 .or. nel > 9) stop 'shap2: only for Q8 & Q9 elements!'

s2 = (1.d0-s*s)*0.5d0
t2 = (1.d0-t*t)*0.5d0

do i = 5,9
  do j = 1,3
    shp(j,i) = 0.0d0
  end do
end do

!     Midside nodes (serendipity)

shp(1,5) = -s*(1.d0-t)
shp(2,5) = -s2
shp(3,5) = s2*(1.d0-t)
      
shp(1,6) = t2
shp(2,6) = -t*(1.d0+s)
shp(3,6) = t2*(1.d0+s)

shp(1,7) = -s*(1.d0+t)
shp(2,7) = s2
shp(3,7) = s2*(1.d0+t)

shp(1,8) = -t2
shp(2,8) = -t*(1.d0-s)
shp(3,8) = t2*(1.d0-s)

!     Interior node (lagrangian)
shp(1,9) = -4.d0*s*t2
shp(2,9) = -4.d0*t*s2
shp(3,9) =  4.d0*s2*t2

      
!     Correct edge nodes for interior node (lagrangian)
if (nel == 9) then

	do j= 1,3
		do i = 1,4
			shp(j,i) = shp(j,i) - 0.25d0*shp(j,9)
		end do
		do i = 5,8
			shp(j,i) = shp(j,i) - 0.5d0*shp(j,9)
		end do
	end do
endif

!     Correct corner nodes for presence of midside nodes
k = 8
do i = 1,4
	l = i + 4
	do j = 1,3
		shp(j,i) = shp(j,i) - 0.5d0*(shp(j,k)+shp(j,l))
	end do
	k = l
end do


end subroutine shap2






subroutine tint2d(l,lint,el)

!-----------------------------------------------------------------
!      Purpose: Set gauss points and weights for triangular elements

!      Inputs:
!         l       - Number of gauss points indicator

!      Outputs:
!         lint    - Total number of points
!         el(4,*) - Area coordinate points and weights for quadrature
!-----------------------------------------------------------------

implicit none

integer :: l, lint
real :: el(4,*)
real :: r0,r1,r2,one3,two3,one6


one3 = 1.0/3.0
two3 = 2.0*one3
one6 = 1.0/6.0

!     1-point gauss integration

if(l.eq.1) then
        el(1,1) = one3
        el(2,1) = one3
        el(3,1) = one3
        el(4,1) = 1.d0
        lint    = 1

!     3-point integration: mid-edge points

elseif(l.eq.3) then
        el(1,1) = 0.d0
        el(2,1) = 0.5d0
        el(3,1) = 0.5d0
        el(4,1) = one3

        el(1,2) = 0.5d0
        el(2,2) = 0.d0
        el(3,2) = 0.5d0
        el(4,2) = one3

        el(1,3) = 0.5d0
        el(2,3) = 0.5d0
        el(3,3) = 0.d0
        el(4,3) = one3

        lint    = 3

!     3-point integration: interior points

elseif(l.eq.-3) then

        el(1,1) = two3
        el(2,1) = one6
        el(3,1) = one6
        el(4,1) = one3

        el(1,2) = one6
        el(2,2) = two3
        el(3,2) = one6
        el(4,2) = one3

        el(1,3) = one6
        el(2,3) = one6
        el(3,3) = two3
        el(4,3) = one3

        lint    = 3

!     4-point gauss integration: NOT RECOMMENDED DUE TO NEGATIVE WEIGHT

elseif(l.eq.4) then
        el(1,1) =  one3
        el(2,1) =  one3
        el(3,1) =  one3
        el(4,1) = -27.d0/48.d0

        el(1,2) =  0.6d0
        el(2,2) =  0.2d0
        el(3,2) =  0.2d0
        el(4,2) =  25.d0/48.d0

        el(1,3) =  0.2d0
        el(2,3) =  0.6d0
        el(3,3) =  0.2d0
        el(4,3) =  el(4,2)

        el(1,4) =  0.2d0
        el(2,4) =  0.2d0
        el(3,4) =  0.6d0
        el(4,4) =  el(4,2)

        lint    =  4

!     6-point nodal integration

elseif(l.eq.6) then

        el(1,1) =  1.0d0
        el(2,1) =  0.0d0
        el(3,1) =  0.0d0
        el(4,1) =  one6

        el(1,2) =  0.0d0
        el(2,2) =  1.0d0
        el(3,2) =  0.0d0
        el(4,2) =  one6

        el(1,3) =  0.0d0
        el(2,3) =  0.0d0
        el(3,3) =  1.0d0
        el(4,3) =  one6

        el(1,4) =  0.5d0
        el(2,4) =  0.5d0
        el(3,4) =  0.0d0
        el(4,4) =  one6

        el(1,5) =  0.0d0
        el(2,5) =  0.5d0
        el(3,5) =  0.5d0
        el(4,5) =  one6

        el(1,6) =  0.5d0
        el(2,6) =  0.0d0
        el(3,6) =  0.5d0
        el(4,6) =  one6

        lint    =  6

!     6-point order 4 formula

elseif(l.eq.-6) then

        el(1,1) = 0.816847572980459d0
        el(2,1) = 0.091576213509771d0
        el(3,1) = 0.091576213509771d0
        el(4,1) = 0.109951743655322d0

        el(1,2) = 0.091576213509771d0
        el(2,2) = 0.816847572980459d0
        el(3,2) = 0.091576213509771d0
        el(4,2) = 0.109951743655322d0

        el(2,3) = 0.091576213509771d0
        el(1,3) = 0.091576213509771d0
        el(3,3) = 0.816847572980459d0
        el(4,3) = 0.109951743655322d0

        el(1,4) = 0.108103018168070d0
        el(2,4) = 0.445948490915965d0
        el(3,4) = 0.445948490915965d0
        el(4,4) = 0.223381589678011d0

        el(1,5) = 0.445948490915965d0
        el(2,5) = 0.108103018168070d0
        el(3,5) = 0.445948490915965d0
        el(4,5) = 0.223381589678011d0

        el(1,6) = 0.445948490915965d0
        el(2,6) = 0.445948490915965d0
        el(3,6) = 0.108103018168070d0
        el(4,6) = 0.223381589678011d0

        lint    = 6

!     7-point gauss integration

elseif(l.eq.7) then
        r0      =  sqrt(15.0d0)
        r1      =  3.d0/7.d0
        r2      =  (r0 + r0)/21.d0

        el(1,1) =  one3
        el(2,1) =  el(1,1)
        el(3,1) =  el(1,1)
        el(4,1) =  0.225d0

        el(1,2) =  r1 + r2
        el(2,2) =  0.5d0 - 0.5d0*el(1,2)
        el(3,2) =  el(2,2)
        el(4,2) =  (155.d0 - r0)/1200.d0

        el(1,3) =  el(2,2)
        el(2,3) =  el(1,2)
        el(3,3) =  el(2,2)
        el(4,3) =  el(4,2)

        el(1,4) =  el(2,2)
        el(2,4) =  el(2,2)
        el(3,4) =  el(1,2)
        el(4,4) =  el(4,2)

        el(1,5) =  r1 - r2
        el(2,5) =  0.5d0 - 0.5d0*el(1,5)
        el(3,5) =  el(2,5)
        el(4,5) =  (155.d0 + r0)/1200.d0

        el(1,6) =  el(2,5)
        el(2,6) =  el(1,5)
        el(3,6) =  el(2,5)
        el(4,6) =  el(4,5)

        el(1,7) =  el(2,5)
        el(2,7) =  el(2,5)
        el(3,7) =  el(1,5)
        el(4,7) =  el(4,5)

        lint    =  7

!     7-point nodal integration

elseif(l.eq.-7) then

        el(1,1) =  1.0d0
        el(2,1) =  0.0d0
        el(3,1) =  0.0d0
        el(4,1) =  0.05d0

        el(1,2) =  0.0d0
        el(2,2) =  1.0d0
        el(3,2) =  0.0d0
        el(4,2) =  0.05d0

        el(1,3) =  0.0d0
        el(2,3) =  0.0d0
        el(3,3) =  1.0d0
        el(4,3) =  0.05d0

        el(1,4) =  0.5d0
        el(2,4) =  0.5d0
        el(3,4) =  0.0d0
        el(4,4) =  2.0d0/15.0d0

        el(1,5) =  0.0d0
        el(2,5) =  0.5d0
        el(3,5) =  0.5d0
        el(4,5) =  el(4,4)

        el(1,6) =  0.5d0
        el(2,6) =  0.0d0
        el(3,6) =  0.5d0
        el(4,6) =  el(4,4)

        el(1,7) =  one3
        el(2,7) =  one3
        el(3,7) =  one3
        el(4,7) =  0.45d0

        lint    =  7
else
	stop 'tint2d: Number of quadrature points is invalid!'
endif

end subroutine tint2d





subroutine trishp(el,xl,ndm,iord, xsj,shp)

!-----------------------------------------------------------------
!      Purpose: Triangular shape function routine

!       Type:  |iord| =  1:  Linear  three-node
!              |iord| =  2:  Quadratic six-node
!              |iord| =  3:  Quadratic seven-node
!              |iord| =  4:  Quadratic + 3 bubbles (Zienkiewicz/Lefebre)
!              |iord| = 10:  Cubic 10-node

!               iord  > 0:  Mid-side and center node are  global  coords.
!               iord  < 0:  Mid-side and center node heirarchical coords.

!      Inputs:
!         el(3)     - Area coordinates for point
!         xl(ndm,*) - Nodal coordinates for element
!         ndm       - Spatial dimension of mesh
!         iord      - Order of shape functions (see above)

!      Outputs:
!         xsj       - Jacobian determinant at point
!         shp(3,*)  - Shape functions and derivatives at point
!-----------------------------------------------------------------

implicit none

integer :: ndm,iord, i
real :: xsj,xsjr, fel1,fel2,fel3, fel12,fel23,fel31
real :: x1,x2,x3,x4,x5,x6,x7, y1,y2,y3,y4,y5,y6,y7
real :: one3, four9
real :: el(3),xl(ndm,*), shp(3,*)


if(abs(iord) > 9) stop 'trishp: too many iord!'

one3 = 1.0/3.0
four9 = 4.0/9/0

!       Form Jacobian terms

x1 = xl(1,1)
x2 = xl(1,2)
x3 = xl(1,3)

y1 = xl(2,1)
y2 = xl(2,2)
y3 = xl(2,3)

if(abs(iord).gt.1) then

          fel1 = 4.d0*el(1)
          fel2 = 4.d0*el(2)
          fel3 = 4.d0*el(3)

          x4   = xl(1,4)
          x5   = xl(1,5)
          x6   = xl(1,6)

          y4   = xl(2,4)
          y5   = xl(2,5)
          y6   = xl(2,6)

!         Form shape functions in total coordinates

          if(iord.gt.0) then
            x4 = x4 - 0.5d0*(x1 + x2)
            x5 = x5 - 0.5d0*(x2 + x3)
            x6 = x6 - 0.5d0*(x3 + x1)
            x7 = one3*(x1 + x2 + x3) + four9*(x4 + x5 + x6)

            y4 = y4 - 0.5d0*(y1 + y2)
            y5 = y5 - 0.5d0*(y2 + y3)
            y6 = y6 - 0.5d0*(y3 + y1)
            y7 = one3*(y1 + y2 + y3) + four9*(y4 + y5 + y6)
          endif

          x1   = x1 + x4*fel2 + x6*fel3
          x2   = x2 + x5*fel3 + x4*fel1
          x3   = x3 + x6*fel1 + x5*fel2

          y1   = y1 + y4*fel2 + y6*fel3
          y2   = y2 + y5*fel3 + y4*fel1
          y3   = y3 + y6*fel1 + y5*fel2

          fel12 = el(1)*el(2)
          fel23 = el(2)*el(3)
          fel31 = el(3)*el(1)

          if(iord.eq.3) then

            fel12 = 27.d0*fel12
            fel23 = 27.d0*fel23
            fel31 = 27.d0*fel31

            x7    = xl(1,7) - x7
            x1    = x1 + x7*fel23
            x2    = x2 + x7*fel31
            x3    = x3 + x7*fel12

            y7    = xl(2,7) - y7
            y1    = y1 + y7*fel23
            y2    = y2 + y7*fel31
            y3    = y3 + y7*fel12

          endif
endif

xsj  = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
xsjr = 1.0d0
if(xsj.ne.0.0d0) then
  xsjr = 1.0d0/xsj
endif
xsj  = 0.5d0*xsj

!       Specify shape functions and their derivatives

shp(1,1) = (y2-y3)*xsjr
shp(2,1) = (x3-x2)*xsjr
shp(3,1) = el(1)

shp(1,2) = (y3-y1)*xsjr
shp(2,2) = (x1-x3)*xsjr
shp(3,2) = el(2)

shp(1,3) = (y1-y2)*xsjr
shp(2,3) = (x2-x1)*xsjr
shp(3,3) = el(3)


!       Add quadratic hierarchical functions

if(abs(iord).gt.1) then

          shp(1,4) = shp(1,1)*fel2 + shp(1,2)*fel1
          shp(2,4) = shp(2,1)*fel2 + shp(2,2)*fel1
          shp(3,4) = el(1)*fel2

          shp(1,5) = shp(1,2)*fel3 + shp(1,3)*fel2
          shp(2,5) = shp(2,2)*fel3 + shp(2,3)*fel2
          shp(3,5) = el(2)*fel3

          shp(1,6) = shp(1,3)*fel1 + shp(1,1)*fel3
          shp(2,6) = shp(2,3)*fel1 + shp(2,1)*fel3
          shp(3,6) = el(3)*fel1

!         Bubble at baricenter

          if(abs(iord).eq.3) then
            shp(1,7) = shp(1,1)*fel23 + shp(1,2)*fel31 + shp(1,3)*fel12
            shp(2,7) = shp(2,1)*fel23 + shp(2,2)*fel31 + shp(2,3)*fel12
            shp(3,7) = fel12*el(3)
          elseif(abs(iord).eq.4) then
            shp(1,7) = (shp(1,1)*fel23*2.d0+shp(1,2)*fel31+shp(1,3)*fel12)*el(1)
            shp(2,7) = (shp(2,1)*fel23*2.d0+shp(2,2)*fel31+shp(2,3)*fel12)*el(1)
            shp(3,7) = fel12*fel31

            shp(1,8) = (shp(1,1)*fel23+shp(1,2)*fel31*2.d0+shp(1,3)*fel12)*el(2)
            shp(2,8) = (shp(2,1)*fel23+shp(2,2)*fel31*2.d0+shp(2,3)*fel12)*el(2)
            shp(3,8) = fel12*fel23

            shp(1,9) = (shp(1,1)*fel23+shp(1,2)*fel31+shp(1,3)*fel12*2.d0)*el(3)
            shp(2,9) = (shp(2,1)*fel23+shp(2,2)*fel31+shp(2,3)*fel12*2.d0)*el(3)
            shp(3,9) = fel31*fel23
          endif

!         Modify shape functions for mid-side and interior values

          if(iord.gt.1) then

!           Modify vertex and mid-side values for bubble

            if(iord.eq.3) then
              do i = 1,3
                shp(i,1) = shp(i,1) -  one3*shp(i,7)
                shp(i,2) = shp(i,2) -  one3*shp(i,7)
                shp(i,3) = shp(i,3) -  one3*shp(i,7)
                shp(i,4) = shp(i,4) - four9*shp(i,7)
                shp(i,5) = shp(i,5) - four9*shp(i,7)
                shp(i,6) = shp(i,6) - four9*shp(i,7)
              end do ! i
            endif

!           Modify vertex shape functions for mid-side values

            do i = 1,3
              shp(i,1) = shp(i,1) - 0.5d0*(shp(i,4) + shp(i,6))
              shp(i,2) = shp(i,2) - 0.5d0*(shp(i,5) + shp(i,4))
              shp(i,3) = shp(i,3) - 0.5d0*(shp(i,6) + shp(i,5))
            end do ! i

          endif
          
endif


end subroutine trishp





subroutine bbar2s(phi,shp,dvol,lint,nel,npm,hh,irad,theta,bbar)

!-----------------------------------------------------------------
!     Purpose: Compute mixed formulation for the volumetric response.
!              4-node and 9-node case for small deformation problem.

!     Inputs:
!        phi(6,*)      - Stress functions
!        shp(3,9,*)    - Shape functions and derivatives
!        vol(*)        - Volume elements
!        lint          - Number of quadrature points
!        nel           - Number of nodes on element (should be 4 or 9)
!        npm           - Number of pressure modes
!        irad(*)       - Inverse radius (or zero) at quadrature points
!        theta(3,*)    - Volumetric strain from displacements

!     Outputs:
!        hh(3,3)       - Volume/pressure shape integrals (inverse)
!        theta(3,*)    - Mixed volumetric strain
!        bbar(2,16,*)  - Mixed volumetric derivatives of shape functions.
!                        (N.B. Includes axisymmetric part using irad(*))
!-----------------------------------------------------------------
  
implicit none

integer :: lint,   nel,  npm,   i,  j,  k, l
real :: shp(3,16,*),  dvol(*)   ,  ht(6,2),   h1
real :: gg(6,2,16) ,  hh(6,6)   ,  hv(6,2)   ,  hg(6,2,16)
real :: irad(*)    ,  phi(6,25) ,  theta(3,*),  bbar(2,16,*)


!     Constant pressure elements

if (npm == 1) then

        do j = 1,nel
          bbar(1,j,1) = 0.0d0
          bbar(2,j,1) = 0.0d0
        end do ! j
        hh(1,1) = 0.d0
        ht(1,1) = 0.0d0
        ht(1,2) = 0.0d0

        do l = 1,lint

!         H-array and G-array

          hh(1,1) = hh(1,1) + dvol(l)
          ht(1,1) = ht(1,1) + theta(1,l)*dvol(l)
          ht(1,2) = ht(1,2) + theta(2,l)*dvol(l)

!         G-array

          do j = 1,nel
            bbar(1,j,1) = bbar(1,j,1) + (shp(3,j,l) * irad(l) + shp(1,j,l))* dvol(l)
            bbar(2,j,1) = bbar(2,j,1) +  shp(2,j,l) * dvol(l)
          end do ! j
        end do ! l

!       Average Jacobian

        hh(1,1)    = 1.d0 / hh(1,1)

!       Small deformation case

        theta(1,1) = hh(1,1)*ht(1,1)
        theta(2,1) = hh(1,1)*ht(1,2)
        theta(3,1) = theta(1,1) - theta(2,1)

!       Modify bbar for B-bar type computations

        do j = 1,nel
          do i = 1,2
            bbar(i,j,1) = bbar(i,j,1)*hh(1,1)
          end do ! i
        end do ! j

!       Copy for other quadrature points

        do l = 2,lint
          theta(1,l) = theta(1,1)
          theta(2,l) = theta(2,1)
          theta(3,l) = theta(3,1)
          do j = 1,nel
            bbar(1,j,l) = bbar(1,j,1)
            bbar(2,j,l) = bbar(2,j,1)
          end do ! j
        end do ! l


!     Higher order elements: npm = 3 (8 or 9 nodes); npm = 6 (16 nodes)
else

        do i = 1,npm
          do j = 1,nel
            gg(i,1,j) = 0.0d0
            gg(i,2,j) = 0.0d0
          end do ! j
          do j = 1,npm
            hh(i,j) = 0.0d0
          end do ! j
          ht(i,1) = 0.0d0
          ht(i,2) = 0.0d0
        end do ! i

!       Quadrature loop

        do l = 1,lint
          do j = 1,npm

!           H-array

            h1      = phi(j,l) * dvol(l)
            ht(j,1) = ht(j,1)  + theta(1,l)*h1
            ht(j,2) = ht(j,2)  + theta(2,l)*h1
            do i = 1,npm
              hh(i,j)   = hh(i,j)   + phi(i,l)*h1
            end do ! i

!           G-array

            do i = 1,nel
              gg(j,1,i) = gg(j,1,i) + (shp(1,i,l) + shp(3,i,l)*irad(l))*h1
              gg(j,2,i) = gg(j,2,i) +  shp(2,i,l)*h1
            end do ! i
          end do ! j

        end do ! l

!       Invert H-array

        call invert(hh,npm,6)

        do j = 1,nel
          do i = 1,npm
            hg(i,1,j) = 0.0d0
            hg(i,2,j) = 0.0d0
            do k = 1,npm
              hg(i,1,j) = hg(i,1,j) + hh(i,k)*gg(k,1,j)
              hg(i,2,j) = hg(i,2,j) + hh(i,k)*gg(k,2,j)
            end do ! k
          end do ! i
        end do ! j

        do j = 1,2
          do i = 1,npm
            hv(i,j) = 0.0d0
            do k = 1,npm
              hv(i,j) = hv(i,j) + hh(i,k)*ht(k,j)
            end do ! k
          end do ! i
        end do ! j

        do l = 1,lint
          theta(1,l) = hv(1,1)
          theta(2,l) = hv(1,2)
          do k = 2,npm
            theta(1,l) = theta(1,l) + phi(k,l)*hv(k,1)
            theta(2,l) = theta(2,l) + phi(k,l)*hv(k,2)
          end do ! k
        end do ! l

        do l = 1,lint
          theta(3,l) = theta(1,l) - theta(2,l)
          do j = 1,nel
            bbar(1,j,l) = hg(1,1,j)
            bbar(2,j,l) = hg(1,2,j)
            do k = 2,npm
              bbar(1,j,l) = bbar(1,j,l) + phi(k,l)*hg(k,1,j)
              bbar(2,j,l) = bbar(2,j,l) + phi(k,l)*hg(k,2,j)
            end do ! k
          end do ! j
        end do ! l

endif

end subroutine bbar2s





subroutine dmatmx (aa,dd)

!-----------------------------------------------------------------
!     Purpose: Project 6x6 AA-matrix onto a 7x7 DD-matrix for mixed model
!              (N.B. Jacobian ratio is in dvol)
!                     | D_11   D_12 |
!                DD = |             |
!                     | D_21   D_22 |
!              where:
!                D_11  = 6 x 6 Deviatoric part of matrix
!                D_12  = 6 x 1 Coupling   part of matrix
!                D_21  = 1 x 6 Coupling   part of matrix
!                D_22  = 1 x 1 Volumetric part of matrix

!     Inputs:
!        aa(6,6)   - Material tangent matrix (based on F)

!     Outputs:
!        dd(7,7)   - Mixed material tangent for stiffness computations.
!-----------------------------------------------------------------

implicit none

integer :: i, j
real :: aa(6,6), dd(7,7)


!     Load moduli from constitution

      do i = 1,6
        do j = 1,6
          dd(i,j) = aa(i,j)
        end do ! j
      end do ! i

!     Compute left and right multiples with trace

      do i = 1,6
        dd(i,7) = (aa(i,1) + aa(i,2) + aa(i,3))/3.0
        dd(7,i) = (aa(1,i) + aa(2,i) + aa(3,i))/3.0
      end do ! i

!     Convert upper 6 x 6 to a deviatoric D_11

      do i = 1,6
        do j = 1,3
          dd(i,j) = dd(i,j) - dd(i,7)
          dd(j,i) = dd(j,i) - dd(7,i)
        end do ! j
      end do ! i

!     Form last term, D_22

      dd(7,7) = (dd(1,7) + dd(2,7) + dd(3,7))/3.0

!     Final update to form D_12 and D_21

      do i = 1,3
        dd(i,7) = dd(i,7) - dd(7,7)
        dd(7,i) = dd(7,i) - dd(7,7)
        do j = 1,3
          dd(i,j) = dd(i,j) + dd(7,7)
        end do ! j
      end do ! i

end subroutine dmatmx




subroutine strn2m(shp,xl,ul,theta,irad,ndm,ndf,nel,eps)

!-----------------------------------------------------------------
!      Purpose: Compute mixed strain for near incompressible formulation.

!      Inputs:
!        shp(3,nen,*)  = Shape functions
!        xl(ndm,nen)   = Nodal coordinates
!        ul(ndf,nen,*) = Nodal solution parameters
!        theta         = Volume change (mixed form)
!        irad          = Inverse radius (or zero for plane).
!        ndm           = Spatial dimension of mesh
!        ndf           = DOF/node (max)
!        nel           = Number nodes on element (4 or 9)

!      Outputs:
!        eps(*)        = Mixed strain at point
!-----------------------------------------------------------------

implicit none

integer :: ndm,ndf,nel,k
real :: irad,theta,dtheta
real :: shp(3,*), xl(ndm,*), ul(ndf,*), eps(*)


!     Compute strain tensor for constitutive equations

      do k = 1,6
        eps(k) = 0.0d0
      end do ! k


      do k = 1,nel
        eps(1)  = eps(1)  + shp(1,k)*ul(1,k)
        eps(2)  = eps(2)  + shp(2,k)*ul(2,k)
        eps(3)  = eps(3)  + shp(3,k)*ul(1,k)
        eps(4)  = eps(4)  + shp(1,k)*ul(2,k) + shp(2,k)*ul(1,k)
      end do ! k

!     Modification for plane/axisymmetry

      eps(3) = eps(3)*irad
!     Correct strains and incremental strains for mixed formulation

      dtheta = (theta - eps(1) - eps(2) - eps(3))/3.0
      eps(1) = eps(1) + dtheta
      eps(2) = eps(2) + dtheta
      eps(3) = eps(3) + dtheta
      
!     write(*,*) 'eps(1:2)&eps(3)-hoop strain: ',eps(1),eps(2), eps(3)

end subroutine strn2m





subroutine char_length(xl, chleng)

! 2D - 4 nodes elements only

real, intent(in) :: xl(2,4)
real, intent(out) :: chleng
real :: lx, ly

lx = 0.5e0*(sqrt((xl(1,1)-xl(1,4))**2+(xl(2,1)-xl(2,4))**2) &
	+sqrt((xl(1,2)-xl(1,3))**2+(xl(2,2)-xl(2,3))**2))
ly = 0.5e0*(sqrt((xl(1,1)-xl(1,2))**2+(xl(2,1)-xl(2,2))**2) &
	+sqrt((xl(1,4)-xl(1,3))**2+(xl(2,4)-xl(2,3))**2))
	
chleng = 0.5e0*lx + 0.5e0*ly

end subroutine char_length





subroutine invert(a,nmax,ndm)

!-----------------------------------------------------------------
!      Purpose: Invert small square matrix

!      Inputs:
!         a(ndm,*) - Matrix to be inverted
!         nmax     - Size of upper submatrix to invert
!         ndm      - Dimension of array

!      Outputs:
!         a(ndm,*) - Submatrix replaces original terms, others not
!                    changed
!-----------------------------------------------------------------

implicit none

integer :: i,j,n,ndm,nmax
real :: d, a(ndm,*)


      do n = 1,nmax
        if(a(n,n).ne.0.0d0) then
          d = 1.d0/a(n,n)
          do j = 1,nmax
            a(n,j) = -a(n,j)*d
          end do ! j

          do i = 1,nmax
            if(n.ne.i) then
              do j = 1,nmax
                if(n.ne.j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
              end do ! j
            endif
            a(i,n) = a(i,n)*d
          end do ! i
          a(n,n) = d
        else
          write(*,*) ' *WARNING* Zero pivot in INVERT'
        endif
      end do ! n

end subroutine invert