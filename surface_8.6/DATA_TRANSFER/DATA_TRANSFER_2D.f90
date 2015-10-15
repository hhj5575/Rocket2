MODULE DATA_TRANSFER_MOD_2D

IMPLICIT NONE; PRIVATE
PUBLIC INTERPOLATION

CONTAINS

SUBROUTINE INTERPOLATION(snNode,tnNode,sx,sy,sz,tx,ty,tz,interpolation_scheme)

IMPLICIT NONE
INTEGER :: snNode, tnNode,i,j
INTEGER :: cnnode, qnnode,index
real*8 sx(snNode), sy(snNode), sz(snNode)
real*8 tx(tnNode), ty(tnNode), tz(tnNode)
real*8, dimension(:) , allocatable :: cx,cy,cx1,cy1,cx2,cy2,cz
real*8, dimension(:) , allocatable :: qx, qy, load


integer shape_order, gauss_order
integer differ
integer temp_num
integer interpolation_scheme
integer poly_order

real*8, dimension(:,:), allocatable :: mass
real*8 alpha, beta,w_sobolev
integer, dimension(:), allocatable :: ipoint,jpoint
real*8 tolerance
integer buffer


cnnode = snNode + tnNode
allocate(cx(cnnode), cy(cnnode))
allocate(cx1(cnnode),cy1(cnnode),cx2(cnnode),cy2(cnnode),cz(cnnode))
allocate(ipoint(cnnode),jpoint(cnnode))
allocate(load(tnNode))

! Initialization
cx  = 0.0d0 ; cy  = 0.0d0 ; cz  = 0.0d0
cx1 = 0.0d0 ; cy1 = 0.0d0 ; cx2 = 0.0d0; cy2 = 0.0d0
ipoint = 0  ; jpoint = 0  ; load = 0.0d0;

buffer = tnNode

! Coefficient
alpha = 0.5d0
beta = 1.0d0-alpha
w_sobolev = 0.0d0

gauss_order = 9
shape_order = 1
poly_order = 1

tolerance = 1e-12

differ = 0
if(snNode.eq.tnNode) then
  do i = 1, snNode
      if(sx(i).NE.tx(i).and.sy(i).NE.ty(i)) differ = differ+1
  end do
end if

if(differ.ne.0.and.snNode.eq.tnNode) then

  cnnode = snNode
  do i = 1, cnnode
      cx(i) = sx(i)
      cy(i) = sy(i)
  end do

else

  ! Source mesh leading, trailing / Target mesh leading, trailing
  call reference_interface(snNode,tnNode,cnnode,sx,sy,tx,ty,cx1,cy1,ipoint,temp_num)
  call reference_interface(tnNode,snNode,cnnode,tx,ty,sx,sy,cx2,cy2,jpoint,temp_num)

  cnnode = temp_num

  if(interpolation_scheme.eq.1) then

    do i=1, cnnode
      cx(i) = alpha*cx1(i) + beta*cx2(i)
      cy(i) = alpha*cy1(i) + beta*cy2(i)
    end do

  else

    do i=1, cnnode
          cx(i) = cx1(i)
          cy(i) = cy1(i)
    end do
    goto 100

  end if

end if

do i = 1, cnnode-1
  if(i.eq.cnnode) exit
    if((dabs(cx(i+1)-cx(i)).LT.tolerance).and.(dabs(cy(i+1)-cy(i)).LT.tolerance)) then
      do j = i+1, cnnode-1
        cx(j) = cx(j+1)
        cy(j) = cy(j+1)
        ipoint(j) = ipoint(j+1)
        jpoint(j) = jpoint(j+1)
      end do
    cnnode = cnnode-1
    ipoint(i) = 1
    jpoint(i) = 1
  end if
end do

ipoint(1) = 1
jpoint(1) = 1
ipoint(cnnode) = 1
jpoint(cnnode) = 1

qnnode = (cnnode-1)*gauss_order
allocate(qx(qnnode), qy(qnnode))
qx = 0.d0 ; qy = 0.0d0

call quadrature_point(cnnode,qnnode,cx,cy,qx,qy,gauss_order)

allocate(mass(tnNode,shape_order+1))
mass = 0.0d0

do i = 1, tnNode
  load(i) = 0.0d0
  do j = 1, shape_order+1
    mass(i,j) = 0.0d0
  end do
end do

call Sobolev_minimization_Matrix(gauss_order,shape_order,snNode,cnnode,tnNode,qnnode,&
                                sx,sy,sz,tx,ty,cx,cy,qx,qy,ipoint,jpoint,load,mass,w_sobolev)

call band_solver(tnNode,shape_order+1,tz,mass,load)

tz(1) = sz(1)
tz(tnNode) = sz(snNode)

deallocate(qx, qy, mass)


100 continue

if(interpolation_scheme.eq.2) then
  call Polynomial_interpolation(poly_order,snNode,tnNode,cnnode,ipoint,jpoint,cx,cy,sz,tz,index)
else if(interpolation_scheme.eq.3) then
  call TPS(snNode,tnNode,cnnode,ipoint,jpoint,cx,cy,sz,tz,index)
end if

tnNode = buffer

tz(1) = sz(1)
tz(tnNode) = sz(snNode)

deallocate(cx, cy, cx1, cy1, cx2, cy2, cz, ipoint, jpoint, load)

END SUBROUTINE INTERPOLATION

subroutine reference_interface(snNode,tnNode,cnnode,sx,sy,tx,ty,cx,cy,point,temp_num)

implicit none
integer snNode, tnNode, cnnode,i,j

double precision sx(snNode), sy(snNode)
double precision tx(tnNode), ty(tnNode)
double precision cx(cnnode), cy(cnnode)
integer point(cnnode), temp_num

double precision sx1, sx2, sx3, sy1, sy2, sy3
double precision sxvec1, syvec1, sdis
double precision sxvec2, syvec2
double precision sx_unit1, sy_unit1, sx_unit2, sy_unit2
double precision sx_norm, sy_norm

double precision tx1, tx2, ty1, ty2
double precision txvec, tyvec, tdis
double precision tx_unit, ty_unit

double precision a1, a2, b1, b2
double precision xcoord, ycoord, eps, eps2, xmin, ymin, xmax, ymax
double precision dist

integer temp_index
integer index, previous

eps = 1e-16
eps2 = 1e-7
previous   = 1
temp_index = 1

index = 1
cx(1) = tx(1)
cy(1) = ty(1)

do i = 1, tnNode-1

  tx1 = tx(i)
  ty1 = ty(i)

  tx2 = tx(i+1)
  ty2 = ty(i+1)

  txvec = tx2-tx1
  tyvec = ty2-ty1
  tdis = dsqrt((txvec**2)+(tyvec**2))

  tx_unit = txvec/tdis
  ty_unit = tyvec/tdis

  a1 = ty_unit/(tx_unit+eps)
  b1 = ty1 - a1*tx1

  xmin = dmin1(tx1, tx2)
  xmax = dmax1(tx1, tx2)

  ymin = dmin1(ty1, ty2)
  ymax = dmax1(ty1, ty2)


  do j = previous+1, snNode-1

    sx1 = sx(j-1)
    sy1 = sy(j-1)

    sx2 = sx(j)
    sy2 = sy(j)

    sx3 = sx(j+1)
    sy3 = sy(j+1)

    sxvec1 = sx2-sx1
    syvec1 = sy2-sy1
    sdis = dsqrt(sxvec1**2 + syvec1**2)
    sx_unit1 = sxvec1/sdis
    sy_unit1 = syvec1/sdis

    sxvec2 = sx3-sx2
    syvec2 = sy3-sy2
    sdis = dsqrt(sxvec2**2 + syvec2**2)
    sx_unit2 = sxvec2/sdis
    sy_unit2 = syvec2/sdis

    sx_norm = -0.5d0*(sy_unit1+sy_unit2)
    sy_norm =  0.5d0*(sx_unit1+sx_unit2)

    a2 = sy_norm/(sx_norm+eps)
    b2 = sy2 - a2*sx2


    if(a1.eq.a2) then
      goto 10
    else
      xcoord = (b2-b1)/(a1-a2)
      ycoord = a1*xcoord + b1
    end if

    if(dabs(sy_norm-1.0d0).LE.eps.or.dabs(sy_norm+1.0d0).LE.eps) then
      xcoord = sx2
      ycoord = a1*xcoord + b1
    end if

    if(dabs(sx_norm-1.0d0).LE.eps.or.dabs(sx_norm+1.0d0).LE.eps) then
      ycoord = sy2
      xcoord = (ycoord - b1)/(a1+eps)
    end if

    if(dabs(tx_unit-1.0d0).LE.eps.or.dabs(tx_unit+1.0d0).LE.eps) then
      ycoord = ty2
      xcoord = (ycoord - b2)/(a2+eps)
    end if

    if(dabs(ty_unit-1.0d0).LE.eps.or.dabs(ty_unit+1.0d0).LE.eps) then
      xcoord = tx2
      ycoord = a2*xcoord + b2
    end if

    if(xcoord.GT.xmax+eps2.or.xcoord.LT.xmin-eps2) goto 10
    if(ycoord.GT.ymax+eps2.or.ycoord.LT.ymin-eps2) goto 10

    dist = dsqrt((sx2-xcoord)**2 + (sy2-ycoord)**2)
    ! larger then maximum distance between source and target
    if(dist.GT.0.001d0) goto 10

    temp_index = j

    index = index+1
    cx(index) = xcoord
    cy(index) = ycoord
    point(index) = 1


10  continue

  end do

  index = index+1
  cx(index) = tx2
  cy(index) = ty2

  previous = temp_index

end do

temp_num = index

end subroutine reference_interface

Subroutine quadrature_point(cnnode,qnnode,cx,cy,qx,qy,kind)

implicit none
integer cnnode, qnnode,i,j
real*8 cx(cnnode), cy(cnnode)
real*8 qx(qnnode), qy(qnnode)
integer kind

    integer index
real(8) cx1, cx2, cy1, cy2, cdis
real*8 cxvec, cyvec
real*8 cx_unit, cy_unit
real*8 xx
    real*8 xi(9,9), w(9,9)

call gauss_point(xi,w)

index = 0

do i = 1, cnnode-1

  cx1 = cx(i)
  cy1 = cy(i)
  cx2 = cx(i+1)
  cy2 = cy(i+1)

  cxvec = cx2-cx1
  cyvec = cy2-cy1

  cdis = dsqrt((cxvec**2) + (cyvec**2))

  cx_unit = cxvec/cdis
  cy_unit = cyvec/cdis

  do j = 1, kind
    xx = 0.5d0*(1+xi(j,kind))*cdis
    index = index+1
    qx(index) = cx1 + xx*cx_unit
    qy(index) = cy1 + xx*cy_unit
  end do

end do

End Subroutine quadrature_point

subroutine sobolev_minimization_matrix(kind,order,snNode,cnnode,tnNode,qnnode,&
           sx,sy,sz,tx,ty,cx,cy,qx,qy,ipoint,jpoint,load,global_mass,w_sobolev)

integer kind, order, snNode, cnnode, tnNode, qnnode
integer ipoint(cnnode),jpoint(cnnode)
real*8 sx(snNode), sy(snNode), sz(snNode)
real*8 tx(tnNode), ty(tnNode)
real*8 qx(qnnode),qy(qnnode)
real*8 cx(cnnode), cy(cnnode)
real*8 load(tnNode)
real*8 global_mass(tnNode,order+1)
real*8 w_sobolev


real*8 xi(9,9), w(9,9)
real*8 psi(order+1), dpsi(order+1)
real*8 spsi(order+1), sdpsi(order+1)
real*8 tpsi(order+1), tdpsi(order+1)
integer snode(order+1), tnode(order+1)
real*8 temp(order+1)
real*8 mass(order+1,order+1)

integer quad_index, sele_num, tele_num

integer sflag, source_index
integer tflag, target_index

real*8 sx1, sx2, sy1, sy2
real*8 sxvec, syvec, sdis
real*8 sx_unit, sy_unit

real*8 tx1, tx2, ty1, ty2
real*8 txvec, tyvec, tdis
real*8 tx_unit, ty_unit

real*8 sqx1, qx2, sqy1, qy2, tqx1, tqy1

real*8 qxvec, qyvec
real*8 innerp

real*8 tempx, tempy
real*8 app_val, dapp_val
real*8 xxi
integer i, j, k, l , m
reaL(8) :: c_delta

call Gauss_point(xi,w)

quad_index = 0
sele_num = 1
tele_num = 1
source_index = 0
target_index = 0

do i = 1, cnnode-1

  if(ipoint(i).eq.1) source_index = source_index+1
  if(jpoint(i).eq.1) target_index = target_index+1

  do j = 1, order+1
    snode(j) = (order*(sele_num-1))+j
    tnode(j) = (order*(tele_num-1))+j
  end do

  do j = 1, order+1
      temp(j) = 0.0d0
  end do

  do j = 1, order+1
    do k = 1, order+1
      mass(j,k) = 0.0d0
    end do
  end do

  c_delta = dsqrt((cx(i+1)-cx(i))**2 + (cy(i+1)-cy(i))**2)

  sx1 = sx(source_index)
  sy1 = sy(source_index)

  if (source_index.eq.snNode) then
    sx1 = sx(source_index-1)
    sy1 = sy(source_index-1)
    sx2 = sx(source_index)
    sy2 = sy(source_index)
  else
    sx2 = sx(source_index+1)
    sy2 = sy(source_index+1)
  endif

  sxvec = sx2-sx1
  syvec = sy2-sy1
  sdis = dsqrt((sxvec**2) + (syvec**2))

  sx_unit = sxvec / sdis
  sy_unit = syvec / sdis

  tx1 = tx(target_index)
  ty1 = ty(target_index)

  if (target_index.eq.tnNode) then
    tx1 = tx(target_index-1)
    ty1 = ty(target_index-1)
    tx2 = tx(target_index)
    ty2 = ty(target_index)
  else
    tx2 = tx(target_index+1)
    ty2 = ty(target_index+1)
  endif

  txvec = tx2-tx1
  tyvec = ty2-ty1
  tdis = dsqrt((txvec**2) + (tyvec**2))

  tx_unit = txvec / tdis
  ty_unit = tyvec / tdis

  sqx1 = sx1
  sqy1 = sy1

  tqx1 = tx1
  tqy1 = ty1

  do j = 1, kind

    quad_index = quad_index+1

    ! source mesh·ÎÀÇ orthogonal projection
    qx2 = qx(quad_index)
    qy2 = qy(quad_index)

    qxvec = qx2-sqx1
    qyvec = qy2-sqy1

    innerp = (sx_unit*qxvec) + (sy_unit*qyvec)

    tempx = sx1 + innerp*sx_unit
    tempy = sy1 + innerp*sy_unit

    ! source mesh·ÎÀÇ orthogonal projection ÁŸ·á
    ! source mesh¿¡Œ­ÀÇ approximation °ª
    call coordinate_transformation(sele_num,order,snNode,source_index,sx,sy,tempx,tempy,xxi)
    call Shape_function(xxi,order,spsi,sdpsi)
    call legendre_interpolation(1,source_index,snNode,sz,xxi,app_val,dapp_val)

    ! source mesh¿¡Œ­ÀÇ approximation °ª ±žÇÏ±â ÁŸ·á
    ! target mesh·ÎÀÇ orthogonal projection

    qxvec = qx2-tqx1
    qyvec = qy2-tqy1

    innerp = (tx_unit*qxvec) + (ty_unit*qyvec)

    tempx = tx1 + innerp*tx_unit
    tempy = ty1 + innerp*ty_unit

    ! target mesh·ÎÀÇ orthogonal projection ÁŸ·á
    call coordinate_transformation(tele_num,order,tnNode,target_index,tx,ty,tempx,tempy,xxi)
    call Shape_function(xxi,order,tpsi,tdpsi)
    call shape_function(xi(j,kind),order,psi,dpsi)

    do l = 1, order+1
      do m = 1, order+1
        mass(l,m) = mass(l,m)+(tpsi(l)*tpsi(m)+w_sobolev*tdpsi(l)*tdpsi(m))*c_delta*w(j,kind)
      end do
    end do

    do k = 1, order+1
      temp(k) = temp(k)+((tpsi(k)*app_val)+(tdpsi(k)*dapp_val*w_sobolev))*c_delta*w(j,kind)
    end do

  end do

  call assembly(tnode(1),tnNode,order,mass,global_mass)

  do j = 1, order+1
    load(tnode(j)) = load(tnode(j)) + temp(j)
  end do

  sflag = mod(source_index,order)
  tflag = mod(target_index,order)

  if((sflag.eq.0).and.(ipoint(i+1).eq.1)) sele_num=sele_num+1
  if((tflag.eq.0).and.(jpoint(i+1).eq.1)) tele_num=tele_num+1

end do

end subroutine sobolev_minimization_matrix

Subroutine TPS(snNode,tnNode,cnnode,ipoint,jpoint,cx,cy,sz,tz,index)

implicit none
integer snNode, tnNode, cnnode,index,i,j
real*8 cx(cnnode), cy(cnnode), sz(snNode),tz(tnNode+2)
integer ipoint(cnnode), jpoint(cnnode)

integer sindex(snNode), tindex(tnNode)
integer temp_index1, temp_index2

real*8 A(snNode+3,snNode+3), B(snNode+3)
real*8 results(snNode+3)
real*8 eps
real*8 dist

eps = 1e-5

temp_index1 = 0
temp_index2 = 0
do i = 1, cnnode
  if(ipoint(i).eq.1) then
    temp_index1 = temp_index1 + 1
    sindex(temp_index1) = i
  end if

  if(jpoint(i).eq.1) then
    temp_index2 = temp_index2 + 1
    tindex(temp_index2) = i
  end if
  continue
end do

sindex(1) = 1
sindex(snNode) = cnnode

if(index.eq.1) then
  do i = 1, snNode
    sindex(i) = i
  end do
  do i = 1, tnNode
    tindex(i) = i+1
  end do
end if

do i  = 1, snNode+3
  do j = 1, snNode+3
    A(i,j) = 0.0d0
  end do
end do



do i = 1, snNode

  do j = i, snNode
    dist = (cx(sindex(i))-cx(sindex(j)))**2 + (cy(sindex(i))-cy(sindex(j)))**2
    A(i,j-i+1) = dist*dlog(dist+eps)
    continue
  end do

  A(i,1) = eps
  A(i,snNode+2-i) = 1.0d0
  A(i,snNode+3-i) = cx(sindex(i))
  A(i,snNode+4-i) = cy(sindex(i))

  B(i) = sz(i)

end do

B(snNode+1) = 0.0d0
B(snNode+2) = 0.0d0
B(snNode+3) = 0.0d0

call Band_Solver(snNode+3,snNode+3,results,A,B)

call TPS_interpolation(snNode+3,tnNode+2,cnnode,snNode,sindex,tindex,cx,cy,results,tz)

end subroutine TPS

Subroutine TPS_interpolation(Fnum,tnNode,cnnode,snNode,sindex,tindex,cx,cy,F,tz)

implicit none
integer Fnum, tnNode,cnnode,snNode,i,j
integer sindex(snNode), tindex(tnNode-2)
real*8 cx(cnnode), cy(cnnode),eps,temp
real*8 F(Fnum), tz(tnNode)
real*8 R

eps = 1d-5

continue

do j = 1, tnNode-2

  temp = 0.d0

  do i = 1, snNode
    R = (cx(tindex(j))-cx(sindex(i)))**2+(cy(tindex(j))-cy(sindex(i)))**2
    temp = temp + F(i)*R*dlog(R+eps)
  end do

  tz(j) = F(snNode+1) + (F(snNode+2)*cx(tindex(j))) + (F(snNode+3)*cy(tindex(j))) + temp

end do

end subroutine TPS_interpolation

Subroutine Shape_function(xi,m,psi,dpsi)

!	Calculate the Legendre polynomial shape function
!
!	xi : coordinate at the transformed coordinate
!	m : Order of legendre polynomial function
!	psi : calculation results of function
!	dpsi : calcuation results of derivative of function

implicit none
real*8 xi
integer m,j
real*8 psi(10), dpsi(10)
real*8, dimension(:), allocatable :: xxi
real*8 delta_xi
real*8, dimension(:), allocatable :: num,dinum
real*8, dimension(:), allocatable :: d_num
integer i

allocate(xxi(m+1))
allocate(d_num(m+1), num(m+1), dinum(m+1))

delta_xi = 2.0/m

xxi(1) = -1
xxi(m+1) = 1

do i = 2, m
  xxi(i) = xxi(i-1)+delta_xi
end do

do i = 1, m+1
  num(i) = 1
  dinum(i) = 1
  d_num(i) = 0
end do

do i = 1, m+1
  do j = 1, m+1
    if(i.NE.j) then
      num(i) = num(i)*(xi-xxi(j))
      dinum(i) = dinum(i)*(xxi(i)-xxi(j))
    end if
  end do
  do j =1, m+1
    if((i.NE.j).and.xi.NE.xxi(j)) then
    d_num(i) = d_num(i) + num(i)/(xi-xxi(j))
    end if
  end do
  psi(i) = num(i)/dinum(i)
  dpsi(i) = d_num(i)/dinum(i)
end do

deallocate(xxi)
deallocate(d_num, num, dinum)

end subroutine Shape_function

! Execuate the legendre interpolation via source x and y
Subroutine Legendre_interpolation(order,ele_num,nnode,z,xxi,results,deri)

implicit none
integer order, ele_num, nnode,l
real*8 z(nnode)
real*8 xxi
real*8 results,deri

real*8 psi(order+1), dpsi(order+1)
integer node(order+1)

do l = 1, order+1
    node(l) = (order*(ele_num-1))+l
end do

call Shape_function(xxi,order,psi,dpsi)

results = 0.0d0
deri = 0.0d0
do l = 1, order+1
  results = results + psi(l)*z(node(l))
  deri = deri + dpsi(l)*z(node(l))
end do

end Subroutine Legendre_interpolation

Subroutine Polynomial_interpolation(order,snNode,tnNode,cnnode,ipoint,jpoint,cx,cy,sz,tz,index)

implicit none
integer order, snNode, tnNode, cnnode,index,i,j
real*8 cx(cnnode), cy(cnnode), sz(snNode),tz(tnNode)
integer ipoint(cnnode), jpoint(cnnode)

integer sindex(snNode), tindex(tnNode)
integer temp_index1, temp_index2
integer total_ele
integer node(order+1)
integer object(tnNode), object_num, tnNode2
real*8 ele_z(order+1)
real*8 xxi
real*8 temp3

temp_index1 = 1
temp_index2 = 0
do i = 2, cnnode-1
  if(ipoint(i).eq.1) then
    temp_index1 = temp_index1 + 1
    sindex(temp_index1) = i
  end if
  if(jpoint(i).eq.1) then
    temp_index2 = temp_index2 + 1
    tindex(temp_index2) = i
  end if
  continue
end do
tnNode2 = temp_index2
sindex(1) = 1
sindex(snNode) = cnnode

if(index.eq.1) then
  do i = 1, snNode
    sindex(i) = i
  end do
  do i = 1, tnNode
    tindex(i) = i+1
  end do
end if

total_ele = (snNode-1)/order

temp_index1 = 2
do i = 1, total_ele

  do j = 1, order+1
    node(j) = sindex((order*(i-1))+j)
    ele_z(j) = sz((order*(i-1))+j)
  end do

  object_num = 0
  do j = 1, tnNode2
    if(tindex(j).GE.node(1).and.tindex(j).LE.node(order+1)) then
      object_num = object_num+1
      object(object_num) = tindex(j)
    end if
  end do

  do j = 1, object_num
    call coordinate_transformation2(node,order,cnnode,object(j),cx,cy,cx(object(j)),cy(object(j)),xxi)
    call Legendre_interpolation(order,i,snNode,sz,xxi,tz(temp_index1),temp3)
    temp_index1 = temp_index1+1
    continue
  end do

  continue
  do j = 1, tnNode2-object_num
    tindex(j) = tindex(j+object_num)
  end do
  tnNode2 = tnNode2 - object_num

end do

End subroutine Polynomial_interpolation

Subroutine Gauss_point(xi,w)

implicit none
real*8 xi(9,9), w(9,9)

!	Gauss point 1
xi(1,1)=0.D0
w(1,1)=2.D0

!	Gauss point 2
xi(1,2)=-1.D0/dsqrt(3.D0)
xi(2,2)=-xi(1,2)
w(1,2)=1.D0
w(2,2)=w(1,2)

!	Gauss point 3
xi(1,3)=-dsqrt(3.D0/5.D0)
xi(2,3)=0.D0
xi(3,3)=-xi(1,3)
w(1,3)=5.D0/9.D0
w(2,3)=8.D0/9.D0
w(3,3)=w(1,3)

!	Gauss point 4
xi(1,4)=-0.861136311594053
xi(2,4)=-0.339981043584856
xi(3,4)=-xi(2,4)
xi(4,4)=-xi(1,4)
w(1,4)=0.347854845137454
w(2,4)=0.652145154862546
w(3,4)=w(2,4)
w(4,4)=w(1,4)

!	Gauss point 5
xi(1,5)=-0.906179845938664
xi(2,5)=-0.538469310105683
xi(3,5)=0.D0
xi(4,5)=-xi(2,5)
xi(5,5)=-xi(1,5)
w(1,5)=0.236926885056189
w(2,5)=0.478628670499366
w(3,5)=0.568888888888889
w(4,5)=w(2,5)
w(5,5)=w(1,5)

!	Gauss Point 6
xi(1,6) = -0.932469514
xi(2,6) = -0.661209387
xi(3,6) = -0.238619186
xi(4,6) = -xi(3,6)
xi(5,6) = -xi(2,6)
xi(6,6) = -xi(1,6)
w(1,6) = 0.171324492
w(2,6) = 0.360761573
w(3,6) = 0.467913935
w(4,6) = w(3,6)
w(5,6) = w(2,6)
w(6,6) = w(1,6)

!	Gauss point 7
xi(1,7) = -0.949107912342758
xi(2,7) = -0.741531185599394
xi(3,7) = -0.405845151377397
xi(4,7) = 0.D0
xi(5,7) = -xi(3,7)
xi(6,7) = -xi(2,7)
xi(7,7) = -xi(1,7)
w(1,7) = 0.129484966168871
w(2,7) = 0.279705391489279
w(3,7) = 0.381830050505119
w(4,7) = 0.417959183673469
w(5,7) = w(3,7)
w(6,7) = w(2,7)
w(7,7) = w(1,7)

!	Gauss point 8
xi(1,8) = -0.960289857
xi(2,8) = -0.796666478
xi(3,8) = -0.525532410
xi(4,8) = -0.183434642
xi(5,8) = -xi(4,8)
xi(6,8) = -xi(3,8)
xi(7,8) = -xi(2,8)
xi(8,8) = -xi(1,8)
w(1,8) = 0.101228536
w(2,8) = 0.222381034
w(3,8) = 0.313706646
w(4,8) = 0.362683783
w(5,8) = w(4,8)
w(6,8) = w(3,8)
w(7,8) = w(2,8)
w(8,8) = w(1,8)

!	Gauss Point 9
xi(1,9) = -0.968160239507626
xi(2,9) = -0.836031107326636
xi(3,9) = -0.61337143270059
xi(4,9) = -0.324253423403809
xi(5,9) = 0
xi(6,9) = -xi(4,9)
xi(7,9) = -xi(3,9)
xi(8,9) = -xi(2,9)
xi(9,9) = -xi(1,9)
w(1,9) = 0.081274388361573
w(2,9) = 0.180648160694854
w(3,9) = 0.260610696402936
w(4,9) = 0.312347077040003
w(5,9) = 0.33023935500126
w(6,9) = w(4,9)
w(7,9) = w(3,9)
w(8,9) = w(2,9)
w(9,9) = w(1,9)

end Subroutine Gauss_point

Subroutine coordinate_transformation(ele_num,order,nnode,npoint,x,y,target_x,target_y,xxi)

implicit none
integer ele_num, order, nnode,l,m
integer npoint
real*8 x(nnode), y(nnode)
real*8 target_x, target_y
real*8 xxi

integer node(order+1)
real*8 total_dis, delta, dis

do l = 1, order+1
  node(l) = (order*(ele_num-1))+l
end do

total_dis = 0.0d0
do l =1, order
  delta = dsqrt((x(node(l+1))-x(node(l)))**2+(y(node(l+1))-y(node(l)))**2)
  total_dis = total_dis + delta
end do

dis = 0.0d0

if(npoint.ne.node(1)) then
  do m = 1, order
    if(node(m).eq.npoint) exit
    dis = dis + dsqrt((x(node(m+1))-x(node(m)))**2 +(y(node(m+1))-y(node(m)))**2)
  end do
end if

dis = dis + dsqrt((target_x-x(npoint))**2 +(target_y-y(npoint))**2)

xxi = ((2.0*dis)/total_dis)-1

end subroutine coordinate_transformation

Subroutine coordinate_transformation2(node,order,nnode,tindex,x,y,target_x,target_y,xxi)
implicit none
integer order, nnode, tindex
integer l, m

real*8 x(nnode), y(nnode)
real*8 target_x, target_y
real*8 xxi

integer node(order+1), temp_index
real*8 total_dis, delta, dis

total_dis = 0.0d0
do l =1, order
    delta = dsqrt((x(node(l+1))-x(node(l)))**2+(y(node(l+1))-y(node(l)))**2)
    total_dis = total_dis + delta
end do

dis = 0.0d0

do m = 1, order
  if(node(m+1).GE.tindex) then
    temp_index = m
    exit
  end if
  dis = dis + dsqrt((x(node(m+1))-x(node(m)))**2+(y(node(m+1))-y(node(m)))**2)
end do

dis = dis + dsqrt((target_x-x(node(temp_index)))**2+(target_y-y(node(temp_index)))**2)

xxi = ((2.0*dis)/total_dis)-1

end subroutine coordinate_transformation2

Subroutine Band_Solver(n, band_size,results,mass,load)

implicit none
integer n, band_size,i
real*8 results(n)
real*8 mass(n, band_size)
real*8 load(n)

do i = 1, n
  results(i) = 0.0d0
end do

call TRIB(n, band_size,mass)
call RHSB(n, band_size,results,mass,load)

end Subroutine Band_Solver

Subroutine TRIB(n, band_size,mass)

implicit none
integer n, band_size,i,j,k
real*8 mass(n, band_size)

integer m1,k1
real*8 sum

do i = 2, n
  m1 = min0(band_size-1,n-i+1)
  do j =1, m1
    sum = 0
    k1 = min0(i-1,band_size-j)
    do k = 1, k1
      sum = sum+mass(i-k,k+1)*mass(i-k,j+k)/mass(i-k,1)
    end do
    mass(i,j) = mass(i,j)-sum
  end do
end do

End Subroutine TRIB

Subroutine RHSB(n, band_size,results,mass,load)

implicit none
integer n, band_size,j,k
real*8 results(n)
real*8 mass(n, band_size)
real*8 load(n)

integer i,np1,k1,j1,j2,mm
real*8 sum

np1 = n+1

do i = 2, n
  sum = 0
  k1 = min0(band_size-1,i-1)
  do k = 1, k1
    sum = sum+mass(i-k,k+1)/mass(i-k,1)*load(i-k)
  end do
  load(i) = load(i)-sum
end do

results(n) = load(n) / mass(n,1)

do k = 2, n
  i = np1-k
  j1 = i+1
  j2 = min0(n,i+band_size-1)
  sum = 0
  do j = j1, j2
    mm = j-j1+2
    sum = sum+results(j)*mass(i,mm)
  end do
  results(i) = (load(i) - sum) / mass(i,1)
end do

END Subroutine RHSB

subroutine Assembly(location,tnNode,order,mass,global_mass)

implicit none
integer location, tnNode, order
real*8 mass(order+1,order+1)
real*8 global_mass(tnNode, order+1)

integer a, b
integer index1, index2, index3

index1 = location
index2 = 1

do a = 1, order+1
  index3 = a
  do b = 1, order + 2 - a
    global_mass(index1,b) = global_mass(index1,b)+ mass(index2,index3)
    index3 = index3+1
  end do
  index1 = index1 + 1
  index2 = index2 + 1
end do

end subroutine Assembly

END MODULE DATA_TRANSFER_MOD_2D

MODULE SRMS_DATATRANS_MOD_2D
!---------------------------------------------------------------------------------!
! Description:                                                                    !
!   Data Transfer Using Common-refinement                                         !
!---------------------------------------------------------------------------------!
!                                                    Written by: MDH & KDH & LCS  !
!---------------------------------------------------------------------------------!
! History:                                                                        !
! Version   Date      Comment                                                     !
! -------   ----      -------                                                     !
! 1.0     14/11/03    Original Code                                               !
! 2.0     15/02/02    DATA_TRANSFER2 subroutine added for unstructured data struc !
!---------------------------------------------------------------------------------!
IMPLICIT NONE; PRIVATE
PUBLIC DATA_TRANSFER2
PUBLIC DATA_TRANSFER

CONTAINS

SUBROUTINE DATA_TRANSFER2(snNode, sxyz, snFace, sf2n, sdata, tnNode, txyz, tnFace, tf2n, tdata)

USE DATA_TRANSFER_MOD_2D
! In or Out Variables
INTEGER, INTENT(in   ) :: snNode
REAL(8), INTENT(in   ) :: sxyz(2,snNode)
INTEGER, INTENT(in   ) :: snFace
INTEGER, INTENT(in   ) :: sf2n(2,snFace)
REAL(8), INTENT(in   ) :: sdata(snNode)
INTEGER, INTENT(in   ) :: tnNode
REAL(8), INTENT(in   ) :: txyz(2,tnNode)
INTEGER, INTENT(in   ) :: tnFace
INTEGER, INTENT(in   ) :: tf2n(2,tnFace)
REAL(8), INTENT(  out) :: tdata(tnNode)
! Local Variables
INTEGER :: psNode(snNode)
INTEGER :: ptNode(tnNode)
REAL(8) :: sxyz2(2,snNode)
REAL(8) :: sdata2(snNode)
REAL(8) :: txyz2(2,tnNode)
REAL(8) :: tdata2(tnNode)

! Re-arrange node array
CALL REARRANGE(snNode,snFace,sf2n,psNode)
CALL REARRANGE(tnNode,tnFace,tf2n,ptNode)

! Set source coordi & data 
sxyz2(:,:) = sxyz(:,psNode(:))
txyz2(:,:) = txyz(:,ptNode(:))
sdata2(:) = sdata(psNode(:))


! Data transfer using common-refinement
CALL DATA_TRANSFER(snNode, sxyz2, sdata2, tnNode, txyz2, tdata2)

! Set target data
tdata(ptNode(:)) = tdata2(:)

CONTAINS

SUBROUTINE REARRANGE(nNode, nFace, f2n, pNode)

! In or Out Variables
INTEGER, INTENT(in   ) :: nNode
INTEGER, INTENT(in   ) :: nFace
INTEGER, INTENT(in   ) :: f2n(2,nFace)
INTEGER, INTENT(  out) :: pNode(nNode) 
! Local Variables
INTEGER :: Face, tempFace
INTEGER :: n1, n2, n3, n4
INTEGER :: fpoint(nFace)

tempFace = 1
fpoint(1) = 1
DO WHILE(tempFace<nFace)
  DO Face = 2, nFace
    n1 = f2n(1,Face)
    n2 = f2n(2,Face)
    n3 = f2n(1,fpoint(1))
    n4 = f2n(2,fpoint(tempFace))
    IF( n2 == n3 ) THEN
      tempFace = tempFace + 1 
      fpoint(2:tempFace) = fpoint(1:tempFace-1)
      fpoint(1) = Face
      EXIT
    ELSE IF( n1 == n4 ) THEN
      tempFace = tempFace + 1
      fpoint(tempFace) = Face
      EXIT
    END IF
  END DO
END DO

DO Face = 1, nFace
  pNode(Face) = f2n(1,fpoint(Face))
END DO
pNode(nFace+1) = f2n(2,fpoint(nFace))

END SUBROUTINE REARRANGE

END SUBROUTINE DATA_TRANSFER2

SUBROUTINE DATA_TRANSFER(snNode, sxyz, sdata, tnNode, txyz, tdata)

USE DATA_TRANSFER_MOD_2D
! In or Out Variables
INTEGER, INTENT(in   ) :: snNode
REAL(8), INTENT(in   ) :: sxyz(2,snNode)
REAL(8), INTENT(in   ) :: sdata(snNode)
INTEGER, INTENT(in   ) :: tnNode
REAL(8), INTENT(in   ) :: txyz(2,tnNode)
REAL(8), INTENT(  out) :: tdata(tnNode)
! Local Variables
REAL(8), ALLOCATABLE :: nsxyz(:,:)
REAL(8), ALLOCATABLE :: ntxyz(:,:)

! Normalize Source / Target Mesh
ALLOCATE( nsxyz(2,snNode) )
ALLOCATE( ntxyz(2,tnNode) )
CALL LENGTH_NORMAL(snNode, sxyz, nsxyz)
CALL LENGTH_NORMAL(tnNode, txyz, ntxyz)

! Data Transfer
CALL INTERPOLATION(snNode, tnNode, nsxyz(1,:), nsxyz(2,:), sdata(:), ntxyz(1,:), ntxyz(2,:), tdata(:), 1)

CONTAINS

SUBROUTINE LENGTH_NORMAL(nNode, xyz, nxyz)

! In or Out Variables
INTEGER, INTENT(in   ) :: nNode
REAL(8), INTENT(in   ) :: xyz(2,nNode)
REAL(8), INTENT(  out) :: nxyz(2,nNode)
! Local Variables
INTEGER :: Node
REAL(8) :: Length

nxyz(1,1) = 0d0
DO Node = 2, nNode
  Length = DSQRT( SUM( ( xyz(:,Node) - xyz(:,Node-1) )**2 ) )
  nxyz(1,Node) = nxyz(1,Node-1) + Length
END DO

Length = nxyz(1,nNode)
DO Node = 1, nNode
  nxyz(1,Node) = nxyz(1,Node) / Length * 100d0
  nxyz(2,Node) = 0d0
END DO

END SUBROUTINE LENGTH_NORMAL

END SUBROUTINE DATA_TRANSFER

END MODULE SRMS_DATATRANS_MOD_2D
