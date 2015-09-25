!-----------------------------------------------------------------!
!                                                                 !
!    SETUP POINT LIST AND EDGE DATA STRUCTURE FOR BOUNDARY NODES  !
!                                                                 !
!-----------------------------------------------------------------!

subroutine initbd
    
    use limitvars
    use pntvars
    use edgvars
    use trivars
    use lnkvars
    use bdyvars
    use typvars
    
    implicit none
    integer :: l, i, it, n, np1, nm1, l1, l2
    integer :: nstrt, inext, ltyp
    real(8) :: xmin, ymin, xmax, ymax
    real(8) :: dx1, dy1, ds1, dx2, dy2, ds2
    real(8) :: dxn, dyn, dsn, ds
    real(8) :: relax, rn, fac
    real(8) :: chord(mcmp), chmax
    
    
    nedge = 0
    nnode = 0
    nbnd1 = 0
    ktyp = 0
    
    if( ncmp .eq. 0 ) goto 45
    
!   POINT AND EDGES FOR INTERNAL BOUNDARY
    
    do 40 l=1,ncmp
        
    ktyp=ktyp+1
    nnode=nnode+1
    nstrt=nnode
    x(1,nnode)=xint(l,1)
    x(2,nnode)=yint(l,1)
    xmin=x(1,nnode) 
    ymin=x(2,nnode) 
    xmax=x(1,nnode) 
    ymax=x(2,nnode) 
    ityp(nnode)=ktyp
    ipoint(nnode)=nnode
    
    dx1=xint(l,1)-xint(l,nint(l))
    dy1=yint(l,1)-yint(l,nint(l))
    ds1=dsqrt(dx1*dx1+dy1*dy1)
    dx2=xint(l,2)-xint(l,1)
    dy2=yint(l,2)-yint(l,1)
    ds2=dsqrt(dx2*dx2+dy2*dy2)
    dxn=dy1/ds1+dy2/ds2
    dyn=-dx1/ds1-dx2/ds2
    dsn=dsqrt(dxn*dxn+dyn*dyn)
    xnrm(1,nnode)=dxn/dsn
    xnrm(2,nnode)=dyn/dsn

    do 10 i=2,nint(l)
        
    nnode=nnode+1
    x(1,nnode)=xint(l,i)
    x(2,nnode)=yint(l,i)
    xmin=dmin1(xmin,x(1,nnode))
    ymin=dmin1(ymin,x(2,nnode))
    xmax=dmax1(xmax,x(1,nnode))
    ymax=dmax1(ymax,x(2,nnode))
    
    ityp(nnode)=ktyp
    inext=mod(i,nint(l))+1
    
    dx1=xint(l,i)-xint(l,i-1)
    dy1=yint(l,i)-yint(l,i-1)
    ds1=dsqrt(dx1*dx1+dy1*dy1)
    dx2=xint(l,inext)-xint(l,i)
    dy2=yint(l,inext)-yint(l,i)
    ds2=dsqrt(dx2*dx2+dy2*dy2)
    dxn=dy1/ds1+dy2/ds2
    dyn=-dx1/ds1-dx2/ds2     
    dsn=dsqrt(dxn*dxn+dyn*dyn)
    xnrm(1,nnode)=dxn/dsn
    xnrm(2,nnode)=dyn/dsn

    nedge=nedge+1
    nbh(1,nedge)=nnode+1
    if(i.eq.nint(l)) nbh(1,nedge)=nnode-nint(l)+1
    nbh(2,nedge)=nnode-2
    if(i.eq.2) nbh(2,nedge)=nint(l)+nnode-2
    ndg(1,nedge)=nnode
    ndg(2,nedge)=nnode-1
    ndg(3,nedge)=-1
    ndg(4,nedge)=0
    nfrt(nedge)=nedge
    ipoint(nnode)=nedge+1
    if(i.eq.nint(l)) ipoint(nnode)=0
    npoint(nedge)=0
    if(i.eq.2) npoint(nedge)=nint(l)+nedge-1
    
10  continue
    
    nedge=nedge+1
    nbh(1,nedge)=nnode-nint(l)+2
    nbh(2,nedge)=nnode-1
    ndg(1,nedge)=nnode-nint(l)+1
    ndg(2,nedge)=nnode
    ndg(3,nedge)=-1
    ndg(4,nedge)=0
    npoint(nedge)=0
    nfrt(nedge)=nedge

!   DETERMINE CHORD FOR ELEMENT AND SET STEP LENGTH
    
    chord(ktyp)=dsqrt((xmax-xmin)**2+(ymax-ymin)**2)
    
!   SMOOTH NORMALS FOR ELEMENT
    
    relax=0.5d0
    do 20 it=1,100 
    do 20 n=nstrt,nnode
    np1=n+1
    if(n.eq.nnode) np1=nstrt
    nm1=n-1
    if(n.eq.nstrt) nm1=nnode
    xnrm(1,n)=(1d0-relax)*xnrm(1,n) +0.5d0*relax*(xnrm(1,np1)+xnrm(1,nm1))
    xnrm(2,n)=(1d0-relax)*xnrm(2,n) +0.5d0*relax*(xnrm(2,np1)+xnrm(2,nm1))
20  continue
    
    do 30 n=nstrt,nnode
    rn=1d0/dsqrt(xnrm(1,n)**2+xnrm(2,n)**2)
    xnrm(1,n)=rn*xnrm(1,n)
    xnrm(2,n)=rn*xnrm(2,n)
30  continue
    
40  continue
    
    nbnd1=nnode
    
    chmax=chord(1)
    do 42 i=1,ktyp
    chmax=dmax1(chord(i),chmax)
42  continue
    
    do 43 i=1,ktyp
    chord(i)=chord(i)/chmax
43  continue
    
!   POINTS AND EDGES FOR EXTERNAL BOUNDARY
    
45  ktyp=ktyp+1
    nnode=nnode+1
    x(1,nnode)=xout(1)
    x(2,nnode)=yout(1)
    ityp(nnode)=ktyp
    ipoint(nnode)=nedge+1
    
    dx1=xout(1)-xout(nout)
    dy1=yout(1)-yout(nout)
    ds1=dsqrt(dx1*dx1+dy1*dy1)
    dx2=xout(2)-xout(1)
    dy2=yout(2)-yout(1)
    ds2=dsqrt(dx2*dx2+dy2*dy2)
    dxn=dy1/ds1+dy2/ds2
    dyn=-dx1/ds1-dx2/ds2
    dsn=dsqrt(dxn*dxn+dyn*dyn)
    xnrm(1,nnode)=-dxn/dsn
    xnrm(2,nnode)=-dyn/dsn
    
    do 50 i=2,nout
    
    nnode=nnode+1
    x(1,nnode)=xout(i)
    x(2,nnode)=yout(i)
    ityp(nnode)=ktyp
    inext=mod(i,nout)+1
    
    dx1=xout(i)-xout(i-1)
    dy1=yout(i)-yout(i-1)
    ds1=dsqrt(dx1*dx1+dy1*dy1)
    dx2=xout(inext)-xout(i)
    dy2=yout(inext)-yout(i)
    ds2=dsqrt(dx2*dx2+dy2*dy2)
    dxn=dy1/ds1+dy2/ds2
    dyn=-dx1/ds1-dx2/ds2
    dsn=dsqrt(dxn*dxn+dyn*dyn)
    xnrm(1,nnode)=-dxn/dsn
    xnrm(2,nnode)=-dyn/dsn
    
    nedge=nedge+1
    nbh(1,nedge)=nnode-2
    if(i.eq.2) nbh(1,nedge)=nout+nnode-2
    nbh(2,nedge)=nnode+1
    if(i.eq.nout) nbh(2,nedge)=nnode-nout+1
    ndg(1,nedge)=nnode-1
    ndg(2,nedge)=nnode  
    ndg(3,nedge)=-1   
    ndg(4,nedge)=0    
    nfrt(nedge)=nedge
    ipoint(nnode)=nedge+1
    if(i.eq.nout) ipoint(nnode)=0
    npoint(nedge)=0
    if(i.eq.2) npoint(nedge)=nedge+nout-1
    
50  continue
    
    nedge=nedge+1
    nbh(1,nedge)=nnode-1
    nbh(2,nedge)=nnode-nout+2
    ndg(1,nedge)=nnode          
    ndg(2,nedge)=nnode-nout+1   
    ndg(3,nedge)=-1             
    ndg(4,nedge)=0              
    nfrt(nedge)=nedge
    npoint(nedge)=0
    nfront=nedge
    chord(ktyp)=1d0
    nbnd2=nnode

!   INITIALIZE POINT DENSITY FUNCTION
    
    do 60 n=1,nnode
    dens(n)=0d0
60  continue
    
    do 70 l=1,nedge
    l1=ndg(1,l)
    l2=ndg(2,l)
    ds=dsqrt((x(1,l1)-x(1,l2))**2+(x(2,l1)-x(2,l2))**2)
    fac=0.58d0 ! fac 0.58? why?
    ltyp=ityp(l1)
    fac=fac*chord(ltyp)
    dens(l1)=dens(l1)+0.5d0*fac*ds
    dens(l2)=dens(l2)+0.5d0*fac*ds
70  continue

    return
    
end subroutine