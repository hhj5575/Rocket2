!---------------------------------------------------------------------!
!                                                                     !
!                    ADD TRIANGLE LCEL TO ORDERED LIST                !
!                                                                     !
!---------------------------------------------------------------------!

subroutine putin(lcel)

    use pntvars
    use trivars
    use lnkvars
    use ordvars
    
    implicit none
    integer, intent(in) :: lcel
    integer :: ii, m1, m2, m3, ltest, n1, n2, n3
    integer :: ia, ib, imid, i, ixch, lxch
    real(8) :: rcel, rtest
    
    m1 = ndc(1,lcel)
    m2 = ndc(2,lcel)
    m3 = ndc(3,lcel)
    rcel = dens(m1) + dens(m2) + dens(m3)
    ii = 1
      
    if (nwait.eq.0) goto 25
    
    ltest = listc(1)
    n1 = ndc(1,ltest)
    n2 = ndc(2,ltest)
    n3 = ndc(3,ltest)
    rtest = dens(n1) + dens(n2) + dens(n3)
    if (rc(lcel)*rtest.gt.rc(ltest)*rcel) goto 17
    
    ii = nwait + 1
    ltest = listc(nwait)
    n1 = ndc(1,ltest)
    n2 = ndc(2,ltest)
    n3 = ndc(3,ltest)
    rtest = dens(n1) + dens(n2) + dens(n3)
    if (rc(lcel)*rtest.le.rc(ltest)*rcel) goto 25
    
    ia = 1
    ib = nwait
5   imid = (ia+ib)/2
    ltest = listc(imid)
    n1 = ndc(1,ltest)
    n2 = ndc(2,ltest)
    n3 = ndc(3,ltest)
    rtest = dens(n1) + dens(n2) + dens(n3)
    if (rc(lcel)*rtest.le.rc(ltest)*rcel) goto 10
    ib = imid
    goto 15    
10  ia = imid
15  if(ib.gt.ia+1) goto 5
    ii = ia
    if (ia.ne.ib) ii=ib
      
17  do 20 i=ii,nwait
    ixch = nwait+ii-i
    lxch = listc(ixch)
    listc(ixch+1) = lxch
    iorder(lxch) = ixch + 1  
20  continue
    
25  listc(ii) = lcel
    iorder(lcel) = ii
    nwait = nwait + 1
    return
    
end subroutine
