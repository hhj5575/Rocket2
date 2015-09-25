!---------------------------------------------------------------------!
!                                                                     !
!                REMOVE TRIANGLE LCEL FROM ORDERED LIST               !
!                                                                     !
!---------------------------------------------------------------------!

subroutine pulout(lcel)

    use pntvars
    use trivars
    use ordvars
    
    implicit none
    integer, intent(in) :: lcel
    integer :: ibeg, nwtm1, i, lnext
    
    ibeg = iorder(lcel)
    if(ibeg.eq.nwait) goto 20
    nwtm1 = nwait - 1
    do 10 i=ibeg,nwtm1
    lnext = listc(i+1)
    iorder(lnext) = i
    listc(i) = lnext
10  continue
20  nwait = nwait - 1
    return
      
end subroutine