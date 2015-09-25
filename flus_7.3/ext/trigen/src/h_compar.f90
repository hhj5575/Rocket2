!---------------------------------------------------------------------!
!                                                                     !
!    CHECK WHETHER NEW EDGE IS ALREADY IN THE EDGE LIST               !
!                                                                     !
!---------------------------------------------------------------------!

subroutine compar( l1, l2, nold, iedg )

    use pntvars
    use edgvars
    use lnkvars
    
    implicit none
    integer, intent(in) :: l1, l2
    integer, intent(inout) :: nold, iedg
    integer :: i
    
    i          = ipoint(l1)
    iedg       = i
    nold       = 0
10  if (i.eq.0) return
    if (l2.eq.max0(ndg(1,i),ndg(2,i))) go to 20
    iedg       = i
    i          = npoint(i)
    go to 10
20  nold       = i
    return

end subroutine