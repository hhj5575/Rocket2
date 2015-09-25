!-----------------------------------------------------------------!
!                                                                 !
!               ARRAY DIMENSION ALLOCATE AND DEALLOCATE           !
!                                                                 !
!-----------------------------------------------------------------!
    
subroutine allocvars
    
    use limitvars
    use pntvars
    use edgvars
    use trivars
    use cndvars
    use lnkvars
    use frtvars
    use bdyvars
    use qudvars
    use typvars
    use nrmvars
    use ordvars
    use rebvars
    
    implicit none
    
!   pntvars
    allocate( x(2,mnode) )
    
!   edgvars
    allocate( ndg(4,medge), nfrt(medge) )
    
!   trivars
    allocate( ndc(3,mcell), nbh(3,mcell) )
    allocate( xc(mcell), yc(mcell), rc(mcell) )
    
!   cndvars
    allocate( ntest(mtest) )
    allocate( xcen(mtest), ycen(mtest) )

!   lnkvars
    allocate( ipoint(mnode), npoint(medge), ncelpt(mnode) )
    allocate( dens(mnode) )
    
!   frtvars
    allocate( nfill(medge) )
    
!   bdyvars
    allocate( xint(mcmp,mbdy), yint(mcmp,mbdy) )
    allocate( xout(mbdy), yout(mbdy) )
    allocate( patch(mbdy) )
    allocate( nint(mcmp) )
    
!   qudvars
    allocate( xquad(2,mquad), yquad(2,mquad), nquad(7,mquad) )
    
!   typvars
    allocate( ityp(mnode) )
    allocate( xnrm(2, mnode) )
    
!   nrmvars
    allocate( nstop(mcell) )
    
!   ordvars
    allocate( iorder(mcell), listc(mcell) )
    
!   rebvars
    allocate( nacpt(mcell) )
    
end subroutine

subroutine dealloc

    use pntvars
    use edgvars
    use trivars
    use cndvars
    use lnkvars
    use frtvars
    use bdyvars
    use qudvars
    use typvars
    use nrmvars
    use ordvars
    use rebvars
    
    implicit none
    
!   pntvars
    deallocate( x )
    
!   edgvars
    deallocate( ndg, nfrt )
    
!   trivars
    deallocate( ndc, nbh )
    deallocate( xc, yc, rc )
    
!   cndvars
    deallocate( ntest )
    deallocate( xcen, ycen )

!   lnkvars
    deallocate( ipoint, npoint, ncelpt )
    deallocate( dens )
    
!   frtvars
    deallocate( nfill )
    
!   bdyvars
    deallocate( xint, yint )
    deallocate( xout, yout )
    deallocate( patch )
    deallocate( nint )
    
!   qudvars
    deallocate( xquad, yquad, nquad )
    
!   typvars
    deallocate( ityp )
    deallocate( xnrm )
    
!   nrmvars
    deallocate( nstop )
    
!   ordvars
    deallocate( iorder, listc )
    
!   rebvars
    deallocate( nacpt )
    
end subroutine