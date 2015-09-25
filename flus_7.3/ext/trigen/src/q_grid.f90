!---------------------------------------------------------------------!
!                                                                     !
!             WRITE OUT FILE WITH LIST OF POINTS AND TRIANGLES.       !
!                                                                     !
!---------------------------------------------------------------------!
subroutine grid(nNode1, nCell1, xyz1, c2n1)

    use pntvars
    use trivars
    use bdyvars
    implicit none
    INTEGER, INTENT(  out) :: nNode1
    INTEGER, INTENT(  out) :: nCell1
    REAL(8), INTENT(  out) :: xyz1(2,*)
    INTEGER, INTENT(  out) :: c2n1(3,*)
    integer :: i, ncnt
    
    ncnt = 0
    do 30 i=1,ncell
    if(ndc(1,i).eq.-1) goto 30
    ncnt = ncnt + 1
    ndc(1,ncnt) = ndc(1,i)
    ndc(2,ncnt) = ndc(2,i)
    ndc(3,ncnt) = ndc(3,i)
30  continue
    ncell = ncnt
    
    nNode1 = nnode
    xyz1(:,1:nnode) = x(:,1:nnode)
    nCell1 = ncell
    c2n1(:,1:ncell) = ndc(:,1:ncell)
    
!    open(10,file='./output/mesh.plt')
!    write(10,*) 'TITLE = "Mesh"'
!    write(10,*) 'VARIABLES= "X", "Y"'
!    write(10,*) 'ZONE F=FEPOINT, N=',nnode,'E=',ncell,'ET=TRIANGLE'
!    do 10 i=1,nnode
!10  write(10,*) x(1,i),x(2,i)
!    do 20 i=1,ncell
!20  write(10,105) ndc(1,i),ndc(2,i),ndc(3,i)
!    close(10)
    
105 format(3i10)
600 format(/5x,'triangulation complete'//5x,'number of points = ',i5//5x,'number of cells  = ',i5)
610 format(/5x,'viscous layer completed')
      
end subroutine

subroutine midtest

    use pntvars
    use trivars
    use bptvars
    
    implicit none
    integer :: i, ncell2
    character(16) :: filename
    
    ncell2 = 0
    do 30 i = 1,ncell
        if( ndc(1,i) .eq. -1 ) goto 30
        ncell2 = ncell2 + 1
30  end do
    
    write(filename,"('midtest',i5.5,'.plt')") it
    open(1,file='./output/'//filename)
    write(1,*) 'TITLE = "test"'
    write(1,*) 'VARIABLES= "X", "Y"'
    write(1,*) 'ZONE F=FEPOINT, N=',nnode,',E=',ncell2,', ET=TRIANGLE'
    
    do i = 1,nnode
        write(1,*) x(1,i), x(2,i)
    end do
    
    do 20 i = 1,ncell
        if( ndc(1,i) .eq. -1 ) goto 20
        write(1,*) ndc(1,i), ndc(2,i), ndc(3,i)
20  end do
    
10  close(1)
    
end subroutine
