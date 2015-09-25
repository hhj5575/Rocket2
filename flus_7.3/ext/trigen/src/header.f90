module limitvars
    integer, parameter :: mcmp = 5
    integer, parameter :: mbdy = 10000
    integer, parameter :: mnode = 500000
    integer, parameter :: medge = 3*mnode
    integer, parameter :: mcell = 2*mnode
    integer, parameter :: mtest = 10000
    integer, parameter :: mquad = 300000
    integer, parameter :: msrch = 5000
end module

module pntvars
    real(8), dimension(:,:), allocatable :: x
    integer :: nnode
end module

module edgvars
    integer, dimension(:,:), allocatable :: ndg
    integer, dimension(:), allocatable :: nfrt
    integer :: nedge, nfront, nrept
end module

module trivars
    integer, dimension(:,:), allocatable :: ndc, nbh
    real(8), dimension(:), allocatable :: xc, yc, rc
    integer :: ncell
end module

module cndvars
    integer, dimension(:), allocatable :: ntest
    real(8), dimension(:), allocatable :: xcen, ycen
    integer :: itest
end module

module lnkvars
    integer, dimension(:), allocatable :: ipoint, npoint, ncelpt
    real(8), dimension(:), allocatable :: dens
end module

module frtvars
    integer, dimension(:), allocatable :: nfill
    integer :: newfr
end module

module bdyvars
    real(8), dimension(:,:), allocatable :: xint, yint
    real(8), dimension(:), allocatable :: xout, yout
    integer, dimension(:), allocatable :: nint
    integer, dimension(:), allocatable :: patch
    integer :: ncmp, nout
end module

module qudvars
    real(8), dimension(:,:), allocatable :: xquad, yquad
    integer, dimension(:,:), allocatable :: nquad
    integer :: iquad
end module

module typvars
    integer, dimension(:), allocatable :: ityp
    real(8), dimension(:,:), allocatable :: xnrm
    integer :: nbnd1, nbnd2, ktyp
end module

module ctlvars
    integer :: ncode
end module

module bptvars
    integer :: it
end module

module nrmvars
    integer, dimension(:), allocatable :: nstop
    integer :: lcel
end module

module ordvars
    integer, dimension(:), allocatable :: iorder, listc
    integer :: nwait
end module

module rebvars
    integer, dimension(:), allocatable :: nacpt
end module
