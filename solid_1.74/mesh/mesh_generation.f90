!****************************************************************************
!
!  PROGRAM : paving
!
!  PURPOSE : 2D Mesh Generation
!
!****************************************************************************
!include 'mesh_module.f90'

subroutine mesh_gen_main(cs, rn, ngroup, ng, tol_nn, tol_ne, ori_nn, new_node, ori_ne, new_ele, flag)

use ale_domain

implicit none

integer, intent(in) :: cs, rn, ng, ngroup
integer, intent(inout) :: tol_nn, tol_ne, ori_nn, ori_ne, new_ele(4,ori_ne*3), flag
real, intent(inout) :: new_node(2,ori_nn*3)

real,allocatable :: node(:,:), b(:,:), dd(:), p(:,:), div_b(:,:)
real,allocatable :: temp_node(:,:), edge(:,:), ori_edge(:,:)
integer, allocatable :: ele(:,:), db(:,:), dbs(:), line(:), b_num(:)
integer, allocatable :: temp_ele(:,:), temp_db(:,:)

integer :: pre_nn, pre_ne

real :: df, d, ori_d, mean_dis, dis
integer :: i, j, status, pn, nn, bn, anes, dn, whole_i, temp_bs, ne, count
integer :: ins_num, first_nn, enn, on, op_num, cn1, cn2, iter, max_iter, gen_type
logical :: ins, iter_check
!==============================================================================
enn = 0;  op_num = 0;  ins = .FALSE.;  ins_num = 0;  on = 1;  max_iter = 13
if (flag == 5) then
    open (Unit=30, File='partition_check.plt', STATUS='old', ACTION='read', IOSTAT=status)
    read(30,*) nn
    allocate (p(2,nn), line(nn), edge(2,nn), ori_edge(2,nn))
    do i = 1, nn
        read (30,*) p(:,i)
    enddo
    close(30)
    gen_type = 3
    ori_edge = 0.0
    ori_d = 0.0
else
    !call calc_ori_d(rn, ng, ori_d)
    cn1 = mesh_set(rn)%corn_group(ng)
    if (ng == ngroup) then
        cn2 = mesh_set(rn)%cn
    else
        cn2 = mesh_set(rn)%corn_group(ng+1) - 1
    endif

    pn = cn2 - cn1 + 1
    allocate ( p(2,pn), line(pn), edge(2,pn), ori_edge(2,pn) )
   
    !write (*,*) '>> Original edge'
    count = 0
    do i = cn1, cn2
        count = count + 1
        line(count) = 1
        !p(:,count) = surface_set(ng)%coord(:,mesh_set(rn)%corn(i))
        p(:,count) = mesh_set(rn)%corn_coord(:,i)
        ori_edge(:,count) = mesh_set(rn)%edge(:,i)
        !write (*,'(A,I4,4F12.7)') 'edge:', i, ori_edge(:,count), p(:,count)
    enddo
    gen_type = flag
    
    ! check corner
    open (Unit=30, File='./output/solid/remesh/check_corner.plt', STATUS='replace', ACTION='write', IOSTAT=status)
    do i = 1, count
        Write (30,*) p(:,i)
    enddo
    close(30)
    
    edge = ori_edge
    mean_dis = mesh_set(rn)%ave_d(ng)
    d = mean_dis
    write (*,*) 'mean_dis, ng, flag:', mean_dis, ng, flag
    if (flag >= 3) then
        gen_type = 2
        flag = 1
    endif
endif

do iter = 1, max_iter
    iter_check = .TRUE.
    write (*,'(I3,A)') iter, " iteration."
    if (gen_type == 1) then
        allocate(div_b(2,pn*500))
       
        if (iter == 1) then
            edge = ori_edge
        else
            !count = iter/2
            !edge = ori_edge * (1.0 + ((-1)**iter)*(count-1)*0.01)
            !d = ori_d * (1.05 - (iter-1)*0.01)
            !edge = ori_edge * (1.07 - (iter-1)*0.01)    ! temporory
            !d = ori_d * (1.07 - (iter-1)*0.01)    ! temporory
            do i = 1, pn
                do j = 1, 2
                    if (edge(j,i) <= mean_dis) then
                        edge(j,i) = edge(j,i) * 1.02
                    else
                        edge(j,i) = edge(j,i) * 0.98
                    endif
                enddo
            enddo
        endif
        !edge = 0.0035                           ! temporory 
        !edge = edge * (1.04 - (iter-1)*0.02)    ! temporory
        call divide_point(pn, p, d, line, bn, div_b, edge)
        nn = bn
        flag = 1
    elseif (gen_type == 2) then
        nn = surface_set(ng)%nn
        allocate (div_b(2,nn*2))
        div_b(:,1:nn) = surface_set(ng)%coord(:,1:nn)
        gen_type = 1
    elseif (gen_type == 3) then
        allocate (div_b(2,nn*2))
        div_b(:,1:nn) = p(:,1:nn)
    endif
    
    anes = nn * 100;  bn = nn * 10
    ne = 0;  first_nn = nn;  enn = nn
    allocate(dbs(nn), db(nn,bn))
    allocate(node(2,anes), ele(5,anes))
    
    node = 0.0;  ele = 0
    node(:,1:nn) = div_b(:,1:nn)
    
    db = 0;  dbs = 0
    dbs(1) = nn
    do i = 1, nn
        db(1,i) = i
    enddo

    !d = temp_d
    call partition(nn, anes, node, bn, db, dbs, dn, first_nn)

    do i = 1, dn
        if ( 2 * (dbs(i) / 2) /= dbs(i) ) then
            iter_check = .FALSE.
            write (*,*) 'mesh_gen: partition - false'
            allocate(dd(dn))
            goto 100
        endif
    enddo
    
    allocate( dd(dn), temp_db(dn,nn) )
    dd = 0.0
    temp_db = db(1:dn,1:nn)
    call calc_d(anes, nn, node, dn, dbs(1:dn), temp_db, dd)

    deallocate ( temp_db )

    do whole_i = 1, dn
	    temp_bs = dbs(whole_i)
        allocate ( b_num(temp_bs) )
	    allocate ( b(2,temp_bs) )
	    do j = 1, temp_bs
		    b_num(j) = db(whole_i,j)
		    b(:,j) = node(:,db(whole_i,j))
	    enddo
	    df = dd(whole_i)
	    call mesh_gen( temp_bs, b_num, b, df, anes, nn, node, ne, ele, whole_i, iter_check )
	    deallocate ( b_num, b )
        
	    if ( iter_check == .FALSE. ) goto 100
    enddo
    
    pre_nn = first_nn
    pre_ne = 1
    
    call combines_two_ele(anes, pre_nn, nn, ne, node, ele, iter_check)
    if ( iter_check == .FALSE. ) then
        write (*,*) 'mesh_gen: combines_two_ele - false'
        goto 100
    endif
    
    call divide_bound_node(anes, pre_nn, nn, ne, node, ele, iter_check)
    if ( iter_check == .FALSE. ) goto 100
    
    call divide_over_ele(anes, nn, ne, node, ele, pre_nn, pre_ne, iter_check)
    if ( iter_check == .FALSE. ) then
        write (*,*) 'mesh_gen: divide_over_ele - false'
        goto 100
    endif
    
    call combines_two_ele(anes, pre_nn, nn, ne, node, ele, iter_check)
    if ( iter_check == .FALSE. ) then
        write (*,*) 'mesh_gen: combines_two_ele - false'
        goto 100
    endif

    call check_ele_angle(anes, ne, node, ele, first_nn, pre_ne)
    
    call inner_soft_object(anes, pre_nn, pre_ne, nn, ne, node, ele, iter_check)

    100 continue
    if ( iter_check ) then
        exit
    else
        call write_tecplot( nn, node(:,1:nn), ne, ele(2:5,1:ne), rn, iter, 6 )
        write (*,*) 'Fail generating mesh at ', iter, ' iteration!!'
        deallocate ( dbs, db, node, ele, dd, div_b )
    endif
enddo
deallocate (p, line, div_b, edge, ori_edge)

if (flag /= 5) then
    surface_set(ng)%adjn = pre_nn
    allocate (surface_set(ng)%adj(pre_nn))
    do i = 1, pre_nn
        surface_set(ng)%adj(i) = i + tol_nn
    enddo
endif

allocate (temp_node(2,nn), temp_ele(4,ne))
temp_node = node(:,1:nn)
temp_ele = ele(2:5,1:ne)
call write_tecplot( nn, temp_node, ne, temp_ele, rn, ng, 2 )

if (flag /= 5) then
    surface_set(ng)%num_nodes = tol_nn + 1
    surface_set(ng)%num_elements = tol_ne + 1
    new_node(:,1+tol_nn:nn+tol_nn) = node(:,1:nn)
    new_ele(:,1+tol_ne:ne+tol_ne) = ele(2:5,1:ne)+tol_nn
    tol_nn = tol_nn + nn
    tol_ne = tol_ne + ne
    write (*,'(2(A,I7))') 'Total number of node :', tol_nn, ', Total number of element :', tol_ne
endif 
if (flag == 4) stop

deallocate (db, dbs, dd)
deallocate (node, ele)
deallocate (temp_node, temp_ele)

end subroutine mesh_gen_main
!==============================================================================

subroutine meshless(rn, ng, tol_nn, tol_ne, ori_nn, new_node, ori_ne, new_ele)
                    
use ale_domain
use physical_domain

implicit none

integer, intent(in) :: rn, ng, ori_nn, ori_ne
integer, intent(inout) :: tol_nn, tol_ne, new_ele(4,ori_ne*2)
real, intent(inout) :: new_node(2,ori_nn*2)

integer :: i, pre_nn, nn, ne
integer :: ngroup, nn1, nn2, ne1, ne2, ns1, ns2
integer, allocatable :: ele(:,:)
real, allocatable :: node(:,:)

ngroup = region(rn)%num_groups
nn1 = region(rn)%node_group(ng)
ne1 = region(rn)%element_group(ng)
ns1 = region(rn)%skin_group_pointer(ng)
if (ng < ngroup) then
    nn2 = region(rn)%node_group(ng+1) - 1
    ne2 = region(rn)%element_group(ng+1) - 1
    ns2 = region(rn)%skin_group_pointer(ng+1) - 1
else
    nn2 = region(rn)%num_nodes
    ne2 = region(rn)%num_elements
    ns2 = region(rn)%num_skin
endif
pre_nn = nn1 - 1     
nn = nn2 - nn1 + 1
ne = ne2 - ne1 + 1
surface_set(ng)%adjn = ns2 - ns1 + 1
allocate ( node(2,nn), ele(4,ne), surface_set(ng)%adj(surface_set(ng)%adjn) )
if (rn == 1) then
    node = mesh_set(rn)%coord(:,nn1:nn2)
else
    node = region(rn)%node_coord(:,nn1:nn2)
endif
ele = region(rn)%element(4:7,ne1:ne2) - pre_nn
surface_set(ng)%adj = region(rn)%skin_nodes(ns1:ns2) - pre_nn + tol_nn

call write_tecplot( nn, node, ne, ele, rn, ng, 2 ) 
do i = 1, nn
    !write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i+tol_nn, ', ', node(1,i), ', ', node(2,i)
    new_node(:,tol_nn+i) = node(:,i)
enddo

do i = 1, ne
    !write (ele_unit_num,'(I5,A,I5,A,I5,A,I5,A,I5)') i+tol_ne,', ',ele(1,i)+tol_nn,', ',ele(2,i)+tol_nn,', ',ele(3,i)+tol_nn,', ',ele(4,i)+tol_nn
    new_ele(:,tol_ne+i) = ele(:,i)+tol_nn
enddo

surface_set(ng)%num_nodes = tol_nn + 1
surface_set(ng)%num_elements = tol_ne + 1
!new_node(:,1+tol_nn:nn+tol_nn) = node(:,1:nn)
!new_ele(:,1+tol_ne:ne+tol_ne) = ele(:,1:ne)
tol_nn = tol_nn + nn
tol_ne = tol_ne + ne
write (*,*) 'tol_nn :', tol_nn, ', tol_ne :', tol_ne

deallocate ( node, ele ) 

end subroutine meshless
!==============================================================================
subroutine write_tecplot( nn, node, ne, ele, rn, ng, num ) 

implicit none 

integer, intent(in) :: nn , ne, ele(1:4,1:ne), rn, num, ng
real, intent(in) :: node(1:2,1:nn)

integer :: status, i
character(len=50) :: fn

if ( num == 1 ) then
    fn = "./output/solid/remesh/2D_rempre00_00.plt"
elseif ( num == 2 ) then
    fn = "./output/solid/remesh/2D_remesh00_00.plt"
elseif ( num == 3 ) then
    fn = "./output/solid/move/2D_movepr00_00.plt"
elseif ( num == 4 ) then
    fn = "./output/solid/move/2D_movemd00_00.plt"
elseif ( num == 5 ) then
    fn = "./output/solid/move/2D_moveaf00_00.plt"
elseif ( num == 6 ) then
    fn = "./output/solid/remesh/2D_fail00_00.plt"
elseif ( num == 7 ) then
    fn = "./output/solid/move/2D_surf00_000000.plt"
endif

if (num <= 3) then
    write ( fn(32:33), '(I2.2)' ) rn
    write ( fn(35:36), '(I2.2)' ) ng
elseif (num == 7) then
    write ( fn(28:29), '(I2.2)' ) rn
    write ( fn(31:36), '(I6.6)' ) ng
else
    write ( fn(30:31), '(I2.2)' ) rn
    write ( fn(33:34), '(I2.2)' ) ng
endif
open (Unit=20, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)

!fn = "remesh00.txt"
!write ( fn(7:8), '(I2.2)' ) rn
!open (Unit=21, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)


Write (20,'(A,A,A)') 'TITLE="' ,fn(10:20) ,'"'
Write (20,*) 'VARIABLES="x", "y"'
if (num == 7) then
    Write (20,'(A,I5,A,I5, A)') 'ZONE N =', nn , ', E =', nn, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
else
    Write (20,'(A,I5,A,I5, A)') 'ZONE N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
endif
!write (21,*) nn, ne

do i = 1, nn
	write (20,*) node(:,i)
!	write (21,*) node(:,i)
enddo
if (num == 7) then
    do i = 1, nn
        if (i == nn) then
            write (20,*) i, 1
        else
	        write (20,*) i, i+1
        endif
    enddo
else
    do i = 1, ne
	    write (20,*) ele(:,i)
    !	write (21,*) ele(:,i)
    enddo
endif
close(unit=20)
!close(unit=21)

end subroutine 
!==============================================================================

!==============================================================================
subroutine write_elem( filename, on, ele_num) 

implicit none 

character(len=20), intent(in) :: filename 
integer, intent(in) :: on
integer, intent(in) :: ele_num(on)

integer :: status, i
character(len=30) :: fn

fn = ".\ele_num\"
write ( fn(11:30), '(A)' ) filename
open (Unit=30, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)

Write (30,'(A,I3)') 'Number of objects = ', on
do i = 1, on
	write (30,*) i, ele_num(i)
enddo

close(unit=30)

end subroutine
!==============================================================================
subroutine write_input(nn, ne, coord, ele, node_unit_num, ele_unit_num, flag)

implicit none

integer, intent(in) :: nn, ne, node_unit_num, ele_unit_num, ele(4,ne*2), flag
real, intent(in) :: coord(2,nn*2)

! flag - 1 : not seperation
!      - 2 : seperation
integer :: i

if ( flag == 2 ) then
    do i = 1, nn
        write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i, ', ', coord(1,i), ', ', coord(2,i)
    enddo
    write (node_unit_num,*)
    do i = 1, ne
        write (ele_unit_num,'(I5,A,I5,A,I5,A,I5,A,I5)') i,', ',ele(1,i),', ',ele(2,i),', ',ele(3,i),', ',ele(4,i)
    enddo
    write (ele_unit_num,*)
endif

!call write_tecplot( nn, coord, ne, ele, rn, 1, 2 ) 

end subroutine
!==============================================================================
subroutine write_tecplot_3D(nn, node, ne, ele, rn, ng, num) 

implicit none 

integer, intent(in) :: nn , ne, ele(8,ne), rn, num, ng
real, intent(in) :: node(3,nn)

integer :: status, i
character(len=50) :: fn

if ( num == 1 ) then
    fn = "./output/solid/remesh/3D_rempre00_00.plt"
elseif ( num == 2 ) then
    fn = "./output/solid/remesh/3D_remesh00_00.plt"
elseif ( num == 3 ) then
    fn = "./output/solid/move/3D_movepr00_00.plt"
elseif ( num == 4 ) then
    fn = "./output/solid/move/3D_movemd00_00.plt"
elseif ( num == 5 ) then
    fn = "./output/solid/move/3D_moveaf00_00.plt"
elseif ( num == 6 ) then
    fn = "./output/solid/remesh/3D_fail00_00.plt"
elseif ( num == 7 ) then
    fn = "./output/solid/check_model00_00000.plt"
endif

if (num <= 2) then
    write ( fn(32:33), '(I2.2)' ) rn
    write ( fn(35:36), '(I2.2)' ) ng
elseif (num == 7) then
    write ( fn(27:28), '(I2.2)' ) rn
    write ( fn(30:34), '(I5.5)' ) ng
else
    write ( fn(30:31), '(I2.2)' ) rn
    write ( fn(33:34), '(I2.2)' ) ng
endif
open (Unit=20, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)

!fn = "remesh00.txt"
!write ( fn(7:8), '(I2.2)' ) rn
!open (Unit=21, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)

Write (20,'(A,A,A)') 'TITLE="' ,fn(10:20) ,'"'
Write (20,*) 'VARIABLES="x", "y", "z"'
Write (20,'(A,I5,A,I5, A)') 'ZONE N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
!write (21,*) nn, ne

do i = 1, nn
	write (20,*) node(:,i)
!	write (21,*) node(:,i)
enddo
do i = 1, ne
	write (20,'(8(I7,2X))') ele(:,i)
!	write (21,*) ele(:,i)
enddo
close(unit=20)
!close(unit=21)

end subroutine write_tecplot_3D
!==============================================================================