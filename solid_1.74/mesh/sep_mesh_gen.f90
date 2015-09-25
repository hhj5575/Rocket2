!****************************************************************************
!
!  PROGRAM : paving
!
!  PURPOSE : 2D Mesh Generation
!
!****************************************************************************
subroutine sep_mesh_gen(crn, bs, ori_div_b, ori_d, tol_nn, tol_ne, ori_nn, ori_ne, new_node, new_ele, adjn)

use ale_domain

implicit none

!integer, intent(inout) :: crn, cnn, cne, bs, temp_nn, temp_ne, temp_ele(4,cne), tol_nn, tol_ne
!real, intent(inout) :: div_b(2,bs), d, temp_coord(2,cnn)
integer, intent(in) :: crn, adjn, ori_nn, ori_ne
integer, intent(inout) :: bs, tol_nn, tol_ne, new_ele(4,ori_ne*4)
real, intent(inout) :: ori_div_b(2,bs), ori_d, new_node(2,ori_nn*4)

real,allocatable :: node(:,:), b(:,:), dd(:), corner_node(:,:), p(:,:), op(:,:), div_b(:,:), edge(:,:)
integer, allocatable :: ele(:,:), db(:,:), dbs(:), opn(:), line(:), b_num(:), temp_corn(:,:), temp_db(:,:)

! 각 영역에 대한 element정보를 넣기 위한 변수
integer :: pre_nn, pre_ne, adva

real :: temp_d, tl, len, cp(2), df, ang, cri, d
integer :: i, j, k, status, pn, ln, nn, bn, anes, dn, whole_i, temp_bs, ne, count, nl, num
integer :: ins_num, first_nn, choice, enn, object, on, op_num, on_i, bound, mn, temp_cn
integer :: up, dp, thick, inout, sp, lp, max_iter, iter
integer :: this(3)
logical :: ins, vect, check, iter_check
character(len=5) :: keyword
character(len=21) :: point_type

!==============================================================================
max_iter = 10
allocate ( p(2,bs), line(bs), div_b(2,bs*500), edge(2,pn) )
p = ori_div_b
ori_d = ori_d * 0.9
line = 1
edge(1,:) = ori_d
edge(2,:) = 1.0

do iter = 1, max_iter
    write (*,'(I3,A)') iter, " iteration."
    iter_check = .TRUE.
    
    d = ori_d * (1 - (iter)*0.05)
    if (iter == 1) then
        bn = bs
        div_b(:,1:bs) = p(:,1:bs)
    else
        call divide_point(pn, p, d, line, bn, div_b, edge)
    endif
    
    nn = bn
    anes = nn * 100;  bn = nn * 10
    ne = 0;  first_nn = nn;  enn = nn
    allocate( dbs(nn), db(nn,bn) )
    allocate( node(2,anes), ele(5,anes) )

    node = 0.0;  ele = 0

    filename = 'check_sep_bound00.plt'
    write ( filename(16:17), '(I2)' ) crn
    open (Unit=30, File=filename, STATUS='replace', ACTION='write', IOSTAT=status)
    node(:,1:nn) = div_b(:,1:nn)
    do i = 1, nn
        write (30,*) node(:,i)
    enddo
    close(30)

    db = 0;  dbs = 0
    dbs(1) = nn
    do i = 1, nn
        db(1,i) = i
    enddo

    !d = temp_d
    call partition(nn, anes, node, bn, db, dbs, dn, first_nn )
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
	    call mesh_gen(temp_bs, b_num, b, df, anes, nn, node, ne, ele, whole_i, iter_check)
        deallocate ( b )
        if ( iter_check == .FALSE. ) goto 100
    enddo

    pre_nn = first_nn
    pre_ne = 1

    ! 2개의 Element로 결합된 절점의 소멸
    call combines_two_ele(anes, pre_nn, nn, ne, node, ele, iter_check)
    if ( iter_check == .FALSE. ) then
        write (*,*) 'mesh_gen: combines_two_ele - false'
        goto 100
    endif
    
    ! 6개 이상 모인 절점의 분할 
    call divide_over_ele(anes, nn, ne, node, ele, pre_nn, pre_ne, iter_check)
    if ( iter_check == .FALSE. ) then
        write (*,*) 'mesh_gen: divide_over_ele - false'
        goto 100
    endif
    
    ! 2개의 Element로 결합된 절점의 소멸
    call combines_two_ele(anes, pre_nn, nn, ne, node, ele, iter_check)
    if ( iter_check == .FALSE. ) then
        write (*,*) 'mesh_gen: combines_two_ele - false'
        goto 100
    endif
    
    ! 180각도 수정
    call check_ele_angle(anes, ne, node, ele, first_nn, pre_ne)

    ! Soften inner nodes
    call inner_soft_object(anes, pre_nn, pre_ne, nn, ne, node, ele, iter_check)
    if ( iter_check == .FALSE. ) then
        write (*,*) 'mesh_gen: inner_soft_object - false'
        goto 100
    endif
    
    100 continue
    if ( iter_check ) then
        exit
    else
        !call write_tecplot( nn, node(:,1:nn), ne, ele(2:5,1:ne), rn, ng, 5 )
        deallocate ( dbs, db, node, ele, b_num, dd )
        write (*,*) 'Fail generating mesh at ', iter, ' iteration!!'
    endif
enddo
deallocate ( p, line, div_b, edge )

! Write the mesh data into mdomain
call sep_write_tecplot( nn, node(:,1:nn), ne, ele(:,1:ne), crn )

allocate (surface_set(crn)%adj(pre_nn))
surface_set(crn)%adjn = pre_nn
do i = 1, pre_nn
    surface_set(crn)%adj(i) = i + tol_nn
enddo

do i = 1, nn
    !write (node_unit_num,'(I5,A,F14.7,A,F14.7)') i+tol_nn, ', ', node(1,i), ', ', node(2,i)
    new_node(:,tol_nn+i) = node(:,i)
enddo

do i = 1, ne
    !write (ele_unit_num,'(I5,A,I5,A,I5,A,I5,A,I5)') i+tol_ne,', ',ele(2,i)+tol_nn,', ',ele(3,i)+tol_nn,', ',ele(4,i)+tol_nn,', ',ele(5,i)+tol_nn
    new_ele(:,tol_ne+i) = ele(2:5,i)+tol_nn
enddo

tol_nn = tol_nn + nn
tol_ne = tol_ne + ne

deallocate (db, dbs, b_num, dd)
deallocate (node, ele) 

end subroutine sep_mesh_gen
!=================================================================================
subroutine sep_write_tecplot( nn, node, ne, ele, rn ) 

implicit none 

integer, intent(in) :: nn , ne, ele(5,ne), rn
real, intent(in) :: node(2,nn)

integer :: status, i
character(len=30) :: fn

fn = "2D_sep_remesh00.plt"
write ( fn(14:15), '(I2.2)' ) rn
open (Unit=20, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)

Write (20,'(A,A,A)') 'TITLE="' ,fn(10:20) ,'"'
Write (20,*) 'VARIABLES="x", "y"'
Write (20,'(A,I5,A,I5, A)') 'ZONE N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'

do i = 1, nn
	write(20,*) node(:,i)
enddo
do i = 1, ne
    write (20,* ) ele(2:5,i)
enddo
close(unit=20)

end subroutine sep_write_tecplot