subroutine set_fsi_inte(tol_rn, flag, option, current_step)

!option: 0 - initialize
!        1 - reset(n_remesh = 0)
!        2 - reset(n_remesh = 1)

use surface_domain
implicit none 

integer, intent(in) :: tol_rn, flag, option, current_step

integer :: rn, sp, ln, lp, fn, ng, seq, temp_adj, max_adjn
integer :: i, j, num, temp_j, max_num, temp_group, next_adj, cp, tadjn
integer :: temp(2)
real :: temp_node(2,3), dis(3), tol, min_dis, temp_dis
logical :: check
character(len=50) :: filename

integer :: case_adjn
integer, allocatable :: temp_inte(:,:), case_adj(:), inte_info(:)
integer, allocatable :: except(:), prop_sp(:), prop_adj(:,:), prop_adjn(:), tadj(:)
real, allocatable :: node(:,:,:), inte_coord(:,:)
!----------------------------------------------------------------------------
rn = fsi_inte_set%fsi_inte_point(1)
sp = fsi_inte_set%fsi_inte_point(2)
ln = fsi_inte_set%fsi_inte_point(3)
lp = fsi_inte_set%fsi_inte_point(4)

if (tol_rn >= 2) then
    case_adjn = region(2)%num_skin
    allocate ( case_adj(case_adjn) )
    case_adj = region(2)%skin_nodes
    tol = fsi_inte_set%tol
endif

ng = region(1)%num_groups
allocate ( prop_adjn(ng), prop_sp(ng), except(ng) )
do i = 1, ng
    call set_skin(1, i, prop_sp(i), prop_adjn(i))
enddo
max_adjn = maxval(prop_adjn(:))
allocate ( prop_adj(ng, max_adjn) )
do i = 1, ng
    prop_adj(i,1:prop_adjn(i)) = region(1)%skin_nodes(prop_sp(i):prop_sp(i)+prop_adjn(i)-1)
enddo

except = 0
num = 0
max_num = 0
do i = 1, tol_rn
    num = num + region(i)%num_nodes
    if ( max_num < region(i)%num_nodes ) max_num = region(i)%num_nodes
enddo
allocate ( temp_inte(2,num), node(tol_rn,2,max_num) )
temp_inte = 0

node = 0.0
do i = 1, tol_rn
    if (i == 1 .AND. flag == 2) then
        node(i,:,1:mesh_set(i)%nn) = mesh_set(i)%coord(:,1:mesh_set(i)%nn)
    else
        node(i,:,1:region(i)%num_nodes) = region(i)%node_coord(1:2,1:region(i)%num_nodes)
    endif
enddo
temp_group = 0
if (rn == 1) then
    fn = 0
    seq = 0
    do i = 1, ng
        tadjn = prop_adjn(i)
        allocate (tadj(tadjn))
        tadj = prop_adj(i,1:tadjn)
        call find_adj_seq(sp, tadjn, tadj, seq)
        deallocate (tadj)
        if (seq /= 0) then
            temp_group = i
            exit
        endif
    enddo
elseif (rn == 2) then
    call find_adj_seq(sp, case_adjn, case_adj, seq)
    fn = 1
    temp_inte(:,fn) = (/ rn, sp /)
endif
temp_adj = seq

do
    if ( rn == 2 ) then  !<--- case(domain 2)
        temp_node(:,1) = node(rn,:,case_adj(temp_adj))
        if ( temp_adj == case_adjn ) then
            next_adj = 1
        else
            next_adj = temp_adj + 1
        endif
        temp_node(:,2) = node(rn,:,case_adj(next_adj))
    
        call calc_len( temp_node(:,1), temp_node(:,2), dis(1) )
        check = .TRUE.
        do i = 1, ng
            if ( except(i) == 0 ) then
                min_dis = 10e8
                temp_j = 0
                do j = prop_adjn(i), 1, -1 
                    temp_node(:,3) = node(1,:,prop_adj(i,j))
                    call calc_len( temp_node(:,1), temp_node(:,3), dis(2) )
                    call calc_len( temp_node(:,2), temp_node(:,3), dis(3) )
                    if ( abs(dis(1)-dis(2)-dis(3)) < tol/100. ) then
                        if (min_dis > dis(2)) then
                            temp_j = j
                            temp_dis = dis(3)
                            min_dis = dis(2)
                        endif
                    endif
                enddo
                !inte_set%cont_slp(j,i)
                ! j: 1 = 시작점, 2 = 끝점
                ! i: 그룹 번호
                !sp = inte_set%cont_slp(j,2)       
                !call find_adj_seq(sp, prop_adjn, prop_adj, seq)
                !temp_node(:,3) = node(1,:,prop_adj(2,i))
                !call calc_len( temp_node(:,1), temp_node(:,3), dis(2) )
                !call calc_len( temp_node(:,2), temp_node(:,3), dis(3) )
                !if ( abs(dis(1)-dis(2)-dis(3)) < tol/1000. ) then
                !    if (min_dis > dis(2)) then
                !        temp_j = j
                !        temp_dis = dis(3)
                !        min_dis = dis(2)
                !    endif
                !endif
                
                if (temp_j /= 0 ) then 
!                    if ( temp_dis < tol/1000. ) then ! last point at case(temp_node(2)) correspond to point of propellant
!                        fn = fn + 1
!                        temp_inte(:,fn) = (/ rn, case_adj(next_adj) /)
!                    endif
                    fn = fn + 1
                    rn = 1
                    temp_inte(:,fn) = (/ rn, prop_adj(i,temp_j) /)
                    if (temp_j == prop_adjn(i)) then
                        temp_adj = 1
                    else
                        temp_adj = temp_j + 1
                    endif
                    temp_group = i
                    check = .FALSE.
                    exit
                endif
            endif
        enddo
        if ( check ) then
            fn = fn + 1
            temp_inte(:,fn) = (/ rn, case_adj(next_adj) /)
            temp_adj = next_adj
            if ( rn == ln .AND. case_adj(temp_adj) == lp ) exit  ! end of main loof
        endif
    elseif ( rn == 1 ) then   !<--- prop(domain 1)
        temp_node(:,1) = node(rn,:,prop_adj(temp_group,temp_adj))
        check = .TRUE.
        if (tol_rn >= 2) then
            do i = case_adjn, 1, -1
                temp_node(:,2) = node(2,:,case_adj(i))
                if ( i == 1 ) then
                    temp_node(:,3) = node(2,:,case_adj(case_adjn))
                else 
                    temp_node(:,3) = node(2,:,case_adj(i-1))
                endif
                call calc_len( temp_node(:,2), temp_node(:,3), dis(1) )
                call calc_len( temp_node(:,1), temp_node(:,2), dis(2) )
                call calc_len( temp_node(:,1), temp_node(:,3), dis(3) )
                !if ( abs(dis(1)-dis(2)-dis(3)) < tol/1000. ) then
                if ( abs(dis(1)-dis(2)-dis(3)) < dis(1)/100. ) then
                    fn = fn + 1
                    temp_inte(:,fn) = (/ rn, prop_adj(temp_group,temp_adj) /)
                    rn = 2
                    if (i == 1) then
                        temp_adj = case_adjn
                    else
                        temp_adj = i - 1
                    endif
                    check = .FALSE.
    !                if ( dis(2) < tol/1000. ) then ! start point at case(temp_node(1)) correspond to point of propellant
    !                    fn = fn + 1
    !                    temp_inte(:,fn) = (/ rn, case_adj(temp_adj) /)
    !                endif
                    except(temp_group) = 1
                    exit
                endif
            enddo
        endif
        if ( check ) then
             fn = fn + 1
             temp_inte(:,fn) = (/ rn, prop_adj(temp_group,temp_adj) /)
             if ( rn == ln .AND. prop_adj(temp_group,temp_adj) == lp ) exit  ! end of main roof
             if ( temp_adj == prop_adjn(temp_group) ) then
                temp_adj = 1
             else
                temp_adj = temp_adj + 1
             endif
        endif                
    endif
enddo

allocate ( inte_coord(2,fn) )
do i = 1, fn
    rn = temp_inte(1,i)
    cp = temp_inte(2,i)

    !inte_coord(:,i) = region(rn)%node_coord(1:2,cp)
    inte_coord(:,i) = node(rn,:,cp)
enddo
!if (option == 1) then
!    filename = "./output/solid/interface/inte_coord000000.plt"
!    write ( filename(36:41), '(I6.6)' ) current_step
!    open (Unit=30, File=filename, STATUS='replace', ACTION='write')
!    do i = 1, fn
!        write (30,*) inte_coord(:,i)
!    enddo
!    close(30)
!endif

if ( option /= 0 ) then
    if (fsi_inte_set%fsi_ni /= fn) then
        if (allocated(fsi_inte_set%old_fsi_inte)) deallocate(fsi_inte_set%old_fsi_inte)
        allocate (fsi_inte_set%old_fsi_inte(2,fsi_inte_set%fsi_ni))
        fsi_inte_set%old_fsi_inte = fsi_inte_set%fsi_inte_node
        do i = 1, fn
            check = .true.
            do j = 1, fsi_inte_set%fsi_ni
                temp = fsi_inte_set%fsi_inte_node(:,j)
                if (temp_inte(1,i) == temp(1) .and. temp_inte(2,i) == temp(2)) then
                    check = .false.
                    exit
                endif
            enddo
            if (check) temp_inte(1,i) = 3
        enddo
    endif    
    deallocate (fsi_inte_set%fsi_inte_node, fsi_inte_set%B_flag, fsi_inte_set%ni_flag)
    deallocate (fsi_inte_set%B_rate, fsi_inte_set%fsi_coord, fsi_inte_set%new_fsi_node)
endif
allocate (fsi_inte_set%fsi_inte_node(2,fn), fsi_inte_set%ni_flag(fn))
allocate (fsi_inte_set%fsi_coord(2,fn), fsi_inte_set%new_fsi_node(fn))
allocate (fsi_inte_set%B_rate(fn-1), fsi_inte_set%B_flag(fn-1))
fsi_inte_set%old_fsi_ni = fsi_inte_set%fsi_ni
fsi_inte_set%fsi_ni = fn
fsi_inte_set%fsi_inte_node = temp_inte(:,1:fn)
fsi_inte_set%fsi_coord = inte_coord

fsi_node_change = .FALSE.
if (option == 0) then
    fsi_inte_set%new_fsi_node = .FALSE.
else
    if (fsi_inte_set%old_fsi_ni /= fsi_inte_set%fsi_ni) then
        fsi_node_change = .TRUE.
        call check_new_fsi_node
    else 
        fsi_inte_set%new_fsi_node = .FALSE.
    endif
endif

if ( flag == 0 ) then
    allocate (inte_info(fn))
    inte_info = temp_inte(1,1:fn)
    call fsi_load_set(1, fn, inte_info)
    deallocate (inte_info)
endif

deallocate ( temp_inte, node )
deallocate(inte_coord, prop_adj, prop_adjn, prop_sp, except)
if ( tol_rn >= 2 ) deallocate(case_adj)

end subroutine set_fsi_inte