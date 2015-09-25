subroutine post_build_interface(current_step)

use physical_domain
use coarse_domain
use surface_domain

implicit none
integer, intent(in) :: current_step

integer, allocatable :: adj(:), case_adj(:)
real, allocatable :: prop_coord(:,:), case_coord(:,:)
real, allocatable :: contact_pt_to(:,:,:), contact_pt_from(:,:,:)

! contact_pt_from : propellent
! contact_pt_to : case

integer :: i, j, k, adjn, case_nn, prop_nn, case_adjn, count, cn, target_node, seq, temp_seq, rn, case_cont_num
integer :: id, reg_from, reg_to, contact_type, num_contact_group, sp, set_number, contact_dim, num_points
integer :: temp(3)
real :: tol, cont_tol, s, cp(2), cont_coord(4,2), temp_coord(2)
character(len=5) :: keyword
logical :: near_check
logical, allocatable :: check(:)
! backup file data
integer :: cont_unit_num, status
character(len=40) :: bak_cont_file
logical :: flag

flag = .false.
if (mod(current_step,abs(output_interval2)) == 0) then
    flag = .true.
    bak_cont_file = './restart/solid0000000.inp.cont'
    write (bak_cont_file(16:22), '(I7.7)') current_step
    cont_unit_num = 5102
    open(unit=cont_unit_num, file=bak_cont_file , STATUS='replace', IOSTAT=status)
endif

!call free_coarse_domain
!call create_contact(inte_set%cont_num)

tol = mesh_set(2)%ave_d(1)*10.0
cont_coord = 0.0
prop_nn = region(1)%num_nodes
allocate (prop_coord(2,prop_nn))
prop_coord = region(1)%node_coord
case_nn = region(2)%num_nodes
call set_skin(2, 1, sp, case_adjn)
allocate ( case_coord(2,case_nn), case_adj(case_adjn))
case_coord = region(2)%node_coord
case_adj = region(2)%skin_nodes(sp:sp+case_adjn-1) 

rewind (unit_num_cont)
do cn = 1, inte_set%cont_num
    read(unit_num_cont,*) keyword, id, reg_from, reg_to, contact_type, contact_dim, num_contact_group
    if (contact_dim == 3) then
        num_points = 4
    else
        num_points = 2
    endif
    if (flag) write(cont_unit_num,'(A,5(I2,A),I2)') '*cont, ', id, ', ', reg_from, ', ', reg_to, ', ', contact_type, ', ', contact_dim, ', ', num_contact_group
    set_number = num_disp_sets + 2*id - 1
    call free_constr(set_number)
    set_number = num_disp_sets + 2*id
    call free_constr(set_number)
    call free_skin_contact_nodes(reg_to)
    call free_skin_contact_nodes(reg_from)
    allocate(contact_pt_from(2,2,num_contact_group) )
    allocate(contact_pt_to(2,2,num_contact_group) )
    
    if ( contact_type == 2 ) then
        if (reg_from == 2) then
            do i = 1, num_contact_group
                contact_pt_from(:,1,i) = region(reg_from)%node_coord(:,inte_set%cont_slp(1,i))
                contact_pt_from(:,2,i) = region(reg_from)%node_coord(:,inte_set%cont_slp(2,i))
                contact_pt_to(:,1,i) = region(reg_to)%node_coord(:,inte_set%cont_slp(3,i))
                contact_pt_to(:,2,i) = region(reg_to)%node_coord(:,inte_set%cont_slp(4,i))
                if (flag) then
                    write(cont_unit_num,*) contact_pt_from(:,1,i), contact_pt_from(:,2,i)
                    write(cont_unit_num,*) contact_pt_to(:,1,i), contact_pt_to(:,2,i)
                endif
            enddo
        elseif (reg_from == 1) then
            do i = 1, num_contact_group
                temp_coord = region(1)%node_coord(:,inte_set%cont_slp(1,i))  ! start point
                target_node = inte_set%cont_slp(4,i)
                !write (*,*) 'target_node:', target_node
                call find_adj_seq(target_node, case_adjn, case_adj, seq)
                temp_seq = seq
                do j = 1, 3
                    temp(1:2) = (/ seq, seq-1 /)
                    if (seq == 1) temp(2) = case_adjn
                    seq = temp(2)
                    temp(1:2) = case_adj(temp(1:2))
                
                    call near_point(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, cp, tol, near_check)
                    if (near_check) then
                        call projection_s(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, s)
                        !write (*,'(A,3F12.7)') 's: ', s, case_coord(:,temp(2))
                        if (s >= 0.00 .AND. s <= 1.00) then
                            cont_coord(3:4,1) = temp_coord
                            cont_coord(3:4,2) = case_coord(:,temp(1))
                            inte_set%cont_slp(4,i) = temp(1)
                            exit
                        endif
                    endif
                enddo
            
                if (near_check == .false.) then
                    write (*,*) i, ' start check point is out of 3 edges'
                    if (temp_seq == case_adjn) then
                        temp(1:2) = (/ 1, temp_seq /)
                    else
                        temp(1:2) = (/ temp_seq+1, temp_seq /)
                    endif
                    temp(1:2) = case_adj(temp(1:2))
                    call near_point(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, cp, tol, near_check)
                    if (near_check) then
                        call projection_s(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, s)
                        !write (*,'(A,3F12.7)') 's: ', s, case_coord(:,temp(2))
                        if (s >= 0.00 .AND. s <= 1.00) then
                            cont_coord(3:4,1) = temp_coord
                            cont_coord(3:4,2) = case_coord(:,temp(1))
                            inte_set%cont_slp(4,i) = temp(1)
                        endif
                    else
                        write (*,'(A,I,I,L, 3F12.7)') 'j, check, s: ', i, j, near_check, s, case_coord(:,temp_seq)
                        stop 'post_build_interface: cannot find contact target node from start ball node!!'
                    endif
                endif
                    
                temp_coord = region(1)%node_coord(:,inte_set%cont_slp(2,i))  ! last point
                target_node = inte_set%cont_slp(3,i)
                call find_adj_seq(target_node, case_adjn, case_adj, seq)
                temp_seq = seq
                do j = 1, 3
                    temp(1:2) = (/ seq, seq+1 /)
                    if (seq == case_adjn) temp(2) = 1
                    seq = temp(2)
                    temp(1:2) = case_adj(temp(1:2))
                
                    call near_point(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, cp, tol, near_check)
                    if (near_check) then
                        call projection_s(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, s)
                        !write (*,'(A,3F12.7)') 's: ', s, case_coord(:,temp(2))
                        if (s >= 0.00 .AND. s <= 1.00) then
                            cont_coord(1:2,1) = temp_coord
                            cont_coord(1:2,2) = case_coord(:,temp(1))
                            inte_set%cont_slp(3,i) = temp(1)
                            exit
                        endif
                    endif
                enddo        
            
                if (near_check == .false.) then
                    write (*,*) i, ' end check point is out of 3 edges'
                    if (temp_seq == 1) then
                        temp(1:2) = (/ case_adjn, temp_seq /)
                    else
                        temp(1:2) = (/ temp_seq-1, temp_seq /)
                    endif
                    temp(1:2) = case_adj(temp(1:2))
                    call near_point(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, cp, tol, near_check)
                    if (near_check) then
                        call projection_s(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, s)
                        if (s >= 0.00 .AND. s <= 1.00) then
                            !write (*,'(A,3F12.7)') 's: ', s, case_coord(:,temp(2))
                            cont_coord(1:2,1) = temp_coord
                            cont_coord(1:2,2) = case_coord(:,temp(1))
                            inte_set%cont_slp(3,i) = temp(1)
                        endif
                    else
                        write (*,'(A,I,I,L, 3F12.7)') 'j, check, s: ', i, j, near_check, s, case_coord(:,temp_seq)
                        stop 'post_build_interface: cannot find contact target node from end ball node!!'
                    endif
                endif
            
                contact_pt_from(:,1,i) = cont_coord(3:4,1)
                contact_pt_from(:,2,i) = cont_coord(1:2,1)
                contact_pt_to(:,1,i) = cont_coord(1:2,2)
                contact_pt_to(:,2,i) = cont_coord(3:4,2)
                if (flag) then
                    write(cont_unit_num,*) contact_pt_from(:,1,i), contact_pt_from(:,2,i)
                    write(cont_unit_num,*) contact_pt_to(:,1,i), contact_pt_to(:,2,i)
                endif
                !deallocate ( adj, check )
            enddo
        endif
        !call reset_contact(id, contact_type)
        do i = 1, num_contact_group
            write (*,'(A,I3,A)') '---------------', i, ' Group contact points --------------'
            write (*,'(4(F12.7,2X))') contact_pt_from(:,1,i), contact_pt_from(:,2,i)
            write (*,'(4(F12.7,2X))') contact_pt_to(:,1,i), contact_pt_to(:,2,i)
        enddo
        call write_contact(2, id, reg_from, reg_to, contact_type, contact_pt_from, contact_pt_to, contact_dim, num_points, num_contact_group)
    else
        do i = 1, num_contact_group
            read(unit_num_cont,*) cont_coord(:,1)
            read(unit_num_cont,*) cont_coord(:,2)
            contact_pt_from(:,1,i) = cont_coord(1:2,1)
            contact_pt_from(:,2,i) = cont_coord(3:4,1)
            contact_pt_to(:,1,i) = cont_coord(1:2,2)
            contact_pt_to(:,2,i) = cont_coord(3:4,2)
            if (flag) then
                write(cont_unit_num,*) cont_coord(:,1)
                write(cont_unit_num,*) cont_coord(:,2)
            endif
        enddo
        call write_contact(2, id, reg_from, reg_to, contact_type, contact_pt_from, contact_pt_to, contact_dim, num_points, num_contact_group)
    endif
enddo
rewind (unit_num_cont)
if (flag) close(cont_unit_num)

deallocate (prop_coord, case_coord, case_adj)

end subroutine post_build_interface

!=========================================================================================
