subroutine move_surface(num_rn, crn, nn, node, output, mp)

use surface_domain

implicit none

integer, intent(in) :: num_rn, crn, nn
integer, intent(inout) :: output, mp(2)
real, intent(inout) :: node(2,nn)

integer :: i, j, k, cp, sp, mn, ng, rn, count
integer :: adjn, seq, cg, pos
integer :: this(3), fb(2)
real :: off, ang, dis, vel, tol_vel, tolerance
real :: vec(2,2), temp_coord(2,4), np(2)
integer, allocatable :: adj(:), move_nodes(:), temp_seq(:)
real, allocatable :: temp_node(:,:)
logical :: check, cont_check, change

allocate ( temp_node(2,nn), move_nodes(nn) )
temp_node = node
change = .FALSE.
tol_vel = 0.0
mn = 0
ng = region(crn)%num_groups
if (num_rn >= 2) then
    if (ng == 1) then
        !tolerance = 0.4*mesh_set(2)%ave_d(1)
        tolerance = 0.2*mesh_set(crn)%ave_d(1) 
    else
        !tolerance = 0.15*mesh_set(2)%ave_d(1)
        tolerance = 0.15*mesh_set(crn)%ave_d(1) 
    endif
else
    tolerance = 0.5*mesh_set(crn)%ave_d(1) 
endif

do i = 1, fsi_inte_set%fsi_ni
    if ( fsi_inte_set%B_flag(i) == 1 ) then
        
        ! move boundary
        vel = fsi_inte_set%B_rate(i)
        tol_vel = tol_vel + vel
        rn = fsi_inte_set%fsi_inte_node(1,i)
        cp = fsi_inte_set%fsi_inte_node(2,i)
        if (rn == 1) then
            mn = mn + 1
            move_nodes(mn) = cp
        endif
        
        if (vel > mesh_set(rn)%min_d(1)*10e-9) then
            call find_group(rn, cp, cg)
            call set_skin(rn, cg, sp, adjn)
            allocate( adj(adjn) )
            adj = region(rn)%skin_nodes(sp:sp+adjn-1) 
            call find_adj_seq(cp, adjn, adj, seq)
        
            if ( seq == 1 ) then
                this(:) = (/ adj(adjn), adj(seq), adj(seq+1) /)
            elseif ( seq == adjn ) then
                this(:) = (/ adj(seq-1), adj(seq), adj(1) /)
            else
                this(:) = (/ adj(seq-1), adj(seq), adj(seq+1) /)
            endif
            call find_pos(rn, this, pos)
            deallocate ( adj )
        
            call calc_angle(node(:,this(1)), node(:,this(2)), node(:,this(3)), ang)
            vec(:,1) = node(:,this(1)) - node(:,this(2))
            vec(:,2) = node(:,this(3)) - node(:,this(2))

            !write (*,*) temp_node(this(2),:)
            if ( pos == 1 ) then ! <- start interface node
                if ( abs(ang-180.) > 0.1 ) then
                    call calc_len(node(:,this(1)), node(:,this(2)), dis)
                    off = vel/sin(ang*pi/180.)
                    vec(:,1) = (vec(:,1)*off)/dis
                    temp_node(:,this(2)) = node(:,this(2)) + vec(:,1)
                    mp(1) = this(2)
                else
                    write (*,*) "Error : Angle at SP is 180 deg!!(", this(:), ")"
                    stop
                endif
            elseif ( pos == 3 ) then ! <- Last interface node
                if ( abs(ang-180.) > 0.1 ) then
                    call calc_len(node(:,this(3)), node(:,this(2)), dis)
                    off = vel/sin(ang*pi/180.)
                    vec(:,2) = (vec(:,2)*off)/dis
                    temp_node(:,this(2)) = node(:,this(2)) + vec(:,2)
                    mp(2) = this(2)
                else
                    write (*,*) "Error : Angle at LP is 180 deg!!(", this(:), ")"
                    stop
                endif
            else
                call calc_len(node(:,this(3)), node(:,this(2)), dis)
                off = vel/sin(ang*pi/360)
                !if ( ang > 190. ) off = off * (1 - 0.5*(sin((ang-180.)*pi/360.)))
                vec(:,2) = (vec(:,2)*off)/dis
                temp_node(1,this(2)) = node(1,this(2))+(cos(ang*pi/360.)*vec(1,2)-sin(ang*pi/360.)*vec(2,2))
                temp_node(2,this(2)) = node(2,this(2))+(sin(ang*pi/360.)*vec(1,2)+cos(ang*pi/360.)*vec(2,2))
            endif
        endif
    endif
enddo
node = temp_node

if (tol_vel > mesh_set(1)%min_d(1)*10e-8) then
    move_flag = 1
    do i = 1, ng
        cont_check = .FALSE.
        call set_skin(crn, i, sp, adjn)
        allocate( adj(adjn), temp_seq(adjn) )
        adj = region(crn)%skin_nodes(sp:sp+adjn-1) 

        do j = adjn, 1, -1
            call move_point_check(mn, move_nodes(1:mn), adj(j), check)
            
            if (check) then
                temp_coord(:,3) = node(:,adj(j))
                count = 0
                fb(1) = adjn
                if (j-8 < 1) fb(1) = adjn+j-8
                do k = j+7, fb(1)
                    count = count + 1
                    temp_seq(count) = k
                enddo
                fb(2) = 1
                if (j+7 > adjn) fb(2) = j+7-adjn
                do k = fb(2), j-8
                    count = count + 1
                    temp_seq(count) = k
                enddo

                do k = 1, count
                    cp = temp_seq(k)
                    if ( cp == adjn ) then
                        this(1:2) = (/ cp, 1 /)
                    else
                        this(1:2) = (/ cp, cp+1 /)
                    endif
                    !if ( (this(1)>j+6 .OR. this(1)<j-6) .AND. (this(2)>j+6 .OR. this(2)<j-6) ) then
                        temp_coord(:,1) = node(:,adj(this(1)))
                        temp_coord(:,2) = node(:,adj(this(2)))
                        
                        call near_point(temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), np, tolerance, cont_check)
                        !if ( ng == 1 ) then
                        !    call near_point(temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), np, 0.35*mesh_set(crn)%ave_d(i), cont_check)
                        !else
                        !    call near_point(temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), np, 0.2*mesh_set(crn)%ave_d(i), cont_check)
                        !endif
                        if ( cont_check ) then
                            call calc_len( temp_coord(:,3), np, dis)
                            sep_seq(:) = (/ this(1), this(2), j /)
                            sep_num(:) = (/ adj(this(1)), adj(this(2)), adj(j) /)
                            sep_point = np
                            output = 1
                            write (*,*) 'dis, tol_dis:', dis, 0.35*mesh_set(crn)%ave_d(i)
                            write (*,*) 'np:', np
                            write (*,*) 'sep_seq:', sep_seq
                            write (*,*) 'sep_num:', sep_num
                            exit
                        endif
                    !endif
                enddo
            endif
            if ( cont_check ) exit
        enddo
        deallocate ( adj, temp_seq ) 
        if ( cont_check ) exit
    enddo
else
    move_flag = 2
endif

deallocate ( temp_node )

end subroutine move_surface
!==========================================================================
subroutine find_group(rn, cp, cg)

use physical_domain

implicit none

integer, intent(in) :: rn, cp
integer, intent(inout) :: cg

integer :: i, ngroup

ngroup = region(rn)%num_groups

cg = 0
do i = 1, ngroup
    if (i == ngroup) then
        if ( region(rn)%node_group(i) <= cp .AND. cp <= region(rn)%num_nodes ) cg = i
    else
        if ( region(rn)%node_group(i) <= cp .AND. cp <= region(rn)%node_group(i+1)-1 ) cg = i
    endif
enddo

if (cg == 0) stop 'find_group(set_skin): Cannot find the group of cp!!!'

end subroutine find_group

!==========================================================================
subroutine find_pos(rn, this, pos)

use fsi_domain
!pos = 1 : start point
!      2 : mid point
!      3 : last point
implicit none

integer, intent(in) :: rn, this(3)
integer, intent(out) :: pos

integer :: i, ni
logical :: check(2)

check = .FALSE.
ni = fsi_inte_set%fsi_ni
do i = 1, ni
    if ( fsi_inte_set%fsi_inte_node(1,i) == rn ) then
        if ( fsi_inte_set%fsi_inte_node(2,i) == this(1) ) then
            check(1) = .TRUE.
        elseif ( fsi_inte_set%fsi_inte_node(2,i) == this(3) ) then
            check(2) = .TRUE.
        endif
        if ( check(1) .AND. check(2) ) exit
    endif
enddo

if ( check(1) .AND. check(2) == .FALSE. ) then
    pos = 3
elseif (check(1) == .FALSE. .AND. check(2) ) then
    pos = 1
else
    pos = 2
endif

end subroutine find_pos

subroutine move_point_check(mn, move_nodes, adj, check)

implicit none

integer, intent(in) :: mn, move_nodes(mn), adj
logical, intent(out) :: check

integer :: i

check = .false.
do i = 1, mn
    if (move_nodes(i) == adj) then
        check = .true.
        exit
    endif
enddo

end subroutine move_point_check


