subroutine Press2Load(Dim, press, PROPEL_POINT_NUM, PROPEL_NBFACE, PROPEL_POINT, PROPEL_F2N, PROPEL_PLOAD, &
                                  CASE_POINT_NUM, CASE_NBFACE, CASE_POINT, CASE_F2N, CASE_PLOAD)
!use fsi_domain
implicit none

integer, intent(in) :: Dim, PROPEL_POINT_NUM, CASE_POINT_NUM, PROPEL_NBFACE, CASE_NBFACE
integer, intent(in) :: PROPEL_F2N((Dim-1)*2,PROPEL_NBFACE), CASE_F2N((Dim-1)*2,CASE_NBFACE)
real(8), intent(   in) :: press, PROPEL_POINT(Dim,PROPEL_POINT_NUM), CASE_POINT(Dim,CASE_POINT_NUM)
real(8), intent(inout) :: PROPEL_PLOAD(Dim,PROPEL_POINT_NUM), CASE_PLOAD(Dim,CASE_POINT_NUM)

integer :: i, j, k, q, count, adjn, num, cp, seq, flag, nBNode, tol_rn, rn, elem_num, face_num, pre_num
integer :: case_adjn, ngroup, max_adjn, inte_rn(2)
real(8) :: dis, temp_node(2,2), temp_load(2), vec(2), temp_vec(2), temp_val
integer, allocatable :: case_adj(:), prop_sp(:), prop_adjn(:), prop_adj(:,:)
integer, allocatable :: inte_load_num(:,:), load_point(:,:), load_seq(:)
real(8), allocatable :: inte_load(:), load_value(:), pload(:)
real(8), parameter :: pi = 3.141592653589793238462643383279502884197169399375
logical :: check

! for 3D
integer :: press_node(PROPEL_POINT_NUM), csn(2), cs_nodes(9,PROPEL_POINT_NUM), end_nodes(9), seq3(2), this(3)
real(8) :: coord_3D(3,4), load_3D(3,4), temp_face(PROPEL_NBFACE), cons_nodes(200)
integer, allocatable :: inte_node(:)
!-------------------------------------------------------------------------------
if (Dim == 2) then
    !tol_rn = num_of_regions
    !nBNode = fsi_inte_set%fsi_ni
    !allocate (inte_load_num(3,nBNode*2), inte_load(nBNode*2), pload(nBNode), load_seq(nBNode))
    !inte_load_num = 0
    !inte_load = 0.0
    !
    !ngroup = region(1)%num_groups
    !allocate (prop_adjn(ngroup), prop_sp(ngroup))
    !do i = 1, ngroup
    !    call set_skin(1, i, prop_sp(i), prop_adjn(i))
    !enddo
    !max_adjn = maxval(prop_adjn(:))
    !allocate (prop_adj(ngroup, max_adjn))
    !do i = 1, ngroup
    !    prop_adj(i,1:prop_adjn(i)) = region(1)%skin_nodes(prop_sp(i):prop_sp(i)+prop_adjn(i)-1)
    !enddo
    !
    !if (tol_rn >= 2) then
    !    call set_skin(2, 1, num, case_adjn)
    !    allocate (case_adj(case_adjn))
    !    case_adj = region(2)%skin_nodes(num:num+case_adjn-1)
    !endif
    !
    !pload = 0.0
    !do i = 1, nBNode
    !    rn = fsi_inte_set%fsi_inte_node(1,i)
    !    cp = fsi_inte_set%fsi_inte_node(2,i)
    !    if (tol_rn < rn) then
    !        fsi_inte_set%fsi_inte_node(1,i) = 2
    !        rn = 2
    !    endif
    !    if (rn == 1) then
    !        num = 0
    !        seq = 0
    !        do j = ngroup, 1, -1
    !            do k = 1, prop_adjn(j)
    !                if ( prop_adj(j,k) == cp ) then
    !                    seq = k
    !                    exit
    !                endif
    !            enddo
    !            if (seq /= 0) then
    !                load_seq(i) = prop_adjn(j)-seq+1+num
    !                pload(i) = press
    !                exit
    !            endif
    !            num = num + prop_adjn(j)
    !        enddo
    !        if (seq == 0) then
    !            write (*,*) 'RN, CP: ', rn, cp
    !            STOP 'Error : find_adj_seq: cannot find sequence of node(cp) in adj!!!'
    !        endif
    !    elseif (rn == 2) then
    !        call find_adj_seq(cp, case_adjn, case_adj, seq)
    !        pload(i) = press
    !        load_seq(i) = CASE_POINT_NUM-seq+1
    !    endif
    !enddo
    !
    !do i = 1, nBNode
    !    temp_load = 0.0
    !    flag = 0
    !    do j = 1, 2
    !        check = .TRUE.
    !        if (j == 1 .AND. i /= 1) then
    !            inte_rn(2) = fsi_inte_set%fsi_inte_node(1,i-1)
    !            cp = fsi_inte_set%fsi_inte_node(2,i-1)
    !            temp_node(:,1) = region(inte_rn(2))%node_coord(:,cp) + history(inte_rn(2))%data(region(inte_rn(2))%node(4:5,cp),1)
    !            inte_rn(1) = fsi_inte_set%fsi_inte_node(1,i)
    !            cp = fsi_inte_set%fsi_inte_node(2,i)
    !            temp_node(:,2) = region(inte_rn(1))%node_coord(:,cp) + history(inte_rn(1))%data(region(inte_rn(1))%node(4:5,cp),1)
    !            temp_val = (pload(i)+(pload(i)+pload(i-1))*0.5)*0.5
    !        elseif (j == 2 .AND. i /= nBNode) then
    !            inte_rn(1) = fsi_inte_set%fsi_inte_node(1,i)
    !            cp = fsi_inte_set%fsi_inte_node(2,i)
    !            temp_node(:,1) = region(inte_rn(1))%node_coord(:,cp) + history(inte_rn(1))%data(region(inte_rn(1))%node(4:5,cp),1)
    !            inte_rn(2) = fsi_inte_set%fsi_inte_node(1,i+1)
    !            cp = fsi_inte_set%fsi_inte_node(2,i+1)
    !            temp_node(:,2) = region(inte_rn(2))%node_coord(:,cp) + history(inte_rn(2))%data(region(inte_rn(2))%node(4:5,cp),1)
    !            temp_val = (pload(i)+(pload(i)+pload(i+1))*0.5)*0.5
    !        else
    !            check = .FALSE.
    !        endif
    !    
    !        if (check) then
    !            call calc_len(temp_node(:,1), temp_node(:,2), dis)
    !            vec(:) = ((temp_node(:,2)-temp_node(:,1))/dis)
    !            temp_vec(1) = -vec(2)
    !            temp_vec(2) = vec(1)
    !            temp_load = temp_load + temp_vec*dis*0.5*temp_val
    !        endif
    !    enddo
    !    if (flag == 1) then
    !        temp_load = temp_load * 2.0
    !    endif
    !
    !    if (fsi_inte_set%fsi_inte_node(1,i) == 1) then
    !        PROPEL_PLOAD(1,load_seq(i)) = temp_load(2)
    !        PROPEL_PLOAD(2,load_seq(i)) = temp_load(1)
    !    else
    !        CASE_PLOAD(1,load_seq(i)) = temp_load(2)
    !        CASE_PLOAD(2,load_seq(i)) = temp_load(1)
    !    endif
    !enddo
    !
    !deallocate (inte_load_num, inte_load)
    !deallocate (prop_adjn, prop_sp, prop_adj)
    !if ( tol_rn >= 2 ) deallocate(case_adj)
elseif (Dim == 3) then
    !do i = 1, fsi_inte_set%fsi_ni
    !    rn = fsi_inte_set%fsi_inte_node(1,i)
    !    elem_num = fsi_inte_set%fsi_inte_node(2,i)
    !    face_num = fsi_inte_set%fsi_inte_node(3,i)
    !    !coord_3D = region(rn)%node_coord(:,region(rn)%bound_face(:,face_num))
    !    call press2load_3D(coord_3D, press, load_3D)
    !enddo
    !temp_face(1:16) = (/ 22, 24, 26, 28, 80, 81, 82, 83, 108, 109, 110, 111, 146, 148, 150, 152 /)
    !do i = 1, 16
    !   coord_3D = PROPEL_POINT(:,PROPEL_F2N(:,temp_face(i)))
    !   call press2load_3D(coord_3D, press, load_3D)
    !   PROPEL_PLOAD(:,PROPEL_F2N(:,temp_face(i))) = PROPEL_PLOAD(:,PROPEL_F2N(:,temp_face(i))) + load_3D
    !enddo
    
    csn = (/ 5, 1 /)
    cs_nodes(1:5,1) = (/ 125, 309, 393, 601, 602 /)
    end_nodes(1:5) = (/ 202, 341, 425, 710, 712 /)
    cs_nodes(1,2) = 127
    !csn = (/ 9, 1 /)
    !cs_nodes(:,1) = (/ 1401, 1867, 2074, 2281, 2488, 2695, 2902, 3109, 3751 /)
    !end_nodes = (/ 1370, 1845, 2052, 2259, 2466, 2673, 2880, 3087, 3721 /)
    !cs_nodes(1,2) = 1536
    pre_num = 0
    find_interface_node: do 
        csn(2) = csn(2) + 1
        ! find first node(cs_nodes(1,csn(2))
        if (csn(2) /= 2) then
            this = (/ cs_nodes(1,csn(2)-1), cs_nodes(2,csn(2)-1), pre_num /)
            do j = 1, PROPEL_NBFACE
                call find_count(2, this(1:2), 4, PROPEL_F2N(:,j), count)
                if (count == 2) then
                    call find_count(3, this, 4, PROPEL_F2N(:,j), count)
                    if (count == 2) then
                        do k = 1, 4
                            if (PROPEL_F2N(k,j) == this(1)) then
                                if (k == 1) then
                                    seq3 = (/ 2, 4 /)
                                elseif (k == 4) then
                                    seq3 = (/ 3, 1 /)
                                else
                                    seq3 = (/ k-1, k+1 /)
                                endif
                                if (PROPEL_F2N(seq3(1),j) == this(2)) then
                                    cs_nodes(1,csn(2)) = PROPEL_F2N(seq3(2),j)
                                else
                                    cs_nodes(1,csn(2)) = PROPEL_F2N(seq3(1),j)
                                endif
                                exit ! k
                            endif
                        enddo
                        exit ! j
                    endif
                endif
            enddo
        endif

        ! find other nodes(cs_nodes(2:csn(1),csn(2))
        do i = 2, csn(1) 
            this = (/ cs_nodes(i-1,csn(2)-1), cs_nodes(i,csn(2)-1), cs_nodes(i-1,csn(2)) /)
            do j = 1, PROPEL_NBFACE
                call find_count(3, this, 4, PROPEL_F2N(:,j), count)
                if (count == 3) then
                    do k = 1, 4
                        num = PROPEL_F2N(k,j)
                        do q = 1, 3
                            if (PROPEL_F2N(k,j) == this(q)) then
                                num = 0
                                exit  ! q
                            endif
                        enddo
                        if (num /= 0) then
                            cs_nodes(i,csn(2)) = num
                            exit ! k
                        endif
                    enddo
                    exit ! j
                endif
            enddo
        enddo
        pre_num = cs_nodes(1,csn(2)-1)
            
        if (cs_nodes(csn(1),csn(2)) == end_nodes(csn(1))) exit find_interface_node
    enddo find_interface_node
    allocate (inte_node(PROPEL_POINT_NUM))
    inte_node = 0
    do i = 1, csn(1)
        do j = 1, csn(2)
            inte_node(cs_nodes(i,j)) = 1
        enddo
    enddo
    
    temp_face = 0
    do i = 1, PROPEL_NBFACE
        num = sum(inte_node(PROPEL_F2N(:,i)))
        if (num == 4) then
            temp_face(i) = 1
        endif
    enddo
    PROPEL_PLOAD = 0.0
    do i = 1, PROPEL_NBFACE
        if (temp_face(i) == 1) then
            coord_3D = PROPEL_POINT(:,PROPEL_F2N(:,i))
            call press2load_3D(coord_3D, press, load_3D)
            PROPEL_PLOAD(:,PROPEL_F2N(:,i)) = PROPEL_PLOAD(:,PROPEL_F2N(:,i)) + load_3D
        endif
    enddo
    do i = 1, csn(2)
        PROPEL_PLOAD(2,cs_nodes(5,i)) = 0.0
    enddo
    deallocate(inte_node)
endif

end subroutine Press2Load


subroutine press2load_3D(coord, press, load)

implicit none 

real(8), intent(in) :: coord(3,4), press
real(8), intent(inout) :: load(3,4)

integer :: i, j
integer :: this(3,2)
real(8) :: pl(4), cp(3), area(4), s, tol_area, dis(3), temp_coord(3,4), tol_load(3)

call calc_plane(coord(:,1), coord(:,2), coord(:,3), pl)
cp = (coord(:,1)+coord(:,2)+coord(:,3)+coord(:,4))*0.25

this(:,1) = (/ 1, 2, 3 /)
this(:,2) = (/ 3, 4, 1 /)

tol_area = 0.0
do i = 1, 2
    call calc_length_3D(coord(:,this(1,i)), coord(:,this(2,i)), dis(1))
    call calc_length_3D(coord(:,this(2,i)), coord(:,this(3,i)), dis(2))
    call calc_length_3D(coord(:,this(3,i)), coord(:,this(1,i)), dis(3))
    s = (dis(1)+dis(2)+dis(3))/2.
    if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
        tol_area = tol_area + sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
    endif
enddo

temp_coord(:,4) = cp
do i = 1, 4
    temp_coord(:,2) = coord(:,i)
    if (i == 1) then
        temp_coord(:,1) = (coord(:,1)+coord(:,4))*0.5
        temp_coord(:,3) = (coord(:,1)+coord(:,2))*0.5
    elseif (i == 4) then
        temp_coord(:,1) = (coord(:,i-1)+coord(:,i))*0.5
        temp_coord(:,3) = (coord(:,i)+coord(:,1))*0.5
    else
        temp_coord(:,1) = (coord(:,i-1)+coord(:,i))*0.5
        temp_coord(:,3) = (coord(:,i)+coord(:,i+1))*0.5
    endif
    
    area(i) = 0.0
    do j = 1, 2
        call calc_length_3D(temp_coord(:,this(1,j)), temp_coord(:,this(2,j)), dis(1))
        call calc_length_3D(temp_coord(:,this(2,j)), temp_coord(:,this(3,j)), dis(2))
        call calc_length_3D(temp_coord(:,this(3,j)), temp_coord(:,this(1,j)), dis(3))
        s = (dis(1)+dis(2)+dis(3))/2.
        if (s*(s-dis(1))*(s-dis(2))*(s-dis(3)) > 0.0) then
            area(i) = area(i) + sqrt(s*(s-dis(1))*(s-dis(2))*(s-dis(3)))
        endif
    enddo
enddo
!write (*,*) 'pl:', pl(1:3)
tol_load = press*tol_area*(-pl(1:3))
do i = 1, 4
    load(:,i) = (area(i)/tol_area)*tol_load
    !write (*,'(I5,A,4(F12.7,2X))') i, ':', load(:,i), area(i)
enddo

end subroutine press2load_3D