module surface_domain 

use fsi_domain

implicit none
save

type mesh_type       ! for ale 
    integer :: cn = 0, nn, ne, convex_cons_num, line_cons_num
    real    :: cri_d
    integer, allocatable :: convex_boundary(:,:), line_boundary(:,:), boundary(:)
    integer ,allocatable :: corn(:), corn_group(:), flag(:), ele(:,:)
    real,    allocatable :: convex_cons(:,:), line_cons(:,:), edge(:,:), corn_coord(:,:)
    real,    allocatable :: ave_d(:), min_d(:), coord(:,:), nndis(:), ele_area(:)
    
    ! temporary for 3D
    integer, allocatable :: sub_num(:,:)
    integer :: bn, bfn, inn
    integer, allocatable :: bound_node(:), surf2ori(:), mat_array(:)
    integer, allocatable :: element(:,:), bound_face(:,:), bound_lv(:)
    integer, allocatable :: inner_node(:), inn_nci_num(:), inn_nci(:,:)
    real,    allocatable :: node_coord(:,:)
end type

type surface_type    ! for remeshing
    integer :: nn = 0, ne = 0, num_nodes = 0, num_elements = 0, adjn = 0
    integer :: part_ne = 0, part_nn = 0, part_ni = 0, part_adjn = 0, part_inte_seq(2)
    logical :: part_remesh, part_first
    integer, allocatable :: adj(:)
    integer ,allocatable :: part_ele(:), part_nodes(:), part_inte_nodes(:), part_adj(:)
    real,    allocatable :: coord(:,:), part_inte_coord(:,:)
    ! for 3D
    integer :: tri_nn, tri_ne, tri_np2, ridge_group_num
    integer, allocatable :: tri_elem(:,:), tri_p2(:)
    integer, allocatable :: ridge_num(:), ridge(:,:), ridge_conn(:)
    real,    allocatable :: tri_coord(:,:)
end type

type inte_type
    integer :: sn = 0, cn = 0, ce = 0
    integer :: cont_num = 0, cont_node(4), abla_node(4), surf_num(2) = 0, cont_region(2) = 0
    real    :: cont_off, move_vel
    integer, allocatable :: cont_slp(:,:), s_nodes(:), c_nodes(:), c_elem(:)
    integer, allocatable :: surf_nodes(:,:), surf_nodes_info(:,:)
    integer, allocatable :: start_nodes(:), end_nodes(:)
end type

type move_type
    ! for 3D
    integer :: p2_csn = 0, p2cn = 0, match_bfn = 0
    integer, allocatable :: csn(:,:), cs_nodes1(:,:,:), cs_nodes2(:,:,:), prop2case(:,:), p2c_nn(:)
    integer, allocatable :: match_b_face(:,:)
    real, allocatable :: seep_dis(:,:)
end type

type revol_type
    ! for 3D
    integer :: revol_info(3)
    integer, allocatable :: revol_nodes(:,:)
end type

type(inte_type) :: inte_set
type(mesh_type), allocatable :: mesh_set(:)
type(surface_type), allocatable :: surface_set(:)
type(move_type), allocatable :: move_set(:)
type(revol_type), allocatable :: revol_set(:)

!output variables
integer :: out_elem_region, out_elem_interval, out_elem_str_out_cnr
integer :: mesh_set_num = 0, sep_num(3) = 0, sep_seq(3) = 0
real :: sep_point(2) = 0.0,  cri_tol_length = 0.0
integer, allocatable :: remesh_flag(:)

!remesh variables
real :: accumulate_move_dis
integer :: move_flag, abla_flag
integer :: revol_region_num = 0

contains

!===========================================================================
Subroutine create_mesh_set(unit_num_elem, num_domain, num_materials)

implicit none

integer, intent(in) :: unit_num_elem, num_domain, num_materials

integer :: i, j, k, rn, ng, adjn, sp, nn, count, line_num, first_line_num, num_dim, num_group
integer :: this(3)
integer, allocatable :: adj(:), temp_corn(:,:), temp_cn(:)
real :: ang, dis, tol_dis, min_dis, cri_dis, first_line_dis, line_dis(2)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375
real, allocatable :: node(:,:), temp_edge(:,:)
! for 3D
integer :: bn, bfn, elem_num, face_num, inn, cn, ne, sub_nn, sub_ne, rn_nn, rn_ne
integer, allocatable :: elem(:,:), conn(:), rn_elem(:,:), rn_conn(:)
real, allocatable :: coord(:,:), rn_coord(:,:)
logical :: check

num_dim = region(1)%num_dim
if (num_dim == 2) then
    allocate(mesh_set(num_domain))
    do rn = 1, num_domain
        ng = region(rn)%num_groups
        nn = region(rn)%num_nodes
        
        if (rn == 1) then
            allocate (surface_set(ng))
            surface_set(ng)%part_remesh = .FALSE. ! Don`t use part remeshing routine
        endif
        allocate ( mesh_set(rn)%corn_group(ng) )
        allocate ( mesh_set(rn)%ave_d(ng), mesh_set(rn)%min_d(ng), mesh_set(rn)%flag(ng) )
        allocate ( mesh_set(rn)%nndis(region(rn)%num_skin) )
        allocate ( temp_corn(ng,nn), temp_cn(ng), node(2,nn), temp_edge(2,nn) )
        
        node = region(rn)%node_coord
        temp_corn = 0
        temp_edge = 0.0
        cri_dis = 10e8
        do i = 1, ng
            call set_skin(rn, i, sp, adjn)
            allocate ( adj(adjn) )
            adj = region(rn)%skin_nodes(sp:sp+adjn-1)
            temp_cn(i) = 0
            line_num = 0
            line_dis = 0.0
            !temp_cn(i) = 1
            !temp_corn(i,1) = adj(1)
            if ( i == 1 ) then
                mesh_set(rn)%corn_group(i) = 1
            else
                mesh_set(rn)%corn_group(i) = temp_cn(i-1)+1
            endif
            tol_dis = 0.
            min_dis = 10e8
            do j = 1, adjn
                if (j == 1) then
                    this(:) = (/ adj(adjn), adj(j), adj(j+1) /)
                elseif ( j == adjn ) then
                    this(:) = (/ adj(j-1), adj(j), adj(1) /)
                else
                    this(:) = (/ adj(j-1), adj(j), adj(j+1) /)
                endif
            
                call calc_angle(node(:,this(1)), node(:,this(2)), node(:,this(3)), ang) 
                call calc_len(node(:,this(2)), node(:,this(3)), dis)
                mesh_set(rn)%nndis(sp+j-1) = dis
                tol_dis = tol_dis + dis
                if (min_dis > dis) min_dis = dis

                if (ang < 179.0) dis = dis * sin(ang*pi/180.)
                if (cri_dis > dis) cri_dis = dis
            
                if ( ang < 179.0 .OR. ang > 181.0 ) then
                    temp_cn(i) = temp_cn(i) + 1
                    temp_corn(i,temp_cn(i)) = this(2)
                
                    ! extract edge information
                    if (line_num == 0) then
                        line_num = 1
                        call calc_len(node(:,this(2)), node(:,this(3)), line_dis(1))
                        call calc_len(node(:,this(2)), node(:,this(1)), first_line_dis)
                        first_line_num = 0
                    else
                        if (temp_cn(i) == 1) then
                            call calc_len(node(:,this(2)), node(:,this(1)), first_line_dis)
                            first_line_num = line_num
                        else
                            call calc_len(node(:,this(2)), node(:,this(1)), line_dis(2))
                            count = mesh_set(rn)%corn_group(i)+temp_cn(i)-2 
                            temp_edge(1,count) = line_dis(1)
                            temp_edge(2,count) = line_dis(2)
                        endif
                        call calc_len(node(:,this(2)), node(:,this(3)), line_dis(1))
                        line_num = 1
                    endif
                else
                    line_num = line_num + 1        
                endif
            enddo
            mesh_set(rn)%ave_d(i) = tol_dis / (adjn-1)
            mesh_set(rn)%min_d(i) = min_dis
            
            ! extract last edge information
            count = mesh_set(rn)%corn_group(i)+temp_cn(i)-1
            line_num = line_num + first_line_num
            temp_edge(1,count) = line_dis(1)
            temp_edge(2,count) = first_line_dis
            deallocate ( adj )
        enddo
        mesh_set(rn)%cri_d = cri_dis
        mesh_set(rn)%cn = sum(temp_cn)
        allocate ( mesh_set(rn)%corn(mesh_set(rn)%cn), mesh_set(rn)%edge(2,mesh_set(rn)%cn) )
        allocate ( mesh_set(rn)%corn_coord(2,mesh_set(rn)%cn) )
        mesh_set(rn)%corn = 0
        count = 0
        do i = 1, ng
            do j = 1, temp_cn(i)
                count = count + 1
                mesh_set(rn)%corn(count) = temp_corn(i,j)
            enddo
        enddo
        mesh_set(rn)%edge = temp_edge(:,1:mesh_set(rn)%cn)
        deallocate ( temp_corn, temp_cn, node, temp_edge )
    enddo
elseif (num_dim == 3) then
    ! 1. divide two region(prop / etc)
    allocate(mesh_set(num_domain))
    ng = sub_mesh_region
    
    nn = region(1)%num_nodes
    ne = region(1)%num_elements
    allocate (surface_set(ng))
    allocate(mesh_set(1)%sub_num(ng+1,5))
    mesh_set(1)%sub_num(1,:) = 1
    surface_set(ng)%part_remesh = .FALSE. ! Don`t use part remeshing routine
    
    rewind(unit=unit_num_elem)
    allocate (elem(8,ne), coord(3,nn), conn(nn))
    allocate (rn_elem(8,ne), rn_coord(3,nn), rn_conn(nn))
    do i = 1, divided_mesh_region
        if (i == 1) then
            num_group = ng
        elseif (i == 2) then
            num_group = 1
        endif
        rn_nn = 0
        rn_ne = 0
        do j = 1, num_group
            coord = region(1)%node_coord
            elem = region(1)%element(4:11,:)
            call read_mat_element(unit_num_elem, num_materials, sub_mesh_region, nn, ne, coord, elem, conn, sub_nn, sub_ne, i, j)
            rn_coord(:,rn_nn+1:rn_nn+sub_nn) = coord(:,1:sub_nn)
            rn_elem(:,rn_ne+1:rn_ne+sub_ne) = elem(:,1:sub_ne)+rn_nn
            rn_conn(rn_nn+1:rn_nn+sub_nn) = conn(1:sub_nn)
            rn_nn = rn_nn + sub_nn
            rn_ne = rn_ne + sub_ne
            if (i == 1) then
                mesh_set(i)%sub_num(j+1,1) = rn_nn+1
                mesh_set(i)%sub_num(j+1,2) = rn_ne+1
            endif
        enddo
        mesh_set(i)%nn = rn_nn
        mesh_set(i)%ne = rn_ne
        allocate(mesh_set(i)%node_coord(3,rn_nn), mesh_set(i)%element(8,rn_ne))
        mesh_set(i)%node_coord = rn_coord(:,1:rn_nn)
        mesh_set(i)%element = rn_elem(:,1:rn_ne)
        
        allocate(mesh_set(i)%surf2ori(rn_nn))
        mesh_set(i)%surf2ori = rn_conn(1:rn_nn)
        
        allocate (mesh_set(i)%ave_d(num_group), mesh_set(i)%min_d(num_group), mesh_set(i)%flag(num_group))
        mesh_set(i)%ave_d = 0.0
        mesh_set(i)%min_d = 0.0
        mesh_set(i)%flag = 0
        call get_boundary_info_3D(0, i, num_group)
    enddo
    deallocate (elem, coord, conn)
    deallocate (rn_elem, rn_coord, rn_conn)
endif

allocate(remesh_flag(num_domain))

end subroutine create_mesh_set
!===========================================================================
subroutine set_smoothing_nodes(rn)

implicit none
integer, intent(in) :: rn

integer :: roof, i, j, k, fn, ne, nn, se, ce, cp, line_lv
integer ::temp_nn, temp_ne, num, next_np, temp_np, adjn, count
integer, allocatable :: inte_nodes(:,:), ele(:,:), node_conn(:,:), node_cn(:)
integer, allocatable :: temp_nodes(:), temp_ele(:), temp_point(:), next_point(:)
logical :: check

if (rn == 1) then
    fn = fsi_inte_set%fsi_ni
    allocate(inte_nodes(2,fn))
    inte_nodes = fsi_inte_set%fsi_inte_node
elseif (rn == 2) then
    call set_skin(rn, 1, num, adjn)
    allocate(inte_nodes(2,adjn))
    check = .false.
    fn = 0
    do i = num, num+adjn-1
        if (region(rn)%skin_nodes(i) == inte_set%abla_node(2)) then
            fn = fn + 1
            inte_nodes(:,fn) = (/ rn, region(rn)%skin_nodes(i) /)
            check = .true.
        elseif (region(2)%skin_nodes(i) == inte_set%abla_node(4) .AND. check == .TRUE. ) then
            fn = fn + 1
            inte_nodes(:,fn) = (/ rn, region(rn)%skin_nodes(i) /)
            check = .false.
            exit
        elseif (check) then
            fn = fn + 1
            inte_nodes(:,fn) = (/ rn, region(rn)%skin_nodes(i) /)
        endif
    enddo
    if (check) then
        do i = num, num+adjn-1
            if (region(rn)%skin_nodes(i) == inte_set%abla_node(4)) then
                fn = fn + 1
                inte_nodes(:,fn) = (/ rn, region(rn)%skin_nodes(i) /)
                exit
            else
                fn = fn + 1
                inte_nodes(:,fn) = (/ rn, region(rn)%skin_nodes(i) /)
            endif
        enddo
    endif
endif

ne = region(rn)%num_elements
nn = region(rn)%num_nodes
allocate(ele(4,ne), node_conn(6,nn), node_cn(nn))
allocate(temp_nodes(nn), temp_ele(ne), temp_point(nn), next_point(nn))
ele = region(rn)%element(4:7,:)
temp_ne = 0 ;  temp_ele = 0
temp_nn = 0 ;  temp_nodes = 0
node_cn = 0 ;  node_conn = 0
do i = 1, ne
    do j = 1, 4
        node_cn(ele(j,i)) = node_cn(ele(j,i)) + 1
        node_conn(node_cn(ele(j,i)), ele(j,i)) = i
    enddo
enddo

temp_np = 0
do i = 1, fn
    if (inte_nodes(1,i) == rn) then
        temp_np = temp_np + 1 
        temp_point(temp_np) = inte_nodes(2,i)
    endif
enddo
temp_nn = temp_np
temp_nodes(1:temp_np) = temp_point(1:temp_np)
se = 1

if (rn == 1) then
    if ( ne < 500 ) then
        line_lv = 3
    elseif ( ne >= 500 .AND. ne < 1000 ) then
        line_lv = 4
    else
        line_lv = 5
    endif
elseif (rn == 2) then
    line_lv = 3
endif  
    
do roof = 1, line_lv
    do i = 1, temp_np
        cp = temp_point(i)
        num = node_cn(cp)
        do j = 1, num
            check = .TRUE.
            do k = 1, temp_ne
                if ( node_conn(j,cp) == temp_ele(k) ) then
                    check = .FALSE.
                    exit
                endif
            enddo
            if ( check ) then
                temp_ne = temp_ne + 1
                temp_ele(temp_ne) = node_conn(j,cp)
            endif
        enddo
    enddo
    
    next_np = 0 ;  next_point = 0
    do i = se, temp_ne
        ce = temp_ele(i)
        do j = 1, 4
            check = .TRUE.
            do k = 1, temp_np
                if ( ele(j,ce) == temp_point(k) ) then
                    check = .FALSE.
                    exit
                endif
            enddo
            if ( check ) then
                do k = 1, next_np
                    if ( ele(j,ce) == next_point(k) ) then
                        check = .FALSE.
                        exit
                    endif
                enddo
            endif
            if ( check ) then
                next_np = next_np + 1
                next_point(next_np) = ele(j,ce)
            endif
        enddo
    enddo
    
    if (roof == line_lv) exit
    temp_nodes(temp_nn+1:temp_nn+next_np) = next_point(1:next_np)
    temp_nn = temp_nn + next_np   
    se = temp_ne + 1
    temp_np = next_np
    temp_point = next_point
enddo

if (rn == 1) then
    inte_set%sn = temp_nn
    allocate(inte_set%s_nodes(temp_nn))
    inte_set%s_nodes = temp_nodes(1:temp_nn)
elseif (rn == 2) then
    inte_set%cn = temp_nn
    inte_set%ce = temp_ne
    allocate(inte_set%c_nodes(temp_nn), inte_set%c_elem(temp_ne))
    inte_set%c_nodes = temp_nodes(1:temp_nn)
    inte_set%c_elem = temp_ele(1:temp_ne)
endif

deallocate(inte_nodes, ele, node_conn, node_cn)
deallocate(temp_nodes, temp_ele, temp_point, next_point)

end subroutine set_smoothing_nodes
!===========================================================================
subroutine calc_ori_d(rn, ng, ori_d)

implicit none

integer, intent(in) :: rn, ng
real, intent(inout) :: ori_d

integer :: i, j, adjn, sp, nn, num
integer, allocatable :: adj(:)
real :: tol_dis, temp_dis
real, allocatable :: node(:,:), adj_dis(:)

call set_skin(rn, ng, sp, adjn)
allocate ( adj(adjn), adj_dis(adjn) )
adj = region(rn)%skin_nodes(sp:sp+adjn-1)
nn = region(rn)%num_nodes
allocate (node(2,nn))
node = region(rn)%node_coord
tol_dis = 0.0

do i = 1, adjn
    if (i == 1) then
        call calc_len(node(:,adj(adjn)), node(:,adj(1)), adj_dis(i))
    else
        call calc_len(node(:,adj(i-1)), node(:,adj(i)), adj_dis(i))
    endif
    tol_dis = tol_dis + adj_dis(i)
enddo
do i = 1, adjn-1
    do j = i+1, adjn
        if (adj_dis(i) > adj_dis(j)) then
            temp_dis = adj_dis(i)
            adj_dis(i) = adj_dis(j)
            adj_dis(j) = temp_dis
        endif
    enddo
enddo
num = int(adjn*0.2)+1
ori_d = sum(adj_dis(1:num))/real(num)

deallocate(adj, node)

end subroutine calc_ori_d



subroutine update_mesh_info(DIm, rn, nn, ne, coord, element, flag)

implicit none

integer, intent(in) :: Dim, rn, nn, ne, element((DIm-1)*4,ne), flag
real, intent(in) :: coord(Dim,nn)
integer :: i, j, k, adj(nn+1), adjn(4*nn), nodes(nn)

mesh_set(rn)%nn = nn
mesh_set(rn)%ne = ne

if (Dim == 2) then
    if (allocated(mesh_set(rn)%coord)) deallocate(mesh_set(rn)%coord)
    if (allocated(mesh_set(rn)%ele)) deallocate(mesh_set(rn)%ele)
    allocate(mesh_set(rn)%coord(Dim,nn), mesh_set(rn)%ele((Dim-1)*4,ne))
    mesh_set(rn)%coord = coord(:,1:nn)
    mesh_set(rn)%ele = element(:,1:ne)
    if (flag == 2) then
        if (allocated(mesh_set(rn)%boundary)) deallocate(mesh_set(rn)%boundary)
        allocate (mesh_set(rn)%boundary(nn))
        call METIS_MeshTONodal(ne, nn, element, 4, 1, adj, adjn)

        nodes = 0
        do i = 1, ne
            do j = 1, 4
                nodes(element(j,i)) = nodes(element(j,i)) + 1
            enddo
        enddo
        mesh_set(rn)%boundary = 1
        !write (*,*) 'nn, ne:', nn, ne
        do i = 1, nn  ! over all nodes
          k = adj(i+1) - adj(i)
          !write (*,*) i, k, nodes(i)
          if (k > nodes(i)) then
            mesh_set(rn)%boundary(i) = 0
          endif
        enddo
    endif
elseif (Dim == 3) then
    if (allocated(mesh_set(rn)%node_coord)) deallocate(mesh_set(rn)%node_coord)
    allocate(mesh_set(rn)%node_coord(Dim,nn))
    mesh_set(rn)%node_coord = coord(:,1:nn)
    if (flag == 2) then
        if (allocated(mesh_set(rn)%element)) deallocate(mesh_set(rn)%element)
        allocate(mesh_set(rn)%element((Dim-1)*4,ne))
        mesh_set(rn)%element = element(:,1:ne)
    endif
endif

end subroutine  update_mesh_info



subroutine free_mesh_info(rn)

integer, intent(in) :: rn

deallocate(mesh_set(rn)%coord, mesh_set(rn)%ele)

end subroutine free_mesh_info
!===========================================================================

subroutine change_corner_node(ng, adjn, del_pn, del_p)

implicit none

integer, intent(in) :: ng, del_pn, adjn
integer, intent(inout) :: del_p(adjn)

integer :: i, j, k, num_groups, temp_p, temp_pn
integer :: seq(del_pn), temp_del_p(del_pn)

temp_pn = 0
do i = 1, del_pn
    do j = 1, mesh_set(1)%cn
        if ( del_p(i) == mesh_set(1)%corn(j) ) then
            temp_pn = temp_pn + 1
            seq(temp_pn) = j
            temp_del_p(temp_pn) = del_p(i)
            exit
        endif
    enddo
enddo

do i = 1, temp_pn
    do j = i+1, temp_pn
        if ( temp_del_p(i) < temp_del_p(j) ) then
            temp_p = temp_del_p(i)
            temp_del_p(i) = temp_del_p(j)
            temp_del_p(j) = temp_p
        endif
    enddo
enddo

num_groups = region(1)%num_groups
do i = 1, temp_pn
    do j = 1, mesh_set(1)%cn
        if ( temp_del_p(i) == mesh_set(1)%corn(j) ) then
            do k = j, mesh_set(1)%cn-1
                mesh_set(1)%corn(k) = mesh_set(1)%corn(k+1)
            enddo
            mesh_set(1)%cn = mesh_set(1)%cn - 1
            do k = ng+1, num_groups
                mesh_set(1)%corn_group(k) = mesh_set(1)%corn_group(k) -1
            enddo
            exit
        endif
    enddo
enddo

end subroutine change_corner_node





subroutine adjust_corner(rn, ng, nn, coord, adjn, adj, ngroups, flag)

implicit none

integer, intent(in) :: rn, ng, nn, adjn, adj(adjn), ngroups, flag
real, intent(in) :: coord(2,nn)

integer :: i, j, k, min_j, temp_cn, add_cn, tol_cn, cnum, cn1, cn2, sp, ep, count, num, seq, target_adjn, target_cn
integer :: temp_corn(adjn), tol_corn_group(ngroups), temp_corn_seq(adjn), tol_corn(adjn), this(3), slp(2), cont_node(adjn)
integer :: prop_cont_slp(2)
real :: cp(2), p(2), tol_edge(2,adjn), tol_coord(2,adjn), edge_dis(2)
real :: temp_coord(2,adjn), temp_edge(2,adjn), add_coord(2,adjn), add_edge(2,adjn)
real :: angle, dis, min_dis, r, a0, an, seg_dis, scale_factor, sn, first_dis, large_val
logical :: check, cont_check
integer, allocatable :: ori_corn(:), check_corn(:), target_adj(:), target_cp(:)
real, allocatable :: corn_coord(:,:), ori_edge(:,:)

temp_cn = 0;  temp_corn = 0
temp_edge = 0.0
large_val = 10e8

! find original corner coord
if (flag == 1) then
    cn1 = mesh_set(rn)%corn_group(ng)
    if ( ng == region(rn)%num_groups ) then
        cn2 = mesh_set(rn)%cn
    else
        cn2 = mesh_set(rn)%corn_group(ng+1) - 1
    endif
    temp_cn = cn2 - cn1 + 1
    allocate (ori_corn(temp_cn), corn_coord(2,temp_cn), ori_edge(2,temp_cn))

    count = temp_cn
    ori_corn = mesh_set(rn)%corn(cn1:cn2)
    ori_edge = mesh_set(rn)%edge(:,cn1:cn2)
endif
! end finding original corner coord


! find sequence of slp
if (inte_set%cont_region(1) == 1) then
    prop_cont_slp = inte_set%cont_slp(1:2,ng)
else
    prop_cont_slp = inte_set%cont_slp(3:4,ng)
endif
if (flag == 1) then
    do i = 1, 2
        do j = 1, adjn
            if (adj(j) == prop_cont_slp(i)) then
                slp(i) = j
                exit
            endif
        enddo
    enddo
elseif (flag == 2) then
    ! find the contact points
    call set_skin(2, 1, sp, target_adjn)
    allocate (target_adj(target_adjn), target_cp(target_adjn))
    target_adj = region(2)%skin_nodes(sp:sp+target_adjn-1)
    sp = inte_set%cont_node(2);  ep = inte_set%cont_node(4)
    call find_contact_point(sp, ep, target_adjn, target_adj, flag, target_cn, target_cp)
    
    cont_node = 0
    do i = 1, adjn
        cp = coord(:,adj(i))
        do j = 1, target_cn-1
            call calc_len(region(2)%node_coord(:,target_cp(j)), region(2)%node_coord(:,target_cp(j+1)), dis)
            call near_point(region(2)%node_coord(:,target_cp(j)), region(2)%node_coord(:,target_cp(j+1)), cp, p, dis*0.01, cont_check)
            if (cont_check) then
                cont_node(i) = 1
                exit
            endif
        enddo
    enddo
    deallocate(target_adj, target_cp)
    
    do i = 1, 2
        cp = region(rn)%node_coord(:,prop_cont_slp(i))
        min_dis = 10e8
        do j = 1, adjn
            if (cont_node(j) == 1) then
                call calc_len(coord(:,adj(j)), cp, dis)
                if (dis < min_dis) then
                    min_j = j
                    min_dis = dis
                endif
            endif
        enddo
        slp(i) = min_j
    enddo
endif
write (*,*) 'ball_sp, lp:', slp, adj(slp)

temp_cn = 0
do i = 1, adjn
    if ( i == adjn ) then
        this = (/ adj(i-1), adj(i), adj(1) /)
    elseif ( i == 1 ) then
        this = (/ adj(adjn), adj(i), adj(i+1) /)
    else
        this = (/ adj(i-1), adj(i), adj(i+1) /)
    endif

    call calc_angle(coord(:,this(1)), coord(:,this(2)), coord(:,this(3)), angle)
    
    if ( angle < 179. .OR. angle > 181. ) then
        temp_cn = temp_cn + 1
        temp_corn(temp_cn) = this(2)
        temp_corn_seq(temp_cn) = i
        temp_coord(:,temp_cn) = coord(:,this(2))
        if (flag == 2) then
            call calc_len(coord(:,this(1)), coord(:,this(2)), edge_dis(1))
            call calc_len(coord(:,this(2)), coord(:,this(3)), edge_dis(2))
            if (temp_cn == 1) then
                first_dis = edge_dis(1)
                temp_edge(1,temp_cn) = edge_dis(2)
            else
                temp_edge(2,temp_cn-1) = edge_dis(1)
                temp_edge(1,temp_cn) = edge_dis(2)
            endif
        endif
    endif
enddo
if (flag == 2) temp_edge(2,temp_cn) = first_dis

!write (*,*) 'temp_corn:', temp_corn(1:temp_cn)
!write (*,*) 'temp_corn_seq:', temp_corn_seq(1:temp_cn)
if (flag == 1) then
    allocate (check_corn(temp_cn))
    do i = 1, count
        this(1:2) = (/ i, i+1 /)
        if (i == count) this(2) = 1
        sp = ori_corn(this(1))
        ep = ori_corn(this(2))
        a0 = ori_edge(1,i)
        an = ori_edge(2,i)
        if (sp == ep) then
            write (*,*) 'sp = ep'
            write (*,*) 'sp, ep:', sp, ep
            stop
        endif
        call find_adj_seq(sp, adjn, adj, seq)
        sp = seq
        call find_adj_seq(ep, adjn, adj, seq)
        ep = seq
    
        call calc_len(coord(:,adj(sp)), coord(:,adj(ep)), dis)
        r = ((an-a0)/dis) + 1.0
        cp = coord(:,adj(sp))
    
        num = 0
        if (sp > ep) then
            do j = sp, adjn
                do k = 1, temp_cn
                    if (j == temp_corn_seq(k)) then
                        num = num + 1
                        check_corn(num) = k
                        exit
                    endif
                enddo
            enddo
            do j = 1, ep
                do k = 1, temp_cn
                    if (j == temp_corn_seq(k)) then
                        num = num + 1
                        check_corn(num) = k
                        exit
                    endif
                enddo
            enddo
        else
            do j = sp, ep
                do k = 1, temp_cn
                    if (j == temp_corn_seq(k)) then
                        num = num + 1
                        check_corn(num) = k
                        exit
                    endif
                enddo
            enddo
        endif
        !write (*,*) 'sp, ep:', sp, ep
        !write (*,*) 'check_corn:', check_corn(1:num)
        k = 0
        do j = 1, num
            if (temp_corn_seq(check_corn(j)) == sp) then
                temp_edge(1,check_corn(j)) = a0
            else
                this(3) = temp_corn_seq(check_corn(j))
                call calc_len(cp, coord(:,adj(this(3))), seg_dis)
                sn = 0.0
                find_k: do
                    k = k + 1
                    sn = sn + a0*(r**(k-1))
                    if (sn*1.01 >= seg_dis) then
                        scale_factor = seg_dis/sn
                        exit find_k
                    endif
                enddo find_k
                if (check_corn(j) == 1) then
                    temp_edge(2,temp_cn) = a0*(r**(k-1))*scale_factor
                    temp_edge(1,temp_cn) = temp_edge(1,temp_cn)*scale_factor
                else
                    temp_edge(2,check_corn(j)-1) = a0*(r**(k-1))*scale_factor
                    temp_edge(1,check_corn(j)-1) = temp_edge(1,check_corn(j)-1)*scale_factor
                endif
            
                if (j /= num .or. temp_corn_seq(check_corn(j)) /= ep) then
                    temp_edge(1,check_corn(j)) = a0*(r**k)
                endif
                cp = coord(:,adj(this(3)))
            endif
        enddo
    enddo

    !write (*,*) 'original edge info: ', count
    !do i = 1, count
    !    write (*,*) i, ': ', ori_edge(:,i)
    !enddo
    deallocate (ori_corn, corn_coord, ori_edge)
    deallocate (check_corn)
endif

! Add nodes in the contact region to edge info 
check = .false.
add_coord(:,1:temp_cn) = temp_coord(:,1:temp_cn)
add_edge(:,1:temp_cn) = temp_edge(:,1:temp_cn)
add_cn = temp_cn
temp_cn = 0
add_corner: do i = 1, add_cn
    if (temp_corn_seq(i) == slp(2)) then
        temp_cn = temp_cn + 1
        temp_edge(:,temp_cn) = add_edge(:,i)
        temp_coord(:,temp_cn) = add_coord(:,i)
        check = .false.
        do j = i+1, add_cn
            if (temp_corn_seq(j) == slp(1)) then
                check = .true.
                exit add_corner
            else
                temp_cn = temp_cn + 1
                temp_edge(:,temp_cn) = add_edge(:,j)
                temp_coord(:,temp_cn) = add_coord(:,j)
            endif
        enddo
        do j = 1, i-1
            if (temp_corn_seq(j) == slp(1)) then
                check = .true.
                exit add_corner
            else
                temp_cn = temp_cn + 1
                temp_edge(:,temp_cn) = add_edge(:,j)
                temp_coord(:,temp_cn) = add_coord(:,j)
            endif
        enddo
        if (check == .false.) stop 'Error(adjust_corner): cannot find slp(1)(start point)'
    endif
enddo add_corner
if (slp(1) < slp(2)) then
    do i = slp(1), slp(2)-1
        temp_cn = temp_cn + 1
        temp_edge(:,temp_cn) = large_val
        temp_coord(:,temp_cn) = coord(:,adj(i))
    enddo
else
    do i = slp(1), adjn
        temp_cn = temp_cn + 1
        temp_edge(:,temp_cn) = large_val
        temp_coord(:,temp_cn) = coord(:,adj(i))
    enddo
    do i = 1, slp(2)-1
        temp_cn = temp_cn + 1
        temp_edge(:,temp_cn) = large_val
        temp_coord(:,temp_cn) = coord(:,adj(i))
    enddo
endif

!write (*,*) 'new edge info: ', temp_cn
!do i = 1, temp_cn
!    write (*,*) i, ': ', temp_edge(:,i)
!    write (*,*) i, ': ', temp_coord(:,i)
!enddo

tol_cn = 0;  tol_corn_group = 0
tol_corn = 0
tol_edge = 0.0
tol_coord = 0.0

do i = 1, ngroups
    if ( ng == i ) then
        tol_corn(tol_cn+1:tol_cn+temp_cn) = temp_corn(1:temp_cn)
        tol_edge(:,tol_cn+1:tol_cn+temp_cn) = temp_edge(:,1:temp_cn)
        tol_coord(:,tol_cn+1:tol_cn+temp_cn) = temp_coord(:,1:temp_cn)
        tol_corn_group(i) = tol_cn + 1
        tol_cn = tol_cn + temp_cn
    else
        cn1 = mesh_set(1)%corn_group(i)
        if ( i == ngroups ) then
            cn2 = mesh_set(rn)%cn
        else
            cn2 = mesh_set(rn)%corn_group(i+1) - 1
        endif
        cnum = cn2 - cn1 + 1
        tol_corn(tol_cn+1:tol_cn+cnum) = mesh_set(rn)%corn(cn1:cn2)
        tol_edge(:,tol_cn+1:tol_cn+cnum) = mesh_set(rn)%edge(:,cn1:cn2)
        tol_coord(:,tol_cn+1:tol_cn+cnum) = mesh_set(rn)%corn_coord(:,cn1:cn2)
        tol_corn_group(i) = tol_cn + 1
        tol_cn = tol_cn + cnum
    endif
enddo

deallocate ( mesh_set(rn)%corn, mesh_set(rn)%edge, mesh_set(rn)%corn_coord )
mesh_set(rn)%cn = tol_cn
allocate ( mesh_set(rn)%corn(mesh_set(rn)%cn) )
allocate ( mesh_set(rn)%edge(2,mesh_set(rn)%cn) )
allocate ( mesh_set(rn)%corn_coord(2,mesh_set(rn)%cn) )
mesh_set(rn)%corn = tol_corn(1:tol_cn)
mesh_set(rn)%edge = tol_edge(:,1:tol_cn)
mesh_set(rn)%corn_coord = tol_coord(:,1:tol_cn)
mesh_set(rn)%corn_group = tol_corn_group
write (*,*) ' > corn_group:', mesh_set(rn)%corn_group
write (*,*) ' > corn:', mesh_set(rn)%corn
write (*,*) '>> Complete adjust corner'

end subroutine adjust_corner


!Remesh check routine start ==================================================

subroutine check_boundary_length( rn, ng, nn, coord, adjn, adj, del_pn, del_p)

implicit none

integer, intent(in) :: rn, ng, nn, adjn, adj(adjn)
integer, intent(inout) :: del_pn, del_p(adjn)
real, intent(inout) :: coord(2,nn)

integer :: i, j, k, temp_nn, check_num, num_rn, temp_adjn, ngroups, sp
integer :: this(2), temp(2), inte_node(2), temp_adj(adjn)
real :: tol, dis, tol_dis, angle(2), vec(2,2), temp_coord(2,4), constr_node(2), cp(2), np(2)
logical :: check, del_point_check, constr_check(2), near_check

integer, allocatable :: corner_adj(:)
real, allocatable :: temp_node(:,:)
!--------------------------------------------------------------------------------
num_rn = num_of_regions
if (num_rn == 1) then
    !tol = 0.1*mesh_set(rn)%ave_d(ng)
    tol = 0.30
else
    !tol = 0.5*mesh_set(rn)%ave_d(ng)
    tol = 0.30
endif
sp = region(rn)%skin_group_pointer(ng)

temp_nn = nn
allocate ( temp_node(2,temp_nn) )
temp_node(:,1:temp_nn) = coord(:,1:temp_nn)
check_num = 0

do j = 1, adjn
    if ( check_num /= adj(j) ) then
        tol_dis = mesh_set(rn)%nndis(sp-1+j)
        if ( j == adjn ) then
            this = (/ j, 1 /)
            !tol_dis = (mesh_set(rn)%nndis(sp-2+j)+mesh_set(rn)%nndis(sp-1+j)+mesh_set(rn)%nndis(sp))/3.0
        elseif ( j == 1 ) then
            this = (/ j, j+1 /)
            !tol_dis = (mesh_set(rn)%nndis(sp-1+adjn)+mesh_set(rn)%nndis(sp)+mesh_set(rn)%nndis(sp+1))/3.0
        else 
            this = (/ j, j+1 /)
            !tol_dis = (mesh_set(rn)%nndis(sp-2+j)+mesh_set(rn)%nndis(sp-1+j)+mesh_set(rn)%nndis(sp+j))/3.0
        endif
        
        del_point_check = .TRUE.
        do i = 1, 2
            do k = 1, del_pn
                if (adj(this(i)) == del_p(k)) then
                    del_point_check = .FALSE.
                    exit                    
                endif
            enddo
            if (del_point_check == .FALSE.) exit
        enddo                    
                
        call calc_len(temp_node(:,adj(this(1))), temp_node(:,adj(this(2))), dis)
        if ( tol_dis*tol > dis .AND. del_point_check) then
            !write (*,*) "Short length!!!"
            !write (*,*) temp_node(:,adj(this(1)))
            !write (*,*) temp_node(:,adj(this(2)))
            !write (*,*) adj(this(1:2)), dis
    
            inte_node = 0
            do i = 1, fsi_inte_set%fsi_ni
                if ( rn == fsi_inte_set%fsi_inte_node(1,i) .AND. adj(this(1)) == fsi_inte_set%fsi_inte_node(2,i) ) inte_node(1) = 1
                if ( rn == fsi_inte_set%fsi_inte_node(1,i) .AND. adj(this(2)) == fsi_inte_set%fsi_inte_node(2,i) ) inte_node(2) = 1
                if ( sum(inte_node) == 2 ) exit
            enddo
                
            constr_node = 0
            if ( sum(inte_node) /= 1 ) then
                constr_check = region(rn)%dof_activity_flag(region(rn)%node(4:5,adj(this(1))))
                if (constr_check(1) == .FALSE. .OR. constr_check(2) == .FALSE.) constr_node(1) = 1
                constr_check = region(rn)%dof_activity_flag(region(rn)%node(4:5,adj(this(2))))
                if (constr_check(1) == .FALSE. .OR. constr_check(2) == .FALSE.) constr_node(2) = 1
            endif
                
            if ( sum(inte_node) == 1 ) then
                del_pn = del_pn + 1
                if ( inte_node(1) == 1 ) then
                    del_p(del_pn) = adj(this(2))
                else
                    del_p(del_pn) = adj(this(1))
                endif
            elseif ( sum(constr_node) == 1) then
                del_pn = del_pn + 1
                if ( constr_node(1) == 1 ) then
                    del_p(del_pn) = adj(this(2))
                else
                    del_p(del_pn) = adj(this(1))
                endif
            else
                if ( this(1) == 1 ) then
                    temp(1:2) = (/ adj(adjn), adj(this(2)+1) /)
                elseif ( this(2) == adjn ) then
                    temp(1:2) = (/ adj(this(1)-1), adj(1) /)
                elseif ( this(1) == adjn ) then
                    temp(1:2) = (/ adj(this(1)-1), adj(2) /)
                else
                    temp(1:2) = (/ adj(this(1)-1), adj(this(2)+1) /)
                endif
                    
                call calc_angle(temp_node(:,temp(1)), temp_node(:,adj(this(1))), temp_node(:,adj(this(2))), angle(1))
                call calc_angle(temp_node(:,adj(this(1))), temp_node(:,adj(this(2))), temp_node(:,temp(2)), angle(2))
                    
                if ((angle(1) > 170. .AND. angle(1) < 190.) .AND. (angle(2) > 170. .AND. angle(2) < 190.)) then !.OR. (angle(1)+angle(2) < 270)) then
                    if ( j /= adjn ) then
                        coord(:,adj(this(1))) = 0.5*(temp_node(:,adj(this(1)))+temp_node(:,adj(this(2))))
                    endif
                    write (*,*) adj(this(1)), " and", adj(this(2)), "Both Angles are 170~190!!"
                elseif ( angle(1) < 120. .AND. angle(2) < 120. ) then
                    write (*,*) region(rn)%dof_activity_flag(region(rn)%node(4:5,adj(this(1))))
                    write (*,*) region(rn)%dof_activity_flag(region(rn)%node(4:5,adj(this(2))))
                    if ( j /= adjn ) then
                        coord(:,adj(this(1))) = 0.5*(temp_node(:,adj(this(1)))+temp_node(:,adj(this(2))))
                    endif
                    write (*,*) adj(this(1)), " and", adj(this(2)), "Both Angles are small(<110)!!"
                else
                    vec(:,1) = temp_node(:,temp(1)) - temp_node(:,adj(this(1)))
                    vec(:,2) = temp_node(:,temp(2)) - temp_node(:,adj(this(2)))
                    temp_coord(:,1) = temp_node(:,adj(this(1))) - (100*vec(:,1))
                    temp_coord(:,2) = temp_node(:,adj(this(1))) + (100*vec(:,1))
                    temp_coord(:,3) = temp_node(:,adj(this(2))) - (100*vec(:,2))
                    temp_coord(:,4) = temp_node(:,adj(this(2))) + (100*vec(:,2))
                    call calc_cross(temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), temp_coord(:,4), check)
                    if (check == .false.) then
                        call calc_cross_p(temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), temp_coord(:,4), cp) !coord(:,adj(this(2))))
                        call near_point(temp_node(:,adj(this(1))), temp_node(:,adj(this(2))), cp, np, dis*0.01, near_check)
                        if (near_check) then
                            coord(:,adj(this(2))) = cp
                        else
                            coord(:,adj(this(2))) = (temp_node(:,adj(this(1)))+temp_node(:,adj(this(2))))*0.5
                        endif
                        if (j /= adjn) then
                            coord(:,adj(this(1))) = coord(:,adj(this(2)))
                        endif
                    !if ( check == .FALSE. .AND. j == adjn ) then
                    !    call calc_cross_p(temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), temp_coord(:,4), cp) !coord(:,adj(this(2))))
                    !    
                    !elseif ( check == .FALSE. .AND. j /= adjn ) then
                    !    call calc_cross_p(temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), temp_coord(:,4), cp) !coord(:,adj(this(1))))
                    !    coord(:,adj(this(2))) = coord(:,adj(this(1)))
                    !    !write (*,*) temp(1), adj(this(1:2)), temp(2)
                    !    !do i = 1, mesh_set(rn)%cn
                    !    !    if (mesh_set(1)%corn(i) == adj(this(2))) then
                    !    !        write (*,*) "convex node : node number change!!"
                    !    !        write (*,*) adj(this(1)), ' -> ', adj(this(2))
                    !    !        this(2) = this(1)
                    !    !        exit
                    !    !    endif
                    !    !enddo
                        if ( (angle(2) < 179.) .AND. angle(1) > 179.) then
                            write (*,*) "convex node : node number change!!"
                            write (*,*) adj(this(1)), ' -> ', adj(this(2))
                            this(2) = this(1)
                        endif
                    else
                        write (*,*) "Can't find cross point!!(check_boundary_length)"
                        stop
                    endif
                    
                    write (*,*) adj(this(1)), " and", adj(this(2)), "Both Angles are not 180!!"
                endif
                del_pn = del_pn + 1
                if ( j == adjn ) then
                    del_p(del_pn) = adj(this(1))
                else
                    del_p(del_pn) = adj(this(2))
                endif
                check_num = del_p(del_pn)
            endif                
        endif
    endif
enddo
deallocate (temp_node)

if ( del_pn /= 0 ) then
    remesh_flag(rn) = 2
    !mesh_set(rn)%flag(ng) = 3
    mesh_set(rn)%flag(ng) = 1
    surface_set(ng)%coord = coord
    ngroups = region(rn)%num_groups
    temp_adjn = adjn
    temp_adj = adj
    
    do i = del_pn, 1, -1
        do j = 1, temp_adjn 
            if (del_p(i) == temp_adj(j)) then
                temp_adjn = temp_adjn - 1
                do k = j, temp_adjn
                    temp_adj(k) = temp_adj(k+1)
                enddo
                exit
            endif
        enddo
    enddo
    allocate (corner_adj(temp_adjn))
    corner_adj = temp_adj(1:temp_adjn)
    call adjust_corner(rn, ng, nn, coord, temp_adjn, corner_adj, ngroups, 1)  
    deallocate (corner_adj)
    !sp = del_p(1)
    !call check_part_element(rn, ng, temp_adjn, corner_adj, 12, sp)
endif

end subroutine check_boundary_length




subroutine check_boundary_angle(rn, nn, coord, adjn, adj, del_pn, del_p)


implicit none

integer, intent(in) :: rn, nn, adjn, adj(adjn)
integer, intent(inout) :: del_pn, del_p(adjn)
real, intent(inout) :: coord(2,nn)

integer :: i, j, check_num, num_rn, temp_nn
integer :: this(3), inte_node(2)
real :: tol, ang
logical :: last_node

real, allocatable :: temp_node(:,:)

num_rn = num_of_regions
if (num_rn == 1) then
    tol = 25.
else
    tol = 22.
endif

temp_nn = nn
allocate ( temp_node(2,temp_nn) )
temp_node(:,1:temp_nn) = coord(:,1:temp_nn)

check_num = 0
last_node = .FALSE.
do i = 1, adjn
    if ( check_num /= adj(i) ) then
        if (i == 1) then
            this = (/ adjn, 1,  2/)
        elseif (i == adjn) then
            this = (/ i-1, i, 1 /)
        else
            this = (/ i-1, i, i+1 /)
        endif
        
        call calc_angle(temp_node(:,adj(this(1))), temp_node(:,adj(this(2))), temp_node(:,adj(this(3))), ang)
        if (tol > ang .OR. tol > abs(360-ang)) then
            write (*,*) "small angle!!!"
            write (*,*) temp_node(:,adj(this(1)))
            write (*,*) temp_node(:,adj(this(2)))
            write (*,*) temp_node(:,adj(this(3)))
            write (*,*) adj(this(1:3)), ang
            
            inte_node = 0
            do j = 1, fsi_inte_set%fsi_ni
                if ( rn == fsi_inte_set%fsi_inte_node(1,j) .AND. adj(this(1)) == fsi_inte_set%fsi_inte_node(2,j) ) inte_node(1) = 1
                if ( rn == fsi_inte_set%fsi_inte_node(1,j) .AND. adj(this(3)) == fsi_inte_set%fsi_inte_node(2,j) ) inte_node(2) = 1
                if ( sum(inte_node) == 2 ) exit
            enddo
            
            del_p(del_pn+1) = adj(this(1))
            del_p(del_pn+2) = adj(this(3))
            del_pn = del_pn + 2
            check_num = adj(this(3))
            if (j == 1) last_node = .TRUE.
            if (j == adjn-1 .AND. last_node) check_num = adj(adjn)

            if (sum(inte_node) == 1) then
                if (inte_node(1) == 1) then
                    coord(:,adj(this(2))) = (temp_node(:,adj(this(2)))+temp_node(:,adj(this(3))))*0.5
                else
                    coord(:,adj(this(2))) = (temp_node(:,adj(this(2)))+temp_node(:,adj(this(1))))*0.5
                endif              
            elseif (sum(inte_node) == 2) then
                coord(:,adj(this(2))) = (temp_node(:,adj(this(1)))+temp_node(:,adj(this(2)))+temp_node(:,adj(this(2))))/3.0
            endif
        endif
    endif
enddo

end subroutine check_boundary_angle


subroutine sharp_check(rn, ng, nn, coord, adjn, adj)

implicit none

integer, intent(in) :: rn, ng, nn, adjn, adj(adjn)
real, intent(inout) :: coord(2,nn)

integer :: j, adj_sp
integer :: temp_adjn, check_pn, ngroups, num, flag, count
integer :: this(3), temp(3), temp_adj(adjn)
integer, allocatable :: corner_adj(:)
real :: tol, angle, cri_angle(2)
real :: temp_coord(2,3), prop_coord(2,nn), dis(2), check_dis(3)
logical :: check
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

! save information of case
ngroups = region(rn)%num_groups
tol = 0.3
cri_angle = (/ 22., 22. /)
adj_sp = region(rn)%skin_group_pointer(ng)
prop_coord = coord
temp_adjn = adjn
temp_adj = adj
check_pn = inte_set%cont_slp(2,ng)
if (ng == ngroups .AND. fsi_inte_set%fsi_inte_point(3) == 1) check_pn = fsi_inte_set%fsi_inte_point(4)

! last node
num = 0
do 
    num = num + 1
    check = .TRUE.
    call find_adj_seq(check_pn, temp_adjn, temp_adj, temp(2))
    temp(1) = temp(2) - 1
    temp(3) = temp(2) + 1
    if (temp(2) == temp_adjn) then
        temp(3) = 1
    elseif ( temp(2) == 1 ) then
        temp(1) = temp_adjn
    endif
    
    if (num == 1) then
        check_dis(2) = mesh_set(rn)%nndis(adj_sp-1+temp(2))
        check_dis(1) = mesh_set(rn)%nndis(adj_sp-1+temp(1))
        check_dis(3) = (check_dis(1)+check_dis(2))*0.5
    endif
    
    this = (/ temp_adj(temp(1)), temp_adj(temp(2)), temp_adj(temp(3)) /)
    temp_coord(:,1) = coord(:,this(1))
    temp_coord(:,2) = coord(:,this(2))
    temp_coord(:,3) = coord(:,this(3))

    call calc_angle( temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), angle)
    call calc_len(temp_coord(:,1), temp_coord(:,2), dis(1))
    call calc_len(temp_coord(:,2), temp_coord(:,3), dis(2))
    flag = 0
    if ( check_dis(3)*tol > dis(2)) then
        flag = 1
    endif
    if ( check_dis(3)*tol > dis(1)) then
        if (flag == 1) then
            flag = 3                          
        else
            flag = 2                                              
        endif
    elseif ( check_dis(3)*tol > sin(angle*pi/180.)*dis(2)) then  
        if (flag == 1) then
            flag = 3
        else
            flag = 1
        endif
    !elseif (num > 1 .AND. angle < cri_angle(2)) then
    !    flag = 1
    endif
    
    ! flag: 1 = (1) node delete
    !       2 = (3) node delete
    !       3 = (1), (3) node delete
    if (flag /= 0) then
        write (*,*) 'Shart_check(Case-1:Last point), flag:', flag
        check = .false.
        check_pn = this(2)
        mesh_set(rn)%flag(ng) = 1
        
        if (flag == 1) then
            do j = temp(3), temp_adjn-1
                temp_adj(j) = temp_adj(j+1)
            enddo
            temp_adjn = temp_adjn - 1
            check_pn = this(2)
        elseif (flag == 2) then
            do j = temp(1), temp_adjn-1
                temp_adj(j) = temp_adj(j+1)
            enddo
            temp_adjn = temp_adjn - 1
        elseif (flag == 3) then
            if (temp(2) == temp_adjn) then
                count = 1
            else
                count = 2
            endif
            do j = temp(3), temp_adjn-count
                temp_adj(j) = temp_adj(j+count)
            enddo
            temp_adjn = temp_adjn - 2
            coord(:,this(3)) = temp_coord(:,2)
        endif
    endif
    
    if ( check .OR. temp_adjn <= 10 ) exit
enddo

!fsi_inte_set%fsi_inte_coord(:,2) = coord(:,check_pn)
if (ng == ngroups .AND. fsi_inte_set%fsi_inte_point(3) == 1) fsi_inte_set%fsi_inte_point(4) = check_pn
check_pn = inte_set%cont_slp(1,ng)
if (ng == 1.AND. fsi_inte_set%fsi_inte_point(1) == 1) check_pn = fsi_inte_set%fsi_inte_point(2)
! start node
num = 0
do 
    num = num + 1
    check = .TRUE.
    call find_adj_seq(check_pn, temp_adjn, temp_adj, temp(2))
    temp(1) = temp(2) - 1
    temp(3) = temp(2) + 1
    if (temp(2) == temp_adjn) then
        temp(3) = 1
    elseif ( temp(2) == 1 ) then
        temp(1) = temp_adjn
    endif
    if (num == 1) then
        check_dis(2) = mesh_set(rn)%nndis(adj_sp-1+temp(2))
        check_dis(1) = mesh_set(rn)%nndis(adj_sp-1+temp(1))
        check_dis(3) = (check_dis(1)+check_dis(2))*0.5
    endif

    this = (/ temp_adj(temp(1)), temp_adj(temp(2)), temp_adj(temp(3)) /)
    !temp_coord(:,1) = prop_coord(:,this(1))
    !temp_coord(:,2) = prop_coord(:,this(2))
    !temp_coord(:,3) = prop_coord(:,this(3))
    temp_coord(:,1) = coord(:,this(1))
    temp_coord(:,2) = coord(:,this(2))
    temp_coord(:,3) = coord(:,this(3))
    
    call calc_angle( temp_coord(:,1), temp_coord(:,2), temp_coord(:,3), angle)
    call calc_len(temp_coord(:,1), temp_coord(:,2), dis(1))
    call calc_len(temp_coord(:,2), temp_coord(:,3), dis(2))
    flag = 0
    if ( check_dis(3)*tol > dis(1)) then
        write (*,*) 'sharp_check: case1'
        flag = 1
    endif
    if ( check_dis(3)*tol > dis(2)) then
        if (flag == 1) then
            write (*,*) 'sharp_check: case4'
            flag = 3                          
        else
            write (*,*) 'sharp_check: case2'
            flag = 2                                              
        endif
    elseif ( check_dis(3)*tol > sin(angle*pi/180.)*dis(1)) then  
        if (flag == 1) then
            write (*,*) 'sharp_check: case5'
            flag = 3
        else
            write (*,*) 'sharp_check: case3'
            flag = 1
        endif
    !elseif (num > 1 .AND. angle < cri_angle(2)) then
    !    write (*,*) 'case6'
    !    flag = 1
    endif
    
    ! flag: 1 = (3) node delete
    !       2 = (1) node delete
    !       3 = (1), (3) node delete
    if (flag /= 0) then
        write (*,*) 'Shart_check(Case-2:Start point), flag:', flag
        write (*,*) 'temp:', temp(1:3)
        write (*,*) 'this:', this(1:3)
        check = .false.
        check_pn = this(2)
        mesh_set(rn)%flag(ng) = 1
        
        if (flag == 1) then
            do j = temp(1), temp_adjn-1
                temp_adj(j) = temp_adj(j+1)
            enddo
            temp_adjn = temp_adjn - 1
            check_pn = this(2)
        elseif (flag == 2) then
            do j = temp(3), temp_adjn-1
                temp_adj(j) = temp_adj(j+1)
            enddo
            temp_adjn = temp_adjn - 1
        elseif (flag == 3) then
            if (temp(1) == temp_adjn) then
                count = 1
            else
                count = 2
            endif
            do j = temp(3), temp_adjn-count
                temp_adj(j) = temp_adj(j+count)
            enddo
            temp_adjn = temp_adjn - 2
            coord(:,this(3)) = temp_coord(:,2)
        endif
    endif
    
    if ( check .OR. temp_adjn <= 10 ) exit
enddo
!fsi_inte_set%fsi_inte_coord(:,1) = coord(:,check_pn)
if (ng == 1.AND. fsi_inte_set%fsi_inte_point(1) == 1) fsi_inte_set%fsi_inte_point(2) = check_pn

if ( mesh_set(rn)%flag(ng) == 1) then ! .AND. surface_set(ng)%part_remesh == .FALSE.) then
    remesh_flag(rn) = 2
    ngroups = region(rn)%num_groups
    surface_set(ng)%coord = coord
    
    allocate (corner_adj(temp_adjn))
    corner_adj = temp_adj(1:temp_adjn)
    call adjust_corner(rn, ng, nn, coord, temp_adjn, corner_adj, ngroups, 1)
    deallocate (corner_adj)
!elseif ( surface_set(ng)%part_remesh ) then
!    allocate (corner_adj(temp_adjn))
!    corner_adj = temp_adj(1:temp_adjn)
!    call check_part_element(rn, ng, temp_adjn, corner_adj, 20, sp)
!    remesh_flag(rn) = 2
!    mesh_set(rn)%flag(ng) = 3
!    surface_set(ng)%coord = coord
endif

end subroutine sharp_check

subroutine check_non_matching_disp(ngroup, check, current_step)

implicit none

integer, intent(in) :: ngroup, current_step
logical, intent(inout) :: check(ngroup)

integer :: i, j, k, case_nn, prop_nn, ng, output_flag
integer :: case_adjn, prop_adjn, sp, lp, cn, count, seq, pos(2)
real :: mp(2), np(2), target_coord(2,2), cont_off, dist, dis(3), vec(2), tol
logical :: cont_check(2), odd_check
integer, allocatable :: case_adj(:), prop_adj(:), cp(:), this(:), cont_node(:)
real, allocatable :: case_coord(:,:), prop_coord(:,:), case_disp(:,:), prop_disp(:,:), temp_coord(:,:)
character(len=60) :: fn

output_flag = 1
prop_nn = region(1)%num_nodes
case_nn = region(2)%num_nodes
allocate (prop_coord(2,prop_nn), prop_disp(2,prop_nn))
allocate (case_coord(2,case_nn), case_disp(2,case_nn))
do i=1,prop_nn
    prop_coord(:,i) = region(1)%node_coord(:,i)
    prop_disp(:,i) = history(1)%data(region(1)%node(4:5,i), 4)
enddo
do i=1,case_nn
    case_coord(:,i) = region(2)%node_coord(:,i)
    case_disp(:,i) = history(2)%data(region(2)%node(4:5,i), 4)
enddo

call set_skin(2, 1, sp, case_adjn)
allocate (case_adj(case_adjn), cp(case_adjn))
case_adj = region(2)%skin_nodes(sp:sp+case_adjn-1) 
sp = inte_set%cont_node(2);  lp = inte_set%cont_node(4)
call find_contact_point(sp, lp, case_adjn, case_adj, 2, cn, cp)

check = .false.
do ng = 1, ngroup
    call set_skin(1, ng, sp, prop_adjn)
    allocate (prop_adj(prop_adjn), this(prop_adjn))
    prop_adj = region(1)%skin_nodes(sp:sp+prop_adjn-1)
    pos = inte_set%cont_slp(:,ng)
    do i = 1, 2
        call find_adj_seq(pos(i), prop_adjn, prop_adj, seq)
        pos(i) = seq
    enddo
    count = 0
    if ( pos(1) > pos(2) ) then
        do i = pos(1), prop_adjn
            count = count + 1
            this(count) = i
        enddo
        do i = 1, pos(2)
            count = count + 1
            this(count) = i
        enddo
    else
        do i = pos(1), pos(2)
            count = count + 1
            this(count) = i
        enddo
    endif

    do j = 2, count-1
        mp(:) = prop_coord(:,prop_adj(this(j)))
        np = 0.0
        do k = 1, cn-1
            target_coord(:,1) = case_coord(:,cp(k))
            target_coord(:,2) = case_coord(:,cp(k+1))
            call calc_len(target_coord(:,1), target_coord(:,2), cont_off)
            call near_point(target_coord(:,1), target_coord(:,2), mp, np, cont_off*0.01, cont_check(1))
            if (cont_check(1)) then
                target_coord(:,1) = case_coord(:,cp(k)) + case_disp(:,cp(k))
                target_coord(:,2) = case_coord(:,cp(k+1)) + case_disp(:,cp(k+1))
                mp = prop_coord(:,prop_adj(this(j))) + prop_disp(:,prop_adj(this(j)))
                call dist_between_point_line(target_coord(:,1), target_coord(:,2), mp, cont_off*0.01, dist)
                call calc_len(mp(:), prop_coord(:,prop_adj(this(j-1))), dis(1))
                call calc_len(mp(:), prop_coord(:,prop_adj(this(j+1))), dis(2))
                tol = min(dis(1), dis(2))
                !write (*,*) prop_adj(this(j)), 'dist:', dist, cont_off*0.25
                if (dist > tol*0.2) then
                    write (*,*) prop_adj(this(j)), 'dist:', dist, tol*0.2
                    check(ng) = .true.
                endif
                exit
            endif
        enddo  !(k)
        if (check(ng)) exit
    enddo
    
    if (check(ng)) then
        allocate(temp_coord(2,prop_adjn*2))
        temp_coord = 0.0
        count = 0
        if ( pos(2) > pos(1) ) then
            do i = pos(2), prop_adjn
                count = count + 1
                temp_coord(:,count) = prop_coord(:,prop_adj(i))
            enddo
            do i = 1, pos(1)
                count = count + 1
                temp_coord(:,count) = prop_coord(:,prop_adj(i))
            enddo
        else
            do i = pos(2), pos(1)
                count = count + 1
                temp_coord(:,count) = prop_coord(:,prop_adj(i))
            enddo
        endif
        lp = count
        mp(:) = temp_coord(:,count)
        do i = cn, 2, -1
            target_coord(:,1) = case_coord(:,cp(i))
            target_coord(:,2) = case_coord(:,cp(i-1))
            call calc_len(target_coord(:,1), target_coord(:,2), cont_off)
            call near_point(target_coord(:,1), target_coord(:,2), mp, np, cont_off*0.01, cont_check(1))
            if (cont_check(1)) then
                count = count + 1
                temp_coord(:,count) = target_coord(:,2)
                mp = temp_coord(:,1)
                do j = i-1, 2, -1
                    target_coord(:,1) = case_coord(:,cp(j))
                    target_coord(:,2) = case_coord(:,cp(j-1))
                    call calc_len(target_coord(:,1), target_coord(:,2), cont_off)
                    call near_point(target_coord(:,1), target_coord(:,2), mp, np, cont_off*0.01, cont_check(2))
                    if (cont_check(2)) then
                        exit !(j)
                    else
                        count = count + 1
                        temp_coord(:,count) = target_coord(:,2)
                    endif
                enddo
                exit
            endif
        enddo  
        call calc_len(temp_coord(:,1), temp_coord(:,2), dis(1))
        call calc_len(temp_coord(:,1), temp_coord(:,count), dis(2))
        call calc_len(temp_coord(:,1), temp_coord(:,count-1), dis(3))
        if (dis(1) > dis(2)+dis(3)*0.25) then
            count = count - 1
        else
            temp_coord(:,count) = (temp_coord(:,count-1)+temp_coord(:,1))*0.5
        endif
        call calc_len(temp_coord(:,lp), temp_coord(:,lp-1), dis(1))
        call calc_len(temp_coord(:,lp), temp_coord(:,lp+1), dis(2))
        if (lp+1 /= count) then
            call calc_len(temp_coord(:,lp), temp_coord(:,lp+2), dis(3))
            if (dis(1) > dis(2)+dis(3)*0.25) then
                do i = lp+1, count-1
                    temp_coord(:,i) = temp_coord(:,i+1)
                enddo
                count = count - 1
            else
                temp_coord(:,lp+1) = (temp_coord(:,lp+2)+temp_coord(:,lp))*0.5
            endif
        endif
        deallocate(prop_adj)
        allocate(prop_adj(count), cont_node(count))
        prop_adjn = count
        do i = 1, count
            prop_adj(i) = i
        enddo
        cont_node = 0
        cont_node(1) = 1
        do i = lp, count
            cont_node(i) = 1
        enddo
        call odd_to_even_node(prop_adjn, temp_coord(:,1:prop_adjn+1), cont_node, odd_check)
        deallocate (cont_node)
        
        if (output_flag == 1) then
            fn = './output/solid/remesh/check_rearrange0000000_2.plt'
            write (fn(38:44), '(I7.7)' ) current_step
            open(unit=30, file=fn, status='REPLACE', action='WRITE')
            Write (30,'(A)') 'TITLE="Check random edge"'
            Write (30,*) 'VARIABLES="x", "y"'
            Write (30,'(A,I5,A,I5, A)') 'ZONE N =', prop_adjn , ', E =', prop_adjn, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
            do i = 1, prop_adjn
                write (30,*) temp_coord(:,i)
            enddo
            do i = 1, prop_adjn
                if (i == prop_adjn) then
                    write (30,*) i, 1
                else
                    write (30,*) i, i+1
                endif
            enddo
            close(30)
        endif
        call adjust_corner(1, ng, count, temp_coord(:,1:count), prop_adjn, prop_adj, ngroup, 2)
        if (allocated(surface_set(ng)%coord)) deallocate(surface_set(ng)%coord)
        allocate(surface_set(ng)%coord(2,prop_adjn))
        surface_set(ng)%nn = prop_adjn
        surface_set(ng)%coord(:,1:prop_adjn) = temp_coord(:,1:prop_adjn)
        
        mesh_set(1)%flag(ng) = 3
        deallocate(temp_coord)
    endif
    deallocate(prop_adj, this)
enddo

deallocate (prop_coord, prop_disp, case_coord, case_disp)
deallocate (case_adj, cp)

end subroutine check_non_matching_disp

!Remesh check routine end ====================================================

!Contact routine start =========================================================

subroutine adjust_contact(prop_nn, prop_node)

implicit none

integer, intent(in) :: prop_nn
real, intent(inout) :: prop_node(2,prop_nn)

integer :: i, target_nn, ngroup
integer :: case_adjn, sp, prop_adjn, ball_nn
real :: cri_ang(2)
integer, allocatable :: case_adj(:)
real, allocatable :: target_node(:,:), ball_node(:,:)

ball_nn = region(1)%num_nodes
target_nn = region(2)%num_nodes
allocate ( target_node(2,target_nn), ball_node(2,ball_nn) )
do i=1,target_nn
    target_node(1,i) =  region(2)%node_coord(1,i) ! + history(2)%data(  region(2)%node(4,i) , 4 )
    target_node(2,i) =  region(2)%node_coord(2,i) ! + history(2)%data(  region(2)%node(5,i) , 4 )
enddo
ngroup = region(1)%num_groups

cri_ang(:) = (/ 160., 210. /)

call set_skin(2, 1, sp, case_adjn)
allocate (case_adj(case_adjn))
case_adj = region(2)%skin_nodes(sp:sp+case_adjn-1) 

! prop_node -> ball_node
!do i = 1, region(1)%num_skin
!    ball_node(:,region(1)%skin_nodes(i)) = prop_node(:,i)
!enddo
ball_node = prop_node

call move_contact_point(ngroup, case_adjn, case_adj, target_nn, target_node, ball_nn, ball_node, cri_ang, 1)

! ball_node -> prop_node
!do i = 1, region(1)%num_skin
!    prop_node(:,i) = ball_node(:,region(1)%skin_nodes(i))
!enddo
prop_node = ball_node

deallocate(target_node, ball_node, case_adj)

end subroutine adjust_contact


subroutine remesh_contact(ng, prop_nn, prop_coord)

implicit none

integer, intent(inout) :: ng, prop_nn
real, intent(inout) :: prop_coord(2,prop_nn)

integer, allocatable :: case_adj(:)
real, allocatable :: case_coord(:,:)

integer :: i, ngroup, case_nn, case_adjn, num
real :: cri_ang(2)

cri_ang(:) = (/ 160., 210. /)
call set_skin(2, 1, num, case_adjn)
allocate (case_adj(case_adjn))
case_adj = region(2)%skin_nodes(num:num+case_adjn-1) 

case_nn = region(2)%num_nodes
allocate ( case_coord(2,case_nn) )
do i=1,case_nn
    case_coord(1,i) =  region(2)%node_coord(1,i) ! + history(2)%data(  region(2)%node(4,i) , 4 )
    case_coord(2,i) =  region(2)%node_coord(2,i) ! + history(2)%data(  region(2)%node(5,i) , 4 )
enddo
!ngroup = region(1)%num_groups
ngroup = ng

call move_contact_point(ngroup, case_adjn, case_adj, case_nn, case_coord, prop_nn, prop_coord, cri_ang, 2)

deallocate(case_adj, case_coord)

end subroutine remesh_contact


subroutine move_contact_point(ngroup, target_adjn, target_adj, target_nn, target_node, ball_nn, ball_node, cri_ang, flag)

!flag: 1 - contact after move
!      2 - contact after remesh
implicit none

integer, intent(in) :: flag, ngroup, target_adjn, target_nn, ball_nn, target_adj(target_adjn)
real, intent(in) :: target_node(2, target_nn), cri_ang(2)
real, intent(inout) :: ball_node(2, ball_nn)

integer :: i, j, k, temp_j, sp, lp, temp_k, count, ng, cn, adjn, seq
integer :: pos(2), prop_cont_slp(2)
integer, allocatable :: cp(:), adj(:), this(:)
real :: ang, min_dis, max_dis, cont_off, np_dis, temp_cont_off
real :: np(2), mp(2), dis(2), temp_np(2), temp_p(2)
logical :: check, cont_check

! find the contact points
sp = inte_set%cont_node(2);  lp = inte_set%cont_node(4)
allocate(cp(target_adjn))
!call find_contact_point(sp, lp, target_adjn, target_adj, flag, cn, cp)
call find_contact_point(sp, lp, target_adjn, target_adj, 2, cn, cp)   ! temporory
max_dis = 0.0
do i = 1, cn-1
    call calc_len(target_node(:,cp(i)), target_node(:,cp(i+1)), dis(1))
    if (max_dis < dis(1)) then
        max_dis = dis(1)
    endif
enddo

if (region(1)%num_groups /= ngroup) then
    call delete_cont_slp(region(1)%num_groups, ngroup, ball_nn, ball_node)
endif

do ng = 1, ngroup
    if (inte_set%cont_region(1) == 1) then
        prop_cont_slp = inte_set%cont_slp(1:2,ng)
    else
        prop_cont_slp = inte_set%cont_slp(3:4,ng)
    endif
    !if ( mesh_set(1)%flag(ng) /= 4 ) then
        ! check the range of contact boundary
        if ( flag == 1 ) then
            call set_skin(1, ng, sp, adjn)
            allocate ( adj(adjn), this(adjn) )
            adj = region(1)%skin_nodes(sp:sp+adjn-1) 
            pos = inte_set%cont_slp(:,ng)
            do i = 1, 2
                call find_adj_seq(pos(i), adjn, adj, seq)
                pos(i) = seq
            enddo
        elseif ( flag == 2 ) then
            adjn = surface_set(ng)%adjn
            allocate (adj(adjn), this(adjn))
            adj = surface_set(ng)%adj(1:adjn)
            pos = 0
            do i = 1, 2
                temp_p = region(1)%node_coord(:,prop_cont_slp(i))
                min_dis = 10e8
                do j = 1, adjn
                    call calc_len(ball_node(:,adj(j)), temp_p, dis(1))
                    if ( dis(1) < 10e-6 ) then
                        pos(i) = j
                        exit
                    elseif ( dis(1) < min_dis ) then   
                        do k = 1, cn-1
                            call calc_len(target_node(:,cp(k)), target_node(:,cp(k+1)), cont_off)
                            call near_point(target_node(:,cp(k)), target_node(:,cp(k+1)), ball_node(:,adj(j)), np, cont_off*0.1, cont_check)
                            if ( cont_check ) then
                                temp_j = j
                                min_dis = dis(1)
                            endif
                        enddo
                    endif
                enddo
                if ( pos(i) == 0 ) pos(i) = temp_j
            enddo
        endif
        count = 0
        if ( pos(1) > pos(2) ) then
            do i = pos(1), adjn
                count = count + 1
                this(count) = i
            enddo
            do i = 1, pos(2)
                count = count + 1
                this(count) = i
            enddo
        else
            do i = pos(1), pos(2)
                count = count + 1
                this(count) = i
            enddo
        endif

        do j = 1, count
            mp(:) = ball_node(:,adj(this(j)))
            np = 0.0
            temp_k = 0
            min_dis = 10e8
            do k = 1, cn-1
                call calc_len(target_node(:,cp(k)), target_node(:,cp(k+1)), cont_off)
                call near_point(target_node(:,cp(k)), target_node(:,cp(k+1)), mp, np, cont_off*0.1, cont_check)
                if ( cont_check ) then
                    call calc_len(mp, np, np_dis)
                    if ( min_dis > np_dis ) then
                        temp_np = np
                        temp_k = k
                        min_dis = np_dis
                        temp_cont_off = cont_off
                    endif
                endif
            enddo
            
            if ( temp_k == 0 ) then
                write (*,*) 'Contact(move_contact_point): Cannot find contact node!!'
                do k = 1, cn-1
                    call near_point(target_node(:,cp(k)), target_node(:,cp(k+1)), mp, np, max_dis, cont_check)
                    if ( cont_check ) then
                        call calc_len(mp, np, np_dis)
                        if ( min_dis > np_dis ) then
                            temp_np = np
                            temp_k = k
                            min_dis = np_dis
                            temp_cont_off = cont_off
                        endif
                    endif
                enddo
                if ( temp_k == 0 ) then
                    write (*,*) 'sp, lp : ', prop_cont_slp, cp(1), cp(cn)
                    write (*,*) 'ng, j, adj : ', ng, this(j), adj(this(j))
                    stop 
                endif
            endif
            k = temp_k
            ball_node(:,adj(this(j))) = temp_np
            if ( j /= 1 .AND. j /= count ) then
                call calc_len(target_node(:,cp(k)), ball_node(:,adj(this(j))), dis(1))
                call calc_len(target_node(:,cp(k+1)), ball_node(:,adj(this(j))), dis(2))    
                if ( dis(1) <= temp_cont_off*0.005 .AND. k /= 1 ) then
                    call calc_angle(target_node(:,cp(k-1)), target_node(:,cp(k)), target_node(:,cp(k+1)), ang)
                    if ( ang <= cri_ang(1) .OR. ang >= cri_ang(2) ) ball_node(:,adj(this(j))) = target_node(:,cp(k))
                elseif ( dis(2) <= temp_cont_off*0.005 .AND. k /= cn-1 ) then
                    call calc_angle(target_node(:,cp(k)), target_node(:,cp(k+1)), target_node(:,cp(k+2)), ang)
                    if ( ang <= cri_ang(1) .OR. ang >= cri_ang(2) ) ball_node(:,adj(this(j))) = target_node(:,cp(k+1))
                endif
            endif
        enddo
        deallocate ( adj, this )
    !endif
enddo
deallocate(cp)

end subroutine move_contact_point


subroutine delete_cont_slp(ori_ngroup, ngroup, ball_nn, ball_node)

implicit none

integer, intent(in) :: ori_ngroup, ngroup, ball_nn
real, intent(in) :: ball_node(2,ball_nn)

integer :: i, j, k, adjn, num, temp_slp(4,ori_ngroup), group_chk(ngroup)
real :: dis, temp_p(2), tol
logical :: check
integer, allocatable :: adj(:)

do i = 1, ori_ngroup
    temp_slp(:,i) = inte_set%cont_slp(:,i)
enddo

num = 0
group_chk = 0
do i = 1, ori_ngroup
    tol = mesh_set(1)%min_d(i)
    temp_p = region(1)%node_coord(:,temp_slp(1,i))
    
    do j = 1, ngroup
        if (group_chk(j) == 0) then
            adjn = surface_set(j)%adjn
            allocate (adj(adjn))
            adj = surface_set(j)%adj(1:adjn)
            
            check = .false.
            do k = 1, adjn
                call calc_len(temp_p, ball_node(:,adj(k)), dis)
                if (dis < tol*0.1) then
                    num = num + 1
                    inte_set%cont_slp(:,j) = temp_slp(:,i)
                    group_chk(j) = 1
                    check = .true.
                    exit
                endif
            enddo
            if (check) exit
        endif
    enddo
enddo
        
end subroutine delete_cont_slp

!Contact routine end ===========================================================

!Coordinate interporation ========================================================
subroutine coord_interporation(Dim, rn, ng, nn, coord, thin_chk)

implicit none

integer, intent(in) :: Dim, rn, ng, nn
real, intent(inout) :: coord(Dim,nn)
logical, intent(in) :: thin_chk(ng)

integer :: i, sn, adjn, sp
integer, allocatable :: s_nodes(:), adj(:)
real, allocatable :: nndis(:)

if (Dim == 2) then
    !sn = inte_set%sn
    !allocate(s_nodes(sn))
    !s_nodes = inte_set%s_nodes
    sn = nn
    allocate(s_nodes(sn))
    do i = 1, nn
        s_nodes(i) = i
    enddo
    !do i = 1, ng
    !    !if (thin_chk(i)) then
    !        call set_skin(rn, i, sp, adjn)
    !        allocate(adj(adjn), nndis(adjn))
    !        adj = region(rn)%skin_nodes(sp:sp+adjn-1) 
    !        nndis = mesh_set(rn)%nndis(sp:sp+adjn-1)
    !
    !        call bound_smooth(adjn, adj, nndis, sn, s_nodes, nn, coord)
    !        deallocate ( adj, nndis )
    !    !endif
    !enddo
    call inner_smooth(rn, sn, s_nodes, nn, coord)
    deallocate(s_nodes)
elseif (Dim == 3) then
    !call surface_smooth_3D(rn, nn, coord)
    call inner_smooth_3D(rn, nn, coord)
endif

end subroutine coord_interporation
!Coordinate interporation end ===================================================

subroutine get_quad_coord( coord )

implicit none

real, intent(inout) ::  coord(2,4)
real :: quad(2,4)
real :: point 
real :: position(2,4)
real :: xi , eta 
integer :: i 

point = 1.0 / sqrt(3.0) 
position(1, :) = (/ -point, point, point ,-point /)
position(2, :) = (/ -point, -point, point ,point /)

do i=1,4
        xi = position(1, i)
        eta = position(2, i)

        quad(1,i) = ( (1-xi)*(1-eta)*coord(1,1) + (1+xi)*(1-eta)*coord(1,2)+ &
       (1+xi)*(1+eta)*coord(1,3)+(1-xi)*(1+eta)*coord(1,4) )*0.25
        quad(2,i) = ((1-xi)*(1-eta)*coord(2,1) + (1+xi)*(1-eta)*coord(2,2)+ &
       (1+xi)*(1+eta)*coord(2,3) + (1-xi)*(1+eta)*coord(2,4))*0.25
enddo

coord = quad

end subroutine 


subroutine check_contact_node(ng, prop_nn, prop_coord, unit_num_cont, cont_unit_num, flag)

implicit none

integer, intent(in) :: ng, prop_nn, unit_num_cont, cont_unit_num, flag
real, intent(in) :: prop_coord(2,prop_nn)

integer, allocatable :: adj(:), case_adj(:)
real, allocatable :: case_coord(:,:)

integer :: i, j, k, adjn, case_nn, case_adjn, count, cn, sp
integer :: temp(3), val(6)
real :: tol, dis, s, cp(2), temp_coord(2), temp_cont(4)
character(len=5) :: keyword
logical :: near_check, flag_check 
real, allocatable :: cont_coord(:,:)
logical, allocatable :: check(:)

case_nn = region(2)%num_nodes
call set_skin(2, 1, sp, case_adjn)
allocate ( case_coord(2,case_nn), case_adj(case_adjn) )
case_coord = region(2)%node_coord
case_adj = region(2)%skin_nodes(sp:sp+case_adjn-1) 
!tol = 0.02 * 0.03
tol = mesh_set(2)%ave_d(1) * 0.03

rewind (unit_num_cont)
do cn = 1, inte_set%cont_num
    read(unit_num_cont,*) keyword, val(:)
    count = 0
    if (flag == 1) then
        val(5) = val(5) + 1
    !elseif ( val(4) == 2 ) then
    !    do i = 1, val(5)
    !        if ( mesh_set(1)%flag(i) == 4 ) then
    !            count = count + 1
    !        endif
    !    enddo
    endif
    write (cont_unit_num,*) keyword, val(1:5), ng
    if ( val(4) == 2 ) then
        !ng = val(5)-count
        allocate (cont_coord(4,ng*2))
        cont_coord = 0.0
        do i = 1, ng
            flag_check = .TRUE.
            !if ( flag == 0 ) then
            !    if ( mesh_set(1)%flag(i) == 4 ) flag_check = .FALSE.
            !endif
            
            if ( flag_check ) then
                adjn = surface_set(i)%adjn       
                allocate ( adj(adjn), check(adjn) )
                adj = surface_set(i)%adj
                check = .FALSE.
                do j = 1, adjn
                    do k = case_adjn, 1, -1
                        if ( k == case_adjn ) then
                            temp(1:2) = (/ case_adj(k), case_adj(1) /)
                        else
                            temp(1:2) = (/ case_adj(k), case_adj(k+1) /)
                        endif
                        
                        call calc_len(case_coord(:,temp(1)), case_coord(:,temp(2)), dis)
                        call near_point(case_coord(:,temp(1)), case_coord(:,temp(2)), prop_coord(:,adj(j)), cp, dis*0.01, check(j))
                        if ( check(j) ) exit
                    enddo
                enddo
                
                count = 0
                do j = 1, adjn
                    temp = (/ j-1, j, j+1 /)
                    if ( j == 1 ) temp(1) = adjn
                    if ( j == adjn ) temp(3) = 1
                    
                    if ( check(temp(1)) == .FALSE. .AND. check(temp(2)) .AND. check(temp(3)) ) then  ! start point
                        count = count + 1
                        cont_coord(1:2,(i-1)*2+2) = prop_coord(:,adj(j))
                    elseif ( check(temp(1)) .AND. check(temp(2)) .AND. check(temp(3)) == .FALSE. ) then  ! last point
                        count = count + 1
                        cont_coord(3:4,(i-1)*2+2) = prop_coord(:,adj(j))
                    endif
                    if ( count == 2 ) exit
                enddo
                      
                do j = 1, 2
                    temp_coord = cont_coord(2*j-1:2*j,(i-1)*2+2)
                    do k = case_adjn, 1, -1
                        if ( k == case_adjn ) then
                            temp(1:2) = (/ case_adj(1), case_adj(k) /)
                        else
                            temp(1:2) = (/ case_adj(k+1), case_adj(k) /)
                        endif
                        call calc_len(case_coord(:,temp(1)), case_coord(:,temp(2)), dis)
                        call near_point(case_coord(:,temp(1)), case_coord(:,temp(2)), temp_coord, cp, dis*0.05, near_check)
                        
                        if ( near_check ) then
                            !if ( j == 1 ) then
                            !    call calc_len(case_coord(:,temp(1)), cp, dis)
                            !    if ( dis < tol*0.01) then
                            !        if ( k == case_adjn ) then
                            !            cont_coord(3:4,(i-1)*2+2) = case_coord(:,case_adj(2))
                            !        elseif ( k == case_adjn-1 ) then
                            !            cont_coord(3:4,(i-1)*2+2) = case_coord(:,case_adj(1))
                            !        else
                            !            cont_coord(3:4,(i-1)*2+2) = case_coord(:,case_adj(k+2))
                            !        endif
                            !    else
                            !        cont_coord(3:4,(i-1)*2+2) = case_coord(:,temp(1))
                            !    endif
                            !else
                            !    call calc_len(case_coord(:,temp(2)), cp, dis)
                            !    if ( dis < tol*0.01) then
                            !        if ( k == 1 ) then
                            !            cont_coord(1:2,(i-1)*2+2) = case_coord(:,case_adj(adjn))
                            !        else
                            !            cont_coord(1:2,(i-1)*2+2) = case_coord(:,case_adj(k-1))
                            !        endif
                            !    else
                            !        cont_coord(1:2,(i-1)*2+2) = case_coord(:,temp(2))
                            !    endif                                  
                            !endif
                            !exit
                            call projection_s_cont(case_coord(:,temp(1)), case_coord(:,temp(2)), cp, s)
                            if (s >= 0.0 .AND. s <= 1.0) then
                                if ( j == 1 ) then
                                    cont_coord(3:4,(i-1)*2+1) = case_coord(:,temp(1))
                                else
                                    cont_coord(1:2,(i-1)*2+1) = case_coord(:,temp(2))
                                endif
                                exit
                            endif
                        endif
                    enddo
                enddo
                deallocate ( adj, check )
            endif
        enddo
        if (ng == 2) then
            call calc_len(cont_coord(3:4,1), cont_coord(1:2,3), dis)
            if (dis < tol * 10e-6) then
                if (flag == 1) then
                    do k = 1, case_adjn
                        call calc_len(cont_coord(1:2,3), case_coord(:,case_adj(k)), dis)
                        if (dis <  tol * 10e-6) then
                            if (k == 1) then
                                cont_coord(1:2,3) = case_coord(:,case_adj(case_adjn))
                            else
                                cont_coord(1:2,3) = case_coord(:,case_adj(k-1))
                            endif
                            exit
                        endif
                    enddo
                else
                    read(unit_num_cont,*) temp_cont(:)
                    read(unit_num_cont,*) temp_cont(1:2), cont_coord(3:4,1)
                    read(unit_num_cont,*) temp_cont(:)
                    read(unit_num_cont,*) cont_coord(1:2,3), temp_cont(1:2)
                endif
            endif
        endif
        do i = 1, ng*2
            write (cont_unit_num,'(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)') cont_coord(:,i)
        enddo
        deallocate (cont_coord)
    else
        read(unit_num_cont,*) cont_coord(:,1)
        read(unit_num_cont,*) cont_coord(:,2)
        write (cont_unit_num,'(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)') cont_coord(:,1)
        write (cont_unit_num,'(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)') cont_coord(:,2)
    endif
enddo
rewind (unit_num_cont)

end subroutine check_contact_node

!=============================================================================================

subroutine write_inte_file(del_flag, num_rn, prop_nn, prop_coord, case_nn, case_coord, inte_unit_num)

implicit none

integer, intent(in) :: del_flag, num_rn, prop_nn, case_nn, inte_unit_num
real, intent(in) :: prop_coord(2,prop_nn), case_coord(2,case_nn)

integer :: i, j, temp_j, temp_inte(4)
real :: tol, dis, min_dis, temp_coord(2)

tol = 10e-6

temp_inte = 0
do i = 1, 2
    if (fsi_inte_set%fsi_inte_point(i*2-1) == 1) then
        temp_coord = fsi_inte_set%fsi_inte_coord(:,i)
        do j = 1, prop_nn
            call calc_len(prop_coord(:,j), temp_coord, dis)
            if (dis < tol) then
                temp_inte(i*2-1:i*2) = (/ 1, j /)
                exit
            endif
        enddo
        if (temp_inte(i*2-1) == 0) then
            if (del_flag == 1) then
                min_dis = 10e8
                do j = 1, case_nn
                    call calc_len(temp_coord, case_coord(:,j), dis)
                    if (min_dis > dis) then
                        temp_j = j
                        min_dis = dis
                    endif
                enddo
                temp_inte(i*2-1:i*2) = (/ num_rn, temp_j /)
            else
                write (*,*) 'surface_domain(write_inte_file) : cannot find inte_point!!!'
                write (*,*) 'coord : ', temp_coord
                stop
            endif
        endif
    else
        temp_inte(i*2-1:i*2) = fsi_inte_set%fsi_inte_point(i*2-1:i*2)
    endif
enddo
!write (*,*) fsi_inte_set%fsi_inte_point
!write (*,*) temp_inte
write (inte_unit_num,*) '*inte'
write (inte_unit_num,*) temp_inte
write (inte_unit_num,*) ' '
if (inte_set%cont_node(1) /= 0) then
    write (inte_unit_num,*) '*cont'
    write (inte_unit_num,*) inte_set%cont_node
    write (inte_unit_num,*) ' '
endif
if (inte_set%abla_node(1) /= 0) then
    write (inte_unit_num,*) '*abla'
    write (inte_unit_num,*) inte_set%abla_node
    write (inte_unit_num,*) ' '
endif
write (inte_unit_num,*) '*end'


end subroutine write_inte_file          
    

!=============================================================================================

subroutine build_cont_slp(unit_num, number_contacts)

implicit none

integer, intent(in) :: unit_num, number_contacts

character(len=5) :: keyword
integer :: id, reg_from, dim, reg_to, contact_type, i, j, k, status, num_contact_group = 1
real :: dis, tol
real, allocatable :: contact_pt_from(:,:,:), contact_pt_to(:,:,:)

! contact_pt_from : case
! contact_pt_to : propellant
tol = minval(mesh_set(1)%min_d(:))
write (*,*) 'tol:', tol
if (number_contacts > 0) then
  rewind(unit=unit_num)
  
  do  
    read(unit_num, *, iostat=status) keyword

    if (status == -1 .or. keyword == '*end ') then
      write(*,*) '*** EOF in contact info file ! -----------------------------'
      EXIT
    elseif (keyword == '*cont') then
      backspace(unit_num)
      read(unit_num,*) keyword, id, reg_from, reg_to, contact_type, dim, num_contact_group
      inte_set%cont_region = (/ reg_from, reg_to /)
      if (id > number_contacts) then
        STOP 'build_interface: contact id number error!'
      else
        allocate(contact_pt_from(2,2,num_contact_group))
        allocate(contact_pt_to(2,2,num_contact_group))
      endif
    
      do i = 1, num_contact_group
        if (contact_type == 2) then
            read (unit_num, *) contact_pt_from(1,1,i), contact_pt_from(2,1,i), contact_pt_from(1,2,i), contact_pt_from(2,2,i)
            do j = 1, 2
                do k = 1, region(reg_from)%num_nodes
                    call calc_len(contact_pt_from(:,j,i), region(reg_from)%node_coord(:,k), dis)
                    if ( dis < tol*10e-3 ) then
                        inte_set%cont_slp(j,i) = k
                        !inte_set%cont_slp(j+2,i) = k
                        exit
                    endif
                enddo
            enddo
            
            read (unit_num, *) contact_pt_to(1,1,i), contact_pt_to(2,1,i), contact_pt_to(1,2,i), contact_pt_to(2,2,i)
            do j = 1, 2
                do k = 1, region(reg_to)%num_nodes
                    call calc_len(contact_pt_to(:,j,i), region(reg_to)%node_coord(:,k), dis) 
                    if ( dis < tol*10e-3 ) then
                        inte_set%cont_slp(j+2,i) = k
                        !inte_set%cont_slp(j,i) = k
                        exit
                    endif
                enddo
            enddo
        endif
        write (*,'(A,4I8)') 'contact domain, nodes:', inte_set%cont_region(1), inte_set%cont_slp(1:2,i)
        write (*,'(A,4I8)') 'contact domain, nodes:', inte_set%cont_region(2), inte_set%cont_slp(3:4,i)
      end do
      deallocate(contact_pt_from, contact_pt_to)
    endif
  enddo
endif

end subroutine build_cont_slp


subroutine check_cons_condition(rn, unit_num)

implicit none

integer, intent(in) :: rn, unit_num

integer :: i, j, temp_i, count, sp, constr_pt, ndof, status, adjn, number, nn
integer :: convex_cons_num, line_cons_num, adj_count, num(2)
integer, allocatable :: temp_constr(:), constr_node(:), boundary(:,:), temp_boundary(:,:)
integer, allocatable :: adj(:), temp_adj(:), line_boundary(:,:), convex_boundary(:,:)
real :: temp_coord(2)
real, allocatable :: coord(:,:), line_cons(:,:), convex_cons(:,:)
character(len=5) :: keyword
logical :: check, corner_check, line_check

ndof = region(rn)%num_node_dof	
constr_pt = 0
rewind(unit_num)
do
    read(unit_num, *, iostat=status) keyword

    if (keyword == '*doma') then
        backspace(unit_num)
        read(unit_num,*) keyword, number
        if (number == 1) then
            do
                read(unit_num, *, iostat=status) keyword
	            if (keyword == '*cons') then
        		    backspace(unit_num)
        		    read(unit_num,*) keyword, constr_pt
        		    allocate(boundary(ndof,constr_pt), constr_node(constr_pt))
        		    do i= 1, constr_pt
			            read (unit_num, *) constr_node(i), (boundary(j,i), j= 1, ndof)
		            enddo
                    exit        		    
        		elseif (keyword == '*doma' .OR. keyword == '*mate' .OR. status == -1) then
        		    exit
        		endif
            enddo
            exit
        endif
    elseif (keyword == '*mate' .OR. status == -1) then
        exit
    endif
enddo

if (constr_pt /= 0) then
    call set_skin(rn, 1, sp, adjn)
    allocate (adj(adjn))
    adj = region(rn)%skin_nodes(sp:sp+adjn-1)
    
    nn = region(rn)%num_nodes
    allocate(coord(2,nn))
    coord = region(rn)%node_coord(1:2, :)
    
    allocate(temp_constr(constr_pt), temp_boundary(ndof,constr_pt), temp_adj(adjn))
    temp_constr = 0 
    do i = 1, adjn
        do j = 1, constr_pt
            if (adj(i) == constr_node(j)) then
                call find_corner_node(rn, constr_node(j), check)
                if (check) then
                    temp_constr(1) = constr_node(j)
                    temp_boundary(:,1) = boundary(:,j)
                    temp_adj(1) = adj(i)
                    temp_i = i
                    exit
                endif
            endif
        enddo
        if (temp_constr(1) /= 0) exit
    enddo
    
    count = 1
    adj_count = 1
    do i = temp_i+1, adjn
        do j = 1, constr_pt
            if (adj(i) == constr_node(j)) then
                count = count + 1
                temp_constr(count) = constr_node(j)
                temp_boundary(:,count) = boundary(:,j)
                exit
            endif
        enddo
        adj_count = adj_count + 1
        temp_adj(adj_count) = adj(i)
    enddo
    do i = 1, temp_i-1
        do j = 1, constr_pt
            if (adj(i) == constr_node(j)) then
                count = count + 1
                temp_constr(count) = constr_node(j)
                temp_boundary(:,count) = boundary(:,j)
                exit
            endif
        enddo
        adj_count = adj_count + 1
        temp_adj(adj_count) = adj(i)
    enddo
    constr_node = temp_constr
    boundary = temp_boundary
    adj = temp_adj
    deallocate (temp_constr, temp_boundary, temp_adj)
    
    convex_cons_num = 0
    line_cons_num = 0
    
    line_check = .FALSE.
    num = 1
    allocate (convex_cons(2,constr_pt), convex_boundary(2,constr_pt))
    allocate (line_cons(4,constr_pt), line_boundary(2,constr_pt))
    convex_cons = 0
    line_cons =0.0
    line_boundary = 0
    do 
        if (constr_node(num(1)) == adj(num(2))) then
            call find_corner_node(rn, constr_node(num(1)), corner_check)
            if (corner_check) then
                convex_cons_num = convex_cons_num + 1
                convex_cons(:,convex_cons_num) = coord(:,constr_node(num(1)))
                convex_boundary(:,convex_cons_num) = boundary(:,num(1))
                if (line_check == .FALSE.) then ! cons pt가 corner pt이면서 line cons의 시작점 후보
                    temp_coord = coord(:,constr_node(num(1)))
                else                       ! cons pt가 corner pt이면서 line cons의 마지막 지점 & 시작점 후보
                    line_cons(3:4,line_cons_num) = coord(:,constr_node(num(1))) 
                    temp_coord = coord(:,constr_node(num(1)))
                    line_check = .FALSE.
                endif                    
            elseif (corner_check == .FALSE.) then ! cons pt가 cons line인 지점 
                if (line_check == .FALSE.) then ! line cons의 두번째지점
                    line_cons_num = line_cons_num + 1
                    line_cons(1:2,line_cons_num) = temp_coord
                    line_boundary(:,line_cons_num) = boundary(:,num(1))
                    line_check = .TRUE.
                endif
            endif
            num(1) = num(1) + 1
            num(2) = num(2) + 1
        else
            line_check = .FALSE.
            num(2) = num(2) + 1
        endif
        
        if (num(1) > constr_pt .OR. num(2) > adjn) then
            if (line_check) then
                line_cons(3:4,line_cons_num) = coord(:,constr_node(1))
            endif
            exit
        endif
    enddo
    
    allocate (mesh_set(1)%convex_cons(2,convex_cons_num), mesh_set(1)%convex_boundary(2,convex_cons_num))
    allocate (mesh_set(1)%line_cons(4,line_cons_num), mesh_set(1)%line_boundary(2,line_cons_num))
    mesh_set(1)%convex_cons_num = convex_cons_num
    mesh_set(1)%line_cons_num = line_cons_num
    
    do i = 1, convex_cons_num
        !write (*,'(2F12.7,2X,2I5)') convex_cons(:,i), convex_boundary(:,i)
        mesh_set(1)%convex_cons(:,i) = convex_cons(:,i)
        mesh_set(1)%convex_boundary(:,i) = convex_boundary(:,i)
    enddo
    do i = 1, line_cons_num
        !write (*,'(4F12.7,2X,2I5)') line_cons(:,i), line_boundary(:,i)
        mesh_set(1)%line_cons(:,i) = line_cons(:,i)
        mesh_set(1)%line_boundary(:,i) = line_boundary(:,i)
    enddo
    
    deallocate (convex_cons, convex_boundary)
    deallocate (line_cons, line_boundary)
endif

deallocate (adj, coord)

end subroutine check_cons_condition




subroutine check_part_element(rn, ng, temp_adjn, temp_adj, limit_num, sp)

implicit none

integer, intent(in) :: rn, ng, limit_num, sp, temp_adjn
integer, intent(in) :: temp_adj(temp_adjn)

integer :: i, j, k, num, temp_ne, temp_nn, ne, nn, cp, se, ce
integer :: temp_np, next_np, adjn, ori_adjn
integer :: temp_ele(50), temp_nodes(100), temp_point(50), next_point(50)
integer :: inte_point(2), inte_seq(2), part_seq
integer, allocatable :: ele(:,:), adj(:), ori_adj(:)
integer, allocatable :: node_conn(:,:), node_cn(:), elem(:,:)
logical :: check

ne = region(rn)%num_elements
nn = region(rn)%num_nodes
allocate(ele(4,ne), node_conn(6,nn), node_cn(nn))
ele = region(rn)%element(4:7,:)
call set_skin(rn, ng, num, ori_adjn)
allocate ( ori_adj(ori_adjn) )
ori_adj = region(rn)%skin_nodes(sp:sp+ori_adjn-1)
        
temp_ne = 0 ;  temp_ele = 0
temp_nn = 0 ;  temp_nodes = 0
node_cn = 0 ;  node_conn = 0
do i = 1, ne
    do j = 1, 4
        node_cn(ele(j,i)) = node_cn(ele(j,i)) + 1
        node_conn(node_cn(ele(j,i)), ele(j,i)) = i
    enddo
enddo
  
temp_np = 1 ;  temp_point(1) = sp
temp_nn = 1 ;  temp_nodes(1) = sp
se = 1
! extract elements(minimum number over 20)
do
    do i = 1, temp_np
        cp = temp_point(i)
        num = node_cn(cp)
        do j = 1, num
            check = .TRUE.
            do k = 1, temp_ne
                if ( node_conn(j,cp) == temp_ele(k) ) then
                    check = .FALSE.
                    exit
                endif
            enddo
            if ( check ) then
                temp_ne = temp_ne + 1
                temp_ele(temp_ne) = node_conn(j,cp)
            endif
        enddo
    enddo
        
    next_np = 0 ;  next_point = 0
    do i = se, temp_ne
        ce = temp_ele(i)
        do j = 1, 4
            check = .TRUE.
            do k = 1, temp_np
                if ( ele(j,ce) == temp_point(k) ) then
                    check = .FALSE.
                    exit
                endif
            enddo
            if ( check ) then
                do k = 1, next_np
                    if ( ele(j,ce) == next_point(k) ) then
                        check = .FALSE.
                        exit
                    endif
                enddo
            endif
            if ( check ) then
                next_np = next_np + 1
                next_point(next_np) = ele(j,ce)
            endif
        enddo
    enddo
    temp_nodes(temp_nn+1:temp_nn+next_np) = next_point(1:next_np)
    temp_nn = temp_nn + next_np   
    if ( temp_ne > limit_num ) exit
        
    se = temp_ne + 1
    temp_np = next_np
    temp_point = next_point
enddo
!write (*,*) 'temp_ne(',temp_ne, '):', temp_ele(1:temp_ne)
!write (*,*) 'temp_nn(',temp_nn, '):', temp_nodes(1:temp_nn)
allocate(surface_set(ng)%part_ele(temp_ne))
allocate(surface_set(ng)%part_nodes(temp_nn))
surface_set(ng)%part_ne = temp_ne
surface_set(ng)%part_ele = temp_ele(1:temp_ne)
surface_set(ng)%part_nn = temp_nn
surface_set(ng)%part_nodes = temp_nodes(1:temp_nn)
    
allocate(elem(4,temp_ne))
allocate(adj(temp_nn))
elem = ele(:,temp_ele(1:temp_ne))
! extract skin of part elements
call write_adj(temp_nn, temp_nodes, temp_ne, elem, sp, adjn, adj)
deallocate(ele, node_conn, node_cn)
 
! check including first point in adj
surface_set(ng)%part_first = .FALSE.
do i = 1, adjn
    if ( adj(i) == region(rn)%node_group(ng) ) then
        surface_set(ng)%part_first = .TRUE.
        exit
    endif
enddo
            
! find start, last point at part interface
if (surface_set(ng)%part_first) then
    num = 0
    do i = 1, adjn
        check = .FALSE.
        do j = 1, ori_adjn
            if (adj(i) == ori_adj(j)) then
                check = .TRUE.
                exit
            endif
        enddo
        if (check == .FALSE. .AND. num == 0) then              
            num = num+1
            inte_point(num) = adj(i-1)
        elseif (check .AND. num == 1) then
            num = num+1
            inte_point(num) = adj(i)
            exit
        endif
    enddo
else
    num = 3
    do i = 1, ori_adjn
        check = .FALSE.
        do j = 1, adjn
            if (adj(j) == ori_adj(i)) then
                check = .TRUE.
                exit
            endif
        enddo
        if (check .AND. num == 3) then
            num = num-1
            inte_point(num) = ori_adj(i)
        elseif ( check == .FALSE. .AND. num == 2) then    
            num = num - 1
            inte_point(num) = ori_adj(i-1)
            exit
        endif
    enddo
endif

!num = 0
!do i = 1, temp_adjn
!    do j = 1, next_np
!        if ( temp_adj(i) == next_point(j) ) then
!            num = num + 1
!            inte_point(num) = temp_adj(i)
!            if ( num == 1 ) then
!                temp_i = i
!            elseif ( num == 2 ) then
!                if ( i-temp_i == 1 ) then
!                    num = num - 1
!                    inte_point(num) = temp_adj(i)
!                    temp_i = i
!                    do k = j, next_np
!                        next_point(k-1) = next_point(k)
!                    enddo
!                    next_np = next_np - 1
!                endif
!            endif
!            exit
!        endif
!    enddo
!    if ( num == 2 ) exit
!enddo
allocate(surface_set(ng)%part_inte_nodes(next_np))
do i = 1, adjn
    if ( adj(i) == inte_point(1) ) then
        part_seq = i
        exit
    endif
enddo
!write (*,*) 'inte_point:', inte_point
! find skin of part elements by using temp_adj
temp_np = 0 ;  temp_point = 0
! 1) find skin in temp_adj
do i = 1, temp_adjn
    if ( temp_adj(i) == inte_point(1) ) then
        exit
    else
        temp_np = temp_np + 1
        temp_point(temp_np) = temp_adj(i)
    endif
enddo         
! 2) find skin in part_inte_adj
num = 0
do i = part_seq, adjn
    temp_np = temp_np + 1
    temp_point(temp_np) = adj(i)
    num = num + 1
    surface_set(ng)%part_inte_nodes(num) = adj(i)
    if ( adj(i) == inte_point(2) ) then
        inte_seq(2) = temp_np
        exit
    endif
enddo
surface_set(ng)%part_ni = num
check = .FALSE.
! 3) find skin in temp_adj
do i = 1, temp_adjn
    if ( check ) then
        temp_np = temp_np + 1
        temp_point(temp_np) = temp_adj(i) 
    else
        if ( temp_adj(i) == inte_point(2) ) check = .TRUE.
    endif
enddo 
do i = 1, temp_np
    if ( temp_point(i) == inte_point(1) ) then
        inte_seq(1) = i
    elseif ( temp_point(i) == inte_point(2) ) then
        inte_seq(2) = i
        exit
    endif
enddo
allocate(surface_set(ng)%part_adj(temp_np))
surface_set(ng)%part_adjn = temp_np
surface_set(ng)%part_adj = temp_point(1:temp_np)
surface_set(ng)%part_inte_seq = inte_seq

deallocate(ori_adj, elem, adj)

end subroutine check_part_element


subroutine update_coarse_disp

implicit none

integer :: num_rn, rn, sn, i, m1, num
real :: disp

if (is_coarse_problem) then
    num_rn = num_of_regions
    do rn = 1, num_rn
        m1 = region(rn)%num_constr
        sn = num_disp_sets+rn
        do i = 1, m1
            num = region(rn)%disp_bc(i)
            disp = history(rn)%data(num,1)
            call update_constr(sn, num, disp)
        enddo
    enddo
endif

end subroutine update_coarse_disp




subroutine find_corner_node(rn, cp, check)

implicit none

integer, intent(in) :: rn, cp
logical, intent(out) :: check

integer :: cn1, cn2, pn, i

cn1 = mesh_set(rn)%corn_group(1)
if ( 1 == region(rn)%num_groups ) then
    cn2 = mesh_set(rn)%cn
else
    cn2 = mesh_set(rn)%corn_group(2) - 1
endif
pn = cn2 - cn1 + 1

check = .FALSE.
do i = cn1, cn2
    if (mesh_set(rn)%corn(i) == cp) then
        check = .TRUE.
        exit
    endif
enddo

end subroutine find_corner_node


subroutine update_corner_info(rn)

implicit none

integer, intent(in) :: rn

integer :: i, j, k, nn, ng, sp, cn1, cn2, temp_cn, adjn, memory
real :: dis(2)
integer, allocatable :: adj(:)
real, allocatable :: node(:,:)

nn = region(rn)%num_nodes
ng = region(rn)%num_groups
allocate(node(2,nn))
node = region(rn)%node_coord

do i = 1, ng
    cn1 = mesh_set(rn)%corn_group(ng)
    if ( ng == region(rn)%num_groups ) then
        cn2 = mesh_set(rn)%cn
    else
        cn2 = mesh_set(rn)%corn_group(ng+1) - 1
    endif
    temp_cn = cn2 - cn1 + 1
    
    call set_skin(1, i, sp, adjn)
    allocate (adj(adjn))
    adj = region(1)%skin_nodes(sp:sp+adjn-1)
    
    memory = 1
    do j = cn1, cn2
        do k = memory, adjn
            if (mesh_set(rn)%corn(j) == adj(k)) then
                if (k == 1) then
                    call calc_len(node(:,adj(k)), node(:,adj(adjn)), dis(1))
                    call calc_len(node(:,adj(k)), node(:,adj(k+1)), dis(2))
                elseif (k == adjn) then
                    call calc_len(node(:,adj(k)), node(:,adj(k-1)), dis(1))
                    call calc_len(node(:,adj(k)), node(:,adj(1)), dis(2))
                else
                    call calc_len(node(:,adj(k)), node(:,adj(k-1)), dis(1))
                    call calc_len(node(:,adj(k)), node(:,adj(k+1)), dis(2))
                endif
                
                if (j == cn1) then
                    mesh_set(rn)%edge(2,cn2) = dis(1)
                    mesh_set(rn)%edge(1,j) = dis(2)
                else
                    mesh_set(rn)%edge(2,j-1) = dis(1)
                    mesh_set(rn)%edge(1,j) = dis(2)
                endif
                memory = k+1
            endif
        enddo
    enddo
    deallocate (adj)
enddo
mesh_set(rn)%corn_coord = node(:,mesh_set(rn)%corn)

end subroutine update_corner_info


subroutine free_surface_domain

if (allocated(inte_set%cont_slp)) deallocate(inte_set%cont_slp)
if (allocated(inte_set%s_nodes)) deallocate(inte_set%s_nodes)
if (allocated(inte_set%c_nodes)) deallocate(inte_set%c_nodes, inte_set%c_elem)
if (allocated(move_set)) deallocate(move_set)
deallocate(mesh_set)
deallocate(surface_set)
deallocate(remesh_flag)
inte_set%surf_num = 0
if (allocated(inte_set%surf_nodes)) then
    deallocate (inte_set%surf_nodes, inte_set%surf_nodes_info)
    deallocate (inte_set%start_nodes, inte_set%end_nodes)
endif

end subroutine free_surface_domain

!============================================================================

subroutine set_domain_3D(current_domain, max_nodes, max_elements)

use file_info_house

implicit none

integer, intent(in) :: current_domain, max_nodes, max_elements
integer :: status, i, num, count, val
character(len=5) :: keyword

mesh_set(current_domain)%nn = max_nodes
mesh_set(current_domain)%ne = max_elements
allocate (mesh_set(current_domain)%node_coord(3,max_nodes))
allocate (mesh_set(current_domain)%element(8,max_elements))

! read node
read_nodes: do
    read(unit_num_node, *, iostat=status) keyword
    if (status == -1 .or. keyword == '*end ') then
		EXIT
    elseif (keyword == '*node') then
        backspace(unit_num_node)
        read(unit_num_node, *, iostat=status) keyword, val
        if (val == current_domain) then
            count = 0
            do i = 1, max_nodes
                count = count + 1
                read(unit_num_node, *, iostat=status) num, mesh_set(current_domain)%node_coord(:,count)
            enddo
            exit read_nodes
        endif
    endif
enddo read_nodes
            
! read element
read_elements: do
    read(unit_num_elem, *, iostat=status) keyword
    if (status == -1 .or. keyword == '*end ') then
		EXIT
    elseif (keyword == '*elem') then
        backspace(unit_num_elem)
        read(unit_num_elem, *, iostat=status) keyword, val
        if (val == current_domain) then
            count = 0
            do i = 1, max_elements
                count = count + 1
                read(unit_num_elem, *, iostat=status) num, mesh_set(current_domain)%element(:,count)
            enddo
            exit read_elements
        endif
    endif
enddo read_elements

end subroutine set_domain_3D


end module surface_domain