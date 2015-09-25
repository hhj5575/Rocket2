subroutine remesh_3D_main(ng, cur_ng, temp_arr, tri_nn, tri_ne, tri_coord, tri_elem, tri_np2, tri_p2, &
                          ridge_group_num, ridge_num, ridge, ridge_conn, new_nn, new_ne, new_coord, new_elem)
    
use remesh_domain
implicit none

integer, intent(in) :: ng, cur_ng, temp_arr, tri_nn, tri_ne, tri_elem(3,tri_ne), tri_np2, tri_p2(tri_np2)
integer, intent(in) :: ridge_group_num, ridge_num(1:ridge_group_num), ridge(500,ridge_group_num), ridge_conn(ridge_group_num)
real,    intent(in) :: tri_coord(3,tri_nn)
integer, intent(inout) :: new_nn, new_ne, new_elem(8,temp_arr)
real,    intent(inout) :: new_coord(3,temp_arr)

integer :: nn, ne, ln, ned, ele_nodes, flag, np2, fsi_flag
integer, allocatable :: elem(:,:), temp_line(:,:), line(:,:), edge(:), p2(:,:), cen(:), control_edge(:,:)
real, allocatable :: node(:,:), line_ang(:)

integer :: i, j, k, q, temp_k, count, tln, remesh_file_num, status, a, tcen, n, add_nn
integer :: this(2), in(2)
real :: cri_ang
real :: p(3,3), pl(4), add_coord(3,20)
logical :: check, existence
character(len=50) :: filename

if (allocated(Tri%node)) call free_geo_domain(1)
allocate(Tri%node(3,tri_nn), Tri%elem(4,tri_ne))
Tri%nn = tri_nn
Tri%ne = tri_ne
Tri%node = tri_coord
Tri%elem = tri_elem

write (*,*) 'tri_nn:', tri_nn, ', tri_ne:', tri_ne
filename = 'check_tri_surface.plt'
call tecplot_surface(filename, tri_nn, tri_ne, tri_coord, tri_elem)
write (*,*) 'tri_np2:,', tri_np2
write (*,*) tri_p2
!call alloc_mesh_info(tri_nn, tri_ne, tri_coord, tri_elem)

allocate(temp_line(4,tri_ne*3))
tln = 0;   temp_line = 0
do i = 1, tri_ne
    do j = 1, 3
        if (j == 3) then
            this = (/ tri_elem(3,i), tri_elem(1,i) /)
        else
            this = (/ tri_elem(j,i), tri_elem(j+1,i) /)
        endif
        
        temp_k = 0
        do k = 1, tln
            count = 0
            do q = 1, 2
                if (temp_line(q,k) == this(1)) count = count + 1
                if (temp_line(q,k) == this(2)) count = count + 1
            enddo
            if (count == 2) then
                temp_k = k
                exit
            endif
        enddo
        if (temp_k /= 0) then
            temp_line(4,temp_k) = i
        else
            tln = tln + 1
            temp_line(1:2,tln) = this
            temp_line(3,tln) = i
        endif
    enddo
enddo
ln = tln
allocate (line(4,ln), line_ang(ln), edge(ln), p2(2,ln)) 
allocate (cen(ln), control_edge(ln,200))
line = 0;  line_ang = 0.0;  edge = 0;  p2 = 0;  ned = 0
tcen = 0;  cen = 0;  control_edge = 0
line = temp_line(:,1:ln)
deallocate (temp_line)

! fsi_flag = 1: make boundary data
!            2: get boundary data from surface module
fsi_flag = 2
if (fsi_flag == 1) then
    line_ang= 0.0
    ele_nodes = 3
    call calc_line_angle(ele_nodes, tri_nn, tri_ne, ln, tri_coord, tri_elem, line, line_ang)
    write (*,*) '>> complete calculation line angle'
    if (tri_np2 > 0) then
        write (*,*) 'np2:', tri_np2
        write (*,*) tri_p2(1:tri_np2)
    endif

    cri_ang = 24.
    do i = 1, 11
        ned = 0;  np2 = 0;  p2 = 0
        call remesh_find_edge_bak(ln, line, line_ang, ned, np2, edge, p2, tri_np2, tri_p2, cri_ang, check)
        !call remesh_find_edge(ln, line, line_ang, tri_nn, tri_coord, tri_ne, tri_elem, ned, tri_np2, edge, tri_p2, cri_ang, check)
        if (check) then
            write (*,*) '>> complete find edge line'
            exit
        else
            write (*,*) 'fail finding edge when criteria angle is', cri_ang
            if (i == 11) stop
            cri_ang = 30. - 2.*i
        endif
    enddo
elseif (fsi_flag == 2) then
    np2 = tri_np2
    do i = 1, np2
        p2(1,i) = tri_p2(i)
    enddo
    tcen = ridge_group_num
    allocate (Hexa(cur_ng)%b_edge_conn(tcen))
    do i = 1, ridge_group_num
        cen(i) = ridge_num(i)
        control_edge(1:cen(i),i) = ridge(1:cen(i),i)
        Hexa(cur_ng)%b_edge_conn(i) = ridge_conn(i)
    enddo
    !write (*,*) 'ln:', ln
    !write (*,*) 'np2:', np2
    !write (*,'(100I8,1X)') p2(1,1:np2)
    !write (*,*) tcen
    !write (*,'(100I8,1X)') cen(1:tcen)
    !do i = 1, tcen
    !    write (*,'(/,(10I9))') control_edge(1:cen(i),i)
    !enddo
endif

call add_data_for_remesh(ng, np2, p2(:,1:np2), add_nn, add_coord)
write (*,*) '>> remeshing: start setting boundary representaion'
call set_b_rep(ng, cur_ng, tri_nn, tri_coord, tri_ne, tri_elem, add_nn, add_coord, ln, line, line_ang, np2, p2, ned, edge, tcen, cen, control_edge, fsi_flag)
write (*,*) '>> remeshing: complete setting boundary representaion'
deallocate (line, line_ang, edge, p2)
deallocate (cen, control_edge)

call make_hexa_mesh(ng, cur_ng, add_nn, add_coord, fsi_flag)
write (*,*) '>> remeshing: complete making hexa mesh'

new_coord(:,new_nn+1:new_nn+Hexa(cur_ng)%nn) = Hexa(cur_ng)%node(:,1:Hexa(cur_ng)%nn)
new_elem(:,new_ne+1:new_ne+Hexa(cur_ng)%ne) = Hexa(cur_ng)%elem(:,1:Hexa(cur_ng)%ne)+new_nn
new_nn = new_nn+Hexa(cur_ng)%nn
new_ne = new_ne+Hexa(cur_ng)%ne

end subroutine remesh_3D_main
!=====================================================================================




