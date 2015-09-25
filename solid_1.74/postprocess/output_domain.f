module output_domain 
    
implicit none
save
    
type mat_mesh
    integer :: nn, ne, rn, rn_nn, mat_num, ele_num
    integer, allocatable :: elem(:,:), node(:), conn(:), inv_conn(:), elem_num(:)
end type

type(mat_mesh), allocatable :: mat_mesh_set(:)

integer :: mat_mesh_set_num

contains
!====================================================================================
subroutine set_mat_mesh(unit_num_elem, unit_num_node, num_rn, num_materials)

implicit none

integer, intent(in) :: unit_num_elem, unit_num_node, num_rn, num_materials

integer :: i, j, rn, status, mat_num, num, mat_nn, mat_ne, temp_ne, temp_num, ele_num
integer :: elem_nodes
integer, allocatable :: node_num(:), temp_node_num(:), nn(:), temp_elem(:,:), temp_elem_num(:)
!real, allocatable :: coord(:,:,:)
character(len=5) :: keyword
character(len=3) :: elem_str

allocate (nn(num_rn))
rewind(unit=unit_num_node)
do i = 1, num_rn
    nn(i) = 0
    read(unit_num_node,*, iostat=status) keyword
    if (keyword == '*node') then
        do
            read(unit_num_node,*, iostat=status) keyword
            backspace(unit_num_node)
            if (status == -1 .or. keyword == '*node' .or. keyword == '*end ') then
                exit
            else
                read(unit_num_node,*, iostat=status) num
                nn(i) = nn(i) + 1
            endif
        enddo
    endif
enddo

allocate(mat_mesh_set(num_materials))
mat_mesh_set_num = num_materials
rewind(unit=unit_num_elem)
temp_ne = 0
do
    read(unit_num_elem,*, iostat=status) keyword
    temp_ne = temp_ne + 1
    if (status == -1) exit
enddo
allocate(temp_elem_num(temp_ne))
rewind(unit=unit_num_elem)

i = 0
find_mat_elem: do
    read(unit_num_elem,*, iostat=status) keyword
    if (keyword == '*elem' .or. keyword == '*spri') then
        backspace(unit_num_elem)
        read(unit_num_elem,*, iostat=status) keyword, rn, ele_num, mat_num
        i = i + 1
        mat_mesh_set(i)%mat_num = mat_num
        mat_mesh_set(i)%ele_num = ele_num
        mat_ne = 0
        
        if (ele_num >= 400 .or. ele_num <= 200) then
            elem_nodes = 2
        else
            write (elem_str,'(I3)') ele_num
            if (elem_str(3:3) == '1' .or. elem_str(3:3) == '2') then
                elem_nodes = 4
            elseif (elem_str(3:3) == '3') then
                elem_nodes = 8
            endif
        endif
        allocate (temp_elem(elem_nodes,temp_ne))
        do
            read(unit_num_elem,*, iostat=status) keyword
            backspace(unit_num_elem)
            if (status == -1 .or. keyword == '*elem' .or. keyword == '*end ' .or. keyword == '*spri') then
                if (keyword == '*elem' .or. keyword == '*spri') backspace(unit_num_elem)
                EXIT   ! exit from do-while loop
            else
                mat_ne = mat_ne + 1
                read(unit_num_elem,*, iostat=status) temp_elem_num(mat_ne), (temp_elem(j,mat_ne), j=1,elem_nodes)
            endif
        enddo
        
        mat_mesh_set(i)%rn = rn
        mat_mesh_set(i)%ne = mat_ne
        allocate (mat_mesh_set(i)%elem(elem_nodes,mat_ne), mat_mesh_set(i)%elem_num(mat_ne))
        allocate (node_num(mat_ne*elem_nodes), temp_node_num(mat_ne*elem_nodes))
        
        do j = 1, mat_ne
            mat_mesh_set(i)%elem_num(j) = temp_elem_num(j)
            node_num((j-1)*elem_nodes+1:j*elem_nodes) = temp_elem(1:elem_nodes,j)
        enddo
        call sortc(node_num, mat_ne*elem_nodes)

        temp_num = 0
        mat_nn = 0
        do j = 1, mat_ne*elem_nodes
            if (temp_num /= node_num(j)) then
                mat_nn = mat_nn + 1
                temp_node_num(mat_nn) = node_num(j)
                temp_num = node_num(j)
            endif
        enddo
        mat_mesh_set(i)%nn = mat_nn
        mat_mesh_set(i)%rn_nn = nn(rn)
        !allocate (mat_mesh_set(i)%node(mat_nn), mat_mesh_set(i)%coord(2,mat_nn))
        allocate (mat_mesh_set(i)%conn(nn(rn)), mat_mesh_set(i)%inv_conn(mat_nn))
        !mat_mesh_set(i)%node = temp_node_num(1:mat_nn)
                
        mat_mesh_set(i)%conn = 0
        mat_mesh_set(i)%inv_conn = 0
        do j = 1, mat_nn
            !mat_mesh_set(i)%coord(:,j) = coord(:,temp_node_num(j),rn)
            mat_mesh_set(i)%conn(temp_node_num(j)) = j
            mat_mesh_set(i)%inv_conn(j) = temp_node_num(j)
        enddo
                
        do j = 1, mat_ne
            mat_mesh_set(i)%elem(1:elem_nodes,j) = mat_mesh_set(i)%conn(temp_elem(1:elem_nodes,j))
        enddo
                
        !write (*,'(A,5I5)') 'mat_num, ele_num, rn, mat_nn, mat_ne:', ele_num, mat_num, rn, mat_nn, mat_ne
        !write (*,*) 'mat_elem:', temp_elem(1:elem_nodes,mat_ne)
        deallocate (node_num, temp_node_num, temp_elem)
    elseif (status == -1) then
        exit find_mat_elem        
    endif
enddo find_mat_elem
deallocate(temp_elem_num)

end subroutine set_mat_mesh

subroutine rearrange_mesh(nn, coord, ne, elem, mat_num, mat_nn, mat_coord, mat_ne, mat_elem, conn)

implicit none

integer, intent(in) :: nn, ne, elem(4,ne), mat_num
real, intent(in) :: coord(2,nn)
integer, intent(inout) :: mat_nn, mat_ne, mat_elem(4,mat_ne), conn(nn)
real, intent(inout) :: mat_coord(2,mat_nn)

integer :: i
integer :: node_num(mat_nn), temp_ele(4,mat_ne)

node_num = mat_mesh_set(mat_num)%node
temp_ele = mat_mesh_set(mat_num)%elem
conn = 0
do i = 1, mat_nn
    mat_coord(:,i) = coord(:,node_num(i))
    conn(node_num(i)) = i
enddo
do i = 1, mat_ne
    mat_elem(:,i) = conn(elem(:,i))
enddo

end subroutine rearrange_mesh


subroutine free_output_domain

deallocate(mat_mesh_set)

end subroutine free_output_domain


end module