subroutine nodal_interpolation2(Dim, mat_num, mat_ne, cur_time, num_var, quad_data, tec_unit_num, ier, flag )

use time_domain
use physical_domain
use system_domain
use output_domain
    
implicit none
INCLUDE 'tecio.f90'

integer, intent(in) :: Dim, mat_num, mat_ne, num_var, tec_unit_num, flag
integer, intent(inout) :: ier
real, intent(in) :: cur_time, quad_data(num_var,mat_ne*(Dim-1)*4)

integer :: i, j, nn, ne, rn_nn, rn, cp, ele_node
integer ::  el_data((Dim-1)*4)
real :: temp, coe(8,3), u_cartesian(Dim)

integer, allocatable :: node_count(:), elem(:,:), conn(:), node_table(:)
real, allocatable :: variable(:,:), new_var(:,:), poly_mat(:,:)
real*4, allocatable :: coord(:,:), node_value(:,:), disp(:,:), temp_val(:)

INTEGER, allocatable :: ValueLocation(:)
INTEGER :: disdouble
CHARACTER(len=8) :: zone_name
CHARACTER(1) :: NulChar = CHAR(0)
INTEGER :: Null(*)
POINTER    (nullptr,Null)

!======================================================================
nn = mat_mesh_set(mat_num)%nn
ne = mat_mesh_set(mat_num)%ne
rn = mat_mesh_set(mat_num)%rn
rn_nn = mat_mesh_set(mat_num)%rn_nn 
allocate (coord(Dim,nn), elem((Dim-1)*4,ne), conn(nn), node_table(nn))
allocate (disp(Dim,nn), temp_val(nn))
conn = mat_mesh_set(mat_num)%inv_conn
elem = mat_mesh_set(mat_num)%elem
coord = region(rn)%node_coord(:,conn)
if (allocated(sdomain(rn)%u_cartesian)) then
    do i = 1, nn
        !disp(:,i) = history(rn)%data(region(rn)%node(4:3+Dim,conn(i)),1) 
        disp(:,i) = sdomain(rn)%u_cartesian(region(rn)%node(4:Dim+3,conn(i)))
        !disp(:,i) = sdomain(rn)%u_cartesian(sdomain(rn)%dof_index(region(rn)%node(4:Dim+3,conn(i))))
        !disp(:,i) = 0.0
    enddo
else
    disp = 0.0
endif
do i = 1, nn
    node_table(region(rn)%node(1,i)) = i
enddo

if (Dim == 2) then
    allocate (poly_mat(4,4), ValueLocation(14))
    temp = 0.5*1.732050808
    poly_mat(:,1) = (/1+temp, -0.5, 1-temp, -0.5/)
    poly_mat(:,2) = (/-0.5, 1+temp, -0.5, 1-temp/)
    poly_mat(:,3) = (/1-temp, -0.5, 1+temp, -0.5/)
    poly_mat(:,4) = (/-0.5, 1-temp, -0.5, 1+temp/)
    ele_node = 4
elseif (Dim == 3) then
    allocate(poly_mat(8,8), ValueLocation(21))
    temp = sqrt(3.0)
    coe(:,1) = (/ -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0 /)
    coe(:,2) = (/ -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0 /)
    coe(:,3) = (/ -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 /)
    do i = 1, 8
        do j = 1, 8
            poly_mat(j,i) = 0.125*(1+coe(i,1)*coe(j,1)*temp)*(1+coe(i,2)*coe(j,2)*temp)*(1+coe(i,3)*coe(j,3)*temp)
        enddo
    enddo
    ele_node = 8
endif
allocate (node_value(nn,num_var),node_count(nn))
allocate (variable(ele_node,num_var), new_var(ele_node,num_var))
node_count = 0
node_value = 0.0
    
do i = 1, ne
    el_data = elem(:,i)
    node_count(el_data) = node_count(el_data)+1
    do j = 1, ele_node
		variable(j,:) = quad_data(:,ele_node*(i-1)+j)
    enddo
    new_var = matmul(poly_mat, variable)

    node_value(el_data,:) = node_value(el_data,:) + new_var
enddo

do i = 1, nn
    if ( node_count(i) == 0 ) then
        stop 'node_interporation error : count of node is 0'
    else
        node_value(i,:) = node_value(i,:) / float(node_count(i))
    endif
enddo

515 format (2F12.7,2X,100ES15.7)
516 format (8(I10,1X))
if (flag == 1) then
    if (Dim == 2) then
        Write (tec_unit_num,'(A,E14.7,A,I1,A,I5,A,I5, A)') 'ZONE, SolutionTime = ', cur_time, ', T = "Mat_num', mat_num, '", N, =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
        do i = 1, nn
            cp = node_table(i)
            write(tec_unit_num,515) coord(2,cp), coord(1,cp), history(rn)%data(region(rn)%node(5,conn(cp)),1), history(rn)%data(region(rn)%node(4,conn(cp)),1), ( node_value(cp,j), j=1,num_var )
        enddo
    elseif (Dim == 3) then
        Write (tec_unit_num,'(A,E14.7,A,I1,A,I5,A,I5, A)') 'ZONE, SolutionTime = ', cur_time, ', T = "Mat_num', mat_num, '", N, =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
        do i = 1, nn
            cp = node_table(i)
            !history_data = history(rn)%data(region(rn)%node(4:3+Dim,conn(cp)),1)
            !u_cartesian = sdomain(rn)%u_cartesian(sdomain(rn)%dof_index(region(rn)%node(4:Dim+3,conn(cp))))
            u_cartesian = sdomain(rn)%u_cartesian(region(rn)%node(4:Dim+3,conn(cp)))
            write(tec_unit_num,515) coord(:,cp), u_cartesian, ( node_value(cp,j), j=1,num_var )
        enddo
    endif
    
    do i = 1, ne
        write(tec_unit_num,516) elem(:,i)
    enddo
elseif (flag == 2) then
    ! teczne112
    nullptr = 0
    zone_name = 'mat_num0'
    write (zone_name(8:8), '(I1.1)') mat_num
    if (Dim == 2) then
        ValueLocation = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
        ier = teczne112(Trim(zone_name)//nulchar,3,nn,ne,1,0,0,0,cur_time,1,0,1,0,0,0,0,0,Null,ValueLocation,Null,0)
    elseif (Dim == 3) then
        ValueLocation = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
        ier = teczne112(Trim(zone_name)//nulchar,5,nn,ne,1,0,0,0,cur_time,1,0,1,0,0,0,0,0,Null,ValueLocation,Null,0)
    endif
    
    ! tecdat112
    Disdouble = 0
    if (Dim == 2) then
        temp_val = coord(2,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
        temp_val = coord(1,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
        temp_val = disp(2,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
        temp_val = disp(1,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
    elseif (Dim == 3) then
        ! tecdat112
        temp_val = coord(1,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
        temp_val = coord(2,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
        temp_val = coord(3,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
        temp_val = disp(1,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
        temp_val = disp(2,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
        temp_val = disp(3,node_table(1:nn))
        ier = tecdat112(nn,temp_val,disdouble)
    endif
    do i = 1, num_var
        temp_val = node_value(node_table(1:nn),i)
        ier = tecdat112(nn,temp_val,disdouble)
    enddo
    ! tecnod112
    ier = tecnod112(elem)
endif

deallocate (coord, elem, disp, conn, poly_mat, ValueLocation, node_table)
deallocate (node_value, node_count, variable, new_var, temp_val)

end subroutine nodal_interpolation2