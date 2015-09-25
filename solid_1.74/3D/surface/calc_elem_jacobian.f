subroutine calc_elem_jacobian(nn, ne, coord, elem, det_jacobi)
    
implicit none

integer, intent(in) :: nn, ne, elem(8,ne)
real, intent(in) :: coord(3,nn)
real, intent(inout) :: det_jacobi

integer :: i, j, k, q, sfn(8,3), this(2), mesh_node_num(nn)
real :: M33DET, det_j, jacobi(3,3), ele_coord(3,8)
real :: dis(3), tol_dis, node_det_j(nn)

sfn(:,1) = (/ -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0 /)
sfn(:,2) = (/ -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0 /)
sfn(:,3) = (/ -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 /)
mesh_node_num = 0
!node_det_j = 10e8
node_det_j = 0.0
det_jacobi = 1.0
do i = 1, ne
    ele_coord = coord(:,elem(:,i))
    jacobi = 0.0
    do j = 1, 3
        do k = 1, 3
            do q = 1, 8
                jacobi(j,k) = jacobi(j,k) + sfn(q,k)*ele_coord(j,q)
            enddo
        enddo
    enddo
    jacobi = jacobi/8.0
    det_j = M33DET(jacobi)
    tol_dis = 0.0
    do j = 1, 4
        this = (/ j, j+1 /)
        if (j == 4) this = (/ j, 1 /)
        call calc_length_3D(ele_coord(:,this(1)), ele_coord(:,this(2)), dis(1))
        call calc_length_3D(ele_coord(:,this(1)+4), ele_coord(:,this(2)+4), dis(2))
        call calc_length_3D(ele_coord(:,this(1)), ele_coord(:,this(1)+4), dis(3))
        tol_dis = tol_dis + sum(dis)
    enddo
    det_j = abs(det_j * (2.0/(tol_dis/12.0))**3)
    do j = 1, 8
        mesh_node_num(elem(j,i)) = mesh_node_num(elem(j,i)) + 1
        node_det_j(elem(j,i)) = node_det_j(elem(j,i)) + det_j
    enddo
    if (det_jacobi > det_j) det_jacobi = det_j
enddo
do i = 1, nn
    if (mesh_node_num(i) == 0) then
        write (*,'(A,I7,A)') 'Error(calc_ele_jacobian): Not include ', i, ' node in mesh'
        stop
    endif
    node_det_j(i) = node_det_j(i)/float(mesh_node_num(i))
    !if (node_det_j(i) > 1.0) write (*,'(A,I6,A,F12.7,A,F12.7)') 'jacobian determinant of ', i, ' element: ', det_j, ', mean_dis:', tol_dis/12.0
enddo

call check_mesh_quality(nn, coord, ne, elem, node_det_j)


end subroutine calc_elem_jacobian
    
    
FUNCTION M33DET (A) RESULT (DET)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
DOUBLE PRECISION :: DET

DET =   A(1,1)*A(2,2)*A(3,3)  &
    - A(1,1)*A(2,3)*A(3,2)  &
    - A(1,2)*A(2,1)*A(3,3)  &
    + A(1,2)*A(2,3)*A(3,1)  &
    + A(1,3)*A(2,1)*A(3,2)  &
    - A(1,3)*A(2,2)*A(3,1)

RETURN

END FUNCTION M33DET


subroutine check_mesh_quality(nn, node, ne, ele, det_j)

implicit none

integer :: nn, ne, ele(8,ne)
real :: node(3,nn), det_j(nn)

integer :: i, status

open (Unit=37, File='check_mesh_quality.plt', STATUS='replace', ACTION='write', IOSTAT=status)
Write (37,'(A)') 'TITLE="check_proj_num!!!"'
Write (37,*) 'VARIABLES="x", "y", "z", "det_jacobian"'
Write (37,'(A,I6,A,I6, A)') 'ZONE N=', nn , ', E=', ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'        
do i = 1, nn
	write(37,'(3(F12.7,2X),3x,F12.7)') node(:,i), det_j(i)
enddo
do i = 1, ne
    write (37,'(8(I7,2X))') ele(:,i)
enddo
close(37)

end subroutine check_mesh_quality
