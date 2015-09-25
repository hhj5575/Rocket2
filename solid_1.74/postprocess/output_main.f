subroutine output_interpolation(Dim, unit_out_num, tol_rn, mode, output_region2, cur_time, current_i, output_interval, flag)

! flag = 1: ASCII Data
!        2: Binary Data
    
use physical_domain
use output_domain

implicit none
INCLUDE 'tecio.f90'

integer, intent(in) :: Dim, unit_out_num, tol_rn, mode, output_region2, current_i, output_interval, flag
real, intent(in) :: cur_time

real, allocatable :: quad_data(:,:), temp_quad_data(:,:)

integer :: rn, ne, tol_nn, tol_ne,  num_var, vec_unit_num, tec_unit_num, mat, ele_num
integer :: max_mat_num, i, error, status, count, num, st_i, ed_i, cur_ne, mat_ne
character(len=40) :: fn
character(len=200) :: all_var
character(len=8) :: a(8)
character(len=8), allocatable :: var(:)
CHARACTER(200) :: VarsName
CHARACTER(1) :: NulChar = CHAR(0)
INTEGER :: Filetype, Debug, Visdouble
INTEGER :: ier

if (mod(current_i,abs(output_interval)) == 0) then  !----------------------------------------------

vec_unit_num = 113
tec_unit_num = 114

tol_nn = 0;  tol_ne = 0
rewind(unit_out_num)
do rn = 1, tol_rn
    if ( output_region2 < 1 .OR. output_region2 == rn ) then
        ne = region(rn)%num_elements
        tol_nn = tol_nn + region(rn)%num_nodes
        tol_ne = tol_ne + ne
    endif
enddo

!write (*,*) 'Number of node = ', nn, ', Number of element = ', ne
!allocate ( coord(2,tol_nn), elem(4,tol_ne), material(tol_ne) )
!allocate( bound_node(tol_nn), skin_info(tol_nn), bound_normal(2,tol_nn), skin_normal(2,tol_nn) )

if ( mode == 0 ) then
	st_i = 0; ed_i = 0 
elseif ( mode == 1 ) then
	st_i = 1; ed_i = max_mat_num
else
	st_i = 0; ed_i = max_mat_num
endif

read (unit_out_num,*,iostat=status) a(1:8), num_var
allocate ( var(num_var) ) 
read (unit_out_num,*,iostat=status) a(1:3), var(:)
count = 1
all_var = '"'
do i = 1, num_var
    num = len_trim(var(i))
    write (all_var(count+1:count+num), '(A)') var(i)
    if ( i /= num_var ) then
        write (all_var(count+num+1:count+num+4), '(A)') '", "'
    else
        write (all_var(count+num+1:count+num+1), '(A)') '"'
    endif
    count = count + num + 4
enddo

!if ( mat_i == st_i ) then
fn = "./output/solid/tec/tec_eldata0000000.plt"
write (fn(30:36), '(I7.7)' ) current_i
if (flag == 1) then
    open (Unit=tec_unit_num, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)
	Write (tec_unit_num,'(A)') 'TITLE="Xfemout_quad_data"'
    if (Dim == 2) then
	    Write (tec_unit_num,'(A,A)') 'VARIABLES= "x", "y", "dx", "dy", ', all_var
    elseif (Dim == 3) then
        Write (tec_unit_num,'(A,A)') 'VARIABLES= "x", "y", "z", "dx", "dy", "dz", ', all_var
    endif
elseif (flag == 2) then
    ! Call tecini112
    if (Dim == 2) then
        VarsName = "x y dx dy stress11 stress22 stress33 stress12 strain11 strain22 strain33 strain12 prncpl_1 prncpl_2"
    elseif (Dim == 3) then
        VarsName = "x y z dx dy dz stress11 stress22 stress33 stress12 stress23 stress31 strain11 strain22 strain33 strain12 strain23 strain31 prncpl_1 prncpl_2 prncpl_3"
    endif
    Filetype = 0
    Debug = 0
    Visdouble = 0
    ier = tecini112("Result"//nulchar,TRIM(VarsName)//nulchar,TRIM(fn)//nulchar,".",FileType,Debug,Visdouble)
endif

num = 0 
cur_ne = 0
allocate (temp_quad_data(num_var,tol_ne*(Dim-1)*4))
call ext_quad_data(Dim, unit_out_num, num_var, tol_ne, temp_quad_data)
do i = 1, mat_mesh_set_num
    if ( output_region2 < 1 .OR. output_region2 == mat ) then
        mat_ne = mat_mesh_set(i)%ne
        ele_num = mat_mesh_set(i)%ele_num
        if (ele_num > 200 .and. ele_num < 400) then
            allocate (quad_data(num_var,mat_ne*(Dim-1)*4), stat=error)
            quad_data = temp_quad_data(:,cur_ne+1:cur_ne+mat_ne*(Dim-1)*4)
        
            call nodal_interpolation2(Dim, i, mat_ne, cur_time, num_var, quad_data, tec_unit_num, ier, flag)
        
            cur_ne = cur_ne + mat_ne*(Dim-1)*4
            deallocate (quad_data)
        endif
    endif
end do
rewind(unit_out_num) 
if (flag == 1) then
    close(tec_unit_num)
elseif (flag == 2) then
    ier = tecend112()
endif

deallocate (var, temp_quad_data) 
    
endif  !-------------------------------------------------------------------------------------------



end subroutine output_interpolation