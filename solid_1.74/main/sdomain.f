module system_domain

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (since October 28, 2006)
! Modification History:	December 05, 2006 (Version 1.0 beta)
!                       August 20, 2010 (Major modification: Symmetric SparseBlas-DSS)
!                       Feb, 2013 (Domain decomposition)
!                       April, 2013 (CSR-format Unsymmetric SparseBlas)
!                       March, 2014: element-level poor condition EXIT

use physical_domain
use time_domain

implicit none
save

type sdomain_type
  integer :: problem_size
  real :: max_stiff
  logical :: system_changed =.FALSE.
  real, allocatable :: u(:), p(:), p0(:), p_int(:), del_u(:), u_physical(:), u_cartesian(:)
  real, allocatable :: sparse_k(:), k(:,:), k_bb(:,:), p_b(:)
  real, allocatable :: bc_sparse_k(:)
  integer, allocatable :: g_ix(:)
  integer, allocatable :: row_index(:), col_addr(:), dof_index(:) ! CSR format
  integer, allocatable :: bc_row_index(:), bc_col_addr(:)    ! CSR format
  logical, allocatable :: row_existence(:), bc_row_existence(:)
  logical, allocatable :: p0_trans_flag(:)
end type	
  
integer, protected :: num_of_sdomains
integer, protected :: problem_type = 1  ! written in 'event_control'
logical, protected :: system_sparse =.FALSE., system_symmetricity =.FALSE., coarse_schur =.FALSE., coarse_jacobi =.FALSE.
logical, protected :: nonlinear_geometric_load = .FALSE.     ! written in 'femula_foro'
type (sdomain_type), allocatable :: sdomain(:) 


contains   
!==============================================================================


subroutine create_sdomain(number)

implicit none

integer, intent(in) :: number

allocate(sdomain(number))   
num_of_sdomains = number

end subroutine create_sdomain





subroutine create_full_system(rn)

! rn: region number

implicit none

integer, intent(in) :: rn
integer ::  m, m1, m2, size


m1 = region(rn)%num_dofs      ! # of active dof
m2 = region(rn)%num_constr    ! # of constrained dof
size = m1 + m2
sdomain(rn)%problem_size = size
sdomain(rn)%system_changed =.FALSE.

allocate(sdomain(rn)%dof_index(size))
allocate(sdomain(rn)%u(size))
allocate(sdomain(rn)%p(size))
allocate(sdomain(rn)%p0(size))
allocate(sdomain(rn)%p0_trans_flag(size))
allocate(sdomain(rn)%p_int(size))
allocate(sdomain(rn)%del_u(size))
allocate(sdomain(rn)%u_physical(size))
allocate(sdomain(rn)%u_cartesian(size))
allocate(sdomain(rn)%k(size,size))

sdomain(rn)%dof_index = region(rn)%dof_to_array
sdomain(rn)%u = 0.0
sdomain(rn)%del_u = 0.0
sdomain(rn)%p0_trans_flag = .FALSE.  ! not coordinate transformed (for MPC) yet

if (is_coarse_problem) then
  m = 2*region(rn)%num_skin_contact_nodes   ! 2 dof per node
  allocate(sdomain(rn)%k_bb(m,m))
  allocate(sdomain(rn)%p_b(m))
  allocate(sdomain(rn)%g_ix(m))
endif

end subroutine create_full_system





subroutine create_sparse_system(rn)

! rn: region number

implicit none

integer, intent(in) :: rn
integer :: m, m1, m2, size

m = region(rn)%num_nodes
m1 = region(rn)%num_dofs      ! # of active dof
m2 = region(rn)%num_constr    ! # of constrained dof
size = m1 + m2
sdomain(rn)%problem_size = size
sdomain(rn)%system_changed =.FALSE.

allocate(sdomain(rn)%dof_index(size))
allocate(sdomain(rn)%u(size))
allocate(sdomain(rn)%p(size))
allocate(sdomain(rn)%p0(size))
allocate(sdomain(rn)%p0_trans_flag(size))
allocate(sdomain(rn)%p_int(size))
allocate(sdomain(rn)%del_u(size))
allocate(sdomain(rn)%u_physical(size))
allocate(sdomain(rn)%u_cartesian(size))

sdomain(rn)%dof_index = region(rn)%dof_to_array
sdomain(rn)%u = 0.0
sdomain(rn)%del_u = 0.0
sdomain(rn)%p0_trans_flag = .FALSE.  ! not coordinate transformed (for MPC) yet

allocate(sdomain(rn)%bc_row_index(size+1)) 
allocate(sdomain(rn)%bc_row_existence(size))

call set_sparse(rn)

if (is_coarse_problem) then
  m = 2*region(rn)%num_skin_contact_nodes   ! 2 dof per node
  allocate(sdomain(rn)%k_bb(m,m))
  allocate(sdomain(rn)%p_b(m))
  allocate(sdomain(rn)%g_ix(m))
endif

end subroutine create_sparse_system







subroutine set_sparse(rn)

implicit none

integer, intent(in) :: rn 
integer :: i, j, k, nn, nnsp, ndof, nnd, nix, m1, size, node, num1, num2, st, ed, number, number0, bandwidth, ix_size
real :: time1, time2
integer,allocatable :: ix(:), inode(:,:), numberboard(:,:), boardsize(:), boardsize0(:), compact_row(:), board_row(:)


call cpu_time(time1)
write(*,*) '### set sparse starts --->'

nn = region(rn)%num_elements
nnsp = region(rn)%num_springs
ndof = region(rn)%num_node_dof
nnd = region(rn)%nodes_element
nix = ndof*nnd

m1 = region(rn)%num_dofs
size = sdomain(rn)%problem_size
bandwidth = MAXVAL(region(rn)%node(2,:))
write(*,*) '   dof bandwith =', bandwidth

if (allocated(sdomain(rn)%row_existence)) then
  deallocate(sdomain(rn)%row_index)
  deallocate(sdomain(rn)%row_existence)
endif
allocate(sdomain(rn)%row_index(m1+1))
allocate(sdomain(rn)%row_existence(m1))
allocate(ix(nix))
allocate(inode(ndof,nnd))

allocate(numberboard(size,bandwidth))
allocate(boardsize(size))
allocate(boardsize0(size))

numberboard = 0

boardsize = 0
boardsize0 = 0


do i = 1, nn  ! over all 2D/3D solid elements
  if (region(rn)%element_participation(i)) then
    do j = 1, nnd
      node = region(rn)%element(3+j,i)    ! connectivity nodes have been converted to internal numbers
      inode(1:ndof,j) = sdomain(rn)%dof_index(region(rn)%node(4:(ndof+3),node))
    end do
    ix = pack(inode, .TRUE.)
    ix_size = nnd*ndof
    forall (k = 1:ix_size)
      numberboard(ix(k),boardsize(ix(k))+1:boardsize(ix(k))+ix_size) = ix
      boardsize(ix(k)) = boardsize(ix(k)) + ix_size
    end forall
  endif
end do

! for spring elements -------
do i = 1, nnsp
  if (region(rn)%spring_participation(i)) then
    do j = 1, 2
      node = region(rn)%spring(3+j,i)    ! connectivity nodes have been converted to internal numbers
      inode(1:ndof,j) = sdomain(rn)%dof_index(region(rn)%node(4:(ndof+3),node))
    end do
    ix = pack(inode, .TRUE.)
    ix_size = 2*ndof
    forall (k = 1:ix_size)
      numberboard(ix(k),boardsize(ix(k))+1:boardsize(ix(k))+ix_size) = ix
      boardsize(ix(k)) = boardsize(ix(k)) + ix_size
    end forall
  endif
end do
!----------------------------

num1 = 1
num2 = 1
sdomain(rn)%row_index(1) = 1
sdomain(rn)%bc_row_index(1) = 1
sdomain(rn)%row_existence = .TRUE.
sdomain(rn)%bc_row_existence = .TRUE.

allocate(compact_row(bandwidth))
allocate(board_row(bandwidth))
compact_row = 0
board_row = 0

do i = 1, m1
  board_row = numberboard(i,:)

  call quicksort(1, board_row, m1, boardsize(i), compact_row, number0, number)

  num1 = num1 + number
  sdomain(rn)%row_index(i+1) = num1
  if (number == 0) sdomain(rn)%row_existence(i) = .FALSE.

  boardsize(i) = number
  boardsize0(i) = number0

  number = number0 - number
  num2 = num2 + number
  sdomain(rn)%bc_row_index(i+1) = num2
  if (number == 0) sdomain(rn)%bc_row_existence(i) = .FALSE.

  numberboard(i,:) = compact_row
end do

do i = m1+1, size   ! for K21 & K22 partitions
  board_row = numberboard(i,:)

  call quicksort(1, board_row, size, boardsize(i), compact_row, number0, number)

  num2 = num2 + number  ! number0 = number in rows belonging to [m1+1, size]
  sdomain(rn)%bc_row_index(i+1) = num2
  if (number == 0) sdomain(rn)%bc_row_existence(i) = .FALSE.

  boardsize(i) = 0
  boardsize0(i) = number0
  numberboard(i,:) = compact_row

end do

if (allocated(sdomain(rn)%col_addr)) then
  deallocate(sdomain(rn)%col_addr)
  deallocate(sdomain(rn)%sparse_k)
  deallocate(sdomain(rn)%bc_col_addr)
  deallocate(sdomain(rn)%bc_sparse_k)
endif  
  
num1 = num1 - 1
num2 = num2 - 1
allocate(sdomain(rn)%col_addr(num1))
allocate(sdomain(rn)%bc_col_addr(num2))

do i = 1, m1
  j = boardsize(i)
  if (sdomain(rn)%row_existence(i)) then
    st = sdomain(rn)%row_index(i)
    ed = sdomain(rn)%row_index(i+1) - 1
    sdomain(rn)%col_addr(st:ed) = numberboard(i,1:j)
  endif
  
  if (sdomain(rn)%bc_row_existence(i)) then
  	st = sdomain(rn)%bc_row_index(i)
    ed = sdomain(rn)%bc_row_index(i+1) - 1
    sdomain(rn)%bc_col_addr(st:ed) = numberboard(i,j+1:boardsize0(i))
  endif
end do

do i = m1+1, size
  j = boardsize0(i)

  if (sdomain(rn)%bc_row_existence(i)) then
  	st = sdomain(rn)%bc_row_index(i)
    ed = sdomain(rn)%bc_row_index(i+1) - 1
    sdomain(rn)%bc_col_addr(st:ed) = numberboard(i,1:j)
  endif
end do

deallocate(ix, inode, numberboard, boardsize, compact_row, board_row)

allocate(sdomain(rn)%sparse_k(num1))
allocate(sdomain(rn)%bc_sparse_k(num2))

call cpu_time(time2)
write(*,*) '### set sparse time: ', time2-time1
write(*,*) '### set sparse ends!'

end subroutine set_sparse







subroutine update_system_changed(action, rn, toggle)

! action = 0: reflect region(rn)%dof_changed (toggle becomes dummy)
!        = 1: write toggle into sdomain(rn)%system_changed

implicit none

integer, intent(in) :: action, rn
logical, intent(in) :: toggle

select case(action)
case(0)
  if (region(rn)%dof_changed) sdomain(rn)%system_changed = .TRUE.
case(1)
  sdomain(rn)%system_changed = toggle
  if (.NOT. toggle) call write_dof_changed(rn, .FALSE.)
case default
  STOP 'update_system_changed: wrong action number!'
end select

end subroutine update_system_changed






subroutine write_problem_type(type_number)

! write the protected variable, problem_type
! mainly used by 'eventcontrol'

implicit none

integer, intent(in) :: type_number

problem_type = type_number

end subroutine write_problem_type






subroutine write_system_type(sparse, symm)

! write the protected variables, system_sparse & system_symmetricity
! mainly used by 'eventcontrol'

implicit none

logical, intent(in) :: sparse, symm

system_sparse = sparse
if (.NOT. system_sparse) system_symmetricity = symm ! sparse option has priority

end subroutine write_system_type





subroutine write_domain_decomposition_type(flag)

! write the protected variables, system_sparse & system_symmetricity
! mainly used by 'eventcontrol'
! flag = 1: Schur
!      = 2: Jacobi

implicit none

integer, intent(in) :: flag


coarse_jacobi = .FALSE.

if (flag == 0) then
  coarse_schur = .FALSE.
elseif (flag == 1) then
  coarse_schur = .TRUE.
elseif (flag == 2) then
  coarse_schur = .TRUE.
  coarse_jacobi = .TRUE.
endif 

end subroutine write_domain_decomposition_type





subroutine write_load_geometric_nonlinearity

! write the protected variables, nonlinear_geometric_load
! called by 'femula_foro'

implicit none

integer :: i

nonlinear_geometric_load = .FALSE.
do i = 1, num_load_sets   ! num_load_sets: p_domain module global variable
  if ((load_set(i)%load_type == 2).OR.(load_set(i)%load_type == 4)) nonlinear_geometric_load = .TRUE.
enddo

end subroutine write_load_geometric_nonlinearity





subroutine linear_solver(rn, error, debug_target)

! rn: region number

implicit none

integer, intent(in) :: rn, debug_target(2)
real, intent(out) :: error
integer :: m1, m2, ass_flag, iter

m1 = region(rn)%num_dofs     ! # of active dof
m2 = region(rn)%num_constr   ! # of constrained dof

if (m1 < 1) STOP 'linear_solver: active dof is zero or less!'

call prescribed_u(rn, m2, 0)

iter = 0  ! no iteration here

if (system_sparse) then
  call assemble_sparse_k(rn, 1, ass_flag, iter, debug_target)
  call sp_solve(rn, m1)
else
  call assemble_k(rn, 1, ass_flag, iter, debug_target)
  call general_solve(rn, m1)
endif

error = 0.0
sdomain(rn)%p0_trans_flag = .FALSE.

end subroutine linear_solver






subroutine newton_solver(rn, max_itr, tolerance, error, status, debug_target)

! rn: region number

implicit none

integer, intent(in) :: rn, max_itr, debug_target(2)
real, intent(in) :: tolerance
real, intent(out) :: error
integer, intent(inout) :: status    ! converge_status = {-1, 1, 2, 3 ...}

integer :: m1, m2, iter, size, k11_size, ass_flag
real :: tol, tol_factor, error_old=0.0, DNRM2


m1 = region(rn)%num_dofs     ! # of active dof
m2 = region(rn)%num_constr   ! # of inactive dof
size = sdomain(rn)%problem_size

if (m1 < 1) STOP 'newton_solver: active dof is zero or less!'

tol_factor = 1.0
iter = 0

call prescribed_u(rn, m2, 0)
call assemble_p(rn)

do  !!!!!!!!!!!!!!!!!!!!!!

  call update_system_changed(0, rn, .TRUE.) ! dummy value: .TRUE.

  if (system_sparse) then

    if (sdomain(rn)%system_changed) call set_sparse(rn)

    call assemble_sparse_k(rn, 1, ass_flag, iter, debug_target)

  else
    call assemble_k(rn, 1, ass_flag, iter, debug_target)
  endif

  if (ass_flag < 0) then
     write(*,*) '>>>>>>>>>>> newton_solver: Element-level poor condition!'
     status = -1
     EXIT   !---------------------------------------------------------------------------->>>>>
  endif

  call update_system_changed(1, rn, .FALSE.)

  if (iter == 0) then
    error =  DNRM2(size, sdomain(rn)%p, 1) 
    if (abs(error) < tolerance) then
        status = 1
        EXIT   !------------------------------------------------->>>>>
    endif
    tol = max(tolerance, error*tolerance)
  else
    if (iter == 1) then
      call energy_norm(sdomain(rn)%p, sdomain(rn)%u, size, error)
      tol = max(tolerance, error*tolerance)
    endif

    call energy_norm(sdomain(rn)%p, sdomain(rn)%del_u, size, error)
  endif

  write(*,'(A,I6,A,2ES16.8)') '>>> iter #:', iter, '  error & tol: ', error*tolerance*tol_factor/tol, tolerance*tol_factor

  if (error < tol) then

    if ((iter < MIN(5, max_itr - 1)) .AND. (1.0e2*error_old**2 > error)) then
      status = MAX(0, status) + 1
      write(*,*) '>>>>>>>>>>> newton_solver: Good convergency!', status
    else
      status = 1
    endif

    EXIT    !---------------------------------------------------------------------------->>>>>
  endif

  if ((iter > 2) .AND. (error_old/error < 0.1)) then
     write(*,*) '>>>>>>>>>>> newton_solver: Poor convergency!'
     status = -1
     EXIT   !---------------------------------------------------------------------------->>>>>  
  elseif (iter == max_itr) then
    tol_factor = 1.0/sqrt(abs(tolerance)) 
    tol = tol_factor*tol  ! down the original tolerance
  elseif (iter  > max_itr) then
     write(*,*) '>>>>>>>>>>> newton_solver: Exceed maximum iteration!'
     status = -1
     EXIT   !---------------------------------------------------------------------------->>>>>
  endif

  if (system_sparse) then  
    call sp_solve(rn, m1)
  else
    call general_solve(rn, m1)
  endif

  iter = iter + 1
  error_old = error

end do  !!!!!!!!!!!!!!!!!!!!!!

sdomain(rn)%p0_trans_flag = .FALSE.

end subroutine newton_solver






subroutine general_solve(rn, m1)

! solve_flag = 1, 3: Symmetric system
!            = 2, 4: General system
! m1:  # of active dof
! m2: # of constrained dof

implicit none

integer, intent(in) :: rn, m1
real :: k11(m1,m1), pp(m1), uu(m1), time1, time2


k11 = sdomain(rn)%k(1:m1,1:m1)
pp = sdomain(rn)%p(1:m1)
uu = 0.0

! call cpu_time(time1)
! write(*,*) 'k11 =', k11
! write(*,'(A/,(12E11.3))') 'p =', sdomain(rn)%p

if (system_symmetricity) then
  call cholesky(k11, pp, uu, m1, 1)
else    
  call general_lu(k11, pp, uu, m1, 1)
endif

! call cpu_time(time2)
! write(*,*) '########## full-matrix solving cpu time:', time2-time1

call update_u(rn, uu, m1)

end subroutine general_solve






subroutine sp_solve(rn, m1)

! rn: region number
! m1: size of active_dof

use mkl_dss   ! MKL header module included in 'xfem_math.f'

implicit none

integer, intent(in) :: rn, m1
integer :: k11_size, perm(1), error
real :: pp(m1), uu(m1), time1, time2
type (MKL_DSS_HANDLE) :: handle 


k11_size = sdomain(rn)%row_index(m1+1) - 1
perm(1) = 0


pp = sdomain(rn)%p(1:m1)
uu = 0.0

! call cpu_time(time1)

error = dss_create(handle, MKL_DSS_DEFAULTS)
error = dss_define_structure(handle, MKL_DSS_SYMMETRIC_STRUCTURE, sdomain(rn)%row_index, m1, m1, sdomain(rn)%col_addr, k11_size)
error = dss_reorder(handle, MKL_DSS_METIS_OPENMP_ORDER, perm)
error = dss_factor(handle, MKL_DSS_DEFAULTS, sdomain(rn)%sparse_k)
error = dss_solve(handle, MKL_DSS_DEFAULTS, pp, 1, uu)
error = dss_delete(handle, MKL_DSS_DEFAULTS) 

! call cpu_time(time2)
! write(*,*) '########## sparse solving cpu time', time2-time1

call update_u(rn, uu, m1)

end subroutine sp_solve







subroutine prescribed_u(rn, m2, iter)

! 1. Prescribed Displacement for dispalcement control 
! 2. Acceleration and Velocity initialization at each timestep

! working in single region (rn)
! upload u from 'time_domain' data considering restart procedure 
! disp_factor: disp load factor - global module var from 'time_domain'

! iter = 0: initial iteration

implicit none

integer, intent(in) :: rn, m2, iter
integer :: i, ix1, ix2, constr_index(m2), init_flag
real :: disp_prescribed(m2)


sdomain(rn)%dof_index = region(rn)%dof_to_array

sdomain(rn)%u(sdomain(rn)%dof_index) = history(rn)%data(:,1)  

ix1 = sdomain(rn)%problem_size - m2 + 1
ix2 = sdomain(rn)%problem_size
constr_index = region(rn)%disp_bc

sdomain(rn)%del_u = 0.0
do i = 1, num_disp_sets
	if (disp_set(i)%region_info == rn) then		
		disp_prescribed = disp_factor(i)*disp_set(i)%constr(constr_index)
    sdomain(rn)%del_u(ix1:ix2) = sdomain(rn)%del_u(ix1:ix2) + disp_prescribed	
	endif

enddo

do i = num_disp_sets+1, total_disp_sets ! no disp_factor multiplication for contacts
	if (disp_set(i)%region_info == rn) then		
		disp_prescribed = disp_set(i)%constr(constr_index) 
    sdomain(rn)%del_u(ix1:ix2) = sdomain(rn)%del_u(ix1:ix2) + disp_prescribed	
	endif

enddo

sdomain(rn)%del_u(ix1:ix2) = sdomain(rn)%del_u(ix1:ix2) - sdomain(rn)%u(ix1:ix2)
sdomain(rn)%u = sdomain(rn)%u + sdomain(rn)%del_u

sdomain(rn)%u_physical = sdomain(rn)%u(sdomain(rn)%dof_index)

if (iter == 0) then
  init_flag = 0
else
  init_flag = 1
endif

call time_integration(init_flag, rn, problem_type, sdomain(rn)%u_physical)

end subroutine prescribed_u






subroutine update_u(rn, uu, m1)

! rn: region number
! m1: size of active_dof
! u, del_u: module global variables

implicit none

integer, intent(in) :: rn, m1
real, intent(in) :: uu(:)


sdomain(rn)%del_u = 0.0
sdomain(rn)%del_u(1:m1) = uu
sdomain(rn)%u = sdomain(rn)%u + sdomain(rn)%del_u
sdomain(rn)%u_physical = sdomain(rn)%u(sdomain(rn)%dof_index)

call time_integration(1, rn, problem_type, sdomain(rn)%u_physical)



end subroutine update_u






subroutine assemble_k(rn, flag, ass_flag, iter, debug_target)

! rn: region number
! flag : not used currently

implicit none

integer, intent(in) :: rn, flag, iter, debug_target(2)
integer, intent(out) :: ass_flag

integer :: i, j, k, m, m1, node, ndof, nnd, nix, nix2, action, nquad, n_space, n_out
integer :: num_elements, num_springs
integer, allocatable :: ix(:), inode(:,:), nodes_elem(:)
integer :: mpc_nodes(27), n_nodes, st_addr, ed_addr, plane_num
real :: time1, time2, time3, theta
real, allocatable :: ke(:,:), pe(:), ue(:)
real, allocatable :: k_ig(:,:), matrix(:,:), k11(:,:), k_gi(:,:), trans(:,:), p0_trans(:)
logical, allocatable :: node_flag(:)
logical :: trans_flag = .FALSE.

! call cpu_time(time1)

sdomain(rn)%k = 0.0
sdomain(rn)%p_int = 0.0

if (problem_type >= 3) then
  action = 2
else
  action = 1
endif
ass_flag = 1

num_elements = region(rn)%num_elements
num_springs = region(rn)%num_springs

ndof = region(rn)%num_node_dof
n_space = region(rn)%num_dim
n_out = region(rn)%size_str_output

if (num_elements > 0) then
  nnd = region(rn)%nodes_element
  nix = ndof*nnd
  nquad = region(rn)%num_quad_pts

  allocate(ke(nix,nix))
  allocate(pe(nix))
  allocate(ue(nix))
  allocate(inode(ndof,nnd))
  allocate(ix(nix))
  allocate(nodes_elem(nnd))
  allocate(trans(nix,nix))
  allocate(node_flag(nnd))
  allocate(p0_trans(nix))
endif

do i = 1, num_elements
  if (region(rn)%element_participation(i)) then
    m = region(rn)%element(2,i)   ! m: element type number

    do j = 1, nnd
      node = region(rn)%element(3+j,i)
      nodes_elem(j) = node
      inode(1:ndof,j) = sdomain(rn)%dof_index(region(rn)%node(4:(ndof+3),node))
    end do
    ix = pack(inode, .TRUE.)

    trans_flag = .FALSE.
    if (region(rn)%mpc_existence) then
      if (region(rn)%mpc_element_group(i) > 0) then
        trans_flag = .TRUE.
        theta = region(rn)%mpc_theta(region(rn)%mpc_element_group(i))
        plane_num = region(rn)%mpc_plane(region(rn)%mpc_element_group(i))

        st_addr = region(rn)%mpc_element_addr(i)
        ed_addr = region(rn)%mpc_element_addr(i+1) - 1
        n_nodes = region(rn)%mpc_element_addr(i+1) - st_addr
        mpc_nodes(1:n_nodes) = region(rn)%mpc_element_nodes(st_addr:ed_addr)

        node_flag = .FALSE.

        do k = 1, nnd
          node_flag(k) = ANY(nodes_elem(k) == mpc_nodes(1:n_nodes))
        end do

        call mpc_trans(nnd, ndof, node_flag, plane_num, theta, nix, trans)
!        write(*,'(A/,(24F5.2))') 'tranformation matrix = ', ((trans(j,k), k=1,nix),j=1,nix)
      endif
    endif

    call element_interface(action, iter, rn, i, m, pe, ke, ndof, nix, nnd, nquad, n_space, n_out, .TRUE., trans, ue, trans_flag)

    if (trans_flag) then
      ke = MATMUL(ke, trans)
      trans = TRANSPOSE(trans)
      ke = MATMUL(trans, ke)
      pe = MATMUL(trans, pe)
      if (iter == 0) then
        p0_trans = MATMUL(trans, sdomain(rn)%p0(ix))
        do k = 1, nix
          if (.NOT. sdomain(rn)%p0_trans_flag(ix(k))) then
            sdomain(rn)%p0(ix(k)) = p0_trans(k)
            sdomain(rn)%p0_trans_flag(ix(k)) = .TRUE.
          endif
        end do
      endif
    endif

    if ((rn == debug_target(1)) .AND. (i == debug_target(2))) then
      call matrix_output_debug(rn, i, iter, nix, nix, ke)
    endif


    if (action < 0) then  ! Element-level poor condition exit
      ass_flag = -1
      RETURN   !---------------------------------------------------------------------------->>>>>
    endif
  
    sdomain(rn)%k(ix,ix) = sdomain(rn)%k(ix,ix) + ke
    sdomain(rn)%p_int(ix) = sdomain(rn)%p_int(ix) + pe
    sdomain(rn)%u_cartesian(ix) = ue

  endif
end do

if (num_elements > 0) deallocate(ke, pe, ue, inode, ix, nodes_elem, trans, node_flag, p0_trans)

! for spring elements ------------

if (num_springs > 0) then
  nnd = 2
  nix = ndof*nnd
  nquad = 1

  allocate(ke(nix,nix))
  allocate(pe(nix))
  allocate(ue(nix))
  allocate(inode(ndof,nnd))
  allocate(ix(nix))
  allocate(trans(nix,nix))
endif

do i = 1, num_springs
  if (region(rn)%spring_participation(i)) then  
    m = region(rn)%spring(2,i)   ! m: element type number
    do j = 1, 2
      node = region(rn)%spring(3+j,i)
      inode(1:ndof,j) = sdomain(rn)%dof_index(region(rn)%node(4:(ndof+3),node))
    end do
    ix = pack(inode, .TRUE.)

    call element_interface(action, iter, rn, i, m, pe, ke, ndof, nix, nnd, nquad, n_space, n_out, .FALSE., trans, ue, .FALSE.)

    if (action < 0) then  ! Element-level poor condition exit
      ass_flag = -1
      RETURN   !---------------------------------------------------------------------------->>>>>
    endif
    sdomain(rn)%k(ix,ix) = sdomain(rn)%k(ix,ix) + ke
    sdomain(rn)%p_int(ix) = sdomain(rn)%p_int(ix) + pe
    sdomain(rn)%u_cartesian(ix) = ue
  endif
end do

if (num_springs > 0) deallocate(ke, pe, ue, inode, ix, trans)
!--------------------------------


!write(*,*) 'k',rn, sdomain(rn)%k
sdomain(rn)%p = sdomain(rn)%p_int + sdomain(rn)%p0
sdomain(rn)%u_cartesian = sdomain(rn)%u_cartesian(sdomain(rn)%dof_index)
!write(*,'(A/,(12E11.3))') 'u_cartesian =', sdomain(rn)%u_cartesian

! call cpu_time(time2)

!----------------------------------------------------------------------------
if (is_coarse_problem .AND. coarse_schur) then

  m = region(rn)%num_skin_contact_nodes

  nix = ndof*m

  allocate(inode(ndof,m))
  allocate(ix(nix))
  if (allocated(sdomain(rn)%g_ix)) then
    deallocate(sdomain(rn)%g_ix, sdomain(rn)%K_bb, sdomain(rn)%p_b)
    allocate(sdomain(rn)%k_bb(nix,nix))
    allocate(sdomain(rn)%p_b(nix))
    allocate(sdomain(rn)%g_ix(nix))
  endif

	do j = 1, m
	  node = region(rn)%skin_contact_nodes(j)
	  inode(1:ndof,j) = sdomain(rn)%dof_index(region(rn)%node(4:(ndof+3),node))
	end do
	ix = pack(inode, .TRUE.)  

  sdomain(rn)%g_ix = ix
  sdomain(rn)%k_bb = sdomain(rn)%k(ix,ix)
  sdomain(rn)%p_b = sdomain(rn)%p(ix)
  
  if (.NOT. coarse_jacobi) then
    m1 = region(rn)%num_dofs
    nix2 = nix + 1
  
    allocate(k_ig(m1,nix2))
    allocate(k_gi(nix,m1))
    allocate(k11(m1,m1))
    allocate(matrix(m1,nix2))
  
    k11 = sdomain(rn)%k(1:m1,1:m1)
    k_ig(1:m1,1:nix) = sdomain(rn)%k(1:m1,ix)
    k_ig(1:m1,nix2) = sdomain(rn)%p(1:m1)
    k_gi = sdomain(rn)%k(ix,1:m1)
    call general_lu(k11, k_ig, matrix, m1, nix2)
  
    sdomain(rn)%k_bb = sdomain(rn)%k_bb - MATMUL(k_gi, matrix(:,1:nix))      
    sdomain(rn)%p_b = sdomain(rn)%p_b - MATMUL(k_gi, matrix(:,nix2))            
    deallocate(k_ig, k11, k_gi, matrix)
  endif
  
  deallocate(ix, inode)
  
endif
! call cpu_time(time3)

! write(*,*) '########## coarse assemble full-matrix time', time3-time1
end subroutine assemble_k






subroutine assemble_sparse_k(rn, flag, ass_flag, iter, debug_target) 

! rn: region number
! flag = 0: local assemble only

use mkl_dss   ! MKL header module included in 'xfem_math.f'

implicit none

integer, intent(in) :: rn, flag, iter, debug_target(2)
integer, intent(out) :: ass_flag

integer :: i, j, k, l, m, node, ndof, nnd, nix, nix2, action, nquad, n_space, n_out, m1, size, nshape(2)
integer :: num_elements, num_springs
integer, allocatable :: ix(:), inode(:,:), full_to_sparse(:), jx(:), nodes_elem(:)
integer :: mpc_nodes(27), n_nodes, st_addr, ed_addr, plane_num
integer :: k11_size, error, perm(1)
real :: time1, time2, time3, theta
real, allocatable :: ke(:,:), pe(:), ue(:), trans(:,:), p0_trans(:)
real, allocatable :: k_gi(:,:), k_ig(:,:), matrix(:,:), vector(:), vvv(:),vvv2(:)
logical, allocatable :: node_flag(:)
logical :: trans_flag = .FALSE.

type (MKL_DSS_HANDLE) :: handle 


! call cpu_time(time1)
        
m1 = region(rn)%num_dofs
size = sdomain(rn)%problem_size
allocate(full_to_sparse(size))

sdomain(rn)%sparse_k = 0.0 
sdomain(rn)%bc_sparse_k = 0.0
sdomain(rn)%p_int = 0.0

if (problem_type >= 3) then
  action = 2
else
  action = 1
endif
ass_flag = 1

num_elements = region(rn)%num_elements
num_springs = region(rn)%num_springs

ndof = region(rn)%num_node_dof
n_space = region(rn)%num_dim
n_out = region(rn)%size_str_output

if (num_elements > 0) then
  nnd = region(rn)%nodes_element
  nix = ndof*nnd
  nquad = region(rn)%num_quad_pts

  allocate(ke(nix,nix))
  allocate(pe(nix))
  allocate(ue(nix))
  allocate(inode(ndof,nnd))
  allocate(ix(nix))
  allocate(jx(nix))
  allocate(nodes_elem(nnd))
  allocate(trans(nix,nix))
  allocate(node_flag(nnd))
  allocate(p0_trans(nix))
endif

do i = 1, num_elements
  if (region(rn)%element_participation(i)) then
    m = region(rn)%element(2,i)   ! m: element type number
    do j = 1, nnd
      node = region(rn)%element(3+j,i)
      nodes_elem(j) = node
      inode(1:ndof,j) = sdomain(rn)%dof_index(region(rn)%node(4:(ndof+3),node))
    end do
    ix = pack(inode, .true.)

    trans_flag = .FALSE.
    if (region(rn)%mpc_existence) then
      if (region(rn)%mpc_element_group(i) > 0) then
        trans_flag = .TRUE.
        theta = region(rn)%mpc_theta(region(rn)%mpc_element_group(i))
        plane_num = region(rn)%mpc_plane(region(rn)%mpc_element_group(i))

        st_addr = region(rn)%mpc_element_addr(i)
        ed_addr = region(rn)%mpc_element_addr(i+1) - 1
        n_nodes = region(rn)%mpc_element_addr(i+1) - st_addr
        mpc_nodes(1:n_nodes) = region(rn)%mpc_element_nodes(st_addr:ed_addr)

        node_flag = .FALSE.

        do k = 1, nnd
          node_flag(k) = ANY(nodes_elem(k) == mpc_nodes(1:n_nodes))
        end do

        call mpc_trans(nnd, ndof, node_flag, plane_num, theta, nix, trans)
!        write(*,'(A/,(24F5.2))') 'tranformation matrix = ', ((trans(j,k), k=1,nix),j=1,nix)
      endif
    endif

    call element_interface(action, iter, rn, i, m, pe, ke, ndof, nix, nnd, nquad, n_space, n_out, .TRUE., trans, ue, trans_flag)

    if (trans_flag) then
      ke = MATMUL(ke, trans)
      trans = TRANSPOSE(trans)
      ke = MATMUL(trans, ke)
      pe = MATMUL(trans, pe)
      if (iter == 0) then
        p0_trans = MATMUL(trans, sdomain(rn)%p0(ix))
        do k = 1, nix
          if (.NOT. sdomain(rn)%p0_trans_flag(ix(k))) then
            sdomain(rn)%p0(ix(k)) = p0_trans(k)
            sdomain(rn)%p0_trans_flag(ix(k)) = .TRUE.
          endif
        end do
      endif
    endif

    if ((rn == debug_target(1)) .AND. (i == debug_target(2))) then
      call matrix_output_debug(rn, i, iter, nix, nix, ke)
    endif

    
    if (action < 0) then  ! Element-level poor condition exit
      ass_flag = -1
      write(*,*) ' assemble_sparse_k: return & exit due to element-level poor condition'
      write(*,'(A,2I9/,A/,3(8I8/))') '   element (internal & user): ', i, region(rn)%element(1,i), '   nodes (user): ', (region(rn)%node(1,nodes_elem(j)), j=1,nnd)
      RETURN   !---------------------------------------------------------------------------->>>>>
    endif
  
    do j = 1, nix   ! over row number
      if (ix(j) > m1) then
        call sparse_address(rn, ix(j), full_to_sparse, size, 0)
        jx = full_to_sparse(ix)
        sdomain(rn)%bc_sparse_k(jx) = sdomain(rn)%bc_sparse_k(jx) + ke(j,:)
      else
        call sparse_address(rn, ix(j), full_to_sparse, size, 1) 
        jx = full_to_sparse(ix)

        do l = 1, nix   ! over column number
          if (ix(l) > m1) then 
            sdomain(rn)%bc_sparse_k(jx(l)) = sdomain(rn)%bc_sparse_k(jx(l)) + ke(j,l)
          else
            sdomain(rn)%sparse_k(jx(l)) = sdomain(rn)%sparse_k(jx(l)) + ke(j,l)      
          endif
        end do
      endif
    end do
    sdomain(rn)%p_int(ix) = sdomain(rn)%p_int(ix) + pe
  endif
  sdomain(rn)%u_cartesian(ix) = ue
end do

if (num_elements > 0) deallocate(ke, pe, ue, inode, ix, jx, nodes_elem, trans, node_flag,  p0_trans)


! for spring elements ------------

if (num_springs > 0) then
  nnd = 2
  nix = ndof*nnd
  nquad = 1

  allocate(ke(nix,nix))
  allocate(pe(nix))
  allocate(ue(nix))
  allocate(inode(ndof,nnd))
  allocate(ix(nix))
  allocate(jx(nix))
  allocate(trans(nix,nix))
endif

do i = 1, num_springs
  if (region(rn)%spring_participation(i)) then 
    m = region(rn)%spring(2,i)   ! m: element type number
    do j = 1, 2
      node = region(rn)%spring(3+j,i)
      inode(1:ndof,j) = sdomain(rn)%dof_index(region(rn)%node(4:(ndof+3),node))
    end do
    ix = pack(inode, .TRUE.)

    call element_interface(action, iter, rn, i, m, pe, ke, ndof, nix, nnd, nquad, n_space, n_out, .FALSE.,trans, ue, trans_flag)
    if (action < 0) then  ! Element-level poor condition exit
      ass_flag = -1
      write(*,*) ' assemble_sparse_k: return & exit due to element-level poor condition'
      RETURN   !---------------------------------------------------------------------------->>>>>
    endif
  
    do j = 1, nix   ! over row number
      if (ix(j) > m1) then
        call sparse_address(rn, ix(j), full_to_sparse, size, 0) 
        jx = full_to_sparse(ix)
        sdomain(rn)%bc_sparse_k(jx) = sdomain(rn)%bc_sparse_k(jx) + ke(j,:)
      else
        call sparse_address(rn, ix(j), full_to_sparse, size, 1) 
        jx = full_to_sparse(ix)
        do l = 1, nix   ! over column number
          if (ix(l) > m1) then 
            sdomain(rn)%bc_sparse_k(jx(l)) = sdomain(rn)%bc_sparse_k(jx(l)) + ke(j,l)
          else
            sdomain(rn)%sparse_k(jx(l)) = sdomain(rn)%sparse_k(jx(l)) + ke(j,l)      
          endif
        end do
      endif
    end do
    sdomain(rn)%p_int(ix) = sdomain(rn)%p_int(ix) + pe
  endif
  sdomain(rn)%u_cartesian(ix) = ue
end do

if (num_springs > 0) deallocate(ke, pe, ue, inode, ix, jx, trans)
!---------------------------------

deallocate(full_to_sparse)

!write(*,*) 'bc_sparse_k',rn, sdomain(rn)%bc_sparse_k

sdomain(rn)%p = sdomain(rn)%p_int + sdomain(rn)%p0
sdomain(rn)%u_cartesian = sdomain(rn)%u_cartesian(sdomain(rn)%dof_index)
!write(*,'(A/,(12E11.3))') 'u_cartesian =', sdomain(rn)%u_cartesian


! call cpu_time(time2)

!----------------------------------------------------------------------------
if (is_coarse_problem .AND. coarse_schur) then

  m = region(rn)%num_skin_contact_nodes
  nix = ndof*m

  allocate(inode(ndof,m))
  allocate(ix(nix))  
  allocate(vector(size))
  if (allocated(sdomain(rn)%g_ix)) then
    deallocate(sdomain(rn)%g_ix, sdomain(rn)%k_bb, sdomain(rn)%p_b)
    allocate(sdomain(rn)%k_bb(nix,nix))
    allocate(sdomain(rn)%p_b(nix))
    allocate(sdomain(rn)%g_ix(nix))
  endif

	do j = 1, m
	  node = region(rn)%skin_contact_nodes(j)
	  inode(1:ndof,j) = sdomain(rn)%dof_index(region(rn)%node(4:(ndof+3),node))
	end do
	ix = pack(inode, .TRUE.)  

  sdomain(rn)%g_ix = ix
  sdomain(rn)%p_b = sdomain(rn)%p(ix)

  if (.NOT. coarse_jacobi) then ! Schur complement
  
    nix2 = nix + 1
    m1 = region(rn)%num_dofs
    k11_size = sdomain(rn)%row_index(m1+1) - 1

    allocate(k_ig(m1,nix2))  
    allocate(k_gi(nix,m1))

    do i = 1, nix
      j = ix(i)
      call bc_sparse_to_full(rn, j, vector, size)
      sdomain(rn)%k_bb(i,:) = vector(ix)
      k_gi(i,:) = vector(1:m1)
    end do  
  
    do i = 1, m1
      call bc_sparse_to_full(rn, i, vector, size)
      k_ig(i,1:nix) = vector(ix)
    end do  
    k_ig(1:m1,nix2) = sdomain(rn)%p(1:m1)      

    allocate(vvv(m1*nix2))
    vvv = pack(k_ig, .TRUE.)
    
    nshape(1) = m1
    nshape(2) = nix2

    deallocate(k_ig)
    allocate(vvv2(m1*nix2))
    allocate(matrix(m1,nix2))  ! this one must exist just after vvv2 for good performance

    error = dss_create(handle, MKL_DSS_DEFAULTS)
    error = dss_define_structure(handle, MKL_DSS_SYMMETRIC_STRUCTURE, sdomain(rn)%row_index, m1, m1, sdomain(rn)%col_addr, k11_size)
    error = dss_reorder(handle, MKL_DSS_METIS_OPENMP_ORDER, perm)
    error = dss_factor(handle, MKL_DSS_DEFAULTS, sdomain(rn)%sparse_k)
    error = dss_solve(handle, MKL_DSS_DEFAULTS, vvv, nix2, vvv2)
    matrix = reshape(vvv2, nshape)  ! this one must exist here for good performance
    error = dss_delete(handle, MKL_DSS_DEFAULTS) 
    
    deallocate(vvv)   ! this one must exist here for good performance
    deallocate(vvv2)  ! this one must exist here for good performance
    
    sdomain(rn)%k_bb = sdomain(rn)%k_bb - MATMUL(k_gi, matrix(:,1:nix))                
    sdomain(rn)%p_b = sdomain(rn)%p_b - MATMUL(k_gi, matrix(:,nix2))   
        
    deallocate(matrix, k_gi)     
  
  elseif (coarse_jacobi) then  ! Jacobi
  
    nix2 = nix + 1
    m1 = region(rn)%num_dofs

    allocate(k_ig(m1,nix2))  
    allocate(k_gi(nix,m1))
    allocate(matrix(m1,nix2)) 
    
    do i = 1, nix
      j = ix(i)
      call bc_sparse_to_full(rn, j, vector, size)
      sdomain(rn)%k_bb(i,:) = vector(ix)
      k_gi(i,:) = vector(1:m1)
    end do
    
    do i = 1, m1
      call bc_sparse_to_full(rn, i, vector, size)
      k_ig(i,1:nix) = vector(ix)
    end do  
    k_ig(1:m1,nix2) = sdomain(rn)%p(1:m1)      
    
    do i = 1, m1
        call sparse_to_full(rn, i, vector, size)
        matrix(i,:) = k_ig(i,:)/vector(i)
    end do
           
    sdomain(rn)%k_bb = sdomain(rn)%k_bb - MATMUL(k_gi, matrix(:,1:nix))                
    sdomain(rn)%p_b = sdomain(rn)%p_b - MATMUL(k_gi, matrix(:,nix2))  
     
    deallocate(k_ig)       
    deallocate(matrix, k_gi)            
    
  endif     
    
  deallocate(ix, inode, vector)
endif
! call cpu_time(time3)

! write(*,*) '########## coarse assemble sparse time', time3-time1

end subroutine assemble_sparse_k







subroutine mpc_trans(nnd, ndof, node_flag, plane_num, theta, nix, trans)

implicit none

logical, intent(in) :: node_flag(nnd)
integer,intent(in) :: nnd, ndof, plane_num, nix
real, intent(in) :: theta
real,intent(out) :: trans(nix,nix)

integer :: i, k
real :: c, s


trans = 0.0
do i = 1, nix
  trans(i,i) = 1.0
end do

c = COS(theta)
s = SIN(theta)

if (ndof < 3) then
  do i = 1, nnd
    if (node_flag(i)) then
      k = (i-1)*ndof + 1
      trans(k,k) = c
      trans(k+1,k+1) = c
      trans(k,k+1) = -s
      trans(k+1,k) = s
    endif
  end do
else
  select case(plane_num)
  case(1) ! x-y plane
    do i = 1, nnd
      if (node_flag(i)) then
        k = (i-1)*ndof + 1
        trans(k,k) = c
        trans(k+1,k+1) = c
        trans(k,k+1) = -s
        trans(k+1,k) = s
      endif
    end do
  case(2) ! y-z plane
    do i = 1, nnd
      if (node_flag(i)) then
        k = (i-1)*ndof + 2
        trans(k,k) = c
        trans(k+1,k+1) = c
        trans(k,k+1) = -s
        trans(k+1,k) = s
      endif
    end do
  case(3) ! z-x plane
    do i = 1, nnd
      if (node_flag(i)) then
        k = (i-1)*ndof + 1
        trans(k,k) = c
        trans(k+2,k+2) = c
        trans(k,k+2) = s
        trans(k+2,k) = -s
      endif
    end do
  case default
    STOP 'mpc_trans: invalid MPC plane number!'
  end select
endif


end subroutine mpc_trans







subroutine sparse_address(rn, row_num, sparse_addr, dim, flag)

! CSR format
! row_num: row number
! dim: problem size
! flag = 1: from sparse_K and bc_sparse_K
!      = 0: from bc_sparse_K only
!      = 2: from sparse_K only
! sparse_addr(dim): full-size addresse vector (negative if lower) in sparse K matrix

implicit none

integer,intent(in) :: rn, row_num, dim, flag
integer,intent(out) :: sparse_addr(dim)
integer :: i, j, k, st, ed

sparse_addr = 0

if (flag >= 1) then
  if (sdomain(rn)%row_existence(row_num)) then
    st = sdomain(rn)%row_index(row_num)
    ed = sdomain(rn)%row_index(row_num+1) - 1
    do i = st, ed
      sparse_addr(sdomain(rn)%col_addr(i)) = i
    end do
  endif
endif

if (flag <= 1) then
  if (sdomain(rn)%bc_row_existence(row_num)) then
    st = sdomain(rn)%bc_row_index(row_num)
    ed = sdomain(rn)%bc_row_index(row_num+1) - 1
    do i = st, ed
      sparse_addr(sdomain(rn)%bc_col_addr(i)) = i
    end do
  endif
endif

end subroutine sparse_address






subroutine bc_sparse_to_full(rn, row_num, vector, dim)

! New version: April, 2013
!
! CSR format
! row_num: row number
! dim: problem size

! vector(dim): full-size row vector

implicit none

integer, intent(in) :: rn, row_num, dim
real, intent(out) :: vector(dim)
integer :: st, ed

vector = 0.0

if (sdomain(rn)%bc_row_existence(row_num)) then
  st = sdomain(rn)%bc_row_index(row_num)
  ed = sdomain(rn)%bc_row_index(row_num+1) - 1

  vector(sdomain(rn)%bc_col_addr(st:ed)) = sdomain(rn)%bc_sparse_k(st:ed)
endif

end subroutine bc_sparse_to_full





subroutine sparse_to_full(rn, row_num, vector, dim)

! New version: April, 2013
!
! CSR format
! row_num: row number
! dim: problem size

! vector(dim): full-size row vector

implicit none

integer, intent(in) :: rn, row_num, dim
real, intent(out) :: vector(dim)
integer :: st, ed

vector = 0.0

if (sdomain(rn)%row_existence(row_num)) then
  st = sdomain(rn)%row_index(row_num)
  ed = sdomain(rn)%row_index(row_num+1) - 1

  vector(sdomain(rn)%col_addr(st:ed)) = sdomain(rn)%sparse_k(st:ed)
endif

end subroutine sparse_to_full





subroutine checkboard_sparse_to_full(rn, row_num, checkboard, dim)

! New version: Sept, 2015
!
! CSR format
! row_num: row number
! dim: problem size

! checkboard(dim): full-size logical row vector

implicit none

integer, intent(in) :: rn, row_num, dim
logical, intent(out) :: checkboard(dim)
integer :: st, ed


checkboard = .FALSE.

if (sdomain(rn)%row_existence(row_num)) then
  st = sdomain(rn)%row_index(row_num)
  ed = sdomain(rn)%row_index(row_num+1) - 1

  checkboard(sdomain(rn)%col_addr(st:ed)) = .TRUE.
endif

if (sdomain(rn)%bc_row_existence(row_num)) then
  st = sdomain(rn)%bc_row_index(row_num)
  ed = sdomain(rn)%bc_row_index(row_num+1) - 1

  checkboard(sdomain(rn)%bc_col_addr(st:ed)) = .TRUE.
endif

end subroutine checkboard_sparse_to_full





subroutine find_pop_up(rn, tol, flag)

! find extreme values in p vector
! rn: region number

implicit none

integer, intent(in) :: rn
real, intent(in) :: tol
logical, intent(out) :: flag
integer :: i, size

size = sdomain(rn)%problem_size
flag = .FALSE.

do i = 1, size
  if (sdomain(rn)%p_int(i) > tol) then
    write(*,*) '!!!!! find DOF number and p_int value:', i, sdomain(rn)%p_int(i)
    flag = .TRUE.
  endif
    if (sdomain(rn)%p0(i) > tol) then
    write(*,*) '!!!!! find DOF number and p0 value:', i, sdomain(rn)%p0(i)
    flag = .TRUE.
  endif
  if (sdomain(rn)%p(i) > tol) then
    write(*,*) '!!!!! find DOF number and p value:', i, sdomain(rn)%p(i)
    flag = .TRUE.
  endif
end do

end subroutine find_pop_up




subroutine pressure_load_update(rn)

! Last modification: Aug, 2015 by Jeeho Lee
! 2D & 3D version
! 2D: over all groups in a region
! rn: region number

implicit none

integer, intent(in) :: rn
integer :: i, j, ndim, ngroup, num_skin, n1, n2, ix(3), ix_monitor(3)
integer :: num_surfaces, face_number, elem_number, face_nodes(4)
integer, allocatable :: the_dof(:)
real :: min_len, corners(3,4)
integer, allocatable :: skin_nodes(:)
real, allocatable :: uuu(:,:), trans(:,:,:)


ndim = region(rn)%num_dim


if (ndim == 2) then
ngroup = region(rn)%num_groups
 
do i = 1, ngroup
  if (i < ngroup) then  
    num_skin = region(rn)%skin_group_pointer(i+1) - region(rn)%skin_group_pointer(i)
  else
    num_skin = region(rn)%num_skin + 1 - region(rn)%skin_group_pointer(i)  
  endif
  
  allocate(skin_nodes(num_skin))
  allocate(uuu(ndim, num_skin))
  
  n1 = region(rn)%skin_group_pointer(i) ! first node
  n2 = n1 + num_skin - 1            ! last node
  skin_nodes = region(rn)%skin_nodes(n1:n2)
  
  do j = 1, num_skin
    ix(1:ndim) = sdomain(rn)%dof_index(region(rn)%node(4:(ndim+3),skin_nodes(j)))
  	call extract_ue(rn, uuu(:,j), ix, ndim)  
  end do
  
  call update_skin_normal(1, rn, i, ndim, num_skin, uuu, min_len)

  do j = 1, num_load_sets  ! global variable in physical_domain
    if (load_set(j)%surface_pressure_existence) call surface_load(1, j, rn)
  end do  ! j

  deallocate(skin_nodes, uuu)

end do     
            
elseif (ndim == 3) then
  num_surfaces = region(rn)%num_skin_surfaces

  allocate(trans(3,3,num_surfaces))
  allocate(uuu(ndim,4))

  do i = 1, num_surfaces

    face_number = region(rn)%skin_surfaces(1,i)
    elem_number = region(rn)%skin_surfaces(2,i)
    face_nodes = region(rn)%element_face(1:4,face_number,elem_number)

    do j = 1, 4
      ix(1:3) = sdomain(rn)%dof_index(region(rn)%node(4:6,face_nodes(j)))
      call extract_ue(rn, uuu(:,j), ix, ndim)
    end do  ! j

    corners = region(rn)%node_coord(1:3,face_nodes)
    corners =   corners + uuu

    call surface_trans_matrix(corners, trans(:,:,i))

  end do  ! i

  call update_surface_normal(rn, num_surfaces, trans(3,:,:))

  do i = 1, num_load_sets  ! global variable in physical_domain
    if (load_set(i)%surface_pressure_existence) call surface_load(1, i, rn)
  end do  ! i

  deallocate(trans, uuu)

endif


end subroutine pressure_load_update






subroutine assemble_p(rn)

! rn: region number

implicit none

integer, intent(in) :: rn
integer :: i, j, k, m, node, radius_dof, ld_type_num
integer, allocatable :: the_dof(:)
real, allocatable :: p_value(:)
real :: radius, rfactor
logical :: flag


sdomain(rn)%p0 = 0.0

if (nonlinear_geometric_load) then    ! nonlinear geometric update
  call pressure_load_update(rn)
endif

do i = 1, num_load_sets
	if (load_set(i)%region_info == rn) then
    ld_type_num = load_set(i)%load_type
    
! load_type_number == 2: axisymmetric (nonlinear)
!                  == 3: axisymmetric loads (linear)    
! load_factor: load factor - global module var from 'time_domain'

		m = load_set(i)%num_load
		allocate(the_dof(m))
		allocate(p_value(m))    
		the_dof = sdomain(rn)%dof_index(load_set(i)%load_bc)
    

    if (ld_type_num == 2) then   ! nonlinear geometric axi-symmetric case
      rfactor = 2.0*SIN(0.5)
      do j = 1, m
        node = load_set(i)%load_node(j)
        radius = region(rn)%node_coord(1,node)
        radius_dof = region(rn)%node(4,node)  ! assume radius-direction is 1st dof direction
        radius = radius + sdomain(rn)%u(sdomain(rn)%dof_index(radius_dof))   ! current deformed radius coordinate
        radius = rfactor*radius
        p_value(j) = load_factor(i)*load_set(i)%load(j)*radius 
        
!        if (node == 2000) write(*,*) '********** the skin node=', radius, load_factor(i), p_value(j)
          
      end do

    elseif (ld_type_num == 3) then   ! axi-symmetric case
      rfactor = 2.0*SIN(0.5)
      do j = 1, m
        node = load_set(i)%load_node(j)
        radius = region(rn)%node_coord(1,node)
        radius = rfactor*radius
        p_value(j) = load_factor(i)*load_set(i)%load(j)*radius 
      end do

    else

      p_value = load_factor(i)*load_set(i)%load
      
    endif

    do j = 1, m   ! make this sequential 'do' to prevent errors from the same 'the_dof' numbers
      sdomain(rn)%p0(the_dof(j)) = sdomain(rn)%p0(the_dof(j)) + p_value(j)
    end do

		deallocate(the_dof, p_value)
	endif
  

end do

do i = num_load_sets+1, total_load_sets
  if (load_set(i)%region_info == rn) then
    m = load_set(i)%num_load
    allocate(the_dof(m))
    allocate(p_value(m))
    the_dof = sdomain(rn)%dof_index(load_set(i)%load_bc)
!    ld_type_num = load_set(i)%load_type

    p_value = load_set(i)%load

    do j = 1, m   ! make this sequential 'do' to prevent errors from the same 'the_dof' numbers
      sdomain(rn)%p0(the_dof(j)) = sdomain(rn)%p0(the_dof(j)) + p_value(j)
    end do

    deallocate(the_dof, p_value)
  endif
end do

! write(*,*) 'load vector:', sdomain(rn)%p0

end subroutine assemble_p







subroutine revise_local_r(id, rn1, rn2, del_u_g, del_u_g2, n1, n2)

! rn: region number
! n1, n2: size of active_dof
! u, del_u: module global variables

implicit none

integer, intent(in) :: id, rn1, rn2, n1, n2
real, intent(in) :: del_u_g(n1), del_u_g2(n2)
integer :: m1, m2, ix1(n1), ix2(n2), size, sparse_size
real, allocatable :: vector(:)

m1 = region(rn1)%num_dofs
ix1 = sdomain(rn1)%g_ix
m2 = region(rn2)%num_dofs
ix2 = sdomain(rn2)%g_ix

if (system_sparse) then
  allocate(vector(m1))
  size = sdomain(rn1)%problem_size
  sparse_size = sdomain(rn1)%bc_row_index(size+1) - 1
  call mul_sparse_matrix_mask(rn1, size, m1, 1, sparse_size, n1, ix1, del_u_g, vector)
  sdomain(rn1)%p(1:m1) = sdomain(rn1)%p(1:m1) - vector
  deallocate(vector)

  allocate(vector(m2))
  size = sdomain(rn2)%problem_size
  sparse_size = sdomain(rn2)%bc_row_index(size+1) - 1
  call mul_sparse_matrix_mask(rn2, size, m2, 1, sparse_size, n2, ix2, del_u_g2, vector)
  sdomain(rn2)%p(1:m2) = sdomain(rn2)%p(1:m2) - vector
  deallocate(vector)

else
  sdomain(rn1)%p(1:m1) = sdomain(rn1)%p(1:m1) - MATMUL(sdomain(rn1)%k(1:m1,ix1), del_u_g)
  sdomain(rn2)%p(1:m2) = sdomain(rn2)%p(1:m2) - MATMUL(sdomain(rn2)%k(1:m2,ix2), del_u_g2)
endif

end subroutine revise_local_r






subroutine mul_sparse_matrix_mask(rn, m, l, n, ndim, n_mask, mask, matrix, result)

! R = S[:,mask]A
! sparse(S): sparse matrix (m rows) in CSR format
! matrix(A): a matrix (n_mask-by-n)
! mask(n_mask): mask index array for columns of S
! result matrix(R): l-by-n
! ndim: size(length) of row_info & sparse arrays

implicit none

integer,intent(in) :: rn, m, l, n, ndim, n_mask, mask(n_mask)
real,intent(in) ::  matrix(n_mask,n)
real,intent(out) :: result(l,n)

integer :: i, j
real, allocatable :: vector(:)

result  = 0.0
allocate(vector(m))

do j = 1, n
  do i = 1, l
    call bc_sparse_to_full(rn, i, vector, m)
    result(i,j) = DOT_PRODUCT(vector(mask),matrix(:,j))
  end do
end do 

deallocate(vector)

end subroutine mul_sparse_matrix_mask





subroutine vector_norm(vect,ndim,norm)

implicit none
integer, intent(in) :: ndim
real, intent(in) :: vect(ndim)
real, intent(out) :: norm
real :: DNRM2


norm =  DNRM2(ndim, vect, 1)

end subroutine vector_norm





subroutine energy_norm(vect,du,ndim,norm)

implicit none
integer, intent(in) :: ndim
real, intent(in) :: vect(ndim), du(ndim)
real, intent(out) :: norm
integer :: i
real :: sum

sum = 0.0
do i = 1, ndim
  sum = sum + abs(vect(i)*du(i))
!  if (sum > 1.0e20) then
!     write(*,*) 'energy_norm: too big!: i, vect & du: ', i, vect(i), du(i)
!     stop
!  endif
end do

norm = sum

end subroutine energy_norm





subroutine extract_ue(rn, ue, ix, ndim)

implicit none
integer, intent(in) :: rn, ix(ndim), ndim
real, intent(out) :: ue(ndim)

ue = sdomain(rn)%u(ix)

end subroutine extract_ue





subroutine read_u(rn, nodes, u_output, n, ndim)

! nodes: user node numbers
! n: number of dofs at a node
! ndim: number of nodes to read

implicit none

integer, intent(in) :: rn, nodes(ndim), n, ndim
real, intent(out) :: u_output(n,ndim)
integer :: dofs(n,ndim)
integer :: i

call read_dof_number(rn, nodes, dofs, n, ndim)  ! find dofs using user node numbers

do i = 1, ndim
   u_output(1:n,i) = sdomain(rn)%u_physical(dofs(1:n,i))
!	u_output(1:n,i) = sdomain(rn)%u(sdomain(rn)%dof_index(dofs(1:n,i)))
end do

end subroutine read_u





subroutine read_r(rn, nodes, p_output, n, ndim)

! read the internal forces from p_int
! n: number of dofs at a node
! ndim: number of nodes to read

implicit none

integer, intent(in) :: rn, nodes(ndim), n, ndim
real, intent(out) :: p_output(n,ndim)
integer :: dofs(n,ndim)
integer :: i


call read_dof_number(rn, nodes, dofs, n, ndim)  ! find dofs using user node numbers

do i = 1, ndim
	p_output(1:n,i) = sdomain(rn)%p_int(region(rn)%dof_to_array(dofs(1:n,i)))
end do

end subroutine read_r





subroutine read_reaction(rn, nodes, p_output, n, ndim, count, rct_nodes)

! read the internal forces from p_int
! n: number of dofs at a node
! ndim: number of nodes to read

implicit none

integer, intent(in) :: rn, nodes(ndim), n, ndim
real, intent(out) :: p_output(n,ndim)
integer, intent(out) :: count, rct_nodes(ndim)
integer :: dofs(n,ndim)
integer :: i, j


call read_dof_number(rn, nodes, dofs, n, ndim)  ! find dofs using user node numbers

count = 0
do i = 1, ndim
  if (.NOT. ALL(region(rn)%dof_activity_flag(dofs(1:n,i)))) then
    count = count + 1
    rct_nodes(count) = i
    p_output(1:n,count) = sdomain(rn)%p(region(rn)%dof_to_array(dofs(1:n,i)))
    do j = 1, n
      if (region(rn)%dof_activity_flag(dofs(j,i))) p_output(j,count) = 0.0
    end do
  endif
end do

end subroutine read_reaction






logical function is_my_node(rn, node)

integer,intent(in) :: rn, node
integer,allocatable ::  node_info(:)
integer :: num, i

num = region(rn)%num_nodes
allocate(node_info(num))
node_info = region(rn)%node(1,:)

is_my_node = .FALSE.

do i = 1, num
	if (node_info(i) == node) then
		is_my_node = .TRUE.
		exit
	endif
end do

end function is_my_node






subroutine node_output(unit_num, key, rn, output_rn, node)

! key: length = 6   ('output_cnt' in femula_prima)
! node: output nodal point in user node number
! called by: feap_prima

implicit none

integer, intent(in) :: unit_num, rn, output_rn, node
character(len=*), intent(in) :: key
integer :: m, n, i, j, n_reaction
integer, allocatable :: nodes(:), rct_nodes(:)
real, allocatable :: results(:,:), results2(:,:), sum(:)

select case(key)
case ('displa') ! displacements
  if (rn == output_rn) then
    m = region(rn)%num_node_dof
    n = 1
    allocate(nodes(n))
    allocate(results(m,n))
	
    nodes(1) = node

    call read_u(rn, nodes, results, m, n)

    do i = 1, n
!      write(*,100)  ' * Displacement result at region # =', rn, ' /node # =', node, (results(j,i), j=1,m)
      write(unit_num,101)  time, (results(j,i), j=1,m)
    end do
!    write(*,*) ' ----------------------------------------------------------'
    deallocate(results)
  endif
  
case ('reacti') ! reactions
  if (rn == output_rn) then
    m = region(rn)%num_node_dof
    n = 1
    allocate(nodes(n))
    allocate(results(m,n))
	
    nodes(1) = node

    call read_r(rn, nodes, results, m, n)
    results = -results
    do i = 1, n
!      write(*,100)  ' * Reaction result at region # =', rn, ' /node # =', node, (results(j,i), j=1,m)
      write(unit_num,101)  time, (results(j,i), j=1,m)
    end do
!    write(*,*) ' ----------------------------------------------------------'
    deallocate(results)
  endif

case ('disrct','disrea','dsprct','dsprea') ! displacements & reactions
  if (rn == output_rn) then
    m = region(rn)%num_node_dof
    n = 1
    allocate(nodes(n))
    allocate(results(m,n))
    allocate(results2(m,n))
	
    nodes(1) = node

    call read_u(rn, nodes, results, m, n)
    call read_r(rn, nodes, results2, m, n)
    results2 = -results2
    do i = 1, n
      write(unit_num,101)  time, (results(j,i), j=1,m), (results2(j,i), j=1,m)
    end do
!    write(*,*) ' ----------------------------------------------------------'
    deallocate(results, results2)
  endif
  
case ('rctall','reaall') ! all reactions
        m = region(rn)%num_node_dof
        n = region(rn)%num_nodes
        allocate(nodes(n), rct_nodes(n))
        allocate(results(m,n), sum(m))

        nodes = region(rn)%node(1,:)

        call read_reaction(rn, nodes, results, m, n, n_reaction, rct_nodes)
        results = -results
!        write(*,'(A,I3)')  ' * Reaction result at region # = ', rn
				write(unit_num,'(A,ES15.7)') 'time = ', time
        write(unit_num,'(A,I5)') 'region = ', rn
        write(unit_num,'(A)')   'node #  X coordinate Y coordinate Z coordinate'
        sum = 0.0
        do i = 1, n_reaction
          write(unit_num,102)   rct_nodes(i), (results(j,i), j=1,m)
          sum = sum + results(:,i)
        end do
        write(unit_num,'(A,3ES13.5)') ' sum =', (sum(j), j=1,m)
        write(unit_num,'(/)')
!        write(*,*) '>>>> All reaction values are written for this time step!'   
!        write(*,*) ' ----------------------------------------------------------'
				
!				rewind(unit=unit_num)
				
        deallocate(nodes, rct_nodes, results, sum)

  
case ('dspall','disall') ! all displacements
        m = region(rn)%num_node_dof
        n = region(rn)%num_nodes
        allocate(nodes(n))
        allocate(results(m,n))

        nodes = region(rn)%node(1,:)

        call read_u(rn, nodes, results, m, n)
	    
!        write(*,'(A,I3)')  ' * Displacement result at region # = ', rn
				write(unit_num,'(A,ES15.7)') 'time = ', time
        write(unit_num,'(A,I5)') 'region = ', rn
        write(unit_num,'(A)')   'node #  X coordinate Y coordinate Z coordinate'

        do i = 1, n
          write(unit_num,102)   nodes(i), (results(j,i), j=1,m)
        end do
        write(unit_num,'(/)')
!        write(*,*) '>>>> All displacement values are written for this time step!'   
!        write(*,*) ' ----------------------------------------------------------'
				
!				rewind(unit=unit_num)
				
        deallocate(nodes, results)
        
case default
        write(*,*) '>>>>>>>>> Nothing is writeen ! key= ', key
end select

100 format (/,A,I3,A,I6,/,6ES13.5)
101 format (7ES13.5)
102 format (I6,6ES13.5)

end subroutine node_output


!==============================================================================

subroutine update_motion_3D(rn, dim,  position , val)

implicit none

integer,intent(in) :: rn, dim
integer,intent(in) :: position(dim)
real, intent(in) :: val(dim,3)

sdomain(rn)%u_cartesian(position) = val(:,1)

end subroutine update_motion_3D

subroutine free_sdomain

implicit none

deallocate ( sdomain )
problem_type = 1                ! written in 'event_control'
system_sparse =.FALSE.; system_symmetricity =.FALSE.
coarse_schur =.FALSE.; coarse_jacobi =.FALSE.

end subroutine

!==============================================================================

end module system_domain
