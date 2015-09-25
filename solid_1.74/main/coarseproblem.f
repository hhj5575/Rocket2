module coarse_domain

! Copyright(2012~) Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (February 20, 2012)

use file_info_house
use physical_domain
use time_domain
use system_domain

implicit none
save

type interface_type
  integer :: contact_type, contact_region(2), num_contact_nodes(2)  ! 1: ball region, 2: target region
  integer :: num_dofs, dof_offset, contact_dim, num_contact_group
  integer, allocatable :: contact_skin_ball(:), contact_skin_target(:) ! contact node numbers
  integer, allocatable :: contact_target_partner(:,:) ! contact ball face info for each target node: 3D only
  integer, allocatable :: contact_group_end_addr1(:), contact_group_end_addr2(:)
  logical(1), allocatable :: contact_ball_island(:) ! island (isolated) ball surface = .TRUE.

  integer, allocatable :: col_index(:), row_addr(:) ! CSC format for Trans matrix
  real, allocatable :: sparse_Trans(:), u_g(:)
  logical(1), allocatable :: relation(:,:)

  integer, allocatable :: global_row_index(:), global_col_addr(:) ! CSR format
  real, allocatable :: global_K(:), global_p(:)
  integer, allocatable :: ix_to_g_pointer(:)
  logical(1), allocatable :: checkboard_gg(:,:), checkboard_gi(:,:), checkboard_ig(:,:)

  integer, allocatable :: partner(:,:), old_number(:)
  real, allocatable :: partner_d(:,:)
  logical(1), allocatable :: attach_info(:), nonoverlap_dof(:)

end type interface_type

type(interface_type), protected, allocatable :: contact(:)
integer, protected :: number_of_contacts = 0
real, parameter :: d_tolernace = 1.0e-2

! parameter constants from 'physical_domain': PI, Numerical_Zero

contains
!==============================================================================


subroutine create_contact(n)

implicit none

integer, intent(in) :: n
integer :: state_flag, stat
character(len=5) :: keyword
character(len=3) :: toggle

number_of_contacts = n
write(*,'(/A,I5/)') '>>> Number of contacts:', n
if (n > 0) then
  allocate(contact(n))

  is_coarse_problem = .TRUE.  ! global variable in physical_domain
  rewind(unit_number)
  call find_starword(unit_number, 'schu', state_flag)
  if (state_flag < 0) toggle = 'off'    ! default
  
  read(unit_number, *, iostat=stat) keyword, toggle
  
  if (toggle(1:2) == 'on') then
    call write_domain_decomposition_type(1)
    write(*,*) '  > Schur complement: ON'
  elseif (toggle == 'off') then
    call write_domain_decomposition_type(0)
    write(*,*) '  > Schur complement: OFF'
  elseif (toggle == 'jac') then
    call write_domain_decomposition_type(2)
    write(*,*) '  > Jacobi preconditioner: ON'
  endif
  write(*,*)
  rewind(unit_number)
endif

end subroutine create_contact






subroutine reset_contact(id, contact_type)

implicit none

integer, intent(in) :: id, contact_type


deallocate(contact(id)%contact_skin_ball, contact(id)%contact_skin_target)
deallocate(contact(id)%contact_group_end_addr1, contact(id)%contact_group_end_addr2)
if (contact_type == 2) then
  deallocate(contact(id)%relation, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, contact(id)%u_g)
elseif (contact_type == 3) then
  deallocate(contact(id)%partner, contact(id)%partner_d, contact(id)%attach_info)
endif

end subroutine reset_contact






subroutine write_contact(action, id, reg_from, reg_to, contact_type, contact_pt_from, contact_pt_to, contact_dim, num_points, num_contact_group)

! action = 1: create & build all types of contacts
!        = 2: rebuild only for contact_type = 2 (when nodal numbers are NOT changed)

implicit none

integer, intent(in) :: action, id, reg_from, reg_to, contact_type, contact_dim, num_points, num_contact_group
real, intent(in) :: contact_pt_from(contact_dim,num_points,num_contact_group), contact_pt_to(contact_dim,num_points,num_contact_group)

integer, allocatable :: snode(:), bound_node_1(:), bound_node_2(:), skin_node(:), node_dofs(:,:), moving_dofs(:), old_bound_node(:), free_node(:)
integer :: i, j, count, n_skin, number, num1, num2, set_number, n_dofs, start_point, end_point, shape(2)
integer :: group1, group2, size_of_col, size_of_row, old_number
integer :: num_dim
real :: tol, corners(3,4)
logical :: bingo


if (contact_type /= 2) STOP 'write_contact: currently only contact_type 2 is valid!'

contact(id)%contact_type = contact_type
contact(id)%contact_region(1) = reg_from
contact(id)%contact_region(2) = reg_to
contact(id)%contact_dim = contact_dim
contact(id)%num_contact_group = num_contact_group

if (contact_dim /= region(reg_from)%num_dim) then
  write(*,*) 'write_contact: contact dimension is not same wuth region dim!', contact_dim, region(reg_from)%num_dim
  STOP
endif

if (action == 2) then
  old_number = contact(id)%num_contact_nodes(2)
  allocate(old_bound_node(old_number))
  old_bound_node = contact(id)%contact_skin_target(:)
  call reset_contact(id, contact_type)    
endif  ! action == 2


if (contact_dim == 2) then  ! 2D -------------------------------------------------------------------

! ball (region_from) surface -----------------------
  n_skin = region(reg_from)%num_skin
  allocate(bound_node_1(n_skin))
  allocate(skin_node(n_skin))
  allocate(contact(id)%contact_group_end_addr1(num_contact_group))
  
  tol = 0.1*region(reg_from)%skin_length ! tolerance based on minimum distance between two skin nodes

  num1 = 0  
  do i = 1, num_contact_group  
    call find_surface_node(reg_from, n_skin, contact_pt_from(:,:,i), tol, group1, start_point, end_point)
 
    call find_skin(reg_from, group1, start_point, end_point, number, bound_node_1, skin_node, num1)
    num1 = num1 + number
    contact(id)%contact_group_end_addr1(i) = num1
  end do
  
  allocate(contact(id)%contact_skin_ball(num1))

  contact(id)%contact_skin_ball(:) = bound_node_1(1:num1)
  contact(id)%num_contact_nodes(1) = num1
  
  if (contact_type == 2) call write_skin_contact(reg_from, bound_node_1(1:num1), num1) ! physical_domain

  deallocate(skin_node)


! target (region_to) surface -------------------------
  n_skin = region(reg_to)%num_skin
  allocate(bound_node_2(n_skin))
  allocate(skin_node(n_skin))
  allocate(contact(id)%contact_group_end_addr2(num_contact_group))

  tol = 0.1*region(reg_to)%skin_length  ! tolerance based on minimum distance between two skin nodes
  num2 = 0
  do i = 1, num_contact_group 
    call find_surface_node(reg_to, n_skin, contact_pt_to(:,:,i), tol, group1, start_point, end_point)
    
    call find_skin(reg_to, group1, start_point, end_point, number, bound_node_2, skin_node, num2)
    num2 = num2 + number  
    contact(id)%contact_group_end_addr2(i) = num2
  end do

  allocate(contact(id)%contact_skin_target(num2))
  
  contact(id)%contact_skin_target(:) = bound_node_2(1:num2)
  contact(id)%num_contact_nodes(2) = num2

  if (contact_type == 2) call write_skin_contact(reg_to, bound_node_2(1:num2), num2) ! physical_domain

  deallocate(skin_node)


! -----------------------------------------------------
! allocate & initialize the relation matrices
  if (contact_type == 2) then
    size_of_col = 2*num1
    size_of_row = 2*num2
    allocate(contact(id)%u_g(size_of_col))   ! 2 dof per node
    contact(id)%u_g = 0.0
  
    allocate(contact(id)%relation(2*num2,2*num1))  ! 2 dof per node
    allocate(contact(id)%col_index(size_of_col+1))

    write(*,*) 'num1, num2 =', num1, num2
    contact(id)%relation = .FALSE.

    call set_contact_relation_2d(id)

  else
    STOP 'write_contact: wrong contact_type number!'
  endif


elseif (contact_dim == 3) then  ! 3D ---------------------------------------------------------------

  call set_contact_group(reg_from, num_contact_group)
  call set_contact_group(reg_to, num_contact_group)

! ball (region_from) surface ----------------------------
  n_skin = 4*region(reg_from)%num_skin_surfaces ! assume 4 nodes per one face
  allocate(bound_node_1(n_skin))
  allocate(contact(id)%contact_group_end_addr1(num_contact_group))

  num1 = 0
  do i = 1, num_contact_group
    corners = contact_pt_from(1:3,1:4,i)
    call find_contact_surface(reg_from, i, corners)

    call extract_nodes_from_faces(reg_from, i, n_skin, bound_node_1, number)

    num1 = num1 + number
    contact(id)%contact_group_end_addr1(i) = num1

  end do

  allocate(contact(id)%contact_skin_ball(num1))
  contact(id)%contact_skin_ball(:) = bound_node_1(1:num1)
  contact(id)%num_contact_nodes(1) = num1

  if (contact_type == 2) call write_skin_contact(reg_from, bound_node_1(1:num1), num1) ! physical_domain

! target (region_to) surface ------------------------------
  n_skin = 4*region(reg_to)%num_skin_surfaces ! assume 4 nodes per one face
  allocate(bound_node_2(n_skin))
  allocate(contact(id)%contact_group_end_addr2(num_contact_group))

  num2 = 0
  do i = 1, num_contact_group
    corners = contact_pt_to(1:3,1:4,i)
    call find_contact_surface(reg_to, i, corners)

    call extract_nodes_from_faces(reg_to, i, n_skin, bound_node_2, number)

    num2 = num2 + number
    contact(id)%contact_group_end_addr2(i) = num2

  end do

  allocate(contact(id)%contact_skin_target(num2))
  contact(id)%contact_skin_target(:) = bound_node_2(1:num2)
  contact(id)%num_contact_nodes(2) = num2

  allocate(contact(id)%contact_target_partner(2,num2))
  contact(id)%contact_target_partner = 0

  if (contact_type == 2) call write_skin_contact(reg_to, bound_node_2(1:num2), num2) ! physical_domain

! -----------------------------------------------------
! allocate & initialize the relation matrices
  if (contact_type == 2) then
    size_of_col = 3*num1
    size_of_row = 3*num2
    allocate(contact(id)%u_g(size_of_col))   ! 3 dof per node
    contact(id)%u_g = 0.0
  
    allocate(contact(id)%relation(size_of_row,size_of_col))  ! 3 dof per node
    allocate(contact(id)%col_index(size_of_col+1))

    write(*,*) '3D: num1, num2 =', num1, num2
    contact(id)%relation = .FALSE.

    call set_contact_relation_3d(id)

  else
    STOP 'write_contact: wrong contact_type number!'
  endif


endif  !--------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
if ((contact_type == 2) .AND. coarse_schur) then 

! for region_from (contact slave) -----------------------------------------------
  n_dofs = region(reg_from)%num_node_dof 
  allocate(node_dofs(n_dofs,num1))
  call read_dof_number(reg_from, bound_node_1(1:num1), node_dofs, n_dofs, num1)
  shape(1) = n_dofs*num1    ! size of the column vector
  shape(2) = 1
  allocate(moving_dofs(shape(1)))

  moving_dofs = 0
  count = 1
  do i = 1, num1
    do j = 1, n_dofs
      moving_dofs(count) = node_dofs(j,i)
      count = count + 1
    end do
  end do

  call dof_transformer(reg_from, -1, moving_dofs, shape(1))   ! module 'physical_domain'
  deallocate(node_dofs, moving_dofs)

  set_number = num_disp_sets + 2*id - 1  ! num_disp_sets: global variable in 'physical_domain'
  call set_disp_bc(set_number, reg_from)

! for region_to (contact master) -----------------------------------------------
  n_dofs = region(reg_to)%num_node_dof 
  allocate(node_dofs(n_dofs,num2))

  call read_dof_number(reg_to, bound_node_2(1:num2), node_dofs, n_dofs, num2)
  shape(1) = n_dofs*num2    ! size of the column vector
  shape(2) = 1

  allocate(moving_dofs(shape(1)))

  moving_dofs = 0
  count = 1
  do i = 1, num2
    do j = 1, n_dofs
      moving_dofs(count) = node_dofs(j,i)
      count = count + 1
    end do
  end do

  call dof_transformer(reg_to, -1, moving_dofs, shape(1))   ! module 'physical_domain'
  deallocate(node_dofs, moving_dofs)

  if (action == 2) then ! rebuild  --------------------------------------------------
    allocate(free_node(old_number))
    count = 0
    do i = 1, old_number
      bingo = .FALSE.
      do j = 1, num2
        if (old_bound_node(i) == bound_node_2(j)) then
          bingo = .TRUE.
          EXIT
        endif
      end do
      if (bingo == .FALSE.) then
        count = count + 1
        free_node(count) = old_bound_node(i)
      endif
    end do
    
    if (count > 0) then ! if there are freed nodes in 'region_to' contact part
      num2 = count
      allocate(node_dofs(n_dofs,num2))

      call read_dof_number(reg_to, free_node(1:num2), node_dofs, n_dofs, num2)
      shape(1) = n_dofs*num2    ! size of the column vector
      shape(2) = 1

      allocate(moving_dofs(shape(1)))

      moving_dofs = 0
      count = 1
      do i = 1, num2
        do j = 1, n_dofs
          moving_dofs(count) = node_dofs(j,i)
          count = count + 1
        end do
      end do

      call dof_transformer(reg_to, 1, moving_dofs, shape(1))   ! module 'physical_domain'
      deallocate(node_dofs, moving_dofs)
    endif ! count > 0
    
    deallocate(old_bound_node, free_node)    
  endif ! action == 2

  set_number = num_disp_sets + 2*id  ! num_disp_sets: global variable in 'physical_domain'
  call set_disp_bc(set_number, reg_to)

endif

!---------------------------------------------------------------------------------------------------

deallocate(bound_node_1, bound_node_2)


end subroutine write_contact







subroutine set_global_sparse(id)


implicit none

integer, intent(in) :: id

integer :: rn0, rn1, rn2, m, m0, m1, m2, nix(2), nix1, nix2, node1, i, j, k, l, number, old_number, st, ed, st0, ed0, global_num_dofs, node_dofs, pb_size
integer, allocatable :: inode(:,:), ix(:), ix_to_g_pointer(:)
logical, allocatable :: mask(:), checkboard(:)


do i = 1, 2 ! over two contact regions
  rn0 = contact(id)%contact_region(i)
  m0 = region(rn0)%num_dofs

  node_dofs = region(rn0)%num_node_dof

  nix(i) = node_dofs * contact(id)%num_contact_nodes(i)
  m = contact(id)%num_contact_nodes(i)
  allocate(inode(node_dofs,m))
  allocate(ix(nix(i)))  

  if (allocated(sdomain(rn0)%g_ix)) then
    deallocate(sdomain(rn0)%g_ix)
    allocate(sdomain(rn0)%g_ix(nix(i)))
  endif

	do j = 1, m
	  node1 = region(rn0)%skin_contact_nodes(j)
	  inode(1:node_dofs,j) = sdomain(rn0)%dof_index(region(rn0)%node(4:(node_dofs+3),node1))
	end do
	ix = pack(inode, .TRUE.)  

  sdomain(rn0)%g_ix = ix

  if (i == 1)  then
    contact(id)%dof_offset = m0
    allocate(ix_to_g_pointer(m0))
    ix_to_g_pointer = 0
    do j = 1, nix(i)
      ix_to_g_pointer(ix(j)) = j
    end do

  elseif (i == 2) then
    if (m0 <= nix(i)) STOP 'set_coarse_sparse: dofs must larger than contact dofs!'
    allocate(contact(id)%nonoverlap_dof(m0))

    contact(id)%nonoverlap_dof = .TRUE.
    contact(id)%nonoverlap_dof(ix) = .FALSE.
    allocate(contact(id)%old_number(m0 - nix(i)))
    l = 1
    do j = 1, m0
      if (contact(id)%nonoverlap_dof(j)) then
        contact(id)%old_number(l) = j
        l = l + 1
      endif
    end do
   endif

  deallocate(inode, ix)
end do  ! i

if (allocated(contact(id)%global_row_index)) then
  deallocate(contact(id)%global_row_index)
endif 


! checkboard update considering contact nodes in region_from (soft region)
rn1 = contact(id)%contact_region(1)
rn2 = contact(id)%contact_region(2)
nix1 = nix(1)
nix2 = nix(2)
m1 = region(rn1)%num_dofs
m2 = region(rn2)%num_dofs - nix2
m0 = m1 + 1

global_num_dofs = m1 + m2
contact(id)%num_dofs = global_num_dofs

write(*,*) '###### global number of dofs =', global_num_dofs
allocate(contact(id)%checkboard_gg(nix1,nix1))
allocate(contact(id)%checkboard_gi(nix1,m2))
allocate(contact(id)%checkboard_ig(m2,nix1))


contact(id)%checkboard_gg = .FALSE.
contact(id)%checkboard_gi = .FALSE.
contact(id)%checkboard_ig = .FALSE.

pb_size = sdomain(rn2)%problem_size
allocate(checkboard(pb_size))

do i = 1, nix1
  do j = 1, nix1
    if (ANY(contact(id)%relation(:,i)) .AND. ANY(contact(id)%relation(:,j))) then
      contact(id)%checkboard_gg(i,j) = .TRUE.
    endif
  end do
!  contact(id)%checkboard_gg((i+1):nix1,i) = contact(id)%checkboard_gg(i,(i+1):nix1) ! symmetric part
  
  do j = 1, m2
    k = contact(id)%old_number(j)

    call checkboard_sparse_to_full(rn2, k, checkboard, pb_size)
    contact(id)%checkboard_gi(i,j) = ANY(contact(id)%relation(:,i) .AND. checkboard(sdomain(rn2)%g_ix))
    contact(id)%checkboard_ig(j,i) = ANY(checkboard(sdomain(rn2)%g_ix) .AND. contact(id)%relation(:,i))
  end do
end do

deallocate(checkboard)

! write global row_index

if (allocated(contact(id)%global_row_index)) then
  deallocate(contact(id)%global_row_index)
  deallocate(contact(id)%global_col_addr)
  deallocate(contact(id)%global_K)
endif 
allocate(contact(id)%global_row_index(global_num_dofs+1))
allocate(mask(contact(id)%num_dofs))

pb_size = sdomain(rn1)%problem_size
allocate(checkboard(pb_size))

number = 1
contact(id)%global_row_index(1) = 1
do i = 1, m1
  if (ix_to_g_pointer(i) == 0) then ! purely-free dofs
    number = number + sdomain(rn1)%row_index(i+1) - sdomain(rn1)%row_index(i)
  else  ! contact-points dofs
    mask = .FALSE.
    k = ix_to_g_pointer(i)

    call checkboard_sparse_to_full(rn1, i, checkboard, pb_size)
    mask(1:m1) = checkboard(1:m1)
    mask(sdomain(rn1)%g_ix) = mask(sdomain(rn1)%g_ix) .OR. contact(id)%checkboard_gg(k,:)
    mask(m0:global_num_dofs) = contact(id)%checkboard_gi(k,:)

    number = number + COUNT(mask)
  endif
  contact(id)%global_row_index(i+1) = number
end do

deallocate(checkboard)
pb_size = sdomain(rn2)%problem_size
allocate(checkboard(pb_size))

do i = 1, m2
  k = contact(id)%old_number(i)

  call checkboard_sparse_to_full(rn2, k, checkboard, pb_size)
  number = number + COUNT(checkboard(contact(id)%old_number))
  number = number + COUNT(contact(id)%checkboard_ig(i,:))
    
  contact(id)%global_row_index(m1+i+1) = number
end do


! write global col_addr
number = number - 1
allocate(contact(id)%global_col_addr(number))
allocate(contact(id)%global_K(number))
contact(id)%global_K = 0.0  ! initialize just in case

deallocate(checkboard)
pb_size = sdomain(rn1)%problem_size
allocate(checkboard(pb_size))

do i = 1, m1
  if (ix_to_g_pointer(i) == 0) then ! purely-free dofs  
    if (sdomain(rn1)%row_existence(i)) then
      st0 = sdomain(rn1)%row_index(i)
      ed0 = sdomain(rn1)%row_index(i+1) - 1
      st = contact(id)%global_row_index(i)
      ed = contact(id)%global_row_index(i+1) - 1
      contact(id)%global_col_addr(st:ed) = sdomain(rn1)%col_addr(st0:ed0)
    endif
  else  ! contact-points dofs
    k = ix_to_g_pointer(i)

    call checkboard_sparse_to_full(rn1, i, checkboard, pb_size)
    st = contact(id)%global_row_index(i)
    mask(1:m1) = checkboard(1:m1)
    mask(sdomain(rn1)%g_ix) = mask(sdomain(rn1)%g_ix) .OR. contact(id)%checkboard_gg(k,:)
    mask(m0:global_num_dofs) = contact(id)%checkboard_gi(k,:)
    number = 0
    do j = 1, global_num_dofs
      if (mask(j)) then
        contact(id)%global_col_addr(st+number) = j
        number = number + 1
      endif
    end do  ! j
  endif
end do  ! i


deallocate(checkboard)
pb_size = sdomain(rn2)%problem_size
allocate(checkboard(pb_size))

do i = 1, m2
  k = contact(id)%old_number(i)

  call checkboard_sparse_to_full(rn2, k, checkboard, pb_size)
  st = contact(id)%global_row_index(m1+i)
  mask(1:m1) = .FALSE.
  mask(sdomain(rn1)%g_ix) = contact(id)%checkboard_ig(i,:)
  mask(m0:global_num_dofs) = checkboard(contact(id)%old_number)

  number = 0
  do j = 1, global_num_dofs
    if (mask(j)) then
      contact(id)%global_col_addr(st+number) = j
      number = number + 1
    endif
  end do  ! j  
end do  ! i

allocate(contact(id)%ix_to_g_pointer(m1))
contact(id)%ix_to_g_pointer = ix_to_g_pointer

deallocate(ix_to_g_pointer, mask)
deallocate(checkboard)

allocate(contact(id)%global_p(global_num_dofs))

! write(*,*) 'row:', contact(id)%global_row_index
! write(*,*) 'col:', contact(id)%global_col_addr


end subroutine set_global_sparse







subroutine global_sp_assemble(id)


implicit none

integer, intent(in) :: id

integer :: rn1, rn2, node_dofs, m, m22, m1, m2, nix1, nix2, i, j, k, number, old_number, ndim, st, ed, st0, ed0
real, allocatable :: vect(:), vect2(:), vect3(:), new_vect(:), p(:), k_bb(:,:), k_gg(:,:), k_gi(:,:), k_gi2(:,:), k_ig(:)


rn1 = contact(id)%contact_region(1)
rn2 = contact(id)%contact_region(2)

node_dofs = region(rn1)%num_node_dof
nix1 = node_dofs*contact(id)%num_contact_nodes(1)

node_dofs = region(rn2)%num_node_dof
nix2 = node_dofs*contact(id)%num_contact_nodes(2)
m1 = region(rn1)%num_dofs
m2 = region(rn2)%num_dofs 
m22 = m2 - nix2

allocate(k_gi(nix2,m22))
allocate(p(nix1))
allocate(vect2(m2))

do i = 1, nix2
  j = sdomain(rn2)%g_ix(i)
  call sparse_to_full(rn2, j, vect2, m2)
  k_gi(i,:) = PACK(vect2, contact(id)%nonoverlap_dof)
end do

ndim = contact(id)%col_index(nix1+1) - 1  ! Trans sparse length
allocate(k_gi2(nix1,m22))

! write(*,*) 'g_ix1 =', sdomain(rn1)%g_ix
! write(*,*) 'g_ix1 =', sdomain(rn2)%g_ix
! write(*,*) 'col_index =', contact(id)%col_index
! write(*,*) 'row_addr =', contact(id)%row_addr
! write(*,*) 'sparse_Trans =', contact(id)%sparse_Trans
k_gi2 = 0.0
call mul_tr_sparse_matrix(nix2, nix1, m22, ndim, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, k_gi, k_gi2)
! write(*,*) 'k_gi =', k_gi
! write(*,*) 'k_gi2 =', k_gi2
deallocate(k_gi)

allocate(k_gg(nix1,nix1))
allocate(k_bb(nix2,nix2))

do i = 1, nix2
  j = sdomain(rn2)%g_ix(i)
  call sparse_to_full(rn2, j, vect2, m2)
  k_bb(i,:) = vect2(sdomain(rn2)%g_ix)
end do  

k_gg = 0.0
call triple_mul_sparse(nix2, nix1, ndim, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, k_bb, k_gg)
deallocate(k_bb)

allocate(vect(m1+m22))

do i = 1, m1
  if (contact(id)%ix_to_g_pointer(i) == 0) then ! purely-free dofs
    if (sdomain(rn1)%row_existence(i)) then
      st0 = sdomain(rn1)%row_index(i)
      ed0 = sdomain(rn1)%row_index(i+1) - 1
      st = contact(id)%global_row_index(i)
      ed = contact(id)%global_row_index(i+1) - 1
      contact(id)%global_K(st:ed) = sdomain(rn1)%sparse_k(st0:ed0)
    endif
  else  ! contact-points dofs
    st = contact(id)%global_row_index(i)
    ed = contact(id)%global_row_index(i+1) - 1
    
    if (ed >= st) then
      k = contact(id)%ix_to_g_pointer(i)    
      vect = 0.0

      call sparse_to_full(rn1, i, vect(1:m1), m1)

      vect(sdomain(rn1)%g_ix) = vect(sdomain(rn1)%g_ix) + k_gg(k,:)
      vect(m1+1:m1+m22) = k_gi2(k,:)

      contact(id)%global_K(st:ed) = vect(contact(id)%global_col_addr(st:ed))
    endif
  endif
end do  ! i

deallocate(vect, k_gi2)

allocate(vect(m1))
allocate(vect3(nix2))
allocate(new_vect(m22))
allocate(k_ig(nix1))

do i = 1, m22
  st = contact(id)%global_row_index(m1+i)
  ed = contact(id)%global_row_index(m1+i+1) - 1
  if (ed >= st) then
    k = contact(id)%old_number(i)
    call sparse_to_full(rn2, k, vect2, m2)
    vect3 = vect2(sdomain(rn2)%g_ix)
    call mul_matrix_sparse(1, nix2, nix1, ndim, vect3, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, k_ig)
    number = COUNT(contact(id)%checkboard_ig(i,:))
      vect = 0.0
      vect(sdomain(rn1)%g_ix) = k_ig
      contact(id)%global_K(st:(st+number-1)) = vect(contact(id)%global_col_addr(st:(st+number-1)))

    if (ed >= (st+number)) then
      new_vect = PACK(vect2, contact(id)%nonoverlap_dof)
      contact(id)%global_K((st+number):ed) = new_vect(contact(id)%global_col_addr((st+number):ed) - m1)
    endif
  endif
end do  ! i

! assemble global load vector

contact(id)%global_p(1:m1) = sdomain(rn1)%p(1:m1)
p = 0.0
vect3 = sdomain(rn2)%p(sdomain(rn2)%g_ix)
call mul_tr_sparse_matrix(nix2, nix1, 1, ndim, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, vect3, p)
contact(id)%global_p(sdomain(rn1)%g_ix) = sdomain(rn1)%p(sdomain(rn1)%g_ix) + p
contact(id)%global_p(m1+1:m1+m22) = PACK(sdomain(rn2)%p, contact(id)%nonoverlap_dof)

deallocate(vect, vect2, vect3, new_vect, k_gg, k_ig)


end subroutine global_sp_assemble







subroutine update_global_u(id, u, size)

implicit none

integer, intent(in) :: id, size
real, intent(in) :: u(size)
integer :: rn0, rn1, rn2, i, m0, m1, m22, ndim, nix1, nix2, node_dofs
real, allocatable :: region_u(:), vect1(:), vect2(:)


rn0 = contact(id)%contact_region(1)
node_dofs = region(rn0)%num_node_dof
nix1 = node_dofs * contact(id)%num_contact_nodes(1)

rn0 = contact(id)%contact_region(2)
node_dofs = region(rn0)%num_node_dof
nix2 = node_dofs * contact(id)%num_contact_nodes(2)

allocate(vect1(nix1))
allocate(vect2(nix2))

do i = 1, 2
  rn0 = contact(id)%contact_region(i)
  m0 = region(rn0)%num_dofs

  allocate(region_u(m0))
  if (i == 1) then
    region_u = u(1:m0)
    m1 = m0
    rn1 = rn0
  elseif (i == 2) then
    rn2 = rn0
    ndim = contact(id)%col_index(nix1+1) - 1  ! Trans sparse length
    region_u(contact(id)%old_number) = u(m1+1:size)
    vect1 = u(sdomain(rn1)%g_ix)

    call mul_sparse_matrix(nix2, nix1, 1, ndim, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, vect1, vect2)
    region_u(sdomain(rn2)%g_ix) = vect2
  endif

  call update_u(rn0, region_u, m0)
  deallocate(region_u)
end do

deallocate(vect1, vect2)


end subroutine update_global_u








subroutine global_newton_solver(num_rn, max_itr, tolerance, total_error, status)

! rn: region number

implicit none

integer, intent(in) :: num_rn, max_itr
real, intent(in) :: tolerance
real, intent(out) :: total_error
integer, intent(inout) :: status    ! converge_status = {-1, 1, 2, 3 ...}

integer :: i, m1, m2, iter, size, size1, size2, rn, id, num1, contact_rn, ass_flag, contact_dim
integer :: debug_target(2)
real :: tol, tol_factor, error2, error, error_old= 0.0, tol_rn(num_rn), DNRM2
logical :: existence


if (.NOT. system_sparse) STOP 'global_newton_solver: currently only for sparse solver!'

tol_factor = 1.0
error2 = 0.0
tol = 0.0

iter = 0

iteration: do  !!!!!!!!!!!!!!!!!!!!!!

  total_error = 0.0
  !write(*,'(/A,I6,A,2ES16.8)') '>>> iter #:', iter, '  ---------------------------------------------'
  !write(*,'(A,I6)') ' % local problem: number of regions =', num_rn
  
  do  rn = 1, num_rn
  
    !write(*,'(A,I5)') '  >> region:', rn
    m1 = region(rn)%num_dofs     ! # of active dof
    m2 = region(rn)%num_constr   ! # of inactive dof
    size = sdomain(rn)%problem_size
    !write(*,'(A,2I8)') '     number of active dof & constrained dof =', m1, m2
    if (m1 < 1) STOP 'global_newton_solver: active dof is zero or less!'

    if (iter == 0) then
      call prescribed_u(rn, m2, iter)
      call assemble_p(rn)
    endif
    
    call update_system_changed(0, rn, .TRUE.) ! dummy value: .TRUE.
    
    if (sdomain(rn)%system_changed) call set_sparse(rn)
    call assemble_sparse_k(rn, 1, ass_flag, iter, debug_target)


    if (ass_flag < 0) then
      write(*,*) '>>>>>>>>>>> newton_solver: Element-level poor condition! @ rn=', rn
      status = -1
      EXIT iteration  !------------------------------------------------------------------------>>>>>
    endif  
        
    call update_system_changed(1, rn, .FALSE.)

    if (iter == 0) then
      error =  DNRM2(size, sdomain(rn)%p, 1)
      tol_rn(rn) = MAX(tolerance, error*tolerance)
    else
      if (iter == 1) then
        call energy_norm(sdomain(rn)%p, sdomain(rn)%u, size, error)
        tol_rn(rn) = MAX(tolerance, error*tolerance)
      endif
      call energy_norm(sdomain(rn)%p, sdomain(rn)%del_u, size, error)
    endif  

    !write(*,'(A,2ES16.8)') '   %% local problem error & tol:', error*tolerance*tol_factor/tol_rn(rn), tol_rn(rn)*tol_factor

    tol = tol + tol_rn(rn)    
    total_error = total_error + error
    
  end do  ! rn

!--------------------------------------------------------------------------------------------------

  tol = tol/num_rn
  total_error = total_error/num_rn
  !write(*,'(A,I6,A,2ES16.8)') '>>> iter #:', iter, '    error & tol: ', total_error*tolerance*tol_factor/tol, tolerance*tol_factor
  !write(*,'(A,I6,A,2ES16.8)') '>>> iter #:', iter, ' error_old & error: ', (10.0**(6-iter))*error_old**2, total_error
  if (total_error < tol) then
    if ((iter < MIN(5, max_itr - 1)) .AND. (1.0e2*error_old**2 > total_error)) then
      status = MAX(0, status) + 1
    !  write(*,*) '>>>>>>>>>>> newton_solver: Good convergency!', status
    else
      status = 1
    endif
    EXIT    !---------------------------------------------------------------------------->>>>>
  endif
  
  if ((iter > 2) .AND. (error_old/total_error < 0.1)) then
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
    
  do id = 1, number_of_contacts
    if (contact(id)%contact_type == 2) then

      call global_sp_assemble(id)
      call global_sp_solve(id, contact(id)%num_dofs, error2)

    endif
  end do
       
  iter = iter + 1  
  error_old = total_error
  if (iter == 1) tol = 0.0
  
end do iteration  !!!!!!!!!!!!!!!!!!!!!!

do rn = 1, num_rn
  sdomain(rn)%p0_trans_flag = .FALSE.
end do

!--------------------------------------------------------------------------------------------------

end subroutine global_newton_solver








subroutine coarse_newton_solver(num_rn, max_itr, tolerance, total_error, status)

! rn: region number

implicit none

integer, intent(in) :: num_rn, max_itr
real, intent(in) :: tolerance
real, intent(out) :: total_error
integer, intent(inout) :: status    ! converge_status = {-1, 1, 2, 3 ...}

integer :: i, m1, m2, iter, size, size1, size2, k11_size, rn, id, num1, contact_rn, coarse_flag, ass_flag, contact_dim
integer :: debug_target(2)
real :: tol, tol_factor, error2, error, error_old(num_rn), tol_rn(num_rn), tol_coarse, DNRM2
logical :: convergence(num_rn), last_call, existence, rn_goodness(num_rn)


tol_factor = 1.0
convergence = .FALSE.
error2 = 0.0
error_old = 0.0
last_call = .FALSE.
coarse_flag = 1
tol_coarse = tolerance
rn_goodness = .FALSE.

iter = 0

iteration: do  !!!!!!!!!!!!!!!!!!!!!!

  total_error = 0.0
  if (iter < 2) tol_rn = 0.0
  write(*,'(/A,I6,A,2ES16.8)') '>>> iter #:', iter, '  ---------------------------------------------'
  write(*,'(A,I6)') ' % local problem: number of regions =', num_rn
  
  do  rn = 1, num_rn
  
    write(*,'(A,I5)') '  >> region:', rn
    m1 = region(rn)%num_dofs     ! # of active dof
    m2 = region(rn)%num_constr   ! # of inactive dof
    size = sdomain(rn)%problem_size
    write(*,'(A,2I8)') '     number of active dof & constrained dof =', m1, m2
    if (m1 < 1) STOP 'coarse_newton_solver: active dof is zero or less!'

    call prescribed_u(rn, m2, iter) ! considering prescribed displacement along contact boundaries

    if (iter == 0) then
      call assemble_p(rn)
    endif
    
    call update_system_changed(0, rn, .TRUE.) ! dummy value: .TRUE.
    
    if (system_sparse) then
      if (sdomain(rn)%system_changed) call set_sparse(rn)
      call assemble_sparse_k(rn, coarse_flag, ass_flag, iter, debug_target)
    else
      call assemble_k(rn, coarse_flag, ass_flag, iter, debug_target)
    endif

    if (ass_flag < 0) then
      write(*,*) '>>>>>>>>>>> newton_solver: Element-level poor condition! @ rn=', rn
      status = -1
      EXIT iteration  !------------------------------------------------------------------------>>>>>
    endif  
        
    call update_system_changed(1, rn, .FALSE.)

    if (iter == 0) then
      error =  DNRM2(size, sdomain(rn)%p, 1)
      tol_rn(rn) = max(tolerance, error*tolerance)
      tol_coarse = tol_coarse + tol_rn(rn)
    else
      if (iter == 1) then
        call energy_norm(sdomain(rn)%p, sdomain(rn)%u, size, error)
        tol_rn(rn) = max(tolerance, error*tolerance)
      endif
      call energy_norm(sdomain(rn)%p, sdomain(rn)%del_u, size, error)
    endif  

    write(*,'(A,2ES16.8)') '  %%  local problem error & tol:', error*tolerance*tol_factor/tol_rn(rn), tol_rn(rn)*tol_factor


    tol = tol_rn(rn)

    if (error < tol) then
      convergence(rn) = .TRUE.
      if ((iter < MIN(5, max_itr - 1)) .AND. (1.0e2*error_old(rn)**2 > error)) then
        rn_goodness(rn) = .TRUE.
      endif
    else
      if ((iter > 2) .AND. (error_old(rn)/error < 0.1)) then
        write(*,*) '>>>>>>>>>>> newton_solver: Poor convergency! @ rn=', rn
        status = -1
        EXIT iteration  !---------------------------------------------------------------------->>>>>
      endif
      convergence(rn) = .FALSE.
    endif
    
    total_error = total_error + error
    error_old(rn) = error

  end do  ! rn

!--------------------------------------------------------------------------------------------------

  if (iter == 0) then
    tol_coarse = sqrt(tol_coarse)

  elseif (iter > 0) then

    if (iter == 1) then
      tol_coarse = max(sqrt(abs(tolerance)), error2*sqrt(abs(tolerance)))
    endif

    write(*,'(/A,2ES16.8)') '  %% coarse problem energy norm & tol:', error2, tol_coarse

    if (error2 < tol_coarse) then
      last_call = .TRUE.
      coarse_flag = 0
      write(*,'(A)') '  %% coarse problem boundary is locked!'
    endif
  
    if (last_call .AND. ALL(convergence)) then
      if (ALL(rn_goodness)) then
        status = MAX(0, status) + 1
        write(*,*) '>>>>>>>>>>> newton_solver: Good convergency!', status
      else
        status = 1
      endif
      EXIT  ! converged! ---------------------------------------------------------------------->>>>>
    endif
  endif
  
  
  if (last_call) then
    if (iter == max_itr) then
      tol_factor = 1.0/sqrt(abs(tolerance))
      do rn = 1, num_rn
        tol_rn(rn) = tol_factor*tol_rn(rn)  ! down the original tolerance
      end do
    elseif (iter > max_itr) then
      status = -1
      EXIT   !--------------------------------------------------------------------------------->>>>>
!      STOP 'coarse_newton_solver: exceed maximum iteration!'
    endif
  elseif (iter > max_itr) then
    status = -1
    EXIT   !----------------------------------------------------------------------------------->>>>>
!    STOP 'coarse_newton_solver: exceed maximum iteration!'
  else  ! '.not. last_call' case

    error2 = 0.0
    do id = 1, number_of_contacts
      if (contact(id)%contact_type == 2) then
        contact_dim = contact(id)%contact_dim
        size1 = contact_dim * contact(id)%num_contact_nodes(1)    ! contact_dim-dof per node
        size2 = contact_dim * contact(id)%num_contact_nodes(2)    ! contact_dim-dof per node
        
        call coarse_handler(id, size1, size2, error)

        error2 = error2 + error
      endif
    end do    
    do rn = 1, num_rn    
      m1 = region(rn)%num_dofs     ! # of active dof

      if (system_sparse) then  
        call sp_solve(rn, m1)      
      else
        call general_solve(rn, m1)
      endif
    end do

  endif
  
  iter = iter + 1  
  if (iter == 1) tol = 0.0
  
end do iteration  !!!!!!!!!!!!!!!!!!!!!!

do rn = 1, num_rn
  sdomain(rn)%p0_trans_flag = .FALSE.
end do

!--------------------------------------------------------------------------------------------------


end subroutine coarse_newton_solver








subroutine global_sp_solve(id, size, error_norm)


! m1: size of active_dof

use mkl_dss   ! MKL header module included in 'xfem_math.f'

implicit none

integer, intent(in) :: id, size
real, intent(out) :: error_norm
integer :: k11_size, perm(1), error
real :: pp(size), uu(size), time1, time2
type (MKL_DSS_HANDLE) :: handle 


k11_size = contact(id)%global_row_index(size+1) - 1
perm(1) = 0

pp = contact(id)%global_p
uu = 0.0

! write(*,*) ' p =', pp

! call cpu_time(time1)

error = dss_create(handle, MKL_DSS_DEFAULTS)
error = dss_define_structure(handle, MKL_DSS_SYMMETRIC_STRUCTURE, contact(id)%global_row_index, size, size, contact(id)%global_col_addr, k11_size)

error = dss_reorder(handle, MKL_DSS_METIS_OPENMP_ORDER, perm)
error = dss_factor(handle, MKL_DSS_DEFAULTS, contact(id)%global_K)
error = dss_solve(handle, MKL_DSS_DEFAULTS, pp, 1, uu)
error = dss_delete(handle, MKL_DSS_DEFAULTS) 

! call cpu_time(time2)
! write(*,*) '########## sparse solving cpu time', time2-time1

! write(*,*) 'del u =', uu

call update_global_u(id, uu, size)

call energy_norm(pp, uu, size, error_norm)

!write(*,*) '   %% error_norm =', error_norm

end subroutine global_sp_solve







subroutine coarse_handler(id, m1, m2, error)

! solve_flag = 1, 3: Symmetric system
!            = 2, 4: General system

implicit none

integer, intent(in) :: id, m1, m2
real, intent(out) :: error
integer :: reg_ball, reg_target, contact_type, ndim
real :: k_gg(m1,m1), p_g(m1), del_u_g(m1), del_u_g2(m2)


contact_type = contact(id)%contact_type
reg_ball = contact(id)%contact_region(1)
reg_target = contact(id)%contact_region(2)

contact(id)%u_g = sdomain(reg_ball)%u(sdomain(reg_ball)%g_ix)

call assemble_coarse(id, reg_ball, reg_target, k_gg, p_g, m1, m2)
call general_lu(k_gg, p_g, del_u_g, m1, 1)

! write(*,*) ' coarse del_u:', del_u_g

ndim = contact(id)%col_index(m1+1) - 1
call mul_sparse_matrix(m2, m1, 1, ndim, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, del_u_g, del_u_g2)

if (coarse_schur) call revise_local_r(id, reg_ball, reg_target, del_u_g, del_u_g2, m1, m2)

call update_coarse_u(id, reg_ball, reg_target, del_u_g, m1, m2)

call energy_norm(p_g, del_u_g, m1, error)

end subroutine coarse_handler







subroutine assemble_coarse(id, rn1, rn2, k_gg, p_g, m1, m2)

! rn: region number

implicit none 

integer, intent(in) :: id, rn1, rn2, m1, m2
real, intent(out) :: k_gg(m1,m1), p_g(m1)
integer :: i, j, ndim, n1, n2, ix(m1), ix2(m2), psize1, psize2
real :: mat(m1,m2), mat2(m2,m1), p_g2(m2)
real, allocatable :: k_gi(:,:), vector(:)


k_gg = sdomain(rn1)%k_bb
ndim = contact(id)%col_index(m1+1) - 1

call triple_mul_sparse(m2, m1, ndim, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, sdomain(rn2)%k_bb, k_gg)

if (coarse_schur) then   ! coarse_schur: from 'system domain'
  p_g = sdomain(rn1)%p_b
  p_g2 = sdomain(rn2)%p_b
else
  n1 = region(rn1)%num_dofs
  n2 = region(rn2)%num_dofs 
  psize1 = sdomain(rn1)%problem_size
  psize2 = sdomain(rn2)%problem_size
  ix = sdomain(rn1)%g_ix
  ix2 = sdomain(rn2)%g_ix

  if (system_sparse) then
    allocate(k_gi(m1,n1), vector(psize1))
    do i = 1, m1
      j = ix(i)
      call bc_sparse_to_full(rn1, j, vector, psize1)
      k_gi(i,:) = vector(1:n1)
    end do
    p_g = sdomain(rn1)%p_b - MATMUL(k_gi, sdomain(rn1)%del_u(1:n1))  
    deallocate(k_gi, vector)
  
    allocate(k_gi(m2,n2), vector(psize2)) 
    do i = 1, m2
      j = ix2(i)
      call bc_sparse_to_full(rn2, j, vector, psize2)
      k_gi(i,:) = vector(1:n2)
    end do
    p_g2 = sdomain(rn2)%p_b - MATMUL(k_gi, sdomain(rn2)%del_u(1:n2))
    deallocate(k_gi, vector)
  else
    p_g = sdomain(rn1)%p_b - MATMUL(sdomain(rn1)%k(ix,1:n1), sdomain(rn1)%del_u(1:n1))
    p_g2 = sdomain(rn2)%p_b - MATMUL(sdomain(rn2)%k(ix2,1:n2), sdomain(rn2)%del_u(1:n2))
  endif
endif


call mul_tr_sparse_matrix(m2, m1, 1, ndim, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, p_g2, p_g)

! write(*,*) ' coarse p_g:', p_g

end subroutine assemble_coarse







subroutine update_coarse_u(id, rn1, rn2, del_u_g, n1, n2)

! rn: region number
! n1, n2: size of active_dof
! u, del_u: module global variables

implicit none

integer, intent(in) :: id, rn1, rn2, n1, n2
real, intent(in) :: del_u_g(n1)
integer :: sn, ix1(n1), ix2(n2), ndim
real :: vect(n2), u_g(n1)


u_g = contact(id)%u_g + del_u_g
sn = num_disp_sets + 2*id - 1  ! num_disp_sets: global variable in 'physical_domain'
ix1 = sdomain(rn1)%g_ix

!write(*,*) 'rn1:', u_g
call update_disp_bc(sn, rn1, ix1, u_g, n1, 2)

sn = num_disp_sets + 2*id  ! num_disp_sets: global variable in 'physical_domain'
ix2 = sdomain(rn2)%g_ix

ndim = contact(id)%col_index(n1+1) - 1
call mul_sparse_matrix(n2, n1, 1, ndim, contact(id)%col_index, contact(id)%row_addr, contact(id)%sparse_Trans, u_g, vect)

call update_disp_bc(sn, rn2, ix2, vect, n2, 2)

contact(id)%u_g = u_g

! write(*,*) 'u_g:', u_g

end subroutine update_coarse_u







subroutine set_contact_relation_3d(id)

! July, 2015 by Jeeho Lee
!   3D version

implicit none

integer, intent(in) :: id

integer :: ball_rn, target_rn, n_skin_ball, n_skin_ball_st, n_skin_target, n_skin_target_st
integer :: i, j, k, n, l, jj, target_node
integer :: number, num_contact_group
integer :: n_faces, face_nodes(4), elem_num, face_num
real :: s_face(4), target_point(3), corners(3,4), areas(4)
integer, allocatable :: contact_face(:)
real, allocatable :: matrix(:,:), contact_trans(:,:,:)
logical :: in_or_out, skin_surface_trans_flag


ball_rn = contact(id)%contact_region(1)
target_rn = contact(id)%contact_region(2)
num_contact_group = contact(id)%num_contact_group

n_skin_ball = contact(id)%contact_group_end_addr1(num_contact_group)
n_skin_target = contact(id)%contact_group_end_addr2(num_contact_group)

allocate(matrix(3*n_skin_target,3*n_skin_ball))  ! 3 dof per node ----------------------------------
matrix = 0.0    ! zero-initialize just in case

do k = 1, num_contact_group

  n_skin_ball = contact(id)%contact_group_end_addr1(k)
  if (k == 1) then
    n_skin_ball_st = 1
  else
    n_skin_ball_st = contact(id)%contact_group_end_addr1(k-1) + 1
  endif
  
  n_skin_target = contact(id)%contact_group_end_addr2(k)
  if (k == 1) then
    n_skin_target_st = 1
  else
    n_skin_target_st = contact(id)%contact_group_end_addr2(k-1) + 1
  endif

  n_faces = COUNT(region(ball_rn)%contact_face_participation(:,k))
  allocate(contact_face(n_faces))

  n = 0
  do i = 1, region(ball_rn)%num_skin_surfaces
    if (region(ball_rn)%contact_face_participation(i,k)) then
      n = n + 1
      contact_face(n) = i
    endif
  end do  ! i

  skin_surface_trans_flag = .FALSE.
  allocate(contact_trans(3,3,n_faces))

  do i = n_skin_target_st, n_skin_target  ! over target nodes
    target_node = contact(id)%contact_skin_target(i)
    target_point = region(target_rn)%node_coord(1:3,target_node)

    do j = 1, n_faces ! over ball surface faces considering 'num_contact_group'
      elem_num = region(ball_rn)%skin_surfaces(2,contact_face(j))
      face_num = region(ball_rn)%skin_surfaces(1,contact_face(j))
      face_nodes(1:4) = region(ball_rn)%element_face(1:4,face_num,elem_num)
      corners = region(ball_rn)%node_coord(1:3,face_nodes)

      if (.NOT. skin_surface_trans_flag) call surface_trans_matrix(corners, contact_trans(:,:,j))

      call point_position(ball_rn, corners, contact_trans(:,:,j), target_point, s_face, areas, in_or_out)

      if (in_or_out) then
        contact(id)%contact_target_partner(1,i) = face_num
        contact(id)%contact_target_partner(2,i) = elem_num

!       expansion from node to DOF (3 dof per node assumed)
        do n = 1, 4
          jj = 0
          do l = n_skin_ball_st, n_skin_ball
            if (contact(id)%contact_skin_ball(l)==face_nodes(n))  jj = l
          end do  ! l
          if (jj > 0) then
            contact(id)%relation(3*i-2,3*jj-2) = .TRUE.
            contact(id)%relation(3*i-1,3*jj-1) = .TRUE.
            contact(id)%relation(3*i,3*jj)     = .TRUE.

            matrix(3*i-2,3*jj-2) = s_face(n)
            matrix(3*i-1,3*jj-1) = s_face(n)
            matrix(3*i,3*jj)     = s_face(n)
          endif
        end do
                        
        EXIT  ! out of loop for j  --------------------------------------->>
      endif
    
    end do  ! j  

    skin_surface_trans_flag = .TRUE.

  end do  ! i

  deallocate(contact_face, contact_trans)

!  assign s_factor for outsider_nodes ----------------------------------------------------
!  call add_outsider_3d(id, k, matrix)

end do ! k


! compress 'matrix' array into sparse form ---------------------------------------------------------

number = 1
contact(id)%col_index(1) = 1
do j = 1, 3*n_skin_ball ! over column
  number = number + COUNT(contact(id)%relation(:,j))
  contact(id)%col_index(j+1) = number  
end do  ! j

number = number - 1
allocate(contact(id)%row_addr(number))
allocate(contact(id)%sparse_Trans(number))  

do j = 1, 3*n_skin_ball ! over column
  number = contact(id)%col_index(j)
  do i = 1, 3*n_skin_target ! over row
    if (contact(id)%relation(i,j)) then
      contact(id)%row_addr(number) = i
      contact(id)%sparse_Trans(number) = matrix(i,j)
      number = number + 1
    endif
  end do  ! i
end do  ! j

!write(*,'(A/,(12L8/))') 'relation =', TRANSPOSE(contact(id)%relation)
!write(*,'(A/,(12F8.2/))') 'matrix =', TRANSPOSE(matrix)
deallocate(matrix)


end subroutine set_contact_relation_3d







subroutine add_outsider_3d(id, contact_group, matrix)

! August, 2015 by Jeeho Lee
!   3D version

implicit none

integer, intent(in) :: id, contact_group
real, intent(inout) :: matrix(:,:)

integer :: ball_rn, target_rn, n_skin_ball, n_skin_ball_st, n_skin_target, n_skin_target_st
integer :: i, j, k, n, l, jj, ll, target_node
integer :: n_faces, face_nodes(4), elem_num, face_num, num_nodes, order(4), face_nodes0(4)
integer :: outsider, outsider_right, outsider_left, outsider_opp, ball_elem_num, ball_face_num, ball_face_nodes(4)
real :: s_face(4), corners(3,4), trans(3,3), ball_points(3,4), outsider_coord(3), outsider_s, areas(4)
real :: DNRM2, distance(4), vector(3)
integer, allocatable :: contact_face(:), node_to_node(:)

logical :: in_or_out, skin_surface_trans_flag, outsider_flag
integer, parameter :: circle_order(0:6) = [4, 1, 2, 3, 4, 1, 2]


k = contact_group

ball_rn = contact(id)%contact_region(1)
target_rn = contact(id)%contact_region(2)
num_nodes = region(target_rn)%num_nodes

n_skin_ball = contact(id)%contact_group_end_addr1(k)
if (k == 1) then
    n_skin_ball_st = 1
else
    n_skin_ball_st = contact(id)%contact_group_end_addr1(k-1) + 1
endif
  
n_skin_target = contact(id)%contact_group_end_addr2(k)
if (k == 1) then
    n_skin_target_st = 1
else
    n_skin_target_st = contact(id)%contact_group_end_addr2(k-1) + 1
endif

allocate(node_to_node(num_nodes))
node_to_node = 0

do i = n_skin_target_st, n_skin_target  ! over target nodes
  target_node = contact(id)%contact_skin_target(i)
  node_to_node(target_node) = i ! store a contact node number in the target region node number array
end do  ! i

n_faces = COUNT(region(target_rn)%contact_face_participation(:,k))
allocate(contact_face(n_faces))

n = 0
do i = 1, region(ball_rn)%num_skin_surfaces
  if (region(target_rn)%contact_face_participation(i,k)) then
    n = n + 1
    contact_face(n) = i
  endif
end do  ! i

do i = 1, n_faces ! over all target surfaces
  elem_num = region(target_rn)%skin_surfaces(2,contact_face(i))
  face_num = region(target_rn)%skin_surfaces(1,contact_face(i))
  face_nodes0(1:4) = region(target_rn)%element_face(1:4,face_num,elem_num)
  face_nodes = node_to_node(face_nodes0)
  outsider_flag = ANY(contact(id)%contact_target_partner(2,face_nodes(1:4)) < 1)

  if (outsider_flag) then
!    write(*,'(A,4I7)') '   >> outsider nodes (pdomain node number system) =', face_nodes0
!    write(*,'(A,4I7/)') '   >> outsider nodes (contact node number system) =', face_nodes
!    corner points of a considered target surface
    corners = region(target_rn)%node_coord(1:3,region(target_rn)%element_face(1:4,face_num,elem_num))
    trans = 0.0
    call surface_trans_matrix(corners, trans)

    do j = 1, 4 ! over all target nodes of a outsider surface
      outsider = face_nodes(j)
      outsider_right = face_nodes(circle_order(j+1))
      outsider_left  = face_nodes(circle_order(j-1))
      outsider_opp   = face_nodes(circle_order(j+2))

      if (contact(id)%contact_target_partner(2,outsider) < 1) then

        if (contact(id)%contact_target_partner(2,outsider_right) > 0) then
          n = outsider_right
        elseif (contact(id)%contact_target_partner(2,outsider_left) > 0) then
          n = outsider_left
        elseif (contact(id)%contact_target_partner(2,outsider_opp) > 0) then
          n = outsider_opp
        else
          n = 0
          EXIT  !------------------------------------>>>
        endif

!         write(*,'(A,2I7)') '   >> outsider and friend nodes:', outsider, n

        ball_elem_num = contact(id)%contact_target_partner(2,face_nodes(n))
        ball_face_num = contact(id)%contact_target_partner(1,face_nodes(n))
        ball_face_nodes(1:4) = region(ball_rn)%element_face(1:4,ball_face_num,ball_elem_num)
        ball_points = region(ball_rn)%node_coord(1:3,ball_face_nodes)

!         write(*,'(A,4I7)') '   >> outsider ball face nodes:', ball_face_nodes

        outsider_coord = region(target_rn)%node_coord(1:3,face_nodes0(j))

        do l = 1, 4
          vector = ball_points(1:3,l)-outsider_coord(1:3)
          distance(l) = DNRM2(3,vector,1)
        end do  ! l

!         write(*,'(A,4(3F10.3/))') 'ball_points =', ball_points
!         write(*,'(A,4(3F10.3/))') 'outsider_coord =', outsider_coord
!         write(*,'(A,4F10.3)') 'distance(1:4) =', distance

        call ranking_by_value(1, 4, distance, order)

        do l = 1, 4 ! over all ball nodes just in case that one is outside of the target surface
!           write(*,'(A,5I7/)') '   >> closest ball surface node:', order(l), order

          call point_position(target_rn, corners, trans, ball_points(1:3,order(l)), s_face, areas, in_or_out)
          if (in_or_out) then
            outsider_s = s_face(j)
            ll = order(l)
            EXIT !----------------------------------->>>
          endif
        end do  ! l

        jj = 0
        do l = n_skin_ball_st, n_skin_ball
          if (contact(id)%contact_skin_ball(l)==ball_face_nodes(ll))  jj = l
        end do  ! l
        if (jj < 1) STOP 'set_contact_relation_3d: fail to find ball node in outsider procedure!'

        if (outsider_s > matrix(3*outsider,3*jj)) then
          contact(id)%relation(3*outsider-2,3*jj-2) = .TRUE.
          contact(id)%relation(3*outsider-1,3*jj-1) = .TRUE.
          contact(id)%relation(3*outsider,3*jj)     = .TRUE.

          matrix(3*outsider-2,3*jj-2) = outsider_s
          matrix(3*outsider-1,3*jj-1) = outsider_s
          matrix(3*outsider,3*jj)     = outsider_s

          contact(id)%contact_target_partner(1,outsider) = -ball_face_num
          contact(id)%contact_target_partner(2,outsider) = -ball_elem_num
        endif

      endif

    end do  ! j
  endif ! outsider_flag
end do  ! i

deallocate(contact_face, node_to_node)


end subroutine add_outsider_3d







subroutine set_contact_relation_2d(id)

! Jan, 2013 by Jeeho Lee
! November, 2014:
!   numerical_zero: parameter from physical_domain
!   Sparse data format: CSC
!   Fast & efficient algorithm version

implicit none

integer, intent(in) :: id

integer :: ball_rn, target_rn, n_skin_ball, n_skin_ball_st, n_skin_target, n_skin_target_st
integer :: i, j, k, node_a, node_b, target_node
integer :: j_memory, number, num_contact_group
real :: s, l, d, left, right
real, allocatable :: matrix(:,:)


ball_rn = contact(id)%contact_region(1)
target_rn = contact(id)%contact_region(2)
num_contact_group = contact(id)%num_contact_group

n_skin_ball = contact(id)%contact_group_end_addr1(num_contact_group)
n_skin_target = contact(id)%contact_group_end_addr2(num_contact_group)

allocate(matrix(2*n_skin_target,2*n_skin_ball))  ! 2 dof per node ----------------------------------
matrix = 0.0    ! zero-initialize just in case

do k = 1, num_contact_group

  n_skin_ball = contact(id)%contact_group_end_addr1(k)
  if (k == 1) then
    n_skin_ball_st = 1
  else
    n_skin_ball_st = contact(id)%contact_group_end_addr1(k-1) + 1
  endif
  
  n_skin_target = contact(id)%contact_group_end_addr2(k)
  if (k == 1) then
    n_skin_target_st = 1
  else
    n_skin_target_st = contact(id)%contact_group_end_addr2(k-1) + 1
  endif

! contact(id)%relation = .FALSE.   <--- already initialized outside

  j_memory = n_skin_ball - 1

  do i = n_skin_target_st, n_skin_target
    target_node = contact(id)%contact_skin_target(i)

    do j = j_memory, n_skin_ball_st, -1

      node_a = contact(id)%contact_skin_ball(j)
      node_b = contact(id)%contact_skin_ball(j+1)
  
      call line_position(ball_rn, node_a, node_b, target_rn, target_node, s, l, d)
      
      if (s >= 0.0 .AND. s <= 1.0) then
        if (d <= 0.1*l) then
          left = 1.0 - s    ! influence factor (1-s)
          right = s         ! influence factor 
!      simple expansion from node to DOF (2 dof per node assumed)
          contact(id)%relation(2*i-1,2*j-1) = .TRUE.
          contact(id)%relation(2*i-1,2*j+1) = .TRUE.
          contact(id)%relation(2*i,2*j)     = .TRUE.
          contact(id)%relation(2*i,2*j+2)   = .TRUE.
          matrix(2*i-1,2*j-1) = left
          matrix(2*i-1,2*j+1) = right
          matrix(2*i,2*j)     = left
          matrix(2*i,2*j+2)   = right
                        
          j_memory = j  ! to avoid unnecessary trials
          EXIT  ! out of loop for j  --------------------------------------->>
        endif
      endif
    
    end do  ! j  
  end do  ! i

!  assign s_factor for outsider_nodes ----------------------------------------------------
!  call add_outsider_2d(id, k, matrix)
end do ! k


! compress 'matrix' array into sparse form ---------------------------------------------------------

number = 1
contact(id)%col_index(1) = 1
do j = 1, 2*n_skin_ball ! over column
  number = number + COUNT(contact(id)%relation(:,j))
  contact(id)%col_index(j+1) = number  
end do  ! j

number = number - 1
allocate(contact(id)%row_addr(number))
allocate(contact(id)%sparse_Trans(number))  

do j = 1, 2*n_skin_ball ! over column
  number = contact(id)%col_index(j)
  do i = 1, 2*n_skin_target ! over row
    if (contact(id)%relation(i,j)) then
      contact(id)%row_addr(number) = i
      contact(id)%sparse_Trans(number) = matrix(i,j)
      number = number + 1
    endif
  end do  ! i
end do  ! j

!write(*,'(A/,(6L8/))') 'relation =', TRANSPOSE(contact(id)%relation)
!write(*,'(A/,(6F8.2/))') 'matrix =', TRANSPOSE(matrix)
deallocate(matrix)


end subroutine set_contact_relation_2d








subroutine line_position(surface_region, node_a, node_b, point_region, point_node, s, l, d)

! Revised: March 14, 2013
! node_a, node_b: nodes on surface_region
! point_node: point to be measured in point_region
! s: proportional length of point_node along the given line (0 <= s <= 1) in [node_a, node_b]
! l: length of face (node_b, node_a)
! d: projectional foot length

implicit none

integer, intent(in) :: surface_region, node_a, node_b, point_region, point_node
real, intent(out) :: s, l, d
real :: point_a(2), point_b(2), test(2), l1, l2, calc_length, vector1(2)


point_a = region(surface_region)%node_coord(:,node_a)
point_b = region(surface_region)%node_coord(:,node_b)
test = region(point_region)%node_coord(:,point_node)

vector1 = point_b - point_a
if (DOT_PRODUCT(vector1, vector1) < 1.0e-10) then
  write(*,*) 'surface_region, node_a, node_b, point_region, point_node =', surface_region, node_a, node_b, point_region, point_node
  STOP 'line_position: same points!'
endif

call projection_s(point_a, point_b, test, s)
l = calc_length(point_b, point_a)
l1 = s*l
l2 = calc_length(test, point_a)
d = sqrt(abs(l2*l2 - l1*l1))

end subroutine line_position







subroutine add_outsider_2d(id, contact_group, matrix)

! August, 2015 by Jeeho Lee
!   2D version

implicit none

integer, intent(in) :: id, contact_group
real, intent(inout) :: matrix(:,:)

integer :: ball_rn, target_rn, n_skin_ball, n_skin_ball_st, n_skin_target, n_skin_target_st
integer :: i, j, k, node_a, node_b, target_node
real :: s, l, d, left, right
logical :: outsider_flag

k = contact_group

ball_rn = contact(id)%contact_region(1)
target_rn = contact(id)%contact_region(2)

n_skin_ball = contact(id)%contact_group_end_addr1(k)
if (k == 1) then
    n_skin_ball_st = 1
else
    n_skin_ball_st = contact(id)%contact_group_end_addr1(k-1) + 1
endif
  
n_skin_target = contact(id)%contact_group_end_addr2(k)
if (k == 1) then
  n_skin_target_st = 1
else
  n_skin_target_st = contact(id)%contact_group_end_addr2(k-1) + 1
endif


do i = n_skin_target_st, n_skin_target

  outsider_flag = .NOT.(ANY(contact(id)%relation(2*i,:)))

  if (outsider_flag) then
    target_node = contact(id)%contact_skin_target(i)

    do j = (n_skin_ball - 1), n_skin_ball_st, -1
      node_a = contact(id)%contact_skin_ball(j)
      node_b = contact(id)%contact_skin_ball(j+1)
  
      call line_position(ball_rn, node_a, node_b, target_rn, target_node, s, l, d)
      
      if ((i == n_skin_target_st) .AND. (s >= 1.0)) then
!        left = MAX((2.0 - s), 0.0)    ! influence factor
        left = 1.0
        contact(id)%relation(2*i-1,2*j+1) = .TRUE.
        contact(id)%relation(2*i,2*j+2)   = .TRUE.
        matrix(2*i-1,2*j+1) = left
        matrix(2*i,2*j+2)   = left
                        
        EXIT  ! out of loop for j  --------------------------------------->>
      elseif ((i == n_skin_target) .AND. (j == n_skin_ball_st) .AND. (s <= 0.0)) then
!        right = MAX((1.0 + s), 0.0)    ! influence factor     
        right = 1.0
        contact(id)%relation(2*i-1,2*j-1) = .TRUE.
        contact(id)%relation(2*i,2*j)     = .TRUE.
        matrix(2*i-1,2*j-1) = right
        matrix(2*i,2*j)     = right
          
        EXIT  ! out of loop for j  --------------------------------------->>
      endif
    end do  ! j

  endif ! outsider_flag

end do  ! i


end subroutine add_outsider_2d


!==============================================================================

subroutine free_coarse_domain

implicit none 

deallocate ( contact )

end subroutine 

!===============================================================================

end module coarse_domain
