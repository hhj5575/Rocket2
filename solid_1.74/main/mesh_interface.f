subroutine mesh_interface(num_domains, num_contacts, num_materials, remesh)

! External subroutine

use system_info_house
use file_info_house
use ale_domain
use remesh_domain

implicit none

integer, intent(in) :: num_domains, num_contacts, num_materials, remesh
integer :: i, status, ng, max_dim
character(len=5) :: keyword
logical :: existence

max_dim = maxval(region(1:num_domains)%num_dim)

inquire (file=inte_file, EXIST=existence)
if (existence .and. max_dim == 3) then
    do 
        read(unit_num_inte,*) keyword
        backspace(unit_num_inte)
        if (keyword == '*divd') then
            read(unit_num_inte,*) keyword, divided_mesh_region
        elseif (keyword == '*subd') then
            read(unit_num_inte,*) keyword, sub_mesh_region
        elseif (keyword == '*rety') then
            read(unit_num_inte,*) keyword, remesh_type_number
        elseif (keyword == '*redi') then
            read(unit_num_inte,*) keyword, remesh_dis
        elseif (keyword == '*p2cn') then
            exit
        endif
    enddo
else
    divided_mesh_region = 1
    sub_mesh_region = 1
endif
! interface file read
allocate(ale_set(num_domains))
if (max_dim == 2) then
    call create_mesh_set(unit_num_elem, num_domains, num_materials)
elseif (max_dim == 3) then
    call create_mesh_set(unit_num_elem, divided_mesh_region, num_materials)
    call create_remesh_domain(sub_mesh_region, 3)
endif

do i = 1, num_domains
    !if (region(i)%num_dim == 2) then
        call check_ale_set(region(i)%num_dim, i)
    !endif
!    !if (i==1) call check_cons_condition(i, unit_number)
enddo
write (*,*) '>>> end ale_set'

if (max_dim == 2) then
    if (existence) then
        inte_set%cont_node = 0
        do
          read(unit_num_inte, *, iostat=status) keyword
          write(*,*) ' read keyword-1:   ', keyword

          if (status == -1 .or. keyword == '*end ') then
            EXIT   ! exit from do-while if meets EOF or step command
          elseif (keyword == '*inte') then
            read(unit_num_inte, *, iostat=status) fsi_inte_set%fsi_inte_point(:)
          elseif (keyword == '*cont') then
            read(unit_num_inte, *, iostat=status) inte_set%cont_node(:)
          elseif (keyword == '*abla') then  
            ablation_flag = .true.
            read(unit_num_inte, *, iostat=status) inte_set%abla_node(:)
          else
            write(*,*) 'buildmodel: build mesh model error-invalid or illegal command line ->  ', keyword
            STOP 'buildmodel'
          endif
        enddo
        write (*,*) '>>> end check_bound'
    
        ng = region(1)%num_groups
        allocate(inte_set%cont_slp(4,ng))
        inte_set%cont_slp = 0
        if (num_domains == 1) then
            if (ng == 1) then
                inte_set%cont_slp(1,ng) = fsi_inte_set%fsi_inte_point(2)
                inte_set%cont_slp(2,ng) = fsi_inte_set%fsi_inte_point(4)
            else
                stop 'STOP(mesh_interface) : num_domain = 1, over num_group 1'
            endif
        else
            call build_cont_slp(unit_num_cont, num_contacts)
        endif
        if (sum(inte_set%cont_node) == 0) inte_set%cont_node = fsi_inte_set%fsi_inte_point
    endif

    if (remesh == 0) cri_tol_length = 5.0*mesh_set(1)%ave_d(1)
    inte_set%cont_num = num_contacts
    if (num_domains>=2) fsi_inte_set%tol = mesh_set(2)%ave_d(1)

    call set_fsi_inte(num_domains, 0, 0, 0)

    call set_smoothing_nodes(1)
    if (ablation_flag) call set_smoothing_nodes(2)
elseif (max_dim == 3) then
    if (existence) then
        call set_ale_region(unit_num_inte, max_dim)
        if (remesh_type_number == 2) call ale_region_smoothing2(mesh_set(1)%nn, mesh_set(1)%node_coord, 1)
    endif
endif

end subroutine mesh_interface
