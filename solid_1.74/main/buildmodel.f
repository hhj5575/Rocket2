subroutine build_model(remesh)

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (October 28, 2006)
! Modification History: December 05, 2006 (by Jeeho Lee)
!                       October 25, 2009
!                       August 20, 2010

use file_info_house
use physical_domain
use materials
! use ale_domain
use output_domain

implicit none

integer, intent(inout) :: remesh

integer :: unit_num, num_domains, num_materials, num_load_bd_sets, num_disp_bd_sets, num_contacts
integer :: status, var(10), max_node_dof, max_dim, max_nodes, max_elements, nodes_elem, num_quad_pts, num_groups=1
integer :: node_addr, elem_addr
integer :: i, j, k, l, n, ig, count, number, place, old, current_domain, node_count, elem_count, element_type
integer :: spring_count
integer :: mat_model_num, material_num, rn, group_number=1
integer :: node_input_flag, input_flag, material_count, load_type_number, set_number=1
integer :: max_max_dim = 1, stress_dim = 4

real :: area, density, consist_mass_ratio, value , time_a, time_b, rtc, damping_factor(2), ref_temperature = 0.0
real :: points(3,4)

character(len=5) :: keyword
character(len=20) :: load_factor_file, gr_file_name
character(len=60) :: title

logical :: flag_bc, flag_load

integer, allocatable :: connect(:), boundary(:,:)
integer, allocatable :: load_point(:,:)
real, allocatable :: field(:), load_value(:), time(:), temp_load(:)

!===============================================================================

!if (remesh == 0) then
  unit_num = unit_number
!else
!  unit_num = unit_num_remesh
!endif

time_a = RTC()
gr_motion_flags = .FALSE.

do
	read(unit_num, *, iostat=status) title
	if (title(1:1) == '!') then
		write(*,'(/A,A60/)') ' >> FEMULA: start reading input file: ', title(1:60)
	elseif (status == -1) then
		STOP 'build_model: End of File error from input file!'
	else
		backspace(unit_num)
		read(unit_num, *) keyword, num_domains, num_materials, num_load_bd_sets, num_disp_bd_sets, num_contacts
		write(*,*) ' read keyword:   ', keyword, num_domains, num_materials, num_load_sets, num_disp_sets, num_contacts
		EXIT
	endif
end do

if ((num_load_bd_sets + num_disp_bd_sets) < 1) STOP 'build_model: (num_load_sets + num_disp_sets) must not less than 1!'

! create physical domains
call create_domain(num_domains)
write(*,*) 'FEMULA: completed: domain space creation'

call create_boundary_set(num_load_bd_sets, num_disp_bd_sets, num_contacts)
write(*,*) 'FEMULA: completed: boundary (load & disp) set creation'

call create_material(num_materials)
write(*,*) 'FEMULA: completed: material space creation'

var = 0

do    ! do over all domain regions =================================================================
	read(unit_num, *, iostat=status) keyword
	write(*,*) ' read keyword-1:   ', keyword

	if (status == -1 .or. keyword == '*end ') then
		EXIT   ! exit from do-while if meets EOF or contact or step command
	elseif (keyword == '*doma') then
		backspace(unit_num)
		read(unit_num, *, iostat=status) keyword, current_domain
		if (current_domain < 1 .OR. current_domain > num_domains) then
      write(*,*) 'build_model: invalid domain number! current_domain, num_domains:', current_domain, num_domains
			STOP
    endif
    read(unit_num,*) max_node_dof, max_dim, max_nodes, max_elements, nodes_elem, num_quad_pts, num_groups
    rewind(unit=unit_num_node)
    rewind(unit=unit_num_elem)

    if (max_dim < 1) then
      STOP 'build_model: max_dim is less than 1!'
    else
      max_max_dim = MAX(max_max_dim, max_dim)
    endif

    call set_domain(current_domain,max_node_dof,max_dim,max_nodes,max_elements,nodes_elem,num_quad_pts,num_groups)
    write(*,*) '      set_domain completed ! region number:', current_domain
    if (num_groups > 1) then
      do ig = 1, num_groups
        read(unit_num,*) group_number, node_addr, elem_addr
        if (group_number < 0 .OR. group_number > num_groups) STOP 'buildmodel: input group_number is invalid!'
        call write_subregion_group(0, current_domain, group_number, node_addr, elem_addr)
      enddo
    endif        
    write(*,*) '      set_domain completed ! number of subregion groups:', num_groups

    read(unit_num, *, iostat=status) keyword  ! check if there are spring elements
    backspace(unit_num)
    if (keyword == '*spri') then
      read(unit_num, *, iostat=status) keyword, max_elements  ! max number of spring elements
      write(*,*) 'build_model: max number of spring elements =', max_elements
      call set_spring_elements(current_domain, max_elements)
    endif
     
    call build_nodes(unit_num_node, current_domain)
    write(*,*) '      build_nodes completed !'
    
    elem_count = 0
    spring_count = 0
    do
      read(unit_num_elem, *, iostat=status) keyword
!      write(*,*) 'read keyword in element input file: ', keyword
      
      if (status == -1 .or. keyword == '*end ') then
				EXIT   ! exit from do-while loop
        
      elseif (keyword == '*elem') then
        backspace(unit_num_elem)
        read(unit_num_elem,*) keyword, rn, element_type, material_num, area, density, consist_mass_ratio, damping_factor(1),damping_factor(2)
        if (num_domains < rn) then
          write(*,*) '* build_model: Invalid domain region number in element input file!', rn
          EXIT    ! skip reading element file due to invalid domain number
        endif
        if (rn == current_domain) then
          write(*,*) '-----------------------------------------------------------------'
          write(*,*) ' read element file for Domain Region #:', rn
          write(*,*) '   Element Type #:', element_type
          write(*,*) '   Material Type #:', material_num
          write(*,*) '   Area:    ', area 
          write(*,*) '   Density: ', density
          write(*,*) '   Consistent Mass Ratio:', consist_mass_ratio
          write(*,'(A,E12.5,A,E12.5)') '    Rayleigh damping (M & K proportional): alpha =', &
                    damping_factor(1), '  beta =', damping_factor(2)
          write(*,*) '-----------------------------------------------------------------'
        
          if (material_num > num_materials) STOP 'build_model: invalid material number!'
          if ((consist_mass_ratio < 0.0).or.(consist_mass_ratio > 1.0)) STOP 'consist_mass_ratio input data error!'
    
          call material_counter(1, material_num, material_count)    ! from module 'materials'
          call build_elements(unit_num_elem, current_domain, element_type, material_num,area, &
                             density, consist_mass_ratio, damping_factor, material_count,elem_count,.TRUE.)
          call material_counter(2, material_num, material_count)
        endif
        
      elseif (keyword == '*spri') then  
        backspace(unit_num_elem)
        read(unit_num_elem,*) keyword, rn, element_type, material_num, area

        if (num_domains < rn) then
          write(*,*) '* build_model: Invalid domain region number in element input file!', rn
          EXIT    ! skip reading element file due to invalid domain number
        endif
        if (rn == current_domain) then
          write(*,*) '-----------------------------------------------------------------'
          write(*,*) ' read element file for Domain Region #:', rn
          write(*,*) '   Spring Element Type #:', element_type
          write(*,*) '   Material Type #:', material_num
          write(*,*) '   Area:    ', area 
          write(*,*) '-----------------------------------------------------------------'
        
          if (material_num > num_materials) STOP 'build_model: invalid material number!'
    
          call material_counter(1, material_num, material_count)    ! from module 'materials'
          
          write(*,*) '****** spring material_count: in =', material_count
          call build_elements(unit_num_elem, current_domain, element_type, material_num,area, &
                             density, consist_mass_ratio, damping_factor, material_count,spring_count,.FALSE.)
          write(*,*) '****** spring material_count: out =', material_count
          call material_counter(2, material_num, material_count)
        endif                
      endif
      
    end do
    
    if (region(current_domain)%num_elements > 0) then
      if (region(current_domain)%num_dim == 2) then
        call write_skin(unit_num_elem, current_domain)
        write(*,*) '>> FEMULA: write_skin completed !'
      elseif (region(current_domain)%num_dim == 3) then
        call write_surface(current_domain)
        write(*,*) '>> FEMULA: write_surface completed !'
      endif
    endif


  ! Read boundary conditions

    write(*,'(/A)') '>> FEMULA: read_boundary starts !'

    call read_boundary(unit_num, current_domain) 
    write(*,'(A/)') '>> FEMULA: read_boundary completed !'
    
    call assign_global_dof(current_domain)
    write(*,*) '>> FEMULA: completed: DOF assignment !'


    input_flag = 2	
		
 
 		write(*,*) '-------- completed: set domain for # ', current_domain
    write(*,*) '------------------------------------------------------------------'
		
  ! Initialize counters and flags in a region
    flag_bc = .FALSE.
    flag_load = .FALSE.

    write(*,'(/A)') '>> FEMULA: read_loads starts !'

    call read_loads(unit_num, current_domain, flag_bc, flag_load, input_flag, gr_unit_numbers, gr_motion_flags)
    write(*,'(A/)') '>> FEMULA: read_loads completed !'

    do i = 1, num_load_sets  ! global variable in physical_domain
      if (load_set(i)%surface_pressure_existence) call surface_load(0, i, current_domain)
    end do  ! i
    write(*,'(A/)') '>> FEMULA: surface_load completed !'


!   =============================================================

	elseif (keyword == '*mate') then
		backspace(unit_num)
		read(unit_num, *, iostat=status) keyword, material_num, mat_model_num, number
    allocate(field(number))
		n = MOD(number, 10)
		count = 1
		do i = 1, number/10                                 ! each line: 10 input data
			read(unit_num,*) (field(j), j=count, count+9)   ! limitation: field(50) -> 5 input lines
			count = count + 10
		end do
		if (n > 0) then
		  read(unit_num,*) (field(i), i=count, count+n-1)  ! read last line if less than 10 data exist there
		endif
		count = count + n - 1
		
		write(*,'(A,/,(10ES10.3))') ' material data: ', (field(j), j=1, count)

		if (number /= count) STOP 'buildmodel: number of model input data inconsistence!'
		
		call set_material(material_num, mat_model_num, number, field)

    deallocate(field)

	elseif (keyword == '*fact') then
		backspace(unit_num)
		read(unit_num, *, iostat=status) keyword, load_factor_file
		call loadfactor(unit_num_factor,load_factor_file)
		
	elseif (keyword == '') then
		EXIT
	else
    write(*,*) 'buildmodel: build mesh model error-invalid or illegal command line ->  ', keyword
		STOP 'buildmodel'
	endif
end do    ! do over all domain regions =============================================================


write (*,*) '>> FEM_GXL: mesh_interface'
call mesh_interface(num_domains, num_contacts, num_materials, remesh)

if (num_contacts > 0) then
  call build_interface(1, unit_num_cont, num_contacts)    ! external subroutine
  write(*,*) '>> FEM_GXL: contact info build completed! #of contacts =', num_contacts
endif
    
write(*,*) '------------------------------------------------------------'
write(*,*) 'FEMULA: completed: read all data'

if (max_max_dim < 3) then
  stress_dim = 4
else
  stress_dim = 6
endif

write(*,*) '  ****** maximum of max_dim =', max_max_dim
write(*,*) '  ******         stress_dim =', stress_dim

do i = 1, num_materials
  mat_model_num = material(i)%model_num
  write(*,*) 'Mat model # :', mat_model_num
  call set_material_db(num_materials,i,mat_model_num)   ! set material db in physical domain
  if (mat_model_num == 21) then  
    call material_initialize_chsv1_2d(0,i,mat_model_num)
  else
    call material_initialize(0,i,mat_model_num,stress_dim)  ! in materials domain
  endif
  call allocate_states(.FALSE.,i)                     ! in materials domain
end do

write(*,*) '  ** about to : set_mat_mesh'
call set_mat_mesh(unit_num_elem, unit_num_node, num_domains, num_materials)
	
write(*,*) '------------------------------------------------------------'
write(*,*) 'FEMULA: completed: all buildmodel procedures'	
	
! write(*,*) 'Load Factor DB Table: ', load_factor_db(1:2,1:2)

end subroutine build_model