subroutine fem_gxl_initialize(Dim, action, solid_flag, remesh)

use file_info_house

implicit none

integer, intent(in) :: Dim, action, solid_flag
integer, intent(inout) :: remesh

logical :: existence


select case(action)
case (0)  ! open major files
  output_node_file = './output/solid/xfemout_node'
  output_element_file = './output/solid/xfemout_elem'

  open(unit=unit_number, file=filename, status='OLD', action='READ')
  open(unit=out1_unit_number, file=output_node_file, status='REPLACE', action='READWRITE')
  open(unit=out2_unit_number, file=output_element_file, status='REPLACE', action='READWRITE')
  main_input_file = filename
  write(*,*) '>>FEM_GXL: open input file: Sucessfully performed!'
  
case (1)  ! initialize & build
  if (solid_flag == 0) then
    call add_char(filename, node_file, 1)
    call add_char(filename, element_file, 2)
    call add_char(filename, contact_file, 3)  
    call add_char(filename, inte_file, 4)
      
  elseif (solid_flag == 1) then
    remesh = 0
    call add_char(filename, node_file, 1)
    call add_char(filename, element_file, 2)
    call add_char(filename, contact_file, 3)
    call add_char(filename, inte_file, 4)

  elseif (solid_flag == 5) then
    !rewind(unit_number)
    close(unit_number)
    close(unit_num_node)
    close(unit_num_elem)
    close(unit_num_cont)
    close(unit_num_inte)

    remesh_file = './remesh/remesh000.main'
    node_file = './remesh/remesh000.node' 
    element_file = './remesh/remesh000.elem'
    inte_file = './remesh/remesh000.inte'
    contact_file = './remesh/remesh000.cont'
    
    write(remesh_file(16:18), '(I3.3)') remesh
    write(node_file(16:18), '(I3.3)') remesh
    write(element_file(16:18), '(I3.3)') remesh
    write(inte_file(16:18), '(I3.3)') remesh
    write(contact_file(16:18), '(I3.3)') remesh
    
    open(unit=unit_number, file=remesh_file, status='OLD', action='READ')
    
  else
    STOP 'fem_gxl_initialize: invalid solid_flag!'
  endif


  inquire (file=node_file, EXIST=existence)
  if (existence) then
    open(unit=unit_num_node, file=node_file, status='OLD', action='READ')
    write(*,*) '>>FEM_GXL: open node data file: Sucessfully performed!'
  else
    write(*,*) '>>FEM_GXL: open node data file: File does NOT exist!!!'
    STOP
  endif

  inquire (file=element_file, EXIST=existence)
  if (existence) then
    open(unit=unit_num_elem, file=element_file, status='OLD', action='READ')
    write(*,*) '>>FEM_GXL: open element data file: Sucessfully performed!'
  else
    write(*,*) '>>FEM_GXL: open element data file: File does NOT exist!!!'
    STOP
  endif

  inquire (file=contact_file, EXIST=existence)
  if (existence) then
    open(unit=unit_num_cont, file=contact_file, status='OLD', action='READ')
    write(*,*) '>>FEM_GXL: open contact info file: Sucessfully performed!'
  else
    write(*,*) '>>FEM_GXL: open contact info file: File does NOT exist!!!'
  endif

  inquire (file=inte_file, EXIST=existence)
  if (existence) then
    open(unit=unit_num_inte, file=inte_file, status='OLD', action='READ')
    write(*,*) '>>FEM_GXL: open interface data file: Sucessfully performed!'
  else
    write(*,*) '>>FEM_GXL: open interface data file: File does NOT exist!!!'
  endif
  
  gr_unit_numbers(1) = gr1_unit_number
  gr_unit_numbers(2) = gr2_unit_number
  gr_unit_numbers(3) = gr3_unit_number

  !if (Dim == 2) then
    call build_model(remesh)
  !elseif (Dim == 3) then
  !  call build_model_3D(remesh)
  !endif
  write(*,*) '>> FEM_GXL: build_model: Sucessfully performed!'

case default

	STOP 'fem_gxl_initialize: invalid action number!'
end select


end subroutine fem_gxl_initialize
