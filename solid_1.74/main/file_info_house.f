module file_info_house

implicit none
save

integer, parameter :: unit_number = 101, unit_num_node = 102, unit_num_elem = 103, unit_num_cont = 104
integer, parameter :: unit_num_factor = 105
integer, parameter :: unit_num_inte=106, unit_num_remesh = 107 
integer, parameter :: out1_unit_number = 110, out2_unit_number = 111, out3_unit_number = 112
integer, parameter :: gr1_unit_number = 121, gr2_unit_number = 122, gr3_unit_number = 123
character(len=35) :: filename
character(len=35) :: main_input_file, node_file, element_file, contact_file
character(len=35) :: inte_file, remesh_file
character(len=35) :: output_node_file, output_element_file
integer :: output_node, output_region1, output_region2, output_interval = 1, output_interval2 = 10
integer :: debug_out_cnr = 0, debug_region, debug_elem, debug_timestep
integer :: str_out_cnr = 0
character(len=6) :: output_cnt
integer :: gr_unit_numbers(3)
logical :: gr_motion_flags(3), element_output_flag = .FALSE.

integer, parameter :: unit_num_remesh_energy = 131
character(len=40) :: remesh_energy_file
logical :: remesh_energy_flag = .FALSE.

end module file_info_house