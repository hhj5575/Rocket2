module system_info_house

implicit none
save

integer, parameter :: system_max_itr = 20
real, parameter :: system_tolerance = 1.0e-12

logical :: fsi_flag = .TRUE., solver_flag = .TRUE., ablation_flag = .FALSE.
integer :: max_itr
real :: tolerance
character(len=6) :: solver_kind = 'linear'
logical :: auto_time_stepping = .FALSE.

real, allocatable :: loadf(:), dispf(:)

end module system_info_house
