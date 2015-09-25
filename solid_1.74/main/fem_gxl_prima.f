subroutine fem_gxl_prima(Dim, step_number, converge_status, error, last_step)

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (October 28, 2006)
! Modification History: December 05, 2006
!                       October 31, 2009
!                       May, 2013
!                       March, 2014
!
! One timestep computation

use file_info_house
use system_info_house
use system_domain
use materials
use coarse_domain
! use ale_domain

implicit none

integer, intent(in) :: Dim, step_number
integer, intent(inout) :: converge_status
real, intent(inout) :: error
logical, intent(in) :: last_step

integer :: num_rn
integer :: rn
integer :: contact_id, contact_flag, ld_factor_type
integer :: converge_status0, converge_status1
integer :: debug_target(2)
real :: dummy


!if (Dim == 2) then
    !write(*,'(/A,I6,5X,A)') ' >> FEM_GXL: computation start for time step #: ', step_number, solver_kind

    num_rn = num_of_regions ! global variable from physical_domain module

    if (fsi_flag) then
      ld_factor_type = 0
    else
      ld_factor_type = 1
    endif

    call update_time(ld_factor_type, num_rn, problem_type, step_number, loadf, num_load_sets, dispf, num_disp_sets) ! time_domain & system_domain, and system_control modules

    call update_ground_accel(gr_motion_flags, gr_unit_numbers)
    call backup_disp_bc(1) 
    
    if (step_number == debug_timestep) then  ! 'debug_timestep', 'debug_region', 'debug_elem' from module 'file_info_house'
      debug_target(1) = debug_region
      debug_target(2) = debug_elem 
    else
      debug_target(1) = 0
      debug_target(2) = 0
    endif
    
    if (solver_flag) then
        if (is_coarse_problem) then
          if (coarse_schur) then
            call coarse_newton_solver(num_rn, max_itr, tolerance, error, converge_status)
          else
            call global_newton_solver(num_rn, max_itr, tolerance, error, converge_status)
          endif
        else
          converge_status0 = converge_status
          converge_status1 = 100
          do rn = 1, num_rn
            converge_status = converge_status0    ! for rn > 1 case
            select case(solver_kind)
            case('linear')
              call linear_solver(rn, error, debug_target)
            case('newton')
              call newton_solver(rn, max_itr, tolerance, error, converge_status, debug_target)
              if (converge_status < 0) then   ! Exit due to failed/poor convergency
                converge_status1 = converge_status
                EXIT    !---------------------------------------------------------------------------->>>>>
              endif
            case default
              STOP 'fem_gxl_prima: undefined solver type'
            end select
            converge_status1 = MIN(converge_status, converge_status1) ! minimum over rn
          end do
  
          converge_status = converge_status1
        endif
        !write(*,*) 'converge_status', converge_status
    else 
        converge_status = 1
        error = 0.0
    endif

    if (converge_status > 0) then
    ! output control variables are stored in 'file_info_house' module
      do rn = 1, num_rn
        call node_output(out1_unit_number, output_cnt, rn, output_region1, output_node)
        if (element_output_flag .and. last_step) call element_output(out2_unit_number, str_out_cnr, rn, output_region2, step_number, output_interval)        
      end do
        
      !write(*,'(A,I6,A,ES12.5,A,ES12.5/)') ' >>FEM_GXL: done @ step # =', step_number, '  *time=', time, '  *error norm=', error

      if (element_output_flag.AND. output_region2 == 0 .and. last_step) call output_interpolation(Dim, out2_unit_number, num_rn, 0, output_region2, time, step_number, output_interval, 2)
      
      call material_history_update(.FALSE., 1)  ! update all material state variables ('materials' module)
    elseif (converge_status < 0) then
      call rewind_time(ld_factor_type, num_rn, step_number, loadf, num_load_sets, dispf, num_disp_sets) 
      call backup_disp_bc(-1)
    endif
!elseif (Dim == 3) then
!    call update_3D_time
!    call output_interpolation(Dim, out2_unit_number, num_rn, 0, output_region2, time, step_number, output_interval, 2)
!endif

end subroutine fem_gxl_prima
