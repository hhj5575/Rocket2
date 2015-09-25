subroutine recover_step(Dim, remesh, num_rn, step_num, delta_t, flag)

use ale_domain
! External subroutine

implicit none

integer, intent(in) :: Dim, remesh, num_rn, flag  ! flag = 1: not remesh, 2: remesh
integer, intent(out) :: step_num
real, intent(in) :: delta_t
integer :: status
character(len=5) :: keyword

call restore_states(Dim, flag, num_rn, remesh, step_num)
if (is_coarse_problem .and. coarse_schur) call update_coarse_disp
call update_init_time(step_num, delta_t)

end subroutine recover_step
