subroutine contact_control(rn, id, flag)

! This unit is originally coded by Jeeho Lee (December 5, 2012)
! Modification History:		 
! External subroutine

use coarse_domain

implicit none

integer, intent(in) :: rn, id, flag
integer :: reg_ball, reg_target, contact_type


contact_type = contact(id)%contact_type
reg_ball = contact(id)%contact_region(1)
reg_target = contact(id)%contact_region(2)

if (rn == reg_ball) then  ! update reaction-loads in target region
  call force_charge(id, reg_target, reg_ball)
elseif (rn == reg_target) then  ! update constrained disp in target region
  call disp_charge(id, reg_ball, reg_target) 
endif


end subroutine contact_control
