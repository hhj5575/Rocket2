module time_domain

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (October 28, 2006)
! Modification History:		December, 2006 (by Jeeho Lee)
!                         March, 2014


implicit none
save

integer, parameter :: max_num_factor = 10
type history_type
	real, allocatable :: data(:,:)
end type

type (history_type), protected, allocatable :: history(:) 

integer, protected :: num_region_history = 1, history_array_flag = 0, timestep = 0, intgrn_method = 1 
integer, protected :: timestep_order = 0
real, protected :: time = 0.0, delt = 0.0, load_factor(max_num_factor) = 0.0, disp_factor(max_num_factor) = 0.0
real, protected :: intgrn_para(3) = 0.0, vel_para, acc_para, dis_para, vel_a_n_para, vel_a_para
real, protected :: alpha, beta, gamma	
real, protected :: ground_accel_data(3) = 0.0
logical, allocatable :: timestep_order_toggle(:)

contains
!==============================================================================

subroutine create_history(num_rn)


implicit none

integer, intent(in) :: num_rn
integer :: r


if (history_array_flag == 0) then
  allocate(history(num_rn))
  num_region_history = num_rn
  history_array_flag = 1
endif

end subroutine create_history






subroutine initialize(problem_type, rn, rn_ndof)

! called by 'femula_foro'

implicit none

integer, intent(in) :: problem_type, rn, rn_ndof

! data (column #)
!    1: displacement at time n+1 (under computing)
!    2: velocity
!    3: acceleration
!    4: displacement at time n (converged at the previous timestep)
!    5: velocity at time n (converged at the previous timestep)
!    6: acceleration at time n (converged at the previous timestep)


allocate(history(rn)%data(rn_ndof, 6))
history(rn)%data = 0.0


end subroutine initialize






subroutine update_time(load_factor_type, num_rn, problem_type, count, lf, lfsize, df, dfsize)

implicit none

integer, intent(in) :: load_factor_type, num_rn, problem_type, count, lfsize, dfsize
real, intent(in) :: lf(lfsize), df(dfsize)
integer :: i


if ((lfsize > max_num_factor).OR.(dfsize > max_num_factor)) STOP 'update_time: load set size exceeds maximum!'
if (count > 0) then
  timestep = timestep + 1
  time = time + delt

! save the disp, vel, acc data for time n
  do i = 1, num_rn  ! loop over regions
    history(i)%data(:,4) = history(i)%data(:,1)
    history(i)%data(:,5) = history(i)%data(:,2)
    history(i)%data(:,6) = history(i)%data(:,3)
  end do
elseif (count == 0) then  
  timestep = 0
endif

if (load_factor_type == 0) then ! no load/disp factor multiplication
  if (lfsize > 0) then
    load_factor(1:lfsize) = 1.0
  endif
  if (dfsize > 0) then
    disp_factor(1:dfsize) = 1.0
  endif
else
  if (lfsize > 0) then
    load_factor(1:lfsize) = load_factor(1:lfsize) + delt*lf(1:lfsize)
  endif
  if (dfsize > 0) then
    disp_factor(1:dfsize) = disp_factor(1:dfsize) + delt*df(1:dfsize)
  endif
endif

if (problem_type >= 3) then
  select case(intgrn_method)
  case(1)   ! Newmark method
    call newmark
  case default
    stop 'update_time: Invalid time integration method name!'
  end select
endif

end subroutine update_time







subroutine rewind_time(load_factor_type, num_rn, count, lf, lfsize, df, dfsize)

implicit none

integer, intent(in) :: load_factor_type, num_rn, count, lfsize, dfsize
real, intent(in) :: lf(lfsize), df(dfsize)
integer :: i


if ((lfsize > max_num_factor).OR.(dfsize > max_num_factor)) STOP 'update_time: load set size exceeds maximum!'
if (count > 0) then
  timestep = timestep - 1
  time = time - delt

! back to the disp, vel, acc data for time n
  do i = 1, num_rn  ! loop over regions
    history(i)%data(:,1) = history(i)%data(:,4)
    history(i)%data(:,2) = history(i)%data(:,5)
    history(i)%data(:,3) = history(i)%data(:,6)
  end do
elseif (count == 0) then  
  timestep = 0
endif

if (load_factor_type == 0) then ! no load/disp factor multiplication
!
else
  if (lfsize > 0) then
    load_factor(1:lfsize) = load_factor(1:lfsize) - delt*lf(1:lfsize)
  endif
  if (dfsize > 0) then
    disp_factor(1:dfsize) = disp_factor(1:dfsize) - delt*df(1:dfsize)
  endif
endif


end subroutine rewind_time



subroutine timestep_control(flag, converge_status, min_order, max_order, delta_t0)

implicit none

integer, intent(in) :: flag, converge_status, min_order, max_order
real, intent(in) :: delta_t0
integer, parameter :: converge_status_promotion = 2


if (flag == 0) then ! initialize 'timestep_order' & 'delt'
  timestep_order = 0
elseif (flag == 1) then
  if (converge_status > 0) then
    if (converge_status > converge_status_promotion) then
      timestep_order = MIN(max_order, timestep_order + 1)
      write(*,*) '***** step order is changed to:', timestep_order        
    endif
  elseif (timestep_order < min_order) then
    STOP 'timestep_control: Fail to converge!'    !---------------------------------------------->>>
  else
    timestep_order = timestep_order - 1
    write(*,*) '***** step order is changed to:', timestep_order        
  endif
endif
delt = delta_t0*(2.0**timestep_order)
if (flag == 2) delt = delta_t0

end subroutine timestep_control



subroutine timestep_control_bak(flag, converge_status, min_order, delta_t0)

implicit none

integer, intent(in) :: flag, converge_status, min_order
real, intent(in) :: delta_t0
integer, parameter :: converge_status_promotion = 2 ! decide when the timestep enlarged to accelerate the analysis execution


if (flag == 0) then ! initialize 'timestep_order' & 'delt'
  timestep_order = 0
  if (allocated(timestep_order_toggle)) deallocate(timestep_order_toggle)
  allocate(timestep_order_toggle(0:ABS(min_order)))
  timestep_order_toggle = .FALSE.
else
  if (converge_status > 0) then
    if (timestep_order < 0) then
      if ((converge_status > converge_status_promotion).AND.timestep_order_toggle(ABS(timestep_order))) then
        timestep_order = timestep_order + 1
        timestep_order_toggle(ABS(timestep_order)) = .NOT.timestep_order_toggle(ABS(timestep_order))
        write(*,*) '***** step order is changed to:', timestep_order        
      endif
    endif
  elseif (timestep_order == min_order) then
    STOP 'timestep_control: Fail to converge!'    !---------------------------------------------->>>
  else
    timestep_order = timestep_order - 1
    timestep_order_toggle(ABS(timestep_order)) = .TRUE.
    write(*,*) '***** step order is changed to:', timestep_order        
  endif
endif
      
delt = delta_t0*(2.0**timestep_order)

if (timestep_order < 0) then
  timestep_order_toggle(ABS(timestep_order)) = .NOT.timestep_order_toggle(ABS(timestep_order))
  if (timestep_order_toggle(ABS(timestep_order))) timestep_order_toggle(ABS(timestep_order)-1) = .NOT.timestep_order_toggle(ABS(timestep_order)-1)
endif

end subroutine timestep_control_bak





subroutine update_ground_accel(flags, file_numbers)

implicit none

logical, intent(in) :: flags(3)
integer, intent(in) :: file_numbers(3)
integer :: i


do i = 1, 3
  if (flags(i)) then
    call time_table(file_numbers(i), time, ground_accel_data(i))
  endif
end do

end subroutine update_ground_accel







subroutine time_table(fnum, t, value)

! time table reader: forward & backward time step reading for fast searching in a pre-opened file
! two column data table (1st field: time, 2nd: data value)
! linear interpolation

implicit none

integer, intent(in) :: fnum
real, intent(in) :: t
real, intent(out) :: value
		
integer :: status
real :: ttble, oldtt, ug, oldug
real :: oldtime, factor
 
           
do
	read(fnum, *, iostat=status) ttble, ug
	if (status ==	-1) then
		write(*,*) 'time_table: warning - current time exceeds upper time table range! 0.0 is assigned.'
		value = 0.0
		EXIT
	elseif (t == ttble) then
		value = ug
		EXIT
	elseif (t < ttble) then   ! forward reading
		backspace(fnum)
		do
			backspace(fnum)

			read(fnum, *, iostat=status) oldtt, oldug      ! backward reading

			if (t >= oldtt) then
				EXIT
			endif
			backspace(fnum)
			ttble = oldtt
			ug = oldug
		end do	
		
      factor = (t - ttble)/(ttble - oldtt)
      value = factor*(ug - oldug) + ug
			EXIT
	endif
end do

end subroutine time_table






subroutine time_integration(flag, rn, problem_type, u_physical)

!   flag = 0: initial iteration step
!        = 1: 2nd or further more iteration steps

implicit none

integer, intent(in) :: flag, problem_type, rn
real, intent(inout) :: u_physical(:)    ! u_physical: dof numbers are arranged for physical domain

integer :: ndim
real, allocatable :: acc(:), vel(:)
real :: gt1, gt2


if (problem_type >= 3) then   ! Dynamic problems

  if (timestep == 0) then
    if (flag == 1) then
      history(rn)%data(:,3) = u_physical   ! initial acceleration computed based on initial force
      u_physical = 0.0
    endif		
  else	
    ndim = size(u_physical)
    allocate(acc(ndim))
    allocate(vel(ndim))
    acc = history(rn)%data(:,6)
    vel = history(rn)%data(:,5)

! Displacement vector update			
    history(rn)%data(:,1) = u_physical
  
! Acceleration vector update
		if (flag == 0) then   ! when initial iteration
			history(rn)%data(:,3) = vel_para*vel + acc_para*acc
		else
			history(rn)%data(:,3) = vel_para*vel + acc_para*acc + dis_para*(u_physical - history(rn)%data(:,4))
		endif

! Velocity vector update
    history(rn)%data(:,2) = vel + vel_a_n_para*acc + vel_a_para*history(rn)%data(:,3)
		
		deallocate(acc, vel)
  endif
  
else     ! Static or Quasi-static problems
  history(rn)%data(:,1) = u_physical
endif  


end subroutine time_integration





subroutine intgrn_para_setup(method, v1, v2, v3)

implicit none

character(len=5), intent(in) :: method
real, intent(in) :: v1, v2, v3


write(*,*) 'Time integration method: ', method
select case(method)
case('newma')
  intgrn_method = 1
  beta = v1   ! v1: beta
  gamma = v2  ! v2: gamma
  alpha = 0.0
case default
  stop 'intgrn_para_setup: Invalid time integration method name!'
end select

end subroutine





subroutine newmark

! timestep = 0: initial step
!         >= 1: normal step

implicit none

! write(*,*) 'beta, delt', beta, delt

if (timestep == 0) then
  intgrn_para(1) = 0.0   ! coefficient for tangent stiffness matrix
  intgrn_para(2) = 0.0   ! coefficient for tangent damping matrix
  intgrn_para(3) = 1.0   ! coefficient for mass matrix
else
  intgrn_para(1) = 1.0   ! coefficient for tangent stiffness matrix
  intgrn_para(2) = gamma/(beta*delt)   ! coefficient for tangent damping matrix
  intgrn_para(3) = 1.0/(beta*delt*delt)   ! coefficient for mass matrix
endif

vel_para = -1.0/(beta*delt)   ! for acceleration computation
acc_para = 1.0 - 0.5/beta     ! for acceleration computation
dis_para =  intgrn_para(3)    ! for acceleration computation

vel_a_n_para = (1.0-gamma)*delt   ! for velocity computation
vel_a_para = gamma*delt           ! for velocity computation

end subroutine newmark


!==============================================================================

subroutine update_init_time(step_number, dt)

implicit none

integer, intent(in) :: step_number
real, intent(in) :: dt

time = step_number * dt

end subroutine update_init_time


subroutine update_motion(rn, dim,  position , val)

implicit none

integer,intent(in) :: rn, dim
integer,intent(in) :: position(dim)
real, intent(in) :: val(dim,3)

history(rn)%data(position,1) = val(:,1)
history(rn)%data(position,2) = val(:,2)
history(rn)%data(position,3) = val(:,3)

end subroutine update_motion


subroutine free_history

implicit none 

deallocate( history )

history_array_flag = 0
num_region_history = 1; intgrn_method = 1
intgrn_para(3) = 0.0
ground_accel_data(3) = 0.0

end subroutine 

!==============================================================================
end module time_domain
