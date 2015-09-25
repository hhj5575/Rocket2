module materials

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (October 28, 2006)
! Modification History:		December 05, 2006 (by Jeeho Lee)
!                         March 2014
!

!   state1(:,:) - previous step data
!   state2(:,:) - CURRENT step data
!   state3(:,:) - back-up for previous step data for rewinding

implicit none
save

type material_type
  logical :: flag
  integer :: num_row_state = 0, num_col_state = 0
  integer :: size_para, size_d, model_num
  real, allocatable :: para(:)
	real, allocatable :: d(:), state1(:,:), state2(:,:), state3(:,:)
end type material_type

type(material_type), protected, allocatable :: material(:)
integer :: number_of_materials

contains
!==============================================================================


subroutine create_material(number)

implicit none

integer, intent(in) :: number


allocate(material(number)) 
material(number)%flag = .FALSE.
number_of_materials = number

end subroutine create_material





subroutine set_d_array(switch,mat_num,size_d,d_array)

implicit none

integer, intent(in) :: switch,mat_num
integer, intent(inout) :: size_d
real, intent(inout) :: d_array(size_d)

if (switch == 1) then
	material(mat_num)%size_d = size_d
	allocate(material(mat_num)%d(size_d))

	material(mat_num)%d = d_array
	material(mat_num)%flag = .TRUE.
else
  d_array = 0.0
	d_array = material(mat_num)%d  
endif

end subroutine set_d_array






subroutine material_counter(switch, mat_num, count)

! read/record to count the number of quad pts per each material across domains (regions)
! switch = 1: read,  = 2: record

implicit none

integer, intent(in) :: switch, mat_num
integer, intent(inout) :: count


select case (switch)
case(1)
  count = material(mat_num)%num_col_state
case(2)
  material(mat_num)%num_col_state = count
end select

end subroutine material_counter




subroutine number_material_state_var(mat_num, num_state_var)

implicit none

integer, intent(in) :: mat_num, num_state_var


material(mat_num)%num_row_state = num_state_var

end subroutine number_material_state_var





subroutine allocate_states(tm_flag, mat_num)

! tm_flag: time machine flag (.TRUE. if time machine is on)

implicit none

logical, intent(in) :: tm_flag
integer, intent(in) :: mat_num
integer :: m, n

m = material(mat_num)%num_row_state


if (m >= 1) then
  n = material(mat_num)%num_col_state
  allocate(material(mat_num)%state1(m,n))
  allocate(material(mat_num)%state2(m,n))
  material(mat_num)%state1 = 0.0
  material(mat_num)%state2 = 0.0
  if (tm_flag) then
    allocate(material(mat_num)%state3(m,n))
    material(mat_num)%state3 = 0.0
  endif
endif

write(*,*) '******* allocate_states: m, n =', m, n

end subroutine allocate_states





subroutine material_history(flag, mat_num, point_num, mat_state, ndim)

! flag = 0: download from the current timestep result (state2)
!      = 1: download from the previous timestep result (state1)
!      = 2: upload to the new timestep array (state2)
!      = 3: update (state2 -> state1)


implicit none

integer, intent(in) :: flag, mat_num, point_num, ndim
real :: mat_state(ndim)

if (flag == 0) then
  mat_state = material(mat_num)%state2(1:ndim,point_num)
elseif (flag == 1) then
  mat_state = material(mat_num)%state1(1:ndim,point_num)
elseif (flag == 2) then
  material(mat_num)%state2(1:ndim,point_num) = mat_state
elseif (flag == 3) then
  material(mat_num)%state1(1:ndim,point_num) =  material(mat_num)%state2(1:ndim,point_num)
else
  STOP 'material_history: flag number is not valid!'
endif

end subroutine material_history





subroutine material_history_update(tm_flag, update_flag)

! Update material state variables after the iteration is converged at a given time step
! tm_flag: time machine flag (.TRUE. if time machine is on)
!         = -1: rewind (valid only when tm_flag is 'ON')
! update_flag = 1: normal update
!             = 2: backward rewind

!   state1(:,:) - previous step data
!   state2(:,:) - CURRENT step data
!   state3(:,:) - back-up for previous step data for rewinding

implicit none

logical, intent(in) :: tm_flag
integer, intent(in) :: update_flag

integer :: i, m


if (tm_flag) then
  if (update_flag > 0) then
    do i = 1, number_of_materials
      m = material(i)%num_row_state
      if (m >= 1) then
        material(i)%state1 =  material(i)%state2
      endif
    enddo
  else  ! rewind mode
    do i = 1, number_of_materials
      m = material(i)%num_row_state
      if (m >= 1) then
        material(i)%state1 =  material(i)%state3
      endif
    enddo  
  endif
else
  if (update_flag <= 0) STOP 'material_history_update: rewind is impossible if time machine is off!'
  do i = 1, number_of_materials
    m = material(i)%num_row_state
    if (m >= 1) then
      material(i)%state1 =  material(i)%state2
    endif
  enddo
endif

end subroutine material_history_update





subroutine set_material(n, mat_model_num, size, values)

implicit none

integer, intent(in) :: n, mat_model_num, size
real, intent(in) :: values(size)
integer :: size2

size2 = size
if (size2 < 3) size2 = 3
allocate(material(n)%para(size2))

material(n)%model_num = mat_model_num
material(n)%size_para = size
material(n)%para = values

end subroutine set_material
   




subroutine material_data(mat_num, para, ndim)

implicit none

integer, intent(in) :: mat_num, ndim
real, intent(out) :: para(ndim)


para(1:ndim) = material(mat_num)%para(1:ndim)

end subroutine material_data





subroutine print_material(mat_num,sign)

implicit none

integer, intent(in) :: mat_num,sign


write(*,*) '******* print_material', sign
write(*,*) '*****  material d size: ', material(mat_num)%para, material(mat_num)%size_d

end subroutine print_material


!==============================================================================

subroutine get_variable( mat_num , point_num , size, value )

implicit none

integer, intent(in) :: mat_num, point_num , size

real, intent(out)  :: value(size)


if ( size == material(mat_num)%num_row_state ) then
    value = material(mat_num)%state1(:,point_num)
else
    STOP "SIZE is wrong"
endif

end subroutine


subroutine update_states(mat_num, csize, rsize, val )

implicit none

integer, intent(in) :: mat_num, csize, rsize
real, intent(in) :: val( csize ,rsize )

if ( ( material(mat_num)%num_row_state .EQ. rsize ) .OR. (material(mat_num)%num_col_state .EQ. csize ) ) then
    STOP "State variable update FAILD!! "
endif
material(mat_num)%state1 = val

end subroutine update_states

subroutine update_states2(mat_num, csize, rsize, val )

implicit none

integer, intent(in) :: mat_num, csize, rsize
real, intent(in) :: val( csize ,rsize )

if ( ( material(mat_num)%num_row_state .EQ. rsize ) .OR. (material(mat_num)%num_col_state .EQ. csize ) ) then
    STOP "State variable update FAILD!! "
endif
material(mat_num)%state2 = val

end subroutine update_states2

subroutine update_states21(mat_num)

implicit none

integer,  intent(in) :: mat_num
integer :: m

m = material(mat_num)%num_row_state
if (m >= 1) then
    material(mat_num)%state1 =  material(mat_num)%state2
endif

end subroutine update_states21
  
subroutine free_material

implicit none 

deallocate( material )
number_of_materials = 0 

end subroutine 

!===========================================================================

end module materials


