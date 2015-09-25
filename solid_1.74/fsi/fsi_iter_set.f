subroutine fsi_iter_set(n_remesh, max_interface_pt, ni, pload, inte_coord, inc_pload, iter, delta_t, temp_delta_t)
    
use fsi_domain
implicit none

integer, intent(in) :: n_remesh, max_interface_pt, ni
integer, intent(inout) :: iter
real, intent(in) :: delta_t, inte_coord(2,max_interface_pt), pload(max_interface_pt)
real, intent(inout) :: temp_delta_t, inc_pload(max_interface_pt)

integer :: i
real :: max_val, pre_max_pload, cri_pload

!cri_pload = 101325*10.0  ! 10 atmosphere
!inc_pload = 0.0
!iter = 0
!if (pre_interface_nn == 0) then
!    iter = 1
!    inc_pload(1:ni) = pload(1:ni)
!else
!    if (n_remesh /= 0) then
!        call pload_interpolation(max_interface_pt, pre_interface_nn, pre_pload, pre_inte_coord, ni, inte_coord)
!    endif
!
!    max_val = 0.0
!    pre_max_pload = maxval(pre_pload(1:ni))
!    do i = 1, 2*(ni-1)
!        inc_pload(i) = pload(i) - pre_pload(i)
!        if (abs(inc_pload(i)) > max_val) max_val = abs(inc_pload(i))
!    enddo
!    write (*,*) 'inc max pload:', max_val
!    if (maxval(pload(1:ni)) > cri_pload .and. max_val > cri_pload .and. pre_max_pload > cri_pload) then
!        max_val = max_val/pre_max_pload
!        iter = int(max_val/0.1)
!    else
!        iter = 1
!    endif
!    if (iter == 0) iter = 1
!    if (iter <= 2) then
!        if (n_remesh /= 0) iter = 3
!    endif
!endif
!
!write (*,*)
!write (*,*) 'iter:', iter
!inc_pload(1:ni) = inc_pload(1:2*(ni-1))/float(iter)
!temp_delta_t = delta_t/float(iter)


end subroutine fsi_iter_set


subroutine difference_pload(Dim, rn, nn, pload, inc_pload)
    
use fsi_domain
implicit none

integer, intent(in) :: Dim, rn, nn
real, intent(in) :: pload(Dim,nn)
real, intent(inout) :: inc_pload(Dim,nn)


if (rn == 1) then
    inc_pload = pload - fsi_inte_set%prop_pre_pload
else
    inc_pload = pload - fsi_inte_set%case_pre_pload
endif

end subroutine difference_pload