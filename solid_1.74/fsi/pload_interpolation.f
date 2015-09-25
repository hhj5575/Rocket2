subroutine pload_interpolation(max_interface_pt, old_ni, old_pload, old_inte_coord, new_ni, new_inte_coord)

implicit none

integer, intent(in) :: max_interface_pt
integer, intent(inout) :: old_ni
real, intent(inout) :: old_pload(max_interface_pt), old_inte_coord(2,max_interface_pt)
integer, intent(in) :: new_ni
real, intent(in) :: new_inte_coord(2,max_interface_pt)

integer :: i, j, memory_j
real :: dis, tol_dis, val, ratio, p(2,2), sum_old_load, sum_new_load
real :: new_pload(new_ni), old_arr(old_ni), new_arr(new_ni)

tol_dis = 0.0
old_arr(1) = 0.0
sum_old_load = 0.0
do i = 1, old_ni-1
    p(:,1) = old_inte_coord(:,i)
    p(:,2) = old_inte_coord(:,i+1)
    call calc_len(p(:,1), p(:,2), dis)
    tol_dis = tol_dis + dis
    old_arr(i+1) = tol_dis
    sum_old_load = sum_old_load + sum(old_pload(i:i+1))*0.5*dis
enddo
tol_dis = 0.0
new_arr(1) = 0.0
do i = 1, new_ni-1
    p(:,1) = new_inte_coord(:,i)
    p(:,2) = new_inte_coord(:,i+1)
    call calc_len(p(:,1), p(:,2), dis)
    tol_dis = tol_dis + dis
    new_arr(i+1) = tol_dis
enddo
ratio = old_arr(old_ni)/new_arr(new_ni)
new_arr = new_arr*ratio
new_pload(1) = old_pload(1)
new_pload(new_ni) = old_pload(old_ni)

memory_j = 1
do i = 2, new_ni-1
    do j = memory_j, old_ni-1
        if (new_arr(i) <= old_arr(j+1)) then
            memory_j = j
            exit
        endif
    enddo
    
    val = old_pload(memory_j+1)-old_pload(memory_j)
    val = val/(old_arr(memory_j+1) - old_arr(memory_j))
    val = val*(new_arr(i)-old_arr(memory_j))
    
    new_pload(i) = old_pload(memory_j) + val
enddo

sum_new_load = 0.0
do i = 1, new_ni-1
    p(:,1) = new_inte_coord(:,i)
    p(:,2) = new_inte_coord(:,i+1)
    call calc_len(p(:,1), p(:,2), dis)
    sum_new_load = sum_new_load + sum(new_pload(i:i+1))*0.5*dis
enddo

write (*,*) 'old_ni, new_ni:', old_ni, new_ni
write (*,*) 'old sum load:', sum_old_load
write (*,*) 'new sum load:', sum_new_load
!write (*,'(10F12.5)') new_pload(1:10)

old_ni = new_ni
old_pload(1:new_ni) = new_pload(1:new_ni)
old_inte_coord(:,1:new_ni) = new_inte_coord(:,1:new_ni)

end subroutine pload_interpolation