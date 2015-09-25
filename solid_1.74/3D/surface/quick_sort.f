module qsort_c_module

implicit none
public :: QsortC
private :: Partition_Q

contains

recursive subroutine QsortC(A, col)
  integer, intent(inout), dimension(:,:) :: A
  integer, intent(in) :: col
  integer :: iq

  if(size(A(1,:)) > 1) then
     call Partition_Q(A, iq, col)
     call QsortC(A(:,:iq-1), col)
     call QsortC(A(:,iq:), col)
  endif
end subroutine QsortC

subroutine Partition_Q(A, marker, col)
  integer, intent(inout), dimension(:,:) :: A
  integer, intent(in) :: col
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp(col)
  real :: x      ! pivot point
  x = A(1,1)
  i= 0
  j= size(A(1,:)) + 1

  do
     j = j-1
     do
        if (A(1,j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(1,i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(:,i)
        A(:,i) = A(:,j)
        A(:,j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
end subroutine Partition_Q

end module qsort_c_module