subroutine build_interface(action, unit_num, number_contacts)

! This unit is originally coded by Jeeho Lee (February 23, 2012)
! Modification History:		 
! External subroutine
!
! action = 1: create & build all types of contacts
!        = 2: rebuild only for contact_type = 2

use coarse_domain

implicit none

integer, intent(in) :: action, unit_num, number_contacts
character(len=5) :: keyword
integer :: i, j, k, id, reg_from, reg_to, contact_type, status, num_contact_group = 1, contact_dim, num_points
real, allocatable :: contact_pt_from(:,:,:), contact_pt_to(:,:,:)


if (action == 1) call create_contact(number_contacts)

if (number_contacts > 0) then
  rewind(unit=unit_num)
  
  do  
    read(unit_num, *, iostat=status) keyword

    if (status == -1 .or. keyword == '*end ') then
      write(*,*) '*** EOF in contact info file ! -----------------------------'
      EXIT
    elseif (keyword == '*cont') then
      backspace(unit_num)
      read(unit_num,*) keyword, id, reg_from, reg_to, contact_type, contact_dim, num_contact_group

      if ((contact_dim < 2) .OR. (contact_dim > 3)) STOP 'build_interface: contact dimension must be in the range of [2, 3]!'

      if (id > number_contacts) then
        STOP 'build_interface: contact id number error!'
      else
        write(*,'(A/,6I10)') 'read contact id, region from, region to, contact type, contact_dim, num_contact_group:', id, reg_from, reg_to, contact_type, contact_dim, num_contact_group

        if (contact_dim == 3) then
          num_points = 4
        else
          num_points = 2
        endif

        allocate(contact_pt_from(contact_dim,num_points,num_contact_group))
        allocate(contact_pt_to(contact_dim,num_points,num_contact_group))
      endif


      do i = 1, num_contact_group
        read (unit_num, *) ((contact_pt_from(j,k,i), j=1,contact_dim), k=1,num_points)
        read (unit_num, *) ((contact_pt_to(j,k,i), j=1,contact_dim), k=1,num_points)
      end do
      if (action == 1) then
        call write_contact(action, id, reg_from, reg_to, contact_type, contact_pt_from, contact_pt_to, contact_dim, num_points, num_contact_group)
      elseif ((action == 2) .AND. (contact_type == 2)) then
        call write_contact(action, id, reg_from, reg_to, contact_type, contact_pt_from, contact_pt_to, contact_dim, num_points, num_contact_group)
      endif
    endif
    deallocate(contact_pt_from, contact_pt_to)
  enddo

endif

end subroutine build_interface
