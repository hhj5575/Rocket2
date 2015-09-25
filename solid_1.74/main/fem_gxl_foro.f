subroutine fem_gxl_foro(Dim, action, n_remesh, remesh, init_step, delta_t, flag)

! Copyright(2006~) LASSCOM: Large Structural System Computing Lab, DGU, Seoul
! This unit is originally coded by Jeeho Lee (October 28, 2006)
! Modification History: December 05, 2006
!                       October 31, 2009
!                       May 2013

use file_info_house
use system_info_house
use system_domain
!use coarse_domain
use ale_domain

implicit none

integer, intent(in) :: Dim, action, n_remesh, flag  ! flag = 1: not remesh, 2: remesh
integer, intent(inout) :: remesh
integer, intent(inout) :: init_step
real, intent(in) :: delta_t

integer :: num_rn, step_flag, id
integer :: j, rn, rn_ndof, status
logical :: step_virgin = .TRUE.

integer :: unit_num_node_bak, unit_num_cont_bak
character(len=40) :: bak_node_file, bak_cont_file
logical :: existence
! step_flag = 
!      -1: EOF or stop
!				0: initial stage (problem type: static, quasi-static, dynamic, nonlinearity)
!				1: preparing stage (solver type, sparse, symmetricity options)
!				2: preparing result output & post-processing output
!				3: running solver stage (steps)

step_flag = 0
num_rn = num_of_regions ! global variable from physical_domain module
write(*,*) '>>FEM_GXL: start pre-execution procedure.'

init_step = 0
remesh = 0
do  ! loop over command lines reading
	if (step_flag < 0) then
		EXIT  ! finish pre-execution procedure
	endif
    
    call read_pre_step(unit_number, step_flag, remesh) ! <-------------------------
  
end do ! loop over command lines reading

if (step_flag == -3) then ! when next execution is '*step' 

    if (action == 0) step_virgin = .TRUE.  ! at present

    if (step_virgin) then  
        call create_sdomain(num_rn)
        
        if (solver_flag) then
            if (system_sparse) then
                do j = 1, num_rn
                    call create_sparse_system(j)
                end do
            else
                do j = 1, num_rn
                    call create_full_system(j)
                end do
            endif      
        endif

        id = 1  ! temporary setting

        if (is_coarse_problem .AND.(.NOT. coarse_schur)) call set_global_sparse(id)

        call create_history(num_rn)
        write(*,*) ' >> FEM_GXL: create_history DONE.'
  
        !if (init_step == 0) then
            call write_init_outputfile(Dim)
            if (Dim == 3) call write_bound_outputfile(Dim, sub_mesh_region, n_remesh)
            !if (problem_type < 3) then ! 'system_domain': static problem case
            !  init_step = 1
            !endif
        !endif

        do rn = 1, num_rn
            rn_ndof = region(rn)%num_dofs + region(rn)%num_constr
            write(*,*) 'initialize history DB: size =', rn_ndof
            call initialize(problem_type, rn, rn_ndof)  ! time_domain via system_domain
            call update_system_changed(1, rn, .FALSE.)
        end do
        write(*,*) ' >>FEM_GXL: initialize DONE.'
        
        step_virgin = .FALSE.
    endif  	
    
    ! restart routine
    if (init_step > 0 .and. solver_flag) then
        call recover_step(Dim, remesh, num_rn, init_step, delta_t, flag)
    endif
endif

! set 'nonlinear_geometric_load': s_domain module protected global variable
call write_load_geometric_nonlinearity

allocate(loadf(num_load_sets))  ! 'system_info_house'
allocate(dispf(num_disp_sets))  ! 'system_info_house'   
loadf = 0.0
dispf = 0.0

write (*,*) '>>> back up data of main file'
!if (remesh == 0) then
    call back_up_main_file(Dim, unit_number)
!else
!    call back_up_main_file(Dim, unit_num_remesh)
!endif
write (*,*) 'init_step:', init_step

contains
!==============================================================================

subroutine read_pre_step(unit_num, step_flag, remesh)

! global variables from subroutine event_control: unit_num

implicit none

integer, intent(in) :: unit_num
integer, intent(inout) :: step_flag, remesh
integer :: flag, stat, problem_type_number
real :: v1, v2, v3
character(len=5) :: keyword
character(len=5) :: problem, method, option
character(len=3) :: toggle

if (step_flag == 0) then
	do
		call find_star(unit_num, flag)
		
		if (flag < 0) EXIT   ! exit from do-while: EOF case
		
		read(unit_num, *, iostat=stat) keyword
		write(*,*) ' > read keyword-3:   ', step_flag, keyword
		if (keyword == '*star') then
			backspace(unit_num)
			read(unit_num, *, iostat=stat) keyword, problem
			
			select case(problem)
            case('none') 
                problem_type_number = 0
                solver_flag = .FALSE.
			case('stati')
				problem_type_number = 1
            case('quasi')
				problem_type_number = 2
			case('dynam', 'lindy')
				problem_type_number = 3	
			case('nlndy','nldyn')
				problem_type_number = 4
			case default
				stop 'read_step: not a valid problem type!'
			end select
						
            call write_problem_type(problem_type_number)  ! write on 'problem_type' in 'system_domain'
      
			if (problem_type_number >= 3) then
            read(unit_num, *, iostat=stat) method, v1, v2, v3
            backspace(unit_num)
            call intgrn_para_setup(method, v1, v2, v3)
      endif
			step_flag = 1
			EXIT   ! exit from do-while to read the next step line
		endif		
	end do
endif

if (step_flag > 0) then
	
	call find_star(unit_num, flag)
	if (flag > 0) then   ! otherwise: EOF case
		
		read(unit_num, *, iostat=stat) keyword
		write(*,*) ' > read keyword-4:   ', step_flag, keyword		
		select case(keyword)
		case('*solv')
			backspace(unit_num)
			read(unit_num, *, iostat=stat) keyword, solver_kind, tolerance, max_itr
			step_flag = 1
            write(*,*) '  > Solver kind:  ', solver_kind
            write(*,*) '  > tolerance & maximum iteration:  ', tolerance, max_itr
		case('*spar')
			backspace(unit_num)
			read(unit_num, *, iostat=stat) keyword
			step_flag = 1
            call write_system_type(.TRUE., .FALSE.)
            write(*,*) '  > Sparse matrix system: ON'
		case('*symm')
			backspace(unit_num)
			read(unit_num, *, iostat=stat) keyword
			step_flag = 1
            call write_system_type(.FALSE., .TRUE.)
            write(*,*) '  > Symmetic matrix system: ON'
        case('*schu')
			backspace(unit_num)
			read(unit_num, *, iostat=stat) keyword, toggle
			step_flag = 1
            if (toggle(1:2) == 'on') then
                call write_domain_decomposition_type(1)
                write(*,*) '  > Schur complement: ON'
            elseif (toggle == 'off') then
                call write_domain_decomposition_type(0)
                write(*,*) '  > Schur complement: OFF'
            elseif (toggle == 'jac') then
                call write_domain_decomposition_type(2)
                write(*,*) '  > Jacobi preconditioner: ON'
            endif
 		case('*hydr')
			backspace(unit_num)
			read(unit_num, *, iostat=stat) keyword, option
			if (option(1:3) == 'fsi') then
				fsi_flag = .TRUE.
				write(*,*) '  > Hydro dynamic:  FSI on!'
			endif
			step_flag = 1
        case('*outn','*outp','*outr')
			backspace(unit_num)
			read(unit_num, *, iostat=stat) keyword, output_cnt, output_region1, output_node
            if (output_region1 > num_rn) stop 'eventcontrol: invalid region number for nodal output!'
			step_flag = 2
		case('*oute')
			backspace(unit_num)
			read(unit_num, *, iostat=stat) keyword, str_out_cnr, output_region2, output_interval, output_interval2
            if (output_region2 > num_rn) stop 'eventcontrol: invalid region number for element output!'
            element_output_flag = .TRUE.
			step_flag = 2
        case('*auto')
            auto_time_stepping = .TRUE.
            step_flag = 2  
        case('*step')
            backspace(unit_num)
			step_flag = -3
		case('*stop')
			step_flag = -1
        ! add keyword for fsi
        case('*init')
            backspace(unit_num)
            read(unit_num, *, iostat=stat) keyword, init_step
            write (*,*) keyword, init_step
            step_flag = 2
        case('*reme')
            backspace(unit_num)
            read(unit_num, *, iostat=stat) keyword, remesh
            step_flag = 2
		case default
			STOP 'read_step: invalid pre-execution command input!'
		end select
  else
    step_flag = -1  ! EOF case
	endif
endif

end subroutine read_pre_step

!==============================================================================

subroutine write_init_outputfile(Dim)

use output_domain

implicit none

integer, intent(in) :: Dim

integer :: i, mat_num, nn, ne, rn, tec_unit_num, ele_node, ele_num
integer, allocatable :: elem(:,:), conn(:)
real, allocatable :: coord(:,:)
real, allocatable :: zero_val(:)
character(len=40) :: fn

if (Dim == 2) then
    allocate(zero_val(12))
    zero_val = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
elseif (Dim == 3) then
    allocate(zero_val(13))
    zero_val = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
endif

fn = "./output/solid/tec/tec_eldata0000000.plt"
tec_unit_num = 114
open (Unit=tec_unit_num, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)
Write (tec_unit_num,'(A)') 'TITLE="Xfemout_quad_data"'
if (Dim == 2) then
    Write (tec_unit_num,'(A)') 'VARIABLES= "x", "y", "dx", "dy", "stress11", "stress22", "stress33", "stress12", "strain11", "strain22", "strain33", "strain12", "prncpl_1", "prncpl_2"'
elseif (Dim == 3) then
    Write (tec_unit_num,'(A,A)') 'VARIABLES= "x", "y", "z", "dx", "dy", "dz", "stress11", "stress22", "stress33", "stress12", "strain11", "strain22", "strain33", "strain12", "prncpl_1", "prncpl_2"'
endif
do mat_num = 1, mat_mesh_set_num
    nn = mat_mesh_set(mat_num)%nn
    ne = mat_mesh_set(mat_num)%ne
    rn = mat_mesh_set(mat_num)%rn
    ele_num = mat_mesh_set(mat_num)%ele_num
    if (ele_num > 200 .and. ele_num < 400) then
        ele_node = (Dim-1)*4
        allocate (coord(Dim,nn), elem(ele_node,ne), conn(nn))
        conn = mat_mesh_set(mat_num)%inv_conn
        elem = mat_mesh_set(mat_num)%elem
        coord = region(rn)%node_coord(:,conn)
    
        if (Dim == 2) then
            Write (tec_unit_num,'(A, I1,A,I5,A,I5, A)') 'ZONE SolutionTime = 0.0, T = "Mat_num', mat_num, '", N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
            do i = 1, nn
                write(tec_unit_num,'(100(ES15.7,2x))') coord(2,i), coord(1,i), zero_val
            enddo
        elseif (Dim == 3) then
            Write (tec_unit_num,'(A, I1,A,I5,A,I5, A)') 'ZONE SolutionTime = 0.0, T = "Mat_num', mat_num, '", N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEBRICK'
            do i = 1, nn
                write(tec_unit_num,'(100(ES15.7,2x))') coord(:,i), zero_val
            enddo
        endif
    
        do i = 1, ne
            write(tec_unit_num,*) elem(:,i)
        enddo
        deallocate(coord, elem, conn)
    elseif (ele_num == 112 .or. ele_num == 113) then
        ele_node = 2
        allocate (coord(Dim,nn), elem(ele_node,ne), conn(nn))
        conn = mat_mesh_set(mat_num)%inv_conn
        elem = mat_mesh_set(mat_num)%elem
        coord = region(rn)%node_coord(:,conn)
        if (Dim == 2) then
            Write (tec_unit_num,'(A, I1,A,I5,A,I5, A)') 'ZONE SolutionTime = 0.0, T = "Mat_num', mat_num, '", N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
            do i = 1, nn
                write(tec_unit_num,'(100(ES15.7,2x))') coord(2,i), coord(1,i), zero_val
            enddo
        elseif (Dim == 3) then
            Write (tec_unit_num,'(A, I1,A,I5,A,I5, A)') 'ZONE SolutionTime = 0.0, T = "Mat_num', mat_num, '", N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FELINESEG'
            do i = 1, nn
                write(tec_unit_num,'(100(ES15.7,2x))') coord(:,i), zero_val
            enddo
        endif
        do i = 1, ne
            write(tec_unit_num,*) elem(:,i)
        enddo
    endif
enddo
close(tec_unit_num)

end subroutine write_init_outputfile

!==============================================================================
subroutine write_bound_outputfile(Dim, num_group, n_remesh)

use output_domain

implicit none

integer, intent(in) :: Dim, num_group, n_remesh

integer :: i, nn, ne, rn, tec_unit_num, status, elem(4)
integer, allocatable :: couple(:)
real, allocatable :: bn_coord(:,:)
character(len=40) :: fn

fn = "./output/solid/Init_bound_face.plt"
tec_unit_num = 114
open (Unit=tec_unit_num, File=fn, STATUS='replace', ACTION='write', IOSTAT=status)
Write (tec_unit_num,'(A)') 'TITLE="Boundary face"'
if (Dim == 3) then
    Write (tec_unit_num,'(A,A)') 'VARIABLES= "x", "y", "z"'
    do rn = 1, divided_mesh_region
        nn = mesh_set(rn)%bn
        ne = mesh_set(rn)%bfn
        allocate(couple(mesh_set(rn)%nn), bn_coord(3,nn))
        couple = 0
        Write (tec_unit_num,'(A, I1,A,I5,A,I5, A)') 'ZONE T = "Mat_num', rn, '", N =', nn , ', E =', ne, ', DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
        do i = 1, nn
            bn_coord(:,i) = mesh_set(rn)%node_coord(:,mesh_set(rn)%bound_node(i))
	        write(tec_unit_num,*) bn_coord(:,i) 
	        couple(mesh_set(rn)%bound_node(i)) = i
        enddo
        do i = 1, ne
            !elem = couple(mesh_set(rn)%bound_face(:,i))
            elem = mesh_set(rn)%bound_face(:,i)
            write (tec_unit_num,'(4I7)') elem
        enddo
        deallocate(couple)
        !if (rn == 1 .and. n_remesh /= 0) then
        !    !call revise_add_data_for_remesh(nn, bn_coord, mesh_set(rn)%bound_lv)
        !    call revise_b_edge_info(num_group, nn, bn_coord)
        !endif
        deallocate (bn_coord)
    enddo
    !call free_geo_domain(2)
    !call create_remesh_domain(num_group, 2)
endif
close(tec_unit_num)

end subroutine write_bound_outputfile


end subroutine fem_gxl_foro
