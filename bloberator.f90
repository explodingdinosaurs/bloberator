!######################################################################
! bloberator
!
! usage: bloberator input.in
!
!
! This program constructs pqr files of two "rigid body" proteins
! by rotation and translation to produce various spacial orientations. 
!
! Uses md_mod.f90 
! 
!
!
!
! Michael Thomas, June 2016
!######################################################################

  
!-------------------------------------------------------------------------------
!   PROGRAMMMMMMMMMMMMMMMMMM! 
!-------------------------------------------------------------------------------
  program bloberator
    
    use md_mod
    use bloberator_mod
    
    implicit none

      ! Global Variables
    integer, parameter :: CONF_UNIT = 80
    integer, parameter :: PQR_UNIT1 = 81, PQR_UNIT2 = 82
    integer, parameter :: PQR_OUTUNIT = 90
    integer, parameter :: LOG_OUTUNIT = 91
    ! Be careful with the following parameters. This is probably the
    ! main source of memory usage
    
    character*80 :: conf_filename, pqr_filename1, pqr_filename2
    character*80 :: output_pqr_filename, log_filename, min_filename
    real :: cutoff_dist, grid_spacing, overlap_dist
    double precision :: theta, phi, gamma
    integer, dimension(:), allocatable :: resid1, resid2
    real, dimension(:,:,:), allocatable :: unique_coords
    real, dimension(:,:), allocatable :: coords1_origin, coords1, tcoords1, coords2, tcoords2
    real, dimension(:,:), allocatable :: coords2_origin, coords2_gridstart
    real, dimension(:,:), allocatable :: coords2_t, coords2_rt, coords2_rrt, coords2_rrrt
    real, dimension(:,:), allocatable :: exp_constraints
    real, dimension(1,3) :: overlap_coords
    real, dimension(:), allocatable :: o1, b1, o2, b2
    character(4), dimension(:), allocatable :: atomtype1, atom1, atomtype2, atom2
    character(4), dimension(:), allocatable :: restype1, chain1, restype2, chain2
    integer :: n_atoms1, n_atoms2, n_grid_points
    integer :: n_theta_rotations, n_phi_rotations, n_gamma_rotations
    integer :: n_unique_coords
    real, dimension(3) :: CTER_C_coords, NTER_N_coords
    real, dimension(3) :: small_TER_coords, big_TER_coords, fuck_coords
    integer :: i, j, k, a, b, c, y
    logical :: is_coords_valid, within_cutoff, is_unique, no_overlap, valid_exp_constraints
    integer :: total_counter, valid_counter
    real :: min_dist
    double precision :: min_theta, min_phi, min_gamma
    character(1) :: small_n_or_c, big_n_or_c
    integer :: n_coords_arrays_in_check_storage
    logical, dimension(:,:,:), allocatable :: unique 
    double precision, dimension(:,:,:,:,:), allocatable :: array_of_rot_arrays
    real :: check_start, check_end
    real :: total_start, total_end

    call cpu_time(total_start)
    
    ! Open, read and close config file
    call process_commandline_conf_file(conf_filename)
    call open_input_file_stop_if_not_found(CONF_UNIT, conf_filename)
    call read_conf_file(CONF_UNIT, pqr_filename1, pqr_filename2, &
         cutoff_dist, grid_spacing, theta, phi, gamma, overlap_dist,&
         n_coords_arrays_in_check_storage, exp_constraints)
    call close_file(CONF_UNIT)

    ! Open log file and write header
    call generate_log_filename(pqr_filename1, pqr_filename2, log_filename)
    call open_output_file(LOG_OUTUNIT, log_filename)
    call write_log_file_header(LOG_OUTUNIT, conf_filename, pqr_filename1, pqr_filename2, cutoff_dist, &
         grid_spacing, theta, phi, gamma, overlap_dist, exp_constraints)
    
    !open the two pqr files and read in variables, then close.
    ! pqr1 is protein with CTER connected to linker
    ! pqr2 is protein with NTER connected to linker
    call open_input_file_stop_if_not_found(PQR_UNIT1, pqr_filename1)
    call open_input_file_stop_if_not_found(PQR_UNIT2, pqr_filename2)
    call read_pqr_file(PQR_UNIT1, atomtype1, coords1, atom1, restype1, resid1, chain1, n_atoms1, o1, b1)
    call read_pqr_file(PQR_UNIT2, atomtype2, coords2, atom2, restype2, resid2, chain2, n_atoms2, o2, b2)
    call close_file(PQR_UNIT1)
    call close_file(PQR_UNIT2)

    ! Determine smallest protein. If nter is smallest, swap all pqr variables
    call get_small_protein(coords1, coords1, small_n_or_c, big_n_or_c)
    call reassign_pqr_variables(small_n_or_c,&
         atomtype1, coords1, atom1, restype1, resid1, chain1, n_atoms1, o1, b1, &
      atomtype2, coords2, atom2, restype2, resid2, chain2, n_atoms2, o2, b2)

    ! Initiate check coord array
    allocate(unique_coords(n_atoms2, 3, n_coords_arrays_in_check_storage))
    unique_coords = 0



    ! Move pqr1 protein (reference protein) so that it's CTER is at the origin
    call get_TER_coords(coords1, resid1, atomtype1, big_n_or_c, big_TER_coords)
    call translate_coords(coords1, -1*big_TER_coords, coords1_origin)

    ! Move pqr2 protein so that it's NTER is at the origin too
    call get_TER_coords(coords2, resid2, atomtype2, small_n_or_c, small_TER_coords)
    call translate_coords(coords2, -1*small_TER_coords, coords2_origin)

    ! Generate logfile name
    call generate_log_filename(pqr_filename1,pqr_filename2,log_filename)
    
    ! Now we move pqr2 protein to the start of the grid

    call get_TER_coords(coords2_origin, resid2, atomtype2, small_n_or_c, small_TER_coords)
    call translate_coords(coords2_origin, [-1*cutoff_dist, -1*cutoff_dist, -1*cutoff_dist], coords2_gridstart)
    call get_TER_coords(coords2_gridstart, resid2, atomtype2, small_n_or_c, small_TER_coords)
   
    ! Now cycle pqr2 protein through the grid, and at each gridpoint perform rotations about NTER_N_coords
    ! Run checks and if it passes, write new coordinates to file.
    n_grid_points = int((2*cutoff_dist)/grid_spacing)
    n_theta_rotations = int((2*PI)/theta)
    n_phi_rotations = int((2*PI)/phi)
    n_gamma_rotations = int((2*PI)/(2*gamma)) 

    total_counter = 0
    valid_counter = 0
    min_dist = 2*cutoff_dist
    min_theta = 0
    min_phi = 0
    min_gamma = 0

    print*, "------Checking translation matrix uniqueness-----"

    ! Do a translation matrix check here. This will replace uniqueness check
    call translation_matrix_check(theta, phi, gamma, n_theta_rotations, n_phi_rotations, n_gamma_rotations, &
         unique, array_of_rot_arrays)

    ! Start doing awesome stuff
    do i=0, n_grid_points
       ! Print of % progres
       print*, int(100*real(i*n_grid_points**2 * n_theta_rotations * n_phi_rotations * n_gamma_rotations)/&
            real(n_grid_points**3 * n_theta_rotations * n_phi_rotations * n_gamma_rotations)), "% Complete"
              
       do j=0, n_grid_points
          do k=0, n_grid_points

             call cpu_time(check_start)
             call translate_coords(coords2_gridstart, [i*grid_spacing, j*grid_spacing, k*grid_spacing], coords2_t) 

             call get_TER_coords(coords2_t, resid2, atomtype2, small_n_or_c, fuck_coords)
             
             ! Get small protein TER c or n coords
             call get_TER_coords(coords2_t, resid2, atomtype2, small_n_or_c, small_TER_coords)

             ! Get large protein TER c or n coords (should be 0 0 0)
             call get_TER_coords(coords1_origin, resid1, atomtype1, big_n_or_c, big_TER_coords)

             overlap_coords(1,:) = small_TER_coords
             
             ! Check to see if distance is valid, saves going around angle loops
             ! See if small TER is overlapping big protein.
             ! If not valid, exit out of k loop.
             call within_cutoff_distance(big_TER_coords, small_TER_coords, cutoff_dist, within_cutoff)

             if (within_cutoff .eqv. .FALSE.) then
                cycle
             end if
             
             ! Check to see if small_TER_atom overlaps with any big_TER protein
             call is_no_overlap(coords1, overlap_coords, overlap_dist, no_overlap)
             
             if (no_overlap .eqv. .FALSE.) then
                cycle
             end if
             
             ! Rotate around new grid position
             do a=1, n_theta_rotations
                do b=1, n_phi_rotations
                   do c=1, n_gamma_rotations
                      total_counter = total_counter + 1
                      
                      ! Use matrix uniqueness test 
                      if (unique(a,b,c) .eqv. .FALSE.) then
                         cycle
                      end if

                      call rotate_coords_with_rot_mat(coords2_t, small_TER_coords, array_of_rot_arrays(a,b,c,:,:), coords2_rrrt)
                      
                      ! Perform expirement constraint check
                      call within_exp_constraints(exp_constraints, coords1_origin, coords2_rrrt, &
                           small_n_or_c, valid_exp_constraints)
                      
                      if (valid_exp_constraints .eqv. .FALSE.) then 
                         cycle
                      end if
                                            
                      ! Peform overlap check
                      call is_no_overlap(coords1_origin, coords2_rrrt, overlap_dist, no_overlap)
                           
                      if (no_overlap .eqv. .FALSE.) then   
                         cycle
                      end if
                                            
                      valid_counter = valid_counter + 1
                            
                      ! Generate and open output pqr file
                      call generate_output_filename( (/ (i*grid_spacing)-cutoff_dist, (j*grid_spacing)-cutoff_dist,&
                           (k*grid_spacing)-cutoff_dist /), (a-1)*theta, (b-1)*phi, (c-1)*gamma, output_pqr_filename)
                      call open_output_file(PQR_OUTUNIT, output_pqr_filename)
                            
                      ! Check min distance and update if min dist
                      call get_TER_coords(coords2_rrrt, resid2, atomtype2, small_n_or_c, small_TER_coords)
                      call get_TER_coords(coords1_origin, resid1, atomtype1, big_n_or_c, big_TER_coords)
                      
                      if (distance(big_TER_coords, small_TER_coords) .lt. min_dist) then
                         min_dist = distance(big_TER_coords, small_TER_coords)
                         min_theta = (a-1)*theta
                         min_phi = (b-1)*phi
                         min_gamma = (c-1)*gamma
                         min_filename = output_pqr_filename
                      end if
                            
                      ! Have to write both proteins to pqr
                      ! Write NTER first, then CTER (NTER is where the probe attaches - so the actual
                      ! CTER is the where the probe attaches to the NTER, and the actual NTER is where
                      ! the probe attaches to the CTER. Confusing hey. 
                      
                      if (big_n_or_c == 'c') then
                         call write_to_output_pqr(PQR_OUTUNIT, 1, atomtype1, coords1_origin, atom1,&
                              restype1, resid1, chain1, n_atoms1, o1, b1)
                         call write_to_output_pqr(PQR_OUTUNIT, n_atoms1+1, atomtype2, coords2_rrrt, atom2,&
                              restype2, resid2, chain2, n_atoms2, o2, b2)
                      elseif (small_n_or_c == 'c') then
                         call write_to_output_pqr(PQR_OUTUNIT, 1, atomtype2, coords2_rrrt, atom2,&
                              restype2, resid2, chain2, n_atoms2, o2, b2)
                         call write_to_output_pqr(PQR_OUTUNIT, n_atoms2+1, atomtype1, coords1_origin, atom1,&
                              restype1, resid1, chain1, n_atoms1, o1, b1)
                      else
                         write(ERROR_UNIT,*) "Definition of NTER and CTER is fucked up"
                         stop
                      end if
                      
                      call close_file(PQR_OUTUNIT)
                      
                      ! Write to log file
                      call write_valid_to_log_file(LOG_OUTUNIT, distance(small_TER_coords, big_TER_coords),&
                           (a-1)*theta, (b-1)*phi, (c-1)*gamma, output_pqr_filename)
                      
                      
                   end do
                end do
             end do
             call cpu_time(check_end)
             print*, "Checked grid index:", i, j, k, "    Time taken:", check_end-check_start, "seconds"
          end do
       end do
    end do

    call write_log_file_summary(LOG_OUTUNIT, valid_counter, total_counter)
    write(LOG_OUTUNIT,*) "Closest fit:"
    call write_valid_to_log_file(LOG_OUTUNIT, min_dist, min_theta, min_phi, min_gamma, min_filename)

    close(LOG_OUTUNIT)

    call cpu_time(total_end)
    print*, "Total time:", total_end-total_start, "seconds"
             
  end program bloberator
