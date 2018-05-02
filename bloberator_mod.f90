!-------------------------------------------------------------------------------
!  - bloberator_mod.f90 is the library for bloberator.f90
!  - it requires the md_mod library for compilation
!  - Michael Thomas Apr 2018
!-------------------------------------------------------------------------------


module bloberator_mod
  use md_mod
  implicit none
  
contains

!-------------------------------------------------------------------------------
!   Read basename from commandline
!-------------------------------------------------------------------------------
    subroutine process_commandline_conf_file(input_filename)
      character*80, intent(out) :: input_filename
      character(len=80) arg

      ! Read input_filename
      call get_command_argument(1, input_filename)
      if (len_trim(input_filename) == 0) call print_usage_then_stop

    end subroutine process_commandline_conf_file


  
!-------------------------------------------------------------------------------
!   Print out the correct usage of commandline args and then stop
!-------------------------------------------------------------------------------
  subroutine print_usage_then_stop
    write(ERROR_UNIT,*) "Incorrect number of arguments. Correct Usage:"
    write(ERROR_UNIT,*) "bloberator <input_file>"
    write(ERROR_UNIT,*) "Example:"
    write(ERROR_UNIT,*) "bloberator input.in"
    stop
  end subroutine print_usage_then_stop
  

!-------------------------------------------------------------------------------
!   Read the config input file. Input file should have the following format:
!   - pqr file 1, this should be the protein to which the linker is connected     
!     to the CTER
!   - pqr file 2, this is the protein to with the linker is connected to the 
!     NTER
!   - max length of the linker (also called the cutoff_distance
!   - grid spacing for translation (in angstrom)
!   - theta - angle by which protein is rotated in x axis (input file is degrees, 
!     subroutine outputs radians for internal goodness )   
!   - phi - angle by which protein is rotated in y axis (input file in degress,
!     subroutine outputs in radians)
!   - gamma - angle by which protein is rotated in y axis (input file in degress,
!     subroutine outputs in radians)    
!   - overlap detection cutoff (in angstrom)
!   - experimental constrains. integer. if 0, stop. If n != 0 then read next n
!     lines, index1, index2, center and tolerance    
!-------------------------------------------------------------------------------
    subroutine read_conf_file(conf_file_unit, pqr_filename1, pqr_filename2, &
         cutoff_dist, grid_spacing, theta, phi, gamma, overlap_dist,&
         n_coords_arrays_in_check_storage, exp_constraints)

      integer, intent(in) :: conf_file_unit
      character*80, intent(out) :: pqr_filename1, pqr_filename2
      real, intent(out) :: cutoff_dist, overlap_dist, grid_spacing      
      double precision, intent(out) :: theta, phi, gamma
      real, dimension(:,:), allocatable, intent(out) :: exp_constraints
      integer, intent(out) :: n_coords_arrays_in_check_storage
      
      integer :: i, n_exp_constraints
      double precision :: theta_deg, phi_deg, gamma_deg
      
      read(conf_file_unit,*) pqr_filename1
      read(conf_file_unit,*) pqr_filename2
      read(conf_file_unit,*) cutoff_dist
      read(conf_file_unit,*) grid_spacing
      read(conf_file_unit,*) theta_deg
      read(conf_file_unit,*) phi_deg
      read(conf_file_unit,*) gamma_deg
      read(conf_file_unit,*) overlap_dist
      read(conf_file_unit,*) n_coords_arrays_in_check_storage
      read(conf_file_unit,*) n_exp_constraints

      ! Is there are experimental constraints?
      ! If no.. create null array..???
      if (n_exp_constraints == 0 ) then
         ! Set null array?
         allocate(exp_constraints(0,0))
         ! If there are, allocate array
      else
         allocate(exp_constraints(n_exp_constraints,4))

         ! Read experimental constraints
         do i=1, n_exp_constraints
            read(conf_file_unit,*) exp_constraints(i,1), exp_constraints(i,2),&
                 exp_constraints(i,3), exp_constraints(i,4)
         end do
      end if
      
      ! Convert theta and phi from degrees to radians
      theta = (theta_deg * PI)/180
      phi = (phi_deg * PI)/180
      gamma = (gamma_deg * PI)/180
      
    end subroutine read_conf_file

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------
!   get_small_protein: see if CTER or NTER protein is smallest.
!   Output n if nter is smallest, c if cter is smallest   
!-------------------------------------------------------------------------------   

    subroutine get_small_protein(coords1, coords2, small_n_or_c, big_n_or_c)
      real, dimension(:,:), intent(in) :: coords1, coords2
      character(1), intent(out) :: small_n_or_c, big_n_or_c

      if (size(coords1,1) .le. size(coords2,1)) then
         small_n_or_c = "c"
         big_n_or_c = "n"
      else
         small_n_or_c = "n"
         big_n_or_c = "c"
      end if

    end subroutine get_small_protein


!-------------------------------------------------------------------------------
!   reassign_pqr_variables: if NTER is smallest, swap all variables around - 
!   e.g coords1 -> coords2, coords2 -> coords1
!   if CTER is smallest do nothing.    
!-------------------------------------------------------------------------------   

    subroutine reassign_pqr_variables(small_n_or_c,&
         atomtype1, coords1, atom1, restype1, resid1, chain1, n_atoms1, o1, b1, &
      atomtype2, coords2, atom2, restype2, resid2, chain2, n_atoms2, o2, b2)

      character(1), intent(in) :: small_n_or_c
      integer, dimension(:), allocatable, intent(inout) :: resid1, resid2 
      real, dimension(:,:), allocatable, intent(inout) :: coords1, coords2
      real, dimension(:), allocatable, intent(inout) :: o1, b1, o2, b2
      character(4), dimension(:), allocatable, intent(inout) :: atomtype1, atom1, restype1, chain1
      character(4), dimension(:), allocatable, intent(inout) :: atomtype2, atom2, restype2, chain2
      integer, intent(inout) :: n_atoms1, n_atoms2
      
      integer, dimension(size(resid1,1)) :: resid1t 
      real, dimension(size(coords1,1),size(coords1,2)) :: coords1t
      real, dimension(size(o1,1)) :: o1t, b1t
      character(4), dimension(size(atomtype1,1)) :: atomtype1t, atom1t, restype1t, chain1t

      integer, dimension(size(resid2,1)) :: resid2t 
      real, dimension(size(coords2,1),size(coords2,2)) :: coords2t
      real, dimension(size(o2,1)) :: o2t, b2t
      character(4), dimension(size(atomtype2,1)) :: atomtype2t, atom2t, restype2t, chain2t

      
      ! If NTER is smallest, then swap everything. If CTER is smallest, do nothing - program
      ! assumes that input has smallest CTER

      if (small_n_or_c == "c") then
         resid1t = resid1
         coords1t = coords1
         o1t = o1
         b1t = b1
         atomtype1t = atomtype1
         atom1t = atom1
         restype1t = restype1
         chain1t = chain1 

         resid2t = resid2
         coords2t = coords2
         o2t = o2
         b2t = b2
         atomtype2t = atomtype2
         atom2t = atom2
         restype2t = restype2
         chain2t = chain2 

         ! deallocate all arrays
         deallocate(resid1, coords1, o1, b1, atomtype1, atom1, restype1, chain1,&
              resid2, coords2, o2, b2, atomtype2, atom2, restype2, chain2)
         
         ! allocate new arrays
         allocate(resid1(size(resid2t,1)))
         allocate(coords1(size(coords2t,1),size(coords2t,2)))
         allocate(o1(size(o2t,1)))
         allocate(b1(size(b2t,1)))
         allocate(atomtype1(size(atomtype2t,1)))
         allocate(atom1(size(atom2t,1)))
         allocate(restype1(size(restype2t,1)))
         allocate(chain1(size(chain2t,1)))

         allocate(resid2(size(resid1t,1)))
         allocate(coords2(size(coords1t,1),size(coords1t,2)))
         allocate(o2(size(o1t,1)))
         allocate(b2(size(b1t,1)))
         allocate(atomtype2(size(atomtype1t,1)))
         allocate(atom2(size(atom1t,1)))
         allocate(restype2(size(restype1t,1)))
         allocate(chain2(size(chain1t,1)))

         ! Swap values
         resid1 = resid2t
         coords1 = coords2t
         o1 = o2t
         b1 = b2t
         atomtype1 = atomtype2t
         atom1 = atom2t
         restype1 = restype2t
         chain1 = chain2t

         resid2 = resid1t
         coords2 = coords1t
         o2 = o1t
         b2 = b1t
         atomtype2 = atomtype1t
         atom2 = atom1t
         restype2 = restype1t
         chain2 = chain1t   

         n_atoms1 = size(coords1,1)
         n_atoms2 = size(coords2,1)
         
      end if
         
      
      
    end subroutine reassign_pqr_variables


!-------------------------------------------------------------------------------
!   translation_matrix_check: 
!   
!-------------------------------------------------------------------------------
    subroutine translation_matrix_check(theta, phi, gamma, n_theta_rotations, n_phi_rotations, n_gamma_rotations, &
         unique, array_of_rot_arrays)
      
       double precision, intent(in) :: theta, phi, gamma
       integer, intent(in) :: n_theta_rotations, n_phi_rotations, n_gamma_rotations 
       logical, dimension(:,:,:), allocatable, intent(out) :: unique
       double precision, dimension(:,:,:,:,:), allocatable, intent(out) :: array_of_rot_arrays
       
       integer :: a, b, c, d, e, f
       double precision, dimension(3,3) :: theta_rot_mat, phi_rot_mat, gamma_rot_mat
       logical :: equal

       allocate(unique(n_theta_rotations,n_phi_rotations,n_gamma_rotations))
       allocate(array_of_rot_arrays(n_theta_rotations,n_phi_rotations,n_gamma_rotations,3,3))
       
       unique=.TRUE. 
      
       ! generate a 3x3 matrix with all the rotation matrices as entries.. so a (n_a, n_b, n_c, 3, 3) array
       do a=1, n_theta_rotations
          do b=1, n_phi_rotations
             do c=1, n_gamma_rotations
                call rotate((a-1)*theta, 1, theta_rot_mat)
                call rotate((b-1)*phi, 2, phi_rot_mat)
                call rotate((c-1)*gamma, 3, gamma_rot_mat)
                array_of_rot_arrays(a,b,c,:,:) = matmul(matmul(theta_rot_mat, phi_rot_mat), gamma_rot_mat)
             end do
          end do
       end do

       ! Check if they're equal
       do a=1, n_theta_rotations
          do b=1, n_phi_rotations
             
             do c=1, n_gamma_rotations
                if (unique(a,b,c) .eqv. .TRUE.) then 
                   
                   do d=1, n_theta_rotations
                      if (d .ne. a) then 
                         
                         do e=1, n_phi_rotations
                            if (e .ne. b) then 
                               
                               do f=1, n_gamma_rotations
                                  if (f .ne. c) then 
                                    
                                     call are_matrices_equal(array_of_rot_arrays(a,b,c,:,:), array_of_rot_arrays(d,e,f,:,:), equal) 
                                     if (equal .eqv. .TRUE.) then
                                        unique(d,e,f) = .FALSE.
                                     end if

                                  end if
                               end do
                               
                            end if
                         end do
                         
                      end if
                   end do
                   
                end if
             end do
             
          end do
       end do
          
     end subroutine translation_matrix_check

!-------------------------------------------------------------------------------
!   is_no_overlap: read as is_no_overlap?
!   checks across all coordinates of each protein to see if any distance 
!   between atoms is less than some cutoff distance.
!   If there is no overlap, return true
!   if overlap, return false
!   
!   
!-------------------------------------------------------------------------------

    subroutine is_no_overlap(coords1, coords2, overlap_dist, no_overlap)
      real, dimension(:,:), intent(in) :: coords1, coords2
      real, intent(in) :: overlap_dist
      logical, intent(out) :: no_overlap

      integer :: i,j
      real :: overlap_check_start, overlap_check_finish


      call cpu_time(overlap_check_start)
      no_overlap = .TRUE. 

      ! Loop over all coords in protein 1, then protein 2
      !$omp parallel do private(j)
      loop1: do i=1, size(coords1,1)
         do j=1, size(coords2,1)
            
            ! Check if distance between atoms in the two proteins is less than
            ! cutoff distance. If yes, is_no_overlap=false, exit from loops.
            ! If no overlap, is_no_overlap stays as true. 
            if (distance(coords1(i,:),coords2(j,:)) .lt. overlap_dist) then 
               no_overlap = .FALSE.
               exit loop1
            end if

         end do
       end do loop1
       !$omp end parallel do

      call cpu_time(overlap_check_finish)
!      print*, "Overlap check time: ", overlap_check_finish - overlap_check_start

    end subroutine is_no_overlap

!-------------------------------------------------------------------------------
!   within_cutoff_distance: read as within_cutoff_distance?
!   checks to see if NTER N atom and CTER C atom are within cutoff distance
!   of one another
!   Within cutoff, return true
!   not within cutoff, return false
!   
!   
!-------------------------------------------------------------------------------

    subroutine within_cutoff_distance(atom_coord1, atom_coord2, cutoff_dist, within_cutoff)
      real, dimension(3) :: atom_coord1, atom_coord2
      real, intent(in) :: cutoff_dist
      logical, intent(out) :: within_cutoff

      within_cutoff = .TRUE.

      if (distance(atom_coord1, atom_coord2) .gt. cutoff_dist) then
         within_cutoff = .FALSE.
      end if

    end subroutine within_cutoff_distance

!-------------------------------------------------------------------------------
!   within_exp_constraints: read as within_exp_constraints?
!   checks to see if indexes are within distance +/-t olerance
!   Within distance, return true
!   not within distance, return false
!   
!   
!-------------------------------------------------------------------------------

    subroutine within_exp_constraints(exp_constraints, coords1, coords2, small_n_or_c, valid_exp_constraints)
      real, dimension(:,:), intent(in) :: exp_constraints
      real, dimension(:,:), intent(in) :: coords1, coords2
      logical, intent(out) :: valid_exp_constraints
      character*1, intent(in) :: small_n_or_c
      
      integer :: index1, index2, i
      real :: min, max, dist
      real, dimension(3) :: atomic_coords1, atomic_coords2
      
      valid_exp_constraints = .TRUE.

      do i=1, size(exp_constraints,1)

         ! WARNING
         ! I THINK THE BELOW IS RIGHT, BUT IM  NOT ENTIRELY SURE
         ! IM TOO TIRED TO MAKE SENSE OF IT
         ! Make sure that atomic indexes point at correct protein
         ! if CTER is smaller protein, leave contraint setup as is in input file
         if (small_n_or_c == "n") then 
            index1 = nint(exp_constraints(i,1))
            index2 = nint(exp_constraints(i,2))
         ! if NTER is smaller protein, swap indexes around, as atom coords are swapped elsewhere.  
         elseif (small_n_or_c == "c") then
            index2 = nint(exp_constraints(i,1))
            index1 = nint(exp_constraints(i,2))
         else
            write(ERROR_UNIT,*) "ERROR: Somethings funking in subroutine: within_exp_constraints in md_mod.f90"
         end if
         
         atomic_coords1 = coords1(index1,:)
         atomic_coords2 = coords2(index2,:)
         dist = distance(atomic_coords1,atomic_coords2)
         min = exp_constraints(i,3)-exp_constraints(i,4)
         max = exp_constraints(i,3)+exp_constraints(i,4)
         
         if ((dist .lt. min) .or. (dist .gt. max)) then 
            valid_exp_constraints = .FALSE.
            exit
         end if
         
      end do


    end subroutine within_exp_constraints


!-------------------------------------------------------------------------------
!   is_coords_unique:
!   Different combinations of angles can produce the same coordinates. is_coords_unique
!   checks current set of coords against an array of saved unique coordinate arrays 
!   to see if it mathces any. If no match, current array is unique, and output  
!   is TRUE. If there is a match, array is not unique, output is FALSE.
!   
!-------------------------------------------------------------------------------

    subroutine is_coords_unique(coords, unique_coords, n_unique_coords, is_unique)
      real, dimension(:,:), intent(in) :: coords
      real, dimension(:,:,:), intent(inout) :: unique_coords
      integer, intent(inout) :: n_unique_coords
      logical, intent(out) :: is_unique

      ! unique_coords(n_atoms,3,10)
      
      integer :: i, j, k, max_check
      real :: tol 
      logical, dimension(:), allocatable :: is_coords_set_unique

      is_unique = .FALSE.
      ! Because of stupid real maths I need to add a tolerance in here
      ! Values of xyz coordinates can be 'identical' but differ in the
      ! third decimal place by 1.. so 0.001. But I'm going to be a bit
      ! hamfisted and increased the tolerance a bit more, in case there
      ! are another other larger straglers.
      tol = 0.01

      ! Set n_unique_coords to 0 on first call
      if ( all(unique_coords(:,:,1)==0) ) then
         n_unique_coords = 0
      end if

      ! Set the number of set of coords to check over
      ! We should only store a certain number of coord arrays in
      ! unique_coords.. I'm worried about memory usage. This is
      ! defined in the program with an allocate statement. If
      ! the number of unique coordinates found is less than this
      ! we only need to check over this number, as checking over the
      ! empty (0) arrays will show everything up to size(unique_coords,3)
      ! is unique. And we can't have that. Can we now? 
      if ( n_unique_coords .lt. size(unique_coords,3) ) then
         max_check = n_unique_coords
      else
         max_check = size(unique_coords,3)
      end if

      ! is_coords_set_unique is an 1-D logical array of size max_check
      ! (number of arrays to check over). An entry is true if the current
      ! tested coords are unique against the ith set of coords in
      ! unique_coords. It enters False if current coords are identical
      ! to ith set of coords. All entries need to be TRUE for the current
      ! set of coords to be declared unique.
      allocate(is_coords_set_unique(max_check))
      
      ! If this is the first time, just add coords into unique_coords
      if (n_unique_coords == 0) then
         unique_coords(:,:,1) = coords
         n_unique_coords = n_unique_coords + 1
         is_unique = .TRUE.
      else
         ! Check over last 10 unique_coords
         do i=1, max_check
            is_coords_set_unique(i) = .FALSE.
            
            loop1: do j=1, size(coords,1)

               ! If current xyz match reference xyz, then cool.. go to next xyz
               ! until we run out of atomsto look at. If all match, then 
               ! is_coords_set_unique(i) remains false.
               if ( (abs(coords(j,1) - unique_coords(j,1,i)) .lt. tol) .and. &
                    (abs(coords(j,2) - unique_coords(j,2,i)) .lt. tol) .and. &
                    (abs(coords(j,3) - unique_coords(j,3,i)) .lt. tol) ) then 
                  continue
               else
                  ! If one of the current xyz coords do not match ref xyz coords,
                  ! then current coords are unique. Set is_coords_set_unique(i)
                  ! to true
                  is_coords_set_unique(i) = .TRUE.
                  exit loop1
               end if
            end do loop1
         end do

         ! Check over all the TRUE/FALSE entries for each saved array to see if new coords
         ! is unique

         ! If current coords are unique against last $max_check coords then set is_unique
         ! to true and update unique coords array
         if ( all(is_coords_set_unique .eqv. .TRUE.) ) then
            is_unique=.TRUE.
            n_unique_coords = n_unique_coords + 1

            ! Shift unique coords along 1 in unique_coords array
            do k=size(unique_coords,3)-1,1,-1
               unique_coords(:,:,k+1) = unique_coords(:,:,k)
            end do
            ! Add current coords into unique coords array
            unique_coords(:,:,1) = coords
         else
            ! If current coords matches at least one of previous coords,
            ! is_unique is false
            is_unique=.FALSE.
         end if
      end if

    end subroutine is_coords_unique

    
!-------------------------------------------------------------------------------
!   GENERATE_PROTEIN_COMBINATION: takes coordinates, performs translation and
!   rotation returns new coordinates. 
!
!
!-------------------------------------------------------------------------------
    subroutine generate_protein_coords(coords, resid, atomtype, translate_vector, rot_point, theta, phi, rrt_coords)
      ! Theta is rotation in x, phi is rotation in y
      real, dimension(:,:), intent(in) :: coords
      integer, dimension(:), intent(in) :: resid
      character(4), dimension(:), intent(in) :: atomtype
      real, dimension(3), intent(in) :: translate_vector, rot_point
      double precision, intent(in) :: theta, phi
      real, dimension(:,:), allocatable, intent(out) :: rrt_coords
      
      real, dimension(:,:), allocatable :: t_coords, rt_coords

      ! Assume protein being moved around is second protein, to which the linker is
      ! connected to the NTER. 
      call translate_coords(coords, translate_vector, t_coords)
      call rotate_coords(t_coords, rot_point, theta, 1, rt_coords)
      call rotate_coords(rt_coords, rot_point, phi, 2, rrt_coords)

    end subroutine generate_protein_coords


!-------------------------------------------------------------------------------
!   check_coords_are_valid: uses distance and overlap checks to see if new
!   coords are valid and suitable for output. Returns true if valid coordinates
!   returns false if not valid coordinates.     
!
!   Assumes coords are protein with linker connecting to NTER (second protein)
!   Assumes ref_coords are protein with linker connecting to CTER (first protein)    
!-------------------------------------------------------------------------------  

    subroutine check_coords_are_valid(coords, resid, atomtype, ref_coords, ref_resid, ref_atomtype,&
      exp_constraints, cutoff_dist, overlap_dist, small_n_or_c, big_n_or_c, is_coords_valid)
      ! Ref is large protein
      
      real, dimension(:,:), intent(in) :: coords, ref_coords, exp_constraints
      integer, dimension(:), intent(in) :: resid, ref_resid
      character(4), dimension(:), intent(in) :: atomtype, ref_atomtype
      logical, intent(out) :: is_coords_valid
      real, intent(in) :: cutoff_dist, overlap_dist
      character(1), intent(in) :: small_n_or_c, big_n_or_c
      
      real, dimension(3) :: small_TER_coords, big_TER_coords
      logical :: within_cutoff, no_overlap, valid_exp_constraints

      ! Have to demonstrate new coordinates are valid
      is_coords_valid = .FALSE.
     
      ! Get link
      call get_TER_coords(coords, resid, atomtype, small_n_or_c, small_TER_coords)
      call get_TER_coords(ref_coords, ref_resid, ref_atomtype, big_n_or_c, big_TER_coords)

      ! Only print file if distance between linker bases are within cutoff
      ! and if there is no overlap of the two proteins

      ! Run checks
      call within_cutoff_distance(small_TER_coords, big_TER_coords, cutoff_dist, within_cutoff)     
      call within_exp_constraints(exp_constraints, coords, ref_coords, small_n_or_c, valid_exp_constraints)
      call is_no_overlap(coords, ref_coords, overlap_dist, no_overlap)


      ! See if validation critera is met. If met, return True. 
      if ((within_cutoff .eqv. .TRUE.) .and. (no_overlap .eqv. .TRUE.)&
           .and. (valid_exp_constraints .eqv. .TRUE.)) then 
         is_coords_valid = .TRUE.
      end if

      end subroutine check_coords_are_valid


!-------------------------------------------------------------------------------
!   generate_output_filename: takes translation vector and angles 
!   to create a unique name for the output pqr file. 
!-------------------------------------------------------------------------------
    subroutine generate_output_filename(trans_vector, theta, phi, gamma, output_filename)
      real, dimension(3), intent(in) :: trans_vector
      double precision, intent(in) :: theta, phi, gamma
      character*80, intent(out) :: output_filename

      character*80 :: ctv1, ctv2, ctv3, ctheta, cphi, cgamma, cdist
      real :: theta_deg, phi_deg, gamma_deg, dist
      
      ! Round trans_vector, theta and phi
      ! to one decimal place
      ! r = rounded

      theta_deg = (theta * 180)/PI
      phi_deg = (phi * 180)/PI
      gamma_deg = (gamma * 180)/PI

      dist = mag(trans_vector)
      
      write(cdist,'(F6.1)') dist
      write(ctv1,'(F6.1)') trans_vector(1)
      write(ctv2,'(F6.1)') trans_vector(2)
      write(ctv3,'(F6.1)') trans_vector(3) 
      write(ctheta,'(F6.1)') theta_deg
      write(cphi,'(F6.1)')  phi_deg
      write(cgamma,'(F6.1)')  gamma_deg 
    
      ! Combine to form output filename with .pqr extension
      output_filename =  trim(adjustl(cdist)) // '_' // &
           trim(adjustl(ctv1)) // '_' // trim(adjustl(ctv2)) // '_' // trim(adjustl(ctv3)) //&
           '_' // trim(adjustl(ctheta)) // '_' // trim(adjustl(cphi)) &
           // '_' // trim(adjustl(cphi)) //'.pdb'

    end subroutine generate_output_filename


!-------------------------------------------------------------------------------
!   generate_logfile_filename: takes pqr file names, smooshes them together
!   makes a new filename
!-------------------------------------------------------------------------------
    subroutine generate_log_filename(pqr_filename1, pqr_filename2, log_filename)
      character*80, intent(in) :: pqr_filename1, pqr_filename2
      character*80, intent(out) :: log_filename

      character*80 :: prefix1, prefix2 
      integer :: pos1, pos2

      ! Get position of full stop of each pqr file
      pos1 = scan(trim(pqr_filename1),".", BACK= .true.)
      pos2 = scan(trim(pqr_filename2),".", BACK= .true.)

      ! Define prefix of filename as everything before full stop
      if ( pos1 > 0 ) prefix1 = pqr_filename1(1:pos1-1)
      if ( pos2 > 0 ) prefix2 = pqr_filename2(1:pos2-1)
      
      ! Combine prefixes to form output filename with .log extension
      log_filename = trim(adjustl(prefix1)) // '_' // trim(adjustl(prefix2)) // '.log'
      
    end subroutine generate_log_filename


!-------------------------------------------------------------------------------
!   write_valid_log_file: write stats about conformations found to be valid
!    
!-------------------------------------------------------------------------------
    subroutine write_log_file_header(log_unit, conf_filename, pqr_filename1, pqr_filename2, cutoff_dist, &
         grid_spacing, theta, phi, gamma, overlap_dist, exp_constraints)

      integer, intent(in) :: log_unit
      character*80, intent(in) :: conf_filename, pqr_filename1, pqr_filename2
      real, intent(in) :: cutoff_dist, overlap_dist, grid_spacing      
      double precision, intent(in) :: theta, phi, gamma
      real, dimension(:,:), intent(in) :: exp_constraints
      
      double precision :: theta_deg, phi_deg, gamma_deg
      integer :: i

      theta_deg = (theta * 180)/PI
      phi_deg = (phi * 180)/PI
      gamma_deg = (gamma * 180)/PI
      
      write(log_unit,*) "######################### BloBerator #########################"
      write(log_unit,*) "Written by Michael Thomas for Macka"
      write(log_unit,*) "June 2016"
      write(log_unit,*) "Australian National University"
      write(log_unit,*) "Questions?: email michael.thomas@anu.edu.au"
      write(log_unit,*) ""
      write(log_unit,*) ""
      write(log_unit,*) "###################### Input Paramaters ######################"
      write(log_unit,*) "Input file: ", trim(adjustl(conf_filename))
      write(log_unit,*) "Protein with CTER linker: ", trim(adjustl(pqr_filename1))
      write(log_unit,*) "Protein with NTER linker: ",trim(adjustl(pqr_filename2))
      write(log_unit,*) "Linker length:", cutoff_dist, "Angstroms"
      write(log_unit,*) "Grid spacing:", grid_spacing, "Angstroms"
      write(log_unit,*) "Angle of rotation about x axis: ", theta_deg, "degrees"
      write(log_unit,*) "Angle of rotation about y axis: ", phi_deg, "degrees"
      write(log_unit,*) "Angle of rotation about z axis: ", gamma_deg, "degrees"
      write(log_unit,*) "vDW distance for overlap criteria: ", overlap_dist, "Angstroms"
      write(log_unit,*) "Number of experimental constraints: ", size(exp_constraints,1)
      if (size(exp_constraints,1) .gt. 0) then
         write(log_unit,*) "Atom index 1     Atom index 2   Constraint Distance   Constraint Tolerance"
         do i=1, size(exp_constraints,1)
            write(log_unit,*) exp_constraints(i,:)
         end do
      end if
      write(log_unit,*) ""
      write(log_unit,*) ""
      write(log_unit,*) "################ Valid Protein Conformations ################"
      write(log_unit,*) "Distance (A)       Theta (degrees)           Phi (degre&
           es)             Gamma(degrees)          Coordinate Filename"
    end subroutine write_log_file_header
    
    
    
!-------------------------------------------------------------------------------
!   write_valid_log_file: write stats about conformations found to be valid
!    
!-------------------------------------------------------------------------------
    subroutine write_valid_to_log_file(log_unit, dist, theta, phi, gamma, output_pqr_filename)

      integer, intent(in) :: log_unit
      real, intent(in) :: dist
      double precision, intent(in) :: theta, phi, gamma 
      character*80, intent(in) :: output_pqr_filename
      
      double precision :: theta_deg, phi_deg, gamma_deg
      
      theta_deg = (theta * 180)/PI
      phi_deg = (phi * 180)/PI
      gamma_deg = (gamma * 180)/PI
      
      write(log_unit,*) dist, theta_deg, phi_deg, gamma_deg, trim(adjustl(output_pqr_filename))
      
    end subroutine write_valid_to_log_file

!-------------------------------------------------------------------------------
!   write_log_file_summary: print some stats to end of log file
!    
!-------------------------------------------------------------------------------
    subroutine write_log_file_summary(log_unit, n_valid_configs, n_total_configs)

      integer, intent(in) :: n_valid_configs, n_total_configs, log_unit
      character*80 :: c_nvc, c_ntc, c_per

      write(c_nvc,*) n_valid_configs
      write(c_ntc,*) n_total_configs
      write(c_per,*) 100*(real(n_valid_configs)/real(n_total_configs))

      write(log_unit,*) trim(adjustl(c_nvc)), " valid protein configurations"
      write(log_unit,*) trim(adjustl(c_ntc)), " total configurations tested"
      write(log_unit,*) trim(adjustl(c_per)), "%"
      
    end subroutine write_log_file_summary


end module bloberator_mod
