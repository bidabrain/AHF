!Simple program to read all gadget files in a directory and print out their header information

! Compile with gfortran -cpp gadgetheader.f90 -o gadgetheader -fdefault-integer-8

! Usage ./gadgetheader /eg/dir/ramses2gadget_181.

program gadgetheader

	implicit none

	logical :: verbose = .true.
	logical :: ramses2gadget = .true.

	integer*4 :: nfiles, ifile
	character(len=128) :: suffix
	character(len=128) :: fname, target

	! Gadget Header Block:

	integer*4 :: header_length_ghead      ! header begins with an integer sizeof(header)
	integer(kind=4) :: npart_ghead(6)           ! number of particles of each type in this file
	real*8 :: mass_ghead(6)               ! mass of particles of each type
	real*8 :: time_ghead                  ! expansion factor
	real*8 :: redshift_ghead              ! redshift of snapshot file
	integer*4 :: flag_sfr_ghead           ! flags whether the simulation was including star formation 
	integer*4 :: flag_feedback_ghead      ! flags whether feedback was included (obsolete) 
	integer*4 :: npartTotal_ghead(6)      ! total number of particles of each type in this snapshot. 
	integer*4 :: npartTotal_ghead_8byte(6)
	integer*4 :: flag_cooling_ghead       ! flags whether cooling was included  
	integer*4 :: num_files_ghead          ! number of files in multi-file snapshot 
	real*8 :: BoxSize_ghead               ! box-size of simulation in case periodic boundaries were used. Could be kpc/h or Mpc/h set with gadget_length_flag
	real*8 :: Omega0_ghead                ! matter density in units of critical density 
	real*8 :: OmegaLambda_ghead           ! cosmological constant parameter 
	real*8 :: HubbleParam_ghead           ! Hubble parameter in units of 100 km/sec/Mpc 
	integer*4 :: flag_stellarage_ghead    ! flags whether the file contains formation times of star particles 
	integer*4 :: flag_metals_ghead                ! flags whether the file contains metallicity values for gas and star particles
	integer*4 :: npartTotalHighWord_ghead(6)      ! High word of the total number of particles of each type !!!! CHECK THIS - unsigned!!!!!!
	integer*8 :: npartTotalHighWord_ghead_8byte(6)
	integer*4 :: flag_entropy_instead_u_ghead     ! flags that IC-file contains entropy instead of u 
	character(60) :: fill_ghead           ! fill rest of header to 256 bytes
	integer*4 :: part_data_size

	! Read program arguments
	call read_args(fname)

	do ifile = 0, 63

	! Loop over all the files
	!do ifile = 0, nfiles

		! Write filenumber to string
		write(suffix, '(i2)') ifile

		suffix = adjustl(suffix)

		! hardcode prefix for now
		target = trim(fname)//trim(suffix)
		!fname = trim('62.5_dm_020.')//trim(suffix)

		! Open the file    TODO PERFORM CHECKS
		call open_gadget_file(target)

		! Read the header
		call read_gadget_header

		! Some verbose logging
		if (verbose) write(*,*) 'Reading header for file ', target

		! Print the header to the console
		call print_gadget_header

	enddo

	contains

		subroutine open_gadget_file(fname)

			character(len=128) fname

			open (unit=10,file=fname,access='stream',status='old')

		end subroutine open_gadget_file

		subroutine read_gadget_header

			integer*4 :: i, int_4byte
			integer*8 :: int_8byte
			character(4) :: dummy_char

			if (ramses2gadget) then

				! Use dummy_char to match up with ramses2gadget.f90
				read(10) dummy_char,dummy_char,header_length_ghead,dummy_char,dummy_char,npart_ghead,mass_ghead,time_ghead,&
					&redshift_ghead,flag_sfr_ghead,flag_feedback_ghead,npartTotal_ghead,&
					&flag_cooling_ghead,num_files_ghead,BoxSize_ghead,Omega0_ghead,OmegaLambda_ghead,&
					&HubbleParam_ghead,flag_stellarage_ghead,flag_metals_ghead,&
					&npartTotalHighWord_ghead,flag_entropy_instead_u_ghead, fill_ghead

			else

				! Read as done in GADGET2CUBEP3M
				read(10) header_length_ghead,npart_ghead,mass_ghead,time_ghead,&
			        &redshift_ghead,flag_sfr_ghead,flag_feedback_ghead,npartTotal_ghead,&
			        &flag_cooling_ghead,num_files_ghead,BoxSize_ghead,Omega0_ghead,OmegaLambda_ghead,&
			        &HubbleParam_ghead,flag_stellarage_ghead,flag_metals_ghead,&
			        &npartTotalHighWord_ghead,flag_entropy_instead_u_ghead, fill_ghead

			endif

			  ! convert unsigned integers into 8 byte fortran integers:

			do i = 1, 6

				int_4byte = npartTotal_ghead(i)
				int_8byte = int(int_4byte,8)

				if(int_8byte .lt. 0) int_8byte = 4294967296 - abs(int_8byte)

				npartTotal_ghead_8byte(i) = int_8byte

			enddo

			do i = 1, 6

				int_4byte = npartTotalHighWord_ghead(i)

				int_8byte = int(int_4byte,8)

				if(int_8byte .lt. 0) int_8byte = 4294967296 - abs(int_8byte)

				npartTotalHighWord_ghead_8byte(i) = int_8byte

			enddo

		end subroutine read_gadget_header

		subroutine print_gadget_header

			print*, 'Gadget Header Size =    ', header_length_ghead
			print*, 'Gas =                   ', npart_ghead(1)
			print*, 'DM =                    ', npart_ghead(2)  
			print*, 'Disk =                  ', npart_ghead(3)
			print*, 'Bulge =                 ', npart_ghead(4)  
			print*, 'Stars =                 ', npart_ghead(5)
			print*, 'Bndry =                 ', npart_ghead(6)  
			print*
			print*, 'Particle Masses (in Gadget internal units):'
			print*, 'Gas =                   ', mass_ghead(1)
			print*, 'DM =                    ', mass_ghead(2)  
			print*, 'Disk =                  ', mass_ghead(3)
			print*, 'Bulge =                 ', mass_ghead(4)  
			print*, 'Stars =                 ', mass_ghead(5)
			print*, 'Bndry =                 ', mass_ghead(6)  
			print*
			print*, 'Expansion Factor:       ', time_ghead
			print*, 'Redshift:               ', redshift_ghead
			print*, 'Omega M0:               ', Omega0_ghead
			print*, 'Omega Lambda0:          ', OmegaLambda_ghead
			print*, 'Boxsize (Mpc/h):        ', BoxSize_ghead
			print*, 'H_0 (km/s):             ', HubbleParam_ghead
			print* 
			print*, 'Number of files for ICs:', num_files_ghead
			print* 
			print*, 'Total Particles in Simulation'
			print*, 'Gas =                   ', npartTotal_ghead_8byte(1)
			print*, 'DM =                    ', npartTotal_ghead_8byte(2)  
			print*, 'Disk =                  ', npartTotal_ghead_8byte(3)
			print*, 'Bulge =                 ', npartTotal_ghead_8byte(4)  
			print*, 'Stars =                 ', npartTotal_ghead_8byte(5)
			print*, 'Bndry =                 ', npartTotal_ghead_8byte(6)  
			print*
			print*, 'Star Formation Flag:    ', flag_sfr_ghead
			print*, 'Feedback Flag:          ', flag_feedback_ghead
			print*, 'Cooling Flag:           ', flag_cooling_ghead
			print*, 'Steller Age Flag:       ', flag_stellarage_ghead
			print*, 'Flag Metals:            ', flag_metals_ghead
			print*, 'Entropy Flag:           ', flag_entropy_instead_u_ghead
			print* 
			print*, 'Total Particles High Word'
			print*, 'Gas =                   ', npartTotalHighWord_ghead_8byte(1)
			print*, 'DM =                    ', npartTotalHighWord_ghead_8byte(2)  
			print*, 'Disk =                  ', npartTotalHighWord_ghead_8byte(3)
			print*, 'Bulge =                 ', npartTotalHighWord_ghead_8byte(4)  
			print*, 'Stars =                 ', npartTotalHighWord_ghead_8byte(5)
			print*, 'Bndry =                 ', npartTotalHighWord_ghead_8byte(6)  
			print*

		end subroutine print_gadget_header

end program gadgetheader

! Subroutine to read command line arguments and set variables
subroutine read_args(fname)

	implicit none

	character(len=128)				:: arg1
	character(len=128), intent(out) :: fname

	call getarg(1, arg1)
	fname = trim(arg1)

	return

end subroutine read_args