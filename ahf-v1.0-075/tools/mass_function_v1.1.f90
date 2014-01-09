	
	! Produces mass functions based on spherical overdensity (SOD) and
	! friend-of-friends (FOF) halo finders and also produces 
	! Press-Schechter and Sheth-Tormen curves for comparison
  ! v1.1 added cubep3m PID support
	
	! COMPILE WITH: 
	
	! (willmb) gfortran -cpp mass_function_v1.1.f90 ST_PS_files/deltac.f ST_PS_files/sorthalf.f ST_PS_files/spline.f ST_PS_files/splint.f ST_PS_files/growth.f -o mf -DGFORTRAN
	! (ranger) ifort -fpp mass_function_v1.1.f90 ST_PS_files/sorthalf.f ST_PS_files/deltac.f ST_PS_files/spline.f ST_PS_files/splint.f ST_PS_files/growth.f -o mf -DBINARY

	! Compilation options: -DBINARY (for intel compilers)
	!                      -DGFORTRAN (for gnu fortran compiler)


 
	implicit none
  
  logical :: analytic_only = .true. ! set to true for just analytic predictions for this z and boxsize
  integer :: min_part_cutoff = 20 ! particle cut-off for mass function data
  integer :: ahf_version = 1 ! 0 for old, 1 for new
  integer :: ahf_substructure = 0 ! 0 for exclude, 1 for include (only valid with new version)
  integer :: ahf_overdensity = 0 ! 0 for not included in file string, otherwise choose which overdensity the halos have
  integer :: volume_correction = 0 ! 0 for no, 1 for yes (no need for boxes over 1Gpc)

  ! Simulation details

  real*8 :: z = 0.000 ! redshift of data - can generalise this to an input table in the future
  real :: box_size = 1000
  integer :: nc_sim = 6912
  integer :: nodes_dim_sim = 12
  integer :: cpm_version = 1 ! 0 = Old version, 1 = with PIDs, 2 = new version no PIDs
  integer*4 :: num_files_ahf_sim = 0 ! Number of files in the ahf output  
  integer*4 :: num_files_fof_chunk_sim = 0 ! Number of files in the fof output  
  
  ! Cosmology flag
    
  integer :: cosmology = 2! add different cosmologies here currently we have 0 and 1 (both WMAP 5), and 2 (planck)
  
	real*8 :: omegabh2,omega0,lambda0,h,an,tcmb,s8, omega_z
  

  ! where are the halo catalogues?

	character(*), parameter :: input_path ='/scratch/01454/ww60/cubepm_100909_12_3456_1Gpc/results/'


  ! where should we write the datafiles?

	character(*), parameter :: output_path ='./'


  character(len=512) :: cpm_file,ahf_file,fof_file,fofcor_file,ahf_dfile,cpm_dfile,fof_dfile,fofcor_dfile,afile

  ! do not alter these numbers!

	real*8 :: mass_bin_fac_coarse = 9.399 ! factor for scaling coarse mass bins
	real*8 ::	m_min = 76400 ! this gives a 400 particle lower end to a mass bin in the 20Mpc, 114 Mpc and 425 Mpc boxes
	integer,parameter :: bin_dim = 59
	integer,parameter :: dim = 100
	integer,parameter :: bin_dim_ps=bin_dim*10
	integer,parameter :: dim_ps=dim*10
  
  real :: box
  real*8 :: sigma_400_part, sigma_box, sigma_test
  integer*4 :: num_files_ahf, num_files_fof_chunk, sigma_limit

! Unit conversion variables
	
	real*8,parameter :: M_sun = 1.989e+33, cm_to_Mpc=3.086e+24
	real*8,parameter :: pi=3.1415927, grav_const=6.672e-8 !cm^3/g/s^2
	
! General variables

	integer :: i, j, k, m,nn(3), ll, nbf, nbo, lob, simulation, nc, nodes_dim
	real*8 :: omegab, deltac, growth, mass8, sigma8, length_scale
	real*8 :: rho_crit_0, rho_crit_z, rho_crit_0_MpcM_sun, rho_sub
	real*8 :: rho_bar, deltac_z, nu, nu_old,a, m_ave

! File variables 

	character*456 line ! header of 1st halo file
	character*200 filename, file, redshift, box_string, part_num_string, overdensity_string
	character(len=512) :: ifile
	
	integer*4 :: file_end_int
	integer*4 total_part
	
	character(512) :: nfile1, mfile1

! cubepm halo data variables
	integer, parameter :: nhv=17 ! number of halo quantities in each catalog entry
	integer*4 :: n_halo
	real(4), dimension(22) :: halo_input_buffer
  integer, parameter :: N_p = 50
  integer(8), dimension(N_p) :: pid_halo
  real(4), dimension(6,N_p) :: xv_halo
  real(4), dimension(6) :: I_ij  
	
  
  
	! Histogram data: number of bin, number of halos in bin fof, number of halos
	! in bin sod, bin lower mass, bin upper mass, mid-point bin mass
	
	real*8, dimension (dim) :: data_mass_ave, data_mass_high, data_mass_low
	real*8, dimension (dim) :: data_mass_cubep3m,data_mass_fof,data_mass_ahf,data_mass_fof_cor
	real*8, dimension (dim) :: data_cubep3m, data_fof_chunk, data_fof_chunk_cor, data_ahf
	real*8, dimension (dim) :: data_ps, data_st, data_jenkins, data_warren, data_reed_2003, data_reed_2007, data_tinker
	real*8, dimension (dim) :: total_ahf_chunk, total_cubep3m, total_fof_chunk,total_fof_chunk_cor ! total mass of halos in a given bin - used for averaging mass
  real*8, dimension (dim) :: error_cubep3m_low, error_fof_chunk_low, error_fof_chunk_cor_low, error_ahf_low
  real*8, dimension (dim) :: error_cubep3m_high, error_fof_chunk_high, error_fof_chunk_cor_high, error_ahf_high
  real*8, dimension (dim_ps) :: data_ps_del,data_st_del,data_jen_del,data_war_del,data_tin_del,data_r03_del,data_r07_del
  real*8, dimension (dim_ps) :: data_ang_del,data_cro_del,data_cou_del
  real*8, dimension (dim) :: data_ps_int,data_st_int,data_jen_int,data_war_int,data_tin_int,data_r03_int,data_r07_int
  real*8, dimension (dim) :: data_ang_int,data_cro_int,data_cou_int
	

  
! cubepm particle data variables

	integer*4 :: np_local

	real*8 :: halo_mass

! AHF file variables

	integer*4 :: ahf_halos_in_file
	integer*4 :: ahf_total_halos
  real*4 :: ahf_dummy
  real*4, dimension(3) :: ahf_pos

	integer(8) :: npart, nbins, ID, HostHalo, NumSubStruct, ID_count
  real(8) :: pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir
  real(8) :: R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE
  real(8) :: Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx, M_vir
  real(8) :: Ecy,Ecz,ovdens,fMhires,Ekin,Epot,SurfP,Phi0,cNFW

! FOF_chunk file variables
	
	integer*4	:: Ngroups_chunk, Ngroups_old_chunk, TotNgroups_chunk, N_FOF_CHUNK

	integer*4 :: Length_chunk
	real*4 :: Length_chunk_corrected
	integer*4  Offset_chunk
	real*4 :: Group_Mass_chunk
	real*4, dimension(3) :: CM_chunk, vel_chunk

! PS and ST variables

	real*8 :: m_grid, m_box, mass1
	real*4 :: dummy
	real*8 :: n_fine(dim_ps+1), mms_fine(0:dim_ps+1), mms_fine_ave(0:dim_ps+1)
	real*8 :: slope(dim_ps+1),slope_z(dim_ps+1),sigma(dim_ps),sigma_z(dim_ps)
	real*8, allocatable :: sigma_local(:),mass_local_del(:),sigma_local_del(:)
	real*8, allocatable :: mms_fine_sorted_by_sigma(:),sigma_local_sorted_by_mass(:)
	real*8 :: sigma_bin_ave_del(dim_ps),sigma_bin_ave_del_z(dim_ps)
	real*8 :: sigma_bin_ave(dim),sigma_bin_ave_z(dim)
  real*8 :: sigma_mass_ave(dim), sigma_mass_low(dim), sigma_mass_high(dim)
  real*8 :: slope_coarse(dim+1),slope_z_coarse(dim+1),sigma_z_coarse(dim+1), sigma_coarse(dim+1)
	real*8 :: n(dim+1), n_simul(dim+1)
  integer,parameter :: max_size = 1000
  real*8 :: mass(max_size), sigma1(max_size), sigma_cobe(max_size), mass_cobe(max_size)
  real*8 :: mass_2(max_size), sigma1_2(max_size)
	real*8 :: table_flat(1000,2),table_open(1000,2)
	real*8,parameter ::  alphapar=0.707, ppar=0.3,  Apar=0.322
	real*8 :: n_ST_fine(dim_ps+1), n_ST(dim+1) 
	real*8 :: f_ST(dim_ps),f_PS(dim_ps)
	
! Jenkins variables
	
	real*8 :: f_jenkins(dim_ps)
	real*8 :: n_jenkins(dim_ps)
	real*8 :: f_ps_st(dim_ps)

! Reed_2003 variables

	real*8 :: f_reed_2003(dim_ps)
	real*8 :: n_reed_2003(dim_ps)

! Reed 2007 variables

	real*8 :: f_reed_2007(dim_ps)
	real*8 :: n_reed_2007(dim_ps)
	real*8,parameter :: c_par = 1.08
	real*8 :: n_eff(dim_ps), G_1(dim_ps), G_2(dim_ps)
	
! Warren variables

	real*8 :: f_warren(dim_ps)
	real*8 :: n_warren(dim_ps)

! Tinker variables

	real*8 :: c_tink
	real*8 :: big_a_tink
	real*8 :: small_a_tink
	real*8 :: b_tink
	real*8 :: alpha_tink
	
	real*8 :: f_tinker(dim_ps)
	real*8 :: n_tinker(dim_ps)
  
! Angulo variables

	real*8 :: f_angulo(dim_ps)  
	real*8 :: n_angulo(dim_ps) 

	real*8, dimension (dim) :: data_angulo
  
  
! Crocce Variables

	real*8 :: f_crocce(dim_ps)  
	real*8 :: n_crocce(dim_ps) 

	real*8, dimension (dim) :: data_crocce
  
  real*8 :: large_a_crocce,a_crocce, b_crocce,c_crocce

!courtin variables

  real*8 :: f_courtin(dim_ps)
  real*8 :: n_courtin(dim_ps)  

	real*8, dimension (dim) :: data_courtin
  
  real*8 :: large_a_courtin = 0.348
  real*8 :: small_a_courtin = 0.695
  real*8 :: p_courtin = 0.1
  
  
! Watson variables

  real*8 :: f_watson_fof_uni(dim_ps)
  real*8 :: n_watson_fof_uni(dim_ps)  

	real*8, dimension (dim) :: data_watson_fof_uni
  
  real*8 :: a_watson_fof_uni = 0.194!0.282244!0.19544
  real*8 :: b_watson_fof_uni = 2.267!1.40631!1.98919
  real*8 :: c_watson_fof_uni = 1.805!2.163!1.74004
  real*8 :: d_watson_fof_uni = 1.287!1.2096!1.15987
  
 
 
  real*8 :: f_watson_cpm_ahf(dim_ps)
  real*8 :: n_watson_cpm_ahf(dim_ps)  

	real*8, dimension (dim) :: data_watson_cpm_ahf
  
  real*8 :: a_watson_cpm_ahf 
  real*8 :: b_watson_cpm_ahf
  real*8 :: c_watson_cpm_ahf
  real*8 :: d_watson_cpm_ahf 
  

  real*8, dimension (dim_ps) :: data_watson_fof_uni_del, data_watson_cpm_ahf_del
  real*8, dimension (dim) ::  data_watson_fof_uni_int, data_watson_cpm_ahf_int

        
  common /interp/ table_flat, nbo, nbf

!------------------------------------------------------------------------------------------------------------------------!

! configure cosmology

  if(cosmology .eq. 0) then
  
  omegabh2=0.02156
  omega0=0.27
  lambda0=0.73
  h=0.7
  an=0.96d0
  tcmb=2.726
  s8=0.8
  
  endif
  
  if(cosmology .eq. 1) then
  
  omegabh2=0.02156
  omega0=0.279
  lambda0=0.721
  h=0.701
  an=0.96d0
  tcmb=2.726
  s8=0.807
  
  endif

  if(cosmology .eq. 2) then
  
  omegabh2=0.022140
  omega0= 0.307115
  lambda0=1.0-omega0
  h=0.6777
  an=0.9611d0
  tcmb=2.725
  s8=0.8288
  
  endif

  omega_z = omega0/(omega0+(1.0/(1+z))**3.0*(1-omega0))


! configure output files:

  write(redshift,'(f8.3)') z
  redshift=adjustl(redshift)

  write(box_string,'(i4)') int(box_size)
  box_string=adjustl(box_string)

  write(part_num_string,'(i6)') int(nc_sim / 2.0)
  part_num_string=adjustl(part_num_string)

  print*, "Redshift: ",redshift(1:len_trim(redshift))
  print*, "Box size (Mpc/h): ",box_string(1:len_trim(box_string))  
  print*, "Particles ^ 1/3: ",part_num_string(1:len_trim(part_num_string))    

  cpm_file = output_path//redshift(1:len_trim(redshift))//'_cpm_mf_'//&
            &part_num_string(1:len_trim(part_num_string))//'_'//&
            &box_string(1:len_trim(box_string))//'Mpc.dat'

  ahf_file = output_path//redshift(1:len_trim(redshift))//'_ahf_mf_'//&
            &part_num_string(1:len_trim(part_num_string))//'_'//&
            &box_string(1:len_trim(box_string))//'Mpc.dat'

  fof_file = output_path//redshift(1:len_trim(redshift))//'_fof_mf_'//&
            &part_num_string(1:len_trim(part_num_string))//'_'//&
            &box_string(1:len_trim(box_string))//'Mpc.dat'

  fofcor_file = output_path//redshift(1:len_trim(redshift))//'_fofcor_mf_'//&
            &part_num_string(1:len_trim(part_num_string))//'_'//&
            &box_string(1:len_trim(box_string))//'Mpc.dat'

  ahf_dfile = output_path//redshift(1:len_trim(redshift))//'_ahf_residuals_'//&
            &part_num_string(1:len_trim(part_num_string))//'_'//&
            &box_string(1:len_trim(box_string))//'Mpc.dat'

  cpm_dfile = output_path//redshift(1:len_trim(redshift))//'_cpm_residuals_'//&
            &part_num_string(1:len_trim(part_num_string))//'_'//&
            &box_string(1:len_trim(box_string))//'Mpc.dat'

  fof_dfile = output_path//redshift(1:len_trim(redshift))//'_fof_residuals_'//&
            &part_num_string(1:len_trim(part_num_string))//'_'//&
            &box_string(1:len_trim(box_string))//'Mpc.dat'

  fofcor_dfile = output_path//redshift(1:len_trim(redshift))//'_fof_cor_residuals_'//&
            &part_num_string(1:len_trim(part_num_string))//'_'//&
            &box_string(1:len_trim(box_string))//'Mpc.dat'

	  afile = output_path//redshift(1:len_trim(redshift))//'_analytic_mass_function_'//&
	    &part_num_string(1:len_trim(part_num_string))//'_'//&
               &box_string(1:len_trim(box_string))//'Mpc.dat'



  print*, cpm_file(1:len_trim(cpm_file))
  print*, ahf_file(1:len_trim(ahf_file))
  print*, fof_file(1:len_trim(fof_file))
  print*, fofcor_file(1:len_trim(fofcor_file))  
  print*, afile(1:len_trim(afile))
  print*, cpm_dfile(1:len_trim(cpm_dfile))
  print*, ahf_dfile(1:len_trim(ahf_dfile))    
  print*, fof_dfile(1:len_trim(fof_dfile))
  print*, fofcor_dfile(1:len_trim(fofcor_dfile))  

! Calculate cosmological parameters

  nc = nc_sim
  nodes_dim = nodes_dim_sim
  box = real(box_size)
  num_files_ahf = num_files_ahf_sim
  num_files_fof_chunk = num_files_fof_chunk_sim
  
	omegab = omegabh2/h**2
	rho_crit_0=1.88d-29*h**2!g/cm^3
	rho_crit_0_MpcM_sun=rho_crit_0*cm_to_Mpc**3/M_sun
	rho_bar = rho_crit_0_MpcM_sun*omega0
  
	m_box = rho_bar*(box/h)**3
	m_grid = m_box/real(nc)**3

	print*,'m_box=',m_box
	print*,'m_grid=',m_grid
	print*,'m_min = ',m_min
	print*,'rho_bar = ',rho_bar

  total_ahf_chunk	= 0
  total_fof_chunk	= 0
  total_cubep3m	= 0            
		
  ! set bins for data:
		
	do i=1,dim
			
		data_fof_chunk(i) = 0    
		data_fof_chunk_cor(i) = 0    
		data_cubep3m(i) = 0
		data_ahf(i) = 0
		data_ps(i) = 0	
		data_st(i) = 0				
		data_jenkins(i) = 0
		data_warren(i) = 0
		data_reed_2003(i) = 0
		data_reed_2007(i) = 0	
		data_tinker(i) = 0				
		data_angulo(i) = 0	    		
		data_crocce(i) = 0	
		data_courtin(i) = 0	
		data_watson_fof_uni(i) = 0	
												
		data_mass_high(i) = m_min*10.**(mass_bin_fac_coarse*(real(i))/bin_dim) !mass bins for plotting
			
		if(i==1) then
			data_mass_low(i) = m_min-1
		else
			data_mass_low(i) = data_mass_high(i-1)		
		endif
		
		data_mass_ave(i) = data_mass_low(i) + (data_mass_high(i) - data_mass_low(i))/2
    
		
	enddo


!---------------------------------------------------------------------------------------------

! Now produce Press-Schechter and Sheth-Tormen approximations
	
	! Read interpolation table
  
200	open(unit=1,file='ST_PS_files/deltac_flat.dat',status='old')
	nbf=0
	do i=1,1000
		read(1,*,end=6) dummy, table_flat(i,1), table_flat(i,2)
		nbf=nbf+1
	enddo
6	close(unit=1)
  
  if(cosmology .eq. 0)	open(20,file='ST_PS_files/ms.z0_wmap3plus_high_a') ! contains mass scale vs. sigma
  if(cosmology .eq. 1)	open(20,file='ST_PS_files/ms.z0_wmap3plus_high_b') ! contains mass scale vs. sigma
  if(cosmology .eq. 2)  open(20,file='ST_PS_files/planck_ms.dat')
	do i=1,max_size
		read(20,*) mass(i),sigma1(i)
		mass(i)=mass(i)/h ! since mass is in 1/h units in sigma file from cmbfast
	end do

	close(20)
	

	if(cosmology .eq. 0) open(30,file='ST_PS_files/ms.z0_wmap3plus_high_a_2') ! contains mass scale vs. sigma
	if(cosmology .eq. 1) open(30,file='ST_PS_files/ms.z0_wmap3plus_high_b_2') ! contains mass scale vs. sigma
        if(cosmology .eq. 2)  open(30,file='ST_PS_files/planck_ms.dat_2')

	do i=1,max_size
		read(30,*) sigma1_2(i),mass_2(i)
		mass_2(i)=mass_2(i)/h ! since mass is in 1/h units in sigma file from cmbfast
	end do

	close(30)

	call spline(mass,sigma1,max_size,1d31,1d31,sigma_cobe)
	call spline(sigma1_2,mass_2,max_size,1d31,1d31,mass_cobe)  

	! Check the normalization of sigma - Normalization length scale in Mpc.

	length_scale=8/h
	mass8 = 4.*pi*omega0*h**2/3*(length_scale/1.533e-04)**3
	call splint(mass,sigma1,sigma_cobe,max_size,mass8,sigma8)

	sigma1=sigma1*s8/sigma8

	print*,'check, sigma8=',sigma8,'mass8=',mass8,'length=',length_scale

	call splint(mass,sigma1,sigma_cobe,max_size,mass8,sigma8)

	! mass bins for PS integrals

  mms_fine = 0

print*	

	do i=1,dim_ps
    
    
		mms_fine(i)=m_min*10.**(mass_bin_fac_coarse*(real(i-1))/bin_dim_ps)
  
  end do
    
	do i=1,dim_ps
      
      if(i .eq. 1) mms_fine_ave(i) = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) mms_fine_ave(i) = (mms_fine(i)+mms_fine(i-1))/2.
      
  end do    

	do i=1,dim_ps      
    call splint(mass,sigma1,sigma_cobe,max_size,mms_fine_ave(i),sigma(i))
  end do  
  
  
  do i=1,dim+1
		call splint(mass,sigma1,sigma_cobe,max_size,data_mass_ave(i),sigma_coarse(i))
	end do

  
  ! calculate overdensity for this redshift:
  
  deltac_z = deltac(omega0,lambda0,z)

  slope = 0
  slope_z = 0
  slope_coarse = 0
  slope_z_coarse = 0

  


	do i=2,dim_ps
	if(sigma(i) .gt. 0) slope(i) = -(log(sigma(i))-log(sigma(i-1)))/(log(mms_fine(i))-log(mms_fine(i-1)))
	end do
  
  ! adjust sigma for redshift
  
 	do i=1, dim_ps
		sigma_z(i) = sigma(i) / growth(omega0,lambda0,z) 													     
	end	do

	do i=1, dim_ps		
		slope_z(i) = slope(i)!-(log(sigma_z(i))-log(sigma_z(i-1)))/(log(mms_fine(i))-log(mms_fine(i-1)))

	enddo
  
	do i=1, dim
		sigma_z_coarse(i) = sigma_coarse(i) / growth(omega0,lambda0,z) 
	end	do

  do i=2,dim
		if(sigma_coarse(i) .gt. 0) slope_coarse(i) = -(log(sigma_coarse(i))-log(sigma_coarse(i-1)))/&
                      &(log(data_mass_low(i))-log(data_mass_low(i-1)))
	end do

	do i=1, dim	
    slope_z_coarse(i) = slope_coarse(i) !-(log(sigma_z_coarse(i))-log(sigma_z_coarse(i-1)))/&
                                        !&(log(mms_fine(i))-log(mms_fine(i-1)))
	enddo

	do j=2,dim_ps
    
     nu = deltac_z/sigma_z(j)

     m_ave = (mms_fine(j)+mms_fine(j-1))/2.

    if(sigma_z(j) .gt. 0) then
    
		n_fine(j)=sqrt(2./pi)*rho_bar/m_ave**2*nu*exp(-nu**2/2.)*slope(j)!m_sun^-1 Mpc^-3
		
    n_ST_fine(j)=Apar*(1.0+(sqrt(alphapar)*nu)**(-2.*ppar))*sqrt(2./pi)*rho_bar/m_ave**2* &
		nu*sqrt(alphapar)*exp(-alphapar*nu**2.0/2.)*slope(j)!m_sun^-1 Mpc^-3

		f_ST(j) = n_ST_fine(j) * m_ave**2 / slope(j) / rho_bar
		f_PS(j) = n_fine(j) * m_ave**2 / slope(j) / rho_bar
    
    else
    
    n_fine(j) = 0
    n_ST_fine(j) = 0
    f_ST(j) = 0
    f_PS(j) =0    
    
    endif
                    
  end do
  
    

	! do number of halos integrals over coarse mass bins


	do i=1,dim
        data_ps(i) = 0.0d0
        data_st(i) = 0.0d0
	end do


	
	do j=2,dim
		do i=2,dim_ps       
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_ps(j) = data_ps(j)+(n_fine(i)+n_fine(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
				data_st(j) = data_st(j)+(n_ST_fine(i)+n_ST_fine(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do
		
        data_ps(j)=data_ps(j)/(data_mass_high(j)-data_mass_low(j)) 
        data_st(j)=data_st(j)/(data_mass_high(j)-data_mass_low(j))
	  
	end do
	

	print*,'done PS and ST'
  

	
! Now do Jenkins Approximation

  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.
		
    if(sigma_z(i) .gt. 0) then
    
		f_jenkins(i) = 0.315*exp(-(abs(log(1/sigma_z(i))+0.61))**3.8)
		
    else
    
    f_jenkins(i) = 0
    
    endif
    
		n_jenkins(i) = f_jenkins(i) * rho_bar * slope_z(i) / m_ave**2
    

		
		
	enddo
	
! aggregate into coarser bins:
		
	do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_jenkins(j) = data_jenkins(j)+(n_jenkins(i)+n_jenkins(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do
	

		data_jenkins(j)=data_jenkins(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM



	end do

  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.
    
    if(sigma_z(i) .gt. 0) then
		
		f_warren(i) = 0.7234*(sigma_z(i)**(-1.625) + 0.2538)*exp(-1.1982/sigma_z(i)**2)
		
    else
    
    f_warren(i) = 0
    
    endif
    
		n_warren(i) = f_warren(i) * rho_bar * slope_z(i) / m_ave**2
	
	enddo
	

do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_warren(j) = data_warren(j)+(n_warren(i)+n_warren(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do

	data_warren(j)=data_warren(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
		
	end do

! Tinker

 alpha_tink = 10**(-(0.75/log10(200.0/75.0))**1.2)


 big_a_tink = 0.186*(1.0+z)**(-0.14)
 small_a_tink = 1.47*(1.0+z)**(-0.06)
 c_tink = 1.19
 b_tink = 2.57*(1.0+z)**(-alpha_tink)



  
  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.
    
    if(sigma_z(i) .gt. 0) then
		
		f_tinker(i) = big_a_tink*((sigma_z(i)/b_tink)**(-small_a_tink)+1)*exp(-c_tink/sigma_z(i)**2)
		
    else
    
    f_tinker(i) = 0
    
    endif
    
		n_tinker(i) = f_tinker(i) * rho_bar * slope_z(i) / m_ave**2
	
	enddo
	

do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_tinker(j) = data_tinker(j)+(n_tinker(i)+n_tinker(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do

	data_tinker(j)=data_tinker(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
		
	end do


  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.
		
    if(sigma_z(i) .gt. 0) then
    
		f_reed_2003(i) = f_ST(i)*exp(-0.7/(sigma_z(i)*cosh(2*sigma_z(i))**5))
		
    else
    
    f_reed_2003(i) = 0
    
    endif
    
		n_reed_2003(i) = f_reed_2003(i) * rho_bar * slope_z(i) / m_ave**2
	
	enddo


do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_reed_2003(j) = data_reed_2003(j)+(n_reed_2003(i)+n_reed_2003(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do

	data_reed_2003(j)=data_reed_2003(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
		
	end do

  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.
                
    G_1(i) = exp(-(((log(1/sigma_z(i))-0.4)**2.0)/(2.0*0.6**2.0)))
                
    G_2(i) = exp(-(((log(1/sigma_z(i))-0.75)**2.0)/(2.0*0.2**2.0)))
                
    n_eff(i) = 6.*slope(i)  - 3.

    if(sigma_z(i) .gt. 0) then
                
    f_reed_2007(i) = Apar*sqrt(2.*alphapar/pi)*(1+((sigma_z(i)**2)/(alphapar*(deltac_z)**2))**ppar+0.6*G_1(i)+0.4*G_2(i))&
    *(deltac_z/sigma_z(i))*exp(-(c_par*alphapar*(deltac_z**2/sigma_z(i)**2)/2.) - 0.03*((deltac_z/sigma_z(i))**0.6)/(n_eff(i)+3)**2)
                
    else

    f_reed_2007(i) = 0

    endif

                n_reed_2007(i) = f_reed_2007(i) * rho_bar * slope_z(i) / m_ave**2
        
        enddo


do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_reed_2007(j) = data_reed_2007(j)+(n_reed_2007(i)+n_reed_2007(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do

	data_reed_2007(j)=data_reed_2007(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
		
	end do
  

  ! Angulo:
  
  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.0
    
    if(sigma_z(i) .gt. 0) then
		
		f_angulo(i) = 0.201*((2.08/sigma_z(i))**(1.7)+1.0)*exp(-1.172/sigma_z(i)**2)
		
    else
    
    f_angulo(i) = 0
    
    endif
    
		n_angulo(i) = f_angulo(i) * rho_bar * slope_z(i) / m_ave**2

	enddo
	

  do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_angulo(j) = data_angulo(j)+(n_angulo(i)+n_angulo(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do

	data_angulo(j)=data_angulo(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
		
	end do

! Crocce:
  
    large_a_crocce = 0.58*(1.0+z)**(-0.13)
    a_crocce = 1.37*(1.0+z)**(-0.15)
    b_crocce = 0.3*(1.0+z)**(-0.084)
    c_crocce = 1.036*(1.0+z)**(-0.024)

  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.
		
    if(sigma_z(i) .gt. 0) then
    
		f_crocce(i) = large_a_crocce*(sigma_z(i)**(-a_crocce) + b_crocce)*exp(-c_crocce/(sigma_z(i)**2.0))
		
    else
    
    f_crocce(i) = 0
    
    endif
    
		n_crocce(i) = f_crocce(i) * rho_bar * slope_z(i) / m_ave**2

	enddo
	

  do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_crocce(j) = data_crocce(j)+(n_crocce(i)+n_crocce(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do

	data_crocce(j)=data_crocce(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
		
	end do

! Courtin:
  
  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.
		
    if(sigma_z(i) .gt. 0) then 
           
		f_courtin(i) = large_a_courtin*(2.0*small_a_courtin/pi)**0.5 * &
    & (deltac_z/sigma_z(i))*(1 + (deltac_z/(sigma_z(i)*small_a_courtin**0.5))**(-2.0*p_courtin))*&
    & exp(-deltac_z**2.0*small_a_courtin/(2.0*sigma_z(i)**2))
		
    else
    
    f_courtin(i) = 0
    
    endif
    
    
		n_courtin(i) = f_courtin(i) * rho_bar * slope_z(i) / m_ave**2

	enddo
	

  do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_courtin(j) = data_courtin(j)+(n_courtin(i)+n_courtin(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do

	data_courtin(j)=data_courtin(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
		
	end do

! Watson:
  
     a_watson_cpm_ahf = omega_z*(1.097*(1.0+z)**(-3.216)+0.074)
     b_watson_cpm_ahf = omega_z*(3.136*(1.0+z)**(-3.058)+2.349)
     c_watson_cpm_ahf = omega_z*(5.907*(1.0+z)**(-3.599)+2.344)
     d_watson_cpm_ahf = 1.318

print*, omega_z
  
  
  do i=1, dim_ps
	
      if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
      if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.

		
    if(sigma_z(i) .gt. 0) then 
           
		f_watson_fof_uni(i) = a_watson_fof_uni*((b_watson_fof_uni/sigma_z(i))**c_watson_fof_uni + 1.0)&
    &*exp(-d_watson_fof_uni/sigma_z(i)**2)
		f_watson_cpm_ahf(i) = a_watson_cpm_ahf*((b_watson_cpm_ahf/sigma_z(i))**c_watson_cpm_ahf + 1.0)&
    &*exp(-d_watson_cpm_ahf/sigma_z(i)**2)		    

    else
    
    f_watson_fof_uni(i) = 0
    f_watson_cpm_ahf(i) = 0
    
    endif
    
    
		n_watson_fof_uni(i) = f_watson_fof_uni(i) * rho_bar * slope_z(i) / m_ave**2
		n_watson_cpm_ahf(i) = f_watson_cpm_ahf(i) * rho_bar * slope_z(i) / m_ave**2

	enddo
	

  do j=2,dim
		do i=2,dim_ps
			if(mms_fine(i)>data_mass_low(j).and.mms_fine(i)<=data_mass_high(j))then
				data_watson_fof_uni(j) = data_watson_fof_uni(j)+(n_watson_fof_uni(i)+n_watson_fof_uni(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
				data_watson_cpm_ahf(j) = data_watson_cpm_ahf(j)+(n_watson_cpm_ahf(i)+n_watson_cpm_ahf(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
			end if
		end do

  	data_watson_fof_uni(j)=data_watson_fof_uni(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
  	data_watson_cpm_ahf(j)=data_watson_cpm_ahf(j)/(data_mass_high(j)-data_mass_low(j)) ! dividing dn by dM
		
	end do


  ! convert to f(sigma_coarse)
  
  do i = 1, dim
  
    if(slope_coarse(i) .ne. 0) then
  
		data_ps(i) = data_ps(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_st(i) = data_st(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_jenkins(i) = data_jenkins(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_warren(i) = data_warren(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_reed_2003(i) = data_reed_2003(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_reed_2007(i) = data_reed_2007(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_tinker(i) = data_tinker(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_angulo(i) = data_angulo(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_crocce(i) = data_crocce(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_courtin(i) = data_courtin(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_watson_fof_uni(i) = data_watson_fof_uni(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
		data_watson_cpm_ahf(i) = data_watson_cpm_ahf(i) * data_mass_ave(i)**2 / rho_bar / slope_coarse(i)
  
    else
    
		data_ps(i) = 0
		data_st(i) = 0
		data_jenkins(i) = 0
		data_warren(i) = 0
		data_reed_2003(i) = 0
		data_reed_2007(i) = 0
		data_tinker(i) = 0
		data_angulo(i) = 0
		data_crocce(i) = 0
		data_courtin(i) = 0
		data_watson_fof_uni(i) = 0
		data_watson_cpm_ahf(i) = 0
    
    endif
  
  enddo
  
  


    open (unit=19,file=afile,form='formatted')
    
    write(19,'(a185)') '#Mass / Msun  #ln(sigma^-1) #Press-Schechter #Sheth-Torman &
    & #Jenkins  #Warren   #Reed (2003) #Reed (2007)  #Tinker   #Angulo    #Crocce   #Courtin   &
    &#Watson_ahf_z0    #Watson_cpm_ahf'
	   
     	   do i=2, dim_ps
         
          if(i .eq. 1) m_ave = (mms_fine(i)+m_min)/2.0
          if(i .gt. 1) m_ave = (mms_fine(i)+mms_fine(i-1))/2.
 
     
          write(19, '(2es18.10,12f15.10)') mms_fine_ave(i),log(1.0/(sigma_z(i))),&
          &log(f_PS(i)), log(f_ST(i)), log(f_jenkins(i)), log(f_warren(i)), &
          &log(f_reed_2003(i)), log(f_reed_2007(i)), log(f_tinker(i)), log(f_angulo(i)),&
	  &log(f_crocce(i)), log(f_courtin(i)), log(f_watson_fof_uni(i)), log(f_watson_cpm_ahf(i))
 

         enddo
     
    close(19)
    
    ! calculate spline so we can interpolate sigma values across the mass range:

    call spline(mms_fine(1:dim_ps),sigma,dim_ps,1d31,1d31,sigma_bin_ave_del)
    call spline(mms_fine(1:dim_ps),sigma_z,dim_ps,1d31,1d31,sigma_bin_ave_del_z)



    call splint(mass,sigma1,sigma_cobe,max_size,m_box,sigma_box)
    
    sigma_box = sigma_box!/growth(omega0,lambda0,z)
  
    print*, 'sigma (squared) of box mass = ',sigma_box,'(',sigma_box**2.0,')'
    
!  call splint(mass,sigma1,sigma_cobe,max_size,m_grid*400.0*8.0,sigma_400_part)
  
!  sigma_400_part = sigma_400_part / growth(omega0,lambda0,z)  
  
!  print*, '400 particle cut-off (old) is at ln(sigma^-1) = ',log(sigma_400_part**(-1.0))
  
!  sigma_400_part = ((sigma_400_part*growth(omega0,lambda0,z))**2.0 + sigma_box**2.0)**0.5
  
!  sigma_400_part = sigma_400_part / growth(omega0,lambda0,z)  
  
!  print*, '400 particle cut-off (new) is at ln(sigma^-1) = ',log(sigma_400_part**(-1.0))  

  ! We have now calculated analytic predictions for a global mass function.
  ! We need to now adjust our sigma values to reflect that we are dealing with
  ! a finite volume, not an infinite one

  ! correct sigmas for finite volumes:

  if(volume_correction .eq. 1) then

      do i = 1, dim_ps

      if((sigma(i)**2.0 - sigma_box**2.0) .lt. 0 .and. &
      &(sigma(i-1)**2.0 - sigma_box**2.0) .gt. 0) then 
      
        sigma_limit = i-1
        goto 100
  
      endif
      enddo
  
100  print*, 'sigma limit = ', sigma_limit
     print*, 'sigma (global) = ',sigma(sigma_limit)
     print*
  
     allocate(sigma_local(sigma_limit))
     allocate(sigma_local_sorted_by_mass(sigma_limit))
     allocate(sigma_local_del(sigma_limit))
     allocate(mass_local_del(sigma_limit))
     allocate(mms_fine_sorted_by_sigma(sigma_limit))
        
  print*      
  do i=1,sigma_limit
       sigma_local(i) = (sigma(i)**2.0 - sigma_box**2.0)**0.5 
	end do

    
   ! sort the sigma_local and mms_fine arrays
    
    mms_fine_sorted_by_sigma = mms_fine(1:sigma_limit)
    sigma_local_sorted_by_mass = sigma_local
    
    call ssort(sigma_local, mms_fine_sorted_by_sigma, sigma_limit, 2)

   ! calculate spline to allow us to interpolate mass values for the local sigma
  
  	call spline(sigma_local,mms_fine_sorted_by_sigma,sigma_limit,1d31,1d31,mass_local_del)
  	call spline(mms_fine(1:sigma_limit),sigma_local_sorted_by_mass,sigma_limit,1d31,1d31,sigma_local_del)
    

    do i = 1, dim
    
    call splint(mms_fine(1:sigma_limit),sigma_local_sorted_by_mass,sigma_local_del,sigma_limit,&
    &data_mass_low(i),sigma_mass_low(i))
    call splint(mms_fine(1:sigma_limit),sigma_local_sorted_by_mass,sigma_local_del,sigma_limit,&
    &data_mass_high(i),sigma_mass_high(i))
    call splint(mms_fine(1:sigma_limit),sigma_local_sorted_by_mass,sigma_local_del,sigma_limit,&
    &data_mass_ave(i),sigma_mass_ave(i))
        
    enddo
    
    
  endif
  
  
  

!---------------------------------------------------------------------------------------------


! Now read the cubep3m SOD halofinder data
 
 
    if(analytic_only .neqv. .true.) then


      do i=0, nodes_dim**3-1
  
        write(file,'(i4)') i
          file=adjustl(file)
				
        ifile=input_path//redshift(1:len_trim(redshift))//'halo'//file(1:len_trim(file))//'.dat'	
 
	        
#ifdef BINARY
        open (unit=13,file=ifile,form='binary')    
#else
        open (unit=13,file=ifile,access='stream')
#endif



        print*, 'Halo file "', ifile(1:len_trim(ifile)), '" successfully opened'

        read(13) n_halo
        print*,'Number of Halos: ', n_halo
		
        do k=1, n_halo
	
          if(cpm_version .eq. 0)  then
        
            read(13, end=122) halo_input_buffer(1:17)
      
            do j=1,dim
              if(halo_input_buffer(15)*m_grid .gt. data_mass_low(j) .and. &
                & halo_input_buffer(15)*m_grid .le. data_mass_high(j) .and. &
                  halo_input_buffer(15)/8.0 .gt. min_part_cutoff) then
                    data_cubep3m(j) = data_cubep3m(j) + 1.0
                    total_cubep3m(j) = total_cubep3m(j) + halo_input_buffer(15)*m_grid                
              endif
            enddo	
        

          elseif(cpm_version .eq. 1) then

            read(13, end=122) halo_input_buffer, pid_halo, xv_halo
      
             do j=1,dim
               if(halo_input_buffer(17)*m_grid .gt. data_mass_low(j) .and. &
                & halo_input_buffer(17)*m_grid .le. data_mass_high(j) .and. &
                  halo_input_buffer(17)/8.0 .gt. min_part_cutoff) then
                    data_cubep3m(j) = data_cubep3m(j) + 1.0
                    total_cubep3m(j) = total_cubep3m(j) + halo_input_buffer(17)*m_grid                
                endif
             enddo	

          elseif(cpm_version .eq. 2) then
        
            read(13, end=122) halo_input_buffer,I_ij
      
            do j=1,dim
              if(halo_input_buffer(17)*m_grid .gt. data_mass_low(j) .and. &
                & halo_input_buffer(17)*m_grid .le. data_mass_high(j) .and. &
                  halo_input_buffer(17)/8.0 .gt. min_part_cutoff) then
                    data_cubep3m(j) = data_cubep3m(j) + 1.0
                    total_cubep3m(j) = total_cubep3m(j) + halo_input_buffer(17)*m_grid                
              endif
            enddo	

          endif
        
        enddo

122     close(13)

      enddo
      
      
      do i=1,dim

        data_mass_cubep3m(i) = total_cubep3m(i) / data_cubep3m(i)

      enddo


    ! calculate errors, NB must do this before volume correction!

  do i = 1, dim
    
    error_cubep3m_low(i) = (real(data_cubep3m(i))+0.25)**0.5 - 0.5
    error_cubep3m_high(i) = (real(data_cubep3m(i))+0.25)**0.5 + 0.5
    
  enddo


  ! interpolate the sigma values to get the value of sigma at the bin average:
  
  if(volume_correction .eq. 0) then
  
  do i=1,dim+1
		call splint(mms_fine(1:dim_ps),sigma,sigma_bin_ave_del,dim_ps,data_mass_cubep3m(i),sigma_bin_ave(i))
		call splint(mms_fine(1:dim_ps),sigma_z,sigma_bin_ave_del_z,dim_ps,data_mass_cubep3m(i),sigma_bin_ave_z(i))
	end do
  
  endif
  
  if(volume_correction .eq. 1) then
  
  do i = 1, dim
  call splint(mms_fine(1:sigma_limit),sigma_local_sorted_by_mass,sigma_local_del,sigma_limit,&
  &data_mass_cubep3m(i),sigma_bin_ave(i))
  
  enddo
  
  
    ! convert halo counts and masses to account for finite volume
    
  do i = 1, dim
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_bin_ave(i),mass1)
    data_cubep3m(i) = data_cubep3m(i)*data_mass_cubep3m(i)/mass1
    error_cubep3m_low(i) = error_cubep3m_low(i)*data_mass_cubep3m(i)/mass1
    error_cubep3m_high(i) = error_cubep3m_high(i)*data_mass_cubep3m(i)/mass1

    
    data_mass_cubep3m(i) = mass1


  enddo
  
  ! convert bin masses to account for finite volume
  
  do i = 1, dim
  
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_low(i),mass1)
    data_mass_low(i) = mass1
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_high(i),mass1)
    data_mass_high(i) = mass1
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_ave(i),mass1)
    data_mass_ave(i) = mass1
  
  enddo
  
  endif
  
  	! Convert Counts of halos into log(M^2/rho_bar * dN/dM) and normalise to unit volume:

	do i=1,dim
  
		data_cubep3m(i) = (data_cubep3m(i) / (box/h)**3 * data_mass_cubep3m(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))

		error_cubep3m_low(i) = (error_cubep3m_low(i) / (box/h)**3 * data_mass_cubep3m(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))

		error_cubep3m_high(i) = (error_cubep3m_high(i) / (box/h)**3 * data_mass_cubep3m(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))

	enddo
  
  
  ! reset bins for data:
		
	do i=1,dim
															
		data_mass_high(i) = m_min*10.**(mass_bin_fac_coarse*(real(i))/bin_dim) !mass bins for plotting
			
		if(i==1) then
			data_mass_low(i) = m_min-1
		else
			data_mass_low(i) = data_mass_high(i-1)		
		endif
		
		data_mass_ave(i) = data_mass_low(i) + (data_mass_high(i) - data_mass_low(i))/2
    
		
	enddo

  
  
  
  
  ! convert to f(sigma_coarse)
  
  do i = 1, dim
    
		data_cubep3m(i) = data_cubep3m(i) / slope_coarse(i)
	
		error_cubep3m_low(i) = error_cubep3m_low(i) / slope_coarse(i)		
		error_cubep3m_high(i) = error_cubep3m_high(i) / slope_coarse(i)		
    
    sigma_bin_ave_z(i) = sigma_bin_ave(i) / growth(omega0,lambda0,z)

  enddo


	! Print data into a file:
	
    open (unit=45,file=cpm_file,form='formatted')

    write(45,'(a218)') '#Mass/Msun (bin ave) #ln(sigma^-1) (bin ave)  #Mass/Msun (bin centre) #ln(sigma^-1) (Bin centre)&
      &  #particles (bin ave) #particles (bin centre) #particles (bin low cut) #ln(f(sigma)    #error low        #error high'
	   
     	   do i=1, dim   
     
          write(45,'(4es20.6,3i20,3f25.6)') data_mass_cubep3m(i), log(1.0/sigma_bin_ave_z(i)),&
           &data_mass_ave(i), log(1.0/sigma_z_coarse(i)),&
           &int(data_mass_cubep3m(i)/m_grid/8.0),int(data_mass_ave(i)/m_grid/8.0),int(data_mass_low(i)/m_grid/8.0),&
           &log(data_cubep3m(i)),log(data_cubep3m(i)-error_cubep3m_low(i)),&
           &log(data_cubep3m(i)+error_cubep3m_high(i))
     
         enddo
     
	  close(45)
    


    ! Now calculate residuals. For this we need to interpolate the values of the analytic 
    ! functions at the precise values of sigma for which we have numerical data
    
    call spline(mms_fine_ave,f_PS,dim_ps,1d31,1d31,data_ps_del)
    call spline(mms_fine_ave,f_ST,dim_ps,1d31,1d31,data_st_del)
    call spline(mms_fine_ave,f_jenkins,dim_ps,1d31,1d31,data_jen_del)
    call spline(mms_fine_ave,f_warren,dim_ps,1d31,1d31,data_war_del)
    call spline(mms_fine_ave,f_tinker,dim_ps,1d31,1d31,data_tin_del)
    call spline(mms_fine_ave,f_reed_2003,dim_ps,1d31,1d31,data_r03_del)
    call spline(mms_fine_ave,f_reed_2007,dim_ps,1d31,1d31,data_r07_del)
    call spline(mms_fine_ave,f_angulo,dim_ps,1d31,1d31,data_ang_del)                    
    call spline(mms_fine_ave,f_crocce,dim_ps,1d31,1d31,data_cro_del)
    call spline(mms_fine_ave,f_courtin,dim_ps,1d31,1d31,data_cou_del)    
    call spline(mms_fine_ave,f_watson_fof_uni,dim_ps,1d31,1d31,data_watson_fof_uni_del)    
    call spline(mms_fine_ave,f_watson_cpm_ahf,dim_ps,1d31,1d31,data_watson_cpm_ahf_del)    

  	do i=1,dim
    
  		call splint(mms_fine_ave,f_PS,data_ps_del,dim_ps,data_mass_cubep3m(i),data_ps_int(i))
  		call splint(mms_fine_ave,f_ST,data_st_del,dim_ps,data_mass_cubep3m(i),data_st_int(i))	  
  		call splint(mms_fine_ave,f_jenkins,data_jen_del,dim_ps,data_mass_cubep3m(i),data_jen_int(i))
  		call splint(mms_fine_ave,f_warren,data_war_del,dim_ps,data_mass_cubep3m(i),data_war_int(i))
      call splint(mms_fine_ave,f_tinker,data_tin_del,dim_ps,data_mass_cubep3m(i),data_tin_int(i))
      call splint(mms_fine_ave,f_reed_2003,data_r03_del,dim_ps,data_mass_cubep3m(i),data_r03_int(i))
      call splint(mms_fine_ave,f_reed_2007,data_r07_del,dim_ps,data_mass_cubep3m(i),data_r07_int(i))
      call splint(mms_fine_ave,f_angulo,data_ang_del,dim_ps,data_mass_cubep3m(i),data_ang_int(i))
      call splint(mms_fine_ave,f_crocce,data_cro_del,dim_ps,data_mass_cubep3m(i),data_cro_int(i))
      call splint(mms_fine_ave,f_courtin,data_cou_del,dim_ps,data_mass_cubep3m(i),data_cou_int(i))
      call splint(mms_fine_ave,f_watson_fof_uni,data_watson_fof_uni_del,dim_ps,data_mass_cubep3m(i),data_watson_fof_uni_int(i))
      call splint(mms_fine_ave,f_watson_cpm_ahf,data_watson_cpm_ahf_del,dim_ps,data_mass_cubep3m(i),data_watson_cpm_ahf_int(i))
                                                                                                                  
    end do
    
    open (unit=46,file=cpm_dfile,form='formatted')        
    
    
    write(46,'(a1000)') '#Mass/Msun (bin ave)   #ln(sigma^-1) #Mass/Msun (bin centre) #ln(sigma^-1) (Bin centre)&
    & #particles (bin ave) #particles (bin centre) #particles (bin low cut) &
    &   #(data-PS)/PS             #error PS low&
    &            #error PS high           #(data-ST)/ST             #error ST low          #error ST high &
    &         #(data)/Jen            #error Jen low           # error Jen high &
    &         #(data)/War         #error War low           # error War high &
    &         #(data)/R03         #error R03 low           # error R03 high &
    &         #(data)/War         #error R07 low           # error R07 high &    
    &         #(data)/Tin          #error Tin low           # error Tin high    &
    &         #(data)/Ang         #error Ang low           # error Ang high      &
    &         #(data)/Cro         #error Cro low           # error Cro high      &
    &         #(data)/Cou         #error Cou low           # error Cou high      &
    &         #(data)/Wat_fof_uni         #error Wat_fof_uni low           # error Wat_fof_uni high          '    
    
   do i = 1, dim
    
      write(46, '(4es18.10,3i20,36f15.10)') data_mass_cubep3m(i), log(1.0/sigma_bin_ave_z(i)),&
      &data_mass_ave(i), log(1.0/sigma_z_coarse(i)),&
      &int(data_mass_cubep3m(i)/m_grid/8.0),int(data_mass_ave(i)/m_grid/8.0),int(data_mass_low(i)/m_grid/8.0),&
      &(data_cubep3m(i))/data_ps_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_ps_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_ps_int(i),&
      &(data_cubep3m(i))/data_st_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_st_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_st_int(i),&
      &(data_cubep3m(i))/data_jen_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_jen_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_jen_int(i),&
      &(data_cubep3m(i))/data_war_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_war_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_war_int(i),&
      &(data_cubep3m(i))/data_r03_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_r03_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_r03_int(i),&
      &(data_cubep3m(i))/data_r07_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_r07_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_r07_int(i),&
      &(data_cubep3m(i))/data_tin_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_tin_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_tin_int(i),&
      &(data_cubep3m(i))/data_ang_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_ang_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_ang_int(i),&
      &(data_cubep3m(i))/data_cro_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_cro_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_cro_int(i),&
      &(data_cubep3m(i))/data_cou_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_cou_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_cou_int(i),&
      &(data_cubep3m(i))/data_watson_fof_uni_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_watson_fof_uni_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_watson_fof_uni_int(i),&
      &(data_cubep3m(i))/data_watson_cpm_ahf_int(i),&
      &(data_cubep3m(i)-error_cubep3m_low(i))/data_watson_cpm_ahf_int(i),&
      &(data_cubep3m(i)+error_cubep3m_high(i))/data_watson_cpm_ahf_int(i)

    
    enddo

  close(46)
  
!---------------------------------------------------------------------------------------------  
  
! Now read the AHF data
    
      do i=0, num_files_ahf-1

        write(file,'(i4)') i
        file=adjustl(file)
				
        if(ahf_overdensity .eq. 0) then
        ifile=input_path//redshift(1:len_trim(redshift))//'_AHF_halos_chunk_'&
              &//file(1:len_trim(file))//'.dat'
        endif
        
        if(ahf_overdensity .ne. 0) then
        
        write(overdensity_string,'(i5)') ahf_overdensity
        overdensity_string=adjustl(overdensity_string)
        
        ifile=input_path//redshift(1:len_trim(redshift))//'_AHF_halos_chunk_'&
              &//file(1:len_trim(file))//'_'//overdensity_string(1:len_trim(overdensity_string))//'.dat'
        
        endif

        open (unit=23,file=ifile,form='formatted')
      
        print*, ifile(1:len_trim(ifile)),' opened successfully'
		
        ahf_halos_in_file = 0
		
        k = 0
			
        do ! This is unbounded as it should exit when the file has been read	
			
      
          if(ahf_version .eq. 0) read(23,*, end=123) npart

          if(ahf_version .eq. 1) read(23,*, end=123) ID,HostHalo,NumSubStruct,M_vir,npart!,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir,  &  
!      &  R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx,& 
!      &  Ecy,Ecz,ovdens,nbins,fMhires,Ekin,Epot,SurfP,Phi0,cNFW 

          halo_mass = npart*m_grid*8 ! multiply by 8 to convert from grid mass to particle mass
          ahf_halos_in_file = ahf_halos_in_file + 1
			
          do j=1,dim
            if(halo_mass .gt. data_mass_low(j)) then
              if(halo_mass .le. data_mass_high(j) .and. &
                  npart .gt. min_part_cutoff) then
                  
                    if(ahf_version .eq. 0) then
                    
                    data_ahf(j) = data_ahf(j) + 1.0
                    total_ahf_chunk(j) = total_ahf_chunk(j) + halo_mass                                  
                    
                    endif
                    
                    if(ahf_version .eq. 1 .and. ahf_substructure .eq. 1) then
                    
                      data_ahf(j) = data_ahf(j) + 1.0
                      total_ahf_chunk(j) = total_ahf_chunk(j) + halo_mass              
                    
                    endif
                    
                    if(ahf_version .eq. 1 .and. ahf_substructure .eq. 0 &
                        & .and. HostHalo .eq. -1) then                    
                    
                      data_ahf(j) = data_ahf(j) + 1.0
                      total_ahf_chunk(j) = total_ahf_chunk(j) + halo_mass              

                    endif

              endif
            endif
          enddo

          k = k + 1
      
        enddo
		
        123 close(23)
		
        print*, 'halos in file = ', ahf_halos_in_file
        
      enddo

      do i=1,dim

        data_mass_ahf(i) = total_ahf_chunk(i) / data_ahf(i)

      enddo

    ! calculate errors, NB must do this before volume correction!

  do i = 1, dim
  
    error_ahf_low(i) = (real(data_ahf(i))+0.25)**0.5 - 0.5
    error_ahf_high(i) = (real(data_ahf(i))+0.25)**0.5 + 0.5 
    
  enddo


  ! interpolate the sigma values to get the value of sigma at the bin average:
  
  if(volume_correction .eq. 0) then
  
  do i=1,dim+1
		call splint(mms_fine(1:dim_ps),sigma,sigma_bin_ave_del,dim_ps,data_mass_ahf(i),sigma_bin_ave(i))
		call splint(mms_fine(1:dim_ps),sigma_z,sigma_bin_ave_del_z,dim_ps,data_mass_ahf(i),sigma_bin_ave_z(i))
	end do
  
  endif
  
  if(volume_correction .eq. 1) then
  
  do i = 1, dim
  call splint(mms_fine(1:sigma_limit),sigma_local_sorted_by_mass,sigma_local_del,sigma_limit,&
  &data_mass_ahf(i),sigma_bin_ave(i))
  enddo
  
    ! convert halo counts to account for finite volume
    
  do i = 1, dim
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_bin_ave(i),mass1)
    data_ahf(i) = data_ahf(i)*data_mass_ahf(i)/mass1
    error_ahf_low(i) = error_ahf_low(i)*data_mass_ahf(i)/mass1
    error_ahf_high(i) = error_ahf_high(i)*data_mass_ahf(i)/mass1    
    

    
    data_mass_ahf(i) = mass1

  enddo

  ! convert bin masses to account for finite volume
  
  do i = 1, dim
  
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_low(i),mass1)
    data_mass_low(i) = mass1
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_high(i),mass1)
    data_mass_high(i) = mass1
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_ave(i),mass1)
    data_mass_ave(i) = mass1
  
  enddo

  endif

  
  	! Convert Counts of halos into log(M^2/rho_bar * dN/dM) and normalise to unit volume:

	do i=1,dim
  
		data_ahf(i) = (data_ahf(i) / (box/h)**3 * data_mass_ahf(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))
		
		error_ahf_low(i) = (error_ahf_low(i) / (box/h)**3 * data_mass_ahf(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))

		error_ahf_high(i) = (error_ahf_high(i) / (box/h)**3 * data_mass_ahf(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))                            
                
	enddo



  ! reset bins for data:
		
	do i=1,dim
															
		data_mass_high(i) = m_min*10.**(mass_bin_fac_coarse*(real(i))/bin_dim) !mass bins for plotting
			
		if(i==1) then
			data_mass_low(i) = m_min-1
		else
			data_mass_low(i) = data_mass_high(i-1)		
		endif
		
		data_mass_ave(i) = data_mass_low(i) + (data_mass_high(i) - data_mass_low(i))/2
    
		
	enddo

  
  ! convert to f(sigma_coarse)
  
  do i = 1, dim
  
		data_ahf(i) = data_ahf(i) / slope_coarse(i)
		error_ahf_low(i) = error_ahf_low(i) / slope_coarse(i)
		error_ahf_high(i) = error_ahf_high(i) / slope_coarse(i)
    
    sigma_bin_ave_z(i) = sigma_bin_ave(i) / growth(omega0,lambda0,z) 													     

  enddo



	! Print data into a file:
	
    open (unit=45,file=ahf_file,form='formatted')

    write(45,'(a218)') '#Mass/Msun (bin ave) #ln(sigma^-1) (bin ave)  #Mass/Msun (bin centre) #ln(sigma^-1) (Bin centre)&
      &  #particles (bin ave) #particles (bin centre) #particles (bin low cut) #ln(f(sigma)    #error low        #error high'
	   
     	   do i=1, dim   
     
          write(45,'(4es20.6,3i20,3f25.6)') data_mass_ahf(i), log(1.0/sigma_bin_ave_z(i)),&
           &data_mass_ave(i), log(1.0/sigma_z_coarse(i)),&
           &int(int(data_mass_ahf(i)/m_grid/8.0)),int(data_mass_ave(i)/m_grid/8.0),int(data_mass_low(i)/m_grid/8.0),&
           &log(data_ahf(i)),log(data_ahf(i)-error_ahf_low(i)),&
           &log(data_ahf(i)+error_ahf_high(i))
     
         enddo
     
	  close(45)
    


    ! Now calculate residuals. For this we need to interpolate the values of the analytic 
    ! functions at the precise values of sigma for which we have numerical data
    
  	do i=1,dim
                                                        
  		call splint(mms_fine_ave,f_PS,data_ps_del,dim_ps,data_mass_ahf(i),data_ps_int(i))
  		call splint(mms_fine_ave,f_ST,data_st_del,dim_ps,data_mass_ahf(i),data_st_int(i))	  
  		call splint(mms_fine_ave,f_jenkins,data_jen_del,dim_ps,data_mass_ahf(i),data_jen_int(i))
  		call splint(mms_fine_ave,f_warren,data_war_del,dim_ps,data_mass_ahf(i),data_war_int(i))
      call splint(mms_fine_ave,f_tinker,data_tin_del,dim_ps,data_mass_ahf(i),data_tin_int(i))
      call splint(mms_fine_ave,f_reed_2003,data_r03_del,dim_ps,data_mass_ahf(i),data_r03_int(i))
      call splint(mms_fine_ave,f_reed_2007,data_r07_del,dim_ps,data_mass_ahf(i),data_r07_int(i))
      call splint(mms_fine_ave,f_angulo,data_ang_del,dim_ps,data_mass_ahf(i),data_ang_int(i))
      call splint(mms_fine_ave,f_crocce,data_cro_del,dim_ps,data_mass_ahf(i),data_cro_int(i))
      call splint(mms_fine_ave,f_courtin,data_cou_del,dim_ps,data_mass_ahf(i),data_cou_int(i))
      call splint(mms_fine_ave,f_watson_fof_uni,data_watson_fof_uni_del,dim_ps,data_mass_ahf(i),data_watson_fof_uni_int(i))
      call splint(mms_fine_ave,f_watson_cpm_ahf,data_watson_cpm_ahf_del,dim_ps,data_mass_ahf(i),data_watson_cpm_ahf_int(i))
                                                    
    end do
    
    open (unit=46,file=ahf_dfile,form='formatted')        
    
    
    write(46,'(a1000)') '#Mass/Msun (bin ave)   #ln(sigma^-1) #Mass/Msun (bin centre) #ln(sigma^-1) (Bin centre)&
    & #particles (bin ave) #particles (bin centre) #particles (bin low cut) &
    &   #(data-PS)/PS             #error PS low&
    &            #error PS high           #(data-ST)/ST             #error ST low          #error ST high &
    &         #(data)/Jen            #error Jen low           # error Jen high &
    &         #(data)/War         #error War low           # error War high &
    &         #(data)/R03         #error R03 low           # error R03 high &
    &         #(data)/War         #error R07 low           # error R07 high &    
    &         #(data)/Tin          #error Tin low           # error Tin high    &
    &         #(data)/Ang         #error Ang low           # error Ang high      &
    &         #(data)/Cro         #error Cro low           # error Cro high      &
    &         #(data)/Cou         #error Cou low           # error Cou high      &
    &         #(data)/Wat_fof_uni         #error Wat_fof_uni low           # error Wat_fof_uni high          '    
    
   do i = 1, dim
    
      write(46, '(4es18.10,3i20,36f15.10)') data_mass_ahf(i), log(1.0/sigma_bin_ave_z(i)),&
      &data_mass_ave(i), log(1.0/sigma_z_coarse(i)),&
      &int(data_mass_ahf(i)/m_grid/8.0),int(data_mass_ave(i)/m_grid/8.0),int(data_mass_low(i)/m_grid/8.0),&
      &(data_ahf(i))/data_ps_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_ps_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_ps_int(i),&
      &(data_ahf(i))/data_st_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_st_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_st_int(i),&
      &(data_ahf(i))/data_jen_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_jen_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_jen_int(i),&
      &(data_ahf(i))/data_war_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_war_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_war_int(i),&
      &(data_ahf(i))/data_r03_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_r03_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_r03_int(i),&
      &(data_ahf(i))/data_r07_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_r07_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_r07_int(i),&
      &(data_ahf(i))/data_tin_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_tin_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_tin_int(i),&
      &(data_ahf(i))/data_ang_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_ang_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_ang_int(i),&
      &(data_ahf(i))/data_cro_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_cro_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_cro_int(i),&
      &(data_ahf(i))/data_cou_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_cou_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_cou_int(i),&
      &(data_ahf(i))/data_watson_fof_uni_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_watson_fof_uni_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_watson_fof_uni_int(i),&
      &(data_ahf(i))/data_watson_cpm_ahf_int(i),&
      &(data_ahf(i)-error_ahf_low(i))/data_watson_cpm_ahf_int(i),&
      &(data_ahf(i)+error_ahf_high(i))/data_watson_cpm_ahf_int(i)


    
    enddo

  close(46)
  


!----------------------------------------------------------------------------------------------------------------------------


! Read the FOF_chunk data

	 
	do file_end_int = 1, num_files_fof_chunk
	
	write(file,'(i4)') file_end_int-1
  file=adjustl(file)


      
  ifile=input_path//redshift(1:len_trim(redshift))//'_FOF_halos_chunk_'//file(1:len_trim(file))//'.dat'
          

#ifdef GFORTRAN		           
		open (unit=10,file=ifile,access='stream')
#endif
#ifdef BINARY
		open (unit=10,file=ifile,form='binary')
#endif        


	!read fof file header
        
	read(10) Ngroups_chunk
         	   
	print*,'Ngroups = ',Ngroups_chunk

	if(Ngroups_chunk /= 0) then
	
	! read file data
  
  do i = 1, Ngroups_chunk
  
	read(10) Length_chunk, Offset_chunk, Group_Mass_chunk, CM_chunk(:), vel_chunk(:)
  
  Length_chunk_corrected = real(Length_chunk)*(1.0-1.0/(real(Length_chunk)**0.6))
  

	! Put data from the first file into histogram table
			do j=1, dim
				if(real(Length_chunk*m_grid*8.0) .gt. data_mass_low(j) .and. &
					&real(Length_chunk*m_grid*8.0)	.le. data_mass_high(j)) then
					data_fof_chunk(j) = data_fof_chunk(j) + 1.0
          total_fof_chunk(j) = total_fof_chunk(j) +Length_chunk*m_grid*8.0
				endif
			enddo

			do j=1, dim
				if(real(Length_chunk_corrected*m_grid*8.0) .gt. data_mass_low(j) .and. &
					&real(Length_chunk_corrected*m_grid*8.0)	.le. data_mass_high(j)) then
					data_fof_chunk_cor(j) = data_fof_chunk_cor(j) + 1.0
          total_fof_chunk_cor(j) = total_fof_chunk_cor(j)+Length_chunk_corrected*m_grid*8.0          
				endif
			enddo

  enddo

	endif

	close(10)

  enddo
      
      do i=1,dim

          data_mass_fof(i) = total_fof_chunk(i) / data_fof_chunk(i)
          data_mass_fof_cor(i) = total_fof_chunk_cor(i) / data_fof_chunk_cor(i)

      enddo      
  
      ! calculate errors, NB must do this before volume correction!

  do i = 1, dim
  
    error_fof_chunk_low(i) = (real(data_fof_chunk(i))+0.25)**0.5 - 0.5
    error_fof_chunk_high(i) = (real(data_fof_chunk(i))+0.25)**0.5 + 0.5
    error_fof_chunk_cor_low(i) = (real(data_fof_chunk_cor(i))+0.25)**0.5 - 0.5
    error_fof_chunk_cor_high(i) = (real(data_fof_chunk_cor(i))+0.25)**0.5 + 0.5
    
  enddo

  ! interpolate the sigma values to get the value of sigma at the bin average:
  
  if(volume_correction .eq. 0) then
  
  do i=1,dim+1
		call splint(mms_fine(1:dim_ps),sigma,sigma_bin_ave_del,dim_ps,data_mass_fof(i),sigma_bin_ave(i))
		call splint(mms_fine(1:dim_ps),sigma_z,sigma_bin_ave_del_z,dim_ps,data_mass_fof(i),sigma_bin_ave_z(i))
	end do
  
  endif
  
  if(volume_correction .eq. 1) then
  do i = 1, dim
  call splint(mms_fine(1:sigma_limit),sigma_local_sorted_by_mass,sigma_local_del,sigma_limit,&
  &data_mass_fof(i),sigma_bin_ave(i))
  enddo
    ! convert halo counts to account for finite volume
    
  do i = 1, dim
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_bin_ave(i),mass1)
    data_fof_chunk(i) = data_fof_chunk_cor(i)*data_mass_fof(i)/mass1
    error_fof_chunk_low(i) = error_fof_chunk_low(i)*data_mass_fof(i)/mass1
    error_fof_chunk_high(i) = error_fof_chunk_high(i)*data_mass_fof(i)/mass1


    
    data_mass_fof(i) = mass1

  enddo

  ! convert bin masses to account for finite volume
  
  do i = 1, dim
  
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_low(i),mass1)
    data_mass_low(i) = mass1
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_high(i),mass1)
    data_mass_high(i) = mass1
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_ave(i),mass1)
    data_mass_ave(i) = mass1
  
  enddo

  endif

  	! Convert Counts of halos into log(M^2/rho_bar * dN/dM) and normalise to unit volume:

	do i=1,dim
  
		data_fof_chunk(i) = (data_fof_chunk(i) / (box/h)**3 * data_mass_fof(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))                

		error_fof_chunk_low(i) = (error_fof_chunk_low(i) / (box/h)**3 * data_mass_fof(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))                

		error_fof_chunk_high(i) = (error_fof_chunk_high(i) / (box/h)**3 * data_mass_fof(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))                
                
	enddo
  
  
    ! reset bins for data:
		
	do i=1,dim
															
		data_mass_high(i) = m_min*10.**(mass_bin_fac_coarse*(real(i))/bin_dim) !mass bins for plotting
			
		if(i==1) then
			data_mass_low(i) = m_min-1
		else
			data_mass_low(i) = data_mass_high(i-1)		
		endif
		
		data_mass_ave(i) = data_mass_low(i) + (data_mass_high(i) - data_mass_low(i))/2
    
		
	enddo

  
  
  ! convert to f(sigma_coarse)
  
  do i = 1, dim
  
		data_fof_chunk(i) = data_fof_chunk(i) / slope_coarse(i)    
		error_fof_chunk_low(i) = error_fof_chunk_low(i) / slope_coarse(i)
		error_fof_chunk_high(i) = error_fof_chunk_high(i) / slope_coarse(i)
    
    sigma_bin_ave_z(i) = sigma_bin_ave(i) / growth(omega0,lambda0,z) 													     

  enddo

	! Print data into a file:
	
    open (unit=45,file=fof_file,form='formatted')

    write(45,'(a218)') '#Mass/Msun (bin ave) #ln(sigma^-1) (bin ave)  #Mass/Msun (bin centre) #ln(sigma^-1) (Bin centre)&
      &  #particles (bin ave) #particles (bin centre) #particles (bin low cut) #ln(f(sigma)    #error low        #error high'
	   
     	   do i=1, dim   
     
          write(45,'(4es20.6,3i20,3f25.6)') data_mass_fof(i), log(1.0/sigma_bin_ave_z(i)),&
           &data_mass_ave(i), log(1.0/sigma_z_coarse(i)),&
           &int(data_mass_fof(i)/m_grid/8.0),int(data_mass_ave(i)/m_grid/8.0),int(data_mass_low(i)/m_grid/8.0),&
           &log(data_fof_chunk(i)),log(data_fof_chunk(i)-error_fof_chunk_low(i)),&
           &log(data_fof_chunk(i)+error_fof_chunk_high(i))
     
         enddo
     
	  close(45)
    


    ! Now calculate residuals. For this we need to interpolate the values of the analytic 
    ! functions at the precise values of sigma for which we have numerical data
    

  	do i=1,dim
    
  		call splint(mms_fine_ave,f_PS,data_ps_del,dim_ps,data_mass_fof(i),data_ps_int(i))
  		call splint(mms_fine_ave,f_ST,data_st_del,dim_ps,data_mass_fof(i),data_st_int(i))	  
  		call splint(mms_fine_ave,f_jenkins,data_jen_del,dim_ps,data_mass_fof(i),data_jen_int(i))
  		call splint(mms_fine_ave,f_warren,data_war_del,dim_ps,data_mass_fof(i),data_war_int(i))
      call splint(mms_fine_ave,f_tinker,data_tin_del,dim_ps,data_mass_fof(i),data_tin_int(i))
      call splint(mms_fine_ave,f_reed_2003,data_r03_del,dim_ps,data_mass_fof(i),data_r03_int(i))
      call splint(mms_fine_ave,f_reed_2007,data_r07_del,dim_ps,data_mass_fof(i),data_r07_int(i))
      call splint(mms_fine_ave,f_angulo,data_ang_del,dim_ps,data_mass_fof(i),data_ang_int(i))
      call splint(mms_fine_ave,f_crocce,data_cro_del,dim_ps,data_mass_fof(i),data_cro_int(i))
      call splint(mms_fine_ave,f_courtin,data_cou_del,dim_ps,data_mass_fof(i),data_cou_int(i))
      call splint(mms_fine_ave,f_watson_fof_uni,data_watson_fof_uni_del,dim_ps,data_mass_fof(i),data_watson_fof_uni_int(i))
      call splint(mms_fine_ave,f_watson_cpm_ahf,data_watson_cpm_ahf_del,dim_ps,data_mass_fof(i),data_watson_cpm_ahf_int(i))
                                                        
    end do
    
    open (unit=46,file=fof_dfile,form='formatted')        
    
    
    write(46,'(a1000)') '#Mass/Msun (bin ave)   #ln(sigma^-1) #Mass/Msun (bin centre) #ln(sigma^-1) (Bin centre)&
    & #particles (bin ave) #particles (bin centre) #particles (bin low cut) &
    &   #(data-PS)/PS             #error PS low&
    &            #error PS high           #(data-ST)/ST             #error ST low          #error ST high &
    &         #(data)/Jen            #error Jen low           # error Jen high &
    &         #(data)/War         #error War low           # error War high &
    &         #(data)/R03         #error R03 low           # error R03 high &
    &         #(data)/War         #error R07 low           # error R07 high &    
    &         #(data)/Tin          #error Tin low           # error Tin high    &
    &         #(data)/Ang         #error Ang low           # error Ang high      &
    &         #(data)/Cro         #error Cro low           # error Cro high      &
    &         #(data)/Cou         #error Cou low           # error Cou high      &
    &         #(data)/Wat_fof_uni         #error Wat_fof_uni low           # error Wat_fof_uni high          '    
    
   do i = 1, dim
    
      write(46, '(4es18.10,3i20,36f15.10)') data_mass_fof(i), log(1.0/sigma_bin_ave_z(i)),&
      &data_mass_ave(i), log(1.0/sigma_z_coarse(i)),&
      &int(data_mass_fof(i)/m_grid/8.0),int(data_mass_ave(i)/m_grid/8.0),int(data_mass_low(i)/m_grid/8.0),&
      &(data_fof_chunk(i))/data_ps_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_ps_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_ps_int(i),&
      &(data_fof_chunk(i))/data_st_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_st_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_st_int(i),&
      &(data_fof_chunk(i))/data_jen_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_jen_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_jen_int(i),&
      &(data_fof_chunk(i))/data_war_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_war_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_war_int(i),&
      &(data_fof_chunk(i))/data_r03_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_r03_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_r03_int(i),&
      &(data_fof_chunk(i))/data_r07_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_r07_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_r07_int(i),&
      &(data_fof_chunk(i))/data_tin_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_tin_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_tin_int(i),&
      &(data_fof_chunk(i))/data_ang_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_ang_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_ang_int(i),&
      &(data_fof_chunk(i))/data_cro_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_cro_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_cro_int(i),&
      &(data_fof_chunk(i))/data_cou_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_cou_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_cou_int(i),&
      &(data_fof_chunk(i))/data_watson_fof_uni_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_watson_fof_uni_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_watson_fof_uni_int(i),&
      &(data_fof_chunk(i))/data_watson_cpm_ahf_int(i),&
      &(data_fof_chunk(i)-error_fof_chunk_low(i))/data_watson_cpm_ahf_int(i),&
      &(data_fof_chunk(i)+error_fof_chunk_high(i))/data_watson_cpm_ahf_int(i)

    
    enddo

  close(46)
  

  ! interpolate the sigma values to get the value of sigma at the bin average:
  
  if(volume_correction .eq. 0) then
  
  do i=1,dim+1
		call splint(mms_fine(1:dim_ps),sigma,sigma_bin_ave_del,dim_ps,data_mass_fof_cor(i),sigma_bin_ave(i))
		call splint(mms_fine(1:dim_ps),sigma_z,sigma_bin_ave_del_z,dim_ps,data_mass_fof_cor(i),sigma_bin_ave_z(i))
	end do
  
  endif
  
  if(volume_correction .eq. 1) then
  do i = 1, dim
  call splint(mms_fine(1:sigma_limit),sigma_local_sorted_by_mass,sigma_local_del,sigma_limit,&
  &data_mass_fof_cor(i),sigma_bin_ave(i))
  enddo
 
  
    ! convert halo counts to account for finite volume
    
  do i = 1, dim
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_bin_ave(i),mass1)
    data_fof_chunk_cor(i) = data_fof_chunk_cor(i)*data_mass_fof_cor(i)/mass1
    error_fof_chunk_cor_low(i) = error_fof_chunk_cor_low(i)*data_mass_fof_cor(i)/mass1
    error_fof_chunk_cor_high(i) = error_fof_chunk_cor_high(i)*data_mass_fof_cor(i)/mass1

 

    data_mass_fof_cor(i) = mass1

  enddo

  ! convert bin masses to account for finite volume
  
  do i = 1, dim
  
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_low(i),mass1)
    data_mass_low(i) = mass1
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_high(i),mass1)
    data_mass_high(i) = mass1
    call splint(sigma1_2,mass_2,mass_cobe,max_size,sigma_mass_ave(i),mass1)
    data_mass_ave(i) = mass1
  
  enddo


  endif
    
  
  	! Convert Counts of halos into log(M^2/rho_bar * dN/dM) and normalise to unit volume:

	do i=1,dim
  
		data_fof_chunk_cor(i) = (data_fof_chunk_cor(i) / (box/h)**3 * data_mass_fof_cor(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))                

		error_fof_chunk_cor_low(i) = (error_fof_chunk_cor_low(i) / (box/h)**3 * data_mass_fof_cor(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))                
                
		error_fof_chunk_cor_high(i) = (error_fof_chunk_cor_high(i) / (box/h)**3 * data_mass_fof_cor(i)**2 / rho_bar / &
								& (data_mass_high(i) - data_mass_low(i)))                                
                
	enddo
  
      ! reset bins for data:
		
	do i=1,dim
															
		data_mass_high(i) = m_min*10.**(mass_bin_fac_coarse*(real(i))/bin_dim) !mass bins for plotting
			
		if(i==1) then
			data_mass_low(i) = m_min-1
		else
			data_mass_low(i) = data_mass_high(i-1)		
		endif
		
		data_mass_ave(i) = data_mass_low(i) + (data_mass_high(i) - data_mass_low(i))/2
    
		
	enddo

  
  
  ! convert to f(sigma_coarse)
  
  do i = 1, dim
  
		data_fof_chunk_cor(i) = data_fof_chunk_cor(i) / slope_coarse(i)
    error_fof_chunk_cor_low(i) = error_fof_chunk_cor_low(i) / slope_coarse(i)
		error_fof_chunk_cor_high(i) = error_fof_chunk_cor_high(i) / slope_coarse(i)
    
    sigma_bin_ave_z(i) = sigma_bin_ave(i) / growth(omega0,lambda0,z) 													     

  enddo

  
  
  
	! Print data into a file:
	
    open (unit=45,file=fofcor_file,form='formatted')

	
    write(45,'(a218)') '#Mass/Msun (bin ave) #ln(sigma^-1) (bin ave)  #Mass/Msun (bin centre) #ln(sigma^-1) (Bin centre)&
      &  #particles (bin ave) #particles (bin centre) #particles (bin low cut) #ln(f(sigma)    #error low        #error high'
	   
     	   do i=1, dim   
     
          write(45,'(4es20.6,3i20,3f25.6)') data_mass_fof_cor(i), log(1.0/sigma_bin_ave_z(i)),&
           &data_mass_ave(i), log(1.0/sigma_z_coarse(i)),&
           &int(data_mass_fof_cor(i)/m_grid/8.0),int(data_mass_ave(i)/m_grid/8.0),int(data_mass_low(i)/m_grid/8.0),&
           &log(data_fof_chunk_cor(i)),log(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i)),&
           &log(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))
     
         enddo
     
	  close(45)
    


    ! Now calculate residuals. For this we need to interpolate the values of the analytic 
    ! functions at the precise values of sigma for which we have numerical data
    
  	do i=1,dim
    
    
  		call splint(mms_fine_ave,f_PS,data_ps_del,dim_ps,data_mass_fof_cor(i),data_ps_int(i))
  		call splint(mms_fine_ave,f_ST,data_st_del,dim_ps,data_mass_fof_cor(i),data_st_int(i))	  
  		call splint(mms_fine_ave,f_jenkins,data_jen_del,dim_ps,data_mass_fof_cor(i),data_jen_int(i))
  		call splint(mms_fine_ave,f_warren,data_war_del,dim_ps,data_mass_fof_cor(i),data_war_int(i))
      call splint(mms_fine_ave,f_tinker,data_tin_del,dim_ps,data_mass_fof_cor(i),data_tin_int(i))
      call splint(mms_fine_ave,f_reed_2003,data_r03_del,dim_ps,data_mass_fof_cor(i),data_r03_int(i))
      call splint(mms_fine_ave,f_reed_2007,data_r07_del,dim_ps,data_mass_fof_cor(i),data_r07_int(i))
      call splint(mms_fine_ave,f_angulo,data_ang_del,dim_ps,data_mass_fof_cor(i),data_ang_int(i))
      call splint(mms_fine_ave,f_crocce,data_cro_del,dim_ps,data_mass_fof_cor(i),data_cro_int(i))
      call splint(mms_fine_ave,f_courtin,data_cou_del,dim_ps,data_mass_fof_cor(i),data_cou_int(i))
      call splint(mms_fine_ave,f_watson_fof_uni,data_watson_fof_uni_del,dim_ps,data_mass_fof_cor(i),data_watson_fof_uni_int(i))
      call splint(mms_fine_ave,f_watson_cpm_ahf,data_watson_cpm_ahf_del,dim_ps,data_mass_fof_cor(i),data_watson_cpm_ahf_int(i))
    
                                                    
    end do
    
    open (unit=46,file=fofcor_dfile,form='formatted')        
    
    
    write(46,'(a1000)') '#Mass/Msun (bin ave)   #ln(sigma^-1) #Mass/Msun (bin centre) #ln(sigma^-1) (Bin centre)&
    & #particles (bin ave) #particles (bin centre) #particles (bin low cut) &
    &   #(data-PS)/PS             #error PS low&
    &            #error PS high           #(data-ST)/ST             #error ST low          #error ST high &
    &         #(data)/Jen            #error Jen low           # error Jen high &
    &         #(data)/War         #error War low           # error War high &
    &         #(data)/R03         #error R03 low           # error R03 high &
    &         #(data)/War         #error R07 low           # error R07 high &    
    &         #(data)/Tin          #error Tin low           # error Tin high    &
    &         #(data)/Ang         #error Ang low           # error Ang high      &
    &         #(data)/Cro         #error Cro low           # error Cro high      &
    &         #(data)/Cou         #error Cou low           # error Cou high      &
    &         #(data)/Wat_fof_uni         #error Wat_fof_uni low           # error Wat_fof_uni high          '        

   do i = 1, dim
    
      write(46, '(4es18.10,3i20,36f15.10)') data_mass_fof_cor(i), log(1.0/sigma_bin_ave_z(i)),&
      &data_mass_ave(i), log(1.0/sigma_z_coarse(i)),&
      &int(data_mass_fof_cor(i)/m_grid/8.0),int(data_mass_ave(i)/m_grid/8.0),int(data_mass_low(i)/m_grid/8.0),&
      &(data_fof_chunk_cor(i))/data_ps_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_ps_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_ps_int(i),&
      &(data_fof_chunk_cor(i))/data_st_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_st_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_st_int(i),&
      &(data_fof_chunk_cor(i))/data_jen_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_jen_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_jen_int(i),&
      &(data_fof_chunk_cor(i))/data_war_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_war_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_war_int(i),&
      &(data_fof_chunk_cor(i))/data_r03_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_r03_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_r03_int(i),&
      &(data_fof_chunk_cor(i))/data_r07_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_r07_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_r07_int(i),&
      &(data_fof_chunk_cor(i))/data_tin_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_tin_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_tin_int(i),&
      &(data_fof_chunk_cor(i))/data_ang_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_ang_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_ang_int(i),&
      &(data_fof_chunk_cor(i))/data_cro_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_cro_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_cro_int(i),&
      &(data_fof_chunk_cor(i))/data_cou_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_cou_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_cou_int(i),&
      &(data_fof_chunk_cor(i))/data_watson_fof_uni_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_watson_fof_uni_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_watson_fof_uni_int(i),&
      &(data_fof_chunk_cor(i))/data_watson_cpm_ahf_int(i),&
      &(data_fof_chunk_cor(i)-error_fof_chunk_cor_low(i))/data_watson_cpm_ahf_int(i),&
      &(data_fof_chunk_cor(i)+error_fof_chunk_cor_high(i))/data_watson_cpm_ahf_int(i)
      

    enddo

  close(46)
  
  endif ! analytic only
  
  
if(volume_correction .eq. 1) then  

     deallocate(sigma_local)
     deallocate(sigma_local_sorted_by_mass)
     deallocate(sigma_local_del)
     deallocate(mass_local_del)
     deallocate(mms_fine_sorted_by_sigma)

endif

!----------------------------------------------------------------------------------------------------------------------------


end program

