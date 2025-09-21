
## Namelist data for module: axisym_toroid_processor_m

### namelist /axisym_toroid_processor_list/

Number of k vectors to plot for each ray in graphics
<br>    integer :: num_plot_k_vectors

Scale plot k vectors to kmax on the ray, (True, False)
<br>    character(len = 5) :: scale_k_vec = 'True'

Base length for drawing k vectors as fraction of fig size
<br>    real(KIND=rkind) :: k_vec_base_length

Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
<br>    character(len = 5) :: set_XY_lim = 'True'

Number of plasma boundary points to calculate
<br>	integer :: n_boundary_points

R,Z boundary points
<br>    real(KIND=rkind), allocatable :: R_boundary(:), Z_boundary(:)

Number of points in R,Z grid for normalized psi and eq_RZ_grid netCDF files
<br>	integer ::  N_pointsR_eq, N_pointsZ_eq

Number of points in psiN grid for 1D profiles
<br>	integer ::  n_psiN

Number of points in rho grid for 1D profiles
<br>	integer ::  n_rho

Tolerance for finding R(psi) grid by bisection
<br>    real(KIND=rkind) :: bisection_eps

Flags determining what to do besides plot rays.  Can be overridden (i.e. turned off)
in namelist file
<br>	logical :: calculate_dep_profiles = .true. ! Calculate all depositon profiles

<br>	logical :: write_dep_profiles = .true. ! Write depositon profile netCDF file

<br>	logical :: calculate_ray_diag = .true. ! Write detailed ray diagnostic netCDF file

<br>	logical :: write_contour_data = .true.  ! Write data needed to plot contours netCDF file

<br>	logical :: write_eq_RZ_grid_data = .true.  ! Write data for equilibrium on RZ grid netCDF file

<br>	logical :: write_eq_radial_profile_data = .true.  ! Write data for radial profiles netCDF file


## Namelist data for module: deposition_profiles_m

### namelist /template_list/

Number of bins to distribute the profile
<br>    integer :: n_bins = 0

Switch whether to write list-directed ASCII file
<br>    logical :: write_results_list_directed = .true.


## Namelist data for module: mirror_processor_m

### namelist /mirror_processor_list/

Number of k vectors to plot for each ray in graphics
<br>    integer :: num_plot_k_vectors

Scale plot k vectors to kmax on the ray, (True, False)
<br>    character(len = 5) :: scale_k_vec = 'True'

Base length for drawing k vectors as fraction of fig size
<br>    real(KIND=rkind) :: k_vec_base_length

Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
<br>    character(len = 5) :: set_XY_lim = 'True'

Number of plasma boundary points to calculate
<br>	integer :: n_boundary_points

R,Z boundary points
<br>    real(KIND=rkind), allocatable :: R_boundary(:), Z_boundary(:)

Number of points in X,Z grid for normalized Aphi and eq_XZ_grid netCDF files
<br>	integer ::  N_pointsX_eq, N_pointsZ_eq

Number of points in AphiN grid for 1D profiles
<br>	integer ::  n_AphiN

Number of points in rho grid for 1D profiles
<br>	integer ::  n_rho

Tolerance for finding R(Aphi) grid by bisection
<br>    real(KIND=rkind) :: bisection_eps

Reference z location at which to evaluate radial equilibrium profiles
<br>    real(KIND=rkind) :: z_reference

<br>	character(len=25) :: char_z_ref

Radius of LUFS at z = z_reference
<br>    real(KIND=rkind) :: r_LUFS_at_z_ref

Flags determining what to do besides plot rays.  Can be overridden (i.e. turned off)
in namelist file
<br>	logical :: calculate_dep_profiles = .true. ! Calculate all depositon profiles

<br>	logical :: write_dep_profiles = .true. ! Write depositon profile netCDF file

<br>	logical :: calculate_ray_diag = .true. ! Write detailed ray diagnostic netCDF file

<br>	logical :: write_contour_data = .true.  ! Write data needed to plot contours netCDF file

<br>	logical :: write_eq_XZ_grid_data = .true.  ! Write data for equilibrium on RZ grid netCDF file

<br>	logical :: write_eq_radial_profile_data = .true.  ! Write data for radial profiles netCDF file

<br>	logical :: do_OX_conv_analysis = .false.  ! Special case OX conversion so, default = false


## Namelist data for module: post_processing_m

### namelist /post_process_list/

Switch to select specific post processor
<br>    character(len=80) :: processor = ''

Selector for ray data input mode -> LD, NC, or ASCII (obsolete but still here)
<br>    character(len=80) :: ray_data_input_mode = ''


## Namelist data for module: slab_processor_m

### namelist /slab_processor_list/

Number of k vectors to plot for each ray in graphics
<br>    integer :: num_plot_k_vectors

Scale plot k vectors to kmax on the ray, (True, False) N.B. character not logical -> python
<br>    character(len = 5) :: scale_k_vec = 'True'

Base length for drawing k vectors as fraction of fig size
<br>    real(KIND=rkind) :: k_vec_base_length

Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
<br>    character(len = 5) :: set_XY_lim = 'True'

number of x points in equilibrium profile plots
<br>	integer :: n_X = 51 ! Default can be reset on input

Flags determining what to do besides plot rays.  Can be overridden (i.e. turned off)
in namelist file
<br>	logical :: calculate_dep_profiles = .true. ! Calculate all depositon profiles

<br>	logical :: write_dep_profiles = .true. ! Write depositon profile netCDF file

<br>	logical :: calculate_ray_diag = .true. ! Write detailed ray diagnostic netCDF file

<br>	logical :: write_eq_X_profile_data = .true.  ! Write data for radial profiles netCDF file


## Namelist data for module: solovev_processor_m

### namelist /solovev_processor_list/

<br>    character(len=80) :: processor = ''

Number of k vectors to plot for each ray in graphics
<br>    integer :: num_plot_k_vectors

Scale plot k vectors to kmax on the ray, (True, False)
<br>    character(len = 5) :: scale_k_vec = 'True'

Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
<br>    character(len = 5) :: set_XY_lim = 'True'

