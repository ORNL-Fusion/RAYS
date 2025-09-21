
## Namelist data for module: axisym_toroid_eq_m

### namelist /axisym_toroid_eq_list/

data for magnetics
<br>    character(len=60) :: magnetics_model

Geometry data
Magnetic axis
<br>    real(KIND=rkind) :: r_axis, z_axis

data for bounding box
<br>    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax

data for plasma boundary
<br>    real(KIND=rkind) :: inner_bound, outer_bound, upper_bound, lower_bound

Maximum value of psi considered to be inside of plasma.  Defaults to 1.0 but can be
reset in namelist.
<br>    real(KIND=rkind) :: plasma_psi_limit = one

Data for density
<br>    character(len=60) :: density_prof_model

<br>    real(KIND=rkind) :: alphan1 ! parameters for parabolic model

<br>    real(KIND=rkind) :: alphan2

Density outside psi = 1 as a fraction of ne0, defaults to 0. but can be set in namelist
<br>    real(KIND=rkind) :: d_scrape_off = zero

Data for temperature
<br>    character(len=60), allocatable :: temperature_prof_model(:)

Parabolic model parameters
<br>    real(KIND=rkind), allocatable :: alphat1(:) ! Can be dfferent for different species

<br>    real(KIND=rkind), allocatable :: alphat2(:) ! Can be dfferent for different species

Temperature outside psi = 1 as a fraction of Te0, defaults to 0. but can be set in namelist
<br>    real(KIND=rkind) :: T_scrape_off = zero

<br>	integer :: i, is


## Namelist data for module: axisym_toroid_ray_init_R_Z_nphi_ntheta_m

### namelist /axisym_toroid_ray_init_R_Z_nphi_ntheta_list/

Data for initial launch position
N.B. R_launch <--> major radius, Z is relative to equatorial plane
integer:: n_R_launch = 1, n_Z_launch = 1

<br>    real(KIND=rkind) ::  R_launch0 = 0., Z_launch0 = 0.

Data initial launch refractive index: theta = poloidal direction, phi = toroidal direction
integer:: n_rindex_theta = 1

<br>    real(KIND=rkind) ::  rindex_theta0 = 0., delta_rindex_theta = 0.

integer:: n_rindex_phi = 1

<br>    real(KIND=rkind) ::  rindex_phi0 = 0., delta_rindex_phi = 0.


## Namelist data for module: damping_m

### namelist /damping_list/

Switch to select which model to use for damping calculation
damping_model = 'no_damp" do not calculate damping
damping_model = 'poynting' use multicomponent Poyntings theorem. (Not yet)
damping_model = 'fund_ECH' use simple weak damping approximation for fundamental ECH
<br>    character(len=60) :: damping_model

Multi species damping.  Only meaningful if damping_model /= 'no_damp"
If .true. integrate damping by individual species as well as total damping
If .false. integrate damping only total damping
<br>    logical :: multi_spec_damping = .false.

Ray is considered totally damped if damping > total_damping_limit. Can reset in input
<br>    real(KIND=rkind) :: total_damping_limit = 0.99


## Namelist data for module: density_spline_interp_m

### namelist /density_spline_interp_list/

<br>    character (len = 60) :: spline_density_model ! Not used, YET

<br>	integer :: ngrid ! Actual number of points to be splined. <= n_grid_max

<br>	real(KIND=rkind) ::  ne_in(n_grid_max) ! Values on grid (grid assumed uniform 0 to 1)

Stuff for 1D spline profiles
type(cube_spline_function_1D) :: ne_profile_N  ! ne profile normalized to 1. on axis

<br> 	character (len = 80) :: profile_name = 'ne_profile'


## Namelist data for module: eqdsk_magnetics_lin_interp_m

### namelist /eqdsk_magnetics_lin_interp_list/

Name of input eqdsk file
<br>    character (len = 100) :: eqdsk_file_name


## Namelist data for module: eqdsk_magnetics_spline_interp_m

### namelist /eqdsk_magnetics_spline_interp_list/

Name of input eqdsk file
<br>    character (len = 100) :: eqdsk_file_name


## Namelist data for module: equilibrium_m

### namelist /equilib_model/

Switch to select specific equilibrium model.
<br>    character(len = 15) :: equilib_model


## Namelist data for module: file_input_ray_init_m

### namelist /file_input_ray_init_list/

Number of initial condition sets to be read in from namelist
integer:: n_rays_in

Initial positions and directions to be read in from namelist file.
N.B. The ordering of indices in rvec_in(nray_max,3), rindex_vec_in(nray_max,3) is
opposite of that for rvec0(3, nray), rindex_vec0(3, nray)!!!  This is to make it
simpler to put data into the namelist file i.e. (X, Y, Z) on one line of the
namelist.
<br>    real(KIND=rkind), allocatable :: rvec_in(:,:), rindex_vec_in(:,:)

<br>    real(KIND=rkind), allocatable :: ray_pwr_wt_in(:)


## Namelist data for module: mirror_magnetics_spline_interp_m

### namelist /mirror_magnetics_spline_interp_list/

Output netCDF file used by multiple_mirror_eq_m
<br>    character (len = 100) :: mirror_field_NC_file


## Namelist data for module: multiple_mirror_eq_m

### namelist /multiple_mirror_eq_list/

<br>    character(len=60) :: magnetics_model

Maximum value of normalized Aphi considered to be inside of plasma.  Defaults to 1.0
but can be reset in namelist.
<br>    real(KIND=rkind) :: plasma_AphiN_limit = one

Data for density
<br>    character(len=60) :: density_prof_model

<br>    real(KIND=rkind) :: alphan1, alphan2 ! parameters for parabolic model

<br>    real(KIND=rkind) :: AphiN0_d, delta_d ! parameters for hyperbolic model

Density outside plasma_AphiN_limit = 1 as a fraction of ne0. Defaults to 0. but can
be reset in namelist
<br>    real(KIND=rkind) :: d_scrape_off = zero

Data for temperature
<br>    character(len=60), allocatable :: temperature_prof_model(:)

Parabolic model parameters
<br>    real(KIND=rkind), allocatable :: alphat1(:), alphat2(:) ! Can be dfferent for different species

<br>    real(KIND=rkind), allocatable :: Aphin0_t(:),delta_t(:) ! Can be dfferent for different species

Temperature outside rho = 1 as a fraction of Te0, defaults to 0. but can be set in namelist
<br>    real(KIND=rkind) :: T_scrape_off = zero


## Namelist data for module: ode_m

### namelist /ode_list/

Switch to select ODE solvers for ray tracing.  Presently supported solvers are:
ode_solver= SG_ode: subroutine ODE developed by L. F. Shampine and M. K. Gordon.
ode_solver = RK4_ode: Simple Runge-Kutta 4th order integrator.
<br>    character (len = 15) :: ode_solver_name

Name of routine used to calculate RHS derivatives for ray equations. Presently
supported routines are:
ray_deriv_name = cold
ray_deriv_name = numerical
<br>    character(len=60) :: ray_deriv_name

Maximum length of ray
<br>     real(KIND=rkind) :: s_max

ODE step size.
<br>    real(KIND=rkind) :: ds

Maximum no. of steps allowed.
<br>    integer :: nstep_max

Namelist

## Namelist data for module: one_ray_init_XYZ_k_direction_m

### namelist /one_ray_init_XYZ_k_direction_list/

Data for initial launch position
<br>    real(KIND=rkind) ::  X = zero, Y = zero, Z = zero

Data for initial launch direction
<br>    real(KIND=rkind) ::  nX = zero, nY = zero, nZ = zero

Switch to turn off solution of disp rel.  Just use the input nX,nY,nZ as is
<br>	logical :: use_this_n_vec = .false.


## Namelist data for module: openmp_m

### namelist /openmp_list/

Number of OMP threads to use
<br>    integer :: num_threads


## Namelist data for module: ray_init_m

### namelist /ray_init_list/

Name of the specific model used for ray initialization
<br>    character(len=60) :: ray_init_model

Maximum number of rays allowed.
<br>    integer :: nray_max


## Namelist data for module: ray_results_m

### namelist /ray_results_list/

<br>    logical :: write_results_list_directed = .false.

<br>    logical :: write_results_netCDF = .false.


## Namelist data for module: rf_m

### namelist /rf_list/

Name of dispersion model used in ray tracing
<br>    character(len=60) :: ray_dispersion_model

RF in Hz,
<br>    real(KIND=rkind) :: frf

A switch to select which root from the dispersion to be used for
k initialization: i.e., fast wave, slow wave, IBW, or KAW, etc.
Cold plasma gives quadratic for nx^2(nz).  So can choose wave mode = plus, minus
fast or slow, where fast is the smaller of |nx|
<br>    character(len=60) :: wave_mode

A switch to determine whether the rays are initialized to propagate into
or out of the plasma.  Typically k_init_sign = 1 for waves to propagate inward
(unless the wave is a backward wave).
integer:: k0_sign

ray_param = 'arcl': default. Integrate with respect to the arclength along the ray.
ray_param = 'time': Integrate with respect to time along the ray.
<br>    character(len = 4) :: ray_param = 'arcl'

Maximum allowable residual of dispersion function. Used in checksave()
<br>    real(KIND=rkind) :: dispersion_resid_limit


## Namelist data for module: SG_ode_m

### namelist /axisym_toroid_eq_list/

Target relative and absolute error tolerances for the ODE solver.
<br>    real(KIND=rkind) :: rel_err0, abs_err0

Total ODE error limit abs(rel_err)+abs(abs_err) above which to bail.
<br>    real(KIND=rkind) :: SG_error_limit = 0.1  ! Default


## Namelist data for module: simple_slab_ray_init_m

### namelist /simple_slab_ray_init_list/

integer:: n_x_launch = 1

<br>    real(KIND=rkind) ::  x_launch0 = zero, dx_launch = zero

integer:: n_y_launch = 1

<br>    real(KIND=rkind) ::  y_launch0 = zero, dy_launch = zero

integer:: n_z_launch = 1

<br>    real(KIND=rkind) ::  z_launch0 = zero, dz_launch = zero

integer:: n_ky_launch, n_kz_launch

<br>    real(KIND=rkind) ::  rindex_y0, delta_rindex_y0, rindex_z0, delta_rindex_z0


## Namelist data for module: slab_eq_m

### namelist /slab_eq_list/

Geometry data
data for bounding box (meters)
<br>    real(KIND=rkind) :: xmin, xmax, ymin, ymax, zmin, zmax

location of x = 0 for tokamak-like models
<br>    real(KIND=rkind) :: rmaj

minor radius-like scale length for parabolic profiles or Gaussian
<br>    real(KIND=rkind) :: rmin

center position of Linear_2 profiles
<br>    real(KIND=rkind) :: x0

data for slab magnetics
Model names for Bx, By, Bz
<br>    character(len=12) :: bx_prof_model, by_prof_model, bz_prof_model

Magnetic field in Tesla at x,y,z = 0
<br>    real(KIND=rkind) :: bx0, by0, bz0

Parameters for linear models of Bz and By shear
<br>    real(KIND=rkind) :: LBy_shear_scale, LBz_scale

Slope for Linear_2 model
<br>    real(KIND=rkind) :: dBzdx

data for slab density
Model name for density profiles
<br>    character(len=60) :: dens_prof_model

Parameters for linear model of ne
<br>    real(KIND=rkind) :: Ln_scale

Slope for Linear_2 model
<br>    real(KIND=rkind) :: dndx

Parameters for parabolic model
<br>    real(KIND=rkind) :: alphan1

<br>    real(KIND=rkind) :: alphan2

Density outside x = rmin as a fraction of ne0, defaults to 0. but can be set in namelist
<br>    real(KIND=rkind) :: n_min = zero

data for slab temperature
Model name for temperature profiles
<br>    character(len=20), allocatable :: t_prof_model(:)

Parameters for linear models of Te
<br>    real(KIND=rkind) :: LT_scale

Parameters for Linear_2 models
<br>    real(KIND=rkind) :: dtdx

Parameters for parabolic model
<br>    real(KIND=rkind), allocatable :: alphat1(:)

<br>   real(KIND=rkind), allocatable :: alphat2(:)

Temperature outside x = rmin as a fraction of T0s, defaults to 0. but can be set in namelist
<br>   real(KIND=rkind), allocatable :: T_min(:)


## Namelist data for module: solovev_eq_m

### namelist /solovev_eq_list/

data for magnetics
<br>    real(KIND=rkind) :: rmaj, kappa, bphi0, iota0

<br>    real(KIND=rkind) :: outer_bound

Flux function psi at plasma boundary
<br>    real(KIND=rkind) :: psiB

data for density and temperature
<br>    character(len=60) :: dens_prof_model

<br>    real(KIND=rkind) :: alphan1

<br>    real(KIND=rkind) :: alphan2

<br>    character(len=20), allocatable :: t_prof_model(:)

<br>    real(KIND=rkind), allocatable :: alphat1(:)

<br>    real(KIND=rkind), allocatable :: alphat2(:)

data for bounding box
<br>    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax


## Namelist data for module: solovev_magnetics_m

### namelist /solovev_magnetics_list/

data for magnetics
<br>    real(KIND=rkind) :: rmaj, kappa, bphi0, iota0

<br>    real(KIND=rkind) :: inner_bound, outer_bound, vert_bound, r_Zmax, zmax_sq

<br>    real(KIND=rkind) :: outer_boundary

data for bounding box
<br>    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax


## Namelist data for module: solovev_ray_init_nphi_ntheta_m

### namelist /solovev_ray_init_nphi_ktheta_list/

Data for initial launch position
N.B. r_launch <--> minor radius, theta is relative to equatorial plane
integer:: n_r_launch = 1

<br>    real(KIND=rkind) ::  r_launch0 = 0., dr_launch = 0.

integer:: n_theta_launch = 1

<br>    real(KIND=rkind) ::  theta_launch0 = 0., dtheta_launch = 0.

Data initial launch refractive index: theta = poloidal direction, phi = toroidal direction
integer:: n_rindex_theta = 1

<br>    real(KIND=rkind) ::  rindex_theta0 = 0., delta_rindex_theta = 0.

integer:: n_rindex_phi = 1

<br>    real(KIND=rkind) ::  rindex_phi0 = 0., delta_rindex_phi = 0.


## Namelist data for module: species_m

### namelist /species_list/

Electron density at reference point (e.g. at magnetic axis or peak electron density)
The profiles generated in the various equilibrium modules are normalized to one at
the reference location (i.e. where ne = n0s(0) = n0).  The ion densities are specified
as a fraction of electron density, eta(i).
<br>    real(KIND=rkind) :: n0

is = species number
is=0 is reserved for electrons and the rest for ions.
qs: charge of species is
ms: mass of species is
eta: concentration as fraction of electron density
t0s_eV: temperature in eV.  Entered in eV in namelist file for convenience
tseps_eV: temperature in eV.  Entered in eV in namelist file for convenience
alfas: T_perp/T_paral
v0s: parallel drift velocity
nus: collision frequency i.e. nu/omega
<br>    character(len=12), dimension(0:nspec0) :: spec_name = ''

<br>    real(KIND=rkind), dimension(0:nspec0) :: qs = 0.

<br>    real(KIND=rkind), dimension(0:nspec0) :: ms = 0.

<br>    real(KIND=rkind), dimension(0:nspec0) :: eta = 0.

<br>    real(KIND=rkind), dimension(0:nspec0) :: nseps = 0.

<br>    real(KIND=rkind), dimension(0:nspec0) :: t0s_eV = 0.

<br>    real(KIND=rkind), dimension(0:nspec0) :: tseps_eV = 0.

<br>    real(KIND=rkind), dimension(0:nspec0) :: alfas = 0.

<br>    real(KIND=rkind), dimension(0:nspec0) :: v0s = 0.

<br>    real(KIND=rkind), dimension(0:nspec0) :: nus = 0.

An array indicating which plasma dispersion model is to be used for each species
spec_model(is) = 'cold' susceptibility model is cold plasma
spec_model(is) = 'bessel' susceptibility model is full besssel function. Not yet.
<br>    character(len=12) :: spec_model(0:nspec0) = ''

Criterion for checking charge neutrality.  Default here can be reset on input.
<br>    real(KIND=rkind) :: neutrality = 1.e-10


## Namelist data for module: temperature_spline_interp_m

### namelist /temperature_spline_interp_list/

<br>    character (len = 60) :: spline_Te_model, spline_Ti_model ! Not used, YET

<br>	integer :: ngrid ! Actual number of points to be splined. <= n_grid_max

<br>	real(KIND=rkind) ::  Te_in(n_grid_max) ! Values on grid (grid assumed uniform 0 to 1)

<br>	real(KIND=rkind) ::  Ti_in(n_grid_max) ! Values on grid (grid assumed uniform 0 to 1)

