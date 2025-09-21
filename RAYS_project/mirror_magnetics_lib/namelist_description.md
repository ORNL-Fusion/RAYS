
## Namelist data for module: mirror_magnetics_m

### namelist /mirror_magnetics_list/

<br>    character(len=80) :: coil_data_file = ''

<br>    character(len=80) :: current_data_file = ''

<br>    character(len=80) :: case_name = ''

<br>    character(len=80) :: base_file_name = ''

<br>    character(len=100) :: NC_file_name = ''

<br>	integer :: n_r, n_z ! Size of R,Z grid

<br>    real(KIND=rkind)  :: r_min = zero ! Can be reset in namelist

<br>    real(KIND=rkind)  :: r_max, z_min, z_max  ! Boundary of R,Z grid

<br>	integer :: n_coils

Location of scrape-off point of last flux surface
<br>    real(KIND=rkind) :: r_LUFS, z_LUFS


### namelist /coil_data_list/

<br>    character(len=80) :: coil_set_name = ''

<br>    real(KIND=rkind), allocatable :: inner_radius(:) ! Inner radius of coil conductors

<br>    real(KIND=rkind), allocatable :: outer_radius(:) ! Outer radius of coil conductors

<br>    real(KIND=rkind), allocatable :: z_width(:) ! z width of coils

<br>    real(KIND=rkind), allocatable :: z_center(:) ! z of coil center

<br>    integer(KIND=rkind), allocatable :: n_turns(:) ! number of turns in coil

<br>	integer, allocatable :: n_r_layers(:) ! number of turn layers in radius

<br>	integer, allocatable :: n_z_slices(:) ! number of turn slices in z


### namelist /current_data_list/

<br>    character(len=80) :: current_set_name = ''

<br>    real(KIND=rkind), allocatable :: I_coil(:) ! current in coil

