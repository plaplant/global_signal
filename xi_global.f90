module xi_global
  ! Intel
  use IFPORT
  use OMP_LIB


  ! FGSL
  use FGSL


  ! Default
  implicit none


!------------------------------------------------------------------------------!


  ! Parameters


  ! Physical constants
  real(8), parameter :: pi = acos(-1D0)
  real(8), parameter :: c_cgs = 2.99792458D10
  real(8), parameter :: c_mks = c_cgs/1D2


  ! Computation
  integer(4), parameter :: Ncpu = 16


  ! HEALPix parameters
  integer(4), parameter :: N_side  = 128
  integer(4), parameter :: LMAX    = 3*N_side - 1
  integer(4), parameter :: N_pix   = 12*N_side**2


  ! Rotation angle values
  logical, parameter :: rotate_maps = .true.
  real(8), parameter :: r_psi       = 0D0
  real(8), parameter :: r_theta     = -((90 + 30.72D0) * (pi/180D0))
  real(8), parameter :: r_phi       = 0D0


  ! Baseline angle values
  real(8),    parameter :: b_theta = pi/2
  real(8),    parameter :: phi_min = 0
  real(8),    parameter :: phi_max = 2*pi
  real(8),    parameter :: d_phi   = pi/12
  integer(4), parameter :: N_phi   = int((phi_max-phi_min)/d_phi)


  ! Data parameters
  integer(4), parameter :: N_freq = 201
  real(8),    parameter :: nu_min = 100
  real(8),    parameter :: d_nu   = 0.5
  integer(4), parameter :: N_bl   = 8
  real(8),    parameter, dimension(N_bl) :: r_bl = (/ 14.6D0, 25.3D0, 29.2D0, &
       38.6D0, 43.8D0, 50.6D0, 52.6D0, 58.4D0 /)


  ! File locations
  character(*), parameter :: infile = '/data4/paper/plaplant/beams/'&
       //'HERA_ijones.hdf5'
  character(*), parameter :: outdir = 'Output/HERA/'


  ! V_xx vs V_I
  logical, parameter :: make_stokes_vis = .false.


!------------------------------------------------------------------------------!


  ! Variables


  ! FGSL
  type(fgsl_error_handler_t) :: default_errh


  ! IO
  logical :: no_maps


!------------------------------------------------------------------------------!


  ! Arrays


  ! Input
  complex(8), dimension(:,:,:,:), allocatable :: all_maps
  real(8),    dimension(:,:,:,:), allocatable :: mueller_maps
  complex(8), dimension(:,:,:),   allocatable :: vis_maps


  ! Output
  real(8),    dimension(N_phi)               :: phi_vals
  real(8),    dimension(N_freq)              :: nu_vals
  real(8),    dimension(N_bl)                :: tauh_vals
  complex(8), dimension(N_freq,N_phi,4,N_bl) :: xi_nu


end module xi_global
