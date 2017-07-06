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


  ! Computation
  integer(4), parameter :: Ncpu = 8


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
  real(8),    parameter :: TAUH   = 50D-9


  ! File locations
  character(*), parameter :: infile = '/data4/paper/plaplant/beams/'&
       //'HERA_ijones.hdf5'
  character(*), parameter :: outdir = 'Output/'


!------------------------------------------------------------------------------!


  ! Variables


  ! FGSL
  type(fgsl_error_handler_t) :: default_errh


!------------------------------------------------------------------------------!


  ! Arrays


  ! Input
  complex(8), dimension(:,:,:,:), allocatable :: all_maps
  real(8),    dimension(:,:),     allocatable :: ii_maps


  ! Output
  real(8),    dimension(N_phi)        :: phi_vals
  real(8),    dimension(N_freq)       :: nu_vals
  complex(8), dimension(N_freq,N_phi) :: xi_nu


end module xi_global
