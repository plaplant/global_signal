module xi_global
  ! Intel
  use IFPORT
  use OMP_LIB


  ! Default
  implicit none


!------------------------------------------------------------------------------!


  ! Parameters


  ! Physical constants
  real(8), parameter :: pi = acos(-1D0)


  ! Computation
  integer(4), parameter :: Ncpu = 16


  ! HEALPix parameters
  integer(4), parameter :: N_side = 128
  integer(4), parameter :: LMAX   = 3*N_side - 1
  integer(4), parameter :: N_pix  = 12*N_side**2


  ! Data parameters
  integer(4), parameter :: N_freq = 201
  real(8),    parameter :: nu_min = 100
  real(8),    parameter :: d_nu   = 0.5
  real(8),    parameter :: TAUH   = 50D-9


  ! File locations
  character(*), parameter :: infile = '/data4/paper/plaplant/beams/'&
       'HERA_ijones.hdf5'


!------------------------------------------------------------------------------!


  ! Arrays


  ! Input
  real(8), dimension(:,:,:,:), allocatable :: all_maps
  real(8), dimension(:,:),     allocatable :: ii_maps


  ! Output
  complex(8), dimension(2,N_freq) :: xi_nu


end module xi_global
