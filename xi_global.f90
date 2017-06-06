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


  ! File locations
  character(*), parameter :: infile = '/data4/paper/plaplant/beams/'&
       'HERA_ijones.hdf5'


!------------------------------------------------------------------------------!


  ! Arrays


  ! Output
  complex(8), dimension(2,N_freq) :: xi_nu


end module xi_global
