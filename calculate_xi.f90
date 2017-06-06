program main
  ! Intel
  use IFPORT
  use OMP_LIB


  ! Global
  use xi_global


  ! Tools
  use io_tools
  use xi_tools


  ! Do work
  call OMP_SET_NUM_THREADS(Ncpu)

  call read_data
  call calc_II_map
  call compute_xi
  call write_xi


end program main
