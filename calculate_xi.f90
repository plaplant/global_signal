program main
  ! Intel
  use IFPORT
  use OMP_LIB


  ! FGSL
  use FGSL


  ! Global
  use xi_global


  ! Tools
  use io_tools
  use xi_tools


  ! Initialize
  call OMP_SET_NUM_THREADS(Ncpu)
  default_errh = fgsl_set_error_handler_off()
  call set_baseline


  ! Do work
  call init_maps
  call read_data
  call calc_II_map(all_maps, ii_maps)
  call compute_xi(ii_maps, xi_nu)
  call write_xi


end program main
