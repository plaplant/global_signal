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


  ! Do work
  call init_maps
  call read_data
  if (make_stokes_vis) then
     call read_mueller_maps
     if (no_maps) then
        call calc_mueller_maps(all_maps, mueller_maps)
        call write_mueller_maps
     endif
     call compute_xi_mueller(mueller_maps, xi_nu)
  else
     call read_vis_maps
     if (no_maps) then
        call calc_vis_maps(all_maps, vis_maps)
        call write_vis_maps
     endif
     call compute_xi_vis(vis_maps, xi_nu)
  endif
  call write_xi


end program main
