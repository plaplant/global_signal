module io_tools
  ! Intel
  use IFPORT
  use OMP_LIB

  ! HDF5
  use ISO_C_BINDING
  use HDF5


  ! Global
  use xi_global


  ! Default
  implicit none


contains


  subroutine read_data
    ! Default
    implicit none


    ! Local variables
    character(100) :: obj_name
    integer(4)     :: error


    ! HDF5 variables
    integer(HID_T)  :: file_id,dset_id
    integer(HID_T)  :: dtype_id
    integer(SIZE_T) :: offset,o0
    type(C_PTR)     :: f_ptr


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Open HDF5 file
    call h5open_f(error)
    call h5fopen_f(infile, H5F_ACC_RDONLY_F, file_id, error)
    write(obj_name,'(a)') '/Data/MuellerMatrices'
    call h5dopen_f(file_id, obj_name, dset_id, error)

    ! Make a "complex derived type"
    offset = h5offsetof(c_loc(all_maps(1,1,1,1)), c_loc(all_maps(2,1,1,1)))
    o0     = 0
    call h5tcreate_f(H5T_COMPOUND_F, offset, dtype_id, error)
    call h5tinsert_f(dtype_id, "r", o0, H5T_NATIVE_DOUBLE, error)
    call h5tinsert_f(dtype_id, "i", offset/2, H5T_NATIVE_DOUBLE, error)

    ! read in data
    f_ptr = c_loc(all_maps)
    call h5dread_f(dset_id, dtype_id, f_ptr, error)

    ! Close it down
    call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called read data'
    return
  end subroutine read_data


!------------------------------------------------------------------------------!


  subroutine write_xi
    ! Default
    implicit none


    ! Local variables
    character(100)  :: fn,rot_dir
    type(C_PTR)     :: f_ptr
    integer(4)      :: error
    integer(HID_T)  :: file_id,dset_id
    integer(HID_T)  :: dspace_id,dtype_id
    integer(HID_T)  :: xi_id,attr_id,header_id
    integer(SIZE_T) :: hint,offset,c_size
    integer(HSIZE_T), dimension(1) :: adims
    integer(HSIZE_T), dimension(4) :: ddims


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Open file and write out
    if (rotate_maps) then
       rot_dir = 'beam_zenith'
    else
       rot_dir = 'beam_default'
    endif

    ! Make filename
    fn = outdir//trim(rot_dir)//'/xi_nu_phi.hdf5'

    write(*,*) "Writing ",trim(fn)

    call h5open_f(error)
    call h5fcreate_f(trim(fn), H5F_ACC_TRUNC_F, file_id, error)

    ! Create header
    hint = 0
    call h5gcreate_f(file_id, "/Header", header_id, error, size_hint=hint)

    ! Write attributes
    adims = (/ 1 /)

    ! Baseline theta value
    call h5screate_f(H5S_SCALAR_F, dset_id, error)
    call h5acreate_f(header_id, "b_theta", H5T_NATIVE_DOUBLE, dset_id, &
         attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, b_theta, adims, error)
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dset_id, error)

    ! Frequency information
    call h5screate_f(H5S_SCALAR_F, dset_id, error)
    call h5acreate_f(header_id, "N_freq", H5T_NATIVE_INTEGER, dset_id, &
         attr_id, error)
    f_ptr = c_loc(N_freq)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, f_ptr, error)
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dset_id, error)

    call h5screate_f(H5S_SCALAR_F, dset_id, error)
    call h5acreate_f(header_id, "nu_min", H5T_NATIVE_DOUBLE, dset_id, &
         attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, nu_min, adims, error)
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dset_id, error)

    call h5screate_f(H5S_SCALAR_F, dset_id, error)
    call h5acreate_f(header_id, "d_nu", H5T_NATIVE_DOUBLE, dset_id, &
         attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, d_nu, adims, error)
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dset_id, error)

    ! Write dataset
    hint = 0
    call h5gcreate_f(file_id, "/Data", xi_id, error, size_hint=hint)

    ! Write out phi values
    adims = (/ N_phi /)
    call h5screate_simple_f(1, adims, dspace_id, error)
    call h5dcreate_f(xi_id, "phi", H5T_NATIVE_DOUBLE, dspace_id, dset_id, &
         error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, phi_vals, adims, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    ! Write out frequencies
    adims = (/ N_freq /)
    call h5screate_simple_f(1, adims, dspace_id, error)
    call h5dcreate_f(xi_id, "nu", H5T_NATIVE_DOUBLE, dspace_id, dset_id, &
         error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, nu_vals, adims, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    ! Write out tauh values
    adims = (/ N_bl /)
    call h5screate_simple_f(1, adims, dspace_id, error)
    call h5dcreate_f(xi_id, "tauh", H5T_NATIVE_DOUBLE, dspace_id, dset_id, &
         error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tauh_vals, adims, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    ! Make custom complex type
    c_size = h5offsetof(c_loc(xi_nu(1,1,1,1)), c_loc(xi_nu(2,1,1,1)))
    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, c_size, dtype_id, error)
    call h5tinsert_f(dtype_id, "r", offset, H5T_NATIVE_DOUBLE, error)
    offset = offset + c_size/2
    call h5tinsert_f(dtype_id, "i", offset, H5T_NATIVE_DOUBLE, error)

    ! Write out the data
    ddims = (/ N_freq, N_phi, 4, N_bl /)
    call h5screate_simple_f(4, ddims, dspace_id, error)
    call h5dcreate_f(xi_id, "xi", dtype_id, dspace_id, dset_id, error)
    f_ptr = c_loc(xi_nu)
    call h5dwrite_f(dset_id, dtype_id, f_ptr, error)
    call h5tclose_f(dtype_id, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)


    ! Close it down
    call h5gclose_f(xi_id, error)
    call h5gclose_f(header_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called write xi'
    return
  end subroutine write_xi


end module io_tools
