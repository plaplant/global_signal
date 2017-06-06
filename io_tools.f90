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
    character(100) :: fn
    integer(4)     :: i


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Open file and write out
    fn = 'xi_nu_fortran.txt'
    open(11,file=fn)
    write(11,'(a,a15,2a16)') '#','nu [MHz]','re(Xi)','im(Xi)'
    do i=1,N_freq
       write(11,'(3es16.8)') dble(xi_nu(1,i)/1D6), dble(xi_nu(2,i)), &
            dimag(xi_nu(2,i))
    enddo
    close(11)


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called write xi'
    return
  end subroutine write_xi


end module io_tools
