program main
  ! Intel
  use IFPORT
  use OMP_LIB


  ! Global
  use xi_global


  ! Tools
  use alm_tools
  use fgsl
  use io_tools
  use xi_tools


  ! Default
  implicit none


  ! Do work
  if (.false.) then
     call init_maps
     call read_data
     call calc_II_map(all_maps, ii_maps)
     call test_data
  endif
  call test_special_functions


contains


  subroutine test_data
    ! Default
    implicit none


    ! Local variables
    integer(4)     :: i,j,k,l
    real(8)        :: x
    complex(8)     :: z
    character(100) :: fn


    ! Local arrays
    real(8),    dimension(N_pix)           :: mp
    complex(8), dimension(1,0:LMAX,0:LMAX) :: alm


    ! Read raw data from a few indices to compare with python
    i = 1
    j = 1
    k = 147713
    l = 151
    z = all_maps(i,j,k,l)
    write(*,'(a,i1,a,i1,a,i6,a,i3,a,2es16.8)') &
         "input_maps(",i,",",j,",",k,",",l,"): ",z

    i = 1
    j = 2
    k = 147715
    l = 126
    z = all_maps(i,j,k,l)
    write(*,'(a,i1,a,i1,a,i6,a,i3,a,2es16.8)') &
         "input_maps(",i,",",j,",",k,",",l,"): ",z

    ! Compare I -> I' maps
    k = 147713
    l = 151
    x = ii_maps(k,l)
    write(*,'(a,i6,a,i3,a,es16.8)') &
         "ii_maps(",k,",",l,"): ",x

    k = 147678
    l = 101
    x = ii_maps(k,l)
    write(*,'(a,i6,a,i3,a,es16.8)') &
         "ii_maps(",k,",",l,"): ",x

    ! Compare a_lm's of HEALPix transform
    mp = ii_maps(:,98)

    ! Save for comparing with python
    fn = 'test_map.txt'
    open(11,file=fn)
    do i=1,N_pix
       write(11,'(es20.12)') mp(i)
    enddo
    close(11)

    i = 0
    j = 0
    z = alm(1,i,j)
    write(*,'(a,i3,a,i3,a,2es16.8)') "alm(1,",i,",",j,"): ",z
    call map2alm(N_side,LMAX,LMAX,mp,alm)
    i = 7
    j = 3
    z = alm(1,i,j)
    write(*,'(a,i3,a,i3,a,2es16.8)') "alm(1,",i,",",j,"): ",z
    call map2alm(N_side,LMAX,LMAX,mp,alm)
    i = 128
    j = 99
    z = alm(1,i,j)
    write(*,'(a,i3,a,i3,a,2es16.8)') "alm(1,",i,",",j,"): ",z

  end subroutine test_data


  subroutine test_special_functions
    ! Default
    implicit none

    ! Local variables
    integer(4) :: l,m
    real(8)    :: jl,theta,phi,x
    complex(8) :: y

    ! Compare spherical harmonic values
    l     = 3
    m     = -2
    theta = 0.2D0
    phi   = 0.1D0
    y     = Ylm(l,m,theta,phi)
    write(*,'(a,i1,a,i1,a,f5.2,a,f5.2,a,2es16.8)') &
         "Ylm(",l,",",m,",",theta,",",phi,"):     ",y
    l     = 194
    m     = -85
    theta = 2.6D0
    phi   = 4.9D0
    y     = Ylm(l,m,theta,phi)
    write(*,'(a,i3,a,i3,a,f5.2,a,f5.2,a,2es16.8)') &
         "Ylm(",l,",",m,",",theta,",",phi,"): ",y
    l  = 12
    x  = 10.4
    jl = fgsl_sf_bessel_jsl(l,x)
    write(*,'(a,i2,a,f5.2,a,es16.8)') "jl(",l,",",x,"): ",jl


  end subroutine test_special_functions


end program main
