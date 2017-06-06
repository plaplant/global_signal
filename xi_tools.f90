module xi_tools
  ! Intel
  use IFPORT
  use OMP_LIB


  ! HEALPix
  use alm_tools


  ! FGSL
  use fgsl


  ! Global
  use xi_global


contains


  subroutine calc_II_map(input_map, output_map)
    ! Default
    implicit none


    ! Subroutine arguments
    real(8), dimension(:,:,:,:), intent(in)  :: input_map
    real(8), dimension(:,:),     intent(out) :: output_map


    ! Local variables
    integer(4) :: i,j,ipix,inu


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Allocate output array
    allocate(output_map(Npix,Nfreq))
    output_map = 0

    ! Compute entries for output map in parallel
    !$omp parallel do           &
    !$omp default(shared)       &
    !$omp private(i,j,ipix,inu)
    do inu=1,N_freq
       do ipix=1,N_pix
          do j=1,2
             do i=1,2
                output_map(ipix,inu) = output_map(ipix,inu) &
                     + abs(input_map(i,j,ipix,inu))**2
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do

    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called calc II map'
    return
  end subroutine calc_II_map


!------------------------------------------------------------------------------!


  subroutine compute_xi(input_map,xi_nu)
    ! Default
    implicit none


    ! Subroutine arguments
    real(8),    dimension(:,:), intent(in)  :: input_map
    complex(8), dimension(:,:), intent(out) :: xi_nu


    ! Local variables
    integer(4) :: inu,il,im
    real(8)    :: arg,jl
    complex(8) :: prefac,y,a,xi


    ! Local arrays
    complex(8), dimension(1,0:LMAX,0:LMAX) :: alm
    real(8),    dimension(N_pix)           :: mp


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Initialize frequency array
    do inu=1,N_freq
       xi_nu(1,inu) = Nu_min + (inu-1)*d_nu
    enddo

    ! Convert from MHz to Hz
    xi_nu = xi_nu*1D6

    ! Loop over frequencies
    do inu=1,N_freq
       ! Get frequency in Hz
       nu = xi_nu(1,inu)

       ! Get healpix map and convert to a_lm
       mp = input_map(:,inu)
       call map2alm(N_side,LMAX,LMAX,mp,alm)

       ! Compute xi value
       xi = 0
       !$omp parallel do         &
       !$omp default(shared)     &
       !$omp private(il,im,idx)  &
       !$omp private(arg,jl,y,a) &
       !$omp reduction(+:xi)
       do il=0,LMAX
          arg = 2*pi*TAUH*nu
          jl  = fgsl_sf_bessel_jsl(il, arg)

          do im=0,LMAX
             prefac = (0D0, 1D0)**(3*il + 2*im)
             y      = Ylm(il,im,0D0,0D0)
             a      = alm(1,il,im)

             xi = xi + prefac*alm*jl*y
             if (im > 0) then
                ! for negative m, we have
                ! a_{l, -m} = a_{l, m}^*
                ! Y_l^{-m}  = (-1)**m * (Y_l^m)^*
                prefac = prefac * (-1)**im
                xi     = xi + conjg(alm)*jl*conjg(y)
             endif
          enddo
       enddo
       !$omp end parallel do

       ! Save value
       xi_nu(2,inu) = xi
    enddo


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called compute xi'
    return
  end subroutine compute_xi


!------------------------------------------------------------------------------!


  function Ylm(l, m, theta, phi)
    ! Default
    implicit none

    ! Function arguments
    integer(4), intent(in) :: l,m
    real(8),    intent(in) :: theta,phi
    complex(8)             :: Ylm

    ! Local variables
    real(8) :: x
    x   = fgsl_sf_legendre_sphplm(l, abs(m), cos(theta))
    Ylm = (-1)**m * x * exp((0D0,1D0)*abs(m)*phi)
    if (m < 0) then
       ! use symmetry to compute -m values
       Ylm = (-1)**m * conjg(Ylm)
    endif

    return
  end function Ylm


end module xi_tools
