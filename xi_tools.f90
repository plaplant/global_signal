module xi_tools
  ! Intel
  use IFPORT
  use OMP_LIB


  ! HEALPix
  use alm_tools


  ! FGSL
  use fgsl
  use, intrinsic :: iso_c_binding


  ! Global
  use xi_global


contains


  subroutine init_maps
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k,l


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Allocate
    allocate(all_maps(2,2,N_pix,N_freq))
    allocate(ii_maps(N_pix,N_freq))

    ! First touch in parallel
    !$omp parallel         &
    !$omp default(shared)  &
    !$omp private(i,j,k,l)
    !$omp do
    do l=1,N_freq
       do k=1,N_pix
          do j=1,2
             do i=1,2
                all_maps(i,j,k,l) = 0
             enddo
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do j=1,N_freq
       do i=1,N_pix
          ii_maps(i,j) = 0
       enddo
    enddo
    !$omp end do
    !$omp end parallel


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called init maps'
    return
  end subroutine init_maps


!------------------------------------------------------------------------------!


  subroutine calc_II_map(input_map, output_map)
    ! Default
    implicit none


    ! Subroutine arguments
    complex(8), dimension(:,:,:,:), intent(in)  :: input_map
    real(8),    dimension(:,:),     intent(out) :: output_map


    ! Local variables
    integer(4) :: i,j,ipix,inu


    ! Local arrays
    real(8),    dimension(N_pix)           :: mp
    complex(8), dimension(1,0:LMAX,0:LMAX) :: alm


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Initialize output array
    output_map = 0

    ! Compute entries for output map in parallel
    !$omp parallel do           &
    !$omp default(shared)       &
    !$omp private(i,j,ipix,inu) &
    !$omp private(mp,alm)
    do inu=1,N_freq
       do ipix=1,N_pix
          do j=1,2
             do i=1,2
                output_map(ipix,inu) = output_map(ipix,inu) &
                     + abs(input_map(i,j,ipix,inu))**2
             enddo
          enddo
       enddo

       if (rotate_maps) then
          ! Also rotate map
          mp = output_map(:,inu)
          call map2alm(N_side,LMAX,LMAX,mp,alm)
          call rotate_alm(LMAX,alm,r_psi,r_theta,r_phi)
          call alm2map(N_side,LMAX,LMAX,alm,mp)
          output_map(:,inu) = mp
       endif
    enddo
    !$omp end parallel do


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called calc II map'
    return
  end subroutine calc_II_map


!------------------------------------------------------------------------------!


  subroutine compute_xi(input_map,xi)
    ! Default
    implicit none


    ! Subroutine arguments
    real(8),    dimension(:,:), intent(in)  :: input_map
    complex(8), dimension(:,:), intent(out) :: xi


    ! Local variables
    integer(4) :: inu,il,im,iphi
    real(8)    :: arg,jl,nu,theta,phi
    complex(8) :: prefac,y,a,x


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
       nu_vals(inu) = Nu_min + (inu-1)*d_nu
    enddo

    ! Convert from MHz to Hz
    nu_vals = nu_vals*1D6

    ! Loop over phi values
    do iphi=1,N_phi
       theta = b_theta
       phi   = phi_min + (iphi-1)*d_phi
       phi_vals(iphi) = phi

       ! Write out as a progress report
       write(*,*) "theta, phi: ", theta, phi

       ! Loop over frequencies
       do inu=1,N_freq
          ! Get frequency in Hz
          nu  = nu_vals(inu)
          arg = 2*pi*TAUH*nu

          ! Get healpix map and convert to a_lm
          mp = input_map(:,inu)
          call map2alm(N_side,LMAX,LMAX,mp,alm)

          ! Compute xi value
          x = 0
          !$omp parallel do         &
          !$omp default(shared)     &
          !$omp private(il,jl,im)   &
          !$omp private(y,a,prefac) &
          !$omp reduction(+:x)
          do il=0,LMAX
             jl  = fgsl_sf_bessel_jsl(il, arg)

             do im=0,il
                prefac = (0D0, 1D0)**(3*il + 2*im)
                y      = Ylm(il,im,theta,phi)
                a      = alm(1,il,im)

                x = x + prefac*a*jl*y
                if (im > 0) then
                   ! for negative m, we have
                   ! a_{l, -m} = a_{l, m}^*
                   ! Y_l^{-m}  = (-1)**m * (Y_l^m)^*
                   prefac = prefac * (-1)**im
                   x      = x + prefac*conjg(a)*jl*conjg(y)
                endif
             enddo
          enddo
          !$omp end parallel do

          ! Save value
          xi(inu,iphi) = x
       enddo
    enddo

    ! Add overall normalization
    xi = xi*sqrt(4*pi)


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
    real(8)    :: x
    integer(4) :: am

    am  = abs(m)
    x   = fgsl_sf_legendre_sphplm(l, am, cos(theta))
    Ylm = (-1)**m * x * exp((0D0,1D0)*am*phi)
    if (m < 0) then
       ! use symmetry to compute -m values
       Ylm = (-1)**m * conjg(Ylm)
    endif

    return
  end function Ylm


end module xi_tools
