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
    if (make_stokes_vis) then
       allocate(mueller_maps(N_pix,N_freq,4,4))
    else
       allocate(vis_maps(N_pix,N_freq,4))
    endif

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
    if (make_stokes_vis) then
       !$omp do
       do l=1,4
          do k=1,4
             do j=1,N_freq
                do i=1,N_pix
                   mueller_maps(i,j,k,l) = 0
                enddo
             enddo
          enddo
       enddo
       !$omp end do
    else
       !$omp do
       do k=1,4
          do j=1,N_freq
             do i=1,N_pix
                vis_maps(i,j,k) = 0
             enddo
          enddo
       enddo
       !$omp end do
    endif
    !$omp end parallel


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called init maps'
    return
  end subroutine init_maps


!------------------------------------------------------------------------------!


  subroutine calc_mueller_maps(input_map, output_map)
    ! Default
    implicit none


    ! Subroutine arguments
    complex(8), dimension(:,:,:,:), intent(in)  :: input_map
    real(8),    dimension(:,:,:,:), intent(out) :: output_map


    ! Local variables
    integer(4) :: i,j,ipix,inu


    ! Local arrays
    real(8),    dimension(N_pix)           :: mp
    complex(8), dimension(2,2)             :: Pi,Pj,J_j
    complex(8), dimension(1,0:LMAX,0:LMAX) :: alm


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Initialize output array
    output_map = 0

    ! Compute Mueller matrix entries
    do j=1,4
       Pj = stokes_matrix(j)
       do i=1,4
          Pi = stokes_matrix(i)
          ! Parallelize over frequency
          !$omp parallel do           &
          !$omp default(shared)       &
          !$omp private(ipix,inu,J_j) &
          !$omp private(mp,alm)
          do inu=1,N_freq
             do ipix=1,N_pix
                J_j = input_map(:,:,ipix,inu)
                ! Equivalent to numpy np.einsum('ab,bc,cd,ad',Pi,J,Pj,J.conj())
                ! Also need to take real part and divide by 2
                output_map(ipix,inu,i,j) = &
                     dble(sum(matmul(matmul(Pi,J_j),Pj)*conjg(J_j)))/2
             enddo

             if (rotate_maps) then
                ! Also rotate map
                mp = output_map(:,inu,i,j)
                call map2alm(N_side,LMAX,LMAX,mp,alm)
                call rotate_alm(LMAX,alm,r_psi,r_theta,r_phi)
                call alm2map(N_side,LMAX,LMAX,alm,mp)
                output_map(:,inu,i,j) = mp
             endif
          enddo
          !$omp end parallel do
       enddo
    enddo


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called calc Mueller maps'
    return


  contains


    pure function stokes_matrix(i)
      ! Default
      implicit none

      ! Function arguments
      integer(4), intent(in)     :: i
      complex(8), dimension(2,2) :: stokes_matrix

      select case(i)
      case(1)
         ! Stokes I
         stokes_matrix(1,1) = (1D0, 0D0)
         stokes_matrix(2,1) = (0D0, 0D0)
         stokes_matrix(1,2) = (0D0, 0D0)
         stokes_matrix(2,2) = (1D0, 0D0)
      case(2)
         ! Stokes Q
         stokes_matrix(1,1) = ( 1D0, 0D0)
         stokes_matrix(2,1) = ( 0D0, 0D0)
         stokes_matrix(1,2) = ( 0D0, 0D0)
         stokes_matrix(2,2) = (-1D0, 0D0)
      case(3)
         ! Stokes U
         stokes_matrix(1,1) = (0D0, 0D0)
         stokes_matrix(2,1) = (1D0, 0D0)
         stokes_matrix(1,2) = (1D0, 0D0)
         stokes_matrix(2,2) = (0D0, 0D0)
      case(4)
         ! Stokes V
         stokes_matrix(1,1) = (0D0,  0D0)
         stokes_matrix(2,1) = (0D0,  1D0)
         stokes_matrix(1,2) = (0D0, -1D0)
         stokes_matrix(2,2) = (0D0,  0D0)
      case default
         ! Garbage
         stokes_matrix = -1
      end select

      return
    end function stokes_matrix


  end subroutine calc_mueller_maps


!------------------------------------------------------------------------------!


  subroutine calc_vis_maps(input_map, output_map)
    ! Default
    implicit none


    ! Subroutine arguments
    complex(8), dimension(:,:,:,:), intent(in)  :: input_map
    complex(8), dimension(:,:,:),   intent(out) :: output_map


    ! Local variables
    integer(4) :: i,j,ipix,inu


    ! Local arrays
    real(8),    dimension(N_pix)           :: mp1,mp2
    complex(8), dimension(2,2)             :: J_i,J_j
    complex(8), dimension(1,0:LMAX,0:LMAX) :: alm


    ! Timing variables
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
    tr1 = omp_get_wtime()


    ! Initialize output array
    output_map = 0

    ! Compute visibility map elements
    ! Assumes there is no Stokes component on the sky except I

    ! Parallelize over frequency
    !$omp parallel do               &
    !$omp default(shared)           &
    !$omp private(ipix,inu,mp1,mp2) &
    !$omp private(J_i,J_j,alm,i,j)
    do inu=1,N_freq
       do ipix=1,N_pix
          J_i = input_map(:,:,ipix,inu)

          ! Compute conjugate transpose
          J_j = transpose(conjg(J_i))

          ! Multiply matrices and unpack elements of visibility
          ! Want to compute J * J^\dagger
          J_i = matmul(J_i, J_j)
          output_map(ipix,inu,1) = J_i(1,1)
          output_map(ipix,inu,2) = J_i(2,1)
          output_map(ipix,inu,3) = J_i(1,2)
          output_map(ipix,inu,4) = J_i(2,2)
       enddo

       if (rotate_maps) then
          ! Also rotate map
          ! Rotate each component separately
          do i=1,4
             ! Rotate real and complex separately
             mp1 = dble(output_map(:,inu,i))
             mp2 = dimag(output_map(:,inu,i))

             ! Real part
             call map2alm(N_side,LMAX,LMAX,mp1,alm)
             call rotate_alm(LMAX,alm,r_psi,r_theta,r_phi)
             call alm2map(N_side,LMAX,LMAX,alm,mp1)

             ! Imaginary part
             call map2alm(N_side,LMAX,LMAX,mp2,alm)
             call rotate_alm(LMAX,alm,r_psi,r_theta,r_phi)
             call alm2map(N_side,LMAX,LMAX,alm,mp2)

             ! Pack it up together
             output_map(:,inu,i) = cmplx(mp1, mp2)
          enddo
       endif
    enddo
    !$omp end parallel do


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called calc vis maps'
    return
  end subroutine calc_vis_maps


!------------------------------------------------------------------------------!


  subroutine compute_xi_mueller(input_maps,xi)
    ! Default
    implicit none


    ! Subroutine arguments
    real(8),    dimension(:,:,:,:), intent(in)  :: input_maps
    complex(8), dimension(:,:,:,:), intent(out) :: xi


    ! Local variables
    integer(4) :: inu,il,im,iphi,itau,ibl,istokes
    real(8)    :: arg,jl,nu,theta,phi,tauh
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

    ! Initialize tauh array
    do itau=1,N_bl
       tauh_vals(itau) = r_bl(itau)/c_mks
    enddo

    ! Initialize phi_vals array
    do iphi=1,N_phi
       phi_vals(iphi) = phi_min + (iphi-1)*d_phi
    enddo

    ! Loop over Stokes parameters
    do istokes=1,4
       ! Progress report
       write(*,*) "Stokes ",istokes

       ! Loop over frequency to get alm's of Stokes map
       do inu=1,N_freq
          ! Get frequency in Hz
          nu  = nu_vals(inu)

          ! Get healpix map and convert to a_lm
          ! We assume leakage is dominant over intrinsic polarization signal,
          !   e.g., I -> Q' >> Q -> Q'
          mp = input_maps(:,inu,istokes,1)
          call map2alm(N_side,LMAX,LMAX,mp,alm)

          ! Loop over phi values
          do iphi=1,N_phi
             theta = b_theta
             phi   = phi_vals(iphi)

             ! Loop over baseline values
             do itau=1,N_bl
                ! Get tauh in s
                tauh = tauh_vals(itau)

                ! Compute argument of spherical Bessel function
                arg  = 2*pi*tauh*nu

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
                         ! a_{l, -m} = (-1)**m * a_{l, m}^*
                         ! Y_l^{-m}  = (-1)**m * (Y_l^m)^*
                         x = x + prefac*conjg(a)*jl*conjg(y)
                      endif
                   enddo
                enddo
                !$omp end parallel do

                ! Save value
                xi(inu,iphi,istokes,itau) = x
             enddo
          enddo
       enddo
    enddo

    ! Add overall normalization
    xi = xi*sqrt(4*pi)

    ! Convert frequencies back to MHz
    nu_vals = nu_vals/1D6


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called compute Xi Mueller'
    return
  end subroutine compute_xi_mueller


!------------------------------------------------------------------------------!


  subroutine compute_xi_vis(input_maps,xi)
    ! Default
    implicit none


    ! Subroutine arguments
    complex(8), dimension(:,:,:),   intent(in)  :: input_maps
    complex(8), dimension(:,:,:,:), intent(out) :: xi


    ! Local variables
    integer(4) :: inu,il,im,iphi,itau,ibl,istokes
    real(8)    :: arg,jl,nu,theta,phi,tauh
    complex(8) :: prefac,y,a,x


    ! Local arrays
    complex(8), dimension(1,0:LMAX,0:LMAX) :: alm1,alm2
    real(8),    dimension(N_pix)           :: mp1,mp2


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

    ! Initialize tauh array
    do itau=1,N_bl
       tauh_vals(itau) = r_bl(itau)/c_mks
    enddo

    ! Initialize phi_vals array
    do iphi=1,N_phi
       phi_vals(iphi) = phi_min + (iphi-1)*d_phi
    enddo

    ! Loop over visibility parameters
    do istokes=1,4
       write(*,*) "Stokes ",istokes

       ! Loop over frequency to get alm's of Stokes map
       do inu=1,N_freq
          ! Get frequency in Hz
          nu  = nu_vals(inu)

          ! Get healpix map and convert to a_lm
          ! We assume leakage is dominant over intrinsic polarization signal,
          !   e.g., I -> Q' >> Q -> Q'
          mp1 = dble(input_maps(:,inu,istokes))
          call map2alm(N_side,LMAX,LMAX,mp1,alm1)

          ! Off-diagonal terms are in general complex
          if (istokes==2 .or. istokes==3) then
             mp2 = dimag(input_maps(:,inu,istokes))
             call map2alm(N_side,LMAX,LMAX,mp2,alm2)
          endif

          ! Loop over phi values
          do iphi=1,N_phi
             theta = b_theta
             phi   = phi_vals(iphi)

             ! Loop over baseline values
             do itau=1,N_bl
                ! Get tauh in s
                tauh = tauh_vals(itau)

                ! Compute argument of spherical Bessel function
                arg  = 2*pi*tauh*nu

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
                      a      = alm1(1,il,im)

                      x = x + prefac*a*jl*y
                      if (im > 0) then
                         ! for negative m, we have
                         ! a_{l, -m} = (-1)**m * a_{l, m}^*
                         ! Y_l^{-m}  = (-1)**m * (Y_l^m)^*
                         x = x + prefac*conjg(a)*jl*conjg(y)
                      endif

                      ! Handle off-diagonal elements too
                      if (istokes==2 .or. istokes==3) then
                         a = (0D0,1D0)*alm2(1,il,im)
                         x = x + prefac*a*jl*y
                         if (im > 0) then
                            x = x + prefac*conjg(a)*jl*conjg(y)
                         endif
                      endif
                   enddo
                enddo
                !$omp end parallel do

                ! Save value
                xi(inu,iphi,istokes,itau) = x
             enddo
          enddo
       enddo
    enddo

    ! Add overall normalization
    xi = xi*sqrt(4*pi)

    ! Convert frequencies back to MHz
    nu_vals = nu_vals/1D6

    ! Convert baseline separation to ns
    tauh_vals = tauh_vals*1D9


    tr2 = omp_get_wtime()
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called compute Xi vis'
    return
  end subroutine compute_xi_vis


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
