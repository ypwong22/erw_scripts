! gfortran main.F90 -I /sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/openblas-0.3.23-ejxcjjxiy43ruq6kjat43xptaqscotvn/include/ -L /sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/openblas-0.3.23-ejxcjjxiy43ruq6kjat43xptaqscotvn/lib/ -lopenblas

module global_parameters
    implicit none
    integer, parameter :: mixing_layer = 6
    ! change to uniform 0.1m grid, 10 layers
    real(8), parameter :: zsoi(7) = (/0.0071, 0.0279, 0.0623, 0.1189, 0.2122, 0.3661, 0.6198/)
    real(8), parameter :: zisoi(7) = (/0.0175, 0.0451, 0.0906, 0.1655, 0.2891, 0.4929, 0.8289/)
    real(8), parameter :: dzsoi_decomp(7) = (/0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360/)
end module global_parameters

program main
  use global_parameters
  implicit none

  real(8) :: conc_trcr(1:mixing_layer)
  real(8) :: adv_flux(1:mixing_layer+1)
  real(8) :: diffus(1:mixing_layer)
  real(8) :: source(1:mixing_layer)
  real(8) :: surf_bc
  real(8) :: dtime
  real(8) :: vwc(1:mixing_layer)
  real(8) :: conc_change_rate(1:mixing_layer)
  integer :: i, j, unit_num

  vwc(1:mixing_layer) = 0.25

  conc_trcr(1:mixing_layer) = 1e-5 * vwc(1) * 1e3 ! mol/L water => mol/m3 soil
  adv_flux(1:mixing_layer+1) = -1e-10 ! m/s
  diffus(1:mixing_layer) = 7.93e-10
  source(1:mixing_layer) = 0.
  dtime = 3600.
  surf_bc = 1e-3 * 1e3 ! mol/L water => mol/m3 soil

  open(newunit=unit_num, file='output.csv', status='replace', action='write')

  write(unit_num, '(I0,",",1pe13.4,",",1pe13.4,",",1pe13.4,",",1pe13.4,",",1pe13.4,",",1pe13.4)') 0, conc_trcr

  do i = 1,2678400
    call advection_diffusion(conc_trcr(1:mixing_layer),adv_flux(1:mixing_layer+1), &
      diffus(1:mixing_layer),source(1:mixing_layer),surf_bc,dtime, & 
      vwc(1:mixing_layer),conc_change_rate(1:mixing_layer))

    do j = 1,mixing_layer
      conc_trcr(j) = conc_change_rate(j) * dtime + conc_trcr(j)
    end do

    if (mod(i, 3600) == 1) then
      write(unit_num, '(I0,",",1pe13.4,",",1pe13.4,",",1pe13.4,",",1pe13.4,",",1pe13.4,",",1pe13.4)') i, conc_trcr
    end if
  end do

  close(unit_num)

contains

  subroutine advection_diffusion(conc_trcr,adv_flux,diffus,source,surf_bc,dtime,vwc,conc_change_rate)
    ! From B. Sulman; edited layer depth; soil bulk concentration can use g/m3
    ! 
    ! Advection and diffusion for a single tracer in one column given diffusion coefficient, flow, and source-sink terms
    ! Based on SoilLittVertTranspMod, which implements S. V. Patankar, Numerical Heat Transfer and luid Flow, Series in Computational Methods in Mechanics and Thermal Sciences, Hemisphere Publishing Corp., 1980. Chapter 5
    ! Not sure if this belongs here or somewhere else. Is it bad to do this in the EMI subroutine?

    real(8), intent(in) :: conc_trcr(1:mixing_layer)  ! Bulk concentration (e.g. mol/m3)
    real(8), intent(in) :: adv_flux(1:mixing_layer+1) ! (m/s), vertical into layer (down is negative)
    real(8), intent(in) :: diffus(1:mixing_layer)  ! diffusivity (m2/s)
    real(8), intent(in) :: source(1:mixing_layer)  ! Source term (mol/m3/s)

    real(8), intent(in) :: surf_bc                 ! Surface boundary layer concentration (for infiltration) (e.g. mol/m3 [water])
    real(8), intent(in) :: dtime                   ! Time step (s)
    real(8), intent(in) :: vwc(1:mixing_layer)     ! Volumetric soil moisture in layer (m3/m3)
    real(8), intent(out):: conc_change_rate(1:mixing_layer) ! Bulk concentration (e.g. mol/m3/s)

    ! Local variables
    real(8) :: aaa                          ! "A" function in Patankar
    real(8) :: pe                           ! Pe for "A" function in Patankar
    real(8) :: w_m1, w_p1                   ! Weights for calculating harmonic mean of diffusivity
    real(8) :: d_m1, d_p1                   ! Harmonic mean of diffusivity
    real(8) :: vwc_m1, vwc_p1               ! Harmonic mean of soil moisture
    real(8) :: a_tri(0:mixing_layer+1)      ! "a" vector for tridiagonal matrix
    real(8) :: b_tri(0:mixing_layer+1)      ! "b" vector for tridiagonal matrix
    real(8) :: c_tri(0:mixing_layer+1)      ! "c" vector for tridiagonal matrix
    real(8) :: r_tri(0:mixing_layer+1)      ! "r" vector for tridiagonal solution
    real(8) :: d_p1_zp1(1:mixing_layer+1)   ! diffusivity/delta_z for next j  (set to zero for no diffusion)
    real(8) :: d_m1_zm1(1:mixing_layer+1)   ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(8) :: f_p1(1:mixing_layer+1)       ! water flux for next j
    real(8) :: f_m1(1:mixing_layer+1)       ! water flux for previous j
    real(8) :: pe_p1(1:mixing_layer+1)      ! Peclet # for next j
    real(8) :: pe_m1(1:mixing_layer+1)      ! Peclet # for previous j
    real(8) :: dz_node(1:mixing_layer+1)    ! difference between nodes
    real(8) :: a_p_0
    real(8) :: conc_after(0:mixing_layer+1)

    integer :: j, info

    ! Statement function
    aaa (pe) = max (0., (1. - 0.1 * abs(pe))**5)  ! "A" function from Patankar, Table 5.2, pg 95

    ! Set the distance between the node and the one ABOVE it   
    dz_node(1) = zsoi(1)
    do j = 2,mixing_layer+1
      dz_node(j)= zsoi(j) - zsoi(j-1)
    enddo

    !print *, 'adv_flux',adv_flux(1:mixing_layer+1)
    !print *, 'diffus',diffus(1:mixing_layer)
    !print *, 'source',source(1:mixing_layer)

    ! Calculate the D and F terms in the Patankar algorithm
    ! d: diffusivity
    ! f: flow
    ! m: layer above ("E" East in Pantakar)
    ! p: layer below ("W" West in Pantakar; positive flow is W->E)
    ! pe: Peclet number (ratio of convection to diffusion)
    do j = 1,mixing_layer
      if (j == 1) then
        d_m1_zm1(j) = 0. ! no diffusion between j = 1 and atmosphere
        w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
        if ( diffus(j+1) > 0. .and. diffus(j) > 0.) then
          d_p1 = 1. / ((1. - w_p1) / diffus(j) + w_p1 / diffus(j+1)) ! Harmonic mean of diffus
        else
          d_p1 = 0.
        endif
        d_p1_zp1(j) = d_p1 / dz_node(j+1) ! Eq. 5.9, D = \Tau / \delta x
        vwc_m1 = vwc(j)
        vwc_p1 = 1. / ((1. - w_p1) / vwc(j) + w_p1 / vwc(j+1))
        f_m1(j) = adv_flux(j) / vwc_m1 ! Include infiltration here
        f_p1(j) = adv_flux(j+1) / vwc_p1
        pe_m1(j) = 0.
        pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
      elseif (j == mixing_layer) then
          ! At the bottom, assume no gradient in d_z (i.e., they're the same)
          w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
          if ( diffus(j) > 0. .and. diffus(j-1) > 0.) then
            d_m1 = 1. / ((1. - w_m1) / diffus(j) + w_m1 / diffus(j-1)) ! Harmonic mean of diffus
          else
            d_m1 = 0.
          endif
          d_m1_zm1(j) = d_m1 / dz_node(j)
          d_p1_zp1(j) = d_m1_zm1(j) ! Set to be the same
          vwc_m1 = 1. / ((1. - w_m1) / vwc(j-1) + w_m1 / vwc(j))
          f_m1(j) = adv_flux(j) / vwc_m1
          ! f_p1(j) = adv_flux(j+1) / vwc(j)
          ! f_p1(j) = 0.
          pe_m1(j) = f_m1(j) / d_m1_zm1(j) ! Peclet #
          pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
      else
          ! Use distance from j-1 node to interface with j divided by distance between nodes
          w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
          if ( diffus(j-1) > 0. .and. diffus(j) > 0.) then
            d_m1 = 1. / ((1. - w_m1) / diffus(j) + w_m1 / diffus(j-1)) ! Harmonic mean of diffus
          else
            d_m1 = 0.
          endif
          w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
          if ( diffus(j+1) > 0. .and. diffus(j) > 0.) then
            d_p1 = 1. / ((1. - w_p1) / diffus(j) + w_p1 / diffus(j+1)) ! Harmonic mean of diffus
          else
            d_p1 = (1. - w_p1) * diffus(j) + w_p1 * diffus(j+1) ! Arithmetic mean of diffus
          endif
          d_m1_zm1(j) = d_m1 / dz_node(j)
          d_p1_zp1(j) = d_p1 / dz_node(j+1)
          vwc_m1 = 1. / ((1. - w_m1) / vwc(j-1) + w_m1 / vwc(j))
          vwc_p1 = 1. / ((1. - w_p1) / vwc(j) + w_p1 / vwc(j+1))
          f_m1(j) = adv_flux(j) / vwc_m1
          f_p1(j) = adv_flux(j+1) / vwc_p1
          pe_m1(j) = f_m1(j) / d_m1_zm1(j) ! Peclet #
          pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
      end if
    enddo ! j; mixing_layer


    ! Calculate the tridiagonal coefficients
    ! Coefficients of tridiagonal problem: a_i*x_(i-1) + b_i*(x_i) + c_i*x_(i+1) = r_i
    ! Here, this is equivalent to Patankar equation 5.56 and 5.57 (but in one dimension):
    ! a_P*phi_P = a_E*phi_E + a_W*phi_W + b [phi is concentration, = x in tridiagonal]. Converting East/West to above/below
    ! -> -a_E*phi_E + a_P*phi_P - a_W+phi_W = b [= sourcesink - d rho * \phi/dt]
    ! -a_tri = a_above = D_above*A(Pe)+max(-F_above,0); D_above=diffus_above/dz
    ! b_tri = a_above+a_below+vwc*dz/dt
    ! -c_tri = D_below*A(Pe)+max(F_below,0); D_below = diffus_below/dz
    ! r_tri = b = source_const*dz + conc*vwc*dz/dt
    do j = 0,mixing_layer +1

      if (j > 0 .and. j < mixing_layer+1) then
          a_p_0 =  dzsoi_decomp(j) / dtime / vwc(j) ! Should this be multiplied by layer water content (for vwc)?
      endif

      if (j == 0) then ! top layer (atmosphere)
          a_tri(j) = 0.
          b_tri(j) = 1.
          c_tri(j) = -1.
          r_tri(j) = 0.
      elseif (j == 1) then
          a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0.)) ! Eqn 5.47 Patankar
          c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0.))
          b_tri(j) = -a_tri(j) - c_tri(j) + a_p_0
          ! r_tri includes infiltration assuming same concentration as top layer. May want to change to either provide upper boundary condition or include in source term
          ! r_tri(j) = source(j) * dzsoi_decomp(j) + (a_p_0 - adv_flux(j)) * conc_trcr(j)
          r_tri(j) = source(j) * dzsoi_decomp(j) + a_p_0 * conc_trcr(j)
          if(adv_flux(j)<0) then ! downward flow (infiltration)
            r_tri(j) = r_tri(j) - adv_flux(j)*surf_bc
            !  write (iulog,*) __LINE__,adv_flux(j),surf_bc,adv_flux(j)*surf_bc
          else ! upward flow to the surface
            r_tri(j) = r_tri(j) - adv_flux(j)*conc_trcr(j)
            ! write (iulog,*) __LINE__,adv_flux(j),conc_trcr(j),adv_flux(j)*conc_trcr(j)
          endif
          
      elseif (j < mixing_layer+1) then
          a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0.)) ! Eqn 5.47 Patankar
          c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0.))
          b_tri(j) = -a_tri(j) - c_tri(j) + a_p_0
          r_tri(j) = source(j) * dzsoi_decomp(j) + a_p_0 * conc_trcr(j) ! Eq. 5.57

      else ! j==mixing_layer+1; 0 concentration gradient at bottom
          a_tri(j) = -1.
          b_tri(j) = 1.
          c_tri(j) = 0. 
          r_tri(j) = 0.
      endif
    enddo ! j; mixing_layer

    ! Solve for the concentration profile for this time step
    ! call Tridiagonal(0, mixing_layer+1, 0, a_tri, b_tri, c_tri, r_tri, conc_after)
    ! This is the LAPACK tridiagonal solver which gave more accurate results in my testing
    call dgtsv( mixing_layer+2, 1, c_tri(0:mixing_layer), b_tri, a_tri(1:mixing_layer+1),  & 
                r_tri, mixing_layer+2, info )

    ! if(info < 0) call endrun(msg='dgtsv error in adv_diff line __LINE__: illegal argument')
    ! if(info > 0) call endrun(msg='dgtsv error in adv_diff line __LINE__: singular matrix')
    do j = 1, mixing_layer
      conc_after(j) = r_tri(j) ! dgtsv solver saves output into the r_tri matrix
      conc_change_rate(j) = (conc_after(j)-conc_trcr(j))/dtime
    end do

  end subroutine advection_diffusion

end program main