! gfortran main.F90 -I /sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/openblas-0.3.23-ejxcjjxiy43ruq6kjat43xptaqscotvn/include/ -L /sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/openblas-0.3.23-ejxcjjxiy43ruq6kjat43xptaqscotvn/lib/ -lopenblas

module global_parameters
    use, intrinsic :: iso_fortran_env, only : r8 => real64

    implicit none

    ! change to uniform 0.1m grid, 10 layers
    integer, parameter :: nlevbed = 10
    integer, parameter :: nlevsoi = 10
    real(r8) :: zsoi(1:nlevsoi) = (/0.007100635_r8, 0.027925_r8, 0.06225858_r8, 0.1188651_r8, &
        0.2121934_r8, 0.3660658_r8, 0.6197585_r8, 1.038027_r8, 1.727635_r8, 2.864607_r8/)
    real(r8) :: zisoi(1:nlevsoi+1) = (/0._r8, 0.0175_r8, 0.0451_r8, 0.0906_r8, 0.1655_r8, &
        0.2891_r8, 0.4929_r8, 0.8289_r8, 1.3828_r8, 2.2961_r8, 3.8019_r8/)

end module global_parameters

program main
  use global_parameters
  use, intrinsic :: iso_fortran_env, only : r8 => real64

  implicit none

  real(r8) :: c_prev(1:nlevsoi)
  real(r8) :: rain_conc
  real(r8) :: q_int(1:nlevsoi+1)
  real(r8) :: theta(1:nlevsoi)
  real(r8) :: watsat(1:nlevsoi)
  real(r8) :: srcsk(1:nlevsoi)
  real(r8) :: d0
  real(r8) :: dt
  real(r8) :: dz(1:nlevsoi)
  real(r8) :: c_next(1:nlevsoi)
  integer  :: i,j,k,unit_num

  c_prev(1:nlevsoi) = 0._r8
  rain_conc = 0._r8
  q_int(1:nlevsoi+1) = 1e-10_r8 ! m/s
  theta(1:nlevsoi) = 0.4_r8
  watsat(1:nlevsoi) = 0.6_r8
  srcsk(1:nlevsoi) = 0._r8
  srcsk(1) = 1e-5_r8
  d0 = 7.93e-10_r8
  dt = 1._r8

  do j = 1, nlevsoi
      dz(j) = zisoi(j+1) - zisoi(j)
  end do

  open(newunit=unit_num, file='output.csv', status='replace', action='write')

  write(unit_num,'(I0,10(",",1PE13.4))')  0, (c_prev(k), k = 1, nlevsoi)

  do i = 1,100
    call advection_diffusion(c_prev, rain_conc, q_int, theta, watsat, srcsk, d0, dt, &
                             dz, nlevbed, c_next)
    write(unit_num,'(I0,10(",",1PE13.4))')  i, (c_next(k), k = 1, nlevsoi)

    c_prev(1:nlevsoi) = c_next(1:nlevsoi)
  end do

  close(unit_num)

contains

  pure function analytical_c(C0, r, k, t1, t2) result(c)
    !------------------------------------------------------------------
    ! Concentration at t2, starting from C0 at t1,  dC/dt = k C + r
    !------------------------------------------------------------------
    real(r8), intent(in) :: C0, r, k, t1, t2
    real(r8)             :: c
    c = (C0 + r/k) * exp(k*(t2 - t1)) - r/k
  end function analytical_c


  pure function analytical_c_int(C0, r, k, t1, t2) result(c_int)
    !------------------------------------------------------------------
    ! Integral of C over [t1, t2] for the same ODE
    !------------------------------------------------------------------
    real(r8), intent(in) :: C0, r, k, t1, t2
    real(r8)             :: c_int
    c_int = (C0 + r/k)/k * (exp(k*(t2 - t1)) - 1.0_r8) - (r/k) * (t2 - t1)
  end function analytical_c_int


  pure function analytical_dt(c1, c2, r, k) result(dt)
    !------------------------------------------------------------------
    ! Time Δt needed for C to evolve from c1 to c2
    !------------------------------------------------------------------
    real(r8), intent(in) :: c1, c2, r, k
    real(r8)             :: dt
    dt = log((c2 + r/k) / (c1 + r/k)) / k
  end function analytical_dt


  subroutine advection_diffusion(c_prev, rain_conc, q_int, theta, watsat, srcsk, &
                                 d0, dt, dz, nlevbed, c_next)
    !------------------------------------------------------------------
    ! Use the explicit solution at each time step
    ! 
    ! If we ignore spatial dependency, the equation becomes very simple
    ! 
    ! Consider a generic cell C_i with adjacent cells
    !     C_{i-1} <--- Δ x_i ---> C_i, D_{eff,i} <--- Δ x_{i+1} ---> C_{i+1}
    ! 
    ! Total mass change within C_i due to outflow
    !     Δ x  * θ * \frac{ dC_i }{ dt }
    ! 
    ! Diffusion to the upper cell, if C_i > C_{i-1}
    !     dup = D_{eff,i} * I[C_i > C_{i-1}] * \frac{C_i - C_{i-1}}{Δ x_i}
    ! 
    ! Diffusion to the lower cell, if C_i > C_{i+1}
    !     dlow = D_{eff,i} * I[C_i > C_{i+1}] * \frac{C_i - C_{i+1}}{Δ x_{i+1}}
    ! 
    ! Advection to upper cell, if q_{in,i} < 0
    !     aup = I[q_{in,i} < 0] * abs( q_{in,i} ) * C_i
    ! 
    ! Advection to lower cell, if q_{out,i} > 0
    !     alow = I[q_{out,i} > 0] * q_{out,i} * C_i
    ! 
    ! Rain inflow
    !     rain = I[q_{in,i} > 0] * q_{in,0} * C_{rain}
    ! 
    ! The relationship between self-change and the four outflow and internal source/sink, 
    !     after simplification, is
    ! 
    !     Δ x * θ * \frac{ dC_i }{ dt } = 
    !       - ( \frac{ D_{eff,i} * I[C_i > C_{i-1}] }{Δ x_i} + 
    !           \frac{ D_{eff,i} * I[C_i > C_{i+1}] }{Δ x_{i+1}} + 
    !           I[q_{in,i} < 0] * abs( q_{in,i} ) + I[q_{out,i} > 0] * q_{out,i} ) * C_i
    !       + ( D_{eff,i} * I[C_i > C_{i-1}] \frac{C_{i-1}}{Δ x_i} )
    !       + ( D_{eff,i} * I[C_i > C_{i+1}] \frac{C_{i+1}}{Δ x_{i+1}} )
    !       + I[q_{in,i} > 0] * q_{in,0} * C_{rain}
    !       + R
    ! 
    !     , which is of the form dC/dt = kC + r, and the analytical solution is just
    !          C = (C0 + r/k)e^{kt} - r/k
    !          on the interval [t, t + Δ t]
    ! 
    ! The equations means when t → ∞, C → -r/k
    ! 
    ! It is guaranteed that k < 0. Therefore, 
    ! 
    ! if r > 0, the concentration will not be negative for large t (we can be sure )
    ! 
    ! If r < 0, the concentration may be negative because the equilibrium value is
    !     negative. This can only happen when r < 0. In this case, the analytical 
    !     form may need to be assessed piecewise to insure non-negativity:
    ! 
    !     Step 1. We evaluate the analytical form at (t + Δ t),
    !             if the resulting C_t' > C_{i-1} and C_t' > C_{i+1}, then C_t = C_t', 
    !             otherwise, we need to go to step 2
    !     Step 2. We assess piecewise: 
    !             first, find the dt' such that C(t + dt') = max(C_{i-1}, C_{i+1})
    !             then, on the interval [t+dt', t + Δ t], set the appropriate I[⋅]
    !                    term to 0, re-evaluate to t + Δ t
    !             third, if the result gives C(t + Δ t) > min(C_{i-1}, C_{i}), stop
    !                    otherwise, need a third piecewise, go to step 3
    !     Step 3. find the dt" such that C(t + dt") = min(C_{i-1}, C_{i+1}), on the
    !             interval [t+dt", t + Δ t], all the diffusion terms need to be 0.
    !             Re-integrate the pure-advection term to t + Δ t
    !------------------------------------------------------------------

    ! Input-output variables
    real(r8), intent(in) :: c_prev(1:nlevsoi) ! start concentration (mol/m3-water)
    real(r8), intent(in) :: rain_conc ! upper boundary condition (rain chemistry) (mol/m3-water)
    real(r8), intent(in) :: q_int(1:nlevsoi+1) ! water flux at grid boundaries, positive downwards (m/s)
    real(r8), intent(in) :: theta(1:nlevsoi) ! soil moisture (m3/m3)
    real(r8), intent(in) :: watsat(1:nlevsoi) ! soil porosity (m3/m3)
    real(r8), intent(in) :: srcsk(1:nlevsoi) ! source/sink strength (mol/m3-soil/s)
    real(r8), intent(in) :: d0 ! diffusion coefficient in water (m2/s)
    real(r8), intent(in) :: dt ! time step size (s)
    real(r8), intent(in) :: dz(1:nlevsoi) ! soil layer thickness (m)
    integer , intent(in) :: nlevbed ! number of hydrologically active layers
    real(r8), intent(out) :: c_next(1:nlevsoi) ! updated concentration (mol/m3-water)

    ! Local variables
    real(r8) :: Deff(1:nlevsoi) ! effective diffusion coefficient in porous media (m2/s)
    real(r8) :: dx(1:nlevsoi+2) ! padded dz; the first & last elements don't matter so long as >0
    real(r8) :: scaler(1:nlevsoi) ! \Delta x * 0 on the LHS
    real(r8) :: sourcesink(1:nlevsoi) ! source sink, converted from mol/m3-soil/s to mol/m2/s; with upper boundary
    integer  :: i, i1, i2, i3, i4 ! indices and boolean indicators
    integer  :: i1_keep, i2_keep  ! side of diffusion
    real(r8) :: k, r ! temporary variables to save the ODE parameters
    real(r8) :: c, c_p, cmax, cmin, c_int, c_int_p, c_int_pp, dt_p, dt_pp ! temporary variables for ODE integration

    ! Local variables, can be converted to diagnostic outputs if needed
    real(r8) :: niter(1:nlevsoi) ! keep track how many piecewise integrations are done
    real(r8) :: dc_up(1:nlevsoi) ! flow upward due to diffusion and advection (mol/m2-ground/s)
    real(r8) :: dc_down(1:nlevsoi) ! flow downward due to diffusion and advection (mol/m2-ground/s)


    Deff = d0 * theta**(10.0_r8/3.0_r8) / watsat**2
    !print *, Deff

    dx(2:nlevsoi+1) = dz(1:nlevsoi)
    dx(1) = 100._r8
    dx(nlevsoi+2) = 100._r8

    scaler = dz * theta
    !print *, scaler

    ! convert sourcesink from mol/m3-soil/s to mol/m2/s
    sourcesink = srcsk * dz
    if (q_int(1) > 0._r8) then
      sourcesink(1) = sourcesink(1) + q_int(1) * rain_conc
    end if

    do i = 1,nlevbed
      if (i == 1) then
        ! top layer: no diffusion to above
        i1 = 0
        i2 = merge(1, 0, c_prev(i) > c_prev(i+1)) ! boolean to integer

      elseif (i < nlevsoi) then
        ! middle layers
        i1 = merge(1, 0, c_prev(i) > c_prev(i-1))
        i2 = merge(1, 0, c_prev(i) > c_prev(i+1))

      else
        ! bottom layer: no diffusion to below
        i1 = merge(1, 0, c_prev(i) > c_prev(i-1))
        i2 = 0
      end if

      i3 = merge(1, 0, q_int(i) < 0)
      i4 = merge(1, 0, q_int(i+1) > 0)

      k = - (i1*Deff(i)/dx(i) + i2*Deff(i)/dx(i+1) + i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)

      if (i == 1) then
        r = (Deff(i)*i2*c_prev(i+1)/dx(i+1) + sourcesink(i)) / scaler(i)
      elseif (i < nlevsoi) then
        r = (Deff(i)*i1*c_prev(i-1)/dx(i) + Deff(i)*i2*c_prev(i+1)/dx(i+1) \
                + sourcesink(i)) / scaler(i)
      else
        r = (Deff(i)*i1*c_prev(i-1)/dx(i) + sourcesink(i)) / scaler(i)
      end if

      !print *, i, r

      ! if there is no flow out of this cell, degrades to linear source
      if (k == 0._r8) then
        c_next(i) = c_prev(i) + r*dt
        dc_up(i) = 0
        dc_down(i) = 0

      else
        ! otherwise, do actual calculations

        ! step 1: find end-of-step solution
        c = analytical_c(c_prev(i), r, k, 0._r8, dt)

        ! top layer
        if (i == 1) then
          cmax = c_prev(i+1)

          ! end-of-step solution works, or there is no diffusion to begin with
          if ((c > cmax) .or. (i2 == 0)) then
              c_next(i) = c
              niter(i) = 1

              ! integration of c over [0, dt]
              c_int = analytical_c_int(c_prev(i), r, k, 0._r8, dt)

              ! use c_int to calculate the fluxes
              dc_up(i) = i3*abs(q_int(i))*c_int
              dc_down(i) = i2*Deff(i)/dx(i+1) * (c_int - c_prev(i+1)*dt) + i4*q_int(i+1)*c_int

          else
              ! step 2: piecewise integrate to cmax, then update k, r and continue
              dt_p = analytical_dt(c_prev(i), cmax, r, k)

              ! integration of c over [0, dt']
              c_int_p = analytical_c_int(c_prev(i), r, k, 0._r8, dt_p)

              ! use c_int to calculate the fluxes during [0, dt']
              dc_up(i) = i3*abs(q_int(i))*c_int_p
              dc_down(i) = i2*Deff(i)/dx(i+1) * (c_int_p - c_prev(i+1)*dt_p) &
                            + i4*q_int(i+1)*c_int_p

              ! update k & r (no more diffusion)
              k = - (i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)
              r = sourcesink(i) / scaler(i)

              ! itegration of c over [dt', dt]
              c_int_pp = analytical_c_int(cmax, r, k, dt_p, dt)

              ! add to the fluxes
              dc_up(i) = dc_up(i) + i3*abs(q_int(i))*c_int_pp
              dc_down(i) = dc_down(i) + i4*q_int(i+1)*c_int_pp

              ! find final value
              c_next(i) = analytical_c(cmax, r, k, dt_p, dt)
              niter(i) = 2

          end if

        ! intermediate layer
        else if (i < nlevsoi) then

          if ((i1 == 0) .and. (i2 == 0)) then
            ! if there is no diffusion to begin with, end of solution works

            ! final value found
            c_next(i) = c
            niter(i) = 1

            ! integration of c over (0, dt)
            c_int = analytical_c_int(c_prev(i), r, k, 0._r8, dt)

            ! calculate the fluxes during (0, dt) (all diffusion)
            dc_up(i) = i3*abs(q_int(i))*c_int
            dc_down(i) = i4*q_int(i+1)*c_int

            !print *, i, c_next(i), dc_up(i), dc_down(i)

          else
            ! if there is only one side diffusion
            if (i1 == 0) then
                cmax = c_prev(i+1) ! diffuse down
                cmin = 0._r8
                i1_keep = 0
                i2_keep = 1

            else if (i2 == 0) then
                cmax = c_prev(i-1) ! diffuse up
                cmin = 0._r8
                i1_keep = 1
                i2_keep = 0

            else
                ! both sides have diffusion, decide which side to begin with

                if (c_prev(i-1) > c_prev(i+1)) then
                    cmax = c_prev(i-1)
                    cmin = c_prev(i+1)
                    i1_keep = 0
                    i2_keep = 1
                else
                    cmax = c_prev(i+1)
                    cmin = c_prev(i-1)
                    i1_keep = 1
                    i2_keep = 0
                end if

            end if

            ! if end of solution works
            if (c > cmax) then

                ! final value found
                c_next(i) = c
                niter(i) = 1

                ! integration of c over [0, dt]
                c_int = analytical_c_int(c_prev(i), r, k, 0._r8, dt)

                ! calculate the fluxes during [0, dt] (all diffusion)
                dc_up(i) = i1*Deff(i)/dx(i) * (c_int - c_prev(i-1)*dt) + i3*abs(q_int(i))*c_int
                dc_down(i) = i2*Deff(i)/dx(i+1) * (c_int - c_prev(i+1)*dt) + i4*q_int(i+1)*c_int

                ! print *, i, c_next(i), dc_up(i), dc_down(i)
                ! print *, i2, Deff(i), dx(i+1), c_int, c_prev(i+1), dt, i4, q_int(i+1), c_int

            else
                ! step 2: piecewise integrate to cmax, then update k, r and continue
                dt_p = analytical_dt(c_prev(i), cmax, r, k)

                ! integration of c over (0, dt')
                c_int_p = analytical_c_int(c_prev(i), r, k, 0._r8, dt_p)

                ! calculate the fluxes during (0, dt') (all diffusion)
                dc_up(i) = i1*Deff(i)/dx(i) * (c_int_p - c_prev(i-1)*dt_p) &
                            + i3*abs(q_int(i))*c_int_p
                dc_down(i) = i2*Deff(i)/dx(i+1) * (c_int_p - c_prev(i+1)*dt_p) &
                              + i4*q_int(i+1)*c_int_p

                ! update k & r (drop one side of diffusion)
                k = - (i1_keep*i1*Deff(i)/dx(i) + i2_keep*i2*Deff(i)/dx(i+1) + &
                       i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)
                r = (Deff(i)*i1_keep*i1*c_prev(i-1)/dx(i) + &
                     Deff(i)*i2_keep*i2*c_prev(i+1)/dx(i+1) + sourcesink(i)) / scaler(i)

                ! continue to end, but we need another check
                c_p = analytical_c(cmax, r, k, dt_p, dt)

                if (c_p > cmin) then

                    ! itegration of c over [dt', dt]
                    c_int_p = analytical_c_int(cmax, r, k, dt_p, dt)

                    ! add the fluxes during (dt', dt) (drop one side of diffusion)
                    dc_up(i) = dc_up(i) &
                        + i1_keep*i1*Deff(i)/dx(i) * (c_int_p - cmax*(dt-dt_p)) &
                        + i3*abs(q_int(i))*c_int_p
                    dc_down(i) = dc_down(i) &
                        + i2_keep*i2*Deff(i)/dx(i+1) * (c_int_p - cmax*(dt-dt_p)) &
                        + i4*q_int(i+1)*c_int_p

                    ! final value at dt
                    c_next(i) = analytical_c(cmax, r, k, dt_p, dt)
                    niter(i) = 2

                else
                    ! step 3: do another piecewise integration

                    ! second stop point
                    dt_pp = analytical_dt(cmax, cmin, r, k)

                    ! itegration of c over [dt', dt"]
                    c_int_pp = analytical_c_int(c_p, r, k, dt_p, dt_pp)

                    ! add the fluxes during [dt', dt"] (drop one side of diffusion)
                    dc_up(i) = dc_up(i) &
                        + i1_keep*i1*Deff(i)/dx(i) * (c_int_pp - cmin*(dt_pp-dt_p)) &
                        + i3*abs(q_int(i))*c_int_pp
                    dc_down(i) = dc_down(i) &
                        + i2_keep*i2*Deff(i)/dx(i+1) * (c_int_pp - cmin*(dt_pp-dt_p)) &
                        + i4*q_int(i+1)*c_int_pp

                    ! update the k & r (drop all diffusion)
                    k = - (i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)
                    r = sourcesink(i) / scaler(i)

                    ! itegration of c over [dt", dt]
                    c_int_pp = analytical_c_int(cmin, r, k, dt_pp, dt)

                    ! add the fluxes during [dt", dt]
                    dc_up(i) = dc_up(i) + i3*abs(q_int(i))*c_int_pp
                    dc_down(i) = dc_down(i) + i4*q_int(i+1)*c_int_pp

                    ! final value
                    c_next(i) = analytical_c(cmin, r, k, dt_pp, dt)
                    niter(i) = 3
              end if
            end if
          end if

          ! print *, i, c_next(i), niter(i)

        ! last soil layer
        else
          cmax = c_prev(i-1)

          ! end-of-step solution works, or there is no diffusion to begin with 
          if ((c > cmax) .or. (i1 == 0)) then
              c_next(i) = c
              niter(i) = 1

              ! integration of c over (0, dt)
              c_int = analytical_c_int(c_prev(i), r, k, 0._r8, dt)

              ! use c_int to calculate the fluxes
              dc_up(i) = i1*Deff(i)/dx(i) * (c_int - c_prev(i-1)*dt) + i3*abs(q_int(i))*c_int
              dc_down(i) = i4*q_int(i+1)*c_int

          else
              ! step 2: piecewise integrate to cmax, then update k, r and continue
              dt_p = analytical_dt(c_prev(i-1), cmax, r, k)

              ! integration of c over [0, dt']
              c_int_p = analytical_c_int(c_prev(i), r, k, 0._r8, dt_p)

              ! calculate the fluxes during [0, dt'] (all diffusion)
              dc_up(i) = i1*Deff(i)/dx(i) * (c_int_p - c_prev(i-1)*dt_p) + i3*abs(q_int(i))*c_int_p
              dc_down(i) = i4*q_int(i+1)*c_int_p

              ! update k & r (no more diffusion)
              k = - (i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)
              r = sourcesink(i) / scaler(i)

              ! itegration of c over [dt', dt]
              c_int_pp = analytical_c_int(cmax, r, k, dt_p, dt)

              ! add to the fluxes
              dc_up(i) = dc_up(i) + i3*abs(q_int(i))*c_int_pp
              dc_down(i) = dc_down(i) + i4*q_int(i+1)*c_int_pp

              ! find final value
              c_next(i) = analytical_c(cmax, r, k, dt_p, dt)
              niter(i) = 2
          end if

        end if ! end distinction of first, last, and intermediate layers

      end if ! end distinction between having diffusion-advection or not

    end do ! end iteration thru soil layers

    ! we still need to catch the case when r is simply too negative
    ! in that case, we really need to reduce r (secondary mineral
    ! precipitation, etc.)
    ! (TBD)

    !print *, 'dc', (c_next - c_prev)*dz*theta
    !print *, 'q_int', q_int
    !print *, 'dc_up', dc_up
    !print *, 'dc_down', dc_down
    !print *, 'dc - outflux', (c_next - c_prev)*dz*theta + (dc_up + dc_down)
    !print *, sum(dc_up + dc_down)
    !print *, sum((c_next - c_prev)*dz*theta)
    !print *, sum((c_next - c_prev)*dz*theta + (dc_up + dc_down))

    ! calculate the net between self-outflow and inflow fluxes
    ! note the inflow fluxes need to be scalerd by soil moisture
    ! to get the correct concentration implications
    c_next(:nlevsoi-1) = c_next(:nlevsoi-1) + dc_up(2:nlevsoi)/theta(:nlevsoi-1)/dz(:nlevsoi-1)
    c_next(2:nlevsoi) = c_next(2:nlevsoi) + dc_down(:nlevsoi-1)/theta(2:nlevsoi)/dz(2:nlevsoi)

  end subroutine advection_diffusion

end program main
