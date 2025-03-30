module my_module
  implicit none
  ! Define a double precision kind for real numbers.
  integer, parameter :: r8 = selected_real_kind(15, 307)
contains

  function u_pdf(u, theta) result(f)
    ! Gamma PDF with k = 0.5 transformed by x = u**2
    !   pdf(x)dx = pdf(u^2) * 2u du
    !   Note: gamma(0.5) = sqrt(pi)
    implicit none
    real(r8), intent(in) :: u, theta
    real(r8) :: f
    f = 2.0_r8 / sqrt(theta * 3.1415926535_r8) * exp(-u**2 / theta)
  end function u_pdf

  function get_ssa(gs) result(a)
    ! Calculate the specific surface area given grain size in um.
    implicit none
    real(r8), intent(in) :: gs
    real(r8) :: a
    a = 69.18_r8 * (gs ** (-1.24_r8))  ! units: m^2 g-1
  end function get_ssa

  subroutine compute_ssa(forc_gra, theta, ssa)
    ! Compute the specific surface area (ssa) by numerically integrating
    ! the gamma PDF weighted by the get_ssa function using Simpson's rule.
    implicit none
    real(r8), intent(in) :: forc_gra, theta
    real(r8), intent(out) :: ssa
    real(r8) :: num, prob
    real(r8) :: u_min, u_max
    real(r8) :: u, du
    integer  :: n_int, ii

    ! Set the number of intervals for Simpson's rule.
    n_int = 1000

    ! Define integration limits based on forc_gra.
    u_min = sqrt(0.01_r8 * forc_gra)
    u_max = sqrt(100.0_r8 * forc_gra)
    du = (u_max - u_min) / n_int

    ! Initialize numerator and denominator with the first endpoint.
    num = get_ssa(u_min**2) * u_pdf(u_min, theta)
    prob = u_pdf(u_min, theta)

    ! Simpson's rule integration loop.
    do ii = 1, n_int-1
       u = u_min + du * ii

       !! print*, "u_pdf(u, theta) ", u_pdf(u, theta)

       if (mod(ii,2) == 1) then
          num = num + 4.0_r8 * get_ssa(u**2) * u_pdf(u, theta)
          prob = prob + 4.0_r8 * u_pdf(u, theta)
       else
          num = num + 2.0_r8 * get_ssa(u**2) * u_pdf(u, theta)
          prob = prob + 2.0_r8 * u_pdf(u, theta)
       end if
    end do

    ! Add the last endpoint.
    num = num + get_ssa(u_max**2) * u_pdf(u_max, theta)
    prob = prob + u_pdf(u_max, theta)

    ! Final Simpson's rule weights.
    num = num * du / 3.0_r8
    prob = prob * du / 3.0_r8

    ! Calculate the specific surface area.
    ssa = num / prob
  end subroutine compute_ssa

end module my_module

program test_ssa
  use my_module
  implicit none
  real(r8) :: forc_gra, theta, ssa

  ! Set test values for forc_gra and theta.
  forc_gra = 107._r8    ! Example grain parameter (can be adjusted)
  theta    = 4.394_r8   ! Example parameter for the gamma PDF

  ! Call the subroutine to compute the specific surface area.
  call compute_ssa(forc_gra, theta, ssa)

  ! Output the computed specific surface area.
  print*, "Computed specific surface area (ssa): ", ssa

end program test_ssa
