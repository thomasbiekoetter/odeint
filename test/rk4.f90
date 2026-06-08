program odeint__test_rk4

  use odeint__config, only : wp
  use odeint__rk4, only : integrate

  implicit none

  real(wp), parameter :: y0(1) = [1.0e0_wp]
  real(wp), parameter :: tstart = 0.0e0_wp
  real(wp), parameter :: tend = 5.0e0_wp
  integer, parameter :: nsteps = 100
  real(wp), allocatable :: t(:)
  real(wp), allocatable :: y(:, :)

  call integrate(exp_decay, y0, tstart, tend, nsteps, t, y)

  write(*, '(A)')         "RK4 integration: dy/dt = -y,  y(0) = 1"
  write(*, '(A, I0)')     "  steps     : ", nsteps
  write(*, '(A, F6.2)')   "  t         : ", t(nsteps)
  write(*, '(A, ES12.5)') "  y(t)      : ", y(nsteps, 1)
  write(*, '(A, ES12.5)') "  exact     : ", exp(-t(nsteps))
  write(*, '(A, ES12.5)') "  error     : ", abs(y(nsteps, 1) - exp(-t(nsteps)))

contains

  ! Example: exponential decay
  function exp_decay(y, t) result(dydt)

      real(wp), intent(in) :: y(:)
      real(wp), intent(in) :: t
      real(wp), allocatable :: dydt(:)

      dydt = -y

  end function exp_decay

end program odeint__test_rk4
