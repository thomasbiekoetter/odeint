program odeint__test_rkck
  use odeint__config, only : wp
  use odeint__integrate, only : integrate

  implicit none

  real(wp), parameter :: y0(1) = [1.0e0_wp]
  real(wp), parameter :: tstart = 0.0e0_wp
  real(wp), parameter :: tend = 5.0e0_wp
  integer, parameter :: nsteps = 25
  real(wp), allocatable :: t(:), y(:, :)
  real(wp) :: err, exact

  ! ---- Same form as the euler/rk4 tests ----
  call integrate(exp_decay, y0, tstart, tend, nsteps, t, y, method="rkck")
  write(*, '(A)')         "RKCK integration: dy/dt = -y,  y(0) = 1"
  write(*, '(A, I0)')     "  steps     : ", nsteps
  write(*, '(A, F6.2)')   "  t         : ", t(nsteps)
  write(*, '(A, ES12.5)') "  y(t)      : ", y(nsteps, 1)
  write(*, '(A, ES12.5)') "  exact     : ", exp(-t(nsteps))
  write(*, '(A, ES12.5)') "  error     : ", abs(y(nsteps, 1) - exp(-t(nsteps)))

  ! ---- Test 1: exponential decay dy/dt = -y, y(0)=1 on [0,5] ----
  call integrate(exp_decay, [1.0e0_wp], 0.0e0_wp, 5.0e0_wp, 6, t, y, &
                 method="rkck", rtol=1.0e-10_wp, atol=1.0e-12_wp)
  exact = exp(-5.0e0_wp)
  err = abs(y(size(t), 1) - exact)
  write(*, '(A)')         "rkck: dy/dt = -y,  y(0)=1"
  write(*, '(A,ES13.6)')  "  y(5)   : ", y(size(t), 1)
  write(*, '(A,ES13.6)')  "  exact  : ", exact
  write(*, '(A,ES13.6)')  "  error  : ", err
  write(*, '(A,6F9.4)')   "  grid   : ", t
  if (err > 1.0e-8_wp) error stop "exp_decay: error above tolerance"

  ! ---- Test 2: harmonic oscillator y'' = -y  ->  [y, v]' = [v, -y] ----
  ! y(0)=1, v(0)=0  =>  y(t)=cos t, v(t)=-sin t.  Integrate to t=10.
  call integrate(sho, [1.0e0_wp, 0.0e0_wp], 0.0e0_wp, 10.0e0_wp, 3, t, y, &
                 method="rkck", rtol=1.0e-9_wp, atol=1.0e-12_wp)
  err = max(abs(y(size(t), 1) - cos(10.0e0_wp)), &
            abs(y(size(t), 2) + sin(10.0e0_wp)))
  write(*, '(A)')         "rkck: harmonic oscillator, t=10"
  write(*, '(A,2ES13.6)') "  y,v    : ", y(size(t), 1), y(size(t), 2)
  write(*, '(A,2ES13.6)') "  exact  : ", cos(10.0e0_wp), -sin(10.0e0_wp)
  write(*, '(A,ES13.6)')  "  error  : ", err
  if (err > 1.0e-6_wp) error stop "sho: error above tolerance"

  ! ---- Test 3: lands exactly on the uniform output grid ----
  call integrate(exp_decay, [1.0e0_wp], 2.0e0_wp, 7.0e0_wp, 11, t, y, method="rkck")
  if (abs(t(1) - 2.0e0_wp) > 1.0e-14_wp) error stop "grid: wrong start"
  if (abs(t(11) - 7.0e0_wp) > 1.0e-13_wp) error stop "grid: wrong end"
  if (abs((t(2) - t(1)) - 0.5e0_wp) > 1.0e-13_wp) error stop "grid: not uniform"
  write(*, '(A)') "rkck: output grid is exactly uniform [OK]"

  write(*, '(A)') "All rkck tests passed."

contains

  function exp_decay(y, t) result(dydt)
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: t
    real(wp), allocatable :: dydt(:)
    dydt = -y
  end function exp_decay

  function sho(y, t) result(dydt)
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: t
    real(wp), allocatable :: dydt(:)
    dydt = [y(2), -y(1)]
  end function sho

end program odeint__test_rkck
