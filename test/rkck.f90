program odeint__test_rkck
  use odeint__config, only : wp
  use odeint__integrate, only : integrate
  use odeint__integrate, only : rkqs_ok
  use odeint__integrate, only : rkqs_err_underflow
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

  implicit none

  real(wp), parameter :: y0(1) = [1.0e0_wp]
  real(wp), parameter :: tstart = 0.0e0_wp
  real(wp), parameter :: tend = 5.0e0_wp
  real(wp), allocatable :: t(:), y(:, :)
  real(wp) :: err, exact
  integer :: ierr

  ! ---- Test 1: exponential decay dy/dt = -y, y(0)=1 on [0,5] ----
  call integrate(exp_decay, [1.0e0_wp], tstart, tend, 6, t, y, &
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

  ! ---- Test 4: backward integration dy/dt = -y, y(0)=1 on [0,-5] ----
  ! Exact solution y(-5) = exp(5).
  call integrate(exp_decay, [1.0e0_wp], 0.0e0_wp, -5.0e0_wp, 6, t, y, &
                 method="rkck", rtol=1.0e-10_wp, atol=1.0e-12_wp)
  exact = exp(5.0e0_wp)
  err = abs(y(size(t), 1) - exact) / exact
  write(*, '(A)')         "rkck: backward, dy/dt = -y,  y(0)=1"
  write(*, '(A,ES13.6)')  "  y(-5)  : ", y(size(t), 1)
  write(*, '(A,ES13.6)')  "  exact  : ", exact
  write(*, '(A,ES13.6)')  "  relerr : ", err
  write(*, '(A,6F9.4)')   "  grid   : ", t
  if (abs(t(size(t)) + 5.0e0_wp) > 1.0e-13_wp) error stop "backward: wrong end"
  if (err > 1.0e-8_wp) error stop "backward: error above tolerance"

  ! ---- Test 5: failure is reported through ierr, not error stop ----
  ! The derivative oscillates far below the step resolution, so the step
  ! size underflows; integrate must return instead of stopping.
  call integrate(wild, [0.0e0_wp], 1.0e0_wp, 2.0e0_wp, 5, t, y, &
                 method="rkck", ierr=ierr)
  write(*, '(A,I0)') "rkck: unresolvable derivative, ierr = ", ierr
  if (ierr /= rkqs_err_underflow) error stop "wild: expected underflow"
  if (y(1, 1) /= 0.0e0_wp) error stop "wild: initial row overwritten"
  if (.not. all(ieee_is_nan(y(2:, 1)))) error stop "wild: rows not NaN"

  ! ---- Test 6: finite-time blow-up dy/dt = y**2, y(0)=1 on [0,2] ----
  ! The exact solution y = 1/(1-t) diverges at t = 1, so the step size
  ! shrinks to zero there. Output grid 0, 2/3, 4/3, 2 (no point on the
  ! singularity itself): the point before t = 1 must be kept.
  call integrate(blowup, [1.0e0_wp], 0.0e0_wp, 2.0e0_wp, 4, t, y, &
                 method="rkck", ierr=ierr)
  write(*, '(A,I0)') "rkck: finite-time blow-up, ierr = ", ierr
  if (ierr /= rkqs_err_underflow) error stop "blowup: expected underflow"
  if (abs(y(2, 1) - 3.0e0_wp) > 1.0e-5_wp)  &
    error stop "blowup: y(2/3) not kept"
  if (.not. all(ieee_is_nan(y(3:, 1)))) error stop "blowup: rows not NaN"

  ! ierr is rkqs_ok on success.
  call integrate(exp_decay, [1.0e0_wp], tstart, tend, 6, t, y, &
                 method="rkck", ierr=ierr)
  if (ierr /= rkqs_ok) error stop "exp_decay: unexpected ierr"

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

  function blowup(y, t) result(dydt)
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: t
    real(wp), allocatable :: dydt(:)
    dydt = y**2
  end function blowup

  function wild(y, t) result(dydt)
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: t
    real(wp), allocatable :: dydt(:)
    dydt = [1.0e280_wp * sin(1.0e20_wp * t)]
  end function wild

end program odeint__test_rkck
