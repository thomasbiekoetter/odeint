program odeint__test_rkqs_vec
  use odeint__config, only : wp
  use odeint__rkck, only : rkqs
  use odeint__rkck, only : rkqs_vec
  use odeint__rkck, only : rkqs_ok
  use odeint__rkck, only : rkqs_err_underflow

  implicit none

  real(wp), allocatable :: y(:), step(:)
  real(wp) :: hdid, hnext, hdid_ref, hnext_ref, hdid_loose, hdid_tight
  integer :: ierr

  ! ---- Test 1: uniform vector tolerances reproduce rkqs exactly ----
  ! One step of the harmonic oscillator; rkqs advances y in place, rkqs_vec
  ! returns the increment. Same tolerances => identical results.
  y = [1.0e0_wp, 0.0e0_wp]
  call rkqs(sho, 0.0e0_wp, y, 0.5e0_wp, 1.0e-12_wp, 1.0e-9_wp,  &
            hdid_ref, hnext_ref)
  call rkqs_vec(sho, 0.0e0_wp, [1.0e0_wp, 0.0e0_wp], 0.5e0_wp,  &
                [1.0e-12_wp, 1.0e-12_wp], [1.0e-9_wp, 1.0e-9_wp],  &
                step, hdid, hnext, ierr)
  write(*, '(A)')         "rkqs_vec: one step of the harmonic oscillator"
  write(*, '(A,2ES13.6)') "  rkqs     : ", y
  write(*, '(A,2ES13.6)') "  rkqs_vec : ", [1.0e0_wp, 0.0e0_wp] + step
  if (ierr /= rkqs_ok) error stop "uniform: unexpected ierr"
  if (any([1.0e0_wp, 0.0e0_wp] + step /= y)) error stop "uniform: step differs"
  if (hdid /= hdid_ref .or. hnext /= hnext_ref)  &
    error stop "uniform: step sizes differ"

  ! ---- Test 2: precomputed derivative gives identical results ----
  call rkqs_vec(sho, 0.0e0_wp, [1.0e0_wp, 0.0e0_wp], 0.5e0_wp,  &
                [1.0e-12_wp, 1.0e-12_wp], [1.0e-9_wp, 1.0e-9_wp],  &
                step, hdid, hnext, ierr,  &
                dydt0=sho([1.0e0_wp, 0.0e0_wp], 0.0e0_wp))
  if (ierr /= rkqs_ok) error stop "dydt0: unexpected ierr"
  if (any([1.0e0_wp, 0.0e0_wp] + step /= y)) error stop "dydt0: step differs"
  if (hdid /= hdid_ref .or. hnext /= hnext_ref)  &
    error stop "dydt0: step sizes differ"
  write(*, '(A)') "rkqs_vec: precomputed dydt0 reproduces the step [OK]"

  ! ---- Test 3: per-component tolerances actually bite ----
  ! Two decoupled decays with very different component scales. Tightening
  ! the absolute tolerance of the small component alone must force a
  ! smaller accepted step.
  call rkqs_vec(two_scales, 0.0e0_wp, [1.0e0_wp, 1.0e-4_wp], 1.0e0_wp,  &
                [1.0e-9_wp, 1.0e-9_wp], [0.0e0_wp, 0.0e0_wp],  &
                step, hdid_loose, hnext, ierr)
  if (ierr /= rkqs_ok) error stop "two_scales loose: unexpected ierr"
  call rkqs_vec(two_scales, 0.0e0_wp, [1.0e0_wp, 1.0e-4_wp], 1.0e0_wp,  &
                [1.0e-9_wp, 1.0e-15_wp], [0.0e0_wp, 0.0e0_wp],  &
                step, hdid_tight, hnext, ierr)
  if (ierr /= rkqs_ok) error stop "two_scales tight: unexpected ierr"
  write(*, '(A)')        "rkqs_vec: per-component absolute tolerances"
  write(*, '(A,ES13.6)') "  hdid, atol = [1e-9, 1e-9]  : ", hdid_loose
  write(*, '(A,ES13.6)') "  hdid, atol = [1e-9, 1e-15] : ", hdid_tight
  if (hdid_tight >= 0.9e0_wp * hdid_loose)  &
    error stop "two_scales: tight tolerance did not shrink the step"

  ! ---- Test 4: failure is reported through ierr, not error stop ----
  ! A derivative oscillating far below the step resolution keeps the error
  ! estimate enormous at every attempted step size, so the controller
  ! shrinks h until t + h == t.
  call rkqs_vec(wild, 1.0e0_wp, [0.0e0_wp], 0.1e0_wp,  &
                [1.0e-9_wp], [1.0e-6_wp], step, hdid, hnext, ierr)
  write(*, '(A,I0)') "rkqs_vec: unresolvable derivative, ierr = ", ierr
  if (ierr /= rkqs_err_underflow)  &
    error stop "wild: expected rkqs_err_underflow"

  write(*, '(A)') "All rkqs_vec tests passed."

contains

  function sho(y, t) result(dydt)
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: t
    real(wp), allocatable :: dydt(:)
    dydt = [y(2), -y(1)]
  end function sho

  function two_scales(y, t) result(dydt)
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: t
    real(wp), allocatable :: dydt(:)
    dydt = [-y(1), -y(2)]
  end function two_scales

  function wild(y, t) result(dydt)
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: t
    real(wp), allocatable :: dydt(:)
    dydt = [1.0e280_wp * sin(1.0e20_wp * t)]
  end function wild

end program odeint__test_rkqs_vec
