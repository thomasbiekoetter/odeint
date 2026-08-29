module odeint__rkck
  use odeint__config, only : wp
  use odeint__interfaces, only : dydt_abstract
  implicit none

  private
  public :: rkck
  public :: rkqs
  public :: rkqs_vec
  public :: rkqs_ok
  public :: rkqs_err_underflow
  public :: rkqs_err_max_tries
  public :: rkqs_err_overflow

  ! Status codes reported by rkqs_vec through `ierr`.
  integer, parameter :: rkqs_ok = 0
  integer, parameter :: rkqs_err_underflow = 1  ! step size rounds down to zero
  integer, parameter :: rkqs_err_max_tries = 2  ! too many step reductions
  integer, parameter :: rkqs_err_overflow = 3   ! step size grows without bound

  ! Cash-Karp parameters
  ! Nodes.
  real(wp), parameter :: a2 = 0.2e0_wp
  real(wp), parameter :: a3 = 0.3e0_wp
  real(wp), parameter :: a4 = 0.6e0_wp
  real(wp), parameter :: a5 = 1.0e0_wp
  real(wp), parameter :: a6 = 0.875e0_wp
  ! Runge-Kutta matrix.
  real(wp), parameter :: b21 = 0.2e0_wp
  real(wp), parameter :: b31 = 3.0e0_wp / 40.0e0_wp
  real(wp), parameter :: b32 = 9.0e0_wp / 40.0e0_wp
  real(wp), parameter :: b41 = 0.3e0_wp
  real(wp), parameter :: b42 = -0.9e0_wp
  real(wp), parameter :: b43 = 1.2e0_wp
  real(wp), parameter :: b51 = -11.0e0_wp / 54.0e0_wp
  real(wp), parameter :: b52 = 2.5e0_wp
  real(wp), parameter :: b53 = -70.0e0_wp / 27.0e0_wp
  real(wp), parameter :: b54 = 35.0e0_wp / 27.0e0_wp
  real(wp), parameter :: b61 = 1631.0e0_wp / 55296.0e0_wp
  real(wp), parameter :: b62 = 175.0e0_wp / 512.0e0_wp
  real(wp), parameter :: b63 = 575.0e0_wp / 13824.0e0_wp
  real(wp), parameter :: b64 = 44275.0e0_wp / 110592.0e0_wp
  real(wp), parameter :: b65 = 253.0e0_wp / 4096.0e0_wp
  ! 5th order weights.
  real(wp), parameter :: c1 = 37.0e0_wp / 378.0e0_wp
  real(wp), parameter :: c3 = 250.0e0_wp / 621.0e0_wp
  real(wp), parameter :: c4 = 125.0e0_wp / 594.0e0_wp
  real(wp), parameter :: c6 = 512.0e0_wp / 1771.0e0_wp
  ! Error weights: (5th order) - (embedded 4th order).
  real(wp), parameter :: dc1 = c1 - 2825.0e0_wp / 27648.0e0_wp
  real(wp), parameter :: dc3 = c3 - 18575.0e0_wp / 48384.0e0_wp
  real(wp), parameter :: dc4 = c4 - 13525.0e0_wp / 55296.0e0_wp
  real(wp), parameter :: dc5 = -277.0e0_wp / 14336.0e0_wp
  real(wp), parameter :: dc6 = c6 - 0.25e0_wp

  ! Step-size controller constants
  real(wp), parameter :: safety = 0.9e0_wp
  real(wp), parameter :: pgrow = -0.2e0_wp    ! 1/5 growth exponent
  real(wp), parameter :: pshrnk = -0.25e0_wp  ! 1/4 shrink exponent
  real(wp), parameter :: errcon = 1.89e-4_wp  ! (5/safety)**(1/pgrow)
  integer, parameter :: max_tries = 100       ! retries before giving up

contains

  ! Single embedded Cash-Karp step.
  ! Returns the 5th-order increment `step` (so y_new = y + step) and the
  ! componentwise error estimate `yerr` = (5th order) - (4th order).
  ! If the derivative at (t, y) is already known it can be passed as
  ! `dydt0`, saving one evaluation of `dydt`.
  subroutine rkck(dydt, t, y, h, step, yerr, dydt0)
    procedure(dydt_abstract) :: dydt
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: h
    real(wp), allocatable, intent(out) :: step(:)
    real(wp), allocatable, intent(out) :: yerr(:)
    real(wp), intent(in), optional :: dydt0(:)

    real(wp), allocatable :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:)

    if (present(dydt0)) then
      k1 = dydt0
    else
      k1 = dydt(y, t)
    end if
    k2 = dydt(y + h * (b21 * k1), t + a2 * h)
    k3 = dydt(y + h * (b31 * k1 + b32 * k2), t + a3 * h)
    k4 = dydt(y + h * (b41 * k1 + b42 * k2 + b43 * k3), t + a4 * h)
    k5 = dydt(y + h * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4), t + a5 * h)
    k6 = dydt(y + h * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5), &
              t + a6 * h)

    step = h * (c1 * k1 + c3 * k3 + c4 * k4 + c6 * k6)
    yerr = h * (dc1 * k1 + dc3 * k3 + dc4 * k4 + dc5 * k5 + dc6 * k6)

  end subroutine rkck

  ! Quality-controlled step. Attempts a step of size `htry` from (t, y),
  ! shrinking until the scaled error is within tolerance, then advances `y`
  ! in place. Reports the step actually taken (`hdid`) and a suggested next
  ! step size (`hnext`). Tolerance per component is max(atol, rtol*|y|), so
  ! the accepted error is below the larger of the absolute and relative bound.
  subroutine rkqs(dydt, t, y, htry, atol, rtol, hdid, hnext)
    procedure(dydt_abstract) :: dydt
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: htry
    real(wp), intent(in) :: atol
    real(wp), intent(in) :: rtol
    real(wp), intent(out) :: hdid
    real(wp), intent(out) :: hnext

    real(wp), allocatable :: step(:)
    integer :: ierr

    call rkqs_vec(dydt, t, y, htry, spread(atol, 1, size(y)),  &
                  spread(rtol, 1, size(y)), step, hdid, hnext, ierr)
    if (ierr == rkqs_err_underflow) error stop "rkqs: step size underflow"
    if (ierr == rkqs_err_overflow) error stop "rkqs: step size overflow"
    if (ierr /= rkqs_ok) error stop "rkqs: too many step reductions"
    y = y + step

  end subroutine rkqs

  ! Quality-controlled step with per-component (vector) tolerances.
  ! Like rkqs, but:
  !   * `atol` and `rtol` are arrays, so each component can carry its own
  !     tolerance (e.g. position and velocity living on different scales);
  !     the tolerance for component i is max(atol(i), rtol(i) * |y(i)|).
  !   * `y` is not advanced in place; the accepted 5th-order increment is
  !     returned in `step` (so y_new = y + step).
  !   * failures are reported through `ierr` (rkqs_ok, rkqs_err_underflow,
  !     rkqs_err_max_tries, rkqs_err_overflow) instead of stopping the
  !     program, so library
  !     callers can handle them. On ierr /= rkqs_ok the outputs `step`,
  !     `hdid` and `hnext` are meaningless.
  !   * if the derivative at (t, y) is already known it can be passed as
  !     `dydt0`, saving one evaluation of `dydt`.
  subroutine rkqs_vec(dydt, t, y, htry, atol, rtol, step, hdid, hnext,  &
                      ierr, dydt0)
    procedure(dydt_abstract) :: dydt
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: htry
    real(wp), intent(in) :: atol(:)
    real(wp), intent(in) :: rtol(:)
    real(wp), allocatable, intent(out) :: step(:)
    real(wp), intent(out) :: hdid
    real(wp), intent(out) :: hnext
    integer, intent(out) :: ierr
    real(wp), intent(in), optional :: dydt0(:)

    real(wp), allocatable :: k1(:), yerr(:), scaler(:)
    real(wp) :: h, htemp, errmax
    integer :: tries

    ierr = rkqs_ok
    hdid = 0.0e0_wp
    hnext = htry

    if (present(dydt0)) then
      k1 = dydt0
    else
      k1 = dydt(y, t)
    end if
    scaler = max(atol, rtol * abs(y))

    h = htry
    do tries = 1, max_tries
      call rkck(dydt, t, y, h, step, yerr, dydt0=k1)
      errmax = maxval(abs(yerr / scaler))
      if (errmax <= 1.0e0_wp) exit
      ! Error too large: shrink, but by no more than a factor of 10.
      htemp = safety * h * errmax ** pshrnk
      h = sign(max(abs(htemp), 0.1e0_wp * abs(h)), h)
      if (t + h == t) then
        ierr = rkqs_err_underflow
        return
      end if
      if (tries == max_tries) then
        ierr = rkqs_err_max_tries
        return
      end if
    end do

    ! Step accepted: propose the next step size (grow by at most factor 5).
    if (errmax > errcon) then
      hnext = safety * h * errmax ** pgrow
    else
      if (abs(h) > 0.2e0_wp * huge(h)) then
        ! 5*h would overflow. This means a caller kept feeding hnext back
        ! in while its stopping condition never triggered (the step size
        ! has grown without bound); report it instead of raising an FPE.
        ierr = rkqs_err_overflow
        return
      end if
      hnext = 5.0e0_wp * h
    end if

    hdid = h

  end subroutine rkqs_vec

end module odeint__rkck
