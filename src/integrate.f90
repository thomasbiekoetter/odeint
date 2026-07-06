module odeint__integrate
  use odeint__config, only : wp
  use odeint__util, only : linspace
  use odeint__interfaces, only : dydt_abstract
  use odeint__interfaces, only : stepper_abstract
  use odeint__interfaces, only : adaptive_stepper_abstract
  use odeint__euler, only : euler
  use odeint__rk4, only : rk4
  use odeint__rkck, only : rkqs
  implicit none

  private
  public :: integrate

contains

  subroutine integrate(dydt, y0, tstart, tend, nsteps, t, y, method, atol, rtol)
    procedure(dydt_abstract) :: dydt
    real(wp), intent(in) :: y0(:)
    real(wp), intent(in) :: tstart
    real(wp), intent(in) :: tend
    integer, intent(in) :: nsteps
    real(wp), allocatable, intent(out) :: t(:)
    real(wp), allocatable, intent(out) :: y(:, :)
    character(len=*), intent(in), optional :: method
    real(wp), intent(in), optional :: atol
    real(wp), intent(in), optional :: rtol

    character(:), allocatable :: method_
    real(wp) :: atol_
    real(wp) :: rtol_

    if (present(method)) then
      method_ = method
    else
      method_ = "rk4"
    end if

    atol_ = 1.0e-9_wp
    if (present(atol)) atol_ = atol
    rtol_ = 1.0e-6_wp
    if (present(rtol)) rtol_ = rtol

    select case (method_)
      case ("euler")
        call run(dydt, euler, y0, tstart, tend, nsteps, t, y)
      case ("rk4")
        call run(dydt, rk4, y0, tstart, tend, nsteps, t, y)
      case ("rkck")
        call run_adaptive(  &
          dydt, rkqs, y0, tstart, tend, nsteps, t, y, atol_, rtol_)
      case default
        error stop "integrate: unknown method '" // method_ // "'"
    end select

  end subroutine integrate

  subroutine run(dydt, step, y0, tstart, tend, nsteps, t, y)
    procedure(dydt_abstract) :: dydt
    procedure(stepper_abstract) :: step
    real(wp), intent(in) :: y0(:)
    real(wp), intent(in) :: tstart
    real(wp), intent(in) :: tend
    integer, intent(in) :: nsteps
    real(wp), allocatable, intent(out) :: t(:)
    real(wp), allocatable, intent(out) :: y(:, :)

    integer :: n
    integer :: i
    real(wp) :: h

    n = size(y0)
    allocate(y(nsteps, n))
    y(1, :) = y0

    t = linspace(tstart, tend, nsteps)
    h = t(2) - t(1)

    do i = 1, nsteps - 1
      y(i + 1, :) = y(i, :) + step(dydt, t(i), y(i, :), h)
    end do

  end subroutine run

  ! Takes any quality-controlled step `step_qc` (rkqs for
  ! Cash-Karp, or a future Dormand-Prince/Fehlberg equivalent) and substeps
  ! with error control, but stores the solution on the same uniform grid as
  ! `run`: it lands exactly on each output point t(i+1). The number of
  ! internal steps is chosen by error control; `nsteps` only sets how many
  ! output samples are returned.
  subroutine run_adaptive(  &
      dydt, step_qc, y0, tstart, tend, nsteps, t, y, atol, rtol)
    procedure(dydt_abstract) :: dydt
    procedure(adaptive_stepper_abstract) :: step_qc
    real(wp), intent(in) :: y0(:)
    real(wp), intent(in) :: tstart
    real(wp), intent(in) :: tend
    integer, intent(in) :: nsteps
    real(wp), allocatable, intent(out) :: t(:)
    real(wp), allocatable, intent(out) :: y(:, :)
    real(wp), intent(in) :: atol
    real(wp), intent(in) :: rtol

    integer :: n
    integer :: i
    real(wp) :: tcur
    real(wp) :: h        ! carried "natural" step size across output points
    real(wp) :: hclamp
    real(wp) :: hdid
    real(wp) :: hnext
    real(wp), allocatable :: ycur(:)

    n = size(y0)
    allocate(y(nsteps, n))
    y(1, :) = y0

    t = linspace(tstart, tend, nsteps)

    ! Initial step guess: the output spacing. step_qc adapts from here.
    h = t(2) - t(1)
    ycur = y0

    do i = 1, nsteps - 1
      tcur = t(i)
      do while (tcur < t(i + 1))
        if (tcur + h >= t(i + 1)) then
          ! Clamp the final substep to land exactly on the output point.
          hclamp = t(i + 1) - tcur
          call step_qc(dydt, tcur, ycur, hclamp, atol, rtol, hdid, hnext)
          tcur = tcur + hdid
          ! If the step had to shrink for accuracy we did not reach the point
          ! yet, so adopt hnext. If the clamped step was accepted whole, keep
          ! the carried `h` so the next interval is not stuck at a tiny step.
          if (hdid < hclamp) h = hnext
        else
          call step_qc(dydt, tcur, ycur, h, atol, rtol, hdid, hnext)
          tcur = tcur + hdid
          h = hnext
        end if
      end do
      y(i + 1, :) = ycur
    end do

  end subroutine run_adaptive

end module odeint__integrate
