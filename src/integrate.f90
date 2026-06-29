module odeint__integrate
  use odeint__config, only : wp
  use odeint__util, only : linspace
  use odeint__interfaces, only : dydt_abstract
  use odeint__interfaces, only : stepper_abstract
  use odeint__euler, only : euler
  use odeint__rk4, only : rk4
  implicit none

  private
  public :: integrate

contains

  subroutine integrate(dydt, y0, tstart, tend, nsteps, t, y, method)
    procedure(dydt_abstract) :: dydt
    real(wp), intent(in) :: y0(:)
    real(wp), intent(in) :: tstart
    real(wp), intent(in) :: tend
    integer, intent(in) :: nsteps
    real(wp), allocatable, intent(out) :: t(:)
    real(wp), allocatable, intent(out) :: y(:, :)
    character(len=*), intent(in), optional :: method

    character(:), allocatable :: method_

    if (present(method)) then
      method_ = method
    else
      method_ = "rk4"
    end if

    select case (method_)
      case ("euler")
        call run(dydt, euler, y0, tstart, tend, nsteps, t, y)
      case ("rk4")
        call run(dydt, rk4, y0, tstart, tend, nsteps, t, y)
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

end module odeint__integrate
