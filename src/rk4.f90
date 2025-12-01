module odeint__rk4

  use odeint__config, only : wp
  use odeint__util, only : linspace

  implicit none

  private :: rk4
  public :: integrate

  abstract interface
    function dydt_abstract(y, t) result(dydt)
      import :: wp
      implicit none
      real(wp), intent(in) :: y(:)
      real(wp), intent(in) :: t
      real(wp), allocatable :: dydt(:)
    end function dydt_abstract
  end interface

contains

  subroutine integrate(dydt, y0, tstart, tend, nsteps, t, y)

    procedure(dydt_abstract) :: dydt
    real(wp), intent(in) :: y0(:)
    real(wp), intent(in) :: tstart
    real(wp), intent(in) :: tend
    integer, intent(in) :: nsteps
    real(wp), allocatable, intent(out) :: t(:)
    real(wp), allocatable, intent(out) :: y(:, :)

    integer :: n
    integer :: i
    real(wp) :: h
    real(wp), allocatable :: step(:)

    n = size(y0)
    allocate(y(nsteps, n))
    y(1, :) = y0

    t = linspace(tstart, tend, nsteps)
    h = t(2) - t(1)

    do i = 1, nsteps - 1
      step = rk4(dydt, t(i), y(i, :), h)
      y(i + 1, :) = y(i, :) + step
    end do

  end subroutine integrate

  function rk4(dydt, t, y, h) result(step)

    procedure(dydt_abstract) :: dydt
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: h
    real(wp), allocatable :: step(:)

    integer :: n
    real(wp), allocatable :: k1(:)
    real(wp), allocatable :: k2(:)
    real(wp), allocatable :: k3(:)
    real(wp), allocatable :: k4(:)
    real(wp), allocatable :: yt(:)

    n = size(y)
    allocate(k1(n))
    allocate(k2(n))
    allocate(k3(n))
    allocate(k4(n))
    allocate(yt(n))
    allocate(step(n))

    k1 = dydt(y, t)
    yt = y + 0.5e0_wp * h * k1
    k2 = dydt(yt, t + 0.5e0_wp * h)
    yt = y + 0.5e0_wp * h * k2
    k3 = dydt(yt, t + 0.5e0_wp * h)
    yt = y + h * k3
    k4 = dydt(yt, t + h)
    step = h * (k1 + 2.0e0_wp * k2 + 2.0e0_wp * k3 + k4) / 6.0e0_wp

  end function

end module odeint__rk4
