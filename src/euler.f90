module odeint__euler

  use odeint__config, only : wp
  use odeint__util, only : linspace

  implicit none

  private :: linspace

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

    n = size(y0)
    allocate(y(nsteps, n))
    y(1, :) = y0

    t = linspace(tstart, tend, nsteps)
    h = t(2) - t(1)

    do i = 1, nsteps - 1
      y(i + 1, :) = y(i, :) + h * dydt(y(i, :), t(i))
    end do

  end subroutine integrate

end module odeint__euler
