module odeint__rk4
  use odeint__config, only : wp
  use odeint__util, only : linspace
  use odeint__interfaces, only : dydt_abstract
  implicit none

  public :: rk4

contains

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

  end function rk4

end module odeint__rk4
