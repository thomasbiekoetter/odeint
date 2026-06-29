module odeint__euler
  use odeint__config, only : wp
  use odeint__util, only : linspace
  use odeint__interfaces, only : dydt_abstract
  implicit none

  public :: euler

contains

  function euler(dydt, t, y, h) result(step)
    procedure(dydt_abstract) :: dydt
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: h
    real(wp), allocatable :: step(:)

    step = h * dydt(y, t)

  end function euler

end module odeint__euler
