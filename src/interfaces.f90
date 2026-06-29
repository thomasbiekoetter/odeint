module odeint__interfaces
  use odeint__config, only : wp
  implicit none

  abstract interface
    function dydt_abstract(y, t) result(dydt)
      import :: wp
      implicit none
      real(wp), intent(in) :: y(:)
      real(wp), intent(in) :: t
      real(wp), allocatable :: dydt(:)
    end function dydt_abstract
  end interface

  abstract interface
    function stepper_abstract(dydt, t, y, h) result(step)
      import :: wp, dydt_abstract
      implicit none
      procedure(dydt_abstract) :: dydt
      real(wp), intent(in) :: y(:)
      real(wp), intent(in) :: t
      real(wp), intent(in) :: h
      real(wp), allocatable :: step(:)
    end function stepper_abstract
  end interface

end module odeint__interfaces
