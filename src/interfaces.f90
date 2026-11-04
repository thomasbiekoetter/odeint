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

  ! Quality-controlled adaptive step. Attempts a step of size `htry` from
  ! (t, y); on success advances `y` in place, reports the step actually taken
  ! (`hdid`) and a suggested next step size (`hnext`). Error is controlled to
  ! within max(atol, rtol*|y|) componentwise. If `ierr` is present, failures
  ! are reported through it (leaving `y` unchanged) instead of stopping the
  ! program. Any adaptive method (Cash-Karp, Dormand-Prince, Fehlberg, ...)
  ! can implement this and be driven by run_adaptive.
  abstract interface
    subroutine adaptive_stepper_abstract(  &
        dydt, t, y, htry, atol, rtol, &
        hdid, hnext, ierr)
      import :: wp, dydt_abstract
      implicit none
      procedure(dydt_abstract) :: dydt
      real(wp), intent(in) :: t
      real(wp), intent(inout) :: y(:)
      real(wp), intent(in) :: htry
      real(wp), intent(in) :: atol
      real(wp), intent(in) :: rtol
      real(wp), intent(out) :: hdid
      real(wp), intent(out) :: hnext
      integer, intent(out), optional :: ierr
    end subroutine adaptive_stepper_abstract
  end interface

end module odeint__interfaces
