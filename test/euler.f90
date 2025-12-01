program odeint__test_euler

  use odeint__config, only : wp
  use odeint__euler, only : integrate

  implicit none

  real(wp), parameter :: y0(1) = [1.0e0_wp]
  real(wp), parameter :: tstart = 0.0e0_wp
  real(wp), parameter :: tend = 5.0e0_wp
  integer, parameter :: nsteps = 10000
  real(wp), allocatable :: t(:)
  real(wp), allocatable :: y(:, :)

  call integrate(exp_decay, y0, tstart, tend, nsteps, t, y)

  write(*,*) t(nsteps), y(nsteps, :)

contains

  ! Example: exponential decay
  function exp_decay(y, t) result(dydt)

      real(wp), intent(in) :: y(:)
      real(wp), intent(in) :: t
      real(wp), allocatable :: dydt(:)

      dydt = -y

  end function exp_decay

end program odeint__test_euler
