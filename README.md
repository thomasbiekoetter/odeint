# odeint

A Fortran library for integrating ordinary differential equations (ODEs), built with the [Fortran Package Manager (fpm)](https://fpm.fortran-lang.org/). Provides Euler and fourth-order Runge-Kutta (RK4) integrators with a uniform interface.

## Requirements

- [fpm](https://fpm.fortran-lang.org/) >= 0.13
- A Fortran compiler (gfortran >= 10 or Intel ifx/ifort)

## Installation

Clone the repository and build with fpm:

```sh
git clone https://gitlab.com/thomas.biekoetter/odeint.git
cd odeint
fpm build --profile='release'
```

Run the bundled tests to verify the build:

```sh
fpm test
```

Each test integrates exponential decay (dy/dt = -y, y(0) = 1) and prints the numerical solution, exact value, and absolute error at t = 5:

```
RK4 integration: dy/dt = -y,  y(0) = 1
  steps     : 100
  t         :   5.00
  y(t)      :  6.73795E-03
  exact     :  6.73795E-03
  error     :  1.90522E-09
```

### Quadruple precision

By default the library uses double precision (`real64`). To build with quadruple precision, pass the preprocessor flag:

```sh
fpm build --flag "-DQUAD"
fpm test  --flag "-DQUAD"
```

## Using odeint in another fpm project

Add odeint as a dependency in your project's `fpm.toml`:

```toml
[dependencies]
odeint = { git = "https://gitlab.com/thomas.biekoetter/odeint.git", profile="release" }
```

Then `use` whichever integrator module you need:

```fortran
use odeint__euler, only : integrate   ! Euler method
use odeint__rk4,   only : integrate   ! Runge-Kutta 4
use odeint__config, only : wp         ! working precision kind
```

## User interface

Both integrators expose the same `integrate` subroutine, so you can switch between methods by changing a single `use` statement.

### `integrate` subroutine

```fortran
subroutine integrate(dydt, y0, tstart, tend, nsteps, t, y)
```

| Argument | Intent | Type | Description |
|---|---|---|---|
| `dydt` | in | `procedure` | Right-hand side function f(y, t) |
| `y0(:)` | in | `real(wp)` | Initial state vector |
| `tstart` | in | `real(wp)` | Start time |
| `tend` | in | `real(wp)` | End time |
| `nsteps` | in | `integer` | Number of time steps (including t0) |
| `t(:)` | out | `real(wp), allocatable` | Time grid, length `nsteps` |
| `y(:,:)` | out | `real(wp), allocatable` | Solution, shape `(nsteps, size(y0))` |

The time grid is uniformly spaced: `t(1) = tstart`, `t(nsteps) = tend`, step size `h = (tend - tstart) / (nsteps - 1)`.

The solution array is indexed as `y(i, j)`, where `i` is the time index and `j` is the component index.

### `dydt` function signature

The right-hand side function must match this interface:

```fortran
function dydt(y, t) result(dydt_out)
  use odeint__config, only : wp
  real(wp), intent(in) :: y(:)      ! current state vector
  real(wp), intent(in) :: t         ! current time
  real(wp), allocatable :: dydt_out(:)  ! time derivative
end function
```

### Working precision

The kind parameter `wp` is exported from `odeint__config`. All real literals in your `dydt` function and initial conditions should use this kind to avoid silent precision mismatches:

```fortran
use odeint__config, only : wp
real(wp) :: x = 1.0e0_wp
```

## Example: exponential decay

The RK4 test program (`test/rk4.f90`) integrates dy/dt = -y with y(0) = 1 from t = 0 to t = 5:

```fortran
program odeint__test_rk4

  use odeint__config, only : wp
  use odeint__rk4, only : integrate

  implicit none

  real(wp), parameter :: y0(1) = [1.0e0_wp]
  real(wp), parameter :: tstart = 0.0e0_wp
  real(wp), parameter :: tend = 5.0e0_wp
  integer, parameter :: nsteps = 100
  real(wp), allocatable :: t(:)
  real(wp), allocatable :: y(:, :)

  call integrate(exp_decay, y0, tstart, tend, nsteps, t, y)

  write(*, '(A)')         "RK4 integration: dy/dt = -y,  y(0) = 1"
  write(*, '(A, I0)')     "  steps     : ", nsteps
  write(*, '(A, F6.2)')   "  t         : ", t(nsteps)
  write(*, '(A, ES12.5)') "  y(t)      : ", y(nsteps, 1)
  write(*, '(A, ES12.5)') "  exact     : ", exp(-t(nsteps))
  write(*, '(A, ES12.5)') "  error     : ", abs(y(nsteps, 1) - exp(-t(nsteps)))

contains

  function exp_decay(y, t) result(dydt)

      real(wp), intent(in) :: y(:)
      real(wp), intent(in) :: t
      real(wp), allocatable :: dydt(:)

      dydt = -y

  end function exp_decay

end program odeint__test_rk4
```

To solve a system of equations, use a state vector with more than one component and return a `dydt` array of the same length.

## Choosing a method

| Method | Module | Order | Recommended use |
|---|---|---|---|
| Euler | `odeint__euler` | 1 | Quick prototyping; requires many steps for accuracy |
| RK4 | `odeint__rk4` | 4 | General purpose; accurate with far fewer steps |

RK4 makes four function evaluations per step but converges as O(h⁴), so it typically needs far fewer steps than Euler to achieve the same accuracy. Prefer RK4 unless you have a specific reason to use Euler.

## License

This project is licensed under the [GNU Affero General Public License v3.0](LICENSE) (AGPL-3.0). Derivative works and projects that use this library over a network must also be released under the AGPL-3.0.
