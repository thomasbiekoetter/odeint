# odeint

A Fortran library for integrating ordinary differential equations (ODEs), built with the [Fortran Package Manager (fpm)](https://fpm.fortran-lang.org/). Provides fixed-step Euler and fourth-order Runge-Kutta (RK4) integrators, plus an adaptive Cash-Karp Runge-Kutta method (RK45) with automatic step-size control, all behind a uniform interface.

## Requirements

- [fpm](https://fpm.fortran-lang.org/) >= 0.13
- A Fortran compiler (gfortran >= 10, Intel ifx/ifort or flang)

## Installation

Clone the repository and build with fpm:

```sh
git clone https://github.com/thomasbiekoetter/odeint.git
cd odeint
fpm build --profile='release'
```

Run the bundled tests to verify the build:

```sh
fpm test --profile='release'
```

Each test integrates exponential decay (dy/dt = -y, y(0) = 1) and prints the numerical solution, exact value, and absolute error at t = 5. For example, the `rk4` test gives:

```
RK4 integration: dy/dt = -y,  y(0) = 1
  steps     : 100
  t         :   5.00
  y(t)      :  6.73795E-03
  exact     :  6.73795E-03
  error     :  1.90522E-09
```

### Quadruple precision

By default the library uses double precision (`real64`). To build with quadruple precision, pass the preprocessor flag (not yet supported by flang):

```sh
fpm build --flag "-DQUAD" --profile='release'
fpm test  --flag "-DQUAD" --profile='release'
```

## Using odeint in another fpm project

Add odeint as a dependency in your project's `fpm.toml`:

```toml
[dependencies]
odeint = { git = "https://github.com/thomasbiekoetter/odeint.git", profile="release" }
```

Then `use` the integrator module, plus the working precision kind:
```fortran
use odeint__integrate, only : integrate   ! generic integrator
use odeint__config,    only : wp          ! working precision kind
```

Select the method with the optional `method` argument; `"euler"`, `"rk4"`, or
`"rkck"`. If omitted, it defaults to `"rkck"`:
```fortran
call integrate(dydt, y0, tstart, tend, nsteps, t, y, method="euler") ! Euler (fixed step)
call integrate(dydt, y0, tstart, tend, nsteps, t, y, method="rk4")   ! Runge-Kutta 4 (default, fixed step)
call integrate(dydt, y0, tstart, tend, nsteps, t, y, method="rkck")  ! adaptive Cash-Karp RK45
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
| `method` | in (optional) | `character(len=*)` | Integration method (`"euler"`, `"rk4"`, or `"rkck"`; default `"rk4"`) |
| `atol` | in (optional) | `real(wp)` | Absolute error tolerance (adaptive methods only; default `1e-9`) |
| `rtol` | in (optional) | `real(wp)` | Relative error tolerance (adaptive methods only; default `1e-6`) |

The time grid is uniformly spaced: `t(1) = tstart`, `t(nsteps) = tend`, step size `h = (tend - tstart) / (nsteps - 1)`. This holds for every method — the adaptive method also returns its solution on this uniform grid (see [Adaptive step size and error control](#adaptive-step-size-and-error-control)).

The solution array is indexed as `y(i, j)`, where `i` is the time index and `j` is the component index.

### Adaptive step size and error control

`"rkck"` is an adaptive Cash-Karp Runge-Kutta method — an embedded 4th/5th order pair. Rather than taking one fixed step per output interval, it substeps internally with automatic step-size control, shrinking the step where the solution changes quickly and growing it where the solution is smooth.

For adaptive methods, `nsteps` sets only the number of *output* samples returned: the solver still lands exactly on the uniform grid `t`, but the number of internal steps taken between output points is chosen automatically to meet the requested accuracy. (For the fixed-step methods `"euler"` and `"rk4"`, `nsteps` is the number of steps actually taken.)

Accuracy is set by the optional `atol` and `rtol` arguments. The per-component error is kept below `max(atol, rtol*|y|)`, so `atol` dominates near zero crossings and `rtol` dominates for large-magnitude components. Tighter tolerances produce smaller internal steps and more function evaluations:

```fortran
call integrate(  &
  dydt, y0, tstart, tend, nsteps, t, y, &
   method="rkck", rtol=1.0e-9_wp, atol=1.0e-12_wp)
```

Both tolerances are ignored by the fixed-step methods.

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
  use odeint__integrate, only : integrate
  implicit none

  real(wp), parameter :: y0(1) = [1.0e0_wp]
  real(wp), parameter :: tstart = 0.0e0_wp
  real(wp), parameter :: tend = 5.0e0_wp
  integer, parameter :: nsteps = 100
  real(wp), allocatable :: t(:)
  real(wp), allocatable :: y(:, :)

  call integrate(exp_decay, y0, tstart, tend, nsteps, t, y, method="rk4")

  write(*, '(A)')         "RK4 integration: dy/dt = -y,  y(0) = 1"
  write(*, '(A, I0)')     "  steps     : ", nsteps
  write(*, '(A, F6.2)')   "  t         : ", t(nsteps)
  write(*, '(A, ES12.5)') "  y(t)      : ", y(nsteps, 1)
  write(*, '(A, ES12.5)') "  exact     : ", exp(-t(nsteps))
  write(*, '(A, ES12.5)') "  error     : ", abs(y(nsteps, 1) - exp(-t(nsteps)))

contains

  ! Example: exponential decay
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

| Method | Module | Order | Step size | Recommended use |
|---|---|---|---|---|
| Euler | `odeint__euler` | 1 | fixed | Quick prototyping; requires many steps for accuracy |
| RK4 | `odeint__rk4` | 4 | fixed | General purpose; accurate with far fewer steps |
| Cash-Karp (RK45) | `odeint__rkck` | 5(4) | adaptive | Error-controlled; efficient when the solution has both fast and slow regions |

RK4 makes four function evaluations per step but converges as O(h⁴), so it typically needs far fewer steps than Euler to achieve the same accuracy. Prefer RK4 unless you have a specific reason to use Euler.

The adaptive Cash-Karp method (`"rkck"`) controls its own error instead of relying on a preset step count. It is the best choice when you need a target accuracy without hand-tuning `nsteps`, or when the solution has both fast and slow regions. It spends internal steps only where they are needed. Each accepted step costs six function evaluations, but it typically takes far fewer steps than a fixed-step method at the same accuracy. Set the accuracy with `atol`/`rtol` rather than with `nsteps` (see [Adaptive step size and error control](#adaptive-step-size-and-error-control)).

## License

This project is licensed under the [GNU Affero General Public License v3.0](LICENSE) (AGPL-3.0). Derivative works and projects that use this library over a network must also be released under the AGPL-3.0.
