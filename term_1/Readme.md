I created generalized cauchy problem solver launched with json configuration file which describes:
- value type to use in calculations (float, double, long double)
- summation algorithm (naive, kahan),
- the problem (simple harmonic oscillator right now) with its initial value and parameters,
- method (euler, heun, runge_kutta),
- constraints to stop iterations (constraint on the amount of iterations, on a relative deviation of coordinates from analytical solution, on a relative deviation from invariant of the problem aka energy).
