#include "harmonic.hpp"

using std::complex;
using std::vector;

#define PI 3.141592653589793


static size_t factorial(size_t n) {
    return (n == 0)? 1: n*factorial(n - 1);
}

/* Recursive Hermite polynomials.

See Shankar, pg. 195, 7.3.21 to obtain the base cases,
then 7.3.35 for the recursive relation itself.
*/
static double hermite(size_t n, double x) {
    if (n == 0)
        return 1.0;
    else if (n == 1)
        return 2.0*x;
    else
        return 2.0*x*hermite(n-1, x) - 2.0*(n-1)*hermite(n-2, x);
}

/* Energy eigenstates for the harmonic oscillator referenced from
Shankar, pg 195, eq 7.3.22.*/
complex<double> stationary_state(
    size_t n, double x, double t,
    double m, double omega, double hbar
) {
    complex<double> i (0.0, 1.0);
    double norm_factor = pow(
        m*omega/(PI*hbar*pow(2.0, 2.0*n)*pow(factorial(n), 2.0)), 0.25);
    return norm_factor*(
        exp(-0.5*m*omega*x*x/hbar - i*omega*(n + 0.5)*t)
        *hermite(n, x*sqrt(m*omega/hbar)));
}

std::complex<double> stationary_states_combination(
    std::vector<std::complex<double>> coefficients,
    double x, double t,
    double m, double omega, double hbar
) {
    std::complex<double> psi_t = 0.0;
    for (size_t n = 0; n < coefficients.size(); n++)
        psi_t += coefficients[n]*stationary_state(n, x, t, m, omega, hbar);
    return psi_t;
}

/* Formula for the harmonic oscillator coherent state from
Shankar, pg. 610, equation 21.1.132 in exercise 21.2.18. */
complex<double> coherent_state(
    double x, double t, 
    double x0, double p0,
    double m, double omega, double hbar
) {
    complex<double> i (0.0, 1.0);
    complex<double> z = 
        sqrt((m*omega)/(2.0*hbar))*x0
        + i*sqrt(1.0/(2.0*m*omega*hbar))*p0; 
    z *= exp(-i*omega*t);
    return std::pow(m*omega/(PI*hbar), 0.25)*(
        std::exp(-real(z)*z
            - m*omega/(2.0*hbar)*x*x 
            + sqrt(2.0*m*omega/hbar)*z*x));
}
