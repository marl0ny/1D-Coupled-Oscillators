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
    // 1/(4*sigma^2) = (m*omega)/(2*hbar)
    // 1/sigma^2 = 2*m*omega/hbar
    // sigma^2 = hbar/(2*m*omega)
    // sigma = sqrt(hbar/(2*m*omega))
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

static std::complex<double> squeezed_state_phase_factor(
    double x, double t, 
    double x0, double p0, double sigma0,
    double m, double omega, double hbar
) {
    double c = cos(omega*t);
    double s = sin(omega*t);
    double c2 = c*c;
    double c3 = c*c*c;
    double s2 = s*s;
    double hbar2 = hbar*hbar;
    double x2 = x*x;
    double x02 = x0*x0;
    double m2 = m*m;
    double m3 = m*m*m;
    double sigma4 = sigma0*sigma0*sigma0*sigma0;
    double p02 = p0*p0;
    double omega2 = omega*omega;
    double omega3 = omega*omega*omega;
    double omega_t = omega*t;
    double phase = (
        2.0*c3*m3*omega3*sigma4*x2 
        + c*m*omega*(hbar2*s2*x2 + hbar2*s2*x02 
                     - 4.0*m2*omega2*sigma4*x2 + 8.0*m*omega*p0*s*sigma4*x 
                     - 4.0*p02*s2*sigma4)/2.0
        + hbar2*s2*x0*(-m*omega*x + p0*s)
        )/(hbar*s*(4.0*c2*m2*omega2*sigma4 + hbar2*s2));
    double eps = 1e-60;
    if (abs(remainder(omega_t, (2.0*PI))) < eps)
        return exp(std::complex<double>(0.0, p0*x/hbar))
            *std::complex<double>(1.0, -1.0)/sqrt(2.0);
    else if (abs(remainder(omega_t, (PI))) < eps)
        return exp(std::complex<double>(0.0, -p0*x/hbar))
            *std::complex<double>(1.0, -1.0)/sqrt(2.0);
    return exp(std::complex<double>(0.0, phase))
        *std::complex<double>(1.0, -1.0)/sqrt(2.0);
}

/*
If the initial wave function in a harmonic oscillator is a Gaussian at t=0,
get it for all subsequent times by integrating this initial wave function 
with the harmonic oscillator propagator.

See equation 7.3.28 on pg 196 chapter 7 of Shankar for the formula for the
harmonic oscillator propagator. In order to evaluate the integral of the
initial Gaussian wave function with the propagator which also happens to take
the form of a Gaussian, equations A.2.4 - A.2.5 in page 660 of the appendix
of Shankar was consulted.
*/
std::complex<double> squeezed_state(
    double x, double t, 
    double x0, double p0, double sigma0,
    double m, double omega, double hbar, bool omit_phase
) {
    double c = cos(t*omega), s = sin(t*omega);
    double hbar2 = hbar*hbar;
    double m2 = m*m;
    double omega2 = omega*omega;
    double sigma4 = sigma0*sigma0*sigma0*sigma0;
    double x_t = x0*c + p0*s/(m*omega);
    double s_t = sqrt(hbar2*s*s + 4.0*m2*omega2*sigma4*c*c)
        /(2.0*m*omega*sigma0);
    std::complex<double> phase_factor = (omit_phase)?
        std::complex<double>(1.0, 0.0):
        squeezed_state_phase_factor(
            x, t, x0, p0, sigma0, m, omega, hbar);
    return exp(-0.25*pow((x - x_t)/s_t, 2.0))/sqrt(s_t*sqrt(2.0*PI))
        *phase_factor;
}
