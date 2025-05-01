#include <complex>
#include <vector>

#ifndef _HARMONIC_
#define _HARMONIC_


/* Energy eigenstates for the harmonic oscillator referenced from
Shankar, pg 195, eq 7.3.22.*/
std::complex<double> stationary_state(
    size_t n, double x, double t,
    double m, double omega, double hbar);

std::complex<double> stationary_states_combination(
    std::vector<std::complex<double>> coefficients,
    double x, double t,
    double m, double omega, double hbar
);

/* Formula for the harmonic oscillator coherent state from
Shankar, pg. 610, equation 21.1.132 in exercise 21.2.18. */
std::complex<double> coherent_state(
    double x, double t, 
    double x0, double p0,
    double m, double omega, double hbar);

std::complex<double> coherent_states_combination(
    double x, double t, 
    std::vector<std::complex<double>> coefficients,
    std::vector<double> x0, std::vector<double> p0,
    double m, double omega, double hbar);


#endif

