"""If the initial wave function in a harmonic oscillator is a Gaussian at t=0,
get it for all subsequent times by integrating this initial wave function 
with the harmonic oscillator propagator.

See equation 7.3.28 on pg 196 chapter 7 of Shankar for the formula for the
harmonic oscillator propagator. In order to evaluate the integral of the
initial Gaussian wave function with the propagator which also happens to take
the form of a Gaussian, equations A.2.4 - A.2.5 in page 660 of the appendix
of Shankar was consulted.
"""
from sympy import Symbol
from sympy import exp, tan, sqrt, I, sin, cos, pi, conjugate, solve, diff
from sympy import limit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re


def gauss_int(a, b):
    return exp(b**2/(4*a))*sqrt(pi/a)

# def print_symbolic_squeezed_wave_func():
#     x = Symbol('x')
#     sigma = Symbol('sigma')
#     x0 = Symbol('x0')
#     p = Symbol('p')
#     m = Symbol('m')
#     hbar = Symbol('hbar')
#     omega = Symbol('omega')
#     T = Symbol('T')
#     int_val = gauss_int(
#         1/(4*sigma)**2 - (I*m*omega/hbar)/(2*tan(omega*T)),
#         x0/(2*sigma**2) + I*p/hbar - (I*m*omega/hbar)*2*x/(2*sin(omega*T))
#     )
#     const = sqrt(m*omega/(2*pi*I*hbar*sin(omega*T)))\
#         *exp((I*m*omega/hbar)*x**2/(2*tan(omega*T))
#             - x0**2/(4*sigma**2))\
#         *(1.0/sqrt(sigma*sqrt(2*pi)))
#     expr = int_val*const
#     print(expr)
#     # print(limit(expr.subs(omega, 1), T, 0))
#     # sqrt(m)*sqrt(-I)*exp(-x**2/(16*sigma**2) + x*x0/(2*sigma**2) 
#     # - x0**2/(4*sigma**2) + I*p*x/hbar)/(sqrt(hbar)*sqrt(-I*m/hbar))
#     # print(limit(expr.subs(omega, 1), T, pi/2))
#     # 2*sqrt(2)*sqrt(m)*sqrt(-I)*exp(3*x0**2/(4*sigma**2) 
#     # - 4*I*m*x*x0/hbar + 4*I*p*x0/hbar - 4*m**2*sigma**2*x**2/hbar**2
#     #  + 8*m*p*sigma**2*x/hbar**2 
#     # - 4*p**2*sigma**2/hbar**2)/(sqrt(hbar)*sqrt(sigma**(-2)))
#     # print(limit(expr.subs(omega, 1), T, pi))
#     # sqrt(I)*sqrt(m)*exp(-x**2/(16*sigma**2) - x*x0/(2*sigma**2)
#     #  - x0**2/(4*sigma**2) - I*p*x/hbar)/(sqrt(hbar)*sqrt(-I*m/hbar))
#     # print(limit(expr.subs(omega, 1), T, -pi/2))
#     # 2*sqrt(2)*sqrt(I)*sqrt(m)*exp(3*x0**2/(4*sigma**2)
#     #  + 4*I*m*x*x0/hbar + 4*I*p*x0/hbar - 4*m**2*sigma**2*x**2/hbar**2
#     #  - 8*m*p*sigma**2*x/hbar**2 
#     # - 4*p**2*sigma**2/hbar**2)/(sqrt(hbar)*sqrt(sigma**(-2)))

def print_squeezed_state_relations():
    x = Symbol('x', real=True)
    sigma = Symbol('sigma', real=True)
    x0 = Symbol('x0', real=True)
    p0 = Symbol('p0', real=True)
    m = Symbol('m', real=True)
    hbar = Symbol('hbar', real=True)
    omega = Symbol('omega', real=True)
    T = Symbol('T', real=True)
    b = x0/(2*sigma**2) + I*p0/hbar - (I*m*omega/hbar)*2*x/(2*sin(omega*T))
    a = 1/(4*sigma**2) - (I*m*omega/hbar)/(2*tan(omega*T))
    e_factor = b**2/(4*a) \
        + (I*m*omega/hbar)*x**2/(2*tan(omega*T)) - x0**2/(4*sigma**2)
    e_factor_re_part = e_factor/2 + conjugate(e_factor/2)
    res = e_factor_re_part.expand().simplify()
    a = -diff(diff(res, x), x)/2
    b = diff(res, x).subs(x, 0)
    x0_t = (b/(2*a)).expand().simplify()
    sigma_t = 1/sqrt(4*a.expand().simplify())
    print('x(t) = ', x0_t)
    print('sigma(t) = ', sigma_t)
    e_factor_im_part = (e_factor/2 
                        - conjugate(e_factor/2)).expand().simplify()
    omega_t = Symbol('omega_t')
    e_factor_im_part = e_factor_im_part.subs(T*omega, omega_t)
    print('i phase(t=0) = ', limit(e_factor_im_part, omega_t, 0))
    print('i phase(t=pi/2) = ', limit(e_factor_im_part, omega_t, pi/2))
    print('i phase(t=pi) = ',limit(e_factor_im_part, omega_t, pi))
    print('i phase(t=3pi/2) = ',limit(e_factor_im_part, omega_t, 3*pi/2))
    c, s = Symbol('c'), Symbol('s')
    e_factor_im_part = e_factor_im_part.subs(sin(omega_t), s)
    e_factor_im_part = e_factor_im_part.subs(cos(omega_t), c)
    e_factor_im_part = e_factor_im_part.subs(tan(omega_t), s/c)
    print(re.sub('\*\*', '', str(e_factor_im_part.simplify())))
    # print('i phase(t) = ', e_factor_im_part.simplify())
    # e_factor_im_part_res_den = (hbar*(hbar**2*tan(T*omega)**2 
    #                                   + 4*m**2*omega**2*sigma**4
    #                                   )*sin(T*omega)**2*tan(T*omega))
    # e_factor_im_part_res_num = 1j*(
    #     hbar**2*m*omega*x**2*sin(T*omega)**2*tan(T*omega)**2/2 
    #     - hbar**2*m*omega*x*x0*sin(T*omega)*tan(T*omega)**3 
    #     + hbar**2*m*omega*x0**2*sin(T*omega)**2*tan(T*omega)**2/2 
    #     + hbar**2*p0*x0*sin(T*omega)**2*tan(T*omega)**3 
    #     + 2*m**3*omega**3*sigma**4*x**2*sin(T*omega)**2 
    #     - 2*m**3*omega**3*sigma**4*x**2*tan(T*omega)**2 
    #     + 4*m**2*omega**2*p0*sigma**4*x*sin(T*omega)*tan(T*omega)**2 
    #     - 2*m*omega*p0**2*sigma**4*sin(T*omega)**2*tan(T*omega)**2)
    # print((e_factor_im_part_res_num/(sin(T*omega)**2*tan(T*omega))).expand())
    # print(e_factor_im_part_res_den)
    # e_factor_im_part_res = 1j*(
    #     hbar**2*m*omega*x**2*sin(T*omega)**2*tan(T*omega)**2/2 
    #     - hbar**2*m*omega*x*x0*sin(T*omega)*tan(T*omega)**3 
    #     + hbar**2*m*omega*x0**2*sin(T*omega)**2*tan(T*omega)**2/2 
    #     + hbar**2*p0*x0*sin(T*omega)**2*tan(T*omega)**3 
    #     + 2*m**3*omega**3*sigma**4*x**2*sin(T*omega)**2 
    #     - 2*m**3*omega**3*sigma**4*x**2*tan(T*omega)**2 
    #     + 4*m**2*omega**2*p0*sigma**4*x*sin(T*omega)*tan(T*omega)**2 
    #     - 2*m*omega*p0**2*sigma**4*sin(T*omega)**2*tan(T*omega)**2)\
    # /(hbar*(hbar**2*tan(T*omega)**2 + 4*m**2*omega**2*sigma**4)
    #   *sin(T*omega)**2*tan(T*omega))

print_squeezed_state_relations()

def squeezed_state_phase_factor(x, t, x0, p0, m, omega, sigma, hbar):
    c = np.cos(omega*t)
    s = np.sin(omega*t)
    c2 = c*c
    c3 = c*c*c
    s2 = s*s
    hbar2 = hbar*hbar
    x2 = x*x
    x02 = x0*x0
    m2 = m*m
    m3 = m*m*m
    sigma4 = sigma*sigma*sigma*sigma
    p02 = p0*p0
    omega2, omega3 = omega**2, omega**3
    omega_t = omega*t
    phase = (
        2*c3*m3*omega3*sigma4*x2 
        + c*m*omega*(hbar2*s2*x2 + hbar2*s2*x02 
                     - 4*m2*omega2*sigma4*x2 + 8*m*omega*p0*s*sigma4*x 
                     - 4*p02*s2*sigma4)/2 
        + hbar2*s2*x0*(-m*omega*x + p0*s)
        )/(hbar*s*(4*c2*m2*omega2*sigma4 + hbar2*s2))
    eps = 1e-60
    if abs(omega_t % (2.0*pi)) < eps:
        return np.exp(1j*p0*x/hbar)*(1.0 - 1.0j)/np.sqrt(2.0)
    elif abs(omega_t % pi) < eps:
        return np.exp(-1j*p0*x/hbar)*(1.0 - 1.0j)/np.sqrt(2.0)
    return np.exp(1j*phase)*(1.0 - 1.0j)/np.sqrt(2.0)
    


def squeezed_state(x, t, x0, p0, m, omega, sigma, hbar):
    x_t = x0*np.cos(t*omega) + p0*np.sin(t*omega)/(m*omega)
    s_t = np.sqrt(
        hbar**2*np.sin(t*omega)**2 
        + 4*m**2*omega**2*sigma**4*np.cos(t*omega)**2)/(2.0*m*omega*sigma)
    phase_factor = squeezed_state_phase_factor(
        x, t, x0, p0, m, omega, sigma, hbar)
    return np.exp(-(x - x_t)**2/(4*s_t**2))/np.sqrt(s_t*np.sqrt(2*np.pi))\
        *phase_factor


ani = {'t': 0.0}
X = np.linspace(-10.0, 10.0, 1000)
Y0 = squeezed_state(
    X, t=0.0, x0=4.0, p0=0.0, m=1.0, omega=1.0, sigma=0.5, hbar=1.0)
plt.plot(X, np.abs(Y0))
plt.plot(X, np.exp(-0.25*(X - 4.0)**2/0.5**2)/np.sqrt(0.5*np.sqrt(2.0*np.pi))
         , linestyle='--')
plt.show()
plt.close()
# import sys; sys.exit()


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
re_line, = ax.plot(X, np.real(Y0))
im_line, = ax.plot(X, np.imag(Y0))
ab_line, = ax.plot(X, np.abs(Y0), color='black')
ax.set_ylim(-1.0, 1.0)

def animation_func(*_):
    ani['t'] += 0.01
    y = squeezed_state(
        X, t=ani['t'], x0=5.0, p0=0.0, m=1.0, omega=1.0, sigma=1.0, hbar=1.0)
    re_line.set_ydata(np.real(y))
    im_line.set_ydata(np.imag(y))
    ab_line.set_ydata(np.abs(y))
    return re_line, im_line, ab_line


ani['animation'] = animation.FuncAnimation(
    fig, animation_func, blit=True, interval=1000.0/60.0)
plt.show()
plt.close()



