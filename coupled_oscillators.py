"""Simple script intended to sanity check if the equations of motion of the
classical coupled 1D harmonic oscillator system in terms of their real 
positions agree with the equations derived in normal coordinates.

The leftmost and rightmost oscillators are clamped in place so that their
positions are always set to zero; with these boundary conditions the 
discrete sine transformation can be used to get the oscillator positions
in normal coordinates.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fft import dst


N = 100  # Number of oscillators that are free to move
MASS = 1.0  # Mass of each oscillator
OMEGA = 1.0
# Angular frequency of each of the normal modes
OMEGA_DIAG = 2.0*np.sin(0.5*np.pi*(np.arange(N) + 1)/(N + 1))
TOTAL_STEPS = 200  # Total number of time steps
TIME_SLICES = np.linspace(0.0, 40.0, TOTAL_STEPS)
# Matrix that expresses the couplings between each of the oscillators.
# It is used to compute the forces, where
# \frac{d\textbf{p}}{dt} = MASS\OMEGA^2COUPLINGS\textbf{x}
COUPLINGS = np.array([(1.0 if abs(i - j) == 1 else -2.0) 
                      if abs(i - j) < 2 else 0.0 
                      for i in range(N) for j in range(N)]).reshape(N, N)


def eom(xp, *_):
    """Equations of motion for the coupled harmonic oscillator system.
    """
    x, p = xp[0: N], xp[N: 2*N]
    return np.hstack([p/MASS, MASS*OMEGA**2*COUPLINGS @ x])


def eom_diag(xp, *_):
    """Equations of motions of the coupled harmonic oscillator system
    in normal coordinates. In these coordinates the oscillators are
    decoupled from each other.
    """
    x, p = xp[0: N], xp[N: 2*N]
    return np.hstack([p/MASS, -MASS*OMEGA_DIAG**2*x])


if __name__ == '__main__':
    index = np.arange(N)

    # Initial positions and momenta of the coupled oscillators
    x, p = np.exp(-0.5*(index - N/2)**2/(N/30)**2), np.zeros([N])
    xp = np.hstack([x, p])

    # Numerically find the positions and momenta
    # at a later time using Scipy's odeint.
    xp_t = odeint(eom, xp, TIME_SLICES)

    # Get the positions and momenta of the coupled oscillators
    # in normal coordinates, where in this coordinate system
    # the oscillators are decoupled from each other.
    xp_diag = np.hstack([dst(x, type=1, norm='ortho'), 
                        dst(p, type=1, norm='ortho')])
    # Numerically find the normal coordinate positions and momenta
    # at later times by using odeint.
    xp_diag_t = odeint(eom_diag, xp_diag, TIME_SLICES)
    # From the positions and momenta of the coupled oscillators
    # that were computed at later times in normal coordinates,
    # transform them back to the original coordinate system.
    xt_dst = np.array([dst(xp_diag_t[i, 0:N], type=1, norm='ortho')
                       for i in range(TOTAL_STEPS)])
    plt.xlim(0, N-1)
    plt.xlabel("Oscillator index")
    plt.ylabel("Amplitude")
    plt.plot(xt_dst[0, 0:N], 
             label="Oscillator positions at t=0")
    plt.plot(xt_dst[TOTAL_STEPS - 1, 0:N],
             label=f"Numerically computed"
             + f" positions at t={TIME_SLICES[TOTAL_STEPS-1]}")
    plt.plot(xp_t[0, 0:N], linestyle='--', 
             label="Oscillator positions at t=0")
    plt.plot(xp_t[TOTAL_STEPS - 1, 0:N], linestyle='--',
             label=f"Numerically computed transformed normal coordinates"
             + f" positions at t={TIME_SLICES[TOTAL_STEPS-1]}")
    plt.legend()

    plt.show()
