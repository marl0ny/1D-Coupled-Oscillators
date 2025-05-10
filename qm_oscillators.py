"""
Simulating coupled 1D Quantum Harmonic Oscillators.

This is based on the following article:

Scott C. Johnson and Thomas D. Gutierrez, 
"Visualizing the phonon wave function",
American Journal of Physics 70, 227-237 (2002) 
https://doi.org/10.1119/1.1446858 

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.fft import dst


N = 100
MASS = 1.0
OMEGA = 1.0
OMEGA_DIAG = 2.0*np.sin(0.5*np.pi*(np.arange(N) + 1)/(N + 1))


def coherent_state(x, t, x0, p0, m, omega, hbar):
    """Formula for the harmonic oscillator coherent state from
    Shankar, pg. 610, equation 21.1.132 in exercise 21.2.18.
    """
    z_real = np.sqrt((m*omega)/(2.0*hbar))*x0
    z_imag = np.sqrt(1.0/(2.0*m*omega*hbar))*p0
    z = z_real + 1j*z_imag
    z = z*np.exp(-1j*omega*t)
    # z*z = z_real^2 - z_imag^2 + 2j*z_real*z_imag
    # a = m*omega/hbar
    # b = 2.0*np.sqrt(2.0*m*omega/hbar)*z_real
    # b2_a = (8.0*z_real*m*omega/hbar)/(m*omega/hbar)
    norm_factor = (
        (m*omega/(np.pi*hbar))**0.25)
        # *np.exp(-abs(z)**2/2.0))
    return norm_factor*(
        np.exp(-np.real(z)*z
               -m*omega/(2.0*hbar)*x**2
               +np.sqrt(2.0*m*omega/hbar)*z*x)
    )

def coherent_state_prod(x, t, x0, p0, m, omega, hbar):
    prod = 1.0
    for i in range(N):
        prod *= np.abs(coherent_state(x[i], t, x0[i], p0[i], 
                                      m, omega[i], hbar))
    return np.abs(prod)**2


index = np.arange(N)
initial_pos = 0.0*index # 5.0*np.exp(-0.5*(index - N/2)**2/(N/30)**2)
initial_diag_pos = dst(initial_pos, type=1, norm='ortho')


amplitude = coherent_state_prod(
    x=initial_diag_pos, t=0.0,
    x0=initial_diag_pos, p0=np.zeros([N]),
    m=1.0, omega=OMEGA_DIAG, hbar=1.0
)

def metropolis(x0, delta, dist_func, steps, kw):
    """
    See pg. 429 to 430 of "An Introduction to Computer Simulation Methods"
    by Harvey Gould, Jan Tobochnik, and Wolfgang Christian for a general 
    outline of the Metropolis algorithm.

    Gould H., Tobochnik J., Christian W., "Numerical and Monte Carlo Methods,"
    in <i>An Introduction to Computer Simulation Methods</i>,
    2016, ch 11., pg 406-444.
    """
    shape = tuple([steps] + list(x0.shape))
    vals = np.zeros(shape)
    x_curr = np.copy(x0)
    for i in range(steps):
        vals[i] = x_curr
        x_next = x_curr + delta*(2.0*np.random.rand(*shape[1:]) - 1.0)/2.0
        prob_curr = dist_func(x_curr, **kw)
        prob_next = dist_func(x_next, **kw)
        if (prob_next >= prob_curr or 
            np.random.rand() <= prob_next/prob_curr):
            x_curr = x_next
    return vals


print(amplitude)
print(1.0/np.sqrt(OMEGA_DIAG))


configs_diag = metropolis(initial_diag_pos, 1.0/np.sqrt(OMEGA_DIAG),
                          coherent_state_prod, 1000,
                          kw=dict(t=0.0, 
                                  x0=initial_diag_pos, p0=np.zeros([N]),
                                  m=1.0, omega=OMEGA_DIAG, hbar=1.0))
# configs_diag = configs_diag[500:,:]

plt.plot(configs_diag.T, color='black', alpha=0.1)
plt.plot(initial_diag_pos, color='blue', linestyle='--')
plt.plot(initial_diag_pos - 1.0/np.sqrt(OMEGA_DIAG),
         color='orange', linestyle='--')
plt.plot(initial_diag_pos + 1.0/np.sqrt(OMEGA_DIAG),
         color='orange', linestyle='--')
plt.show()
plt.close()

initial_pos = dst(initial_diag_pos, norm='ortho', type=1)
configs = np.array([
    dst(configs_diag[i, :], norm='ortho', type=1)
    for i in range(100)
])
plt.plot(initial_pos)
plt.plot(configs.T, color='black', alpha=0.1)
plt.show()
plt.close()


def show_1d_coherent_state_example():
    x = np.linspace(-10.0, 10.0, 1000)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    line, = ax.plot(x, 
        np.abs(coherent_state(x, 0.0, x0=5.0, p0=10.0, 
                              m=1.0, omega=OMEGA, hbar=1.0)
            + coherent_state(x, 0.0, x0=-2.0, p0=10.0, 
                              m=1.0, omega=OMEGA, hbar=1.0)))
    animation_data = {'t': 0.0}
    def animation_func(*_):
        animation_data['t'] += 0.01
        t = animation_data['t']
        psi = coherent_state(
            x, t, x0=5.0, p0=10.0, 
            m=1.0, omega=OMEGA, hbar=1.0)
        psi += coherent_state(
            x, t, x0=-2.0, p0=10.0, 
            m=1.0, omega=OMEGA, hbar=1.0)
        line.set_ydata(abs(psi))
        return (line, )
    animation_data['ani'] = animation.FuncAnimation(
        fig, animation_func, blit=True, interval=1000.0/60.0)
    plt.show()


show_1d_coherent_state_example()

