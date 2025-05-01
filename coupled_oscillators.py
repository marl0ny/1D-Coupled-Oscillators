import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fft import dst

N = 100
MASS = 1.0
OMEGA = 1.0
OMEGA_DIAG = 2.0*np.sin(0.5*np.pi*(np.arange(N) + 1)/(N + 1))
TOTAL_STEPS = 200
TIME_SLICES = np.linspace(0.0, 40.0, TOTAL_STEPS)


def get_force_matrix(n):
    force_mat = np.zeros([n, n])
    for i in range(n):
        force_mat[i, i] = -2.0
        if i < n-1:
            force_mat[i, i+1] = 1.0
        if i > 0:
            force_mat[i, i-1] = 1.0
    return force_mat


# FORCE_MAT = np.array([(1.0 if abs(i - j) == 1 else -2.0) 
#                       if abs(i - j) < 2 else 0.0
#                       for i in range(N) for j in range(N)]).reshape(N, N)
FORCE_MAT = get_force_matrix(N)


def eom(xp, *_):
    x, p = xp[0: N], xp[N: 2*N]
    return np.hstack([p/MASS, MASS*OMEGA**2*FORCE_MAT @ x])


def eom_diag(xp, *_):
    x, p = xp[0: N], xp[N: 2*N]
    return np.hstack([p/MASS, -MASS*OMEGA_DIAG**2*x])


if __name__ == '__main__':
    index = np.arange(N)
    x = np.sin((index + 1)*np.pi*7/(N + 1))
    x = np.exp(-0.5*(index - N/2)**2/(N/30)**2)
    p = np.zeros([N])
    xp = np.hstack([x, p])
    xp_t = odeint(eom, xp, TIME_SLICES)

    xp_diag = np.hstack([dst(x, type=1, norm='ortho'), 
                        dst(p, type=1, norm='ortho')])
    xp_diag = odeint(eom_diag, xp_diag, TIME_SLICES)

    xt_dst = np.array([dst(xp_diag[i, 0:N], type=1, norm='ortho')
                    for i in range(TOTAL_STEPS)])

    plt.plot(xt_dst[0, 0:N])
    plt.plot(xt_dst[TOTAL_STEPS - 1, 0:N])
    plt.plot(xp_t[0, 0:N], linestyle='--')
    plt.plot(xp_t[TOTAL_STEPS - 1, 0:N], linestyle='--')

    # for i in range(TOTAL_STEPS):
    #     alpha = float(i)/TOTAL_STEPS
    #     plt.plot(xp[i, 0:N].T, color='blue', alpha=alpha)

    plt.show()
