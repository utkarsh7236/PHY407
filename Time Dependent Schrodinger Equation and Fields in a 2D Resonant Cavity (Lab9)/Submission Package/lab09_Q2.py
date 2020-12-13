"""
author: Aslesha
"""

import numpy as np
from numpy import sin, pi
from lab09_Q2_functions import dXYt2, idXYt2
import matplotlib.pyplot as plt
from tqdm import tqdm


# Set given conditions
tau = 0.01
T = 20
N = int(T / tau)  # T = N*tau
Lx = Ly = J0 = m = n = c = 1
P = 32

# imlement discretization
t = np.arange(0, N) * tau
ax = Lx / P
ay = Ly / P
Xp = np.arange(0, P) * ax
Yp = np.arange(0, P) * ay

# get the coordinates at each point
xv, yv = np.meshgrid(Xp, Yp)


def perform_fourier_evolution(Hx, Hy, Ez, Jz):
    """
    Take the fourier transform of Hx, Hy, Jz, and Ez. Evolve the Fourier 
    coefficients using eqns. (11) and reconstruct Hx,Hy, and Ez via a 2D 
    inverse Fourier transform. Return the reconstructed Hx, Hy and Ez.
    """
    # Take the Fourier transforms of Hx, Hy, Jz, and Ez (eqns. 10). 
    Hx_fourier = dXYt2(Hx, 'S', 'C')
    Hy_fourier = dXYt2(Hy, 'C', 'S')
    Jz_fourier = dXYt2(Jz, 'S', 'S')
    Ez_fourier = dXYt2(Ez, 'S', 'S')

    # Evolve the Fourier coefficients using eqns. (11)
    Dx = pi * c * tau / (2 * Lx)
    Dy = pi * c * tau / (2 * Ly)

    # p and q derived from xv and yv
    p = xv / ax
    q = yv / ay

    # eqn 11a
    E_numerator = (1 - p ** 2 * Dx ** 2 - q ** 2 * Dy ** 2) * Ez_fourier + \
                  2 * q * Dy * Hx_fourier + 2 * p * Dx * Hy_fourier + tau * Jz_fourier
    E_denominator = (1 + p ** 2 * Dx ** 2 + q ** 2 * Dy ** 2)
    E_new = E_numerator / E_denominator

    # eqn 11b
    X_new = Hx_fourier - q * Dy * (E_new + Ez_fourier)

    # eqn 11c
    Y_new = Hy_fourier - p * Dx * (E_new + Ez_fourier)

    # Reconstruct Hx,Hy, and Ez via a 2D inverse Fourier transform
    Hx_reconstructed = idXYt2(X_new, 'S', 'C')
    Hy_reconstructed = idXYt2(Y_new, 'C', 'S')
    Ez_reconstructed = idXYt2(E_new, 'S', 'S')

    return Hx_reconstructed, Hy_reconstructed, Ez_reconstructed


def main(omega_list, question='d'):
    """
    Main program to calculate the electromagnetic field and plot the traces.
    omega (List) : list of omega to be used
    question (str) : Question number. Either 'c', 'd' or 'e'.
    """

    # to store max Ez(x=0.5,y=0.5,t) for each omega
    Ez_max_arr = []

    for omega in tqdm(omega_list):
        # use eqn 8 to generate the current pattern J_z(x,y,t)
        Jz = J0 * sin((m * pi * xv) / Lx) * sin((n * pi * yv) / Ly) * sin(omega * t[0])

        # initialize Hx, Hy, Ez to 0
        Hx = np.zeros(xv.shape)
        Hy = np.zeros(xv.shape)
        Ez = np.zeros(xv.shape)

        # array to store Hx(x=0.5, y=0.0) for plotting
        Hx_plot_arr = []
        # find the index
        indx_Hx = np.argwhere((xv == 0.5) & (yv == 0.0))
        # array to store Hy(x=0.0, y=0.5) for plotting
        Hy_plot_arr = []
        # find the index
        indx_Hy = np.argwhere((xv == 0.0) & (yv == 0.5))
        # array to store Ez(x=0.5, y=0.5) for plotting
        Ez_plot_arr = []
        # find the index
        indx_Ez = np.argwhere((xv == 0.5) & (yv == 0.5))

        # loop over time
        for tn in t:
            # call perform_fourier_evolution to get updated Hx, Hy and Hz
            Hx, Hy, Ez = perform_fourier_evolution(Hx, Hy, Ez, Jz)

            # update J_z(x,y,t)
            Jz = J0 * sin((m * pi * xv) / Lx) * sin((n * pi * yv) / Ly) * sin(omega * tn)

            # append the results for plotting
            # use flipped indices coz of meshgrid conventons
            Hx_plot_arr.append(Hx[indx_Hx[0, 1], indx_Hx[0, 0]])
            Hy_plot_arr.append(Hy[indx_Hy[0, 1], indx_Hy[0, 0]])
            Ez_plot_arr.append(Ez[indx_Ez[0, 1], indx_Ez[0, 0]])

        # store the maximum amplitude of Ez(x=0.5,y=0.5,t) for omega
        Ez_max_arr.append(max(np.abs(Ez_plot_arr)))

    # plot the figure
    plt.figure(figsize=(8, 6))

    if question == 'c' or question == 'e':
        # plot Hx, Hy and Ez for c and e
        plt.plot(t, Hx_plot_arr, '-', label="$H_x(x=0.5, y=0)$")
        plt.plot(t, Hy_plot_arr, '--', label="$H_y(x=0, y=0.5)$")
        plt.plot(t, Ez_plot_arr, '-', label="$E_z(x=0.5, y=0.5)$")
        plt.xlabel("Time (s)", fontsize=14)
        plt.ylabel("$H_x, H_y, E_z$", fontsize=14)
        plt.grid()
        if question == 'c':
            plt.title("Traces of $H_x$, $H_y$, $E_z$ for $\omega=3.75$", fontsize=14)
        else:
            plt.title("Traces of $H_x$, $H_y$, $E_z$ for $\omega=\omega_0^{1,1}$", fontsize=14)
    else:
        # plot amplitude of max Ez for d
        normal_f = np.pi * c * np.sqrt((n * Lx) ** (-2) + (m * Ly) ** (-2))  # eqn 21
        plt.plot(omega_list, Ez_max_arr, "-*", label="$max(|E_z(x=0.5, y=0.5, t)|)$")
        # plot a vertical line at normal_f
        plt.axvline(x=normal_f, color='r', label="Normal frequency")
        plt.xlabel("$\omega [s^{-1}]$", fontsize=14)
        plt.ylabel("$max(|E_z|)$", fontsize=14)
        plt.title("Max amplitude of $E_z$ vs $\omega$", fontsize=14)

    # common plotting functions
    plt.legend()
    plt.savefig("q2{}.pdf".format(question))


# run the program
# Q2c
# main([3.75], 'c')

# Q2d
omega_arr = np.linspace(0, 9)
# main(omega_arr, 'd')

# Q2e
normal_omega = [np.pi * c * np.sqrt((n*Lx)**(-2) + (m*Ly)**(-2))]
main(normal_omega, 'e')
