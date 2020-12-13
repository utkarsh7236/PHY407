import numpy as np
from SolveLinear import PartialPivot
import matplotlib.pyplot as plt


# constants for the program
R1 = R3 = R5 = 1000 # in ohm
R2 = R4 = R6 = 2000 # in ohm
C1 = 1*10**(-6) # in F
C2 = 0.5 * 10**(-6) # in F
x_plus = 3 # in V
omega = 1000 # rad/s


def run_experiment(R6_or_replacement, name):
    """
    Run experiment for RC and RLC circuits.
    R6_or_replacement : R6 value or its replacement L value
    name : name of the circuit
    """
    # matrix for the equations given in pg.3 
    A = np.array([[ 1/R1 + 1/R4 + 1j * omega * C1, - 1j * omega * C1, 0  ], 
         [ - 1j * omega * C1, 1/R2 + 1/R5 + 1j * omega * C1 + 1j * omega * C2, -1j * omega * C2], 
         [ 0, - 1j * omega * C2,  1/R3 + 1/R6_or_replacement + 1j * omega * C2 ]], 
                 dtype=np.complex)
    # vector of the system to be solved
    v = np.array([x_plus/R1, x_plus/R2, x_plus/R3], dtype=np.complex)
    
    # solve for x = [x1, x2, x3] using Partial Pivoting
    x = PartialPivot(A, v)
    
    # calculate V1, V2, V3 for t=0
    t = 0
    V1 = x[0] * np.exp(1j * omega * t) 
    V2 = x[1] * np.exp(1j * omega * t) 
    V3 = x[2] * np.exp(1j * omega * t) 
    
    # print amplitudes of voltages
    print("\t\t### For {} ###".format(name))
    print("The amplitudes of the three voltages at t=0 are given below:")
    print("\t|V1| = {}".format(abs(V1)))
    print("\t|V2| = {}".format(abs(V2)))
    print("\t|V3| = {}".format(abs(V3)))
    # print phases for x1, x2, x3
    print("Their phases in degrees at t=0 is:")
    print("\tphase of x1 = {}".format(np.angle(x[0], deg=True)))
    print("\tphase of x2 = {}".format(np.angle(x[1], deg=True)))
    print("\tphase of x3 = {}".format(np.angle(x[2], deg=True)))
    plt.savefig("Q1ci.pdf")
    
    
    # Get period
    T = 2 * np.pi / omega
    # time list for 3 periods 
    time = np.linspace(0, 2*T, 100)
    # get real part of all three voltages for time frame
    V1_real = (x[0] * np.exp(1j * omega * time)).real
    V2_real = (x[1] * np.exp(1j * omega * time)).real
    V3_real = (x[2] * np.exp(1j * omega * time)).real 
    
    
    # plot real voltage against time
    plt.figure()
    plt.plot(time, V1_real, label="V1 real")
    plt.plot(time, V2_real, label="V2 real")
    plt.plot(time, V3_real, label="V3 real")
    plt.legend()
    plt.grid()
    plt.xlabel("Time (s)")
    plt.ylabel("Real Voltage (V)")
    plt.title("Real Voltage for three periods of oscillations for {}".format(name))
    plt.savefig("Q1cii.pdf")


if __name__ == "__main__":
    # run experiment for RC circuit with R6 value
    run_experiment(R6, "RC circuit")
    # Inductor of 2H inductance
    L = R6/omega 
    # run experiment for RLC circuit with L value
    run_experiment(1j * omega * L, "RLC circuit")

