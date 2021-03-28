import os
import numpy as np
import matplotlib.pyplot as plt

# Simulations for the evolution of parasitic viral strains
# Author: Mehdi Munim
#

# System to be visualized numerically
# Inspired from "The Equations of Life" (M. Nowak)
# S is the host susceptible population
# I1, I2 are the strain 1 or strain 2 infected
# R is the "removed" individuals (dead population)
# dS/ dt = a - (b1*I1 + b2*I2)*S - b*S
# dI1/dt  = b1*I1*S - b*I1 - a1*I1
# dI2/dt  = b2*I2*S - b*I2 - a2*I2
# dR/dt = a1*I1 + a2*I2

# Global variables

# Time
step = 1/12  # one month
span = 10  # years
t = np.arange(1, span, step)

# Virus virulence of strains 1 and 2
a1 = 4.1
a2 = a1/4


# Virus transmissiblity of strains 1 and 2
b1 = 0.00001
b2 = 0.00001

# Host recovery rates of strains 1 and 2
v1 = 0.0001
v2 = 0.0001

# Host birth and death rates
a = 1
b = 0.00001

# Initial population infected by 1 and 2
I01 = 100
I02 = 100

# Host population
N = 100000
print(b1*N)

# Take-over rate in the model of Levin and Pimentel (1981)
s = 0.0001

# R0


def R0(beta, alpha):
    return beta*a/((b + alpha)*b)

# Defining functional second member


def fS(I1, I2, S):
    return a - b*S - S*(b1*I1 + b2*I2)


def fI1(I1, S):
    return I1*(b1*S - b - a1)


def fI2(I2, S):
    return I2*(b2*S - b - a2)


def fR(I1, I2):
    return a1*I1 + a2*I2

# Using euler function


def euler(t, N, I01, I02):
    """
    Apply the euler method for differential equations to find S, I1, I2 and R of our systems
    for a defined time array and given initial conditions
    """
    # number of points
    n = len(t)

    # time step
    h = t[1] - t[0]

    # Initializing arrays
    S = np.zeros(n)
    I1 = np.zeros(n)
    I2 = np.zeros(n)
    R = np.zeros(n)

    # Giving initial values
    S[0] = N
    I1[0] = I01
    I2[0] = I02
    # R[0] = 0

    for i in range(1, n):
        # Getting previous values
        Sin = S[i-1]
        I1in = I1[i-1]
        I2in = I2[i-1]
        Rin = R[i-1]

        # Getting the derivative of each function
        Sout = fS(I1in, I2in, Sin)
        I1out = fI1(I1in, Sin)
        I2out = fI2(I2in, Sin)
        Rout = fR(I1in, I2in)

        # Getting the next value
        S[i] = Sin + h*Sout
        I1[i] = I1in + h*I1out
        I2[i] = I2in + h*I2out
        R[i] = Rin + h*Rout

    return S, I1, I2, R


def main():
    (S, I1, I2, R) = euler(t, N, I01, I02)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    ax1.plot(t, S, c="blue")
    ax1.set_title("Susceptible population", fontsize=10)

    ax3.plot(t, R, c="red")
    ax3.set_title("Deceased caused by the virus", fontsize=10)

    ax2.plot(t, I1, c="green")
    ax2.set_title(
        "Infected by strain 1 (R0 = {:.2f})".format(R0(b1, a1)), fontsize=10)

    ax4.plot(t, I2, c="yellow")
    ax4.set_title(
        "Infected by strain 2 (R0 = {:.2f})".format(R0(b2, a2)), fontsize=10)

    fig.tight_layout()
    plt.savefig("a1={:.2e} a2={:.2e} beta={:.2e}.png".format(
        a1, a2, b1), dpi=500)

    plt.show()


main()
exit(0)
