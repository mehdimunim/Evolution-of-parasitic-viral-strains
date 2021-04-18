import os
import numpy as np
import matplotlib.pyplot as plt

# Simulations for the evolution of parasitic viral strains
# Superinfection
# Author: Mehdi Munim

# System to be visualized numerically
# Modeling the superinfection model proposed in "The Equations of Life" (M. Nowak)
# S is the host susceptible population
# I corresponds the strain-infected populations for different strains
# R is the "removed" individuals (dead population)

# Time parameters
step = 1/24
span = 500  # days
t = np.arange(1, 1 + span, step)


class Virus:
    def __init__(self, a, b, s):
        # virulences of the strains
        # ordered by increasing values
        self.a = a.copy()

        # transmissibilities of the strains
        self.b = b.copy()

        # take-over factor
        self.s = s

    def get_virulence(self):
        return self.a.copy()

    def get_transmissibility(self):
        return self.b.copy()

    def get_factor(self):
        return self.s

    def get_strain_number(self):
        return len(self.a)


def R0(transmissibility, virulence, birth, death):
    return transmissibility*birth/((death + virulence)*death)


class Host:
    def __init__(self, population, birth, death, initial_infected):

        # Initial population
        self.population = population

        # Birth rate
        self.birth = birth

        # Death rate
        self.death = death

        # Initial number of infected for all viral strains
        self.initial_infected = initial_infected.copy()

    def get_birth(self):
        return self.birth

    def get_death(self):
        return self.death

    def get_population(self):
        return self.population

    def get_initial_infected(self):
        return self.initial_infected.copy()


def fS(host, virus, S, I_list):
    """
    Growth rate of the susceptible population
    """
    a = virus.get_virulence()
    b = virus.get_transmissibility()

    birth = host.get_birth()
    death = host.get_death()

    result = birth - death*S

    result -= S*np.dot(b, I_list)
    return result


def fI(host, virus, S, I_list, strain):
    """
    Growth rate of the strain-infected population
    """

    a = virus.get_virulence()
    b = virus.get_transmissibility()
    s = virus.get_factor()

    birth = host.get_birth()
    death = host.get_death()

    I_list_before = I_list[:strain]
    I_list_after = I_list[strain+1:]

    result = b[strain]*S - birth - a[strain]

    result += s*b[strain]*np.sum(I_list_before)

    result -= s*np.dot(b[strain+1:], I_list_after)

    return I_list[strain]*result


def fR(virus, S, I_list):
    """
    Growth rate of the dead population
    """
    result = 0

    a = virus.get_virulence()

    result += np.dot(a, I_list)

    return S*result


def euler(host, virus, t):
    """
    Apply the euler method for differential equations to find S, I1, I2 and R of our systems
    for a defined time array and given initial conditions
    """
    # number of points
    n = len(t)

    # time step
    h = t[1] - t[0]

    # strain number
    total_strain = virus.get_strain_number()

    # Initializing arrays
    S = np.zeros(n)
    I_array = np.zeros((n, total_strain))
    R = np.zeros(n)

    # Giving initial values
    S[0] = 1  # the susceptible populations are normalized to 1
    I_array[0, :] = host.get_initial_infected()
    # R[0] = 0

    for i in range(1, n):
        # Getting previous values
        Sin = S[i-1]
        I_list_in = I_array[i-1, :].copy()
        Rin = R[i-1]

        # Getting the derivative of each function
        Sout = fS(host, virus, Sin, I_list_in)
        I_list_out = np.array([fI(host, virus, Sin, I_list_in, strain)
                               for strain in range(total_strain)])
        Rout = fR(virus, Sin, I_list_in)

        # Getting the next value
        S[i] = Sin + h*Sout
        I_array[i, :] = I_list_in + h*I_list_out
        R[i] = Rin + h*Rout

    return S, I_array, R


def main():
    # Number of points
    n = len(t)

    # Defining host and virus parameters
    birth = 1
    death = 1
    total_strain = 100
    s = 1
    N = 10**(3)
    # a few infected at the beginning
    initial_infected = np.array([2/N]*total_strain)

    # random virulence between 0 and 5
    a = 5*np.random.random_sample((total_strain, ))
    a = np.sort(a)
    b = np.array([8*virulence/(1+virulence) for virulence in a])

    host = Host(N, birth, death, initial_infected)
    virus = Virus(a, b, s)

    # applying Euler method
    (S, I_array, R) = euler(host, virus, t)

    # plotting abundance in function of virulence
    # proportions of each strain at the end
    sum_abundance = np.sum(I_array[n-1, :])
    abundance = I_array[n - 1, :]/sum_abundance

    fig, ax1 = plt.subplots()
    ax1.set_ylabel("abundance")
    ax1.set_xlabel("virulence")
    cmap = plt.get_cmap("gist_rainbow")
    def rescale(y): return (y - np.min(y)) / (np.max(y) - np.min(y))

    ax1.bar(a, abundance, width=0.02, color=cmap(rescale(b)))

    # Plotting R0
    ax2 = ax1.twinx()
    R0s = [R0(transmissibility, virulence, birth, death)
           for transmissibility, virulence in zip(b, a)]
    ax2.plot(a, R0s, c="black", linewidth=2)

    plt.title("s= {}".format(s))

    # saving and showing figure
    plt.savefig("superinfection with {} strains, {:e} points and s={}.png".format(
        virus.get_strain_number(), n, virus.get_factor()))
    plt.show()

    # Plotting R0 in function of virulence


main()
exit(0)
