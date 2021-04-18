# Evolution-of-parasitic-viral-strains
Simulations done for the English Tiny Thesis.
Mainly inspired from *The Equations of Life* by Martin Nowak.

## Presentation

My tiny thesis was aiming at studying the evolution of a virus, and its impact on populations. 

As a context, I assumed that a range of strains infects a population at t = 0.  I looked at which strain will survive in the end. With two strains, I also looked at the impact on populations. 

## Basic notions

* Virulence: how the virus is dangerous, often corresponds to the death rate.
* Transmissibility: how the virus is contagious, the rate at which susceptible became infected.
* R0: the expected cases provoked by patient 0 in a whole new susceptible population 

## Code

- simulation.py: first approach with two strains. I looked at which one will dominate, i.e died last. I also plotted the infected, susceptible, deceased evolutions against time.
- superinfection.py: considered a range of strains and looked which ones survive after a long time. 

## Outputs

Outputs can be found in outputs folder.



