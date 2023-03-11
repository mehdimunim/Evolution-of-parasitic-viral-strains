# Evolution of parasitic viral strains

## Presentation

> "Tiny Thesis" project for my engineering school. Loosely based on chapter "Evolution of Virulence" in "Evolutionary Dynamics" (Martin Nowak).

This work aims at modeling the evolution of a virus. We assumed that the virus exists in several strains with different virulences. At t0 , a range of strains infects a population. Through both mathematics and simulations, we look at which strains remain in the end. A strong emphasis is put in the impact on host populations.  

## Basic notions

* Virulence: how the virus is dangerous, often corresponds to the death rate.

* Transmissibility: how the virus is contagious, the rate at which susceptible became infected.

* R_0: the expected cases provoked by patient 0 in a whole new susceptible population. 

* SIR: Susceptible-Infected-Removed.
 
## Code

The code works in Windows with python 3.8. Requirements: matplotlib 3.2, numpy 1.18.

- simulation.py: A first approach with two strains. I looked at which one will dominate, i.e died last. I also plotted the infected, susceptible, deceased evolutions against time.

- superinfection.py: I considered a range of strains and looked which ones survive after a long time. 

## Document

My document explaining my reflections can be found [here](documents/Mehdi%20MUNIM%20Tiny%20Thesis.pdf). My presentation in PPTX is [there](documents/MUNIM.pptx).

## Source

[1] R.M Anderson and R.M. May. “Coevolution of hosts and parasites”. In: Parasitology (1982).

[2] B. R. Levin Antia R. and R.M. May. “Within-host population dynamics and the evolution and maintenance of microparasite virulence”. In: The American Naturalist (1994).

[3] Khan MN Billah MA Miah MM. “Reproductive number of coronavirus: A systematic
review and meta-analysis based on global level evidence”. In: Plos ONE (2020).

[4] F.N. Ratcliffe Frank Fenner. Myxomatosis. Cambridge University Press, 1965.

[5] W. O. Kermack and A. G. McKendrick. “A contribution to the mathematical theory of epidemics”. In: Proceedings of the Royal Society of London (1927).

[6] Peter J. Kerr. “Myxomatosis in Australia and Europe: A model for emerging infectious diseases”. In: Ecosystem Sciences (2012).

[7] Pimentel D. Levin S. “Selection of intermediate rates of increase in parasite-host systems”. In: American Naturalist (1981).

[8] Martin A. Nowak. Evolutionary dynamics: exploring the equations of life. Harvard University Press, 2006

## Contact

Created by Mehdi MUNIM, mehdi.munim@gmail.com.

## License

MIT License 

Copyright c 2021 - Mehdi MUNIM.
