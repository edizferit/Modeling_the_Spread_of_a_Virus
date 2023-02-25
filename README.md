# Modeling the Spread of a Virus

A project for Boğaziçi University's **ME 303: Computer Applications in Mechanical Engineering** class.

The spread of coronavirus in a population of 10000 persons where at first 10 people are exposed is mathematically modeled with first order differential equations and simulated for 100 days.

* Euler, 4th Order Runge Kutta, and ode45 solvers are used and the results are compared with respect to accuracy and computation time.
* The effect of the change of each parameter is investigated.


### Parameters that affect spread of the virus:
1. Encounters per day (c)
2. Transmission probability per encounter (B)
3. Rate at which infected person becomes infectious per day (a)
3. Rate at which infected person becomes symptomatic per day (g)
4. Recovery rate (w)
5. Hospital beds (AB)

### Sub groups of population:
1. Susceptibles (S)
2. Exposed (E)
3. Infected (I)
4. Medically symptomatic (M)
5. Recovered (R)

### Governing differential equations:
* $ dS/dt = -cB(I/(S+E+I+M+R))S $
* $ dE/dt = cB(I/(S+E+I+M+R))S - aE$
* $ dI/dt = aE - gI $
* $ dM/dt = gI - wM $
* $ dR/dt = wM $

### Solution plot:

<p align="center">
  <img src="https://github.com/edizferit/Modeling_the_Spread_of_a_Virus/blob/main/figures/preview.jpg?raw=true" width="50%">
</p>


