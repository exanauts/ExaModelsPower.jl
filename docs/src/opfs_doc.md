# Optimal Power Flow Formulation

The AC optimal power flow (ACOPF) formulation is a problem characterizing the nonlinear optimization problem which independent service operators must solve to route power through transmission systems. The OPF has also been built upon, leading to even more complex models that account for other considerations such as multiperiodicity and reliability. This documentation provides a brief outline of the key features of each version of the OPF modeled in ExaModelsPower.jl.

## Static OPF
The static OPF is the most basic version of the OPF, implemented over just a single time period. While some simplifications can be made to the model to make it easier to solve, the version implemented in ExaModelsPower.jl is the full, unrelaxed version. The focal constraint in the static OPF is the power balance - this constraint requires that the power generated at each bus is matched by the power consumed, shunt losses, and net power flow from the bus to others via transmission lines. 

A central complexity of the ACOPF is calculating the power transfer along transmission lines, which requires the following nonlinear constraint in polar coordinates, where G and B are parameters of the line, n and m are bus indeces, and V and θ represent voltage and phase angle. 

$$$
P_n = \sum_{mk} V_n V_m \left( G_{nmk} \cos \theta_{nm} + B_{nmk} \sin \theta_{nm} \right)
$$$

Due to the computational challenges introduced by these equations, ExaModelsPower.jl implements the OPF in both polar and rectangular coordinates to provide flexibility in the case that a certain algorithm may be more capable of performing computations for certain forms of nonlinear equations. 

A more complete description of the static OPF can be found [here](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=16f56c56e83bb2e2fa98a1400dcf9367f04783b6).

## Multi-period OPF
The multi-period OPF (MPOPF) is a simple variant of the static OPF, but it can introduce powerful results. The MPOPF expands variables across time periods, requiring the same balance as in the static OPF to be balanced for each time period. This allows for the evaluation of real-world problems that involve fluctuating demand over time. The MPOPF also introduces an additional constraint to restrict how quickly generators can ramp up or down. 
$$$
-\Delta P \leq P_{g,t} - P_{g,t-1} \leq \Delta P
$$$

## MPOPF with storage
The MPOPF also enables the introduction of storage devices to the model. We model storage devices as discussed in __[Geth, Coffrin, Fobes 2020](https://arxiv.org/pdf/2004.14768)__ [2]. Transfer of power to and from storage devices along with associated losses are now included in the overall power balance for each bus. 

One challenge of implementing storage is the complementarity constraint that prevents storage devices from charging and discharging simultaneously. In order to not introduce integer variables, the constraint is implemented as so:
$$$
P_{d} \cdot P_{c} = 0
$$$

This formulation can lead to errors arising from numerical imprecision. As a result, ExaModelsPower.jl provides several options to avoid the explicit implementation of this complementarity constraint. 