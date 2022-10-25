# **Argon Molecular Dynamics**

![GitHub last commit](https://img.shields.io/github/last-commit/mateuszk098/Argon-Molecular-Dynamics)

**That is a simple molecular dynamics simulation for one type of atom (especially Argon) interacting with van der Waals' forces. One may observe the solid &rarr; gas phase transition and investigate thermodynamical properties that characterise the gas state.**

## **What is Molecular Dynamics?**
**Referring to Wikipedia:**

_**"Molecular dynamics (MD) is a computer simulation method for analyzing the physical movements of atoms and molecules. The atoms and molecules are allowed to interact for a fixed period of time, giving a view of the dynamic "evolution" of the system. In the most common version, the trajectories of atoms and molecules are determined by numerically solving Newton's equations of motion for a system of interacting particles, where forces between the particles and their potential energies are often calculated using interatomic potentials or molecular mechanics force fields."**_

**More about [Molecular Dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics#:~:text=Molecular%20dynamics%20(MD)%20is%20a,%22evolution%22%20of%20the%20system.).**


## **What this offers?**
- **Simulation of dynamics for an any ideal gas which satisfies the ideal gas law and the Maxwell-Boltzmann distribution.**
- **Simulation of melting the crystals.**
- **Calculations of fundamental properties of the system.**
- **Observations of simulations may be provided by an open-source viewer [Jmol](http://jmol.sourceforge.net/).**

**More about [Maxwell-Boltzmann Distribution](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution) and [Ideal Gas Law](https://en.wikipedia.org/wiki/Ideal_gas_law).**

## **Units used in the system:**
- **Length: 1 nm (10<sup>-9</sup> m).**
- **Time: 1 ps (10<sup>-12</sup> s).**
- **Mass: 1 u (1.66 x 10<sup>-27</sup> kg).**
- **Temperature: 1 K.**

## **How to use?**
**The primary program control is done by setting the system parameters.**
**Parameters to set in program:**

- **n - Number of atoms along the crystal edge (default 6).**
- **m - Atomic mass (default 40.0 - Argon).**
- **e - Minimum of the potential (default 1.0).**
- **R - Interatomic distance for which occurs minimum of the potential (default 0.38).**
- **k - Boltzmann constant (default 8.31e-3).**
- **f - Elastic coefficient of sphere which confines molecules (default 1e4).**
- **L - Radius of sphere which confines atoms (default 6.0).**
- **a - Interatomic distance (default 0.38).**
- **T0  - Initial temperature (default 1e4).**
- **tau - Simulation time step (default 1e-3).**
- **So - Initial number of steps for thermalization of the system (default 5000).**
- **Sd - Number of steps for mainly simulation (default 50000).**
- **Sout - Interval with which information about the system are saved (default 500).**
- **Sxyz - Interval with which positions of the molecules are saved (default 500).**
---

**C++ code to set in main file:**
```c++
// Usage: ./main <1> <2> <3> <4> <5> <6> <7>
// Where:
// <1> - input file with parameters in `Config` folder e.g. parameters.txt
// <2> - output file with initial positions to save in `Out` folder e.g. r0_init.txt
// <3> - output file with initial H, T, P to save in `Out` folder e.g. htp_init.txt
// <4> - output file with initial momenta to save in `Out` folder e.g. p0_init.txt
// <5> - output file with positions from the whole simulation to save in `Out` folder e.g. rt_sim.txt
// <6> - output file with H, T and P from the whole simulation to save in `Out` folder e.g. htp_sim.txt
// <7> - output file with initial momentum histogram to save in `Out` folder e.g. hist.txt

// Create object first.
Argon *A = new Argon;

// Call function `setParameters()` is optional.
// If you do not give file with own parameters, then simulation suppose default values.
A->setParameters(argv[1]);

// Call function `checkParameters()` is optional.
// It is only to information if system is properly set.
A->checkParameters();

// Call function `initialState()` is required if you want to get to simulation.
A->initialState(argv[2], argv[3], argv[4]);

// Get absolute values of momenta, its size and calculated temperature
// is required if you want to calculate statistics.
// You may call this function after `initialState()` or after `simulateDynamics()`.
usint N;
double *pAbs, T, k, m;
std::tie(pAbs, N, T, k, m) = A->getMomentumAbs();

// Call function `simulateDynamics()` is optional.
// But obviously it is the core of entertainment and playing with the system.
// That At the end of the simulation, the program checks if the ideal gas law is 
// fulfilled (It is if the value is around 1). Moreover, while the whole simulation,
// the total energy should be constant.
A->simulateDynamics(argv[5], argv[6]);

// Do not forget to release memory
delete A;

// Calculate statistics from Maxwell-Boltzmann distribution is optional.
// That provides calculation of most probable momentum, mean momentum,
// mean square momentum and kinetic energy. 
Stats *S = new Stats;
S->setInputFromArgon(pAbs, N, T, k, m);
S->evaluateHist(argv[7]);
delete S;
```

## **Example results:**

**Dynamics of 216 molecules at temperature 10<sup>4</sup> K:** | **Crystal of 15625 molecules at temperature 10<sup>2</sup> K:**
:-------------------------------------------------:|:-------------------------------------------------:
<img src="https://github.com/mateuszk098/Argon-Molecular-Dynamics/blob/master/Images/ArgonGasState.gif" width="478"/> | <img src="https://github.com/mateuszk098/Argon-Molecular-Dynamics/blob/master/Images/15625molecules100K.png" width="478"/>

<p align="center">
  **Momentum distribution fulfils Maxwell distribution:**
  <img src="https://github.com/mateuszk098/Argon-Molecular-Dynamics/blob/master/Images/15625molecules_hist.png" width="478"/>
</p>
