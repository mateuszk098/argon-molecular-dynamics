/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Simulation of Argon Gas. That is a molecular dynamics      *
 * simulation for a gas interacting with van der Waals        *
 * forces. We observe a phase transition (solid -> gas)       *
 * in the simulation and study some thermodynamic properties. *
 * Author: Mateusz Kowalczyk                                  *
 * Date: 12.06.2022                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **/

// Compile this: c++ @flags.inp
// Run this: ./main parameters.txt r0_init.txt p0_init.txt htp_init.txt rt_sim.txt htp_sim.txt hist.txt
// Or compile and run:
// c++ @flags.inp && ./main parameters.txt r0_init.txt p0_init.txt htp_init.txt rt_sim.txt htp_sim.txt hist.txt

#include "argon.h"
#include "stats.h"
#include <iostream>
#include <chrono>

int main(int argc, char *argv[])
{
    if (argc < 8)
    {
        std::cerr << "Usage: ./main <1> <2> <3> <4> <5> <6> <7>\n";
        std::cerr << "Where:\n";
        std::cerr << "<1> - input file with parameters in `Config` folder e.g. parameters.txt\n";
        std::cerr << "<2> - output file with initial positions to save in `Out` folder e.g. r0_init.txt\n";
        std::cerr << "<3> - output file with initial H, T, P to save in `Out` folder e.g. htp_init.txt\n";
        std::cerr << "<4> - output file with initial momenta to save in `Out` folder e.g. p0_init.txt\n";
        std::cerr << "<5> - output file with positions from the whole simulation to save in `Out` folder e.g. rt_sim.txt\n";
        std::cerr << "<6> - output file with H, T and P from the whole simulation to save in `Out` folder e.g. htp_sim.txt\n";
        std::cerr << "<7> - output file with initial momentum histogram to save in `Out` folder e.g. hist.txt\n";
        exit(1);
    }

    std::chrono::high_resolution_clock::time_point tp = std::chrono::high_resolution_clock::now();
    // --------------------------------------
    Argon *A = new Argon;

    // Call function `setParameters()` is optional.
    // If you do not give file with own parameters, then simulation suppose default values.
    A->setParameters(argv[1]);

    // Call function `checkParameters()` is otpional.
    // It is only to information if system is properly set.
    A->checkParameters();

    // Call function `initialState()` is required if you want to get to simulation.
    A->initialState(argv[2], argv[3], argv[4]);

    // Get absolute values of momenta, its size and calculated temperature
    // is required if you want to calculate statistics.
    usint N;
    double *pAbs, T, k, m;
    std::tie(pAbs, N, T, k, m) = A->getMomentumAbs();

    // Call function `simulateDynamics()` is optional.
    // But obviously it is the core of entertainment and playing with the system.
    A->simulateDynamics(argv[5], argv[6]);

    // Do not forget to release memory
    delete A;
    // --------------------------------------
    std::chrono::high_resolution_clock::time_point tk = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double(tk - tp);
    std::cout << "`main()` >: Argon execution time on CPU: " << ms_double.count() << " ms.\n";

    // Calculate statistics from Maxwell-Boltzmann distribution
    Stats *S = new Stats;
    S->setInputFromArgon(pAbs, N, T, k, m);
    S->evaluateHist(argv[7]);

    return EXIT_SUCCESS;
}