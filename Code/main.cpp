// Compile this:
// c++ @flags.inp && ./main parameters.txt r0_data.txt p0_data.txt

/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Simulation of Argon Gas. That is a molecular dynamics      *
 * simulation for a gas interacting with van der Waals        *
 * forces. We observe a phase transition (solid -> gas)       *
 * in the simulation and study some thermodynamic properties. *
 * Author: Mateusz Kowalczyk                                  *
 * Date: 12.06.2022                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **/

#include "argon.h"
#include "stats.h"
#include <iostream>
#include <chrono>

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage: ./main <parameters_input> <r0_output> <p0_output>";
        exit(1);
    }

    double *pAbs;
    usint N;

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
    std::tie(pAbs, N) = A->getMomentumAbs();
    // Call function `simulateDynamics()` is optional.
    // But obviously it is the core of entertainment and playing with the system.
    A->simulateDynamics(argv[5], argv[6]);
    delete A;
    // --------------------------------------
    std::chrono::high_resolution_clock::time_point tk = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double(tk - tp);
    std::cout << "Execution time on CPU: " << ms_double.count() << " ms.";

    // for (int i = 0; i < N; i++)
    //     std::cout << *(pAbs + i) << '\n';

    Stats R(0., 30., 25);
    R.printStats(pAbs, N);

    return EXIT_SUCCESS;
}