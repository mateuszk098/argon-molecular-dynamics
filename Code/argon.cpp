#define _USE_MATH_DEFINES
#include "argon.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <iomanip>

/**************************************************************************************
 * Default constructor initializes exmaple parameters and memory to store informations
 * about the system. It also prints appropriate messages.
 * @return Nothing to return.
 *************************************************************************************/
Argon::Argon() noexcept : n(6), So(5000), Sd(50000), Sout(500), Sxyz(500), m(40.), e(1.),
                          R(0.38), k(8.31e-3), f(1e4), L(6.), a(0.38), T0(1e4), tau(1e-3),
                          initialStateCheck(false), mt(std::mt19937(time(nullptr)))
{
    N = n * n * n; // System is defined as 3D
    K = 3;

    std::cout << "`Argon()` :> Initialized parameters to default values." << '\n';
    std::cout << "`Argon()` :> Set pseudo-random number generator std::mt19937." << '\n';

    // Allocate memory and immediately set the values
    b0 = new double[K]{a, 0., 0.};
    b1 = new double[K]{a * 0.5, a * sqrt(3.) * 0.5, 0.};
    b2 = new double[K]{a * 0.5, a * sqrt(3.) / 6., a * sqrt(6.) / 3.};

    // Note the parenthesis at the end of new. These caused the allocated memory's
    // value to be set to zero (value-initialize)
    p = new double[K]();
    pAbs = new double[N]();
    Vs = new double[N]();

    r0 = new double *[N]();
    p0 = new double *[N]();
    Vp = new double *[N]();
    Fs = new double *[N]();
    Fi = new double *[N]();

    Fp = new double **[N]();

    for (usint i = 0; i < N; i++)
    {
        r0[i] = new double[K]();
        p0[i] = new double[K]();
        Vp[i] = new double[N]();
        Fs[i] = new double[K]();
        Fi[i] = new double[K]();

        Fp[i] = new double *[N]();

        for (usint j = 0; j < N; j++)
            Fp[i][j] = new double[K]();
    }

    std::cout << "`Argon()` :> Allocated memory for buffer.\n\n";
}

/**************************************************************************************
 * Destructor frees up buffer memory and print appropriate message about it.
 * @return Nothing to return.
 **************************************************************************************/
Argon::~Argon() noexcept
{
    delete[] b0;
    delete[] b1;
    delete[] b2;

    delete[] p;
    delete[] pAbs;
    delete[] Vs;

    for (usint i = 0; i < N; i++)
    {
        delete[] r0[i];
        delete[] p0[i];
        delete[] Vp[i];
        delete[] Fs[i];
        delete[] Fi[i];

        for (usint j = 0; j < N; j++)
            delete[] Fp[i][j];

        delete[] Fp[i];
    }

    delete[] r0;
    delete[] p0;
    delete[] Vp;
    delete[] Fs;
    delete[] Fi;

    delete[] Fp;

    std::cout << "`~Argon()` :> Memory released.\n\n";
}

/**************************************************************************************
 * This function reads parameters from input file and sets appropriate variables.
 * Then it reallocates required memory for buffer. Moreover it provides exception
 * handling for invalid parameters and files. If something has gone wrong, then
 * the function set default parameters and print appropriate message.
 * @param char* filename with parameters to set.
 * @return Nothing to return.
 *************************************************************************************/
void Argon::setParameters(const char *filename)
{
    std::ifstream input;
    std::string tmp;

    try
    {
        input.open("../Config/" + std::string(filename), std::ios::in);

        if (input.fail())
            throw std::ifstream::failure("Exception opening/reading input file");

        if (fileIsEmpty(input))
            throw std::ifstream::failure("Exception input file is empty");

        input >> tmp >> n >> tmp >> m >> tmp >> e >> tmp >> R >> tmp >> k >> tmp >> f >> tmp >> L >> tmp >> a;
        input >> tmp >> T0 >> tmp >> tau >> tmp >> So >> tmp >> Sd >> tmp >> Sout >> tmp >> Sxyz;

        if (n < 1 || n > 25)
            throw std::invalid_argument("Invalid argument: n. Must be between 1 and 25.");
        if (m < 0.)
            throw std::invalid_argument("Invalid argument: m. Must be positive.");
        if (e < 0.)
            throw std::invalid_argument("Invalid argument: e. Must be positive.");
        if (R < 0.)
            throw std::invalid_argument("Invalid argument: R. Must be positive.");
        if (k < 0. || k > 1.)
            throw std::invalid_argument("Invalid argument: k. Must be between 0 and 1.");
        if (f < 0.)
            throw std::invalid_argument("Invalid argument: f. Must be positive.");
        if (L < 1.22 * (n - 1) * a)
            throw std::invalid_argument("Invalid argument: L. Must be greater than 1.22(n-1)a.");
        if (a < 0.)
            throw std::invalid_argument("Invalid argument: a. Must be positive.");
        if (T0 < 0.)
            throw std::invalid_argument("Invalid argument: T0. Must be positive.");
        if (tau < 0. || tau > 1e-2)
            throw std::invalid_argument("Invalid argument: tau. Must be between 0 and 1e-2.");
        if (So < 0 || So > Sd)
            throw std::invalid_argument("Invalid argument: So. Must be between 0 and Sd.");
        if (Sd < 0)
            throw std::invalid_argument("Invalid argument: Sd. Must be positive.");
        if (Sout < 0 || Sout > Sd)
            throw std::invalid_argument("Invalid argument: Sout. Must be between 0 and Sd.");
        if (Sxyz < 0 || Sxyz > Sd)
            throw std::invalid_argument("Invalid argument: Sxyz. Must be between 0 and Sd.");

        std::cout << "`setParameters()` :> Successfully set parameters from ../Config/" << filename << '\n';

        // N and K are still the same so I can carefully release the memory
        delete[] b0;
        delete[] b1;
        delete[] b2;

        delete[] p;
        delete[] pAbs;
        delete[] Vs;

        for (usint i = 0; i < N; i++)
        {
            delete[] r0[i];
            delete[] p0[i];
            delete[] Vp[i];
            delete[] Fs[i];
            delete[] Fi[i];

            for (usint j = 0; j < N; j++)
                delete[] Fp[i][j];

            delete[] Fp[i];
        }

        delete[] r0;
        delete[] p0;
        delete[] Vp;
        delete[] Fs;
        delete[] Fi;

        delete[] Fp;

        // Here I can set the new values
        N = n * n * n;
        K = 3;

        b0 = new double[K]{a, 0., 0.};
        b1 = new double[K]{a * 0.5, a * sqrt(3.) * 0.5, 0.};
        b2 = new double[K]{a * 0.5, a * sqrt(3.) / 6., a * sqrt(6.) / 3.};

        p = new double[K]();
        pAbs = new double[N]();
        Vs = new double[N]();

        r0 = new double *[N]();
        p0 = new double *[N]();
        Vp = new double *[N]();
        Fs = new double *[N]();
        Fi = new double *[N]();

        Fp = new double **[N]();

        for (usint i = 0; i < N; i++)
        {
            r0[i] = new double[K]();
            p0[i] = new double[K]();
            Vp[i] = new double[N]();
            Fs[i] = new double[K]();
            Fi[i] = new double[K]();

            Fp[i] = new double *[N]();

            for (usint j = 0; j < N; j++)
                Fp[i][j] = new double[K]();
        }

        input.close();
        std::cout << "`setParameters()` :> Successfully reallocated memory for new parameters.\n\n";
    }
    catch (const std::invalid_argument &error)
    {
        // Notice I do not need to reallocate memory because of default parameters
        // and buffer have the same sizes
        n = 6;
        So = 5000;
        Sd = 50000;
        Sout = 500;
        Sxyz = 500;
        m = 40.;
        R = 0.38;
        e = 1.;
        k = 8.31e-3;
        f = 1e4;
        L = 6.;
        a = 0.38;
        T0 = 1e4;
        tau = 1e-3;

        input.close();
        std::cerr << "`setParameters()` :> Exception while setting parameters from ../Config/" << filename << '\n';
        std::cerr << "`setParameters()` :> " << error.what() << '\n';
        std::cerr << "`setParameters()` :> Values are set to default now.\n\n";
    }
    catch (const std::ifstream::failure &error)
    {
        n = 6;
        So = 5000;
        Sd = 50000;
        Sout = 500;
        Sxyz = 500;
        m = 40.;
        R = 0.38;
        e = 1.;
        k = 8.31e-3;
        f = 1e4;
        L = 6.;
        a = 0.38;
        T0 = 1e4;
        tau = 1e-3;

        std::cerr << "`setParameters()` :> " << error.what() << '\n';
        std::cerr << "`setParameters()` :> Values are set to default now.\n\n";
    }
}

/**************************************************************************************
 * This function prints all currently set parameters. Notice that the section
 * parameters is only for the information of printed parameters.
 * This function DOES NOT accept parameters.
 * @param usint n     // Number of atoms along the crystal edge
 * @param usint So    // Thermalisation steps
 * @param usint Sout  // Save informations about the system every \p`Sout` steps
 * @param usint Sxyz  // Save positions of atoms every `Sxyz` steps
 * @param double m    // Mass of the single atom
 * @param double e    // Minimum of the potential
 * @param double R    // Interatomic distance for which occurs minimum of the potential
 * @param double k    // Boltzmann constant
 * @param double f    // Elastic coefficient
 * @param double L    // Radius of sphere which confines atoms
 * @param double a    // Interatomic distance
 * @param double T0   // Initial temperature
 * @param double tau  // Integration step
 * @return Nothing to return.
 **************************************************************************************/
void Argon::checkParameters() const noexcept
{
    std::cout << "`checkParameters()` :> Currently set parameters." << '\n';
    std::cout << "`checkParameters()` :> n:        " << n << '\n';
    std::cout << "`checkParameters()` :> m:        " << m << '\n';
    std::cout << "`checkParameters()` :> e:        " << e << '\n';
    std::cout << "`checkParameters()` :> R:        " << R << '\n';
    std::cout << "`checkParameters()` :> k:        " << k << '\n';
    std::cout << "`checkParameters()` :> f:        " << f << '\n';
    std::cout << "`checkParameters()` :> L:        " << L << '\n';
    std::cout << "`checkParameters()` :> a:        " << a << '\n';
    std::cout << "`checkParameters()` :> T0:       " << T0 << '\n';
    std::cout << "`checkParameters()` :> tau:      " << tau << '\n';
    std::cout << "`checkParameters()` :> So:       " << So << '\n';
    std::cout << "`checkParameters()` :> Sd:       " << Sd << '\n';
    std::cout << "`checkParameters()` :> Sout:     " << Sout << '\n';
    std::cout << "`checkParameters()` :> Sxyz:     " << Sxyz << '\n';
    std::cout << "`checkParameters()` :> End of parameters.\n\n";
}

/**************************************************************************************
 * This function calculates initial positions, momenta, forces and potentials acting on
 * atoms in initial state. Primarily we calculate there the trapping potentials,
 * repulsion related to sphere walls, total forces impact to particles, van der Waals
 * interactions, interaction forces between atoms and total potential. Next that
 * function calculates total energy (Hamiltonian), initial real temperature (T) and
 * initial pressure related to sphere walls. Finally it saves all needed informations
 * to the given files.
 * @param char* filename where to save initial positions,
 * @param char* filename where to save initial momenta,
 * @param char* filename where to save initial H, T and P.
 * @return Nothing to return.
 **************************************************************************************/
void Argon::initialState(const char *rFilename, const char *pFilename, const char *htpFilename) noexcept
{
    // Calculate initial positions of atoms (5)
    for (usint i_0 = 0; i_0 < n; i_0++)
    {
        for (usint i_1 = 0; i_1 < n; i_1++)
        {
            for (usint i_2 = 0; i_2 < n; i_2++)
            {
                for (usint j = 0; j < K; j++)
                {
                    usint i = i_0 + i_1 * n + i_2 * n * n;
                    r0[i][j] = (i_0 - 0.5 * (n - 1)) * b0[j] + (i_1 - 0.5 * (n - 1)) * b1[j] + (i_2 - 0.5 * (n - 1)) * b2[j];
                }
            }
        }
    }

    // Calculate initial momenta of atoms (7)
    for (usint i = 0; i < N; i++)
    {
        for (usint j = 0; j < K; j++)
        {
            // Random number between 0 and 1
            double number = static_cast<double>(mt()) / static_cast<double>(mt.max());
            // Random sign of momentum variable
            char sign = mt() % 2;

            // Calculate momentum
            if (sign == 0)
                p0[i][j] = -sqrt(-0.5 * k * T0 * log(number) * 2. * m); // (7)
            else if (sign == 1)
                p0[i][j] = +sqrt(-0.5 * k * T0 * log(number) * 2. * m); // (7)

            // Accumulate momenta
            p[j] += p0[i][j];
        }
    }

    // Eliminate the centre of mass movement (8) and immediately calculate absolute values
    for (usint i = 0; i < N; i++)
    {
        for (usint j = 0; j < K; j++)
        {
            p0[i][j] = p0[i][j] - (p[j] / N);
            pAbs[i] += p0[i][j] * p0[i][j];
        }

        pAbs[i] = sqrt(pAbs[i]);
    }

    // Calculate initial forces and potentials affecting to atoms
    V = 0.; // Total potential energy
    // Forces and potentials loop
    for (usint i = 0; i < N; i++)
    {
        // Absolute value of r_i -> |r_i|
        double r_i = sqrt(r0[i][0] * r0[i][0] + r0[i][1] * r0[i][1] + r0[i][2] * r0[i][2]);

        // (10)
        if (r_i < L)
            Vs[i] = 0.;
        else if (r_i >= L)
            Vs[i] = 0.5 * f * (r_i - L) * (r_i - L);

        // Accumulate potential related to sphere walls to total potential
        V += Vs[i];

        // (14)
        for (usint j = 0; j < K; j++)
        {
            if (r_i < L)
                Fs[i][j] = 0.;
            else if (r_i >= L)
                Fs[i][j] = f * (L - r_i) * r0[i][j] / r_i;

            // Accumulate repulsive forces related to sphere walls to total forces
            Fi[i][j] = Fs[i][j];
        }

        // (9) and (13)
        for (usint j = 0; j < i; j++)
        {
            // Absolute value of r_i - r_j -> |r_i - r_j|
            double r_ij = sqrt((r0[i][0] - r0[j][0]) * (r0[i][0] - r0[j][0]) + (r0[i][1] - r0[j][1]) * (r0[i][1] - r0[j][1]) +
                               (r0[i][2] - r0[j][2]) * (r0[i][2] - r0[j][2]));

            // Local variables to evaluate powers -> huge increase of performance (instead of calculate with common pow())
            double y = (R / r_ij) * (R / r_ij);
            double x = y * y * y;
            // (9)
            Vp[i][j] = e * x * (x - 2.);

            for (usint k = 0; k < K; k++)
            {
                // Fp is not required as 3D array, we may use just ordinary variable but with Fp is more evident
                // what is happens here
                Fp[i][j][k] = 12. * e * x * (x - 1.) * (r0[i][k] - r0[j][k]) / (r_ij * r_ij);

                // Symmetry of forces matrix (only one triangular matrix needs to be calculated) -> increase performance
                Fi[i][k] += Fp[i][j][k];
                Fi[j][k] -= Fp[i][j][k];
            }

            // Accumulate van der Waals potentials
            V += Vp[i][j];
        }
    }

    // Prepare to accumulate physical parameters at initial time
    H = V; // At this moment Hamiltonian is just total potential
    T = 0.;
    P = 0.;
    Ek = 0.;

    for (usint i = 0; i < N; i++)
    {
        // Local variable to increase performance;
        Ek = (p0[i][0] * p0[i][0] + p0[i][1] * p0[i][1] + p0[i][2] * p0[i][2]) / (2. * m);

        // Accumulate physical parameters
        H += Ek;
        T += 2. / (3. * N * k) * Ek;
        P += sqrt(Fs[i][0] * Fs[i][0] + Fs[i][1] * Fs[i][1] + Fs[i][2] * Fs[i][2]) / (4. * M_PI * L * L);
    }

    initialStateCheck = true;
    saveInitialState(rFilename, pFilename, htpFilename);
    std::cout << "`initialState()` :> Successfully calculated and saved initial state.\n\n";
}

/**************************************************************************************
 * This function carries out the dynamics of the whole simulation. Primarily it
 * calculates required positions, momenta, forces and potentials acting on atoms at
 * given moment in time. Similarly as in the initial state, we calculate here
 * the trapping potentials, repulsion related to sphere walls, total forces impact to
 * particles, van der Waals interactions, interaction forces between atoms and total
 * potential. Next that function calculates total energy (Hamiltonian), real temperature
 * (T) and pressure (P) related to sphere walls.
 * @param char* filename where to save current positions,
 * @param char* filename where to save current H, T and P.
 * @return Nothing to return.
 *************************************************************************************/
void Argon::simulateDynamics(const char *rFilename, const char *htpFilename) noexcept
{
    if (initialStateCheck == false)
    {
        std::cerr << "`simulateDynamics()` :> Error - calculate initial state before!\n\n";
        return;
    }

    std::cout << "`simulateDynamics()` :> System is ready to simulation.\n\n";

    std::ofstream ofileRt("../Out/" + std::string(rFilename), std::ios::out);
    std::ofstream ofileHtp("../Out/" + std::string(htpFilename), std::ios::out);

    ofileRt << std::fixed << std::setprecision(5);
    ofileHtp << std::fixed << std::setprecision(5);

    // Save initial positions and initial H, T and P
    saveCurrentPositions(ofileRt);
    saveCurrentHTP(0., ofileHtp);

    // Current information to track simulation
    printCurrentInfo(0.);

    // At this point, these values are not computed
    Hmean = 0.;
    Tmean = 0.;
    Pmean = 0.;
    Vol = 4. / 3. * M_PI * L * L * L;

    // Informations print interval
    uint infoOut = Sd / 10;

    // Simulation loop
    for (uint s = 1; s <= So + Sd; s++)
    {
        // Print current informations
        if (s % infoOut == 0)
        {
            printCurrentInfo(s * tau);
        }

        // Calculate auxiliary momenta (18a) and positions (18b)
        for (usint i = 0; i < N; i++)
        {
            for (usint j = 0; j < K; j++)
            {
                p0[i][j] = p0[i][j] + 0.5 * Fi[i][j] * tau;
                r0[i][j] = r0[i][j] + p0[i][j] * tau / m;
            }
        }

        // In every step set total potential to zero (IMPORTANT!)
        V = 0.;

        // Dynamics loop
        for (usint i = 0; i < N; i++)
        {
            // Absolute value of r_i -> |r_i|
            double r_i = sqrt(r0[i][0] * r0[i][0] + r0[i][1] * r0[i][1] + r0[i][2] * r0[i][2]);

            // (10)
            if (r_i < L)
                Vs[i] = 0.;
            else if (r_i >= L)
                Vs[i] = 0.5 * f * (r_i - L) * (r_i - L);

            // Accumulate potential related to sphere walls to total potential
            V += Vs[i];

            // (14)
            for (usint j = 0; j < K; j++)
            {
                if (r_i < L)
                    Fs[i][j] = 0.;
                else if (r_i >= L)
                    Fs[i][j] = f * (L - r_i) * r0[i][j] / r_i;

                // Accumulate repulsive forces related to sphere walls to total forces
                Fi[i][j] = Fs[i][j];
            }

            for (usint j = 0; j < i; j++)
            {
                // Absolute value of r_i - r_j -> |r_i - r_j|
                double r_ij = sqrt((r0[i][0] - r0[j][0]) * (r0[i][0] - r0[j][0]) + (r0[i][1] - r0[j][1]) * (r0[i][1] - r0[j][1]) +
                                   (r0[i][2] - r0[j][2]) * (r0[i][2] - r0[j][2]));

                // Local variables to evaluate powers -> huge increase of performance (instead of calculate with common pow())
                double y = (R / r_ij) * (R / r_ij);
                double x = y * y * y;
                // (9)
                Vp[i][j] = e * x * (x - 2.);

                for (usint k = 0; k < K; k++)
                {
                    Fp[i][j][k] = 12. * e * x * (x - 1.) * (r0[i][k] - r0[j][k]) / (r_ij * r_ij);

                    // Symmetry of forces matrix (only one triangular matrix needs to be calculated) -> increase performance
                    Fi[i][k] += Fp[i][j][k];
                    Fi[j][k] -= Fp[i][j][k];
                }

                // Accumulate van der Waals potentials
                V += Vp[i][j];
            }
        } // End of dynamics loop

        // Calculate momenta (18c) and absolute values
        for (usint i = 0; i < N; i++)
        {
            for (usint j = 0; j < K; j++)
            {
                p0[i][j] = p0[i][j] + 0.5 * Fi[i][j] * tau;
                pAbs[i] += p0[i][j] * p0[i][j];
            }

            pAbs[i] = sqrt(pAbs[i]);
        }

        calculateCurrentHTP();

        // Save temporary positions at given time
        if (s % Sxyz == 0)
        {
            saveCurrentPositions(ofileRt);
        }

        // Save temporary H, T and P at given time
        if (s % Sout == 0)
        {
            saveCurrentHTP(s * tau, ofileHtp);
        }

        // Accumulate mean values only when thermalisation is done
        if (s >= So)
        {
            Tmean += T;
            Pmean += P;
            Hmean += H;
        }
    } // End of simulation loop

    printCurrentInfo((So + Sd) * tau); // Latest step

    // Average the cumulative values
    Hmean /= Sd;
    Tmean /= Sd;
    Pmean /= Sd;

    // Ideal gas law -> PV = NkT Volume not potential :)
    IdealGas = (N * k * Tmean) / (Pmean * Vol);

    // Chemical potential from microcanonical ensemble
    u = k * T * log(Vol / N * (4. * M_PI * m * Hmean) / (3. * N) * sqrt((4. * M_PI * m * Hmean) / (3. * N)));

    std::cout << "Mean Total Energy:        " << Hmean << '\n';
    std::cout << "Mean Temperature:         " << Tmean << '\n';
    std::cout << "Mean Pressure:            " << Pmean << '\n';
    std::cout << "Ideal Gas Law:            " << IdealGas << '\n';
    std::cout << "Mean Chemical Potential:  " << u << '\n';
    std::cout << '\n';

    ofileRt.close();
    ofileHtp.close();
}

/**************************************************************************************
 * This function calculates absolute value of momentum for every particle.
 * @return std::tuple<double *, usint, double, double, double> - where the first
 * parameter is a pointer to array with the absolute momentum values, the second is the
 * size of this array, third is the temperature related to this calculated state, fourth
 * is the Boltzmann constant and the fifth is the particle mass.
 *************************************************************************************/
std::tuple<double *, usint, double, double, double> Argon::getMomentumAbs() const noexcept
{
    double *pAbsToReturn = new double[N]();

    for (usint i = 0; i < N; i++)
    {
        pAbsToReturn[i] = pAbs[i];
    }

    return std::make_tuple(pAbsToReturn, N, T, k, m);
}

/**************************************************************************************
 * This function calculates current Hamiltonian, Temperature and Pressure of the
 * system. It uses current momenta and sphere repulsion to this.
 * @return Calculates Hamiltonian, Temperature and Pressure of the system.
 *************************************************************************************/
void Argon::calculateCurrentHTP() noexcept
{
    // Prepare to accumulate physical parameters at initial time
    H = V; // At this moment Hamiltonian is just total potential
    T = 0.;
    P = 0.;
    Ek = 0.;

    for (usint i = 0; i < N; i++)
    {
        // Local variable to increase performance;
        Ek = (p0[i][0] * p0[i][0] + p0[i][1] * p0[i][1] + p0[i][2] * p0[i][2]) / (2. * m);

        // Accumulate physical parameters
        H += Ek;
        T += 2. / (3. * N * k) * Ek;
        P += sqrt(Fs[i][0] * Fs[i][0] + Fs[i][1] * Fs[i][1] + Fs[i][2] * Fs[i][2]) / (4. * M_PI * L * L);
    }
}

/**************************************************************************************
 * Writes to file current time, Hamiltonian, Temperature and Pressure of the system.
 * @param double current time,
 * @param ofstream file where to save current H, T and P.
 * @return Set subsequent lines with current time, H, T and P in the given file.
 *************************************************************************************/
inline void Argon::saveCurrentHTP(const double &time, std::ofstream &ofileHtp) noexcept
{
    ofileHtp << time << '\t' << H << '\t' << T << '\t' << P << '\n';
}

/**************************************************************************************
 * Writes to file current positions of the atoms.
 * @param ofstream file where to save current positions.
 * @return Set subsequent positions of particles in the given file.
 *************************************************************************************/
void Argon::saveCurrentPositions(std::ofstream &ofileRt) noexcept
{
    // Number of atoms to read by Jmol
    ofileRt << N;
    ofileRt << "\n\n";

    for (usint i = 0; i < N; i++)
    {
        ofileRt << "AR\t";

        for (usint j = 0; j < K; j++)
            ofileRt << r0[i][j] << '\t';
        ofileRt << '\n';
    }

    ofileRt << '\n';
}

/**************************************************************************************
 * This function saves initial state to the given files.
 * @param char* filename where to save initial positions,
 * @param char* filename where to save initial momenta,
 * @param char* filename where to save initial H, T and P.
 * @return Set appropriate positions, momenta and H, T, P in given files.
 *************************************************************************************/
void Argon::saveInitialState(const char *rFilename, const char *pFilename, const char *htpFilename) const noexcept
{
    std::ofstream rOut("../Out/" + std::string(rFilename), std::ios::out);
    std::ofstream pOut("../Out/" + std::string(pFilename), std::ios::out);
    std::ofstream htpOut("../Out/" + std::string(htpFilename), std::ios::out);

    rOut << std::fixed << std::setprecision(5);
    pOut << std::fixed << std::setprecision(5);

    // Number of atoms in header to read by Jmol
    rOut << N << "\n\n";
    pOut << N << '\t' << T << '\t' << k << '\t' << m << "\n\n";

    for (usint i = 0; i < N; i++)
    {
        // AR beacuse of we analyse Argon gas
        rOut << "AR\t";
        pOut << "AR\t";

        for (usint j = 0; j < K; j++)
        {
            rOut << r0[i][j] << '\t';
            pOut << p0[i][j] << '\t';
        }

        pOut << pAbs[i] << '\n';
        rOut << '\n';
    }

    // Header with physical parameters
    htpOut << std::fixed << std::setprecision(5);
    htpOut << "t (ps)\t";
    htpOut << "H (kJ/mol)\t";
    htpOut << "T (K)\t";
    htpOut << "P (atm)\n";
    htpOut << 0 << '\t' << H << '\t' << T << '\t' << P << '\n';

    rOut.close();
    pOut.close();
    htpOut.close();
}

/**************************************************************************************
 * Checks if the input file is empty.
 * @param ifstream input.
 * @return True if file is empty, otherwise false.
 *************************************************************************************/
inline bool Argon::fileIsEmpty(std::ifstream &input) const noexcept
{
    return input.peek() == std::ifstream::traits_type::eof();
}

/**************************************************************************************
 * Print current informations about the system while simulation is in progress.
 * @param double current time.
 * @return Prints current informations.
 *************************************************************************************/
void Argon::printCurrentInfo(const double &time) const noexcept
{
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Current Time:             " << time << '\n';
    std::cout << "Current Total Energy:     " << H << '\n';
    std::cout << "Current Total Potential:  " << V << '\n';
    std::cout << "Current Temperature:      " << T << '\n';
    std::cout << "Current Pressure:         " << P << '\n';
    std::cout << '\n';
}
