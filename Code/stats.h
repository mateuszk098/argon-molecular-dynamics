#ifndef STATS_H
#define STATS_H
#include <string>
typedef unsigned short int usint;

class Stats
{
private:
    /// Variables related to histogram printing
    double low;          ///< Minimum value of the histogram range
    double up;           ///< Maximum value of the histogram range
    usint underflow;     ///< Number of samples under minimum value `low`
    usint overflow;      ///< Number of samples overmaximum value `up`
    usint bins;          ///< Number of bins in the histogram
    usint maxStarsIndex; ///< Index of `binRanges` where is the most counts
    usint maxStars;      ///< Number of max counts in the `binRanges`
    double *binRanges;   ///< 1D array with ranges of the histogram
    std::string *stars;  ///< 1D array with the stars '*' to visualise histogram

    /// Variables related to histogram statistics
    double distributionMean;   ///< Mean value of the distribution
    double distributionMeanSq; ///< Mean squared value of the distribution
    double distributionSigma;  ///< Standard deviation of the distribution

    /// Variables related to the system
    double *pAbs; ///< Absolute value of the particles momentum evaluated by Argon library
    usint N;      ///< Number of particles
    double T;     ///< Temperature evaluated by Argon library
    double k;     ///< Boltzmann constant
    double m;     ///< Mass of the single particle

    /// Variables related to Maxwell-Boltzmann statistics
    double pProEmp;    ///< Most probable momentum obtained empirically from argon library
    double pMeanEmp;   ///< Mean momentum obtained empirically from argon library
    double pMeanSqEmp; ///< Mean square momentum obtained empirically from argon library
    double pProMB;     ///< Most probable momentum obtained analytically from Maxwell-Boltzmann distribution
    double pMeanMB;    ///< Mean momentum obtained analytically from Maxwell-Boltzmann distribution
    double pMeanSqMB;  ///< Mean square momentum obtained analytically from Maxwell-Boltzmann distribution
    double EkEmp;      ///< Kinetic energy from ordinary Newton formula
    double EkMB;       ///< Kinetic energy from kinetic theory of gases

public:
    Stats();
    Stats(const double &Low, const double &Up, const usint &Bins);
    Stats(const Stats &H);
    ~Stats();

    void setStats(const double &Low, const double &Up, const usint &Bins);
    void evaluateHist(const char *histFilename);
    void calculateStats();
    void setInputFromArgon(const double *pAbs, const usint &N, const double &T, const double &K, const double &M);
};

#endif // STATS_H