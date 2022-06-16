#ifndef STATS_H
#define STATS_H
#include <string>
typedef unsigned short int usint;

class Stats
{
private:
    double low;
    double up;
    usint underflow;
    usint overflow;
    usint bins;
    double *binRanges;
    std::string *stars;
    usint maxStarsIndex;
    usint maxStars;

    double distributionMean;
    double distributionMeanSq;
    double distributionSigma;

    void calculateHist();

    double *pAbs;
    usint N;
    double T;
    double k;
    double m;

    double pProEmp;    ///< Most probable momentum obtained empirically from argon library
    double pMeanEmp;   ///< Mean momentum obtained empirically from argon library
    double pMeanSqEmp; ///< Mean square momentum obtained empirically from argon library
    double pProMB;     ///< Most probable momentum obtained analytically from Maxwell-Boltzmann distribution
    double pMeanMB;    ///< Mean momentum obtained analytically from Maxwell-Boltzmann distribution
    double pMeanSqMB;  ///< Mean square momentum obtained analytically from Maxwell-Boltzmann distribution

public:
    Stats();
    Stats(const double &Low, const double &Up, const usint &Bins);
    Stats(const Stats &H);
    ~Stats();

    void setStats(const double &Low, const double &Up, const usint &Bins);
    void evaluateHist();
    void calculateStats();
    void setInputFromArgon(const double *pAbs, const usint &N, const double &T, const double &K, const double &M);
};

#endif // STATS_H