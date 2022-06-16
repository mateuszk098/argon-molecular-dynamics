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

    double distributionMean;
    double distributionSigma;

    void calculate();

    double *pAbs;
    usint N;
    double T;

public:
    Stats();
    Stats(const double &Low, const double &Up, const usint &Bins);
    Stats(const Stats &H);
    ~Stats();

    void setStats(const double &Low, const double &Up, const usint &Bins);
    void printStats();
    void setInputFromArgon(const double *pAbs, const usint &N, const double &T);
};

#endif // STATS_H