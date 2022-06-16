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

    void calculateHist();

    double *pAbs;
    usint N;
    double T;
    double k;
    double m;

public:
    Stats();
    Stats(const double &Low, const double &Up, const usint &Bins);
    Stats(const Stats &H);
    ~Stats();

    void setStats(const double &Low, const double &Up, const usint &Bins);
    void evaluateHist();
    void setInputFromArgon(const double *pAbs, const usint &N, const double &T, const double &K, const double &M);
};

#endif // STATS_H