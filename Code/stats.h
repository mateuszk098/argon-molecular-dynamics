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

    void calculate(const double *Values, const usint &Size);

public:
    Stats();
    Stats(const double &Low, const double &Up, const usint &Bins);
    Stats(const Stats &H);
    ~Stats();

    void setStats(const double &Low, const double &Up, const usint &Bins);
    void printStats(const double *Values, const usint &Size);
};

#endif // STATS_H