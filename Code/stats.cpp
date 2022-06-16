#include "stats.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

Stats::Stats() : low(0), up(0), underflow(0), overflow(0), bins(0)
{
    binRanges = new double[0];
    stars = new std::string[0];
}

Stats::Stats(const double &Low, const double &Up, const usint &Bins) : low(Low), up(Up), underflow(0), overflow(0), bins(Bins)
{
    binRanges = new double[bins + 1];
    stars = new std::string[bins];
}

Stats::Stats(const Stats &H) : low(H.low), up(H.up), underflow(H.underflow), overflow(H.overflow), bins(H.bins)
{
    binRanges = new double[bins + 1];
    stars = new std::string[bins];
}

Stats::~Stats()
{
    delete[] binRanges;
    delete[] stars;
}

void Stats::setStats(const double &Low, const double &Up, const usint &Bins)
{
    low = Low;
    up = Low;
    bins = Bins;

    delete[] binRanges;
    delete[] stars;

    binRanges = new double[bins + 1];
    stars = new std::string[bins];
}

void Stats::printStats()
{
    calculate();

    std::cout << '\n';
    std::cout << "Bins:         " << bins << '\n';
    std::cout << "Low:          " << low << '\n';
    std::cout << "Up:           " << up << '\n';
    std::cout << "Underflow:    " << underflow << '\n';
    std::cout << "Overflow:     " << overflow << '\n';
    std::cout << "Mean:         " << distributionMean << '\n';
    std::cout << "StdDev:       " << distributionSigma << '\n';
    std::cout << '\n';
    std::cout << std::fixed << std::showpos << std::setprecision(3);

    std::stringstream ss;
    ss << std::fixed << std::showpos << std::setprecision(3);
    ss << "(" << *(binRanges + bins - 1) << "; " << *(binRanges + bins) << "]: ";
    usint ssSize = ss.str().length();
    ss.str("");

    for (int i = 0; i < bins; i++)
    {
        ss << "(" << *(binRanges + i) << "; " << *(binRanges + i + 1) << "]: ";
        // std::cout << ssSize - ss.str().length() << '\n';
        std::cout << std::setw(ssSize - ss.str().length() + 1) << "(" << *(binRanges + i) << "; " << *(binRanges + i + 1) << "]: ";
        std::cout << std::setw(3) << (*(stars + i)).length() << "  " << *(stars + i) << '\n';
        ss.str("");
    }

    std::cout << sqrt(2. * 8.31e-3 * T / (40.)) * 40.;
}

void Stats::calculate()
{
    distributionMean = 0.;
    distributionSigma = 0.;

    const double histRange = abs(low - up);
    const double binRange = histRange / bins;

    *(binRanges) = low;
    *(binRanges + bins) = up;

    for (int i = 1; i < bins; i++)
        *(binRanges + i) = *(binRanges + i - 1) + binRange;

    for (int i = 0; i < N; i++)
    {
        if (*(pAbs + i) < low)
            ++underflow;
        else if (*(pAbs + i) > up)
            ++overflow;
        else
        {
            for (int j = 0; j < bins; j++)
            {
                if (*(pAbs + i) > *(binRanges + j) && (*(pAbs + i) <= *(binRanges + j + 1)))
                {
                    *(stars + j) += '*';
                    break;
                }
            }
        }

        distributionMean += *(pAbs + i);
    }

    distributionMean /= N;

    for (int i = 0; i < N; i++)
        distributionSigma += (*(pAbs + i) - distributionMean) * (*(pAbs + i) - distributionMean);

    distributionSigma = sqrt(distributionSigma / N);
}

void Stats::setInputFromArgon(const double *pAbsArgon, const usint &NArgon, const double &TArgon)
{
    N = NArgon;
    T = TArgon;

    //delete[] pAbs;
    pAbs = new double[N];

    for (usint i = 0; i < N; i++)
        pAbs[i] = pAbsArgon[i];
}