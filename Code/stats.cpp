#define _USE_MATH_DEFINES
#include "stats.h"
#include <cmath>
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

void Stats::evaluateHist()
{
    // Evaluate histogram related to particles momentum
    distributionMean = 0.;
    distributionSigma = 0.;

    // Calculate range of histogram and its bins
    const double histRange = abs(low - up);
    const double binRange = histRange / bins;

    // Set bin ranges
    binRanges[0] = low;
    binRanges[bins] = up;

    for (int i = 1; i < bins; i++)
        binRanges[i] = binRanges[i - 1] + binRange;

    // Count how many particles are in the specific momentum range
    for (int i = 0; i < N; i++)
    {
        if (pAbs[i] < low)
        {
            ++underflow;
        }
        else if (pAbs[i] > up)
        {
            ++overflow;
        }
        else
        {
            for (int j = 0; j < bins; j++)
            {
                if (pAbs[i] > binRanges[j] && pAbs[i] <= binRanges[j + 1])
                {
                    stars[j] += '*';
                    break;
                }
            }
        }

        // Accumulate momentum to calculate distribution mean value
        distributionMean += pAbs[i];
    }

    distributionMean /= N;

    // Calculate standard deviation of momentum
    for (int i = 0; i < N; i++)
        distributionSigma += (pAbs[i] - distributionMean) * (pAbs[i] - distributionMean);

    distributionSigma = sqrt(distributionSigma / N);

    // Print evaluated histogram with its basic statistics
    std::cout << '\n';
    std::cout << "Bins:         " << bins << '\n';
    std::cout << "Low:          " << low << '\n';
    std::cout << "Up:           " << up << '\n';
    std::cout << "Underflow:    " << underflow << '\n';
    std::cout << "Overflow:     " << overflow << '\n';
    std::cout << "Mean:         " << distributionMean << '\n';
    std::cout << "StdDev:       " << distributionSigma << '\n';
    std::cout << '\n';

    // Print options
    std::cout << std::fixed << std::showpos << std::setprecision(3);

    // Create caption with bin range e.g (+0.123; +0.456]:
    std::stringstream binRangeS;
    binRangeS << std::fixed << std::showpos << std::setprecision(3);
    binRangeS << "(" << binRanges[bins - 1] << "; " << binRanges[bins] << "]: ";
    // Calculate maximum length of binRangeS string
    usint captionLen = binRangeS.str().length();
    binRangeS.str("");

    // Print histogram of particles momentum
    for (int i = 0; i < bins; i++)
    {
        binRangeS << "(" << binRanges[i] << "; " << binRanges[i + 1] << "]: ";
        std::cout << std::setw(captionLen - binRangeS.str().length() + 1) << "(" << binRanges[i] << "; " << binRanges[i + 1] << "]: ";
        std::cout << std::setw(3) << stars[i].length() << ' ' << stars[i] << '\n';
        binRangeS.str("");
    }

    // std::cout << sqrt(8. * k * T / (M_PI * m)) * m;
}

void Stats::setInputFromArgon(const double *pAbsArgon, const usint &NArgon, const double &TArgon, const double &KArgon, const double &MArgon)
{
    N = NArgon;
    T = TArgon;
    k = KArgon;
    m = MArgon;

    // delete[] pAbs;
    pAbs = new double[N];

    for (usint i = 0; i < N; i++)
        pAbs[i] = pAbsArgon[i];
}