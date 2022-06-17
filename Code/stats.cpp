#define _USE_MATH_DEFINES
#include "stats.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <fstream>

Stats::Stats() noexcept : low(0.), up(0.), underflow(0), overflow(0), bins(0), maxStarsIndex(0), maxStars(0),
                          distributionMean(0.), distributionMeanSq(0.), distributionSigma(0.), N(0), T(0.), k(0.),
                          m(0.), pProEmp(0.), pMeanEmp(0.), pMeanSqEmp(0.), pProMB(0.), pMeanMB(0.), pMeanSqMB(0.),
                          EkEmp(0.), EkMB(0.)
{
    binRanges = new double[bins + 1]();
    stars = new std::string[bins]();
    pAbs = new double[N]();
}

Stats::~Stats() noexcept
{
    delete[] binRanges;
    delete[] stars;
    delete[] pAbs;
}

void Stats::evaluateHist(const char *histFilename)
{
    // Calculate range of histogram and its bins
    const double histRange = std::abs(low - up);
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
        distributionMeanSq += pAbs[i] * pAbs[i];
    }

    distributionMean /= N;
    distributionMeanSq = std::sqrt(distributionMeanSq / N);

    // Calculate standard deviation of momentum
    for (int i = 0; i < N; i++)
        distributionSigma += (pAbs[i] - distributionMean) * (pAbs[i] - distributionMean);

    distributionSigma = std::sqrt(distributionSigma / N);

    std::ofstream histOfile("../Out/" + std::string(histFilename), std::ios::out);
    histOfile << std::fixed << std::setprecision(5); // << std::showpos

    // Print evaluated histogram with its basic statistics
    histOfile << "Momentum Distribution:\n\n";
    histOfile << "Bins:         " << bins << '\n';
    histOfile << "Low:          " << low << '\n';
    histOfile << "Up:           " << up << '\n';
    histOfile << "Underflow:    " << underflow << '\n';
    histOfile << "Overflow:     " << overflow << '\n';
    histOfile << "Mean:         " << distributionMean << '\n';
    histOfile << "StdDev:       " << distributionSigma << '\n';
    histOfile << '\n';

    // Create caption with bin range e.g (+0.123; +0.456]:
    std::stringstream binRangeS;
    binRangeS << std::fixed << std::setprecision(2); // << std::showpos
    binRangeS << "(" << binRanges[bins - 1] << "; " << binRanges[bins] << "]: ";
    // Calculate maximum length of binRangeS string
    usint captionLen = binRangeS.str().length();
    binRangeS.str("");

    // Print histogram of particles momentum
    for (int i = 0; i < bins; i++)
    {
        binRangeS << "(" << binRanges[i] << "; " << binRanges[i + 1] << "]: ";
        histOfile << std::setw(captionLen - binRangeS.str().length() + 1) << "(" << binRanges[i] << "; " << binRanges[i + 1] << "]: ";
        histOfile << std::setw(3) << stars[i].length() << ' ' << stars[i] << '\n';

        if (stars[i].length() > maxStars)
        {
            maxStarsIndex = i;
            maxStars = stars[i].length();
        }

        binRangeS.str("");
    }

    pProMB = std::sqrt(2. * k * T * m);
    pMeanMB = std::sqrt(8. * k * T * m / M_PI);
    pMeanSqMB = std::sqrt(3. * k * T * m);
    EkMB = 3. / 2. * k * T;

    for (usint i = 0; i < N; i++)
    {
        if (pAbs[i] > binRanges[maxStarsIndex] && pAbs[i] <= binRanges[maxStarsIndex + 1])
        {
            pProEmp += pAbs[i];
        }
    }

    pProEmp /= maxStars;
    pMeanEmp = distributionMean;
    pMeanSqEmp = distributionMeanSq;
    EkEmp = pMeanEmp * pMeanEmp / (2. * m);

    histOfile << std::fixed << std::setprecision(5);
    histOfile << '\n';
    histOfile << "pProMB:       " << pProMB << '\t' << "Most probable momentum obtained analytically.\n";
    histOfile << "pProEmp:      " << pProEmp << '\t' << "Most probable momentum obtained numerically.\n";
    histOfile << "Error (%):    " << std::abs(pProEmp - pProMB) / pProMB * 100. << '\n';
    histOfile << '\n';
    histOfile << "pMeanMB:      " << pMeanMB << '\t' << "Mean momentum obtained analytically.\n";
    histOfile << "pMeanEmp:     " << pMeanEmp << '\t' << "Mean momentum obtained numerically.\n";
    histOfile << "Error (%):    " << std::abs(pMeanEmp - pMeanMB) / pMeanMB * 100. << '\n';
    histOfile << '\n';
    histOfile << "pMeanSqMB:    " << pMeanSqMB << '\t' << "Mean square momentum obtained analytically.\n";
    histOfile << "pMeanSqEmp:   " << pMeanSqEmp << '\t' << "Mean square momentum obtained numerically.\n";
    histOfile << "Error (%):    " << std::abs(pMeanSqEmp - pMeanSqMB) / pMeanSqMB * 100. << '\n';
    histOfile << '\n';
    histOfile << "EkMB:         " << EkMB << '\t' << "Mean kinetic energy obtained analytically.\n";
    histOfile << "EkEmp:        " << EkEmp << '\t' << "Mean kinetic energy obtained numerically.\n";
    histOfile << "Error (%):    " << std::abs(EkEmp - EkMB) / EkMB * 100. << '\n';
    histOfile << '\n';

    histOfile.close();
}

void Stats::setInputFromArgon(const double *pAbsArgon, const usint &NArgon, const double &TArgon, const double &KArgon, const double &MArgon)
{
    N = NArgon;
    T = TArgon;
    k = KArgon;
    m = MArgon;

    delete[] pAbs;
    pAbs = new double[N];

    low = pAbsArgon[0];

    for (usint i = 0; i < N; i++)
    {
        pAbs[i] = pAbsArgon[i];

        if (pAbs[i] > up)
        {
            up = pAbs[i];
        }
        else if (pAbs[i] < low)
        {
            low = pAbs[i];
        }
    }

    low = floor(low);
    up = ceil(up);
    bins = 25; // static_cast<usint>(up - low);

    delete[] binRanges;
    binRanges = new double[bins + 1]();

    delete[] stars;
    stars = new std::string[bins]();
}