#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>

#include <Rcpp.h>

#include "utils.h"
#include "utils-spectra.h"

namespace {

void specApply(Rcpp::List spectra, double rtMin, double rtMax, double mzMin, double mzMax,
               std::function<void(double, const Rcpp::NumericVector &, const Rcpp::NumericVector &)> func)
{
    const Rcpp::DataFrame specHeader = Rcpp::as<Rcpp::DataFrame>(spectra["header"]);
    const Rcpp::NumericVector hdRetTimes = specHeader["retentionTime"];
    const Rcpp::NumericVector hdMSLevels = specHeader["msLevel"];
    const Rcpp::NumericVector hdSeqNums = specHeader["seqNum"];
    const Rcpp::List specList = spectra["spectra"];
    const int specCount = hdRetTimes.size();
        
    for (int spi=0; spi<specCount; ++spi)
    {
        const double specRt = hdRetTimes[spi];
        
        if (!numberWithin(specRt, rtMin, rtMax, 1E-4) || hdMSLevels[spi] != 1)
            continue;
        
        const Rcpp::DataFrame peaklist = Rcpp::as<Rcpp::DataFrame>(specList[hdSeqNums[spi] - 1]); // -1: R 1 based index
        func(specRt, peaklist["mz"], peaklist["intensity"]);
    }
}

struct BinnedSpectrum
{
    std::vector<double> mzs, intsLeft, intsRight;
};
    

}


// [[Rcpp::export]]
Rcpp::NumericVector loadEICIntensities(Rcpp::List spectra, Rcpp::DataFrame featList, Rcpp::NumericVector rtWindow)
{
    const Rcpp::NumericVector rets = featList["ret"];
    const Rcpp::NumericVector mzmins = featList["mzmin"];
    const Rcpp::NumericVector mzmaxs = featList["mzmax"];
    
    const double rtWin = Rcpp::as<double>(rtWindow);
    const int ftCount = rets.length();
    
    Rcpp::NumericVector intensities(ftCount);
    for (int fti=0; fti<ftCount; ++fti)
    {
        const double featRt = rets[fti];
        const double mzMin = mzmins[fti], mzMax = mzmaxs[fti];
        double closestRTDiff = -1, closestInt = 0;
        
        specApply(spectra, featRt - rtWin, featRt + rtWin, mzMin, mzMax,
                  [&](double specRt, const Rcpp::NumericVector &peakMZs, const Rcpp::NumericVector &peakInts)
                  {
                      const double rtdiff = fabs(featRt - specRt);
                      if (closestRTDiff == -1)
                          closestRTDiff = rtdiff;
                      else if (closestRTDiff < rtdiff)
                          return;
                      
                      const int pCount = peakMZs.length();
                      
                      double totInt = 0;
                      for (int pi=0; pi<pCount; ++pi)
                      {
                          if (numberWithin(peakMZs[pi], mzMin, mzMax, 1E-8))
                              totInt += peakInts[pi];
                      }
                      
                      closestRTDiff = rtdiff;
                      closestInt = totInt;
                  });
        
        intensities[fti] = closestInt;
    }
    
    return intensities;
}

// [[Rcpp::export]]
Rcpp::List loadEICs(Rcpp::List spectra, Rcpp::List rtRanges, Rcpp::List mzRanges)
{
    const int EICCount = rtRanges.length();
    Rcpp::List ret(EICCount);
    for (int eici=0; eici<EICCount; ++eici)
    {
        const Rcpp::NumericVector rtr = rtRanges[eici], mzr = mzRanges[eici];
        const double mzMin = mzr[0], mzMax = mzr[1];
        
        std::vector<double> EICTimes;
        std::vector<double> EICIntensities;
        specApply(spectra, rtr[0], rtr[1], mzMin, mzMax,
                  [&](double specRt, const Rcpp::NumericVector &peakMZs, const Rcpp::NumericVector &peakInts)
                  {
                      const int pCount = peakMZs.length();
                      double totInt = 0;
                      for (int pi=0; pi<pCount; ++pi)
                      {
                          if (numberWithin(peakMZs[pi], mzMin, mzMax, 1E-8))
                              totInt += peakInts[pi];
                      }
                      
                      EICTimes.push_back(specRt);
                      EICIntensities.push_back(totInt);
                  });

        // compress data by removing any zero intensity datapoints that are inbetween two other zero intensitiy points.
        for (size_t ind=0; ind<(EICTimes.size()-2); )
        {
            if (EICIntensities[ind+2] != 0)
                ind += 3;
            else if (EICIntensities[ind+1] != 0)
                ind += 2;
            else if (EICIntensities[ind] != 0)
                ++ind;
            else // all zero
            {
                EICTimes.erase(EICTimes.begin() + (ind + 1));
                EICIntensities.erase(EICIntensities.begin() + (ind + 1));
            }
        }
        
        ret[eici] = Rcpp::DataFrame::create(Rcpp::Named("time") = EICTimes,
                                            Rcpp::Named("intensity") = EICIntensities);
    }
    
    return ret;
}

BinnedSpectrum binSpectra(const Spectrum &specLeft, Spectrum specRight,
                          const std::string &shift, double precDiff, double mzWindow)
{
    // assumptions: specs are ordered on mz
    
    if (shift == "precursor")
    {
        for (double &m : specRight.mzs)
            m -= precDiff;
        // NOTE: negative m/z values will be skipped below
    }
    else if (shift == "both") // UNDONE: other name for "both"?
    {
        // first bin as normal (recursive call)
        const BinnedSpectrum binNone = binSpectra(specLeft, specRight, "none", precDiff, mzWindow);
        Spectrum specLeftUn, specRightUn;
        BinnedSpectrum binOverlap;
        for (size_t i=0; i<binNone.mzs.size(); ++i)
        {
            const double m = binNone.mzs[i];
            if (binNone.intsLeft[i] == 0)
            {
                specRightUn.mzs.push_back(m);
                specRightUn.intensities.push_back(binNone.intsRight[i]);
            }
            else if (binNone.intsRight[i] == 0)
            {
                specLeftUn.mzs.push_back(m);
                specLeftUn.intensities.push_back(binNone.intsLeft[i]);
            }
            else
            {
                binOverlap.mzs.push_back(m);
                binOverlap.intsLeft.push_back(binNone.intsLeft[i]);
                binOverlap.intsRight.push_back(binNone.intsRight[i]);
            }
        }
        
        // bin missing with shift
        BinnedSpectrum binShift = binSpectra(specLeftUn, specRightUn, "precursor", precDiff, mzWindow);
        
        // merge both: add missing from binNone
        binShift.mzs.insert(binShift.mzs.end(), binOverlap.mzs.begin(), binOverlap.mzs.end());
        binShift.intsLeft.insert(binShift.intsLeft.end(), binOverlap.intsLeft.begin(), binOverlap.intsLeft.end());
        binShift.intsRight.insert(binShift.intsRight.end(), binOverlap.intsRight.begin(), binOverlap.intsRight.end());
        
        // UNDONE: sort?
        
        return binShift;
    }
    
    BinnedSpectrum ret;
    std::vector<size_t> usedRightInds;
    size_t lastRightInd = 0;

    for (size_t i=0; i<specLeft.mzs.size(); ++i)
    {
        const double leftMZ = specLeft.mzs[i];
        double rightMZ = 0, rightInt;
        bool foundRight = false;
        
        while (lastRightInd < specRight.mzs.size())
        {
            const double rmz = specRight.mzs[lastRightInd], rmzmin = rmz - mzWindow, rmzmax = rmz + mzWindow;
            
            if (rmz <= 0) // may be (below) zero due to precursor shift
            {
                ++lastRightInd;
                continue;
            }
            
            if (leftMZ < rmzmin)
                break; // surpassed range for left
            
            if (leftMZ >= rmzmin && leftMZ <= rmzmax)
            {
                // overlap
                rightMZ = rmz; rightInt = specRight.intensities[lastRightInd];
                foundRight = true;
                usedRightInds.push_back(lastRightInd);
            }
            
            ++lastRightInd;
            
            if (foundRight)
                break; // done or surpassed range
        }
        
        if (foundRight)
        {
            ret.mzs.push_back((leftMZ + rightMZ) / 2.0);
            ret.intsRight.push_back(rightInt);
        }
        else
        {
            ret.mzs.push_back(leftMZ);
            ret.intsRight.push_back(0);
        }
        ret.intsLeft.push_back(specLeft.intensities[i]);
    }
    
    // add missing from right
    for (size_t j=0; j<specRight.mzs.size(); ++j)
    {
        if (std::find(usedRightInds.begin(), usedRightInds.end(), j) == usedRightInds.end() &&
            specRight.mzs[j] > 0)
        {
            ret.mzs.push_back(specRight.mzs[j]);
            ret.intsLeft.push_back(0);
            ret.intsRight.push_back(specRight.intensities[j]);
        }
    }
    
    // UNDONE: sort?
    
    return ret;
}

// UNDONE: remove
// [[Rcpp::export]]
Rcpp::DataFrame binSpecCPP(Rcpp::DataFrame sp1, Rcpp::DataFrame sp2, Rcpp::CharacterVector shift,
                           Rcpp::NumericVector mzWindow)
{
    Spectrum specLeft{ sp1["mz"], sp1["intensity"] };
    Spectrum specRight{ sp2["mz"], sp2["intensity"] };
    
    normalizeNums(specLeft.intensities); normalizeNums(specRight.intensities);
    
    // figure out precursor masses
    const std::vector<int> isPrecLeft = sp1["precursor"];
    const std::vector<int> isPrecRight = sp2["precursor"];
    
    double precDiff = 0.0;
    auto itl = std::find(isPrecLeft.begin(), isPrecLeft.end(), TRUE);
    auto itr = std::find(isPrecRight.begin(), isPrecRight.end(), TRUE);
    if (itl == isPrecLeft.end() || itr == isPrecRight.end())
    {
        // UNDONE
        // Rcpp::stop("Cannot shift spectra: one or both lack precursor ion!");
    }
    else
        precDiff = specRight.mzs[itr - isPrecRight.begin()] - specLeft.mzs[itl - isPrecLeft.begin()];
    
    BinnedSpectrum binnedSpec = binSpectra(specLeft, specRight, Rcpp::as<std::string>(shift), precDiff, Rcpp::as<double>(mzWindow));
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = binnedSpec.mzs,
                                   Rcpp::Named("intensity_1") = binnedSpec.intsLeft,
                                   Rcpp::Named("intensity_2") = binnedSpec.intsRight);
}

double doCalcSpecSimilarity(Spectrum sp1, Spectrum sp2, const std::string &method,
                            const std::string &shift, double precDiff,
                            double mzWeight, double intWeight, double mzWindow)
{
    normalizeNums(sp1.intensities); normalizeNums(sp2.intensities);
    
    BinnedSpectrum binnedSpec = binSpectra(sp1, sp2, shift, precDiff, mzWindow);
    
    // UNDONE: pearsons/spearman? needs sorting?
    if (method == "cosine")
    {
        std::vector<double> u, v;
        for (size_t i=0; i<binnedSpec.mzs.size(); ++i)
        {
            const double m = std::pow(binnedSpec.mzs[i], mzWeight);
            u.push_back(m * std::pow(binnedSpec.intsLeft[i], intWeight));
            v.push_back(m * std::pow(binnedSpec.intsRight[i], intWeight));
        }
        
        const double dp = std::inner_product(u.begin(), u.end(), v.begin(), 0.0);
        // sqrt(sum(u^2)) * sqrt(sum(v^2)))
        double divu = 0.0, divv = 0.0;
        for (size_t i=0; i<u.size(); ++i)
        {
            divu += (u[i] * u[i]);
            divv += (v[i] * v[i]);
        }
        
        if (divu == 0.0 || divv == 0.0) // lack of any overlap
            return 0.0;
        
        const double div = std::sqrt(divu) * std::sqrt(divv);
        
        return dp / div;
    }
    else if (method == "jaccard")
    {
        // binnedPL[intensity_1 != 0 & intensity_2 != 0, .N] / nrow(binnedPL)
        int both = 0;
        for (size_t i=0; i<binnedSpec.mzs.size(); ++i)
        {
            if (binnedSpec.intsLeft[i] != 0 && binnedSpec.intsRight[i] != 0)
                ++both;
        }
        return static_cast<double>(both) / static_cast<double>(binnedSpec.mzs.size());
    }
    
    return NA_REAL; // shouldn't be here
}

// [[Rcpp::export]]
Rcpp::NumericVector calcSpecSimilarity(Rcpp::DataFrame sp1, Rcpp::DataFrame sp2, Rcpp::CharacterVector method,
                                       Rcpp::CharacterVector shift, Rcpp::NumericVector precDiff,
                                       Rcpp::NumericVector mzWeight, Rcpp::NumericVector intWeight, Rcpp::NumericVector mzWindow)
{
    Spectrum specLeft{ sp1["mz"], sp1["intensity"] };
    Spectrum specRight{ sp2["mz"], sp2["intensity"] };
    
    return Rcpp::NumericVector::create(doCalcSpecSimilarity(Spectrum{ sp1["mz"], sp1["intensity"] }, Spectrum{ sp2["mz"], sp2["intensity"] },
                                                            Rcpp::as<std::string>(method), Rcpp::as<std::string>(shift), Rcpp::as<double>(precDiff),
                                                            Rcpp::as<double>(mzWeight), Rcpp::as<double>(intWeight), Rcpp::as<double>(mzWindow)));
}
