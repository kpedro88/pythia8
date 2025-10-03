// main369.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: top; toponium; openmp;

// Plot ttbar system properties contrasting six different scenarios:
// 0 : pure Born baseline.
// 1 : Coulomb enhancement.
// 2 : Green's function purely for E > 0.
// 3 : Green's function, only for the E < 0 part mirrored to E > 0.
//     Scenario 2 can be added to 3 to give a "complete" scenario.
// 4 : as 3, but in addition the top mass is reduced to give better <E>.
// 5 : full Green's for E > and E < 0, the "best bet" scenario.
// Plots all scenarios in one frame, which may be too much,
// but also shows how to select and combine as desired.
// Runs in parallel if Pythia configured with --with-openmp.
// The same analysis is also provided with main370 but using PythiaParallel.

// Based on threshold factors according to
// V. Fadin,  V. Khoze and T. Sjostrand, Z. Phys. C48 (1990) 613.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

// OpenMP includes and user-specifiable number of threads requested.
// (Will be capped by system maximum. Use -1 to use all available.)
#ifdef OPENMP
#include <omp.h>
const int nThreadsRequested = -1;
#endif

//============================================================================

int main() {

  // Restrict to threshold region only or not.
  bool thresholdOnly = true;

  // Full event generation or only cross-section-oriented studies.
  bool fullEvents = false;

  // Number of events. Fewer in scenarios 3 and 4, since only mirrored part.
  int nEvent       = 1000000;
  int nEventMirror = 0.1 * nEvent;

  // LHC collision energy.
  double eCM = 13000.;

  // Example of possible setup values. Reduced top mass in scenario 4.
  double mt            = 172.5;
  double mtLow         = mt - 3.8;
  double thrWidth      = 10.;
  int    alphasOrder   = 2;
  double alphasValue   = 0.118;
  double ggSingletFrac = 2. / 7.;
  double qqSingletFrac = 0.;

  // Histograms. (Last slot used for E < 0 part of scenario 5.)
  bool separateFrames = true;
  Hist mThrOrig[7], mThrOrigZoom[7],  mHatLow[7];
  for (int scenario = 0; scenario < 7; ++scenario) {
    mThrOrig[scenario].book("original threshold", 100, -20., 80.);
    mThrOrigZoom[scenario].book("original threshold, zoom", 100, -10., 30.);
    mHatLow[scenario].book("ttbar invariant mass, low", 100, 300., 400.);
  }

  // Compare the six scenarios.
  // Automatically run in parallel if Pythia configured with flag
  // --with-openmp.
#ifdef OPENMP
  int nThreadsSys = omp_get_max_threads();
  int nThreadsMax = min(6, nThreadsSys);
  int nThreads    = nThreadsRequested;
  if (nThreads <= 0 || nThreads > nThreadsMax) nThreads = nThreadsMax;
  cout << " OMP parralelisation using "<<nThreads << " threads."
       << " System maximum = "<<nThreadsSys<<"." << endl;
#pragma omp parallel for num_threads(nThreads)
#endif
  for (int scenario = 0; scenario < 6; ++scenario) {

    // Generator.
    Pythia pythia;

    // Feed in desired values.
    int topModel = (scenario < 4) ? scenario : scenario - 1;
    pythia.settings.mode("TopThreshold:model", topModel);
    double mtNow = (scenario == 4) ? mtLow :  mt;
    pythia.particleData.m0(6, mtNow);
    pythia.settings.parm("TopThreshold:width", thrWidth);
    pythia.settings.mode("TopThreshold:alphasOrder", alphasOrder);
    pythia.settings.parm("TopThreshold:alphasValue", alphasValue);
    pythia.settings.parm("TopThreshold:ggSingletFrac", ggSingletFrac);
    pythia.settings.parm("TopThreshold:qqSingletFrac", qqSingletFrac);

    // Process selection and collision energy.
    pythia.readString("Top:gg2ttbar = on");
    pythia.readString("Top:qqbar2ttbar = on");
    pythia.settings.parm("Beams:eCM", eCM);

    // Restrict to threshold region.
    if (thresholdOnly) {
      pythia.readString("PhaseSpace:mHatMin = 300.");
      pythia.readString("PhaseSpace:mHatMax = 400.");
    }
    else pythia.readString("PhaseSpace:mHatMin = 200.");

    // Switch off (most) code parts not relevant here. Reduce printout.
    if (!fullEvents) {
      pythia.readString("PartonLevel:ISR = off");
      pythia.readString("PartonLevel:FSR = off");
      pythia.readString("PartonLevel:MPI = off");
      pythia.readString("HadronLevel:all = off");
    }
    pythia.readString("Next:numberCount = 100000");

    // If Pythia fails to initialize, abort this scenario.
    // Use critical to avoid thread clashes when initializing in parallel.
    bool ok = true;
#ifdef OPENMP
#pragma omp critical
#endif
    {
      cout << " =============================================================="
        "=========\n Initializing PYTHIA for scenario " << scenario
           << " with topModel " << topModel << " and top mass " << mtNow
           << " GeV" << endl;
      ok = pythia.init();
      cout << " Completed initialization of PYTHIA for scenario " << scenario
           << "\n ============================================================"
        "===========\n";
    }
    if (!ok) continue;

    // Begin event loop.
    int nBelow = 0, nBelowgg = 0, nBelowqq = 0, nAbovegg = 0, nAboveqq = 0;
    double eBelow = 0.;
    int nEventNow = (scenario == 3 || scenario == 4) ? nEventMirror : nEvent;
    for (int iEvent = 0; iEvent < nEventNow; ++iEvent) {

      // Generate events.
      if (!pythia.next()) continue;

      // Original threshold value and other properties.
      double eThr = pythia.info.toponiumE;
      if (scenario == 4) eThr -= 2. * (mt - mtLow);
      double mHat = pythia.info.mHat();

      // Histogram threshold value and ttbar invariant mass.
      // Use critical to avoid thread clashes when filling histogram
      // in parallel.
#ifdef OPENMP
#pragma omp critical
#endif
      {
        mThrOrig[scenario].fill( eThr);
        mThrOrigZoom[scenario].fill( eThr);
        mHatLow[scenario].fill( mHat);
        if (scenario == 5 && eThr < 0.) {
          mThrOrig[6].fill( eThr);
          mThrOrigZoom[6].fill( eThr);
          mHatLow[6].fill( mHat);
        }
      }

      // Statistics for below threshold cross sections.
      if (scenario == 3 || scenario == 4 || (scenario == 5 && eThr < 0.)) {
        ++nBelow;
        if (pythia.info.code() == 601) ++nBelowgg;
        else ++nBelowqq;
        eBelow += eThr;
      } else {
        if (pythia.info.code() == 601) ++nAbovegg;
        else ++nAboveqq;
      }

      // End of event loop.
    }

    // Normalization and statistics for this scenario.
    // Use critical to avoid thread clashes.
#ifdef OPENMP
#pragma omp critical
#endif
    {
      cout << "\n ===== End-of-run statistics for scenario " << scenario
           << " with topModel " << topModel << " and top mass " << mtNow
           << " GeV =====" << endl;

      // Normalize histogram to cross section, in pb/GeV.
      double sigmapb = 1e9 * pythia.info.sigmaGen() / nEventNow;
      mThrOrig[scenario]     *= sigmapb;
      mThrOrigZoom[scenario] *= 2.5 * sigmapb;
      mHatLow[scenario]      *= sigmapb;
      if (scenario == 5) {
        mThrOrig[6]          *= sigmapb;
        mThrOrigZoom[6]      *= 2.5 * sigmapb;
        mHatLow[6]           *= sigmapb;
      }

      // Statistics and cross sections.
      pythia.stat();
      double sigAbv   = sigmapb * (nEventNow - nBelow);
      double sigAbvgg = sigmapb * nAbovegg;
      double sigAbvqq = sigmapb * nAboveqq;
      double sigBel   = sigmapb * nBelow;
      double sigBelgg = sigmapb * nBelowgg;
      double sigBelqq = sigmapb * nBelowqq;
      cout << "\n sigma above threshold = " << fixed << setprecision(3)
           << sigAbv << " pb" << endl << "   whereof gg = " << sigAbvgg
           << " and qq = "<< sigAbvqq << endl << " sigma below threshold = "
           << sigBel << " pb" << endl << "   whereof gg = " << sigBelgg
           << " and qq = " << sigBelqq << endl;
      eBelow /= max(1, nBelow);
      cout << " average energy for below-threshold part = " << eBelow
           << " GeV" << endl;
    }
    // End loop over threshold cases.
  }

  // No more parallelism beyond this point.
#ifdef OPENMP
#pragma omp barrier
#endif

  // Histograms with pyplot. Common notation.
  HistPlot hpl("plot369");
  string symbol[6] = {"h,black", "h,olive", "h,blue", "h,orange", "h,magenta",
                      "h,red"};
  string legend[6] = {"Born", "Coulomb", "Green's, $E > 0$",
      "Green's, $E < 0$ mirrored", "Green's, $E < 0$ mirrored + shifted",
      "Green's, all $E$"};

  // The formal E_thr value above/below threshold, key for reweighting.
  hpl.frame("fig369", "Energy above or below threshold",
    "$E$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  for (int scenario = 0; scenario < 6; ++scenario)
    hpl.add( mThrOrig[scenario], symbol[scenario], legend[scenario]);
  hpl.plot();
  hpl.frame("", "Energy above or below threshold",
    "$E$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  for (int scenario = 0; scenario < 6; ++scenario)
    hpl.add( mThrOrigZoom[scenario], symbol[scenario], legend[scenario]);
  hpl.plot();

  // Invariant mass of the ttbar pair.
  hpl.frame("", "invariant mass of ttbar pair, near threshold",
    "$m(t+tbar) (GeV)$", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  for (int scenario = 0; scenario < 6; ++scenario)
    hpl.add( mHatLow[scenario], symbol[scenario], legend[scenario]);
  hpl.plot();

  // Zooming in on three specific threshold scenarios (plots for talk).
  string legMirror = (separateFrames) ? "fig369mirror" : "";
  hpl.frame(legMirror, "", "$E$ (GeV)",
    "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)", 4.8, 4.8);
  hpl.add( mThrOrigZoom[0], "h,black", "Born cross section");
  hpl.add( mThrOrigZoom[2], "h,blue", "Green's above threshold");
  hpl.add( mThrOrigZoom[3], "h,magenta", "Green's below threshold mirrored");
  hpl.add( mThrOrigZoom[2] + mThrOrigZoom[3], "h,red", "Green's in total");
  hpl.plot(-10., 30., 0., 3.5);
  string legShift = (separateFrames) ? "fig369shift" : "";
  hpl.frame(legShift, "", "$E$ (GeV)",
    "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)", 4.8, 4.8);
  hpl.add( mThrOrigZoom[0], "h,black", "Born cross section");
  hpl.add( mThrOrigZoom[2], "h,blue", "Green's above threshold");
  hpl.add( mThrOrigZoom[4], "h,magenta",
    "Green's below threshold shifted + mirrored");
  hpl.add( mThrOrigZoom[2] + mThrOrigZoom[4], "h,red", "Green's in total");
  hpl.plot(-10., 30., 0., 3.5);
  string legFull = (separateFrames) ? "fig369full" : "";
  hpl.frame(legFull, "", "$E$ (GeV)",
    "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)", 4.8, 4.8);
  hpl.add( mThrOrigZoom[0], "h,black", "Born cross section");
  hpl.add( mThrOrigZoom[5] - mThrOrigZoom[6], "h,blue",
    "Green's above threshold");
  hpl.add( mThrOrigZoom[6], "h,magenta", "Green's below threshold");
  hpl.plot(-10., 30., 0., 3.5);

  // Done.
  return 0;
}
