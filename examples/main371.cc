// main371.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: top; toponium

// Study ttbar production above and below threshold,
// starting out from the equations in
// V. Fadin,  V. Khoze and T. Sjostrand, Z. Phys. C48 (1990) 613,
// with some further assumptions in the below-threshold simulation.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//============================================================================

int main() {

  // Restrict to threshold region only or not.
  bool thresholdOnly = true;

  // Full event generation or only cross-section-oriented studies.
  bool fullEvents = false;

  // Number of events.
  int nEvent = 1000000;

  // LHC collision energy.
  double eCM = 13000.;

  // Example of possible setup values.
  // topModel = 0 is Born, = 1 is simple Coulomb,
  // = 2 is Green's above threshold, = 3 adds below-threshold mirrored,
  // = 4 includes above-and-below-threshold parts (most plausible scenario).
  int    topModel       = 4;
  double mt             = 172.5;
  double thresholdWidth = 10.;
  int    alphasOrder    = 2;
  double alphasValue    = 0.118;
  double ggSingletFrac  = 2./7.;
  double qqSingletFrac  = 0.;

  // Generator.
  Pythia pythia;

  // Feed in desired values.
  pythia.settings.mode("TopThreshold:model", topModel);
  pythia.particleData.m0(6, mt);
  pythia.settings.parm("TopThreshold:width", thresholdWidth);
  pythia.settings.mode("TopThreshold:alphasOrder", alphasOrder);
  pythia.settings.parm("TopThreshold:alphasValue", alphasValue);
  pythia.settings.parm("TopThreshold:ggSingletFrac", ggSingletFrac);
  pythia.settings.parm("TopThreshold:qqSingletFrac", qqSingletFrac);

  // Process selection. LHC at 13 TeV.
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

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Histograms.
  Hist mHatOri("original threshold", 100, -20., 80.);
  Hist mHatOriLow("original threshold", 100, -10., 30.);
  Hist mHatThr("corrected threshold", 100, -20., 80.);
  Hist mHatLow("ttbar invariant mass, low", 100, 300., 400.);
  Hist mHatAll("ttbar invariant mass, all", 100, 200., 1200.);
  Hist mTopMax("max(m_t, m_tbar) mass distribution", 100, 100., 200.);
  Hist mTopMin("min(m_t, m_tbar) mass distribution", 100, 100., 200.);
  Hist betaPair("beta factor of pair", 100, 0., 1.);
  Hist mHatOriBT("original threshold, below", 100, -20., 80.);
  Hist mHatOriBTLow("original threshold, below", 100, -10., 30.);
  Hist mHatThrBT("corrected threshold, below", 100, -10., 30.);
  Hist mHatLowBT("ttbar invariant mass, low, below", 100, 300., 400.);
  Hist mHatAllBT("ttbar invariant mass, all, below", 100, 200., 1200.);
  Hist mTopMaxBT("max(m_t, m_tbar) mass distribution, below", 100, 100., 200.);
  Hist mTopMinBT("min(m_t, m_tbar) mass distribution, below", 100, 100., 200.);
  Hist betaPairBT("beta factor of pair, below", 100, 0., 1.);
  Hist mShiftBT("shift in top mass, after - before", 100, -45., 5.);

  // Begin event loop.
  int nBelow = 0, nBelowgg = 0, nBelowqq = 0, nAbovegg = 0, nAboveqq = 0;
  double eBelow = 0.;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events.
    if (!pythia.next()) continue;

    // Original threshold value and other properties.
    double eThr = pythia.info.toponiumE;
    double mt1T = pythia.info.toponiumm3;
    double mt2T = pythia.info.toponiumm4;
    double mHat = pythia.info.mHat();
    double mt1  = pythia.info.m3Hat();
    double mt2  = pythia.info.m4Hat();
    double rat1 = pow2( mt1 / mHat);
    double rat2 = pow2( mt2 / mHat);
    double beta = sqrtpos( pow2(1 - rat1 - rat2) - 4. * rat1 * rat2);

    // Histogram for all events.
    mHatOri.fill( eThr);
    mHatOriLow.fill( eThr);
    mHatThr.fill( mHat - mt1 - mt2);
    mHatLow.fill( mHat);
    mHatAll.fill( mHat);
    mTopMax.fill( max(mt1,mt2) );
    mTopMin.fill( min(mt1,mt2) );
    betaPair.fill( beta);

    // Histogram for below-threshold events.
    if (eThr < 0.) {
      ++nBelow;
      if (pythia.info.code() == 601) ++nBelowgg;
      else ++nBelowqq;
      eBelow += eThr;
      mHatOriBT.fill( eThr);
      mHatOriBTLow.fill( eThr);
      mHatThrBT.fill( mHat - mt1 - mt2);
      mHatLowBT.fill( mHat);
      mHatAllBT.fill( mHat);
      mTopMaxBT.fill( max(mt1,mt2) );
      mTopMinBT.fill( min(mt1,mt2) );
      betaPairBT.fill( beta);
      mShiftBT.fill( mt1 - mt1T);
      mShiftBT.fill( mt2 - mt2T);
    } else {
      if (pythia.info.code() == 601) ++nAbovegg;
      else ++nAboveqq;
    }

  // End of event loop.
  }

  // Normalize histogram to cross section, in pb/GeV.
  double sigmapb = 1e9 * pythia.info.sigmaGen() / nEvent;
  mHatOri      *= sigmapb;
  mHatOriLow   *= 2.5 * sigmapb;
  mHatThr      *= sigmapb;
  mHatLow      *= sigmapb;
  mHatAll      *= 0.1 * sigmapb;
  mTopMax      *= sigmapb;
  mTopMin      *= sigmapb;
  betaPair     *= 100. * sigmapb;
  mHatOriBT    *= sigmapb;
  mHatOriBTLow *= 2.5 * sigmapb;
  mHatThrBT    *= sigmapb;
  mHatLowBT    *= sigmapb;
  mHatAllBT    *= 0.1 * sigmapb;
  mTopMaxBT    *= sigmapb;
  mTopMinBT    *= sigmapb;
  betaPairBT   *= 100. * sigmapb;
  mShiftBT     *= 2. * sigmapb;

  // Statistics and cross sections
  pythia.stat();
  double sigAbv   = sigmapb * (nEvent - nBelow);
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

  // Histograms with pyplot.
  HistPlot hpl("plot371");

  // Spectra relative to event-by-event threshold, including below.
  hpl.frame("fig371", "Energy above/below threshold, original tops",
    "$E$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}E$ (pb/GeV)");
  hpl.add( mHatOriBT, "h,red", "below-threshold part");
  hpl.add( mHatOri - mHatOriBT, "h,black", "above-threshold part");
  //hpl.plot(-20., 80., 0., 5.);
  hpl.plot();

  // Ditto, but with new t/tbar masses for below-threshold part.
  hpl.frame("", "Energy above threshold, new top masses where necessary",
    "$E'$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}E'$ (pb/GeV)");
  hpl.add( mHatThrBT, "h,red", "below-threshold part");
  hpl.add( mHatThr, "h,black", "all");
  //hpl.plot(-10., 40., 0., 5.);
  hpl.plot();

  // Pair mass spectra close to threshold.
  hpl.frame("", "Invariant mass of ttbar pair, near threshold",
    "$m(t+tbar) (GeV)$", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  hpl.add( mHatLowBT, "h,red", "below-threshold part");
  hpl.add( mHatLow, "h,black", "all");
  //hpl.plot(300., 400., 0., 3.5);
  hpl.plot();

  // Pair mass spectra over full range.
  hpl.frame("", "Invariant mass of ttbar pair, full range",
    "$m(t+tbar) (GeV)$", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  hpl.add( mHatAllBT, "h,red", "below-threshold part");
  hpl.add( mHatAll, "h,black", "all");
  hpl.plot();

  // Top/antitop mass spectra.
  hpl.frame("", "Larger of top and antitop masses",
    "max($m(t), m(tbar)$) (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ "
    "(pb/GeV)");
  hpl.add( mTopMaxBT, "h,red", "below-threshold part");
  hpl.add( mTopMax, "h,black", "all");
  hpl.plot(100., 200., 0.01, 100., true, false);

  hpl.frame("", "Smaller of top and antitop masses",
    "min($m(t), m(tbar)$) (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ "
    "(pb/GeV)");
  hpl.add( mTopMinBT, "h,red", "below-threshold part");
  hpl.add( mTopMin, "h,black", "all");
  hpl.plot(100., 200., 0.01, 100., true, false);

  // Top/antitop separation.
  hpl.frame("", "Beta factor in t-tbar pair",
    "$\\beta$", "$\\mathrm{d}\\sigma/\\mathrm{d}\\beta$ (pb)");
  hpl.add( betaPairBT, "h,red", "below-threshold part");
  hpl.add( betaPair, "h,black", "all");
  hpl.plot();

  // Redefinition change of top masses in below-threshold pairs.
  hpl.frame("", "Shift of below-threshold top masses", "$\\Delta m$ (GeV)",
    "$\\mathrm{d}\\sigma/\\mathrm{d}(\\Delta m)$ (pb/GeV)");
  hpl.add( mShiftBT, "h,red", "below-threshold part");
  hpl.plot();

  // For inclusion in talk.
  hpl.frame("fig371mass", "", "$\\hat{m}$ (GeV)",
    "$\\mathrm{d}\\sigma/\\mathrm{d}\\hat{m}$ (pb/GeV)", 4.8, 4.8);
  hpl.add( mHatLow, "h,black", "all");
  hpl.add( mHatLowBT, "h,red", "below threshold");
  hpl.plot();
  hpl.frame("fig371BreitWigners", "", "$m$ (GeV)",
    "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)", 4.8, 4.8);
  hpl.add( mTopMax, "h,black", "max($m_{t1},m_{t2}$) all");
  hpl.add( mTopMin, "h,blue", "min($m_{t1},m_{t2}$) all");
  hpl.add( mTopMaxBT, "h,red", "max($m_{t1},m_{t2}$) below threshold");
  hpl.add( mTopMinBT, "h,magenta", "min($m_{t1},m_{t2}$) below threshold");
  hpl.plot(140., 200., 0.01, 1000., true, false);

  // Done.
  return 0;
}
