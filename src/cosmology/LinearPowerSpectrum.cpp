// $Id: LinearPowerSpectrum.cpp,v 1.4 2006-07-07 18:26:45 garyb Exp $
// Default routines for calculating desired quantities from a linear PS
#include "PowerSpectra.h"
#include "Simpson.h"
#include "Std.h"
#include <cmath>

using namespace cosmology;

// Integrands to match Rob's "wint" integrals.
// Do integrals over ln K here
class LPSIntegrand {
public:
  enum FunctionType {SIGSQ, D1, D2};
  double operator()(double lnK) const {
    double kR2 = exp(2*lnK)*Rsq;
    double Dsq_dlnk = lps.DeltaSq(lnK);
    switch (ftype) {
    case SIGSQ:
      return Dsq_dlnk * exp(-kR2);
    case D1:
      return Dsq_dlnk * exp(-kR2) * 2. * kR2;
    case D2:
      return Dsq_dlnk * exp(-kR2) * 4. * kR2 * (1.-kR2);
    default:
      throw MyException("Bad ftype in LPSIntegrand");
    }
  }
  LPSIntegrand(const LinearPowerSpectrum& lps_, FunctionType f_,
	       double k0): lps(lps_), ftype(f_), Rsq(1./(k0*k0)) {}
private:
  const LinearPowerSpectrum& lps;
  FunctionType ftype;
  double Rsq;
};

class TophatIntegrand {
public:
  double operator()(double k) const {
    double kR = k*R;
    double W = 3*(sin(kR)-kR*cos(kR))*pow(kR,-3.);
    double Dsq_dlnk = lps.DeltaSq(log(k));
    return Dsq_dlnk * W * W / k;
  }
  TophatIntegrand(const LinearPowerSpectrum& lps_, double R_): lps(lps_),
							       R(R_) {}
private:
  const LinearPowerSpectrum& lps;
  double R;
};

double
LinearPowerSpectrum::sigX(double x) const {
  TophatIntegrand th(*this, x);
  return sqrt(Simp1d(th, 1e-4, 20./x));
}

void
LinearPowerSpectrum::integrals(double k0, 
			       double& sigSq, double& neff, double& C) const {
  // Bounds of integrations
  double lnk0 = MIN(log(k0)-3., log(1e-4));  
  double lnk1 = MAX(log(k0)+3., log(1e+4));
  LPSIntegrand lp0((*this), LPSIntegrand::SIGSQ, k0);
  LPSIntegrand lp1((*this), LPSIntegrand::D1, k0);
  LPSIntegrand lp2((*this), LPSIntegrand::D2, k0);
  sigSq = Simp1d(lp0,lnk0, lnk1, 1.e-5);
  double d1 = Simp1d(lp1,lnk0, lnk1, 1.e-5);
  double d2 = Simp1d(lp2,lnk0, lnk1, 1.e-5);
  neff = d1 / sigSq - 3.;
  C = d1*d1/(sigSq*sigSq) + d2/sigSq;
}

void
LinearPowerSpectrum::extendTables(double lnK) const {
  const double tableStep=0.05;
  const double tableSafetyMargin=0.31;
  int evalCount=0;
  const int maxEvals=100;

  if (CTable.size()==0) {
    // Make initial entry in table right here:
    double sigSq, neff, C;
    integrals(exp(lnK), sigSq, neff, C);
    sigSqTable.addEntry(lnK, log(sigSq));
    neffTable.addEntry(lnK, neff);
    CTable.addEntry(lnK, C);
  }
  // March table down in lnK
  double tableMin = CTable.argMin();
  while (lnK < tableMin + tableSafetyMargin) {
    if (++evalCount > maxEvals) throw MyException("Too many extendTables steps");
    tableMin -= tableStep;
    double sigSq, neff, C;
    integrals(exp(tableMin), sigSq, neff, C);
    sigSqTable.addEntry(tableMin, log(sigSq));
    neffTable.addEntry(tableMin, neff);
    CTable.addEntry(tableMin, C);
  }
  double tableMax = CTable.argMax();
  while (lnK > tableMax - tableSafetyMargin) {
    if (++evalCount > maxEvals) throw MyException("Too many extendTables steps");
    tableMax += tableStep;
    double sigSq, neff, C;
    integrals(exp(tableMax), sigSq, neff, C);
    sigSqTable.addEntry(tableMax, log(sigSq));
    neffTable.addEntry(tableMax, neff);
    CTable.addEntry(tableMax, C);
  }
}

LinearPowerSpectrum::LinearPowerSpectrum(): sigSqTable(gtable::Table<>::spline),
					    neffTable(gtable::Table<>::spline),
					    CTable(gtable::Table<>::spline) {}

  
double
LinearPowerSpectrum::sigSq(double lnK) const {
  extendTables(lnK);
  return exp(sigSqTable(lnK));
}
double
LinearPowerSpectrum::neff(double lnK) const {
  extendTables(lnK);
  return neffTable(lnK);
}
double
LinearPowerSpectrum::C(double lnK) const {
  extendTables(lnK);
  return CTable(lnK);
}

double
LinearPowerSpectrum::lnkNL(double G) const {
  // Find the scale at which spectrum is nonlinear when growth factor
  // rises to G, i.e. sigSq(kNL)*G*G=1.;
  const double tolerance=1.e-5; // Tolerance (in ln(sigSq)) for solution.
  const double lnKStart = 0.;	// Start the search at this value
  const double lnKStep = 1.;	// initial lnK step size in soln search

  int evalCount=0;
  const int maxEvals=50;

  double lG2=2*log(G);

  /*cerr << "Looking for kNL at G " << G << endl;*/
  // Note that sigSq has to be an INCREASING function of k, so search is easier

  // First march in some direction at lnKStep until root is bracketed
  double lnKupper, lnKlower, lnKmid;
  double dsupper, dslower, dsmid;
  double testds = log(sigSq(lnKStart))+lG2;
  if (fabs(testds)<tolerance) return lnKStart;

  lnKlower = lnKupper = lnKStart;
  dslower = dsupper = testds;

  if (testds<0.) {
    while (dsupper < 0.) {
      /*cerr << "Going up: lnKu " << lnKupper << " ds " << dsupper << endl;*/
      lnKlower = lnKupper;
      dslower = dsupper;
      lnKupper = lnKlower + lnKStep;
      dsupper =  log(sigSq(lnKupper))+lG2;
      if (++evalCount > maxEvals) throw MyException("Too many kNL steps");
    }
  } else {
    while (dslower > 0.) {
      /*cerr << "Going down: lnKl " << lnKlower << " ds " << dslower << endl;*/
      lnKupper = lnKlower;
      dsupper = dslower;
      lnKlower = lnKupper - lnKStep;
      dslower = log(sigSq(lnKlower))+lG2;
      if (++evalCount > maxEvals) throw MyException("Too many kNL steps");
    }
  }

  if (fabs(dsupper)<tolerance) return lnKupper;
  if (fabs(dslower)<tolerance) return lnKlower;
  
  // Now proceed with bisection
  while (true) {
    lnKmid = 0.5*(lnKupper + lnKlower);
    dsmid = log(sigSq(lnKmid))+lG2;
    /*cerr << "Bisecting: lnKm " << lnKmid << " ds " << dsmid << endl;*/
    if (fabs(dsmid)<tolerance)
      return lnKmid;
    else if (dsmid<0.)
      lnKlower = lnKmid;
    else if (dsmid>0.)
      lnKupper = lnKmid;
    if (++evalCount > maxEvals) throw MyException("Too many kNL steps");
  }
}
