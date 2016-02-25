// $Id: GrowthFunction.cpp,v 1.7 2006-07-13 19:26:09 garyb Exp $

#include "GrowthFunction.h"
#include "odeint.h"
#include "AstronomicalConstants.h"

using namespace cosmology;
using namespace ode;

class GFDerivsH: public ODEFunction {
public:
  GFDerivsH(const Cosmology& c_): c(c_) {}
  int order() const {return 2;}
  DVector operator()(double lna, const DVector& y) const {
    double z = exp(-lna) - 1.;
    DVector ret(2);
    ret[0] = y[1];
    ret[1] = 1.5 * c.OmegaM() * pow(c.H(z),-2.) * exp(-3.*lna) +
      (c.q(z)-1.) * y[1] - y[1]*y[1];
    return ret;
  }
private:
  const Cosmology& c;
};

/*
class DofA;

class GFDerivsChi: public ODEFunction {
public:
  GFDerivsChi(const DofA& c_): c(c_) {}
  int order() const {return 2;}
  DVector operator()(double lna, const DVector& y) const {
    DVector ret(2);
    ret[0] = y[1];
    double chi = c.chi(lna);
    double d1 = c.d1(lna);
    double d2 = c.d2(lna);
    ret[1] = 1.5 * omegaM * d1 * d1 *exp(-lna)/(chi*chi) - 
      (d1 + d2/d1 -1.) * y[1] - y[1]*y[1];
    return ret;
  }
private:
  const DofA& c;
  double omegaM;
};
*/
void
GrowthFunction::buildTables(double aMax) {
  lnGTable.clear();
  dlnGdlnATable.clear();

  // start integration at z=10 if w1 is in place ???
  const double lnaStart= c.usingWa() ? -log(1.+RecombinationRedshift) : -1.;
  //const double lnaStart=log(1e-7);
  const double lna = log(aMax);
  const double lnaStep = 0.05;	// suggested step

  //  const double eps=1e-6;
  const double eps=1e-6;
  const double minstep=1e-5;

  int numberOfSteps = static_cast<int> (ceil(lna-lnaStart)/lnaStep);
  DVector xp(numberOfSteps);
  DMatrix yp(2, numberOfSteps);

  DVector initialConditions(2);
  initialConditions[0] = 0.;	// Normalize growth function to 1 at rec
  initialConditions[1] = 1.;	// dlogG/dlog a=1. at matter-dominated

  GFDerivsH gfd(c);
  int nok, nbad, kount=0;
  odeint(initialConditions,    
	 xp,
	 yp,
	 lnaStep,
	 kount,
	 lnaStart,
	 lna,
	 eps,
	 lnaStep,
	 minstep,
	 nok,
	 nbad,
	 gfd);

  // Now stuff all the saved ODE points into the two tables
  for (int i=0; i<kount; i++) {
    lnGTable.addEntry(xp[i], yp(0,i));
    dlnGdlnATable.addEntry(xp[i], yp(1,i));
  }
  lnGRecombination = lnGTable(log(1+RecombinationRedshift));
};
double
GrowthFunction::lnG(double lna) const {
  return lnGTable(lna);
}
double
GrowthFunction::dlnGdlnA(double lna) const {
  return dlnGdlnATable(lna);
}

double
GrowthFunction::valueAtZ(double z) const {
  const double zRecombination=1088.;
  const double lnaRec=-log(1.+zRecombination);
  const double lna = -log(1.+z);
  DVector initialConditions(2);
  initialConditions[0] = 0.;	// Normalize growth function to 1 at rec
  initialConditions[1] = 1.;	// dlogG/dlog a=1. at matter-dominated

  GFDerivsH gfd(c);
  DVector xp(0);
  DMatrix yp(0,0);
  int nok, nbad, kount=0;
  odeint(initialConditions,
	 xp,
	 yp,
	 1.,
	 kount,
	 lnaRec,
	 lna,
	 1.e-6,
	 0.2,
	 0.01,
	 nok,
	 nbad,
	 gfd);
  /*
  cerr << "nok/bad " << nok << "/" << nbad
       << " dlnG/dlna " << initialConditions[1]
       << " G/a " << exp(initialConditions[0]) * exp(lnaRec-lna)
       << endl;*/
  return initialConditions[0];
}
