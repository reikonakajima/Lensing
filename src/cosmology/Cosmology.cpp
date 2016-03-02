// $Id: Cosmology.cpp,v 1.19 2006-08-11 13:43:41 garyb Exp $
// Cosmological distances and times, routines
#include "Cosmology.h"
#include "AstronomicalConstants.h"
#include <iostream>

using namespace cosmology;

// The "dark energy" contribution to rho is 
// rho = omega_X * (1+z)^(3(1+w)), and can put
// w = w + w' z if desired; 
// now w = w0 + w1(1-a) is the default (useWa=true), Linder's model.
// w = P/rho = -1 for cosmological constant
// For now I'm not going to deal with w'.

static double aRatio=1.05;	//Factor for tabulation intervals of a
static double Tolerance=1.e-5;	//fractional tolerance for each integration
static double minA=1e-10;	//smallest a to search for dominance

/// Change integrands!
/// add derivative integrals
/// and vector input constructor

// class CosmoFunction contains all the useful functions
void
Cosmology::CosmoFunction::setup() {
  Om = c->Om;
  Ol = c->Ol;
  w0 = c->w0;
  w1 = c->w1;
  Or = c->Or;
  Ok = 1 - Om - Ol - Or;
  useWa = c->useWa;
  if (useWa) 
    ww = 4-3*(1.+w0+w1);
  else
    ww = 4-3*(1.+w0-w1);
}

inline
double Cosmology::CosmoFunction::expfunc(double a) const {
  return useWa ? exp(-3.*w1*(1.-a)) : exp(3.*w1*(1./a-1.));
}

inline
double Cosmology::CosmoFunction::dLogExpfuncdA(double a) const {
  return useWa ? 3.*w1 : -3.*w1/(a*a);
}

double 
Cosmology::CosmoFunction::operator()(double a) const {
  // First calculate H^2 a^4 function from Friedman eqn:
  double hh = Or + a*(Om+a*Ok) + Ol*pow(a,ww)*expfunc(a);
  double qq;
  switch (ftype) {
  case Cosmology::H2:
    return hh;
  case Cosmology::Q:
    qq=a*(Om+2.*a*Ok) + Ol*pow(a,ww)*expfunc(a)*
      (ww + a*dLogExpfuncdA(a));
    return 1.-0.5*qq/hh;
  case Cosmology::H2Abs:
    // A version of H2 that does not allow negative Omegas to cancel.
    return fabs(Or) + a*(fabs(Om)+a*fabs(Ok)) + 
      fabs(Ol)*pow(a,ww)*expfunc(a);
  case Cosmology::DC:
    return 1./sqrt(hh);
  case Cosmology::T:
    return a/sqrt(hh);
  case Cosmology::MassDeriv:
    return -0.5* pow(hh,-1.5) * (1.-a)*a; 
  case Cosmology::VacDeriv:
    return -0.5* pow(hh,-1.5) * (pow(a,ww)*expfunc(a) - a*a);
  case Cosmology::RadDeriv:
    return -0.5* pow(hh,-1.5) * (1.-a*a); 
  case Cosmology::W0Deriv: 
    return 1.5* pow(hh,-1.5) * pow(a,ww)*expfunc(a)*Ol
      * log(a);
  case Cosmology::W1Deriv:
    return -1.5* pow(hh,-1.5) * pow(a,ww)* expfunc(a) * Ol *
      (useWa? (-log(a) - (1.-a))
       : (log(a) + 1./a-1.) );
  default:
    throw CosmologyError("Unknown case in CosmoFunction op()");
  }
}

// Utility function used for many quantities
double
Cosmology::RskR(double d) const {
  if (invR0*d < 0.01) {
    // close to flat
    return d*(1+Ok*d*d/6.);
  } else if (Ok<0.) {
    /* closed */
    return sin(d*invR0)/invR0;
  } else if (Ok>0.) {
    /* open */
    return sinh(d*invR0)/invR0;
  } else {
    throw CosmologyError("Should not reach this in RskR");
  }
}

// Convert z to a, check for out of bounds
double
Cosmology::getA(double z) const {
  if (z<0.) throw CosmologyError("Negative redshift");
  double a=1./(1.+z);
  if (badW1 && a<aBounce) 
    throw CosmologyError("Cannot use w1>0 for high redshifts");
  if (doesBounce && a<aBounce) 
    throw CosmologyError("Redshift not attained in bouncing cosmology");
  return a;
}

double
Cosmology::Dc(double zsrc, double zobs) const {
  return getIntegral(zsrc,DC) - getIntegral(zobs,DC);
}

//proper-motion dist
double 
Cosmology::Dpm(double zsrc, double zobs) const {
  return RskR(getIntegral(zsrc,DC) - getIntegral(zobs,DC))/(1+zobs);
}

//angular diam dist
double
Cosmology::DA(double zsrc, double zobs) const {
  return RskR(getIntegral(zsrc,DC) - getIntegral(zobs,DC))/(1+zsrc);
}

//comoving angular diam dist
double
Cosmology::cDA(double zsrc, double zobs) const {
  return RskR(getIntegral(zsrc,DC) - getIntegral(zobs,DC));
}

//luminosity dist
double
Cosmology::DL(double zsrc, double zobs) const {
  return RskR(getIntegral(zsrc,DC) - getIntegral(zobs,DC))*(1+zsrc) /
    ((1+zobs)*(1+zobs));
}

//Hubble constant
double
Cosmology::H(double z) const {
  func.setMode(H2);
  double a=1./(1+z);
  return sqrt(func(a))/(a*a);
}

//Deceleration parameter
double
Cosmology::q(double z) const {
  func.setMode(Q);
  double a=1./(1+z);
  return func(a);
}

// dV/dz, comoving volume per unit solid angle per z
double 
Cosmology::dVdz(double zsrc) const {
  return pow(DA(zsrc),2.)/sqrt(func(1./(1+zsrc), H2));
}

// Lensing shear/mag strength, Dol Dls / Dos.
// Critical density is c^2/ (4 pi G LensShear c/H), which is
// 0.115 h g/cm^2 / LensShear
double
Cosmology::LensShear(double zsrc, double zlens) const {
  if (zsrc <= zlens || zsrc<=0.) return 0.;
  else return DA(zlens)*DA(zsrc,zlens)/DA(zsrc);
}

//Distance modulus for fixed rest-frame luminosity
double
Cosmology::mu(double zsrc) const {
  return 5.*log10(RskR(getIntegral(zsrc,DC))*HubbleDistance/(10.*Parsec)) 
    +5.*log10(1+zsrc);
}

double
Cosmology::DcDeriv(Parameter p,
		   double zsrc,
		   double zobs) const {
  FunctionType f;
  bool secondTerm=false;	//Need curvature deriv terms?
  
  switch (p) {
  case Omatter:
    f = MassDeriv;
    break;
  case Ovacuum:
    f = VacDeriv;
    break;
  case Oradiation:
    f = RadDeriv;
    break;
  case W:
    f = W0Deriv;
    break;
  case Wprime:
    f = W1Deriv;
    break;
  default:
    throw CosmologyError("Bad parameter for muDeriv");
  }
  
  return getIntegral(zsrc, f) - getIntegral(zobs, f);
}

double
Cosmology::muDeriv(Parameter p,
		   double zsrc) const {
  FunctionType f;
  bool secondTerm=false;	//Need curvature deriv terms?
  
  switch (p) {
  case Omatter:
    f = MassDeriv;
    secondTerm = true;
    break;
  case Ovacuum:
    f = VacDeriv;
    secondTerm = true;
    break;
  case Oradiation:
    f = RadDeriv;
    secondTerm = true;
    break;
  case W:
    f = W0Deriv;
    secondTerm = false;
    break;
  case Wprime:
    f = W1Deriv;
    secondTerm = false;
    break;
  default:
    throw CosmologyError("Bad parameter for muDeriv");
  }
  
  double dc = getIntegral(zsrc, DC);
  double result;

  // First term is common to all derivatives, depend upon geometry:
  // This is |kappa| S'(|kappa| I ) / S(|kappa| I) in Dragan
  if (invR0*dc < 0.01) {
    // Nearly flat
    result = (1. + Ok * dc * dc / 3.) / dc;
  } else if (Ok > 0.) { 
    // Open Universe:
    result = invR0 / tanh(invR0*dc);
  } else if (Ok < 0.) {
    // Closed Universe:
    result = invR0 / tan(invR0*dc);
  } else {
    throw CosmologyError("Should never get here #1");
  }
  result *= getIntegral(zsrc, f);	//times deriv of Dc w.r.t. parameter

  // Second term is for derivative of curvature w.r.t. parameter
  if (secondTerm) {
    if (invR0*dc < 0.01) {
      // Nearly flat
      result -= dc * dc / 6.;
    } else if (Ok > 0.) { 
      // Open Universe:
      result += (1. - dc* invR0 / tanh(invR0*dc)) / (2.*Ok);
    } else if (Ok < 0.) {
      // Closed Universe:
      result += (1. - dc* invR0 / tan(invR0*dc)) / (2.*Ok);
    } else {
      throw CosmologyError("Should never get here #2");
    }
  }
  result *= 5/log(10.);
  return result;
}

//comoving Horizon radius
double 
Cosmology::Horizon(double zobs) const {
  if (badW1) throw CosmologyError("Cannot use w1>0 to horizon.");
  if (doesBounce) throw CosmologyError("No horizon in bouncing cosmology");
  if (earlyIndex>=2.) throw CosmologyError("Divergent Horizon");
  return getIntegral(1./minA,DC) - getIntegral(zobs,DC);
}

//Age of Universe at chosen time
double 
Cosmology::Age(double zobs) const {
  if (badW1) throw CosmologyError("Cannot use w1>0 to large z for age");
  if (doesBounce) throw CosmologyError("No age for bouncing cosmology");
  if (earlyIndex>=4.) throw CosmologyError("Divergent Age");
  return getIntegral(1./minA,T) - getIntegral(zobs,T);
}

double 
Cosmology::LookbackTime(double zsrc, double zobs) const {
  return getIntegral(zsrc,T) - getIntegral(zobs,T);
}

Cosmology::~Cosmology() {
    for (int i=0; i<integrals.size(); i++) delete integrals[i];
}
// Create all the derived quantities from cosmology specs
void
Cosmology::setup() {
  func.reset();
  Ok = 1. - Om - Ol - Or;
  if (Ok<0.) {
    /* closed */
    invR0 = sqrt(-Ok);
  } else if (Ok==0.) {
    /* flat */
    invR0 = 0.;
  } else if (Ok>0.) {
    /* open */
    invR0 = sqrt(Ok);
  }

  // flush the integral tables
  for (int i=0; i<integrals.size(); i++) delete integrals[i];
  integrals.clear();

  // Determine the component which dominates energy at early times
  // The index is that for the H2 function - H^2 a^4.
  earlyOmega = 0.;
  earlyIndex = 99.;	//Should always be overridden by one of below
  // unless there is just cosmo with non-zero w1, which has no good
  // asymptotic behavior.
  if (Or>0.)
    if (0. < earlyIndex) {
      earlyIndex = 0.;
      earlyOmega = Or;
    }
  if (Om>0.)
    if (1. < earlyIndex) {
      earlyIndex = 1.;
      earlyOmega = Om;
    }
  if (Ok!=0.)
    if (2. < earlyIndex) {
      earlyIndex = 2.;
      earlyOmega = Ok;
    }
  if (Ol!=0. && w1==0.)
    if (1-3*w0 < earlyIndex) {
      earlyIndex = 1-3*w0;
      earlyOmega = Ol;
    }
  if (Ol!=0. && useWa && w1!=0.)
    if (1.-3*(w0+w1) < earlyIndex) {
      earlyIndex = 1.-3*(w0+w1); // ?? check this
      earlyOmega = Ol;
    }
  if (earlyOmega<0.) 
    throw CosmologyError("Should not have negative Omega at early times!");
  
  // Starting from present day, run back in time to see either
  // where there's a bounce, or where single-component model suffices.
  doesBounce = false;
  badW1 = false;
  earlyA = 0.;
  for (double a=1.; a>minA; a/=aRatio) {
    double h=func(a,H2);
    if (h<=0.) {
      // We have bounced.  Set the bounce redshift to the
      // previous step.  Accurate to lnStep of course.
      doesBounce = true;
      aBounce = a*aRatio;
      break;
    }
    if (Ol!=0. && w1>0. && ( !useWa && ( w0+w1*(1/a-1) > 2.)) ) {
      // Positive w1 is giving nonsensical behavior at this z or above
      // ??? check for this problem with Wa formulation
      badW1 = true;
      aBounce = a*aRatio;
      break;
    }
    if ( fabs( 1. - 
	       earlyOmega*pow(a,earlyIndex)
	       / func(a,H2Abs)) < Tolerance) {
      // Now have dominant component, analytic beyond this a.
      earlyA = a;
      break;
    }
    if (a<=minA) throw CosmologyError("Cannot find z for bounce/dominance");
  }
}


//Does this bounce?
bool 
Cosmology::Bounce(double *zBounce) const {
  if (doesBounce && zBounce!=0)
    *zBounce = 1./aBounce - 1.;
  return doesBounce;
}

// This routine gives an integral of any of the functions that you
// want.  Points along the integral stored as tables to save
// replication, so routine manages these tables.
double
Cosmology::getIntegral(double z, FunctionType ftype) const {
  double a = getA(z);	//This will throw exception if z is too high

  
  // First extend table vector to this function
  while (integrals.size() < ftype+1) {
    integrals.push_back(0);
  }
  // Make a table for this function if there is none
  if (!integrals[ftype]) {
    integrals[ftype] = new gtable::Table<>(gtable::Table<>::spline);
    // Set up tables with a few very low-z points
    gtable::Table<>& tab= *(integrals[ftype]);
    tab.addEntry(1., 0.);	//z=0;
    double da=sqrt(Tolerance)/2.;
    func.setMode(ftype);
    double value=Simp1d(func, 1.-da, 1., Tolerance);
    tab.addEntry(1.-da, value);
    da /= 2.;
    value=Simp1d(func, 1.-da, 1., Tolerance);
    tab.addEntry(1.-da, value);
  }

  gtable::Table<>& tab=*(integrals[ftype]);
  double result=0.;

  // See if we are going to extrapolate to early times with dominant
  // component.
  if (a<earlyA) {
    // Add on the analytic part of integral
    double delta;
    switch (ftype) {
    case DC:
      if (earlyIndex==2.) delta = log(earlyA / a);
      else delta = (pow(earlyA, 1-earlyIndex/2.)-pow(a, 1-earlyIndex/2.))
	     / (1-earlyIndex/2.);
      result += delta / sqrt(earlyOmega);
      break;
    case T:
      if (earlyIndex==4.) delta = log(earlyA / a);
      else delta = (pow(earlyA, 2-earlyIndex/2.)-pow(a, 2-earlyIndex/2.))
	     / (2-earlyIndex/2.);
      result += delta / sqrt(earlyOmega);
      break;
    default:
      // ??? currently won't be doing this for derivatives
      throw CosmologyError("Only set up to extrapolate Dc and T to early times");
    }
  }


  //value to look up in table:
  double lookupa = MAX(a,earlyA);

  // Extend the integral table to this value if needed.
  // The spline table should go a bit past our point, or to z=1:
  double tableA = MIN(lookupa/aRatio/aRatio,0.5);
  if (doesBounce) tableA=MAX(tableA,aBounce);	//don't exceed bounce z
  if (tableA < tab.argMin()) {
    // build table up to this value
    double cumulant=tab(tab.argMin());
    func.setMode(ftype);
    for (double a1=tab.argMin(); a1>tableA; a1/=aRatio) {
      double a2 = a1/aRatio;
      cumulant += Simp1d(func, a2, a1, Tolerance);
      tab.addEntry(a2, cumulant);
    }
  }

  // Look up our value in the table.
  result += tab(lookupa);
  return result;
}
