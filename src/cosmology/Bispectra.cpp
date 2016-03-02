// $Id: Bispectra.cpp,v 1.10 2006-08-08 20:24:04 garyb Exp $
#define TAKADA

#include "Bispectra.h"

#ifdef TAKADA
#include "Cosmology.h"
#include "GrowthFunction.h"
#endif

using namespace cosmology;

// Step sizes for derivative calculation:
const double domegaM=0.01;
const double domegaB=0.001;
const double dns = 0.01;
const double dlnG = 0.02;
const double dlnAs = dlnG;
const double dlnk = 0.02;

const double d2p = 0.01;  // ??? this is ok for OmegaM, ugly for other cases!

#ifdef TAKADA
double savedSig8;
#endif

void
PTBispectrum::setPerturbation(double z, PerturbationType pt,
			      bool ifUp, int iParam) const {
  savedZ = z;
#ifdef TAKADA
  /**/static Cosmology c(0.27, 1.-0.27);
  /**/static GrowthFunction gf(c);
  savedSig8 = 0.9*gf(savedZ)/gf(0.); /**/
#endif
  pType = pt;
  perturbUp = ifUp;
  ip = iParam;
  switch (pType) {
  case PNone:
    dp = 0.;
    break;
  case PLnK:
    dp = dlnk;
    break;
  case PLnG:
    dp = dlnG;
    break;
  case P2p:
    dp = d2p;
    break;
  case PParam:
    {
      const FisherParameter* fp=dependsOn(1., 1., 0., savedZ)[iParam];
      if ( dynamic_cast<const omegaMParam*> (fp) )
	dp = domegaM;
      else if ( dynamic_cast<const omegaBParam*> (fp) )
	dp = domegaB;
      else if ( dynamic_cast<const nsParam*> (fp) )
	dp = dns;
      else if ( dynamic_cast<const lnAParam*>  (fp) )
	dp = dlnAs;
      else 
	throw MyException("Unknown parameter type in PTBispectrum::setFunction");
    }
    break;
  default:
    throw MyException("Unreachable line in PTBispectrum::setFunction");
  }
}

double 
PTBispectrum::PNL(double k) const {
  double increment = perturbUp? +0.5*dp : -0.5*dp;
  switch (pType) {
  case PNone:
    return nlps(k, savedZ);
  case PLnK:
    return nlps(k,savedZ)*exp(increment*nlps.dlnPdlnk(k,savedZ));
  case PLnG:
    return nlps(k,savedZ)*exp(increment*nlps.dlnPdlnG(k,savedZ));
  case P2p:
    return nlps(k,savedZ)*exp(increment*nlps.dlnPd2p(k,savedZ));
  case PParam:
    return nlps(k,savedZ)*exp(increment*nlps.dlnPdParam(k,savedZ,ip));
  default:
    throw MyException("Unreachable line in PTBispectrum::P");
  }
}

double
PTBispectrum::F(double k1, double k2, double costheta3, 
		double a1, double b1, double c1,
		double a2, double b2, double c2) const {
  return (5./7.)*a1*a2 
    + 0.5*(k1/k2 + k2/k1)*costheta3*b1*b2
    + (2./7.)*costheta3*costheta3*c1*c2;
};

double 
PTBispectrum::B(double k1, double k2, double costheta12) const {

  double k3 = sqrt(k1*k1+k2*k2+2*k1*k2*costheta12);
  double P1 = PNL(k1);
  double P2 = PNL(k2);
  double P3 = PNL(k3);

  double costheta23 = (k1*k1 - k2*k2 - k3*k3)/(2*k2*k3);
  double costheta31 = (k2*k2 - k3*k3 - k1*k1)/(2*k3*k1);

  double a1,b1,c1;
  double a2,b2,c2;
  double a3,b3,c3;
  abc(a1,b1,c1,
      a2,b2,c2,
      a3,b3,c3,
      k1,k2,k3);
  
  /**cerr << "State " << pType << " " << perturbUp
	   << " k's: " << k1 << ", " << k2 << ", " <<k3 
	   << " abc: "
	   << a1 << " " << b1 << " " << c1 << endl;//*/


  return 2.* (F(k1,k2,costheta12,a1,b1,c1,a2,b2,c2)*P1*P2 +
	      F(k2,k3,costheta23,a2,b2,c2,a3,b3,c3)*P2*P3 +
	      F(k3,k1,costheta31,a3,b3,c3,a1,b1,c1)*P3*P1);
}

double 
PTBispectrum::operator()(double k1, double k2, 
			 double costheta, double z) const {
  setPerturbation(z, PNone);
  return B(k1,k2,costheta);
}

double 
PTBispectrum::dBdlnk(double k1, double k2, 
			 double costheta, double z) const
{
  setPerturbation(z,PNone); // Just handle the derivative here
  const double kfac=exp(dlnk/2.);
  return ( B(k1*kfac, k2*kfac, costheta)
	   - B(k1/kfac, k2/kfac, costheta)) / dlnk;
}

double 
PTBispectrum::dBdlnG(double k1, double k2, 
			 double costheta, double z) const
{
  setPerturbation(z, PLnG, true);
  double b = B(k1,k2,costheta);
  setPerturbation(z, PLnG, false);
  return (b - B(k1,k2,costheta)) / dp;
}

double 
PTBispectrum::dBd2p(double k1, double k2, 
		    double costheta, double z) const
{
  if (!has2p()) return 0.;
  setPerturbation(z, P2p, true);
  double b = B(k1,k2,costheta);
  setPerturbation(z, P2p, false);
  return (b - B(k1,k2,costheta)) / dp;
}

const ParameterVector& 
PTBispectrum::dependsOn(double k1, double k2, 
			double costheta, double z) const {
  return nlps.dependsOn(k1, z);
}

double
PTBispectrum::dBdParam(double k1, double k2, 
		       double costheta, double z,
		       int iParam) const
{
  setPerturbation(z, PParam, true, iParam);
  double b = B(k1,k2,costheta);
  setPerturbation(z, PParam, false, iParam);
  double b2 = B(k1,k2,costheta);
  //  cerr << "b+ " << b << " b- " << b2 << " dp " << dp << endl;
  return (b-b2)/dp;
}

//////////////////////////////////////////////////////////////
// Stuff to implement the HEPT abc's, needing values from the
// linear power spectrum each time.


// Function to produce linear power spectrum or perturbations to it
// for derivatives.  Class variable determines which function
// is being returned.
double
HEPTBispectrum::lnDSqFunc(double lnK) const {
  double increment = perturbUp? +0.5*dp : -0.5*dp;
  double k = exp(lnK);
  double lplin = nlps.lnLinearDeltaSq(k, savedZ);
  switch (pType) {
  case PNone:
    return lplin;
  case PLnK:
    return nlps.lnLinearDeltaSq(k*exp(increment), savedZ);
  case PLnG:
    return lplin + 2*increment;
  case P2p:
    return lplin;	//2nd parameter does not affect linear spectrum
  case PParam:
    return lplin + (increment*nlps.dlnLinearPdParam(k,savedZ,ip));
  default:
    throw MyException("Unreachable line in PTBispectrum::lnDSqFunc");
  }
}

double
HEPTBispectrum::neff(double lnK) const {
  return (lnDSqFunc(lnK+dlnk/2.) -
	  lnDSqFunc(lnK-dlnk/2.) ) / dlnk - 3.;
}

double
HEPTBispectrum::kNL() const {
  // Find kNL for the saved redshift.  This kNL is defined as simply
  // Delta-squared(k)=1.
  // ??? could cache last result
  double lastLnK = savedZ-2.; // A rough initial guess
  double lastDsq = lnDSqFunc(lastLnK);

  int nsteps=0;
  const int MaxSteps = 10;
  const double lnKTolerance = 0.001;
  const double dSqTolerance = 0.001;


  while (nsteps<MaxSteps) {
    double nextLnK = lastLnK - lastDsq/(neff(lastLnK)+3.);
    double nextDsq = lnDSqFunc(nextLnK);
    if (abs(nextDsq)<dSqTolerance &&
	abs(nextLnK-lastLnK) < lnKTolerance) {
      //**/ cerr << "nsteps " << nsteps << endl;
      return exp(nextLnK);
    }
    lastLnK = nextLnK;
    lastDsq = nextDsq;
    nsteps++;
  }

  // Get here if no solution:
  throw MyException("Too many steps for kNL in HEPTBispectrum");
}

void
HEPTBispectrum::abc(double& a1, double& b1, double& c1,
		    double& a2, double& b2, double& c2,
		    double& a3, double& b3, double& c3,
		    double k1, double k2, double k3) const {
  double knl=kNL();
  SCabc(a1,b1,c1,k1,knl,neff(log(k1)));
  SCabc(a2,b2,c2,k2,knl,neff(log(k2)));
  SCabc(a3,b3,c3,k3,knl,neff(log(k3)));
  //**/cerr << " kNL " << knl << " neff " << neff(log(k1))  << endl;
}

void 
HEPTBispectrum::SCabc(double& a, double& b, double& c, 
		      double k, double kNL, double neff) const {
  // From Scoccimarro & Couchman (1991) via Takada & Jain, but I'm using
  // the SC version that does NOT depend upon sigma8, since it was in
  // 8h^-1 Mpc tophat, and I don't want h in here:

#ifdef TAKADA
  /* Here are the values that Masahiro used - these need a factor
  ** of sigma8^-0.2 in the expression for a.*/
  const double a1=0.25;
  const double a2=3.5;
  const double a3=2.;
  const double a4=1.;
  const double a5=2.;

  double q = k / kNL;
  double Q3 = (4-pow(2.,neff)) /
    (1 + pow(2., neff+1));


  a = ( 1 + pow(savedSig8,-0.2)*sqrt(0.7*Q3)*pow(a1*q, neff+a2)) /
    (1 + pow(a1*q, neff+a2));
#else
  // Here are the kludge factors they give:
  const double a1=0.25;
  const double a2=4.;
  const double a3=1.5;
  const double a4=1.;
  const double a5=1.5; 

  double q = k / kNL;
  double Q3 = (4-pow(2.,neff)) /
    (1 + pow(2., neff+1));


  a = ( 1 + sqrt(0.7*Q3)*pow(a1*q, neff+a2)) /
    (1 + pow(a1*q, neff+a2));
#endif

  b = (1 + 0.2*a3*(neff+3.)*pow(q,neff+3.)) / 
    ( 1 + pow(q,neff+3.5));

  c = ( 1 + 4.5*a4*pow(a5*q,neff+3.)/(1.5+pow(neff+3,4.))) /
    ( 1 + pow(a5*q,neff+3.5));


  /*cerr << "kNL: " << kNL
	   << " neff: " << neff
	   << " a,b,c: " << a << " " << b << " " << c
	   << endl; //**/
}

 
