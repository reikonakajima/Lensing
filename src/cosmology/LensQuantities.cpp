// $Id: LensQuantities.cpp,v 1.7 2006-08-08 00:27:04 garyb Exp $
#include "LensQuantities.h"
#include "Simpson.h"
#include <cmath>

using namespace cosmology;

// Fractional tolerance for all integrations:
const double Tolerance=1e-4;

// The functions that convert correlation function into MapSq:
double Tplus(double x) {
  double xx=x*x;
  return (1 - 2*asin(x/2)/PI)*6*(2-15*xx)/5
    + x*sqrt(4-xx)*(120 + xx*(2320+xx*(-754+xx*(132-9*xx))))/(100*PI);
}
double Tminus(double x) {
  return 192*pow(x,3.)*pow( 1-x*x/4,3.5)/(35*PI);
}


double
ZDist::integral(double z0,
		double z1) const {
  return Simp1d(*this, z0, z1,Tolerance);
}

class WeightIntegrand {
public:
  WeightIntegrand(const ZDist& zd_, const Cosmology& c_,
		  double zl_): zd(zd_), c(c_), zl(zl_) {}
  double operator()(double zs) const {
    return zd(zs) * c.cDA(zs,zl) / c.cDA(zs);
  }
private:
  const ZDist& zd;
  const Cosmology& c;
  const double zl;
};

double
ZDistLensWeight::operator()(double zl) const {
  if (zl>=zd.zMax()) return 0.;
  WeightIntegrand I(zd,c,zl);
  return Simp1d(I, MAX(zl,zd.zMin()), zd.zMax(), Tolerance)/nTot;
}

class PkIntegrand {
public:
  PkIntegrand(const LensWeight& lw1_,
	      const LensWeight& lw2_,
	      const FiducialCosmology& fc_,
	      const NonLinearPowerSpectrum& nlps_,
	      double l_): 
    lw1(lw1_), lw2(lw2_), fc(fc_), nlps(nlps_), l(l_) {}
  double operator()(double zl) const {
    double ww = lw1(zl) * lw2(zl);
    if (ww<=0.) return 0.;
    /*cerr << "PkI at zl=" << zl
	     << " nlps " << nlps(l/fc.cDA(zl)/HubbleLengthMpc,zl) 
	     << " wt " << w
	     << endl;*/
    return ww*nlps(l/fc.cDA(zl)/HubbleLengthMpc,zl) 
      * (1+zl)*(1+zl) / fc.H(zl);
  }
private:
  const LensWeight& lw1;
  const LensWeight& lw2;
  const FiducialCosmology& fc;
  const NonLinearPowerSpectrum& nlps;
  const double l;
};


double
Pkappa::valueAt(double l) const {
  double minZ=0.001; // ??? kludge
  PkIntegrand I(lw1,lw2,fc,nlps,l);
  return Simp1d(I, minZ, MIN(lw1.zsMax(), lw2.zsMax()), Tolerance)
    * ((9/4.)*fc.omegaM()*fc.omegaM()*pow(HubbleLengthMpc,-3.));
}

double
Pkappa::operator()(double l) const {
  // Maintain a lookup table in (log l, log Pk) space.
  const double logLStep=0.05;	// Step size for table entries
  const int nBuffer=4;		// Number of table points needed above/below
  
  // Extend table as needed:

  double nearLogL = logLStep * floor(log(l)/logLStep + 0.5);

  for (int i=-nBuffer; i<=nBuffer; i++) {
    double logL= nearLogL + i*logLStep;
    if (pktable.size()<=0 
	|| logL<pktable.argMin()
	|| logL>pktable.argMax()) {
      // Add new point to the table:
      double y = log(valueAt(exp(logL)));
      pktable.addEntry(logL, y);
    }
  }

  return exp(pktable(log(l)));
}

class MapSqIntegrand {
public:
  MapSqIntegrand(double theta_, const Pkappa& pk_): theta(theta_),
						    pk(pk_) {}
  double operator()(double l) const {
    double lt = l*theta;
    if (lt<=0.) return 0.;
    double j=jn(4,lt);
    return 576*j*j*pow(lt,-4.)*pk(l)*l;
  }
private:
  const double theta;
  const Pkappa& pk;
};

double
MapSq(double theta, const Pkappa& pk) {
  MapSqIntegrand I(theta,pk);
  return Simp1d(I, 0., 20./theta, Tolerance) / (2*PI);;
}

class XiPlusIntegrand {
public:
  XiPlusIntegrand(double theta_, const Pkappa& pk_): theta(theta_),
						     pk(pk_) {}
  double operator()(double l) const {
    double lt = l*theta;
    if (lt<=0.) return 0.;
    return j0(lt)*pk(l)*l;
  }
private:
  const double theta;
  const Pkappa& pk;
};

double
XiPlus(double theta, const Pkappa& pk) {
  XiPlusIntegrand I(theta,pk);
  return Simp1d(I, 0., 20./theta, Tolerance) / (2*PI);;
}

class XiMinusIntegrand {
public:
  XiMinusIntegrand(double theta_, const Pkappa& pk_): theta(theta_),
						     pk(pk_) {}
  double operator()(double l) const {
    double lt = l*theta;
    if (lt<=0.) return 0.;
    return jn(4,lt)*pk(l)*l;
  }
private:
  const double theta;
  const Pkappa& pk;
};

double
XiMinus(double theta, const Pkappa& pk) {
  XiMinusIntegrand I(theta,pk);
  return Simp1d(I, 0., 20./theta, Tolerance) / (2*PI);;
}

class BkIntegrand {
public:
  BkIntegrand(const LensWeight& lw1_,
	      const LensWeight& lw2_,
	      const LensWeight& lw3_,
	      const FiducialCosmology& fc_,
	      const Bispectrum& b3d_,
	      double l1_,
	      double l2_,
	      double costheta): 
    lw1(lw1_), lw2(lw2_), lw3(lw3_), fc(fc_), b3d(b3d_),
  l1(l1_), l2(l2_), ct(costheta) {}
  double operator()(double zl) const {
    double w1 = lw1(zl);
    if (w1<=0.) return 0.;
    double w2 = lw2(zl);
    if (w2<=0.) return 0.;
    double w3 = lw3(zl);
    if (w3<=0.) return 0.;

    double dA = fc.cDA(zl) * HubbleLengthMpc; // Put in Mpc
    double H = fc.H(zl) / HubbleLengthMpc; // Put in  inv Mpc

    /*cerr << "zl " << zl
	     << " w " << w1
	     << " b " << b3d(l1/dA, l2/dA, ct, zl) 
	     << " HdA " << dA * H
	     << endl; //**/
    return w1*w2*w3*pow(1+zl, 3.) * b3d(l1/dA, l2/dA, ct, zl) / (dA*H);
  }
private:
  const LensWeight& lw1;
  const LensWeight& lw2;
  const LensWeight& lw3;
  const FiducialCosmology& fc;
  const Bispectrum& b3d;
  const double l1;
  const double l2;
  const double ct;
};

double
Bkappa::operator()(double l1, double l2, double costheta) const {
  double minZ=0.001; // ??? kludge
  double maxZ = MAX(lw1.zsMax(), lw2.zsMax());
  maxZ = MAX(maxZ, lw3.zsMax());
  BkIntegrand I(lw1,lw2,lw3,fc,b3d,l1,l2,costheta);
  return Simp1d(I, minZ, maxZ,Tolerance)
    * pow(1.5*fc.omegaM(), 3.) * pow(HubbleLengthMpc,-6.);
}





