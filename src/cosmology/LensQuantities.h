// $Id: LensQuantities.h,v 1.6 2006-08-07 13:22:00 garyb Exp $
// Classes to give standard lensing power spectra, etc.
#include "FisherMatrices.h"
#include "FiducialCosmology.h"
#include "PowerSpectra.h"
#include "Bispectra.h"
#include "Table.h"
#include "AstronomicalConstants.h"

#ifndef LENSQUANTITIES_H
#define LENSQUANTITIES_H

namespace cosmology {

  // Abstract base class for dN/dz 
  class ZDist {
  public:
    virtual double operator()(double z) const=0;
    virtual double integral(double z0, double z1) const;
    virtual double integral() const {return integral(zMin(),zMax());}
    virtual double zMin() const=0;
    virtual double zMax() const=0;
    virtual ~ZDist() {}
  };

  // and 3 particular implementations 
  class ZDistExp: public ZDist {
  private:
    double z0;
    double norm;
    double zmin;
    double zmax;
    static const double DefaultZMax=5.;
  public:
    ZDistExp(double zmed, double density, 
	     //***	     double zmin_=0., double zmax_=DefaultZMax): z0(zmed/1.414), 
	     double zmin_=0., double zmax_=5.): z0(zmed/1.414), 
      norm(density*pow(z0,-3.)/0.667), zmin(zmin_), zmax(zmax_) {}
    double operator()(double z) const {
      return z*z*exp(-pow(z/z0,1.5))*norm;
    }
    double zMin() const {return zmin;}
    double zMax() const {return zmax;}
  };

  class ZDistTakada: public ZDist {
  private:
    double z0;
    double norm;
    double zmin;
    double zmax;
    static const double DefaultZMax=5.;
  public:
    ZDistTakada(double zmed, double density, 
		double zmin_=0., double zmax_=DefaultZMax): z0(zmed/3.), 
      norm(density*0.5*pow(z0,-3.)), zmin(zmin_), zmax(zmax_) {}
    double operator()(double z) const {
      return z*z*exp(-z/z0)*norm;
    }
    double zMin() const {return zmin;}
    double zMax() const {return zmax;}
  };

  /* class ZDist... : public ZDist
  // Make one with arbitrary pow(z,alpha)*exp(-pow(z/z0,beta)).
  // Note the normalization becomes
  // z0^{1+alpha}*Gamma( (1+alpha)/beta ) / beta.
  // Need to find median numerically. */

  class ZDistTable: public ZDist {
  private:
    Table<> logNTab;
    double zmin;
    double zmax;
  public:
    ZDistTable(istream& is, double zmin=-1., double zmax=-1.): logNTab(is) {
      if (zmin<0.) zmin=logNTab.argMin();
      if (zmax<0.) zmax=logNTab.argMax();
    }
    // Table was log10 of number per square degree
    double operator()(double z) const {
      return pow(10.,logNTab(z))/DEGREE/DEGREE;
    }
    double zMin() const {return zmin;}
    double zMax() const {return zmax;}
  };

  // I'm NOT writing a delta-function ZDist since it's singular

  // Abstract class for the lensing weight function in projection
  class LensWeight {
  public:
    virtual double operator()(double zl) const =0;
    virtual double zsMin() const =0;	//Least distant source
    virtual double zsMax() const =0;	//Most distant source
    virtual double shearNoise() const=0;  //White-noise level for these sources
    virtual ~LensWeight() {}
  };

  // And two derived implementations:
  class SingleSourceLensWeight: public LensWeight {
  public:
    SingleSourceLensWeight(double zs_, const Cosmology& c_,
			   double density=10/ARCMIN/ARCMIN,
			   double sigmaGamma=0.3): 
      zs(zs_), c(c_), noise(sigmaGamma*sigmaGamma/density) {}
    virtual double operator()(double zl) const {
      return zs>zl ? ( c.cDA(zs,zl) / c.cDA(zs)) : 0.;
    }
    virtual double zsMax() const {return zs;}
    virtual double zsMin() const {return zs;}
    virtual double shearNoise() const {return noise;}
  private:
    double zs;
    const Cosmology& c;
    const double noise;
  };

  class ZDistLensWeight: public LensWeight {
  public:
    ZDistLensWeight(const ZDist& zd_, const Cosmology& c_,
		    double sigmaGamma=0.3): zd(zd_), c(c_) {
      nTot = zd.integral(); 
      noise = sigmaGamma*sigmaGamma / nTot;
    }
    virtual double operator()(double zl) const;
    virtual double zsMin() const {return zd.zMin();}
    virtual double zsMax() const {return zd.zMax();}
    double shearNoise() const {return noise;}
  private:
    const ZDist& zd;
    const Cosmology& c;
    double nTot;
    double noise;
  };

  // Integrate the power spectrum down the line of sight:
  class Pkappa {
  public:
    Pkappa(const LensWeight& lw_,
	   const FiducialCosmology& fc_,
	   const NonLinearPowerSpectrum& nlps_): lw1(lw_), lw2(lw_), 
      fc(fc_), nlps(nlps_),
      pktable(Table<>::spline) {}
    Pkappa(const LensWeight& lw1_,
	   const LensWeight& lw2_,
	   const FiducialCosmology& fc_,
	   const NonLinearPowerSpectrum& nlps_): lw1(lw1_), lw2(lw2_),
      fc(fc_), nlps(nlps_),
      pktable(Table<>::spline) {}
    double operator()(double l) const;
    double valueAt(double l) const;	//compute without lookup table
  private:
    const LensWeight& lw1;
    const LensWeight& lw2;
    const FiducialCosmology& fc;
    const NonLinearPowerSpectrum& nlps;
    mutable Table<> pktable;
  };

  // Get Map and correlation functions from power spectrum:
  double
    MapSq(double theta, const Pkappa& pk);
  double
    XiPlus(double theta, const Pkappa& pk);
  double
    XiMinus(double theta, const Pkappa& pk);

  // The functions that convert correlation function into MapSq:
  extern
  double Tplus(double x);
  extern
  double Tminus(double x);

  // Integrate the power spectrum down the line of sight:
  class Bkappa {
  public:
    Bkappa(const LensWeight& lw1_,
	   const LensWeight& lw2_,
	   const LensWeight& lw3_,
	   const FiducialCosmology& fc_,
	   const Bispectrum& b3d_): 
      lw1(lw1_), lw2(lw2_), lw3(lw3_), fc(fc_), b3d(b3d_) {}
    double operator()(double l1, double l2, double costheta) const;
  private:
    const LensWeight& lw1;
    const LensWeight& lw2;
    const LensWeight& lw3;
    const FiducialCosmology& fc;
    const Bispectrum& b3d;
  };

} // namespace cosmology

#endif // LENSQUANTITIES_H
