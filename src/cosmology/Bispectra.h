// $Id: Bispectra.h,v 1.4 2006-07-24 22:16:40 garyb Exp $
// Classes to calculate bispectra and related quantitities
// Units of wavevectors are always comoving Mpc^{-1}.

#include "PowerSpectra.h"

#ifndef BISPECTRA_H
#define BISPECTRA_H

namespace cosmology {

  class Bispectrum {
  public:
    virtual double operator()(double k1, double k2, 
			      double costheta, double z) const=0;
    // Derivative if all k's are scaled by some factor:
    virtual double dBdlnk(double k1, double k2, 
			      double costheta, double z) const=0;
    virtual double dBdlnG(double k1, double k2, 
			      double costheta, double z) const=0;

    // Is there a 2nd parameter for non-linear growth?
    virtual bool has2p() const {return false;}
    virtual double dBd2p(double k1, double k2, 
			 double costheta, double z) const { return 0.;}
    // Which other parameter id's have partials for this?
    virtual const ParameterVector& dependsOn(double k1, double k2, 
			      double costheta, double z) const=0;
    virtual double dBdParam(double k1, double k2, 
			    double costheta, double z,
			    int iParam) const=0;
    virtual ~Bispectrum() {};
  };

  class PTBispectrum: public Bispectrum {
  public:
    PTBispectrum(const NonLinearPowerSpectrum& nlps_): nlps(nlps_)  {}
    double operator()(double k1, double k2, 
		      double costheta, double z) const;
    double dBdlnk(double k1, double k2, 
		    double costheta, double z) const;
    double dBdlnG(double k1, double k2, 
		    double costheta, double z) const;
    virtual bool has2p() const {return nlps.has2p();}
    virtual double dBd2p(double k1, double k2, 
			 double costheta, double z) const;
    
    // Which other parameter id's have partials for this?
    const
    ParameterVector& dependsOn(double k1, double k2, 
			       double costheta, double z) const;
    double dBdParam(double k1, double k2, 
		    double costheta, double z,
		    int iParam) const;
    virtual ~PTBispectrum() {}
  protected:
    const NonLinearPowerSpectrum& nlps;
    // In calculating derivatives, the following tell what kind we're doing
    // at the moment:
    enum PerturbationType {PNone, PLnK, PLnG, P2p, PParam};
    mutable double savedZ;
    mutable PerturbationType pType;
    mutable int ip;
    mutable bool perturbUp;
    mutable double dp;	// size of perturbation to parameter
    void setPerturbation(double z, PerturbationType pt, 
			 bool ifUp=true,
			 int iParam=-1) const;

    // Get the abc values, using saved redshift & perturbation
    virtual void abc(double& a1, double& b1, double& c1,
		     double& a2, double& b2, double& c2,
		     double& a3, double& b3, double& c3,
		     double k1, double k2, double k3) const =0;

  private:
    double F(double k1, double k2, double costheta3, 
	     double a1, double b1, double c1,
	     double a2, double b2, double c2) const;
    // These functions are calculated with the current perturbation type & z:
    double B(double k1, double k2, double costheta3) const;
    double PNL(double k) const;
  };

  class HEPTBispectrum: public PTBispectrum {
  public:
    HEPTBispectrum(const NonLinearPowerSpectrum& nlps_): 
      PTBispectrum(nlps_) {}
  protected:
    void abc(double& a1, double& b1, double& c1,
	     double& a2, double& b2, double& c2,
	     double& a3, double& b3, double& c3,
	     double k1, double k2, double k3) const;
  private:
    // Helper functions to query the linear power spectrum for the HEPT a,b,c's:
    double kNL() const;
    double lnDSqFunc(double lnK) const;
    double neff(double lnK) const;
    void SCabc(double& a, double& b, double& c, 
	       double k, double kNL, double neff) const;
  };

  // Tree-level perturbation theory
  class TreePTBispectrum: public PTBispectrum {
  public:
    TreePTBispectrum(const NonLinearPowerSpectrum& nlps_): 
      PTBispectrum(nlps_) {}
  protected:
    void abc(double& a1, double& b1, double& c1,
	     double& a2, double& b2, double& c2,
	     double& a3, double& b3, double& c3,
	     double k1, double k2, double k3) const {
      a1=b1=c1=a2=b2=c2=a3=b3=c3=1.;
    };
  };
}  // namespace cosmology
#endif  // BISPECTRA_H
