// $Id: PowerSpectra.h,v 1.7 2007-11-26 03:30:39 garyb Exp $
// Cosmological Power Spectra.
// Note units of k are inverse Mpc - NOT h Mpc^-1 - for everything.
#ifndef POWERSPECTRA_H
#define POWERSPECTRA_H
#include <cmath>
#include "GTable.h"
#include "FisherMatrices.h"
#include "AstronomicalConstants.h"

namespace cosmology {
  // Abstract base classes for linear and non-linear power spectra
  // Convention in these classes is that k is in inverse Mpc (no h's).

  class LinearPowerSpectrum {
  public:
    LinearPowerSpectrum();
    virtual double DeltaSq(double lnK) const =0;
    virtual double sigSq(double lnK) const;
    virtual double neff(double lnK) const;
    virtual double C(double lnK) const;	// Curvature
    virtual double lnkNL(double G) const; // solve for sigSq(kNL)=1.
    virtual double sigX(double x) const; // RMS in tophats of radius x Mpc
    virtual ~LinearPowerSpectrum() {}
  private:
    void integrals(double k0, 
		   double& sigSq, double& neff, double& C) const;
    // All of the integral quantities are cached in spline tables:
    void extendTables(double lnK) const;
    mutable gtable::Table<> sigSqTable;
    mutable gtable::Table<> neffTable;
    mutable gtable::Table<> CTable;
  };

  class NonLinearPowerSpectrum {
  public:
    virtual double operator()(double k, double z) const =0;
    virtual double dlnPdlnk(double k, double z) const =0;
    virtual double dlnPdlnG(double k, double z) const=0;
    virtual double dlnPdz(double k, double z) const=0;
    // Which other parameter id's have partials for this?
    virtual const ParameterVector& dependsOn(double k, double z) const =0;
    virtual double dlnPdParam(double k, double z, int iParam) const =0;

    // Is there a 2nd parameter besides growth factor that governs
    // the non-linear behavior?
    virtual bool has2p() const {return false;}
    // If so, there's a derivative for it:
    virtual double dlnPd2p(double k, double z) const {return 0.;}

    // Sometimes we need to know the underlying linear model:
    virtual double lnLinearDeltaSq(double k, double z) const =0;
    virtual double dlnLinearPdParam(double k, double z, int iParam) const =0;

    virtual ~NonLinearPowerSpectrum() {}
  };


  // Abstract base class for transfer functions
  class TransferFunction {
  public:
    virtual double operator()(double k) const=0;
    virtual ~TransferFunction() {};
  };

  // The Bond-Efstathiou transfer function.  Not used right now.
  class TransferFunctionBE: public TransferFunction {
  private:
    double omegaM;
  public:
    TransferFunctionBE(double omegaM_):omegaM(omegaM_) {}
    double operator()(double k) const {
      double q=k/omegaM;
      return pow(1+pow(6.4*q+pow(3.0*q,1.5) + 1.7*1.7*q*q ,1.13), -1.0/1.13);
    }
  };

  // The Eisenstein-Hu transfer function.
  class TransferFunctionEH: public TransferFunction {
  public:
    TransferFunctionEH(double omegaM, double omegaB);
    double operator()(double k) const;
  private:
    /* The following are set in TFmdm_set_cosm() */
    double  alpha_gamma;	/* sqrt(alpha_nu) */
    double  alpha_nu;	/* The small-scale suppression */
    double  beta_c;	/* The correction to the log in the small-scale */
    double  num_degen_hdm;/* Number of degenerate massive neutrino species */
    double  f_baryon;	/* Baryon fraction */
    double  f_bnu;	/* Baryon + Massive Neutrino fraction */
    double  f_cb;		/* Baryon + CDM fraction */
    double  f_cdm;	/* CDM fraction */
    double  f_hdm;	/* Massive Neutrino fraction */
    double  growth_k0;	/* D_1(z) -- the growth function as k->0 */
    double  growth_to_z0;	/* D_1(z)/D_1(0) -- the growth relative to z=0 */
    double  k_equality;	/* The comoving wave number of the horizon at equality*/
    double  obhh;		/* Omega_baryon * hubble^2 */
    double  omega_curv;	/* = 1 - omega_matter - omega_lambda */
    double  omega_lambda_z; /* Omega_lambda at the given redshift */
    double  omega_matter_z;/* Omega_matter at the given redshift */
    double  omhh;		/* Omega_matter * hubble^2 */
    double  onhh;		/* Omega_hdm * hubble^2 */
    double  p_c;		/* The correction to the exponent before drag epoch */
    double  p_cb;		/* The correction to the exponent after drag epoch */
    double  sound_horizon_fit;  /* The sound horizon at the drag epoch */
    double  theta_cmb;	/* The temperature of the CMB, in units of 2.7 K */
    double  y_drag;	/* Ratio of z_equality to z_drag */
    double  z_drag;	/* Redshift of the drag epoch */
    double  z_equality;	/* Redshift of matter-radiation equality */
  };

  // A linear spectrum with power-law primordial spectrum
  const double DefaultDeltaSqCurvature=pow(5.07e-5,2.); // from HuJain
  const double DefaultKNorm=0.05;

  class PowerLawPS: public LinearPowerSpectrum {
  public:
    PowerLawPS(const TransferFunction* tf_,
	       double omegaM_,
	       double ns_=1.,
	       double DeltaSqCurvature=DefaultDeltaSqCurvature,
	       double kNorm=DefaultKNorm): tf(tf_), ns(ns_) {
      normalize(kNorm, DeltaSqCurvature, omegaM_);
    }
    // As a default, we use the Eisenstein-Hu transfer function
    PowerLawPS(double omegaM_, double omegaB_, 
	       double ns_=1.,
	       double DeltaSqCurvature=DefaultDeltaSqCurvature,
	       double kNorm=DefaultKNorm): 
      tf(new TransferFunctionEH(omegaM_, omegaB_)), ns(ns_) {
	normalize(kNorm, DeltaSqCurvature, omegaM_);
      }
    ~PowerLawPS() {delete tf;}
    double DeltaSq(double lnK) const {
      return primordial(exp(lnK))*pow((*tf)(exp(lnK)),2.);
    }
  private:
    // Hide copy constructor and assignment to avoid trouble with
    // the pointer to TransferFunction;
    PowerLawPS(const PowerLawPS& rhs) {}
    void operator=(const PowerLawPS& rhs) {}

    const TransferFunction* tf;
    double norm;
    double ns;

    // The normalization will set linear power to be correctly
    // normalized at high (matter-dominated) redshift, if the multiplying
    // growth function is G=a.
    // Assuming that curvature fluctuations at kNorm in matter-dominated period
    // are specified.
    void normalize(double kNorm, double DeltaSqCurvature, double omegaM) {
      norm = DeltaSqCurvature * pow(kNorm, -ns-3.) * 
	pow(omegaM, -2.) * (4./25.)* pow(kNorm*HubbleLengthMpc, 4.);
    }
    double primordial(double k) const {return norm * pow(k, ns+3.);}
  };
    

  } // namespace cosmology

#endif //POWERSPECTRA_H
