// $Id: SmithPS.h,v 1.8 2007-11-26 03:30:39 garyb Exp $
#ifndef SMITHPS_H
#define SMITHPS_H
#include "AstronomicalConstants.h"
#include "FiducialCosmology.h"
#include "GrowthFunction.h"
#include "PowerSpectra.h"

namespace cosmology {

  class SmithPS: public NonLinearPowerSpectrum {
  public:
    // Defaults in constructor are WMAP values; the last flag is set if
    // normalization is to primordial spectrum (T=1)
    SmithPS(FiducialCosmology& fc_, double ns=1., 
	    double norm=DefaultDeltaSqCurvature,
	    double kNorm=DefaultKNorm);
    // operation() returns P(k), not Delta^2(k):
    virtual double operator()(double k, double z) const;
    virtual double DeltaSq(double k, double z) const;
    virtual double dlnPdlnk(double k, double z) const;
    virtual double dlnPdlnG(double k, double z) const;
    virtual double dlnPdz(double k, double z) const;
    // Remember that all partials (except lnG) are for fixed growth factor,
    // with the exception of dlnPdz.
    // Which other parameter id's have partials for this?
    const ParameterVector& dependsOn(double k, double z) const {
      return pv;
    }
    double dlnPdParam(double k, double z, int iParam) const;

    // OmegaM(z) is a 2nd parameter for Smith model:
    double dlnPdOM(double k, double z) const;
    virtual bool has2p() const {return true;}
    virtual double dlnPd2p(double k, double z) const {
      return dlnPdOM(k,z);
    }

    double lnLinearDeltaSq(double k, double z) const;
    double dlnLinearPdParam(double k, double z, int iParam) const;
  private:
    FiducialCosmology& fc;
    const GrowthFunction gf;
    const PowerLawPS lps;
    const PowerLawPS lpsDoMp;
    const PowerLawPS lpsDoMm;
    const PowerLawPS lpsDoBp;
    const PowerLawPS lpsDoBm;
    const PowerLawPS lpsDnsp;
    const PowerLawPS lpsDnsm;
    // Multiply growth function by this to make it equal to a in matter dom.
    double gfnorm;	
    double halofit()  const;
    void haloPrep(double k, double z) const;
    mutable double y;
    mutable double lnkNL;
    mutable double plin;
    mutable double neff;
    mutable double cur;
    mutable double OmegaM;
    ParameterVector pv;
    lnAParam lnAp;
    nsParam nsp;
  };

} // namespace cosmology

#endif // SMITHPS_H
