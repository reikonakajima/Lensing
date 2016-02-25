// $Id: FiducialCosmology.h,v 1.4 2006-08-07 13:22:00 garyb Exp $
// Subroutines for lensing cross-correlation study
#ifndef FIDUCIALCOSMOLOGY_H
#define FIDUCIALCOSMOLOGY_H

#include <sstream>
#include "Std.h"
#include "FisherMatrices.h"
#include "Cosmology.h"

namespace cosmology {

// I'm going to assume flat fiducial
// This differs from Cosmology in that units are absolute, not H-relative.
class FiducialCosmology {
private:
  const double h;
  const double oM;
  const double oB;
  Cosmology c;
  omegaMParam oMP;
  omegaKParam oKP;
  omegaBParam oBP;
public:
  FiducialCosmology(double omegaM=0.146, double omegaB=0.024, double h_=0.725,
		    double w0=-1., double wa=0.,
		    bool useWa=true,
		    double omegaQ=-1.):
    h(h_), oM(omegaM), oB(omegaB),
    c(omegaM/h/h, 
      omegaQ>=0. ? omegaQ/h/h : 1.-omegaM/h/h, 
      w0, wa, useWa) {}
    // with radiation???  , 5e-5/h/h) {}
  const Cosmology& cosmology() const {return c;}
  double omegaM() const {return oM;}
  double omegaB() const {return oB;}
  double cDA(double z) const {return c.cDA(z)/h;}
  double cDA(double z1, double z2) const {return c.cDA(z1,z2)/h;}
  double Dc(double z) const {return c.Dc(z)/h;}
  double H(double z) const {return c.H(z)*h;}
  double OmegaM(double z) const {
    return oM*pow(1+z,3.)*pow(H(z),-2.);
  }
  FisherParameter* oMParam() {return &oMP;}
  FisherParameter* oKParam() {return &oKP;}
  FisherParameter* oBParam() {return &oBP;}
  ParameterVector getPV() {
    ParameterVector pv; 
    pv.push_back(&oMP);
    pv.push_back(&oKP);
    pv.push_back(&oBP);
    return pv;
  }
};

} // namespace cosmology

#endif // FIDUCIALCOSMOLOGY_H
