// $Id: ExpansionProjection.h,v 1.9 2006-07-24 22:16:40 garyb Exp $
// Generate derivatives of growth/distances w.r.t. cosmological param set
// *** Some of these derivatives may assume that the fiducial is flat
#ifndef EXPANSIONPROJECTION_H
#define EXPANSIONPROJECTION_H
#include "GrowthFunction.h"
#include "Cosmology.h"
#include "Matrix.h"
#include "AstronomicalConstants.h"

namespace cosmology {
  class ExpansionProjection {
  private:
    // Fiducials:
    const Cosmology& c;
    GrowthFunction g;
  
    DVector delta; // Steps for differentials:
    double h2;	 // h-squared

    // Differentials
    vector<Cosmology*> dcp;
    vector<Cosmology*> dcm;
    vector<GrowthFunction*> dgp;
    vector<GrowthFunction*> dgm;
  public:
    static const int DoM=0;
    static const int DoK=1;
    static const int DoQ=2;
    static const int Dw0=3;
    static const int Dwa=4;
    static const int Dlna=5;
    static const int DSIZE=6;
    ExpansionProjection(const Cosmology& c_,
			double omegaM): c(c_),
      g(c),
      delta(DSIZE),
      h2(omegaM / c.OmegaM()),
      dcp(DSIZE), dcm(DSIZE),
      dgp(DSIZE), dgm(DSIZE) {
      delta[DoM]=0.002;
      delta[DoK]=0.002;
      delta[DoQ]=0.002;
      delta[Dw0]=0.02;
      delta[Dwa]=0.02;
      delta[Dlna]=0.02;

      // ??? will need to change below if radiation is to be in growth

      const double posd=+0.5;
      const double negd=-0.5;
     // oM differentials:
      double deltaOM = delta[DoM] * (1-c.OmegaM()) / h2;
      double deltaOL = delta[DoM] * (-c.OmegaLambda()) / h2;
      dcp[DoM] = new Cosmology(c.OmegaM() + posd*deltaOM,
			      c.OmegaLambda() + posd*deltaOL,
			      c.getw0(), c.getwa(), c.usingWa());
      dcm[DoM] = new Cosmology(c.OmegaM() + negd*deltaOM,
			      c.OmegaLambda() + negd*deltaOL,
			      c.getw0(), c.getwa(), c.usingWa());
      // oK differentials:
      deltaOM = delta[DoK] * (-c.OmegaM()) / h2;
      deltaOL = delta[DoK] * (-c.OmegaLambda()) / h2;
      dcp[DoK] = new Cosmology(c.OmegaM() + posd*deltaOM,
			      c.OmegaLambda() + posd*deltaOL,
			      c.getw0(), c.getwa(), c.usingWa());
      dcm[DoK] = new Cosmology(c.OmegaM() + negd*deltaOM,
			      c.OmegaLambda() + negd*deltaOL,
			      c.getw0(), c.getwa(), c.usingWa());

      // oQ differentials:
      deltaOM = delta[DoQ] * (-c.OmegaM()) / h2;
      deltaOL = delta[DoQ] * (1-c.OmegaLambda()) / h2;
      dcp[DoQ] = new Cosmology(c.OmegaM() + posd*deltaOM,
			      c.OmegaLambda() + posd*deltaOL,
			      c.getw0(), c.getwa(), c.usingWa());
      dcm[DoQ] = new Cosmology(c.OmegaM() + negd*deltaOM,
			      c.OmegaLambda() + negd*deltaOL,
			      c.getw0(), c.getwa(), c.usingWa());

      // w differentials:
      dcp[Dw0] = new Cosmology(c.OmegaM(), c.OmegaLambda(),
			      c.getw0() + posd*delta[Dw0], 
			       c.getwa(), c.usingWa());
      dcm[Dw0] = new Cosmology(c.OmegaM(), c.OmegaLambda(),
			      c.getw0() + negd*delta[Dw0], 
			       c.getwa(), c.usingWa());

      dcp[Dwa] = new Cosmology(c.OmegaM(), c.OmegaLambda(),
			      c.getw0(), c.getwa() + posd*delta[Dwa], 
			       c.usingWa());
      dcm[Dwa] = new Cosmology(c.OmegaM(), c.OmegaLambda(),
			      c.getw0(), c.getwa() + negd*delta[Dwa], 
			      c.usingWa());

      for (int i=0; i<DSIZE-1; i++) {
	dgp[i] = new GrowthFunction(*dcp[i]);
	dgm[i] = new GrowthFunction(*dcm[i]);
      }
    }
    ~ExpansionProjection() {
      for (int i=0; i<DSIZE; i++) {
	if (i==Dlna) continue;
	delete dcp[i]; delete dgp[i];
	delete dcm[i]; delete dgm[i];}
    }

    DVector dlnDdp(double z) {
      DVector ret(DSIZE);
      for (int i=0; i<DSIZE; i++)
	if (i==Dlna) {
	  ret[i] = -(1+z)/(c.H(z)*c.cDA(z)); // Flat assumed, needs C_k
	                                     // more generally.  ???
	} else {
	  ret[i] = (log(dcp[i]->cDA(z)) - log(dcm[i]->cDA(z))) / delta[i];
	}
      ret[DoM] -= 0.5/h2; // value of h is altered by omegaM, omegaK
      ret[DoK] -= 0.5/h2;
      ret[DoQ] -= 0.5/h2;
      return ret;
    }
    DVector dlnGdp(double z) {
      DVector ret(DSIZE);
      const double lnaRecombination = -log(1+RecombinationRedshift);
      const double lna=-log(1+z);
      for (int i=0; i<DSIZE; i++)
	if (i==Dlna) {
	  ret[i] = g.dlnGdlnA(lna);
	} else {
	  ret[i] = ( dgp[i]->lnG(lna) - dgp[i]->lnG(lnaRecombination) -
		     (dgm[i]->lnG(lna) - dgm[i]->lnG(lnaRecombination)) )
	    / delta[i];
	}
      return ret;
    }
    DVector dlnHdp(double z) {
      DVector ret(DSIZE);
      for (int i=0; i<DSIZE; i++)
	if (i==Dlna) {
	  double a1 = exp(delta[i])/(1+z);
	  double z1 = 1./a1 - 1.;
	  ret[i] = (log(c.H(z1)) - log(c.H(z))) / delta[i];
	} else {
	  ret[i] = (log(dcp[i]->H(z)) - log(dcm[i]->H(z))) / delta[i];
	}
      ret[DoM] += 0.5/h2; // value of h is altered by omegaM, omegaK
      ret[DoK] += 0.5/h2;
      ret[DoQ] += 0.5/h2;
      return ret;
    }

    // Derivative of OmegaM(z) w.r.t. parameters
    DVector dOMdp(double z) {
      DVector ret(DSIZE);
      for (int i=0; i<DSIZE; i++)
	if (i==Dlna) {
	  double a1 = exp(delta[i])/(1+z);
	  double z1 = 1./a1 - 1.;
	  ret[i] = c.OmegaM()*
	    (pow(1+z1,3.)*pow(c.H(z1),-2.) -
	     pow(1+z,3.)*pow(c.H(z),-2.) ) / delta[i];
	} else {
	  ret[i] = pow(1+z,3.)*
	    (dcp[i]->OmegaM()*pow(dcp[i]->H(z),-2.) -
	     dcm[i]->OmegaM()*pow(dcm[i]->H(z),-2.) ) / delta[i];
	}
      return ret;
    }
  };

} // namespace cosmology

#endif // EXPANSIONPROJECTION_H

  
