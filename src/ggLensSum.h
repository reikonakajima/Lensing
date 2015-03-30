#ifndef GGLENSSUM_H
#define GGLENSSUM_H
#include <iostream>
#include "Std.h"
//#include "RedshiftObjects.h"
//#include "LensObjects.h"
//#include "SourceObjects.h"
#include "Bounds.h"

//
// structure to keep photcat signal tally
//
class ggLensSum {
 public:
 ggLensSum() :
  lenscounts(0), paircounts(0), weights(0.), wsq(0.), DeltaSigma_t(0.), DeltaSigma_s(0.),
    responsivity(0.), variance_t(0.), variance_s(0.), invsigmacrit(0.), lum(0.), z(0.),
    zwidth(0.), stmass(0.), amag(0.) {}
  inline void addLensCounts(int lcount = 1) { lenscounts += lcount; return; }
  inline void addPairCounts(int pcount = 1) { paircounts += pcount; return; }
  inline void addWeight(double weight) { weights += weight; return; }
  inline void addWeightSq(double _wsq) { wsq += _wsq; return; }
  inline void addResponsivity(double resp) { responsivity += resp; return; }
  inline void addDeltaSigma_t(double delsig) { DeltaSigma_t += delsig; return; }
  inline void addDeltaSigma_s(double delsig) { DeltaSigma_s += delsig; return; }
  inline void addVariance_t(double var) { variance_t += var; return; }
  inline void addVariance_s(double var) { variance_s += var; return; }
  inline void addInvSigmaCrit(double invsc) {invsigmacrit += invsc; return; }
  inline void addLum(double lum_) { lum += lum_; return; }
  inline void addZ(double z_) { z += z_; return; }
  inline void addZWidth(double zw) { zwidth += zw; return; }
  inline void addStellarMass(double sm) { stmass += sm; return; }
  inline void addApparentMag(double am) { amag += am; return; }
  inline int    getLensCounts() const { return lenscounts; }
  inline int    getPairCounts() const { return paircounts; }
  inline double getWeights() const { return weights; }
  inline double getWeightSq() const { return wsq; }
  inline double getResponsivity() const { return responsivity; }
  inline double getDeltaSigma_t() const { return DeltaSigma_t; }
  inline double getDeltaSigma_s() const { return DeltaSigma_s; }
  inline double getVariance_t() const { return variance_t; }
  inline double getVariance_s() const { return variance_s; }
  inline double getInvSigmaCrit() const { return invsigmacrit; }
  inline double getLum() const { return lum; }
  inline double getZ() const { return z; }
  inline double getZWidth() const { return zwidth; }
  inline double getStellarMass() const { return stmass; }
  inline double getApparentMag() const { return amag; }
 private:
  int    lenscounts;   // the number of lens contributing to this sum
  int    paircounts;   // numb   (pair counts)
  double weights;      // wtotw  (total weight)
  double wsq;          // w^2
  double responsivity; // sshup  (sum for shear responsivity)
  double DeltaSigma_t; // sigup  (weighted tangent shear -- weight should include SigmaCrit)
  double DeltaSigma_s; //        (weighted cross shear -- weight should include SigmaCrit)
  double variance_t;   // variance, tangential
  double variance_s;   // variance, cross
  double invsigmacrit; // sumscinv2
  double lum;          // luminosity
  double z;            // redshift (median or otherwise)
  double zwidth;       // redshift uncertainty
  double stmass;       // stellar mass
  double amag;         // apparent magnitude
};


class LensCounts {
 public:
  LensCounts() : lenscounts(0) {}
  inline void addLensCounts(int lcount = 1) { lenscounts += lcount; return; }
  inline int  getLensCounts() const { return lenscounts; }
 private:
  int lenscounts;   //        (lens counts)
};


#endif // GGLENSSUM_H
