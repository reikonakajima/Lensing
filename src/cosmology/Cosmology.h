// $Id: Cosmology.h,v 1.17 2006-07-07 18:26:45 garyb Exp $
// Cosmological distances and times
// NOTE:  Because distances & ages are tabulated relative to present day,
//  there may be roundoff error if this routine is used to measure
//  distances/times between early epochs.  Not intended for this use.
//  Could alter the program to use dominant-component analysis, but
//  probably not worth it since the constant-Omega assumptions here will 
//  not be valid before e+-e- annihilation anyway.

#ifndef COSMOLOGY_H
#define COSMOLOGY_H
#include "Std.h"
#include "Simpson.h"
#include "GTable.h"

namespace cosmology {

// An exception class:
class CosmologyError: public MyException {
public:
  CosmologyError(const string &m=""):
    MyException("Cosmology Error: " + m) {};
};

// All distances are in units of c/H0, which is unspecified.  Time
// in units of 1/H0.

class Cosmology {
public:
  Cosmology(double o_m=1., double o_l=0., double ww0=-1., 
	    double ww1=0.,  bool useWa_=true, double o_r=0.):
    Om(o_m), Or(o_r), Ol(o_l), w0(ww0), w1(ww1), func(this),
    useWa(useWa_) {setup();}
  Cosmology(const vector<double>& vecCP);
  // copy constructor and op=
  // Note that all tables, etc., start afresh after copy.
  Cosmology(const Cosmology& rhs): Om(rhs.Om), Or(rhs.Or),
    Ol(rhs.Ol), w0(rhs.w0), w1(rhs.w1), func(this),
    useWa(rhs.useWa) {setup();}
  void operator=(const Cosmology& rhs) {Om=rhs.Om; Or=rhs.Or;
    Ol=rhs.Ol; w0=rhs.w0; w1=rhs.w1; useWa=rhs.useWa; setup();}
  // destructor must kill all the stored integral tables
  ~Cosmology() ;

  enum   Parameter {Omatter=0, Ovacuum, W, Wprime, Oradiation, CPMAX};

  double Dc(double zsrc, double zobs=0.) const ; //comoving distance
  double DA(double zsrc, double zobs=0.) const ; //angular diam dist
  double cDA(double zsrc, double zobs=0.) const ; //comoving angular diam dist
  double DL(double zsrc, double zobs=0.) const ; //luminosity dist
  double Dpm(double zsrc, double zobs=0.) const ; //proper-motion dist
  double Horizon(double zobs=0.) const ; //comoving Horizon radius
  double LookbackTime(double zsrc, double zobs=0.) const ;
  double Age(double zobs=0.) const; //Time since a=0
  double dVdz(double zsrc) const;
  bool Bounce(double *zBounce=0) const;	//Does this bounce?
  double LensShear(double zsrc, double zlens) const; //inv. critical density
  double q0() const { return 0.5*(Om+2.*Or+(1+3*w0)*Ol);}
  double mu(double zsrc) const ; //Distance modulus for fixed rest filter
  double muDeriv(Parameter p, double zsrc) const;  //deriv w.r.t. params
  double DcDeriv(Parameter p, double zsrc, double zobs=0.) const ;
  double H(double z) const;	//Hubble constant, units H0.
  double q(double z) const;	//q(z)= -a''*a/a'*a'
  double OmegaM() const {return Om;}
  double OmegaLambda() const {return Ol;}
  double getw0() const {return w0;}
  double getwa() const {return w1;}
  double usingWa() const {return useWa;}
private:
  double Om;	//Omega_matter, present day
  double Or;	//Omega_radiation
  double Ol;	//Omega_Lambda (or quintessence)
  double w0;
  double w1;
  bool   useWa; // use w_a formula (otherwise w_1)
  void setup();	// configure all of the quantities below

  mutable double Ok;	//Omega in Curvature 
  mutable double invR0;	// |R0|^{-1/2}
  mutable bool doesBounce;	// true if this cosmology bounces
  mutable bool badW1;		// true if w' gives ridiculous high-z 
  mutable double aBounce;	// a value at bounce
  mutable double earlyIndex;	// Index that dominates early behavior
  mutable double earlyOmega;	// Omega for early dominant component
  mutable double earlyA;	// a below which there is dominance
  double RskR(double d) const ; // Return R0 S_k(r/R0).
  double getA(double z) const;  // get a=1/(1+z), check bounds/bounce

  enum FunctionType {DC=0, T, 
		     MassDeriv, VacDeriv, RadDeriv, W0Deriv, W1Deriv,
		     H2, H2Abs, Q};
  // Get one of the cosmology functions from lookup tables:
  double getIntegral(double z, FunctionType t) const; 
  // Stored tables of the various functional integrals
  mutable vector< gtable::Table<>* > integrals;

  // This class will be used as integrand and knows all the key
  // formulae:
  class CosmoFunction {
  public:
    CosmoFunction(const Cosmology* _c): c(_c) {setup();}
    void reset() {setup();}
    void setMode(const Cosmology::FunctionType _t) {ftype=_t;}
    double operator()(double a, Cosmology::FunctionType _t) {
      ftype=_t; return operator()(a);}
    double operator()(double a) const;
  private:
    void setup();
    const Cosmology* c;
    double Om, Ol, w0, w1, Or;
    double ww, Ok;
    bool useWa;
    double expfunc(double a) const;
    double dLogExpfuncdA(double a) const;
    Cosmology::FunctionType ftype;
  };
  friend class CosmoFunction;
  mutable CosmoFunction func;
};

} // namespace cosmology
#endif

