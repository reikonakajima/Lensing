// $Id: SmithPS.cpp,v 1.10 2007-11-26 03:30:39 garyb Exp $
#include "SmithPS.h"

using namespace cosmology;

// Smith et al non-linear power spectrum formula
const double domegaM=0.01;
const double domegaB=0.001;
const double dns = 0.02;
const double dlnG = 0.02;
const double dlnk = 0.02;
const double dlna = 0.02;

const double dOM = 0.01;

const double posd = +0.5;
const double negd = -0.5;

// Order of parameter storage:
const int oMIndex=0;
const int nsIndex=1;
const int AsIndex=2;
const int oBIndex=3;
const int paramSize=4;

SmithPS::SmithPS(FiducialCosmology& fc_, double ns, 
		 double norm,
		 double knorm):
  fc(fc_), gf(fc.cosmology()), 
  lps(fc.omegaM(), fc.omegaB(), ns, norm,knorm), 
  lpsDoMp(fc.omegaM()+posd*domegaM, fc.omegaB(), ns, norm,knorm),
  lpsDoMm(fc.omegaM()+negd*domegaM, fc.omegaB(), ns, norm,knorm),
  lpsDoBp(fc.omegaM(), fc.omegaB()+posd*domegaB, ns, norm,knorm),
  lpsDoBm(fc.omegaM(), fc.omegaB()+negd*domegaB, ns, norm,knorm),
  lpsDnsp(fc.omegaM(), fc.omegaB(), ns + posd*dns, norm,knorm),
  lpsDnsm(fc.omegaM(), fc.omegaB(), ns + negd*dns, norm,knorm),
  gfnorm( 1./(1+RecombinationRedshift)/gf(RecombinationRedshift)),
  pv(paramSize) {
  // Build the parameter vector:
  pv[oMIndex] = fc.oMParam();
  pv[nsIndex] = &nsp;

  pv[AsIndex] = &lnAp;
  pv[oBIndex] = fc.oBParam();
}

void
SmithPS::haloPrep(double k, double z) const {
  // Find the growth factor at this redshift
  double growth = gfnorm*gf(z);
  // Find kNL for this redshift
  /*cerr << "haloprep, growth = " << growth << endl;*/
  double lnK = log(k);
  lnkNL = lps.lnkNL(growth);
  y= exp(lnK-lnkNL);
  plin = lps.DeltaSq(lnK) * growth * growth;
  neff = lps.neff(lnkNL);
  cur = lps.C(lnkNL);
  /*cerr << "growth " << growth
	   << " kNL " << exp(lnkNL)
	   << " sigsq " << lps.sigSq(lnkNL)*growth*growth
	   << " plin " << plin
	   << " neff " << neff
	   << " cur " << cur
	   << " sig8 " << lps.sigX(8./0.72)*growth
	   << endl;*/
  // Here is my cheap way of putting a value of Omega_M into
  // Rob's formula for an arbitrary cosmology:
  //***  double q = fc.cosmology().q(z);
  //***  OmegaM = 2.*(1+q)/3.;
  OmegaM = fc.omegaM()*pow(1+z,3.)*pow(fc.H(z),-2.);
}

// Call to the formula
double
SmithPS::DeltaSq(double k, double z) const {
  haloPrep(k,z);
  /*cerr << "kNL " << exp(lnkNL)
	   << " neff " << neff
	   << " cur " << cur
	   << " plin " << plin
	   << " OmegaM " << OmegaM
	   << endl;*/
  return halofit();
}

double 
SmithPS::lnLinearDeltaSq(double k, double z) const {
  double growth = gfnorm*gf(z);
  double lnK = log(k);
  return log(lps.DeltaSq(lnK) * growth * growth);
}

double
SmithPS::dlnLinearPdParam(double k, double z, int iParam) const {
  DVector dp(paramSize);
  double growth = gfnorm*gf(z);
  double lnK = log(k);
  double plin= lps.DeltaSq(lnK);
  if (iParam==oMIndex)
    return log(lpsDoMp.DeltaSq(lnK) / lpsDoMm.DeltaSq(lnK)) / domegaM;
  else if (iParam==oBIndex)
    return log(lpsDoBp.DeltaSq(lnK) / lpsDoBm.DeltaSq(lnK)) / domegaB;
  else if (iParam == nsIndex)
    return log(lpsDnsp.DeltaSq(lnK) / lpsDnsm.DeltaSq(lnK)) / dns;
  else if (iParam == AsIndex)
    return 1.;
  else
    throw MyException("Request for unknown parameter in"
		      " SmithPS::dlnLinearPdParam");
}


// Turn DeltaSq into P(k) for op()
double
SmithPS::operator()(double k, double z) const {
  return DeltaSq(k,z) * pow(k, -3.) * (2. * PI * PI);
}

// Experimental derivative for the "hidden" parameter 
double
SmithPS::dlnPdOM(double k, double z) const {
  haloPrep(k,z);
  double lnK = log(k);
  OmegaM += posd*dOM;
  double lpp = log(halofit());
  OmegaM -= posd*dOM;
  OmegaM += negd*dOM;
  double lpm = log(halofit());
  return (lpp-lpm)/dOM;
}

double
SmithPS::dlnPdlnk(double k, double z) const {
  haloPrep(k,z);
  double lnK = log(k);
  double growth = gfnorm * gf(z);

  y= exp(posd*dlnk + lnK-lnkNL);
  plin = lps.DeltaSq(lnK+posd*dlnk) * growth * growth;
  double lpp=log(halofit());

  y= exp(negd*dlnk + lnK-lnkNL);
  plin = lps.DeltaSq(lnK+negd*dlnk) * growth * growth;
  double lpm=log(halofit());

  return (lpp-lpm) / dlnk
    - 3.; // Remember that halofit gives DeltaSq=k^3 P(k)
}

double
SmithPS::dlnPdlnG(double k, double z) const {
  haloPrep(k,z);

  double growth = gfnorm * gf(z) * exp(posd*dlnG);
  double lnK = log(k);
  lnkNL = lps.lnkNL(growth);
  y = exp(lnK-lnkNL);
  plin = lps.DeltaSq(lnK) * growth * growth;
  neff = lps.neff(lnkNL);
  cur = lps.C(lnkNL);
  double lpp=log(halofit());

  growth = gfnorm * gf(z) * exp(negd*dlnG);
  lnkNL = lps.lnkNL(growth);
  y = exp(lnK-lnkNL);
  plin = lps.DeltaSq(lnK) * growth * growth;
  neff = lps.neff(lnkNL);
  cur = lps.C(lnkNL);
  double lpm=log(halofit());

  return (lpp-lpm) / dlnG;
}

double
SmithPS::dlnPdz(double k, double z) const {
  haloPrep(k,z);
  double dz = posd*dlna*(1+z);
  haloPrep(k,z+dz);
  double lpp=log(halofit());

  dz = negd*dlna*(1+z);
  haloPrep(k,z+dz);
  double lpm=log(halofit());
  return (lpp-lpm) / dz;
}


double
SmithPS::dlnPdParam(double k, double z, int iParam) const {
  if (iParam==AsIndex)
    return 0.5*dlnPdlnG(k,z);

  haloPrep(k,z);
  double lnK = log(k);
  double growth = gfnorm * gf(z);

  if (iParam == oMIndex) {
    lnkNL = lpsDoMp.lnkNL(growth);
    y = exp(lnK-lnkNL);
    plin = lpsDoMp.DeltaSq(lnK) * growth * growth;
    neff = lpsDoMp.neff(lnkNL);
    cur = lpsDoMp.C(lnkNL);
    double lpp = log(halofit());

    lnkNL = lpsDoMm.lnkNL(growth);
    y = exp(lnK-lnkNL);
    plin = lpsDoMm.DeltaSq(lnK) * growth * growth;
    neff = lpsDoMm.neff(lnkNL);
    cur = lpsDoMm.C(lnkNL);
    double lpm = log(halofit());

    return (lpp-lpm) / domegaM;

  } else if (iParam == oBIndex) {
    lnkNL = lpsDoBp.lnkNL(growth);
    y = exp(lnK-lnkNL);
    plin = lpsDoBp.DeltaSq(lnK) * growth * growth;
    neff = lpsDoBp.neff(lnkNL);
    cur = lpsDoBp.C(lnkNL);
    double lpp = log(halofit());

    lnkNL = lpsDoBm.lnkNL(growth);
    y = exp(lnK-lnkNL);
    plin = lpsDoBm.DeltaSq(lnK) * growth * growth;
    neff = lpsDoBm.neff(lnkNL);
    cur = lpsDoBm.C(lnkNL);
    double lpm = log(halofit());

    return (lpp-lpm) / domegaB;

  } else if (iParam == nsIndex) {
    lnkNL = lpsDnsp.lnkNL(growth);
    y = exp(lnK-lnkNL);
    plin = lpsDnsp.DeltaSq(lnK) * growth * growth;
    neff = lpsDnsp.neff(lnkNL);
    cur = lpsDnsp.C(lnkNL);
    double lpp = log(halofit());

    lnkNL = lpsDnsm.lnkNL(growth);
    y = exp(lnK-lnkNL);
    plin = lpsDnsm.DeltaSq(lnK) * growth * growth;
    neff = lpsDnsm.neff(lnkNL);
    cur = lpsDnsm.C(lnkNL);
    double lpm = log(halofit());

    return (lpp-lpm) / dns;
  }
}

// The formula itself: it's a function of y=k/knl, then neff and cur at
// knl, and the linear power.  Cosmological parameter comes in too.
double
SmithPS::halofit() const 
{
  double gam=0.86485+0.2989*neff+0.1631*cur;
  double a=1.4861+1.83693*neff+1.67618*neff*neff+0.7940*neff*neff*neff+
	0.1670756*neff*neff*neff*neff-0.620695*cur;
  a=pow(10.,a);
  double b=pow(10.,0.9463+0.9466*neff+0.3084*neff*neff-0.940*cur);
  double c=pow(10.,-0.2807+0.6669*neff+0.3214*neff*neff-0.0793*cur);
  double xmu=pow(10.,-3.54419+0.19086*neff);
  double xnu=pow(10.,0.95897+1.2857*neff);
  double alpha=1.38848+0.3701*neff-0.1452*neff*neff;
  double beta=0.8291+0.9854*neff+0.3400*neff*neff;

  double f1a=pow(OmegaM,-0.0732);
  double f2a=pow(OmegaM,-0.1423);
  double f3a=pow(OmegaM, 0.0725);
  double f1b=pow(OmegaM,-0.0307);
  double f2b=pow(OmegaM,-0.0585);
  double f3b=pow(OmegaM, 0.0743);

  //  double frac=om_v/(1.-OmegaM);
  const double frac=1.; // I'm assuming a nearly-flat case, since
  // intermediate cases are not well modelled anyway.
  double f1=frac*f1b + (1-frac)*f1a;
  double f2=frac*f2b + (1-frac)*f2a;
  double f3=frac*f3b + (1-frac)*f3a;

  /*cerr << "*Old halofit: lnknl, neff, C: " << lnkNL
	   << " " << neff << " " << cur << endl;
  cerr << "       OM, plin: " << OmegaM << " " << plin
       << endl; //**/

  double ph=a*pow(y,f1*3.)/(1+b*pow(y,f2)+pow(f3*c*y,3-gam));
  ph=ph/(1+xmu/y+xnu/(y*y));
  double pq=plin*pow(1+plin,beta)*exp(-y/4.0-y*y/8.0)/(1+plin*alpha);

  return pq+ph;
}

