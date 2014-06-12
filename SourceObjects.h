//
// SourceObjects.h
//
#ifndef SOURCEOBJECTS_H
#define SOURCEOBJECTS_H

#include <iostream>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include <list>
#include <vector>
#include "Bounds.h"
using std::list;


class SourceObjectsError : public MyException {
 public:
 SourceObjectsError(const string& m="") : 
  MyException("SourceObjectsError: "+m) {}
};


//
// Base Source Object Class
//
class SourceObject {

 public:
  SourceObject(ifstream& ifs);

  double getRA() const { return ra; }
  double getDec() const { return dec; }
  /*
  float getTrueZ() const { return trueZ; }
  float getSpectroZ() const { return trueZ; }
  float getPhotoZ() const { return trueZ; }  // this is not a joke.  
  float getRMag() const { return rmag; }

  void setTrueZ(float truez_) { trueZ = truez_; }
  */
  Shear getShear() const {return Shear(-e1,e2);}  // parity change
  float getE1() const {return -e1;}               // parity change
  float getE2() const {return e2;}
  float getESq() const {return e1*e1 + e2*e2;}
  float getShapeError() const  { return 2.*shapeerr; }  // ask rachel...
  float getResolution() const  { return resr; }
  float getERms() const  { return eRMS; }

  double getLensingWeight() const { if (wt < 0) setLensingWeight(); return wt; }
  void setLensingWeight() const {
    setVars();
    double invshapeweight = (vare + varSN);
    wt = 1. / invshapeweight;
    return;
  }

  double getResponsivity(double et) { if (responsiv < -1.) setResponsivity(et); return responsiv; }
  void setResponsivity(double et) {
    setVars();
    double k0 = varSN*vare / (varSN+vare);
    double k1 = varSN / (varSN+vare);
    k1 *= k1;
    responsiv = getLensingWeight() * (1. - k0 - k1*et*et);
    return;
  }

  /*
  float getRResolution() const  { return resr; }
  float getIResolution() const  { return resi; }
  float getZLRG() const { return zlrg; }
  bool  isLRG() const { return islrg; }

  void setShapeError(float serr) { shapeerr = serr; }
  void setResolution(float res_r) { resr = res_r; }
  void setIResolution(float res_i) { resi = res_i; }
  void setERms(float e_rms) { eRMS = e_rms; }

  void printLine(ostream& os) const;
  void printLineInBinary(ofstream& ofs) const;
  */
  // for use with Mesh object
  double getX() const { return ra; }
  double getY() const { return dec; }
  double getZ() const { return 0.; }

 protected:
  double ra, dec;  // position
  float e1, e2;    // measured shape
  float shapeerr;  // shape measurement error
  float eRMS;      // shape noise (calculated from rmag)
  mutable double vare;     // shapeerr^2
  mutable double varSN;    // eRMS^2
  void setVars() const {
    if (vare < 0 || varSN < 0) {
      vare = shapeerr * shapeerr;
      varSN = eRMS * eRMS;
    }
    return;
  }
  mutable double wt;         // calculated weight
  mutable double responsiv;  // responsivity

  float rmag;      // r band magnitude

  float zlrg;      // LRG photometric redshift
  bool islrg;      // is this an LRG?
  float resr, resi;   // resolution in r/i bands
  mutable float trueZ;   // redshift, if known...
};




//
// Base Source Object List Class
//
class SourceObjectList : public list<SourceObject*> {
 public:
  SourceObjectList(ifstream& ifs) {};
  void sortByRA(); 
  void sortByDec();
  void setBounds();
  Bounds<double> getBounds() { if (!bounds) setBounds(); return bounds; }
  vector<SourceObject*> getVectorForm();
 protected:
  bool    isOriginalList;
  mutable Bounds<double> bounds;
};




#endif // SOURCEOBJECTS_H

