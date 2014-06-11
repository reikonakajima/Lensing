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
  SourceObject() {} 

  double getRA() const { return ra; }
  double getDec() const { return dec; }
  /*
  float getTrueZ() const { return trueZ; }
  float getSpectroZ() const { return trueZ; }
  float getPhotoZ() const { return trueZ; }  // this is not a joke.  
  float getRMag() const { return rmag; }

  void setTrueZ(float truez_) { trueZ = truez_; }

  Shear getShear() const {return Shear(-e1,e2);}  // parity change
  float getE1() const {return -e1;}               // parity change
  float getE2() const {return e2;}
  float getESq() const {return e1*e1 + e2*e2;}
  float getShapeError() const  { return 2.*shapeerr; }  // ask rachel...
  float getRResolution() const  { return resr; }
  float getIResolution() const  { return resi; }
  float getERms() const  { return eRMS; }
  float getZLRG() const { return zlrg; }
  bool  isLRG() const { return islrg; }
  */
  void setShapeError(float serr) { shapeerr = serr; }
  void setRResolution(float res_r) { resr = res_r; }
  void setIResolution(float res_i) { resi = res_i; }
  void setERms(float e_rms) { eRMS = e_rms; }
  /*
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

