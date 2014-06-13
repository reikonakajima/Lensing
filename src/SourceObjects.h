//
// SourceObjects.h
//
#ifndef SOURCEOBJECTS_H
#define SOURCEOBJECTS_H

#include <iostream>
#include <vector>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include "Bounds.h"
using std::vector;


class SourceObjectsError : public MyException {
 public:
 SourceObjectsError(const string& m="") : 
  MyException("SourceObjectsError: "+m) {}
};


/*
 * class SourceObject
 *
 * This is the base class for source objects, and contains minimum information needed to calculate
 * lensing signal around the lens.
 *
 */

class SourceObject {

 public:
  SourceObject(ifstream& ifs);

  double getRA() const { return ra; }
  double getDec() const { return dec; }

  Shear getShear() const {return Shear(e1,e2);}  // parity change
  float getE1() const {return e1;}               // parity change
  float getE2() const {return e2;}
  float getESq() const {return e1*e1 + e2*e2;}
  float getShapeError() const  { return 2.*shapeerr; }  // ask rachel...
  float getResolution() const  { return res; }
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

  // for use with Mesh object
  double getX() const { return ra; }
  double getY() const { return dec; }
  double getZ() const { return 0.; }

 protected:
  double ra, dec;  // position
  float e1, e2;    // measured shape
  float shapeerr;  // shape measurement error
  float eRMS;      // shape noise (calculated from rmag)
  float res;       // resolution

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
};



/*
 * class SourceObjectList<>
 *
 * This is the base class for a list of source objects, and contains minimum information
 * needed to calculate a tangential shear signal around any given lens.
 *
 */

template <class ObjPtr>
class SourceObjectList {
 public:
  SourceObjectList(ifstream& ifs) {};

  Bounds<double> getBounds() { if (!bounds) findBounds(); return bounds; }
  void findBounds();

  typename vector<ObjPtr>::iterator begin() { return source_list.begin(); }
  typename vector<ObjPtr>::iterator end() { return source_list.end(); }
  int size() { return source_list.size(); }

  void sortByRA(); 
  void sortByDec();

  // for use with Mesh object
  const vector<ObjPtr>& getVectorForm() const { return source_list; }

 protected:
  SourceObjectList() {}

  vector<ObjPtr> source_list;
  mutable Bounds<double> bounds;

  static bool Compare_Source_RA(ObjPtr lhs, ObjPtr rhs) {
      return lhs->getRA() < rhs->getRA(); // sort in increasing order
  }
  static bool Compare_Source_Dec(ObjPtr lhs, ObjPtr rhs) {
      return lhs->getDec() < rhs->getDec(); // sort in increasing order
  }

};




#endif // SOURCEOBJECTS_H

