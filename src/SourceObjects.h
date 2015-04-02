//
// SourceObjects.h
//
#ifndef SOURCEOBJECTS_H
#define SOURCEOBJECTS_H

#include <iostream>
#include <vector>
#include <valarray>
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
  SourceObject();
  SourceObject(long int _id, double _ra, double _dec, double _s1, double _s2,
	       float _zB, double _wt=1., bool _shapeIsReducedShear=true) :
  id(_id), ra(_ra), dec(_dec), zB(_zB), wt(_wt), responsiv(-99.),
  inputIsReducedShear(_shapeIsReducedShear)
    {
      if (inputIsReducedShear) {
	s.setG1G2(_s1,_s2);
	g1=_s1, g2=_s2;
      }
      else {
	s.setE1E2(_s1,_s2);
	s.getG1G2(g1,g2);
      }
    }

  long int getId() const { return id; }

  double getRA() const { return ra; }
  double getDec() const { return dec; }

  double getE1() const { return s.getE1(); }
  double getE2() const { return s.getE2(); }
  double getESq() const { return s.getESq(); }
  double getG1() const { return g1; }
  double getG2() const { return g2; }
  Shear  getShear() const { return s; }

  float  getRedshift() const { return zB; }
  std::valarray<float>& getPz() {
    if (pz.size()==0)
      throw SourceObjectsError("getPz() not implemented");
    else
      return pz;
  }

  double getWeight() const { return wt; }
  void setWeight(double _wt) const { wt = _wt; return; }
  /*
  double getShapeError() const  { return 2.*shapeerr; }  // ask rachel...
  double getResolution() const  { return res; }
  double getERms() const  { return eRMS; }

  double getLensingWeight() const { if (wt < 0) setLensingWeight(); return wt; }
  void setLensingWeight() const {
    setVars();
    double invshapeweight = (vare + varSN);
    wt = 1. / invshapeweight;
    return;
  }
  */
  double getResponsivity(double et) const {
    if (responsiv < -1.) setResponsivity(et);
    return responsiv;
  }
  void setResponsivity(double et) const {
    /*
    setVars();
    double k0 = varSN*vare / (varSN+vare);
    double k1 = varSN / (varSN+vare);
    k1 *= k1;
    */
    double k0 = 0.;
    double k1 = 1.;
    responsiv = getWeight() * (1. - k0 - k1*et*et);
    return;
  }

  // for use with Mesh object
  double getX() const { return getRA(); }
  double getY() const { return getDec(); }
  double getZ() const { return 0.; }

 protected:
  long int id;
  double ra, dec;  // position
  Shear s;         // shear saves e1,e2 values
  double g1, g2;   // reduced shear
  mutable double wt;         // calculated weight
  bool inputIsReducedShear;
  float zB;             // z_B of photo-z
  std::valarray<float> pz;  // p(z) of photo-z
/*
  double shapeerr;  // shape measurement error
  double eRMS;      // shape noise (calculated from rmag)
  double res;       // resolution

  mutable double vare;     // shapeerr^2
  mutable double varSN;    // eRMS^2
  void setVars() const {
    if (vare < 0 || varSN < 0) {
      vare = shapeerr * shapeerr;
      varSN = eRMS * eRMS;
    }
    return;
  }
*/
  void setVars() const {
  }
  mutable double responsiv;  // responsivity
};



/*
 * class SourceObjectList<>
 *
 * This is the base class for a list of source objects, and contains minimum information
 * needed to calculate a tangential shear signal around any given lens.
 *
 * Note that the constructor is not implemented; this is a virutual base class at the moment.
 */

template <class ObjPtr>
class SourceObjectList {
 public:

  /* Currently, the destructor is not properly implemented (the following will give run-time errors)
   * The proper implementation will keep track of pointer counts, and will delete only when
   * the pointer count goes to zero.
  ~SourceObjectList() {
    for (int i=0; i<size(); ++i) delete source_list[i];
  }
   */

  Bounds<double> getBounds() {
    if (size() <= 0) throw SourceObjectsError("bounds undefined for empty source list");
    if (!bounds) findBounds(); return bounds;
  }
  void findBounds();

  typename vector<ObjPtr>::iterator begin() { return source_list.begin(); }
  typename vector<ObjPtr>::iterator end() { return source_list.end(); }
  int size() { return source_list.size(); }

  void sortByRA(); 
  void sortByDec();

  std::valarray<float>& getPzBins() { 
    if (pzbins.size() == 0)
      throw SourceObjectsError("getPzBins() not implemented");
    else
      return pzbins;
  }

  // for use with Mesh object
  const vector<ObjPtr>& getVectorForm() const { return source_list; }

 protected:
  SourceObjectList() {}

  std::valarray<float> pzbins;      // p(z) redshifts values
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

