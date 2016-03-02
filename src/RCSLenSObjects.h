//
// RCSLenSObjects.h
//
#ifndef RCSLENSOBJECTS_H
#define RCSLENSOBJECTS_H

#include <iostream>
#include <vector>
#include <CCfits/CCfits>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include "Bounds.h"
#include "SourceObjects.h"
using std::vector;


class RCSLenSObjectsError : public MyException {
 public:
 RCSLenSObjectsError(const string& m="") : 
  MyException("RCSLenSObjectsError: "+m) {}
};

/*
 * class RCSLenSObject
 *
 * This is the RCSLenS object class, containing minimum information needed to calculate
 * the tangential shear around the lens for systematic test purposes (i.e., no redshift info). 
 *
 */

class RCSLenSObject : public SourceObject {

 public:
  
  static const int NUM_SHEAR=4;   // number of blindings

  RCSLenSObject() {}
  RCSLenSObject(long int _id, double ra, double dec, float _mag, float _xpos, float _ypos,
		float _fwhm_image,
		double _g1_A, double _g2_A, double _g1_B, double _g2_B,
		double _g1_C, double _g2_C, double _g1_D, double _g2_D,
		float sn_ratio, double _wt=1.) :
  SourceObject(_id, ra, dec, 99., 99., _wt),
    mag(_mag), xpos(_xpos), ypos(_ypos), fwhm(_fwhm_image), sn(sn_ratio) {
    shear[0] = Shear().setG1G2(_g1_A, -_g2_A);  // lensfit flips g2 sign
    shear[1] = Shear().setG1G2(_g1_B, -_g2_B);
    shear[2] = Shear().setG1G2(_g1_C, -_g2_C);
    shear[3] = Shear().setG1G2(_g1_D, -_g2_D);
    this->setShearAndG1G2(shear[0], _g1_A, -_g2_A);  // copy the _A shear into "main shape"
  }

  void usePixCoord(bool _usePix) { usePixelCoords = _usePix; return; }
  bool isPixCoordUsed() const { return usePixelCoords; }

  double getRA() const { if (usePixelCoords) return xpos; else return SourceObject::ra; }
  double getDec() const { if (usePixelCoords) return ypos; else return SourceObject::dec; }
  float getSNratio() const { return sn; }

  // for use with Mesh object: important to override the SourceObject::getX() and getY()!!
  double getX() const { return getRA(); }
  double getY() const { return getDec(); }

  void printLine(ostream& os) const;

  // these static members will be set when the object list is generated
  static bool usePixelCoords; // to keep track of which coordinates we use

 private:
  float mag;
  float xpos, ypos;
  float fwhm;
  Shear shear[NUM_SHEAR];
  float sn;

  // make sure the shear is set in the base SourceObject class
  void setShearAndG1G2(Shear s_new, double _g1, double _g2) {
    SourceObject::s = s_new;
    SourceObject::g1 = _g1;
    SourceObject::g2 = _g2;
  }

};


class RCSLenSObjectList : public SourceObjectList<RCSLenSObject*> {

 public:
  RCSLenSObjectList(const string fits_filename);
  void usePixelCoord(bool _usePix) { source_list[0]->usePixCoord(_usePix); return; }
  void setShearIndex(int _index) {
    checkShearIndex(_index);
    throw RCSLenSObjectsError("update base SourceObject shears (s,g1,g2) for ABCD selection");
    return;
  }
  int getShearIndex() const { return index; }

 private:
  mutable int index;  // to keep track of which shear to return

  // which of the ABCD shears should be used?
  void checkShearIndex(int i) const {
    if (i < 0 || i >= RCSLenSObject::NUM_SHEAR) throw RCSLenSObjectsError("Wrong shear index");
  }
};

#endif // RCSLENSOBJECTS_H

