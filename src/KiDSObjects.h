//
// KiDSObjects.h
//
#ifndef KIDSOBJECTS_H
#define KIDSOBJECTS_H

#include <iostream>
#include <vector>
#include <CCfits/CCfits>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include "Bounds.h"
#include "SourceObjects.h"
using std::vector;


class KiDSObjectsError : public MyException {
 public:
 KiDSObjectsError(const string& m="") : 
  MyException("KiDSObjectsError: "+m) {}
};


/*
 * class KiDSObject
 *
 * This is the KiDS object class, containing minimum information needed to calculate
 * the tangential shear around the lens for systematic test purposes (i.e., no redshift info). 
 *
 *
 * class KiDSObjectList
 *
 * Upon generation, lensfit weight == 0 objects are automatically excluded from this list.
 *
 *
 * TODO:
 *   read Massimo's paper, see if any additional KiDS catalog info is necessary
 *   add applyMask() function
 *   correct KiDS catalog entries (such as magnitudes)
 *   add photometric redshift information
 *
 */

class KiDSObject : public SourceObject {

 public:
  
  static const int NUM_SHEAR=4;  // number of blindings

  KiDSObject() {}
  KiDSObject(long int _id, double ra, double dec, float _mag, float _xpos, float _ypos,
	     float _fwhm_image,
	     double _g1_A, double _g2_A, double _g1_B, double _g2_B,
	     double _g1_C, double _g2_C, double _g1_D, double _g2_D,
	     float sn_ratio, double _zB, valarray<float> _pz_full, int _mask, double _wt) :
  SourceObject(_id, ra, dec, 99., 99., _wt),  // temporarily fill in g1 and g2 in base source object
    mag(_mag), xpos(_xpos), ypos(_ypos), fwhm(_fwhm_image), sn(sn_ratio),
    zB(_zB), pz(_pz_full), mask(_mask) {
    shear[0] = Shear().setG1G2(_g1_A, -_g2_A);  // ra runs in negative direction,
    shear[1] = Shear().setG1G2(_g1_B, -_g2_B);  // lensfit flips g2 sign
    shear[2] = Shear().setG1G2(_g1_C, -_g2_C);  // presumably because it fits in pixel space
    shear[3] = Shear().setG1G2(_g1_D, -_g2_D);
    this->setShearAndG1G2(shear[0], _g1_A, -_g2_A);  // copy the _A shear into "main shape"
  }

  // do we want to use pixel coordinates, instead of ra/dec?
  void usePixCoord(bool _usePix) { usePixelCoords = _usePix; return; }
  bool isPixCoordUsed() const { return usePixelCoords; }

  double getRA() const { if (usePixelCoords) return xpos; else return SourceObject::ra; }
  double getDec() const { if (usePixelCoords) return ypos; else return SourceObject::dec; }
  float getSNratio() const { return sn; }
  int getMask() const { return mask; }

  float getRedshift() const { return zB;}
  valarray<float>& getPz() { return pz; }

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
  float sn;            // Signal-to-noise
  float zB;             // z_B of photo-z
  valarray<float> pz;  // p(z) of photo-z
  int   mask;

  // make sure the shear is set in the base SourceObject class
  void setShearAndG1G2(Shear s_new, double _g1, double _g2) {
    SourceObject::s = s_new;
    SourceObject::g1 = _g1;
    SourceObject::g2 = _g2;
  }
};


class KiDSObjectList : public SourceObjectList<KiDSObject*> {

 public:
  KiDSObjectList(const string fits_filename);
  void usePixelCoord(bool _usePix) { source_list[0]->usePixCoord(_usePix); return; }
  void setShearIndex(int _index) {
    checkShearIndex(_index);
    throw KiDSObjectsError("update base SourceObject shears (s,g1,g2) for ABCD selection");
    return;
  }
  // apply mask such that objects with "MAN_MASK <= mask_thres" are kept, return number of kept obj.
  int applyMask(int mask_thres=0);
  // apply bit mask, so that (MAN_MASK & bitmask)>0 objects are excluded
  int applyBitMask(int bitmask);

  // return the Pz redshift bins corresponding to the p(z) that KiDSObject returns
  valarray<float>& getPzBins() { return pzbins; }

  static const int NUM_PZ_ELEM = 70;  // parameters to set p(z) bins
  static const float DELTA_Z = 0.05;

 private:

  mutable int index;  // to keep track of which shear to return
  valarray<float> pzbins;      // p(z) redshifts values

  // which of the ABCD shears should be used?
  void checkShearIndex(int i) const {
    if (i < 0 || i >= KiDSObject::NUM_SHEAR) throw KiDSObjectsError("Wrong shear index");
  }
};

#endif // KIDSOBJECTS_H

