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

class SourceObject;

class KiDSObject : public SourceObject {

 public:
  static const int NUM_SHEAR=4;  // number of blindings
  
  KiDSObject() {}
  KiDSObject(long int _id, double ra, double dec, float _mag, float _xpos, float _ypos,
		float _fwhm_image,
		double _e1_A, double _e2_A, double _e1_B, double _e2_B,
		double _e1_C, double _e2_C, double _e1_D, double _e2_D,
		float sn_ratio, double _wt=1.) :
  SourceObject(_id, ra, dec, 99., 99., _wt),  // temporarily fill in e1 and e2 in base source object
    mag(_mag), xpos(_xpos), ypos(_ypos), fwhm(_fwhm_image), sn(sn_ratio) {
    shear[0] = Shear().setG1G2(_e1_A, _e2_A);
    shear[1] = Shear().setG1G2(_e1_B, _e2_B);
    shear[2] = Shear().setG1G2(_e1_C, _e2_C);
    shear[3] = Shear().setG1G2(_e1_D, _e2_D);
    SourceObject::e1 = this->getE1();
    SourceObject::e2 = this->getE2();
  }

  // apply mask such that objects with "MAN_MASK > mask_thres" are excluded
  int applyMask(int mask_thres);

  // do we want to use pixel coordinates, instead of ra/dec?
  void usePixCoord(bool _usePix) { usePixelCoords = _usePix; return; }
  bool isPixCoordUsed() const { return usePixelCoords; }

  double getRA() const { if (usePixelCoords) return xpos; else return SourceObject::ra; }
  double getDec() const { if (usePixelCoords) return ypos; else return SourceObject::dec; }

  // for use with Mesh object: important to override the SourceObject::getX() and getY()!!
  double getX() const { return getRA(); }
  double getY() const { return getDec(); }

  void setShearIndex(int i) { index = i; return; }
  int getShearIndex() const { return index; }

  const Shear& getShear() const { checkShearIndex(index); return shear[index]; }
  double getE1() const { checkShearIndex(index); return shear[index].getE1(); }
  double getE2() const { checkShearIndex(index); return shear[index].getE2(); }
  double getESq() const { checkShearIndex(index); return shear[index].getESq(); }

  float getSNratio() const { return sn; }

  void printLine(ostream& os) const;

  // these static members will be set when the object list is generated
  static bool usePixelCoords; // to keep track of which coordinates we use
  static int index;  // to keep track of which shear to return

 private:
  float mag;
  float xpos, ypos;
  float fwhm;
  Shear shear[NUM_SHEAR];
  float sn;

  void checkShearIndex(int i) const { 
    if (i < 0 || i >= NUM_SHEAR) throw KiDSObjectsError("Wrong shear index");
  }
};


class KiDSObjectList : public SourceObjectList<KiDSObject*> {

 public:
  KiDSObjectList(const string fits_filename);
  void usePixelCoord(bool _usePix) { source_list[0]->usePixCoord(_usePix); return; }
  void setShearIndex(int _index) { source_list[0]->setShearIndex(_index); return; }

 private:
};

#endif // KIDSOBJECTS_H

