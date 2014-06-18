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

class SourceObject;

class RCSLenSObject : public SourceObject {

 public:
  static const int NUM_SHEAR=4;
  
  RCSLenSObject() {}
  RCSLenSObject(long int _id, double ra, double dec, float _mag, float _xpos, float _ypos,
		float _fwhm_image,
		float _e1_A, float _e2_A, float _e1_B, float _e2_B,
		float _e1_C, float _e2_C, float _e1_D, float _e2_D,
		float sn_ratio, double _wt=1.) :
  SourceObject(_id, ra, dec, _e1_A, _e2_A, _wt),
      mag(_mag), xpos(_xpos), ypos(_ypos), fwhm(_fwhm_image), sn(sn_ratio) {
      shear[0] = Shear(_e1_A, _e2_A);
      shear[1] = Shear(_e1_B, _e2_B);
      shear[2] = Shear(_e1_C, _e2_C);
      shear[3] = Shear(_e1_D, _e2_D);
  }


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
  float getE1() const { checkShearIndex(index); return shear[index].getE1(); }
  float getE2() const { checkShearIndex(index); return shear[index].getE2(); }
  float getESq() const { checkShearIndex(index); return shear[index].getESq(); }

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
    if (i < 0 || i >= NUM_SHEAR) throw RCSLenSObjectsError("Wrong shear index");
  }
};


class RCSLenSObjectList : public SourceObjectList<RCSLenSObject*> {

 public:
  RCSLenSObjectList(const string fits_filename);
  void usePixelCoord(bool _usePix) { source_list[0]->usePixCoord(_usePix); return; }
  void setShearIndex(int _index) { source_list[0]->setShearIndex(_index); return; }

 private:
};

#endif // RCSLENSOBJECTS_H

