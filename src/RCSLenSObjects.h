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
  RCSLenSObject(double ra, double dec, float _mag, float _xpos, float _ypos, float _fwhm_image,
		float _e1_A, float _e2_A, float _e1_B, float _e2_B,
		float _e1_C, float _e2_C, float _e1_D, float _e2_D,
		float sn_ratio, double _wt=1.) :
  SourceObject(ra, dec, _e1_A, _e2_A, _wt), 
      mag(_mag), xpos(_xpos), ypos(_ypos), fwhm(_fwhm_image), sn(sn_ratio), index(0)  {
      shear[0] = Shear(_e1_A, _e2_A);
      shear[1] = Shear(_e1_B, _e2_B);
      shear[2] = Shear(_e1_C, _e2_C);
      shear[3] = Shear(_e1_D, _e2_D);
  }

  void setShearIndex(int i) { index = i; return; }
  int getShearIndex() const { return index; }

  const Shear& getShear() const { checkShearIndex(index); return shear[index]; }
  float getE1() const { checkShearIndex(index); return shear[index].getE1(); }
  float getE2() const { checkShearIndex(index); return shear[index].getE2(); }
  float getESq() const { checkShearIndex(index); return shear[index].getESq(); }

 private:
  float mag;
  float xpos, ypos;
  float fwhm;
  Shear shear[NUM_SHEAR];
  float sn;

  mutable int index;  // to keep track of which shear to return

  void checkShearIndex(int i) const { 
    if (i < 0 || i >= NUM_SHEAR) throw RCSLenSObjectsError("Wrong shear index");
  }
};



class RCSLenSObjectList : public SourceObjectList<RCSLenSObject*> {

 public:
  RCSLenSObjectList(const string fits_filename);

 private:
};

#endif // RCSLENSOBJECTS_H

