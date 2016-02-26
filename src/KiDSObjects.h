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
 *  *** Upon generation, lensfit weight == 0 objects are automatically excluded from this list.
 *
 */

class KiDSObject;

class KiDSObjectList : public SourceObjectList<KiDSObject*> {

 friend class KiDSObject;

 public:
  enum blind {Z, X, Y};       // Z=(A), X=(B), Y=(C)  index shuffling, where (A)(B)(C) are the original
  static blind blind_index;  // to keep track of which shear to return

  KiDSObjectList(const string fits_filename, int bitmask=0,  // bitmask=0 means *no* masking
		 KiDSObjectList::blind blind_index=KiDSObjectList::X);
  // selection source objects according to the blindings 0, 1, or 2
  void setBlinding(KiDSObjectList::blind _blind_index);
  // apply mask such that objects with "MASK <= mask_thres" are kept, return number of kept obj.
  int applyMask(int mask_thres=0);  // mask_thres=0 means mask *everything* except MASK==0
  // apply bit mask, so that (MASK & bitmask)>0 objects are excluded
  int applyBitMask(int bitmask);
};


class KiDSObject : public SourceObject {

  friend class KiDSObjectList;

 public:
  
  static const int NUM_SHEAR=3;  // number of blindings

  // constructors
  KiDSObject() {}
  KiDSObject(long int _id, double ra, double dec, float _mag, float _xpos, float _ypos,
	     float _fwhm_image,
	     double _g1_A, double _g2_A, double _g1_B, double _g2_B, double _g1_C, double _g2_C,
	     double _wt_A, double _wt_B, double _wt_C,
	     float sn_ratio, double _zB, int _mask, int blind_index) :
  SourceObject(_id, ra, dec, 99., 99., _zB, 0.),  // store junk g1, g2, wt in base source object
    mag(_mag), xpos(_xpos), ypos(_ypos), fwhm(_fwhm_image), sn(sn_ratio), mask(_mask) {

    shear[0] = Shear().setG1G2(_g1_A, -_g2_A);  // ra runs in negative direction,
    shear[1] = Shear().setG1G2(_g1_B, -_g2_B);  // so the g2 sign needs to be flipped
    shear[2] = Shear().setG1G2(_g1_C, -_g2_C);

    g1[0] = _g1_A;
    g1[1] = _g1_B;
    g1[2] = _g1_C;

    g2[0] = -_g2_A;
    g2[1] = -_g2_B;
    g2[2] = -_g2_C;

    weight[0] = _wt_A;
    weight[1] = _wt_B;
    weight[2] = _wt_C;

    // assign blinded shears and weights to SourceObject
    setShearG1G2BiasCorrections(shear[blind_index],
				g1[blind_index], g2[blind_index], weight[blind_index]);
  }

  // override SourceObject::setWeight() function
  void   setWeight() const {
    throw KiDSObjectsError("setWeight() invalidated, use KiDSObjectList::setBlinding() instead");
  }

  // ordinary return value functions
  float getSNratio() const { return sn; }
  int   getMask() const { return mask; }

  // get the array of shapes
  Shear*  getShearArray() {return shear;}
  double* getG1Array() {return g1;}
  double* getG2Array() {return g2;}
  double* getWeightArray() {return weight;}

  // for use with Mesh object: important to override the SourceObject::getX() and getY()!!
  double getX() const { return getRA(); }
  double getY() const { return getDec(); }

  void printLine(ostream& os) const;

 private:
  float mag;
  float xpos, ypos;
  float fwhm;
  Shear shear[NUM_SHEAR];
  double g1[NUM_SHEAR];  // because it's faster to save these values...
  double g2[NUM_SHEAR];  // ...than to have to recalculate them
  double weight[NUM_SHEAR];
  float sn;            // Signal-to-noise
  int   mask;

  // make sure the shear is set in the base SourceObject class
  void setShearG1G2BiasCorrections(Shear s_new, double _g1, double _g2, double _wt) {
    SourceObject::s = s_new;
    SourceObject::g1 = _g1;
    SourceObject::g2 = _g2;
    SourceObject::wt = _wt;
  }
};


#endif // KIDSOBJECTS_H

