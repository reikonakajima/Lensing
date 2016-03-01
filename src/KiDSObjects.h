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
  enum blind {C, A, B};      // 1=(A), 2=(B), 0=(C)  index shuffling, where (A)(B)(C) are the original
  static char blind_str(blind b) {
    switch(b) {
      case A: return 'A'; break;
      case B: return 'B'; break;
      case C: return 'C'; break;
    }
  }
  static int blind_index;  // to keep track of which shear to return

  KiDSObjectList(const string fits_filename,
		 int bitmask=0,       // bitmask=0 means *no* masking
		 int blind_index=0,
		 valarray<float> pz_full=valarray<float>());
  // selection source objects according to the blindings 0, 1, or 2
  void setBlinding(int _blind_index);
  // apply mask such that objects with "MASK <= mask_thres" are kept, return number of kept obj.
  int applyMask(int mask_thres=0);  // mask_thres=0 means mask *everything* except MASK==0
  // apply bit mask, so that (MASK & bitmask)>0 objects are excluded
  int applyBitMask(int bitmask);

  // read in p(z) from a specz_file
  static valarray<float> getPZ(const string specz_fits_filename, float minz, float maxz, int bitmask);

 private:
  // check that the blinding index is within the correct range
  void checkBlinding(int _blind_index);

  static const int NUM_PZ_ELEM = 70;
  static const float DELTA_Z = 0.05;
  static const float INIT_Z = 0.025;
};


class KiDSObject : public SourceObject {

  friend class KiDSObjectList;

 public:
  
  static const int NUM_SHEAR=3;  // number of blindings

  // constructors
  KiDSObject() {}
  KiDSObject(long int _id, double ra, double dec, float _mag, float _xpos, float _ypos,
	     float _fwhm_image, float sn_ratio, 
	     double _g1_A, double _g2_A, double _g1_B, double _g2_B, double _g1_C, double _g2_C,
	     double _wt_A, double _wt_B, double _wt_C,
	     double _zB, valarray<float> _pz_full,
	     int _mask, int blind_index) :
  // assign id, ra, dec, zB to base SourceObject
  SourceObject(_id, ra, dec, 99., 99., _zB, 0.),  // store junk g1, g2, wt in base source object
    // assign mag, xpos, ypos, fwhm, sn, mask into KiDSObject
    mag(_mag), xpos(_xpos), ypos(_ypos), fwhm(_fwhm_image), sn(sn_ratio), mask(_mask) {

    // assign shear and weights to KiDSObject (with blinding-shuffling!)
    shear[KiDSObjectList::A] = Shear().setG1G2(_g1_A, -_g2_A);  // ra runs in negative direction,
    shear[KiDSObjectList::B] = Shear().setG1G2(_g1_B, -_g2_B);  // so the g2 sign needs to be flipped
    shear[KiDSObjectList::C] = Shear().setG1G2(_g1_C, -_g2_C);

    g1[KiDSObjectList::A] = _g1_A;
    g1[KiDSObjectList::B] = _g1_B;
    g1[KiDSObjectList::C] = _g1_C;

    g2[KiDSObjectList::A] = -_g2_A;
    g2[KiDSObjectList::B] = -_g2_B;
    g2[KiDSObjectList::C] = -_g2_C;

    weight[KiDSObjectList::A] = _wt_A;
    weight[KiDSObjectList::B] = _wt_B;
    weight[KiDSObjectList::C] = _wt_C;

    // assign blinded shears and weights to base SourceObject
    setShearG1G2BiasCorrections(shear[blind_index],
				g1[blind_index], g2[blind_index], weight[blind_index]);

    // assign p(z) to baseSourceObject
    SourceObject::pz = _pz_full;
  }

  // override SourceObject::setWeight() function
  void   setWeight() const {
    throw KiDSObjectsError("setWeight() invalidated fore KiDSObjectList (determined by LensFit)");
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

