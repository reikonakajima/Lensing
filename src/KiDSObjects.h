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

class KiDSObject : public SourceObject {

  friend class KiDSObjectList;

 public:
  
  static const int NUM_SHEAR=3;  // number of blindings

  KiDSObject() {}
  KiDSObject(long int _id, double ra, double dec, float _mag, float _xpos, float _ypos,
	     float _fwhm_image,
	     double _g1_A, double _g2_A, double _g1_B, double _g2_B, double _g1_C, double _g2_C,
	     double _wt_A, double _wt_A, double _wt_B, double _wt_B, double _wt_C, double _wt_C,
	     float sn_ratio, double _zB, int _mask,
	     //float _m_corr, 
	     //float _c1_A, float _c2_A, float _c1_B, float _c2_B, float _c1_C, float _c2_C
	     ) :
  SourceObject(_id, ra, dec, 99., 99., _zB, 0.),  // temporary junk g1, g2, wt in base source object
    mag(_mag), xpos(_xpos), ypos(_ypos), fwhm(_fwhm_image), sn(sn_ratio), mask(_mask) {

    shear[0] = Shear().setG1G2(_g1_A, -_g2_A);  // ra runs in negative direction,
    shear[1] = Shear().setG1G2(_g1_B, -_g2_B);  // so the g2 sign needs to be flipped
    shear[2] = Shear().setG1G2(_g1_C, -_g2_C);

    /*
    m_corr = _m_corr;

    c1_corr[0] = _c1_A;
    c1_corr[1] = _c1_B;
    c1_corr[2] = _c1_C;

    c2_corr[0] = -_c2_A;
    c2_corr[1] = -_c2_B;
    c2_corr[2] = -_c2_C;
    */

    // copy the _A info into "main"  // TODO: make this specifiable via parameter (file)
    this->setShearG1G2BiasCorrections(shear[0], _g1_A, -_g2_A, _wt_A);
    //this->setShearG1G2BiasCorrections(shear[1], _g1_B, -_g2_B, _wt_B);
    //this->setShearG1G2BiasCorrections(shear[2], _g1_C, -_g2_C, _wt_C);
  }

  // do we want to use pixel coordinates, instead of ra/dec?
  void usePixCoord(bool _usePix) { usePixelCoords = _usePix; return; }
  bool isPixCoordUsed() const { return usePixelCoords; }

  double getRA() const { if (usePixelCoords) return xpos; else return SourceObject::ra; }
  double getDec() const { if (usePixelCoords) return ypos; else return SourceObject::dec; }
  float getSNratio() const { return sn; }
  int getMask() const { return mask; }

  // get the array of shapes, c1/2 corrections
  Shear* getShearArray() {return shear;}
  //float* getC1Array() {return c1_corr;}
  //float* getC2Array() {return c2_corr;}

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
  //float m_corr;
  //float c1_corr[NUM_SHEAR];
  //float c2_corr[NUM_SHEAR];
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


class KiDSObjectList : public SourceObjectList<KiDSObject*> {

 public:
  KiDSObjectList(const string fits_filename, int bitmask=0);  // bitmask=0 means *no* masking
  void usePixelCoord(bool _usePix) { source_list[0]->usePixCoord(_usePix); return; }
  void setShearIndex(int _index) {
    checkShearIndex(_index);
    throw KiDSObjectsError("update base SourceObject shears (s,g1,g2) for ABCD selection");
    return;
  }
  // apply mask such that objects with "MASK <= mask_thres" are kept, return number of kept obj.
  int applyMask(int mask_thres=0);  // mask_thres=0 means mask *everything* except MASK==0
  // apply bit mask, so that (MASK & bitmask)>0 objects are excluded
  int applyBitMask(int bitmask);
  // selection source objects according to the blindings A, B, or C
  void setBlinding(char blinding);   // THIS NEEDS DEBUGGING!!!  THE SHEAR ARRAY IS NOT WORKING!!

 private:

  mutable int index;  // to keep track of which shear to return

  // which of the ABC shears should be used?
  void checkShearIndex(int i) const {
    if (i < 0 || i >= KiDSObject::NUM_SHEAR) throw KiDSObjectsError("Wrong shear index");
  }
};

#endif // KIDSOBJECTS_H

