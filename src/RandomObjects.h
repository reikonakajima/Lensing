//
// RandomObjects.h
//
#ifndef RANDOMOBJECTS_H
#define RANDOMOBJECTS_H

#include <iostream>
#include <vector>
#include <CCfits/CCfits>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include "Bounds.h"
#include "LensObjects.h"
using std::vector;


class RandomObjectsError : public MyException {
 public:
 RandomObjectsError(const string& m="") : 
  MyException("RandomObjectsError: "+m) {}
};


/*
 * class RandomObject
 *
 * This is the Random object class, containing minimum information needed to calculate
 * the tangential shear around the random lens point
 *
 */

class LensObject;

class RandomObject : public LensObject {

  public:
  RandomObject() {}
  RandomObject(const string buffer);
  RandomObject(long int _id, double _ra, double _dec, float _xpos, float _ypos, int _mask) : 
  LensObject(_id, _ra, _dec, 0., 0.) {  // redshift and mag set to 0
    xpos = _xpos;
    ypos = _ypos;
    mask = _mask;
  }
  double getRedshift() const {
    return 0.34;  // for M7 objects, the <zlens> is 0.34  TO BE CORRECTED
    //throw RandomObjectsError("random objects have no redshifts");
  }
  float  getMag() const { throw RandomObjectsError("random objects have no magnitudes"); }

 private:
  float xpos, ypos;
  int   mask;

};


class RandomObjectList : public LensObjectList<RandomObject*> {

 public:
  RandomObjectList(const string fits_filename, int max_count=-1);

 private:
  // none
};


#endif // RANDOMOBJECTS_H

