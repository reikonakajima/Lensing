//
// StarMaskObjects.h
//
#ifndef STARMASKOBJECTS_H
#define STARMASKOBJECTS_H
#include "LensObjects.h"

class StarMaskObjectsError : public MyException {
 public:
  StarMaskObjectsError(const string& m="") :
  MyException("StarMaskObjectsError: " +m) {}
};


/*
 * class StarMaskObject
 *
 * This class reads a file with star positions (in pixels) as well as halo center (which may be
 * offset from the star position), and sets it up as a "lens" in galaxy-galaxy lensing format.
 * The columns are formatted: mag, x (of halo center), y (of halo center), rhalo, xstar, ystar
 * 
 * [Note: files named <field>.tight.dat have x==xstar, and y==ystar, as they are star mask files.
 * files named <field>.wide.dat have x and y as the halo center, as they are stellar halo masks.
 * The x and y values are stored in the LensObject::RA/Dec variables.]
 * 
 * The purpose is to estimate systematic tangential shear around a star, or about a stellar "halo"
 * (the reflection ghost of a star).
 *
 * For use in systematic test analysis for lensing.
 *
 */

class StarMaskObject : public LensObject {

 public:
    StarMaskObject(const string buffer, int id);
  float getMag() const {return mag;}
  void printLine(ostream& os) const;

 private:
  float mag;
  float xstar;  // star centroid (in pixel units)
  float ystar;
  float rhalo;
};


class StarMaskObjectList : public LensObjectList<StarMaskObject*> {

 public:
  StarMaskObjectList() {}
  StarMaskObjectList(istream& is);

 private:
  // none

};


#endif // LENSOBJECTS_H
