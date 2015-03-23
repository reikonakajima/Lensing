//
// GAMAObjects.h
//
#ifndef GAMAOBJECTS_H
#define GAMAOBJECTS_H
#include "LensObjects.h"

class GAMAObjectsError : public MyException {
 public:
  GAMAObjectsError(const string& m="") :
  MyException("GAMAObjectsError: " +m) {}
};


/*
 * class GAMAObject
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

class GAMAObject : public LensObject {

 public:
  GAMAObject(const string buffer, int id);
  float getMag() const { return mag; }
  int   getType() const { return type; }
  void printLine(ostream& os) const;

 private:
  int   type;
  float mag;
  float xstar;  // star centroid (in pixel units)
  float ystar;
  float rhalo;
};


class GAMAObjectList : public LensObjectList<GAMAObject*> {

 public:
  GAMAObjectList(istream& is);

 private:
  // none

};


#endif // GAMAOBJECTS_H
