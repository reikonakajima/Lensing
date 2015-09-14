//
// GAMARandomObjects.h
//
#ifndef GAMARANDOMOBJECTS_H
#define GAMARANDOMOBJECTS_H

#include <iostream>
#include <vector>
#include <CCfits/CCfits>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include "Bounds.h"
#include "LensObjects.h"
using std::vector;


class GAMARandomObjectsError : public MyException {
 public:
 GAMARandomObjectsError(const string& m="") : 
  MyException("GAMARandomObjectsError: "+m) {}
};


/*
 * class GAMARandomObject
 *
 * This is the Random object class, containing minimum information needed to calculate
 * the tangential shear around the random lens point
 *
 */

class LensObject;

class GAMARandomObject : public LensObject {

  public:
  GAMARandomObject() {}
  GAMARandomObject(const string buffer);
  GAMARandomObject(long int _id, double _ra, double _dec, float _z, float _r_comoving,
                   float _logmstar, float _absmag_r, float _uminusr) :
             LensObject(_id, _ra, _dec, _z, -99), r_comoving(_r_comoving), logmstar(_logmstar),
	     absmag_r(_absmag_r), uminusr(_uminusr) {}
  float  getMag() const {
    throw GAMARandomObjectsError("GAMA random objects have no apparent magnitudes (yet)"); }
  float getLogMStar() const { return logmstar; }

 private:
  float r_comoving;
  float logmstar;
  float absmag_r;
  float uminusr;
};


class GAMARandomObjectList : public LensObjectList<GAMARandomObject*> {

 public:
  // Constructor:  If max_count < 0, use the entire input.
  //               Otherwise, randomize input rows and return the first max_count objects
  GAMARandomObjectList(const string fits_filename, int max_count=-1);
  // apply cuts based on logmstar
  int applyLogMStarCut(float min_logmstar, float max_logmstar);

 private:
  // none
};


#endif // RANDOMOBJECTS_H

