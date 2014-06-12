//
// GGLens.h
//
#ifndef GGLENS_H
#define GGLENS_H
#include <iostream>
#include <list>
#include <vector>
#include "Bounds.h"
#include "Bins.h"
#include "ggLensSum.h"
#include "LensObjects.h"
#include "SourceObjects.h"

using std::list;
using std::istringstream;


class GGLensError : public MyException {
 public:
 GGLensError(const string& m="") :
  MyException("GGLensError: " +m) {}
};

/*
 * GGLensObject
 *
 * This is an AstronomicalObject which keeps track of the shear signal around itself in radial bins.
 * The actual bin range is not specified here; it is specified in the GGLensObjectList.
 * The AstronomicalObject properties should be inherited from AstronomicalObject.   FIXME!! :  TODO
 */
class GGLensObject {
 public:
  GGLensObject(const int num_radial_bins);
  ggLensSum& operator()(int i_rad) {
    if (i_rad < 0 || i_rad >= tangential_shears.size()) {
      throw GGLensError("radial bin index error");
    }
    return tangential_shears[i_rad];
  }

 private:
  vector<ggLensSum> tangential_shears;
  string id;
  double ra, dec;
  float mag;
};


/*
 * GGLensObjectList
 *
 * This is a list of GGLensObjects.
 * The driver will sort through these objects and then stack the shear signals as desired.
 *
 */
class GGLensObjectList {
 public:
  enum geometry { Flat, SphericalSurface };
  GGLensObjectList() {}  // empty list
  GGLensObjectList(LensObjectList lens_list,
		   SourceObjectList source_list,  // FIXME!!  with sourceObjectList
		   GenericBins radial_bin,
		   geometry = Flat,             // FIXME!! Change default to SphericalSurface
		   double mesh_size = 30.);     // FIXME!! mesh_size will default to width / 100
  int size() { return gglens_object_list.size(); }
  /*
  void sortByRA();
  void sortByDec();
  Bounds<double> getBounds() { if (!bounds) findBounds(); return bounds;}
  void findBounds();  // finds the bounds of the objects in this list and saves it
  */
 private:
  list<GGLensObject*> gglens_object_list;
  geometry geom;      // geometry as to detemine distance between two points
  double mesh_size;
  GenericBins radial_bin;  // stores the radial bin edges

  /*
  Bounds<double> bounds;  // run findBounds() to set value
  LensObjectList::iterator searchRA(LensObjectList::iterator first, 
				    LensObjectList::iterator last,
				    const double ra);
  LensObjectList::iterator searchDec(LensObjectList::iterator first, 
				     LensObjectList::iterator last,
				     const double dec);
  */
};


#endif // GGLENS_H
