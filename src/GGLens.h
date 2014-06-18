//
// GGLens.h
//
#ifndef GGLENS_H
#define GGLENS_H
#include <iostream>
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
 * This is an LensObject which keeps track of the shear signal by (specified) radial bins.
 *
 * The actual bin range is specified in the GGLensObjectList.
 * The LensObject properties should be inherited from the LensObject.
 *
 */
template<class lensObjPtr>
class GGLensObject {
 public:
  GGLensObject(lensObjPtr lensptr, const int num_radial_bins) :
    lens_ptr(lensptr), tangential_shears(num_radial_bins) {}
  ggLensSum& operator[](int i_rad) {
    if (i_rad < 0 || i_rad >= tangential_shears.size()) {
      throw GGLensError("radial bin index error");
    }
    return tangential_shears[i_rad];
  }
  lensObjPtr getLensPtr() { return lens_ptr; };

 private:
  lensObjPtr lens_ptr;
  vector<ggLensSum> tangential_shears;
};


/*
 * GGLensObjectList
 *
 * This is a list of GGLensObjects.
 * The driver will sort through these objects and then stack the shear signals as desired.
 *
 */
template<class lensObjPtr, class srcObjPtr>
class GGLensObjectList {
 public:
  enum geometry { Flat, SphericalSurface };
  GGLensObjectList() {}  // empty list
  GGLensObjectList(LensObjectList<lensObjPtr> lens_list,
		   SourceObjectList<srcObjPtr> source_list,
		   GenericBins radial_bin,
		   geometry = Flat,             // FIXME!! Change default to SphericalSurface
		   double mesh_frac = 0.);
  int size() { return gglens_object_list.size(); }
  GGLensObject<lensObjPtr>* operator[](int ilens) { return gglens_object_list[ilens]; }
  /*
  void sortByRA();
  void sortByDec();
  Bounds<double> getBounds() { if (!bounds) findBounds(); return bounds;}
  void findBounds();  // finds the bounds of the objects in this list and saves it
  */
 private:
  vector<GGLensObject<lensObjPtr>*> gglens_object_list;
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
