//
// TreeCorrObjects.h
//
#ifndef TREECORROBJECTS_H
#define TREECORROBJECTS_H

#include <iostream>
#include <vector>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include "Bounds.h"
#include "GGLensData.h"
using std::vector;


class TreeCorrObjectsError : public MyException {
 public:
 TreeCorrObjectsError(const string& m="") : 
  MyException("TreeCorrObjectsError: "+m) {}
};


/*
 * class TreeCorrNGObject
 *
 * This is the TreeCorr object class, containing classes that can ingest TreeCorr output files
 * for NGCorrelation (galaxy-galaxy lensing correlation function).
 *
 */


class TreeCorrNGObject: public GGLensData {

 public:
  
  // constructors
  TreeCorrNGObject() {}
  TreeCorrNGObject(string ascii_filename);
  // ordinary return value functions
  float getTotalNPairs() const {
    float sum = 0.0;
    for (int i=0; i < npairs.size(); ++i) {
      sum += npairs[i];
    }
    return sum;
  }
  vector<float> getR_nom() { return R_nom; }

 private:

  vector<float> R_nom;  // nominal radial bin center
  vector<float> weight;
  vector<float> npairs;

};


#endif // TREECORROBJECTS_H

