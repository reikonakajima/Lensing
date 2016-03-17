//
// TreeCorrObjects.h
//
#ifndef TREECORROBJECTS_H
#define TREECORROBJECTS_H

#include <iostream>
#include <vector>
#include <CCfits/CCfits>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include "Bounds.h"
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


class TreeCorrNGObject {

 public:
  
  // constructors
  TreeCorrNGObject() {}
  TreeCorrNGObject(string ascii_filename);
  // ordinary return value functions
  int  getTotalNPairs() const { return 0; }
  vector<float> getMeanR() { return meanR; }

 private:

  vector<float> meanR;
  vector<float> gamT;
  vector<float> gamX;
  vector<float> sigma;
  vector<float> weight;
  vector<int>   npairs;

};


#endif // TREECORROBJECTS_H

