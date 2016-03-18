//
// GGLensData.h
//
#ifndef GGLENSDATA_H
#define GGLENSDATA_H

#include <iostream>
#include <vector>
#include "Std.h"
#include "StringStuff.h"
#include "Shear.h"
#include "Bounds.h"
using std::vector;


class GGLensDataError : public MyException {
 public:
 GGLensDataError(const string& m="") : 
  MyException("GGLensDataError: "+m) {}
};


//
// GGLensData: base class for any GGLens data output
//
class GGLensData {

 public:

  // constructor
  GGLensData() {}
  // unit conversions for meanR
  void rescaleMeanR(float scale);
  // interpolation functions
  void getValuesAt(vector<float> radial_vals,
		   vector<float>& signalT, vector<float>& signalX, vector<float>& var);
  // ordinary return value functions
  int getRBinSize() const { return meanR.size(); }
  vector<float> getMeanR() { return meanR; }

 protected:

  // minimum data to describe a GGLensData
  vector<float> meanR;
  vector<float> gamT;
  vector<float> gamX;
  vector<float> sigma;
};


#endif // GGLENSDATA_H
