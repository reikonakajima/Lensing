//
// GGLens.h
//
#ifndef GGLENS_H
#define GGLENS_H
#include <iostream>
#include <list>
#include <vector>
#include "Bounds.h"
using std::list;
using std::istringstream;


class GGLensError : public MyException {
 public:
 GGLensError(const string& m="") :
  MyException("GGLensError: " +m) {}
};


class GGLensObject {
 public:
  GGLensObject(const string buffer);

 private:
  string id;
  double ra, dec;
  float mag;
};


class GGLensObjectList : public list<GGLensObject*> {
 public:
  GGLensObjectList() {}  // empty list
  /*
  GGLensObjectList(istream& is);
  void sortByRA();
  void sortByDec();
  LensObjectList cullByRA(double minra, double maxra);
  LensObjectList cullByDec(double mindec, double maxdec);
  Bounds<double> getBounds() { if (!bounds) findBounds(); return bounds;}
  void findBounds();  // finds the bounds of the objects in this list and saves it
  vector<LensObject*> getVectorForm();
  */
 private:
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
