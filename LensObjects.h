//
// LensObjects.h
//

#ifndef LENSOBJECTS_H
#define LENSOBJECTS_H
#include <iostream>
#include <list>
#include <vector>
#include "Bounds.h"
using std::list;
using std::istringstream;


class LensObjectsError : public MyException {
 public:
 LensObjectsError(const string& m="") :
  MyException("LensObjectsError: " +m) {}
};


/*
 * class LensObject
 *
 * This is the base class for lens objects, and contains minimum information needed to calculate
 * lensing signal around the lens.
 *
 */

class LensObject {
 public:
  LensObject(const string buffer);

  string getId() const {return id;}
  double getRA() const {return ra;}
  double getDec() const {return dec;}
  Position<double> getRADec() const {return Position<double>(ra,dec);}

  void printLine(ostream& os) const;

  // for use with Mesh.h
  double getX() const {return ra;}
  double getY() const {return dec;}
  double getZ() const {return 0.;}

 protected:
  LensObject() {}
  string id;
  double ra, dec;
};


/*
 * class LensObjectList<>
 *
 * This is the base class for a list of lens objects, and contains minimum information
 * needed to calculate lensing signal around the lens.
 *
 */

template <class ObjPtr>
class LensObjectList {
 public:
  LensObjectList() {}  // empty list
  LensObjectList(istream& is);

  Bounds<double> getBounds() { if (!bounds) findBounds(); return bounds;}
  void findBounds();  // finds the bounds of the objects in this list and saves it

  typename vector<ObjPtr>::iterator begin() { return lens_list.begin(); }
  typename vector<ObjPtr>::iterator end() { return lens_list.end(); }
  int size() { return lens_list.size(); }
  
  void sortByRA();
  void sortByDec();
  LensObjectList cullByRA(double minra, double maxra);
  LensObjectList cullByDec(double mindec, double maxdec);

 protected:

  vector<ObjPtr> lens_list;
  Bounds<double> bounds;  // run findBounds() to set value
  typename vector<ObjPtr>::iterator searchRA(typename vector<ObjPtr>::iterator first,
					     typename vector<ObjPtr>::iterator last,
					     const double ra);
  typename vector<ObjPtr>::iterator searchDec(typename vector<ObjPtr>::iterator first,
					      typename vector<ObjPtr>::iterator last,
					      const double dec);
  static bool Compare_Source_RA(ObjPtr lhs, ObjPtr rhs) {
      return lhs->getRA() < rhs->getRA(); // sort in increasing order
  }
  static bool Compare_Source_Dec(ObjPtr lhs, ObjPtr rhs) {
      return lhs->getDec() < rhs->getDec(); // sort in increasing order
  }


};


#endif // LENSOBJECTS_H
