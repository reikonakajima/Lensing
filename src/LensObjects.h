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
  LensObject(long int _id, double _ra, double _dec, double _z, float _mag) {
    id = _id;
    ra = _ra;
    dec = _dec;
    z = _z;
    mag = _mag;
  }

  long int getId() const {return id;}
  double getRA() const {return ra;}
  double getDec() const {return dec;}
  double getRedshift() const {return z;}
  float  getMag() const {return mag;}
  Position<double> getRADec() const {return Position<double>(ra,dec);}

  void printLine(ostream& os) const;

  // for use with Mesh.h
  double getX() const {return ra;}
  double getY() const {return dec;}
  double getZ() const {return 0.;}

 protected:
  LensObject() {}
  long int id;
  double ra, dec;
  double z;
  float  mag;
};


/*
 * class LensObjectList<>
 *
 * This is the base class for a list of lens objects, and contains minimum information
 * needed to calculate the tangential shear signal around the lenses.
 *
 */

template <class ObjPtr>
class LensObjectList {
 public:
  /* Currently, the destructor is not properly implemented (the following will give run-time errors)
   * The proper implementation will keep track of pointer counts, and will delete only when
   * the counter goes to zero.
  ~LensObjectList() {
    for (int i=0; i<size(); ++i) delete lens_list[i];
  }
   */

  Bounds<double> getBounds() {
    if (size() <= 0) throw LensObjectsError("bounds undefined for empty lens list");
    if (!bounds) findBounds(); return bounds;
  }
  void findBounds();  // finds the bounds of the objects in this list and saves it

  typename vector<ObjPtr>::iterator begin() { return lens_list.begin(); }
  typename vector<ObjPtr>::iterator end() { return lens_list.end(); }
  int size() { return lens_list.size(); }
  
  void sortByRA();
  void sortByDec();
  LensObjectList cullByRA(double minra, double maxra);
  LensObjectList cullByDec(double mindec, double maxdec);

 protected:
  LensObjectList() {}  // hide default constructor: must be defined in derived classes

  vector<ObjPtr> lens_list;
  mutable Bounds<double> bounds;
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
