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

 private:
  string id;
  double ra, dec;
};


class LensObjectList : public list<LensObject*> {
 public:
  LensObjectList() {}  // empty list
  LensObjectList(istream& is);
  void sortByRA();
  void sortByDec();
  LensObjectList cullByRA(double minra, double maxra);
  LensObjectList cullByDec(double mindec, double maxdec);
  Bounds<double> getBounds() { if (!bounds) findBounds(); return bounds;}
  void findBounds();  // finds the bounds of the objects in this list and saves it
  vector<LensObject*> getVectorForm();
 private:
  Bounds<double> bounds;  // run findBounds() to set value
  LensObjectList::iterator searchRA(LensObjectList::iterator first, 
				    LensObjectList::iterator last,
				    const double ra);
  LensObjectList::iterator searchDec(LensObjectList::iterator first, 
				     LensObjectList::iterator last,
				     const double dec);
};


#endif // LENSOBJECTS_H
