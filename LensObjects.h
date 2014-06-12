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


class LensObjectList {
 public:
  LensObjectList() {}  // empty list
  LensObjectList(istream& is);

  Bounds<double> getBounds() { if (!bounds) findBounds(); return bounds;}
  void findBounds();  // finds the bounds of the objects in this list and saves it

  list<LensObject*>::iterator begin() { return lens_list.begin(); }
  list<LensObject*>::iterator end() { return lens_list.end(); }
  int size() { return lens_list.size(); }
  
  void sortByRA();
  void sortByDec();
  LensObjectList cullByRA(double minra, double maxra);
  LensObjectList cullByDec(double mindec, double maxdec);

  vector<LensObject*> getVectorForm();

 private:

  list<LensObject*> lens_list;
  Bounds<double> bounds;  // run findBounds() to set value
  list<LensObject*>::iterator searchRA(list<LensObject*>::iterator first, 
				       list<LensObject*>::iterator last,
				       const double ra);
  list<LensObject*>::iterator searchDec(list<LensObject*>::iterator first, 
					list<LensObject*>::iterator last,
					const double dec);
};


#endif // LENSOBJECTS_H
