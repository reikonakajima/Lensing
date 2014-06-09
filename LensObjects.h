//
// LensObjects.h
//

#ifndef LENSOBJECTS_H
#define LENSOBJECTS_H
#include <iostream>
#include <list>
#include <vector>
//#include "Cosmology.h"
#include "Bounds.h"
//#include "AstronomicalConstants.h"
//#include "Maggies.h"
//using namespace cosmology;
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
  /*
  float getRedshift() const {return z;}
  float getRedshiftError() const {return zerr;}
  float getRedshiftConfidence() const {return zconf;}

  void setRedshift(float z_) {z = z_; return;}
  void setRedshiftError(float zerr_) {zerr = zerr_; return;}
  void setRedshiftConfidence(float zconf_) {zconf = zconf_; return;}

  float getGFracDeV() const {return fracdev_g;}
  float getRFracDeV() const {return fracdev_r;}
  float getIFracDeV() const {return fracdev_i;}
  float getFracDeV() const { return (fracdev_g + fracdev_r + fracdev_i)/3.; }

  float getUMag() const {return umag;}
  float getGMag() const {return gmag;}
  float getRMag() const {return rmag;}
  float getIMag() const {return imag;}
  float getZMag() const {return zmag;}

  float getUMagErr() const {return umagerr;}
  float getGMagErr() const {return gmagerr;}
  float getRMagErr() const {return rmagerr;}
  float getIMagErr() const {return imagerr;}
  float getZMagErr() const {return zmagerr;}

  string getUGRIZString() const;
  string getMaggieString() const;
  string getMaggieInvVarString() const;

  float getDistanceModulus(const Cosmology& c, float hubbleconst, float alternatez = 0.) const;
  float getSimpleRAbsMag(const Cosmology& c,  float hubbleconst, float alternatez = 0.) const 
  { return rmag - this->getDistanceModulus(c, hubbleconst, alternatez); }
  */
  void printLine(ostream& os) const;
  //void printLineWithModifiedRADec(ostream& os, double newra, double newdec) const;

  // for use with Mesh.h
  double getX() const {return ra;}
  double getY() const {return dec;}
  double getZ() const {return 0.;}

 private:
  string id;
  double ra, dec;
  double mag;
  /*
  float z, zerr, zconf;
  float fracdev_g, fracdev_r, fracdev_i;
  float umag, gmag, rmag, imag, zmag;
  float umagerr, gmagerr, rmagerr, imagerr, zmagerr;
  */
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
