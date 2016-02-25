//
// Mesh class 
//   Generates a mesh of pointers for easy access in 2d/3d.  
//   The objects in mesh should have methods getX() getY() and getZ()
//
// 2009.10.01  original code from Martin White 
// 2009.10.05  RN format modified
//             RN MeshClass -> Mesh 
//                instead of float (for position type), use Tpos template
// 2009.10.12  RN give options to use (x,y,z) instead of pos[3]
// 2009.10.26  RN getNearAngleMap(): add spherical angle distances for 2d case of ra/dec
// 2009.11.03  RN explicit instantiation: RADecObjects.h
// 
#ifndef MESH_H
#define MESH_H
#include <list>
#include <vector>
#include <map>
using std::list;
using std::vector;
using std::multimap;
#include "Bounds.h"
#include "AstronomicalConstants.h"


//
// Mesh: specialized here for 3D vectors.
// For 2D objects just set nn3=1.
// Ttype should have a member function getPos(), which returns Tpos[]
//
template <class Ttype, class Tpos = double>   // Ttype = SDSSpzObject*, SourceObject*, etc.
class Mesh {
 public:
  // constructor for [min, max)
  Mesh(int nn1, int nn2, int nn3, vector<Ttype> &P, 
       bool period=false,
       Tpos minx = 0., Tpos maxx = 1., 
       Tpos miny = 0., Tpos maxy = 1., 
       Tpos minz = 0., Tpos maxz = 1.);
  // Returns a list of indices for objects near pos.
  // list<int> nearobj(Tpos pos[], Tpos rmin, Tpos rmax);
  list<int> getNearMeshList(Tpos x, Tpos y, Tpos z, Tpos rmax, Tpos rmin = 0);
  // Returns a map of (dist.sq., obj.index)
  multimap<double,int> getNearMeshMap(Tpos x, Tpos y, Tpos z, Tpos rmax, Tpos rmin = 0);
  // Returns a map of (spherical angle separation, obj.index)
  // *assumes ra,dec,theta are in degrees*
  multimap<double,int> getNearAngleMap(Tpos ra, Tpos dec, Tpos z, Tpos thetamax, Tpos thetamin=0);
  vector<multimap<double,int> > getNearAngleMap(Tpos ra, Tpos dec, Tpos z,
						vector<Tpos> thetabin_edges);
  // Returns a list of distances^2 for objects near pos.
  // list<Tpos> dist2(Tpos pos[], Tpos rmin, Tpos rmax);
  list<Tpos> getDistSqList(Tpos x, Tpos y, Tpos z, Tpos rmax, Tpos rmin = 0);
 private:
  vector<int>	head,next,nm;
  vector<Ttype>	*dat;
  bool isPeriodic;
  Tpos periodicX(Tpos x);      // Wraps x periodically in the range [xmin,xmax).
  Tpos periodicY(Tpos y);      // Wraps y periodically in the range [ymin,ymax).
  Tpos periodicZ(Tpos z);      // Wraps z periodically in the range [zmin,zmax).
  Tpos diffperiodicX(Tpos x);  // Periodic differences.
  Tpos diffperiodicY(Tpos y);  // Periodic differences.
  Tpos diffperiodicZ(Tpos z);  // Periodic differences.
  // Returns the distance^2 between points a and b.
  //Tpos  distance2(Tpos a[], Tpos b[]); 
  Tpos	distance2(Tpos ax, Tpos ay, Tpos az, Tpos bx, Tpos by, Tpos bz); 
  // Returns a list (vector) of the mesh points within sr of i.
  vector<int> closemeshes(int ix, int iy, int iz,
			  int srx, int sry, int srz);
  // prepares srx, sry, srz given ra/dec/thetamax (in degrees)
  void prepare_angular_indices(Tpos ra, Tpos dec, Tpos thetamax,  /* input */
			       int& ix, int& iy, int& iz,    /* output */
			       int& srx, int& sry, int& srz  /* output */ );
  Tpos calculate_angular_separation(Tpos ra0, Tpos sin_dec0, Tpos cos_dec0, Tpos ra1, Tpos dec1);
  Tpos xmin, xmax, dx;
  Tpos ymin, ymax, dy;
  Tpos zmin, zmax, dz;
  Tpos xmid, ymid, zmid;  // for use in diffperiodic
};


/*  // remove when using g++ (this works with icc), and add to the file used instead
#include "SourceObjects.h"
template class Mesh<SourceObject*, double>;  // explicit instantiation
#include "lensObjects.h"
template class Mesh<LensObject*, double>;    // explicit instantiation
#include "RedshiftObjects.h"
template class Mesh<SDSSpzObject*, double>;  // explicit instantiation
#include "RADecObjects.h"
template class Mesh<RADecEntry*, double>;    // explicit instantiation
#include "ggLensObjects.h"
template class Mesh<SDSSpzSourceObject*, double>; // explicit instantiation
template class Mesh<SDSSpzLensObject*, double>;   // explicit instantiation
#include "bgCatObjects.h"
template class Mesh<bgCatObject*, double>;     // explicit instantiation
#include "photCatObjects.h"
template class Mesh<photCatObject*, double>;     // explicit instantiation
#include "SDSSLensObject.h"
template class Mesh<SDSSLensObject*, double>;  // explicit instantiation
*/

/// remove this line if you are compiling Mesh, and plan to link with Mesh.o
//#include "Mesh.incl"

#endif // MESH_H
