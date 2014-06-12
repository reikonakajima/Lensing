#include <algorithm>
#include "LensObjects.h"
#include "StringStuff.h"
using namespace std;

LensObject::LensObject(const string buffer) {
  istringstream iss(buffer);
  if (!(iss >> id >> ra >> dec)) {
    cerr << "## " << buffer << endl;
    throw LensObjectsError("error reading LensObject");
  }
}


void 
LensObject::printLine(ostream& os) const{
  os << fixed << setprecision(0)
     << id << " "
     << setprecision(5) << setw(10)
     << ra << " " << setw(10) << dec
     << endl;
  return;
}

template <class ObjPtr>
LensObjectList<ObjPtr>::LensObjectList(istream& is) {
  string buffer;
  while (getlineNoComment(is, buffer)) {
    ObjPtr ptr = new LensObject(buffer);
    lens_list.push_back(ptr);
  }
}


template <class ObjPtr>
void 
LensObjectList<ObjPtr>::sortByRA() {
  std::sort(lens_list.begin(), lens_list.end(), LensObjectList<ObjPtr>::Compare_Source_RA);
  return;
}


template <class ObjPtr>
void 
LensObjectList<ObjPtr>::sortByDec() {
  std::sort(lens_list.begin(), lens_list.end(), LensObjectList<ObjPtr>::Compare_Source_Dec);
  return;
}


template <class ObjPtr>
LensObjectList<ObjPtr>
LensObjectList<ObjPtr>::cullByRA(double minra, double maxra) {
  LensObjectList<ObjPtr> culledlist;
  this->sortByRA();
  typename vector<ObjPtr>::iterator i0 = searchRA(lens_list.begin(), lens_list.end(), minra);
  typename vector<ObjPtr>::iterator i1 = searchRA(i0, lens_list.end(), maxra);
  culledlist.lens_list.assign(i0, i1);
  return culledlist;
}


template <class ObjPtr>
LensObjectList<ObjPtr>
LensObjectList<ObjPtr>::cullByDec(double mindec, double maxdec) {
  LensObjectList<ObjPtr> culledlist;
  this->sortByDec();
  typename vector<ObjPtr>::iterator i0 = searchDec(lens_list.begin(), lens_list.end(), mindec);
  typename vector<ObjPtr>::iterator i1 = searchDec(i0, lens_list.end(), maxdec);
  culledlist.lens_list.assign(i0, i1);
  return culledlist;
}


template <class ObjPtr>
typename vector<ObjPtr>::iterator
LensObjectList<ObjPtr>::searchRA(typename vector<ObjPtr>::iterator first,
				 typename vector<ObjPtr>::iterator last,
				 const double ra) {
  typename vector<ObjPtr>::iterator i = first;
  while ( ra > (*i)->getRA() && i != last)
    ++i;
  return i;
}


template <class ObjPtr>
typename vector<ObjPtr>::iterator
LensObjectList<ObjPtr>::searchDec(typename vector<ObjPtr>::iterator first,
				  typename vector<ObjPtr>::iterator last,
				  const double dec) {
  typename vector<ObjPtr>::iterator i = first;
  while ( dec > (*i)->getDec() && i != last)
    ++i;
  return i;
}


template <class ObjPtr>
void 
LensObjectList<ObjPtr>::findBounds() {

  double minra, maxra, mindec, maxdec;

  this->sortByRA();
  typename vector<ObjPtr>::const_iterator i = lens_list.begin();
  minra = (*i)->getRA();
  i = lens_list.end();  --i;
  maxra = (*i)->getRA();

  this->sortByDec();
  i = lens_list.begin();
  mindec = (*i)->getDec();
  i = lens_list.end();  --i;
  maxdec = (*i)->getDec();

  bounds.setXMin(minra);
  bounds.setXMax(maxra);
  bounds.setYMin(mindec);
  bounds.setYMax(maxdec);

  return;
}



//
// explicit instantiations
//

template class LensObjectList<LensObject*>;
