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


LensObjectList::LensObjectList(istream& is) {
  string buffer;
  while (getlineNoComment(is, buffer)) {
    LensObject* ptr = new LensObject(buffer);
    lens_list.push_back(ptr);
  }
}


bool Compare_Source_RA(LensObject* rhs, LensObject* lhs) {
  return rhs->getRA() < lhs->getRA(); // sort in increasing order
}
bool Compare_Source_Dec(LensObject* lhs, LensObject* rhs) {
  return lhs->getDec() < rhs->getDec(); // sort in increasing order
}


void 
LensObjectList::sortByRA() {
  lens_list.sort(Compare_Source_RA);
  return;
}


void 
LensObjectList::sortByDec() {
  lens_list.sort(Compare_Source_Dec);
  return;
}


LensObjectList
LensObjectList::cullByRA(double minra, double maxra) {
  LensObjectList culledlist;
  this->sortByRA();
  list<LensObject*>::iterator i0 = searchRA(lens_list.begin(), lens_list.end(), minra);
  list<LensObject*>::iterator i1 = searchRA(i0, lens_list.end(), maxra);
  culledlist.lens_list.assign(i0, i1);
  return culledlist;
}


LensObjectList
LensObjectList::cullByDec(double mindec, double maxdec) {
  LensObjectList culledlist;
  this->sortByDec();
  list<LensObject*>::iterator i0 = searchDec(lens_list.begin(), lens_list.end(), mindec);
  list<LensObject*>::iterator i1 = searchDec(i0, lens_list.end(), maxdec);
  culledlist.lens_list.assign(i0, i1);
  return culledlist;
}


list<LensObject*>::iterator 
LensObjectList::searchRA(list<LensObject*>::iterator first, list<LensObject*>::iterator last, 
			 const double ra) {
  list<LensObject*>::iterator i = first;
  while ( ra > (*i)->getRA() && i != last)
    ++i;
  return i;
}

list<LensObject*>::iterator 
LensObjectList::searchDec(list<LensObject*>::iterator first, list<LensObject*>::iterator last, 
			  const double dec) {
  list<LensObject*>::iterator i = first;
  while ( dec > (*i)->getDec() && i != last)
    ++i;
  return i;
}


void 
LensObjectList::findBounds() {

  double minra, maxra, mindec, maxdec;

  this->sortByRA();
  list<LensObject*>::const_iterator i = lens_list.begin();
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


vector<LensObject*> 
LensObjectList::getVectorForm() {

  vector<LensObject*> vectorform;
  vectorform.reserve(lens_list.size());
  list<LensObject*>::const_iterator i = lens_list.begin();
  for (; i != lens_list.end(); ++i) {
    vectorform.push_back(*i);
  }
  return vectorform;
}


