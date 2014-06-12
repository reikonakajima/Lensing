#include "LensObjects.h"
#include "StringStuff.h"
using namespace std;

LensObject::LensObject(const string buffer) {
  istringstream iss(buffer);
  if (!(iss >> id >> ra >> dec)) {
    cerr << "## " << buffer << endl;
    throw LensObjectsError("can't read LensObject");
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
    this->push_back(ptr);
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
  this->sort(Compare_Source_RA);
  return;
}


void 
LensObjectList::sortByDec() {
  this->sort(Compare_Source_Dec);
  return;
}


LensObjectList
LensObjectList::cullByRA(double minra, double maxra) {
  LensObjectList culledlist;
  this->sortByRA();
  LensObjectList::iterator i0 = searchRA(this->begin(), this->end(), minra);
  LensObjectList::iterator i1 = searchRA(i0, this->end(), maxra);
  culledlist.assign(i0, i1);
  return culledlist;
}


LensObjectList
LensObjectList::cullByDec(double mindec, double maxdec) {
  LensObjectList culledlist;
  this->sortByDec();
  LensObjectList::iterator i0 = searchDec(this->begin(), this->end(), mindec);
  LensObjectList::iterator i1 = searchDec(i0, this->end(), maxdec);
  culledlist.assign(i0, i1);
  return culledlist;
}


LensObjectList::iterator 
LensObjectList::searchRA(LensObjectList::iterator first, LensObjectList::iterator last, 
			 const double ra) {
  LensObjectList::iterator i = first;
  while ( ra > (*i)->getRA() && i != last)
    ++i;
  return i;
}

LensObjectList::iterator 
LensObjectList::searchDec(LensObjectList::iterator first, LensObjectList::iterator last, 
			  const double dec) {
  LensObjectList::iterator i = first;
  while ( dec > (*i)->getDec() && i != last)
    ++i;
  return i;
}


void 
LensObjectList::findBounds() {

  double minra, maxra, mindec, maxdec;

  this->sortByRA();
  LensObjectList::const_iterator i = this->begin();
  minra = (*i)->getRA();
  i = this->end();  --i;
  maxra = (*i)->getRA();

  this->sortByDec();
  i = this->begin();
  mindec = (*i)->getDec();
  i = this->end();  --i;
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
  vectorform.reserve(this->size());
  LensObjectList::const_iterator i = this->begin();
  for (; i != this->end(); ++i) {
    vectorform.push_back(*i);
  }
  return vectorform;
}


