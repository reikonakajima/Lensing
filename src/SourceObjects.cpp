#include <iostream>
using namespace std;
#include "Std.h"
#include "StringStuff.h"
#include "SourceObjects.h"
#include "AstronomicalConstants.h"



bool Compare_Source_RA(SourceObject* rhs, SourceObject* lhs) {
  return rhs->getRA() < lhs->getRA(); // sort in increasing order
}


void 
SourceObjectList::sortByRA() {
  this->sort(Compare_Source_RA);  // check?
  return;
}

bool Compare_Source_Dec(SourceObject* lhs, SourceObject* rhs) {
  return lhs->getDec() < rhs->getDec(); // sort in increasing order
}

void 
SourceObjectList::sortByDec() {
  this->sort(Compare_Source_Dec);  // check?
  return;
}



void 
SourceObjectList::setBounds() {

  double minra, maxra, mindec, maxdec;

  this->sortByRA();
  SourceObjectList::const_iterator i = this->begin();
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



vector<SourceObject*>
SourceObjectList::getVectorForm() {

  vector<SourceObject*> vectorform;
  SourceObjectList::const_iterator i = this->begin();
  for (; i != this->end(); ++i) {
    vectorform.push_back(*i);
  }
  return vectorform;
}


SourceObject::SourceObject(ifstream& ifs) {

  wt = -1.;          // indicate that the weight has not been set (if < 0)
  responsiv = -10.;  // indicate that the responsivity has not been set (if < -1)

}

/*
void
SourceObject::printLine(ostream& os) const {

  os << fixed << setprecision(5)
     << ra << " " << dec << "  "
    /*
     << angle << " " 
     << setprecision(0)
     << setw(6) << setfill('0')
     << run << " "
     << reduction << " " << camcol << " " 
     << setw(4) << setfill('0') << field << " " 
     << setw(4) << setfill('0') << id << " " 
     << setw(4) << setfill('0') << x_ccd << " " 
     << setw(4) << setfill('0') << y_ccd << " "
    * /
     << setprecision(6)
     << e1 << " " << e2 << "  "
     << shapeerr << " " << eRMS << "  "
     << resr << " " << resi << "  "
     << rmag << "  " 
     << trueZ << " " << zlrg << " " 
     << endl;

  return;
}



void
SourceObject::printLineInBinary(ofstream& ofs) const {

  ofs.write((char*) &(*this), sizeof(SourceObject));

  return;
}


*/
