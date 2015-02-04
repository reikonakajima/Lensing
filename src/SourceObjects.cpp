//
// SourceObjects.cpp
//
#include "SourceObjects.h"
using namespace std;


template <class ObjPtr>
void 
SourceObjectList<ObjPtr>::sortByRA() {
  std::sort(source_list.begin(), source_list.end(), SourceObjectList<ObjPtr>::Compare_Source_RA);
  return;
}


template <class ObjPtr>
void 
SourceObjectList<ObjPtr>::sortByDec() {
  std::sort(source_list.begin(), source_list.end(), SourceObjectList<ObjPtr>::Compare_Source_Dec);
  return;
}


template <class ObjPtr>
void 
SourceObjectList<ObjPtr>::findBounds() {

  double minra, maxra, mindec, maxdec;

  this->sortByRA();
  typename vector<ObjPtr>::const_iterator i = this->begin();
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


SourceObject::SourceObject() {

  wt = -1.;          // indicate that the weight has not been set (if < 0)
  responsiv = -10.;  // indicate that the responsivity has not been set (if < -1)
  //vare = varSN = -1.;  // indicate that these quantities have not been set

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

//
// explicit instantiations
//

template class SourceObjectList<SourceObject*>;
#include "RCSLenSObjects.h"
template class SourceObjectList<RCSLenSObject*>;
#include "KiDSObjects.h"
template class SourceObjectList<KiDSObject*>;
