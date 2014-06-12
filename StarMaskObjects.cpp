//
// StarMaskObjects.cpp
//
#include "StarMaskObjects.h"
#include "StringStuff.h"
using namespace std;


StarMaskObject::StarMaskObject(const string buffer, const int _id) {
  istringstream iss(buffer);
  if (!(iss >> mag >> LensObject::ra >> LensObject::dec >> rhalo >> xstar >> ystar)) {
    cerr << "## " << buffer << endl;
    throw StarMaskObjectsError("error reading StarMaskObject");
  }
  LensObject::id = _id;
}


void 
StarMaskObject::printLine(ostream& os) const{
  os << fixed << setprecision(0)
     << id << " "
     << setprecision(5) << setw(10)
     << ra << " " << setw(10) << dec << " "
     << setprecision(1) << setw(10)
     << mag
     << endl;
  return;
}
