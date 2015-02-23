//
// StarMaskObjects.cpp
//
#include "StarMaskObjects.h"
#include "StringStuff.h"
using namespace std;


StarMaskObject::StarMaskObject(const string buffer, const int _id) {
  istringstream iss(buffer);
  if (!(iss >> type >> mag >> xstar >> ystar >> LensObject::ra >> LensObject::dec >> rhalo)) {
    cerr << "## " << buffer << endl;
    throw StarMaskObjectsError("error reading StarMaskObject");
  }
  LensObject::id = type;
}


void 
StarMaskObject::printLine(ostream& os) const{
  os << fixed << setprecision(0)
     << id << " "
     << setprecision(5) << setw(10)
     << ra << " " << setw(10) << dec << " "
     << setprecision(1) << setw(10)
     << mag << " "
     << rhalo << " "
     << setprecision(5) << setw(10)
     << xstar << " " << ystar << " "
     << endl;
  return;
}


StarMaskObjectList::StarMaskObjectList(istream& is) {
  string buffer;
  int id = 0;
  while (getlineNoComment(is, buffer)) {
    StarMaskObject* ptr = new StarMaskObject(buffer, id);
    lens_list.push_back(ptr);
    ++id;
  }
}
