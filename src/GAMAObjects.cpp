//
// GAMAObjects.cpp
//
#include "GAMAObjects.h"
#include "StringStuff.h"
using namespace std;


GAMAObject::GAMAObject(const string buffer, const int _id) {
  istringstream iss(buffer);
  if (!(iss >> type >> mag >> xstar >> ystar >> LensObject::ra >> LensObject::dec >> rhalo)) {
    cerr << "## " << buffer << endl;
    throw GAMAObjectsError("error reading GAMAObject");
  }
  LensObject::id = _id;
}


void 
GAMAObject::printLine(ostream& os) const{
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


GAMAObjectList::GAMAObjectList(istream& is) {
  string buffer;
  int id = 0;
  while (getlineNoComment(is, buffer)) {
    GAMAObject* ptr = new GAMAObject(buffer, id);
    lens_list.push_back(ptr);
    ++id;
  }
}
