//
// GAMAObjects.cpp
//
#include "GAMAObjects.h"
#include "StringStuff.h"
using namespace std;


GAMAObject::GAMAObject(const string buffer) {
  istringstream iss(buffer);
  if (!(iss >> LensObject::id >> group_id >> LensObject::ra >> LensObject::dec >> LensObject::z
	>> absmag_r >> d_absmag_r >> logmstar >> d_logmstar >> uminusr
	>> d_uminusr >> rpetro >> logmoverl >> d_logmoverl >> loglwage
	>> d_loglwage >> metal >> d_metal >> logtau >> d_logtau
	>> logmremnants >> d_logmremnants >> rankbcg >> Nfof >> Zmax_19P8
	>> Zmax_19P4)) {
    cerr << "## " << buffer << endl;
    throw GAMAObjectsError("error reading GAMAObject");
  }
  LensObject::mag = rpetro;  // OK?
}


void 
GAMAObject::printLine(ostream& os) const{
  os << fixed << setprecision(0)
     << id << " "
     << setprecision(5) << setw(10)
     << ra << " " << setw(10) << dec << " "
     << setprecision(1) << setw(10)
     << z << " "
     << absmag_r << " "
     << d_absmag_r << " "
     << logmstar << " "
     << d_logmstar << " "
     << uminusr << " "
     << d_uminusr << " "
     << endl;
  return;
}


GAMAObjectList::GAMAObjectList(istream& is) {
  string buffer;
  while (getlineNoComment(is, buffer)) {
    GAMAObject* ptr = new GAMAObject(buffer);
    lens_list.push_back(ptr);
  }
}


int
GAMAObjectList::applyLogMStarCut(float min_logmstar, float max_logmstar) {
  for (vector<GAMAObject*>::iterator it = lens_list.begin();
       it != lens_list.end(); /* no increment */) {
    if ( (*it)->getLogMStar() < min_logmstar or (*it)->getLogMStar() >= max_logmstar ) {
      it = lens_list.erase(it);
    } else {
      ++it;
    }
  }
  return lens_list.size();
}
