#include "LensObjects.h"
#include "StringStuff.h"
using namespace std;

LensObject::LensObject(const string buffer) {
  istringstream iss(buffer);
  if (!(iss >> id >> ra >> dec >> mag
	/*
	>> z >> zerr >> zconf
	>> fracdev_g >> fracdev_r >> fracdev_i
	>> umag >> gmag >> rmag >> imag >> zmag
	>> umagerr >> gmagerr >> rmagerr >> imagerr >> zmagerr
	*/
	)) {
    cerr << "## " << buffer << endl;
    throw LensObjectsError("can't read LensObject");
  }
}


void 
LensObject::printLine(ostream& os) const{
  os << fixed << setprecision(0)
     << id << " "
     << setprecision(5) << setw(10)
     << ra << " " << setw(10) << dec << "  "
     << setprecision(4) << setw(6)
     << mag
    /*
     << z << " " << zerr << " " << zconf << "  "
     << fracdev_g << " " << fracdev_r << " " << fracdev_i << "  "
     << umag << " " << gmag << " " << rmag << " " << imag << " " << zmag << "  "
     << umagerr << " " << gmagerr << " " << rmagerr << " " << imagerr << " " << zmagerr << "  "
    */
     << endl;
  return;
}

/*
void 
LensObject::printLineWithModifiedRADec(ostream& os, double newra, double newdec) const {
  os << fixed << setprecision(0)
     << id << " "
     << setprecision(5) << setw(10)
     << newra << " " << setw(10) << newdec << "  "
     << setprecision(4) << setw(6)
     << z << " " << zerr << " " << zconf << "  "
     << fracdev_g << " " << fracdev_r << " " << fracdev_i << "  "
     << umag << " " << gmag << " " << rmag << " " << imag << " " << zmag << "  "
     << umagerr << " " << gmagerr << " " << rmagerr << " " << imagerr << " " << zmagerr << "  "
     << endl;
  return;
}


string 
LensObject::getUGRIZString() const {
  ostringstream oss;
  oss << setprecision(5);
  oss << umag << " "
      << gmag << " "
      << rmag << " "
      << imag << " "
      << zmag;
  return oss.str();
}


string 
LensObject::getMaggieString() const {
  ostringstream oss;
  oss << setprecision(5);
  oss << mag2maggies(umag, sdss_b_u) << " "
      << mag2maggies(gmag, sdss_b_g) << " "
      << mag2maggies(rmag, sdss_b_r) << " "
      << mag2maggies(imag, sdss_b_i) << " "
      << mag2maggies(zmag, sdss_b_z);
  return oss.str();
}


string 
LensObject::getMaggieInvVarString() const {
  ostringstream oss;
  oss << setprecision(5);
  oss << calculateInvVarMaggies(umag, umagerr, sdss_b_u) << " "
      << calculateInvVarMaggies(gmag, gmagerr, sdss_b_g) << " "
      << calculateInvVarMaggies(rmag, rmagerr, sdss_b_r) << " "
      << calculateInvVarMaggies(imag, imagerr, sdss_b_i) << " "
      << calculateInvVarMaggies(zmag, zmagerr, sdss_b_z);
  return oss.str();
}


float 
LensObject::getDistanceModulus(const Cosmology& c, float hubbleconst, float alternatez) const {
  
  float z = this->getRedshift();
  if (z == 0. && alternatez == 0.)
    throw LensObjectsError("getDistanceModulus: redshift not specified");
  else if (z == 0. || alternatez > 0.)
    z = alternatez;
  if (z < 0.)
    return 1000.0;
  double Dlum = c.DL(z) * HubbleLengthMpc;   // Lum Dist in [Mpc/h]
  Dlum *= 1.0e6 * hubbleconst;               // Lum Dist in [pc]
  float  distmod = 5. * (log10(Dlum) - 1.);  // distance modulus
  return distmod;
}
*/

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


