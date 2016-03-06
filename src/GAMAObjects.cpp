//
// GAMAObjects.cpp
//
#include "GAMAObjects.h"
#include "StringStuff.h"
#include <algorithm>  // std::random_shuffle
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
    //<< d_absmag_r << " "
     << logmstar << " "
    //<< d_logmstar << " "
     << uminusr << " "
    //<< d_uminusr << " "
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


GAMAObjectList::GAMAObjectList(const string fits_filename, int max_count) {

  // open FITS file
  const string obj_extension = "JOINED";
  auto_ptr<CCfits::FITS> pInfile(0);

  try {
    // open the fits table and go to the right extension
    pInfile.reset(new CCfits::FITS(fits_filename, CCfits::Read, obj_extension));

  } catch (CCfits::FITS::CantOpen &fitsError) {
      throw GAMAObjectsError(string("Error opening the fits file with message: ") +
			     fitsError.message());
  } catch (CCfits::FITS::NoSuchHDU &fitsError) {
      throw GAMAObjectsError(string("Error going to the HDU: ") + fitsError.message());
  } catch (CCfits::FitsException &fitsError) {
      throw GAMAObjectsError(string("Error in FITS: ") + fitsError.message());
  }

  // goto OBJECTS extension
  CCfits::ExtHDU& table = pInfile->extension(obj_extension);

  // read the following columns (annoyingly, only one column can be read at a time):
  //  id, RA/DEC, Z_1, R_COMOVING, ...

  valarray<double> ra;
  CCfits::Column& column1 = table.column("RA");
  column1.read( ra, 1, column1.rows() );
  valarray<double> dec;
  CCfits::Column& column2 = table.column("DEC");
  column2.read( dec, 1, column2.rows() );

  valarray<float> z;
  CCfits::Column& column3 = table.column("Z_TONRY");
  column3.read( z, 1, column3.rows() );

  valarray<float> logmstar;
  CCfits::Column& column5 = table.column("logmstar");
  column5.read( logmstar, 1, column5.rows() );

  valarray<float> absmag_r;
  CCfits::Column& column6 = table.column("absmag_r");
  column6.read( absmag_r, 1, column6.rows() );

  valarray<float> uminusr;
  CCfits::Column& column7 = table.column("uminusr");
  column7.read( uminusr, 1, column7.rows() );

  valarray<float> r_comoving;
  CCfits::Column& column8 = table.column("R_COMOVING");
  column8.read( r_comoving, 1, column8.rows() );


  // create array of indices 0 through (num_randoms-1)
  int num_randoms = column1.rows();
  int index[num_randoms];
  for (int i=0; i<num_randoms; ++i) {
    index[i] = i;
  }
  // randomize the index array, then truncate at max_count
  if (max_count <= 0) {
    max_count = num_randoms;
  }
  random_shuffle(&index[0], &index[num_randoms]);
  lens_list.reserve(max_count);
  for (long int i=0; i<max_count; ++i) {
    int ii = index[i];  // randomized index
    GAMAObject* ptr = new GAMAObject(ii, ra[ii], dec[ii], z[ii], r_comoving[ii],
				     logmstar[ii], absmag_r[ii], uminusr[ii]);
    lens_list.push_back(ptr);
  }

  return;
}


int
GAMAObjectList::applyLogMStarCut(float min_logmstar, float max_logmstar) {
  for (vector<GAMAObject*>::iterator it = lens_list.begin();
       it != lens_list.end(); /* no increment */) {
    float logmstar = (*it)->getLogMStar();
    if ( (*it)->getLogMStar() < min_logmstar or (*it)->getLogMStar() >= max_logmstar ) {
      it = lens_list.erase(it);
    } else {
      ++it;  // increment here; only if there were no deletion
    }
  }
  LensObjectList::bounds = Bounds<double>();  // undefine ra/dec bounds if list is modified
  return lens_list.size();
}
