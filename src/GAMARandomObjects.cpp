//
// GAMARandomObjects.cpp
//
#include "GAMARandomObjects.h"
#include <algorithm>  // std::random_shuffle
using namespace std;


GAMARandomObjectList::GAMARandomObjectList(const string fits_filename, int max_count) {

/*
 * Content of GAMARandomObject FITS file
 *
  >>> import astropy.io.fits as pf
  >>> hdulist = pf.open('GAMA_rand_Planck14_129p0_0p5.fits')
  >>> hdulist.info()
Filename: GAMA_rand_Planck14_129p0_0p5.fits
No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU       4   ()              
1    RANDOMS_BINARY  BinTableHDU     29   198296R x 7C   [1E, 1E, 1E, 1E, E, E, E]   
>>> hdulist[1].columns
ColDefs(
    name = 'RA'; format = '1E'; unit = 'degrees'
    name = 'DEC'; format = '1E'; unit = 'degrees'
    name = 'Z_1'; format = '1E'
    name = 'R_COMOVING'; format = '1E'; unit = 'Mpc/h (Planck14)'
    name = 'logmstar'; format = 'E'; unit = 'dex'
    name = 'absmag_r'; format = 'E'; unit = 'mag'
    name = 'uminusr'; format = 'E'; unit = 'mag'
)
 *
 */


  // open FITS file
  const string obj_extension = "RANDOMS_BINARY";
  auto_ptr<CCfits::FITS> pInfile(0);

  try {
    // open the fits table and go to the right extension
    pInfile.reset(new CCfits::FITS(fits_filename, CCfits::Read, obj_extension));

  } catch (CCfits::FITS::CantOpen &fitsError) {
      throw GAMARandomObjectsError(string("Error opening the file with message: ")+fitsError.message());
  } catch (CCfits::FITS::NoSuchHDU &fitsError) {
      throw GAMARandomObjectsError(string("Error going to the HDU: ") + fitsError.message());
  } catch (CCfits::FitsException &fitsError) {
      throw GAMARandomObjectsError(string("Error in FITS: ") + fitsError.message());
  }

  // goto OBJECTS extension
  CCfits::ExtHDU& table = pInfile->extension(obj_extension);

  // read the following columns (annoyingly, only one column can be read at a time):
  //  id, RA/DEC, Z_1, FIELD_POS, R_COMOVING

  valarray<double> ra;
  CCfits::Column& column1 = table.column("RA");
  column1.read( ra, 1, column1.rows() );
  valarray<double> dec;
  CCfits::Column& column2 = table.column("DEC");
  column2.read( dec, 1, column2.rows() );

  valarray<float> z;
  CCfits::Column& column3 = table.column("Z_1");
  column3.read( z, 1, column3.rows() );

  valarray<float> r_comoving;
  CCfits::Column& column4 = table.column("R_COMOVING");
  column4.read( r_comoving, 1, column4.rows() );

  valarray<float> logmstar;
  CCfits::Column& column5 = table.column("logmstar");
  column5.read( logmstar, 1, column5.rows() );

  valarray<float> absmag_r;
  CCfits::Column& column6 = table.column("absmag_r");
  column6.read( absmag_r, 1, column6.rows() );

  valarray<float> uminusr;
  CCfits::Column& column7 = table.column("uminusr");
  column7.read( uminusr, 1, column7.rows() );


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
    GAMARandomObject* ptr = new GAMARandomObject(ii, ra[ii], dec[ii], z[ii], r_comoving[ii],
						 logmstar[ii], absmag_r[ii], uminusr[ii]);
    lens_list.push_back(ptr);
  }

  return;
}


int
GAMARandomObjectList::applyLogMStarCut(float min_logmstar, float max_logmstar) {
  for (vector<GAMARandomObject*>::iterator it = lens_list.begin();
       it != lens_list.end(); /* no increment */) {
    float logmstar = (*it)->getLogMStar();
    if ( (*it)->getLogMStar() < min_logmstar or (*it)->getLogMStar() >= max_logmstar ) {
      it = lens_list.erase(it);
    } else {
      ++it;  // increment here; only if there were no deletion
    }
  }
  return lens_list.size();
}


/*
void
RandomObject::printLine(ostream& os) const {
  os << id << " "
     << setprecision(5) << setw(10)
     << ra << " " << setw(10) << dec << " "
     << setprecision(3) << setw(10)
     << this->getE1() << " " << this->getE2() << " "
     << wt << " " << mag << " "
     << xpos << " " << ypos << " "
     << fwhm << " " << sn << " ";
  for (int i=0; i<NUM_SHEAR; ++i) {
    os << shear[i] << " ";
  }
  os << endl;
  return;
}
*/
