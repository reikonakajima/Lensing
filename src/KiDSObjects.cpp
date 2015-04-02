//
// KiDSObjects.cpp
//
#include "KiDSObjects.h"
using namespace std;


// The following variables are static members of the KiDSObject class.
// They need to be defined outside the scope of the class in order to be accessible from outside.
bool KiDSObject::usePixelCoords;


KiDSObjectList::KiDSObjectList(const string fits_filename, int bitmask) {

  // open FITS file
  const string obj_extension = "OBJECTS";
  auto_ptr<CCfits::FITS> pInfile(0);

  try {
    // open the fits table and go to the right extension
    pInfile.reset(new CCfits::FITS(fits_filename,CCfits::Read, obj_extension));

  } catch (CCfits::FITS::CantOpen &fitsError) {
      throw KiDSObjectsError(string("Error opening the file with message: ")+fitsError.message());
  } catch (CCfits::FITS::NoSuchHDU &fitsError) {
      throw KiDSObjectsError(string("Error going to the HDU: ") + fitsError.message());
  } catch (CCfits::FitsException &fitsError) {
      throw KiDSObjectsError(string("Error in FITS: ") + fitsError.message());
  }

  // goto OBJECTS extension
  CCfits::ExtHDU& table = pInfile->extension(obj_extension);

  // read the following columns (annoyingly, only one column can be read at a time):
  //  e1/2_A/B/C/D, ALPHA_J2000/DELTA_J2000, Xpos/Ypos,
  //  MAN_MASK, MAG_GAAP_r_CALIB, MAGERR_GAAP_r,
  //  PZ_full, Z_B, (Z_B_MIN, Z_B_MAX)
  //  FWHM_IMAGE, weight, and assign it to SourceObject

  valarray<float> g1a;
  CCfits::Column& column1 = table.column("e1_A");
  int max_src_count = column1.rows();
  // DEBUG REMOVE
  max_src_count = MIN(50000, max_src_count);  // for debugging purposes, to be removed!

  column1.read( g1a, 1, max_src_count );
  valarray<float> g2a;
  CCfits::Column& column2 = table.column("e2_A");
  column2.read( g2a, 1, max_src_count );
  valarray<float> g1b;
  CCfits::Column& column3 = table.column("e1_B");
  column3.read( g1b, 1, max_src_count );
  valarray<float> g2b;
  CCfits::Column& column4 = table.column("e2_B");
  column4.read( g2b, 1, max_src_count );
  valarray<float> g1c;
  CCfits::Column& column5 = table.column("e1_C");
  column5.read( g1c, 1, max_src_count );
  valarray<float> g2c;
  CCfits::Column& column6 = table.column("e2_C");
  column6.read( g2c, 1, max_src_count );
  valarray<float> g1d;
  CCfits::Column& column7 = table.column("e1_D");
  column7.read( g1d, 1, max_src_count );
  valarray<float> g2d;
  CCfits::Column& column8 = table.column("e2_D");
  column8.read( g2d, 1, max_src_count );

  valarray<double> ra;
  CCfits::Column& column9 = table.column("ALPHA_J2000");
  column9.read( ra, 1, max_src_count );
  valarray<double> dec;
  CCfits::Column& column10 = table.column("DELTA_J2000");
  column10.read( dec, 1, max_src_count );

  valarray<float> xpos;
  CCfits::Column& column11 = table.column("Xpos");
  column11.read( xpos, 1, max_src_count );
  valarray<float> ypos;
  CCfits::Column& column12 = table.column("Ypos");
  column12.read( ypos, 1, max_src_count );

  valarray<int> mask;
  CCfits::Column& column13 = table.column("MAN_MASK");
  column13.read( mask, 1, max_src_count );

  valarray<float> mag;
  CCfits::Column& column14 = table.column("MAG_GAAP_r_CALIB");
  column14.read( mag, 1, max_src_count );
  valarray<float> magerr;
  CCfits::Column& column15 = table.column("MAGERR_GAAP_r");
  column15.read( magerr, 1, max_src_count );

  valarray<float> fwhm;
  CCfits::Column& column16 = table.column("FWHM_IMAGE");
  column16.read( fwhm, 1, max_src_count );

  valarray<float> weight;
  CCfits::Column& column17 = table.column("weight");
  column17.read( weight, 1, max_src_count );

  valarray<float> sn;
  CCfits::Column& column18 = table.column("SNratio");
  column18.read( sn, 1, max_src_count );

  valarray<double> z_B;
  CCfits::Column& column19 = table.column("Z_B");
  column19.read( z_B, 1, max_src_count );

  vector<valarray<float> > pz_full;
  CCfits::Column& column20 = table.column("PZ_full");
  column20.readArrays( pz_full, 1, max_src_count );

  // append objects to this list
  source_list.reserve(max_src_count);
  for (int i=0; i<max_src_count; ++i) {
      if (weight[i] == 0) continue;
      if (mask[i] & bitmask) continue;
      KiDSObject* ptr = new KiDSObject(i, ra[i], dec[i], mag[i], xpos[i], ypos[i], fwhm[i],
				       g1a[i], g2a[i], g1b[i], g2b[i],
				       g1c[i], g2c[i], g1d[i], g2d[i],
				       sn[i], z_B[i], pz_full[i], mask[i], weight[i]);
      source_list.push_back(ptr);
  }

  // set up the p(z) redshift bins
  // "a vector of length 70 giving P(z) at redshifts spanning 0<z<3.5 with dz=0.05"
  float array[NUM_PZ_ELEM];
  for (int i=0; i<NUM_PZ_ELEM; ++i) {
    array[i] = i * DELTA_Z;
  }
  SourceObjectList::pzbins = valarray<float>(array, NUM_PZ_ELEM);

  // The following two commands needs to be set after the KiDSObject list has been filled:
  // - initially set shear to "A" (driver code can change this)
  // setShearIndex(0);  // TODO -- be able to set to any of the ABCD shears!  (initially set to 0)
  // - initially set coordinates to use ra/dec
  //   (driver code must specify if pixel coordinate is to be used)
  usePixelCoord(false);

  return;
}


int
KiDSObjectList::applyMask(int mask_thres) {
  for (vector<KiDSObject*>::iterator it = source_list.begin();
       it != source_list.end(); /* no increment */) {
    if ((*it)->getMask() > mask_thres) {
      it = source_list.erase(it);
    } else {
      ++it;
    }
  }
}


int
KiDSObjectList::applyBitMask(int bitmask){
  for (vector<KiDSObject*>::iterator it = source_list.begin();
       it != source_list.end(); /* no increment */) {
    if ((*it)->getMask() & bitmask) {
      it = source_list.erase(it);
    } else {
      ++it;
    }
  }
}


void
KiDSObject::printLine(ostream& os) const {
  os << id << " "
     << setprecision(5) << setw(10)
     << ra << " " << setw(10) << dec << " "
     << setprecision(3) << setw(10)
     << this->getG1() << " " << this->getG2() << " "
     << wt << " " << mag << " "
     << xpos << " " << ypos << " "
     << fwhm << " " << sn << " ";
  for (int i=0; i<NUM_SHEAR; ++i) {
    os << shear[i] << " ";
  }
  os << endl;
  return;
}
