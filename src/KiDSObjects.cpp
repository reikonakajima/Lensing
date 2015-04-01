//
// KiDSObjects.cpp
//
#include "KiDSObjects.h"
using namespace std;


// The following variables are static members of the KiDSObject class.
// They need to be defined outside the scope of the class in order to be accessible from outside.
bool KiDSObject::usePixelCoords;


KiDSObjectList::KiDSObjectList(const string fits_filename) {

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

  valarray<float> e1a;
  CCfits::Column& column1 = table.column("e1_A");
  column1.read( e1a, 1, column1.rows() );
  valarray<float> e2a;
  CCfits::Column& column2 = table.column("e2_A");
  column2.read( e2a, 1, column2.rows() );
  valarray<float> e1b;
  CCfits::Column& column3 = table.column("e1_B");
  column3.read( e1b, 1, column3.rows() );
  valarray<float> e2b;
  CCfits::Column& column4 = table.column("e2_B");
  column4.read( e2b, 1, column4.rows() );
  valarray<float> e1c;
  CCfits::Column& column5 = table.column("e1_C");
  column5.read( e1c, 1, column5.rows() );
  valarray<float> e2c;
  CCfits::Column& column6 = table.column("e2_C");
  column6.read( e2c, 1, column6.rows() );
  valarray<float> e1d;
  CCfits::Column& column7 = table.column("e1_D");
  column7.read( e1d, 1, column7.rows() );
  valarray<float> e2d;
  CCfits::Column& column8 = table.column("e2_D");
  column8.read( e2d, 1, column8.rows() );

  valarray<double> ra;
  CCfits::Column& column9 = table.column("ALPHA_J2000");
  column9.read( ra, 1, column9.rows() );
  valarray<double> dec;
  CCfits::Column& column10 = table.column("DELTA_J2000");
  column10.read( dec, 1, column10.rows() );

  valarray<float> xpos;
  CCfits::Column& column11 = table.column("Xpos");
  column11.read( xpos, 1, column11.rows() );
  valarray<float> ypos;
  CCfits::Column& column12 = table.column("Ypos");
  column12.read( ypos, 1, column12.rows() );

  valarray<float> mask;
  CCfits::Column& column13 = table.column("MAN_MASK");
  column13.read( mask, 1, column13.rows() );

  valarray<float> mag;
  CCfits::Column& column14 = table.column("MAG_GAAP_r_CALIB");
  column14.read( mag, 1, column14.rows() );
  valarray<float> magerr;
  CCfits::Column& column15 = table.column("MAGERR_GAAP_r");
  column15.read( magerr, 1, column15.rows() );

  valarray<float> fwhm;
  CCfits::Column& column16 = table.column("FWHM_IMAGE");
  column16.read( fwhm, 1, column16.rows() );

  valarray<float> weight;
  CCfits::Column& column17 = table.column("weight");
  column17.read( weight, 1, column17.rows() );

  valarray<float> sn;
  CCfits::Column& column18 = table.column("SNratio");
  column18.read( sn, 1, column18.rows() );

  valarray<double> z_B;
  CCfits::Column& column19 = table.column("Z_B");
  column19.read( z_B, 1, column19.rows() );

  vector<valarray<float> > pz_full;
  CCfits::Column& column20 = table.column("PZ_full");
  column20.readArrays( pz_full, 1, column20.rows() );

  // append objects to this list
  source_list.reserve(column1.rows());
  for (int i=0; i<column1.rows(); ++i) {
      if (weight[i] == 0) continue;
      KiDSObject* ptr = new KiDSObject(i, ra[i], dec[i], mag[i], xpos[i], ypos[i], fwhm[i],
				       e1a[i], e2a[i], e1b[i], e2b[i],
				       e1c[i], e2c[i], e1d[i], e2d[i],
				       sn[i], z_B[i], pz_full[i], weight[i]);
      source_list.push_back(ptr);
  }

  // The following two commands needs to be set after the KiDSObject list has been filled:
  // - initially set shear to "A" (driver code can change this)
  // setShearIndex(0);  // TODO -- be able to set to any of the ABCD shears!  (initially set to 0)
  // - initially set coordinates to use ra/dec
  //   (driver code must specify if pixel coordinate is to be used)
  usePixelCoord(false);

  return;
}


void
KiDSObject::printLine(ostream& os) const {
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
