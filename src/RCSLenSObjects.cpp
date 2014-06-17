//
// RCSLenSObjects.cpp
//
#include "RCSLenSObjects.h"
using namespace std;


RCSLenSObjectList::RCSLenSObjectList(const string fits_filename) {

  // open FITS file
  const string obj_extension = "OBJECTS";
  auto_ptr<CCfits::FITS> pInfile(0);

  try {
    // open the fits table and go to the right extension
    pInfile.reset(new CCfits::FITS(fits_filename,CCfits::Read, obj_extension));

  } catch (CCfits::FITS::CantOpen &fitsError) {
      throw RCSLenSObjectsError(string("Error opening the file with message: ")+fitsError.message());
  } catch (CCfits::FITS::NoSuchHDU &fitsError) {
      throw RCSLenSObjectsError(string("Error going to the HDU: ") + fitsError.message());
  } catch (CCfits::FitsException &fitsError) {
      throw RCSLenSObjectsError(string("Error in FITS: ") + fitsError.message());
  }

  // goto OBJECTS extension
  CCfits::ExtHDU& table = pInfile->extension(obj_extension);

  // read the following columns (annoyingly, only one column can be read at a time):
  //  e1/2_A/B/C/D, ALPHA_J2000/DELTA_J2000, Xpos/Ypos, CANDIDATEMASK, MAG_BEST, MAGERR_BEST,
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
  CCfits::Column& column13 = table.column("MASK");
  column13.read( mask, 1, column13.rows() );

  valarray<float> mag;
  CCfits::Column& column14 = table.column("MAG_BEST");
  column14.read( mag, 1, column14.rows() );
  valarray<float> magerr;
  CCfits::Column& column15 = table.column("MAGERR_BEST");
  column15.read( magerr, 1, column15.rows() );

  valarray<float> fwhm;
  CCfits::Column& column16 = table.column("FWHM_IMAGE");
  column16.read( fwhm, 1, column16.rows() );

  valarray<float> weight;
  CCfits::Column& column17 = table.column("weight");
  column17.read( weight, 1, column17.rows() );

  valarray<float> sn;
  CCfits::Column& column18 = table.column("SNratio");
  column17.read( sn, 1, column18.rows() );

  // DEBUG
  for (int i=0; i<5; ++i) {
      cerr << e1a[i] << " ";
  }
  cerr << endl;

  // append objects to this list
  source_list.reserve(column1.rows());
  for (int i=0; i<column1.rows(); ++i) {
      if (weight[i] == 0) continue;
      RCSLenSObject* ptr = new RCSLenSObject(ra[i], dec[i], mag[i], xpos[i], ypos[i], fwhm[i],
					     e1a[i], e2a[i], e1b[i], e2b[i],
					     e1c[i], e2c[i], e1d[i], e2d[i],
					     sn[i], weight[i]);
      source_list.push_back(ptr);
  }

  return;
}

