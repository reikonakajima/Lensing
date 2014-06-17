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

  // read the following columns
  //  e1/2_A/B/C/D, ALPHA_J2000/DELTA_J2000, Xpos/Ypos, CANDIDATEMASK, MAG_BEST, MAGERR_BEST,
  //  FWHM_IMAGE, weight, and assign it to SourceObject

  vector< valarray<float> > e1a;
  table.column("e1_A").readArrays( e1a, 0, 0 );
  vector< valarray<float> > e2a;
  table.column("e2_A").readArrays( e2a, 0, 0 );
  vector< valarray<float> > e1b;
  table.column("e1_B").readArrays( e1b, 0, 0 );
  vector< valarray<float> > e2b;
  table.column("e2_B").readArrays( e2b, 0, 0 );
  vector< valarray<float> > e1c;
  table.column("e1_C").readArrays( e1c, 0, 0 );
  vector< valarray<float> > e2c;
  table.column("e2_C").readArrays( e2c, 0, 0 );
  vector< valarray<float> > e1d;
  table.column("e1_D").readArrays( e1d, 0, 0 );
  vector< valarray<float> > e2d;
  table.column("e2_D").readArrays( e2d, 0, 0 );

  vector< valarray<double> > ra;
  table.column("ALPHA_J2000").readArrays( ra, 0, 0 );
  vector< valarray<double> > dec;
  table.column("DELTA_J2000").readArrays( dec, 0, 0 );

  vector< valarray<float> > xpos;
  table.column("Xpos").readArrays( xpos, 0, 0 );
  vector< valarray<float> > ypos;
  table.column("Ypos").readArrays( ypos, 0, 0 );

  vector< valarray<float> > mask;
  table.column("CANDIDATEMASK").readArrays( mask, 0, 0 );

  vector< valarray<float> > mag;
  table.column("MAG_BEST").readArrays( mag, 0, 0 );
  vector< valarray<float> > magerr;
  table.column("MAGERR_BEST").readArrays( magerr, 0, 0 );

  vector< valarray<float> > fwhm;
  table.column("FWHM_IMAGE").readArrays( fwhm, 0, 0 );

  vector< valarray<float> > weight;
  table.column("weight").readArrays( weight, 0, 0 );

  // append object to this list
  for (int i=0; i<5; ++i) {
      cout << weight[0][i] << " ";
  }
  cout << endl;

}

