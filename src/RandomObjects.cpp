//
// RandomObjects.cpp
//
#include "RandomObjects.h"
using namespace std;


RandomObjectList::RandomObjectList(const string fits_filename, int max_count) {

  // open FITS file
  const string obj_extension = "OBJECTS";
  auto_ptr<CCfits::FITS> pInfile(0);

  try {
    // open the fits table and go to the right extension
    pInfile.reset(new CCfits::FITS(fits_filename, CCfits::Read, obj_extension));

  } catch (CCfits::FITS::CantOpen &fitsError) {
      throw RandomObjectsError(string("Error opening the file with message: ")+fitsError.message());
  } catch (CCfits::FITS::NoSuchHDU &fitsError) {
      throw RandomObjectsError(string("Error going to the HDU: ") + fitsError.message());
  } catch (CCfits::FitsException &fitsError) {
      throw RandomObjectsError(string("Error in FITS: ") + fitsError.message());
  }

  // goto OBJECTS extension
  CCfits::ExtHDU& table = pInfile->extension(obj_extension);

  // read the following columns (annoyingly, only one column can be read at a time):
  //  ALPHA_J2000/DELTA_J2000, Xpos/Ypos, MASK

  valarray<double> ra;
  CCfits::Column& column1 = table.column("ALPHA_J2000");
  column1.read( ra, 1, column1.rows() );
  valarray<double> dec;
  CCfits::Column& column2 = table.column("DELTA_J2000");
  column2.read( dec, 1, column2.rows() );

  valarray<float> xpos;
  CCfits::Column& column3 = table.column("Xpos");
  column3.read( xpos, 1, column3.rows() );
  valarray<float> ypos;
  CCfits::Column& column4 = table.column("Ypos");
  column4.read( ypos, 1, column4.rows() );

  valarray<float> mask;
  CCfits::Column& column5 = table.column("MASK");
  column5.read( mask, 1, column5.rows() );

  // append objects to this list
  if (max_count < 0)
    max_count = column1.rows();
  else if (column1.rows() < max_count)
    max_count = column1.rows();

  lens_list.reserve(max_count);
  for (long int i=0; i<max_count; ++i) {
      RandomObject* ptr = new RandomObject(i, ra[i], dec[i], xpos[i], ypos[i], mask[i]);
      lens_list.push_back(ptr);
  }

  return;
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
