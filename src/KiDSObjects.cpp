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
  //  e1/2_A/B/C, ALPHA_J2000/DELTA_J2000, Xpos/Ypos,
  //  MAN_MASK, MAG_GAAP_r_CALIB, MAGERR_GAAP_r,
  //  PZ_full, Z_B, (Z_B_MIN, Z_B_MAX)
  //  FWHM_IMAGE, weight, and assign it to SourceObject

  valarray<float> g1a;
  valarray<float> g2a;
  valarray<float> g1b;
  valarray<float> g2b;
  valarray<float> g1c;
  valarray<float> g2c;
  valarray<float> wa;
  valarray<float> wb;
  valarray<float> wc;

  int max_src_count;

  try {
    CCfits::Column& column1 = table.column("e1_A");
    max_src_count = column1.rows();
    column1.read( g1a, 1, max_src_count );
    CCfits::Column& column2 = table.column("e2_A");
    column2.read( g2a, 1, max_src_count );
    CCfits::Column& column3 = table.column("weight_A");
    column3.read( wa, 1, max_src_count );

    CCfits::Column& column4 = table.column("e1_B");
    column4.read( g1b, 1, max_src_count );
    CCfits::Column& column5 = table.column("e2_B");
    column5.read( g2b, 1, max_src_count );
    CCfits::Column& column6 = table.column("weight_B");
    column6.read( wb, 1, max_src_count );

    CCfits::Column& column7 = table.column("e1_C");
    column7.read( g1c, 1, max_src_count );
    CCfits::Column& column8 = table.column("e2_C");
    column8.read( g2c, 1, max_src_count );
    CCfits::Column& column9 = table.column("weight_c");
    column9.read( wc, 1, max_src_count );


  }
  catch (CCfits::Table::NoSuchColumn& m) {
    CCfits::Column& column1 = table.column("e1");
    max_src_count = column1.rows();
    column1.read( g1a, 1, max_src_count );
    g1b = g1c = g1a;
    CCfits::Column& column2 = table.column("e2");
    column2.read( g2a, 1, max_src_count );
    g2b = g2c = g2a;
  }

  valarray<double> ra;
  CCfits::Column& column10 = table.column("ALPHA_J2000");
  column10.read( ra, 1, max_src_count );
  valarray<double> dec;
  CCfits::Column& column11 = table.column("DELTA_J2000");
  column11.read( dec, 1, max_src_count );

  valarray<float> xpos;
  CCfits::Column& column12 = table.column("Xpos_THELI");
  column12.read( xpos, 1, max_src_count );
  valarray<float> ypos;
  CCfits::Column& column13 = table.column("Ypos_THELI");
  column13.read( ypos, 1, max_src_count );

  valarray<int> mask;
  CCfits::Column& column14 = table.column("MASK");
  column14.read( mask, 1, max_src_count );

  valarray<float> mag;
  CCfits::Column& column15 = table.column("MAG_GAAP_r_CALIB");
  column15.read( mag, 1, max_src_count );
  valarray<float> magerr;
  CCfits::Column& column16 = table.column("MAGERR_GAAP_r");
  column16.read( magerr, 1, max_src_count );

  valarray<float> fwhm;
  CCfits::Column& column17 = table.column("FWHM_IMAGE_THELI");
  column17.read( fwhm, 1, max_src_count );

  valarray<float> sn;
  CCfits::Column& column18 = table.column("SNratio");
  column18.read( sn, 1, max_src_count );

  valarray<double> z_B;
  CCfits::Column& column19 = table.column("Z_B");
  column19.read( z_B, 1, max_src_count );

  /*
  valarray<float> c1a;
  valarray<float> c2a;
  valarray<float> c1b;
  valarray<float> c2b;
  valarray<float> c1c;
  valarray<float> c2c;
  valarray<float> m_corr;
  try {
    CCfits::Column& column20 = table.column("c1_A");
    column20.read( c1a, 1, max_src_count );
    CCfits::Column& column21 = table.column("c2_A");
    column21.read( c2a, 1, max_src_count );
    CCfits::Column& column22 = table.column("c1_B");
    column22.read( c1b, 1, max_src_count );
    CCfits::Column& column23 = table.column("c2_B");
    column23.read( c2b, 1, max_src_count );
    CCfits::Column& column24 = table.column("c1_C");
    column24.read( c1c, 1, max_src_count );
    CCfits::Column& column25 = table.column("c2_C");
    column25.read( c2c, 1, max_src_count );
    CCfits::Column& column26 = table.column("m_cor");
    column26.read( m_corr, 1, max_src_count );
  }
  catch (CCfits::Table::NoSuchColumn& m) {
    CCfits::Column& column21 = table.column("c1");
    column21.read( c1a, 1, max_src_count );
    c1b = c1c = c1a;
    CCfits::Column& column22 = table.column("c2");
    column22.read( c2a, 1, max_src_count );
    c2b = c2c = c2a;
    CCfits::Column& column23 = table.column("m_cor");
    column23.read( m_corr, 1, max_src_count );
  }
  */

  // append objects to this list
  source_list.reserve(max_src_count);
  for (int i=0; i<max_src_count; ++i) {
      if (weight[i] == 0) continue;
      if (mask[i] & bitmask) continue;
      KiDSObject* ptr = new KiDSObject(i, ra[i], dec[i], mag[i], xpos[i], ypos[i], fwhm[i],
				       g1a[i], g2a[i], g1b[i], g2b[i],
				       g1c[i], g2c[i],
				       sn[i], z_B[i], pz_full[i], mask[i], weight[i],
				       // m_corr[i], c1a[i], c2a[i], c1b[i], c2b[i], c1c[i], c2c[i]
				       );
      source_list.push_back(ptr);
  }

  // set up the p(z) redshift bins
  // "a vector of length 70 giving P(z) at redshifts spanning 0<z<3.5 with dz=0.05"
  float array[NUM_PZ_ELEM];
  for (int i=0; i<NUM_PZ_ELEM; ++i) {
    array[i] = INIT_Z + i * DELTA_Z;
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
  return source_list.size();
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
  return source_list.size();
}


void
KiDSObjectList::setBlinding(char blinding){
  int blind_index = 0;
  switch (blinding) {
    case 'A':
      blind_index = 0;
      break;
    case 'B':
      blind_index = 1;
      break;
    case 'C':
      blind_index = 2;
      break;
    case 'D':
      blind_index = 3;
      break;
    default:
      throw KiDSObjectsError("setBlinding: wrong blinding specified");
  }
  for (vector<KiDSObject*>::iterator it = source_list.begin();
       it != source_list.end(); ++it) {
    Shear s = (*it)->getShearArray()[blind_index];
    double g1, g2;
    s.setG1G2(g1, g2);
    (*it)->setShearG1G2BiasCorrections(s, g1, g2, (*it)->getM(),
				       (*it)->getC1Array()[blind_index],
				       (*it)->getC2Array()[blind_index]);
  }

  return;
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
