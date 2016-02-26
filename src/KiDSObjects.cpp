//
// KiDSObjects.cpp
//
#include "KiDSObjects.h"
using namespace std;


KiDSObjectList::KiDSObjectList(const string fits_filename, int bitmask,
			       KiDSObjectList::blind blind_index) {
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
  //  Z_B, (Z_B_MIN, Z_B_MAX)
  //  FWHM_IMAGE, weight, and assign it to SourceObject

  valarray<float> g1a;
  valarray<float> g2a;
  valarray<float> g1b;
  valarray<float> g2b;
  valarray<float> g1c;
  valarray<float> g2c;
  valarray<float> wt_a;
  valarray<float> wt_b;
  valarray<float> wt_c;

  int max_src_count;

  try {
    CCfits::Column& column1 = table.column("e1_A");
    max_src_count = column1.rows();
    column1.read( g1a, 1, max_src_count );
    CCfits::Column& column2 = table.column("e2_A");
    column2.read( g2a, 1, max_src_count );
    CCfits::Column& column3 = table.column("weight_A");
    column3.read( wt_a, 1, max_src_count );

    CCfits::Column& column4 = table.column("e1_B");
    column4.read( g1b, 1, max_src_count );
    CCfits::Column& column5 = table.column("e2_B");
    column5.read( g2b, 1, max_src_count );
    CCfits::Column& column6 = table.column("weight_B");
    column6.read( wt_b, 1, max_src_count );

    CCfits::Column& column7 = table.column("e1_C");
    column7.read( g1c, 1, max_src_count );
    CCfits::Column& column8 = table.column("e2_C");
    column8.read( g2c, 1, max_src_count );
    CCfits::Column& column9 = table.column("weight_C");
    column9.read( wt_c, 1, max_src_count );


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


  // append objects to this list
  source_list.reserve(max_src_count);
  for (int i=0; i<max_src_count; ++i) {
      if (wt_a[i] == 0) continue;
      if (mask[i] & bitmask) continue;
      KiDSObject* ptr = new KiDSObject(i, ra[i], dec[i], mag[i], xpos[i], ypos[i], fwhm[i],
				       g1a[i], g2a[i], g1b[i], g2b[i], g1c[i], g2c[i],
				       wt_a[i], wt_b[i], wt_c[i], sn[i], z_B[i], mask[i], blind_index);
      source_list.push_back(ptr);
  }

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
KiDSObjectList::setBlinding(KiDSObjectList::blind blind_index){

  if ((blind_index<0) || (blind_index >=KiDSObject::NUM_SHEAR))
    throw KiDSObjectsError("setBlinding(): wrong blinding specified");

  for (vector<KiDSObject*>::iterator it = source_list.begin();
       it != source_list.end(); ++it) {
    Shear s = (*it)->getShearArray()[blind_index];
    double g1 = (*it)->getG1Array()[blind_index];
    double g2 = (*it)->getG2Array()[blind_index];
    double wt = (*it)->getWeightArray()[blind_index];
    (*it)->setShearG1G2BiasCorrections(s, g1, g2, wt);
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
