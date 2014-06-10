//
// GGLensDriver.cpp   : For testing the GGLens C++ module
//
#include <iostream>
#include "Std.h"
#include "StringStuff.h"
#include "Bounds.h"
#include "LensObjects.h"
//#include "ggLensObjects.h" 
#include "ggLensSum.h"
//#include "photCatObjects.h"
//#include "bgCatObjects.h"
//#include "Cosmology.h"
#include "Bins.h"
//#include "StellarMass.h"
#include "Mesh.h"
//#include "AstronomicalConstants.h"
//#include "SigmaCrit.h"
using std::ostringstream;
using std::setw;
using std::setfill;
using std::fixed;
using std::setprecision;


const string outfprefix = "gglens";
const string configfname = "config.par";
const string suffix = ".dat";
const bool flat = true;

const string usage =
  "\n"
  "GGLensDriver: calculate tangential shear around a point (galaxy-galaxy lensing signal)\n"
  "\n"
  " usage: GGLensDriver <lens_catalog> <source_catalog>\n"
  "  lens_catalog:  lens catalog which containts the columns\n"
  "  \n"
  "  source_catalog:  source catalog, which contains the columns\n"
  "  \n"
  " output #1: file name:\" "+outfprefix+suffix+"\"\n"
  " stdin:  (none)\n"
  " stdout: (none)\n";


int
main(int argc, char* argv[]) {

  try {

    //
    // process arguments
    //
    if (argc != 3) {
      cerr << usage;
      exit(2);
    }
    int iarg = 0;
    const string lens_filename = argv[++iarg];
    const string source_filename = argv[++iarg];
    
    /// open lens file
    ifstream lensf(lens_filename.c_str());
    if (!lensf) 
      throw MyException("lens catalog file " + lens_filename + " not found");

    /// construct source filenames, open files
    ifstream sourcef(source_filename.c_str());
    if (!sourcef) 
      throw MyException("source catalog file " + source_filename + " not found");

    //
    // setup bins (magnitude)
    //

    ifstream radialbinf("r.txt");

    /// setup radial bins (in pixel units)
    double min_theta = 5.0;
    double max_theta = 300.0;
    int rad_nbin = 21;
    if (radialbinf) {
      if (!(radialbinf >> min_theta >> max_theta >> rad_nbin))
	throw MyException("radialbin file type error");
      if (rad_nbin < 2 || min_theta < 0 || max_theta < min_theta)
	throw MyException("radialbin file specification error");
    }
    LogarithmicBins radial_bin(min_theta, max_theta, rad_nbin);

    /// create lens magnitude bins (r-band magnitude)
    int mag_nbins = 3;
    vector<double> mvec(mag_nbins+1);
    for (int i=0; i<=mag_nbins; ++i) {  // TODO:  FIXME  !!!
      mvec[i] = i + 10;
    }
    ArbitraryWidthBins magnitude_bin(mvec);
    const int nmagbin = magnitude_bin.size(); // remove?


    //
    // setup lens/source samples
    //

    // Needs: 
    //   list(==vector):  bounds, 
    //   each object: position, magnitude, (optional: redshift, sed type)
    //                may have multiple bands
    LensObjectList master_lens_list(lensf);    // TODO:  MAKE ME  !!!
    LensObjectList lens_list(lensf);     // TODO:  FIXME  !!!
    // Needs: 
    //   list(==vector):  bounds, 
    //   each object: position, shear, resolution, (optional: redshift, magnitude)
    //                may have multiple bands
    LensObjectList master_source_list(sourcef);  // TODO:  MAKE ME  !!!
    LensObjectList source_list(sourcef);    // TODO:  FIXME  !!!


    //
    // diagnostic error messages
    //
    cerr << "=== GGLensDriver ===" << endl;
    cerr << "lens catalog ...... " << lens_filename << endl;
    cerr << "     count ........ " << lens_list.size() << "/" << master_lens_list.size() << endl;
    cerr << "     bounds ....... " << lens_list.getBounds() << endl;
    //cerr << "     rmag range ... " << lens_list.getMaxMag << endl;

    if (lens_list.size() == 0) {
      return(9);
    }

    cerr << "source catalog .... " << source_filename << endl;
    cerr << "     count ........ " << source_list.size() << "/" << master_source_list.size() << endl;
    cerr << "     bounds ....... " << source_list.getBounds() << endl;

    if (source_list.size() == 0) {
      return(9);
    }

    cerr << "radial bin range .. " << radial_bin[0] << " ... " 
	 << radial_bin[radial_bin.size()-1] << endl;
    cerr << "magnitude bin range " << magnitude_bin[0] << " ... " 
	 << magnitude_bin[magnitude_bin.size()-1] << endl;


    /// create source mesh list
    Bounds<double> srcbounds = source_list.getBounds();
    vector<LensObject*> srcvector = source_list.getVectorForm();   // TODO:  FIXME !!!
    double ramin = srcbounds.getXMin();
    double ramax = srcbounds.getXMax();
    double decmin = srcbounds.getYMin();
    double decmax = srcbounds.getYMax();
    const int meshdimX = static_cast<int>((ramax-ramin) / (10. / 3600.));  // rawidth / 10 as
    const int meshdimY = static_cast<int>((decmax-decmin) / (10. / 3600.));// decwidth / 10 as
    const bool isPeriodic = false;          // we want a non-periodic mesh
    Mesh<LensObject*> srcmesh(meshdimX, meshdimY, 1, srcvector, isPeriodic,
			      ramin, ramax, decmin, decmax);

    //
    // iterate over all lens; generate lensing signal profile for each lens object
    //


    LensObjectList::iterator it = lens_list.begin();
    for (; it != lens_list.end(); ++it) {

      LensObject* lensobj = *it;

      // we're running systematics tests, so we'll be in units of pixels for distances
      double lensra = lensobj->getRA();   // in pixels
      double lensdec = lensobj->getDec(); // in pixels

      // print out magnitude for binning later on
      float mag = lensobj->getMag();

      // DEBUG
      cerr << " === lens radec: " << lensra << ", " << lensdec << endl;
      cerr << "mag bin #: " << mag  << endl;
      

      /// find all possible matching sources (get their indicies of srcvector) in radial bins
      vector<multimap<double, int> > bglist(rad_nbin+1);
      for (int irad = 0; irad <= rad_nbin; ++irad) {
	bglist[irad] = srcmesh.getNearAngleMap(lensra, lensdec, 0., 
					       radial_bin[irad+1], radial_bin[irad]);
	// DEBUG radial_bin
	//if (bglist[irad].size() > 0)
	//  cerr << irad << "th radial bin, bg obj count: " << bglist[irad].size() << endl;  
	//cerr << "limits are " << radial_bin[irad] << " to " << radial_bin[irad+1] << endl;
      }

      //
      // LOOP over sources for this lens object
      //

      int bg_count = 0;  // source object count for this lens object
      int lost_et = 0;
      int bad_src = 0;
      int lost_bgcount = 0;  // not enough BG count

      multimap<double, int>::const_iterator isrc;
      vector<ggLensSum> shear_sum(rad_nbin);  // shear data storage
      for (int irad = 0; irad <= rad_nbin; ++irad) {
	for (isrc = bglist[irad].begin(); isrc != bglist[irad].end(); ++isrc) {
	  
	  /// calculate tangential/skew shear
	  LensObject* srcobj = srcvector[isrc->second];
	  double sra = srcobj->getRA();   // in pixels
	  double sdec = srcobj->getDec(); // in pixels
	  double dra = sra - lensra;
	  double ddec = sdec - lensdec;
	  double theta;
	  if (flat) {
	    theta = hypot(dra, ddec);  // 2d euclidean
	  } else {
	    theta = atan2(cos(lensdec)*sin(sdec)-sin(lensdec)*cos(sdec)*cos(dra),
			  cos(sdec)*sin(dra));      // spherical surface
	  }

	  double ct = cos(theta);
	  double st = sin(theta);
	  double c2t = ct*ct-st*st;
	  double s2t = 2.*ct*st;
	  
	  double et = -c2t * srcobj->getE1() - s2t * srcobj->getE2();
	  double es =  s2t * srcobj->getE1() - c2t * srcobj->getE2();
	  
	  // DEBUG
	  cerr << bg_count << "   "
	       << fixed << setprecision(6) 
	       << lensra << " " << lensdec << "   "
	       << srcobj->getRA() << " " << srcobj->getDec() << "  "
	       << srcobj->getE1() << " " << srcobj->getE2() << "  "
	       << et << " " << es << "   "
	       << theta / DEGREE << " " 
	       << endl;
	  
	  if (isnan(et)) {
	    lost_et++;
	    continue;
	  }
	  if (srcobj->getERms() <= 0. || srcobj->getShapeError() <= 0.) {
	    bad_src++;
	    continue;
	  }
	  
	  /// calculate weight
	  double vare = srcobj->getShapeError() * srcobj->getShapeError(); 
	  /////////////  CHECK  (need to calculate shape noise properly)
	  double varSN = srcobj->getERms() * srcobj->getERms(); 
	  double invshapeweight = (vare + varSN);
	  double weight = 1 / invshapeweight;
	  

	  /// calculate other quantities
	  double k0 = varSN*vare/(varSN+vare);
	  double k1 = varSN/(varSN+vare);
	  k1 *= k1;
	  double weightedsignal_t = et / invshapeweight; // ok  * invsigcrit
	  double weightedsignal_s = es / invshapeweight; // ok  * invsigcrit
	  /////////////  CHECK  (where does this responsivity come from?)
	  double responsiv = weight * (1. - k0 - k1*et*et);           // ok
	  double weightederror_t = weight / invshapeweight * et * et; // ?
	  double weightederror_s = weight / invshapeweight * es * es; // ?
	  double weightedinvsigcrit = weight;  // * invsigcrit
	  double weightedrmag = weight * lensobj->getMag();
	  
	  /// add source object to sm bin 
	  shear_sum[irad].addPairCounts();
	  shear_sum[irad].addWeight(weight);
	  shear_sum[irad].addResponsivity(responsiv);
	  shear_sum[irad].addDeltaSigma_t(weightedsignal_t);
	  shear_sum[irad].addDeltaSigma_s(weightedsignal_s);
	  shear_sum[irad].addError_t(weightederror_t);
	  shear_sum[irad].addError_s(weightederror_s);
	  shear_sum[irad].addApparentMag(weightedrmag);
	} // END source irad 'for' loop
	
	/// count available source objects for this lens
	bg_count += shear_sum[irad].getPairCounts();

      } // [end of first source irad 'for' loop]
      
      /// make sure there are enough bg counts for this object
      if (bg_count < 10) {
	lost_bgcount++;
	continue; 
      }

      //
      // print output
      //
      cout << shear_sum[irad].getPairCounts() << " "
	   << shear_sum[irad].getWeights() << " "
	   << shear_sum[irad].getResponsivity() << " "
	   << shear_sum[irad].getDeltaSigma_t()	 << " "
	   << shear_sum[irad].getDeltaSigma_s()	 << " "
	   << shear_sum[irad].getError_t() << " "
	   << shear_sum[irad].getError_s() << " "
	   << shear_sum[irad].getApparentMag();

    } // [end of lens objects for loop]


    // DEBUG
    cerr << "lost to et: " << lost_et << endl;
    cerr << "lost to source obj: " << bad_src << endl;
    cerr << "lost to bg_count: " << lost_bgcount << endl;


  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
