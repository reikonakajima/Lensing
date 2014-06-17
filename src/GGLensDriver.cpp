//
// GGLensDriver.cpp   : For testing the GGLens C++ module
//
#include <iostream>
#include "Std.h"
#include "StringStuff.h"
#include "Bounds.h"
#include "LensObjects.h"
#include "StarMaskObjects.h"
#include "SourceObjects.h"
#include "RCSLenSObjects.h"
#include "GGLens.h"
#include "Bins.h"
using std::ostringstream;
using std::setw;
using std::setfill;
using std::fixed;
using std::setprecision;


const string outfprefix = "gglens";
const string configfname = "config.par";
const string suffix = ".dat";

const string usage =
  "\n"
  "GGLensDriver: calculate tangential shear around a point (galaxy-galaxy lensing signal)\n"
  "\n"
  " usage: GGLensDriver <lens_catalog> <source_catalog>\n"
  "  lens_catalog:  lens catalog which containts the columns\n"
  "  r.txt:   radial bin info (3 numbers): [min_theta, max_theta, rad_nbin]\n"
  "  \n"
  "  source_catalog:  source catalog, which contains the columns\n"
  "  \n"
  //  " output #1: file name:\" "+outfprefix+suffix+"\"\n"
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
    // setup radial bins (in pixel units)
    //
    ifstream radialbinf("r.txt");
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

    //
    // setup bins (magnitude)
    //

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
    StarMaskObjectList master_lens_list(lensf);
    StarMaskObjectList lens_list(lensf);     // TODO:  FIXME  !!!
    // Needs: 
    //   list(==vector):  bounds, 
    //   each object: position, shear, resolution, (optional: redshift, magnitude)
    //                may have multiple bands
    RCSLenSObjectList master_source_list(source_filename);
    RCSLenSObjectList source_list(source_filename);  // TODO: add any extra cuts


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


    //
    // create GGLensObjectList from lens_list and source_list
    //
    GGLensObjectList<StarMaskObject*, RCSLenSObject*> gglens_list(lens_list, source_list, radial_bin);

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
