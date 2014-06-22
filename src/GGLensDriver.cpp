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
    double max_theta = 400.0;
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
    ifstream magbinf("mag.txt");
    ArbitraryWidthBins magnitude_bin(magbinf);


    //
    // setup lens/source samples
    //

    // Needs: 
    //   list(==vector):  bounds, 
    //   each object: position, magnitude, (optional: redshift, sed type)
    //                may have multiple bands
    StarMaskObjectList master_lens_list(lensf);
    StarMaskObjectList lens_list(master_lens_list);     // TODO: add any extra cuts
    // Needs: 
    //   list(==vector):  bounds, 
    //   each object: position, shear, resolution, (optional: redshift, magnitude)
    //                may have multiple bands
    RCSLenSObjectList master_source_list(source_filename);
    RCSLenSObjectList source_list(master_source_list);  // TODO: add any extra cuts
    source_list.usePixelCoord(true);  // use pixel coordinates


    //
    // diagnostic error messages
    //
    cerr << "=== GGLensDriver ===" << endl;
    cerr << "lens catalog ...... " << lens_filename << endl;
    cerr << "     count ........ " << lens_list.size() << "/" << master_lens_list.size() << endl;
    cerr << "     bounds ....... " << lens_list.getBounds() << endl;

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
	 << radial_bin[radial_bin.binSize()] << endl;
    cerr << "magnitude bin range " << magnitude_bin[0] << " ... " 
	 << magnitude_bin[magnitude_bin.binSize()] << endl;


    //
    // create GGLensObjectList from lens_list and source_list (sums tangential shears for each lens)
    //
    GGLensObjectList<StarMaskObject*, RCSLenSObject*> gglens_list(lens_list, source_list,
								  radial_bin);

    //
    // sort each lens into binned_lists
    //
    vector<GGLensObjectList<StarMaskObject*, RCSLenSObject*> > binned_lists
	= gglens_list.splitList(magnitude_bin.binSize());

    for (int ilens = 0; ilens < gglens_list.size(); ++ilens) {
      int index = magnitude_bin.getIndex(gglens_list[ilens]->getLensPtr()->getMag());
      if (index != -1) {
	binned_lists[index].push_back(gglens_list[ilens]);
      }
    }


    //
    // sum tangential shear according to the given lens binning
    //
    vector<vector<ggLensSum> > radial_shears(magnitude_bin.binSize());
    for (int imag = 0; imag < magnitude_bin.binSize(); ++imag) {
      /// initialize radial_shears
      vector<ggLensSum> temp(radial_bin.binSize());
      radial_shears[imag] = temp;
      /// sum tangential shear quantities over all lenses
      for (int ilens = 0; ilens < binned_lists[imag].size(); ++ilens) {
	for (int irad = 0; irad < radial_bin.binSize(); ++irad) {
	  radial_shears[imag][irad].addPairCounts(
	      (*binned_lists[imag][ilens])[irad].getPairCounts());
	  double w = (*binned_lists[imag][ilens])[irad].getWeights();
	  radial_shears[imag][irad].addWeight(w);
	  radial_shears[imag][irad].addWeightSq(w*w);
	  radial_shears[imag][irad].addResponsivity(
	      (*binned_lists[imag][ilens])[irad].getResponsivity());
	  radial_shears[imag][irad].addDeltaSigma_t(
	      (*binned_lists[imag][ilens])[irad].getDeltaSigma_t());
	  radial_shears[imag][irad].addDeltaSigma_s(
	      (*binned_lists[imag][ilens])[irad].getDeltaSigma_s());
	  radial_shears[imag][irad].addVariance_t(
	      (*binned_lists[imag][ilens])[irad].getVariance_t());
	  radial_shears[imag][irad].addVariance_s(
	      (*binned_lists[imag][ilens])[irad].getVariance_s());
	}
      }
    }


    //
    // provide output per bin
    //
    cout << "#imag irad pairs sum(weights) sum(w^2) sum(responsivity) sum(w*et) sum(w*ex) "
	 << "sum(w*var(et)) sum(w*var(ex))" << endl;

    cout << "#magbins: ";
    for (int imag=0; imag<magnitude_bin.vectorSize(); ++imag) {
      cout << magnitude_bin[imag] << " ";
    }
    cout << endl;

    cout << "#radbins: ";
    for (int irad=0; irad<radial_bin.vectorSize(); ++irad) {
      cout << radial_bin[irad] << " ";
    }
    cout << endl;

    for (int imag=0; imag<magnitude_bin.binSize(); ++imag) {
      for (int irad=0; irad<radial_bin.binSize(); ++irad) {
	if (radial_shears[imag][irad].getPairCounts() == 0)
	  continue;
	cout << imag << " " << irad << "  ";
	cout << radial_shears[imag][irad].getPairCounts() << " "
	     << radial_shears[imag][irad].getWeights() << " "
	     << radial_shears[imag][irad].getWeightSq() << " "
	     << radial_shears[imag][irad].getResponsivity() << " "
	     << radial_shears[imag][irad].getDeltaSigma_t() << " "
	     << radial_shears[imag][irad].getDeltaSigma_s() << " "
	     << radial_shears[imag][irad].getVariance_t() << " "
	     << radial_shears[imag][irad].getVariance_s() << endl;
      }
    }


  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
