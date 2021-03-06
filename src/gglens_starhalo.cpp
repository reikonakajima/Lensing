//
// gglens_starhalo.cpp   : For measuring the GGLens signal around a stellar reflection halo
//
#include <iostream>
#include "Std.h"
#include "StringStuff.h"
#include "Bounds.h"
#include "StarMaskObjects.h"
#include "KiDSObjects.h"
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
  "gglens_starhalo: calculate tangential shear around a given central point (halo center)\n"
  "\n"
  " usage: gglens_starhalo <lens_cat> <source_cat> <specz_src_cat> <radial_bin_info> <outfile_prefix> <blinding>\n"
  "  lens_catalog:    lens catalog which contains the columns\n"
  "  source_catalog:  source catalog, which contains the columns\n"
  "  specz_src_catalog:  source specz catalog, which contains the columns\n"
  "  radial_bin_info: radial bin info (3 numbers, in arcseconds): [min_theta, max_theta, rad_nbin]\n"
  "  blinding:        A, B, C, or D (for 4 different blindings)\n"
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
    if (argc != 7) {
      cerr << usage;
      exit(2);
    }
    int iarg = 0;
    const string lens_filename = argv[++iarg];
    const string source_filename = argv[++iarg];
    const string specz_src_filename = argv[++iarg];
    const string radial_bin_filename = argv[++iarg];
    const string outf_prefix = argv[++iarg];
    const char blinding = argv[++iarg][0];
    
    /// check blinding
    if (blinding != 'A' and blinding != 'B' and blinding != 'C' and blinding != 'D')
      throw MyException("blinding not within specified range");

    /// open lens file
    ifstream lensf(lens_filename.c_str());
    if (!lensf) 
      throw MyException("lens catalog file " + lens_filename + " not found");

    /// construct source filenames, open files
    ifstream sourcef(source_filename.c_str());
    if (!sourcef) 
      throw MyException("source catalog file " + source_filename + " not found");
    sourcef.close();

    //
    // setup radial bins (in arcsec)
    //
    ifstream radialbinf(radial_bin_filename.c_str());
    double min_theta = 5.0;
    double max_theta = 4000.0;
    int rad_nbin = 50;
    if (radialbinf) {
      if (!(radialbinf >> min_theta >> max_theta >> rad_nbin))
	throw MyException("radialbin file type error");
      if (rad_nbin < 2 || min_theta < 0 || max_theta < min_theta)
	throw MyException("radialbin file specification error");
    }
    LogarithmicBins radial_bin_arcsec(min_theta, max_theta, rad_nbin);  // for debug/info
    LogarithmicBins radial_bin(min_theta, max_theta, rad_nbin);  // keep in arcsec

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
    KiDSObjectList master_source_list(source_filename);
    KiDSObjectList source_list(master_source_list);  // TODO: add any extra cuts
    //source_list.usePixelCoord(true);  // must use ra/dec for lens and source
    //source_list.setBlinding(blinding);

    //
    // diagnostic error messages
    //
    cerr << "=== gglens_starhalo ===" << endl;
    cerr << "lens catalog .......... " << lens_filename << endl;
    cerr << "     count ............ " << lens_list.size() << "/"
	 << master_lens_list.size() << endl;
    cerr << "     bounds ........... " << lens_list.getBounds() << endl;

    if (lens_list.size() == 0) {
      cerr << "no lens objects, exiting" << endl;
      return(9);
    }

    cerr << "source catalog ........ " << source_filename << endl;
    cerr << "     count ............ " << source_list.size() << "/"
	 << master_source_list.size() << endl;
    cerr << "     bounds ........... " << source_list.getBounds() << endl;

    if (source_list.size() == 0) {
      cerr << "no source objects, exiting" << endl;
      return(9);
    }

    cerr << "radial bin range (\") .. " << radial_bin_arcsec[0] << " ... "
	 << radial_bin_arcsec[radial_bin_arcsec.binSize()] << endl;

    cerr << "magnitude bin range ... ";
    for (int i=0; i<magnitude_bin.vectorSize(); ++i)  cerr << magnitude_bin[i] << " ";
    cerr << endl;


    //
    // create GGLensObjectList from lens_list and source_list (sums tangential shears for each lens)
    //
    bool radialBinIsMpc = false;
    bool normalizeToSigmaCrit = false;
    GGLensObjectList<StarMaskObject*, KiDSObject*> gglens_list(lens_list, source_list, radial_bin,
							       radialBinIsMpc, normalizeToSigmaCrit);

    //
    // sort each lens into binned_lists
    //

    /// first, sort by halo type
    const int num_i_type = 5;    /// THIS IS A HARD-CODED ITEM; SHOULD BE FIXED
    for (int i_type = 0; i_type < num_i_type; ++i_type) {

      // make output filename
      std::stringstream sstm;
      sstm << outf_prefix << "_type" << i_type << ".dat";
      string out_filename = sstm.str();
      ofstream ofs(out_filename.c_str());

      vector<GGLensObjectList<StarMaskObject*, KiDSObject*> > binned_lists  // vector over <mags>
        = gglens_list.splitList(magnitude_bin.binSize());  // initialize the split lists
      // fill in the split lists
      for (int ilens = 0; ilens < gglens_list.size(); ++ilens) {
	int index = magnitude_bin.getIndex(gglens_list[ilens]->getLensPtr()->getMag());
	long int halo_type = gglens_list[ilens]->getLensPtr()->getType();
	if (index != -1 and halo_type == i_type) {
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
	    radial_shears[imag][irad].addLensCounts();  // increment by one lens
	    radial_shears[imag][irad].addPairCounts(    // increment by source counts
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
      ofs << "#imag irad pairs sum(weights) sum(w^2) sum(responsivity) sum(w*et) sum(w*ex) "
	  << "sum(w*var(et)) sum(w*var(ex)) n_lens" << endl;

      ofs << "#magbins: ";
      for (int imag=0; imag<magnitude_bin.vectorSize(); ++imag) {
	ofs << magnitude_bin[imag] << " ";
      }
      ofs << endl;

      ofs << "#radbins(arcsec): ";
      for (int irad=0; irad<radial_bin_arcsec.vectorSize(); ++irad) {
	ofs << radial_bin_arcsec[irad] << " ";
      }
      ofs << endl;

      for (int imag=0; imag<magnitude_bin.binSize(); ++imag) {
	for (int irad=0; irad<radial_bin.binSize(); ++irad) {
	  if (radial_shears[imag][irad].getPairCounts() == 0)
	    continue;
	  ofs << imag << " " << irad << "  ";
	  ofs << radial_shears[imag][irad].getPairCounts() << " "
	      << radial_shears[imag][irad].getWeights() << " "
	      << radial_shears[imag][irad].getWeightSq() << " "
	      << radial_shears[imag][irad].getResponsivity() << " "
	      << radial_shears[imag][irad].getDeltaSigma_t() << " "
	      << radial_shears[imag][irad].getDeltaSigma_s() << " "
	      << radial_shears[imag][irad].getVariance_t() << " "
	      << radial_shears[imag][irad].getVariance_s() << " "
	      << radial_shears[imag][irad].getLensCounts() << endl;
	}
      }

    } // end for(i_type)


  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
