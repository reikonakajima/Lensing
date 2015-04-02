//
// gglens_gama_kids.cpp   : For measuring the GGLens signal around gama objects
//
#include <iostream>
#include "Std.h"
#include "StringStuff.h"
#include "Bounds.h"
#include "GAMAObjects.h"
#include "KiDSObjects.h"
#include "Cosmology.h"
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
  "gglens_gama_kids: calculate tangential shear around a given central point (halo center)\n"
  "\n"
  " usage: gglens_gama_kids <lens_catalog> <source_catalog> <radial_bin_info> <outfile_prefix>\n"
  "  lens_catalog:    lens catalog which contains the columns\n"
  "  source_catalog:  source catalog, which contains the columns\n"
  "  radial_bin_info: radial bin info (3 numbers, in arcseconds): [min_theta, max_theta, rad_nbin]\n"
  "  \n"
  //  " output #1: file name:\" "+outfprefix+suffix+"\"\n"
  " stdin:  (none)\n"
  " stdout: (none)\n";


// source redshift cuts (based on z_B)
const double MIN_SRC_Z = 0.005;
const double MAX_SRC_Z = 1.2;
const double MIN_LENS_SRC_SEP = 0.15;


int
main(int argc, char* argv[]) {

  try {

    //
    // process arguments
    //
    if (argc != 5) {
      cerr << usage;
      exit(2);
    }
    int iarg = 0;
    const string lens_filename = argv[++iarg];
    const string source_filename = argv[++iarg];
    const string radial_bin_filename = argv[++iarg];
    const string outf_prefix = argv[++iarg];
    
    /// open lens file
    ifstream lensf(lens_filename.c_str());
    if (!lensf) 
      throw MyException("lens catalog file " + lens_filename + " not found");

    /// construct source filenames, open files
    ifstream sourcef(source_filename.c_str());
    if (!sourcef) 
      throw MyException("source catalog file " + source_filename + " not found");

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
    LogarithmicBins radial_bin(min_theta, max_theta, rad_nbin);
    radial_bin.rescale(1./3600.);  // convert arcsec to degree

    //
    // setup bins (magnitude)
    //
    ifstream magbinf("absmag.txt");
    ArbitraryWidthBins magnitude_bin(magbinf);

    //
    // setup bins (magnitude)
    //
    ifstream logmstarf("logmstar.txt");
    ArbitraryWidthBins logmstar_bin(logmstarf);

    //
    // setup lens/source samples
    //
    GAMAObjectList master_lens_list(lensf);
    GAMAObjectList lens_list(master_lens_list);
    lens_list.applyLogMStarCut(logmstar_bin.getMin(), logmstar_bin.getMax());

    int bitmask = 1;  /// TODO: UPDATE BITMASK OPTIONS  2^15 - 1 = 32767
    KiDSObjectList master_source_list(source_filename, bitmask);
    KiDSObjectList source_list(master_source_list);
    source_list.applyRedshiftCut(MIN_SRC_Z, MAX_SRC_Z);

    //
    // diagnostic error messages
    //
    cerr << "=== " << argv[0] << " ===" << endl;
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

    cerr << "radial bin range ...... " << radial_bin_arcsec[0] << " ... "
	 << radial_bin_arcsec[radial_bin_arcsec.binSize()] << " (arcsec)" << endl;
    /*
    cerr << "magnitude bin range ... " << magnitude_bin[0] << " ... "
	 << magnitude_bin[magnitude_bin.binSize()] << endl;
    */
    cerr << "log(mstar) bin range ... ";
    for (int i=0; i<logmstar_bin.vectorSize(); ++i)  cerr << logmstar_bin[i] << " ";
    cerr << endl;


    //
    // create GGLensObjectList from lens_list and source_list (sums tangential shears for each lens)
    //
    double Om=0.27, Ol=1.-Om;
    cosmology::Cosmology cosmo(Om, Ol);
    GGLensObjectList<GAMAObject*, KiDSObject*> gglens_list(lens_list, source_list, radial_bin,
							   cosmo);

    //
    // sort each lens into binned_lists
    //

    // make output filename
    std::stringstream sstm;
    sstm << outf_prefix << ".dat";    // output filename
    string out_filename = sstm.str();
    ofstream ofs(out_filename.c_str());

    vector<GGLensObjectList<GAMAObject*, KiDSObject*> > binned_lists  // vector over <mags>
      = gglens_list.splitList(logmstar_bin.binSize());  // initialize the split lists
    // fill in the split lists
    for (int ilens = 0; ilens < gglens_list.size(); ++ilens) {
      int index = logmstar_bin.getIndex(gglens_list[ilens]->getLensPtr()->getLogMStar());
      if (index != -1) {
	binned_lists[index].push_back(gglens_list[ilens]);
      }
    }


    //
    // sum tangential shear according to the given lens binning
    //
    vector<vector<ggLensSum> > radial_shears(logmstar_bin.binSize(),
					     vector<ggLensSum>(radial_bin.binSize()));
    for (int imag = 0; imag < logmstar_bin.binSize(); ++imag) {
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
    } // end for(imag)


    //
    // provide output per bin
    //
    ofs << "#imag irad pairs sum(weights) sum(w^2) sum(responsivity) sum(w*et) sum(w*ex) "
	<< "sum(w^2*var(et)) sum(w^2*var(ex)) n_lens" << endl;

    ofs << "#magbins: ";
    for (int imag=0; imag<logmstar_bin.vectorSize(); ++imag) {
      ofs << logmstar_bin[imag] << " ";
    }
    ofs << endl;

    ofs << "#radbins(arcsec): ";
    for (int irad=0; irad<radial_bin_arcsec.vectorSize(); ++irad) {
      ofs << radial_bin_arcsec[irad] << " ";
    }
    ofs << endl;

    for (int imag=0; imag<logmstar_bin.binSize(); ++imag) {
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

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
