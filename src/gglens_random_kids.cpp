//
// gglens_random_kids.cpp   : For measuring the GGLens signal around a random objects
//
#include <iostream>
#include "Std.h"
#include "StringStuff.h"
#include "Bounds.h"
#include "RandomObjects.h"
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
  "gglens_random_kids: calculate tangential shear around a given central point (halo center)\n"
  "\n"
  " usage: gglens_random_kids <lens_catalog> <source_catalog> <radial_bin_info> <outfile_prefix> <max_num_randoms>\n"
  "  lens_catalog:    lens catalog which contains the columns\n"
  "  source_catalog:  source catalog, which contains the columns\n"
  "  radial_bin_info: radial bin info (3 numbers, in arcseconds): [min_theta, max_theta, rad_nbin]\n"
  "  max_num_randoms: maximum number of randoms to use to calculate signal\n"
  "  \n"
  //  " output #1: file name:\" "+outfprefix+suffix+"\"\n"
  " stdin:  (none)\n"
  " stdout: (none)\n";

/// source redshift cuts (based on z_B)   TODO: These should be specifiable from a parameter file
const double MIN_SRC_Z = 0.005;
const double MAX_SRC_Z = 1.2;
const double MIN_LENS_SRC_SEP = 0.15;
const double h = 1.0;  // H0 = 100 h km/s/Mpc

int
main(int argc, char* argv[]) {

  try {

    //
    // process arguments
    //
    if (argc != 6) {
      cerr << usage;
      exit(2);
    }
    int iarg = 0;
    const string lens_filename = argv[++iarg];
    const string source_filename = argv[++iarg];
    const string radial_bin_filename = argv[++iarg];
    const string outf_prefix = argv[++iarg];
    const int    max_lens_count = atoi(argv[++iarg]);
    
    /// construct source filenames, open files
    ifstream sourcef(source_filename.c_str());
    if (!sourcef) 
      throw MyException("source catalog file " + source_filename + " not found");

    //
    // setup radial bins (in arcsec)
    //
    ifstream radialbinf(radial_bin_filename.c_str());
    double min_Mpc = 1e-2;  // 10 kpc/h minimum
    double max_Mpc = 10;    // 10 Mpc/h maximum
    int rad_nbin = 16;
    if (radialbinf) {
      if (!(radialbinf >> min_Mpc >> max_Mpc >> rad_nbin))
	throw MyException("radialbin file type error");
      if (rad_nbin < 2 || min_Mpc < 0 || max_Mpc < min_Mpc)
	throw MyException("radialbin file specification error");
    }
    // note that radial_bin will eventually need to be in degrees for use with the Mesh class
    LogarithmicBins radial_bin(min_Mpc, max_Mpc, rad_nbin);


    //
    // setup lens/source samples
    //
    GAMARandomObjectList master_lens_list(lens_filename.c_str(), max_lens_count);
    GAMARandomObjectList lens_list(master_lens_list);     // TODO: add any extra cuts

    KiDSObjectList master_source_list(source_filename);
    KiDSObjectList source_list(master_source_list);  // TODO: add any extra cuts
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

    cerr << "radial bin range ...... " << radial_bin[0] << " ... "
	 << radial_bin[radial_bin.binSize()] << " (Mpc/h)" << endl;


    //
    // create GGLensObjectList from lens_list and source_list (sums tangential shears for each lens)
    //
    double Om=0.315, Ol=1.-Om;
    cosmology::Cosmology cosmo(Om, Ol);
    bool radialBinIsMpc = true;
    bool normalizeToSigmaCrit = true;
    if (normalizeToSigmaCrit) {
      cerr << "cosmology is .......... " << "(Om=" << Om << ", Ol=" << Ol << ")" << endl;
      cerr << "h is .................. " << h << endl;
    }

    GGLensObjectList<GAMARandomObject*, KiDSObject*>
      gglens_list(lens_list, source_list, radial_bin, radialBinIsMpc, normalizeToSigmaCrit,
		  cosmo, h, MIN_LENS_SRC_SEP);

    //
    // make output filename
    //
    std::stringstream sstm;
    sstm << outf_prefix << ".dat";    // output filename
    string out_filename = sstm.str();
    ofstream ofs(out_filename.c_str());


    //
    // sum tangential shear according to the given lens binning
    //
    vector<ggLensSum> radial_shears(vector<ggLensSum>(radial_bin.binSize()));

    /// sum tangential shear quantities over all lenses
    for (int ilens = 0; ilens < gglens_list.size(); ++ilens) {
      for (int irad = 0; irad < radial_bin.binSize(); ++irad) {
	radial_shears[irad].addLensCounts();  // increment by one lens
	radial_shears[irad].addPairCounts(    // increment by source counts
	  (*gglens_list[ilens])[irad].getPairCounts());
	double w = (*gglens_list[ilens])[irad].getWeights();
	radial_shears[irad].addWeight(w);
	radial_shears[irad].addWeightSq(w*w);
	radial_shears[irad].addResponsivity(
	  (*gglens_list[ilens])[irad].getResponsivity());
	radial_shears[irad].addDeltaSigma_t(
	  (*gglens_list[ilens])[irad].getDeltaSigma_t());
	radial_shears[irad].addDeltaSigma_s(
	  (*gglens_list[ilens])[irad].getDeltaSigma_s());
	radial_shears[irad].addVariance_t(
	  (*gglens_list[ilens])[irad].getVariance_t());
	radial_shears[irad].addVariance_s(
	  (*gglens_list[ilens])[irad].getVariance_s());
	radial_shears[irad].addMBias(
	  (*gglens_list[ilens])[irad].getMBias());
      }
    }
    

    //
    // provide output per bin
    //
    ofs << "#imag irad pairs sum(weights) sum(w^2) sum(responsivity) sum(w*et) sum(w*ex) "
	<< "sum(w*var(et)) sum(w*var(ex)) n_lens m_corr" << endl;

    ofs << "#radbins(Mpc/h): ";
    for (int irad=0; irad<radial_bin.vectorSize(); ++irad) {
      ofs << radial_bin[irad] << " ";
    }
    ofs << endl;

    const int imag = 0;  // no magnitude binning in random objects
    for (int irad=0; irad<radial_bin.binSize(); ++irad) {
      if (radial_shears[irad].getPairCounts() == 0)
	continue;
      ofs << imag << " " << irad << "  ";
      ofs << radial_shears[irad].getPairCounts() << " "
	  << radial_shears[irad].getWeights() << " "
	  << radial_shears[irad].getWeightSq() << " "
	  << radial_shears[irad].getResponsivity() << " "
	  << radial_shears[irad].getDeltaSigma_t() << " "
	  << radial_shears[irad].getDeltaSigma_s() << " "
	  << radial_shears[irad].getVariance_t() << " "
	  << radial_shears[irad].getVariance_s() << " "
	  << radial_shears[irad].getLensCounts() << " "
	  << radial_shears[irad].getMBias() << endl;
    }

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
