//
// gglens_gama_kids.cpp   : For measuring the GGLens signal around gama objects
//
#include <iostream>
#include <cstdlib>
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

const string suffix = ".dat";

const string usage =
  "\n"
  "gglens_gama_kids: calculate tangential shear around a given central point (halo center)\n"
  "\n"
  " usage: gglens_gama_kids <lens_cat> <src_cat> <src_specz_cat> <src_bitmask> <src_blind_index>\n"
  "         <radial_bin_info> <stellar_mass_bin_info> <random_shear_file> <outfile_prefix>\n"
  "  lens_cat:        lens catalog which contains the columns\n"
  "  src_cat:         source catalog, which contains the columns\n"
  "  src_specz_cat:   spec_z catalog which determines the p(z) of the source catalog\n"
  "  src_bitmask:     bitmask to apply to the source catalog.  Default = [0x7c14]\n"
  "  src_blind_index: choose blinding for source catalog.  option = [0,1,2]\n"
  "  radial_bin_info: radial bin info (3 numbers, in arcmin): [min_angle, max_angle, rad_nbin]\n"
  "  stellar_mass_bin_info: stellar mass bin info ([bin edges in log(M/Msun)]\n"
  "  random_shear_file: TreeCorr output for randoms, for use in subtracting random signals\n"
  "  outfile_prefix:  prefix for the output file (suffix is "+suffix+")\n"
  "  \n"
  " stdin:  (none)\n"
  " stdout: (none)\n";


/// source redshift cuts (based on z_B)   TODO: These should be specifiable from a parameter file
const double MIN_SRC_ZB = 0.005;
const double MAX_SRC_ZB = 1.2;
const double MIN_LENS_SRC_SEP = 0.15;
const double h = 1.0;  // H0 = 100 h km/s/Mpc


int
main(int argc, char* argv[]) {

  try {

    //
    // process arguments (TODO: convert to option flags)
    //
    if (argc != 10) {
      cerr << usage;
      exit(2);
    }
    int iarg = 0;
    const string lens_filename = argv[++iarg];
    const string source_filename = argv[++iarg];
    const string specz_src_filename = argv[++iarg];
    const string src_bitmask_str = argv[++iarg];         // remove bitmask masked objs.
    const int    src_bitmask = strtol(src_bitmask_str.c_str(), NULL, 0); // convert hex/dec str
    const int    src_blind_index = atoi(argv[++iarg]);   // choose blinding options.
    const string radial_bin_filename = argv[++iarg];     // radial mass bin
    const string sm_bin_filename = argv[++iarg];         // stellar mass bin
    const string random_shear_filename = argv[++iarg];
    const string outf_prefix = argv[++iarg];
    
    /// open lens file
    ifstream lensf(lens_filename.c_str());
    if (!lensf)
      throw MyException("lens catalog file " + lens_filename + " not found");
    lensf.close();

    /// check source filenames
    ifstream sourcef(source_filename.c_str());
    if (!sourcef)
      throw MyException("source catalog file " + source_filename + " not found");
    sourcef.close();
    ifstream speczf(specz_src_filename.c_str());
    if (!speczf)
      throw MyException("source specz catalog file " + specz_src_filename + " not found");
    speczf.close();


    //
    // setup radial bins (in Mpc or Mpc/h, depending on definition of h)
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
    // setup bins (magnitude)
    //
    ifstream logmstarf(sm_bin_filename.c_str());
    ArbitraryWidthBins logmstar_bin(logmstarf);

    //
    // setup lens sample
    //
    GAMAObjectList master_lens_list(lens_filename);
    GAMAObjectList lens_list(master_lens_list);
    lens_list.applyLogMStarCut(logmstar_bin.getMin(), logmstar_bin.getMax());

    //
    // setup source sample
    //
    valarray<float> pz_list = KiDSObjectList::getPZ(specz_src_filename,
						    MIN_SRC_ZB, MAX_SRC_ZB, src_bitmask);
    KiDSObjectList master_source_list(source_filename, src_bitmask, src_blind_index, pz_list);
    KiDSObjectList source_list(master_source_list);
    source_list.applyRedshiftCut(MIN_SRC_ZB, MAX_SRC_ZB);

    //
    // TODO: read in randoms gglens curve (make a new class)
    //
    /*
    TreeCorrNGOutput random_shear(random_shear_filename);
    */

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
    cerr << "     mask ............. " << "0x" << std::hex << source_list.getBitMask();
    cerr << std::dec << " (" << source_list.getBitMask() << ")" << endl;
    cerr << "     blind ID ......... " << source_list.getBlindIndex() << endl;

    if (source_list.size() == 0) {
      cerr << "no source objects, exiting" << endl;
      return(9);
    }

    cerr << "radial bin range ...... " << radial_bin[0] << " to "
	 << radial_bin[radial_bin.binSize()] << " (arcmin)" << endl;

    cerr << "log(mstar) bin range .. ";
    for (int i=0; i<logmstar_bin.vectorSize(); ++i)  cerr << logmstar_bin[i] << " ";
    cerr << endl;

    //
    // TODO: print randoms gglens curve info
    //
    /*
    cerr << "random shear data ..... " << random_shear_filename << endl;
    cerr << "     number of pairs .. " << random_shear.getNPairs() << endl;
    cerr << "radial bin range ...... " << random_shear.getRadialBin()[0] << " to "
	 << random_shear.getRadialBin()[random_shear.getRadialBin().binSize()]
          << " (arcmin)" << endl;
     */

    //
    // create GGLensObjectList from lens_list and source_list (sums tangential shears for each lens)
    //
    double Om=0.315, Ol=1.-Om;
    cosmology::Cosmology cosmo(Om, Ol);  // unused (but unfortunately necessary)
    bool radialBinIsMpc = false;
    bool normalizeToSigmaCrit = false;
    if (normalizeToSigmaCrit) {
      cerr << "cosmology is .......... " << "(Om=" << Om << ", Ol=" << Ol << ")" << endl;
      cerr << "h is .................. " << h << endl;
    }

    //
    // TODO: (1) GGLens constructor takes randoms info as input (2) change order so cosmo comes last
    //

    GGLensObjectList<GAMAObject*, KiDSObject*>
      gglens_list(lens_list, source_list, radial_bin, radialBinIsMpc, normalizeToSigmaCrit,
		  /*random_shear,*/
		  cosmo, h, MIN_LENS_SRC_SEP);

    //
    // sort each lens into binned_lists
    //

    // make output filename
    std::stringstream sstm;
    sstm << outf_prefix << suffix;    // output filename
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
	  radial_shears[imag][irad].addMBias(
	    (*binned_lists[imag][ilens])[irad].getMBias());
	}
      }
    } // end for(imag)


    //
    // provide output per bin
    //
    ofs << "#imag irad pairs sum(weights) sum(w^2) sum(responsivity) sum(w*et) sum(w*ex) "
	<< "sum(w*var(et)) sum(w*var(ex)) n_lens m_corr" << endl;

    ofs << "#magbins: ";
    for (int imag=0; imag<logmstar_bin.vectorSize(); ++imag) {
      ofs << logmstar_bin[imag] << " ";
    }
    ofs << endl;

    ofs << "#radbins(Mpc/h),h="<< h <<": ";
    for (int irad=0; irad<radial_bin.vectorSize(); ++irad) {
      ofs << radial_bin[irad] << " ";
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
	    << radial_shears[imag][irad].getLensCounts() << " "
	    << radial_shears[imag][irad].getMBias() << endl;
      }
    }

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
