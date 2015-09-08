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
  "  radial_bin_info: radial bin info (3 numbers, in Mpc/h): [min_Mpch, max_Mpch, rad_nbin]\n"
  "  \n"
  "  (future TODO:)\n"
  "  (allow stellar mass or luminosity binning as parameter file option)\n"
  "  (allow h = H0/100 km/s/Mpc as parameter file option)\n"
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

    /*/
    // setup bins (magnitude)
    //
    ifstream magbinf("absmag.txt");
    ArbitraryWidthBins magnitude_bin(magbinf);
    /*/

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

    int bitmask = 1;  /// remove bitmask masked objs. TODO: UPDATE BITMASK OPTIONS  2^15 - 1 = 32767
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

    cerr << "radial bin range ...... " << radial_bin[0] << " to "
	 << radial_bin[radial_bin.binSize()] << " (Mpc/h)" << endl;

    cerr << "log(mstar) bin range .. ";
    for (int i=0; i<logmstar_bin.vectorSize(); ++i)  cerr << logmstar_bin[i] << " ";
    cerr << endl;


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

    /*/
    // DEBUG
    cerr << "cosmo test (comoving distanes Dc):" << endl;
    cerr << "z=0.1: " << cosmo.Dc(0.1) * HubbleLengthMpc / h << endl;
    cerr << "z=0.2: " << cosmo.Dc(0.2) * HubbleLengthMpc / h << endl;
    cerr << "z=0.3: " << cosmo.Dc(0.3) * HubbleLengthMpc / h << endl;
    cerr << "z=0.4: " << cosmo.Dc(0.4) * HubbleLengthMpc / h << endl;
    cerr << "z=0.5: " << cosmo.Dc(0.5) * HubbleLengthMpc / h << endl;
    cerr << endl;
    double one_deg_in_radian = 1.0 * DEGREE;
    cerr << "1deg at z=0.1: " << cosmo.Dc(0.1) * HubbleLengthMpc / h * one_deg_in_radian << endl;
    cerr << "1deg at z=0.2: " << cosmo.Dc(0.2) * HubbleLengthMpc / h * one_deg_in_radian << endl;
    cerr << "1deg at z=0.3: " << cosmo.Dc(0.3) * HubbleLengthMpc / h * one_deg_in_radian << endl;
    cerr << "1deg at z=0.4: " << cosmo.Dc(0.4) * HubbleLengthMpc / h * one_deg_in_radian << endl;
    cerr << "1deg at z=0.5: " << cosmo.Dc(0.5) * HubbleLengthMpc / h * one_deg_in_radian << endl;
    cerr << argv[0] << " DEBUG END" << endl;
    exit(1);
    // DEBUG END
    /*/

    GGLensObjectList<GAMAObject*, KiDSObject*>
      gglens_list(lens_list, source_list, radial_bin, radialBinIsMpc, normalizeToSigmaCrit,
		  cosmo, h, MIN_LENS_SRC_SEP);

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
