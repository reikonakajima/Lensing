//
// GGLensDriver.cpp   : For testing the GGLens C++ module
//
#include <iostream>
#include "Std.h"
#include "StringStuff.h"
#include "Bounds.h"
#include "LensObjects.h"
//#include "ggLensObjects.h" 
//#include "ggLensSum.h"
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


const string outfprefix = "gglens.";
const string configfname = "config.par";
const string suffix = ".dat";


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

    /// setup radial bins (in arcseconds)
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

    /*

    //
    // setup summation buffer for nested bins
    //

    /// create summation tally for the nested bins
    const string lumlabel = "l";
    const string morphlabel = "t";
    const string zlabel = "z";
    vector<vector<vector<vector<GGLensSums> > > >tally(nlumbin);  // tally[lum][morph][z][rad]
    for (int l=0; l < nlumbin; ++l) {
      vector<vector<vector<GGLensSums> > > temp1(nmorphbin);
      tally[l] = temp1;
      for (int m=0; m < nmorphbin; ++m) {
	vector<vector<GGLensSums> > temp2(nzbin);
	tally[l][m] = temp2;
	for (int z=0; z < nzbin; ++z) {
	  vector<GGLensSums> temp3(nradbin, GGLensSums(lumlabel));
	  tally[l][m][z] = temp3;
	}
      }
    }
    vector<vector<vector<LensCounts> > > lenscount(nlumbin);  // lenscount[lum][morph][z]
    for (int l=0; l < nlumbin; ++l) {
      vector<vector<LensCounts> > temp1(nmorphbin);
      lenscount[l] = temp1;
      for (int m=0; m < nmorphbin; ++m) {
	vector<LensCounts> temp2(nzbin);
	lenscount[l][m] = temp2;
      }
    }

    cerr << "nlumbin, nmorphbin, nzbin, nradbin: " << nlumbin << ", " << nmorphbin << ", "
	 << nzbin << ", " << nradbin << endl;

    //
    // setup cosmology
    //
    const float Om = 0.27;
    const float Ol = 1. - Om;
    Cosmology lcdm(Om, Ol);
    const float h = 0.742;
    const bool isComoving = true;  // for calculating inverse Sigma_crit


    //
    // setup kCorrect information
    //
    ZebraKCorrectSet kcorr(kcorrectDir, templatecount);
    double kcorr_maxz = kcorr.getMaximumRedshift();

    cerr << "kcorrect dir ...... " << kcorrectDir << endl;
    cerr << "   template count . " << templatecount << endl;
    cerr << "   maximum z ...... " << kcorr_maxz << endl;


    //
    // iterate over all lens-src pair 
    //

    // DEBUG
    int lostz = 0;
    int lostlum = 0;
    int lostmorph = 0;
    int lostet = 0;
    int lostsrc = 0;
    int lostsigcrit = 0;
    int lostbgcount = 0;

    photCatObjectList::iterator it = lenslist.begin();
    for (; it != lenslist.end(); ++it) {

      photCatObject* lensobj = *it;

      double lensra = lensobj->getRA();   // in degrees
      double lensdec = lensobj->getDec(); // in degrees
      double lra = lensra * DEGREE;   // in radians
      double ldec = lensdec * DEGREE; // in radians

      /// calculate binning index (on lens characteristics)
      double zlens = lensobj->getPhotoZ();       // use PhotoZ as zlens
      int imorph = morphologybin.getIndex(lensobj->getSEDType()); 
      int iz = redshiftbin.getIndex(zlens, imorph);
      int ilum = luminositybin.getIndex(lensobj->getAbsoluteMagnitude(ZebraKCorrect::SDSSr, 
								      ZebraKCorrect::SDSSi,
								      lcdm, h, kcorr, filterz),
					iz, imorph);

      /* // DEBUG
      cerr << " === lens radec: " << lensra << ", " << lensdec << endl;
      cerr << "lum, morph, z bin: " << ilum << ", " << imorph << ", " << iz << "   "
	   << lensobj->getAbsoluteRMagnitude(lcdm, h, kcorr, filterz) << ", " 
	   << lensobj->getFracDeV() << ", "
	   << lensobj->getTrueZ() << ", "
	   << lensobj->getPhotoZ() << endl;
      /* /

      /// discard lens object if not in bin
      if (ilum == -1 || imorph == -1 || iz == -1) {
	//cerr << lensra << " " << lensdec  << "   "
	//     << ilum << " " << imorph << " " << iz << endl;
	if (ilum == -1) 
	  lostlum++;
	if (imorph == -1) 
	  lostmorph++;
	if (iz == -1) 
	  lostz++;
	continue;
      }
      
      /* // (stellar mass bin)
      double lensmags[5], lensmagerrs[5];
      lensobj->setSloanMags(lensmags, lensmagerrs);
      double this_stellarmass = calculateStellarMass(lensmags, lensmagerrs);
      int istmass = stellarmassbin.getIndex(this_stellarmass);
      // DEBUG
      //cerr << "stellar mass: " << this_stellarmass << endl;
      /* /

      /// calculate angular separation from this lens
      double angdiamdist = lcdm.DA(zlens) * HubbleLengthMpc;  // angular diameter distance [Mpc/h]
      double angcomvdiamdist = angdiamdist * (1.+zlens);      // corresponding comv distance [Mpc/h]
      /* // DEBUG
      cerr << "z: " << zlens << ", "
           << "DA: " << angdiamdist << " [Mpc/h]" << endl;
      /* /

      /// calculate radial bin tolerances for this lens
      vector<double> angularbin(nradbin + 1);
      // DEBUG radialbin
      //cerr << "angular bins" << endl;
      for (int i = 0; i <= nradbin; ++i) {
	angularbin[i] = radialbin[i] / angcomvdiamdist * 180. / PI;
	// DEBUG angularbin
	//cerr << angularbin[i] << " ";  // DEBUG
      }
      //cerr << endl;  // DEBUG

      /// find all possible matching sources (get their indicies of srcvector) in radial bins
      vector<multimap<double, int> > bglist(nradbin);
      for (int irad = 0; irad < nradbin; ++irad) {
	bglist[irad] = srcmesh.getNearAngleMap(lensra, lensdec, 0., 
					       angularbin[irad+1], angularbin[irad]);
	// DEBUG radialbin
	//if (bglist[irad].size() > 0)
	//  cerr << irad << "th radial bin, bg obj count: " << bglist[irad].size() << endl;  
	//cerr << "limits are " << angularbin[irad] << " to " << angularbin[irad+1] << endl;
      }

      //
      // LOOP over sources for this lens object
      //

      int bgcount = 0;  // source object count for this lens object
      multimap<double, int>::const_iterator isrc;
      vector<GGLensSums> tempobjsum(nradbin, GGLensSums(lumlabel));  // temporary data storage
      for (int irad = 0; irad < nradbin; ++irad) {
	for (isrc = bglist[irad].begin(); isrc != bglist[irad].end(); ++isrc) {
	  
	  /// calculate tangential/skew shear
	  bgCatObject* srcobj = srcvector[isrc->second];
	  double sra = srcobj->getRA() * DEGREE;   // in radians
	  double sdec = srcobj->getDec() * DEGREE; // in radians
	  double dra = sra - lra;
	  /////////////  CHECK  (tangent angle)
	  double theta = atan2(cos(ldec)*sin(sdec)-sin(ldec)*cos(sdec)*cos(dra), 
			       cos(sdec)*sin(dra));   // orientation of src from lens
	  double ct = cos(theta);
	  double st = sin(theta);
	  double c2t = ct*ct-st*st;
	  double s2t = 2.*ct*st;
	  
	  double et = -c2t * srcobj->getE1() - s2t * srcobj->getE2();
	  double es =  s2t * srcobj->getE1() - c2t * srcobj->getE2();
	  
	  /* // DEBUG
	  cerr << bgcount << "   "
	       << fixed << setprecision(6) 
	       << lensra << " " << lensdec << "   "
	       << srcobj->getRA() << " " << srcobj->getDec() << "  "
	       << srcobj->getE1() << " " << srcobj->getE2() << "  "
	       << et << " " << es << "   "
	       << theta / DEGREE << " " 
	       << endl;
	  /* /
	  
	  if (isnan(et)) {
	    lostet++;
	    continue;
	  }
	  if (srcobj->getERms() <= 0. || srcobj->getShapeError() <= 0.) {
	    lostsrc++;
	    continue;
	  }
	  
	  /// calculate invsigcrit
	  double ztol = 0.;
	  double invsigcrit = SigmaCritPrefactor * 
	    getComovingSigmaCritInverse(lcdm, zlens, srcobj->getPhotoZ(), ztol);

	  if (invsigcrit == 0 || isnan(invsigcrit)) {
	    lostsigcrit++;
	    continue;
	  }
	  /// apply correction factor for Sigma_crit inverse (future?)
	  //invsigcrit *= 1.; 
	  
	  /// calculate weight
	  double geomweight = invsigcrit * invsigcrit;
	  double sigcritsq = 1. / geomweight;
	  double vare = srcobj->getShapeError() * srcobj->getShapeError(); 
	  /////////////  CHECK  (need to calculate shape noise properly)
	  double varSN = srcobj->getERms() * srcobj->getERms(); 
	  double invshapeweight = (vare + varSN);
	  double lssweight = 1.;
	  double weight = geomweight / invshapeweight * lssweight; 
	  

	  /// calculate other quantities
	  double k0 = varSN*vare/(varSN+vare);
	  double k1 = varSN/(varSN+vare);
	  k1 *= k1;
	  double weightedsignal_t = et * invsigcrit / invshapeweight; // ok
	  double weightedsignal_s = es * invsigcrit / invshapeweight; // ok
	  /////////////  CHECK  (where does this responsivity come from?)
	  double responsiv = weight * (1. - k0 - k1*et*et);           // ok
	  double weightederror_t = weight / invshapeweight * et * et; // ?
	  double weightederror_s = weight / invshapeweight * es * es; // ?
	  double weightedinvsigcrit = weight * invsigcrit;
	  double weightedlum = weight * lensobj->getAbsoluteMagnitude(ZebraKCorrect::SDSSr, 
								      ZebraKCorrect::SDSSi,
								      lcdm, h, kcorr, filterz);
	  double weightedz = weight * zlens;
	  double weightedzwidth = 0.;
	  double weightedrmag = weight * lensobj->getRMag();
	  
	  /// add source object to sm bin 
	  tempobjsum[irad].addPairCounts();
	  tempobjsum[irad].addWeight(weight);
	  tempobjsum[irad].addResponsivity(responsiv);
	  tempobjsum[irad].addDeltaSigma_t(weightedsignal_t);
	  tempobjsum[irad].addDeltaSigma_s(weightedsignal_s);
	  tempobjsum[irad].addError_t(weightederror_t);
	  tempobjsum[irad].addError_s(weightederror_s);
	  tempobjsum[irad].addInvSigmaCrit(weightedinvsigcrit);
	  tempobjsum[irad].addLum(weightedlum);                
	  tempobjsum[irad].addZ(weightedz);                
	  tempobjsum[irad].addZWidth(weightedzwidth);                
	  tempobjsum[irad].addApparentMag(weightedrmag);                
	}
	
	/// count available source objects for this lens
	bgcount += tempobjsum[irad].getPairCounts();

      } // [end of first source 'for' loop]
      
      /// make sure there are enough bg counts for this object
      if (bgcount < 10) {
	lostbgcount++;
	continue; 
      }

      // count lens by luminosity and morphology
      lenscount[ilum][imorph][iz].addLensCounts();
      
      /// include tempobjsum into the proper GGLensSums bin (redshift, lum...)
      for (int irad = 0; irad < nradbin; ++irad) {
	tally[ilum][imorph][iz][irad].addPairCounts(tempobjsum[irad].getPairCounts());
	tally[ilum][imorph][iz][irad].addWeight(tempobjsum[irad].getWeights());
	tally[ilum][imorph][iz][irad].addResponsivity(tempobjsum[irad].getResponsivity());
	tally[ilum][imorph][iz][irad].addDeltaSigma_t(tempobjsum[irad].getDeltaSigma_t());	
	tally[ilum][imorph][iz][irad].addDeltaSigma_s(tempobjsum[irad].getDeltaSigma_s());	
	tally[ilum][imorph][iz][irad].addError_t(tempobjsum[irad].getError_t());
	tally[ilum][imorph][iz][irad].addError_s(tempobjsum[irad].getError_s());
	tally[ilum][imorph][iz][irad].addInvSigmaCrit(tempobjsum[irad].getInvSigmaCrit());
	tally[ilum][imorph][iz][irad].addLum(tempobjsum[irad].getLum());
	tally[ilum][imorph][iz][irad].addZ(tempobjsum[irad].getZ());
	tally[ilum][imorph][iz][irad].addZWidth(tempobjsum[irad].getZWidth());
	tally[ilum][imorph][iz][irad].addApparentMag(tempobjsum[irad].getApparentMag());
      }      

    } // [end of lens objects for loop]

    //
    // print output
    //

    /// print lens count info
    int totallenscount = 0;
    string lenscountfname;
    lenscountfname = lenscountfprefix + lenstype + "." + patchstr + suffix;
    ofstream lenscountf(lenscountfname.c_str());
    if (!lenscountf)
      throw MyException("cannot open "+lenscountfname);

    //lenscount[ilum][imorph][iz]
    lenscountf << "#patch  type  zbin  ";
    for (int ilum = 0; ilum < nlumbin; ++ilum) {
      lenscountf << setw(6) << "lum" << ilum;
    }
    lenscountf << endl;
    for (int iz = 0; iz < nzbin; ++iz) {
      for (int imorph = 0; imorph < nmorphbin; ++imorph) {
	lenscountf << patchstr << "  " 
		   << "type " << imorph << "  "
		   << "z " << iz << "  ";
	for (int ilum = 0; ilum < nlumbin; ++ilum) {
	  lenscountf << setw(6) << lenscount[ilum][imorph][iz].getLensCounts() << " ";
	  totallenscount += lenscount[ilum][imorph][iz].getLensCounts();
	}
	lenscountf << endl;
      }
    }

    cerr << "used lens count ... " << totallenscount << endl;

    /* // DEBUG
    cerr << "lost to z: " << lostz << endl;
    cerr << "lost to lum: " << lostlum << endl;
    cerr << "lost to morph: " << lostmorph << endl;
    cerr << "lost to et: " << lostet << endl;
    cerr << "lost to source obj: " << lostsrc << endl;
    cerr << "lost to sigcrit: " << lostsigcrit << endl;
    cerr << "lost to bgcount: " << lostbgcount << endl;
    //* /

    /// print the sm bin results in the sm file
    for (int ilum = 0; ilum < nlumbin; ++ilum) {
      for (int imorph = 0; imorph < nmorphbin; ++imorph) {
	for (int iz = 0; iz < nzbin; ++iz) {
	  
	  /// define sm bin output file
	  ostringstream outfilename;
	  outfilename << outfprefix + lenstype + "." + patchstr + "." + lumlabel;
	  outfilename << setw(2) << setfill('0') << ilum << ".";
	  outfilename << morphlabel << imorph << ".";
	  outfilename << zlabel << iz << suffix;
	  ofstream outf(outfilename.str().c_str());
	  if (!outf) 
	    throw MyException("cannot open "+outfilename.str());
	  
	  /// output to stream
	  outf << "# patch " << patchstr << endl;
	  /// tally[lum][morph][z][rad]
	  for (int irad = 0; irad  < nradbin; ++irad ) {
	    outf << irad  << " "                                            //  0 (r index)
		 << tally[ilum][imorph][iz][irad].getPairCounts() << " "    //  1
		 << tally[ilum][imorph][iz][irad].getWeights() << " "       //  2
		 << tally[ilum][imorph][iz][irad].getResponsivity() << " "  //  3
		 << tally[ilum][imorph][iz][irad].getDeltaSigma_t() << " "  //  4
		 << tally[ilum][imorph][iz][irad].getDeltaSigma_s() << " "  //  5
		 << tally[ilum][imorph][iz][irad].getError_t() << " "       //  6
		 << tally[ilum][imorph][iz][irad].getError_s() << " "       //  7
		 << tally[ilum][imorph][iz][irad].getInvSigmaCrit() << " "  //  8
		 << tally[ilum][imorph][iz][irad].getLum() << " "           //  9
		 << tally[ilum][imorph][iz][irad].getZ() << " "             // 10
		 << tally[ilum][imorph][iz][irad].getZWidth() << " "        // 11
		 << tally[ilum][imorph][iz][irad].getApparentMag() << " "   // 12
		 << tally[ilum][imorph][iz][irad].getStellarMass() << " "   // 13
		 << endl;
	  }
	}
      }
    }

    // if there are other bins, add print to here
    
    // print config file
    ofstream configf(configfname.c_str());
    if (!configf)
      throw MyException("cannot open file "+configfname);
    configf << "radial" << endl;
    for (int i = 0; i < nradbin; ++i)
      configf << i << " " << radialbin[i] << " " << radialbin[i+1] << " " << endl;
    configf << lumlabel << endl;
    for (int i = 0; i < nlumbin; ++i)
      configf << i << " " << luminositybin[i] << " " << luminositybin[i+1] << " " << endl;
    configf << morphlabel << endl;
    for (int i = 0; i < nmorphbin; ++i)
      configf << i << " " << morphologybin[i] << " " << morphologybin[i+1] << " " << endl;
    configf << zlabel << endl;
    for (int i = 0; i < nzbin; ++i)
      configf << i << " " << redshiftbin[i] << " " << redshiftbin[i+1] << " " << endl;
	  */

  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
