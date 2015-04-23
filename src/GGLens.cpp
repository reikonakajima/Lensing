//
// GGLens.cpp
//
#include "GGLens.h"
#include "Mesh.h"
#include "Bounds.h"
#include "KiDSObjects.h"
using std::fixed;
using std::setprecision;

const double DEFAULT_MESH_FRAC = 100.;
const int MINIMUM_SOURCE_COUNT_PER_LENS = 1;


template<class lensObjPtr, class srcObjPtr>
GGLensObjectList<lensObjPtr, srcObjPtr>::GGLensObjectList(LensObjectList<lensObjPtr> lens_list,
							  SourceObjectList<srcObjPtr> source_list,
							  GenericBins _radial_bin,
							  bool radialBinInMpc,
							  bool normalizeToSigmaCrit,
							  cosmology::Cosmology cosmo,
							  double h,
							  double min_lens_src_delta_z,
							  geometry _geom,
							  double mesh_frac) :
  radial_bin(_radial_bin) , geom(_geom) {

  //
  // adjust cosmological distance calculations for h
  //
  const double HUBBLE_LENGTH_MPC = HubbleLengthMpc * h;
  const double SIGMA_CRIT_PREFACTOR = SigmaCritPrefactor / h;

  //
  // create source mesh list, for fast iteration over lens-source pairs
  //

  if (mesh_frac == 0.) {
    cerr << "setting mesh_frac to " << DEFAULT_MESH_FRAC << endl;
    mesh_frac = DEFAULT_MESH_FRAC;
  } else if (mesh_frac < 0.) {
    throw GGLensError("mesh_frac set to negative value");
  }

  Bounds<double> srcbounds = source_list.getBounds();
  vector<srcObjPtr> source_vector = source_list.getVectorForm();

  double ramin = srcbounds.getXMin();
  double ramax = srcbounds.getXMax();
  double decmin = srcbounds.getYMin();
  double decmax = srcbounds.getYMax();

  mesh_size = (ramax-ramin) / mesh_frac;
  double temp = (decmax-decmin) / mesh_frac;
  if (mesh_size > temp)
      mesh_size = temp;

  const int meshdimX = static_cast<int>((ramax-ramin) / mesh_size);
  const int meshdimY = static_cast<int>((decmax-decmin) / mesh_size);

  const bool isPeriodic = false;          // we want a non-periodic mesh
  Mesh<srcObjPtr> srcmesh(meshdimX, meshdimY, 1, source_vector, isPeriodic,
			  ramin, ramax, decmin, decmax);

  /// prepare for conversion into DEGREE (still needs to be devided by Angular Diameter Distance)
  if (radialBinInMpc) {
    radial_bin /= DEGREE; // still in Mpc, divided by DEGREE
  } else {
    radial_bin /= 3600.;  // convert from arcsec into degree
  }

  //
  // iterate over all lens; generate lensing signal profile for each lens object
  //

  /// each GGLensObject will have some number of radial bins
  int rad_nbin = radial_bin.binSize();

  // lens objects lost through "not enough BG object count"
  int lost_bgcount = 0;

  typename vector<lensObjPtr>::iterator it = lens_list.begin();
  for (; it != lens_list.end(); ++it) {

    lensObjPtr lensobj = *it;

    double lensra = lensobj->getRA();
    double lensdec = lensobj->getDec();
    double ldec = lensdec * DEGREE;
    double zlens = lensobj->getRedshift();

    /// correct radial binning from Mpc to angular scale (degrees)
    GenericBins angular_radial_bin(radial_bin);  // make a copy
    if (radialBinInMpc) {
      angular_radial_bin /= cosmo.DA(zlens) * HUBBLE_LENGTH_MPC;   // angular_radial_bin in degrees
    }

    if (!srcbounds.includes(Position<double>(lensra, lensdec)))
       continue;

    //
    // collect all matching sources (get their indicies of srcvector) in radial bins
    //
    vector<multimap<double, int> > bglist(rad_nbin);
    for (int irad = 0; irad < rad_nbin; ++irad) {
      if (geom == Flat) {
	throw GGLensError("GGLens should not take Flat geometry");
      } else if (geom == SphericalSurface) {
	bglist[irad] = srcmesh.getNearAngleMap(lensra, lensdec, 0.,
					       angular_radial_bin[irad+1],
					       angular_radial_bin[irad]);
      }
    }

    //
    // LOOP over all sources (by radial bins) for this lens object, accumulate data
    //

    GGLensObject<lensObjPtr> *this_gglens = new GGLensObject<lensObjPtr>(lensobj, rad_nbin);

    int bg_count = 0;      // source object count for this lens object
    int bad_et = 0;
    int bad_src = 0;

    // integrate over redshift probability distribution?
    bool use_pz = true;

    // get the redshift bins for p(z) of source objects
    valarray<float> src_zbins = source_list.getPzBins();

    multimap<double, int>::const_iterator isrc;
    for (int irad = 0; irad < rad_nbin; ++irad) {
      for (isrc = bglist[irad].begin(); isrc != bglist[irad].end(); ++isrc) {

	//
	// calculate tangential/skew shear on Spherical Surface
	//
	srcObjPtr srcobj = source_vector[isrc->second];
	double sra = srcobj->getRA();
	double sdec = srcobj->getDec() * DEGREE;
	double dra = (sra - lensra) * DEGREE;
	double theta;

	if (geom == Flat) {   // 2d euclidean
	  throw GGLensError("GGLens should not take Flat geometry");
	} else if (geom == SphericalSurface) {
	  theta = atan2(cos(ldec)*sin(sdec) - sin(ldec)*cos(sdec)*cos(dra),
			cos(sdec)*sin(dra));      // spherical surface
	}

	double ct = cos(theta);
	double st = sin(theta);
	double c2t = ct*ct-st*st;
	double s2t = 2.*ct*st;

	double et = -c2t * srcobj->getG1() - s2t * srcobj->getG2();
	double es =  s2t * srcobj->getG1() - c2t * srcobj->getG2();

	/// increment available source objects for this lens
	++bg_count;

	if (std::isnan(et)) {
	  bad_et++;
	  continue;  // to the next src obj
	}

	//
	// calculate weights, (responsivities,) shears and shear errors
	// (reject based on redshift)
	//
	double Sigma_crit_inv = 1.; // if we want to stack pure shear signals, keep Sigma_crit=1.
	double obj_weight = srcobj->getWeight();  // weight for the individual src object
	double weight = obj_weight;               // if no geometry/cosmology
	double weight_and_norm = weight;          // if no geometry/cosmology
	double responsiv = 1.;  // responsivity used for e1/e2; irrelevant if using g1/g2
	valarray<float> src_pz; // source redshift PDF

	if (normalizeToSigmaCrit) {
	  // integrate over p(z) to calculate Sigma_crit_inv
	  if (typeid(*srcobj) == typeid(KiDSObject)) {

	    float zsrc = srcobj->getRedshift();  // this is z_B

	    // discard src object that is not far enough behind of the lens
	    if ((zsrc-zlens) < min_lens_src_delta_z)
	      continue;  // to the next src obj

	    if (use_pz) {
	      // first, normalize p(z)
	      src_pz = srcobj->getPz();   // this is p(z), unnormalized
	      double pz_sum = 0;
	      for (int i=0; i<src_pz.size(); ++i) pz_sum += src_pz[i];
	      src_pz /= pz_sum;

	      // integrate over p(z) to get Sigma_crit_inv
	      Sigma_crit_inv = 0;
	      for (int i=0; i<src_pz.size(); ++i) {
		Sigma_crit_inv += src_pz[i] * cosmo.LensShear(src_zbins[i], zlens);
	      }
	    }
	    else {
	      // calculate Sigma_crit using source z_B only
	      Sigma_crit_inv = cosmo.LensShear(zsrc, zlens);
	    }
	    Sigma_crit_inv /= SIGMA_CRIT_PREFACTOR;  // normalize to units of [(h) M_sun pc^-2]

	    // calculate weight with normalization: DeltaSigma = Sum (w*Sigma_crit*gamma_t)
	    weight_and_norm *= Sigma_crit_inv;           // normalization = Sigma_crit
	    weight *= Sigma_crit_inv * Sigma_crit_inv;   // update to include geom_wt=Sigma_crit^-2
	  }
	  else { // FUTURE TODO  // for a class using e1/e2 instead of g1/g2 reduced shear
	    responsiv = srcobj->getResponsivity(et);
	  }
	}  // end: if (normalizeToSigmaCrit)

	// the weighted shears and variance
	double weightedsignal_t = et * weight_and_norm;    // w * w_geom * et * Sigma_crit
	double weightedsignal_s = es * weight_and_norm;    //   = weight_and_norm * et
	double weightedVariance_t = obj_weight * et * et;  // w * w_geom * et^2 * Sigma_crit^2
	double weightedVariance_s = obj_weight * es * es;  //   = w * et^2

	/// add source object to sm bin 
	(*this_gglens)[irad].addPairCounts();
	(*this_gglens)[irad].addWeight(weight);
	(*this_gglens)[irad].addWeightSq(weight*weight);
	(*this_gglens)[irad].addResponsivity(responsiv);
	(*this_gglens)[irad].addDeltaSigma_t(weightedsignal_t);
	(*this_gglens)[irad].addDeltaSigma_s(weightedsignal_s);
	(*this_gglens)[irad].addVariance_t(weightedVariance_t);
	(*this_gglens)[irad].addVariance_s(weightedVariance_s);

      } // END: source object (within bglist[irad]) 'for' loop

    } // END: irad 'for' loop
    
    //
    // reject lenses without enough background source objects
    //
    if (bg_count < MINIMUM_SOURCE_COUNT_PER_LENS) {
      lost_bgcount++;
      delete this_gglens;
      continue;
    }

    //
    // add GGLens object to the list
    //
    gglens_object_list.push_back(this_gglens);

  }  // END: lens_list loop


  /// At this point, the GGLensObjectList should be full
  cerr << "All lens objects loaded and paired with source objects." << endl;
  cerr << "The size of the input LensObjectList is: " << lens_list.size() << endl;
  cerr << "Rejected lenses without enough background source objects: " << lost_bgcount << endl;
  cerr << "The size of this GGLensObjectList is: " << this->size() << endl;
}


template<class lensObjPtr, class srcObjPtr>
vector<GGLensObjectList<lensObjPtr, srcObjPtr> >
GGLensObjectList<lensObjPtr, srcObjPtr>::splitList(const int nSplit) {
  vector<GGLensObjectList<lensObjPtr, srcObjPtr> > split_lists(nSplit);

  // initialize each GGLensObjectList as the same as the parent
  for (int i = 0; i < nSplit; ++i) {
    split_lists[i].geom = this->geom;
    split_lists[i].mesh_size = this->mesh_size;
    split_lists[i].radial_bin = this->radial_bin;
  }
  return split_lists;
}


//
// explicit instantiations
//

#include "RCSLenSObjects.h"
#include "KiDSObjects.h"
#include "StarMaskObjects.h"
#include "GAMAObjects.h"
#include "RandomObjects.h"
template class GGLensObjectList<StarMaskObject*, RCSLenSObject*>;
template class GGLensObjectList<StarMaskObject*, KiDSObject*>;
template class GGLensObjectList<GAMAObject*, KiDSObject*>;
template class GGLensObjectList<RandomObject*, KiDSObject*>;
