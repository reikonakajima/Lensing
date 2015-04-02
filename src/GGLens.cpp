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
							  bool normalizeToSigmaCrit,
							  geometry _geom,
							  double mesh_frac) :
  radial_bin(_radial_bin) , geom(_geom) {

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
					       radial_bin[irad+1], radial_bin[irad]);
      }
    }

    //
    // LOOP over all sources (by radial bins) for this lens object, accumulate data
    //

    GGLensObject<lensObjPtr> *this_gglens = new GGLensObject<lensObjPtr>(lensobj, rad_nbin);

    int bg_count = 0;      // source object count for this lens object
    int bad_et = 0;
    int bad_src = 0;

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
	  continue;
	}

	//
	// calculate weights, responsivities, shears and shear errors
	//
	double Sigma_crit = 1.; // if we want to stack pure shear signals, keep Sigma_crit 1.
	double weight = srcobj->getWeight();  // weight for the individual object
	double responsiv = 1.;  // responsivity for e1/e2; irrelevant if using g1/g2
	valarray<float> src_pz;
	valarray<float> src_zbins;

	// if source is a KiDSObject, then integrate over p(z) to calculate Sigma_crit
	if (normalizeToSigmaCrit) {
	  if (typeid(*srcobj) == typeid(KiDSObject)) {
	    src_pz = srcobj->getPz();
	    src_zbins = source_list.getPzBins();
	    /*/ DEBUG
	    for (int i=0; i<KiDSObjectList::NUM_PZ_ELEM; ++i) {
	      cerr << i << " " << src_zbins[i] << " " << src_pz[i] << endl;
	    }
	    //*/ // end DEBUG
	    //Sigma_crit = cosmo.getSigmaCrit(zlens, src_pz, src_zbins, min_lens_src_sep);
	    //double geom_weight = 1.0/Sigma_crit/Sigma_crit;
	    //double weight *= geom_weight;  // update weight to include geometric weighting
	  }
	  else { // FUTURE TODO  // for a class using e1/e2 instead of g1/g2 reduced shear
	    responsiv = srcobj->getResponsivity(et);
	  }
	}

	// the weighted shears
	double weightedsignal_t = et * weight;
	double weightedsignal_s = es * weight;
	double weightedVariance_t = weight * weight* et * et;
	double weightedVariance_s = weight * weight* es * es;

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
