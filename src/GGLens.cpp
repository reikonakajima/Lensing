//
// GGLens.cpp
//
#include "GGLens.h"
#include "Mesh.h"
#include "Bounds.h"
using std::fixed;
using std::setprecision;

const double DEFAULT_MESH_FRAC = 100.;
const int MINIMUM_SOURCE_COUNT_PER_LENS = 1;


template<class lensObjPtr, class srcObjPtr>
GGLensObjectList<lensObjPtr, srcObjPtr>::GGLensObjectList(LensObjectList<lensObjPtr> lens_list,
							  SourceObjectList<srcObjPtr> source_list,
							  GenericBins _radial_bin,
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
  int rad_nbin = radial_bin.size() - 1;  // radial_bin contains bin edges, so nbin is one less

  typename vector<lensObjPtr>::iterator it = lens_list.begin();
  for (; it != lens_list.end(); ++it) {

    lensObjPtr lensobj = *it;

    double lensra = lensobj->getRA();    // FIXME: if SphericalSurface, convert to radians
    double lensdec = lensobj->getDec();
    if (!srcbounds.includes(Position<double>(lensra, lensdec)))
	continue;

    // DEBUG
    //cerr << " === ";
    //lensobj->printLine(cerr);

    //
    // collect all matching sources (get their indicies of srcvector) in radial bins
    //

    vector<multimap<double, int> > bglist(rad_nbin);
    for (int irad = 0; irad < rad_nbin; ++irad) {
      if (geom == Flat) {                  // FIXME!!  make Mesh class take geometry
	bglist[irad] = srcmesh.getNearMeshMap(lensra, lensdec, 0.,
					      radial_bin[irad+1], radial_bin[irad]);
      } else if (geom == SphericalSurface) {
	bglist[irad] = srcmesh.getNearAngleMap(lensra, lensdec, 0.,
					       radial_bin[irad+1], radial_bin[irad]);
      }
      // DEBUG radial_bin
      if (bglist[irad].size() > 0) {
	  source_vector[bglist[irad].begin()->second]->printLine(cerr);
	  cerr << irad << "th radial bin, bg obj count: " << bglist[irad].size() << endl;
	  cerr << "limits are " << radial_bin[irad] << " to " << radial_bin[irad+1] << endl;
      }
    }

    //
    // LOOP over all sources (by radial bins) for this lens object, accumulate data
    //

    GGLensObject<lensObjPtr> *this_gglens = new GGLensObject<lensObjPtr>(lensobj, rad_nbin);

    int bg_count = 0;      // source object count for this lens object
    int bad_et = 0;
    int bad_src = 0;
    int lost_bgcount = 0;  // lens objects lost through "not enough BG object count"

    multimap<double, int>::const_iterator isrc;
    for (int irad = 0; irad < rad_nbin; ++irad) {
      for (isrc = bglist[irad].begin(); isrc != bglist[irad].end(); ++isrc) {

	//
	// calculate tangential/skew shear
	//
	srcObjPtr srcobj = source_vector[isrc->second];
	double sra = srcobj->getRA();
	double sdec = srcobj->getDec();
	double dra = sra - lensra;
	double ddec = sdec - lensdec;
	double theta;   // is in RADIANS (output of atan2)
	if (geom == Flat) {   // 2d euclidean
	  theta = atan2(ddec, dra);
	} else if (geom == SphericalSurface) {
	  // FIXME!!  convert everything into RADIANs
	  throw GGLensError(" FIXME!!  convert everything into RADIANs");
	  theta = atan2(cos(lensdec)*sin(sdec)-sin(lensdec)*cos(sdec)*cos(dra),
			cos(sdec)*sin(dra));      // spherical surface
	}

	double ct = cos(theta);
	double st = sin(theta);
	double c2t = ct*ct-st*st;
	double s2t = 2.*ct*st;

	double et = -c2t * srcobj->getE1() - s2t * srcobj->getE2();
	double es =  s2t * srcobj->getE1() - c2t * srcobj->getE2();

	/// increment available source objects for this lens
	++bg_count;

	/*/ DEBUG
	cerr << bg_count << "   " << irad << " " << sqrt(isrc->first) << "  "
	     << fixed << setprecision(6)
	     << lensra << " " << lensdec << "   "
	     << fixed << setprecision(3)
	     << srcobj->getId() << " " << srcobj->getSNratio() << " " << srcobj->getWeight() << "  "
	     << fixed << setprecision(6)
	     << srcobj->getRA() << " " << srcobj->getDec() << "  "
	     << srcobj->getE1() << " " << srcobj->getE2() << "  "
	     << et << " " << es << "   "
	     << theta / DEGREE << " "
	     << endl;
        /*/

	if (std::isnan(et)) {
	  bad_et++;
	  continue;
	}
	/*   TODO:  transfer: these will need to be included for SDSS sources...
	if (srcobj->getERms() <= 0. || srcobj->getShapeError() <= 0.) {
	  bad_src++;
	  continue;
	}
	*/

	//
	// calculate weights, responsivities, shears and shear errors
	//
	double weight = srcobj->getWeight();
	double responsiv = srcobj->getResponsivity(et);
	double weightedsignal_t = et * weight;
	double weightedsignal_s = es * weight;
	double weightedVariance_t = weight * weight * et * et;
	double weightedVariance_s = weight * weight * es * es;

	/// add source object to sm bin 
	(*this_gglens)[irad].addPairCounts();
	(*this_gglens)[irad].addWeight(weight);
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

    // DEBUG
    cerr << "breaking (debug)" << endl;
    break;

  }  // END: lens_list loop

  /// At this point, the GGLensObjectList should be full
  cerr << "All lens objects loaded and paired with source objects." << endl;
  cerr << "The size of the input LensObjectList is: " << lens_list.size() << endl;
  cerr << "The size of this GGLensObjectList is: " << this->size() << endl;
}


//
// explicit instantiations
//

#include "RCSLenSObjects.h"
#include "StarMaskObjects.h"
template class GGLensObjectList<StarMaskObject*, RCSLenSObject*>;
