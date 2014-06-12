//
// GGLens.cpp
//
#include "GGLens.h"
#include "Mesh.h"
#include "Bounds.h"
using std::fixed;
using std::setprecision;

const double DEFAULT_MESH_FRAC = 100.;
const int MINIMUM_SOURCE_COUNT_PER_LENS = 10;


GGLensObject::GGLensObject(const int num_radial_bins) :
  tangential_shears(num_radial_bins) {} // initialize the vector object


GGLensObjectList::GGLensObjectList(LensObjectList lens_list,
				   LensObjectList source_list,   // FIXME!!  with sourceObjectList
				   GenericBins _radial_bin,
				   geometry _geom,
				   double _mesh_size) :
  radial_bin(_radial_bin) , geom(_geom), mesh_size(_mesh_size) {

  //
  // create source mesh list, for fast iteration over lens-source pairs
  //

  if (mesh_size == 0.) {
    double mesh_frac = DEFAULT_MESH_FRAC;   // TODO: FIXME!!  calculate
    mesh_size = 30.;
  }

  Bounds<double> srcbounds = source_list.getBounds();
  vector<LensObject*> srcvector = source_list.getVectorForm();   // TODO:  FIXME !!!  vector type
  double ramin = srcbounds.getXMin();
  double ramax = srcbounds.getXMax();
  double decmin = srcbounds.getYMin();
  double decmax = srcbounds.getYMax();
  const int meshdimX = static_cast<int>((ramax-ramin) / mesh_size);
  const int meshdimY = static_cast<int>((decmax-decmin) / mesh_size);
  const bool isPeriodic = false;          // we want a non-periodic mesh
  Mesh<LensObject*> srcmesh(meshdimX, meshdimY, 1, srcvector, isPeriodic,
			    ramin, ramax, decmin, decmax);


  //
  // iterate over all lens; generate lensing signal profile for each lens object
  //

  /// each GGLensObject will have some number of radial bins
  int rad_nbin = radial_bin.size() - 1;  // radial_bin contains bin edges, so nbin is one less

  LensObjectList::iterator it = lens_list.begin();
  for (; it != lens_list.end(); ++it) {

    LensObject* lensobj = *it;

    double lensra = lensobj->getRA();    // FIXME: if SphericalSurface, convert to radians
    double lensdec = lensobj->getDec();

    // DEBUG
    cerr << " === lens radec: " << lensra << ", " << lensdec << endl;


    //
    // find all possible matching sources (get their indicies of srcvector) in radial bins
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
      if (bglist[irad].size() > 0)
        cerr << irad << "th radial bin, bg obj count: " << bglist[irad].size() << endl;
      cerr << "limits are " << radial_bin[irad] << " to " << radial_bin[irad+1] << endl;
    }

    //
    // LOOP over sources for this lens object, accumulate data
    //

    GGLensObject* this_gglens = new GGLensObject(rad_nbin);

    int bg_count = 0;  // source object count for this lens object
    int bad_et = 0;
    int bad_src = 0;
    int lost_bgcount = 0;  // not enough BG count

    multimap<double, int>::const_iterator isrc;
    for (int irad = 0; irad < rad_nbin; ++irad) {
      for (isrc = bglist[irad].begin(); isrc != bglist[irad].end(); ++isrc) {

	//
	// calculate tangential/skew shear
	//
	LensObject* srcobj = srcvector[isrc->second];
	double sra = srcobj->getRA();   // in pixels
	double sdec = srcobj->getDec(); // in pixels
	double dra = sra - lensra;
	double ddec = sdec - lensdec;
	double theta;   // is in RADIANS (output of atan2)
	if (geom == Flat) {   // 2d euclidean
	  theta = atan2(ddec, dra);
	} else if (geom == SphericalSurface) {
	  // FIXME!!  convert everything into RADIANs
	  theta = atan2(cos(lensdec)*sin(sdec)-sin(lensdec)*cos(sdec)*cos(dra),
			cos(sdec)*sin(dra));      // spherical surface
	}

	double ct = cos(theta);
	double st = sin(theta);
	double c2t = ct*ct-st*st;
	double s2t = 2.*ct*st;

	double et = -c2t * srcobj->getE1() - s2t * srcobj->getE2();
	double es =  s2t * srcobj->getE1() - c2t * srcobj->getE2();

	// DEBUG
	cerr << bg_count << "   "
	     << fixed << setprecision(6)
	     << lensra << " " << lensdec << "   "
	     << srcobj->getRA() << " " << srcobj->getDec() << "  "
	     << srcobj->getE1() << " " << srcobj->getE2() << "  "
	     << et << " " << es << "   "
	     << theta / DEGREE << " "
	     << endl;

	if (std::isnan(et)) {
	  bad_et++;
	  continue;
	}
	if (srcobj->getERms() <= 0. || srcobj->getShapeError() <= 0.) {
	  bad_src++;
	  continue;
	}

	//
	// calculate weights
	//
	//// TODO: modularize weights as part of SourceObject
	//// TODO: include invsigcrit
	double vare = srcobj->getShapeError() * srcobj->getShapeError();
	double varSN = srcobj->getERms() * srcobj->getERms();
	double invshapeweight = (vare + varSN);
	double weight = 1 / invshapeweight;

	//
	// calculate other quantities
	//
	double k0 = varSN*vare/(varSN+vare);
	double k1 = varSN/(varSN+vare);
	k1 *= k1;
	double weightedsignal_t = et / invshapeweight; // ok  * invsigcrit
	double weightedsignal_s = es / invshapeweight; // ok  * invsigcrit
	/////////////  CHECK  (where does this responsivity come from?)
	double responsiv = weight * (1. - k0 - k1*et*et);           // ok
	double weightederror_t = weight / invshapeweight * et * et; // ?
	double weightederror_s = weight / invshapeweight * es * es; // ?
	double weightedinvsigcrit = weight;  // * invsigcrit

	/// add source object to sm bin 
	(*this_gglens)(irad).addPairCounts();
	(*this_gglens)(irad).addWeight(weight);
	(*this_gglens)(irad).addResponsivity(responsiv);
	(*this_gglens)(irad).addDeltaSigma_t(weightedsignal_t);
	(*this_gglens)(irad).addDeltaSigma_s(weightedsignal_s);
	(*this_gglens)(irad).addError_t(weightederror_t);
	(*this_gglens)(irad).addError_s(weightederror_s);

      } // END: source object (within bglist[irad]) 'for' loop

	/// count available source objects for this lens
      bg_count += (*this_gglens)(irad).getPairCounts();

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
  cerr << "The size of this GGLensObjectList is: " << this->size() << endl;
}
