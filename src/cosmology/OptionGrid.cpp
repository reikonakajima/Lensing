// Evaluate the PSF size and effective noise level for all SNAP options 
// in 1 filter.
// $Id: OptionGrid.cpp,v 1.2 2007-09-14 01:27:53 garyb Exp $
#include "Psf.h"
#include "Passbands.h"
#include "Solve.h"
#include "StringStuff.h"

#include <list>
using namespace std;
using namespace stringstuff;

class Pair {
public:
  Pair(double x_, double y_): x(x_), y(y_) {}
  double x, y;
  bool operator<(const Pair& rhs) const {return x<rhs.x;}
};


const double pupilObscuration = 0.4;

vector<int> findKeys(string possibleKeys, string keyList) {
  vector<int> idxlist;
  if (possibleKeys.find_first_of(keyList) == string::npos) {
    for (int i=0; i<possibleKeys.size(); i++)
      idxlist.push_back(i);
  } else {
    string::size_type idx=0;
    while ( (idx=possibleKeys.find_first_of(keyList,idx)) != string::npos) {
	idxlist.push_back(idx);
	idx++;
	if (idx >= possibleKeys.size()) break;
    }
  }
  return idxlist;
}

double EEradius(double lOverD, double sigma, double pixsize, 
		double frac=0.5)
{
  //**  double dx=0.006;
  double dx=0.01;

  PsfAiry airy(pupilObscuration, 1./lOverD);
  PsfGaussian gauss(sigma);
  PsfSquare box(pixsize);

  Psf* p1=airy.convolve(gauss);
  Psf* p2=p1->convolve(box);

  double snom=sqrt(p2->AreaSN()/4/PI);

  kTable* kt=p2->getKTable(dx);

  // Some x-space statistics:
  xTable* xt=kt->Transform();
  int N=xt->getN();
  dx = xt->getDx();
  list<Pair> lp;
  double flux=0.;
  for (int i=-N/2; i<N/2; i++)
    for (int j=-N/2; j<N/2; j++) {
      double rad=sqrt(dx*dx*i*i+dx*dx*j*j);
      double val=xt->xval(i,j);
      lp.push_back(Pair(rad,val));
      flux += val;
    }

  // Done with everything now, clean up:
  delete xt; xt=0;
  delete kt; kt=0;
  delete p2; p2=0;
  delete p1; p1=0;

  /*cerr << " N " << N
	   << " dx " << dx
	   << endl;*/
  lp.sort();
  double fsum=0.;
  double eer;
  for (list<Pair>::iterator i=lp.begin();
       i != lp.end();
       ++i) {
    double val = i->y;
    fsum += val;
    if (fsum >= frac*flux)
      return i->x;
  }
}
 
class GeminiProtectedSilver4: public Function1d<> {
  double argMin() const {return 0.32;}
  double argMax() const {return 2.0;}
  double operator()(const double lam) const {
    const double H=0.99;
    const double K=0.009;
    const double L=0.319;
    return pow(H,4.)*exp(-4.*K/(lam-L));
  }
};

int
main(int argc, char *argv[])
{
  vector<double> primaryDiam;
  string primaryDiamKey;
  primaryDiam.push_back(2.0); primaryDiamKey.push_back('A');
  primaryDiam.push_back(1.8); primaryDiamKey.push_back('B');
  primaryDiam.push_back(1.6); primaryDiamKey.push_back('C');
  primaryDiam.push_back(1.4); primaryDiamKey.push_back('a');
  primaryDiam.push_back(1.2); primaryDiamKey.push_back('b');

  vector<double> opticalBlur;
  string opticalBlurKey;
  opticalBlur.push_back(0.020); opticalBlurKey.push_back('E');
  opticalBlur.push_back(0.040); opticalBlurKey.push_back('F');
  opticalBlur.push_back(0.060); opticalBlurKey.push_back('G');
  opticalBlur.push_back(0.090); opticalBlurKey.push_back('e');
  opticalBlur.push_back(0.120); opticalBlurKey.push_back('f');
  opticalBlur.push_back(0.150); opticalBlurKey.push_back('g');

  vector<double> chargeDiffusion;
  string chargeDiffusionKey;
  chargeDiffusion.push_back(5.0); chargeDiffusionKey.push_back('J');
  chargeDiffusion.push_back(4.0); chargeDiffusionKey.push_back('K');
  chargeDiffusion.push_back(3.0); chargeDiffusionKey.push_back('L');

  vector<double> focalRatio;
  string focalRatioKey;
  focalRatio.push_back(9.0); focalRatioKey.push_back('N');
  focalRatio.push_back(11.0); focalRatioKey.push_back('P');
  focalRatio.push_back(13.0); focalRatioKey.push_back('Q');

  vector<double> ccdPixel;
  string ccdPixelKey;
  ccdPixel.push_back(9.0); ccdPixelKey.push_back('S');
  ccdPixel.push_back(10.5); ccdPixelKey.push_back('T');
  ccdPixel.push_back(12.0); ccdPixelKey.push_back('U');

  vector<double> exposureTime;
  string exposureTimeKey;
  exposureTime.push_back(450.); exposureTimeKey.push_back('W');
  exposureTime.push_back(300.); exposureTimeKey.push_back('X');
  exposureTime.push_back(150.); exposureTimeKey.push_back('Y');
  exposureTime.push_back(100.); exposureTimeKey.push_back('w');
  exposureTime.push_back(50.); exposureTimeKey.push_back('x');
  exposureTime.push_back(25.); exposureTimeKey.push_back('y');

  // Other relevant observatory characteristics
  const double exposureOverhead = 30.;	// Read/repoint time between exp.
  const double nirPixel = 18.0;
  const double ccdReadNoise = 4;
  const double nirReadNoise = 9.; // from Bebek's TOTAL 12e in 300s
  const double ccdDarkRate = 50./3600./(10.5*10.5); //e per s per sq micron
  const double nirDarkRate = 0.2/(18.0*18.0); // a conservative spec
  const double filterEfficiency = 0.9;
  const int ccdDitherSteps=4;

  // Focal palane area per NIR filter
  const double nirArea = 12.*pow(2048*nirPixel, 2.); 

  //////////////////////////////////////////////////////////

  // Now get all necessary passband info
  try {
    if (argc<3) {
      cerr << "Usage: OptionGrid <filter number> <filter filename> [keylist]"
	   << endl;
      exit(1);
    }
    string filtname = argv[2];
    Passband pb("@"+filtname);

    string spec;

    string useKeys;
    if(argc>3) useKeys=argv[3];
    vector<int> usePrimaryDiam = findKeys(primaryDiamKey, useKeys);
    vector<int> useOpticalBlur = findKeys(opticalBlurKey, useKeys);
    vector<int> useChargeDiffusion = findKeys(chargeDiffusionKey, useKeys);
    vector<int> useFocalRatio = findKeys(focalRatioKey, useKeys);
    vector<int> useCCDPixel = findKeys(ccdPixelKey, useKeys);
    vector<int> useExposureTime = findKeys(exposureTimeKey, useKeys);

    // Which detector are we using?
    const bool useCCD = pb.Centroid() < 0.96;

    string qename = useCCD ? "ccdQE.txt" : "hgcdteQE.txt";
    Passband qe("@"+qename);
    pb *= qe;

    Passband telescope(new GeminiProtectedSilver4);
    pb *= telescope;

    double lambda = pb.Centroid();

    const double baseMag=25.;
    Constant abBase(ABMagToFlux(baseMag));
    /*    ifstream ifs("zodi.txt");
	  Table<> zodi(ifs);
	  /*** SNAPSim zodi table does not start until 0.4 um */
    ZodiAtPole zodi;


    cout << "# " << stringstuff::taggedCommandLine(argc, argv) << endl;
    cout << "# Code   EE50    mn     m1   skyRate  fsky/yr " << endl;

    // Begin the nested loop series, put things that influence ePSF on the
    // outside.

    bool newPSF = true;
    double ee50 = 0.;
    const double fraction=0.50;

    for (int i0=0; 
	 i0<usePrimaryDiam.size(); 
	 i0++) {
      int iPrimaryDiam = usePrimaryDiam[i0];

      newPSF = true;
      const double D = primaryDiam[iPrimaryDiam];
      const double lOverD = lambda*1e-6 / D / ARCSEC;

      const double area = PI * (D*D/4.) * (1.-pow(pupilObscuration,2.))
	* filterEfficiency;

      for (int i1=0; 
	   i1<useFocalRatio.size(); 
	   i1++) {
	int iFocalRatio=useFocalRatio[i1];
	newPSF = true;
	double focalLength = D * focalRatio[iFocalRatio];

	for (int i2=0; 
	     i2<useOpticalBlur.size(); 
	     i2++) {
	  int iOpticalBlur=useOpticalBlur[i2];
	  newPSF = true;

	  for (int i3=0; 
	       i3<useChargeDiffusion.size(); 
	       i3++) {
	    int iChargeDiffusion=useChargeDiffusion[i3];
	    double diffusion=0.;
	    if (useCCD) {
	      newPSF = true;
	      diffusion = chargeDiffusion[iChargeDiffusion];
	    }

	    double blur = hypot(opticalBlur[iOpticalBlur],
				diffusion*1e-6/focalLength/ARCSEC);

	    for (int i4=0; 
		 i4<useCCDPixel.size(); 
		 i4++) {
	      int iCcdPixel=useCCDPixel[i4];
	      double pixel = useCCD ? ccdPixel[iCcdPixel] : nirPixel;
	      if (useCCD) newPSF = true;
	      pixel *= 1e-6/focalLength/ARCSEC;

	      // Get the ePSF size; all args in ARCSEC
	      if (newPSF) {
		ee50 = EEradius(lOverD, blur, pixel, fraction);
		newPSF = false;
	      }

	      // Get source and sky rates
	      double baseRate = pb.Integrate(abBase)*area;
	      // Sky per pixel per second:
	      /** My zodi function already in arcsec^-2 **
		  double pixArea = pow(pixel*ARCSEC, 2.); */
	      double pixArea = pow(pixel, 2.);
	      double skyRate = pb.Integrate(zodi)*area*pixArea;

	      /** cerr << "ee50: " << ee50
			<< " 25th mag e/s: " << baseRate
			<< " sky e/s/pix: " << skyRate
			<< endl;   //***/

	      for (int i5=0; 
		   i5<useExposureTime.size(); 
		   i5++) {
		int iExposureTime=useExposureTime[i5];
    		double exptime = exposureTime[iExposureTime];

		double fskyPerYear = nirArea*pow(1e-6/focalLength, 2.) 
		  * YEAR / ((exptime+exposureOverhead)*SECOND)
		  / (ccdDitherSteps * 4 * PI);
		// ???? CHANGE FOR 8-FILTER SCHEME:
		int nRead = (useCCD ? 1 : 2) * ccdDitherSteps;
		fskyPerYear *= 0.5; // ??? since 1/2 detector area per vis

		string simCode;
		{
		  ostringstream oss;
		  oss << primaryDiamKey[iPrimaryDiam]
		      << opticalBlurKey[iOpticalBlur]
		      << chargeDiffusionKey[iChargeDiffusion]
		      << focalRatioKey[iFocalRatio]
		      << ccdPixelKey[iCcdPixel]
		      << exposureTimeKey[iExposureTime]
		      << argv[1];
		    simCode = oss.str();
		}

		double m1 = baseMag + 2.5*log10(nRead*exptime*baseRate);
		double darkRate = useCCD ? 
		  ccdDarkRate * pow(ccdPixel[iCcdPixel],2.) : 
		  nirDarkRate * pow(nirPixel,2.);
		double varPerPix = (skyRate + darkRate)*exptime
		  + pow(useCCD ? ccdReadNoise : nirReadNoise , 2.);
		varPerPix *= nRead;
		// Scale to 1 arcsec square and convert to mag:
		double mn = m1 - 1.25*log10(varPerPix / pixel / pixel);

		cout << simCode
		     << " " << setprecision(4) << setw(6) << ee50
		     << " " << setprecision(5) << setw(6) << mn
		     << " " << setprecision(4) << setw(5) << m1
		     << " " << setprecision(3) << setw(6) << skyRate
		     << " " << setprecision(3) << setw(6) << fskyPerYear
		     << endl;
		
	      } // exptime
	    } // ccdpixel
	  } // charge diffusion
	} // optical blur
      } // focal ratio
    } // primary diameter
  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }
}

		
