//
// GGLensData.cpp
//
#include "GGLensData.h"

void
GGLensData::getValuesAt(vector<float> radial_vals,
			vector<float>& signalT, vector<float>& signalX, vector<float>& var) {

  int n = radial_vals.size();
  signalX.resize(n, 0.);
  signalT.resize(n, 0.);
  var.resize(n, 0.);

  for (int i=0; i<n; ++i) {

    // when interpolation not possible (inside the smallest radius)
    if (radial_vals[i] < meanR[0]) {
      // approximate with no bias (signal) and closest value available (var)
      signalT[i] = 0.;
      signalX[i] = 0.;
      var[i] = sigma[0] * sigma[0];
      continue;
    }

    // when interpolation not possible (outside the largest radius)
    else if (radial_vals[i] >= meanR[meanR.size()-1]) {
      // return NaNs to signal that subtraction was unsuccessful
      signalT[i] = 0.0/0.0;
      signalX[i] = 0.0/0.0;
      var[i] = 0.0/0.0;
      continue;
    }

    // interpolation (linear) on the signal
    else {
      for (int j=0; j<meanR.size()-1; ++j) {
	if (radial_vals[i] >= meanR[j]) {
	  float dR = meanR[j+1] - meanR[j];
	  signalT[i] = gamT[j] + dR * (gamT[j+1] - gamT[j]);
	  signalX[i] = gamX[j] + dR * (gamX[j+1] - gamX[j]);
	  var[i] = sigma[j] * sigma[j];  // variance by nearest value, no interpolation
	  break;
	}
      }
    }

  }

  return;
}
