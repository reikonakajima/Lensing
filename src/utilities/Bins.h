#ifndef BINS_H
#define BINS_H
#include <iostream>
#include <vector>
using std::vector;
#include "Std.h"


//
// Bins exceptions
//
class BinsError : public MyException {
 public:
  BinsError(const string& m="") : MyException("BinsError: "+m) {}
};


//
// Bins, base class
//
class MultipleBins;

class GenericBins {
 public:
  GenericBins() : reverseOrder(false) {}
  int getIndex(double val) { return findIndex(val); }
  int binSize() { return nbin; }
  int vectorSize() { return binEdges.size(); }
  double getMin() { return min; }
  double getMax() { return max; }
  void setReverse(bool b) { reverseOrder = b; return; }
  bool isReverse() { return reverseOrder; }
  double operator[](int index) { return binEdges[index]; }
  double operator/=(double _scale) {
    for (int i=0; i<this->vectorSize(); ++i) {
      binEdges[i] /= _scale;
    }
  }
  void rescale(double _scale) {
    for (int i=0; i<this->vectorSize(); ++i) {
      binEdges[i] *= _scale;
    }
    return;
  }
  int trim_high(double _max_val);
  vector<double> getBinEdges() const { return binEdges; }
  vector<float> getCentralValues() {
    vector<float> central_values(nbin, 0);
    for (int i=0; i<nbin; ++i) {
      central_values[i] = (binEdges[i] + binEdges[i+1]) / 2.0;
    }
    return central_values;
  }
 protected:
  vector<double> binEdges;
  int findIndex(double val);
  double min;
  double max;
  double nbin;
  bool   reverseOrder;

  friend class MultipleBins;
};


//
// Log bins
//
class LogarithmicBins : public GenericBins {
 public:
  LogarithmicBins(double min_, double max_, int nbin_);
  double getLogStepSize() { return dlnx; }
  int getIndex(double val);
 private:
  double dlnx;
};


//
// bins with arbitrary width
//
class ArbitraryWidthBins : public GenericBins {
 public:
  ArbitraryWidthBins() {}  
  ArbitraryWidthBins(vector<double> binvec_);
  ArbitraryWidthBins(istream& is);
  void initArbBins(vector<double> binvec_, double& tempmin, double& tempmax);
};


//
// Different bins for different indicies
//
class MultipleBins : public ArbitraryWidthBins {  // base class will contain the nominal binning
 public:
  MultipleBins(istream& is, const int n1, const int n2=1);
  int getIndex(const double val, const int i1, const int i2=0);
 private:
  int nindex1;
  int nindex2;
  vector<ArbitraryWidthBins> actualbins;
  int getN1N2index(const int i1, const int i2) {return (i1 * nindex2) + i2;}
};



#endif // BINS_H
