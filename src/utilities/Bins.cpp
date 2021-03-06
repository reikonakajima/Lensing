#include "Bins.h"
#include "StringStuff.h"
using std::istringstream;

LogarithmicBins::LogarithmicBins(double min_, double max_, int nbin_) {
  min = min_;
  max = max_; 
  nbin = nbin_;
  // check params
  if (min < 0) 
    throw BinsError("invalid Log bin minimum");
  if (max <= 0) 
    throw BinsError("invalid Log bin maximum");
  if (min == 0)
    min = max / nbin / 1000.;
  // calculate & assign bins
  dlnx = log(max/min)/nbin;  
  for (int i = 0; i <= nbin; ++i) 
    binEdges.push_back(min*exp(i*dlnx));
}


int
GenericBins::trim_high(double _max_val) {
  int count = 0;
  if (!reverseOrder) {
    for (int i=this->vectorSize()-1; i >= 0; --i) {
      if (binEdges[i] > _max_val) {
	binEdges.erase(binEdges.begin() + i);
	++count;
      }
    }
    nbin = binEdges.size() - 1;
    min = binEdges[0];
    max = binEdges[nbin-1];
  } else {
    while (*(binEdges.begin()) > _max_val) {
      binEdges.erase(binEdges.begin());
      ++count;
    }
    nbin = binEdges.size() - 1;
    max = binEdges[0];
    min = binEdges[nbin-1];
  }
  return count;
}

int
GenericBins::findIndex(double val) { 
  int index = 0;

  if (!reverseOrder) {
    while (val >= binEdges[index] && val < binEdges[binEdges.size()-1])
      ++index;
  }
  else {
    while (val < binEdges[index] && val >= binEdges[binEdges.size()-1])
      ++index;
  }
  
  if (index == nbin + 1)  // because there are (nbin+1) edges
    return -1;
    
  --index;

  return index;
}


int
LogarithmicBins::getIndex(double val) { 
  int index = 0;
  if (val < binEdges[index])
    return -1;
  while (binEdges[index] < val && index < binEdges.size())
    ++index;
    
  return index - 1;
}



ArbitraryWidthBins::ArbitraryWidthBins(vector<double> binvec_) {
  double tmin, tmax;
  initArbBins(binvec_, tmin, tmax);
  min = tmin;
  max = tmax;
  return;
}


ArbitraryWidthBins::ArbitraryWidthBins(istream& is) {
  // read file
  vector<double> binvec;
  double val;
  double tmin, tmax;
  while (is >> val) {
    binvec.push_back(val);
  }
  initArbBins(binvec, tmin, tmax);
  min = tmin;
  max = tmax;
  return;
}


void
ArbitraryWidthBins::initArbBins(vector<double> binvec_, double& tempmin, double& tempmax) {

  // assign vector
  binEdges.clear();
  binEdges.resize(binvec_.size());
  for (int i = 0; i < binvec_.size(); ++i) {
    binEdges[i] = binvec_[i];
  }
  // set bin parameters
  tempmin = binEdges[0];
  nbin = binvec_.size() - 1;  // because binvec_ defines the bin edges
  tempmax = binEdges[nbin];

  reverseOrder = false;
  if (tempmax < tempmin) {
    SWAP(tempmax, tempmin);
    reverseOrder = true;
  }

  return;
}


MultipleBins::MultipleBins(istream& is, int n1, int n2) {

  string buffer;
  vector<double> initvec;
  double val;
  string junk;
  int nindex = 2;
  if (n2 == 1)
    nindex = 1;

  // read nominal bins (the first line)
  if (!getlineNoComment(is, buffer))
    throw BinsError("MultipleBins: input file has insufficient lines");
  istringstream iss(buffer);
  for (int i=0; i<nindex; ++i) {
    if (!(iss >> junk)) {
      cerr << buffer << endl;
      throw BinsError("MultipleBins: first line has incorrect format: "+buffer);
    }
  }
  while (iss >> val) {
    initvec.push_back(val);
  }
  nbin = initvec.size() - 1;
  if (nbin < 1) {
    cerr << nbin << endl;
    throw BinsError("MultipleBins: input first line has no bins: "+buffer);
  }
  double tmin, tmax;
  this->initArbBins(initvec, tmin, tmax);

  // setup multi-bins for the two indicies
  if (n1 < 1 || n2 < 1)
    throw BinsError("wrong number of dimensions specified");
  nindex1 = n1;
  nindex2 = n2;
  actualbins.clear();
  actualbins.resize(nindex1 * nindex2);
  // read file
  int index1, index2=0;
  for (int i_bin = 0; i_bin < nindex1 * nindex2; ++i_bin) {
    if (!getlineNoComment(is, buffer))
      throw BinsError("input file has insufficient number of lines");
    istringstream iss(buffer);
    if (!(iss >> index1)) {
      cerr << buffer << endl;
      throw BinsError("MultipleBins: input file index and specified format differs in line "
		      +buffer);
    }
    if (nindex2 > 1) {  // because the second index column is not in the input file for such cases
      if (!(iss >> index2)) {
	cerr << buffer << endl;
	throw BinsError("MultipleBins: input file index and specified format differs in line "
			+buffer);
      }
    }
    int binindex = getN1N2index(index1, index2);
    // get the bin for this set of indicies (index1, index2);
    initvec.clear();
    while (iss >> val) {
      initvec.push_back(val);
    }
    if (initvec.size() != nbin+1)
      throw BinsError("input file index and specified format differs in line "+buffer);
    actualbins[binindex].initArbBins(initvec, tmin, tmax);
  }
  return;
}


int 
MultipleBins::getIndex(double val, int i1, int i2) {
  
  // first, check the given external indicies i1 and i2
  if (i1 < 0 || i2 < 0)
    return -1;
  if (i1 >= nindex1 || i2 >= nindex2)
    throw BinsError("getIndex: specified external index out of range");
      
  // then check the binning index
  int binindex = MultipleBins::getN1N2index(i1,i2);

  /*/ debug
  for (int i = 0; i < actualbins[binindex].size()+1; ++i)
    cerr << actualbins[binindex][i] << " ";
  cerr << endl;
  /*/
  return actualbins[binindex].findIndex(val);
}
