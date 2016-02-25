// $Id: FisherMatrices.h,v 1.5 2006-08-08 00:27:04 garyb Exp $
// Classes and subroutines used generically for cosmological Fisher Matrix
// analyses.

#ifndef FISHERMATRICES_H
#define FISHERMATRICES_H

#include <sstream>
#include "Std.h"
#include "Matrix.h"

namespace cosmology {

class FisherParameter {
public:
  // These are "types" of parameters: are they associated
  // with a particular redshift, or redshift & harmonic, or global?
  static const int GLOBAL=0;
  static const int ZBIN=1;
  static const int LZBIN=2;
  FisherParameter(int t_=GLOBAL): type(t_), priorSD(0.) {
    id=idCounter++;
  }
  bool operator==(const FisherParameter& rhs) const {return id==rhs.id;}
  virtual ~FisherParameter() {}
  double getPrior() const {return priorSD;}
  void setPrior(double sd_) {priorSD=sd_;}
  int getType() const {return type;}
  void setType(int t_) {type=t_;}
  int getID() const {return id;}
  virtual string name() const {return "Undefined";}
private:
  FisherParameter(const FisherParameter& rhs) {} // hide copy
  void operator=(const FisherParameter& rhs) {} // and assignment
  double priorSD;  //0 will signify no prior info
  int type;
  static int idCounter;
  int id;
};

class ParameterVector: public vector<FisherParameter*> {
public:
  ParameterVector(int size=0): vector<FisherParameter*>(size) {}
  int findID(int id) const {
    for (int i=0; i<this->size(); i++)
      if ((*this)[i]->getID()==id) return i;
    return -1;
  }
  void deleteMembers() {
    for (int i=0; i<this->size(); i++) {
      delete (*this)[i]; (*this)[i]=0;
    }
  }
  // Add into this vector elements of another that we don't yet have:
  void merge(const ParameterVector& rhs) {
    for (int i=0; i<rhs.size(); i++)
      if (findID(rhs[i]->getID())<0) push_back(rhs[i]);
  }
};

// Manipulations of Fisher matrices: 

extern SqDMatrix
TransformFisher(const SqDMatrix& Fin, const ParameterVector& pvin,
		const DMatrix& project, const ParameterVector& pvout);

extern void 
MarginalizeOver(const SqDMatrix& Fin, const ParameterVector& pvin,
		SqDMatrix& Fout, ParameterVector& pvout,
		vector<int> ids);

// Remove selected parameters from Fisher matrix (i.e. fix the parameter)
extern void 
StrikeOut(const SqDMatrix& Fin, const ParameterVector& pvin,
	  SqDMatrix& Fout, ParameterVector& pvout,
	  vector<int> ids);

// Sum priors for selected parameters into the Fisher matrix
extern void 
AddPriors(SqDMatrix& F, const ParameterVector& pv,
	  vector<int> ids);

extern void
FisherDump(const SqDMatrix& F, const ParameterVector& pv,
	   ostream& os);

SqDMatrix ADBT(const DMatrix& A, const DVector& D, const DMatrix& B);
double  TrADBT(const DMatrix& A, const DVector& D, const DMatrix& B);
double  TrABT(const DMatrix& A, const DMatrix& B);
double  TrAB(const DMatrix& A, const DMatrix& B);



///////////////////////////////////////////////////////////////////////

// Following are specific types of FisherParameters


class lnDParam: public FisherParameter {
public:
  double z;
  lnDParam(double z_=0.,int t_=FisherParameter::ZBIN): FisherParameter(t_),
    z(z_) {}
  ~lnDParam() {}
  string name() const {
    std::ostringstream oss;
    oss << "lnD(" << z << ")";
    return oss.str();
  }
};

class lnHParam: public FisherParameter {
public:
  double z;
  lnHParam(double z_=0.,int t_=FisherParameter::ZBIN): FisherParameter(t_),
    z(z_) {}
  string name() const {
    std::ostringstream oss;
    oss << "lnH(" << z << ")";
    return oss.str();
  }
};

class lnaParam: public FisherParameter {
public:
  double z;
  lnaParam(double z_=0.,int t_=FisherParameter::ZBIN): FisherParameter(t_),
						       z(z_) {}
  string name() const {
    std::ostringstream oss;
    oss << "dln(a)(" << z << ")";
    return oss.str();
  }
};

class lnGParam: public FisherParameter {
public:
  double z;
  lnGParam(double z_=0.,int t_=FisherParameter::ZBIN): FisherParameter(t_),
    z(z_) {}
  ~lnGParam() {}
  string name() const {
    std::ostringstream oss;
    oss << "lnG(" << z << ")";
    return oss.str();
  }
};

// A "second parameter" for non-linear growth
class nl2Param: public FisherParameter {
public:
  double z;
  nl2Param(double z_=0.,int t_=FisherParameter::ZBIN): FisherParameter(t_),
    z(z_) {}
  ~nl2Param() {}
  string name() const {
    std::ostringstream oss;
    oss << "nl2p(" << z << ")";
    return oss.str();
  }
};

class omegaKParam: public FisherParameter {
public:
  omegaKParam(int t_=FisherParameter::GLOBAL): FisherParameter(t_) {}
  string name() const {return "omegaK";}
};

class omegaQParam: public FisherParameter {
public:
  omegaQParam(int t_=FisherParameter::GLOBAL): FisherParameter(t_) {}
  string name() const {return "omegaQ";}
};

class omegaBParam: public FisherParameter {
public:
  omegaBParam(int t_=FisherParameter::GLOBAL): FisherParameter(t_) {}
  string name() const {return "omegaB";}
};

class omegaMParam: public FisherParameter {
public:
  omegaMParam(int t_=FisherParameter::GLOBAL): FisherParameter(t_) {}
  string name() const {return "omegaM";}
};

class w0Param: public FisherParameter {
public:
  w0Param(int t_=FisherParameter::GLOBAL): FisherParameter(t_) {}
  string name() const {return "w0";}
};

class waParam: public FisherParameter {
public:
  waParam(int t_=FisherParameter::GLOBAL): FisherParameter(t_) {}
  string name() const {return "wa";}
};

class lnAParam: public FisherParameter {
public:
  lnAParam(int t_=FisherParameter::GLOBAL): FisherParameter(t_) {}
  string name() const {return "log(As)";}
};

class nsParam: public FisherParameter {
public:
  nsParam(int t_=FisherParameter::GLOBAL): FisherParameter(t_) {}
  string name() const {return "ns";}
};

} // namespace cosmology

#endif	// FISHERMATRICES_H

