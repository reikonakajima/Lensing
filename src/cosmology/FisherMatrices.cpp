// $Id: FisherMatrices.cpp,v 1.4 2006-08-08 00:27:04 garyb Exp $

#include "FisherMatrices.h"
using namespace cosmology;

int
FisherParameter::idCounter=0;

SqDMatrix 
cosmology::ADBT(const DMatrix& A, const DVector& D, const DMatrix& B) {
  Assert(A.getM()==B.getM());
  Assert(A.getN()==D.size());
  Assert(A.getN()==B.getN());
  SqDMatrix ret(A.getM());

  int vsize=D.size();
  for (int i=0; i<ret.getM(); i++)
    for (int j=0; j<ret.getM(); j++) 
      ret(i,j) = A.GetRow(i).dot(B.GetRow(j), D);
  return ret;
}

double  
cosmology::TrADBT(const DMatrix& A, const DVector& D, const DMatrix& B) {
  Assert(A.getM()==B.getM());
  Assert(A.getN()==D.size());
  Assert(A.getN()==B.getN());

  double tr=0.;
  for (int i=0; i<A.getM(); i++)
    tr += A.GetRow(i).dot(B.GetRow(i), D);
  return tr;
}

double  
cosmology::TrABT(const DMatrix& A, const DMatrix& B) {
  Assert (A.getM()==B.getM());
  Assert (A.getN()==B.getN());
  double tr=0;
  for (int i=0; i<A.getM(); i++)
    tr += A.GetRow(i).dot(B.GetRow(i));
  return tr;
}

double  
cosmology::TrAB(const DMatrix& A, const DMatrix& B) {
  Assert (A.getM()==B.getN());
  Assert (A.getN()==B.getM());
  double tr=0;
  for (int i=0; i<A.getM(); i++)
    tr += A.GetRow(i).dot(B.GetCol(i));
  return tr;
}

void 
cosmology::AddPriors(SqDMatrix& F, const ParameterVector& pv,
		     vector<int> ids) {
  for (int i=0; i<ids.size(); i++) {
    int index = pv.findID(ids[i]);
    if (index<0) throw MyException("Can't find parameter in AddPriors");
    double sd = pv[index]->getPrior();
    if (sd>0.) F(index, index) += 1./(sd*sd);
  }
}

void 
cosmology::FisherDump(const SqDMatrix& F, const ParameterVector& pv,
		      ostream& os) {
  Assert(F.getM() == pv.size());
  // First list the parameters and their unmarginalized errors
  for (int i=0; i<pv.size(); i++) {
    double sd=-1.;
    if (F(i,i)>0.) sd = 1./sqrt(F(i,i));
    os << std::setw(2) << i
       << std::setw(20) << pv[i]->name()
       << "  " << sd
       << endl;
  }
  for (int i=0; i<pv.size(); i++)
    for (int j=0; j<pv.size(); j++)
      os << i << " " << j << " " << F(i,j) << endl;
}

    
void 
cosmology::MarginalizeOver(const SqDMatrix& Fin, const ParameterVector& pvin,
			   SqDMatrix& Fout, ParameterVector& pvout,
		vector<int> ids) {

  Assert (Fin.getN() == pvin.size());
  int nCut = ids.size();
  int nKeep = pvin.size() - nCut;
  
  pvout.clear();
  vector<int> indicesCut;
  vector<int> indicesKeep;
  for (int i=0; i<pvin.size(); i++) {
    int id = pvin[i]->getID();
    bool keep=true;
    for (int j=0; j<ids.size(); j++)
      if (id==ids[j]) {
	keep=false;
	break;
      }
    if (keep) {
      indicesKeep.push_back(i);
      pvout.push_back(pvin[i]);
    }
    else indicesCut.push_back(i);
  }
  Assert (indicesCut.size()==nCut); 
  Assert (indicesKeep.size()==nKeep);

  // Resize Fout to nKeep
  Fout.resize(nKeep);
  DMatrix Fkc(nKeep, nCut);
  SqDMatrix Fcc(nCut,0.);

  // Divide original Fisher into four submatrices
  for (int i=0; i<nKeep; i++) {
    int ii = indicesKeep[i];
    for (int j=0; j<nKeep; j++) {
      int jj = indicesKeep[j];
      Fout(i,j) = Fin(ii,jj);
    }
    for (int j=0; j<nCut; j++) {
      int jj = indicesCut[j];
      Fkc(i,j) = Fin(ii,jj);
    }
  }
  for (int i=0; i<nCut; i++) {
    int ii = indicesCut[i];
    for (int j=0; j<nCut; j++) {
      int jj = indicesCut[j];
      Fcc(i,j) = Fin(ii,jj);
    }
  }

  // Do the SVD for Fcc, subtract Fkc Fcc^-1 Fck from Fout:
  mv::SVD ccsvd(Fcc);
  ccsvd.symmetric();
  DMatrix kcv = Fkc * ccsvd.v;
  DMatrix kcu = Fkc * ccsvd.u;

  /*for (int k=0; k<nCut; k++) {
    cerr << "Singular value " << k << " " << ccsvd.w[k] << endl;
    for (int j=0; j<nCut; j++) 
      if (abs(ccsvd.v(j,k)) > 0.02)
	cerr << j << " " << indicesCut[j] << "  " << ccsvd.v(j,k) 
	     << " " << pvin[indicesCut[j]]->name() << endl;
  }*/

  // ??? is this good ???
  const double SMALL=1.e-8;
  // If there's a singular (unconstrained) cut parameter combination,
  // it had better have no influence on the kept parameters:
  for (int k=0; k<nCut; k++) 
    if (ccsvd.w[k]<SMALL)
      for (int i=0; i<nKeep; i++)
	if (false) {
	  //	if (abs(kcv(i,k))>SMALL) {
	  /**/cerr << "kcv(" << i << "," << k << ") = " << kcv(i,k)
		   << " " << pvin[indicesKeep[i]]->name()
		   << " while SV is " << ccsvd.w[k] << endl;
	  for (int j=0; j<nCut; j++) 
	    cerr << j << " " << indicesCut[j] 
		 << " " << pvin[indicesCut[j]]->name() 
		 << "  " << ccsvd.v(j,k) 
		 << "  " << Fkc(i,j)
		 << "  " << ccsvd.v(j,k)*Fkc(i,j)
		 << endl;
	  throw MyException("Singular Fisher marginalization");
	}
  for (int i=0; i<nKeep; i++)
    for (int j=0; j<nKeep; j++) 
      for (int k=0; k<nCut; k++)
	//	if (ccsvd.w[k]>SMALL) Fout(i,j) -= kcu(i,k)*kcv(j,k)/ccsvd.w[k];
	if (ccsvd.w[k]>SMALL) Fout(i,j) -= kcv(i,k)*kcv(j,k)/ccsvd.w[k];
}

void 
cosmology::StrikeOut(const SqDMatrix& Fin, const ParameterVector& pvin,
		     SqDMatrix& Fout, ParameterVector& pvout,
		     vector<int> ids) {

  Assert (Fin.getN() == pvin.size());
  int nCut = ids.size();
  int nKeep = pvin.size() - nCut;
  
  pvout.clear();
  vector<int> indicesKeep;
  for (int i=0; i<pvin.size(); i++) {
    int id = pvin[i]->getID();
    bool keep=true;
    for (int j=0; j<ids.size(); j++)
      if (id==ids[j]) {
	keep=false;
	break;
      }
    if (keep) {
      indicesKeep.push_back(i);
      pvout.push_back(pvin[i]);
    }
  }
  Assert (indicesKeep.size()==nKeep);

  // Resize Fout to nKeep
  Fout.resize(nKeep);

  // Copy over the submatrix
  for (int i=0; i<nKeep; i++) {
    int ii = indicesKeep[i];
    for (int j=0; j<nKeep; j++) {
      int jj = indicesKeep[j];
      Fout(i,j) = Fin(ii,jj);
    }
  }
}

SqDMatrix
cosmology::TransformFisher(const SqDMatrix& Fin, 
			   const ParameterVector& pvin,
			   const DMatrix& project, 
			   const ParameterVector& pvout) {
  Assert (Fin.getN() == pvin.size());
  Assert (project.getN() == pvin.size());
  Assert (project.getM() == pvout.size());
  DMatrix tmp = project * Fin * project.Transpose();
  SqDMatrix Fout(pvout.size());
  for (int i=0; i<pvout.size(); i++)
    for (int j=0; j<pvout.size(); j++)
      Fout(i,j) = tmp(i,j);
  return Fout;
}

