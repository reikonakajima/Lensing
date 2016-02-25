// 	$Id: Matrix.cpp,v 2.25 2007-11-26 19:05:46 garyb Exp $	
// Matrix and vector routines, by Jarvis
#include "Matrix.h"

using std::complex;
using std::ostream;

namespace mv {

  // Declarations of functions (from Numerical Recipes) used herein:
  template <class T>
  void ludcmp(SqMatrix<T>* a, valarray<size_t>* indx, T* d);
  template <class T>
  void lubksb(const SqMatrix<T>& a, const valarray<size_t>& indx, Vector<T>* b);
  void svdcmp(DMatrix* u, SqDMatrix* v, DVector* w);
  template <class T>
  void svbksb(const DMatrix& u, const SqDMatrix& v, const DVector& w,
	      const CSlice_ref<T>& b, Vector<T>* x);
  template <class T>
  void tred2(SqMatrix<T>& a, Vector<T>& d, Vector<T>& e);
  template <class T>
  void tqli(Vector<T>& d, Vector<T>& e, SqMatrix<T>& a);


  template <class T>
  void 
  eigenvectors(const SqMatrix<T>& realSym, 
	       Vector<T>& evals, 
	       SqMatrix<T>& evecs) {
    // ?? verify symmetric matrix?
    evals.Resize(realSym.size());
    evecs=realSym;
    Vector<T> e(realSym.size());
    tred2(evecs,evals,e);
    tqli(evals,e,evecs);
  }

  /////////////////////////////////////////////////////////
  // Following are matrix/vector manipulation building blocks, written
  // in a way that appears to optimize speed on present compilers.
  // For some compilers it is necessary to alter valarray include file
  // to allow taking pointers of const valarray elements.  This is
  // why Matrix.h includes "valarrayKludge.h" instead of <valarray>
  /////////////////////////////////////////////////////////
  template <class T>
  void 
  Slice_ref<T>::swap(Slice_ref<T> rhs) const {
    MVAssert(s.size()==rhs.s.size());
    T* lptr=&(va[s.start()]);
    T* rptr=&(rhs.va[rhs.s.start()]);
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    T dum;
    if (lstride==1 && rstride==1) {
      // stride-1 loop can using SIMD instructions on Pentium.
      for (size_t i=0; i<count; i++) {
	dum = *lptr;
	*lptr = *rptr;
	*rptr = dum;
	lptr++;
	rptr++;
      }
    } else {
      for (size_t i=0; i<count; i++) {
	dum = *lptr;
	*lptr = *rptr;
	*rptr = dum;
	lptr += lstride;
	rptr += rstride;
      }
    }
    return;
  }

#ifdef HAVE_VALARRAY_POINTERS
  template <class T>
  const T CSlice_ref<T>::dot(const CSlice_ref<T>& rhs) const {
    MVAssert(s.size()==rhs.s.size());
    const T* lptr=&(va[s.start()]);
    const T* rptr=&(rhs.va[rhs.s.start()]);
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    T sum=0;
    if (lstride==1 && rstride==1) {
      for (size_t i=0; i<count; i++, lptr++, rptr++) {
      	sum += *lptr * *rptr;
	// Here's a form that Intel will vectorize:
	//for (size_t i=0; i<count; i++) {
      	//sum += lptr[i] * rptr[i];
      } 
    } else {
      // Here's a form that Intel will vectorize:
      //#pragma vector always
      //for (size_t i=0; i<count; i++) {
      //	sum += lptr[i*lstride] * rptr[i*rstride];
      for (size_t i=0; i<count; i++) {
	sum += *lptr * *rptr;
	lptr += lstride;
	rptr += rstride;
      }
    }
    return sum;
  }
  template <class T>
  const T CSlice_ref<T>::dot(const CSlice_ref<T>& rhs1,
			     const CSlice_ref<T>& rhs2) const {
    MVAssert(s.size()==rhs1.s.size());
    MVAssert(s.size()==rhs2.s.size());
    const T* lptr=&(va[s.start()]);
    const T* r1ptr=&(rhs1.va[rhs1.s.start()]);
    const T* r2ptr=&(rhs2.va[rhs2.s.start()]);
    size_t lstride=s.stride();
    size_t r1stride=rhs1.s.stride();
    size_t r2stride=rhs2.s.stride();
    size_t count=s.size();
    T sum=0;
    if (lstride==1 && r1stride==1 && r2stride==1) {
      for (size_t i=0; i<count; i++, lptr++, r1ptr++, r2ptr++) {
      	sum += *lptr * *r1ptr * *r2ptr;
      } 
    } else {
      for (size_t i=0; i<count; i++) {
	sum += *lptr * *r1ptr * *r2ptr;
	lptr += lstride;
	r1ptr += r1stride;
	r2ptr += r2stride;
      }
    }
    return sum;
  }
  template <class T>
  const T CSlice_ref<T>::sum() const {
    const T* lptr=&(va[s.start()]);
    size_t lstride=s.stride();
    size_t count=s.size();
    T sum=0;
    if (lstride==1) {
      // A SIMD-compatible form:
      for (size_t i=0; i<count; i++, lptr++) {
	sum += *lptr;
      }
    } else {
      for (size_t i=0; i<count; i++) {
	sum += *lptr;
	lptr += lstride;
      }
    }
    return sum;
  }
  template <class T>
  void Slice_ref<T>::operator=(const CSlice_ref<T>& rhs) {
    MVAssert(s.size()==rhs.s.size());
    T* lptr=&(va[s.start()]);
    T* rptr=&(rhs.va[rhs.s.start()]);
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    if (lstride==1 && rstride==1) {
      // stride-1 loop can using SIMD instructions on Pentium.
      for (size_t i=0; i<count; i++) {
	*lptr = *rptr;
	lptr++;
	rptr++;
      }
    } else {
      for (size_t i=0; i<count; i++) {
	*lptr = *rptr;
	lptr += lstride;
	rptr += rstride;
      }
    }
    return;
  }

#else
  template <class T>
  const T CSlice_ref<T>::dot(const CSlice_ref<T>& rhs) const {
    MVAssert(s.size()==rhs.s.size());
    size_t lindx=s.start();
    size_t rindx=rhs.s.start();
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    T sum=0;
    if (lstride==1 && rstride==1) {
      for (size_t i=0; i<count; i++) {
	sum += va[lindx]*rhs.va[rindx];
	lindx++;
	rindx++;
      }
    } else {
      for (size_t i=0; i<count; i++) {
	sum += va[lindx]*rhs.va[rindx];
	lindx+=lstride;
	rindx+=rstride;
      }
    }
    return sum;
  }

  template <class T>
  const T CSlice_ref<T>::dot(const CSlice_ref<T>& rhs1,
			     const CSlice_ref<T>& rhs2) const {
    MVAssert(s.size()==rhs1.s.size());
    MVAssert(s.size()==rhs2.s.size());
    size_t lindx=s.start();
    size_t r1indx=rhs1.s.start();
    size_t r2indx=rhs2.s.start();
    size_t lstride=s.stride();
    size_t r1stride=rhs1.s.stride();
    size_t r2stride=rhs2.s.stride();
    size_t count=s.size();
    T sum=0;
    for (size_t i=0; i<count; i++) {
      sum += va[lindx]*rhs1.va[r1indx]*rhs2.va[r2indx];
      lindx+=lstride;
      r1indx+=r1stride;
      r2indx+=r2stride;
    }
    return sum;
  }
  template <class T>
  const T CSlice_ref<T>::sum() const {
    size_t lindx=s.start();
    size_t lstride=s.stride();
    size_t count=s.size();
    T sum=0;
    for (size_t i=0; i<count; i++) {
      sum += va[lindx];
      lindx+=lstride;
    }
    return sum;
  }
  template <class T>
  void Slice_ref<T>::operator=(const CSlice_ref<T>& rhs) {
    MVAssert(s.size()==rhs.s.size());
    size_t lindex=s.start();
    size_t rindex=rhs.s.start();
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    for (size_t i=0; i<count; i++) {
      va[lindex] = rhs.va[rindex];
      lindex += lstride;
      rindex += rstride;
    }
    return;
  }
#endif

  // Rest of Matrix & vector methods:
  template <class T>
  bool Matrix<T>::SVZero() const {return svd && svd->zero;}

  template <class T>
  void Matrix<T>::Write(ostream& fout, double minnonzero) const
  {
    for(size_t i=0;i<GetM();i++) {
      fout << "[ ";
      for(size_t j=0;j<GetN();j++) {
	T temp = Get(i,j);
	//**	if (abs(temp)<minnonzero) fout << "0.0 ";
	//else fout << temp << ' ';
	fout << temp << ' ';
      }
      fout << " ]\n";
    }
  }

  template <class T>
  void Vector<T>::Write(ostream& fout,double minnonzero) const
  {
    for(size_t i=0;i<size();i++) {
      fout << "[ ";
      //**      if (abs(va[i])<minnonzero) fout << "0.0 ";
      //else fout << va[i] << ' ';
      fout << va[i] << ' ';
      fout << " ]\n";
    }
  }

  // ??? this code unchecked
  /**
  template <class T>
  void Matrix<T>::Read(istream& fin, double minnonzero) 
  {
    MVAssert(getM()>0 && getN()>0)
    char c;
    for(size_t i=0;i<GetM();i++) {
      fin >> c;
      for(size_t j=0;j<GetN();j++) {
	T temp;
	//**	if (abs(temp)<minnonzero) fin << "0.0 ";
	//else fin << temp << ' ';
	fin >> temp;
	(*this)(i,j) = temp;
      }
      fin >> c;
    }
  }

  template <class T>
  void Vector<T>::Read(istream& fin,double minnonzero) 
  {
    MVAssert(va.size()>0);
    char c;
    for(size_t i=0;i<size();i++) {
      fin >> c;
      //**      if (abs(va[i])<minnonzero) fin << "0.0 ";
      //else fin << va[i] << ' ';
      fin >> va[i] ;
      fin >> c;
    }
  }
  **/
  template <class T>
  SqMatrix<T> Matrix<T>::InverseATA() const
    // (AtA)^(-1)_ij = Sum_k v_ik v_jk / w_k^2
  {
    MVAssert(svd); // This is only efficient to do if SV decomp already done.
    SqMatrix<T> temp(GetN());
    for(size_t i=0;i<GetN();i++) for(size_t j=0;j<=i;j++) {
      double sum = 0.;
      for(size_t k=0;k<GetN();k++) if (svd->w[k])
	sum += svd->v(i,k)*svd->v(j,k)/pow(svd->w[k],2);
      temp(i,j) = temp(j,i) = sum;
    }
    return temp;
  }

  template<> SqCMatrix CMatrix::InverseATA() const
  {
    MVAssert(false); 
    // This is hard to write, since you have to deal with each component of 
    // the complex numbers as a separate variable.  
    // So, I haven't bothered to write it yet.
    return SqCMatrix(*this);
  }

  template <class T> 
  void Matrix<T>::BackSub(const CSlice_ref<T>& col, Vector<T>* x) const
  {
    MVAssert(col.size() == GetM());
    MVAssert(x->size() == GetN());
    MVAssert(svd);
    svbksb(svd->u,svd->v,svd->w,col,x);
  }

  void SVD::nonNegativeSVs()
  {
    for (int i=0; i<w.size(); i++)
      if (w[i]<0.) {
	// Flip the sign of the SV and the corresponding row of V
	w[i] =-w[i];
	v.GetCol(i) *= -1.;
      }
  }

  bool SVD::symmetric()
  {
    MVAssert(u.GetM()==v.GetM());
    for (int i=0; i<w.size(); i++) {
      double sum = u.GetCol(i).dot(v.GetCol(i));
      if (sum<0.) {
	// Flip the sign of the SV and the corresponding row of V
	v.GetCol(i) *= -1.;
	w[i] =-w[i];
      }
    }
    DMatrix tmp = u - v;
    for (int i=0; i<tmp.GetM(); i++)
      for (int j=0; j<tmp.GetN(); j++)
	if (abs(tmp(i,j))>1e-6) {
	  return false;
	}
    return true;
  }

  template <class T> 
  void SVD::BackSub(const CSlice_ref<T>& col, Vector<T>* x) const
  {
    MVAssert(col.size() == u.GetM());
    MVAssert(x->size() == u.GetN());
    svbksb(u,v,w,col,x);
  }

  // Element of the inverse matrix: ??? Do this efficiently
  double 
  SVD::inverseElement(int i, int j) const {
    int ws = w.size();
    double retval=0.;
    for (int k=0; k<ws; k++)
      if (abs(w[k])>0.)
	retval += v(i,k)*u(j,k)/w[k];
    return retval;
  }

  template<> void CMatrix::BackSub(const CSlice_ref<complex<double> >& col,
				   CVector* x) const
  {
    MVAssert(col.size() == GetM());
    MVAssert(x->size() == GetN());
    MVAssert(svd);
    DVector dcol(2*GetM());
    for(size_t i=0;i<GetM();i++) {
      dcol[2*i] = col[i].real();
      dcol[2*i+1] = col[i].imag();
    }
    DVector dx(2*GetN());
    svbksb(svd->u,svd->v,svd->w,CSlice_ref<double>(dcol),&dx);
    for(size_t i=0;i<GetM();i++) 
      (*x)[i] = complex<double>(dx[2*i],dx[2*i+1]);
  }

  template <class T> 
  void SqMatrix<T>::SqBackSub(Vector<T>* x) const
  {
    MVAssert(x->size() == GetN());
    MVAssert(lud);
    lubksb(lud->lu,lud->indx,x);
  }

  //---------------------------------------------------------------------------

template <class T>
T RCMult(const CSlice_ref<T>& row, const CSlice_ref<T>& col)
{
  //  MVAssert(row.size() == col.size());
  //  T sum=0;
  //  for(size_t i=0;i<row.size();i++) sum += row[i]*col[i];
  //  return sum;
  return row.dot(col);
}

complex<double> RCMult(const CSlice_ref<complex<double> >& row, 
    const CSlice_ref<double>& col)
{
  MVAssert(row.size() == col.size());
  complex<double> sum = 0;
  size_t ss=row.size();
  for(size_t i=0;i<ss;i++) sum += row[i]*col[i];
  return sum;
}

complex<double> RCMult(const CSlice_ref<double>& row,
    const CSlice_ref<complex<double> >& col)
{ return RCMult(col,row); }

template <class T>
Matrix<T>& Matrix<T>::operator*=(const SqMatrix<T>& rhs)
{ 
  MVAssert(GetN() == rhs.GetN());
  Vector<T> temprow(GetN());
  for(size_t i=0;i<GetM();i++) {
    for(size_t j=0;j<GetN();j++) temprow[j] = RCMult(GetRow(i),rhs.GetCol(j));
    SetRow(i,temprow);
  }
  return *this;
}

CMatrix& operator*=(CMatrix& m1, const SqDMatrix& m2)
{
  MVAssert(m1.GetN() == m2.GetN());
  CVector temprow(m1.GetN());
  for(size_t i=0;i<m1.GetM();i++) {
    for(size_t j=0;j<m1.GetN();j++) {
      complex<double> sum = 0;
      for(size_t k=0;k<m1.GetN();k++) sum += m1(i,k) * m2(k,j);
      temprow[j] = sum;
    }
    m1.SetRow(i,temprow);
  }
  return m1;
}

template <class T>
Matrix<T> operator*(const Matrix<T>& m1,const Matrix<T>& m2)
{
  MVAssert(m1.GetN() == m2.GetM());
  Matrix<T> temp(m1.GetM(),m2.GetN());
  //  for(size_t i=0;i<m1.GetM();i++) for(size_t j=0;j<m2.GetN();j++) 
  //    temp(i,j) = RCMult(m1.GetRow(i),m2.GetCol(j));
  for(size_t j=0;j<m2.GetN();j++) {
    valarray<T> tcol(m2.GetCol(j));
    CSlice_ref<T> tslice(tcol, slice(0,tcol.size(),1));
    for(size_t i=0;i<m1.GetM();i++)
      temp(i,j) = m1.GetRow(i).dot(tslice);
  }
  return temp;
}

// Calculate trace of product of two matrices
template <class T>
T TraceAB(const Matrix<T>& A, const Matrix<T>& B) {
  MVAssert(A.getN()==B.getM());
  MVAssert(A.getM()==B.getN());
  T sum=0.;
  for (size_t i=0; i<A.getM(); i++) 
    sum += RCMult(A.GetRow(i),B.GetCol(i));
  return sum;
}

CMatrix operator*(const CMatrix& m1, const DMatrix& m2)
{
  MVAssert(m1.GetN() == m2.GetM());
  CMatrix temp(m1.GetM(),m2.GetN());
  // ?? could speed this up by switching indices.
  for(size_t i=0;i<m1.GetM();i++) for(size_t j=0;j<m2.GetN();j++)
    temp(i,j) = RCMult(m1.GetRow(i),m2.GetCol(j));
  return temp;
}

CMatrix operator*(const DMatrix& m1, const CMatrix& m2)
{
  MVAssert(m1.GetN() == m2.GetM());
  CMatrix temp(m1.GetM(),m2.GetN());
  // ?? could speed this up by switching indices.
  for(size_t i=0;i<m1.GetM();i++) for(size_t j=0;j<m2.GetN();j++)
    temp(i,j) = RCMult(m1.GetRow(i),m2.GetCol(j));
  return temp;
}

template <class T>
Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v)
{
  MVAssert(m.GetN() == v.size());
  Vector<T> temp(m.GetM());
  for(size_t i=0;i<m.GetM();i++)
    temp[i] = RCMult(m.GetRow(i),CSlice_ref<T>(v));
  return temp;
}

template <class T>
Vector<T> operator*(const Vector<T>& v, const Matrix<T>& m)
{
  MVAssert(v.size() == m.GetM());
  Vector<T> temp(m.GetN());
  for(size_t i=0;i<m.GetM();i++)
    Slice_ref<T>(temp).SVadd(v[i], m.GetRow(i));
  return temp;
}

template <class T>
Vector<T>& operator*=(Vector<T>& v, const SqMatrix<T>& m)
{
  MVAssert(v.size() == m.GetM());
  Vector<T> temp = v*m;
  v = temp;
  return v;
}

CVector operator*(const CMatrix& m, const DVector& v)
{
  MVAssert(m.GetN() == v.size());
  CVector temp(m.GetM());
  for(size_t i=0;i<m.GetM();i++)
    temp[i] = RCMult(m.GetRow(i),CSlice_ref<double>(v));
  return temp;
}

CVector operator*(const CVector& v, const DMatrix& m)
{
  MVAssert(v.size() == m.GetM());
  CVector temp(m.GetN());
  for(size_t j=0;j<m.GetN();j++)
    temp[j] = RCMult(CSlice_ref<complex<double> >(v),m.GetCol(j));
  return temp;
}

CVector operator*(const DMatrix& m, const CVector& v)
{
  MVAssert(m.GetN() == v.size());
  CVector temp(m.GetM());
  for(size_t i=0;i<m.GetM();i++)
    temp[i] = RCMult(m.GetRow(i),CSlice_ref<complex<double> >(v));
  return temp;
}

CVector operator*(const DVector& v, const CMatrix& m)
{
  MVAssert(v.size() == m.GetM());
  CVector temp(m.GetN());
  for(size_t j=0;j<m.GetN();j++)
    temp[j] = RCMult(CSlice_ref<double>(v),m.GetCol(j));
  return temp;
}

template <class T>
T operator*(const Vector<T>& v1, const Vector<T>& v2)
{ return RCMult(CSlice_ref<T>(v1),CSlice_ref<T>(v2)); }

complex<double> operator*(const CVector& v1, const DVector& v2)
{ return RCMult(CSlice_ref<complex<double> >(v1),CSlice_ref<double>(v2)); }

complex<double> operator*(const DVector& v1, const CVector& v2)
{ return RCMult(CSlice_ref<double>(v1),CSlice_ref<complex<double> >(v2)); }

template <class T>
SqMatrix<T> TruncateMul(const SqMatrix<T>& m1,const SqMatrix<T>& m2)
{
  size_t order = m1.GetN();
  if (m2.GetN()<order) order=m2.GetN();
  SqMatrix<T> temp(order);
  //for(size_t i=0;i<order;i++) for(size_t j=0;j<order;j++) 
  //  temp(i,j) = RCMult(m1.GetPartRow(i,order),m2.GetPartCol(j,order));
  for(size_t j=0;j<order;j++) {
    valarray<T> tcol(m2.GetPartCol(j,order));
    CSlice_ref<T> tslice(tcol, slice(0,tcol.size(),1));
    for(size_t i=0;i<order;i++)
      temp(i,j) = m1.GetPartRow(i,order).dot(tslice);
  }
  return temp;
}

template <class T>
SqMatrix<T>& TruncateMulEq(SqMatrix<T>& m1,
				const SqMatrix<T>& m2)
{ 
  if (m1.GetN()<=m2.GetN()) {
    size_t order = m1.GetN();
    Vector<T> temprow(order);
    for(size_t i=0;i<order;i++) {
      for(size_t j=0;j<order;j++) temprow[j] = 
		 RCMult(m1.GetRow(i),m2.GetPartCol(j,order));
      m1.SetRow(i,temprow);
    }
  } else {
    size_t order = m2.GetN();
    SqMatrix<T> temp=m1;
    m1.Resize(order,order);
    for(size_t i=0;i<order;i++) {
      for(size_t j=0;j<order;j++) m1(i,j) = 
	  RCMult(temp.GetPartRow(i,order),m2.GetPartCol(j,order));
    }
  }
  return m1;
}

template <class T>
Vector<T> TruncateMul(const Vector<T>& v, const SqMatrix<T>& m)
{
  MVAssert(v.size() <= m.GetM());
  Vector<T> temp(v.size());
  for(size_t j=0;j<v.size();j++)
    temp[j] = RCMult(CSlice_ref<T>(v),m.GetPartCol(j,v.size()));
  return temp;
}

template <class T>
Vector<T> TruncateMul(const SqMatrix<T>& m, const Vector<T>& v)
{
  MVAssert(v.size() <= m.GetM());
  Vector<T> temp(v.size());
  for(size_t j=0;j<v.size();j++)
    temp[j] = RCMult(m.GetPartRow(j,v.size()),CSlice_ref<T>(v));
  return temp;
}

template <class T>
Vector<T>& TruncateMulEq(const SqMatrix<T>& m, Vector<T>& v) {
  //Note this is actually producing M*V, not V*M
  MVAssert(v.size() <= m.GetM());
  Vector<T> temp=v;
  for(size_t j=0;j<v.size();j++)
    v[j] = RCMult(m.GetPartRow(j,v.size()),CSlice_ref<T>(temp));
  return v;
}

template <class T>
Vector<T>& TruncateMulEq(Vector<T>& v, const SqMatrix<T>& m) {
  // This one does v -> v*M
  MVAssert(v.size() <= m.GetM());
  Vector<T> temp=v;
  for(size_t j=0;j<v.size();j++)
    v[j] = RCMult(CSlice_ref<T>(temp),m.GetPartCol(j,v.size()));
  return v;
}
  

template <class T>
Matrix<T> OuterProduct(const Vector<T>& v1, const Vector<T>& v2)
{
  Matrix<T> temp(v1.size(),v1.size());
  for(size_t i=0;i<v1.size();i++)
    for(size_t j=0;j<v2.size();j++)
      temp(i,j) = v1[i]*v2[j];
  return temp;
}

template <class T>
Matrix<T>& Matrix<T>::operator/=(const SqMatrix<T>& rhs)
{
  MVAssert(GetM() == rhs.GetN());
  rhs.SetLUD();
  for(size_t j=0;j<GetN();j++) {
    Vector<T> col(GetCol(j));
    rhs.SqBackSub(&col);
    SetCol(j,col);
  }
  return *this;
}

template <class T>
SqMatrix<T> SqMatrix<T>::Inverse() const
{
  SqMatrix temp(this->GetN());
  temp = T(1);
  temp /= *this;
  return temp;
}

template <class T>
Matrix<T> operator/(const Matrix<T>& m1,const Matrix<T>& m2)
{
  MVAssert(m1.GetM() == m2.GetM());
  Matrix<T> temp(m2.GetN(),m1.GetN());
  m2.SetSVD();
  Vector<T> tempcol(m2.GetN());
  for(size_t j=0;j<m1.GetN();j++) {
    m2.BackSub(m1.GetCol(j),&tempcol);
    temp.SetCol(j,tempcol);
  }
  return temp;
}

CMatrix operator/(const CMatrix& m1, const DMatrix& m2)
{
  MVAssert(m1.GetM() == m2.GetM());
  m2.SetSVD();
  CMatrix temp(m2.GetN(),m1.GetN());
  DVector outrealcol(m2.GetN());
  DVector outimagcol(m2.GetN());
  DVector inrealcol(m1.GetM());
  DVector inimagcol(m1.GetM());
  for(size_t j=0;j<m1.GetN();j++) {
    for(size_t i=0;i<m1.GetM();i++) {
      inrealcol[i] = m1(i,j).real();
      inimagcol[i] = m1(i,j).imag();
    }
    m2.BackSub(Slice_ref<double>(inrealcol),&outrealcol);
    m2.BackSub(Slice_ref<double>(inimagcol),&outrealcol);
    for(size_t i=0;i<m2.GetN();i++) 
      temp.Set(i,j,complex<double>(outrealcol[i],outimagcol[i]));
  }
  return temp;
}

CMatrix operator/(const DMatrix& m1, const CMatrix& m2)
{
  MVAssert(m1.GetM() == m2.GetM());
  m2.SetSVD();
  CMatrix temp(m2.GetN(),m1.GetN());
  CVector incol(m1.GetM());
  CVector outcol(m2.GetN());
  for(size_t j=0;j<m1.GetN();j++) {
    for(size_t i=0;i<m1.GetM();i++) 
      incol[i] = complex<double>(m1(i,j));
    m2.BackSub(Slice_ref<complex<double> >(incol),&outcol);
    temp.SetCol(j,outcol);
  }
  return temp;
}

CMatrix& operator/=(CMatrix& m1, const SqDMatrix& m2)
{
  MVAssert(m1.GetM() == m2.GetM());
  m2.SetLUD();
  DVector tempcolreal(m2.GetN());
  DVector tempcolimag(m2.GetN());
  for(size_t j=0;j<m1.GetN();j++) {
    for(size_t i=0;i<m2.GetN();i++) {
      tempcolreal[i] = m1(i,j).real();
      tempcolimag[i] = m1(i,j).imag();
    }
    m2.SqBackSub(&tempcolreal);
    m2.SqBackSub(&tempcolimag);
    for(size_t i=0;i<m2.GetN();i++) 
      m1.Set(i,j,complex<double>(tempcolreal[i],tempcolimag[i]));
  }
  return m1;
}

CMatrix operator/(const DMatrix& m1, const SqCMatrix& m2)
{
  MVAssert(m1.GetM() == m2.GetM());
  CMatrix temp(m2.GetN(),m1.GetN());
  m2.SetLUD();
  CVector tempcol(m2.GetN());
  for(size_t j=0;j<m1.GetN();j++) {
    for(size_t i=0;i<m2.GetN();i++) 
      tempcol[i] = complex<double>(m1(i,j));
    m2.SqBackSub(&tempcol);
    temp.SetCol(j,tempcol);
  }
  return temp;
}

template <class T>
Vector<T>& Vector<T>::operator/=(const SqMatrix<T>& m)
{
  MVAssert(va.size() == m.GetN());
  m.SetLUD();
  m.SqBackSub(this);
  return *this;
}

template <class T>
Vector<T> operator/(const Vector<T>& v, const Matrix<T>& m)
{
  MVAssert(v.size() == m.GetM());
  m.SetSVD();
  Vector<T> temp(m.GetN());
  m.BackSub(CSlice_ref<T>(v),&temp);
  return temp;
}

CVector operator/(const CVector& v, const DMatrix& m)
{
  MVAssert(v.size() == m.GetM());
  m.SetSVD();
  CVector temp(m.GetN());
  DVector inreal(v.size());
  DVector inimag(v.size());
  DVector outreal(m.GetN());
  DVector outimag(m.GetN());
  for(size_t i=0;i<v.size();i++) {
    inreal[i] = v[i].real();
    inimag[i] = v[i].imag();
  }
  m.BackSub(Slice_ref<double>(inreal),&outreal);
  m.BackSub(Slice_ref<double>(inimag),&outimag);
  for(size_t i=0;i<m.GetN();i++)
    temp.Set(i,complex<double>(outreal[i],outimag[i]));
  return temp;
}

CVector& operator/=(CVector& v, const SqDMatrix& m)
{
  MVAssert(v.size() == m.GetM());
  m.SetLUD();
  DVector tempreal(m.GetN());
  DVector tempimag(m.GetN());
  for(size_t i=0;i<m.GetN();i++) {
    tempreal[i] = v[i].real();
    tempimag[i] = v[i].imag();
  }
  m.SqBackSub(&tempreal);
  m.SqBackSub(&tempimag);
  for(size_t i=0;i<m.GetN();i++)
    v.Set(i,complex<double>(tempreal[i],tempimag[i]));
  return v;
}

template <class T>
LUD<T>::LUD(const SqMatrix<T>& m)
: lu(m),indx(m.GetN()),det(1.)
{
  ludcmp(&lu,&indx,&det);
  for(size_t i=0;i<lu.GetN();i++) det *= lu(i,i);
}

template <class T>
void SqMatrix<T>::SetLUD() const
{
  if (this->isAltered) Changed();
  if (!lud)  lud = new LUD<T>(*this); 
}

template <class T>
SVD::SVD(const Matrix<T>& m) 
: u(m),v(m.GetN()),w(m.GetN()),zero(false)
{
  svdcmp(&u,&v,&w);
}

SVD::SVD(const Matrix<complex<double> >& m)
: u(2*m.GetM(),2*m.GetN()),v(2*m.GetN()),w(2*m.GetN()),zero(false)
{
  for(size_t i=0;i<m.GetM();i++) for(size_t j=0;j<m.GetN();j++) {
    u.Set(2*i,2*j, m.Get(i,j).real());
    u.Set(2*i+1,2*j, m.Get(i,j).imag());
    u.Set(2*i,2*j+1, -m.Get(i,j).imag());
    u.Set(2*i+1,2*j+1, m.Get(i,j).real());
  }
  svdcmp(&u,&v,&w);
}

template <class T>
void Matrix<T>::SetSVD() const
{
  if (isAltered) Changed();
  if (!svd) {
    svd = new SVD(*this);
    SVDThresh(SVDTOL);
  }
}

template <class T>
int Matrix<T>::SVDThresh(double toler,ostream* debugout) const
{
  MVAssert(svd);
  DVector absw(svd->w.size());
  for(size_t i=0;i<absw.size();i++) absw[i] = fabs(svd->w[i]);
  double wmax = absw.GetValArray().max();
  double thresh = wmax*toler;
  if(debugout) (*debugout)<<"SVD wmax = " << wmax ;
  return SVDThreshAbs(thresh, debugout);
}

template <class T>
int Matrix<T>::SVDThreshAbs(double thresh,ostream* debugout) const
{
  MVAssert(svd);
  if(debugout) (*debugout) << "SVD thresh = " << thresh << endl;
  int nzero=0;
  for(size_t i=0;i<svd->w.size();i++) {
    double absw=fabs(svd->w[i]);
    if (absw < thresh) {
      if(debugout) (*debugout)<< "setting w["<< i
			      <<"] = " << svd->w[i]
			      << " to 0.0" << endl;
      svd->w[i] = 0.0;
      svd->zero = true;
      nzero++;
    }
  }
  return nzero;
}

template <class T>
void Matrix<T>::SVDTop(int neigen,ostream* debugout) const
{
  MVAssert(svd);
  MVAssert(neigen > 0)
  MVAssert((size_t) neigen < svd->w.size());
  vector<double> absw(svd->w.size());
  for(size_t i=0;i<absw.size();i++) absw[i] = abs(svd->w[i]);
  vector<double> sortedabsw = absw;
  std::nth_element(sortedabsw.begin(),sortedabsw.begin()+neigen-1,
    sortedabsw.end(),std::greater<double>());
  double thresh = sortedabsw[neigen-1];
  if(debugout) (*debugout)<<"thresh = "<<thresh<<endl;
  for(size_t i=0;i<svd->w.size();i++) {
    if(debugout) (*debugout)<<"absw["<<i<<"] = "<<absw[i]<<endl;
    if (absw[i] < thresh) {
      if(debugout) (*debugout)<<"setting to 0.0\n";
      svd->w[i] = 0.0;
      svd->zero = true;
    }
  }
}


/////////////////////////////////////////////////////////////
// Numerical Recipes routines follow, but they have been altered
// to have value type templated, to use our Matrix/Vector classes,
// and to make use of the above choices for efficient methods of
// doing the innermost loops, e.g. the dot() method.
/////////////////////////////////////////////////////////////


#define TINY 1.0e-20;

template <class T>
void ludcmp(SqMatrix<T>* a,valarray<size_t>* indx,T* d)
{
  MVAssert(indx->size() == a->GetN());
  size_t n = a->GetN();
  DVector vv(n);
  *d=1.0;
  for (size_t i=0;i<n;i++) {
    double big=0.0;
    for (size_t j=0;j<n;j++) {
      double temp = abs((*a)(i,j));
      if (temp > big) big=temp;
    }
    if (big == 0.0) {
	cerr << "mverror(matrix0)\n";
      cerr<<"Singular matrix in routine LUDCMP\n"; throw mverror();
    }
    vv[i]=1.0/big;
  }
  size_t imax = 0;
  for (size_t j=0;j<n;j++) {
    for (size_t i=0;i<j;i++) {
      T sum=(*a)(i,j);
      //      for (size_t k=0;k<i;k++) sum -= (*a)(i,k)*(*a)(k,j);
      if (i>0) sum -= a->partialRow(i,0,i-1).dot(a->partialCol(j,0,i-1));
      (*a)(i,j)=sum;
    }
    double big=0.0;
    for (size_t i=j;i<n;i++) {
      T sum=(*a)(i,j);
      //      for (size_t k=0;k<j;k++)
      //  sum -= (*a)(i,k)*(*a)(k,j);
      if (j>0) sum -= a->partialRow(i,0,j-1).dot(a->partialCol(j,0,j-1));
      (*a)(i,j)=sum;
      double temp = vv[i]*abs(sum);
      if (temp >= big) {
        big=temp;
        imax=i;
      }
    }
    if (j != imax) {
      //      for (size_t k=0;k<n;k++) {
      //  T dum=(*a)(imax,k);
      //  (*a)(imax,k)=(*a)(j,k);
      //  (*a)(j,k)=dum;
      //}
      if (n>0) a->partialRow(imax,0,n-1).swap(a->partialRow(j,0,n-1));
      *d = -(*d);
      vv[imax]=vv[j];
    }
    (*indx)[j]=imax;
    if ((*a)(j,j) == 0.0) (*a)(j,j)=TINY;
    for (size_t i=j+1;i<n;i++) {
      (*a)(i,j) /= (*a)(j,j);
    }
  }
}

#undef TINY

template <class T>
void lubksb(const SqMatrix<T>& a, const valarray<size_t>& indx,
    Vector<T>* b)
{
  size_t n = a.GetN();
  MVAssert(indx.size() == n);
  MVAssert(b->size() == n);
  int ii=-1;
  for (size_t i=0;i<n;i++) {
    size_t ip=indx[i];
    T sum=(*b)[ip];
    (*b)[ip]=(*b)[i];
    if (ii>=0)
      for (size_t j=ii;j<=i-1;j++) sum -= a(i,j)*(*b)[j];
    else if (sum!=0.) ii=i;
    (*b)[i]=sum;
  }
  for (int i=n-1;i>=0;i--) {
    T sum=(*b)[i];
    for (size_t j=i+1;j<n;j++) sum -= a(i,j)*(*b)[j];
    (*b)[i]=sum/a(i,i);
  }
}

#define MMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? abs(a) : -abs(a))

void svdcmp(DMatrix* a, SqDMatrix* v, DVector* w)
{
  const int MaxIterations=40;
  size_t m = a->GetM();
  size_t n = a->GetN();
  MVAssert(w->size() == n);
  MVAssert(v->GetN() == n);
  if (m < n) {
    cerr<<"SVDCMP: You must augment A with extra zero rows"; 
      cerr << "mverror(matrix1)\n";
    throw mverror();
  }

  double at,bt,ct,maxarg1,maxarg2;
  double scale=0.0,anorm=0.0,g=0.0;
  valarray<double> rv1(n);

  for (size_t i=0;i<n;i++) {
    rv1[i]=scale*g;
    g=scale=0.0;
    if (i < m) {
      for (size_t k=i;k<m;k++) scale += abs((*a)(k,i));
      if (scale) {
	double s = 0.0;
        for (size_t k=i;k<m;k++) {
          (*a)(k,i) /= scale;
          s += (*a)(k,i) * (*a)(k,i);
        }
        double f = (*a)(i,i);
        g = -SIGN(sqrt(s),f);
        double h=f*g-s;
        (*a)(i,i)=f-g;
        for (size_t j=i+1;j<n;j++) {
          double s2 = 0.0;
	  //          for (size_t k=i;k<m;k++) s2 += (*a)(k,i)*(*a)(k,j);
	  if (m>0) s2 += a->partialCol(i,i,m-1).dot(a->partialCol(j,i,m-1));
          f=s2/h;
          // for (size_t k=i;k<m;k++) (*a)(k,j) += f*(*a)(k,i); // ???
	  // ???? use AddScalarVector in vectorop template here and below
	  double* lptr=&(*a)(i,j);
	  const double* rptr=&(*a)(i,i);
	  size_t step=a->getN();
	  for (size_t k=i; k<m; k++) {
	    *lptr += f* *rptr;
	    lptr+=step;
	    rptr+=step;
	  }
        }
        for (size_t k=i;k<m;k++) (*a)(k,i) *= scale;
      }
    }
    (*w)[i]=scale*g;
    g=scale=0.0;
    if (i < m && i != n-1) {
      for (size_t k=i+1;k<n;k++) scale += abs((*a)(i,k));
      if (scale) {
        double s = 0.0;
        for (size_t k=i+1;k<n;k++) {
          (*a)(i,k) /= scale;
          s += (*a)(i,k) * (*a)(i,k);
        }
        double f=(*a)(i,i+1);
        g = -SIGN(sqrt(s),f);
        double h=f*g-s;
        (*a)(i,i+1)=f-g;
        for (size_t k=i+1;k<n;k++) rv1[k]=(*a)(i,k)/h;
        for (size_t j=i+1;j<m;j++) {
          double s2 = 0.0;
	  //for (size_t k=i+1;k<n;k++) s2 += (*a)(j,k)*(*a)(i,k);
	  if (n>0) 
	    s2 += a->partialRow(j,i+1,n-1).dot(a->partialRow(i,i+1,n-1));
	  //for (size_t k=i+1;k<n;k++) (*a)(j,k) += s2*rv1[k]; // ???
	  double* lptr=&(*a)(j,i+1);
	  const double* rptr=&rv1[i+1];
	  for (size_t k=i+1; k<n; k++) {
	    *lptr += s2 * *rptr;
	    lptr++;
	    rptr++;
	    }
        }
	//for (size_t k=i+1;k<n;k++) (*a)(i,k) *= scale;
	a->partialRow(i,i+1,n-1)*=scale;
      }
    }
    anorm=MMAX(anorm,(abs((*w)[i])+abs(rv1[i])));
  }

  for (int i=n-1;i>=0;i--) {
    if (g) {
      for (size_t j=i+1;j<n;j++)
        (*v)(j,i)=((*a)(i,j)/(*a)(i,i+1))/g;
      for (size_t j=i+1;j<n;j++) {
        double s2 = 0.0;
	//for (size_t k=i+1;k<n;k++) s2 += (*a)(i,k)*(*v)(k,j);
	 if (n>0)
	   s2 += a->partialRow(i,i+1,n-1).dot(v->partialCol(j,i+1,n-1));
	 // for (size_t k=i+1;k<n;k++) (*v)(k,j) += s2*(*v)(k,i);
	double* lptr=&(*v)(i+1,j);
	const double* rptr=&(*v)(i+1,i);
	size_t step=v->getN();
	for (size_t k=i+1; k<n; k++) {
	  *lptr += s2 * *rptr;
	  lptr+=step;
	  rptr+=step;
	}
      }
    }
    for (size_t j=i+1;j<n;j++) (*v)(i,j)=(*v)(j,i)=0.0;
    (*v)(i,i)=1.0;
    g=rv1[i];
  }

  for (int i=(m<n ? m:n)-1;i>=0;i--) {
    g=(*w)[i];
    for (size_t j=i+1;j<n;j++) (*a)(i,j)=0.0;
    if (g) {
      g=1.0/g;
      for (size_t j=i+1;j<n;j++) {
        double s=0.0;
	//        for (size_t k=i+1;k<m;k++) s += (*a)(k,i)*(*a)(k,j);
	if (m>0)
	  s += a->partialCol(i,i+1,m-1).dot(a->partialCol(j,i+1,m-1));
        double f=(s/(*a)(i,i))*g;
        // for (size_t k=i;k<m;k++) (*a)(k,j) += f*(*a)(k,i);
	double* lptr=&(*a)(i,j);
	const double* rptr=&(*a)(i,i);
	size_t step=a->getN();
	for (size_t k=i; k<m; k++) {
	  *lptr += f* *rptr;
	  lptr+=step;
	  rptr+=step;
	}
      }
      for (size_t j=i;j<m;j++) (*a)(j,i) *= g;
    } else {
      for (size_t j=i;j<m;j++) (*a)(j,i)=0.0;
    }
    ++(*a)(i,i);
  }

  for (int k=n-1;k>=0;k--) {
    for (size_t its=0;its<30;its++) {
      bool flag=true; int ll; int nm=0;
      for (ll=k;ll>=0;ll--) {
        nm=ll-1;
        MVAssert(rv1[0] == 0.0);
        if ((abs(rv1[ll])+anorm) == anorm) {
          flag=false;
          break;
        }
        MVAssert(nm>=0);
        if ((abs((*w)[nm])+anorm) == anorm) break;
      }
      if (flag) {
        MVAssert(ll>0);
        double c=0.0;
        double s=1.0;
        for (int i=ll;i<=k;i++) {
          double f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if (abs(f)+anorm == anorm) break;
          g=(*w)[i];
          double h=hypot(f,g);
          (*w)[i]=h;
          h=1.0/h;
          c=g*h;
          s=(-f*h);
          // for (size_t j=0;j<m;j++) {
          //   double y=(*a)(j,nm);
          //   double z=(*a)(j,i);
          //   (*a)(j,nm)=y*c+z*s;
          //   (*a)(j,i)=z*c-y*s;
          // }
	  double* nmptr=&(*a)(0,nm);
	  double* iptr=&(*a)(0,i);
	  size_t step=a->getN();
	  for (size_t j=0;j<m;j++) {
	    double y=*nmptr;
	    double z=*iptr;
            *nmptr=y*c+z*s;
            *iptr =z*c-y*s;
	    nmptr+=step;
	    iptr +=step;
          }
        }
      }
      double z=(*w)[k];
      if (ll == k) {
        if (z < 0.0) {
          (*w)[k] = -z;
          //for (size_t j=0;j<n;j++) (*v)(j,k)=(-(*v)(j,k));
	  v->partialCol(k,0,n-1)*=-1.;
        }
        break;
      }
      if (its == MaxIterations-1) {
        cerr<< "No convergence in "
	    << MaxIterations 
	    << " SVDCMP iterations; "; 
	cerr << "mverror(matrix3)\n";
        throw mverror(); 
      }
      double x=(*w)[ll];
      nm=k-1;
      double y=(*w)[nm];
      g=rv1[nm];
      double h=rv1[k];
      double f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=hypot(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      double c=1.0,s=1.0;
      for (int j=ll;j<=nm;j++) {
        int i=j+1;
        g=rv1[i];
        y=(*w)[i];
        h=s*g;
        g=c*g;
        z=hypot(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y=y*c;
        // for (size_t jj=0;jj<n;jj++) {
        //   x=(*v)(jj,j);
        //   z=(*v)(jj,i);
        //   (*v)(jj,j)=x*c+z*s;
        //   (*v)(jj,i)=z*c-x*s;
        // }
	double* jptr=&(*v)(0,j);
	double* iptr=&(*v)(0,i);
	size_t step=v->getN();
        for (size_t jj=0;jj<n;jj++) {
          x=*jptr;
          z=*iptr;
          *jptr=x*c+z*s;
          *iptr=z*c-x*s;
	  jptr+=step;
	  iptr+=step;
        }
        z=hypot(f,h);
        (*w)[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=(c*g)+(s*y);
        x=(c*y)-(s*g);
        // for (size_t jj=0;jj<m;jj++) {
        //  y=(*a)(jj,j);
        //  z=(*a)(jj,i);
        //  (*a)(jj,j)=y*c+z*s;
        //  (*a)(jj,i)=z*c-y*s;
        // }
	jptr=&(*a)(0,j);
	iptr=&(*a)(0,i);
	step=a->getN();
        for (size_t jj=0;jj<m;jj++) {
          y=*jptr;
          z=*iptr;
          *jptr=y*c+z*s;
          *iptr=z*c-y*s;
	  jptr+=step;
	  iptr+=step;
        }
      }
      rv1[ll]=0.0;
      MVAssert(rv1[0] == 0.0);
      rv1[k]=f;
      MVAssert(rv1[0] == 0.0);
      (*w)[k]=x;
    }
  }
}


template <class T>
void svbksb(const DMatrix& u, const SqDMatrix& v, const DVector& w,
    const CSlice_ref<T>& b, Vector<T>* x)
{
  size_t m = u.GetM();
  size_t n = u.GetN();
  MVAssert(v.GetN() == n);
  MVAssert(w.size() == n);
  MVAssert(b.size() == m);
  MVAssert(x->size() == n);
  DVector tmp(n);
  for (size_t j=0;j<n;j++) {
    double s=0.0;
    if (w[j] != 0.0) {
      for (size_t i=0;i<m;i++) s += u(i,j)*b[i]; //?? use dot
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (size_t j=0;j<n;j++) {
    double s=0.0;
    for (size_t jj=0;jj<n;jj++) s += v(j,jj)*tmp[jj]; // ?? use dot
    (*x)[j]=s;
  }
}

// Reduction of symmetric matrix to tridiagonal form.
// Numerical recipes routine, altered to work in the context
// of our C++ Matrix objects.
// Matrix a is (symmetric) input matrix, and on ouput has the
// transformation to tridiag form.
// d and e are the diagonal and off-diagonal of the tridiagonal form.

template <class T>
class AddScalarVector {
private:
  T factor;
public:
  AddScalarVector(const T& f): factor(f) {};
  const T operator()(const T a, const T b) {return a + factor*b;}
};

template <class T>
void 
tred2(SqMatrix<T>& a, Vector<T>& d, Vector<T>& e)
{
  size_t k,j,i;
  int l;
  T scale,hh,h,g,f;

  size_t n=a.GetN();
  d.Resize(n);
  e.Resize(n);

  for (i=n-1;i>=1;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<=l;k++)
	scale += abs(a(i,k));
      if (scale == 0.0)
	e[i]=a(i,l);
      else {
	//for (k=0;k<=l;k++) {
	//  a[i,k] /= scale;
	//  h += a[i,k]*a[i,k];
	//}
	a.partialRow(i,0,l)*=1./scale;
	h += a.partialRow(i,0,l).dot(a.partialRow(i,0,l));

	f=a(i,l);
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;
	a(i,l)=f-g;
	f=0.0;
	for (j=0;j<=l;j++) {
	  a(j,i)=a(i,j)/h;
	  g=0.0;
	  //for (k=0;k<=j;k++)
	  //  g += a[j,k]*a[i,k];
	  g += a.partialRow(j,0,j).dot(a.partialRow(i,0,j));
	  //for (k=j+1;k<=l;k++)
	  //  g += a[k,j]*a[i,k];
	  if (l>j) g += a.partialCol(j,j+1,l).dot(a.partialRow(i,j+1,l));

	  e[j]=g/h;
	  f += e[j]*a(i,j);
	}
	hh=f/(h+h);
	for (j=0;j<=l;j++) {
	  f=a(i,j);
	  e[j]=g=e[j]-hh*f;
	  //	  for (k=0;k<=j;k++)
	  //  a[j,k] -= (f*e[k]+g*a[i,k]);
	  a.partialRow(j,0,j).vectorOp(AddScalarVector<T>(-f),
				       e.partialVector(0,j));
	  a.partialRow(j,0,j).vectorOp(AddScalarVector<T>(-g),
				       a.partialRow(i,0,j));

	}
      }
    } else
      e[i]=a(i,l);
    d[i]=h;
  }
  d[0]=0.0;
  e[0]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i,i]; */
  for (i=0;i<n;i++) {
    l=i-1;
    if (i>0 && d[i]) {
      for (j=0;j<=l;j++) {
	g=0.0;
	//for (k=0;k<=l;k++)
	//  g += a[i,k]*a[k,j];
	g += a.partialRow(i,0,l).dot(a.partialCol(j,0,l));
	//for (k=0;k<=l;k++)
	//  a[k,j] -= g*a[k,i];
	a.partialCol(j,0,l).vectorOp(AddScalarVector<T>(-g),
				     a.partialCol(i,0,l));
      }
    }
    d[i]=a(i,i);
    a(i,i)=1.0;
    if (i>0) for (j=0;j<=l;j++) a(j,i)=a(i,j)=0.0;
  }
}

// Find eigenvalues/vectors of tridiagonal symmetric real matrix.
// Adopted for C++ Matrix from Numerical Recipes.
template <class T>
void tqli(Vector<T>& d, Vector<T>& e, SqMatrix<T>& z)
{
  int m,l,iter,i,k;
  T s,r,p,g,f,dd,c,b;

  size_t n=d.size();
  MVAssert (e.size()==n);
  MVAssert (z.GetN()==n);
  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if (abs(e[m])+dd == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) {
	  cerr << "Too many iterations in tqli" << endl;
	  cerr << "mverror(matrix4)\n";
	  throw mverror();
	}
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=hypot(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=hypot(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  //for (k=0;k<n;k++) {
	  //  f=z[k][i+1];
	  //  z[k][i+1]=s*z[k][i]+c*f;
	  //  z[k][i]=c*z[k][i]-s*f;
	  //}
	  T* zptr=&z(0,i+1);
	  for (k=0; k<n; k++) {
	    f=*zptr;
	    *zptr=s*(*(zptr-1)) + c*f;
	    *(zptr-1)=c*(*(zptr-1))-s*f;
	    zptr+=z.GetN();
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}
#undef SIGN
#undef MMAX


#define T double
template void 
tred2<T>(SqMatrix<T>& a, Vector<T>& d, Vector<T>& e);
template void 
tqli<T>(Vector<T>& d, Vector<T>& e, SqMatrix<T>& a);
template void 
eigenvectors(const SqMatrix<T>& realSym, 
	     Vector<T>& evals, 
	     SqMatrix<T>& evecs);

template class Matrix<T>;
template class SqMatrix<T>;
template class Vector<T>;
template class CSlice_ref<T>;
template void SVD::BackSub(const CSlice_ref<T>& col, Vector<T>* x) const;
template SVD::SVD(const Matrix<T>& m);
template Matrix<T> operator*(const Matrix<T>& m1,const Matrix<T>& m2);
template Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v);
template Vector<T> operator*(const Vector<T>& v, const Matrix<T>& m);
template T operator*(const Vector<T>& v1, const Vector<T>& v2);
template Matrix<T> OuterProduct(const Vector<T>& v1, const Vector<T>& v2);
template Matrix<T> operator/(const Matrix<T>& m1,const Matrix<T>& m2);
template Vector<T> operator/(const Vector<T>& v, const Matrix<T>& m);
template T TraceAB(const Matrix<T>& A, const Matrix<T>& B);
#undef T
#define T float
template void 
tred2<T>(SqMatrix<T>& a, Vector<T>& d, Vector<T>& e);
template void 
tqli<T>(Vector<T>& d, Vector<T>& e, SqMatrix<T>& a);
template void 
eigenvectors(const SqMatrix<T>& realSym, 
	     Vector<T>& evals, 
	     SqMatrix<T>& evecs);

template SVD::SVD(const Matrix<T>& m);
template class Matrix<T>;
template class SqMatrix<T>;
template class Vector<T>;
template Matrix<T> operator*(const Matrix<T>& m1,const Matrix<T>& m2);
template Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v);
template Vector<T> operator*(const Vector<T>& v, const Matrix<T>& m);
template T operator*(const Vector<T>& v1, const Vector<T>& v2);
template Matrix<T> OuterProduct(const Vector<T>& v1, const Vector<T>& v2);
template Matrix<T> operator/(const Matrix<T>& m1,const Matrix<T>& m2);
template Vector<T> operator/(const Vector<T>& v, const Matrix<T>& m);
template T TraceAB(const Matrix<T>& A, const Matrix<T>& B);
#undef T
#define T complex<double>
template class mv::Matrix<T>;
template class mv::SqMatrix<T>;
template class mv::Vector<T>;
template Matrix<T> operator*(const Matrix<T>& m1,const Matrix<T>& m2);
template Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v);
template Vector<T> operator*(const Vector<T>& v, const Matrix<T>& m);
template T operator*(const Vector<T>& v1, const Vector<T>& v2);
template Matrix<T> OuterProduct(const Vector<T>& v1, const Vector<T>& v2);
template Matrix<T> operator/(const Matrix<T>& m1,const Matrix<T>& m2);
template Vector<T>& operator*=(Vector<T>& v, const SqMatrix<T>& m);
template Vector<T> operator/(const Vector<T>& v, const Matrix<T>& m);
template SqMatrix<T> TruncateMul(const SqMatrix<T>& m1,const SqMatrix<T>& m2);
template Vector<T> TruncateMul(const Vector<T>& v,const SqMatrix<T>& m);
template Vector<T> TruncateMul(const SqMatrix<T>& m, const Vector<T>& v);
template Vector<T>& TruncateMulEq(Vector<T>& v, const SqMatrix<T>& m);
template Vector<T>& TruncateMulEq(const SqMatrix<T>& m, Vector<T>& v);
template SqMatrix<T>& TruncateMulEq(SqMatrix<T>& m1, const SqMatrix<T>& m2);
template T TraceAB(const Matrix<T>& A, const Matrix<T>& B);
#undef T

} //mv
