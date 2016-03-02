// 	$Id: Matrix.h,v 2.21 2007-07-23 02:34:16 garyb Exp $	
// Matrix and vector routines, by Jarvis.
//---------------------------------------------------------------------------
#ifndef MatrixH
#define MatrixH

// Defines mathematical matrices and vectors
// Matrix<T> is the generalized matrix, with the typedefs:
// DMatrix = Matrix<double>
// CMatrix = Matrix<complex<double> >
//
// SqMatrix<T> is the generalized square matrix
// SqDMatrix = SqMatrix<double>
// SqCMatrix = SqMatrix<complex<double> >
//
// Vector<T> is the generalized vector
// DVector = Vector<double>
// CVector = Vector<complex<double> >
//
// See Stroustrup(2000), section 22.4 for the general structure
// See Numerical Recipes, sections 2.3 and ?? for the division functions

#include "Std.h"
#include <complex>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm>
#include <cstdlib>
using std::abs;
using std::vector;
using std::complex;
using std::cerr;
using std::endl;


// Use a polluted version of valarray that allows pointers to be
// taken from const va's, to help optimize slice loops.
#include "valarrayKludge.h"
using std::valarray;
using std::slice;

namespace mv {


class mverror {};

#ifndef DEBUG
  #define MVAssert(x)
#else
  #define MVAssert(x) if(!(x)) { \
      cerr<<"Error - MVAssert " #x " failed\n"; \
      cerr<<"on line "<<__LINE__<<" in file "<<__FILE__<<endl; \
      throw mverror();}
#endif

const double SVDTOL = 2.22e-16; // = DBL_EPSILON 
// trouble using DBL_EPSILON on sun's, so make it explicit instead

class SVD;
template <class T> class LUD;
template <class T> class Matrix;
template <class T> class SqMatrix;
template <class T> class Vector;


// First a read-only slice class:
template<class T> 
class CSlice_ref {
  template <class U>
  friend class Slice_ref;
private:
  void operator=(const CSlice_ref& cs) {};  //preclude op=
protected:
  const valarray<T>& va; 
  slice s;
  T ref(size_t i) const 
    { MVAssert(i<s.size()); return va[s.start()+i*s.stride()]; }
public:
  CSlice_ref(const valarray<T>& vv, slice ss) : va(vv), s(ss) {}
  CSlice_ref(const Vector<T>& v) : va(v.GetValArray()), s(0,v.size(),1) {}
  T operator[](size_t i) const
    { MVAssert(i<s.size()); return ref(i); } 
  operator valarray<T>() const { return valarray<T>(va[s]); }
  size_t size() const { return s.size(); }
  // Common operations: first is inner product:
  const T dot(const CSlice_ref<T>& rhs) const;
  // Return accumulated product of THREE vectors
  const T dot(const CSlice_ref<T>& rhs1, const CSlice_ref<T>& rhs2) const;
  const T sum() const;
}; // CSlice_ref

// Now a reference to a writable slice
template<class T> 
class Slice_ref: public CSlice_ref<T> {
private:
  valarray<T>& va; 
  slice s;
  T& ref(size_t i) const 
    { MVAssert(i<s.size()); return va[s.start()+i*s.stride()]; }

public:
  Slice_ref(valarray<T>& vv, slice ss) : va(vv), s(ss), 
    CSlice_ref<T>(vv,ss) {}
  Slice_ref(Vector<T>& v) : va(v.GetValArray()), s(0,v.size(),1),
    CSlice_ref<T>(v.GetValArray(), slice(0,v.size(),1)) {}
  T& operator[](size_t i) const
    { MVAssert(i<s.size()); return ref(i); } 
  operator valarray<T>() const { return valarray<T>(va[s]); }
  size_t size() const { return s.size(); }
  // Common operations:
  void swap(Slice_ref<T> rhs) const;
  void operator=(const CSlice_ref<T>& rhs);

  template <class BinOp>
  void scalarOp(BinOp f, const T rhs) {
    T* lptr=&(va[s.start()]);
    size_t lstride=s.stride();
    size_t count=s.size();
    if (lstride==1) {
      for (size_t i=0; i<count; i++) {
	*lptr=f(*lptr,rhs);
	lptr ++;
      }
    } else {
      for (size_t i=0; i<count; i++) {
	*lptr=f(*lptr,rhs);
	lptr += lstride;
      }
    }
  }



  template <class BinOp>
  void vectorOp(BinOp f, CSlice_ref<T> rhs) {
    MVAssert(s.size()==rhs.s.size());
#ifdef HAVE_VALARRAY_POINTERS
    T* lptr=&(va[s.start()]);
    const T* rptr=&(rhs.va[rhs.s.start()]);
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    if (lstride==1 && rstride==1) {
      for (size_t i=0; i<count; i++) {
	*lptr = f(*lptr,*rptr);
	lptr++;
	rptr++;
      }
    } else {
      for (size_t i=0; i<count; i++) {
	*lptr=f(*lptr,*rptr);
	lptr += lstride;
	rptr += rstride;
      }
    }
#else
    size_t lindex=s.start();
    size_t rindex=rhs.s.start();
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    for (size_t i=0; i<count; i++) {
      va[lindex] = f(va[lindex],rhs.va[rindex]);
      lindex += lstride;
      rindex += rstride;
    }
#endif
  }
  void operator+=(const T rhs) {
    scalarOp(std::plus<T>(), rhs);
  }
  void operator*=(const T rhs) {
    scalarOp(std::multiplies<T>(), rhs);
  }
  void operator+=(CSlice_ref<T> rhs) {
    vectorOp(std::plus<T>(), rhs);
  }
  void operator-=(CSlice_ref<T> rhs) {
    vectorOp(std::minus<T>(), rhs);
  }
  void operator*=(CSlice_ref<T> rhs) {
    vectorOp(std::multiplies<T>(), rhs);
  }
  void operator/=(CSlice_ref<T> rhs) {
    vectorOp(std::divides<T>(), rhs);
  }

  // v += scalar*slice
  void SVadd(const T scalar, CSlice_ref<T> rhs) {
    MVAssert(s.size()==rhs.s.size());
#ifdef HAVE_VALARRAY_POINTERS
    T* lptr=&(va[s.start()]);
    const T* rptr=&(rhs.va[rhs.s.start()]);
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    if (lstride==1 && rstride==1) {
      for (size_t i=0; i<count; i++) {
	*lptr += scalar* (*rptr);
	lptr++;
	rptr++;
      }
    } else {
      for (size_t i=0; i<count; i++) {
	*lptr += scalar* (*rptr);
	lptr += lstride;
	rptr += rstride;
      }
    }
#else
    size_t lindex=s.start();
    size_t rindex=rhs.s.start();
    size_t lstride=s.stride();
    size_t rstride=rhs.s.stride();
    size_t count=s.size();
    for (size_t i=0; i<count; i++) {
      va[lindex] += scalar * rhs.va[rindex];
      lindex += lstride;
      rindex += rstride;
    }
#endif
  }
};

//---------------------------------------------------------------------------


template <class T> class Matrix {

public:
  Matrix() : m(0),n(0),va(),svd(0),isAltered(false) {}
  Matrix(size_t mm,size_t nn,T val=0) 
    : m(mm),n(nn),va(val,m*n),svd(0),isAltered(false) {}
  Matrix(const Matrix& rhs)
    : m(rhs.m),n(rhs.n),va(rhs.va),svd(0),isAltered(false) {}
  template <class U> 
    explicit Matrix(const Matrix<U>& rhs)
      : m(rhs.GetM()),n(rhs.GetN()),va(m*n),isAltered(false),svd(0)
    { for(size_t i=0;i<va.size();i++) va[i] = T(rhs.GetValArray()[i]); }
  Matrix(size_t mm,size_t nn,const valarray<T>& vv)
    : m(mm),n(nn),va(vv),svd(0),isAltered(false) 
    { MVAssert(va.size() == m*n); }
  Matrix<T>& operator=(const Matrix<T>& rhs)
    { if (&rhs == this) return *this;
      m=rhs.m;n=rhs.n;va.resize(m*n);
      MVAssert(rhs.m == m && rhs.n == n);
      Changed();
      va = rhs.va; 
      return *this; }
  virtual ~Matrix() { if(svd) delete svd; }
  void Resize(uint mm, uint nn) {resize(mm,nn);}
  void resize(uint mm, uint nn) 
    { m = mm; n = nn; va.resize(m*n); Changed(); }
  Matrix<T>& operator=(const T &value)
    { for(size_t i=0;i<va.size();i++) va[i] = value;  Changed(); return *this;}

  size_t GetM() const {return m;}
  size_t getM() const {return m;}
  size_t GetN() const {return n;}
  size_t getN() const {return n;}
  const valarray<T>& GetValArray() const {return va;}
  valarray<T>& GetValArray() {return va;}

  slice row(size_t i) const {return slice(i*n,n,1); }
  slice col(size_t i) const {return slice(i,m,n); }
  size_t index(size_t i,size_t j) const { return i*n+j; }

  T& operator()(size_t i,size_t j) 
    { MVAssert(i<m && j<n);
      isAltered = true;	// Set a flag for element alteration
      return va[index(i,j)]; }
  T operator()(size_t i,size_t j) const
    { MVAssert(i<m && j<n);
      return va[index(i,j)]; }
  CSlice_ref<T> operator[](size_t i) const
    { MVAssert(i<m);
      return CSlice_ref<T>(va,row(i)); }
  Slice_ref<T> operator[](size_t i) 
    { MVAssert(i<m);
      return Slice_ref<T>(va,row(i)); }
  T Get(size_t i,size_t j) const
    { return (*this)(i,j); }
  CSlice_ref<T> GetRow(size_t i) const
    { MVAssert(i<m);
      return CSlice_ref<T>(va,row(i)); }
  CSlice_ref<T> GetCol(size_t j) const
    { MVAssert(j<n);
      return CSlice_ref<T>(va,col(j)); }
  Slice_ref<T> GetRow(size_t i)
    { MVAssert(i<m);
      return Slice_ref<T>(va,row(i)); }
  Slice_ref<T> GetCol(size_t j)
    { MVAssert(j<n);
      return Slice_ref<T>(va,col(j)); }

  // Parts of rows & columns:
  // get slice for row i from cols jstart to jend
  slice partrow(size_t i, size_t jstart, size_t jend) const {
    return slice(i*n+jstart, jend-jstart+1 ,1); }
  // and col i, rows jstart to jend
  slice partcol(size_t j, size_t istart, size_t iend) const {
    return slice(istart*n+j, iend-istart+1 ,n); }

  // Truncated row/col:
  CSlice_ref<T> GetPartRow(size_t i, size_t k) const
    { MVAssert(i<m && k<=n);
      return CSlice_ref<T>(va,partrow(i,0,k-1)); }
  CSlice_ref<T> GetPartCol(size_t j, size_t k) const
    { MVAssert(j<n && k<=m);
      return CSlice_ref<T>(va,partcol(j,0,k-1)); }

  Slice_ref<T> partialRow(size_t i, size_t jstart, size_t jend){
    MVAssert(i<m && jstart<=jend && jend<n);
      return Slice_ref<T>(va,partrow(i,jstart,jend)); }
  Slice_ref<T> partialCol(size_t j, size_t istart, size_t iend){
    MVAssert(j<n && istart<=iend && iend<m);
      return Slice_ref<T>(va,partcol(j,istart,iend)); }

  void Set(size_t i,size_t j,T t)
    { MVAssert(i<m && j<n);
    isAltered = true;	// Set a flag for element alteration
      va[index(i,j)] = t; }
  void SetRow(size_t i,const valarray<T>& newrow)
    { MVAssert(newrow.size() == n);
      Changed();
      va[row(i)] = newrow; }
  void SetRow(size_t i,const valarray<T>& newrow, size_t start, size_t end)
    { MVAssert(end < newrow.size());
      MVAssert(end+1-start == n);
      Changed();
      va[row(i)] = newrow[slice(start,end+1-start,1)]; }
  void SetCol(size_t j,const valarray<T>& newcol)
    { MVAssert(newcol.size() == m);
      Changed();
      va[col(j)] = newcol; }
  void SetCol(size_t j,const valarray<T>& newcol, size_t start, size_t end)
    { MVAssert(end < newcol.size());
      MVAssert(end+1-start == m);
      Changed();
      va[col(j)] = newcol[slice(start,end+1-start,1)]; }
  void Write(std::ostream& fout,double minnonzero=1.e-30) const;
  void Read(std::istream& fin,double minnonzero=1.e-30);
  void Zero() {zero();}	//backwards compatibility
  void zero()
    { Changed(); va = 0; }

  Matrix& operator+=(const Matrix& rhs)
    { MVAssert(rhs.m == m && rhs.n == n);
      Changed();
      va += rhs.va; 
      return *this; }
  Matrix& operator-=(const Matrix& rhs)
    { MVAssert(rhs.m == m && rhs.n == n);
      Changed();
      va -= rhs.va; 
      return *this; }
  Matrix& operator*=(T x)
    { Changed();
      va *= x; 
      return *this; }
  Matrix& operator/=(T x)
    { Changed();
      va /= x; 
      return *this; }
  Matrix& operator*=(const SqMatrix<T>& rhs);
  Matrix& operator/=(const SqMatrix<T>& rhs);

  Matrix operator-() const {return Matrix(m,n,-va);}
  bool operator==(const Matrix& rhs)
    { if(m != rhs.m || n != rhs.n) return false;
      for(size_t i=0;i<m*n;i++) if (va[i]!= rhs.va[i]) return false;
      return true;  }
  bool operator!=(const Matrix& rhs)
    { return !(*this == rhs); }

  Matrix Transpose() const
    { Matrix temp(n,m);
      for(size_t i=0;i<n;i++) temp.SetRow(i,va[col(i)]);
      return temp; }
  SqMatrix<T> InverseATA() const;

  bool SVZero() const;
  void SetSVD() const;
  // zero out SVD elements smaller than some fraction of the largest.
  // returns # of zero elements in SVD
  int SVDThresh(double toler, std::ostream* debugout=0) const;
  // zero out SVD elements below some absolute threshold
  int SVDThreshAbs(double toler, std::ostream* debugout=0) const;
  void SVDTop(int neigen, std::ostream* debugout=0) const;
  void BackSub(const CSlice_ref<T>& col, Vector<T>* x) const;

protected:

  size_t m,n;
  valarray<T> va;

  virtual void Changed() const {if(svd) delete svd; svd=0; isAltered=false;}

  mutable SVD* svd;
  mutable bool isAltered;

}; // Matrix

//---------------------------------------------------------------------------

template <class T> class SqMatrix : public Matrix<T> {

public:

  SqMatrix() : Matrix<T>(),lud(0) {}
  explicit SqMatrix(size_t nn,T val=0) : Matrix<T>(nn,nn,val),lud(0) {}
  SqMatrix(size_t mm,size_t nn,T val=0) : Matrix<T>(mm,nn,val),lud(0)
    { MVAssert(m==n); }
  SqMatrix(size_t nn, const valarray<T>& vv) : Matrix<T>(nn,nn,vv),lud(0) 
    { MVAssert(vv.size() == nn*nn); }
  explicit SqMatrix(const Matrix<T>& rhs) : Matrix<T>(rhs),lud(0)
    { MVAssert(m==n); }
  SqMatrix(const SqMatrix<T>& rhs) : Matrix<T>(rhs),lud(0) {}
  template <class U>
    explicit SqMatrix(const SqMatrix<U>& rhs) : Matrix<T>(rhs),lud(0)  {}
  Matrix<T>& operator=(const Matrix<T>& rhs) 
    { if (&rhs == this) return *this;
      if (rhs.getM() != rhs.getN()) throw mverror();
      {this->m=rhs.GetM();this->n=rhs.GetN();this->va.resize(this->m*this->n);}
      MVAssert(rhs.GetM() == GetN() && rhs.GetN() == GetN());
      Changed();
      this->va = rhs.GetValArray(); 
      return *this; }
  Matrix<T>& operator=(const SqMatrix<T>& rhs)
    { if (&rhs == this) return *this;
      {this->m=this->n=rhs.GetN(); this->va.resize(this->n*this->n); }
      MVAssert(rhs.GetN() == GetN() && rhs.GetM() == GetM());
      Changed();
      this->va = rhs.GetValArray(); 
      return *this; }
  ~SqMatrix() { if(lud) delete lud; }

  slice diag() { return slice(0,this->n,this->n+1); }
  size_t size() const {return this->GetN();}
  void resize(int newN) {Changed(); 
    this->m=this->n=newN;
    this->va.resize(this->n*this->n);
  }
  SqMatrix<T> operator=(T val)
    { Changed();
      this->va = 0; this->va[diag()] = val; 
      return *this; }
  SqMatrix<T> operator+=(T x)
    { Changed();
      this->va[diag()] += valarray<T>(x,this->n); 
      return *this; }
  SqMatrix<T> operator-=(T x)
    { Changed();
      this->va[diag()] -= valarray<T>(x,this->n); 
      return *this; }
  SqMatrix<T> operator+=(const Matrix<T>& rhs) 
    { Matrix<T>::operator+=(rhs); 
      return *this; }
  SqMatrix<T> operator-=(const Matrix<T>& rhs) 
    { Matrix<T>::operator-=(rhs); 
      return *this; }

  SqMatrix operator-() const {return SqMatrix(this->n,-this->va); }
  SqMatrix Inverse() const;
  SqMatrix Transpose() const
    { SqMatrix temp(this->n);
      for(size_t i=0;i<this->n;i++) temp.SetRow(i,this->va[this->col(i)]);
      return temp; }

  void SetLUD() const;
  void SqBackSub(Vector<T>* x) const;
  T Det() const
    { if (!lud) SetLUD();
      return lud->det; }
  T Tr() const {
    T tr; 
    for(size_t i=0;i<this->n;i++) tr+=this->va[Matrix<T>::index(i,i)];
    return tr;
  }
protected:

  void Changed() const {
    if(lud) delete lud; lud=0; 
    if(this->svd) delete this->svd; this->svd=0; 
    this->isAltered = false;}
  mutable LUD<T>* lud;

}; // SqMatrix

//---------------------------------------------------------------------------

template <class T> class Vector {

public:
  Vector() : va() {}
  explicit Vector(size_t n,T val=0) : va(val,n) {}
  Vector(const valarray<T>& vv) : va(vv) {}
  Vector(const vector<T>& vv) : va(vv.size()) 
    { for(size_t i=0;i<vv.size();i++) va[i] = vv[i]; }
  Vector(const Vector& rhs) : va(rhs.va) {}
  explicit Vector(const CSlice_ref<T> rhs): va(valarray<T> (rhs)) {};
  template <class U>
    explicit Vector(const Vector<U>& rhs)
    : va(rhs.size())
    { for(size_t i=0;i<va.size();i++) va[i] = T(rhs[i]); }
  Vector& operator=(const Vector& rhs)
    { if (&rhs == this) return *this;
      va.resize(rhs.size()); 
      MVAssert(rhs.va.size() == va.size());
      va = rhs.va; 
      return *this; }
  void Resize(int n) { resize(n); }
  void resize(int n) { va.resize(n); }

  size_t size() const {return va.size();}
  const valarray<T>& GetValArray() const {return va;}
  valarray<T>& GetValArray() {return va;}
  operator const valarray<T>&() const {return va;}
  CSlice_ref<T> GetPartVector(size_t i) const
    { MVAssert(i<=va.size());
      return CSlice_ref<T>(va,slice(0,i,1)); }
  CSlice_ref<T> partialVector(size_t istart, size_t iend) const
    { MVAssert(istart<=iend && iend<va.size());
      return CSlice_ref<T>(va,slice(istart,iend-istart+1,1)); }

  Slice_ref<T> GetPartVector(size_t i)
    { MVAssert(i<=va.size());
      return Slice_ref<T>(va,slice(0,i,1)); }
  Slice_ref<T> partialVector(size_t istart, size_t iend)
    { MVAssert(istart<=iend && iend<va.size());
      return Slice_ref<T>(va,slice(istart,iend-istart+1,1)); }

  T& operator[](size_t i) 
    { MVAssert(i<va.size());
      return va[i]; }
  T operator[](size_t i) const 
    { MVAssert(i<va.size());
      return va[i]; }
  
  void Set(size_t i,T t)
    { MVAssert(i<va.size());
      va[i] = t; }
  void Write(std::ostream& fout,double minnonzero=0.) const;
  void Read(std::istream& fin,double minnonzero=0.);
  void Zero() {zero();}	//backwards compatibility
  void zero()
    { va = 0; }

  Vector<T>& operator*=(const Vector& rhs)
    { Slice_ref<T>(*this) *= CSlice_ref<T>(rhs);
      return *this; }
  Vector<T>& operator*=(const CSlice_ref<T> rhs)
    { Slice_ref<T>(*this) *= rhs;
      return *this; }
  Vector<T>& operator+=(const Vector& rhs)
    { Slice_ref<T>(*this) += CSlice_ref<T>(rhs);
      return *this; }
  Vector<T>& operator-=(const Vector& rhs)
    { Slice_ref<T>(*this) -= CSlice_ref<T>(rhs);
      return *this; }
  Vector<T>& operator+=(T x)
    { va += x; 
      return *this; }
  Vector<T>& operator*=(T x)
    { va *= x; 
      return *this; }
  Vector<T>& operator/=(T x)
    { va /= x; 
      return *this; }
  Vector<T>& operator*=(const SqMatrix<T>& rhs);
  Vector<T>& operator/=(const SqMatrix<T>& rhs);
  Vector<T>& operator=(const T &value) {
    size_t ss=va.size();
    for(size_t i=0;i<ss;i++) va[i] = value;  
    return *this;
  }

  const T dot(const Vector<T>& rhs) const {
    return CSlice_ref<T>(*this).dot(CSlice_ref<T>(rhs));
  }
  const T dot(const CSlice_ref<T> rhs) const {
    return CSlice_ref<T>(*this).dot(rhs);
  }

  Vector<T> operator-() const {return Vector<T>(-va);}
  bool operator==(const Vector<T>& rhs) const
    { if (rhs.va.size() != va.size()) return false;
    size_t ss=va.size();
      for(size_t i=0;i<ss;i++) if (va[i] != rhs.va[i]) return false;
      return true; }

private: 
  valarray<T> va;
}; // Vector

//---------------------------------------------------------------------------
struct SVD { 
  Matrix<double> u;
  SqMatrix<double> v;
  Vector<double> w;
  bool zero;
  template <class T>
  SVD(const Matrix<T>& m);
  SVD(const Matrix<complex<double> >& m);
  template <class T>
  void BackSub(const CSlice_ref<T>& col, Vector<T>* x) const;
  double inverseElement(int i, int j) const;
  void nonNegativeSVs();	//Make all SVs non-negative
  void sortSVs();       // sort by SVs
  bool symmetric();	//Test u=v, flipping signs if needed.
};

template <class T>
struct LUD {
  SqMatrix<T> lu;
  valarray<size_t> indx;
  T det;
  LUD(const SqMatrix<T>& m);
};

//---------------------------------------------------------------------------

typedef Matrix<double> DMatrix;
typedef Matrix<complex<double> > CMatrix;
typedef SqMatrix<double> SqDMatrix;
typedef SqMatrix<complex<double> > SqCMatrix;
typedef Vector<double> DVector;
typedef Vector<complex<double> > CVector;

//---------------------------------------------------------------------------

// Define operators for Matrices and Vectors

//
// Output
//

template <class T>
inline std::ostream& operator<<(std::ostream& fout, const Matrix<T>& m)
{ m.Write(fout); return fout;}

template <class T>
inline std::ostream& operator<<(std::ostream& fout, const Vector<T>& v)
{ v.Write(fout); return fout;}

//
// Input
//

template <class T>
inline std::istream& operator>>(std::istream& fin, Matrix<T>& m)
{ m.Read(fin); return fin;}

template <class T>
inline std::istream& operator>>(std::istream& fin, Vector<T>& v)
{ v.Read(fin); return fin;}

//
// Matrix +- Matrix
//

template <class T>
inline Matrix<T> operator+(const Matrix<T>& m1,const Matrix<T>& m2)
{
  Matrix<T> temp = m1;
  temp += m2;
  return temp;
}

inline CMatrix operator+(const DMatrix& m1,const CMatrix& m2)
{
  CMatrix temp(m1);
  temp += m2;
  return temp;
}

inline CMatrix operator+(const CMatrix& m1,const DMatrix& m2)
{ return m2+m1; }

template <class T>
inline Matrix<T> operator-(const Matrix<T>& m1,const Matrix<T>& m2)
{
  Matrix<T> temp = m1;
  temp -= m2;
  return temp;
}

inline CMatrix operator-(const CMatrix& m1,const DMatrix& m2)
{ return -m2 + m1; }

inline CMatrix operator-(const DMatrix& m1,const CMatrix& m2)
{ return -m2 + m1; }

//
// Matrix +- Scalar
//

template <class T>
inline SqMatrix<T> operator+(const SqMatrix<T>& m, T x)
{
  SqMatrix<T> temp = m;
  temp += x;
  return temp;
}

template <class T>
inline SqMatrix<T> operator-(const SqMatrix<T>& m, T x)
{ return m + (-x); }

template <class T>
inline SqMatrix<T> operator+(T x,const SqMatrix<T>& m)
{ return m + x; }

template <class T>
inline SqMatrix<T> operator-(T x,const SqMatrix<T>& m)
{ return -m + x; }

inline SqCMatrix operator+(const SqDMatrix& m, complex<double> x)
{ return SqCMatrix(m) + x; }

inline SqCMatrix operator-(const SqDMatrix& m, complex<double> x)
{ return SqCMatrix(m) - x; }

inline SqCMatrix operator+(complex<double> x,const SqDMatrix& m)
{ return x + SqCMatrix(m); }

inline SqCMatrix operator-(complex<double> x,const SqDMatrix& m)
{ return x - SqCMatrix(m); }


//
// Vector +- Vector
//

template <class T>
inline Vector<T> operator+(const Vector<T>& v1,const Vector<T>& v2)
{
  Vector<T> temp = v1;
  temp += v2;
  return temp;
}

template <class T>
inline Vector<T> operator-(const Vector<T>& v1,const Vector<T>& v2)
{
  Vector<T> temp = v1;
  temp -= v2;
  return temp;
}

inline CVector operator+(const CVector& v1,const DVector& v2)
{ return v1 + CVector(v2); }

inline CVector operator-(const CVector& v1,const DVector& v2)
{ return v1 - CVector(v2); }

inline CVector operator+(const DVector& v1,const CVector& v2)
{ return CVector(v1) + v1; }

inline CVector operator-(const DVector& v1,const CVector& v2)
{ return CVector(v1) - v2; }


//
// Matrix * Matrix
//

template <class T>
Matrix<T> operator*(const Matrix<T>& m1,const Matrix<T>& m2);

CMatrix operator*(const DMatrix& m1,const CMatrix& m2);

CMatrix operator*(const CMatrix& m1,const DMatrix& m2);

template <class T>
inline SqMatrix<T> operator*(const SqMatrix<T>& m1,const SqMatrix<T>& m2)
{
  SqMatrix<T> temp = m1;
  temp *= m2;
  return temp;
}

CMatrix& operator*=(CMatrix& m1,const SqDMatrix& m2);

inline SqCMatrix operator*(const SqDMatrix& m1,const SqCMatrix& m2)
{ return SqCMatrix(m1) * m2; }

inline SqCMatrix operator*(const SqCMatrix& m1,const SqDMatrix& m2)
{
  SqCMatrix temp = m1;
  temp *= m2;
  return temp;
}

// Truncating square matrix/vector multiplications:
template <class T>
SqMatrix<T> TruncateMul(const SqMatrix<T>& m1,const SqMatrix<T>& m2);

template <class T>
Vector<T> TruncateMul(const Vector<T>& v,const SqMatrix<T>& m);

template <class T>
Vector<T> TruncateMul(const SqMatrix<T>& m, const Vector<T>& v);

template <class T>
Vector<T>& TruncateMulEq(Vector<T>& v, const SqMatrix<T>& m);

template <class T>
Vector<T>& TruncateMulEq(const SqMatrix<T>& m, Vector<T>& v);

template <class T>
SqMatrix<T>& TruncateMulEq(SqMatrix<T>& m1, const SqMatrix<T>& m2);

//
// Matrix * Scalar
//

template <class T> 
inline SqMatrix<T> operator*(const SqMatrix<T>& m, T x)
{
  SqMatrix<T> temp = m;
  temp *= x;
  return temp;
}

template <class T> 
inline SqMatrix<T> operator*(T x, const SqMatrix<T>& m)
{ return m*x; }

template <class T> 
inline Matrix<T> operator*(const Matrix<T>& m, T x)
{
  Matrix<T> temp = m;
  temp *= x;
  return temp;
}

template <class T> 
inline Matrix<T> operator*(T x, const Matrix<T>& m)
{ return m*x; }

inline CMatrix operator*(const DMatrix& m,complex<double> x)
{ return CMatrix(m) * x; }

inline CMatrix operator*(complex<double> x,const DMatrix& m)
{ return CMatrix(m) * x; }


//
// Vector * Scalar
//

template <class T>
inline Vector<T> operator*(const Vector<T>& v, T x)
{
  Vector<T> temp = v;
  temp *= x;
  return temp;
}

template <class T> 
inline Vector<T> operator*(T x, const Vector<T>& v)
{ return v*x; }

inline CVector operator*(const DVector& v, complex<double> x)
{ return CVector(v) * x; }

inline CVector operator*(complex<double> x, const DVector& v)
{ return CVector(v) * x; }


//
// Matrix * Vector
//

template <class T>
Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v);

template <class T>
Vector<T> operator*(const Vector<T>& v, const Matrix<T>& m);

CVector operator*(const CMatrix& m, const DVector& v);

CVector operator*(const DMatrix& m, const CVector& v);

CVector operator*(const CVector& v, const DMatrix& m);

CVector operator*(const DVector& v, const CMatrix& m);


//inline CVector& operator*=(const DVector& v, const SqDMatrix& m)
//{ v = v * m; return v; }

//
// Vector * Vector (inner product)
//

template <class T>
T operator*(const Vector<T>& v1, const Vector<T>& v2);

complex<double> operator*(const CVector& v1, const DVector& v2);

complex<double> operator*(const DVector& v1, const CVector& v2);

template <class T>
Matrix<T> OuterProduct(const Vector<T>& v1, const Vector<T>& v2);

inline CMatrix OuterProduct(const CVector& v1, const DVector& v2)
{ return OuterProduct(v1,CVector(v2)); }

inline CMatrix OuterProduct(const DVector& v1, const CVector& v2)
{ return OuterProduct(CVector(v1),v2); }


//
// Matrix / Matrix
//

template <class T>
Matrix<T> operator/(const Matrix<T>& m1,const Matrix<T>& m2);

CMatrix operator/(const CMatrix& m1, const DMatrix& m2);

CMatrix operator/(const DMatrix& m1, const CMatrix& m2);

template <class T>
inline Matrix<T> operator/(const Matrix<T>& m1, const SqMatrix<T>& m2)
{
  Matrix<T> temp = m1;
  temp /= m2;
  return temp;
}

CMatrix& operator/=(CMatrix& m1, const SqDMatrix& m2);

inline CMatrix operator/(const CMatrix& m1, const SqDMatrix& m2)
{
  CMatrix temp = m1;
  temp /= m2;
  return temp;
}

CMatrix operator/(const DMatrix& m1, const SqCMatrix& m2);

template <class T>
inline SqMatrix<T> operator/(const SqMatrix<T>& m1, const SqMatrix<T>& m2)
{
  SqMatrix<T> temp = m1;
  temp /= m2;
  return temp;
}

//
// Matrix / Scalar
//

template <class T>
inline Matrix<T> operator/(const Matrix<T>& m, T x)
{
  Matrix<T> temp = m;
  temp /= x;
  return temp;
}

template <class T>
inline Matrix<T> operator/(T x, const SqMatrix<T>& m)
{ return m.Inverse() * x; }

CMatrix operator/(const DMatrix& m, complex<double> x);
inline CMatrix operator/(complex<double> x, const SqDMatrix& m)
{ return m.Inverse() * x; }


//
// Vector / Scalar
//

template <class T>
inline Vector<T> operator/(const Vector<T>& v, T x)
{
  Vector<T> temp = v;
  temp /= x;
  return temp;
}

inline CVector operator/(const DVector& v, complex<double> x)
{ return CVector(v)/x; }

//
// Vector / Matrix
//

template <class T>
Vector<T> operator/(const Vector<T>& v, const Matrix<T>& m);

CVector operator/(const CVector& v, const DMatrix& m);

inline CVector operator/(const DVector& v, const CMatrix& m)
{ return CVector(v) / m; }

template <class T>
inline Vector<T> operator/(const Vector<T>& v, const SqMatrix<T>& m)
{
  Vector<T> temp = v;
  temp /= m;
  return temp;
}

CVector& operator/=(CVector& v, const SqDMatrix& m);

inline CVector operator/(const CVector& v, const SqDMatrix& m)
{
  CVector temp = v;
  temp /= m;
  return temp;
}

inline CVector operator/(const DVector& v, const SqCMatrix& m)
{ return CVector(v) / m; }

//
// Conjugates, Norm, etc.
//

/*
template <class T>
inline Vector<complex<T> > Conj(const Vector<complex<T> >& v)
{ return Vector<complex<T> >(v.GetValArray().apply(&conj)); }

template <class T>
inline Matrix<complex<T> > Conj(const Matrix<complex<T> >& m)
{ return Matrix<complex<T> >(m.GetM(),m.GetN(),m.GetValArray().apply(&conj)); }
*/
template <class T>
inline double Norm(const Vector<T>& v)
{ return v*v; }

template <class T>
inline double Norm(const Vector<complex<T> >& v)
{
  double sum=0;
  for(size_t i=0;i<v.size();i++) sum += norm(v[i]);
  return sum;
}

template <class T>
T TraceAB(const Matrix<T>& A, const Matrix<T>& B);

// Following assumes that matrix is real & symmetric:
template <class T>
void eigenvectors(const SqMatrix<T>& realSym,
		  Vector<T>& evals, SqMatrix<T>& evecs);

} // namespace mv

using mv::DMatrix;
using mv::CMatrix;
using mv::SqDMatrix;
using mv::SqCMatrix;
using mv::DVector;
using mv::CVector;
using mv::Vector;
using mv::Matrix;
using mv::SqMatrix;

#endif
