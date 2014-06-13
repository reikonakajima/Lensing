//---------------------------------------------------------------------------
#ifndef StdH
#define StdH
//---------------------------------------------------------------------------
// 	$Id: Std.h,v 1.11 2003-03-06 03:29:16 garyb Exp $

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <complex>
#include <vector>
using std::vector;
#include <cstdlib>
#include <cmath>
#include <new>

using std::string;
using std::ofstream;
using std::ifstream;
using std::ostream;
using std::istream;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

// ???using namespace std;

// Useful conventions:
template <class T>
void SWAP(T& a,T& b) {T temp=a; a=b; b=temp;}

template <class T>
T SQR(const T& x) {return x*x;}

template <class T>
const T& MAX(const T& a, const T& b) {return a>b ? a : b;}

template <class T>
const T& MIN(const T& a, const T& b) {return a<b ? a : b;}

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long ulong;
typedef std::complex<double> DComplex;


// ??? replace here with numeric_limits<double>.infinity()
#ifdef HUGE_VAL
const double DOUBLE_INFINITY=HUGE_VAL;
#else
const double DOUBLE_INFINITY=1./0.;
#endif
const double DOUBLE_NEGATIVE_INFINITY=-DOUBLE_INFINITY;
const double DOUBLE_NAN = DOUBLE_NEGATIVE_INFINITY + DOUBLE_INFINITY;
#ifndef PI
#ifdef M_PI
#define PI M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif

// Debugging and exception classes:

extern ostream* dbgout;
#ifdef DEBUGLOGGING
#define dbg if(dbgout) (*dbgout)
#else
#define dbg if(false) (cerr)
#endif

#ifdef ASSERT
  #define Assert(x) \
    { if(!(x)) { \
      cerr << "Error - Assert " #x " failed"<<endl; \
      cerr << "on line "<<__LINE__<<" in file "<<__FILE__<<endl; \
      exit(1);} }
#else
  #define Assert(x)
#endif

inline void error(const char *s,const char *s2 = "")
{
  dbg << "Error: " << s << ' ' << s2 << endl;
  cerr << "Error: " << s << ' ' << s2 << endl;
  exit(1);
}

// Define a base exception class with error message:
class MyException {
 public:
  MyException(const string &m=""): msg(m) {};
  MyException(const char* c): msg(c) {};
  void set(const string &m="") {msg=m;}
  void set(const char* c) {msg=c;}
  void dump(ostream &os) const {os << "Error: " << msg << endl;}
  void quit(const int exit_code=1) {
    // dbg << "Error: " << msg << endl;
    cerr << "Error: " << msg << endl;
    exit(exit_code);
  }
  string msg;
};

// Dump a message in case of memory exhaustion
// main needs to say
// set_new_handler(out_of_memory);
// to use this.
/*void out_of_memory() {
  cerr << "Memory allocation failure." << endl;
  exit(1);
  }*/

// A vector class with compile-optional bounds checking

#ifdef DEBUG
template <class T> 
class myvector : virtual public std::vector<T> {
public:
  myvector() : vector<T>() {}
  explicit myvector(uint n) : vector<T>(n) {}
  myvector(const vector<T>& v) : vector<T>(v) {}
  myvector(uint n, const T& t) : vector<T>(n,t) {}
  template <class InputIterator>
  myvector(InputIterator i1,InputIterator i2) : vector<T>(i1,i2) {}
  virtual ~myvector() {}
  virtual T& operator[](uint n)
  { if(n>=size()) {
      dbg << "bad myvector[] call: n = "<<n<<", size = "<<size()<<endl;
      error("bad myvector[] call");
    }
    else return vector<T>::operator[](n); }
  virtual const T& operator[](uint n) const
    { if(n>=size()) {
      dbg << "bad myvector[] call: n = "<<n<<", size = "<<size()<<endl;
      error("bad myvector[] call");
    }
    return vector<T>::operator[](n); }
};

template<> class myvector<bool> : public std::vector<bool> {
public:
  myvector() : vector<bool>() {}
  explicit myvector(uint n) : vector<bool>(n) {}
  myvector(const vector<bool>& v) : vector<bool>(v) {}
  myvector(uint n, const bool& t) : vector<bool>(n,t) {}
  template <class InputIterator>
    myvector(InputIterator i1,InputIterator i2) : vector<bool>(i1,i2) {}
  ~myvector() {}
};

#else
#define myvector std::vector
#endif

#endif
