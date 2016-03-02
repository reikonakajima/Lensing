// $Id: valarrayKludge.h,v 1.1 2002-12-02 14:23:48 garyb Exp $
// This file is included by Matrix.h instead of <valarray> because
// efficient optimization on compilers requires use of pointer arithmetic
// instead of valarray slices.  So we need a version of valarray that has
// const T& operator[] const   instead of
// T operator[] const.
// *** Note that there better not have been an #include <valarray> before
// *** this!

#ifdef ICC
// Intel compiler:
#include "valarray.icc"
#define HAVE_VALARRAY_POINTERS

#elif  GCC295
// GCC 2.95
#include "valarray.gcc295"
#define HAVE_VALARRAY_POINTERS

#elif  GCC3
// GCC 3.x already does it our way:
#include <valarray>
#define HAVE_VALARRAY_POINTERS

#else
// Something I don't know about
#include <valarray>
#endif
