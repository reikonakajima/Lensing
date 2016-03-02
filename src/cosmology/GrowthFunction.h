// $Id: GrowthFunction.h,v 1.4 2006-07-07 18:26:45 garyb Exp $
// Solve the differential equation for GR growth function
// given a cosmology.
#ifndef GROWTHFUNCTION_H
#define GROWTHFUNCTION_H

#include "Cosmology.h"
#include "GTable.h"

namespace cosmology {

  class GrowthFunction {
  public:
    GrowthFunction(const Cosmology& c_): c(c_), 
      lnGTable(gtable::Table<>::spline),
      dlnGdlnATable(gtable::Table<>::spline) {
      buildTables();
    }
    GrowthFunction(const GrowthFunction& rhs): c(rhs.c), 
      lnGTable(gtable::Table<>::spline),
      dlnGdlnATable(gtable::Table<>::spline) {
      buildTables();  // Make clean copies of tables in a copy.
    }
    double operator()(double z) const {return exp(lnG(-log(1.+z)));}
    double valueAtZ(double z) const;
    double lnG(double lna) const;
    double dlnGdlnA(double lna) const;
  private:
    void operator=(const GrowthFunction& rhs) {} //Hide op=
    const Cosmology& c;
    // The diff eq is solved once at construction, then all values
    // looked up in spline-fit tables.
    gtable::Table<> lnGTable;
    gtable::Table<> dlnGdlnATable;
    void buildTables(double aMax=1.1);
    double lnGRecombination;
  };

}
#endif // GROWTHFUNCTION_H
