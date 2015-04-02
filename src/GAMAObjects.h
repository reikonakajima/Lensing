//
// GAMAObjects.h
//
#ifndef GAMAOBJECTS_H
#define GAMAOBJECTS_H
#include "LensObjects.h"

class GAMAObjectsError : public MyException {
 public:
  GAMAObjectsError(const string& m="") :
  MyException("GAMAObjectsError: " +m) {}
};


/*
 * class GAMAObject
 *
 * This class reads the GAMA data file with galaxy positions as well as halo center (which may be
 * offset from the star position), and sets it up as a "lens" in galaxy-galaxy lensing format.
 *
 * The columns are formatted:
 *   0 -  4 CATAID, GROUPIDA, RA, DEC, Z,
 *   5 -  9 ABSMAG_R, DELABSMAG_R, LOGMSTAR, DELLOGMSTAR, UMINUSR,
 *  10 - 14 DELUMINUSR, RPETRO, LOGMOVERL, DELLOGMOVERL, LOGLWAGE,
 *  15 - 19 DELLOGLWAGE, METAL, DELMETAL, LOGTAU, DELLOGTAU,
 *  20 - 24 LOGMREMNANTS, DELLOGMREMNATS, RANKBCG, NFOF, ZMAX_19P8,
 *  25      ZMAX_19P4
 *
 */

class GAMAObject : public LensObject {

 public:
  GAMAObject(const string buffer);
  float getAbsMagR() const { return absmag_r; }
  float getAbsMagRErr() const { return d_absmag_r; }
  float getLogMStar() const { return logmstar; }
  float getLogMStarErr() const { return d_logmstar; }
  float getRankBCG() const { return rankbcg; }
  int   getNfof() const { return Nfof; }
  void printLine(ostream& os) const;

 private:
  int   group_id;
  float absmag_r, d_absmag_r;
  float logmstar, d_logmstar;
  float uminusr, d_uminusr;
  float rpetro;
  float logmoverl, d_logmoverl;
  float loglwage, d_loglwage;
  float metal, d_metal;
  float logtau, d_logtau;
  float logmremnants, d_logmremnants;
  float rankbcg;
  int   Nfof;
  float Zmax_19P8, Zmax_19P4;
};


class GAMAObjectList : public LensObjectList<GAMAObject*> {

 public:
  GAMAObjectList(istream& is);
  int applyLogMStarCut(float min_logmstar, float max_logmstar);

 private:
  // none

};


#endif // GAMAOBJECTS_H
