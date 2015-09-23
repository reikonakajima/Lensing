//
// GAMAObjects.h
//
#ifndef GAMAOBJECTS_H
#define GAMAOBJECTS_H
#include "LensObjects.h"
#include <CCfits/CCfits>

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
 * The (Edo-generated) ASCII columns are formatted:
 *   0 -  4 CATAID, GROUPIDA, RA, DEC, Z,
 *   5 -  9 ABSMAG_R, DELABSMAG_R, LOGMSTAR, DELLOGMSTAR, UMINUSR,
 *  10 - 14 DELUMINUSR, RPETRO, LOGMOVERL, DELLOGMOVERL, LOGLWAGE,
 *  15 - 19 DELLOGLWAGE, METAL, DELMETAL, LOGTAU, DELLOGTAU,
 *  20 - 24 LOGMREMNANTS, DELLOGMREMNATS, RANKBCG, NFOF, ZMAX_19P8,
 *  25      ZMAX_19P4
 *
 * The FITS table formats have various columns, including
    name = 'CATAID_1'; format = 'J'
    name = 'RA'; format = 'D'; unit = 'deg'
    name = 'DEC'; format = 'D'; unit = 'deg'
    name = 'Z_1'; format = 'E'
    name = 'logmstar'; format = 'E'; unit = 'dex'
    name = 'uminusr'; format = 'E'; unit = 'mag'
    name = 'absmag_r'; format = 'E'; unit = 'mag'
    name = 'R_COMOVING'; format = '1E'; unit = 'Mpc(Planck14)'
 *
 *
 *  Might be nice to retain the following quantities from the original:
 *
    name = 'Rpetro'; format = 'D'; unit = 'mag'
    name = 'RankBCG'; format = 'I'
    name = 'SepBCG'; format = 'D'; unit = 'arcsec'
    name = 'CoSepBCG'; format = 'D'; unit = 'Mpc/h'
    name = 'AngSepBCG'; format = 'D'; unit = 'Mpc/h'
    name = 'RankCen'; format = 'I'
    name = 'SepCen'; format = 'D'; unit = 'arcsec'
    name = 'CoSepCen'; format = 'D'; unit = 'Mpc/h'
    name = 'AngSepCen'; format = 'D'; unit = 'Mpc/h'
 *
 *
 *
 */

class GAMAObject : public LensObject {

 public:
  GAMAObject(const string buffer);
  GAMAObject(long int _id, double _ra, double _dec, float _z, float _r_comoving,
	     float _logmstar, float _absmag_r, float _uminusr) :
             LensObject(_id, _ra, _dec, _z, -99), r_comoving(_r_comoving), logmstar(_logmstar),
	     absmag_r(_absmag_r), uminusr(_uminusr) {}
  float getAbsMagR() const { return absmag_r; }
  //float getAbsMagRErr() const { return d_absmag_r; }
  float getLogMStar() const { return logmstar; }
  //float getLogMStarErr() const { return d_logmstar; }
  //float getRankBCG() const { return rankbcg; }
  //int   getNfof() const { return Nfof; }
  void printLine(ostream& os) const;

 private:
  int   group_id;
  float r_comoving;
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
  // Constructor:  If max_count <= 0, use the entire input.
  //               Otherwise, randomize input rows and return the first max_count objects
  GAMAObjectList(const string fits_filename, int max_count=-1);
  int applyLogMStarCut(float min_logmstar, float max_logmstar);

 private:
  // none

};


#endif // GAMAOBJECTS_H
