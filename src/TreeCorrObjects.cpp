#include "TreeCorrObjects.h"


TreeCorrNGObject::TreeCorrNGObject(string ascii_filename) {

  /// open file
  ifstream is(ascii_filename.c_str());
  if (!is) {
    throw TreeCorrObjectsError("Specified TreeCorrNGObject filename does not exist");
  }

  /// clear database
  R_nom.clear();
  meanR.clear();
  gamT.clear();
  gamX.clear();
  sigma.clear();
  weight.clear();
  npairs.clear();

  /// read file one line at a time
  string buffer;
  while (getlineNoComment(is, buffer)) {
    std::istringstream iss(buffer);
    /// #  R_nom  meanR  meanlogR  gamT  gamX  sigma  weight  npairs
    float f1, f2, f3, f4, f5, f6, f7, f8;
    if (!(iss >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8)) {
      cerr << "##" << buffer << endl;
      throw TreeCorrObjectsError("error reading TreeCorrNGObject file");
    }
    /// assign input to class variables
    R_nom.push_back(f1);
    meanR.push_back(f2);
    gamT.push_back(f4);
    gamX.push_back(f5);
    sigma.push_back(f6);
    weight.push_back(f7);
    npairs.push_back(f8);
  }

  return;
}
