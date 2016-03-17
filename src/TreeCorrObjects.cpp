#include "TreeCorrObjects.h"


TreeCorrNGObject::TreeCorrNGObject(string ascii_filename) {

  /// open file
  ifstream is(ascii_filename.c_str());
  if (!is) {
    throw TreeCorrObjectsError("Specified TreeCorrNGObject filename does not exist");
  }

  /// read file one line at a time
  string buffer;
  while (getlineNoComment(is, buffer)) {
    ;
  }

  /// assign input to class variables
  ;

  return;
}
