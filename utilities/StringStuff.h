// $Id: StringStuff.h,v 2.4 2003-02-18 15:18:17 garyb Exp $
// Common convenience functions with strings
#ifndef STRINGSTUFF_H
#define STRINGSTUFF_H

#include <iostream>
#include <sstream>
#include <string>
using std::string;

namespace stringstuff {
  extern bool isComment(const string& instr);

  extern std::istream& getlineNoComment(std::istream& is, string& s);

  bool nocaseEqual(char c1, char c2);
  bool nocaseEqual(const string& s1, const string& s2);

  // Remove trailing whitespace
  void stripTrailingBlanks(string& s);
  // Remove everything after & including the last "." period
  void stripExtension(string& s);

  // Make a string that holds the current time and the command line
  string taggedCommandLine(int argc, char *argv[]);
} // namespace stringstuff

using namespace stringstuff;

#endif
