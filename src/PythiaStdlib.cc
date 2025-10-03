// PythiaStdlib.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for Pythia utilities.

#include "Pythia8/PythiaStdlib.h"

// Directory creation for POSIX.
#include <sys/stat.h>

namespace Pythia8 {

//==========================================================================

// Convert string to lowercase for case-insensitive comparisons.
// By default remove any initial and trailing blanks or escape characters.

string toLower(const string& name, bool trim) {

  // Copy string without initial and trailing blanks or escape characters.
  string temp = name;
  if (trim) temp = trimString(name);

  // Convert to lowercase letter by letter.
  for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]);
  return temp;

}

//==========================================================================

// Return the given string with any initial and trailing blanks or
// escape characters removed.

string trimString(const string& name) {

  // Copy string without initial and trailing blanks or escape characters.
  if (name.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return "";
  int firstChar = name.find_first_not_of(" \n\t\v\b\r\f\a");
  int lastChar  = name.find_last_not_of(" \n\t\v\b\r\f\a");
  return name.substr( firstChar, lastChar + 1 - firstChar);

}

//==========================================================================

// Convert a double to string using reasonable formatting.

string toString(double val) {

  stringstream ssval;
  bool eform = false;
  if ( val == 0. ) ssval << fixed << setprecision(1);
  else if ( abs(val) < 0.001 ) ssval << scientific << setprecision(4);
  else if ( abs(val) < 0.1 ) ssval << fixed << setprecision(7);
  else if ( abs(val) < 1000. ) ssval << fixed << setprecision(5);
  else if ( abs(val) < 1000000. ) ssval << fixed << setprecision(3);
  else {eform = true; ssval << scientific << setprecision(4);}
  ssval << val;
  string sval = ssval.str();
  if (!eform) sval.erase(sval.find_last_not_of('0') + 1);
  return sval;

}

//==========================================================================

// Split a string by a delimiter.

vector<string> splitString(string val, string delim) {

  vector<string> vectorVal;
  size_t         stringPos(0);
  while (stringPos != string::npos) {
    stringPos = val.find(delim);
    vectorVal.push_back(val.substr(0, stringPos));
    val = val.substr(stringPos + 1);
  }
  return vectorVal;

}

//==========================================================================

} // end namespace Pythia8
