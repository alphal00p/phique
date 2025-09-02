#include <iostream>
#include "polygamma.hpp"

using namespace std;

extern "C" {
  // It is important that the function or subroutine name should be in lower cases
  // for Fortran
  complex<double> polygamma_(const int & n, const complex<double> * z) {
    // call PolyGamma[n,z] function
    return pGamma(n,*z);
  }
}
