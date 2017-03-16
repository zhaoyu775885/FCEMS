#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <mkl.h>
#include "Basics/Constants.h"

using namespace std;
using namespace Constants;

typedef complex<double> DComplex;

typedef DComplex field;

#define Em (Hm*ETA0)
#define Hm 1.0

#define TST 1


#endif
