#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <complex>

#include "math_fns.h"

using namespace std;

double math_fns::rng(int rngnum) {
    srand(time(0) + rngnum);
    double epsilon = (double) rand()/RAND_MAX;
    return epsilon;
}

double math_fns::norm(std::complex<double> vector[2]){
    double nor = 0;
    double norsq = 0;
    for(int i=0; i<2; i++){
        norsq += pow(abs(vector[i]),2); 
    }
    nor = sqrt(norsq);
    return nor;
}