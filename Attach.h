#ifndef ATTACH_H
#define ATTACH_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <ctime>

using namespace std;

static const double PI = 2*acos(0.0);

// We have the following constants
static const double alpha = 4.0;
static const double nb = 3.4;
static const double gammac = 5.0*( pow(10.0, 11.0) );
static const double gamman = 2.0*( pow(10.0, 9.0) );
static const double Nbar_tr = 2.0*( pow(10.0, 24.0) );
static const double eta = pow(10.0, -19.0);
static const double GAMMA = 0.1;
static const double c = 3.0*( pow(10.0, 8.0) );
static double Delta = 0.0;
//double J = 2;

// Which leads to after dimensionalisation
static const double beta = (((2*c*GAMMA*eta)/(nb*gammac))*Nbar_tr)+1;
static const double smallgamma = gammac/(2*gamman);
static const double lambda = 1.0;


// Class B or Class C ?
static bool ClassBC = false;

// Constants for the Actual Equations
static double Lambda = 6.0;
static double Dp = 0.0;
static double Gnc = 0.1;
static double Gpc = 0.1;
static long double dt = 10.0;
static double TotalTime = 10.0;


// Time
static double TransientTime = 4.0;

//Kicks
static double InKick = 6.0;
static double FinKick = 10.0;
static double dKick = 0.1;
static int NKick = 1000;
static double NumPet = 6;
static double Perturb = 0.1;
static double IrunK = 0.0;
static double FrunK = 10.0;
static int FlNum = 1;

// Normalization
static double NormStep = 500;
static double NumOfNorm = 5000;



static double EquScal = 1;
static double KickTime = 20.0;



#include "GramSchmidt.h"
#include "RungeKutta.h"
#include "Lyapunov.h"


#endif
