#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "Attach.h"


class RungeKutta{
public:
	RungeKutta(); // default constructor
	RungeKutta(int NumSteps, double dt);
	RungeKutta(double TimeEvlo, double dt);
	
	vector <vector <double> > RK4(vector<double> (f_yt)(vector <double> vec, double param), vector<double> y0, double K);	
	vector <vector <double> > RK4(vector<double> (f_yt)(vector <double> vec), vector<double> y0);
	
	vector<double> RK4_11(vector <double> (f_yt)(vector<double> vec, double param), vector<double> ENPvec, double K);
	
	vector <double> RK4_11(vector <double> (f_yt)(vector<double> vec), vector<double> ENPvec);
	
	vector<long double> RK4_ld(vector <long double> (f_yt)(vector<long double> vec), vector<long double> ENPvec);
	
private:
	double h; // Step size
	int N; // Number of steps  

	//Slopes
	vector <double> k1;
	vector <double> k2;
	vector <double> k3;
	vector <double> k4;
	vector <double> k5;
	vector <double> k6;
	
	
	vector <long double> k1_l;
	vector <long double> k2_l;
	vector <long double> k3_l;
	vector <long double> k4_l;
	
	
	
	void FillWithVec(vector <double> v1, vector <double> v2);
	
};

#endif