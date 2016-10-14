#ifndef LYAPUNOV_H
#define LYAPUNOV_H

#include "Attach.h"

// This class (will) contain a functional Gram-Schmidt 
// orthogonalization procedure. 

class Lyapunov{
public:
	Lyapunov(); // default constructor
	Lyapunov(int NormalSteps, double TimeStepSize, double TransientTime, double TimeEvolution, vector<double> y0, vector<double> p0);
	Lyapunov(int NormalSteps, double TimeStepSize, double TransientTime, double TimeEvolution, vector<long double> y0, vector<long double> p0);
	
	Lyapunov(int NormalSteps, double TimeStepSize, double TransientTime, double TimeEvolution, vector<double> y0, vector<vector<double> > p0);
	
	vector<double> CalcManyLypunov(vector<double> (*f_yt)(vector <double> vec, double param), double K);
	
	double CalcBigLypunov_Kick_new(vector<double> (*f_yt)(vector<double> vec), vector<double> (*k_yt)(vector<double> vec, double KickSize), double Kicktime, double kicksize);
	
	double CalcBigLypunov_Kick_new_l(vector<long double> (*f_yt)(vector<long double> vec), vector<long double> (*k_yt)(vector<long double> vec, double KickSize), double Kicktime, double kicksize);
	
	
private:
	vector<double> x;
	vector<double> p;
	
	vector<long double> xl;
	vector<long double> pl;
	
	vector<vector<double> > pvol;
	int N; // Dim of system
	int m; // Number of steps before renormalization
	double dt; // Time Steps
	double Tran; // Time to settle into attractor
	double TimeEvlo; // Time to settle into attractor
	double tol;
	
	vector<double> KICK_B(vector<double> ENvec, double kicksize);
	vector<double> KICK_C(vector<double> ENvec, double kicksize);
	
	vector<long double> KICK_B(vector<long double> ENvec, double kicksize);
	vector<long double> KICK_C(vector<long double> ENvec, double kicksize);
	
	vector<long double> KICK_B_Pert(vector<long double> ENvec, double kicksize);
	
	vector <double> normalize(vector<double> v);
	vector <long double> normalize(vector<long double> v);
	
	vector <double> normalize_P_only(vector<double> p);
	
	double GetNorm(vector<double> v);
	long double GetNorm(vector<long double> v);
	double GetVol(vector<vector<long double> > v, int i);
	void LookForConvergence(double m, vector<double> alphas);
};

#endif