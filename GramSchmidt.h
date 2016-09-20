#ifndef GRAMSCHMIDT_H
#define GRAMSCHMIDT_H

#include "Attach.h"

// This class (will) contain a functional Gram-Schmidt 
// orthogonalization procedure. 

class GramSchmidt{
public:
	GramSchmidt(); // default constructor
	GramSchmidt(bool example);
	GramSchmidt(vector< vector<double> > v);
	GramSchmidt(vector< vector<long double> > v);
	
	vector < vector <double> > GSprocess(vector < vector <double> >);
	vector < vector <long double> > GSprocess(vector < vector <long double> >);

private:
	double GetNorm(vector <double>);
	vector <double> projection( vector <double>, vector <double>);
	double inprod(vector <double> a, vector <double> b);
	bool check(vector< vector<double> > v);
	
	long double GetNorm(vector <long double>);
	vector <long double> projection( vector <long double>, vector <long double>);
	long double inprod(vector <long double> a, vector <long double> b);
	bool check(vector< vector<long double> > v);
	
	int Nvec; // The number of vectors to be normalized
	int Ndim; // The dimension of the vectors
	double tol;
	
};

#endif