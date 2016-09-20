#ifndef GRAMSCHMIDT_H
#include "GramSchmidt.h"
#endif

#include "Attach.h"

//Public

GramSchmidt::GramSchmidt(){
	cout << "This is the GS solver.\n";
	tol = pow (1.0, -10.0);
	
}

GramSchmidt::GramSchmidt(bool example){
	// Example as done
	// http://www.math.uconn.edu/~glaz/math2210s14/Section%20Handouts/sec6_4.pdf
	
	if(example == true)
		cout << "Here is an example. \n";
	else
		cout << "Here is an example. \n";
		
	
	Nvec = 3;
	Ndim = 4;
	tol = pow (1.0, -10.0);
	
	vector < vector <double> > sp(Nvec);
	sp[0].push_back(1.0);	sp[1].push_back(1.0); 	sp[2].push_back(1.0);
	sp[0].push_back(2.0);	sp[1].push_back(2.0);	sp[2].push_back(0.0);
	sp[0].push_back(3.0);	sp[1].push_back(0.0);	sp[2].push_back(0.0);
	sp[0].push_back(0.0);	sp[1].push_back(0.0);	sp[2].push_back(1.0);
	
	for(int i = 0; i < Ndim; i++){
		for(int j = 0; j < Nvec; j++){
			cout << " |" << sp[j][i] << "| ";
		}
		cout << endl;
	}
	
	cout << "\n\n";
	
	vector < vector <double> > sspp = GSprocess(sp);
	
	
	for(int i = 0; i < Ndim; i++){
		for(int j = 0; j < Nvec; j++){
			cout << " |" << sspp[j][i] << "| ";
		}
		cout << endl;
	}
	
	bool arewegood = check(sspp);
	
	if(arewegood == false)
		cout << " WE GOT A PROB !!! \n ";
	
}

GramSchmidt::GramSchmidt(vector< vector<double> > v){
	
	//cout << "Gram-Schmidt Orthonormalization.\n";
	
	Nvec = (int)v.size();
	Ndim = (int)v[0].size();
	tol = pow (1.0, -10.0);
	
}

GramSchmidt::GramSchmidt(vector< vector<long double> > v){
	
	//cout << "Gram-Schmidt Orthonormalization.\n";
	
	Nvec = (int)v.size();
	Ndim = (int)v[0].size();
	tol = pow (1.0, -10.0);
	
}

// Private 

double GramSchmidt::GetNorm(vector<double> v){
	
	double norm = 0.0;
	
	for(int i = 0; i < (int)v.size(); i++){
		norm = norm + v[i]*v[i];
	}
	
	norm = sqrt(norm);
	
	return norm;
}

long double GramSchmidt::GetNorm(vector<long double> v){
	
	long double norm = 0.0;
	
	for(int i = 0; i < (int)v.size(); i++){
		norm = norm + v[i]*v[i];
	}
	
	norm = sqrt(norm);
	
	return norm;
}


vector <double> GramSchmidt::projection( vector <double> v, vector <double> u){
	// Projection operator 
	// Proj_u(v) = (<v,u>/<u,u>)u
	
	vector <double> projection;
	double top = 0.0;
	double bottom = 0.0;
	for(int i = 0; i < (int)v.size(); i++){
		top = top + v[i]*u[i];
		bottom = bottom + u[i]*u[i];
	}
	
	double scale = (top/bottom);
	for(int i = 0; i < (int)v.size(); i++){
		projection.push_back(scale*u[i]);
	}
	return projection;
}

vector <long double> GramSchmidt::projection( vector <long double> v, vector <long double> u){
	// Projection operator 
	// Proj_u(v) = (<v,u>/<u,u>)u
	
	vector <long double> projection;
	long double top = 0.0;
	long double bottom = 0.0;
	for(int i = 0; i < (int)v.size(); i++){
		top = top + v[i]*u[i];
		bottom = bottom + u[i]*u[i];
	}
	
	long double scale = (top/bottom);
	for(int i = 0; i < (int)v.size(); i++){
		projection.push_back(scale*u[i]);
	}
	return projection;
}


vector < vector <double> > GramSchmidt::GSprocess(vector < vector <double> > v){
	// Using the modified Gram-Schmidt
	// https://en.wikipedia.org/wiki/Gram–Schmidt_process
	
	
	vector < vector<double> > u(Nvec);
	
	double scale = 1/GetNorm(v[0]);
	for(int i = 0; i < Ndim; i++){
		u[0].push_back(scale*v[0][i]);
	}
	
	
	vector <double> umod(Ndim);
	
	for(int i = 1; i < Nvec; i++){
		
		vector <double> tmp1 = projection(v[i],u[0]);
		for(int j = 0; j < Ndim; j++){
				umod[j] = v[i][j] - tmp1[j];
		}
		
		for(int k = 1; k < i; k++){
			vector <double> tmp2 = projection(umod, u[k]);
			for(int j = 0; j < Ndim; j++){
				umod[j] = umod[j] - tmp2[j];
			}
		}
		
		double scale = 1/GetNorm(umod);
		for(int j = 0; j < Ndim; j++){
			u[i].push_back(scale*umod[j]);
		}	
	}
	
	bool ck = check(u);
	
	if(ck == false){
		cout << "Error !!!" << endl;
	}
	
	return u;
	
}

vector < vector <long double> > GramSchmidt::GSprocess(vector < vector <long double> > v){
	// Using the modified Gram-Schmidt
	// https://en.wikipedia.org/wiki/Gram–Schmidt_process
	
	
	vector < vector<long double> > u(Nvec);
	
	long double scale = 1/GetNorm(v[0]);
	for(int i = 0; i < Ndim; i++){
		u[0].push_back(scale*v[0][i]);
	}
	
	
	vector <long double> umod(Ndim);
	
	for(int i = 1; i < Nvec; i++){
		
		vector <long double> tmp1 = projection(v[i],u[0]);
		for(int j = 0; j < Ndim; j++){
				umod[j] = v[i][j] - tmp1[j];
		}
		
		for(int k = 1; k < i; k++){
			vector <long double> tmp2 = projection(umod, u[k]);
			for(int j = 0; j < Ndim; j++){
				umod[j] = umod[j] - tmp2[j];
			}
		}
		
		long double scale = 1/GetNorm(umod);
		for(int j = 0; j < Ndim; j++){
			u[i].push_back(scale*umod[j]);
		}	
	}
	
	bool ck = check(u);
	
	if(ck == false){
		cout << "Error !!!" << endl;
	}
	
	return u;
	
}

double GramSchmidt::inprod(vector <double> a, vector <double> b){
	double tmp = 0.0;
	for(int i = 0; i < (int)a.size(); i++){
		tmp = tmp + a[i]*b[i];
	}
	return tmp;
}

long double GramSchmidt::inprod(vector <long double> a, vector <long double> b){
	long double tmp = 0.0;
	for(int i = 0; i < (int)a.size(); i++){
		tmp = tmp + a[i]*b[i];
	}
	return tmp;
}
bool GramSchmidt::check(vector< vector<double> > v){
	bool res = true;
	
	for(int i = 0; i < Nvec; i++){
		for(int j = i; j < Nvec; j++){
			if(i == j){
				double val = inprod(v[i],v[j]);
				if( !(val < 1+tol && val > 1-tol) )
					res = false;
			}
			if(i != j){
				double val = inprod(v[i],v[j]);
				if(!(val < tol && val > -tol))
					res = false;
			}
		}
	}
	return res;
}

bool GramSchmidt::check(vector< vector<long double> > v){
	bool res = true;
	
	for(int i = 0; i < Nvec; i++){
		for(int j = i; j < Nvec; j++){
			if(i == j){
				long double val = inprod(v[i],v[j]);
				if( !(val < 1+tol && val > 1-tol) )
					res = false;
			}
			if(i != j){
				long double val = inprod(v[i],v[j]);
				if(!(val < tol && val > -tol))
					res = false;
			}
		}
	}
	return res;
}
















