#ifndef LYAPUNOV_H
#include "Lyapunov.h"
#endif

#include "Attach.h"

//static inline void loadbar(unsigned int x, unsigned int n, unsigned int w );
void OutMatLy(vector <vector <double> > Mat);
void OutMatLy(vector <vector <long double> > Mat);
void OutVecLy(vector <double> Vec);
void OutVecLy(vector <long double> Vec);
void TwoVecoutfile(vector <double> vec1, vector <double> vec2 ,std::string filename);


//Public
Lyapunov::Lyapunov(){
	N = 1;
	m = 100;
	dt = 0.01;
	Tran = 1; 
	TimeEvlo = 5.0;
	tol = 0.00001;
}

Lyapunov::Lyapunov(int NormalSteps, double TimeStepSize, double TransientTime, double TimeEvolution, vector<double> y0, vector<double> p0){
	N = (int)y0.size();
	m =  NormalSteps;
	dt = TimeStepSize;
	Tran = TransientTime;
	TimeEvlo = TimeEvolution;
	tol = 0.0001;
	
	x = y0;
	p = p0;
}

Lyapunov::Lyapunov(int NormalSteps, double TimeStepSize, double TransientTime, double TimeEvolution, vector<double> y0, vector<vector<double> > p0){
	N = (int)y0.size();
	m =  NormalSteps;
	dt = TimeStepSize;
	Tran = TransientTime;
	TimeEvlo = TimeEvolution;
	tol = 0.0001;
	
	x = y0;
	pvol = p0;
}



double Lyapunov::CalcBigLypunov(vector<double> (*f_yt)(vector <double> vec, double param), double K){
	
	int xsize = (int)x.size();
	int psize = (int)p.size();
	vector<double> XP( xsize + psize );
	
	for(int i = 0; i < (int)XP.size(); i++){
		if(i < xsize)
			XP[i] = x[i];
		else
			XP[i] = p[i - xsize];
	}
	
	//OutVecLy(XP);
	
	// Iterate long enough to insure x is in attractor.
	int steps = Tran/dt;
	RungeKutta RunKut(steps, dt);
	int count = 0;
	
	
	for(int i = 0; i < steps; i++){
		XP = RunKut.RK4_11(f_yt, XP, K);
	}
	
	//Renormalize p
	for(int i = xsize; i < xsize+psize; i++){
		p[i - xsize] = XP[i];
	}
	p = normalize(p);
	for(int i = xsize; i < xsize+psize; i++){
		XP[i] = p[i - xsize];
	}
	
	
	steps = (TimeEvlo - Tran)/dt;
	count = 0;
	int NumT = 0;
	vector<double> al; 

	
	for(int i = 0; i <= steps; i++){
		XP = RunKut.RK4_11(f_yt, XP, K);
		
		if((i%m == 0 && i > 1)){
			
			for(int g = xsize; g < xsize+psize; g++){
				p[g - xsize] = XP[g];
			}
			
			NumT = m + NumT;
			count = count + 1;
			
			double tmp1 = 1/((double)(m*dt));
			double tmp2 = GetNorm(p); 
			
			
			al.push_back(tmp1*log(tmp2));
		
			p = normalize(p);
			for(int g = xsize; g < xsize+psize; g++){
				XP[g] = p[g - xsize];
			}
			
		}
		
	}
	
	double LyapExp = 0.0;
	for(int i = 0; i < (int)(al.size()); i++){
		LyapExp = (LyapExp + al[i]);
	}
	
	//LookForConvergence(m,  al);
	
	return (1/(double)(al.size() ))*LyapExp;
	
}

vector<double> Lyapunov::CalcManyLypunov(vector<double> (*f_yt)(vector <double> vec, double param), double K){
	
	int xsize = (int)x.size();
	int psize = (int)pvol.size();
	int NumP = (int)pvol[0].size();
	
	
	vector<vector<double> > XP( NumP, vector<double>(xsize + psize) );
	
	for(int j = 0; j < xsize + psize; j++){
		for(int i = 0; i < NumP; i++){
			if(j < xsize)
				XP[i][j] = x[j];
			else{
				XP[i][j] = pvol[j - xsize][i];
			}
		}	
	}
	
	
	// Iterate long enough to insure x is in attractor.
	int steps = Tran/dt;
	RungeKutta RunKut(steps, dt);
	int count = 0;
	
	for(int j = 0; j < NumP; j++){
		for(int i = 0; i < steps; i++){
			XP[j] = RunKut.RK4_11(f_yt, XP[j], 0.0);
		}
	}
	
	// Need to create a matrix to hold the ps for the G-S
	vector<vector<long double> > GSvec(NumP, vector<long double>(psize));
	
	// First fill GSvec
	for(int j = 0; j < NumP; j++){
		for(int i = 0; i < psize; i++){
			GSvec[j][i] = (long double)XP[j][i + xsize];
		}
	}
	
	
	//Execute the GS process
	GramSchmidt Orth(GSvec);
	GSvec = Orth.GSprocess(GSvec);
	
	
	// Refill XP
	for(int j = 0; j < NumP; j++){
		for(int i = 0; i < psize; i++){
			XP[j][i + xsize] = (double)GSvec[j][i];
		}
	}
	

	steps = (TimeEvlo - Tran)/dt;
	count = 0;
	int NumT = 0;
	vector<vector<double> > al(3); 
	
	for(int i = 0; i <= steps; i++){
	
		for(int j = 0; j < NumP; j++){
			XP[j] = RunKut.RK4_11(f_yt, XP[j], 0.0);
		}
		
		if((i%m == 0 && i > 1)){
			
			// First fill GSvec
			for(int h = 0; h < NumP; h++){
				for(int k = 0; k < psize; k++){
					GSvec[h][k] = (long double)XP[h][k + xsize];
				}
			}
			
			NumT = m + NumT;
			count = count + 1;
			
			double tmp1 = 1/((double)(m*dt));
			
			for(int k = 0; k < NumP; k++ ){
				double tmp2 = GetVol(GSvec, k); 
				al[k].push_back(tmp1*log(tmp2));
			}
			
			GSvec = Orth.GSprocess(GSvec);
			
			for(int h = 0; h < NumP; h++){
				for(int k = 0; k < psize; k++){
					XP[h][k + xsize] = GSvec[h][k];
				}
			}
			
		}
		
		//loadbar(i, steps, 50);
	}
	//cout << endl;
	
	
	vector<double> LyapExp(NumP);
	for(int j = 0; j < NumP; j++){
		for(int i = 0; i < (int)(al[0].size()); i++){
			if(j > 0){
				LyapExp[j] = (LyapExp[j] + al[j][i] - al[j-1][i]);
			}
			else
				LyapExp[j] = (LyapExp[j] + al[j][i]);
		}
		LyapExp[j] = (1/(double)(al[0].size() ))*LyapExp[j];
	}	
	
	//LookForConvergence(m,  al[0]);
		
	return LyapExp;
	
}

double Lyapunov::CalcBigLypunov_Kick(vector<long double> (*f_yt)(vector <long double> vec), double Kicktime, double kicksize){
	
	int xsize = (int)x.size();
	int psize = (int)p.size();
	vector<long double> XP( xsize + psize );
	vector<long double> P(psize);
	
	
	for(int i = 0; i < (int)XP.size(); i++){
		if(i < xsize)
			XP[i] = (long double)x[i];
		else
			XP[i] = (long double)p[i - xsize];
	}
	
	
	// Iterate long enough to insure x is in attractor
	int steps = Tran/dt;
	RungeKutta RunKut(steps, dt);
	int count = 0;
	
	
	for(int i = 0; i < steps; i++){
		XP = RunKut.RK4_ld(f_yt, XP);
	}
	
	//cout << "Amplitude = " << sqrt(XP[0]*XP[0] + XP[1]*XP[1]) << endl;
	//cout << "N = " << XP[4] << endl;
	
	//Renormalize p
	for(int i = xsize; i < xsize+psize; i++){
		P[i - xsize] = XP[i];
	}
	
	
	P = normalize(P);
	
	
	for(int i = xsize; i < xsize+psize; i++){
		XP[i] = P[i - xsize];
	}
	
	steps = (TimeEvlo)/dt;
	
	count = 0;
	int NumT = 0;
	vector<double> al; 
	
	int Kickstep =  (int)(Kicktime/dt);
	
	cout << "Kicktime = " << Kicktime << endl;
	cout << "Kicksize = " << kicksize << endl; 
	
	for(int i = 0; i <= steps; i++){
		XP = RunKut.RK4_ld(f_yt, XP);
		
			
		if( i%Kickstep == 0 && i > 0 ){
			//cout << "hi" << endl;
			if(xsize == 5){
				XP = KICK_C(XP, kicksize);
			}
			else{
				XP = KICK_B(XP, kicksize);
				XP = KICK_B_Pert(XP, kicksize);
			}
		}
		
		/*double phase = atan( XP[1]/XP[0] );
		double Amp = XP[0]*XP[0] + XP[1]*XP[1];
		double ActAmp = Lambda - 1 - Dp*Dp;
		
		if((phase > -0.001 && phase < 0.001 && Amp < ActAmp + 0.001 && Amp > ActAmp - 0.001 )|| i == 794328){
			cout << XP[0] << endl;
			cout << XP[1] << endl;
			cout << XP[2] << endl;
			cout << XP[3] << endl;
			cout << XP[4] << endl;
			
			cout << "Amp = " << Amp << endl;
			cout << "ActAmp = " << ActAmp << endl;
			cout << "Step = " << i << endl;
			
			
			
			cout << "Angle = " << atan( XP[1]/XP[0] ) << endl;
			cout << "Thus there is a kick of " << sin(NumPet*phase) << endl;
			
			int hihih = 0;
			cin >> hihih;
			
		}*/
		
		
		if((i%m == 0 && i > 1)){
			
			
			for(int g = xsize; g < xsize+psize; g++){
				P[g - xsize] = XP[g];
			}
			
			
			NumT = m + NumT;
			count = count + 1;
			
			double tmp1 = 1/((double)(m*dt));
			double tmp2 = GetNorm(P); 
			
			
			al.push_back(tmp1*log(tmp2));
			
			
			P = normalize(P);
			for(int g = xsize; g < xsize+psize; g++){
				XP[g] = P[g - xsize];
			}
			
			
		}
		
	}
	
	double LyapExp = 0.0;
	for(int i = 0; i < (int)(al.size()); i++){
		LyapExp = (LyapExp + al[i]);
	}
	
	LookForConvergence(m,  al);
	
	return (1/(double)(al.size() ))*LyapExp;
	
}

vector<double> Lyapunov::CalcManyLypunov_Kick(vector<double> (*f_yt)(vector <double> vec, double param), double Kicktime){
	
	int xsize = (int)x.size();
	int psize = (int)pvol.size();
	int NumP = (int)pvol[0].size();
	
	
	vector<vector<double> > XP( NumP, vector<double>(xsize + psize) );
	
	for(int j = 0; j < xsize + psize; j++){
		for(int i = 0; i < NumP; i++){
			if(j < xsize)
				XP[i][j] = x[j];
			else{
				XP[i][j] = pvol[j - xsize][i];
			}
		}	
	}
	
	
	// Iterate long enough to insure x is in attractor.
	int steps = Tran/dt;
	RungeKutta RunKut(steps, dt);
	int count = 0;
	
	for(int j = 0; j < NumP; j++){
		for(int i = 0; i < steps; i++){
			XP[j] = RunKut.RK4_11(f_yt, XP[j], 0.0);
		}
	}
	
	// Need to create a matrix to hold the ps for the G-S
	vector<vector<long double> > GSvec(NumP, vector<long double>(psize));
	
	// First fill GSvec
	for(int j = 0; j < NumP; j++){
		for(int i = 0; i < psize; i++){
			GSvec[j][i] = (long double)XP[j][i + xsize];
		}
	}
	
	
	//Execute the GS process
	GramSchmidt Orth(GSvec);
	GSvec = Orth.GSprocess(GSvec);
	
	
	// Refill XP
	for(int j = 0; j < NumP; j++){
		for(int i = 0; i < psize; i++){
			XP[j][i + xsize] = (double)GSvec[j][i];
		}
	}
	

	steps = (TimeEvlo - Tran)/dt;
	count = 0;
	int NumT = 0;
	vector<vector<double> > al(3); 
	
	for(int i = 0; i <= steps; i++){
	
		for(int j = 0; j < NumP; j++){
			XP[j] = RunKut.RK4_11(f_yt, XP[j], 0.0);
		}
		
		if((i%m == 0 && i > 1)){
			
			// First fill GSvec
			for(int h = 0; h < NumP; h++){
				for(int k = 0; k < psize; k++){
					GSvec[h][k] = (long double)XP[h][k + xsize];
				}
			}
			
			NumT = m + NumT;
			count = count + 1;
			
			double tmp1 = 1/((double)(m*dt));
			
			for(int k = 0; k < NumP; k++ ){
				double tmp2 = GetVol(GSvec, k); 
				al[k].push_back(tmp1*log(tmp2));
			}
			
			GSvec = Orth.GSprocess(GSvec);
			
			for(int h = 0; h < NumP; h++){
				for(int k = 0; k < psize; k++){
					XP[h][k + xsize] = GSvec[h][k];
				}
			}
			
		}
		
		//loadbar(i, steps, 50);
	}
	//cout << endl;
	
	
	vector<double> LyapExp(NumP);
	for(int j = 0; j < NumP; j++){
		for(int i = 0; i < (int)(al[0].size()); i++){
			if(j > 0){
				LyapExp[j] = (LyapExp[j] + al[j][i] - al[j-1][i]);
			}
			else
				LyapExp[j] = (LyapExp[j] + al[j][i]);
		}
		LyapExp[j] = (1/(double)(al[0].size() ))*LyapExp[j];
	}	
	
	//LookForConvergence(m,  al[0]);
		
	return LyapExp;
	
}
//****************************************************************************************
//Private 
vector <double> Lyapunov::normalize(vector<double> v){
	
	double norm = 0.0;
	
	for(int i = 0; i < (int)v.size(); i++){
		norm = norm + v[i]*v[i];
	}
	
	norm = sqrt(norm);
	
	vector<double> unitv;
	
	for(int i = 0; i < (int)v.size(); i++){
		unitv.push_back(v[i]/(double)norm);
	}
	
	return unitv;
}

vector <long double> Lyapunov::normalize(vector<long double> v){
	
	long double norm = 0.0;
	
	for(int i = 0; i < (int)v.size(); i++){
		norm = norm + v[i]*v[i];
	}
	
	norm = sqrt(norm);
	
	vector<long double> unitv;
	
	for(int i = 0; i < (int)v.size(); i++){
		unitv.push_back(v[i]/(long double)norm);
	}
	
	return unitv;
}

double Lyapunov::GetNorm(vector<double> v){
	
	double norm = 0.0;
	
	for(int i = 0; i < (int)v.size(); i++){
		norm = norm + v[i]*v[i];
	}
	
	norm = sqrt(norm);
	
	return norm;
}

long double Lyapunov::GetNorm(vector<long double> v){
	
	long double norm = 0.0;
	
	for(int i = 0; i < (int)v.size(); i++){
		norm = norm + v[i]*v[i];
	}
	
	norm = sqrt(norm);
	
	return (double)norm;
}

double Lyapunov::GetVol(vector<vector<long double> > v, int i){

	if(i == 0){
		return GetNorm(v[0]);
	}
	else{
		int lenvecs = v[0].size();
		int NPs = i+1;
		
		// Two matrices.
		vector<vector<long double> > Mat_1(lenvecs, vector<long double> (NPs));
		vector<vector<long double> > Mat_1_Trans(NPs, vector<long double> (lenvecs));
		
		//Fill the matrices
		for(int i = 0; i < NPs; i++){
			for(int j = 0; j < lenvecs; j++){
				Mat_1[j][i] = v[i][j];
				Mat_1_Trans[i][j] = v[i][j];
			}
		}
		
		
		//Multiply the matrices
		vector<vector<long double> > Mul_Mat(NPs, vector<long double> (NPs));
		double tmp = 0.0;
		
		for(int i = 0; i < NPs; i++){
			for(int j = 0; j < NPs; j++){
				
				for(int k = 0; k < lenvecs; k++){
					tmp = tmp + Mat_1_Trans[i][k]*Mat_1[k][j];
				}
				Mul_Mat[i][j] = tmp;
				tmp = 0.0;
	
			}
		}
		
		//Get determinant
		long double det = 0.0;
		
		// 2x2
		if(i == 1){
			
			det = (Mul_Mat[0][0]*Mul_Mat[1][1]) - (Mul_Mat[1][0]*Mul_Mat[0][1]);
			
		}
		
		// 3x3
		if(i == 2){
			
			long double tmp1 = Mul_Mat[0][0]*(Mul_Mat[1][1]*Mul_Mat[2][2] - Mul_Mat[2][1]*Mul_Mat[1][2]);
			long double tmp2 = Mul_Mat[1][0]*(Mul_Mat[0][1]*Mul_Mat[2][2] - Mul_Mat[2][1]*Mul_Mat[0][2]);
			long double tmp3 = Mul_Mat[2][0]*(Mul_Mat[0][1]*Mul_Mat[1][2] - Mul_Mat[1][1]*Mul_Mat[0][2]);
			
			det = tmp1 - tmp2 + tmp3;
		}
		
		return (double)sqrt(fabs(det));
	}
}

void Lyapunov::LookForConvergence(double m, vector<double> alphas){
	// Looks for convergence in Lyapunov exponents 
	
	int NumOfA = (int)alphas.size();
	//cout << "Size of m" << NumOfA << endl;
	vector <double> h(NumOfA);
	vector <double> M(NumOfA);
	int n = 1;
	h[0] = alphas[0];
	M[0] = m;
	for(int i = 1; i < NumOfA; i++){
		n++;
		h[i] = ( 1/(double)n )*( (n - 1)*h[i-1] + alphas[i] );
		M[i] = n*m;
	}
	TwoVecoutfile(M, h, "Check_Conver.txt");
}	

//****************************************************************************************


vector<double> Lyapunov::KICK_C(vector<double> ENvec, double kicksize){
	
	double phase = atan(ENvec[1]/ENvec[0]);
	double radius = sqrt(ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1]);
	
	double Amp = kicksize*sin(NumPet*phase);
	radius = radius + Amp;
	
	vector<double> f(10);
	
	if(ENvec[0] >= 0 && ENvec[1] >= 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
		f[6] = ENvec[6];
		f[7] = ENvec[7];
		f[8] = ENvec[8];
		f[9] = ENvec[9];
		
	}
	else if(ENvec[0] < 0 && ENvec[1] >= 0){
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
		f[6] = ENvec[6];
		f[7] = ENvec[7];
		f[8] = ENvec[8];
		f[9] = ENvec[9];
	}
	else if(ENvec[0] >= 0 && ENvec[1] < 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
		f[6] = ENvec[6];
		f[7] = ENvec[7];
		f[8] = ENvec[8];
		f[9] = ENvec[9];
	}
	else{
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
		f[6] = ENvec[6];
		f[7] = ENvec[7];
		f[8] = ENvec[8];
		f[9] = ENvec[9];
	}
	
	//cout << "phase = " << phase << "  " << ENvec[0] << "   " << ENvec[1] << endl;
	//cout << "PHASE = " << phase << "  " << f[0] << "   " << f[1] << endl;
	
	return f;	
}

vector<double> Lyapunov::KICK_B(vector<double> ENvec, double kicksize){
	
	double phase = atan(ENvec[1]/ENvec[0]);
	double radius = sqrt(ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1]);
	
	double Amp = kicksize*sin(NumPet*phase);
	radius = radius + Amp;
	
	vector<double> f(5);
	
	if(ENvec[0] >= 0 && ENvec[1] >= 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5];
	}
	else if(ENvec[0] < 0 && ENvec[1] >= 0){
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
	}
	else if(ENvec[0] >= 0 && ENvec[1] < 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
	}
	else{
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
	}
	
	//cout << "phase = " << phase << "  " << ENvec[0] << "   " << ENvec[1] << endl;
	//cout << "PHASE = " << phase << "  " << f[0] << "   " << f[1] << endl;
	
	return f;	
}

vector<long double> Lyapunov::KICK_C(vector<long double> ENvec, double kicksize){
	
	long double phase = atan(ENvec[1]/ENvec[0]);
	long double radius = sqrt(ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1]);
	
	long double Amp = kicksize*sin(NumPet*phase);
	radius = radius + Amp;
	
	vector<long double> f(10);
	
	if(ENvec[0] >= 0 && ENvec[1] >= 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
		f[6] = ENvec[6];
		f[7] = ENvec[7];
		f[8] = ENvec[8];
		f[9] = ENvec[9];
		
	}
	else if(ENvec[0] < 0 && ENvec[1] >= 0){
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
		f[6] = ENvec[6];
		f[7] = ENvec[7];
		f[8] = ENvec[8];
		f[9] = ENvec[9];
	}
	else if(ENvec[0] >= 0 && ENvec[1] < 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
		f[6] = ENvec[6];
		f[7] = ENvec[7];
		f[8] = ENvec[8];
		f[9] = ENvec[9];
	}
	else{
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
		f[6] = ENvec[6];
		f[7] = ENvec[7];
		f[8] = ENvec[8];
		f[9] = ENvec[9];
	}
	
	//cout << "phase = " << phase << "  " << ENvec[0] << "   " << ENvec[1] << endl;
	//cout << "PHASE = " << phase << "  " << f[0] << "   " << f[1] << endl;
	
	return f;	
}

vector<long double> Lyapunov::KICK_B(vector<long double> ENvec, double kicksize){
	
	long double phase = atan(ENvec[1]/ENvec[0]);
	long double radius = sqrt(ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1]);
	
	long double Amp = kicksize*sin(NumPet*phase);
	radius = radius + Amp;
	
	vector<long double> f(6);
	
	if(ENvec[0] >= 0 && ENvec[1] >= 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5];
	}
	else if(ENvec[0] < 0 && ENvec[1] >= 0){
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
	}
	else if(ENvec[0] >= 0 && ENvec[1] < 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
	}
	else{
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = ENvec[2]; 
		f[3] = ENvec[3];
		f[4] = ENvec[4];
		f[5] = ENvec[5]; 
	}
	
	//cout << "phase = " << phase << "  " << ENvec[0] << "   " << ENvec[1] << endl;
	//cout << "PHASE = " << phase << "  " << f[0] << "   " << f[1] << endl;
	
	return f;	
}

vector<long double> Lyapunov::KICK_B_Pert(vector<long double> ENvec, double kicksize){
	
	int vs = 3;
	vector<long double> PertVec(vs);
	
	for(int i = 0; i < vs; i++){
		PertVec[i] = ENvec[i] + ENvec[i + vs];
	}
	
	long double phase = atan(PertVec[1]/PertVec[0]);
	long double radius = sqrt(PertVec[0]*PertVec[0] + PertVec[1]*PertVec[1]);
	
	long double Amp = kicksize*sin(NumPet*phase);
	radius = radius + Amp;
	
	vector<long double> f(vs);
	
	if(PertVec[0] >= 0 && PertVec[1] >= 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = PertVec[2]; 
	}
	else if(PertVec[0] < 0 && PertVec[1] >= 0){
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = fabs(radius*sin(phase));
		f[2] = PertVec[2]; 
	}
	else if(PertVec[0] >= 0 && PertVec[1] < 0){
		f[0] = fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = PertVec[2]; 
	}
	else{
		f[0] = -1.0*fabs(radius*cos(phase));
		f[1] = -1.0*fabs(radius*sin(phase));
		f[2] = PertVec[2]; 
	}
	
	for(int i = 0; i < vs; i++){
		ENvec[i + vs] = f[i] - ENvec[i];
	}
	
	
	return ENvec;	
}

//****************************************************************************************
/*static inline void loadbar(unsigned int x, unsigned int n, unsigned int w)
{
    if ( (x != n) && (x % (n/100+1) != 0) ) return;
 
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;
 
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}*/

void OutMatLy(vector <vector <double> > Mat){
	cout << setprecision(60);
	for(int i = 0; i < Mat.size(); i++){
		for(int j = 0; j < Mat[0].size(); j++){
			cout << "  |  " << Mat[i][j];
		}
		cout << "\n" << endl;
	}
}

void OutMatLy(vector <vector <long double> > Mat){
	for(int i = 0; i < Mat.size(); i++){
		for(int j = 0; j < Mat[0].size(); j++){
			cout << "  |  " << Mat[i][j];
		}
		cout << "\n" << endl;
	}
}

void OutVecLy(vector <double> Vec){
		for(int j = 0; j < Vec.size(); j++){
			cout << Vec[j] << endl;
		}
		cout << endl;
}

void OutVecLy(vector <long double> Vec){
		for(int j = 0; j < Vec.size(); j++){
			cout << Vec[j] << endl;
		}
		cout << endl;
}

void TwoVecoutfile(vector <double> vec1, vector <double> vec2 ,std::string filename){
	
	std::ofstream ofile(filename, std::ios::out);
	
	ofile.precision(15);
	
	for(int i = 0; i < (int)vec1.size(); i++)
	{
		ofile << vec1[i] << "\t" << vec2[i] << "\n";
	}
	ofile.close();
}









