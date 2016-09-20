#include "Attach.h"


vector <double> LaserEquation(vector <double> ENvec, double K);
vector <double> JacobianTimesPLaserEquation(vector <double> ENvec, vector <double> p);
vector <double> LaserWithJacobian(vector <double> ENPvec, double K);
vector <double> RosslerEquation(vector <double> ENvec, double K);
void OutMat(vector <vector <double> > Mat);
void OutVec(vector <double> Vec);
void Matoutfile(vector <vector <double> > Mat ,std::string filename);
void Vecoutfile(vector <double> vec ,std::string filename);
void MatToVecOutFiles(vector <vector <double> > Mat ,std::string filename);
static inline void loadbarMain(unsigned int x, unsigned int n, unsigned int w);
static inline void loadbarMain(unsigned int x, unsigned int n, unsigned int w, clock_t start);
vector <double> normalize(vector<double> v);


void SetParameters(string filename);
void ExecuteScaling();

vector<double> ClassCPlaner(vector<double> ENvec, double k);
vector<double> ClassCPlanerWithJacobian(vector<double> ENvec, double k);
vector<double> ClassBPlaner(vector<double> ENvec, double k);
vector<double> ClassBPlanerWithJacobian(vector<double> ENvec, double k);
vector<long double> ClassBPlanerWithJacobian_ld(vector<long double> ENvec);
vector<long double> ClassCPlanerWithJacobian_ld(vector<long double> ENvec);


int main(){

	cout << "Start the clock !!!\n";
	clock_t t = clock();
		
	vector <double> y0(3);
	y0[0] = 0.1;
	y0[1] = 2.1;
	y0[2] = 5.1;
	
	SetParameters("Parameters.txt");
	cout << " dt = " << dt << endl;
	ExecuteScaling();
	cout << " dt = " << dt << endl;
	
//****************************************************************************************

	
	{
		int NumK = (int)((FinKick - InKick)/dKick);
		double TimeEvlo = ((double)NumOfNorm*NormStep*dt);
		
		cout << "Total Number of steps:  " << NumOfNorm*NormStep << endl;

		int vecsize = 3;			
		vector <double> p0(vecsize);				vector <double> y0(vecsize);
		p0[0] = -0.00001;						y0[0] = -2.1;
		p0[1] = -0.00001;						y0[1] = 2.1;
		p0[2] = -0.00001;						y0[2] = 5.1;
//		p0[3] = 0.001;						y0[3] = 2.0;
//		p0[4] = 0.001;						y0[4] = 5.1;
		 
		 
		 
		vector<double> h;
		vector<double> k;
			
		clock_t st = clock();
		for(int i = 0; i < NumK; i++){
			Lyapunov LyapLase(NormStep, dt, TransientTime, TimeEvlo, y0, p0);
			h.push_back(LyapLase.CalcBigLypunov_Kick(ClassBPlanerWithJacobian_ld, InKick + i*dKick, Perturb ));
			k.push_back(InKick + i*dKick);
			
			loadbarMain(i, NumK, 50, st);
			
		}
		
			
		Vecoutfile(h ,"LyapunovSlice.txt");
		Vecoutfile(k ,"Kvals.txt");
	}

//****************************************************************************************
	
	cout << "Stop the clock !!!" << endl;
	t = clock() - t;
	
	cout << " It took me " << t << " clicks ( " << ((float)t)/CLOCKS_PER_SEC << " seconds). " << endl;
	cout << "You are the weakest link—goodbye!" << endl;
	
	cout << "Enter 0 to exit . . . ";
	double kjkjkj = 0;
	cin >> kjkjkj;
	return 0;
}


//****************************************************************************************

vector<double> ClassCPlaner(vector<double> ENvec, double k){
	
	vector<double> f(5);
	k = k+1;
	
	f[0] = - ENvec[0] - Dp*ENvec[1] - ENvec[3];
	f[1] = - ENvec[1] + Dp*ENvec[0] + ENvec[2];
	f[2] = Gpc*( -1*ENvec[2] + Dp*ENvec[3] + ENvec[1]*ENvec[4] );
	f[3] = Gpc*( -1*ENvec[3] - Dp*ENvec[2] - ENvec[0]*ENvec[4] );
	f[4] = Gnc*( Lambda - ENvec[4] + ENvec[3]*ENvec[0] - ENvec[2]*ENvec[1] );
	
	return f;
	
}

vector<double> ClassCPlanerWithJacobian(vector<double> ENvec, double k){
	
	vector<double> f(10);
	k = k+1;
	
	f[0] = - ENvec[0] - Dp*ENvec[1] - ENvec[3];
	f[1] = - ENvec[1] + Dp*ENvec[0] + ENvec[2];
	f[2] = Gpc*( -1*ENvec[2] + Dp*ENvec[3] + ENvec[1]*ENvec[4] );
	f[3] = Gpc*( -1*ENvec[3] - Dp*ENvec[2] - ENvec[0]*ENvec[4] );
	f[4] = Gnc*( Lambda - ENvec[4] + ENvec[3]*ENvec[0] - ENvec[2]*ENvec[1] );
	
	// Jacobian
	
	f[5] = - ENvec[5] - Dp*ENvec[6] - ENvec[8];
	f[6] = Dp*ENvec[5] - ENvec[6] + ENvec[7];
	
	f[7] = Gpc*(   ENvec[4]*ENvec[6] - ENvec[7] + Dp*ENvec[8] + ENvec[1]*ENvec[9] );
	f[8] = Gpc*( - ENvec[4]*ENvec[5] - Dp*ENvec[7] - ENvec[8] - ENvec[0]*ENvec[9] );
	
	f[9] = Gnc*(   ENvec[3]*ENvec[5] - ENvec[2]*ENvec[6] - ENvec[1]*ENvec[7] +
				   ENvec[0]*ENvec[8] - ENvec[9] );
	
	f[0] = EquScal*f[0];
	f[1] = EquScal*f[1];
	f[2] = EquScal*f[2];
	f[3] = EquScal*f[3];
	f[4] = EquScal*f[4];
	
	f[5] = EquScal*f[5];
	f[6] = EquScal*f[6];
	f[7] = EquScal*f[7];
	f[8] = EquScal*f[8];
	f[9] = EquScal*f[9];
	
	
	return f;
	
}

vector<double> ClassBPlaner(vector<double> ENvec, double k){
	
	vector<double> f(3);
	
	f[0] = ( (ENvec[2]/(1 + Dp*Dp)) - 1)*(ENvec[0] + Dp*ENvec[1]) + k;
	f[1] = ( (ENvec[2]/(1 + Dp*Dp)) - 1)*(ENvec[1] - Dp*ENvec[0]);
	f[2] = Gnc*( Lambda - ENvec[2] - ((ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1])*ENvec[2])/(1 + Dp*Dp));
	
	return f;
	
}

vector<double> ClassBPlanerWithJacobian(vector<double> ENvec, double k){
	
	vector<double> f(6);
	
	f[0] = ( (ENvec[2]/(1 + Dp*Dp)) - 1.0)*(ENvec[0] + Dp*ENvec[1]);
	f[1] = ( (ENvec[2]/(1 + Dp*Dp)) - 1.0)*(ENvec[1] - Dp*ENvec[0]);
	f[2] = Gnc*( Lambda - ENvec[2] - ((ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1])*ENvec[2])/(1.0 + Dp*Dp));
	
	k = k+1;
	
	// Just for ease lets define some simple constants
	double alf = (ENvec[2]/(1.0 + Dp*Dp)) - 1.0;
	double Lor = (1.0/(1.0 + Dp*Dp));
	double EEE = ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1];
	
	f[3] = alf*ENvec[3] + Dp*alf*ENvec[4] + (ENvec[0] + Dp*ENvec[1])*Lor*ENvec[5];
	f[4] = -Dp*alf*ENvec[3] + alf*ENvec[4] + (ENvec[1] - Dp*ENvec[0])*Lor*ENvec[5];
	f[5] = -Gnc*( 2.0*ENvec[0]*ENvec[2]*Lor*ENvec[3] + 2.0*ENvec[1]*ENvec[2]*Lor*ENvec[4] + ENvec[5] + EEE*Lor*ENvec[5] );
	
	
	
	f[0] = EquScal*f[0];
	f[1] = EquScal*f[1];
	f[2] = EquScal*f[2];
	
	f[3] = EquScal*f[3];
	f[4] = EquScal*f[4];
	f[5] = EquScal*f[5];
	
	return f;
	
}


vector<long double> ClassCPlanerWithJacobian_ld(vector<long double> ENvec){
	
	vector<long double> f(10);
	
	f[0] = - ENvec[0] - Dp*ENvec[1] - ENvec[3];
	f[1] = - ENvec[1] + Dp*ENvec[0] + ENvec[2];
	f[2] = Gpc*( -1*ENvec[2] + Dp*ENvec[3] + ENvec[1]*ENvec[4] );
	f[3] = Gpc*( -1*ENvec[3] - Dp*ENvec[2] - ENvec[0]*ENvec[4] );
	f[4] = Gnc*( Lambda - ENvec[4] + ENvec[3]*ENvec[0] - ENvec[2]*ENvec[1] );
	
	// Jacobian
	
	f[5] = - ENvec[5] - Dp*ENvec[6] - ENvec[8];
	f[6] = Dp*ENvec[5] - ENvec[6] + ENvec[7];
	
	f[7] = Gpc*(   ENvec[4]*ENvec[6] - ENvec[7] + Dp*ENvec[8] + ENvec[1]*ENvec[9] );
	f[8] = Gpc*( - ENvec[4]*ENvec[5] - Dp*ENvec[7] - ENvec[8] - ENvec[0]*ENvec[9] );
	
	f[9] = Gnc*(   ENvec[3]*ENvec[5] - ENvec[2]*ENvec[6] - ENvec[1]*ENvec[7] +
				   ENvec[0]*ENvec[8] - ENvec[9] );
	
	f[0] = EquScal*f[0];
	f[1] = EquScal*f[1];
	f[2] = EquScal*f[2];
	f[3] = EquScal*f[3];
	f[4] = EquScal*f[4];
	
	f[5] = EquScal*f[5];
	f[6] = EquScal*f[6];
	f[7] = EquScal*f[7];
	f[8] = EquScal*f[8];
	f[9] = EquScal*f[9];
	
	
	return f;
	
}

vector<long double> ClassBPlanerWithJacobian_ld(vector<long double> ENvec){
	
	vector<long double> f(6);
	
	f[0] = ( (ENvec[2]/(1.0 + Dp*Dp)) - 1.0)*(ENvec[0] + Dp*ENvec[1]);
	f[1] = ( (ENvec[2]/(1.0 + Dp*Dp)) - 1.0)*(ENvec[1] - Dp*ENvec[0]);
	f[2] = Gnc*( Lambda - ENvec[2] - ((ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1])*ENvec[2])/(1.0 + Dp*Dp));
	
	
	// Just for ease lets define some simple constants
	double alf = (ENvec[2]/(1 + Dp*Dp)) - 1.0;
	double Lor = (1.0/(1.0 + Dp*Dp));
	double EEE = ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1];
	
	f[3] = alf*ENvec[3] + Dp*alf*ENvec[4] + (ENvec[0] + Dp*ENvec[1])*Lor*ENvec[5];
	f[4] = -Dp*alf*ENvec[3] + alf*ENvec[4] + (ENvec[1] - Dp*ENvec[0])*Lor*ENvec[5];
	f[5] = -Gnc*( 2.0*ENvec[0]*ENvec[2]*Lor*ENvec[3] + 2.0*ENvec[1]*ENvec[2]*Lor*ENvec[4] + ENvec[5] + EEE*Lor*ENvec[5] );
	
	
	f[0] = EquScal*f[0];
	f[1] = EquScal*f[1];
	f[2] = EquScal*f[2];
	
	f[3] = EquScal*f[3];
	f[4] = EquScal*f[4];
	f[5] = EquScal*f[5];
	
	return f;
	
}

//****************************************************************************************

void SetParameters(string filename){
	
	fstream myfile(filename, std::ios_base::in);
	
    myfile >> Lambda >> Dp >> Gnc >> Gpc >> TransientTime >> dt >> InKick >> FinKick 
    	   >> dKick >> Perturb >> NumPet >> NormStep >> NumOfNorm >> EquScal;
    myfile.close();
    
    cout << "     **** Parameters **** " << endl;
    cout << "Lambda               = " << Lambda << endl;
    cout << "Dp                   = " << Dp << endl;
    cout << "Gnc                  = " << Gnc << endl;
    cout << "Gpc                  = " << Gpc << endl;
    cout << "--------------------------------" << endl;
    cout << "Transient Time       = " << TransientTime << endl;
    cout << "dt                   = " << dt << endl;
    cout << "--------------------------------" << endl;
    cout << "Initial Kick         = " << InKick << endl;
    cout << "Final Kick           = " << FinKick << endl;
    cout << "delta Kick           = " << dKick << endl;
    cout << "Kick Size            = " << Perturb << endl;
    cout << "Number of petals     = " << NumPet << endl;
    cout << "--------------------------------" << endl;
    cout << "Steps Till Normalize = " << NormStep << endl;
    cout << "Num of Normalize     = " << NumOfNorm << endl;
    cout << "--------------------------------" << endl;
    cout << "Equation Scaling     = " << EquScal << endl;
    cout << "*********************************" << endl;
    
    
}

void ExecuteScaling(){

	// Scale all the time
	dt = dt / EquScal;
	TransientTime = TransientTime / EquScal;
	
	InKick = InKick / EquScal;
	FinKick = FinKick / EquScal;
	dKick = dKick / EquScal;
	
}

//****************************************************************************************

vector <double> LaserEquation(vector <double> ENvec, double K){
	// Input: Vector ENvec corresponding to y(t) = (Ex, Ey, N) at some time t
	// 		  Varying parameter K.	

	// Output: Vector f(y,t), explained below 
	
	// Solving the following laser equation
	// dE/dt = i*Delta*E + beta*gamma*(1 - i*alpha)*N*E + K
	// dN/dt = lambda - N - (1 + beta*N_n)|E|^2

	// Which when you resolve into the its real and complex parts gives 
	// dEx/dt = - DeltaEy + beta*gamma*N*(Ex + alpha*Ey) + K
	// dEy/dt = DeltaEx + beta*gamma*N*(Ey - alpha*Ex)
	// dN/dt = lambda - N - (1 + beta*N_n)*(Ex^2 + Ey^2)
	
	// No we have an equation dy/dt = f(y,t) where 
	// y = (Ex, Ey, N) and f(y,t) is given by the above equations.
	// Thus we can solve for f(y,t). Where y = Evec.
	
	double f0 = - Delta*ENvec[1] + beta*smallgamma*ENvec[2]*(ENvec[0] + alpha*ENvec[1]) + K;
	double f1 = Delta*ENvec[0] + beta*smallgamma*ENvec[2]*(ENvec[1] - alpha*ENvec[0]);
	double f2 = lambda - ENvec[2] - ((1 + beta*ENvec[2])*(ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1]));

	vector <double> fyt(3);
	fyt[0] = f0;
	fyt[1] = f1;
	fyt[2] = f2;

	return fyt;
}

vector <double> JacobianTimesPLaserEquation(vector <double> ENvec, vector <double> p){
	// Here is the Jacobian of the Laser equation Df(x(t,x0))
	// Input: x at some time t with p
	// Output: The Jacobian Matrix
	
	double N = 3; // Known from Laser Equation
	double M = 3;
	
	vector <vector <double> > Df(N, vector<double> (M));
	
	// Column 0
	Df[0][0] = beta*smallgamma*ENvec[2];
	Df[1][0] = -1*alpha*beta*smallgamma*ENvec[2];
	Df[2][0] = -2*(1 + (beta*ENvec[2]))*ENvec[0];
	
	// Column 1
	Df[0][1] = alpha*beta*smallgamma*ENvec[2];
	Df[1][1] = beta*smallgamma*ENvec[2];
	Df[2][1] = -2*(1 + (beta*ENvec[2]))*ENvec[1];
	
	// Column 2
	Df[0][2] = beta*smallgamma*(ENvec[0] + alpha*ENvec[1]);
	Df[1][2] = beta*smallgamma*(ENvec[1] - alpha*ENvec[0]);
	Df[2][2] = -1 - beta*(ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1]);
	
	vector <double> Dfp(N);
	Dfp[0] = Df[0][0]*p[0] + Df[0][1]*p[1] + Df[0][2]*p[2]; 
	Dfp[1] = Df[1][0]*p[0] + Df[1][1]*p[1] + Df[1][2]*p[2]; 
	Dfp[2] = Df[2][0]*p[0] + Df[2][1]*p[1] + Df[2][2]*p[2]; 
	
	//OutVec(Dfp);
	
	return Dfp;
}

vector <double> LaserWithJacobian(vector <double> ENPvec, double K){

// Input the Vector Corresponding to Electric Field and Number of hole pairs and the displacement vector
// 		 And the Changing parameter

// Output The vector corresponding to f(x,K) and the Jacobian
	
	vector<double> DfD(6);

	DfD[0] = -1.0*Delta*ENPvec[1] + beta*smallgamma*ENPvec[2]*(ENPvec[0] + alpha*ENPvec[1]) + K;
	DfD[1] = Delta*ENPvec[0] + beta*smallgamma*ENPvec[2]*(ENPvec[1] - alpha*ENPvec[0]);
	DfD[2] = lambda - ENPvec[2] - ((1 + beta*ENPvec[2])*(ENPvec[0]*ENPvec[0] + ENPvec[1]*ENPvec[1]));
	
	DfD[3] = beta*smallgamma*ENPvec[2]*ENPvec[3] + (-1.0*Delta + beta*smallgamma*ENPvec[2]*alpha)*ENPvec[4] + beta*smallgamma*(ENPvec[0] + alpha*ENPvec[1])*ENPvec[5];
	DfD[4] = (Delta - alpha*beta*smallgamma*ENPvec[2])*ENPvec[3] + beta*smallgamma*ENPvec[2]*ENPvec[4] + beta*smallgamma*(ENPvec[1] - alpha*ENPvec[0])*ENPvec[5];
	DfD[5] = -2.0*(1 + (beta*ENPvec[2]))*ENPvec[0]*ENPvec[3] - 2.0*(1+(beta*ENPvec[2]))*ENPvec[1]*ENPvec[4] + (-1.0 - beta*(ENPvec[0]*ENPvec[0] + ENPvec[1]*ENPvec[1]))*ENPvec[5];
	
	//OutVec(DfD);
	//cout << K << endl;

	return DfD;

}

vector <double> RosslerEquation(vector <double> ENvec, double K){
	double c = 14.0;
	double b = 0.1;
	
	double f0 = -(ENvec[1] + ENvec[2]);
	double f1 = ENvec[0] + K*ENvec[1];
	double f2 = b + ENvec[2]*(ENvec[0] - c);
	
	vector <double> f(3);
	f[0] = f0;
	f[1] = f1;
	f[2] = f2;
	
	return f;
	
}

void OutMat(vector <vector <double> > Mat){
	for(int i = 0; i < Mat[0].size(); i++){
		for(int j = 0; j < Mat.size(); j++){
			cout << "  " << Mat[j][i];
		}
		cout << "\n" << endl;
	}
}

void OutVec(vector <double> Vec){
		for(int j = 0; j < Vec.size(); j++){
			cout << Vec[j] << endl;;
		}
}

void Matoutfile(vector <vector <double> > Mat ,std::string filename){
	
	std::ofstream ofile(filename, std::ios::out);
	
	ofile.precision(15);
	
	for(int i = 0; i < (int)Mat.size(); i++)
	{
		for(int j = 0; j < (int)Mat[0].size(); j++)
			{ofile << Mat[i][j] << " ";}
		ofile << "\n";	
	}
	ofile.close();
}

void Vecoutfile(vector <double> vec ,std::string filename){
	
	std::ofstream ofile(filename, std::ios::out);
	
	ofile.precision(15);
	
	for(int i = 0; i < (int)vec.size(); i++)
	{
			ofile << vec[i] << "\n";
	}
	ofile.close();
}

void MatToVecOutFiles(vector <vector <double> > Mat ,std::string filename){
	
	ofstream sizefile("NumOfFiles.txt", std::ios::out);
	sizefile << Mat.size() << "\n";	
	sizefile.close();
	
	for(int j = 0; j < Mat.size(); j++){
			
		vector<string> thisname;
  		thisname.push_back(filename);
  		string num = to_string(j);
  		thisname.push_back(num);
  		thisname.push_back(".txt");
  		
  		string thefilename = thisname[0] + thisname[1] + thisname[2]; 
  		
		ofstream ofile(thefilename, std::ios::out);
	
		ofile.precision(15);
	
		for(int i = 0; i < (int)Mat[j].size(); i++)
		{
			ofile << Mat[j][i] << "\n";	
		}
		ofile.close();
	}

}

static inline void loadbarMain(unsigned int x, unsigned int n, unsigned int w){
    if ( (x != n) && (x % (n/100+1) != 0) ) return;
 
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;
 
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}

static inline void loadbarMain(unsigned int x, unsigned int n, unsigned int w, clock_t start){
    if ( (x != n) && (x % (n/100+1) != 0) ) return;
 
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;
    
    clock_t current = clock();
    double tinsec = ((float)current - start)/CLOCKS_PER_SEC;
    double tinmin = tinsec/60.0;
    
    double tin1per = tinmin/(ratio*100.0);
    int trmper = ceil(int((100.0 - (ratio*100.0))*tin1per));
    
    if(trmper < 0){
    	trmper = 0;
    }
 
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "] Estimated Time Remaining: " << trmper << "mins \r" << flush;
}

vector <double> normalize(vector<double> v){
	
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





