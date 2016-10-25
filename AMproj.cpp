#include "Attach.h"

// Helper Functions
void OutMat(vector <vector <double> > Mat);
void OutVec(vector <double> Vec);
void Matoutfile(vector <vector <double> > Mat ,std::string filename);
void Vecoutfile(vector <double> vec ,std::string filename);
void MatToVecOutFiles(vector <vector <double> > Mat ,std::string filename);
static inline void loadbarMain(unsigned int x, unsigned int n, unsigned int w);
static inline void loadbarMain(unsigned int x, unsigned int n, unsigned int w, clock_t start);
vector <double> normalize(vector<double> v);

// Set the parameters - important
void SetParameters(string filename);
void ExecuteScaling();

// Here are the updated systems with forcing included in the Jacobian
// when doing the linearization.
vector<double> ClassCPlanerWithJacobian(vector<double> ENvec);
vector<double> ClassCPlanerWithJacobian_Kicking_Instance(vector<double> ENvec, double Kick_Size);

vector<double> ClassBPlanerWithJacobian(vector<double> ENvec);
vector<double> ClassBPlanerWithJacobian_Kicking_Instance(vector<double> ENvec, double Kick_Size);

vector<long double> ClassBPlanerWithJacobian(vector<long double> ENvec);
vector<long double> ClassBPlanerWithJacobian_Kicking_Instance(vector<long double> ENvec, double Kick_Size);


// Need a get phi function as this can be tricky as atan function wont deal with the full circle
double GetPhi(vector<double> ENvec);
long double GetPhi(vector<long double> ENvec);

int main(){

	cout << "Start the clock !!!\n";
	clock_t t = clock();


	SetParameters("Parameters.txt");
	ExecuteScaling();

//****************************************************************************************


	{
		double TimeEvlo = ((double)NumOfNorm*NormStep*dt);

		cout << "Total Number of steps:  " << NumOfNorm*NormStep << endl;


		// We are excuting on lig Scale
		double kickitlg = log10(InKick);
		double kickfnlg = log10(FinKick);
		double dklg = (kickfnlg - kickitlg)/NKick;

		//Now for the values we deal with
		double runkit = log10(IrunK);
		double runkfn = log10(FrunK);

		int numRuns = (int)((runkfn - runkit)/dklg);

		int vecsize = 5;
		vector <double> p0(vecsize);			vector <double> y0(vecsize);
		p0[0] = 0.001;						y0[0] = 1.0;
		p0[1] = -0.001;						y0[1] = -2.0;
		p0[2] = 0.001;						y0[2] = 3.0;
		p0[3] = 0.001;						y0[3] = 2.0;
		p0[4] = 0.001;						y0[4] = 5.1;

		/*vector <long double> p0(vecsize);			vector <long double> y0(vecsize);
		p0[0] = -0.00001;							y0[0] = -2.1;
		p0[1] = -0.00001;							y0[1] = 1.5;
		p0[2] = -0.00001;							y0[2] = 6.1;*/

		vector<double> h;
		vector<double> k;

		clock_t st = clock();
		for(int i = 0; i < numRuns; i++){
			Lyapunov LyapLase(NormStep, dt, TransientTime, TimeEvlo, y0, p0);
			double t_kick = pow(10, (runkit + i*dklg) );
			//cout << t_kick << endl;
			h.push_back(LyapLase.CalcBigLypunov_Kick_new(ClassCPlanerWithJacobian, ClassCPlanerWithJacobian_Kicking_Instance, t_kick, Perturb ));
			k.push_back(runkit + i*dklg);
			loadbarMain(i, NKick, 50, st);
		}

		string st0 = "LyapunovSlice_C_3_";
		string st1 = "Kvals_C_3_";
		string num = to_string(FlNum);
		string ext = ".txt";

		st0.append(num);
		st0.append(ext);

		st1.append(num);
		st1.append(ext);

		cout << st0 << endl;

		Vecoutfile(h , st0);
		Vecoutfile(k , st1);
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

vector<double> ClassCPlanerWithJacobian(vector<double> ENvec){

	vector<double> f(10);

	// Here, the dx/dt = f(x)
	f[0] = - ENvec[0] - Dp*ENvec[1] - ENvec[3];
	f[1] = - ENvec[1] + Dp*ENvec[0] + ENvec[2];
	f[2] = Gpc*( -1*ENvec[2] + Dp*ENvec[3] + ENvec[1]*ENvec[4] );
	f[3] = Gpc*( -1*ENvec[3] - Dp*ENvec[2] - ENvec[0]*ENvec[4] );
	f[4] = Gnc*( Lambda - ENvec[4] + ENvec[3]*ENvec[0] - ENvec[2]*ENvec[1] );

	// Jacobian
	f[5] = (-1.0)*ENvec[5] + (-1.0*Dp)*ENvec[6] - ENvec[8];
	f[6] = (Dp)*ENvec[5] +  (-1.0)*ENvec[6] + ENvec[7];

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

vector<double> ClassCPlanerWithJacobian_Kicking_Instance(vector<double> ENvec, double Kick_Size){

	vector<double> f(10);

	// From Documentation "Linearization/Linearization.pdf" I am including the functions f_x(T) and f_y(T)
	double phi = GetPhi(ENvec);
	double f_x = Kick_Size*cos(phi)*sin(phi*NumPet);
	double f_y = Kick_Size*sin(phi)*sin(phi*NumPet);

	// Here, the dx/dt = f(x)
	f[0] = ENvec[0] + f_x;
	f[1] = ENvec[1] + f_y;
	f[2] = ENvec[2];
	f[3] = ENvec[3];
	f[4] = ENvec[4];


	// Again from Documentation "Linearization/Linearization.pdf" I am including the necessary differentials of f_x(T) and f_y(T)
	double Dphi_Ex = (-ENvec[1])/( ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1] );
	double Dphi_Ey = ( ENvec[0])/( ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1] );

	double Df_xx = Dphi_Ex*Kick_Size*( NumPet*cos(phi)*cos(phi*NumPet) - sin(phi)*sin(phi*NumPet) );
	double Df_xy = Dphi_Ey*Kick_Size*( NumPet*cos(phi)*cos(phi*NumPet) - sin(phi)*sin(phi*NumPet) );
	double Df_yx = Dphi_Ex*Kick_Size*( NumPet*sin(phi)*cos(phi*NumPet) + cos(phi)*sin(phi*NumPet) );
	double Df_yy = Dphi_Ey*Kick_Size*( NumPet*sin(phi)*cos(phi*NumPet) + cos(phi)*sin(phi*NumPet) );

	// Jacobian
	f[5] = ENvec[5] + Df_xx*ENvec[5] + Df_xy*ENvec[6];
	f[6] = ENvec[6] + Df_yx*ENvec[5] + Df_yy*ENvec[6];
	f[7] = ENvec[7];
	f[8] = ENvec[8];
	f[9] = ENvec[9];


	return f;

}

//****************************************************************************************

vector<double> ClassBPlanerWithJacobian(vector<double> ENvec){

	vector<double> f(6);

	// Here, the dx/dt = f(x)
	f[0] = ( (ENvec[2]/(1.0 + Dp*Dp)) - 1.0)*(ENvec[0] + Dp*ENvec[1]);
	f[1] = ( (ENvec[2]/(1.0 + Dp*Dp)) - 1.0)*(ENvec[1] - Dp*ENvec[0]);
	f[2] = Gnc*( Lambda - ENvec[2] - ((ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1])*ENvec[2])/(1.0 + Dp*Dp));

	// Jacobian
	double alf = (ENvec[2]/(1 + Dp*Dp)) - 1.0;
	double Lor = (1.0/(1.0 + Dp*Dp));
	double EEE = ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1];

	f[3] = (alf)*ENvec[3] + (Dp*alf)*ENvec[4] + (ENvec[0] + Dp*ENvec[1])*Lor*ENvec[5];
	f[4] = (-Dp*alf)*ENvec[3] + (alf)*ENvec[4] + (ENvec[1] - Dp*ENvec[0])*Lor*ENvec[5];
	f[5] = -Gnc*( 2.0*ENvec[0]*ENvec[2]*Lor*ENvec[3] + 2.0*ENvec[1]*ENvec[2]*Lor*ENvec[4] + ENvec[5] + EEE*Lor*ENvec[5] );

	f[0] = EquScal*f[0];
	f[1] = EquScal*f[1];
	f[2] = EquScal*f[2];

	f[3] = EquScal*f[3];
	f[4] = EquScal*f[4];
	f[5] = EquScal*f[5];


	return f;

}

vector<double> ClassBPlanerWithJacobian_Kicking_Instance(vector<double> ENvec, double Kick_Size){

	// From Documentation "Linearization/Linearization.pdf" I am including the functions f_x(T) and f_y(T)
	double phi = GetPhi(ENvec);
	double f_x = Kick_Size*cos(phi)*sin(phi*NumPet);
	double f_y = Kick_Size*sin(phi)*sin(phi*NumPet);

	// Again from Documentation "Linearization/Linearization.pdf" I am including the necessary differentials of f_x(T) and f_y(T)
	double Dphi_Ex = (-ENvec[1])/( ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1] );
	double Dphi_Ey = ( ENvec[0])/( ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1] );

	double Df_xx = Dphi_Ex*Kick_Size*( NumPet*cos(phi)*cos(phi*NumPet) - sin(phi)*sin(phi*NumPet) );
	double Df_xy = Dphi_Ey*Kick_Size*( NumPet*cos(phi)*cos(phi*NumPet) - sin(phi)*sin(phi*NumPet) ); // d/d(Ey) fx(T)
	double Df_yx = Dphi_Ex*Kick_Size*( NumPet*sin(phi)*cos(phi*NumPet) + cos(phi)*sin(phi*NumPet) ); // d/d(Ex) fy(T)
	double Df_yy = Dphi_Ey*Kick_Size*( NumPet*sin(phi)*cos(phi*NumPet) + cos(phi)*sin(phi*NumPet) );


	vector<double> f(6);

	f[0] = ENvec[0] + f_x;
	f[1] = ENvec[1] + f_y;
	f[2] = ENvec[2];

	f[3] = ENvec[3] + ENvec[3]*Df_xx + ENvec[4]*Df_xy;
	f[4] = ENvec[4] + ENvec[3]*Df_yx + ENvec[4]*Df_yy;
	f[5] = ENvec[5];

	return f;

}

//****************************************************************************************

vector<long double> ClassBPlanerWithJacobian(vector<long double> ENvec){

	vector<long double> f(6);

	// Here, the dx/dt = f(x)
	f[0] = ( (ENvec[2]/(1.0 + Dp*Dp)) - 1.0)*(ENvec[0] + Dp*ENvec[1]);
	f[1] = ( (ENvec[2]/(1.0 + Dp*Dp)) - 1.0)*(ENvec[1] - Dp*ENvec[0]);
	f[2] = Gnc*( Lambda - ENvec[2] - ((ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1])*ENvec[2])/(1.0 + Dp*Dp));

	// Jacobian
	long double alf = (ENvec[2]/(1 + Dp*Dp)) - 1.0;
	long double Lor = (1.0/(1.0 + Dp*Dp));
	long double EEE = ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1];

	f[3] = (alf)*ENvec[3] + (Dp*alf)*ENvec[4] + (ENvec[0] + Dp*ENvec[1])*Lor*ENvec[5];
	f[4] = (-Dp*alf)*ENvec[3] + (alf)*ENvec[4] + (ENvec[1] - Dp*ENvec[0])*Lor*ENvec[5];
	f[5] = -Gnc*( 2.0*ENvec[0]*ENvec[2]*Lor*ENvec[3] + 2.0*ENvec[1]*ENvec[2]*Lor*ENvec[4] + ENvec[5] + EEE*Lor*ENvec[5] );

	f[0] = EquScal*f[0];
	f[1] = EquScal*f[1];
	f[2] = EquScal*f[2];

	f[3] = EquScal*f[3];
	f[4] = EquScal*f[4];
	f[5] = EquScal*f[5];


	return f;

}

vector<long double> ClassBPlanerWithJacobian_Kicking_Instance(vector<long double> ENvec, double Kick_Size){

	// From Documentation "Linearization/Linearization.pdf" I am including the functions f_x(T) and f_y(T)
	long double phi = GetPhi(ENvec);
	long double f_x = Kick_Size*cos(phi)*sin(phi*NumPet);
	long double f_y = Kick_Size*sin(phi)*sin(phi*NumPet);

	// Again from Documentation "Linearization/Linearization.pdf" I am including the necessary differentials of f_x(T) and f_y(T)
	long double Dphi_Ex = (-ENvec[1])/( ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1] );
	long double Dphi_Ey = ( ENvec[0])/( ENvec[0]*ENvec[0] + ENvec[1]*ENvec[1] );

	long double Df_xx = Dphi_Ex*Kick_Size*( NumPet*cos(phi)*cos(phi*NumPet) - sin(phi)*sin(phi*NumPet) );
	long double Df_xy = Dphi_Ey*Kick_Size*( NumPet*cos(phi)*cos(phi*NumPet) - sin(phi)*sin(phi*NumPet) ); // d/d(Ey) fx(T)
	long double Df_yx = Dphi_Ex*Kick_Size*( NumPet*sin(phi)*cos(phi*NumPet) + cos(phi)*sin(phi*NumPet) ); // d/d(Ex) fy(T)
	long double Df_yy = Dphi_Ey*Kick_Size*( NumPet*sin(phi)*cos(phi*NumPet) + cos(phi)*sin(phi*NumPet) );


	vector<long double> f(6);

	f[0] = ENvec[0] + f_x;
	f[1] = ENvec[1] + f_y;
	f[2] = ENvec[2];

	f[3] = ENvec[3] + ENvec[3]*Df_xx + ENvec[4]*Df_xy;
	f[4] = ENvec[4] + ENvec[3]*Df_yx + ENvec[4]*Df_yy;
	f[5] = ENvec[5];

	return f;

}

//****************************************************************************************

void SetParameters(string filename){

	fstream myfile(filename.c_str(), std::ios_base::in);

    myfile >> Lambda >> Dp >> Gnc >> Gpc >> TransientTime >> dt >> InKick >> FinKick
    	   >> NKick >> IrunK >> FrunK >> FlNum >> Perturb >> NumPet >> NormStep >> NumOfNorm >> EquScal;
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
	cout << "Now doing log scale:" << endl;
	cout << "Initial Kick nonlog yet = " << InKick << endl;
    cout << "Final Kick nonlog yet   = " << FinKick << endl;
    cout << "Num of diff intvls      = " << NKick << endl;
	cout << "Inital Kick for this run = " << IrunK << endl;
	cout << "Final Kick for this run  = " << FrunK << endl;
	cout << "File Number              = " << FlNum << endl;
    cout << "Kick Size                = " << Perturb << endl;
    cout << "Number of petals         = " << NumPet << endl;
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


	IrunK = IrunK / EquScal;
	FrunK = FrunK / EquScal;
}



//****************************************************************************************

double GetPhi(vector<double> ENvec){

	double phi = fabs(atan(ENvec[1]/ENvec[0]));

	if(ENvec[0] == 0){
		if(ENvec[1] > 0 )
			return PI/2.0;
		else
			return 3*PI/2.0;
	}
	else if( ENvec[0] > 0 && ENvec[1] >= 0 )
		return phi;
	else if( ENvec[0] < 0 && ENvec[1] >= 0 )
		return PI - phi;
	else if( ENvec[0] < 0 && ENvec[1] <= 0 )
		return PI + phi;
	else
		return 2*PI - phi;

}

long double GetPhi(vector<long double> ENvec){

	double phi = fabs(atan(ENvec[1]/ENvec[0]));

	if(ENvec[0] == 0){
		if(ENvec[1] > 0 )
			return PI/2.0;
		else
			return 3*PI/2.0;
	}
	else if( ENvec[0] > 0 && ENvec[1] >= 0 )
		return phi;
	else if( ENvec[0] < 0 && ENvec[1] >= 0 )
		return PI - phi;
	else if( ENvec[0] < 0 && ENvec[1] <= 0 )
		return PI + phi;
	else
		return 2*PI - phi;

}

//****************************************************************************************

void OutMat(vector <vector <double> > Mat){
	for(int i = 0; i < (int)Mat[0].size(); i++){
		for(int j = 0; j < (int)Mat.size(); j++){
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

	std::ofstream ofile(filename.c_str(), std::ios::out);

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

	std::ofstream ofile(filename.c_str(), std::ios::out);

	ofile.precision(15);

	for(int i = 0; i < (int)vec.size(); i++)
	{
			ofile << vec[i] << "\n";
	}
	ofile.close();
}

/*void MatToVecOutFiles(vector <vector <double> > Mat ,std::string filename){

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
*/
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
