#ifndef RUNGEKUTTA_H
#include "RungeKutta.h"
#endif

#include "Attach.h"

void OutMatRK(vector <vector <double> > Mat);
void OutVecRK(vector <double> Vec);
void OutVecRK(vector <long double> Vec);

//Public

RungeKutta::RungeKutta(){
	//cout << "This is the RungeKutta solver.\n";
	h = 0.1;
	N = 1000;
}

RungeKutta::RungeKutta(int NumSteps, double dt){
	//cout << "This is the RungeKutta solver.\n";
	
	h = dt;
	N = NumSteps;

}

RungeKutta::RungeKutta(double TimeEvlo, double dt){
	//cout << "This is the RungeKutta solver.\n";
	
	h = dt;
	N = (int)(TimeEvlo/dt);
}


vector <vector <double> > RungeKutta::RK4(vector <double> (f_yt)(vector<double> vec, double param), vector<double> y0, double K){
	// This is a function that will carry out the RK4 method of solving ODEs 
	// Input: Function containing the f(y,t) of the ODE dy/dt = f(y,t) with some varying parameter
	//        Some initial condition y0
	
	// Output: A matrix Y in which each column
	//		   will contain all of your point calculated via the RK4
	
	int NDim = (int)y0.size();
	
	// Create matrix to hold points 
	vector <vector <double> > Y(N, vector <double>(NDim) );
	
	// Fill the initial column with initial vector
	Y[0] = y0;
	
	// Create two temporary vectors 
	vector <double> TmpVec1(NDim);
	vector <double> TmpVec2(NDim);
	
	// Lets start the RK4
	for(int i = 0; i < N-1; i++){
	
		// Slope 1
		k1 = f_yt(Y[i], K);	
		
		// Slope 2
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = Y[i][j] + 0.5*h*k1[j];
		k2 = f_yt(TmpVec1, K);
		
		
		// Slope 3
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = Y[i][j] + 0.5*h*k2[j];
		k3 = f_yt(TmpVec1, K);
		
		
		// Slope 4
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = Y[i][j] + h*k3[j];
		k4 = f_yt(TmpVec1, K);
		
		
		for(int j = 0; j < NDim; j++){
			TmpVec1[j] = Y[i][j] + (1.0/(double)6.0)*h*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
		}
			
		// Fill the initial column with initial vector
		Y[i+1] = TmpVec1;
	}
	
	return Y;
}
		
vector <vector <double> > RungeKutta::RK4(vector <double> (f_yt)(vector<double> vec), vector<double> y0){
	// This is a function that will carry out the RK4 method of solving ODEs 
	// Input: Function containing the f(y,t) of the ODE dy/dt = f(y,t) with some varying parameter
	//        Some initial condition y0
	
	// Output: A matrix Y in which each column
	//		   will contain all of your point calculated via the RK4
	
	int NDim = (int)y0.size();
	
	// Create matrix to hold points 
	vector <vector <double> > Y(N, vector <double>(NDim) );
	
	// Fill the initial column with initial vector
	Y[0] = y0;
	
	// Create two temporary vectors 
	vector <double> TmpVec1(NDim);
	vector <double> TmpVec2(NDim);
	
	// Lets start the RK4
	for(int i = 0; i < N-1; i++){
	
		// Slope 1
		k1 = f_yt(Y[i]);	
		
		// Slope 2
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = Y[i][j] + 0.5*h*k1[j];
		k2 = f_yt(TmpVec1);
		
		
		// Slope 3
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = Y[i][j] + 0.5*h*k2[j];
		k3 = f_yt(TmpVec1);
		
		
		// Slope 4
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = Y[i][j] + h*k3[j];
		k4 = f_yt(TmpVec1);
		
		
		for(int j = 0; j < NDim; j++){
			TmpVec1[j] = Y[i][j] + (1.0/(double)6.0)*h*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
		}
			
		// Fill the initial column with initial vector
		Y[i+1] = TmpVec1;
	}
	
	return Y;
}



vector <double> RungeKutta::RK4_11(vector <double> (f_yt)(vector<double> vec, double param), vector<double> ENPvec, double K){
	 
	 // This is a function that will carry out the RK4 method of solving ODEs 
	// Input: Function containing the f(y,t) of the ODE dy/dt = f(y,t) with some varying parameter
	//        Some initial condition y0
	
	// Output: A matrix Y in which each column
	//		   will contain all of your point calculated via the RK4
	
	int NDim = (int)ENPvec.size();
	
	// Create two temporary vectors 
	vector <double> TmpVec1(NDim);
	vector <double> TmpVec2(NDim);
	
	for(int i = 0; i < NDim; i++ ){
			TmpVec2[i] = ENPvec[i];
			//cout << "The " << i << " component is " << ENPvec[i] << endl;
	}
	
	
	
	// Lets start the RK4
		// Slope 1
		k1 = f_yt(ENPvec, K);	
		
		// Slope 2
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + 0.5*h*k1[j];
		k2 = f_yt(TmpVec1, K);
		
		
		// Slope 3
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + 0.5*h*k2[j];
		k3 = f_yt(TmpVec1, K);
		
		
		// Slope 4
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + h*k3[j];
		k4 = f_yt(TmpVec1, K);
		
		for(int j = 0; j < NDim; j++){
			TmpVec1[j] = ENPvec[j] + (1.0/(double)6.0)*h*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
		}
	
	//cout<<"here - - - - - - - - - - - - - - - - - - - - - -"<<endl;
	//OutVecRK(TmpVec1);
	
	return TmpVec1;
	
	


}

vector <double> RungeKutta::RK4_11(vector <double> (f_yt)(vector<double> vec, int Current_Step, double Kick_Size, double Time_Between_Kicks), 
										vector<double> ENPvec, int C_Step, double Kick_Size, double Time_Between_Kicks){
	 
	 // This is a function that will carry out the RK4 method of solving ODEs 
	// Input: Function containing the f(y,t) of the ODE dy/dt = f(y,t) with some varying parameter
	//        Some initial condition y0
	
	// Output: A matrix Y in which each column
	//		   will contain all of your point calculated via the RK4
	
	int NDim = (int)ENPvec.size();
	
	// Create two temporary vectors 
	vector <double> TmpVec1(NDim);
	vector <double> TmpVec2(NDim);
	
	for(int i = 0; i < NDim; i++ ){
			TmpVec2[i] = ENPvec[i];
			//cout << "The " << i << " component is " << ENPvec[i] << endl;
	}
	
	
	
	// Lets start the RK4
		// Slope 1
		k1 = f_yt(ENPvec, C_Step, Kick_Size, Time_Between_Kicks);	
		
		// Slope 2
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + 0.5*h*k1[j];
		k2 = f_yt(TmpVec1, C_Step, Kick_Size, Time_Between_Kicks);	
		
		
		// Slope 3
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + 0.5*h*k2[j];
		k3 = f_yt(TmpVec1, C_Step, Kick_Size, Time_Between_Kicks);	
		
		
		// Slope 4
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + h*k3[j];
		k4 = f_yt(TmpVec1, C_Step, Kick_Size, Time_Between_Kicks);	
		
		for(int j = 0; j < NDim; j++){
			TmpVec1[j] = ENPvec[j] + (1.0/(double)6.0)*h*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
		}
	
	//cout<<"here - - - - - - - - - - - - - - - - - - - - - -"<<endl;
	//OutVecRK(TmpVec1);
	
	return TmpVec1;
	
	


}
		
vector<long double> RungeKutta::RK4_ld(vector <long double> (f_yt)(vector<long double> vec), vector<long double> ENPvec){
	 
	// This is a function that will carry out the RK4 method of solving ODEs 
	// Input: Function containing the f(y,t) of the ODE dy/dt = f(y,t) with some varying parameter
	//        Some initial condition y0
	
	// Output: A matrix Y in which each column
	//		   will contain all of your point calculated via the RK4
	
	int NDim = (int)ENPvec.size();
	
	// Create two temporary vectors 
	vector <long double> TmpVec1(NDim);
	vector <long double> TmpVec2(NDim);
	
	for(int i = 0; i < NDim; i++ ){
			TmpVec2[i] = ENPvec[i];
	}
	
	// Lets start the RK4
		// Slope 1
		k1_l = f_yt(ENPvec);	
		
		// Slope 2
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + 0.5*h*k1_l[j];
		k2_l = f_yt(TmpVec1);
		
		
		// Slope 3
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + 0.5*h*k2_l[j];
		k3_l = f_yt(TmpVec1);
		
		
		// Slope 4
		for(int j = 0; j < NDim; j++)
			TmpVec1[j] = ENPvec[j] + h*k3_l[j];
		k4_l = f_yt(TmpVec1);
		
		for(int j = 0; j < NDim; j++){
			TmpVec1[j] = ENPvec[j] + (1.0/(long double)6.0)*h*(k1_l[j] + 2*k2_l[j] + 2*k3_l[j] + k4_l[j]);
		}
	
	//cout<<"here - - - - - - - - - - - - - - - - - - - - - -"<<endl;
	//OutVecRK(TmpVec1);
	
	return TmpVec1;
}
		
//****************************************************************************************
// Private
void RungeKutta::FillWithVec(vector <double> v1, vector <double> v2){
	// Fills vector v1 with vector v2
	
	for(int i = 0; i < (int)v2.size(); i++){
		v1[i] = v2[i];
	}
	
}


//****************************************************************************************
// Helper
void OutMatRK(vector <vector <double> > Mat){
	for(int i = 0; i < Mat[0].size(); i++){
		for(int j = 0; j < Mat[0].size(); j++){
			cout << "  " << Mat[i][j];
		}
		cout << "\n" << endl;
	}
}

void OutVecRK(vector <double> Vec){
		for(int j = 0; j < Vec.size(); j++){
			cout << Vec[j] << endl;
		}
}

void OutVecRK(vector <long double> Vec){
		for(int j = 0; j < Vec.size(); j++){
			cout << Vec[j] << endl;
		}
}




















