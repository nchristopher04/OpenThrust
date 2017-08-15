#include <stdafx.h>
#include <iostream>
#include <cmath>
#include "RPA_to_struct.h"
#include "thermo_functions.h"

using namespace std;
double tankPressure,
double oxyMass;
double At = 10;							// Nozzle throat area
double timeStep=0.1;
float OF = 2.1;							// Oxidizer-fuel ratio assumed constant to start
double Pc, k, R, Tc, nozzleArea,mDot;	// Start values
double time; 


int main() {
	for (int x; x < 1000; x++) { //time steps
		RPALookup(Pc, OF, k, R, Tc);
		mDot = massFlowRate(nozzleArea, Pc, k, R, Tc);
		if (mDot < 0) { throw "massFlowNegative"; break; }
		else { oxyMass -= mDot*timeStep; }
		//compute new chamber pressure
		time = x*timeStep;
		cout << time << oxyMass << Pc;
		if (oxyMass <= 0.01) { break; };
	}

}

double massFlowRate(double nozzleArea, double Pc, double k, const float R, double Tc) {
	double mDot = nozzleArea*Pc*k*sqrt( pow((2 *(k + 1)), ((k + 1) / (k - 1)))) / sqrt(k*R*Tc);
	return mDot;
}

double thrustCoefficient(double Patm, double A2, double Pc) {
	double C12 = 977;//calculate C1,2
	double Cf = C12 - Patm*A2 / (Pc*At);
	return Cf; 
}

void RPALookup(float Pc, double OF, double& k, double& R, double& Tc) {
	//find k, R, Tc from table and return them to the program
	RPA_Table CombustionProps = lookUp(Pc, OF, Create_Table_Array());
	k = CombustionProps.k_value;
	R = CombustionProps.R_value;
	Tc = CombustionProps.Chamber_Temperture;
}