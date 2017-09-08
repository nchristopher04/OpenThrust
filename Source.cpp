#include <stdafx.h>
#include <iostream>
#include <cmath>
#include "RPA_to_struct.h"
#include "thermo_functions.h"

const double PSI_TO_PA = 6894.76; 

using namespace std;
// Will be defined later in program
double oxyMass;							// Initial oxidizer mass [kg]
double Pc;								// Chamber pressure [psi]
double k;								// Heat capacity ratio []
double R;								// Specific gas constant [kJ/kg*k]
double Tc;								// Chamber temperature [k]
double Cf;								// Thrust coefficient []

// Already defined
double At = 0.00153058;					// Nozzle throat area [m^2] (taken from CAD drawing)
double A2 = 0.00724004;					// Nozzle exit area [m^2] (taken from CAD drawing)
double OF = 2.1;						// Oxidizer-fuel ratio assumed constant to start
double timeStep = 0.1;					// [s]
double mDotNozzle, mDotInjector;		// Mass flow rates at the nozzle and the injector [kg/s]
double time[1000], thrust[1000];		// Output arrays that give thrust over time


int main() {

	cout << "Input initial oxidizer mass in [kg]:  ";
	cin >> oxyMass;
	
	for (int x = 0; x < 1000; x++) { //time steps
		
		// Finds all relevant values for thrust
		
		RPALookup(Pc, OF, k, R, Tc);
		mDotNozzle = massFlowRate(At, Pc, k, R, Tc);
		mDotInjector = massFlowRateInjector(mDotNozzle, OF);
		Cf = thrustCoefficient(14.7, A2, Pc);


		if (mDotNozzle < 0) { throw "massFlowNegative"; break; }
		else { oxyMass -= mDotInjector*timeStep; }

		// Creates outputs for each timestep

		time[x] = x*timeStep;
		thrust[x] = At*(Pc*PSI_TO_PA)*Cf;
		cout << "T+" << time << " s =>>> Oxy Mass: " <<  oxyMass << "kg | Chamber Pressure: " <<  Pc << " psi" << endl;
		

		if (oxyMass <= 0.01) { cout << "Empty"; break; };
	}

}

double massFlowRate(double nozzleArea, double Pc, double k, double R, double Tc) {
	// Calculates mass flow rate out of the nozzle
	// Inputs in [m^2], [psi], [], [kJ/kg*K], [k]
	// Outputs in [kg/s]
	Pc = Pc*PSI_TO_PA;
	double mDot = (nozzleArea*Pc*k*sqrt( pow((2 /(k + 1)), ((k + 1) / (k - 1))))) / sqrt(k*R*Tc);
	return mDot;
}

double massFlowRateInjector(double nozzleFlow, double OF_ratio) {
	// Calculates mass flow rate out of the injector
	// Inputs in [kg/s]
	// Outputs in [kg/s]
	double mDotI = nozzleFlow*(OF_ratio / (1 + OF_ratio));
	return mDotI;
}

double thrustCoefficient(double Patm, double A2, double Pc) {
	// Calculates thrust coefficient for different chamber pressures
	// Inputs in [psi], [m^2], [psi]
	// Ouput is unitless
	Patm = Patm*PSI_TO_PA;
	Pc = Pc*PSI_TO_PA;
	double C12 = 977;						//calculate C1,2
	double Cf = C12 - (Patm*A2) / (Pc*At);
	return Cf; 
}

void RPALookup(float Pc, double OF, double& k, double& R, double& Tc) {
	// Find k, R, Tc from table and return them to the program
	// Inputs in [psi], []
	// No output but stores values in variables k, R, Tc
	// Stored in [], [kJ/kg*K], [k]
	RPA_Table CombustionProps = lookUp(Pc, OF, Create_Table_Array());
	k = CombustionProps.k_value;
	R = CombustionProps.R_value;
	Tc = CombustionProps.Chamber_Temperture;
}