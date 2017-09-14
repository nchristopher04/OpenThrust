#include <iostream>
#include <cmath>
#include "RPA_to_struct.h"
#include "thermo_functions.h"
#include <fstream>
#include <stdexcept>
#include "injector_Model.h"
#include "Source.h"

using namespace std;

const double PSI_TO_PA = 6894.76; 

// Will be defined later in program
double oxyMass;							// Initial oxidizer mass [kg]
float Pc;								// Chamber pressure [psi]
double k;								// Heat capacity ratio []
double R;								// Specific gas constant [kJ/kg*k]
double Tc;								// Chamber temperature [k]
double Cf;								// Thrust coefficient []
double T1;								// NOS Tank Temperature
struct options {
	int flowModel;
}MainX;

// Already defined
double At = 0.00153058;					// Nozzle throat area [m^2] (taken from CAD drawing)
double A2 = 0.00724004;					// Nozzle exit area [m^2] (taken from CAD drawing)
float OF = 2.1f;						// Oxidizer-fuel ratio assumed constant to start
double timeStep = 0.1;					// [s]
double mDotNozzle, mDotInjector;		// Mass flow rates at the nozzle and the injector [kg/s]
double time[1000], thrust[1000];		// Output arrays that give thrust over time



template <typename Arg, typename... Args>
void output(ofstream& out, Arg&& arg, Args&&... args)
{
	out << forward<Arg>(arg);
	using expander = int[];
	(void)expander {
		0, (void(out << ',' << std::forward<Args>(args)), 0)...
	};
	out << '\n';
} //this is a variadic print function to output.csv

int main() {
	ofstream simFile("output.csv");
	simFile << "Time (s), Liquid Mass   ,  Chamber Pressure , Thrust , Mass Flow Rate " << '\n'; //setup basic output format
	cout << "Input initial params." << endl;
	cout << "NOS Temperature; Celcius";
	cin >> T1;
	cout << "Abs. Chamber Pressure; PSI";
	cin >> Pc;

	cout << "Input initial oxidizer mass in [kg]:  ";
	cin >> oxyMass;
	
	for (int x = 0; x < 1000; x++) { //time steps
		
		// Finds all relevant values for thrust
		RPALookup(Pc, OF, k, R, Tc);
		mDotNozzle = massFlowRate(At, Pc, k, R, Tc);
		if (MainX.flowModel == 2) {
			int PcRound10 = Pc; //cast chamber pressure to int.
			PcRound10 += 5;
			PcRound10 -= PcRound10 % 10; 
			mDotInjector = injectorModel(T1, PcRound10);
		}
		else mDotInjector = massFlowRateInjector(mDotNozzle, OF);
		Cf = thrustCoefficient(14.7, A2, Pc);

		try {
			if (mDotNozzle < 0 || mDotInjector < 0) { throw "massFlowNegative"; }
		}
		catch (exception& e)
		{
			cout << e.what() << '\n'; //catch exception, display to user and break loop
			break;
		}
		oxyMass -= mDotInjector*timeStep; 

		// Creates outputs for each timestep

		time[x] = x*timeStep;
		thrust[x] = At*(Pc*PSI_TO_PA)*Cf;
		cout << "T+" << time << " s =>>> Oxy Mass: " <<  oxyMass << "kg | Chamber Pressure: " <<  Pc << " psi" << endl;
		output(simFile,time, oxyMass, Pc);//output to csv

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

void RPALookup(float Pc, float OF, double &k, double &R, double &Tc) {
	// Find k, R, Tc from table and return them to the program
	// Inputs in [psi], []
	// No output but stores values in variables k, R, Tc
	// Stored in [], [kJ/kg*K], [k]
	RPA_Table CombustionProps = lookUp(Pc, OF, Create_Table_Array());
	k = CombustionProps.k_value;
	R = CombustionProps.R_value;
	Tc = CombustionProps.Chamber_Temperture;
}

double tankProps(double oxyMass, double Pc, double &Temp, double &TankPressure) {
	
	
	//TankPressure = nox_vp(Temp);
	return 0;
}
//NOS properties from Modelling the Nitrous Run tank Emptying
const float pCrit = 72.51f; /* critical pressure, Bar Abs */
const float rhoCrit = 452.0f; /* critical density, kg/m3 */
const float ZCrit = 0.28f; /* critical compressibility factor */
const float gamma = 1.3f; /* average over subcritical range */

/* Nitrous oxide vapour pressure, Bar */
double nox_vp(double T_Celcius)
{
	const float p[4] = { 1.0f, 1.5f, 2.5f, 5.0f };
	const float b[4] = { -6.71893f, 1.35966f, -1.3779f, -4.051f };
	double Tr = reduced_temperature(T_Celcius);
	float rab = 1.0 - Tr;
	float shona = 0.0;
	for (int dd = 0; dd < 4; dd++)
		shona += b[dd] * pow(rab, p[dd]);
	double bob = pCrit * exp((shona / Tr));
	return(bob);
}