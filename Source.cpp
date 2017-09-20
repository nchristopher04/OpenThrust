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
double Pc, PcNew, PcOld;				// Chamber pressure [psi]
double k;								// Heat capacity ratio []
double R;								// Specific gas constant [kJ/kg*k]
double Tc;								// Chamber temperature [k]
double Cf;								// Thrust coefficient []
double T1;								// NOS Tank Temperature
int PcRound10;							// Casts chamber pressure to integer
double err;								// Used for calculating relative error
struct options {
	int flowModel;
}MainX;

// Already defined
double At = 0.000382646;				// Nozzle throat area [m^2]
double A2 = 0.00181001;					// Nozzle exit area [m^2]
double OF = 2.1;						// Oxidizer-fuel ratio assumed constant to start
double timeStep = 0.1;					// [s]
double mDotNozzle, mDotInjector;		// Mass flow rates at the nozzle and the injector [kg/s]
double time[1000], thrust[1000];		// Output arrays that give thrust over time
Look_Up_Table Table_Array = Create_Table_Array();


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
	simFile << "Time (s), Liquid Mass (kg)  ,  Chamber Pressure (PSI), Thrust (N), Mass Flow Rate (kg/s)" << '\n'; //setup basic output format
	cout << "Input initial params." << endl;
	cout << "Initial NOS Temperature [Celsius]:  ";
	cin >> T1;
	cout << "Abs. Chamber Pressure [PSI]:  ";
	cin >> Pc;

	cout << "Input initial oxidizer mass in [kg]:  ";
	cin >> oxyMass;

	 MainX.flowModel = 2;
	
	for (int x = 0; x < 1000; x++) { //time steps
		
		// Finds all relevant values for thrust
		RPALookup(Pc, OF, k, R, Tc);
		mDotNozzle = massFlowRate(At, Pc, k, R, Tc);
		if (MainX.flowModel == 2) {
			PcOld = Pc;
			for (int i = 0; i < 100; i++) {
				PcRound10 = PcOld; //cast chamber pressure to int.
				PcRound10 += 5;
				PcRound10 -= PcRound10 % 10;
				mDotInjector = injectorModel(T1, PcRound10);
				mDotNozzle = massFlowRateNozzle(mDotInjector, OF);
				PcNew = calcPc(At, mDotNozzle, k, R, Tc);
				err = 100 * (PcOld - PcNew) / PcOld;
				PcOld = PcNew;
				if (err < 5) { Pc = PcNew; break; }
				else if (i == 99) { throw "PressureCalculatorDiverged"; }
			}
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
		cout << "T+" << time[x] << " s =>>> Oxy Mass: " <<  oxyMass << "kg | Chamber Pressure: " <<  Pc << " psi | " << 
			"Injector flow rate: " << mDotInjector << " kg/s" << endl;
		output(simFile,time[x], oxyMass, Pc, thrust[x], mDotInjector);//output to csv

		if (oxyMass <= 0.01) { cout << "Empty"; cin; break; };
	}

}

double massFlowRate(double nozzleThroatArea, double Pc, double k, double R, double Tc) {
	// Calculates mass flow rate out of the nozzle using equations from the paper
	// Inputs in [m^2], [psi], [], [kJ/kg*K], [k]
	// Outputs in [kg/s]
	Pc = Pc*PSI_TO_PA;
	double mDot = (nozzleThroatArea*Pc*k*sqrt( pow((2 /(k + 1)), ((k + 1) / (k - 1))))) / sqrt(k*R*Tc);
	return mDot;
}

double calcPc(double nozzleThroatArea, double mDotNoz, double k, double R, double Tc) {
	// Calculates mass chamber pressure
	// Inputs in [m^2], [kg/s], [], [kJ/kg*K], [k]
	// Outputs in [psi]
	double Pc = (mDotNoz*sqrt(k*R*Tc))/(nozzleThroatArea * k * sqrt(pow((2 / (k + 1)), ((k + 1) / (k - 1)))));
	Pc = Pc / (PSI_TO_PA);
	return Pc;
}

double massFlowRateInjector(double mDotNoz, double OF_ratio) {
	// Calculates mass flow rate out of the injector
	// Inputs in [kg/s]
	// Outputs in [kg/s]
	double mDotI = mDotNoz*(OF_ratio / (1 + OF_ratio));
	return mDotI;
}

double massFlowRateNozzle(double mDotI, double OF_ratio) {
	// Calculates mass flow rate out of the nozzle
	// Inputs in [kg/s]
	// Outputs in [kg/s]
	double mDotNoz = mDotI*((1 + OF_ratio) / OF_ratio);
	return mDotNoz;
}

double thrustCoefficient(double Patm, double A2, double Pc) {
	// Calculates thrust coefficient for different chamber pressures
	// Inputs in [psi], [m^2], [psi]
	// Ouput is unitless
	Patm = Patm*PSI_TO_PA;
	Pc = Pc*PSI_TO_PA;
	double C12 = 2.23;						//calculate C1,2
	double Cf = C12 - (Patm*A2) / (Pc*At);
	return Cf; 
}

void RPALookup(float Pc, double OF, double &k, double &R, double &Tc) {
	// Find k, R, Tc from table and return them to the program
	// Inputs in [psi], []
	// No output but stores values in variables k, R, Tc
	// Stored in [], [kJ/kg*K], [k]
	RPA_Table CombustionProps = lookUp(Pc, OF, Table_Array);
	k = CombustionProps.k_value;
	R = CombustionProps.R_value;
	Tc = CombustionProps.Chamber_Temperture;
}

