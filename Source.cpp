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
//int PcRound10;							// Casts chamber pressure to integer
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
		interpRPAValues(Pc, OF, k, R, Tc);
		mDotNozzle = massFlowRate(At, Pc, k, R, Tc);
		if (MainX.flowModel == 2) {
			PcOld = Pc;
			for (int i = 0; i < 100; i++) {
				//PcRound10 = PcOld; //cast chamber pressure to int.
				//PcRound10 += 5;
				//PcRound10 -= PcRound10 % 10;
				interpRPAValues(PcOld, OF, k, R, Tc);
				mDotInjector = injectorModel(T1, PcOld);
				mDotNozzle = massFlowRateNozzle(mDotInjector, OF);
				PcNew = calcPc(At, mDotNozzle, k, R, Tc);
				err = abs(100 * (PcOld - PcNew) / PcOld);
				PcOld = PcNew;
				if (err < 5) { Pc = PcNew; err = 100; break; }
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

void RPALookup(double Pc, double OF, double &k, double &R, double &Tc) {
	// Find k, R, Tc from table and return them to the program
	// Inputs in [psi], []
	// No output but stores values in variables k, R, Tc
	// Stored in [], [kJ/kg*K], [k]
	RPA_Table CombustionProps = lookUp(Pc, OF, Table_Array);
	k = CombustionProps.k_value;
	R = CombustionProps.R_value;
	Tc = CombustionProps.Chamber_Temperture;
}

double linInterp(double x1, double y1, double x2, double y2, double x) {
	double y;
	// Enter in two coordinates (x1,y1) and (x2,y2)
	// as well as an x value that should be somewhere
	// close or in between
	y = y1 + (x - x1)*(y2 - y1) / (x2 - x1);
	return y;
}

double bilinInterp(double x1, double x2, double y1, double y2,
	double P11, double P12, double P21, double P22,
	double a, double b) {
	// Enter in the values of a function f(x,y)
	// where P1 = f(x1,y1), P2 = f(x1,y2),
	// P3 = f(x2,y1), P4 = f(x2,y2). Also enter in
	// an a and b value to interpolate. The a value
	// should be in between x1 and x2 and the b value
	// should be in betweeen y1 and y2. Returns the
	// interpolated value of f(a,b).
	double H1, H2, H3, H4, H5;
	double P;
	H1 = 1 / ((x2 - x1)*(y2 - y1));
	H2 = P11*(x2 - a)*(y2 - b);
	H3 = P12*(a - x1)*(y2 - b);
	H4 = P21*(x2 - a)*(b - y1);
	H5 = P22*(a - x1)*(b - y1);
	P = H1*(H2 + H3 + H4 + H5);
	return P;
}

void interpRPAValues(double Pc, double OF, double &k, double &R, double &Tc) {
	double OF1Pc1, OF1Pc2, Of2Pc1, Of2Pc2;	// Used to interpolate RPA values
	double OF1, OF2, Pc1, Pc2;				// Used to interpolate RPA values
	double k1, k2, k3, k4;					// Used to interpolate k value
	double R1, R2, R3, R4;					// Used to interpolate R value
	double Tc1, Tc2, Tc3, Tc4;				// Used to interpolate Tc value
	OF1 = floor(OF * 10) / 10;
	OF2 = ceil(OF * 10) / 10;
	Pc1 = floor(Pc / 10) * 10;
	Pc2 = ceil(Pc / 10) * 10;
	if (Pc1 != Pc2 && OF1 != OF2) {
		RPALookup(Pc1, OF1, k1, R1, Tc1);
		RPALookup(Pc1, OF2, k2, R2, Tc2);
		RPALookup(Pc2, OF1, k3, R3, Tc3);
		RPALookup(Pc2, OF2, k4, R4, Tc4);
		k = bilinInterp(Pc1, Pc2, OF1, OF2, k1, k2, k3, k4, Pc, OF);
		R = bilinInterp(Pc1, Pc2, OF1, OF2, R1, R2, R3, R4, Pc, OF);
		Tc = bilinInterp(Pc1, Pc2, OF1, OF2, Tc1, Tc2, Tc3, Tc4, Pc, OF);
	}
	else if (Pc1 != Pc2) {
		RPALookup(Pc1, OF, k1, R1, Tc1);
		RPALookup(Pc2, OF, k2, R2, Tc2);
	}
	else if (OF1 != OF2) {
		RPALookup(Pc, OF1, k3, R3, Tc3);
		RPALookup(Pc, OF2, k4, R4, Tc4);
	}
	else {
		RPALookup(Pc, OF, k, R, Tc);
	}
}