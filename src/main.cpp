#include <iostream>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include "../include/RpaTable.h"
#include "../include/thermo_functions.h"
#include "../include/injector_Model.h"
#include "../include/main.h"
#include "../include/blowdownModel.h"
#include "../include/cfg_file_reader.h"
#include "../include/SolomonModel.h"

using namespace std;

const double PSI_TO_PA = 6894.76; 

// Will be defined later in program
double oxyMass;							// Initial oxidizer mass [kg]
double Pc, PcNew, PcOld;				// Chamber pressure [psi]
double k;								// Heat capacity ratio []
double R;								// Specific gas constant [kJ/kg*k]
double Tc;								// Chamber temperature [k]
double Cf;								// Thrust coefficient []
double Tt, T_Kelvin;					// NOS Tank Temperature
//int PcRound10;						// Casts chamber pressure to integer
double err;								// Used for calculating relative error
double tankPressure;					// Tank absolute pressure [psi]
double ISP;
// Already defined
double mDotNozzle, mDotInjector;		// Mass flow rates at the nozzle and the injector [kg/s]
double mDotInjector_old;				//Prev iteration mass flow rate
double time[1000], thrust[1000];		// Output arrays that give thrust over time

RpaTable RpaTableArray;

void thrustRamp(int dir, double Ramp, double FitThrust, double timestep, int x, ofstream &simfile);

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

int main() 
{
	OptionFileParser UserOptions;
	// Creates parser object, sets its path, and reads the file

	UserOptions.SetPath("./settings.cfg", ":");
	UserOptions.ReadFile();
	RpaTableArray.SetRpaDataFile("RPA_Output_Table.csv");
	RpaTableArray.CreateRpaTable();


	bool SolomonModelFlag = UserOptions.mSolomonFlag;
	double At = UserOptions.mThroatArea;
	double A2 = UserOptions.mExitArea;
	double tankVolume = UserOptions.mOxTankVolume;
	double timeStep = UserOptions.mTimeStep;
	double OF = UserOptions.mOxFuelRatio;
	double RampDownTime = UserOptions.RampDownTime;
	double RampUpTime = UserOptions.RampUpTime;
	
	double liquidMass, vaporizedMass=0;
	ofstream simFile("output.csv");
	simFile << "Time (s), Modeled Liquid Mass  ,  Modeled Chamber Pressure, Modeled Thrust, Modeled Tank Pressure" << '\n'; //setup basic output format
	cout << "Input initial params." << endl;
	bool flag = 1;
	while (flag == 1) {
		cout << "Initial NOS Temperature [Celsius]:  ";
		cin >> Tt;
		T_Kelvin = Tt + 273.15;
		cout << "Abs. Chamber Pressure [PSI]:  ";
		cin >> Pc;
		cout << "Input initial oxidizer mass in [kg]:  ";
		cin >> oxyMass;
		setupTables();
		if (oxyMass > tankVolume*nox_Lrho(Tt)) {
			cout << "Density exceeds sat.Liquid density" << endl;
			cout << "Max expected mass @ T: "<<(tankVolume*nox_Lrho(Tt))<<endl;
			flag = 1;
		}
		else flag = 0;
	}


	////////////////////	
	if (SolomonModelFlag) {
		simFile.close();
		SolomonModel model;
		model.SetInitialProperties(tankVolume, oxyMass, Tt, timeStep);
		model.SetNosDataFile("N20_Neg30_35T.txt");
		model.ImportNosData();
		model.TimeStepLoop();
		cout << "Done!";
		cin >> SolomonModelFlag;
		return 1;
	}
	///////////////////



	liquidMass = oxyMass;//suppress errors from calc below
	tankProps(0, tankVolume, oxyMass, vaporizedMass, liquidMass, T_Kelvin, tankPressure); //calc initial liquid mass and pressure
	cout << "Initial liq Mass:"<<liquidMass;
	for (int x = 0; x < 1000; x++) { //time steps
		
		// Finds all relevant values for thrust
		if (UserOptions.mFlowModel == 2) {
			PcOld = Pc;
			try {
				for (int i = 0; i < 100; i++) {
					interpRPAValues(PcOld, OF, k, R, Tc);
					mDotInjector = interpInjectorModel(Tt, PcOld);
					mDotNozzle = massFlowRateNozzle(mDotInjector, OF);
					PcNew = calcPc(At, mDotNozzle, k, R, Tc);
					err = abs(100 * (PcOld - PcNew) / PcOld);
					PcOld = (PcNew-PcOld)*UserOptions.mConvergenceWeight+PcOld;
					if (err < 4) { Pc = PcNew; err = 100; break; }
					else if (i == 99) { throw runtime_error("PressureCalculatorDiverged"); }
				}
			}
				catch (runtime_error& e)
				{
					cout << e.what() << '\n'; //catch exception, display to user and break loop
				
			}
		}
		else if (UserOptions.mFlowModel==1) { 
			interpRPAValues(Pc, OF, k, R, Tc);
			mDotNozzle = massFlowRate(At, Pc, k, R, Tc);
			mDotInjector = massFlowRateInjector(mDotNozzle, OF); 
		}
		Cf = thrustCoefficient(14.7, A2, Pc, At);

		try { //catch negative flow exception, display to user and break loop
			if (mDotNozzle < 0 || mDotInjector < 0) { throw runtime_error("massFlowNegative"); }
		}
		catch (runtime_error& e)
		{
			cout << e.what() << '\n'; //catch exception, display to user and break loop
			break;
		}
		double deltaM;
		if (UserOptions.mIntegrationType == 1) { 
			deltaM = mDotInjector*timeStep;
		}
		else if (UserOptions.mIntegrationType == 2) {
			deltaM= 0.5 * timeStep * (3.0 * mDotInjector - mDotInjector_old);//addams integration
		}
		oxyMass -= deltaM;
		liquidMass -= deltaM; //assumes pure liquid through injector
		mDotInjector_old = mDotInjector;

		tankProps(timeStep,tankVolume,oxyMass,vaporizedMass,liquidMass,T_Kelvin,tankPressure); //update new Temperature and tank pressure
		
		// Creates outputs for each timestep
		time[x] = x*timeStep;
		thrust[x] = At*(Pc*PSI_TO_PA)*Cf;
		if (x == 1) { thrustRamp(1,RampUpTime, thrust[x], timeStep, x, simFile); }
		ISP += (thrust[x] / (9.81*mDotNozzle));
		cout << "T+" << time[x] << " s =>>> Liquid Mass: " << liquidMass << "kg | Chamber Pressure: " << Pc << " psi | " <<
			"Injector flow rate: " << mDotInjector << " kg/s | Tank Temperature" << T_Kelvin << " K" << endl;
		output(simFile, time[x], liquidMass, Pc, thrust[x], tankPressure);//output to csv
		Tt = T_Kelvin - 273.15;
		if (liquidMass <= 0.05) { cout << "Empty"; thrustRamp(-1, RampDownTime, thrust[x], timeStep, x, simFile); simFile.close(); return 1; };

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

double thrustCoefficient(double Patm, double A2, double Pc, double At) {
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
	RpaTable::RpaDataPoint CombustionProps = RpaTableArray.LookUpRpa(OF, Pc);
	k = CombustionProps.KValue;
	R = CombustionProps.RValue;
	Tc = CombustionProps.ChamberTemperture;
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
		k  = bilinInterp(Pc1, Pc2, OF1, OF2, k1, k2, k3, k4, Pc, OF);
		R  = bilinInterp(Pc1, Pc2, OF1, OF2, R1, R2, R3, R4, Pc, OF);
		Tc = bilinInterp(Pc1, Pc2, OF1, OF2, Tc1, Tc2, Tc3, Tc4, Pc, OF);
	}
	else if (Pc1 != Pc2) {
		RPALookup(Pc1, OF, k1, R1, Tc1);
		RPALookup(Pc2, OF, k2, R2, Tc2);
		k  = linInterp(Pc1, k1, Pc2, k2, Pc);
		R  = linInterp(Pc1, R1, Pc2, R2, Pc);
		Tc = linInterp(Pc1, Tc1, Pc2, Tc2, Pc);
	}
	else if (OF1 != OF2) {
		RPALookup(Pc, OF1, k3, R3, Tc3);
		RPALookup(Pc, OF2, k4, R4, Tc4);
		k = linInterp(OF1, k1, OF2, k2, OF);
		R = linInterp(OF1, R1, OF2, R2, OF);
		Tc = linInterp(OF1, Tc1, OF2, Tc2, OF);
	}
	else {
		RPALookup(Pc, OF, k, R, Tc);
	}
}

double interpInjectorModel(double Tt, double Pc) {
	double Tt1, Tt2, Pc1, Pc2;				// Used to interpolate RPA values
	double mDot1, mDot2, mDot3, mDot4;
	double mDotInjector;
	Tt1 = floor(Tt);
	Tt2 = ceil(Tt);
	Pc1 = floor(Pc / 10) * 10;
	Pc2 = ceil(Pc / 10) * 10;
	if (Pc1 != Pc2 && Tt1 != Tt2) {
		mDot1 = injectorModel(Tt1, Pc1);
		mDot2 = injectorModel(Tt1, Pc2);
		mDot3 = injectorModel(Tt2, Pc1);
		mDot4 = injectorModel(Tt2, Pc2);
		mDotInjector=linInterp(Pc1, linInterp(Tt1, mDot1, Tt2, mDot3, Tt), Pc2, linInterp(Tt1, mDot2, Tt2, mDot4, Tt), Pc);
		//mDotInjector = bilinInterp(Tt1, Tt2, Pc1, Pc2, mDot1, mDot2, mDot3, mDot4, Tt, Pc);
	}
	else if (Tt1 != Tt2) {
		mDot1 = injectorModel(Tt1, Pc);
		mDot2 = injectorModel(Tt2, Pc);
		mDotInjector = linInterp(Tt1, mDot1, Tt2, mDot2, Tt);
	}
	else if (Pc1 != Pc2) {
		mDot1 = injectorModel(Tt, Pc1);
		mDot2 = injectorModel(Tt, Pc2);
		mDotInjector = linInterp(Pc1, mDot1, Pc2, mDot2, Pc);
	}
	else {
		mDotInjector = injectorModel(Tt, Pc);
	}

	return mDotInjector;
}
void thrustRamp(int dir,double Ramp, double FitThrust, double timestep, int x, ofstream &simfile) {
	cout << "Ramping" << dir << endl;
	double A = 1;
	if (dir == -1) { A = FitThrust; }
	int expFactor = exp(log(FitThrust) / Ramp);
	for (int i = 0; i < (floor(Ramp / timestep)); i++) {
		output(simfile, ((x + i)*timestep), 0, 0, (A*pow(expFactor, (i*timestep*dir))));
	}
}
