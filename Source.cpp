#include <head.h>
#include <iostream>
#include <cmath>
#include <utility>
#include <fstream>
#include <stdexcept>
#include "RPA_to_struct.h"
#include "thermo_functions.h"

template <typename Arg, typename... Args>
void output(oftream& simFile, Arg&& arg, Args&&... args) {
	out << forward<Arg>(arg);
	using expander = int[];
	(void)expander {
		0, (void(out << ',' << forward<Args>(args)), 0)...
	};
} //this will be a variadic print function that will output all parameters passed to it to output.csv

int main() {
	ofstream simFile;
	simFile.open("output.csv");
	simFile << "Liquid Mass   ,  Chamber Pressure , Thrust , Mass Flow Rate "; //setup basic output format

	for (int x; x < 1000; x++) { //time steps
		RPALookup(Pc, OF, k, R, Tc);
		mDot = massFlowRate(nozzleArea, Pc, k, R, Tc);
		if (mDot < 0) { 
			throw invalid_argument ("massFlowNegative");
			break; }
		else { oxyMass -= mDot*timeStep; }
		//compute new chamber pressure
		time = x*timeStep;
		output();

		if (oxyMass <= 0.01) {
			cout << "Tank empty";
			break;
		};
	}

}

double massFlowRate(double nozzleArea, double Pc, double k, double R, double Tc) {
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
};

