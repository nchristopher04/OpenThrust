#include <stdafx.h>
#include <iostream>
#include <cmath>
using namespace std;
double tankPressure,
double oxyMass;
double At = 10; //nozzle throat area
double timeStep=0.1;
float OF = 2.1; //Oxidizer-fuel ratio assumed constant to start
int main() {
	for (int x; x < 1000; x++) { //time steps
		RPALookup(Pc);
		mDot = massFlowRate(nozzleArea, Pc, k, R, Tc);
		if (mDot < 0) throw "massFlowNegative"; break
		else	oxyMass -= mDot*timeStep;
		//compute new chamber pressure
		time = x*timeStep;
		cout << time << oxyMass << Pc
			if (oxyMass <= 0.01) { break };
	}

}

double massFlowRate(double nozzleArea, double Pc, float k, const float R, double Tc) {
	mDot = nozzleArea*Pc*k*sqrt(2(k + 1) ^ (k + 1 / (k - 1))) / sqrt(k*R*Tc);
		return mDot;
}

double thrustCoefficient(double Patm, double A2, double Pc) {
	double C12 = 977;//calculate C1,2
	Cf=C12-Patm*A2/(Pc*At)
		return Cf
}

void RPALookup(float Pc, &double k, &double R, &double Tc) {
	//find k, R, Tc from table and return them to the program
}