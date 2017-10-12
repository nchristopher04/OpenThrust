#include <cmath>
#include "blowdownModel.h"
#include "thermo_functions.h"
#include "main.h"
#include <iostream>
#include <stdexcept>
using namespace std;

struct Faults {
	char tempFault;
	bool vaporFault = false;
}blowdownModel;

void tankProps(double timeStep, double tankVolume, double oxyMass, double &vaporizedMass_prev, double &liquidMass_prev, double &T_Kelvin, double &TankPressure) {
	double liquidMass = 0, vaporMass = 0, deltaQ = 0, vaporizedMass = 0;
		double lagged=0.0;
	if (T_Kelvin < (183))
	{
		std::cout << "TempFault L :" << T_Kelvin;
		T_Kelvin = (183); // lower limit -90C
		blowdownModel.tempFault = 'L';
		
	}
	else if (T_Kelvin > 309.15) //upper limit 36C
	{
		std::cout << "TempFault H :" << T_Kelvin;
		T_Kelvin = 309.15;
		blowdownModel.tempFault = 'H';
	}
	deltaQ = vaporizedMass_prev*nox_enthV(T_Kelvin);
	T_Kelvin -= (deltaQ / (oxyMass * nox_Cp(T_Kelvin)*1000));  // define a heat capacity for the whole system

	TankPressure = nox_vp(T_Kelvin);
	double liquidDensity = nox_Lrho(T_Kelvin);
	double vaporDensity = nox_Vrho(T_Kelvin);
	double combinedDensity = (oxyMass / tankVolume);
	double NoxQuality = fluid_quality(combinedDensity, liquidDensity, vaporDensity);
	try{
	liquidMass = (1-NoxQuality)*oxyMass;
	vaporMass = oxyMass - liquidMass;

	vaporizedMass = (liquidMass_prev - liquidMass);
	if (vaporizedMass < 0 || NoxQuality<0) {
		blowdownModel.vaporFault = true;
		throw (invalid_argument("vaporFault"));
	}
	lagged = (timeStep / 0.15) * (vaporizedMass - lagged) + lagged; // 1st-order lag
	vaporizedMass_prev = lagged; //to be used in next iteration
	liquidMass_prev = liquidMass; //set liquidmass_old for next iteration to current, get pushed outside of this function to main iterator
	}
	catch (invalid_argument& e) {
		cout << e.what();
	}
	return;
}

