#include <cmath>
#include "blowdownModel.h"
#include "thermo_functions.h"
#include "Source.h"
using namespace std;

struct Faults {
	char tempFault;
	bool vaportFault = false;
}blowdownModel;

void tankProps(double oxyMass, double &vaporizedMass_prev,double &liquidMass_prev, double &T_Kelvin, double &TankPressure) {
	double liquidMass, vaporMass, deltaQ,vaporizedMass;
	
	deltaQ = vaporizedMass_prev*nox_enthV(T_Kelvin);
	T_Kelvin-=(deltaQ / (liquidMass * nox_Cp(T_Kelvin)));  // can we define a heat capacity for the whole system?
	if (T_Kelvin < (183))
	{
		T_Kelvin = (183); /* lower limit */
		 blowdownModel.tempFault = 'L';
	}
	else if (T_Kelvin >309.15) //upper limit 36C
	{
		T_Kelvin = 309.15;
		blowdownModel.tempFault = 'H';
	}
	
	TankPressure = nox_vp(T_Kelvin);

	liquidMass = (tankVolume - (oxyMass / nox_Vrho(T_Kelvin))) / ((1.0 / nox_Lrho(T_Kelvin)) - (1.0 / nox_Vrho(T_Kelvin)));
	vaporMass = oxyMass - liquidMass;
	
	vaporizedMass = liquidMass_prev - liquidMass;
	if (vaporizedMass < 0) {
		blowdownModel.vaportFault = true;
	}
	double lagged = (timeStep / 0.15) * (vaporizedMass - lagged) + lagged; // 1st-order lag
	vaporizedMass_prev = lagged; //to be used in next iteration
	liquidMass_prev = liquidMass; //set liquidmass_old for next iteration to current, get pushed outside of this function to main iterator
	return;
}

