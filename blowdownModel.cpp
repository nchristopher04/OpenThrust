#include <cmath>
#include "blowdownModel.h"
#include "thermo_functions.h"
#include "Source.h"

struct Faults {
	char tempFault;
	bool vaporFault = false;
}blowdownModel;

void tankProps(double timeStep, double tankVolume, double oxyMass, double &vaporizedMass_prev, double &liquidMass_prev, double &T_Kelvin, double &TankPressure) {
	double liquidMass, vaporMass, deltaQ, vaporizedMass;
		double lagged=0.0;

	deltaQ = vaporizedMass_prev*nox_enthV(T_Kelvin);
	T_Kelvin -= (deltaQ / (oxyMass * nox_Cp(T_Kelvin)));  // define a heat capacity for the whole system
	if (T_Kelvin < (183))
	{
		T_Kelvin = (183); // lower limit -90C
		blowdownModel.tempFault = 'L';
	}
	else if (T_Kelvin > 309.15) //upper limit 36C
	{
		T_Kelvin = 309.15;
		blowdownModel.tempFault = 'H';
	}

	TankPressure = nox_vp(T_Kelvin);

	liquidMass = (tankVolume - (oxyMass / nox_Vrho(T_Kelvin))) / ((1.0 / nox_Lrho(T_Kelvin)) - (1.0 / nox_Vrho(T_Kelvin))); //check this for validity
	vaporMass = oxyMass - liquidMass;

	vaporizedMass = liquidMass_prev - liquidMass;
	if (vaporizedMass < 0) {
		blowdownModel.vaporFault = true;
	}
	lagged = (timeStep / 0.15) * (vaporizedMass - lagged) + lagged; // 1st-order lag
	vaporizedMass_prev = lagged; //to be used in next iteration
	liquidMass_prev = liquidMass; //set liquidmass_old for next iteration to current, get pushed outside of this function to main iterator
	return;
}

