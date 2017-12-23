#include <cmath>
#include "../include/blowdownModel.h"
#include "../include/thermo_functions.h"
#include "../include/main.h"
#include <iostream>
#include <stdexcept>
using namespace std;

void blowdownModel::tankProps(double timeStep, double tankVolume, double oxyMass, double &vaporizedMass_prev, double &liquidMass_prev, double &T_Kelvin, double &TankPressure) {
	double liquidMass = 0, vaporMass = 0, deltaQ = 0, vaporizedMass = 0;
		static double lagged=0.0;
		NIST_Table NoxTable;
		NoxPropertiesT(&NoxTable, T_Kelvin);//compute new nos properties
	if (T_Kelvin < (183))
	{
		std::cout << "TempFault L :" << T_Kelvin;
		T_Kelvin = (183); // lower limit -90C
		fault.tempFault = 'L';
		
	}
	else if (T_Kelvin > 309.15) //upper limit 36C
	{
		std::cout << "TempFault H :" << T_Kelvin;
		T_Kelvin = 309.15;
		fault.tempFault = 'H';
	}
	deltaQ = vaporizedMass_prev*NoxTable.enthV;
	T_Kelvin -= (deltaQ / (liquidMass_prev * NoxTable.Cp*1000));  // define a heat capacity for the whole system

	TankPressure = NoxTable.pVap;
	double combinedDensity = (oxyMass / tankVolume);
	NoxTable.Quality = fluid_quality(combinedDensity, NoxTable.Lrho, NoxTable.Vrho);
	liquidMass = (1-NoxTable.Quality)*oxyMass;
	vaporMass = oxyMass - liquidMass;
	vaporizedMass = liquidMass_prev - liquidMass;
	if ( NoxTable.Quality<0|| NoxTable.Quality>1 || vaporizedMass>oxyMass) {
		fault.vaporFault = true;
		cout << "VaporFault"<< endl;
		system("Pause");
	}
	lagged = (0.005) * (vaporizedMass - lagged) + lagged; // 1st-order lag 
	vaporizedMass_prev = lagged; //to be used in next iteration
	liquidMass_prev = liquidMass; //set liquidmass_old for next iteration to current, get pushed outside of this function to main iterator
	return;
}

double blowdownModel::ventLosses(double Tankpressure, double ambPressure, double ventArea) {
	double mVented;
	return mVented;
}