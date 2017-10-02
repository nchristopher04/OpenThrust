#include <iostream>
#include <cmath>
#include "thermo_functions.h"

const double CELS_TO_KELVIN = 273.15;					// Conversion factor from [C] to [K]
const double KPa_TO_BAR = 100;							//conversion from [kPa] to [bar]
const double CRIT_TEMP_NOS = 309.55;					// Critical Temperature of NOS [K]
const double CRIT_PRES_NOS = 7255;						// Critical Pressure of NOS [kPa]
const double R_CONSTANT = 0.008314;						// Universal Gas Constant [kJ/(mol*K)]
const double MM_NOS = 0.044013;							// Molar Mass of Nitrous Oxide [kg/mol]
const double R_SPEC_NOS = R_CONSTANT / (MM_NOS);		// Specific Gas Constant of Nitrous	Oxide [kJ/(kg*K)]

double volume_tank, mass_combined, temperature_tank;	// Initial conditions

double reduced_temperature(double Temp_Kelvin)
{
	// Returns reduced temperature of NOS
	// Input should be in Celsius
	double Tr = (Temp_Kelvin) / (CRIT_TEMP_NOS);
	return Tr;
}

double reduced_pressure(double pressure_tank)
{
	// Returns reduced pressure of NOS
	// Input should be in kPa
	double Pr = pressure_tank / CRIT_PRES_NOS;
	return Pr;
}

double density_combined(double mass_vapor, double mass_liquid, double volume_tank)
{
	// Assumes uniform mixture of vapor and liquid in ox tank
	// Calculates combined density of oxidizer tank fluids
	// Inputs should be in kg and m^3
	// Output in kg/m^3
	// Equation 3.1 in Soloman
	double density;
	density = (mass_vapor + mass_liquid) / volume_tank;
	return density;
}

double fluid_quality(double density_combined, double density_liquid, double density_vapor)
{
	// Calculates the quality of the fluid
	// Inputs should be in kg/m^3
	// Output is unitless
	// Equation 3.2 in Soloman
	double quality;
	quality = (density_vapor / density_combined)*((density_liquid - density_combined) / (density_liquid - density_vapor));
	return quality;
}

double specific_enthalpy_combined(double specific_enthalpy_vapor, double specific_enthalpy_liquid, double fluid_quality)
{
	// Assumes uniform mixture of vapor and liquid in ox tank
	// Calculates the combined specific enthalpy of oxidizer tank fluids
	// Inputs should be in kJ/kg for enthalpies
	// Outputs in kJ/kg
	double sp_enth_c;
	sp_enth_c = specific_enthalpy_vapor*fluid_quality + (1-fluid_quality)*specific_enthalpy_liquid;
	return sp_enth_c;
}

double specific_entropy_combined(double specific_entropy_vapor, double specific_entropy_liquid, double fluid_quality)
{
	// Assumes uniform mixture of vapor and liquid in ox tank
	// Calculates the combined specific entropy of oxidizer tank fluids
	// Inputs should be in kJ/(kg*K) for entropies
	// Outputs in kJ/(kg*K)
	double sp_entr_c;
	sp_entr_c = specific_entropy_vapor*fluid_quality + (1 - fluid_quality)*specific_entropy_liquid;
	return sp_entr_c;
}

double total_enthalpy(double specific_enthalpy, double mass)
{
	// Calculates total enthalpy
	// Mass input should be in kg
	// Specific Enthalpy input should be in kJ/kg
	// Outputs in kJ
	double tot_enth;
	tot_enth = mass*specific_enthalpy;
	return tot_enth;
}

/* Nitrous oxide vapour pressure, kPa */
double nox_vp(double T_Kelvin)
{
	const float p[4] = { 1.0f, 1.5f, 2.5f, 5.0f };
	const float b[4] = { -6.71893f, 1.35966f, -1.3779f, -4.051f };
	double Tr = reduced_temperature(T_Kelvin);
	double rab = 1.0 - Tr;
	double shona = 0.0;
	for (int dd = 0; dd < 4; dd++)
		shona += b[dd] * pow(rab, p[dd]);
	double  Pvap= (CRIT_PRES_NOS/KPa_TO_BAR) * exp((shona / Tr));
	return(Pvap*KPa_TO_BAR); //should work without the conversions but I'm not sure so I'll leave converting through to bar
}
/* Nitrous liquid Enthalpy (Latent heat) of vaporisation, J/kg */
double nox_enthV(double T_Kelvin)
{
	double Hvap = 0;
	
	double hl = data_grab("Enthalpy (l, kJ/kg)", P2, "P", constP, headP);//liquid enthalpy at 2
	double hg = data_grab("Enthalpy (v, kJ/kg)", P2, "P", constP, headP);//gas enthalpy at 2
	
	Hvap = (hg - hl) * 1000;
	/*
	const float bL[5] = { -200.0f, 116.043f, -917.225f, 794.779f, -589.587f };
	const float bV[5] = { -200.0f, 440.055f, -459.701f, 434.081f, -485.338f };
	double shonaL, shonaV;
	double Tr = reduced_temperature(T_Kelvin);
		double rab;
	rab = 1.0 - Tr;
	shonaL = bL[0];
	shonaV = bV[0];
	for (int dd = 1; dd < 5; dd++)
	{
		shonaL += bL[dd] * pow(rab, (dd / 3.0));  saturated liquid enthalpy 
		shonaV += bV[dd] * pow(rab, (dd / 3.0));  saturated vapour enthalpy 
	}
double Hvap = (shonaV - shonaL) * 1000.0;  net during change from liquid to vapour */
	return(Hvap);
}

double nox_Lrho(double T_Kelvin)
{
	const float b[4] = { 1.72328f, -0.8395f, 0.5106f, -0.10412f };
	double Tr = reduced_temperature(T_Kelvin);
	double rab = (1.0 / Tr) - 1.0;
	double shona = 0.0;
	for (int i = 0; i < 5; i++)
		shona += b[i] * pow(rab, ((i + 1) / 3.0));
	double rho = rhoCrit * exp(shona);
	return(rho);
}
/* Nitrous oxide saturated vapour density, kg/m3 */
double nox_Vrho(double T_Kelvin)
{
	const float b[5] = { -1.009f, -6.28792f, 7.50332f, -7.90463f, 0.629427f };
	double Tr = reduced_temperature(T_Kelvin);
	double rab = (1.0 / Tr) - 1.0;
	double shona = 0.0;
	for (int i = 0; i < 5; i++)
		shona += b[i] * pow(rab, ((i + 1) / 3.0));
	double rho = rhoCrit * exp(shona);
	return(rho);
}

double nox_Cp(double T_Kelvin)
{
	const float b[5] = { 2.49973f, 0.023454f, -3.80136f, 13.0945f, -14.518f };
	double Tr = reduced_temperature(T_Kelvin);
	double rab = 1.0 - Tr;
	double shona = 1.0 + b[1] / rab;
	for (int i = 1; i < 4; i++)
		shona += b[(i + 1)] * pow(rab, i);
	double heatCapacity = b[0] * shona * 1000.0; /* convert from KJ to J */
	return(heatCapacity);
}