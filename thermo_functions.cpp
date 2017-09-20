#include <iostream>
#include <cmath>
#include "thermo_functions.h"

const double CELS_TO_KELVIN = 273.15;					// Conversion factor from [C] to [K]

const double CRIT_TEMP_NOS = 36.4;						// Critical Temperature of NOS [C]
const double CRIT_PRES_NOS = 7255;						// Critical Pressure of NOS [kPa]
const double R_CONSTANT = 0.008314;						// Universal Gas Constant [kJ/(mol*K)]
const double MM_NOS = 0.044013;							// Molar Mass of Nitrous Oxide [kg/mol]
const double R_SPEC_NOS = R_CONSTANT / (MM_NOS);		// Specific Gas Constant of Nitrous	Oxide [kJ/(kg*K)]

double volume_tank, mass_combined, temperature_tank;	// Initial conditions
														//NOS properties from Modelling the Nitrous Run tank Emptying
const float pCrit = 72.51f; /* critical pressure, Bar Abs */
const float rhoCrit = 452.0f; /* critical density, kg/m3 */
const float ZCrit = 0.28f; /* critical compressibility factor */
const float gamma = 1.3f; /* average over subcritical range */


double reduced_temperature(double temperature_tank)
{
	// Returns reduced temperature of NOS
	// Input should be in Celsius
	double Tr = (temperature_tank + CELS_TO_KELVIN) / (CRIT_TEMP_NOS + CELS_TO_KELVIN);
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

/* Nitrous oxide vapour pressure, Bar */
double nox_vp(double T_Celcius)
{
	const float p[4] = { 1.0f, 1.5f, 2.5f, 5.0f };
	const float b[4] = { -6.71893f, 1.35966f, -1.3779f, -4.051f };
	double Tr = reduced_temperature(T_Celcius);
	float rab = 1.0 - Tr;
	float shona = 0.0;
	for (int dd = 0; dd < 4; dd++)
		shona += b[dd] * pow(rab, p[dd]);
	double bob = pCrit * exp((shona / Tr));
	return(bob);
}