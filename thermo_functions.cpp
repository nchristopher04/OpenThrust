#include <iostream>
#include <math.h>
#include "thermo_functions.h"

static double CELS_TO_KELVIN = 273.15;					// Conversion factor from [C] to [K]

static double CRIT_TEMP_NOS = 36.4;						// Critical Temperature of NOS [C]
static double CRIT_PRES_NOS = 7255;						// Critical Pressure of NOS [kPa]
static double R = 0.008314;								// Universal Gas Constant [kJ/(mol*K)]
static double MM_NOS = 0.044013;						// Molar Mass of Nitrous Oxide [kg/mol]
static double R_SPEC_NOS = R / (MM_NOS);				// Specific Gas Constant of Nitrous	Oxide [kJ/(kg*K)]

double volume_tank, mass_combined, temperature_tank;	// Initial conditions

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
