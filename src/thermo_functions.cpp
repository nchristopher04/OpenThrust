#include <iostream>
#include <cmath>
#include <fstream>
#include "../include/thermo_functions.h"
#include "../include/injector_Model.h"
#include <sstream>
#include <string>
using namespace std;

const double CELS_TO_KELVIN = 273.15;					// Conversion factor from [C] to [K]
const double KPA_TO_BAR = 100;							//conversion from [kPa] to [bar]
const double CRIT_TEMP_NOS = 309.55;					// Critical Temperature of NOS [K]
const double CRIT_PRES_NOS = 7255;						// Critical Pressure of NOS [kPa]
const double R_CONSTANT = 0.008314;						// Universal Gas Constant [kJ/(mol*K)]
const double MM_NOS = 0.044013;							// Molar Mass of Nitrous Oxide [kg/mol]
const double R_SPEC_NOS = R_CONSTANT / (MM_NOS);		// Specific Gas Constant of Nitrous	Oxide [kJ/(kg*K)]

double volume_tank, mass_combined, temperature_tank;	// Initial conditions

ifstream Pfile("N20_100_1000PSI.txt");
ifstream Tfile("N20_Neg30_35T.txt");

double constT[30][100], constP[30][100];
string headT[30], headP[30]; //global scope


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
double nox_vp(double T_Celcius)
{
	double pVap= data_grab("Pressure (psia)", T_Celcius, "T", constT, headT);
	return(pVap);
}
/* Nitrous liquid Enthalpy (Latent heat) of vaporisation, J/kg */
double nox_enthV(double T_Celcius)
{
	double Hvap = 0;
	double hl = data_grab("Enthalpy (l, kJ/kg)", T_Celcius, "T", constT, headT);//liquid enthalpy at 2
	double hg = data_grab("Enthalpy (v, kJ/kg)", T_Celcius, "T", constT, headT);//gas enthalpy at 2
	Hvap = (hg - hl) * 1000;
	return(Hvap);
}

double nox_Lrho(double T_Celcius)
{
	double rho = data_grab("Density (l, kg/m3)", T_Celcius, "T", constT, headT);//NIST liquid density
	return(rho);
}
/* Nitrous oxide saturated vapour density, kg/m3 */
double nox_Vrho(double T_Celcius)
{	
	double rho=data_grab("Density (v, kg/m3)", T_Celcius, "T", constT, headT);//NIST VAPOR density
	return(rho);
}

double nox_Cp(double T_Celcius)
{
	double Cp = data_grab("Cp (l, J/g*K)", T_Celcius, "T", constT, headT);
	return(Cp);
}

double linInterp(double x1, double y1, double x2, double y2, double x) {
	double y;
	// Enter in two coordinates (x1,y1) and (x2,y2)
	// as well as an x value that should be somewhere
	// close or in between
	if (x2 == x1) { //prevent division by zero
			return y1;
	}
	else y = y1 + (x - x1)*(y2 - y1) / (x2 - x1);
	return y;
}

void NoxPropertiesT(NIST_Table *props, double T_Kelvin){
	T_Kelvin -= 273.15;
	double T1 = ceil(T_Kelvin);
	double T2 = floor(T_Kelvin);
	props->Cp = linInterp(T1, nox_Cp(T1), T2, nox_Cp(T2), T_Kelvin);
	props->enthV = linInterp(T1, nox_enthV(T1), T2, nox_enthV(T2), T_Kelvin);
	props->Lrho = linInterp(T1, nox_Lrho(T1), T2, nox_Lrho(T2), T_Kelvin);
	props->Vrho = linInterp(T1, nox_Vrho(T1), T2, nox_Vrho(T2), T_Kelvin);
	props->pVap = linInterp(T1, nox_vp(T1), T2, nox_vp(T2), T_Kelvin);
	return;
}
void setupTables() { 
	data_gather(Pfile, constP, headP);
	data_gather(Tfile, constT, headT);
}