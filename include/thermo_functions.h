#ifndef THERMO_FUNCTIONS_H
#define THERMO_FUNCTIONS_H

double reduced_temperature(double);
double reduced_pressure(double);
double density_combined(double, double, double);
double fluid_quality(double, double, double);
double specific_enthalpy_combined(double, double, double);
double specific_entropy_combined(double, double, double);
double total_enthalpy(double, double);
double nox_enthV(double);
double nox_vp(double);
double nox_Lrho(double);
double nox_Vrho(double);
double nox_Cp(double);
double linInterp(double x1, double y1, double x2, double y2, double x);

struct NIST_Table
{
	double pVap;
	double enthV;
	double Lrho;
	double Vrho;					// [kg/m^3]
	double Cp;
	double Quality;
};
void NoxPropertiesT(NIST_Table *props,double T_Kelvin);
void setupTables();
#endif