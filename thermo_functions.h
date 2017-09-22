#ifndef THERMO_FUNCTIONS_H
#define THERMO_FUNCTIONS_H

double reduced_temperature(double);
double reduced_pressure(double);
double density_combined(double, double, double);
double fluid_quality(double, double, double);
double specific_enthalpy_combined(double, double, double);
double specific_entropy_combined(double, double, double);
double total_enthalpy(double, double);
double nox_enthV(double T_Kelvin);
double nox_vp(double T_Kelvin);
double nox_Lrho(double T_Kelvin);
double nox_Vrho(double T_Kelvin);
double nox_Cp(double T_Kelvin);
#endif