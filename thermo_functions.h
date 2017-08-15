#ifndef THERMO_FUNCTIONS_H
#define THERMO_FUNCTIONS_H

double reduced_temperature(double);
double reduced_pressure(double);
double density_combined(double, double, double);
double fluid_quality(double, double, double);
double specific_enthalpy_combined(double, double, double);
double specific_entropy_combined(double, double, double);
double total_enthalpy(double, double);


#endif