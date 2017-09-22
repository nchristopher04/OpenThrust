#ifndef  SOURCE_H
#define SOURCE_H

double massFlowRate(double nozzleThroatArea, double Pc, double k, double R, double Tc);
double calcPc(double nozzleThroatArea, double mDotNoz, double k, double R, double Tc);
double massFlowRateInjector(double mDotNoz, double OF_ratio);
double massFlowRateNozzle(double mDotI, double OF_ratio);
double thrustCoefficient(double Patm, double A2, double Pc);
void RPALookup(float Pc, double OF, double &k, double &R, double &Tc);
double nox_vp(double T_Celcius);

double tankVolume, tankPressure; //these need to be accessible by blowdownModel and Source
const float timeStep = 0.1;				// [s]
#endif // ! SOURCE_H
