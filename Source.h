#ifndef  SOURCE_H
#define SOURCE_H

double massFlowRate(double nozzleThroatArea, double Pc, double k, double R, double Tc);
double calcPc(double nozzleThroatArea, double mDotNoz, double k, double R, double Tc);
double massFlowRateInjector(double mDotNoz, double OF_ratio);
double massFlowRateNozzle(double mDotI, double OF_ratio);
double thrustCoefficient(double Patm, double A2, double Pc);
void RPALookup(float Pc, double OF, double &k, double &R, double &Tc);
double tankProps(double oxyMass, double Pc, double &Temp, double &TankPressure);
double nox_vp(double T_Celcius);

#endif // ! SOURCE_H
