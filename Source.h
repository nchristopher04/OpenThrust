#ifndef  SOURCE_H
#define SOURCE_H

double massFlowRate(double nozzleThroatArea, double Pc, double k, double R, double Tc);
double calcPc(double nozzleThroatArea, double mDotNoz, double k, double R, double Tc);
double massFlowRateInjector(double mDotNoz, double OF_ratio);
double massFlowRateNozzle(double mDotI, double OF_ratio);
double thrustCoefficient(double Patm, double A2, double Pc);
void RPALookup(double Pc, double OF, double &k, double &R, double &Tc);
double nox_vp(double T_Celcius);
double linInterp(double x1, double y1, double x2, double y2, double x);
double bilinInterp(double x1, double x2, double y1, double y2, double P11, double P12, double P21, double P22, double a, double b);
void interpRPAValues(double Pc, double OF, double &k, double &R, double &Tc);
double interpInjectorModel(double Tt, double Pc);


#endif // ! SOURCE_H
