#ifndef  SOURCE_H
#define SOURCE_H

double massFlowRateInjector(double nozzleFlow, double OF_ratio);
double massFlowRate(double nozzleArea, double Pc, double k, double R, double Tc);
void RPALookup(float Pc, float OF, double& k, double& R, double& Tc); //Forward declarations
double thrustCoefficient(double Patm, double A2, double Pc);


#endif // ! SOURCE_H
