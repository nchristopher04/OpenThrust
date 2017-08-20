#pragma once
#ifndef HEAD_H
#define HEAD_H
using namespace std;
double At = 10;							// Nozzle throat area
int timeStep = 10;						//ms
float OF = 2.1;							// Oxidizer-fuel ratio assumed constant to start
double Pc, k, R, Tc, nozzleArea, mDot, oxyMass;	// Initializing vars
int time;
#endif // 