#pragma once
#ifndef BLOWDOWNMODEL_H   // To make sure you don't declare the function more than once by including the header multiple times.
#define BLOWDOWNMODEL_H

double tankProps(double oxyMass, double Pc, double &T_Kelvin, double &TankPressure);
#endif