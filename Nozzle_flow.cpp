/*
 * rocket_sim.cxx
 * 
 * Copyright 2016 Nicholas Christopher <nicholas@Mackenstein>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#define _USE_MATH_DEFINES // define constants for C++
#include <iostream>
#include <cmath>

using namespace std;

int main()
{
	/*
	 * This is a simple program to approximate apogee of our rocket
	 * The equations can be found here: http://www.rocketmime.com/rockets/qref.html
	 * This assumes some things such as constant thrust. a way to make this more accurate would be to use a time dependant eq'n with our actual thrust distribution
	 * also, the coeff of drag is just a guess. We have no data for our actual coeff of drag.
	 * 
	*/
	double g=9.81;//gravity
	double OF=2.1;//mass of oxidizer / mass of fuel grain
	double k=21.5708/16.7097;//cp/cv
	double T=2163;//degrees K
	double P=310;//psi,absolute chamber pressure
	double R=.3661;//KJ/(kg*k)
	double Dt=0.869;//nozzle throat diameter (inch)
	double At=0.00064516*Dt*Dt*0.25*M_PI;//nozzle throat area (m^2)
	R*=1000;
	P*=6894.76;//Pascals

	double numerator=pow(2/(k+1),(k+1)/(k-1));
	double denom=k*R*T;
	double mdotn=At*P*k*sqrt(numerator/denom);
	double mdoti=mdotn*OF/(1+OF);
		
	cout<<"mass flow rate out of rocket: "<<mdotn<<" kg/s\n";
	cout<<"mass flow rate out of injector: "<<mdoti<<" kg/s\n";

	//cout<<"\ntotal distance (feet) "<<yt<<;
	
//	system("PAUSE");
	//cin>>Mt;//so the code doesn't end (couldnt get system("pause") to work
	return 0;
}

