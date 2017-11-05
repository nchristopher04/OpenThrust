#include "../include/SolomonModel.h"
#include "../include/injector_Model.h"
#include "../libs/alglib/ap.h"
#include "../libs/alglib/alglibinternal.h"
#include "../libs/alglib/interpolation.h"
#include "../libs/alglib/solvers.h"
#include "../include/injector_Model.h"

using namespace std;


void SolomonModel::SetInitialProperties(double oxTankVolume, double initialMassOxidizer, double initialTemperature, double timeStep)
{
	// Assume uniform mixture of fluid
	mInitialFluidDensity = (initialMassOxidizer/oxTankVolume);
	mInitialMassOxidizer = initialMassOxidizer;
	mOxTankVolume = oxTankVolume;
	mInitialTemperature = initialTemperature + 273.15;
	mTimeStep = timeStep;

	mFluidDensity = mInitialFluidDensity;
	mFluidMass = mInitialMassOxidizer;
	mNosTemperature = mInitialTemperature;
}

void SolomonModel::SetNosDataFile(const char* temperatureIncrementedFile)
{
	mTemperatureIncrementedFile = temperatureIncrementedFile;
}

void SolomonModel::ImportNosData()
{

	alglib::read_csv(mTemperatureIncrementedFile, '\t', 1 , NistFile);

	int tableLength = NistFile.rows();
	mNistTableLength = tableLength;

	// Temperature Column [K]
	double celsiusToKelvin = 273.15;
	NIST.Temperature.setlength(tableLength);
	for (int i = 0; i < tableLength; i++) {
		NIST.Temperature[i] = NistFile[i][0] + celsiusToKelvin;
	}

	// Pressure Column [psia]
	NIST.Pressure.setlength(tableLength);
	for (int i = 0; i < tableLength; i++) {
		NIST.Pressure[i] = NistFile[i][1];
	}

	// Liquid Density Column [kg/m3]
	NIST.LiqDensity.setlength(tableLength);
	for (int i = 0; i < tableLength; i++) {
		NIST.LiqDensity[i] = NistFile[i][2];
	}

	// Liquid Enthalpy Column [kJ/kg]
	NIST.LiqEnthalpy.setlength(tableLength);
	for (int i = 0; i < tableLength; i++) {
		NIST.LiqEnthalpy[i] = NistFile[i][5];
	}

	// Liquid Entropy Column [J/g*K]
	NIST.LiqEntropy.setlength(tableLength);
	for (int i = 0; i < tableLength; i++) {
		NIST.LiqEntropy[i] = NistFile[i][6];
	}

	// Vapor Density Column [kg/m3]
	NIST.VapDensity.setlength(tableLength);
	for (int i = 0; i < tableLength; i++) {
		NIST.VapDensity[i] = NistFile[i][14];
	}

	// Vapor Enthalpy Column [kJ/kg]
	NIST.VapEnthalpy.setlength(tableLength);
	for (int i = 0; i < tableLength; i++) {
		NIST.VapEnthalpy[i] = NistFile[i][17];
	}

	// Vapor Entropy Column [J/g*K]
	NIST.VapEntropy.setlength(tableLength);
	for (int i = 0; i < tableLength; i++) {
		NIST.VapEntropy[i] = NistFile[i][18];
	}

	mMinTemperature = NIST.Temperature(0);
	mMaxTemperature = NIST.Temperature(tableLength - 1);
}

SolomonModel::NistDataPoint SolomonModel::SetDataPoint(double oxTankTemperature, double effectiveDensity)
{
	NistDataPoint out;
	// Inputs should be in [K] and [kg/m3]
	if (oxTankTemperature < mMinTemperature || oxTankTemperature > mMaxTemperature)
	{
		throw(runtime_error("Temperature out of bounds"));
	}

	// Creates interpolation objects
	alglib::spline1dinterpolant pressureSpline;
	alglib::spline1dinterpolant liqDensitySpline;
	alglib::spline1dinterpolant liqEnthalpySpline;
	alglib::spline1dinterpolant liqEntropySpline;
	alglib::spline1dinterpolant vapDensitySpline;
	alglib::spline1dinterpolant vapEnthalpySpline;
	alglib::spline1dinterpolant vapEntropySpline;

	// Builds the splines into the objects
	alglib::spline1dbuildcubic(NIST.Temperature, NIST.Pressure, pressureSpline);
	alglib::spline1dbuildcubic(NIST.Temperature, NIST.LiqDensity, liqDensitySpline);
	alglib::spline1dbuildcubic(NIST.Temperature, NIST.LiqEnthalpy, liqEnthalpySpline);
	alglib::spline1dbuildcubic(NIST.Temperature, NIST.LiqEntropy, liqEntropySpline);
	alglib::spline1dbuildcubic(NIST.Temperature, NIST.VapDensity, vapDensitySpline);
	alglib::spline1dbuildcubic(NIST.Temperature, NIST.VapEnthalpy, vapEnthalpySpline);
	alglib::spline1dbuildcubic(NIST.Temperature, NIST.VapEntropy, vapEntropySpline);

	// Interpolates using objects and inputted temperatures
	out.Temperature = oxTankTemperature;
	out.Pressure = alglib::spline1dcalc(pressureSpline, oxTankTemperature);
	out.LiqDensity = alglib::spline1dcalc(liqDensitySpline, oxTankTemperature);
	out.LiqEnthalpy = alglib::spline1dcalc(liqEnthalpySpline, oxTankTemperature);
	out.LiqEntropy = alglib::spline1dcalc(liqEntropySpline, oxTankTemperature);
	out.VapDensity = alglib::spline1dcalc(vapDensitySpline, oxTankTemperature);
	out.VapEnthalpy = alglib::spline1dcalc(vapEnthalpySpline, oxTankTemperature);
	out.VapEntropy = alglib::spline1dcalc(vapEntropySpline, oxTankTemperature);

	if (effectiveDensity < 0 || oxTankTemperature < 0) {
		out.Quality = alglib::fp_nan;
		out.Pressure = alglib::fp_nan;
		
		out.FluidDensity = alglib::fp_nan;
		out.FluidSpecificEnthalpy = alglib::fp_nan;
		out.FluidSpecificEntropy = alglib::fp_nan;

		out.LiqDensity = alglib::fp_nan;
		out.LiqEnthalpy = alglib::fp_nan;
		out.LiqEntropy = alglib::fp_nan;
		
		out.VapDensity = alglib::fp_nan;
		out.VapEnthalpy = alglib::fp_nan;
		out.VapEntropy = alglib::fp_nan;
		
		out.State = -1;
	}
	else {
		out.Quality = (out.VapDensity / effectiveDensity) * ((out.LiqDensity - effectiveDensity) / (out.LiqDensity - out.VapDensity));
		out.FluidSpecificEntropy = out.VapEntropy * out.Quality + out.LiqEntropy * (1 - out.Quality);
		out.FluidSpecificEnthalpy= out.VapEnthalpy * out.Quality + out.LiqEnthalpy * (1 - out.Quality);
		out.State = 1;
	}
	return out;
}

void SolomonModel::SetNewState(int i, double time, double mass, double rho, double temp, double tankP,
	double quality, double h, double totH, double mDot, double outP, double st, alglib::real_2d_array &SimState)
{
	SimState[i][0] = time;
	SimState[i][1] = mass;
	SimState[i][2] = rho;
	SimState[i][3] = temp;
	SimState[i][4] = tankP;
	SimState[i][5] = quality;
	SimState[i][6] = h;
	SimState[i][7] = totH;
	SimState[i][8] = mDot;
	SimState[i][9] = outP;
	SimState[i][10] = st;

}

alglib::real_1d_array SolomonModel::pFunc(double temperature, double rho, double P1, double X1)
{
	alglib::real_1d_array pressureQuality;
	double pressure, quality;
	pressureQuality.setlength(2);


	NistDataPoint tableValue = SetDataPoint(temperature, rho);
	pressure = tableValue.Pressure;
	quality = tableValue.Quality;


	pressureQuality[0] = pressure - P1;
	pressureQuality[1] = quality - X1;
	return pressureQuality;
}

/*
alglib::real_1d_array SolomonModel::FindTempRho(double tempGuess, double rhoGuess, alglib::real_1d_array (*func)(double, double, double, double))
{
double temperature, rho;
double xnMinus2 = tempGuess - 10;
double xnMinus1 = tempGuess;
double xn, fnMinus1, fnMinus2;
int maxIterations = 1000;
int i = 0;
double error = 100;
while (error > 0.001) {
fnMinus1 = func(xnMinus1, rhoGuess,
i = i + 1;
// Secant Method
xn = xnMinus1 - fnMinus1*((xnMinus1 - xnMinus2) / (fnMinus1 - fnMinus2));
xnMinus2 = xnMinus1;
xnMinus1 = xn;
error = (abs(xn - xnMinus1) / xn) * 100;
if (i >= maxIterations) { throw(runtime_error("Can't converge")); break; }
}
temperature = xn;

xnMinus2 = rhoGuess + 25;
xnMinus1 = rhoGuess;
i = 0;
error = 100;
while (error > 0.001) {
i = i + 1;
// Secant Method
xn = xnMinus1 - fnMinus1*((xnMinus1 - xnMinus2) / (fnMinus1 - fnMinus2));
xnMinus2 = xnMinus1;
xnMinus1 = xn;
error = (abs(xn - xnMinus1) / xn) * 100;
if (i >= maxIterations) { throw(runtime_error("Can't converge")); break; }
}
}*/

void SolomonModel::TimeStepLoop()
{
	// Volume, initial ox mass, and initial tank temperature should be set using SetInitialProperties
	// This will also set initial density

	double tempGuess, rhoGuess, tankPv;
	NistDataPoint NosTableValues = SetDataPoint(mNosTemperature, mFluidDensity);
	NistDataPoint Upstream, Downstream;
	double Cd = 0.8;
	double Ac = 2.8274333882310001e-05;

	mUpstreamPressure = NosTableValues.Pressure;
	mLiqDensity = NosTableValues.LiqDensity;
	mVapDensity = NosTableValues.VapDensity;
	// paper sets quality differently
	mFluidQuality = NosTableValues.Quality;
	mFluidSpecificEnthalpy = NosTableValues.FluidSpecificEnthalpy;
	mTotalFluidEnthalpy = mFluidSpecificEnthalpy*mFluidMass;
	mState = NosTableValues.State;

	double time = 0;						// [s]
	double mass = mFluidMass;				// [kg]
	double rho = mFluidDensity;				// [kg/m3]
	double temp = mNosTemperature;			// [K]
	double tankP = mUpstreamPressure;		// [psi]
	double quality = mFluidQuality;
	double h = mFluidSpecificEnthalpy;		// [kJ/kg]
	double totH = mTotalFluidEnthalpy;		// [kJ]
	double mDot = 0;						// [kg/s]
	double outP = 14.7;						// psi
	int st = mState;

	alglib::real_2d_array SimState;
	SimState.setlength(1000, 12);
	SetNewState(0, time, mass, rho, temp, tankP, quality, h, totH, mDot, outP, st, SimState);

	for (int i = 1; i < 1000; i++) {
		if (mFluidMass <= 0) { return; }

		/*	guess=[300 300]; %[T,rho]
		pFunc = @(v) [getfield(CO2Props(v(1),v(2)),'P')-P1; getfield(CO2Props(v(1),v(2)),'X') - X1];
		v1 = lsqnonlin(pFunc,guess,0,inf,optimset('Display','off','TolFun',1e-14));
		temp = v1(1); rho = v1(2); */
		/////////////////////////////////////////////////////////////
		alglib::spline1dinterpolant pressureSpline;
		alglib::real_1d_array pressures;
		pressures.setlength(mNistTableLength);
		alglib::spline1dinterpolant qualitySpline;
		alglib::real_1d_array qualities;
		qualities.setlength(mNistTableLength);
		for (int i = 0; i < mNistTableLength; i++)
		{
			pressures[i] = SetDataPoint(NIST.Temperature[i], rho).Pressure;
		}
		alglib::spline1dbuildcubic(pressures, NIST.Temperature, pressureSpline);
		temp = alglib::spline1dcalc(pressureSpline, tankP);
		for (int i = 0; i < mNistTableLength; i++)
		{
			qualities[i] = SetDataPoint(NIST.Temperature[i], rho).Quality;
		}
		/////////////////////////////////////////////////////////////

		Upstream = SetDataPoint(temp, rho);
		tankPv = Upstream.Pressure;
		double rhoL1 = Upstream.LiqDensity;
		h = Upstream.FluidSpecificEnthalpy;
		totH = mass*h;


		/*
		guess=[300 300]; %[T,rho]
		pFunc = @(v) [getfield(CO2Props(v(1),v(2)),'P')-P2, getfield(CO2Props(v(1),v(2)),'h')-h1];
		v2 = lsqnonlin(pFunc,guess,0,inf,optimset('Display','off','TolFun',1e-14));
		T2 = v2(1); rho2 = v2(2);
		*/
		//Downstream=SetDataPoint(T2,rho2);
		outP = Downstream.Pressure;
		double h2 = Downstream.FluidSpecificEnthalpy;


		double k = sqrt((tankP - outP) / (tankPv - outP));
		double W = (1 / (k + 1));
		double mdot_inc = Ac*sqrt(2 * rhoL1*(tankP - outP)*1e6);
		//double mdot_HEM = rho2*Ac*sqrt(2*(h-h2));
		//mDot = Cd*((1-W)*mdot_inc + W*mdot_HEM);

		mDot = injectorModel(temp, outP);
		mass = mass - mDot*mTimeStep;
		totH = totH - h*mTimeStep*mDot;
		rho = mass / mOxTankVolume;
		h = totH / mass;

		/*
		pFunc = @(T_Unknown) getfield(CO2Props(T_Unknown,rho),'h')-h1;
		temp = lsqnonlin(pFunc,300,0,inf,optimset('Display','off','TolFun',1e-14));
		|
		V
		*/
		alglib::spline1dinterpolant enthSpline;
		alglib::real_1d_array enthalpies;
		enthalpies.setlength(mNistTableLength);
		for (int i = 0; i < mNistTableLength; i++)
		{
			enthalpies[i] = SetDataPoint(NIST.Temperature[i], rho).FluidSpecificEnthalpy;
		}
		alglib::spline1dbuildcubic(enthalpies, NIST.Temperature, enthSpline);
		alglib::spline1dcalc(enthSpline, h);

		Downstream = SetDataPoint(temp, rho);
		tankP = Downstream.Pressure;
		quality = Downstream.Quality;
		st = Downstream.State;

		time = time + mTimeStep;
		SetNewState(i, time, mass, rho, temp, tankP, quality, h, totH, mDot, outP, st, SimState);

	}


}
