#include "../include/SolomonModel.h"
#include "../include/injector_Model.h"
#include "../libs/alglib/ap.h"
#include "../libs/alglib/alglibinternal.h"
#include "../libs/alglib/interpolation.h"
#include "../libs/alglib/solvers.h"
#include "../include/injector_Model.h"
#include "../include/main.h"

#include <fstream>
#include <algorithm>

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
		//throw(runtime_error("Temperature out of bounds"));
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
		out.FluidSpecificEnthalpy = out.VapEnthalpy * out.Quality + out.LiqEnthalpy * (1 - out.Quality);
		out.FluidDensity = out.VapDensity * out.Quality + out.LiqDensity * (1 - out.Quality);
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

SolomonModel::NistDataPoint SolomonModel::MatchPressureQuality(double tempGuess, double rhoGuess, double pressureToMatch, double qualityToMatch)
{
	double P = pressureToMatch;
	double X = qualityToMatch;
	double T = tempGuess;
	double r = rhoGuess;
	int k = 0;
	double learnRate = 0.005;
	double h = 0.001;
	double tol = 0.001;
	double relError = 100;

	double TGrad, rGrad, rPrev, TPrev;

	while (relError > tol)
	{

		rPrev = r;
		TPrev = T;

		k = k + 1;
		if (k > 100000) { throw(runtime_error("PresusreQualityMatching: diverged")); break; }
		//cout << k << endl;
		//cout << relError << endl;
		//if (T > mMaxTemperature) { throw(runtime_error("PresusreQualityMatching: high temperature")); break; }
		//if (T < mMinTemperature) { throw(runtime_error("PresusreQualityMatching: low temperature")); break; }

		TGrad = MatchPressureQualityError(T + h, r, P, X);
		TGrad = TGrad - MatchPressureQualityError(T, r, P, X);
		TGrad = TGrad / h;

		rGrad = MatchPressureQualityError(T, r + h, P, X);
		rGrad = rGrad - MatchPressureQualityError(T, r, P, X);
		rGrad = rGrad / h;

		T = T - learnRate * TGrad;
		r = r - learnRate * rGrad;

		relError = max(abs((r - rPrev) / r), abs((T - TPrev) / T));
	}
	NistDataPoint returnPoint = SetDataPoint(T, r);
	return returnPoint;

}

double SolomonModel::MatchPressureQualityError(double temperature, double rho, double pressure, double quality)
{
	NistDataPoint dataPnt = SetDataPoint(temperature, rho);
	double error;
	error = (pressure - dataPnt.Pressure)*(pressure - dataPnt.Pressure);
	error = error + (quality - dataPnt.Quality)*(quality - dataPnt.Quality);
	return error;
}

SolomonModel::NistDataPoint SolomonModel::MatchPressureEnthalpy(double tempGuess, double rhoGuess, double pressureToMatch, double enthalpyToMatch)
{
	double P = pressureToMatch;
	double E = enthalpyToMatch;
	double T = tempGuess;
	double r = rhoGuess;
	int k = 0;
	double learnRate = 0.001;
	double h = 0.001;
	double tol = 0.01;
	double relError = 100;

	double TGrad, rGrad, rPrev, TPrev;

	while (relError > tol)
	{
		TPrev = T;
		rPrev = r;

		k = k + 1;
		if (k > 100000) { throw(runtime_error("PresusreQualityMatching: diverged")); break; }
		//cout << k << endl;
		//cout << relError << endl;
		//if (T > mMaxTemperature) { throw(runtime_error("PresusreQualityMatching: high temperature")); break; }
		//if (T < mMinTemperature) { throw(runtime_error("PresusreQualityMatching: low temperature")); break; }

		TGrad = MatchPressureEnthalpyError(T + h, r, P, E);
		TGrad = TGrad - MatchPressureEnthalpyError(T, r, P, E);
		TGrad = TGrad / h;

		rGrad = MatchPressureEnthalpyError(T, r + h, P, E);
		rGrad = rGrad - MatchPressureEnthalpyError(T, r, P, E);
		rGrad = rGrad / h;

		T = T - learnRate * TGrad;
		r = r - learnRate * rGrad;

		relError = max(abs((r - rPrev) / r),abs((T-TPrev)/T));
	}
	NistDataPoint returnPoint = SetDataPoint(T, r);
	return returnPoint;

}

double SolomonModel::MatchPressureEnthalpyError(double temperature, double rho, double pressure, double enthalpy)
{
	NistDataPoint dataPnt = SetDataPoint(temperature, rho);
	double error;
	error = (pressure - dataPnt.Pressure)*(pressure - dataPnt.Pressure);
	error = error + (enthalpy - dataPnt.FluidSpecificEnthalpy)*(enthalpy - dataPnt.FluidSpecificEnthalpy);
	return error;
}

void SolomonModel::TimeStepLoop()
{
	// Volume, initial ox mass, and initial tank temperature should be set using SetInitialProperties
	// This will also set initial density
	ofstream file;
	remove("output.csv");
	file.open("output.csv", ios::out | ios::trunc);
	double tempGuess, rhoGuess, tankPv;
	NistDataPoint NosTableValues = SetDataPoint(mNosTemperature, mFluidDensity);
	NistDataPoint Upstream, Downstream;
	double Cd = 0.8;
	double Ac = 0.000028274333882310001;

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
	file << "time" << "," << "mass" << "," << "rho" << "," << "temp" << "," << "tankP" << "," << "quality" << "," << "h" << ","
		<< "totH" << "," << "mDot" << "," << "outP" << "," << "st" << "," << "mdot_HEM" << "," << "mdot_inc" << "," << "thrust" << endl;

	for (int i = 1; i < 1000; i++) {
		if (mFluidMass <= 0) { return; }

		Upstream = MatchPressureQuality(temp, rho, tankP, quality);
		temp = Upstream.Temperature;
		rho = Upstream.FluidDensity;
		tankPv = Upstream.Pressure;
		double rhoL1 = Upstream.LiqDensity;
		h = Upstream.FluidSpecificEnthalpy;
		totH = mass*h;



		Downstream = MatchPressureEnthalpy(300, 300, outP, h);

		//outP = Downstream.Pressure;
		double h2 = Downstream.FluidSpecificEnthalpy;
		double rho2 = Downstream.FluidDensity;

		double kinjector = sqrt((tankP - outP) / (tankPv - outP));
		double W = (1 / (kinjector + 1));
		double mdot_inc = Ac*sqrt(2 * rhoL1*abs(tankP - outP)*6894.76);
		double mdot_HEM = rho2*Ac*sqrt(2*abs(h-h2));
		mDot = Cd*((1-W)*mdot_inc + W*mdot_HEM);

		///////////////////////
		//mDot = mdot_HEM;
		double k, R, Tc, Cf;
		double OF = 2.1;
		double nozzleThroatArea = 0.000382646;
		double nozzleExitArea = 0.00181001;
		interpRPAValues(outP, OF, k, R, Tc);
		double mDotNoz = mDot*((1 + OF) / OF);
		double outP = (mDotNoz*sqrt(k*R*Tc)) / (nozzleThroatArea * k * sqrt(pow((2 / (k + 1)), ((k + 1) / (k - 1)))));
		outP = outP / (6894.76);
		Cf = thrustCoefficient(14.7, nozzleExitArea, outP, nozzleThroatArea);
		double thrust = nozzleThroatArea*(outP*6894.76)*Cf;
		///////////////////////

		mass = mass - mDot*mTimeStep;
		totH = totH - h*mTimeStep*mDot;
		rho = mass / mOxTankVolume;
		h = totH / mass;

		alglib::spline1dinterpolant enthSpline;
		alglib::real_1d_array enthalpies;
		alglib::real_1d_array temperatures;
		temperatures.setlength(mNistTableLength);
		enthalpies.setlength(mNistTableLength);
		for (int i = 0; i < mNistTableLength; i++)
		{
			enthalpies[i] = SetDataPoint(NIST.Temperature[i], rho).FluidSpecificEnthalpy;
			temperatures[i] = NIST.Temperature[i];
		}
		try {
			alglib::spline1dbuildcubic(enthalpies, temperatures, enthSpline);
			temp = alglib::spline1dcalc(enthSpline, h);
		}
		catch (alglib::ap_error& e)
		{ 
			break; 
		}
		Upstream = SetDataPoint(temp, rho);
		tankP = Upstream.Pressure;
		quality = Upstream.Quality;
		st = Upstream.State;

		time = time + mTimeStep;
		SetNewState(i, time, mass, rho, temp, tankP, quality, h, totH, mDot, outP, st, SimState);
		file << time << "," << mass << "," << rho << "," << temp << "," << tankP << "," << quality << "," << h << ","
			<< totH << "," << mDot << "," << outP << "," << st << "," << mdot_HEM << "," << mdot_inc << "," << thrust << endl;	
	}


}
