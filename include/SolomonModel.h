#ifndef SOLOMON_MODEL_H
#define SOLOMON_MODEL_H

#include <fstream>
#include <string>
#include "../include/injector_Model.h"
#include "../libs/alglib/ap.h"
#include "../libs/alglib/alglibinternal.h"
#include "../libs/alglib/alglibmisc.h"
using namespace std;

class SolomonModel
{
private:
	// Nomenclature name refers to the name used for these variables in
	// the paper from Solomon et al. under the nomenclature section
	
	// Inputted Values								Nomenclature Name
	double mOxTankVolume;						//	Tank volume
	double mInitialMassOxidizer;				//	Initial fluid mass
	double mInitialTemperature;					//	Initial fluid temperature
	double mInitialFluidDensity;				//	Initial fluid density
	double mTimeStep;							//

	// NOS Properties Data Files and Storage
	const char* mTemperatureIncrementedFile;
	
	alglib::real_2d_array NistFile;

	int mNistTableLength;

	double mMaxTemperature;
	double mMinTemperature;

	// NOS Variables								Nomenclature Name
	int mState;	
	
	double mUpstreamPressure;					//	Pressure
	double mDownstreamPressure;
	double mNosTemperature;						//	

	double mLiqMass;							//
	double mVapMass;							//
	double mFluidMass;							// 
	
	double mLiqDensity;							//	Liquid density
	double mVapDensity;							//	Vapor density
	double mFluidDensity;						// 
	double mDownstreamDensity;
	
	double mLiqSpecificEntropy;					//	Liquid specifc entropy
	double mVapSpecificEntropy;					//	Vapor specifc entropy
	double mFluidSpecificEntropy;				//	Initial fluid specific entropy
	double mEntropyFlowRate;
	double mTotalFluidEntropy;

	double mLiqSpecificEnthalpy;				//	Liquid specific enthalpy
	double mVapSpecificEnthalpy;				//	Vapor specific enthalpy
	double mFluidSpecificEnthalpy;				//	Initial fluid specific enthalpy
	double mDownstreamSpecificEnthalpy;
	double mEnthalpyFlowRate;
	double mTotalFluidEnthalpy;					//	Initial fluid total enthalpy

	double mFluidQuality;						//	Initial fluid quality

	// Orific Variables								Nomenclature Name
	double mOrificeCrossSectionalArea;			// Cross sectional area
	double mOrificeCount;						// 

	double mOrificeDischargeCoefficient;		// Discharge coefficient

	double mNonEquilibriumParameter;

	// Mass Flow Rates
	double mMassFlowRateIncompressible;
	double mMassFlowRateHomogenous;
	double mMassFlowRateDryer;
	double mMassFlowRate;

public:
	
	struct NistData {
		alglib::real_1d_array Temperature;
		alglib::real_1d_array Pressure;
		alglib::real_1d_array LiqDensity;
		alglib::real_1d_array VapDensity;
		alglib::real_1d_array LiqEnthalpy;
		alglib::real_1d_array VapEnthalpy;
		alglib::real_1d_array LiqEntropy;
		alglib::real_1d_array VapEntropy;
	}NIST;

	struct NistDataPoint {
		double Temperature;
		double Pressure;
		double Quality;

		double FluidSpecificEntropy;
		double LiqEntropy;
		double VapEntropy;

		double FluidSpecificEnthalpy;
		double LiqEnthalpy;
		double VapEnthalpy;

		double FluidDensity;
		double LiqDensity;
		double VapDensity;
		
		double State;
	};

	void SetInitialProperties(double initialOxTankVolume, double initialMassOxidizer, double initialTemperature, double timeStep);
	void SetNosDataFile(const char* temperatureIncrementedFile);
	void ImportNosData();
	NistDataPoint SetDataPoint(double oxTankTemperature, double effectiveDensity);

	void SetNewState(int i, double time, double mass, double rho, double temp, double tankP,
		double quality, double h, double totH, double mDot, double outP, double st, alglib::real_2d_array &SimState);
	alglib::real_1d_array pFunc(double mass, double rho, double P1, double X1);
	alglib::real_1d_array SolomonModel::FindTempRho(double, double, alglib::real_1d_array(*func)(double, double, double, double));
	void TimeStepLoop();


};
#endif // !SOLOMON_MODEL_H
