#ifndef RPA_TO_STRUCT_H
#define RPA_TO_STRUCT_H
#include "../libs/alglib/alglibinternal.h"
#include "../libs/alglib/ap.h"
#include "../libs/alglib/interpolation.h"

class RpaTable {
private:
	double mMaxPc;
	double mMinPc;
	double mMaxOF;
	double mMinOF;

	alglib::real_1d_array mOfRatio;
	alglib::real_1d_array mChamberPressure;
	alglib::real_1d_array mNozzleInlet;
	alglib::real_1d_array mNozzleExit;
	alglib::real_1d_array mDensity;
	alglib::real_1d_array mChamberTemperature;
	alglib::real_1d_array mMValue;
	alglib::real_1d_array mGamma;
	alglib::real_1d_array mKValue;
	alglib::real_1d_array mCAsterix;
	alglib::real_1d_array mIsOpt;
	alglib::real_1d_array mIsVac;
	alglib::real_1d_array mCfOpt;
	alglib::real_1d_array mCfVac;
	alglib::real_1d_array mCFactor;
	alglib::real_1d_array mRValue;

	const char* mRpaTablePath;
	int mTableLength;
public:
	struct RpaDataPoint
	{
		double OfRatio;
		double ChamberPressure;	// [psi]
		double NozzleInlet;
		double rho;					// [kg/m^3]
		double NozzleExit;
		double ChamberTemperture;	// [k]
		double MValue;
		double Gamma;
		double KValue;
		double CAsterix;					// [m/s]
		double IsOpt;
		double IsVac;
		double CfOpt;
		double CfVac;
		double CFactor;
		double RValue;
	};

	void SetRpaDataFile(const char* rpaTableFile);
	void CreateRpaTable();
	RpaDataPoint LookUpRpa(double OF, double ChamberPressure);
};

#endif // !RPA_TO_STRUCT_H
