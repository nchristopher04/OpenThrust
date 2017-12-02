#include <iostream>
#include <string>
#include <stdexcept>
#include "../include/RpaTable.h"
#include "../libs/alglib/alglibinternal.h"
#include "../libs/alglib/ap.h"

using namespace std;

const double UNIV_GAS_CONST = 8.314;						// [kJ/kmol*K]

void RpaTable::SetRpaDataFile(const char* rpaTableFile) 
{
	mRpaTablePath = rpaTableFile;
}

void RpaTable::CreateRpaTable()
{
	alglib::real_2d_array RpaTableArray;
	alglib::read_csv(mRpaTablePath, ',', 1, RpaTableArray);
	mTableLength = RpaTableArray.rows();


	mOfRatio.setlength(mTableLength);
	mChamberPressure.setlength(mTableLength);
	mNozzleInlet.setlength(mTableLength);
	mNozzleExit.setlength(mTableLength);
	mDensity.setlength(mTableLength);
	mChamberTemperature.setlength(mTableLength);
	mMValue.setlength(mTableLength);
	mGamma.setlength(mTableLength);
	mKValue.setlength(mTableLength);
	mCAsterix.setlength(mTableLength);
	mIsOpt.setlength(mTableLength);
	mIsVac.setlength(mTableLength);
	mCfOpt.setlength(mTableLength);
	mCfVac.setlength(mTableLength);
	mCFactor.setlength(mTableLength);
	mRValue.setlength(mTableLength);

	for (int i = 0; i < mTableLength; i++) {
		mOfRatio[i] = RpaTableArray[i][0];
	}
	for (int i = 0; i < mTableLength; i++) {
		mChamberPressure[i] = RpaTableArray[i][1];
	}
	for (int i = 0; i < mTableLength; i++) {
		mNozzleInlet[i] = RpaTableArray[i][2];
	}
	for (int i = 0; i < mTableLength; i++) {
		mNozzleExit[i] = RpaTableArray[i][3];
	}
	for (int i = 0; i < mTableLength; i++) {
		mDensity[i] = RpaTableArray[i][4];
	}
	for (int i = 0; i < mTableLength; i++) {
		mChamberTemperature[i] = RpaTableArray[i][5];
	}
	for (int i = 0; i < mTableLength; i++) {
		mMValue[i] = RpaTableArray[i][6];
	}
	for (int i = 0; i < mTableLength; i++) {
		mGamma[i] = RpaTableArray[i][7];
	}
	for (int i = 0; i < mTableLength; i++) {
		mKValue[i] = RpaTableArray[i][8];
	}
	for (int i = 0; i < mTableLength; i++) {
		mCAsterix[i] = RpaTableArray[i][9];
	}
	for (int i = 0; i < mTableLength; i++) {
		mIsOpt[i] = RpaTableArray[i][10];
	}
	for (int i = 0; i < mTableLength; i++) {
		mIsVac[i] = RpaTableArray[i][11];
	}
	for (int i = 0; i < mTableLength; i++) {
		mCfOpt[i] = RpaTableArray[i][12];
	}
	for (int i = 0; i < mTableLength; i++) {
		mCfVac[i] = RpaTableArray[i][13];
	}
	for (int i = 0; i < mTableLength; i++) {
		mCFactor[i] = RpaTableArray[i][14];
	}
	for (int i = 0; i < mTableLength; i++) {
		mRValue[i] = (UNIV_GAS_CONST / mMValue[i])*1000;
	}

	mMinOF = mOfRatio(0);
	mMaxOF = mOfRatio(mTableLength - 1);

	mMinPc = mChamberPressure(0);
	mMaxPc = mChamberPressure(mTableLength - 1);

}

RpaTable::RpaDataPoint RpaTable::LookUpRpa(double OF, double ChamberPressure)
{
	RpaDataPoint out;
	if (OF < mMinOF || OF > mMaxOF || ChamberPressure < mMinPc || ChamberPressure > mMaxPc)
	{
		throw(runtime_error("OF or Chamber Pressure out of bounds"));
	}
	
	for (unsigned int i = 0; i < (mOfRatio.length()); i++)
	{
		if (mOfRatio[i] == OF)
		{
			if (mChamberPressure[i] == ChamberPressure)
			{
				out.OfRatio = mOfRatio[i];
				out.ChamberPressure = mChamberPressure[i];
				out.NozzleInlet = mNozzleInlet[i];
				out.NozzleExit = mNozzleExit[i];
				out.rho = mDensity[i];
				out.ChamberTemperture = mChamberTemperature[i];
				out.MValue = mMValue[i];
				out.Gamma = mGamma[i];
				out.KValue = mKValue[i];
				out.CAsterix = mCAsterix[i];
				out.IsOpt = mIsOpt[i];
				out.IsVac = mIsVac[i];
				out.CfOpt = mCfOpt[i];
				out.CfVac = mCfVac[i];
				out.CFactor = mCFactor[i];
				out.RValue = mRValue[i];
				return out;
			}
		}
	}
}