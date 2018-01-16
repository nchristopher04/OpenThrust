#include <fstream>
#include <string>
#include <map>
#include "../include/cfg_file_reader.h"

using namespace std;


void OptionFileParser::SetPath(string path, string delimiter)
{
	mPath = path;
	mDelimiter = delimiter;
}

void OptionFileParser::ReadFile()
{
	ifstream file;
	string file_line;
	string dict_key, dict_value;
	int delimiter_index;

	file.open(mPath);
	while (file.good())
	{
		// Iterates the file line by line, skipping certain lines
		// and then parses the file by making anything before the
		// delimiter the key and anything after it the value.
		getline(file, file_line);
		// Empty lines are skipped
		if (file_line.empty()) { continue; }
		// Skips lines that start with # (comments)
		if (file_line.find("#") == 0) { continue; }
		delimiter_index = file_line.find(mDelimiter);
		dict_key = file_line.substr(0, delimiter_index);
		dict_value = file_line.substr(delimiter_index + 1, file_line.length());
		mRocketValues[dict_key] = stod(dict_value);
	}

	mOxTankVolume = mRocketValues["oxTankVolume"];
	mTimeStep = mRocketValues["timeStep"];
	mConvergenceWeight = mRocketValues["convergeWeighting"];
	mThroatArea = mRocketValues["throatArea"];
	mExitArea = mRocketValues["exitArea"];
	mOxFuelRatio = mRocketValues["OF"];
	mSolomonFlag = mRocketValues["solomonMode"];
	mRPACf = mRocketValues["RPACf"];
	RampUpTime = mRocketValues["RampUpTime"];
	RampDownTime = mRocketValues["RampDownTime"];

	// To get the integer values from a double, they are rounded to nearest
	// integer value as floating point number (round()), 0.1 is added to ensure 
	// they are above the integer and then they are typecast into int which cuts 
	// off the floating point part.			

	mFlowModel = int(round(mRocketValues["flowModel"] + 0.1));
	mIntegrationType = int(round(mRocketValues["integrationType"]) + 0.1);
	file.close();
}

void OptionFileParser::WriteToFile()
{
	ofstream file;
	file.open(mPath, ofstream::out|ostream::trunc);

	file << "# OpenThrust Settings File" << endl;
	file << "oxTankVolume" << mDelimiter << mOxTankVolume << endl;
	file << "timeStep" << mDelimiter << mTimeStep << endl;
	file << "convergeWeighting" << mDelimiter << mConvergenceWeight << endl;
	file << "throatArea" << mDelimiter << mThroatArea << endl;
	file << "exitArea" << mDelimiter << mExitArea << endl;
	file << "OF" << mDelimiter << mOxFuelRatio << endl;
	file << "flowModel" << mDelimiter << mFlowModel << endl;
	file << "integrationType" << mDelimiter << mIntegrationType << endl;
	file << "RPACf" << mDelimiter << mRPACf << endl;

	file.close();
}