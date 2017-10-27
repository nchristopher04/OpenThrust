#ifndef CFG_FILE_READER
#define CFG_FILE_READER

using namespace std;

#include <map>

class OptionFileParser
{
	private:

		map <string, double> mRocketValues;

		// Options file properties
		string mPath;
		string mDelimiter;

	public:

		// Variables in the options file
		double mOxTankVolume;
		double mTimeStep;
		double mConvergenceWeight;
		double mThroatArea;
		double mExitArea;
		double mOxFuelRatio;
		int mFlowModel;
		int mIntegrationType;

		void SetPath(string path, string delimiter);
		void ReadFile();
};
#endif // !CFG_FILE_READER
