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
		bool mSolomonFlag;
		bool mRPACf;


		// Methods
		void SetPath(string path, string delimiter);
		void ReadFile();
		void WriteToFile();
};
#endif // !CFG_FILE_READER
