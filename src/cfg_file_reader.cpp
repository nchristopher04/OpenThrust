#include <fstream>
#include <string>
#include <map>
#include "../include/cfg_file_reader.h"

using namespace std;

map<string, double> Rocket_Values;

Rocket_Properties Read_File() 
{
	string Path = "./settings.cfg";
	ifstream file;
	string row;
	string delimiter = ":";
	int delimiterIndex;
	string dictKey, dictValue;
	Rocket_Properties Options;

	file.open(Path);
	while (file.good())
	{
		getline(file, row);
		if (row.empty()) { continue; }
		delimiterIndex = row.find(delimiter);
		dictKey = row.substr(0, delimiterIndex);
		dictValue = row.substr(delimiterIndex + 1, row.length());
		Rocket_Values[dictKey] = stod(dictValue);
	}

	Options.convergeWeighting = Rocket_Values["convergeWeighting"];
	Options.oxTankVolume = Rocket_Values["oxTankVolume"];
	Options.throatArea = Rocket_Values["throatArea"];
	Options.exitArea = Rocket_Values["exitArea"];
	Options.timeStep = Rocket_Values["timeStep"];
	Options.OF = Rocket_Values["OF"];

	// To get the integer values from a double, they are rounded to nearest
	// integer value as floating point number (round()), 0.1 is added to ensure 
	// they are above the integer and then they are typecast into int which cuts 
	// off the floating point part.
	Options.integrationType = int(round(Rocket_Values["integrationType"] + 0.1));
	Options.flowModel = int(round(Rocket_Values["flowModel"])+0.1);
	return Options;
}