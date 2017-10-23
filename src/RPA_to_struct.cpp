#include <iostream>
#include <string>
#include <fstream>
#include "../include/RPA_to_struct.h"
#include <stdexcept>

using namespace std;

const int TABLE_LENGTH = 21000;								// Doesn't actually have to be table length, just longer
const double UNIV_GAS_CONST = 8.314;						// [kJ/kmol*K]

Look_Up_Table Create_Table_Array()
{
	ifstream file;
	string row, cell;
	string delimiter = ",";						// Only works for csv files because they are comma seperated
	Look_Up_Table Table;
	RPA_Table single;
	int startS, endS;
	int k = 0;
	double rowArray[15];	
	string tablePath = "./RPA_Output_Table.csv";

	file.open(tablePath);
	Table.RPA_Vector.reserve(TABLE_LENGTH);		// Reserves spots in look up table
	
	while (file.good())
	{
		getline(file, row);
		if (k == 0){ k++; continue; }			// Skip header
		if (row.empty()) { continue; }			// Skip empty rows
		startS = 0;
		for (int i = 0; i < 15; i++)
		{
			endS = row.find(delimiter);
			cell = row.substr(startS, endS);
			row = row.substr((endS+1), row.length());
			rowArray[i] = stod(cell);
		}
		// Input array into structure
		single.OF_Ratio = rowArray[0];								// Oxidizer/Fuel ratio, rounded
		single.Chamber_Pressure = rowArray[1];						// Pressure of the chamber [psi] rounded
		single.Nozzle_inlet = rowArray[2];							// 
		single.Nozzle_exit = rowArray[3];							//
		single.rho = rowArray[4];									// [kg/m^3] Currently has no values
		single.Chamber_Temperture = rowArray[5];					// Temperature of the chamber [K]
		single.M_value = rowArray[6];								// [kg/kmol]
		single.gamma = rowArray[7];									//
		single.k_value = rowArray[8];								//
		single.c = rowArray[9];										// Characteristic velocity [m/s]
		single.Is_opt = rowArray[10];								// Specific impulse, optimal [s]
		single.Is_vac = rowArray[11];								// Specific impulse, vacuum [s]
		single.Cf_opt = rowArray[12];								// Thrust coefficient, optimal 
		single.Cf_vac = rowArray[13];								// Thrust coefficient, vacuum
		single.c_factor = rowArray[14];								//
		single.R_value = (UNIV_GAS_CONST / single.M_value)*1000;	// [kJ/kg*K]

		// Input structures into vector and then into another structure
		Table.RPA_Vector.push_back(single);
	}
	return Table;
}

Limits_Table Find_Limits(Look_Up_Table Table)
{
	Limits_Table Limits;
	unsigned int maxIndex = Table.RPA_Vector.size();
	Limits.mnPc = Table.RPA_Vector[0].Chamber_Pressure;
	Limits.mnOF = Table.RPA_Vector[0].OF_Ratio;
	Limits.mxPc = Table.RPA_Vector[maxIndex-1].Chamber_Pressure;
	Limits.mxOF = Table.RPA_Vector[maxIndex-1].OF_Ratio;
	return Limits;
}

RPA_Table lookUp(double Chamber_Pressure, double OF, Look_Up_Table Table, Limits_Table Limits)
{
	try 
	{
		if (OF < Limits.mnOF || OF > Limits.mxOF) 
		{ 
			throw invalid_argument("Invalid OF at lookup (0.1-3.0)"); 
		}
		else if (Chamber_Pressure < Limits.mnPc || Chamber_Pressure > Limits.mxPc) 
		{ 
			throw invalid_argument("Chamber Pressure out of range (10-700PSI)"); 
		}
	}catch (invalid_argument& e) {cout << e.what() << " at lookup" << endl;}
	try 
	{
		for (unsigned int i = 0; i < (Table.RPA_Vector.size()); i++)
		{
			if (Table.RPA_Vector[i].OF_Ratio == OF)
			{
				if (Table.RPA_Vector[i].Chamber_Pressure == Chamber_Pressure)
				{
					return Table.RPA_Vector[i];
				}
			}
			if (i == (Table.RPA_Vector.size() - 1)) {throw runtime_error("Lookup Failed");}
		}
	}catch (runtime_error& e){cout << e.what() << '\n'; cin >> OF;}
}
		


