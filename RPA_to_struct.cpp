#include <iostream>
#include <string>
#include <fstream>
#include "RPA_to_struct.h"

using namespace std;

const int INTEGER_FILE_LENGTH = 110; // Length of csv table
// File path of csv table 
const string FILE_PATH = "C:/Users/hsaafan/Documents/visual studio 2017/Projects/RPA_to_struct/file.csv";

Look_Up_Table Create_Table_Array()
{
	ifstream file(FILE_PATH);			
	string row;
	string delimiter = ",";				// Only works for csv files because they are comma seperated
	Look_Up_Table Table;
	RPA_Table all_structs[INTEGER_FILE_LENGTH];
	int k = 0;
	string cell;
	while (file.good())
	{
		getline(file, row);
		if (k == 0){ k++; continue; }	// skip header
		if (row.empty()) { continue; }	// skip empty rows
		int startS = 0;
		int endS;

		double rowArray[15];

		for (int i = 0; i < 15; i++)
		{
			endS = row.find(delimiter);
			cell = row.substr(startS, endS);
			row = row.substr((endS+1), row.length());
			rowArray[i] = stod(cell);
		}
		struct RPA_Table single;
		// Input array into structure
		single.OF_Ratio = rowArray[0];
		single.Chamber_Pressure = rowArray[1];
		single.Nozzle_inlet = rowArray[2];
		single.Nozzle_exit = rowArray[3];
		single.rho = rowArray[4];
		single.Chamber_Temperture = rowArray[5];
		single.M_value = rowArray[6];
		single.gamma = rowArray[7];
		single.k_value = rowArray[8];
		single.c = rowArray[9];
		single.Is_opt = rowArray[10];
		single.Is_vac = rowArray[11];
		single.Cf_opt = rowArray[12];
		single.Cf_vac = rowArray[13];
		single.c_factor = rowArray[14];

		// Input structures into array and then into another structure
		// Did this because c++ can't return arrays of structs
		// k-1 to skip header row
		Table.RPA_Array[k-1] = single;
		k++;
	}
	return Table;
}

RPA_Table lookUp(double Chamber_Pressure, double OF_Ratio, Look_Up_Table Table)
{
	Chamber_Pressure = round(Chamber_Pressure);
	for (int i = 0; i < sizeof(Table.RPA_Array); i++) 
	{
		if (Table.RPA_Array[i].OF_Ratio == OF_Ratio && Table.RPA_Array[i].Chamber_Pressure == Chamber_Pressure) 
		{
			return Table.RPA_Array[i];
		}
	}
}
