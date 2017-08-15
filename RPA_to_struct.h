#ifndef RPA_TO_STRUCT_H
#define RPA_TO_STRUCT_H

struct RPA_Table
{
	double OF_Ratio;
	double Chamber_Pressure;	// [psi]
	double Nozzle_inlet;
	double rho;					// [kg/m^3]
	double Nozzle_exit;
	double Chamber_Temperture;	// [k]
	double M_value;
	double gamma;
	double k_value;
	double c;					// [m/s]
	double Is_opt;
	double Is_vac;
	double Cf_opt;
	double Cf_vac;
	double c_factor;
	double R_value;
};

struct Look_Up_Table
{
	RPA_Table RPA_Array[INTEGER_FILE_LENGTH];
};

Look_Up_Table Create_Table_Array();
RPA_Table lookUp(double, double, Look_Up_Table);

#endif // !RPA_TO_STRUCT_H
