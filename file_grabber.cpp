#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>


using namespace std;

double data_grab(string property,double value, string TorP);

int main(){
cout<< data_grab("Enthalpy (l, kJ/kg)",20,"T")<<endl;

return 0;

}

/* Temperature (C)	Pressure (psia)	Density (l, kg/m3)	Volume (l, m3/kg)	Internal Energy (l, kJ/kg)	Enthalpy (l, kJ/kg)	Entropy (l, J/g*K)	Cv (l, J/g*K)	Cp (l, J/g*K)	Sound Spd. (l, m/s)	Joule-Thomson (l, F/psia)	Viscosity (l, uPa*s)	Therm. Cond. (l, W/m*K)	Surf. Tension (l, N/m)	Density (v, kg/m3)	Volume (v, m3/kg)	Internal Energy (v, kJ/kg)	Enthalpy (v, kJ/kg)	Entropy (v, J/g*K)	Cv (v, J/g*K)	Cp (v, J/g*K)	Sound Spd. (v, m/s)	Joule-Thomson (v, F/psia)	Viscosity (v, uPa*s)	Therm. Cond. (v, W/m*K)
*/


double data_grab(string property,double value, string TorP)
{
string file_grab;
int ref_col;

if(TorP=="P")
{
file_grab = "N20_100_1000PSI.txt";
ref_col =1;
}

else
{
file_grab ="N20_Neg30_35T.txt";
ref_col =0;
}

ifstream infile(file_grab);

double a[30][100];
string header[30];

string line;

getline(infile, line);
string delimiter = "\t";
size_t pos = 0;
int e=0;
while ((pos = line.find(delimiter)) != std::string::npos) {
    header[e] = line.substr(0, pos);
   // std::cout << header[e];
    line.erase(0, pos + delimiter.length());
	e++;
}
std::cout << endl;

int col = 0;
int row = 0;
string transfer;
while (getline(infile, line))
{
while (e>=col && (pos = line.find(delimiter)) != std::string::npos) {
    transfer = line.substr(0, pos);
    istringstream isss(transfer);
    isss >> a[col][row];
   // std::cout << a[col][row]<< " ";
    line.erase(0, pos + delimiter.length());
	col++;
}
col=0;

//cout<<a[col][row];
//cout<<endl;
    row++;
}
int find1 =0;
int find2 =0;
while(property != header[find1])
{
find1++;
if (find1>30)
	return 0;
}

while(a[ref_col][find2]<value)
{
find2++;
}
cout<< "value is at: " << a[ref_col][find2]<<endl;
return a[find1][find2];
}
