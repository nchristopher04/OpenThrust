#ifndef INJECTORMODEL_H   // To make sure you don't declare the function more than once by including the header multiple times.
#define INJECTORMODEL_H
using namespace std;

double injectorModel(double T1, double P2F);
double data_grab(string property, double value, string TorP, double(&a)[30][100], string(&header)[30]);
void data_gather(ifstream &infile, double(&a2)[30][100], string(&header2)[30]);
#endif