#ifndef CFG_FILE_READER
#define CFG_FILE_READER

struct Rocket_Properties {
	double oxTankVolume;
	double timeStep;
	int flowModel;
	int integrationType;
	double convergeWeighting;
	double throatArea;
	double exitArea;
	double OF;
};

Rocket_Properties Read_File();


#endif // !CFG_FILE_READER
