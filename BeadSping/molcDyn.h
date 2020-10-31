#ifndef MOLCDYN_H
#define MOLCDYN_H

#include <vector>

using namespace std;

/***INIT.CPP***/
/***Defining Parameter which won't change throughout the program***/
extern long int particles;
extern double boxLength;
extern const double sigma, epsilon, mass,bondLength ,relaxLen, FEND, rc;
extern double dt;
extern double dt2by2, stiffness;
extern double frc, urc;
extern double potentialEnergy, kineticEnergy, totalEnergy;

/***Defining all the dynamic quantites***/
extern vector<double> x;extern vector<double> y;extern vector<double> z;
extern vector<double> vx;extern vector<double> vy;extern vector<double> vz;
extern vector<double> fx;extern vector<double> fy;extern vector<double> fz;
extern vector<double> newForceX, newForceY, newForceZ;


/****/
extern vector<double> xeq;extern vector<double> yeq;extern vector<double> zeq;
extern vector<double> vxeq;extern vector<double> vyeq;extern vector<double> vzeq;
extern vector<double> fxeq;extern vector<double> fyeq;extern vector<double> fzeq;

/***Defining functions to initialise position,force and velocity***/
/**********************************************************************************
	**posInit2D and velocityInit2D are for 2D initialization and initates all z = 0.
	**Test contains code for testing of 2 particles only. Initialises pos, velocity and force all equal to 0.
	  Position of 1st particle is 1.5 units away.
*/

void posInit(int);
void posInit2D(int);
void posInitLJ();
void posLine();
void velocityInit(double, double, int);
void forceInit();
void velocityInit2D(double, double, int);
void test();
void modifyPos(int);
void posEq();
void velocityEq();
void forceEq();

bool excludeVolume(vector<double>& ,vector<double>& ,vector<double>& ,double ,double, double);
double randInRange(double, double );
void initRand(int);			// For a different random number generator everytime
/*--------------------------------------------------------------------------------------------*/

/***MOLCDYN.CPP***/
/***Defining updating functions***/
void updatePos();
void updateVelocity();
void updateForce();
void velocityVerlet(int, int);
void energyPlot(int, int);
void dataCollection(int, int);
void verlet(int, int);
void energyEquil1(int, int);
void energyEquil2(int, int);

/*--------------------------------------------------------------------------------------------*/

/***CALC.CPP***/
/***Defining functions to calculate stuff***/ 
double calcAvg(vector<double>&);
double calcRg(vector<double>&, vector<double>&,vector<double>&, double, double, double);
double endDistance(double, double, double, double, double, double);
double springForce(double, double, double);
double ljForce(double, double, double);
double minDistance(double, double, double);			//For periodic Boundary Conditions
double ppp(double, double);
double calcKE(vector<double>&, vector<double>&,vector<double>&);
double ljPE(double);
double springPE(double);
/*--------------------------------------------------------------------------------------------*/

/***PUBLISH.CPP***/
void writeData();
void trajectoryFile();
void structureFile(vector<double>&, vector<double>& ,vector<double>&);	//Initial 
void trajectoryFile(vector<double>& ,vector<double>& , vector<double>& ,
					vector<double>& ,vector<double>& , vector<double>&, int );			//To print out any vector quantity at a time
void printVector(vector<double>&);
void printVelocity(vector<double>&);
void equilVector(vector<double>&, vector<double>&);
void printConstants();
void finalData();

#endif