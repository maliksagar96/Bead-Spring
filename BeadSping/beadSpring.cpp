#include <iostream>
#include "molcDyn.h"
#include <cmath>
#include <time.h>
#include <vector>
#include <fstream>

using namespace std;

int main () {
	//printConstants();
//	cout<<1/pow(2.5,6);

	long int i,j;

	const int snap = 100;
	const int cycles = 10000;//pow(10,5);
	
	const double vmin = -0.50;
	const double vmax = 0.50;

//	posInit(280);	
	posInitLJ();
	forceInit();
	velocityInit(vmin, vmax, 230);
//	verlet(cycles, snap);
	velocityVerlet(cycles, snap);                                              	
//	energyPlot(cycles, snap);
//	dataCollection(cycles, snap);
	
	return 0;
}

