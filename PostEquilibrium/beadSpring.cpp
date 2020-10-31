#include <iostream>
#include "molcDyn.h"
#include <cmath>
#include <time.h>
#include <vector>
#include <fstream>

using namespace std;

int main () {
	
	long int i,j;
	double vmin = -0.5;
	double vmax = 0.5;
	
	long int cycles = 2000;
	long int snap = 1;
	
//	posInit(2880);	
	posLine();
	forceInit();
	velocityInit(vmin, vmax, 212);
	
	posEq();
	velocityEq();
	forceEq();
	
	energyPlot(cycles, snap);
	
	return 0;
}

