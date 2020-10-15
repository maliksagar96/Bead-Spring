#include <iostream>
#include "molcDyn.h"
#include <cmath>
#include <time.h>
#include <vector>
#include <fstream>

using namespace std;

int main () {
	
	long int i,j;

	const int snap = 250;
	const int cycles = 2*pow(10,4);
	const double vmin = 0.00;
	const double vmax = 0.00;

	test();
/*
	posLine();	
	forceInit();
	velocityInit(vmin, vmax, 230);
*/
	velocityVerlet(cycles, snap);	

//	dataCollection(cycles, snap);
	return 0;
}

