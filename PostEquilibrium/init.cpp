#include <iostream>
#include <time.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include "molcDyn.h"

using namespace std;

long int  particles = 50;
const double bondLength = 1.0, FEND = 0.0;		

const double sigma = 0.8 * bondLength, epsilon = 1.0, mass = 1.0;		//Effectively Epsilon is 1.
const double rc = sigma * pow(2.0,0.16666666);
double dt = 0.0025;
double boxLength = 40.10 * bondLength;

//FEND = in +ve x-direction for particle at (0,0,0) and in -ve X direction for last particle in chain

double dt2by2 = (dt*dt)/(2.0*mass);
double frc = epsilon * (12.0*(pow(sigma, 12)/pow(rc,13))-(6.0*(pow(sigma,6)/pow(rc,7))));
double sig6 = pow(sigma,6)/pow(rc,6);
//double urc = epsilon*sig6*(sig6-1);
double urc = epsilon*((pow(sigma,12.0)/pow(rc,12.0)) - (pow(sigma,6.0)/pow(rc,6.0)));
//double urc = epsilon*((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/(rc*rc*rc*rc*rc*rc*rc*rc*rc*rc*rc*rc)) - (1/(rc*rc*rc*rc*rc*rc)));
double stiffness = 100.0;

double potentialEnergy = 0.0, kineticEnergy = 0.0, totalEnergy = 0.0;

vector<double> x;vector<double> y;vector<double> z;
vector<double> vx;vector<double> vy;vector<double> vz;
vector<double> fx;vector<double> fy;vector<double> fz;
vector<double> newForceX, newForceY, newForceZ;

vector<double> xeq;vector<double> yeq;vector<double> zeq;
vector<double> vxeq;vector<double> vyeq;vector<double> vzeq;
vector<double> fxeq;vector<double> fyeq;vector<double> fzeq;

void posInit(int seed) {
	int i;
	double x_init = boxLength/2.0, y_init = boxLength/2.0, z_init = boxLength/2.0; 
	double phi = 0.0,theta = 0.0, l = bondLength;
	
	x.insert(x.begin(), x_init);
	y.insert(y.begin(), y_init);
	z.insert(z.begin(), z_init);
	
	vector<double> rand1;
	vector<double> rand2;
	
	initRand(seed);
	
	for(i = 1;i<particles;i++) {
		
		theta = 3.14159*randInRange(0,1);		
		phi = 2*3.14159*randInRange(0,1);
		
		x.insert(x.begin() + i,x[i-1] + l*sin(theta)*cos(phi));
		y.insert(y.begin() + i,y[i-1] + l*sin(theta)*sin(phi));
		z.insert(z.begin() + i,z[i-1] + l*cos(theta));
		
		
		if(excludeVolume(x,y,z,x[i],y[i],z[i])) {
			i = i-1;
			x.erase(x.begin() + x.size() - 1);
			y.erase(y.begin() + y.size() - 1);
		}		
	}
	structureFile(x,y,z);	
}

void posEq() {

	ifstream inputFile("pos.dat");

	double x1,y1,z1;

	while (inputFile >> x1 >> y1 >> z1)
	{
		xeq.push_back(x1);
		yeq.push_back(y1);
		zeq.push_back(z1);
	}
	
	equilVector(xeq,x);
	equilVector(yeq,y);
	equilVector(zeq,z);
	/*
	for(int i = 0; i < xeq.size(); i++)
	{
		cout << x[i] << "\t\t" << y[i]<< "\t\t" << z[i] << endl;
	}	*/
}

void velocityEq() {
	ifstream inputFile("vel.dat");

	double x1,y1,z1;

	while (inputFile >> x1 >> y1 >> z1)
	{
		vxeq.push_back(x1);
		vyeq.push_back(y1);
		vzeq.push_back(z1);
	}
	
	equilVector(vxeq,vx);
	equilVector(vyeq,vy);
	equilVector(vzeq,vz);
	/*
	for(int i = 0; i < vxeq.size(); i++)
	{
		cout << vx[i] << "\t\t" << vy[i]<< "\t\t" << vz[i] << endl;
	}*/	
}

void forceEq() {
	ifstream inputFile("forces.dat");

	double x1,y1,z1;

	while (inputFile >> x1 >> y1 >> z1)
	{
		fxeq.push_back(x1);
		fyeq.push_back(y1);
		fzeq.push_back(z1);
	}
	
	equilVector(fxeq,fx);
	equilVector(fyeq,fy);
	equilVector(fzeq,fz);
	
	for(int i = 0; i < fxeq.size(); i++)
	{
		cout << fx[i] << "\t\t" << fy[i]<< "\t\t" << fz[i] << endl;
	}
}

void posInitLJ() {
	int i, maxCount, countx = 1, county = 1, countz = 1;
	double x_init = 1.0, y_init = 1, z_init = 1.0; 
	double sepFactor = sigma * 1.2;
	
	maxCount = (int)(boxLength/sepFactor);
	
	x.insert(x.begin(), x_init);
	y.insert(y.begin(), y_init);
	z.insert(z.begin(), z_init);
	
	for(i = 1;i<particles;i++) {
		x.insert(x.begin()+i, 0);
		y.insert(y.begin()+i, 0);
		z.insert(z.begin()+i, z_init);
	}
	
	for(i = 1;i<particles;i++) {		
		
		if(countx == maxCount) {
			x[i] = x_init;	
			countx = 1;
		}
		else {
			x[i] = x[i-1] + sepFactor;
			countx++;
		}
	}
	
	countx=1;
	
	for(i=1;i<particles;i++) {
		
		
		
		if(countx == maxCount) {
			y[i] = y[i-1] + sepFactor;
			countx = 1;
			county++;
		}
	
	
		else {
			y[i] = y[i-1];
			countx++;
		}
		

		if(county > maxCount){
			y[i] = y_init;
			countx = 1;
			county = 1;
		}		
		
	}
	
	for(i=1;i<particles;i++) {
		
		if(countz == (maxCount*maxCount)) {
			z[i] = z[i-1] + sepFactor;
			countz = 1;
		}
		
		else{
			z[i] = z[i-1];
			countz++;
		}
	}

	structureFile(x,y,z);	
	
}

void velocityInit(double vmin, double vmax, int seed) {
	srand(seed);
	int i;
	double avgVx = 0.0, avgVy = 0.0, avgVz = 0.0;
	vx.insert(vx.begin(), 0.0);
	vy.insert(vy.begin(), 0.0);
	vz.insert(vz.begin(), 0.0);

//Initializing velocity compoenents from vmin to vmax.	
	for(i = 1;i<particles;i++) {
		vx.insert(vx.begin() + i, randInRange(vmin, vmax));
		vy.insert(vy.begin() + i, randInRange(vmin, vmax));
		vz.insert(vz.begin() + i, randInRange(vmin, vmax));	
	}
	
	for(i=0;i<particles;i++) {
			avgVx = avgVx + vx[i];
			avgVy = avgVy + vy[i];
			avgVz = avgVz + vz[i];
	}
	
	avgVx = avgVx/particles;
	avgVy = avgVy/particles;
	avgVz = avgVz/particles;
	
	for(i = 0;i<particles;i++) {
		vx[i] = vx[i] - avgVx;
		vy[i] = vy[i] - avgVy;
		vz[i] = vz[i] - avgVz;
	}
}

void forceInit() {
	int i;
	
//Initializing Forces to 0.
	for(i = 0;i<particles;i++) {
		fx.insert(fx.begin() + i, 0.0);
		fy.insert(fy.begin() + i, 0.0);
		fz.insert(fz.begin() + i, 0.0);	
	}
	
	for(i=0;i<particles;i++) {
		newForceX.insert(newForceX.begin() + i, 0.0);
		newForceY.insert(newForceY.begin() + i, 0.0);
		newForceZ.insert(newForceZ.begin() + i, 0.0);
	}		
	
}


void posInit2D(int seed) {
	int i;
	double x_init = boxLength/2, y_init = boxLength/2, z_init = boxLength/2; 
	double phi = 0.0, l = bondLength;
	//Took kbT = 1, taking stiffness for SAW
	double stiffness = 100/(bondLength*bondLength * pow(particles, 1.2));		
	
	//Initializing x, y and z coordinates together
	x.insert(x.begin(), x_init);
	y.insert(y.begin(), y_init);
	z.insert(z.begin(), z_init);
	
	vector<double> rand1;
	vector<double> rand2;
	
	initRand(seed);
	
	for(i = 1;i<particles;i++) {
			
		phi = 2*3.14159*randInRange(0,1);
		
		x.insert(x.begin() + i,x[i-1] + l*cos(phi));
		y.insert(y.begin() + i,y[i-1] + l*sin(phi));
		z.insert(z.begin() + i,z[i-1] + 0);
		
		
		if(excludeVolume(x,y,z,x[i],y[i],z[i])) {
			i = i-1;
			x.erase(x.begin() + x.size() - 1);
			y.erase(y.begin() + y.size() - 1);
		}
			
	}
	
	//structureFile(x,y,z);
}

void velocityInit2D(double vmin, double vmax, int seed) {
	srand(seed);
	int i;
	double avgVx = 0, avgVy = 0, avgVz = 0;
	vx.insert(vx.begin(), 0.0);
	vy.insert(vy.begin(), 0.0);
	vz.insert(vz.begin(), 0.0);

//Initializing velocity compoenents from vmin to vmax.	
	for(i = 1;i<particles;i++) {
		vx.insert(vx.begin() + i, randInRange(vmin, vmax));
		vy.insert(vy.begin() + i, randInRange(vmin, vmax));
		vz.insert(vz.begin() + i, 0);	
	}
	
	for(i=0;i<particles;i++) {
			avgVx = avgVx + vx[i];
			avgVy = avgVy + vy[i];
	}
	
	avgVx = avgVx/particles;
	avgVy = avgVy/particles;
	
	for(i = 0;i<particles;i++) {
		vx[i] = vx[i] - avgVx;
		vy[i] = vy[i] - avgVy;
	}
}

void test() {
	particles = 2;

	int i;
	double x_init = boxLength/2, y_init = boxLength/2, z_init = boxLength/2; 
	
	double stiffness = 1;
	
	//Initializing x, y and z coordinates together
	x.insert(x.begin(), x_init);
	y.insert(y.begin(), y_init);
	z.insert(z.begin(), z_init);
	
	for(i = 1;i<particles;i++) {
		x.insert(x.begin() + i,x_init + 1.1);
		y.insert(y.begin() + i,y_init);
		z.insert(z.begin() + i,z_init);
	}
	
	vx.insert(vx.begin(), 0.0);
	vy.insert(vy.begin(), 0.0);
	vz.insert(vz.begin(), 0.0);
	
	
	for(i = 1;i<particles;i++) {
		vx.insert(vx.begin() + i,0.0);
		vy.insert(vy.begin() + i,0.0);
		vz.insert(vz.begin() + i,0.0);
	}
	
	//Initializing Forces to 0.
	for(i = 0;i<particles;i++) {
		fx.insert(fx.begin() + i, 0.0);
		fy.insert(fy.begin() + i, 0.0);
		fz.insert(fz.begin() + i, 0.0);	
	}
	
	for(i=0;i<particles;i++) {
		newForceX.insert(newForceX.begin() + i, 0.0);
		newForceY.insert(newForceY.begin() + i, 0.0);
		newForceZ.insert(newForceZ.begin() + i, 0.0);
	}

	structureFile(x,y,z);

}

void posLine() {
	int i;
	boxLength = 2*particles*bondLength;
	double x_init = boxLength/4, y_init = boxLength/2, z_init = boxLength/2; 

	double l = bondLength;
	
	x.insert(x.begin(), x_init);
	y.insert(y.begin(), y_init);
	z.insert(z.begin(), z_init);
	
	for(i = 1;i<particles;i++) {
		
		x.insert(x.begin() + i,x[i-1] + l);
		y.insert(y.begin() + i,y[i-1]);
		z.insert(z.begin() + i,z[i-1]);
		
	}
	
	structureFile(x,y,z);
}

bool excludeVolume(vector<double>& x_,vector<double>& y_,vector<double>& z_,double xi,double yi, double zi) {
	bool vFlag = false;
	double dist;
	for(int i = 0;i<x_.size()-1;i++) {
		dist = ((x_[i] - xi)*(x_[i] - xi)) + ((y_[i] - yi)*(y_[i] - yi) +(z_[i] - zi)*(z_[i] - zi));
		dist = sqrt(dist);
		if(dist <= sigma){
			vFlag = true;
			return vFlag;
		}
	}
	return vFlag;
}	

void initRand(int seed) {
	srand(seed); 
}

double randInRange(double min, double max) {
	return min + (max - min)*(rand()/(double)(RAND_MAX));;
}

//For ensmeble average
void modifyPos(int seed) {
	
	int i;
	double phi = 0.0,theta = 0.0, l = bondLength;
	
	srand(seed);
	
	for(i=0;i<x.size();i++){
		x[i] = 0;
		y[i] = 0;
		z[i] = 0;
	}
	
	srand(seed);
	
	for(i = 1;i<particles;i++) {
		
		theta = 3.14159*randInRange(0.0,1.0);		
		phi = 2.0*3.14159*randInRange(0.0,1.0);
		
		x[i] = x[i-1] + l*sin(theta)*cos(phi);
		y[i] = y[i-1] + l*sin(theta)*sin(phi);
		z[i] = z[i-1] + l*cos(theta);
	}
}




