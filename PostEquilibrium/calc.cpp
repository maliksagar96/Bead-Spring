#include <iostream>
#include <cmath> 
#include <math.h>
#include "molcDyn.h"

using namespace std;

double calcRg(vector<double>& x_,vector<double>& y_,vector<double>& z_,double xavg, double yavg, double zavg) {
	int i;
	double Rg = 0;
	for(i=0;i<x_.size();i++) {
		Rg = Rg + (pow((x_[i] - xavg),2) + pow((y_[i] - yavg),2)+pow((z_[i] - zavg),2));		
	}	
	Rg = Rg/particles;	
	return sqrt(Rg);
}

double endDistance(double x0, double y0, double z0, double xn, double yn, double zn) {
	return sqrt((x0-xn)*(x0-xn) + (y0-yn)*(y0-yn) + (z0-zn)*(z0-zn));
}

double calcAvg(vector<double>& vec) {
	double avg = 0;
	for(int i = 0; i<vec.size();i++) {
		avg = avg + vec[i];
	}
	return (avg/vec.size());
}


//Always force is calculated on xi partilce due to spring between xi and xj particle
double springForce(double xi, double xj,double r) {
	double force;
	//double xDiff = xi-xj;
	double xDiff = minDistance(xi, xj, boxLength);
	//if(xi == xj) return 0.0;
	force = stiffness*(bondLength - abs(r))*(xDiff)/abs(r);
	
	return force ;
}

double ljForce(double xi, double xj, double r) {
	double sig6 = pow(sigma,6)/pow(r,6);
	//double xDiff  = xi-xj;
	double xDiff = minDistance(xi, xj, boxLength);
	double force = ((6*sig6*epsilon*(2*sig6-1))* (xDiff/(r*r))) - ((frc * xDiff/r));
	return force;
}


double springPE(double r) {
	return 0.5*stiffness*(bondLength - r)*(bondLength - r);
}

double ljPE(double r) {
	
	double PE, sig6;
	sig6 = pow(sigma, 6)/pow(r, 6);
	PE = epsilon*sig6*(sig6-1) - (urc) + (r-rc)*frc;      
	return PE;
}

double calcKE(vector<double>& vx_, vector<double>& vy_, vector<double>& vz_) {
	double KE = 0;
	for(int i=0;i<vx_.size();i++) {
		KE = KE + 0.5*mass*((vx_[i]*vx_[i]) + (vy_[i]*vy_[i]) + (vz_[i]*vz_[i])); 
	}
 	return KE;
}

//Tested with input ri=1.0,rj=14.9 with boxSize = 15
double minDistance(double ri, double rj, double boxSize){	
	
	double dist;	
	
	if((ri - rj) > (boxSize /2.0)) {
		dist = ri - rj - boxSize;
		return dist;
	}
    else if ((ri - rj)< (-boxSize /2.0)) {
		dist = ri - rj + boxSize;			
		return dist;
	} 
    else {
		dist = ri - rj;
		return dist;
	}
}

//Periodic boundary Conditions
double ppp(double xi, double boxSize) {
	double xj = fmod(xi, boxSize);
	if(xj >= 0) return xj;
	else return (boxSize + xj);
}

