#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include "molcDyn.h"
#include <iomanip>

using namespace std;


void updatePos(){
	int i;
	for(i=0;i<particles;i++) {
		x[i] = x[i] + (vx[i] * dt) + (fx[i] * dt2by2);
		y[i] = y[i] + (vy[i] * dt) + (fy[i] * dt2by2);
		z[i] = z[i] + (vz[i] * dt) + (fz[i] * dt2by2);
		
		//Periodic Boundary Condition
		x[i] = ppp(x[i], boxLength);
		y[i] = ppp(y[i], boxLength);
		z[i] = ppp(z[i], boxLength);
	
	}
}

void updateForce() {
	int i,j;
	
	double r;
	double xF = 0, yF = 0, zF = 0;
	
	potentialEnergy = 0;
	kineticEnergy = 0;
	totalEnergy = 0;

	for(i=0;i<particles;i++) {
		newForceX[i] = 0.0;
		newForceY[i] = 0.0;
		newForceZ[i] = 0.0;
	}		

	
	//Force Calculation  for 0th particle
	r = sqrt((x[0]-x[1])*(x[0]-x[1])+
			 (y[0]-y[1])*(y[0]-y[1])+
			 (z[0]-z[1])*(z[0]-z[1]));
	
	newForceX[0] = springForce(x[0], x[1], r) + FEND;
	newForceY[0] = springForce(y[0], y[1], r);
	newForceZ[0] = springForce(z[0], z[1], r);
	//Calculating PE for 1st and 2nd Partilce
	potentialEnergy = potentialEnergy + springPE(r);
	
	//Force Calculation for last particle
	r = sqrt((x[particles-1]-x[particles-2])*(x[particles-1]-x[particles-2])+
			 (y[particles-1]-y[particles-2])*(y[particles-1]-y[particles-2])+
			 (z[particles-1]-z[particles-2])*(z[particles-1]-z[particles-2]));
	
	newForceX[particles-1] = -FEND + springForce(x[particles-1], x[particles-2], r);
	newForceY[particles-1] = springForce(y[particles-1], y[particles-2], r);
	newForceZ[particles-1] = springForce(z[particles-1], z[particles-2], r);
	
	for(i = 1;i<(particles - 1);i++) {
		
		r = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+
				 (y[i]-y[i-1])*(y[i]-y[i-1])+
				 (z[i]-z[i-1])*(z[i]-z[i-1]));
		//Force due to previous particle
		newForceX[i] = newForceX[i] + springForce(x[i], x[i-1], r);
		newForceY[i] = newForceY[i] + springForce(y[i], y[i-1], r);
		newForceZ[i] = newForceZ[i] + springForce(z[i], z[i-1], r);
			
		//Force due to next particle
		r =sqrt((x[i]-x[i+1])*(x[i]-x[i+1])+
				(y[i]-y[i+1])*(y[i]-y[i+1])+
				(z[i]-z[i+1])*(z[i]-z[i+1]));
		
		newForceX[i] = newForceX[i] + springForce(x[i], x[i+1], r);
		newForceY[i] = newForceY[i] + springForce(y[i], y[i+1], r);
		newForceZ[i] = newForceZ[i] + springForce(z[i], z[i+1], r);		
		//Calculating PE of all particles except 1st and 2nd
		potentialEnergy = potentialEnergy + springPE(r);
	}
	
	
	
	for(i = 0;i<(particles-1);i++) {
		for(j=(i+1);j<particles;j++) {
			
		r = sqrt(minDistance(x[i], x[j], boxLength)*minDistance(x[i], x[j], boxLength) + minDistance(y[i], y[j], boxLength)*minDistance(y[i], y[j], boxLength) + minDistance(z[i], z[j], boxLength)*minDistance(z[i], z[j], boxLength));
				
	//		r = sqrt((x[i] - x[j])*(x[i] - x[j])+ (y[i] - y[j])*(y[i] - y[j])+(z[i] - z[j])*(z[i] - z[j])); 				 
		
		
			if(r <= rc) {
				
				xF = ljForce(x[i],x[j],r);
				yF = ljForce(y[i],y[j],r);
				zF = ljForce(z[i],z[j],r);
				
				newForceX[i] = newForceX[i] + xF;
				newForceY[i] = newForceY[i] + yF;
				newForceZ[i] = newForceZ[i] + zF;
				
				newForceX[j] = newForceX[j] - xF;
				newForceY[j] = newForceY[j] - yF;
				newForceZ[j] = newForceZ[j] - zF;
				
				potentialEnergy = potentialEnergy + ljPE(r);
			}
		}
	}
	
	
}

void updateVelocity(){
	int i;
	double dtby2 = dt/(mass * 2.0);	
	for(i=0;i<particles;i++) {  
        vx[i] = vx[i] + (fx[i] + newForceX[i]) * dtby2;
        vy[i] = vy[i] + (fy[i] + newForceY[i]) * dtby2;
        vz[i] = vz[i] + (fz[i] + newForceZ[i]) * dtby2;
	}
	
	for(i=0;i<particles;i++) {
		fx[i] = newForceX[i];
		fy[i] = newForceY[i];
		fz[i] = newForceZ[i];
	}

	kineticEnergy = calcKE(vx,vy,vz);
	totalEnergy = kineticEnergy + potentialEnergy;
}

void velocityVerlet(int cycles, int snap) {
	
	long int i,j, rgcount = 0;
	double rg = 0, endDist = 0;
	fstream bs;
	bs.open("bs.trj", ios::out);

	for(j = 0;j<cycles;j++) {
	
		if(((j%snap) == 0)) {
			bs<<"ITEM: TIMESTEP \n"<<j<<endl;
			bs<<"ITEM: NUMBER OF ATOMS\n"<<particles<<endl;
			bs<<"ITEM: BOX BOUNDS pp pp pp"<<endl;	
			bs<<0.0<<"\t"<<boxLength<<endl;
			bs<<0.0<<"\t"<<boxLength<<endl;
			bs<<0.0<<"\t"<<boxLength<<endl;
			bs<<"ITEM: ATOMS id x y z vx vy vz"<<endl;
	
			for(i = 0;i<x.size();i++) {
				bs<<i<<" "<<x[i]<<" "<<y[i]<<" "<<z[i]<<" "<<vx[i]<<" "<<vy[i]<<" "<<vz[i]<<" "<<endl; 
			}		
		cout<<"Completed Simulation = "<<j<<endl;
		
		//printVelocity(vz);
	
		endDist = endDist + endDistance(x[0], y[0], z[0], x[particles -1], y[particles-1], z[particles-1]);	
		rg = rg + calcRg(x,y,z,calcAvg(x), calcAvg(y), calcAvg(z));
		rgcount++;

		}
	
		updatePos();
		updateForce();
		updateVelocity();
	}
	
	rg = rg/rgcount;
	endDist = endDist/rgcount;
	cout<<"Particles = "<<particles<<endl;
	cout<<"Avg end to end Distance = "<<endDist<<endl;	
	cout<<"Avg Rg = "<<rg<<endl;
	
	cout<<"Velocity Verlet Completed"<<endl;
	bs.close();
}

void energyPlot(int cycles, int snap) {
	
	long int i,j;
	fstream energy;
	energy.open("Energy.dat", ios::out);
	
	for(j = 0;j<cycles;j++) {
		
		updatePos();
		updateForce();
		updateVelocity();	
		
		if(((j%snap) == 0)) {
			energy<<setw(6)<<j<<setprecision(10)<<setw(22)<<potentialEnergy<<setw(22)<<kineticEnergy<<setw(22)<<totalEnergy<<endl;	
			cout<<"Completed Simulation = "<<j+snap<<endl;
		}		
	}
	cout<<"Energy Plotting completed"<<endl;
	energy.close();
}

void energyEquil1(int cycles, int snap) {
	
	long int i,j;
	fstream energy1;
	energy1.open("Energyequilibrium.dat", ios::out);
	
	for(j = 0;j<cycles;j++) {
		
		updatePos();
		updateForce();
		updateVelocity();	
		
		if(((j%snap) == 0)) {
			energy1<<setw(6)<<j<<setprecision(10)<<setw(22)<<potentialEnergy<<setw(22)<<kineticEnergy<<setw(22)<<totalEnergy<<endl;
		}
	}
	cout<<"Energy Equilibrium 1 completed"<<endl;
	energy1.close();
}

void energyEquil2(int cycles, int snap) {
	
	long int i,j;
	fstream energy2;
	energy2.open("Energyequilibrium2.dat", ios::out);
	
	for(j = 0;j<cycles;j++) {
		
		updatePos();
		updateForce();
		updateVelocity();	
		
		if(((j%snap) == 0)) {
			energy2<<setw(6)<<j<<setprecision(10)<<setw(22)<<potentialEnergy<<setw(22)<<kineticEnergy<<setw(22)<<totalEnergy<<endl;		
		}
	}
	cout<<"Energy Equilibrium 2completed"<<endl;
	energy2.close();
}

void verlet(int cycles, int snap) {
	
	for(int i = 0;i<cycles;i++) {
		updatePos();
		updateForce();
		updateVelocity();
		if((i%snap) == 0) {
			cout<<"Completed Simulation = "<<i<<endl;
		}
	}
}

void dataCollection(int cycles, int snap) {
	
	long int i,j,k, rgcount = 0;
	double rg = 0,endDist = 0;
	fstream data;
	data.open("rg.dat", ios::out);
	
	for(k = 0;k<10;k++) {
		endDist = 0;
		rg = 0;
		rgcount = 0;
		posLine();
		velocityInit(-0.50,0.50,212);
		forceInit();
		
		cout<<"Equilibrating For K = "<<k<<" Number of particles = "<<particles<<endl;
		
		//DataCollection loop
		for(j = 0;j<cycles;j++) {
		
			updatePos();
			updateForce();
			updateVelocity();
		
			if(((j%snap) == 0)) {
				endDist = endDist + endDistance(x[0], y[0], z[0], x[particles -1], y[particles-1], z[particles-1]);
				rg = rg + calcRg(x,y,z,calcAvg(x), calcAvg(y), calcAvg(z));
				rgcount++;
			}
		}
		rg = rg/rgcount;
		endDist = endDist/rgcount;
		
		cout<<"Equilibration Achieved for K = "<<k<<endl;
		cout<<"Avg Rg = "<<rg<<endl;
		cout<<"End Distance = "<<endDist<<endl;
		data<<particles<<"\t \t"<<rg<<"\t \t"<<endDist<<endl;

		particles = particles + 10;
		
		x.clear();
		y.clear();
		z.clear();
		
		vx.clear();
		vy.clear();
		vz.clear();
		
		fx.clear();
		fy.clear();
		fz.clear();
		newForceX.clear();
		newForceY.clear();
		newForceZ.clear();
		
	}
	
	data.close();	
}