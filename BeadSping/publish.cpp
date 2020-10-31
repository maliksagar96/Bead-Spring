#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "molcDyn.h"

using namespace std;

//Print a vector
void printVector(vector<double>& inVector) {
	for(int i = 0;i<inVector.size();i++) {
		cout<<inVector[i]<<"\n";
	}
}

void printVelocity(vector<double>& vx_) {
	int i;
	double avg = 0;
	for(i=0;i<vx_.size();i++) {
		avg = avg + vx_[i];
	}
	cout<<"Net velocity = "<<avg<<endl;
}

/***Create an OVITO style trajectory file containing positions and velocities***/
void trajectoryFile(vector<double>& x_, vector<double>& y_, vector<double>& z_, 
					vector<double>& vx_, vector<double>& vy_, vector<double>& vz_,
					int loopParameter) {

	int i;
	fstream bs;
	bs.open("bs.trj", ios::out);

	bs<<"ITEM: TIMESTEP \n"<<loopParameter<<endl;
	bs<<"ITEM: NUMBER OF ATOMS\n"<<particles<<endl;
	bs<<"ITEM: BOX BOUNDS pp pp pp"<<endl;	
	bs<<-boxLength<<"\t"<<boxLength<<endl;
	bs<<-boxLength<<"\t"<<boxLength<<endl;
	bs<<-boxLength<<"\t"<<boxLength<<endl;
	bs<<"ITEM: ATOMS id x y z vx vy vz"<<endl;
	
	for(i=0;i<x.size();i++) {
		bs<<setprecision(20)<<i<<" "<<x_[i]<<" "<<y_[i]<<" "<<z_[i]<<" "<<vx_[i]<<" "<<vy_[i]<<" "<<vz_[i]<<" "<<endl; 
	}
	bs.close();
}

/*** To write initial configuration into a structure.dat file ***/
void structureFile(vector<double>& x_,vector<double>& y_, vector<double>& z_) {
	fstream positions;
	positions.open("structure.dat", ios::out);
	positions<<"#Initializing at default lattice positions\n\n";
	positions<<particles<<" atoms\n\n";
	positions<<"1 atom types\n\n";
	positions<<0.0<<"\t"<<setprecision(10)<<boxLength<<" xlo xhi\n";
	positions<<0.0<<"\t"<<boxLength<<" ylo yhi\n";
	positions<<0.0<<"\t"<<boxLength<<" zlo zhi\n\n";
	positions<<"Atoms\n\n";

	for(int i = 0;i<x_.size();i++) {
		positions<<setw(3)<<i+1<<setw(3)<<1<<setw(21)<<x_[i]<<setw(21)<<y_[i]<<setw(21)<<z_[i]<<endl;
	}	
}

void writeData(){
	int a = 10, b = 15, c = 120;
	ofstream dataFile;
	dataFile.open ("example.dat");
	dataFile << "Writing this to a file.\n";
	dataFile << a <<" "<<b<<" "<<c;
	dataFile.close(); 
}

void printConstants() {
	cout<<"Sigma = "<<sigma<<endl;
	cout<<"Rc = "<<rc<<endl;
	cout<<"Urc = "<<urc<<endl;
	cout<<"Frc = "<<frc<<endl;
}
