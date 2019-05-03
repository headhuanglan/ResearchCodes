/*************************************************************************************************************
    SAW: self-avoid random walk polymer chian generator for amorphous polymer
    Copyright (C) 2016 Lan Huang

    headhuanglan@tamu.edu  

    SAW is a free code: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

	The basic algorithm is described by the following paper,SAW use local_density and
	angle between atoms to create a resonable amorphous polymer configuration.SAW modified
	the C code to C++(object oriented),so that the simulation of diffrent chain length distribution
	could be possible.

	Hossain, D., Tschopp, M.A., Ward, D.K., Bouvard, J.L., Wang, P., Horstemeyer, M.F.,
	"Molecular dynamics simulations of deformation mechanisms of amorphous polyethylene," 
	Polymer, 51 (2010) 6071-6083.

	Huang, Lan, et al. "Fracture mechanism of amorphous polymers at strain fields."
	Physical Chemistry Chemical Physics 16.45 (2014): 24892-24898.
 *************************************************************************************************************/
#include "molecule_t.h"
#include "readinput.h"
#include "writeoutput.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "lattice_range_define.h"
using namespace std;


 

//global variable//
double density;
int chain_length[MAX_MOLECULES];
int  number_chains;
double  bond_length;
int  seed;
double mass;
int lattice_length;
int site[MAX_X][MAX_Y][MAX_Z][4]; //define the lattice site and occupied atom influenced local_density
int local_density[MAX_X][MAX_Y][MAX_Z][4];

std::vector<molecule_t> molecule_list;






int main(int argc, char *argv[])
{   
// debug data//
//density=0.5;
//chain_length[0]=1000;chain_length[1]=2000;
//number_chains=2;
//bond_length=1.3;
//seed=22222;
//mass=14.2;
//
//lattice_length=static_cast<int>(pow(((3000)/density/4),1/3.0))+1;  //increment for lattice//
//for (int i=0;i<lattice_length;i++)
//		for(int j=0;j<lattice_length;j++)
//			for(int k=0;k<lattice_length;k++)
//				for(int m=0;m<4;m++)
//				   {
//					 site[i][j][k][m]=0;
//                     local_density[i][j][k][m]=0;
//			     	}
 

//###########################read inputfile###################################################################################################
   
     if(argc==1)
	 {
	   char* InputFileName="input.txt";
	   readinput::ReadInputFile(InputFileName);
	 }
	 else
	 {
	  char* InputFileName=argv[1];
	  readinput::ReadInputFile(InputFileName);
	 }
    
     
 
      
//##########################create molecule####################################################################################################
	   
for (int i=0; i<number_chains ; i++) 
 {	molecule_t molecule(chain_length[i]);
    molecule_list.push_back(molecule);
    //molecule.display();
 }

//##########################write lammpsdata###################################################################################################
writeoutput::WriteLammpsData();
    cout<<"  All done!!!"<<endl;
	cout<<"  tips VMD Command: topo readlammpsdata lammps.data"<<endl;
	cout<<"                    pbc box -on"<<endl;
	cout<<"                    mol modcolor 0 0 ResID"<<endl;
	 
    int a;
	cin>>a;

	return 0;
}