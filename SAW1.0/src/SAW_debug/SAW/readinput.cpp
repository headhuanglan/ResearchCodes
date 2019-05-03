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
#include "readinput.h"

 


 void readinput::ReadInputFile(char* filename)
{
const int LINE_LENGTH=100;char str[LINE_LENGTH];    //temp string  
ifstream   in;
in.open(filename);
cout<<"reading file:"<<filename<<endl;
in>>str;
in>>str;
in>>density;
in>>str;
in>>bond_length;
in>>str;
in>>seed;
in>>str;
in>>mass;
in>>str;
in>>number_chains;
in>>str;

for (int i=0;i<number_chains;i++)
{
int read_chain_length;
in>>read_chain_length;
chain_length[i]=read_chain_length;
}
cout<<" density:"<<density<<"\n number_chains:"<<number_chains<<"\n bond_length:"<<bond_length<<"\n seed:"<<seed<<"\n mass:"<<mass<<endl;
in.close();

int op_sites=0;
for (int j=0;j<number_chains;j++)
{
	op_sites=op_sites+chain_length[j];
	cout<<"chain_length:"<<j+1<<"  "<<chain_length[j]<<endl;
}

               //Initializing  the lattice sites and local_density
lattice_length=static_cast<int>(pow(((op_sites)/density/4),1/3.0))+1;  //increment for lattice//
for (int i=0;i<lattice_length;i++)
		for(int j=0;j<lattice_length;j++)
			for(int k=0;k<lattice_length;k++)
				for(int m=0;m<4;m++)
				   {
					 site[i][j][k][m]=0;
                     local_density[i][j][k][m]=0;
			     	}
 srand (seed);  //Initializaing the random number generator
 }

 