/*************************************************************************************************************
    SAW: self-avoid random walk polymer chian generator for amorphous polymer
    Copyright (C) 2016 Lan Huang

    headhuanglan@tamu.edu  

    SAW is a free code: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

	The basic algrithom is described by the following paper,SAW use local_density and
	angle between atoms to create a resonable amorphous polymer configuration.SAW modified
	the C code to C++(object oriented),so that the simulation of diffrent chain length distribution
	could be possible.

	Hossain, D., Tschopp, M.A., Ward, D.K., Bouvard, J.L., Wang, P., Horstemeyer, M.F.,
	"Molecular dynamics simulations of deformation mechanisms of amorphous polyethylene," 
	Polymer, 51 (2010) 6071-6083.

	Huang, Lan, et al. "Fracture mechanism of amorphous polymers at strain fields."
	Physical Chemistry Chemical Physics 16.45 (2014): 24892-24898.
 *************************************************************************************************************/
#ifndef ATOM_T_H
#define ATOM_T_H
#include <vector>
#include "vector3.h"
#include "vector4.h"
#include "lattice_pointer.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>
using namespace std;


class atom_t
{
public:
	friend class writeoutput;
	friend class molecule_t;
	void display();
	atom_t(int x=-1, int y=-1, int z=-1, int m=-1, int type=1, int q=-1);   
    ~atom_t();
private:
	 void find_neighbor_sites();
	 void clean(); //clean site and local_density 
	 static int atomcounter;
	 int atomid_;//atom id
     int x_;     //atom position defined in FCC by X Y Z M//
	 int y_;
	 int z_;
	 int m_;
     int type_;  //atom type defalt 1//
     int q_;     // charge of atom not used yet//
     
	 int neighbor_sites_lable;    // the  neighbor_sites_lable for which a new atom is created
	 vector4 neighbor_sites[12];  //per atom should have 12 nearest lattice site 
	 vector3 neighbor_sites_vector[12];//store the vectors from the atom_center to neighbor sites 
	 

	 
};
#endif //ATOM_T_H