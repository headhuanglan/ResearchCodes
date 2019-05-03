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
#ifndef MOLECULE_T_H
#define MOLECULE_T_H
#include "atom_t.h"
#include <vector>
#include "lattice_pointer.h"
#include <iostream>
#include "stdlib.h"
using namespace std;



class molecule_t
{
public:
	friend class writeoutput;
	void display();
	molecule_t(int c_length);
	~molecule_t();

private:

	void atoms_gen();  
	static int molcounter;
	int molid_;     //molecular id//
	int c_length_;  //chain length//
    
	int number_bonds_;
	int number_angles_;
	int number_dihdrals_;

    vector<atom_t>   atoms_; //atoms_ array for atom storage//
	

};

#endif //MOLECULE_T_H


