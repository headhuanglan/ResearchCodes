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
#include "writeoutput.h"

void writeoutput::WriteLammpsData()
{
    cout<<"writting LammpsData file..."<<endl;
	cout<<"total molecule number:"<<molecule_list.size()<<endl;
   
	ofstream   out("lammps.data");
   	int atoms(0),bonds(0),angles(0),dihedrals(0),number_molecules(0);
	double lattice;
	lattice=bond_length*sqrt(2.0);
    number_molecules=molecule_list.size();

    for(int i=0;i<number_molecules;i++)
	{
		atoms+=molecule_list[i].c_length_;
		bonds+=molecule_list[i].number_bonds_;
		angles+=molecule_list[i].number_angles_;
		dihedrals+=molecule_list[i].number_dihdrals_;        

	}

	out<<"Amorphous polymer generator by L \n";
	out<<"\n";
	out<<atoms<<"  atoms \n";
	out<<bonds<<"  bonds \n";
	out<<angles<<"  angles \n";
	out<<dihedrals<<"  dihedrals \n";
	out<<"\n";
	out<<"1 atom types \n";
	out<<"1 bond types \n";
	out<<"1 angle types \n";
	out<<"1 dihedral types \n";
	out<<"\n";
	out<<"0.0   "<<lattice*lattice_length<<"  xlo xhi \n";
	out<<"0.0   "<<lattice*lattice_length<<"  ylo yhi \n";
	out<<"0.0   "<<lattice*lattice_length<<"  zlo zhi \n";
	out<<"\n";
	out<<"Masses   \n";
	out<<"\n";
	out<<"1   "<<mass<<"\n";
	out<<"\n";
	out<<"Atoms   \n";
	out<<"\n";
     for(int i=0;i<number_molecules;i++)
	{
		for(int j=0;j<molecule_list[i].c_length_;j++)
		{   
			int molid=molecule_list[i].molid_;
			int atomid=molecule_list[i].atoms_[j].atomid_;
			int atomtype=molecule_list[i].atoms_[j].type_;
			int x=molecule_list[i].atoms_[j].x_;
			int y=molecule_list[i].atoms_[j].y_;
			int z=molecule_list[i].atoms_[j].z_;
			int m=molecule_list[i].atoms_[j].m_;
			vector3 xyz(x,y,z);
			vector3 a; //local FCC coordinate 
             if(m==0)
                {  
					a(0,0,0);
                }
              else if(m==1)
                { 
                  a(0.5,0.5,0.0);
                }
              else if(m==2)
                { 
				  a(0.5,0,0.5);
                }
              else
                { 
                  a(0,0.5,0.5);
                }
             vector3 real_xyz=xyz+a;
			 real_xyz=real_xyz*lattice;
			 //                                          charge                                                    //
	out<<atomid<<"    "<<molid<<"   "<<atomtype<<"   "<<"    0      "<<real_xyz.x<<"   "<<real_xyz.y<<"   "<<real_xyz.z<<"\n";
	
		}
	 }
   out<<"\n"; 
   out<<"Bonds \n";
   out<<"\n";
   int bonds_counter(0);
    for(int i=0;i<number_molecules;i++)
	{
		for(int j=0;j<molecule_list[i].number_bonds_;j++)
		{
			bonds_counter++;
			int ii=molecule_list[i].atoms_[j].atomid_;
			int jj=molecule_list[i].atoms_[j+1].atomid_;

    out<<bonds_counter<<"    1     "<<ii<<"   "<<jj<<"\n";     

		}
	}
   out<<"\n"; 
   out<<"Angles \n";
   out<<"\n";
   int angle_counter(0);
    for(int i=0;i<number_molecules;i++)
	{
		for(int j=0;j<molecule_list[i].number_angles_;j++)
		{
			angle_counter++;
			int iii=molecule_list[i].atoms_[j].atomid_;
			int jjj=molecule_list[i].atoms_[j+1].atomid_;
            int kkk=molecule_list[i].atoms_[j+2].atomid_;
    out<<angle_counter<<"    1     "<<iii<<"   "<<jjj<<"   "<<kkk<<"\n";     

		}
	}
   out<<"\n"; 
   out<<"Dihedrals \n";
   out<<"\n";
   int dihedral_counter(0);
    for(int i=0;i<number_molecules;i++)
	{
		for(int j=0;j<molecule_list[i].number_dihdrals_;j++)
		{
			dihedral_counter++;
			int iiii=molecule_list[i].atoms_[j].atomid_;
			int jjjj=molecule_list[i].atoms_[j+1].atomid_;
            int kkkk=molecule_list[i].atoms_[j+2].atomid_;
			int llll=molecule_list[i].atoms_[j+3].atomid_;
    out<<dihedral_counter<<"    1     "<<iiii<<"   "<<jjjj<<"   "<<kkkk<<"   "<<llll<<"\n";     

		}
	}





 

}
 
