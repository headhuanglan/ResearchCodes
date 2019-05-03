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


int molecule_t::molcounter=0;  


molecule_t::molecule_t(int c_length)
{  
	molcounter++;
    molid_=molcounter;
	c_length_=c_length;

    number_bonds_=c_length_-1;
	number_angles_=c_length_-2;
	if(number_angles_<0){number_angles_=0;}
	number_dihdrals_=c_length_-3;
	if(number_dihdrals_<0){number_dihdrals_=0;}

    atoms_gen();

}


void molecule_t::display()
{

 
cout<<"molid="<<molid_<<" c_length="<<c_length_<<endl;
cout<<c_length_<<" atoms listed in chain"<<endl;
for(unsigned int i=0;i<atoms_.size();i++){atoms_[i].display();} 
}


void molecule_t::atoms_gen()
{
    
   //define the first atom in the molecue chain (abandon the occupied site and the trapped site)
	  int i,j,k,m,reset(0);
	   
   do{
      i=static_cast<int>(lattice_length*((double)rand()/(double)RAND_MAX));     
      j=static_cast<int>(lattice_length*((double)rand()/(double)RAND_MAX));     
      k=static_cast<int>(lattice_length*((double)rand()/(double)RAND_MAX));     
	  m=static_cast<int>(4*((double)rand()/(double)RAND_MAX)); 
      } while(site[i][j][k][m]==1||site[i][j][k][m]==2);

      atom_t atom(i,j,k,m);
	  atoms_.push_back(atom);
  //define the rest atoms in chain //  nic==atom number in chain 
	  for (int nic=1;nic<c_length_;nic++)
	  {
		  
		  int density=0;  //each atom has 12 nearest neighbors, If one site is surrounded by the 12 neighbors, the occupied weight is 12 for this site, the free volume weight is sum(12-local_density), 
		  vector<int> qualified_neighbor_sites_lable;
		  for(int i=0;i<12;i++)  
		  {         
			     //criteria: 1 neighbor list site=0 ; 2. angle of neigbor site is correct//  atoms_.back() return current atom;  atoms_.at(atoms_.size()-2) return the second last atom;
                 //n_s coordinate for current atom
             int n_s_x=atoms_.back().neighbor_sites[i].x;
		     int n_s_y=atoms_.back().neighbor_sites[i].y;
			 int n_s_z=atoms_.back().neighbor_sites[i].z;
		     int n_s_m=atoms_.back().neighbor_sites[i].m;
				 
			  if(site[n_s_x][n_s_y][n_s_z][n_s_m]==0||(nic==c_length_-1&&site[n_s_x][n_s_y][n_s_z][n_s_m]==2))   //unoccupied cites and site==2 for the last chain atom are qualified sites
				 {
					 if(nic==1){qualified_neighbor_sites_lable.push_back(i);}
					 else{
					          int nsl=atoms_.at(atoms_.size()-2).neighbor_sites_lable;
					          vector3 v_old=atoms_.at(atoms_.size()-2).neighbor_sites_vector[nsl];
							  double angle=vector3::angle(v_old,atoms_.back().neighbor_sites_vector[i]);
                                   if ((angle<-60.0&&angle>-150.0)||(angle>60.0&&angle<150.0))
								   {   //cout<<"angle:"<<angle<<endl;//
									   qualified_neighbor_sites_lable.push_back(i);
								   }

				          }  
			     }
				 

		  }//scan the nearest neigbor over
            //if no avaliable nearest neighbor 
         /* for(int i=0;i<qualified_neighbor_sites_lable.size();i++)
		  {cout<<"Chain "<<this->molid_<<"atom"<<nic<<" Qualified "<<qualified_neighbor_sites_lable[i]<<endl;
		  }
		  cout<<"####################"<<endl;*/

           if(qualified_neighbor_sites_lable.size()==0){ 
			                                                  if(reset==1)
                                                                {
                                                                  cout<<"has already been reset and still can not find a neighbor. Chain NO.:"<<molid_<<"  atom NO.:"<<nic<<endl;
                                                                  exit(1);
															     } 
                                                              cout<<"there is not a free neighbor site and must reset"<<endl;
				                                              reset=1;
															  site[atoms_.back().x_][atoms_.back().y_][atoms_.back().z_][atoms_.back().m_]=2;   //this site is trapped, but suitable for the last atom in chain
															  atom_t::atomcounter--;
															  atoms_.back().clean();   //restore the local_density around the atom
															  atoms_.pop_back();      
															  nic=nic-2;
															  continue;
			                                              }
          //randomly choose the atom creation site by the percent of free volume
		   //first let's cal the total free volume around a atom
          for(unsigned int i=0;i<qualified_neighbor_sites_lable.size();i++)
		  {
			  int mark=qualified_neighbor_sites_lable[i];
		  
			  density=density+(12-local_density[atoms_.back().neighbor_sites[mark].x][atoms_.back().neighbor_sites[mark].y][atoms_.back().neighbor_sites[mark].z][atoms_.back().neighbor_sites[mark].m]);
		  }


		  int point=static_cast<int>(density*((double)rand()/(double)RAND_MAX));
          int prob_low=0;
		  int prob_high=(12-local_density[atoms_.back().neighbor_sites[qualified_neighbor_sites_lable[0]].x][atoms_.back().neighbor_sites[qualified_neighbor_sites_lable[0]].y][atoms_.back().neighbor_sites[qualified_neighbor_sites_lable[0]].z][atoms_.back().neighbor_sites[qualified_neighbor_sites_lable[0]].m]);
		  int qualified_neighbor_number=qualified_neighbor_sites_lable.size();
		  for(int nn=0;nn<qualified_neighbor_number;nn++)
		  {
			  int label=qualified_neighbor_sites_lable[nn];

		      if(point>=prob_low&&point<prob_high)
		          {
					//cout<<"qualified neighbor label:"<<label<<"prob_low and high"<<prob_low<<" "<<prob_high<<endl;
					int xx=atoms_.back().neighbor_sites[label].x;
					int yy=atoms_.back().neighbor_sites[label].y;
					int zz=atoms_.back().neighbor_sites[label].z;
					int mm=atoms_.back().neighbor_sites[label].m;
					atoms_.back().neighbor_sites_lable=label;
  				    atom_t atom(xx,yy,zz,mm);
		            atoms_.push_back(atom);
					reset=0;  //if we succesfuly creat the atom, the reset should by restore to 0
				    break;
		          }
			  else
			     { 
					 if(nn==qualified_neighbor_number-1)    //this is alrealy the last avaliable site to creat atom, prob_high is already the higest probility
					 {
						if(reset==1)
                             {              
								 cout<<"has already been reset and still can not find a neighbor. Chain NO.:"<<molid_<<"  atom NO.:"<<nic<<endl;
                              exit(1);
			                  } 
					   cout<<"Cannot find the site to create atom, must reset."<<endl;
					   reset=1;
					   site[atoms_.back().x_][atoms_.back().y_][atoms_.back().z_][atoms_.back().m_]=2;   //this site is trapped, but suitable for the last atom in chain
				       atom_t::atomcounter--;
					   atoms_.back().clean();    
					   atoms_.pop_back();      
					   nic=nic-2;
					   break;
					 }
					 else
					 {
		               prob_low=prob_high;
				       prob_high=prob_high+(12-local_density[atoms_.back().neighbor_sites[qualified_neighbor_sites_lable[nn+1]].x][atoms_.back().neighbor_sites[qualified_neighbor_sites_lable[nn+1]].y][atoms_.back().neighbor_sites[qualified_neighbor_sites_lable[nn+1]].z][atoms_.back().neighbor_sites[qualified_neighbor_sites_lable[nn+1]].m]);
			         } 
			     }
		  }    
            

		   

	  
	  }//chain_creation_over
 




	 
}


molecule_t::~molecule_t()
{
 
}
