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
#include "atom_t.h"


double sqrt_2_2=sqrt(2.0)/2.0;
int atom_t::atomcounter=0;


 


atom_t::atom_t(int x, int y, int z, int m, int type, int q) 
{    
	 atomcounter++;
	 site[x][y][z][m]=1;
	 atomid_=atomcounter;
     x_=x;     //atom position defined in FCC by X Y Z M//
	 y_=y;
	 z_=z;
	 m_=m;
     type_=type;  //atom type//
     q_=q;     // charge of atom//


	 find_neighbor_sites();  //find 12 neighbor sites for this atoms site[x][y][z][m] Then local_density++;
}


 
	

 



void atom_t::find_neighbor_sites()
{

int x=this->x_;
int y=this->y_;
int z=this->z_;
int m=this->m_;
vector3 shift2real(x,y,z); //shift vector from (primitive cell local coordinate 0 0 0) to (real coordinate) 
vector3 a;                //local FCC coordinate 

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
 //find the neighbor sites//
int ii,jj,kk,mm;
vector3 aa;   //define FCC primitive cell coordinate
vector3 va;   //the 
vector3 shift; //shift vector3 for neighbor sites search
vector3 neighbor_site_x_y_z; //store for neighbor x y z data
int neighbor_sites_counter=0;
			  for(ii=0;ii<3;ii++)
			  {
                  for(jj=0;jj<3;jj++)
                    {
                      for(kk=0;kk<3;kk++)
                        {
                          for(mm=0;mm<4;mm++)
                            {
                              if(mm==0)
                                {
                                  aa(0.0,0.0,0.0);
                                }
                              else if(mm==1)
                                { 
                                  aa(0.5,0.5,0.0);
                                }
                              else if(mm==2)
                                { 
                                  aa(0.5,0.0,0.5);
                                }
                              else if(mm==3)
                                { 
                                  aa(0.0,0.5,0.5);
                                }
							  double tempx=ii-1;
							  double tempy=jj-1;
							  double tempz=kk-1;
                              aa=aa+shift(tempx,tempy,tempz); //ii range 0 1 2 ; ii-1 range -1 0 1  search the neighbor sites
							  //cout<<"aa:xyz "<<aa.x<<" "<<aa.y<<" "<<aa.z<<" "<<endl;
							  va=aa-a;
							  double va_mag=va.magnitude();
							  //cout<<"va_amg:"<<va_mag<<endl;
							  
							  if(va_mag<=sqrt_2_2+0.01&&va_mag>=0.1) //within the sqrt_2_2 range is the nearest 12 sites  (exclude the atom itself)
							  {      
							        vector3 v3=(shift2real+shift);
									//lattice pbc conditon//
                                    if(v3.x==lattice_length){v3.x=0;}
	                                if(v3.y==lattice_length){v3.y=0;}
	                                if(v3.z==lattice_length){v3.z=0;}
                                    if(v3.x<0){v3.x=lattice_length-1;}
	                                if(v3.y<0){v3.y=lattice_length-1;}
	                                if(v3.z<0){v3.z=lattice_length-1;}
                                   //lattice pbc conditon//

  						            //if(v3.x!=this->x_&&v3.y!=this->y_&&v3.z!=this->z_&&mm!=this->m_) //exclude the atom itself
									//{
                                         neighbor_site_x_y_z=v3;
									     neighbor_sites[neighbor_sites_counter].x=static_cast<int>(neighbor_site_x_y_z.x);
									     neighbor_sites[neighbor_sites_counter].y=static_cast<int>(neighbor_site_x_y_z.y);
								         neighbor_sites[neighbor_sites_counter].z=static_cast<int>(neighbor_site_x_y_z.z);
									     neighbor_sites[neighbor_sites_counter].m=mm;
										 neighbor_sites_vector[neighbor_sites_counter]=va;
    									 local_density[neighbor_sites[neighbor_sites_counter].x][neighbor_sites[neighbor_sites_counter].y][neighbor_sites[neighbor_sites_counter].z][mm]++;
									     neighbor_sites_counter=neighbor_sites_counter+1;
									//}
							  }
                              
						  }
					  }
				  }

                    
			  }

			  if(neighbor_sites_counter!=12){cout<<"Something is wrong,neighbor sites are not 12."<<endl; exit(1);}
     //cout<<"neighbor site found:"<<neighbor_sites_counter<<endl;
     //for(int i=0;i<neighbor_sites_counter;i++)
	 //cout<<"sites x"<<neighbor_sites[i].x<<" sites y"<<neighbor_sites[i].y<<" sites z"<<neighbor_sites[i].z<<" sites m"<<neighbor_sites[i].m<<endl;

}





void atom_t::display()
{
 cout<<"atomid="<<atomid_<<" x="<<x_<<" y="<<y_<<" z="<<z_<<" m="<<m_<<" type="<<type_<<" q="<<q_<<endl;

}

 

void atom_t::clean()
{
cout<<"due to atom deletion, restore local_density."<<endl;
	for(int i=0;i<12;i++)
		local_density[neighbor_sites[i].x][neighbor_sites[i].y][neighbor_sites[i].z][neighbor_sites[i].m]--;
}

atom_t::~atom_t()
{
	
}
