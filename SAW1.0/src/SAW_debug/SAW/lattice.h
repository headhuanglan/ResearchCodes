#ifndef LATTICE_H
#define LATTICE_H
#include <cmath>

#include <iostream>
using namespace std;

extern double density;
extern int chain_length;
extern int number_chains;
extern int lattice_length;
#define MAX_X         300
#define MAX_Y         300
#define MAX_Z         300

class lattice
{
public:
   
	 
	static void latticeInit();
	static void display();

	static int site_[MAX_X][MAX_Y][MAX_Z][4];
    static int local_density_[MAX_X][MAX_Y][MAX_Z][4];

};
#endif //LATTIC_H


int i,j,k,m;
for( i = 0 ; i <lattice_length ; i ++ ){
	for( j = 0 ; j <lattice_length ; j ++ ){
		for( k = 0 ; k <lattice_length ; k ++ ){
		 for( m = 0 ; m<4; m ++)
			 lattice::site_[lattice_length][lattice_length][lattice_length]=0;
		}
	}
}
