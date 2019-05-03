#ifdef _OPENMP
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#endif
#include <iostream>
#include "MortonKeyCalculator.hh"

MortonKeyCalculator::MortonKeyCalculator(Body *bodies, int N) {
	//Find in parallel the minimum x and y values among bodies and 
	//the range of values (i.e. max(xmax - xmin, ymax - ymin)).

	//You may use stl's min_element and max_element functions if
	//you like.

	//Your code here.

	double * xmin_data;
	double * xmax_data;
	double * ymin_data;
	double * ymax_data;
	int data_length;

#pragma omp parallel shared(bodies,N,data_length,xmin_data,xmax_data,ymin_data,ymax_data)
	{
		int tid = omp_get_thread_num();
		int P = omp_get_num_threads();
#pragma omp master
		{
			data_length = P;
			xmin_data = (double*) malloc(P * sizeof(double));
			xmax_data = (double*) malloc(P * sizeof(double));
			ymin_data = (double*) malloc(P * sizeof(double));
			ymax_data = (double*) malloc(P * sizeof(double));
		}
#pragma omp barrier
		double xmin_local = bodies[0].x();
		double xmax_local = bodies[0].x();
		double ymin_local = bodies[0].y();
		double ymax_local = bodies[0].y();
		for (int i = tid * N / P; i <= (tid + 1) * N / P - 1; i++) {
			//#pragma omp critical
			//std::cout<<tid<<"|"<<i<<std::endl;
			if (bodies[i].x() < xmin_local) {
				xmin_local = bodies[i].x();
			}
			if (bodies[i].x() > xmax_local) {
				xmax_local = bodies[i].x();
			}
			if (bodies[i].y() < ymin_local) {
				ymin_local = bodies[i].y();
			}
			if (bodies[i].y() > ymax_local) {
				ymax_local = bodies[i].y();
			}
		}			//endfor

		xmin_data[tid] = xmin_local;
		xmax_data[tid] = xmax_local;
		ymin_data[tid] = ymin_local;
		ymax_data[tid] = ymax_local;

	}

	xmin = xmin_data[0];
	ymin = ymin_data[0];
	double xmax = xmax_data[0];
	double ymax = ymax_data[0];
	for (int i = 1; i < data_length; i++) {
		if (xmin_data[i] < xmin)
			xmin = xmin_data[i];
		if (ymin_data[i] < ymin)
			ymin = ymin_data[i];
		if (xmax_data[i] > xmax)
			xmax = xmax_data[i];
		if (ymax_data[i] > ymax)
			ymax = ymax_data[i];
	}

	max_range = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);

//	std::cout << "xmin xmax ymin ymax max_range" << "|" << xmin << "|" << xmax
//			<< "|" << ymin << "|" << ymax << "|" << max_range << std::endl;

	free(xmin_data);
	free(xmax_data);
	free(ymin_data);
	free(ymax_data);
}
;

std::ostream& MortonKeyCalculator::printKey(std::ostream &os,
		const Body &b) const {
	uint32_t bx = x(b.x());
	uint32_t by = y(b.y());

	for (int i = 31; i >= 0; i--) {
		os << ((by >> i) % 2);
		os << ((bx >> i) % 2);
	}
	os << std::endl;

	return os;
}



bool MortonKeyCalculator::operator()(const Body &a, const Body &b) const {
	//Return true is a is less than b in a Morton Key ordering; 
	//otherwise, return false.
	uint32_t ax = x(a.x());
	uint32_t ay = y(a.y());
	uint32_t bx = x(b.x());
	uint32_t by = y(b.y());
	int ith_morton_val_a;
	int ith_morton_val_b;
	for (int i = 31; i >= 0; i--) {
		 ith_morton_val_a=(((ay >> i) % 2) << 1) + (ax >> i) % 2;
		 ith_morton_val_b=(((by >> i) % 2) << 1) + (bx >> i) % 2;

       if(ith_morton_val_a<ith_morton_val_b){
    	   return true;
       }else if(ith_morton_val_a==ith_morton_val_b){
    	   continue;
       }else{
    	   return false;
       }


	}

	return false;

}
