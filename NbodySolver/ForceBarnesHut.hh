#ifndef _FORCEBARNESHUT_HH_
#define _FORCEBARNESHUT_HH_

#include "ForceCalculator.hh"
#include "Body.hh"
#include <stdint.h>
class MortonKeyCalculator {
public:

	MortonKeyCalculator(Body *bodies, int N);

	uint32_t x(double xval) const {
		return (uint32_t) (0xffffffff * (xval - xmin) / max_range);
	}
	;

	uint32_t y(double yval) const {
		return (uint32_t) (0xffffffff * (yval - ymin) / max_range);
	}

	double cellWidth(int level) const {
		return max_range / (1 << (level - 1));
	}

	std::ostream& printKey(std::ostream &os, const Body &b) const;

	//Returns true if a is less than b in a MortonKey ordering.
	//Note that the y coordinate is more significant.
	bool operator()(const Body &a, const Body &b) const;

	double getXmin() {
		return xmin ;
	}
	double getYmin() {
		return ymin ;
	}

	int no_of_threads;
private:
	double xmin;
	double ymin;
	double max_range;
};



/*-----------------------------------------*/
class QTree {
public:
	QTree(double xmin_,double ymin_,double width_,QTree* parent_)

	{
		parent=parent_;
		width=width_;
        xmin=xmin_;
		ymin=ymin_;
        sub=NULL;
        m=0.0;
        x=0.0;
        y=0.0;
        isleaf=false;
        body=NULL;
        child0=NULL;
        child1=NULL;
        child2=NULL;
        child3=NULL;
	}

	inline void initsub(){
		this->sub= (class QTree*) calloc(4,sizeof(class QTree));  // pre allocate memory for children
	}

	//QTree(){x=m=y=0;isleaf=false;width=xmin=ymin=0;sub=parent=child0=child1=child2=child3=NULL;body=NULL;}

   ~QTree(){ if (sub!=NULL) free(sub);}

	bool contain(Body* body);
	int  childcontain(Body* body2test); //return which child to insert
	QTree* insert(Body* body2insert);

	QTree* parent;
	QTree* sub;

	QTree* child0;
	QTree* child1;
	QTree* child2;
	QTree* child3;

	Body*  body;

	double width;
    double xmin;
    double ymin;

	bool isleaf;

	double x;
	double y;
	double m;


};

class ForceBarnesHut: public ForceCalculator {

public:
	ForceBarnesHut(Body *body, int N, double theta);

	virtual void operator()(Body *pulled);

 ~ForceBarnesHut() {

	}
	;


QTree* mytree;



private:
	Body *bodies_;
	const int N_;
	double theta_;
};

#endif
