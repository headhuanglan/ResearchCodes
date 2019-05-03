#ifdef _OPENMP
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#endif
#include <iostream>
#include "ForceBarnesHut.hh"
#include "Timer.hh"





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
			no_of_threads = P;
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
		ith_morton_val_a = (((ay >> i) % 2) << 1) + (ax >> i) % 2;
		ith_morton_val_b = (((by >> i) % 2) << 1) + (bx >> i) % 2;

		if (ith_morton_val_a < ith_morton_val_b) {
			return true;
		} else if (ith_morton_val_a == ith_morton_val_b) {
			continue;
		} else {
			return false;
		}

	}

	return false;

}
/*------------------------------------------------------------------------------*/
bool QTree::contain(Body* body) {

	if (body->y() <= this->ymin + this->width && body->y() >= this->ymin) {
		if (body->x() <= xmin + this->width && body->x() >= this->xmin) {
			return true;
		}
	}

	return false;
}
//find which children to insert body return 0 1 2 3
//do a test insert do not do actural insertion
int QTree::childcontain(Body* body2test) {

	double testx = body2test->x();
	double testy = body2test->y();

	if (testy <= ymin + width / 2.0) {
		if (testx <= xmin + width / 2.0) {
			return 0;
		} else {
			return 1;
		}
	} else {
		if (testx <= xmin + width / 2.0) {
			return 2;
		} else {
			return 3;
		}

	}

}
//success return the node inserted  unsuccess return NULL
QTree* QTree::insert(Body* body2insert) {
	if (this->contain(body2insert)) {
		if (this->sub != NULL) {	//this has children

			//update this info
			//		double newm = this->m + body2insert->m();
			//		this->x = (this->x * this->m + body2insert->x() * body2insert->m())
			//				/ newm;
			//		this->y = (this->y * this->m + body2insert->y() * body2insert->m())
			//				/ newm;
			//		this->m = newm;

			int decide = this->childcontain(body2insert);
			switch (decide) {
			case 0: {
				if (child0 == NULL)
					child0 = new (this->sub) QTree(xmin, ymin, width / 2.0,
							this);
				child0->insert(body2insert);
				return child0;
			}
			case 1: {
				if (child1 == NULL)
					child1 = new (this->sub + 1) QTree(xmin + width / 2.0, ymin,
							width / 2.0, this);
				child1->insert(body2insert);
				return child1;

			}
			case 2: {
				if (child2 == NULL)
					child2 = new (this->sub + 2) QTree(xmin, ymin + width / 2.0,
							width / 2.0, this);
				child2->insert(body2insert);
				return child2;

			}
			case 3: {
				if (child3 == NULL)
					child3 = new (this->sub + 3) QTree(xmin + width / 2.0,
							ymin + width / 2.0, width / 2.0, this);
				child3->insert(body2insert);
				return child3;

			}
			default: {
				std::cout << "Error body2insert is not in any of child region"
						<< std::endl;
				exit(1);
			}
			}

		} else if (this->body != NULL) {	//this contain body

			//update this info
			//	double newm = this->m + body2insert->m();
			//	this->x = (this->x * this->m + body2insert->x() * body2insert->m())
			//			/ newm;
			//	this->y = (this->y * this->m + body2insert->y() * body2insert->m())
			//			/ newm;
			//	this->m = newm;

			//move ond to newchild insert new
			this->initsub();
			int decideold2go = this->childcontain(this->body);
			switch (decideold2go) {
			case 0: {
				//placement new
				child0 = new (this->sub) QTree(xmin, ymin, width / 2.0, this);
				child0->insert(this->body);
				this->body = NULL;
				this->isleaf = false;
				break;

			}
			case 1: {
				child1 = new (this->sub + 1) QTree(xmin + width / 2.0, ymin,
						width / 2.0, this);
				child1->insert(this->body);
				this->body = NULL;
				this->isleaf = false;
				break;

			}
			case 2: {
				child2 = new (this->sub + 2) QTree(xmin, ymin + width / 2.0,
						width / 2.0, this);
				child2->insert(this->body);
				this->body = NULL;
				this->isleaf = false;
				break;

			}
			case 3: {
				child3 = new (this->sub + 3) QTree(xmin + width / 2.0,
						ymin + width / 2.0, width / 2.0, this);
				child3->insert(this->body);
				this->body = NULL;
				this->isleaf = false;
				break;

			}
			default: {
				std::cout
						<< "Error after divide old body  is not in any of child region"
						<< std::endl;
				exit(1);
			}
			}	//finsh insert old to child

			//insert new
			int decidenew2go = this->childcontain(body2insert);
			switch (decidenew2go) {
			case 0: {
				if (child0 == NULL)
					child0 = new (this->sub) QTree(xmin, ymin, width / 2.0,
							this);
				child0->insert(body2insert);
				return child0;

			}
			case 1: {
				if (child1 == NULL)
					child1 = new (this->sub + 1) QTree(xmin + width / 2.0, ymin,
							width / 2.0, this);
				child1->insert(body2insert);
				return child1;

			}
			case 2: {
				if (child2 == NULL)
					child2 = new (this->sub + 2) QTree(xmin, ymin + width / 2.0,
							width / 2.0, this);
				child2->insert(body2insert);
				return child2;

			}
			case 3: {
				if (child3 == NULL)
					child3 = new (this->sub + 3) QTree(xmin + width / 2.0,
							ymin + width / 2.0, width / 2.0, this);
				child3->insert(body2insert);
				return child3;

			}
			default: {
				std::cout << "NEWBODYINFO m:" << body2insert->m() << "x:"
						<< body2insert->x() << "y:" << body2insert->y()
						<< " decidenew2go:" << decidenew2go
						<< "this tree info m:" << m << "width" << width
						<< "xmin:" << xmin << "ymin" << ymin
						<< "\nError Cannot insert newbody into any of child region\n"
						<< "This contain body2insert?:"
						<< this->contain(body2insert) << "\nchild who 2 insert?"
						<< this->childcontain(body2insert) << std::endl;
				exit(1);
			}
			}

		} else {
			this->body = body2insert;
			this->isleaf = true;
					this->x = body2insert->x();
					this->y = body2insert->y();
					this->m = body2insert->m();
			return this;
		}
	} else {
		return NULL;
	}

}
/*------------------------------------------------------------------------------*/

void print_mytree(QTree* mytree) {

	if (true) {	//mytree->isleaf){
		std::cout << "mytree width" << mytree->width << "m:" << mytree->m
				<< "x:" << mytree->x << "y:" << mytree->y << std::endl;
	}

	if (mytree->child0 != NULL) {
		std::cout << "child0" << std::endl;
		print_mytree(mytree->child0);
	}
	if (mytree->child1 != NULL) {
		std::cout << "child1" << std::endl;
		print_mytree(mytree->child1);
	}
	if (mytree->child2 != NULL) {
		std::cout << "child2" << std::endl;
		print_mytree(mytree->child2);
	}
	if (mytree->child3 != NULL) {
		std::cout << "child3" << std::endl;
		print_mytree(mytree->child3);
	}

}

/*-----------------------------------------------------------------------------*/

QTree* mergeTrees(QTree* tree, QTree* other) {

	if (tree->isleaf) {
		if (other->isleaf) {	//both leaf
			tree->insert(other->body);
			other->body = NULL;
			return tree;
		} else {	//tree leaf other isnot
			other->insert(tree->body);
			tree->body = NULL;
			return other;

		}
	} else {
		if (other->isleaf) {	//other is leaf this is not
			tree->insert(other->body);
			other->body = NULL;
			return tree;

		} else {	//both not leaf

			double newm = tree->m + other->m;
			tree->x = (tree->x * tree->m + other->x * other->m) / newm;
			tree->y = (tree->y * tree->m + other->y * other->m) / newm;
			tree->m = newm;

#pragma omp parallel
			{
#pragma omp sections nowait
				{
					/*----*/
#pragma omp section
					{
						if (tree->child0 != NULL) {
							if (other->child0 != NULL)
								tree->child0 = mergeTrees(tree->child0,
										other->child0);

						} else {
							if (other->child0 != NULL)
								//treechild0==null  other !=null
								tree->child0 = other->child0;
						}
					}
					/*----*/
#pragma omp section
					{
						if (tree->child1 != NULL) {
							if (other->child1 != NULL)
								tree->child1 = mergeTrees(tree->child1,
										other->child1);

						} else {
							if (other->child1 != NULL)
								//treechild0==null  other !=null
								tree->child1 = other->child1;
						}
					}
					/*----*/
#pragma omp section
					{
						if (tree->child2 != NULL) {
							if (other->child2 != NULL)
								tree->child2 = mergeTrees(tree->child2,
										other->child2);

						} else {
							if (other->child2 != NULL)
								//treechild0==null  other !=null
								tree->child2 = other->child2;
						}
					}
					/*----*/
#pragma omp section
					{
						if (tree->child3 != NULL) {
							if (other->child3 != NULL)
								tree->child3 = mergeTrees(tree->child3,
										other->child3);

						} else {
							if (other->child3 != NULL)
								//treechild0==null  other !=null
								tree->child3 = other->child3;
						}
					}
				}	//endsection
			}	//endparallel
			return tree;

		}

	}

}
void trasverse(QTree* mytree) {
	if (mytree == NULL)
		return;

	if (mytree->child0 != NULL)
		trasverse(mytree->child0);

	if (mytree->child1 != NULL)
		trasverse(mytree->child1);

	if (mytree->child2 != NULL)
		trasverse(mytree->child2);

	if (mytree->child3 != NULL)
		trasverse(mytree->child3);

	if (mytree->isleaf) {
		if (mytree->body != NULL) {
			mytree->m = mytree->body->m();
			mytree->x = mytree->body->x();
			mytree->y = mytree->body->y();
		}

	} else {
		//mytree is middle node

		double newm = 0;
		double newx = 0;
		double newy = 0;
		if (mytree->child0 != NULL) {
			newm += mytree->child0->m;
			newx += mytree->child0->m * mytree->child0->x;
			newy += mytree->child0->m * mytree->child0->y;
		}
		if (mytree->child1 != NULL) {
			newm += mytree->child1->m;
			newx += mytree->child1->m * mytree->child1->x;
			newy += mytree->child1->m * mytree->child1->y;
		}

		if (mytree->child2 != NULL) {
			newm += mytree->child2->m;
			newx += mytree->child2->m * mytree->child2->x;
			newy += mytree->child2->m * mytree->child2->y;
		}

		if (mytree->child3 != NULL) {
			newm += mytree->child3->m;
			newx += mytree->child3->m * mytree->child3->x;
			newy += mytree->child3->m * mytree->child3->y;
		}

		mytree->m = newm;
		mytree->x = newx / newm;
		mytree->y = newy / newm;
	}

}
/*------------------------------------------------------------------------------*/
ForceBarnesHut::ForceBarnesHut(Body *bodies, int N, double theta) :
		bodies_(bodies), N_(N), theta_(theta) {
	//Build the quad tree and do any other needed initialization here.

	MortonKeyCalculator mkcalc(bodies_, N_);
#ifdef _OPENMP
	__gnu_parallel::sort(bodies_, bodies_ + N_, mkcalc);
#else
	std::sort(bodies_, bodies_ + N_, mkcalc);
#endif

	/*---------debug-------*/
	 mytree = new QTree(mkcalc.getXmin(), mkcalc.getYmin(), mkcalc.cellWidth(1),  NULL);
	 //simple insertion
	 //  for (int i = 1; i < N_; i++) {
	 // 		mytree->insert(&bodies_[i]);
	 // 	}


	 //book marking insertion

	 QTree* curnode = mytree->insert(&bodies_[0]);
	 for (int i = 1; i < N_; i++) {

	 if (curnode != NULL)
	 curnode =  curnode->parent;

	 if (curnode == NULL) {
	 curnode = mytree->insert(&bodies_[i]);
	 } else {
	 curnode=curnode->insert(&bodies_[i]);
	 }

	 }

	 trasverse(mytree);

	/* ---*/

	/*--------------parallel subtree building----------*/
/*
	int P = mkcalc.no_of_threads;
	QTree* subtrees[P];
#pragma omp parallel
	{

#pragma omp for nowait
		for (int i = 0; i < P; i++) {
			subtrees[i] = new QTree(mkcalc.getXmin(), mkcalc.getYmin(),
					mkcalc.cellWidth(1), NULL);

			QTree* curnode = subtrees[i]->insert(&bodies_[i * N_ / P]);

			for (int j = i * N_ / P + 1; j <= (i + 1) * N_ / P - 1; j++) { //insert in subtree

				//naive method (subtrees+i)->insert(&bodies_[j]);
				//book marking parent
				if (curnode != NULL)
					curnode = curnode->parent;
				if (curnode == NULL) {
					curnode = subtrees[i]->insert(&bodies_[j]);
				} else {
					curnode = curnode->insert(&bodies_[j]);
				}

			}

		}				//all bodies belong to ith subtree is inserted

		//trasverse to assign middle node x y m info

#pragma omp  for
		for (int i = 0; i < P; i++) {
			trasverse(subtrees[i]);
		}

	}				//endparallel


*/
	/*
	 std::cout << "@@@subtree0" << std::endl;
	 print_mytree(subtrees[0]);
	 std::cout << "@@@subtree1" << std::endl;
	 print_mytree(subtrees[1]);
	 std::cout << "@@@subtree2" << std::endl;
	 print_mytree(subtrees[2]);
	 std::cout << "@@@subtree3" << std::endl;
	 print_mytree(subtrees[3]);
	 */
	/*--------------parallel merge subtree ----------*/
//	std::cout
//			<< "@@@parallel merge subtree@@@p@@@p@@@p@@@p@@@p@@@p@@@p@@@p@@@p@@@p@@@p "
//			<< std::endl;

/*
	int loopcounter = 0;
	for (int m = P; m > 0; m = ceil(m / 2)) {
#pragma omp parallel for shared(P,m,loopcounter)  schedule(dynamic)
		for (int i = 0; i < P; i += 2 * (int(P / m))) {

			if ((i + (1 << loopcounter)) < P) {
				//std::cout << "@@@@@@m(" << m << ")" << "i(" << i << "|||||" << i
				//		<< "-" << i + pow(2, loopcounter) << std::endl;
				int mergeID = i + (1 << loopcounter);
				subtrees[i] = mergeTrees(subtrees[i], subtrees[mergeID]);
				// #pragma omp critical
				//	{std::cout<<"@@@@Mergedtree:"<<i<<"from:"<<i<<"-"<<mergeID<<std::endl;
				//  print_mytree(subtrees[i]);}
			}
		}

		loopcounter++;
	}

	mytree = subtrees[0];
*/
//	 std::cout<<"@@@@@@@@@@@@@@@@@@@mytree@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
 	// print_mytree(mytree);

}

void BH(Body* pulled, QTree* tree, double &theta) {

	double dx = tree->x - pulled->x();
	double dy = tree->y - pulled->y();
	double r2 = dx * dx + dy * dy;
	double d2 = tree->width * tree->width;
	double theta2 = theta * theta;
	if (d2 < theta2 * r2) {
		pulled->accGravityFrom(Body(tree->x, tree->y, 0.0, 0.0, tree->m));
	} else {
		if (tree->sub != NULL) { //tree has children

			if (tree->child0 != NULL)
				BH(pulled, tree->child0, theta);

			if (tree->child1 != NULL)
				BH(pulled, tree->child1, theta);

			if (tree->child2 != NULL)
				BH(pulled, tree->child2, theta);

			if (tree->child3 != NULL)
				BH(pulled, tree->child3, theta);


		} else { //tree do not have children

			pulled->accGravityFrom(Body(tree->x, tree->y, 0.0, 0.0, tree->m));

		}
	}

}

void ForceBarnesHut::operator()(Body *pulled) {

	//Replace the naive code below with an efficient implementation
	//of Barnes-Hut.
	BH(pulled, mytree, theta_);

	//for (int i = 0; i < N_; i++)
	//		pulled->accGravityFrom(bodies_[i]);

}
;
