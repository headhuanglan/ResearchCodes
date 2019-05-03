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
#ifndef VECTOR3_H
#define VECTOR3_H
#include <cmath>
class vector3
{
public:
  vector3(): x(0), y(0), z(0)
    {}

  vector3(double _x, double _y, double _z): x(_x), y(_y), z(_z)
    {}

  inline double magnitude() const
    {
        return sqrt(x * x + y * y + z * z);
    }

  inline void normalize()
    {
        double magnitude = this->magnitude();
        if (magnitude < 1e-7)
        {
            *this = vector3(0, 0, 0);
        }
        else
        {   
			 magnitude=1.0/magnitude;
            *this =vector3(x*magnitude,y*magnitude,z*magnitude);
        }
    }

   inline vector3 normalized() const
    {
        vector3 v(*this);
        v.normalize();
        return v;
    }


  inline double distance(const vector3& other) const
    {
        return (*this - other).magnitude();
    }

  inline double dot(const vector3& other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

   inline vector3 operator+(const vector3& other) const
    {
        return vector3(x + other.x, y + other.y, z + other.z);
    }

    inline vector3 operator-(const vector3& other) const
    {
        return vector3(x - other.x, y - other.y, z - other.z);
    }

    inline vector3 operator*(double s) const
    {
        return vector3(x * s, y * s, z * s);
    }

    vector3 operator()(double a, double b, double c) 
      {   
		  this->x=a;
		  this->y=b;
		  this->z=c;

		  return *this;
      }



    inline static double dot(const vector3& a, const vector3& b)
    {
        return a.dot(b);
    }

    inline static double angle(const vector3& a, const vector3& b)
    {
        return acos(dot(a.normalized(), b.normalized())) * 57.29578;  //57.29578 rad2deg
	}

	

public:
    double x;
    double y;
	double z;

};

#endif //VECTOR3_H
