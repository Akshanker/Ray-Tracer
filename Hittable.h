#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "Vector.h"
#include "Ray.h"
#include<math.h>
struct Hit
{
	float t;
	Vector Ph;
	Vector surnorm;
	Vector revnorm;
	Vector color;
	float spec;
	float kr;
	float ior;

};
class Hittable
{ 

public:

	//Hittable();
	virtual bool intersect(Ray r, float& t) { return 0; }
	virtual Vector getNormal(Vector hitpos) { return Vector(0, 0, 0); }
	virtual Vector getCol() { return Vector(0, 0, 0); }
	virtual float getSpec() { return 0; }
	virtual float getRF() { return 0; }
	virtual float getIOR() { return 0; }
	float getRadius() { return 0; }
	Vector getVert() { return Vector(0, 12, 0); }
	virtual Vector getTex(Ray ray,Vector hitpos) { return Vector(0, 0, 0); }
	virtual std::string getType() { return 0; }
	
};
