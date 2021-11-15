#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "Vector.h"

using namespace std;

class Ray
{
public:
	Vector origin, dir;

	//constructors

	Ray() { };
	Ray(Vector o, Vector d)
	{
		origin = o;
		dir = d;
	}
	Vector getOrigin() { return origin; }
	Vector getDir() { return dir; }

	Vector pos(float k) { return origin + k * dir; }

};