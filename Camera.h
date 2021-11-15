#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "Vector.h"
#include "Ray.h"

#include<math.h>

using namespace std;

class Camera
{
public:
	Vector pe;   //camera position
	Vector v_view;   //horizontal orientation from cam
	Vector v_up;	//verical orientation from cam
	Vector n0, n1, n2;    //unit vectors for screen? n2 is v_view normalized
	Vector pc;    //screen centre position
	Vector p00;  //lower left corner of screen
	//float focl;
	float sx, sy;  //aspect ratio
	//float Xsize, Ysize; //screen size;
	float d;  //camera to screen distance

	float lensX, lensY, focus;
	
	Camera() {}
	Camera(Vector campos, Vector lookat, Vector vertical, float hasp, float vasp, float dist)
	{
		pe = campos;
		v_view = lookat;
		v_up = vertical;
		sx = hasp;
		sy = vasp;
		d = dist;
		n2 = v_view.unitv();
		n0 = (cross(v_view, v_up)).unitv();
		n1 = cross(n0, n2);

		pc = pe + d * n2;
		p00 = pc -  (sx * n0 + sy * n1) / 2;


	}
	
	//Vector origin; () { Vector(0.0f, 0.0f, 0.0f); }
	//float aspectw;
	//float aspecth;
	//Vector origin;

	/*Camera(float aspecth,float aspectw)
	{
		aspecth = aspecth;
		aspectw = aspectw;
		origin = Vector(0.0f, 0.0f, 0.0f);
	}*/
	Ray getRay(float u, float v)
	{
		Vector pp = (p00 + (u * n0 * sx) + (v * n1 * sy));
		pp = (pp - pe).unitv();
		return Ray(pe, pp);

	}

};



