#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "Vector.h"
#include "Ray.h"
#include "Hittable.h"
#include<math.h>
#include "Material.h"
#include "Hit_list.h"

class Triangle :public Hittable
{
public:
	Vector s1;
	Vector s2;
	Vector s3;
	Material* mat;
	//Vector color;
	//float spec;

	Triangle(){}
	Triangle(Vector _s1, Vector _s2, Vector _s3,Material *_mat)
	{
		s1 = _s1;
		s2 = _s2;
		s3 = _s3;
		mat = _mat;
		//color = _color;
		//spec = _spec;
	}
	bool intersect(Ray r, float& t);
	Vector getNormal(Vector hitpos);
	Material getMaterial(Vector hitpos);
	Vector getTex(Ray ray, Vector hitpos);
	Vector getCol() { return Vector(0,8,0); }
	Vector getVert();
	std::string getType();
	float getSpec();
	float getRF();
	float getIOR();
	//float getSpec() { return spec; }

};
Vector Triangle :: getVert() { return s3; }
bool Triangle::intersect(Ray r, float& t)
{
	Vector p1 = s2 - s1;
	Vector p2 = s3 - s1;

	Vector N = cross(p1, p2);// .unitv()) / 2;
	float A = N.length()/2;

	float parll = dot(-N, r.dir);

	float t_min = 0.001, tmax = 1000;

	if (parll > 0.001) {

		Vector p0l0 = r.origin - s1;
		t = -dot(p0l0, -N) / parll;
		if (t<tmax && t>t_min) {

			Vector P = r.origin + t * r.dir;

			Vector C;

			Vector bc = s3 - s2;
			Vector bp = P - s2;
			C = cross(bc, bp);
			float u = (C.length()) / 2;
			u /= A;
			if (dot(N, C) < 0)
			{
				return false;
			}
			Vector ca = s1 - s3;
			Vector cp = P - s3;
			C = cross(ca, cp);
			float v = (C.length()) / 2;
			v /= A;
			if (dot(N, C) < 0) {
				return false;
			}
			float w = 1 - u - v;

			if (w > 0 && w < 1) {
				return true;
			}

		}



	}

	return false;

	/*if (fabs(parll) < 0.000000001)
	{
		return false;
	}
	float d = dot(N, s1-r.origin);

	t = (dot(N, r.origin) + d) / parll;

	if (t < 0)
	{
		return false;
	}
	Vector P = r.origin + t * r.dir;
	Vector C;
	Vector e1 = s2 - s1;
	Vector s1P = P - s1;
	C = cross(e1, s1P);
	if (dot(C, N) < 0)return false;

	Vector e2 = s3 - s2;
	Vector s2P = P - s2;
	C = cross(e2, s2P);
	if (dot(N, C) < 0) return false;

	Vector e3 = s1 - s3;
	Vector s3P = P - s3;
	C = cross(e3, s3P);
	if (dot(N, C) < 0) return false;*/

	//return true;
}
Material Triangle :: getMaterial(Vector hitpos)
{
	return *mat;
}
Vector Triangle::getTex(Ray ray, Vector hitpos)
{
	int channels2, texwid, texht;
	//stbi_set_flip_vertically_on_load(true);
	unsigned char* texture = mat->texmap;// stbi_load("earthmap.jpg", &texwid, &texht, &channels2, 0);
	texwid = mat->width;
	texht = mat->height;

	Vector p1 = s2 - s1;
	Vector p2 = s3 - s2;

	Vector N = cross(p1, p2);
	Vector A = N/ 2;
	float parll = dot(N, ray.dir);
	float d = dot(N, s1 - ray.origin);

	float t = (dot(N, ray.origin) + d) / parll;

	Vector P = ray.origin + t * ray.dir;

	Vector A1 = cross((P - s1), (s3 - P)) / 2;
	Vector A2 = cross((P - s2), (s1 - P)) / 2;

	float u = dot(N, A1) / A.length();
	float v = dot(N, A2) / A.length();

	if (u < 0)
	{
		u = u - floor(u);
		u = 1 - u;
	}
	else
	{
		u = u - floor(u);
	}
	if (v < 0)
	{
		v = v - floor(v);
		v = 1 - v;
	}
	else
	{
		v = v - floor(v);
	}

	float Xtex = u * texwid;
	float Ytex = v * texht;
	int tp = ((int)Ytex * texwid + (int)Xtex) * 3;
	//if(tp>=(texwid*texht*3))
	//cout <<Xtex<<" "<<Ytex<<" "<<tp << endl;
	float r = texture[tp] / 255.0;
	float g = texture[tp + 1] / 255.0;
	float b = texture[tp + 2] / 255.0;

	return Vector(r, g, b);
	//return Vector(0.4, 0.4, 0.7);
}
Vector Triangle :: getNormal(Vector hitpos)
{
	Vector p1 = s2 - s1;
	Vector p2 = s3 - s1;

	Vector N = cross(p1, p2);// .unitv()) / 2;
	
	return N.unitv();
}
std::string Triangle::getType()
{
	return mat->type;
}
float Triangle::getSpec()
{
	return mat->spec;
}
float Triangle::getRF()
{
	return mat->kr;
}
float Triangle::getIOR()
{
	return mat->ior;
}