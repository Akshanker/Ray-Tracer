#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "Vector.h"
#include "Ray.h"
#include "Hittable.h"
#include<math.h>
#include "Material.h"

class Plane: public Hittable
{
public:
	Vector pos;
	Vector normal;
	Material* mat;
	//Vector color;
	//float spec;

	Plane(){}
	Plane(Vector p, Vector n,Material *_mat)
	{
		pos = p;
		normal = n;
		mat = _mat;
		//color = _color;
		//spec = _spec;


	}
	bool intersect(Ray r, float& t);
	Vector getNormal(Vector hitpos);
	Material getMaterial(Vector hitpos);
	Vector getTex(Ray ray, Vector hitpos);
	std::string getType();
	float getSpec();
	float getRF();
	float getIOR();
	//Vector getCol() { return color; }
	//float getSpec() { return spec; }
};
bool Plane::intersect(Ray r, float& t)
{
	normal.normalize();
	Vector npe = r.dir;
	npe = npe.unitv();
	float denom = dot(normal, r.dir);
	if (abs(denom) > 0.00000001)
	{
		Vector dist = pos-r.origin;
		float num = dot(dist, normal);
		//if (num < 0 && denom < 0)
		//{
			t = num / denom;
			if (t >= 0)
			{
				return true;
			}
		//}

		
	}
	return false;
}
Vector Plane :: getTex(Ray ray,Vector hitpos)
{
	int channels2, texwid, texht;
	//stbi_set_flip_vertically_on_load(true);
	unsigned char* texture = mat->texmap;// stbi_load("earthmap.jpg", &texwid, &texht, &channels2, 0);
	texwid = mat->width;
	texht = mat->height;
	float sx = 10, sy = 10;
	//Vector pos
	//Vector pp = (Vector(0,0,25) + (u * n0 * sx) + (v * n1 * sy));

	Vector a = cross(normal, Vector(1, 0, 0));
	Vector btmp = cross(normal, Vector(0, 1, 0));

	Vector maxAB = dot(a, a) < dot(btmp, btmp) ? btmp:a;
	Vector c = cross(normal, Vector(0, 0, 1));

	Vector maxVec = dot(maxAB, maxAB) < dot(c, c) ? c : maxAB;
	maxVec = maxVec.unitv();
	Vector U = maxVec;
	Vector V = cross(normal, U);


	Vector n0 = normal;
	Vector n1, n2;
	n0 = n0.unitv();

	Vector tmp;

	if (abs(n0.y) == 1)
	{
		n1 = Vector(1, 0, 0);
		n1 = n1.unitv();
		n2 = Vector(0, 0, 1);

		n2 = n2.unitv();//Vector(0, 1, 0);//
	}
	else if (abs(n0.x) == 1)
	{
		n1 = Vector(0, 1, 0);
		n1 = n1.unitv();
		n2 = Vector(0, 0, 1);

		n2 = n2.unitv();//Vector(0, 1, 0);//
	}
	else if (abs(n0.z) == 1)
	{
		n1 = Vector(0, 1, 0);
		n1 = n1.unitv();
		n2 = Vector(1, 0, 0);

		n2 = n2.unitv();//Vector(0, 1, 0);//

	}
	Vector pe = Vector(0, -2, 1.2);
	Vector lookat = Vector(0, 0, -1);

	



	float phi = dot(n1, (hitpos - pos));
	float theta = dot(n2, (hitpos - pos));
	float u = phi / sx;//dot(U, hitpos-pos);//
	float v = theta / sy; //dot(V, hitpos-pos);//

	//u = floor((hitpos-pos).x);
	//v = floor((hitpos-pos).z);
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

	float Xtex = u * texwid;//u * texwid;//
	float Ytex = v * texht;// (1 - v)* texht - 0.001;//

	float xx = Xtex - floor(Xtex);
	float yy = Ytex - floor(Ytex);

	int x0 = Xtex;// int(floor(Xtex + 0.5));
	int y0 = Ytex;// int(floor(Ytex + 0.5));
	/*if ((Xtex + Ytex)%2 == 0)
	{
		return Vector(1, 1, 1);
	}
	else
	{
		return Vector(0, 0, 0);
	}*/
	if (x0 < 0)x0 = -x0;
	if (x0 < 0)y0 = -y0;
	if (x0 > texwid)x0 = x0 % (texwid);
	//cout << v << endl;
	if (y0 > texht)y0 = y0 % (texht);

	int tp = (y0 * texwid + x0) * 3;
	//if(tp>=(texwid*texht*3))
	//cout <<Xtex<<" "<<Ytex<<" "<<tp << endl;




	float r = texture[tp] / 255.0;
	float g = texture[tp + 1] / 255.0;
	float b = texture[tp + 2] / 255.0;

	return Vector(r, g, b);
	//return Vector(0.2, 0.5, 0.4);

}
Vector Plane :: getNormal(Vector hitpos)
{


	return -(normal.unitv());
}
std::string Plane::getType()
{
	return mat->type;
}
float Plane::getSpec()
{
	return mat->spec;
}
float Plane::getRF()
{
	return mat->kr;
}
float Plane::getIOR()
{
	return mat->ior;
}