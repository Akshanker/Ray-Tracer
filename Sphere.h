#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "Vector.h"
#include "Ray.h"
#include "Hittable.h"
#include<math.h>
#include "Material.h"
#include <exception>


using namespace std;

class Sphere:public Hittable
{
public:
	Vector centre;
	float radius;
	Material *mat;
	//Vector color;
	//float spec;

	Sphere() {}
	Sphere(Vector _centre, float _radius, Material *_mat)
	{
		centre = _centre;// Vector(0, 0, 0);
		radius = _radius;// 0.5;
		mat=_mat;
		//color = _color;// Vector(0.4, 0.6, 0.3);
		//spec = _spec;
	}
	Vector getCenter()
	{
		return centre;
	}
	bool intersect(Ray r, float& t);
	Vector getNormal(Vector hitpos);
	Material getMaterial(Vector hitpos);
	
	//Vector getCol() { return color; }
	//float getSpec() { return spec; }
	float getRadius() { return radius; }
	Vector getTex(Ray ray,Vector hitpos);
	std::string getType();
	float getSpec();
	float getRF();
	float getIOR();
};
bool Sphere::intersect(Ray r,float& t)
{
	float t0, t1;
	Vector p_dist = r.origin - centre;
	Vector npe = r.dir;
	npe = npe.unitv();

	float a = dot(npe, npe);
	float b = dot(p_dist, npe);
	float c = dot(p_dist, p_dist) - radius * radius;

	float disc = b * b - c;

	if (disc > 0 && b <= 0)
	{
		t0 = (-b + sqrt(disc));
		t1 = (-b - sqrt(disc));

		t = t1;// min(t0, t1);
		return true;
	}
	

	return false;
}
Vector Sphere :: getNormal(Vector hitpos)
{
	Vector sn = (hitpos-centre).unitv();
	return sn;
}

/*Material Sphere :: getMaterial(Vector hitpos)
{
	return *mat;
}*/
Vector Sphere :: getTex(Ray ray,Vector hitpos)
{
	int channels2,texwid,texht;
	//stbi_set_flip_vertically_on_load(true);
	unsigned char* texture = mat->texmap;// stbi_load("earthmap.jpg", &texwid, &texht, &channels2, 0);
	texwid = mat->width;
	texht = mat->height;
	//texture = new unsigned char[height2 * width2 * 3];
	Vector cen = centre;// Vector(0, 0, 20);
	Vector tmpp = (hitpos-cen) / radius;

	//original formula

	//float phi = atan2((tmpp.z), tmpp.x);// (dot(Vector(0, 1, 0), tmpp);// acos(snormal.z);// atan2(hitpos.z, hitpos.x);
	//float theta = asin(tmpp.y);// (dot(Vector(1, 0, 0), tmpp);// acos(snormal.y / sin(phi));// asin(hitpos.y);
	//float v = (theta+M_PI/2)/M_PI;// acos(phi) / M_PI;// (theta + M_PI / 2) / M_PI;
	///float u = 1 - (phi + M_PI) / (2 * M_PI);// acos(theta / sin(M_PI * v)) / 2 * M_PI;
	
	//ergun formula
	float phi= (dot(Vector(0, 1, 0), tmpp));
	float theta = (dot(Vector(1, 0, 0), tmpp));
	float v= acos(phi) / M_PI;
	//float u= acos(theta / sin(M_PI * v)) / 2 * M_PI;
	float u = 1 - (theta + M_PI) / (2 * M_PI);
	float Xtex = u * texwid;
	float Ytex = v * texht;// (1 - v)* texht - 0.001;
	
	int tp =  ((int)Ytex * texwid + (int)Xtex) * 3;// ((int)Ytex * width2 + (int)Xtex) * 3;
	//if(tp>=157941)cout << tp<< endl;
	//cout << acos(phi) << endl;
	//try {
		float r = texture[tp] / 255.0;
		float g = texture[tp + 1] / 255.0;
		float b = texture[tp + 2] / 255.0;
		
	//}
	//catch(exception& e)
	//{
	//	cout << "exception" << tp << endl;
	//}
	return Vector(r, g, b);
	//return Vector(0.4, 0.5, 0.9);

}
std::string Sphere :: getType()
{
	return mat->type;
}
float Sphere::getSpec()
{
	return mat->spec;
}
float Sphere::getRF()
{
	return mat->kr;
}
float Sphere::getIOR()
{
	return mat->ior;
}

