#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include<cmath>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"

#include "Vector.h"
#include "Ray.h"
#include "Camera.h"
#include "Camera2.h"
#include "Sphere.h"
#include"HolSphere.h"
#include "Triangle.h"
#include "Plane.h"
#include "Hittable.h"
#include "Obj.h"
#include "Obj_load.h"
#include "Hit_list.h"


using namespace std;
int width1, height1, channels1,width2, height2, channels2,width3,height3,channels3,width4,height4,channels4, width5, height5, channels5,width7,height7,channels7,width8,height8,channels8;
int width6, height6, channels6;
unsigned char* pixmap;


//int height = 360, width = 400;
int width = 1280;
int height = 720;

float radians(float x)
{
	return x * M_PI / 180;
}
float clamp(float x)
{
	if (x > 255)
	{
		x = 255;
	}
	if (x < 0)
	{
		x = 0;
	}
	return x;
}

float clamp01(float x)
{
	if (x > 1)
	{
		x = 1;
	}
	if (x < 0)
	{
		x = 0;
	}
	return x;
}

float clampp(float min, float max, float x)
{
	if (x < min)x = min;
	if (x > max)x = max;

	return x;
}
Vector dirL = Vector(2, 3.6, 0.3);
//Vector dirL = Vector(-2, 2,-5);
//Vector dirL = Vector(-5, 5, -2);
//Vector dirL = Vector(5, 4, -1);
Vector L = Vector(0,0,1);
// Vector(sin(30 * M_PI / 180), cos(30 * M_PI / 180), -1).unitv();
vector<Hittable*>objects;

vector<Hittable*>objects2;


//Vector col;
unsigned char* texmap6 = stbi_load("sky6.jpg", &width6, &height6, &channels6, 0);
Material* bgmat = new Material(texmap6, width6, height6,"lambert",0,0,0);
Vector getBG(float u, float v, Vector tmpray,Material* bgmat)
{
	int channels2, texwid, texht;
	unsigned char* texture = bgmat->texmap;// stbi_load("earthmap.jpg", &texwid, &texht, &channels2, 0);
	texwid = bgmat->width;
	texht = bgmat->height;
	
	float Xtex = u * texwid;
	float Ytex = v * texht;
	//cout << texwid << endl;
	//cout << u << endl;
	int tp = ((int)Ytex * texwid + (int)Xtex) * 3;// ((int)Ytex * width2 + (int)Xtex) * 3;
	//if(tp>=157941)cout << tp<< endl;
	//cout << acos(phi) << endl;
	//try {
	
	float r = texture[tp] / 255.0;
	float g = texture[tp + 1] / 255.0;
	float b = texture[tp + 2] / 255.0;

	return Vector(r, g, b);
}
/*Vector refract(Vector dir, Vector normal, float ior)
{
	Vector n;
	float etai = 1, etat = ior;

	// clampp(-1, 1, dot(dir, normal));
	
	if (dot(normal, dir) > 0)
	{
		n = -normal;
		etai = ior;
		etat = 1;
	}
	else
	{
		n = normal;

	}
	
	float cosi = dot(normal, dir);

	float eta = etai/etat;

	float k = 1 - eta * eta * (1 - cosi * cosi);

	if (k < 0)k = 0;

	return  eta * (dir - cosi * n) - n * sqrt(k);// (eta * cosi - sqrtf(k))* n;


}*/
Vector refract(Vector& I, Vector& N, const float& ior)
{
	float cosi = clamp(-1, 1, dot(I, N));
	float etai = 1, etat = ior;
	Vector n = N;
	if (cosi < 0) { cosi = -cosi; }
	else { std::swap(etai, etat); n = -N; }
	float eta = etai / etat;
	float k = 1 - eta * eta * (1 - cosi * cosi);
	//Vector refr;
	//if (k < 0) { refr=Vector(0; }
	return k < 0 ? Vector(0,0,0) : (eta * I + (eta * cosi - sqrtf(k)) * n);
}



void fresnel( Vector& I,  Vector& N, const float& ior, float& kr)
{
	float cosi = clamp(-1, 1, dot(I, N));
	float etai = 1, etat = ior;
	if (cosi > 0) { std::swap(etai, etat); }
	// Compute sini using Snell's law
	float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
	// Total internal reflection
	if (sint >= 1) {
		kr = 1;
	}
	else {
		float cost = sqrtf(std::max(0.f, 1 - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
		float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
		kr = (Rs * Rs + Rp * Rp) / 2;
	}
	// As a consequence of the conservation of energy, transmittance is given by:
	// kt = 1 - kr;
}


Vector getColor(Ray r, Hit rec)
{
	
	Vector npe = r.dir;
	
	//Vector ln = (dirL - rec.Ph).unitv(); // wrong - but works
	//Vector ln = (rec.Ph - dirL).unitv(); // right
	Vector ln = L.unitv();

	float sf = (dot(ln, rec.surnorm) + 1) / 2;//max(0,(int)(min((int)dot(ln, sn),1)));
	//if (sf < 0) { sf = 0; }
	float kd = 0.9, ka = 0.1;
	//cout << dot(ln, rec.surnorm) << endl;
	Vector v = (r.origin - rec.Ph).unitv();
	Vector rv = -ln + 2 * (dot(ln, rec.surnorm)) * rec.surnorm;
	Vector hv = (v + rv).unitv();
	float sp = dot(hv, v);

	float diff = (ka + (kd * sf));
	//L = L - dirL;
	//L = L.unitv();
	//if (dot(L, ln) - cos((50 * M_PI) / 180) <= 0)
	//{
	//	diff = 0.3;// (ka + (kd * sf));
	//}
	return ((diff)*rec.color + (1.0 - diff) * Vector(0, 0, 0)) +(rec.spec * pow(sp, 300) * rec.color * Vector(1, 1, 1));// 0.5 * Vector(sn.x + 1, sn.y + 1, sn.z + 1);

}

Vector trace(Ray r,float depth)
{
	Hit rec; 
	Vector pixcol = Vector(1, 0, 0);
	Vector shadowcol;
	//Vector reflcol;
	//Vector bounce = Vector(0, 0, 0);
	Vector bgcol;
	if (depth > 20)return bgcol;
	bool hit = false;
	float t = 0; float tmin = INFINITY, min_index = -1, vis = 1;
	float t1 = 0;
	Vector tmpray = r.dir;
	tmpray.normalize();
	float phi = atan2(-tmpray.z, tmpray.x);
	float theta = asin(tmpray.y);
	float u = 1 - (phi + M_PI) / (2 * M_PI);
	float v = (theta + M_PI / 2) / M_PI;
	//cout << u << endl;
	if (u > 1)
	{
		u = 1;
	}
	if (v > 1)
	{
		v = 1;
	}
	//if (u < 0)
	//{
	//	u = 0;
	//}
	//if (v < 0)
	//{
	//	v = 0;
	//}

	Vector envcol = getBG(u, v, tmpray, bgmat);
	pixcol =  envcol;

	for (int i = 0; i < objects.size(); i++)
	{
		if (objects.at(i)->intersect(r, t))
		{
			if (t < tmin)
			{
				tmin = t;
				min_index = i;
			}
		}
	}
	//if(min_index==2)cout << min_index << endl;
	if (min_index != -1)
	{
		
		
		rec.t = tmin;
		
		rec.Ph = r.origin + (tmin * r.dir);
		rec.surnorm = (objects.at(min_index)->getNormal(rec.Ph));

		//cout << rec.surnorm << endl;
		//rec.revnorm = objects.at(min_index)->getrevNormal(rec.Ph);
		int mode=0;
		if (objects.at(min_index)->getType() == "refl")
			mode = 1;
		if (objects.at(min_index)->getType() == "refr")
			mode = 2;
		if (objects.at(min_index)->getType() == "lambert")
			mode = 6;
		rec.color = objects.at(min_index)->getTex(r, rec.Ph);// textcol(rec.Ph);
		rec.spec = objects.at(min_index)->getSpec();
		rec.kr = objects.at(min_index)->getRF();
		rec.ior = objects.at(min_index)->getIOR();
		//pixcol = getColor(r, rec);
		//float lamb = 0;
		//cout << depth << endl;
		bool outside = dot(r.dir, rec.surnorm) < 0;
		//Reflection ///
		
		
		
		//Vector refrcol = Vector(0, 0, 0);
		//Vector reflcol = Vector(0, 0, 0);
		//if (dot(r.dir,rec.surnorm) > 0) rec.surnorm = -rec.surnorm, inside = true;
		//Vector I= ((rec.Ph - dirL).unitv());
		Vector oppDir = -1.0 * r.dir;
		Vector rrd = -oppDir + 2.0 * (dot(oppDir, rec.surnorm)) * rec.surnorm;
		rrd = rrd.unitv();

		//cout << rrd << endl;
		
		Vector bias = 1.0* rec.surnorm;
		//switch (mode) {
		//case 1:
		if(mode==1)
		    {
			    pixcol = getColor(r, rec);

				Vector srd = L.unitv();// 1.0 * ((rec.Ph - dirL)).unitv();
				//Vector srd= (rec.Ph - dirL).unitv();
				Vector Ph_prime = rec.Ph + 0.0005 * srd;
				Ray shadowray = Ray(Ph_prime, srd);
				for (int j = 0; j < objects.size() && j != min_index; j++)
				{
					if (objects.at(j)->intersect(shadowray, t1))
					{

						pixcol = Vector(0.3, 0.3, 0.3) * pixcol;// 
					}
				}

				bool checkout = dot(rrd, rec.surnorm) > 0;
				Vector newPh = Vector(rec.Ph.x, rec.Ph.y, rec.Ph.z);
				Vector reflp = checkout ? newPh + bias : newPh - bias;//rrd = outside ? rrd : -rrd;// 
				//cout << newPh << endl;
				Ray reflray = Ray(reflp, rrd);
				//cout << reflray.origin << endl;
				//bounce=
				Vector reflcol = trace(reflray, depth + 1);
				float kr = rec.kr;
				//fresnel(r.dir, rec.surnorm, 1.6, kr);
				pixcol = (kr)*reflcol + (1 - kr) * pixcol; // trace(reflray, depth + 1, bounce);
				//return pixcol;
				
			}
			/// /refraction/////


			
		//case 2:
		else if(mode==2)
			{
				//pixcol = getColor(r, rec);
				bool checkout = dot(rrd, rec.surnorm) > 0;
				
				Vector newPh = Vector(rec.Ph.x, rec.Ph.y, rec.Ph.z);
				Vector reflp = checkout ? newPh - bias : newPh + bias;//rrd = outside ? rrd : -rrd;// 
				//cout << newPh << endl;
				Ray reflray = Ray(reflp, rrd);
				//cout << reflray.origin << endl;
				//bounce=
				float kr = rec.kr;
				Vector reflcol = trace(reflray, depth + 1);

				
				Vector tmpn,temph, R;
				float ior = 2.417;
				//Vector refrrd = refract(r.dir, rec.surnorm, 1.45);
				float eta1 = 1.0, eta2 = ior;
				float coseta = dot(rec.surnorm, r.dir);
				if (coseta > 0)
				{
					tmpn = -1.0*rec.surnorm;
					eta1 = ior;
					eta2 = 1.0;
				}
				else
				{
					eta1 = 1.0;
					eta2 = ior;
					coseta = -coseta;
					
				}
				float eta = eta1 / eta2;
				float coseta2 = dot(tmpn, r.dir);
				
				float k = 1.0 - eta * eta * (1.0 - coseta * coseta);
				if (k < 0.0)
				{
					k = 0;// R = rrd;
				}
				k = sqrt(k);
				
				R =  eta* r.dir + (eta * coseta - k) * rec.surnorm;// eta* (r.dir - coseta2 * tmpn) - tmpn * sqrt(k);//
				Vector refrrd = R;
				refrrd = refrrd.unitv();
				//cout << R << endl;
				bool checkoutref = dot(refrrd, rec.surnorm) > 0;
				Vector biasr = 1.4 * refrrd;
				Vector refrp = checkoutref ? rec.Ph - biasr : rec.Ph + biasr;
				Ray refrRay = Ray(refrp, refrrd);


				Vector refrcol = trace(refrRay, depth + 1);
				
				//kr = (pow((eta1 * coseta - eta2 * k) / (eta1 * coseta + eta2 * k), 2.0) + pow((eta2 * coseta - eta1 * k) / (eta1 * k + eta2 * coseta), 2.0)) * 0.5;
				//fresnel(r.dir, rec.surnorm, ior, kr);
				//cout << kr << endl;
				pixcol = ((kr) * reflcol + (1-kr) * refrcol);
				//return pixcol;
				//break;
			}
			//if (!outside) cout << bounce << endl;
			//cout << bounce << endl;
			//pixcol = 0.5*bounce;

			///arealights/

			//for (int j = 0; j < 15; j++)
			//{
				//for (int i = 0; i <15; i++)
				//{

					//float u = (float)i / (float)2;
					//float v = (float)j / (float)2;
					//Vector areal = (dirL + (i * Vector(1,0,0)) + (j * Vector(0,1,0)));
		
			//Vector srd = ((rec.Ph-areal)).unitv();
		//default:
		else
		{
			pixcol = getColor(r, rec);
			//float lamb = 0;
			//for (int j = 0; j < 25; j++)
			//{
				//for (int i = 0; i <25; i++)
				//{

					//float u = (float)i / (float)2;
					//float v = (float)j / (float)2;
					//Vector areal = (dirL + (i * Vector(1,0,0)) + (j * Vector(0,1,0)));
					Vector srd = L.unitv(); // 1.0 * ((rec.Ph - dirL)).unitv();
					//Vector srd= (rec.Ph - dirL).unitv();
					Vector Ph_prime = rec.Ph + 0.0005 * srd;
					Ray shadowray = Ray(Ph_prime, srd);
					for (int j = 0; j < objects.size() && j != min_index; j++)
					{
						if (objects.at(j)->intersect(shadowray, t1))
						{

							pixcol = Vector(0.3, 0.3, 0.3) * pixcol;// 
						}
					}

			//}
		//}
		//lamb /= 225;

		//pixcol = pixcol * (1 - lamb) + Vector(0, 0, 0) * lamb;
			
			//break;
		}
		}
		return pixcol;
	}
	/*else
	{
		Vector tmpray = r.dir;
		tmpray.normalize();
		float phi = atan2(-tmpray.z, tmpray.x);
		float theta = asin(tmpray.y);
		float u = 1 - (phi + M_PI) / (2 * M_PI);
		float v = (theta + M_PI / 2) / M_PI;
		//cout << u << endl;
		if (u > 1)
		{
			u = 1;
		}
		if (v > 1)
		{
			v = 1;
		}
		//if (u < 0)
		//{
		//	u = 0;
		//}
		//if (v < 0)
		//{
		//	v = 0;
		//}

		pixcol =  getBG(u, v, tmpray, bgmat);
		//return bgcol;

	}*/
	
	
	//return pixcol;


void main()
{
	stbi_set_flip_vertically_on_load(true);
	unsigned char* texmap1 = stbi_load("blue.jpg", &width1, &height1, &channels1, 0);
	unsigned char* texmap2 = stbi_load("steel.jpg", &width2, &height2, &channels2, 0);
	unsigned char* texmap3 = stbi_load("carp2.jpg", &width3, &height3, &channels3, 0);
	unsigned char* texmap4 = stbi_load("spek.jpg", &width4, &height4, &channels4, 0);
	unsigned char* texmap5 = stbi_load("jupiter.jpg", &width5, &height5, &channels5, 0);
	unsigned char* texmap7 = stbi_load("fire.jpg", &width7, &height7, &channels7, 0);
	unsigned char* texmap8 = stbi_load("caust.jpg", &width8, &height8, &channels8, 0);

	Material *mat1 = new Material(texmap1, width1, height1,"refr",0.7,0.3,1.1);
	Material* mat2 = new Material(texmap2, width2, height2, "refl",0.7,0.8,0);
	Material* mat3 = new Material(texmap3, width3, height3, "lambert",0,0.0,0);
	Material* mat4 = new Material(texmap4, width4, height4, "refl",0.5,0.6,2.417);
	Material* mat5 = new Material(texmap5, width5, height5, "lambert",0.4,0.7,0);
	Material* mat7 = new Material(texmap7, width7, height7, "refl", 1.2, 0.3, 0);
	Material* mat8 = new Material(texmap8, width8, height8, "lambert", 0, 0.3, 0);



	//for (int frames = 11; frames < 12; frames++)
	//{
		Vector ocent = Vector(0, 0, -2.5);//
		// = mrad.unitv();
		//float mtheta = 0 * M_PI / 180;// acos(dot(Vector(1, 0, 0), mrad));
		float frames = 100;
		float orbit = 2.5;
		float inc = frames;
		float x1 = ocent.x + cos(radians(52 + inc)) * orbit;
		float z1 = ocent.z + sin(radians(52 + inc)) * orbit;

		float x2 = ocent.x + cos(radians(142 + inc)) * orbit;
		float z2 = ocent.z + sin(radians(142 + inc)) * orbit;

		float x3 = ocent.x + cos(radians(232 + inc)) * orbit;
		float z3 = ocent.z + sin(radians(232 + inc)) * orbit;

		float x4 = ocent.x + cos(radians(322 + inc)) * orbit;
		float z4 = ocent.z + sin(radians(322 + inc)) * orbit;

		cout << x1 << endl;
		cout << z1 << endl;
		cout << x2 << endl;
		cout << z2 << endl;
		cout << x3 << endl;
		cout << z3 << endl;
		cout << x4 << endl;
		cout << z4 << endl;


		Plane* pln = new Plane(Vector(0, -0.4, 5), Vector(0, 1, 0), mat3);
		Sphere* sphr = new Sphere(Vector((float)x1, -1, (float)z1), 0.6, mat1);
		Sphere* sphr1 = new Sphere(Vector((float)x3, -1, (float)z3), 0.6, mat2);
		Sphere* sphr2 = new Sphere(Vector((float)x2, -1, (float)z2), 0.6, mat5);
		//Sphere* sphr3 = new Sphere(Vector(x2, -4, z2), 0.8, mat8);
		Sphere* sphr4 = new Sphere(Vector((float)x4, -1, (float)z4), 0.6, mat7);
		//HolSphere* sphr5 = new HolSphere(Vector(0, 0, 20), 150, mat5);
		//objects2.push_back(sphr4);

		//Obj_load("cyll.obj",objects, mat1);

		//Obj_load("torus.obj", objects, mat3);







		/*Plane* pln1 = new Plane(Vector(0, -0.4, -5), Vector(0, 0, 1), mat4);
		Plane* pln2 = new Plane(Vector(-3, -0.4, -1), Vector(-1, 0, 0), mat4);
		Plane* pln3 = new Plane(Vector(3, -0.4, -1), Vector(1, 0, 0), mat4);*/
		//Triangle* tri = new Triangle(Vector(0, 1, 14), Vector(-2, 1, 14), Vector(-2, -1, 14),Vector(1,0,0),0.5);

		/*Triangle* tri1 = new Triangle(Vector(1.319, 3.584, 22.661), Vector(-2.891, 3.584, 19.320), Vector(0.123, 1.265, 15.522),mat1);
		Triangle* tri2 = new Triangle(Vector(0.123, 1.265, 15.522), Vector(4.333, 1.265, 18.863), Vector(1.319, 3.584, 22.661), mat1);
		Triangle* tri3 = new Triangle(Vector(4.333, 1.265, 18.863), Vector(0.123, 1.265, 15.522), Vector(-1.319, -3.584, 17.339), mat1);
		Triangle* tri4 = new Triangle(Vector(-1.319, -3.584, 17.339), Vector(2.891, -3.584, 20.680), Vector(4.333, 1.265, 18.863), mat1);
		Triangle* tri5 = new Triangle(Vector(0.123, 1.265, 15.522), Vector(-1.319, -3.584, 17.339), Vector(-4.333, -1.265, 21.137), mat1);
		Triangle* tri6 = new Triangle(Vector(-4.333, -1.265, 21.137), Vector(-2.891, 3.584, 19.320), Vector(0.123, 1.265, 15.522), mat1);

		Triangle* tri7 = new Triangle(Vector(-0.923, 0.923, -1.923), Vector(-0.923, -0.923, 1.923), Vector(-0.923, -0.923, -0.077), mat1);
		Triangle* tri8 = new Triangle(Vector(6.802, -0.821, 23.427), Vector(7.423, 0.200, 16.9454), Vector(13.363, -0.821, 24.056), mat1);
		Triangle* tri9 = new Triangle(Vector(13.363, -0.821, 24.056), Vector(6.705, 5.691, 24.444), Vector(7.423, 0.200, 16.9454), mat1);*/
		//Obj* mesh = loadPolyMeshFromFile("cubehandle.obj");
		//objects.push_back(pln1);
		//objects.push_back(pln2);
		//objects.push_back(pln3);
		//cout << mesh.nfaces << endl;
		/*objects.push_back(tri1);
		objects.push_back(tri2);
		objects.push_back(tri3);
		objects.push_back(tri4);
		objects.push_back(tri5);
		objects.push_back(tri6);*/
		//objects.push_back(tri7);
		//objects.push_back(tri8);
		//objects.push_back(tri9);
		//objects.push_back(pln);
		//objects.push_back(tri);


		//Obj_load("dmnd2.obj", objects, mat4);
		objects.push_back(sphr);
		objects.push_back(sphr1);
		objects.push_back(sphr2);
		objects.push_back(sphr4);
		//objects.push_back(sphr3);

		objects.push_back(pln);





		//image array//////

		//Camera controls///
		Vector pe = Vector(0, -3, 1.2);
		Vector lookat = Vector(0, 0, -1);
		///image////////
		float aspect = (float(width)) / (float(height));
		float sy = 3;
		float sx = sy * aspect;
		float fov = 90;

		float imageaspratio = width / height;
		pixmap = new unsigned char[height * width * 3];
		float depth = 0;
		///camera orientation////////
		Vector v_view = lookat;
		Vector v_up = Vector(0, 1, 0);
		float d = 3.3;
		float aperture = 0.1;
		float focl = 1.0;
		//Camera cam(pe,v_view,v_up,sx,sy,d);
		//Camra cam2(pe, lookat, v_up, 20, imageaspratio, aperture, d);
		camera2 cam2 = camera2(pe, lookat, v_up, 80, 1.778, 0, 1);
		//cout << cam.p00 << endl;
		//objects////
		Vector color;
		//Sphere sp(Vector(0,0,0),0.5,Vector(0.1,0.6,0.4));
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				int rps = (j * width + i) * 3;
				int gps = rps + 1;
				int bps = rps + 2;
				int M = 10, N = 10;

				/*for (int m = 0; m < M; m++)
				{
					for (int n = 0; n < N; n++)
					{
						float rx = (float)rand() / (float)RAND_MAX;
						float ry = (float)rand() / (float)RAND_MAX;

						float X = i + (n + rx) / N;
						float Y = j + (m + ry) / M;*/

						float s = (float)i / (float)width; //(float)i / width;
						float t = (float)j / (float)height;// (float)j / height - 1;
						Ray r = cam2.get_ray(s, t);
						//cout << r.dir << endl;

						color = trace(r, depth);

						//cout << color<< endl;

						pixmap[rps] = clamp(abs(255 * color.x));
						pixmap[gps] = clamp(abs(255 * color.y));
						pixmap[bps] = clamp(abs(255 * color.z));
					//}

				//}

			}
		}



		/*testing
		Vector col;
		col=getcol(col);
		writeimg(col);
		cout << col.x << endl;



		Vector testy(2, 2, 1);
		Vector t2(0, 8, 6);
		cout << testy.length() << endl;
		//testy.normalize();
		cout << testy.y << endl;

		cout << dot(testy, t2) << endl;
		cout << cross(testy, t2) << endl;

		testy += t2;
		cout << testy << endl;

		testy *= 2;
		cout << testy << endl;

		cout << -testy << endl;
		/*testing*//////////////////
		int framen = 0;
		stringstream st;
		st << "frame" << framen << ".jpg";
		string filename;
		st >> filename;



		//stbi_write_jpg( width1, height1, 3, pixmap3, 100);
		stbi_write_jpg(filename.c_str(), width, height, 3, pixmap, 100);
		objects.clear();

	//}
	cout << "twar";


}

