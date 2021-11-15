#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "Vector.h"
#include "Ray.h"
#include "Hittable.h"
#include<math.h>
#include "Material.h"

class Hit_list : Hittable
{
public:

	Hittable** list;
	int lsize;
	Hit_list(){}
	Hit_list(Hittable** _list, int _lsize)
	{
		list = _list;
		lsize = _lsize;
	}
	bool intersect(Ray r, float& t);
	Vector getNormal(Vector hitpos);


};
bool Hit_list::intersect(Ray r, float& t)
{
	for (int i = 0; i < lsize; i++)
	{
		if (list[i]->intersect(r, t))
		{
			return true;
		}
	}
	return false;
}
Vector Hit_list :: getNormal(Vector hitpos)
{
	for (int i = 0; i < lsize; i++)
	{
		return list[i]->getNormal(hitpos);
	}
}
