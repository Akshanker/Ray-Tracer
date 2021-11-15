#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "Vector.h"
#include "Ray.h"
#include "Hittable.h"
#include<math.h>

class Material
{
public:

	unsigned char* texmap;
	int width;
	int height;
	std::string type;
	float spec;
	float kr;
	float ior;

	Material(){}
	Material(unsigned char* _texmap, int _width, int _height,std::string _type,float _spec,float _kr,float _ior)
	{
		texmap = _texmap;
		width = _width;
		height = _height;
		type = _type;
		spec = _spec;
		kr = _kr;
		ior = _ior;
	}



};