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
#include "Triangle.h"
#include "Obj.h"

using namespace std;
//#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include "extlib\OBJ_Loader.h"
//#include "tiny_obj_loader.h"

vector<Hittable*> Obj_load(string filename, vector<Hittable*>& objects2,  Material *mat)
{
	objl::Loader Loader;
	//vector<Hittable*>objects2;
	// Load .obj File
	bool loadout = Loader.LoadFile(filename);

	bool check = loadout;

	// Check to see if it loaded
	Triangle* tri;
	// If so continue
	if (loadout)
	{

		// Go through each loaded mesh and out its contents
		for (int i = 0; i < Loader.LoadedMeshes.size(); i++)
		{
			objl::Mesh curMesh = Loader.LoadedMeshes[i];
			for (int j = 0; j < curMesh.Vertices.size() - 3; j += 3)
			{
				Vector T_P0, T_P1, T_P2;

				T_P0 = Vector(curMesh.Vertices[j].Position.X, curMesh.Vertices[j].Position.Y, curMesh.Vertices[j].Position.Z);
				
				T_P1 = Vector(curMesh.Vertices[j + 1].Position.X, curMesh.Vertices[j + 1].Position.Y, curMesh.Vertices[j + 1].Position.Z);
				
				T_P2 = Vector(curMesh.Vertices[j + 2].Position.X, curMesh.Vertices[j + 2].Position.Y, curMesh.Vertices[j + 2].Position.Z);
				
				tri = new Triangle(T_P0, T_P1, T_P2, mat);

				objects2.push_back(tri);

			}
		}

	}

	return objects2;
}