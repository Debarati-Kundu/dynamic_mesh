/*
Szymon Rusinkiewicz
Princeton University

Benedict Brown
Katholieke Universiteit Leuven

Debarati Kundu
University of Texas at Austin

dynamic_mesh.cc
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <vector>
using namespace std;


int main(int argc, char *argv[])
{
	// Don't clutter the output
	int i;
	TriMesh::set_verbose(0);
	int mFrames = atoi(argv[2]);
	char buffer[200];

	TriMesh *MeshArr[mFrames];
	
	for(i = 0; i < mFrames; i++)
	{
		sprintf (buffer, "%s0%03d.off", argv[1], i+1);
		MeshArr[i] = TriMesh::read(buffer);

		// Code to check the mesh is deformable, the edge length is varying between the frames	
		MeshArr[i]->need_faces();
		float dist_face1 = dist(MeshArr[i]->vertices[MeshArr[i]->faces[0][1]], MeshArr[i]->vertices[MeshArr[i]->faces[0][2]]);
		float dist_face2 = dist(MeshArr[i]->vertices[MeshArr[i]->faces[0][2]], MeshArr[i]->vertices[MeshArr[i]->faces[0][0]]);
		float dist_face3 = dist(MeshArr[i]->vertices[MeshArr[i]->faces[0][0]], MeshArr[i]->vertices[MeshArr[i]->faces[0][1]]);
		// printf("%d %f %f %f\n", i, dist_face1, dist_face2, dist_face3);
	}
	
	return 0;
}

