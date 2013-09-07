/*
Szymon Rusinkiewicz
Princeton University

Benedict Brown
Katholieke Universiteit Leuven

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
	}
	
//	i = 50;
//	MeshArr[i]->need_faces();
//	printf("done %d\n", MeshArr[i]->faces.size());
//	TriMesh *mesh = TriMesh::read(argv[1]);
	return 0;
}

