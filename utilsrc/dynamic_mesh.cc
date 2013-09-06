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
	TriMesh::set_verbose(0);
	int mPoints = atoi(argv[2]);
	char buffer[200];

	cout << "Number of vertices " << mPoints << endl;
	TriMesh *MeshArr[mPoints];
//	n=sprintf (buffer, "%s00%0d", argv[1],);
//	MeshArr[0] = TriMesh::read(argv[1]);
//	MeshArr[0]->need_faces();
//	printf("size of faces %d\n", (int) MeshArr[0]->faces.size());
//	TriMesh *mesh = TriMesh::read(argv[1]);
	return 0;
}

