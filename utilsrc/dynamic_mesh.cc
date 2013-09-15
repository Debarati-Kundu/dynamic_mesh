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

#define WRITE_TO_FILE
#define SUBTRACT_MEAN
#define DEBUG

int main(int argc, char *argv[])
{
	int i, j, k;
	TriMesh::set_verbose(0);
	int mFrames = atoi(argv[2]);
	int mVertices = 0;
	char buffer[200];

	TriMesh *MeshArr[mFrames];
	
	for(i = 0; i < mFrames; i++)
	{
		sprintf (buffer, "%s0%03d.off", argv[1], i+1);
		MeshArr[i] = TriMesh::read(buffer);

		// Code to check the mesh is deformable, the edge length is varying between the frames (it is a deformable mesh)
		MeshArr[i]->need_faces();
		float dist_face1 = dist(MeshArr[i]->vertices[MeshArr[i]->faces[0][1]], MeshArr[i]->vertices[MeshArr[i]->faces[0][2]]);
		float dist_face2 = dist(MeshArr[i]->vertices[MeshArr[i]->faces[0][2]], MeshArr[i]->vertices[MeshArr[i]->faces[0][0]]);
		float dist_face3 = dist(MeshArr[i]->vertices[MeshArr[i]->faces[0][0]], MeshArr[i]->vertices[MeshArr[i]->faces[0][1]]);

	}
	mVertices = MeshArr[0]->vertices.size();
	vector< vector<float> > trajectories;
	trajectories.resize(mVertices);

	vector<float> mean_vert; // Mean of a vertex accross multiple frames
        mean_vert.assign(3*mFrames, 0);

        for (i = 0; i < mVertices; i++)
        {
              	trajectories[i].resize(mFrames*3);
                for (j = 0; j < mFrames; j++)
                {
		     	mean_vert[j*3+0] += MeshArr[j]->vertices[i][0]/mVertices;
			mean_vert[j*3+1] += MeshArr[j]->vertices[i][1]/mVertices;
			mean_vert[j*3+2] += MeshArr[j]->vertices[i][2]/mVertices;
                }

		for(j = 0; j < mFrames; j++)
		{
			trajectories[i][3*j + 0] = MeshArr[j]->vertices[i][0];
                        trajectories[i][3*j + 1] = MeshArr[j]->vertices[i][1];
                        trajectories[i][3*j + 2] = MeshArr[j]->vertices[i][2];
		}
        }

#ifdef SUBTRACT_MEAN
	for (j = 0; j < mFrames; j++)
	{
		for(i = 0; i < mVertices; i++)
		{
                        trajectories[i][3*j + 0] = MeshArr[j]->vertices[i][0] - mean_vert[j*3+0];
                        trajectories[i][3*j + 1] = MeshArr[j]->vertices[i][1] - mean_vert[j*3+1];
                        trajectories[i][3*j + 2] = MeshArr[j]->vertices[i][2] - mean_vert[j*3+2];
		}
	}
#endif

#ifdef DEBUG
	for (i = 0; i < mFrames; i++)
		cout << "Frame number " << i << endl << "Mean of X coord " << mean_vert[3*i] << endl << "Mean of Y coord " << mean_vert[3*i+1] << endl <<  "Mean of Z coord " << mean_vert[3*i+2] << endl;
#endif

	// Written to file in the format suitable for SVD
#ifdef WRITE_TO_FILE
	FILE *ftraj = fopen("utilsrc/trajectory_matrix.txt", "w");
	for (i = 0; i < mVertices; i++)
	{
		for (j = 0; j < mFrames; j++)
		{
			fprintf(ftraj, "%f ", trajectories[i][3*j + 0]);
			fprintf(ftraj, "%f ", trajectories[i][3*j + 1]);
			fprintf(ftraj, "%f ", trajectories[i][3*j + 2]);
		}
		fprintf(ftraj, "\n");
	}
	fclose(ftraj);
#endif	
	return 0;
}
	
