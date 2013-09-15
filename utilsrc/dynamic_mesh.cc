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

#define MAX_NEIGHBOR 20 // For the time being assume that's the max we can handle.
//#define WRITE_TO_FILE_TRAJ
#define WRITE_TO_FILE_COHERENCE
//#define SUBTRACT_MEAN
//#define DEBUG

float dot_product(vector<float> a, vector<float> b, int n)
{
	float result = 0;
	int ind;
	for(ind = 0; ind < n; ind++)
		result += a[ind]*b[ind];
	return result;
}

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
	for (j = 0; j < mFrames; j++)
		cout << "Frame number " << j << endl << "Mean of X coord " << mean_vert[3*j] << endl << "Mean of Y coord " << mean_vert[3*j+1] << endl <<  "Mean of Z coord " << mean_vert[3*j+2] << endl;
#endif

	// Written to file in the format suitable for SVD
#ifdef WRITE_TO_FILE_TRAJ
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

	// Finding out motion coherency
	vector <vector <float> > motion_vec;
	motion_vec.resize(mVertices);

	for(i = 0; i < mVertices; i++)
		motion_vec[i].resize(3*mFrames);		

	for (j = 0; j < mFrames-1; j++)
	{
                MeshArr[j]->need_neighbors();
		for (i = 0; i < mVertices; i++)
		{
			motion_vec[i][3*j] = MeshArr[j+1]->vertices[i][0] - MeshArr[j]->vertices[i][0];
			motion_vec[i][3*j + 1] = MeshArr[j+1]->vertices[i][1] - MeshArr[j]->vertices[i][1];
			motion_vec[i][3*j + 2] = MeshArr[j+1]->vertices[i][2] - MeshArr[j]->vertices[i][2];
		}
	}

#ifdef WRITE_TO_FILE_COHERENCE
	FILE *fcor = fopen("utilsrc/coherence_matrix.txt", "w");
#endif
	for (j = 0; j < mFrames-1; j++)
	{
		vector< vector <float> > motion_vec_arr;
		motion_vec_arr.resize(mVertices);

		float coherence_matrix[9]; 
		for(i = 0; i < mVertices; i++)
		{
			int num_neighbors = MeshArr[j]->neighbors[i].size();
			vector<float>motion_x;
			vector<float>motion_y;
			vector<float>motion_z;
			motion_x.assign(MAX_NEIGHBOR, 0);
			motion_y.assign(MAX_NEIGHBOR, 0);
			motion_z.assign(MAX_NEIGHBOR, 0);			
			
			for(k = 0; k < num_neighbors; k++)
			{
				motion_x[k] = motion_vec[MeshArr[j]->neighbors[i][k]][3*j];
				motion_y[k] = motion_vec[MeshArr[j]->neighbors[i][k]][3*j + 1];
				motion_z[k] = motion_vec[MeshArr[j]->neighbors[i][k]][3*j + 2];
			} 
			coherence_matrix[0] = dot_product(motion_x, motion_x, num_neighbors);
			coherence_matrix[1] = dot_product(motion_x, motion_y, num_neighbors);
			coherence_matrix[2] = dot_product(motion_x, motion_z, num_neighbors);
		
			coherence_matrix[3] = dot_product(motion_x, motion_y, num_neighbors);
                        coherence_matrix[4] = dot_product(motion_y, motion_y, num_neighbors);
                        coherence_matrix[5] = dot_product(motion_y, motion_z, num_neighbors);

			coherence_matrix[6] = dot_product(motion_x, motion_z, num_neighbors);
                        coherence_matrix[7] = dot_product(motion_z, motion_y, num_neighbors);
                        coherence_matrix[8] = dot_product(motion_z, motion_z, num_neighbors);

#ifdef WRITE_TO_FILE_COHERENCE
			for(k = 0; k < 9; k++)
				fprintf(fcor, "%f ", coherence_matrix[k]);
			fprintf(fcor, "\n");
#endif			
		}
	}
#ifdef WRITE_TO_FILE_COHERENCE
        fclose(fcor);
#endif
	return 0;
}
	
